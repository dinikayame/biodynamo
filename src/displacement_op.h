#ifndef DISPLACEMENT_OP_H_
#define DISPLACEMENT_OP_H_

#include <array>
#include <cmath>
#include <vector>
#include "grid.h"
#include "math_util.h"
#include "param.h"

namespace bdm {

using std::array;

template <typename TSO>
void ApplyBoundingBox(TSO* cell, double lb, double rb) {
  // Need to create a small distance from the positive edge of each dimension;
  // otherwise it will fall out of the boundary of the simulation space
  double eps = 1e-10;
  auto& pos = cell->GetPosition();
  for (int i = 0; i < 3; i++) {
    if (pos[i] < lb) {
      cell->SetCoordinate(i, lb);
    }
    if (pos[i] >= rb) {
      cell->SetCoordinate(i, rb - eps);
    }
  }
}

/// Defines the 3D physical interactions between physical objects
template <typename TGrid = Grid<>>
class DisplacementOp {
 public:
  DisplacementOp() {}
  ~DisplacementOp() {}

  template <typename TContainer>
  void operator()(TContainer* cells, uint16_t type_idx) const {
    std::vector<array<double, 3>> cell_movements;
    cell_movements.reserve(cells->size());

    auto& grid = TGrid::GetInstance();
    auto search_radius = grid.GetLargestObjectSize();
    double squared_radius = search_radius * search_radius;

#pragma omp parallel for shared(grid) firstprivate(squared_radius)
    for (size_t i = 0; i < cells->size(); i++) {
      auto&& cell = (*cells)[i];
      // Basically, the idea is to make the sum of all the forces acting
      // on the Point mass. It is stored in translationForceOnPointMass.
      // There is also a computation of the torque (only applied
      // by the daughter neurites), stored in rotationForce.

      // TODO(roman) : There might be a problem, in the sense that the biology
      // is not applied if the total Force is smaller than adherence.
      // Once, I should look at this more carefully.

      // If we detect enough forces to make us  move, we will re-schedule us
      // setOnTheSchedulerListForPhysicalObjects(false);

      // fixme why? copying
      const auto& tf = cell.GetTractorForce();

      // the 3 types of movement that can occur
      // bool biological_translation = false;
      bool physical_translation = false;
      // bool physical_rotation = false;

      double h = Param::simulation_time_step_;
      std::array<double, 3> movement_at_next_step{0, 0, 0};

      // BIOLOGY :
      // 0) Start with tractor force : What the biology defined as active
      // movement------------
      movement_at_next_step[0] += h * tf[0];
      movement_at_next_step[1] += h * tf[1];
      movement_at_next_step[2] += h * tf[2];

      // PHYSICS
      // the physics force to move the point mass
      std::array<double, 3> translation_force_on_point_mass{0, 0, 0};
      // the physics force to rotate the cell
      // std::array<double, 3> rotation_force { 0, 0, 0 };

      // 1) "artificial force" to maintain the sphere in the ecm simulation
      // boundaries--------
      // 2) Spring force from my neurites (translation and
      // rotation)--------------------------
      // 3) Object avoidance force
      // -----------------------------------------------------------
      //  (We check for every neighbor object if they touch us, i.e. push us
      //  away)

      auto calculate_neighbor_forces = [&](auto&& neighbor,
                                           auto&& neighbor_handle) {
        std::array<double, 3> neighbor_force;
        neighbor.GetForceOn(cell.GetMassLocation(), cell.GetDiameter(),
                            &neighbor_force);
        translation_force_on_point_mass[0] += neighbor_force[0];
        translation_force_on_point_mass[1] += neighbor_force[1];
        translation_force_on_point_mass[2] += neighbor_force[2];
      };

      grid.ForEachNeighborWithinRadius(calculate_neighbor_forces, cell,
                                       SoHandle(type_idx, i), squared_radius);

      // 4) PhysicalBonds
      // How the physics influences the next displacement
      double norm_of_force = std::sqrt(translation_force_on_point_mass[0] *
                                           translation_force_on_point_mass[0] +
                                       translation_force_on_point_mass[1] *
                                           translation_force_on_point_mass[1] +
                                       translation_force_on_point_mass[2] *
                                           translation_force_on_point_mass[2]);

      // is there enough force to :
      //  - make us biologically move (Tractor) :
      //  - break adherence and make us translate ?
      physical_translation = norm_of_force > cell.GetAdherence();

      assert(cell.GetMass() != 0 && "The mass of a cell was found to be zero!");
      double mh = h / cell.GetMass();
      // adding the physics translation (scale by weight) if important enough
      if (physical_translation) {
        // We scale the move with mass and time step
        movement_at_next_step[0] += translation_force_on_point_mass[0] * mh;
        movement_at_next_step[1] += translation_force_on_point_mass[1] * mh;
        movement_at_next_step[2] += translation_force_on_point_mass[2] * mh;

        // Performing the translation itself :

        // but we want to avoid huge jumps in the simulation, so there are
        // maximum distances possible
        if (norm_of_force * mh > Param::simulation_max_displacement_) {
          const auto& norm = Math::Normalize(movement_at_next_step);
          movement_at_next_step[0] =
              norm[0] * Param::simulation_max_displacement_;
          movement_at_next_step[1] =
              norm[1] * Param::simulation_max_displacement_;
          movement_at_next_step[2] =
              norm[2] * Param::simulation_max_displacement_;
        }
      }
      cell_movements[i] = movement_at_next_step;
    }

// Set new positions after all updates have been calculated
// otherwise some cells would see neighbors with already updated positions
// which would lead to inconsistencies
#pragma omp parallel for
    for (size_t i = 0; i < cells->size(); i++) {
      auto&& cell = (*cells)[i];
      cell.UpdateMassLocation(cell_movements[i]);
      if (Param::bound_space_) {
        ApplyBoundingBox(&cell, Param::min_bound_, Param::max_bound_);
        grid.SetDimensionThresholds(Param::min_bound_, Param::max_bound_);
      }
      cell.SetPosition(cell.GetMassLocation());

      // Reset biological movement to 0.
      cell.SetTractorForce({0, 0, 0});
    }
  }
};

/// Keeps the simulation objects contained within the bounds as defined in
/// param.h
class BoundSpace {
 public:
  BoundSpace() {}
  ~BoundSpace() {}

  template <typename TContainer>
  void operator()(TContainer* cells, uint16_t type_idx) const {
    // set new positions after all updates have been calculated
    // otherwise some cells would see neighbors with already updated positions
    // which would lead to inconsistencies
    auto& grid = Grid<>::GetInstance();
#pragma omp parallel for
    for (size_t i = 0; i < cells->size(); i++) {
      auto&& cell = (*cells)[i];
      if (Param::bound_space_) {
        ApplyBoundingBox(&cell, Param::min_bound_, Param::max_bound_);
        grid.SetDimensionThresholds(Param::min_bound_, Param::max_bound_);
      }
      cell.SetPosition(cell.GetMassLocation());

      // Reset biological movement to 0.
      cell.SetTractorForce({0, 0, 0});
    }
  }
};

}  // namespace bdm

#endif  // DISPLACEMENT_OP_H_
