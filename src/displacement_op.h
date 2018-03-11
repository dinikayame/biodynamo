#ifndef DISPLACEMENT_OP_H_
#define DISPLACEMENT_OP_H_

#include "displacement_op_cpu.h"
#include "displacement_op_cuda.h"
#include "displacement_op_opencl.h"
#include "grid.h"
#include "param.h"

namespace bdm {

using std::array;

template <typename TSO>
void ApplyBoundingBox(TSO* cell, double lb, double rb) {
  auto& pos = cell->GetPosition();
  for (int i = 0; i < 3; i++) {
    if (pos[i] < lb) {
      cell->SetCoordinate(i, lb);
    }
    if (pos[i] > rb) {
      cell->SetCoordinate(i, rb);
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
  void operator()(TContainer* cells, uint16_t type_idx) {
    if (Param::use_gpu_) {
      if (Param::use_opencl_) {
        opencl_(cells, type_idx);
      } else {
        cuda_(cells, type_idx);
      }
    } else {
      cpu_(cells, type_idx);
    }
  }

 private:
  DisplacementOpCpu<TGrid> cpu_;
  DisplacementOpCuda<TGrid> cuda_;
  DisplacementOpOpenCL<TGrid> opencl_;
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
#pragma omp parallel for
    for (size_t i = 0; i < cells->size(); i++) {
      auto&& cell = (*cells)[i];
      if (Param::bound_space_) {
        ApplyBoundingBox(&cell, Param::min_bound_, Param::max_bound_);
      }
      cell.SetPosition(cell.GetMassLocation());

      // Reset biological movement to 0.
      cell.SetTractorForce({0, 0, 0});
    }
  }
};

}  // namespace bdm

#endif  // DISPLACEMENT_OP_H_
