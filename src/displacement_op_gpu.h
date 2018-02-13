#ifndef DISPLACEMENT_OP_GPU_H_
#define DISPLACEMENT_OP_GPU_H_

#include "displacement_op.h"
#include "displacement_op_cuda.h"
#include "resource_manager.h"

#include <fstream>

namespace bdm {

using std::array;

/// Defines the 3D physical interactions between physical objects
template <typename TGrid = Grid<>, typename TResourceManager = ResourceManager<>>
class DisplacementOpGpu {
 public:
  DisplacementOpGpu() {}
  ~DisplacementOpGpu() {}

  template <typename TContainer>
  void operator()(TContainer* cells, uint16_t type_idx) const {
    auto& grid = TGrid::GetInstance();

    std::vector<std::array<double,3 >> cell_movements(cells->size());
    std::vector<double> mass(cells->size());
    std::vector<uint32_t> starts;
    std::vector<uint16_t> lengths;
    std::vector<uint32_t> successors(cells->size());
    uint32_t box_length;
    uint32_t N = cells->size();
    std::array<uint32_t, 3> num_boxes_axis;
    std::array<int32_t, 3> grid_dimensions;

    // We need to create a mass vector, because it is not stored by default in
    // a cell container
    cells->FillMassVector(&mass);
    grid.GetSuccessors(&successors);
    grid.GetGPUBoxData(&starts, &lengths);
    grid.GetGridData(&box_length, num_boxes_axis, grid_dimensions);

    std::cout << "Launching CUDA kernel" << std::endl;

    displacement_op_cuda(cells->GetPositionPtr(), cells->GetDiameterPtr(), cells->GetTractorForcePtr(), cells->GetAdherencePtr(), mass.data(), &(Param::simulation_time_step_), &(Param::simulation_max_displacement_), &N, starts.data(), lengths.data(), successors.data(), &box_length, num_boxes_axis.data(), grid_dimensions.data(), cell_movements.data()->data());

    // remove("gpu.txt");
    // std::ofstream ofs("gpu.txt", std::ofstream::out);
    // for (size_t k = 0; k < cell_movements.size(); k++) {
    //   ofs << cell_movements[k][0] << ", " << cell_movements[k][1] << ", " << cell_movements[k][2] << std::endl;
    // }
    // ofs.close();

// set new positions after all updates have been calculated
// otherwise some cells would see neighbors with already updated positions
// which would lead to inconsistencies
#pragma omp parallel for
    for (size_t i = 0; i < cells->size(); i++) {
      auto&& cell = (*cells)[i];
      cell.UpdateMassLocation(cell_movements[i]);
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

#endif  // DISPLACEMENT_OP_GPU_H_
