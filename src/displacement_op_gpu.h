#ifndef DISPLACEMENT_OP_GPU_H_
#define DISPLACEMENT_OP_GPU_H_

#include "displacement_op.h"
#include "gpu/gpu_helper.h"
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
    auto rm = TResourceManager::Get();
    auto context = rm->GetOpenCLContext();
    auto devices = rm->GetOpenCLDeviceList();
    auto programs = rm->GetOpenCLProgramList();

    std::vector<cl_double> mass(cells->size());
    std::vector<array<cl_double, 3>> cell_movements(cells->size());
    std::vector<cl_uint> gpu_starts;
    std::vector<cl_ushort> gpu_lengths;
    std::vector<cl_uint> successors(cells->size());
    cl_uint box_length;
    std::array<cl_uint, 3> num_boxes_axis;
    std::array<cl_int, 3> grid_dimensions;

    // We need to create a mass vector, because it is not stored by default in
    // a cell container
    cells->FillMassVector(&mass);
    grid.GetSuccessors(&successors);
    grid.GetGPUBoxData(&gpu_starts, &gpu_lengths);
    grid.GetGridData(&box_length, num_boxes_axis, grid_dimensions);

    // Create command queue on the GPU device
    cl::CommandQueue queue(*context, (*devices)[0], CL_QUEUE_PROFILING_ENABLE);

    // Allocate GPU buffers
    cl::Buffer positions_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 cells->size() * 3 * sizeof(cl_double), cells->GetPositionPtr());
    cl::Buffer diameters_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 cells->size() * sizeof(cl_double), cells->GetDiameterPtr());
    cl::Buffer tractor_force_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 cells->size() * 3 * sizeof(cl_double), cells->GetTractorForcePtr());
    cl::Buffer adherence_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 cells->size() * sizeof(cl_double), cells->GetAdherencePtr());
    cl::Buffer mass_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 cells->size() * sizeof(cl_double), mass.data());
    cl::Buffer cell_movements_arg(*context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                 cells->size() * 3 * sizeof(cl_double), cell_movements.data()->data());
    cl::Buffer starts_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 gpu_starts.size() * sizeof(cl_uint), gpu_starts.data());
    cl::Buffer lengths_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 gpu_lengths.size() * sizeof(cl_short), gpu_lengths.data());
    cl::Buffer successors_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 successors.size() * sizeof(cl_uint), successors.data());
    cl::Buffer nba_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 3 * sizeof(cl_uint), num_boxes_axis.data());
    cl::Buffer gd_arg(*context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                 3 * sizeof(cl_int), grid_dimensions.data());

    // Create the kernel object from our program
    // TODO(ahmad): generalize the program selection, in case we have more than
    // one. We can maintain an unordered map of programs maybe
    cl::Kernel collide((*programs)[0], "collide");

    // Set kernel parameters
    collide.setArg(0, positions_arg);
    collide.setArg(1, diameters_arg);
    collide.setArg(2, tractor_force_arg);
    collide.setArg(3, adherence_arg);
    collide.setArg(4, mass_arg);
    collide.setArg(5, Param::simulation_time_step_);
    collide.setArg(6, Param::simulation_max_displacement_);

    collide.setArg(7, static_cast<cl_int>(cells->size()));
    collide.setArg(8, starts_arg);
    collide.setArg(9, lengths_arg);
    collide.setArg(10, successors_arg);
    collide.setArg(11, box_length);
    collide.setArg(12, nba_arg);
    collide.setArg(13, gd_arg);
    collide.setArg(14, cell_movements_arg);

    int block_size = 256;

    auto N = cells->size();

    try {
      // std::cout << "Global work size = " << (N + (block_size - (N%block_size))) << std::endl;
      queue.enqueueNDRangeKernel(collide, cl::NullRange, cl::NDRange(N + (block_size - (N%block_size))), cl::NDRange(block_size));
    } catch (const cl::Error &err) {
      std::cerr << "OpenCL error: " << err.what() << "(" << err.err() << ") = " << getErrorString(err.err()) 
                << std::endl;
      throw;
    }

    queue.enqueueReadBuffer(cell_movements_arg, CL_TRUE, 0, cells->size() * 3 * sizeof(cl_double), cell_movements.data()->data());

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
