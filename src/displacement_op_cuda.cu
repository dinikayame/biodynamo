#include "helper_math.h"
#include "displacement_op_cuda.h"
#include "stdio.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ double norm(double3 v) {
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

__device__ int3 get_box_coordinates(double3 pos, int32_t* grid_dimensions, uint32_t box_length) {
  int3 box_coords;
  box_coords.x = (floor(pos.x) - grid_dimensions[0]) / box_length;
  box_coords.y = (floor(pos.y) - grid_dimensions[1]) / box_length;
  box_coords.z = (floor(pos.z) - grid_dimensions[2]) / box_length;
  return box_coords;
}

__device__ uint32_t get_box_id_2(int3 bc, uint32_t* num_boxes_axis) {
  return bc.z * num_boxes_axis[0]*num_boxes_axis[1] + bc.y * num_boxes_axis[0] + bc.x;
}

__device__ uint32_t get_box_id(double3 pos, uint32_t* num_boxes_axis, int32_t* grid_dimensions, uint32_t box_length) {
  int3 box_coords = get_box_coordinates(pos, grid_dimensions, box_length);
  return get_box_id_2(box_coords, num_boxes_axis);
}

__device__ void compute_force(double* positions, double* diameters, uint32_t idx, uint32_t nidx, double3* result) {
  double r1 = 0.5 * diameters[idx];
  double r2 = 0.5 * diameters[nidx];
  // We take virtual bigger radii to have a distant interaction, to get a desired density.
  double additional_radius = 10.0 * 0.15;
  r1 += additional_radius;
  r2 += additional_radius;

  double comp1 = positions[3*idx + 0] - positions[3*nidx + 0];
  double comp2 = positions[3*idx + 1] - positions[3*nidx + 1];
  double comp3 = positions[3*idx + 2] - positions[3*nidx + 2];
  double center_distance = sqrt(comp1 * comp1 + comp2 * comp2 + comp3 * comp3);

  // the overlap distance (how much one penetrates in the other)
  double delta = r1 + r2 - center_distance;

  if (delta < 0) {
    return;
  }

  // to avoid a division by 0 if the centers are (almost) at the same location
  if (center_distance < 0.00000001) {
    result->x += 42.0;
    result->y += 42.0;
    result->z += 42.0;
    return;
  }

  // printf("Colliding cell [%d] and [%d]\n", idx, nidx);
  // printf("Delta for neighbor [%d] = %f\n", nidx, delta);

  // the force itself
  double r = (r1 * r2) / (r1 + r2);
  double gamma = 1; // attraction coeff
  double k = 2;     // repulsion coeff
  double f = k * delta - gamma * sqrt(r * delta);

  double module = f / center_distance;
  result->x += module * comp1;
  result->y += module * comp2;
  result->z += module * comp3;
}

__device__ void default_force(double* positions,
                   double* diameters,
                   uint32_t idx, uint32_t start, uint16_t length,
                   uint32_t* successors,
                   double3* result) {
  uint32_t nidx = start;
  for (uint16_t nb = 0; nb < length; nb++) {
    // implement logic for within radius here
    if (nidx != idx) {
      compute_force(positions, diameters, idx, nidx, result);
    }
    // traverse linked-list
    nidx = successors[nidx];
  }
}

__global__ void collide(
       double* positions,
       double* diameters,
       double* tractor_force,
       double* adherence,
       double* mass,
       double* timestep,
       double* max_displacement,
       uint32_t* N,
       uint32_t* starts,
       uint16_t* lengths,
       uint32_t* successors,
       uint32_t* box_length,
       uint32_t* num_boxes_axis,
       int32_t* grid_dimensions,
       double* result) {
  uint32_t tidx = blockIdx.x * blockDim.x + threadIdx.x;
  if (tidx < N[0]) {
		result[3*tidx + 0] += timestep[0] * tractor_force[3*tidx + 0];
    result[3*tidx + 1] += timestep[0] * tractor_force[3*tidx + 1];
    result[3*tidx + 2] += timestep[0] * tractor_force[3*tidx + 2];
    
    double3 pos;
    pos.x = positions[3*tidx + 0];
    pos.y = positions[3*tidx + 1];
    pos.z = positions[3*tidx + 2];

    double3 collision_force = make_double3(0, 0, 0);

    // Moore neighborhood
    int3 box_coords = get_box_coordinates(pos, grid_dimensions, box_length[0]);
    for (int z = -1; z <= 1; z++) {
      for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
          uint32_t bidx = get_box_id_2(box_coords + make_int3(x, y, z), num_boxes_axis);
          if (lengths[bidx] != 0) {
            default_force(positions, diameters, tidx, starts[bidx], lengths[bidx], successors, &collision_force);
          }
        }
      }
    }

    // Mass needs to non-zero!
    double mh = timestep[0] / mass[tidx];

    if (norm(collision_force) > adherence[tidx]) {
      result[3*tidx + 0] += collision_force.x * mh;
      result[3*tidx + 1] += collision_force.y * mh;
      result[3*tidx + 2] += collision_force.z * mh;

      if (norm(collision_force) * mh > max_displacement[0]) {
        result[3*tidx + 0] = max_displacement[0];
        result[3*tidx + 1] = max_displacement[0];
        result[3*tidx + 2] = max_displacement[0];
      }
    }
  }
}

void displacement_op_cuda(double* positions, double* diameters, double* tractor_force, double* adherence, double* mass, double* timestep, double* max_displacement, uint32_t* N, uint32_t* starts, uint16_t* lengths, uint32_t* successors, uint32_t* box_length, uint32_t* num_boxes_axis, int32_t* grid_dimensions, double* cell_movements) {
	double* d_positions = NULL;
    double* d_diameters = NULL;
    double* d_mass = NULL;
    double* d_timestep = NULL;
    double* d_max_displacement = NULL;
    uint32_t* d_N = NULL;
    double* d_cell_movements = NULL;
    double* d_tractor_force = NULL;
    double* d_adherence = NULL;
    uint32_t* d_starts = NULL;
    uint16_t* d_lengths = NULL;
    uint32_t* d_successors = NULL;
    uint32_t* d_box_length = NULL;
    uint32_t* d_num_boxes_axis = NULL;
    int32_t* d_grid_dimensions = NULL;

    uint32_t num_boxes = num_boxes_axis[0] * num_boxes_axis[1] * num_boxes_axis[2];

    cudaMalloc(&d_positions, 3 * N[0] * sizeof(double));
    cudaMalloc(&d_diameters, N[0] * sizeof(double));
    cudaMalloc(&d_tractor_force, 3 * N[0] * sizeof(double));
    cudaMalloc(&d_adherence, N[0] * sizeof(double));
    cudaMalloc(&d_mass, N[0] * sizeof(double));
    cudaMalloc(&d_timestep, sizeof(double));
    cudaMalloc(&d_max_displacement, sizeof(double));
    cudaMalloc(&d_N, sizeof(uint32_t));
    cudaMalloc(&d_starts, num_boxes * sizeof(uint32_t));
    cudaMalloc(&d_lengths, num_boxes * sizeof(uint16_t));
    cudaMalloc(&d_successors, N[0] * sizeof(uint32_t));
    cudaMalloc(&d_box_length, sizeof(uint32_t));
    cudaMalloc(&d_num_boxes_axis, 3 * sizeof(uint32_t));
    cudaMalloc(&d_grid_dimensions, 3 * sizeof(int32_t));
    cudaMalloc(&d_cell_movements, 3 * N[0] * sizeof(double));

    cudaMemcpy(d_positions, 		positions, 3 * N[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_diameters, 		diameters, N[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tractor_force, 	tractor_force, 3 * N[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_adherence, 		adherence, N[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_mass, 				mass, N[0] * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_timestep, 			timestep, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_max_displacement, 	max_displacement, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_N, 				N, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_starts, 			starts, num_boxes * sizeof(uint32_t), cudaMemcpyHostToDevice);
    gpuErrchk(cudaMemcpy(d_lengths, 			lengths, num_boxes * sizeof(uint16_t), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_successors, 		successors, N[0] * sizeof(uint32_t), cudaMemcpyHostToDevice));
    cudaMemcpy(d_box_length, 		box_length, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_num_boxes_axis, 	num_boxes_axis, 3 * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_grid_dimensions, 	grid_dimensions, 3 * sizeof(uint32_t), cudaMemcpyHostToDevice);

    int blockSize;
    int minGridSize;
    int gridSize;
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, collide, 0, N[0]);
    gridSize = (N[0] + blockSize - 1) / blockSize;

    printf("gridSize = %d  |  blockSize = %d\n", gridSize, blockSize);
    collide<<<gridSize, blockSize>>>(d_positions, d_diameters, d_tractor_force, d_adherence, d_mass, d_timestep, d_max_displacement, d_N, d_starts, d_lengths, d_successors, d_box_length, d_num_boxes_axis, d_grid_dimensions, d_cell_movements);

    cudaDeviceSynchronize();
    cudaMemcpy(cell_movements, d_cell_movements, 3 * N[0] * sizeof(double), cudaMemcpyDeviceToHost);
}
