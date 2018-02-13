#include <math.h>
#include <stdint.h>

void displacement_op_cuda(double* positions, double* diameter, double* tractor_force, double* adherence, double* mass, double* timestep, double* max_displacement, uint32_t* N, uint32_t* starts, uint16_t* lengths, uint32_t* successors, uint32_t* box_length, uint32_t* num_boxes_axis, int32_t* grid_dimensions, double* cell_movements);
