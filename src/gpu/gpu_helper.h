#ifndef GPU_GPU_HELPER_H_
#define GPU_GPU_HELPER_H_

#include <vector>
#include <iostream>

#include "gpu/opencl_kernels.h"
#include "resource_manager.h"

namespace bdm {

template <typename TResourceManager = ResourceManager<>>
static void FindGpuDevices() {
  try {
    // We keep the context and device list in the resource manager to be
    // accessible elsewhere to create command queues and buffers from
    auto rm = TResourceManager::Get();
    cl::Context* context = rm->GetOpenCLContext();
    std::vector<cl::Device>* devices = rm->GetOpenCLDeviceList();

    // Get list of OpenCL platforms.
    std::vector<cl::Platform> platform;
    cl::Platform::get(&platform);

    if (platform.empty()) {
      std::cerr << "OpenCL platforms not found." << std::endl;
    }

    // Lambda to select only the first available GPU device
    auto first_gpu = [&]() {
      return devices->empty();
    };

    // Go over all available platforms and devices until first device is found
    for (auto p = platform.begin(); first_gpu() && p != platform.end();
         p++) {
      std::vector<cl::Device> pldev;

      try {
        p->getDevices(CL_DEVICE_TYPE_GPU, &pldev);

        for (auto d = pldev.begin(); first_gpu() && d != pldev.end(); d++) {
          if (!d->getInfo<CL_DEVICE_AVAILABLE>())
            continue;

          // The OpenCL extension available on this device
          std::string ext = d->getInfo<CL_DEVICE_EXTENSIONS>();

          devices->push_back(*d);
        }
      } catch (...) {
        std::cout << "Found bad platform! Continuing to next one." << std::endl;
        devices->clear();
        continue;
      }
    }

    if (devices->empty()) {
      std::cerr << "No GPU found!" << std::endl;
    }

    *context = cl::Context(*devices);

    std::cout << std::endl << "Found the following GPU: " << (*devices)[0].getInfo<CL_DEVICE_NAME>() << std::endl << std::endl;

    // Create command queue.
    cl::CommandQueue queue(*context, (*devices)[0], CL_QUEUE_PROFILING_ENABLE);
  } catch (const cl::Error &err) {
    std::cerr << "OpenCL error: " << err.what() << "(" << err.err() << ")"
              << std::endl;
  }
}

template <typename TResourceManager = ResourceManager<>>
static void CompileKernels() {
  auto rm = TResourceManager::Get();
  std::vector<cl::Program>* all_programs = rm->GetOpenCLProgramList();
  cl::Context* context = rm->GetOpenCLContext();
  std::vector<cl::Device>* devices = rm->GetOpenCLDeviceList();
  // Compile OpenCL program for found device
  // TODO(ahmad): create more convenient way to compile all OpenCL kernels, by
  // going through a list of header files. Also, create a stringifier that goes
  // from .cl --> .h, since OpenCL kernels must be input as a string here
  cl::Program displacement_op_program(
      *context,
      cl::Program::Sources(1, std::make_pair(displacement_op_kernel, strlen(displacement_op_kernel))));

  all_programs->push_back(displacement_op_program);

  std::cout << "Compiling OpenCL kernel" << std::endl;
  for (auto& prog : *all_programs) {
    try {
        prog.build(*devices);
    } catch (const cl::Error &) {
      std::cerr << "OpenCL compilation error" << std::endl
                << prog.getBuildInfo<CL_PROGRAM_BUILD_LOG>((*devices)[0])
                << std::endl;
    }
  }
}

static const char *getErrorString(cl_int error)
{
switch(error){
    // run-time and JIT compiler errors
    case 0: return "CL_SUCCESS";
    case -1: return "CL_DEVICE_NOT_FOUND";
    case -2: return "CL_DEVICE_NOT_AVAILABLE";
    case -3: return "CL_COMPILER_NOT_AVAILABLE";
    case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5: return "CL_OUT_OF_RESOURCES";
    case -6: return "CL_OUT_OF_HOST_MEMORY";
    case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8: return "CL_MEM_COPY_OVERLAP";
    case -9: return "CL_IMAGE_FORMAT_MISMATCH";
    case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11: return "CL_BUILD_PROGRAM_FAILURE";
    case -12: return "CL_MAP_FAILURE";
    case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15: return "CL_COMPILE_PROGRAM_FAILURE";
    case -16: return "CL_LINKER_NOT_AVAILABLE";
    case -17: return "CL_LINK_PROGRAM_FAILURE";
    case -18: return "CL_DEVICE_PARTITION_FAILED";
    case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
    case -30: return "CL_INVALID_VALUE";
    case -31: return "CL_INVALID_DEVICE_TYPE";
    case -32: return "CL_INVALID_PLATFORM";
    case -33: return "CL_INVALID_DEVICE";
    case -34: return "CL_INVALID_CONTEXT";
    case -35: return "CL_INVALID_QUEUE_PROPERTIES";
    case -36: return "CL_INVALID_COMMAND_QUEUE";
    case -37: return "CL_INVALID_HOST_PTR";
    case -38: return "CL_INVALID_MEM_OBJECT";
    case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40: return "CL_INVALID_IMAGE_SIZE";
    case -41: return "CL_INVALID_SAMPLER";
    case -42: return "CL_INVALID_BINARY";
    case -43: return "CL_INVALID_BUILD_OPTIONS";
    case -44: return "CL_INVALID_PROGRAM";
    case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46: return "CL_INVALID_KERNEL_NAME";
    case -47: return "CL_INVALID_KERNEL_DEFINITION";
    case -48: return "CL_INVALID_KERNEL";
    case -49: return "CL_INVALID_ARG_INDEX";
    case -50: return "CL_INVALID_ARG_VALUE";
    case -51: return "CL_INVALID_ARG_SIZE";
    case -52: return "CL_INVALID_KERNEL_ARGS";
    case -53: return "CL_INVALID_WORK_DIMENSION";
    case -54: return "CL_INVALID_WORK_GROUP_SIZE";
    case -55: return "CL_INVALID_WORK_ITEM_SIZE";
    case -56: return "CL_INVALID_GLOBAL_OFFSET";
    case -57: return "CL_INVALID_EVENT_WAIT_LIST";
    case -58: return "CL_INVALID_EVENT";
    case -59: return "CL_INVALID_OPERATION";
    case -60: return "CL_INVALID_GL_OBJECT";
    case -61: return "CL_INVALID_BUFFER_SIZE";
    case -62: return "CL_INVALID_MIP_LEVEL";
    case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64: return "CL_INVALID_PROPERTY";
    case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66: return "CL_INVALID_COMPILER_OPTIONS";
    case -67: return "CL_INVALID_LINKER_OPTIONS";
    case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
    case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
    case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
    case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
    case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
    case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
    case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default: return "Unknown OpenCL error";
    }
}

}  // namespace bdm

#endif  // GPU_GPU_HELPER_H_
