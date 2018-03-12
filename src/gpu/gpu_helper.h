#ifndef GPU_GPU_HELPER_H_
#define GPU_GPU_HELPER_H_

#include <vector>
#include <iostream>

#ifdef USE_OPENCL
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#endif
#ifdef USE_CUDA
#include "cuda_runtime_api.h"
#endif

#include "gpu/opencl_kernels.h"
#include "param.h"
#include "resource_manager.h"

namespace bdm {

#ifdef USE_OPENCL
static cl_int cl_assert(cl_int const code, char const * const file, int const line, bool const abort);
static const char *getErrorString(cl_int error);
#define cl_ok(err) cl_assert(err, __FILE__, __LINE__, true);
#endif

#ifdef USE_CUDA
static void FindGpuDevicesCuda() {
  int nDevices = 0;

  cudaGetDeviceCount(&nDevices);
  if (nDevices == 0) {
    std::cerr << "No CUDA-compatible GPU found on this machine! Switching to the CPU version..." << std::endl;
    Param::use_gpu_ = false;
    return;
  }

  std::cout << "Found " << nDevices << " CUDA-compatible GPU device(s): " << std::endl;

  for (int i = 0; i < nDevices; i++) {
    cudaSetDevice(i);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    std::cout << "  [" << (i + 1) << "] " << prop.name << std::endl;
  }

  cudaSetDevice(Param::preferred_gpu_);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, Param::preferred_gpu_);
  std::cout << std::endl << "Selected GPU [" << Param::preferred_gpu_ << "]: " << prop.name << std::endl << std::endl;
}
#endif

#ifdef USE_OPENCL
template <typename TResourceManager = ResourceManager<>>
static void FindGpuDevicesOpenCL() {
  try {
    // We keep the context and device list in the resource manager to be
    // accessible elsewhere to create command queues and buffers from
    auto rm = TResourceManager::Get();
    cl::Context* context = rm->GetOpenCLContext();
    cl::CommandQueue* queue = rm->GetOpenCLCommandQueue();
    std::vector<cl::Device>* devices = rm->GetOpenCLDeviceList();

    // Get list of OpenCL platforms.
    std::vector<cl::Platform> platform;
    // If we get stuck here, it might be because the DISPLAY envar is not set.
    // Set it to 0 to avoid getting stuck. It's an AMD specific issue
    cl::Platform::get(&platform);

    if (platform.empty()) {
      std::cerr << "OpenCL platforms not found." << std::endl;
    }

    // Go over all available platforms and devices until first device is found
    for (auto p = platform.begin(); p != platform.end();
         p++) {
      std::vector<cl::Device> pldev;

      try {
        p->getDevices(CL_DEVICE_TYPE_GPU, &pldev);

        for (auto d = pldev.begin(); d != pldev.end(); d++) {
          if (!d->getInfo<CL_DEVICE_AVAILABLE>())
            continue;

          // The OpenCL extension available on this device
          std::string ext = d->getInfo<CL_DEVICE_EXTENSIONS>();

          devices->push_back(*d);
        }
      } catch (...) {
        std::cerr << "Found bad OpenCL platform! Continuing to next one." << std::endl;
        devices->clear();
        continue;
      }
    }

    if (devices->empty()) {
      std::cerr << "No OpenCL-compatible GPU found on this machine! Switching to the CPU version..." << std::endl;
      Param::use_gpu_ = false;
      return;
    }

    *context = cl::Context(*devices);

    std::cout << "Found " << devices->size() << " OpenCL-compatible GPU device(s): " << std::endl << std::endl;

    for (size_t i = 0; i < devices->size(); i++) {
      std::cout << "  [" << (i + 1) << "] " << (*devices)[i].getInfo<CL_DEVICE_NAME>() << std::endl;
    }

    std::cout << std::endl << "Selected GPU [" << Param::preferred_gpu_ << "]: " << (*devices)[Param::preferred_gpu_].getInfo<CL_DEVICE_NAME>() << std::endl;


    // Create command queue for that GPU
    cl_int queue_err;
    *queue = cl::CommandQueue(*context, (*devices)[Param::preferred_gpu_], CL_QUEUE_PROFILING_ENABLE, &queue_err);
    cl_ok(queue_err);

  } catch (const cl::Error &err) {
    std::cerr << "OpenCL error: " << err.what() << "(" << err.err() << ")"
              << std::endl;
  }
}

template <typename TResourceManager = ResourceManager<>>
static void CompileOpenCLKernels() {
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
      cl::Program::Sources(1, std::make_pair(displacement_op_opencl_kernel, strlen(displacement_op_opencl_kernel))));

  all_programs->push_back(displacement_op_program);

  // TODO(ahmad): replace all cout with ROOT logging functions
  std::cout << "Compiling OpenCL kernels..." << std::endl;

  std::string options;
  if (Param::opencl_debug_) {
    std::cout << "Building OpenCL kernels with debugging ON" << std::endl;
    options = "-g -O0";
  }

  for (auto& prog : *all_programs) {
    try {
        prog.build(*devices, options.c_str());
    } catch (const cl::Error &) {
      std::cerr << "OpenCL compilation error" << std::endl
                << prog.getBuildInfo<CL_PROGRAM_BUILD_LOG>((*devices)[0])
                << std::endl;
    }
  }
}
#endif

template <typename TResourceManager = ResourceManager<>>
static void InitializeGPUEnvironment() {
  if (Param::use_opencl_) {
#ifdef USE_OPENCL
    FindGpuDevicesOpenCL<>();
    CompileOpenCLKernels<>();
#else
    std::cout << "You tried to use the GPU (OpenCL) version of BioDynaMo, but no OpenCL installation was detected on this machine.\nSwitching to the CPU version..." << std::endl;
    Param::use_gpu_ = false;
#endif
  } else {
#ifdef USE_CUDA
    FindGpuDevicesCuda();
#else
    std::cout << "You tried to use the GPU (CUDA) version of BioDynaMo, but no CUDA installation was detected on this machine.\nSwitching to the CPU version..." << std::endl;
    Param::use_gpu_ = false;
#endif
  }
}

#ifdef USE_OPENCL
const char *getErrorString(cl_int error) {
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

cl_int cl_assert(cl_int const code, char const * const file, int const line, bool const abort) {
  if (code != CL_SUCCESS)
    {
      char const * const err_str = getErrorString(code);

      fprintf(stderr,
              "\"%s\", line %d: cl_assert (%d) = \"%s\"",
              file,line,code,err_str);

      if (abort)
        {
          // stop profiling and reset device here if necessary
          exit(code);
        }
    }

  return code;
}
#endif

}  // namespace bdm

#endif  // GPU_GPU_HELPER_H_
