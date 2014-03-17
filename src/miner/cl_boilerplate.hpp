#ifndef cl_boilerplate_hpp
#define cl_boilerplate_hpp
#include <vector>
#include <fstream>

#define __CL_ENABLE_EXCEPTIONS
#include "CL/cl.hpp"

#include "bitecoin_log.hpp"

namespace bitecoin {
std::string LoadSource(const char *fileName) {
  // Don't forget to change your_login here
  std::string baseDir = "./";
  if (getenv("HPCE_CL_SRC_DIR")) {
    baseDir = getenv("HPCE_CL_SRC_DIR");
  }

  std::string fullName = baseDir + "/" + fileName;

  std::ifstream src(fullName, std::ios::in | std::ios::binary);
  if (!src.is_open())
    throw std::runtime_error("LoadSource : Couldn't load cl file from '" +
                             fullName + "'.");

  return std::string(
      (std::istreambuf_iterator<char>(src)), // Node the extra brackets.
      std::istreambuf_iterator<char>());
}

void setupOpenCL(std::vector<cl::Platform> &platforms,
                 std::vector<cl::Device> &devices, cl::Device &device,
                 cl::Context &context, std::shared_ptr<ILog> &log) {
  // Initialise OpenCL
  cl::Platform::get(&platforms);
  if (platforms.size() == 0)
    throw std::runtime_error("No OpenCL platforms found.");

  log->Log(Log_Info, "[OpenCL] Found %d platforms", platforms.size());
  for (unsigned i = 0; i < platforms.size(); i++) {
    std::string vendor = platforms[0].getInfo<CL_PLATFORM_VENDOR>();
    log->Log(Log_Info, "[OpenCL]  Platform %i: %s", i, vendor.c_str());
  }

  int selectedPlatform = 0;
  if (getenv("HPCE_SELECT_PLATFORM")) {
    selectedPlatform = atoi(getenv("HPCE_SELECT_PLATFORM"));
  }
  log->Log(Log_Info, "[OpenCL] Choosing platform %d", selectedPlatform);
  cl::Platform platform = platforms.at(selectedPlatform);

  platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
  if (devices.size() == 0) {
    throw std::runtime_error("No opencl devices found.");
  }

  log->Log(Log_Info, "[OpenCL] Found %d devices", devices.size());
  for (unsigned i = 0; i < devices.size(); i++) {
    std::string name = devices[i].getInfo<CL_DEVICE_NAME>();
    log->Log(Log_Info, "[OpenCL]  Device %i: %s", i, name.c_str());
  }

  int selectedDevice = 0;
  if (getenv("HPCE_SELECT_DEVICE")) {
    selectedDevice = atoi(getenv("HPCE_SELECT_DEVICE"));
  }
  log->Log(Log_Info, "[OpenCL] Choosing device %d", selectedDevice);
  device = devices.at(selectedDevice);

  context = cl::Context(devices);

  // Harvesting some information
  size_t value;
  device.getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &value);
  log->Log(Log_Info, "[OpenCL] CL_DEVICE_MAX_COMPUTE_UNITS %d", value);
  device.getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &value);
  log->Log(Log_Info, "[OpenCL] CL_DEVICE_MAX_WORK_GROUP_SIZE %d", value);
  device.getInfo(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, &value);
  log->Log(Log_Info, "[OpenCL] CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS %d", value);
  device.getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &value);
  log->Log(Log_Info, "[OpenCL] CL_DEVICE_MAX_WORK_ITEM_SIZES %d", value);
}
}; // namespace bitecoin
#endif
