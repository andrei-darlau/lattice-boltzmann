cmake_minimum_required(VERSION 3.16)
project(serial_lbm_3d_ogive LANGUAGES CXX)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
add_executable(lbm_3d src/lbm_3d.cpp)
target_compile_features(lbm_3d PRIVATE cxx_std_17)
if(MSVC)
  target_compile_options(lbm_3d PRIVATE /W4 /O2)
else()
  target_compile_options(lbm_3d PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()
enable_testing()
add_test(NAME lbm_3d_smoke COMMAND lbm_3d --nx 28 --ny 28 --nz 64 --steps 20 --every 20 --out smoke_output)

# The solver has no VTK dependency. Enable the reader regression test when
# Python's VTK modules are installed (for example, from ParaView or VTK).
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(COMMAND "${Python3_EXECUTABLE}" -c
    "from vtkmodules.vtkIOXML import vtkXMLImageDataReader; from vtkmodules.vtkFiltersCore import vtkContourFilter"
    RESULT_VARIABLE VTK_PYTHON_IMPORT_RESULT OUTPUT_QUIET ERROR_QUIET)
  if(VTK_PYTHON_IMPORT_RESULT EQUAL 0)
    add_test(NAME lbm_3d_vti_reader
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_vti.py" $<TARGET_FILE:lbm_3d>)
  else()
    message(STATUS "Python VTK modules not found; VTI reader test disabled")
  endif()
endif()
