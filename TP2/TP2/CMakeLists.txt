cmake_minimum_required(VERSION 2.8.5 FATAL_ERROR)
if(POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW) # CMake 3.0
endif()
if(POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW) # CMake 3.1
endif()

PROJECT (Medical)

find_package(VTK COMPONENTS
  vtkFiltersGeometry
  vtkFiltersModeling
  vtkIOImage
  vtkIOXML
  vtkImagingStatistics
  vtkInteractionStyle
  vtkRenderingVolumeOpenGL
  vtkTestingRendering
  vtksys
 )

#FIND_PACKAGE(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(Medical1 Medical1.cxx)
add_executable(Medical4 Medical4.cxx)

target_link_libraries(Medical1 ${VTK_LIBRARIES})
target_link_libraries(Medical4 ${VTK_LIBRARIES})


