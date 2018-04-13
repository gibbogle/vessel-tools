cmake_minimum_required(VERSION 2.6)
 
PROJECT(showpath)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


set(PROJECTNAME "showpath")
ADD_EXECUTABLE(showpath showpath.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(showpath ${VTK_LIBRARIES})
else()
  target_link_libraries(showpath vtkHybrid vtkWidgets)
endif()

ADD_CUSTOM_COMMAND(TARGET ${PROJECTNAME} POST_BUILD COMMAND copy \"$(TargetPath)\" \"../bin/exec\")
