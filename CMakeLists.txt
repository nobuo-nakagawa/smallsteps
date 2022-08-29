cmake_minimum_required(VERSION 3.6)

project(smallsteps)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}/bin)
add_executable(smallsteps ./src/smallsteps.cpp
                      ./src/imgui/imgui.cpp
                      ./src/imgui/imgui_draw.cpp
                      ./src/imgui/imgui_demo.cpp
                      ./src/imgui/imgui_widgets.cpp
                      ./src/imgui/imgui_impl_glut.cpp
                      ./src/imgui/imgui_impl_opengl2.cpp
)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT smallsteps)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
set(CMAKE_SUPPRESS_REGENERATION true)

if (WIN32)
set(OPENGL_LIBRARY_DIR "${PROJECT_SOURCE_DIR}/freeglut/lib/x64")
set(GLUT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}freeglut/include)
endif()

if(MSVC)
add_definitions(-D_UNICODE -DUNICODE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest")

set(CompilerFlags
  CMAKE_CXX_FLAGS
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_RELEASE
)
foreach(CompilerFlag ${CompilerFlags})
  string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
endforeach()

add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PROJECT_SOURCE_DIR}/freeglut/bin/x64/freeglut.dll" "${CMAKE_BINARY_DIR}"
  COMMAND ${CMAKE_COMMAND} -E make_directory    "${PROJECT_SOURCE_DIR}/bin"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PROJECT_SOURCE_DIR}/freeglut/bin/x64/freeglut.dll" "${PROJECT_SOURCE_DIR}/bin"
)
endif(MSVC)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

if (WIN32)
include_directories( ${PROJECT_SOURCE_DIR}/freeglut/include )
else()
include_directories( ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} )
endif()

include_directories( ${PROJECT_SOURCE_DIR}/src/imgui )

target_link_libraries(smallsteps ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
