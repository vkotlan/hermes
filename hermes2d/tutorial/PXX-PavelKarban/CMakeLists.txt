project(PXX-vrtak)

add_executable(${PROJECT_NAME} main.cpp data_table.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
