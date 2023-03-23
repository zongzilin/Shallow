cmake_minimum_required(VERSION 3.10)

project(swsolve)

# Find MPI
find_package(MPI REQUIRED)

# Find BLAS
find_package(BLAS REQUIRED)

# Find Boost program options
find_package(Boost REQUIRED COMPONENTS program_options)

# Set C++ standard to C++11
set(CMAKE_CXX_STANDARD 11)

# Add executable for loop base version
add_executable(swsolve swsolve.cpp cw.h cw_mpi.h)
target_link_libraries(swsolve PUBLIC ${BLAS_LIBRARIES} MPI::MPI_CXX Boost::program_options)

# Add executable for BLAS version
add_executable(swsolve_blas swsolve_blas.cpp cw_mpi.h cw_blas.h)
target_link_libraries(swsolve_blas PUBLIC ${BLAS_LIBRARIES} MPI::MPI_CXX Boost::program_options)
target_include_directories(swsolve_blas PUBLIC ${BLAS_INCLUDE_DIRS})
target_compile_options(swsolve_blas PUBLIC -O3)

# Add build target for BLAS base version
add_custom_target(BLASBASE DEPENDS swsolve_blas)
add_custom_target(LOOPBASE DEPENDS swsolve)

