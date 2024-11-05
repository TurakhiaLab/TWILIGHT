cmake_minimum_required(VERSION 3.10)

project(msa LANGUAGES CXX)

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
	# For CUDA
	enable_language(CUDA)
	set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
	find_package(CUDA 10.0 REQUIRED)
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_86 -lineinfo")
	set(SRC_FILES src/msa-main.cu
				  src/align.cu	
				  src/msa-gpu.cu
				  src/tree.cpp
				  src/setting.cpp
				  src/treePartition.cpp
				  src/TALCO-XDrop.cpp
				  src/util.cpp
				  src/msa-cpu.cpp
				  )
else()
  	message(STATUS "No CUDA support.")
	set(SRC_FILES src/msa-main.cpp
				  src/tree.cpp
				  src/setting.cpp
				  src/treePartition.cpp
				  src/TALCO-XDrop.cpp
				  src/util.cpp
				  src/msa-cpu.cpp
	              )	
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g")


find_package(Boost COMPONENTS program_options REQUIRED)

find_package(ZLIB REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)



# add_executable (msa-gpu
# 	        	src/msa-gpu.cu
# 				src/align.cu
# 		    	src/tree.cpp
# 				src/setting.cu
# 				src/treePartition.cpp
# 				src/TALCO-XDrop.cpp
# 				src/util.cpp
# 				src/msa.cu)
# add_executable (msa-cpu
# 	        	src/msa-cpu.cpp
# 		    	src/tree.cpp
# 				src/setting.cpp
# 				src/treePartition.cpp
# 				src/TALCO-XDrop.cpp
# 				src/util.cpp
# 				src/msa.cpp)

add_executable (msa-accel ${SRC_FILES})
target_link_libraries(msa-accel PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS} stdc++fs)

# target_link_libraries(msa-gpu PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS})
# target_link_libraries(msa-cpu PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS})
