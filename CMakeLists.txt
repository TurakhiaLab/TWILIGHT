cmake_minimum_required(VERSION 3.10)

project(twilight LANGUAGES CXX)

option(USE_CUDA "Build with CUDA" OFF)
option(USE_HIP "Build with HIP" OFF)

include(CheckLanguage)

# Configure TWILIGHT version
set(twilight_VERSION "?.?.?")
execute_process(COMMAND head -n1 ${CMAKE_SOURCE_DIR}/version.txt OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
if ("${VERSION_PATCH}" MATCHES "[0-9]+")
    set(twilight_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()
message(STATUS "TWILIGHT version ${twilight_VERSION}")
configure_file(src/version.hpp.in ${CMAKE_SOURCE_DIR}/src/version.hpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

if(USE_CUDA)
	message(STATUS "Build with CUDA")
	check_language(CUDA)
	if(CMAKE_CUDA_COMPILER)
		enable_language(CUDA)
		message(STATUS "CUDA compiler detected at: ${CMAKE_CUDA_COMPILER}")
		TRY_RUN(
			FIND_GPU_RESULT
			FIND_GPU_COMPILE_RESULT
			${CMAKE_CURRENT_BINARY_DIR}/
			${CMAKE_CURRENT_BINARY_DIR}/../src/cuda/findgpu.cu
			RUN_OUTPUT_VARIABLE GPU_ARCH_OUTPUT
		)
		if("${FIND_GPU_COMPILE_RESULT}" AND ("${FIND_GPU_RESULT}" EQUAL 0))
			message(STATUS "CUDA-capable GPU found. GPU architecture: ${GPU_ARCH_OUTPUT}")
			if(${CMAKE_VERSION} VERSION_LESS "3.18")
				set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_${GPU_ARCH_OUTPUT}")
			else()
				set(CMAKE_CUDA_ARCHITECTURES ${GPU_ARCH_OUTPUT})
			endif()
		else()
			message(STATUS "No CUDA-capable GPU found. Generating binary to support multiple GPU generations.")
			set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo -gencode arch=compute_50,code=sm_50 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_80,code=sm_80")
		endif()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable")
		set(SRC_FILES src/cuda/twilight-main.cu
					  src/tree.cpp
					  src/setting.cpp
					  src/treePartition.cpp
					  src/TALCO-XDrop.cpp
					  src/util.cpp
					  src/msa-cpu.cpp
					  src/cuda/msa-gpu.cu
					  src/cuda/align.cu
		)
	else()
	    message(FATAL_ERROR "No CUDA found")
	endif()
elseif(USE_HIP)
	message(STATUS "Build with HIP")
	find_package(HIP QUIET)
	if(HIP_FOUND)
		enable_language(HIP)
	    message(STATUS "HIP found: ${HIP_VERSION}")
		add_definitions(-D__HIP_PLATFORM_AMD__)
		set(CMAKE_CXX_COMPILER "/opt/rocm-${HIP_COMPILE_VERSION}/bin/hipcc")
		include_directories("/opt/rocm-${HIP_COMPILE_VERSION}/include/")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated")
		set(SRC_FILES src/hip/twilight-main.hip.cpp
			  		  src/tree.cpp
			  		  src/setting.cpp
			  		  src/treePartition.cpp
			  		  src/TALCO-XDrop.cpp
			  		  src/util.cpp
			  		  src/msa-cpu.cpp
			  		  src/hip/msa-gpu.hip.cpp
			  		  src/hip/align.hip.cpp
	    )	
	else()
	    message(FATAL_ERROR "No HIP found")
	endif()
else()
	message(STATUS "Auto-detecting hardware and selecting compute backend.")
	find_package(HIP QUIET)
	check_language(CUDA)
	if(CMAKE_CUDA_COMPILER)
		message(STATUS "CUDA compiler detected at: ${CMAKE_CUDA_COMPILER}")
		enable_language(CUDA)
		TRY_RUN(
    		FIND_GPU_RESULT
			FIND_GPU_COMPILE_RESULT
			${CMAKE_CURRENT_BINARY_DIR}/
			${CMAKE_CURRENT_BINARY_DIR}/../src/cuda/findgpu.cu
			RUN_OUTPUT_VARIABLE GPU_ARCH_OUTPUT
		)
		if("${FIND_GPU_COMPILE_RESULT}" AND ("${FIND_GPU_RESULT}" EQUAL 0))
			message(STATUS "CUDA-capable GPU found. GPU architecture: ${GPU_ARCH_OUTPUT}")
			if(${CMAKE_VERSION} VERSION_LESS "3.18")
				set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_${GPU_ARCH_OUTPUT}")
			else()
				set(CMAKE_CUDA_ARCHITECTURES ${GPU_ARCH_OUTPUT})
			endif()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable")
			set(SRC_FILES src/cuda/twilight-main.cu
					      src/tree.cpp
					      src/setting.cpp
					      src/treePartition.cpp
					      src/TALCO-XDrop.cpp
					      src/util.cpp
					      src/msa-cpu.cpp
						  src/cuda/msa-gpu.cu
					      src/cuda/align.cu
			)
		else()
			message(STATUS "No CUDA-capable GPU found.")
		endif()
		
	elseif(HIP_FOUND)
		enable_language(HIP)	
	    message(STATUS "HIP found: ${HIP_VERSION}")
		add_definitions(-D__HIP_PLATFORM_AMD__)
		set(CMAKE_CXX_COMPILER "/opt/rocm-${HIP_COMPILE_VERSION}/bin/hipcc")
		include_directories("/opt/rocm-${HIP_COMPILE_VERSION}/include/")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated")
		set(SRC_FILES src/hip/twilight-main.hip.cpp
			  		  src/tree.cpp
			  		  src/setting.cpp
			  		  src/treePartition.cpp
			  		  src/TALCO-XDrop.cpp
			  		  src/util.cpp
			  		  src/msa-cpu.cpp
			  		  src/hip/msa-gpu.hip.cpp
			  		  src/hip/align.hip.cpp
	    )
	else()
		message(STATUS "No CUDA or HIP detected — building CPU version.")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated")
		set(SRC_FILES src/twilight-main.cpp
				      src/tree.cpp
				      src/setting.cpp
				      src/treePartition.cpp
				      src/TALCO-XDrop.cpp
				      src/util.cpp
				      src/msa-cpu.cpp
	    )	
	endif()
endif()

find_package(Boost COMPONENTS program_options system filesystem REQUIRED)
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})
find_package(TBB REQUIRED CONFIG)
add_executable (twilight ${SRC_FILES})
target_link_libraries(twilight PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb)
