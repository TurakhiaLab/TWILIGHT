cmake_minimum_required(VERSION 3.10)

project(twilight LANGUAGES C CXX)

option(USE_CUDA "Build with CUDA" OFF)
option(USE_HIP "Build with HIP" OFF)
option(TWILIGHT_BUILD_EXAMPLES "Build TWILIGHT C API examples" ON)

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
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

set(COMMON_CORE_SRC_FILES
    src/phylogeny.cpp
    src/tree.cpp
    src/partitionInfo.cpp
    src/node.cpp
    src/scoring-matrix.cpp
    src/io.cpp
    src/option.cpp
    src/sequencedb.cpp
    src/progressive.cpp
    src/alignment-cpu.cpp
    src/alignment-helper.cpp
    src/pairwise-align.cpp
    src/profile-align.cpp
    src/accurate.cpp
    src/accurate-sampling.cpp
    src/consistency-library.cpp
    src/accurate-kernel.cpp
    src/accurate-plan.cpp
    src/subtree-refine.cpp
    src/score.cpp
    src/TALCO-XDrop.cpp
    src/twilight_api.cpp
)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    message(STATUS "Building TALCO with x86 SIMD")
    add_compile_definitions(TALCO_SIMD)
endif()

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
			set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo")
		else()
			message(STATUS "No CUDA-capable GPU found. Generating binary to support multiple GPU generations.")
			set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86 89 90)
			set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo")
		endif()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -march=native")
		add_definitions(-DCUDA)
		set(MAIN_SRC_FILE src/cuda/twilight-main.cu)
		set(CORE_SRC_FILES
			src/cuda/gpu-info.cu
			src/cuda/alignment-gpu.cu
			src/cuda/device-function.cu
			${COMMON_CORE_SRC_FILES}
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
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated -march=native")
		set(MAIN_SRC_FILE src/hip/twilight-main.hip.cpp)
		set(CORE_SRC_FILES
			src/hip/gpu-info.hip.cpp
			src/hip/alignment-gpu.hip.cpp
			src/hip/device-function.hip.cpp
			${COMMON_CORE_SRC_FILES}
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
			set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -march=native")
			add_definitions(-DCUDA)
			set(MAIN_SRC_FILE src/cuda/twilight-main.cu)
			set(CORE_SRC_FILES
				src/cuda/gpu-info.cu
				src/cuda/alignment-gpu.cu
				src/cuda/device-function.cu
				${COMMON_CORE_SRC_FILES}
			)
		else()
			message(STATUS "No CUDA-capable GPU found.")
		endif()
	endif()
	if(NOT DEFINED CORE_SRC_FILES OR NOT CORE_SRC_FILES)
		if (HIP_FOUND)
			enable_language(HIP)	
	    	message(STATUS "HIP found: ${HIP_VERSION}")
			add_definitions(-D__HIP_PLATFORM_AMD__)
			set(CMAKE_CXX_COMPILER "/opt/rocm-${HIP_COMPILE_VERSION}/bin/hipcc")
			include_directories("/opt/rocm-${HIP_COMPILE_VERSION}/include/")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated -march=native")
			set(MAIN_SRC_FILE src/hip/twilight-main.hip.cpp)
			set(CORE_SRC_FILES
				src/hip/gpu-info.hip.cpp
				src/hip/alignment-gpu.hip.cpp
				src/hip/device-function.hip.cpp
				${COMMON_CORE_SRC_FILES}
			)
		else()
			message(STATUS "No CUDA or HIP detected — building CPU version.")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated -march=native")
			set(MAIN_SRC_FILE src/twilight-main.cpp)
			set(CORE_SRC_FILES ${COMMON_CORE_SRC_FILES})
		endif()
	endif()
endif()

find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})
find_package(TBB REQUIRED CONFIG)

# Core static library (libtwilight.a) exposing the C API (include/twilight.h)
add_library(twilight_lib STATIC ${CORE_SRC_FILES})
set_target_properties(twilight_lib PROPERTIES OUTPUT_NAME twilight)
target_include_directories(twilight_lib
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${ZLIB_INCLUDE_DIR}
)
target_link_libraries(twilight_lib PUBLIC ${ZLIB_LIBRARIES} TBB::tbb)

# CLI executable
add_executable(twilight ${MAIN_SRC_FILE})
target_include_directories(twilight PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(twilight PRIVATE twilight_lib)

# Pure C99 example executable
if(TWILIGHT_BUILD_EXAMPLES)
    add_executable(c_api_example api/c_api_example.c)
    set_target_properties(c_api_example PROPERTIES
        C_STANDARD 99
        C_STANDARD_REQUIRED ON
    )
    target_link_libraries(c_api_example PRIVATE twilight_lib)
endif()
