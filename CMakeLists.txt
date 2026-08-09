cmake_minimum_required(VERSION 3.10)

project(twilight LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
			set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86 89 90)
		endif()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -march=native")
		set(SRC_FILES src/cuda/twilight-main.cu
					  src/cuda/gpu-info.cu
				      src/cuda/alignment-gpu.cu
				      src/cuda/device-function.cu
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
					  src/TALCO-XDrop.cpp
					  
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
		set(SRC_FILES src/hip/twilight-main.hip.cpp
					  src/hip/gpu-info.hip.cpp
				      src/hip/alignment-gpu.hip.cpp
				      src/hip/device-function.hip.cpp
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
					  src/TALCO-XDrop.cpp
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
				set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_${GPU_ARCH_OUTPUT} -G -g -lineinfo")
			else()
				set(CMAKE_CUDA_ARCHITECTURES ${GPU_ARCH_OUTPUT})
			endif()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -march=native")
			set(SRC_FILES src/cuda/twilight-main.cu
					      src/cuda/gpu-info.cu
				          src/cuda/alignment-gpu.cu
				          src/cuda/device-function.cu
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
					      src/TALCO-XDrop.cpp
			)
		else()
			message(STATUS "No CUDA-capable GPU found.")
		endif()
	endif()
	if(NOT DEFINED SRC_FILES OR NOT SRC_FILES)
		if (HIP_FOUND)
			enable_language(HIP)	
	    	message(STATUS "HIP found: ${HIP_VERSION}")
			add_definitions(-D__HIP_PLATFORM_AMD__)
			set(CMAKE_CXX_COMPILER "/opt/rocm-${HIP_COMPILE_VERSION}/bin/hipcc")
			include_directories("/opt/rocm-${HIP_COMPILE_VERSION}/include/")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated -march=native")
			set(SRC_FILES src/hip/twilight-main.hip.cpp
					      src/hip/gpu-info.hip.cpp
				          src/hip/alignment-gpu.hip.cpp
				          src/hip/device-function.hip.cpp
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
					      src/TALCO-XDrop.cpp
	    	)
		else()
			message(STATUS "No CUDA or HIP detected — building CPU version.")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -Wno-self-assign -Wno-deprecated -march=native")
			set(SRC_FILES src/twilight-main.cpp
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
					      src/TALCO-XDrop.cpp
		    )	
		endif()
	endif()
endif()

# add_executable (twilight_mga
# 				src_mga/twilight_mga-main.cpp
# 				src_mga/alignment_new.cpp
# 				src_mga/option.cpp
# 				src_mga/tree.cpp
# 				src_mga/node.cpp
# 				src_mga/block.cpp
# 				src_mga/io.cpp
# 				src_mga/progressive.cpp
# 				src_mga/alignment_kernel_new.cpp
# 				src_mga/blockManager.cpp
# 				src_mga/blockSet.cpp
# 				src_mga/blockMerger.cpp
# 				src_mga/SequenceInfo.cpp
# 				src_mga/selfmapping.cpp
# 				src_mga/validate.cpp
# 				src_mga/global_alignment.cpp
# 				src_mga/realignShallowBlocks.cpp
# 				src_mga/refine.cpp)

add_executable (twilight_mga
				mga/src/twilight_mga-main.cpp
				mga/src/option.cpp
				mga/src/tree.cpp
				mga/src/node.cpp
				mga/src/io.cpp
				mga/src/progressive.cpp
				mga/src/alignment.cpp
				mga/src/coverageTracker.cpp
				mga/src/alignment_util.cpp
				mga/src/minimap2_util.cpp
				mga/src/alignment_kernel.cpp
				mga/src/block.cpp
				mga/src/blockManager.cpp
				mga/src/blockSet.cpp
				mga/src/blockMerger.cpp
				mga/src/merge.cpp
				mga/src/consensus.cpp
				mga/src/Sequence.cpp
				mga/src/validate.cpp
				mga/src/refine.cpp
				mga/src/global_alignment.cpp
				mga/src/self_align.cpp
				)

find_package(Boost COMPONENTS program_options system filesystem REQUIRED)
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR} mga/internal ${CMAKE_CURRENT_SOURCE_DIR}/minimap2)
find_package(TBB REQUIRED CONFIG)
get_target_property(TBB_INC_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
if(TBB_INC_DIRS)
    include_directories(${TBB_INC_DIRS})
endif()
if(TBB_INCLUDE_DIRS)
    include_directories(${TBB_INCLUDE_DIRS})
endif()

# add_executable (twilight ${SRC_FILES})
# target_link_libraries(twilight PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb)
target_link_libraries(twilight_mga PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/minimap2/libminimap2.a
    ${Boost_LIBRARIES} 
    ${ZLIB_LIBRARIES} 
    TBB::tbb
    pthread
    m
)


# =========================================================================
# 🌟 加入 AddressSanitizer (ASan) 進行記憶體除錯
# =========================================================================
# 1. 加入編譯標記 (-g 保留除錯資訊，-fno-omit-frame-pointer 讓 Stack Trace 更完整)
target_compile_options(twilight_mga PRIVATE -fsanitize=address -g -fno-omit-frame-pointer)

# 2. 加入連結標記 (ASan 必須在連結期綁定特殊的記憶體管理函式庫)
target_link_options(twilight_mga PRIVATE -fsanitize=address)