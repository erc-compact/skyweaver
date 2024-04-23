find_package(CUDA REQUIRED)
include_directories(${CUDA_TOOLKIT_INCLUDE})

set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
set(CUDA_PROPAGATE_HOST_FLAGS OFF)

list(APPEND CUDA_NVCC_FLAGS --std=c++${CMAKE_CXX_STANDARD} -Wno-deprecated-gpu-targets --ptxas-options=-v)
list(APPEND CUDA_NVCC_FLAGS_DEBUG -O2 -Xcompiler "-Wextra" --Werror all-warnings -Xcudafe "--diag_suppress=20012")
list(APPEND CUDA_NVCC_FLAGS_PROFILE --generate-line-info)
list(APPEND CUDA_NVCC_FLAGS_RELEASE -O3 -use_fast_math -restrict)
#list(APPEND CUDA_NVCC_FLAGS_RELEASE -gencode=arch=compute_61,code=sm_61) # Titan X Pascal
#list(APPEND CUDA_NVCC_FLAGS_RELEASE -gencode=arch=compute_75,code=sm_75) # GeForce 2080
#list(APPEND CUDA_NVCC_FLAGS_RELEASE -gencode=arch=compute_80,code=sm_80) # A100
#if(CUDA_VERSION GREATER_EQUAL 11.1)
#    message(STATUS "Enabling device specific (arch=86)")
#    list(APPEND CUDA_NVCC_FLAGS_RELEASE -gencode=arch=compute_86,code=sm_86) # GeForce 3090
#endif(CUDA_VERSION GREATER_EQUAL 11.1)
if(CUDA_VERSION GREATER_EQUAL 11.8)
    message(STATUS  "Enabling device specific (arch=89)")
    list(APPEND CUDA_NVCC_FLAGS_RELEASE -gencode=arch=compute_89,code=sm_89) # L40
endif(CUDA_VERSION GREATER_EQUAL 11.8)

# There is some kind of bug here, as setting CMAKE_CUDA_ARCHITECTURES
# does not add compiler flags to the NVCC call, yet omitting this 
# will cause CMAKE to fail.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 89)
endif()

set(CMAKE_CXX_FLAGS "-DENABLE_CUDA ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE_UPPER}")
if(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_NVCC_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE_UPPER STREQUAL "RELEASE")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_NVCC_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE_UPPER STREQUAL "PROFILE")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_NVCC_FLAGS_PROFILE}")
endif()


get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
add_definitions(-DUSE_NVTX)
add_definitions(-DENABLE_CUDA)
find_library(CUDA_NVTOOLSEXT_LIB
    NAMES nvToolsExt
    HINTS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
    list(APPEND CUDA_DEPENDENCY_LIBRARIES ${CUDA_NVTOOLSEXT_LIB})

macro(CUDA_GENERATE_LINK_FILE cuda_link_file cuda_target)

    # Compute the file name of the intermedate link file used for separable
    # compilation.
    CUDA_COMPUTE_SEPARABLE_COMPILATION_OBJECT_FILE_NAME(cuda_link_file ${cuda_target} "${${cuda_target}_SEPARABLE_COMPILATION_OBJECTS}")

    # Add a link phase for the separable compilation if it has been enabled.  If
    # it has been enabled then the ${cuda_target}_SEPARABLE_COMPILATION_OBJECTS
    # variable will have been defined.
    CUDA_LINK_SEPARABLE_COMPILATION_OBJECTS("${cuda_link_file}" ${cuda_target} "${_options}" "${${cuda_target}_SEPARABLE_COMPILATION_OBJECTS}")

endmacro()

#
# @macro CUDA_SUBPACKAGE_COMPILE
# Use to specify that certain objects should be compiled with
# alternative flags (e.g. a specific architecture)
# @example
# cuda_subpackage_compile(${my_cuda_files} OPTIONS "-arch compute_35")
#
macro(CUDA_SUBPACKAGE_COMPILE)
    FILE(APPEND ${SUBPACKAGE_FILENAME}
    "CUDA_ADD_CUDA_INCLUDE_ONCE()\n"
    "cuda_compile(_cuda_objects "
    )
    foreach(arg ${ARGN})
        if(EXISTS "${arg}")
            set(_file "${arg}")
        else()
            set(_file "${CMAKE_CURRENT_SOURCE_DIR}/${arg}")
            if(NOT EXISTS "${_file}")
                set(_file ${arg})
            endif()
        endif()
        FILE(APPEND ${SUBPACKAGE_FILENAME}
            "${_file}\n"
        )
    endforeach(arg ${ARGV})
    FILE(APPEND ${SUBPACKAGE_FILENAME}
    ")\n"
    "list(APPEND lib_obj_cuda \${_cuda_objects})\n"
    )
endmacro(CUDA_SUBPACKAGE_COMPILE)
