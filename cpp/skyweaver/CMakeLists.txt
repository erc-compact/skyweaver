include_directories(${CMAKE_SOURCE_DIR}/cpp/)

set(skyweaver_src
   src/BufferedDispenser.cu
   src/CoherentBeamformer.cu
   src/CoherentDedisperser.cu
   src/DedispersionPlan.cpp
   src/DelayManager.cu
   src/DescribedVector.cpp
   src/FileOutputStream.cpp
   src/Header.cpp
   src/IncoherentBeamformer.cu
   src/IncoherentDedisperser.cu
   src/MultiFileReader.cu
   src/ObservationHeader.cpp
   src/PipelineConfig.cpp
   src/StatisticsCalculator.cu
   src/Timer.cpp
   src/Transposer.cu
   src/WeightsManager.cu
   src/SigprocHeader.cpp

   )

set(skyweaver_inc
    BufferedDispenser.cuh
    CoherentBeamformer.cuh
    CoherentDedisperser.cuh
    DedispersionPlan.hpp
    DelayManager.cuh
    DescribedVector.hpp
    FileOutputStream.hpp
    Header.hpp
    IncoherentBeamformer.cuh
    IncoherentDedisperser.cuh
    MultiFileReader.cuh
    ObservationHeader.hpp
    PipelineConfig.hpp
    StatisticsCalculator.cuh
    Timer.hpp
    Transposer.cuh
    WeightsManager.cuh
    MultiFileWriter.cuh
    SigprocHeader.hpp
    SkyCleaver.hpp
    SkyCleaverConfig.hpp
   )

set(SKYWEAVER_LIBRARIES ${CMAKE_PROJECT_NAME} ${DEPENDENCY_LIBRARIES})
# -- the main library target
add_library(${CMAKE_PROJECT_NAME} ${skyweaver_src})
target_link_libraries(${CMAKE_PROJECT_NAME} 
    PUBLIC 
    ${DEPENDENCY_LIBRARIES} 
    ${CUDA_cufft_LIBRARY} 
    ${CUDA_NVTOOLSEXT_LIB} 
    OpenMP::OpenMP_CXX)
target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: ${CUDA_NVCC_FLAGS}>)

IF (DEFINED SKYWEAVER_NANTENNAS)
    message(STATUS "Compiling for SKYWEAVER_NANTENNAS = " ${SKYWEAVER_NANTENNAS})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC SKYWEAVER_NANTENNAS=${SKYWEAVER_NANTENNAS})
ENDIF (DEFINED SKYWEAVER_NANTENNAS)

IF (DEFINED SKYWEAVER_NBEAMS)
    message(STATUS "Compiling for SKYWEAVER_NBEAMS = " ${SKYWEAVER_NBEAMS})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC SKYWEAVER_NBEAMS=${SKYWEAVER_NBEAMS})
ENDIF (DEFINED SKYWEAVER_NBEAMS)

IF (DEFINED SKYWEAVER_NCHANS)
    message(STATUS "Compiling for SKYWEAVER_NCHANS = " ${SKYWEAVER_NCHANS})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC SKYWEAVER_NCHANS=${SKYWEAVER_NCHANS})
ENDIF (DEFINED SKYWEAVER_NCHANS)

IF (DEFINED SKYWEAVER_IB_SUBTRACTION)
    message(STATUS "Compiling for SKYWEAVER_IB_SUBTRACTION = " ${SKYWEAVER_IB_SUBTRACTION})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC SKYWEAVER_IB_SUBTRACTION=${SKYWEAVER_IB_SUBTRACTION})
ENDIF (DEFINED SKYWEAVER_IB_SUBTRACTION)

IF (DEFINED SKYWEAVER_CB_TSCRUNCH)
    message(STATUS "Compiling for SKYWEAVER_CB_TSCRUNCH = " ${SKYWEAVER_CB_TSCRUNCH})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC SKYWEAVER_CB_TSCRUNCH=${SKYWEAVER_CB_TSCRUNCH})
ENDIF (DEFINED SKYWEAVER_CB_TSCRUNCH)

IF (DEFINED SKYWEAVER_IB_TSCRUNCH)
    message(STATUS "Compiling for SKYWEAVER_IB_TSCRUNCH = " ${SKYWEAVER_IB_TSCRUNCH})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC SKYWEAVER_IB_TSCRUNCH=${SKYWEAVER_IB_TSCRUNCH})
ENDIF (DEFINED SKYWEAVER_IB_TSCRUNCH)

IF (DEFINED SKYWEAVER_CB_FSCRUNCH)
    message(STATUS "Compiling for SKYWEAVER_CB_FSCRUNCH = " ${SKYWEAVER_CB_FSCRUNCH})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC SKYWEAVER_CB_FSCRUNCH=${SKYWEAVER_CB_FSCRUNCH})
ENDIF (DEFINED SKYWEAVER_CB_FSCRUNCH)

IF (DEFINED SKYWEAVER_IB_FSCRUNCH)
    message(STATUS "Compiling for SKYWEAVER_IB_FSCRUNCH = " ${SKYWEAVER_IB_FSCRUNCH})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC SKYWEAVER_IB_FSCRUNCH=${SKYWEAVER_IB_FSCRUNCH})
ENDIF (DEFINED SKYWEAVER_IB_FSCRUNCH)

IF (DEFINED USE_NVTX)
    message(STATUS "Compiling with USE_NVTX = " ${USE_NVTX})
    target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC USE_NVTX=${USE_NVTX})
ENDIF (DEFINED USE_NVTX)

# CLI applications
cuda_add_executable(skyweavercpp src/skyweaver_cli.cu)
target_link_libraries(skyweavercpp 
    ${SKYWEAVER_LIBRARIES} 
    ${DEPENDENCY_LIBRARIES} 
    ${CUDA_cufft_LIBRARY} 
    ${CUDA_NVTOOLSEXT_LIB} 
    OpenMP::OpenMP_CXX)
install(TARGETS skyweavercpp DESTINATION bin)

add_executable(skycleaver src/skycleaver_cli.cpp)
target_link_libraries(skycleaver 
    ${SKYWEAVER_LIBRARIES} 
    ${DEPENDENCY_LIBRARIES} 
    OpenMP::OpenMP_CXX)
install(TARGETS skycleaver DESTINATION bin)

#install (TARGETS ... DESTINATION bin)
#install (TARGETS ${CMAKE_PROJECT_NAME}
#         RUNTIME DESTINATION bin
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib)
install(FILES ${skyweaver_inc} DESTINATION include/skyweaver)

# Install testing tools if testing is enabled
if(ENABLE_TESTING)
    include_directories(${GTEST_INCLUDE_DIR})
    link_directories(${GTEST_LIBRARY_DIR})
    add_subdirectory(test)
endif()

# Add benchmark binaries here to get executed with make benchmark
if(ENABLE_BENCHMARK)
  add_custom_target(benchmark
    COMMAND beamformer_bench --benchmark_counters_tabular=true
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  add_subdirectory(benchmark)
endif()
