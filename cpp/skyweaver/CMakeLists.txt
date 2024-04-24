include_directories(${CMAKE_SOURCE_DIR}/cpp/)

set(skyweaver_src
   src/DelayManager.cu
   )

set(skyweaver_inc
    DelayManager.cuh
   )

set(SKYWEAVER_LIBRARIES ${CMAKE_PROJECT_NAME} ${DEPENDENCY_LIBRARIES})
# -- the main library target
add_library(${CMAKE_PROJECT_NAME} ${skyweaver_src})
target_link_libraries(${CMAKE_PROJECT_NAME} ${DEPENDENCY_LIBRARIES})

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
    COMMAND gspectrometer_benchmark --benchmark_counters_tabular=true
    COMMAND merger_benchmark --benchmark_counters_tabular=true
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()