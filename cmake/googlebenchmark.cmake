if(ENABLE_BENCHMARK AND BUILD_SUBMODULES)
    include(ExternalProject)
    ExternalProject_Add(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.8.3
        PREFIX ${CMAKE_BINARY_DIR}/googlebenchmark
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=/usr/local/ -DBENCHMARK_ENABLE_TESTING=FALSE -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        DEPENDS googletest)
    set(BENCHMARK_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/)
    set(BENCHMARK_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/lib/)
    set(BENCHMARK_LIBRARIES libbenchmark_main.a libbenchmark.a)
elseif(ENABLE_BENCHMARK)
    message(STATUS "ENABLE_BENCHMARK but BUILD_SUBMODULES disabled, not installing googlebenchmark")
endif()
