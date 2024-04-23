if(ENABLE_TESTING)
    include(ExternalProject)
    ExternalProject_Add(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
        PREFIX ${CMAKE_BINARY_DIR}/googletest
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=/usr/local/)
    set(GTEST_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/)
    set(GTEST_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/lib/)
    set(GTEST_LIBRARIES libgtest_main.a libgtest.a libgmock_main.a libgmock.a)
endif()
