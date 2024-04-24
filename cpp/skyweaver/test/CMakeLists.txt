include_directories(${GTEST_INCLUDE_DIR})

link_directories(${GTEST_LIBRARY_DIR})

set(
    gtest_skyweaver_src
    src/DelayManagerTester.cu
    src/gtest_skyweaver.cu
)
cuda_add_executable(gtest_skyweaver ${gtest_skyweaver_src} )
target_link_libraries(gtest_skyweaver ${SKYWEAVER_LIBRARIES} ${GTEST_LIBRARIES})
add_test(NAME gtest_skyweaver 
         COMMAND gtest_skyweaver 
         WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")