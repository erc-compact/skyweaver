if(ENABLE_CUDA)
set(LIBRARIES
    ${DEPENDENCY_LIBRARIES}
    -lboost_program_options
    -lboost_system
    -lpthread
    -lhiredis)

include_directories(..)
include_directories(${PROJECT_BINARY_DIR})


set(fbfuse_src

)

cuda_add_library(${CMAKE_PROJECT_NAME} ${fbfuse_src})

#cuda_add_executable(beamforming_cli src/beamforming_cli.cu)
#target_link_libraries(beamforming_cli ${LIBRARIES})

#cuda_add_executable(weightupdater_cli src/weightupdater_cli.cpp)
#target_link_libraries(weightupdater_cli ${LIBRARIES})

#install(TARGETS weightupdater_cli beamforming_cli DESTINATION bin)

#add_subdirectory(profiling)
#add_subdirectory(test)


endif(ENABLE_CUDA)
