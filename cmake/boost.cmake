# find the usual boost
find_package(Boost COMPONENTS log program_options system REQUIRED)
option(BOOST_ASIO_DEBUG "set to true to enable boost asio handler tracking" OFF)
if(BOOST_ASIO_DEBUG)
    set(BOOST_ASIO_DEBUG true)
    add_definitions(-DBOOST_ASIO_ENABLE_HANDLER_TRACKING)
else(BOOST_ASIO_DEBUG)
    set(BOOST_ASIO_DEBUG false)
endif(BOOST_ASIO_DEBUG)
