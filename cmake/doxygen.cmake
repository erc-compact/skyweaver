
# ------------------------------------------------------------------
# Documentation
# ------------------------------------------------------------------
find_package(Doxygen)
find_program(SPHINX_EXECUTABLE
    NAMES sphinx-build
    DOC "Path to sphinx-build executable")

if(DOXYGEN_FOUND AND SPHINX_EXECUTABLE)
	MESSAGE(STATUS "Found Doxygen and Sphinx to build documentation")
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/DoxygenLayout.xml ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml COPYONLY)
	add_custom_target(doxy ${DOXYGEN_EXECUTABLE}
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR} COMMENT
        "Generating API documentation with Doxygen" VERBATIM)

	include(FindPackageHandleStandardArgs)
	find_package_handle_standard_args(Sphinx
	                                  "Failed to find sphinx-build executable"
	                                  SPHINX_EXECUTABLE)
	set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/doc)
	set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/doc)
	add_custom_target(doc
        ${SPHINX_EXECUTABLE} -b html
            -Dbreathe_projects.psrdada-cpp=${CMAKE_CURRENT_BINARY_DIR}/xml
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Sphinx")
	add_dependencies(doc doxy)
else()
	MESSAGE(STATUS "Doxygen and Sphinx not found: cannot build documentation")
endif(DOXYGEN_FOUND AND SPHINX_EXECUTABLE)
