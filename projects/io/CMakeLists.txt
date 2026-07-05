
#core files
LIST (APPEND io_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/io/private/HepMC3Writer.cxx
    ${PROJECT_SOURCE_DIR}/projects/io/private/HepMC3Reader.cxx
)
add_library(SIREN_io OBJECT ${io_SOURCES})
set_property(TARGET SIREN_io PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(SIREN_io PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/io/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/io/>
)

target_link_libraries(SIREN_io
    PUBLIC
        SIREN_utilities
        SIREN_dataclasses
)

if(HEPMC3_FOUND)
    target_link_libraries(SIREN_io PUBLIC HepMC3::HepMC3)
    target_compile_definitions(SIREN_io PUBLIC SIREN_HAS_HEPMC3=1)
endif()

apply_siren_compile_options(SIREN_io)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/io/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

if(NOT CIBUILDWHEEL)
    if(HEPMC3_FOUND)
        package_add_test(UnitTest_HepMC3Writer ${PROJECT_SOURCE_DIR}/projects/io/private/test/HepMC3Writer_TEST.cxx)
        target_link_libraries(UnitTest_HepMC3Writer HepMC3::HepMC3)
        package_add_test(UnitTest_HepMC3RoundTrip ${PROJECT_SOURCE_DIR}/projects/io/private/test/HepMC3RoundTrip_TEST.cxx)
        target_link_libraries(UnitTest_HepMC3RoundTrip HepMC3::HepMC3)
    endif()
endif()

# The Python module is named `hepmc3` (not `io`) to avoid shadowing the Python
# standard-library `io` module, which breaks importing the siren package.
pybind11_add_module(hepmc3 ${PROJECT_SOURCE_DIR}/projects/io/private/pybindings/io.cxx)
target_link_libraries(hepmc3 PRIVATE SIREN pybind11::embed)
siren_strip_module(hepmc3)
