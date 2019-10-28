cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
cmake_policy(VERSION 3.1.0)

project(lepton_injector VERSION 1.0.0 DESCRIPTION "lepton injector library")
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_STANDARD 99)

INCLUDE (Python)
INCLUDE (HDF5)
INCLUDE (BOOST)

include(GNUInstallDirs)

# some apple thing

#core files 
LIST (append core_SOURCES
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/Controller.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/Coordinates.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/DataWriter.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/EventProps.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/LeptonInjector.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/Particle.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/Random.cxx
    ${CMAKE_SOURCE_DIR}/private/earthmodel-service/EarthModelCalculator.cxx
    ${CMAKE_SOURCE_DIR}/private/earthmodel-service/EarthModelService.cxx
    ${CMAKE_SOURCE_DIR}/private/phys-services/LICrossSection.cxx
)
add_library( lepton_injector SHARED ${core_SOURCES} )



set_target_properties(lepton_injector PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)
configure_file(lepton_injector.pc.in lepton_injector.pc @ONLY)
target_include_directories(lepton_injector PRIVATE .)
install(TARGETS lepton_injector
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/lepton_injector.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


# build and install pyhton module
IF(PYTHON_FOUND)
    MESSAGE("-- building python module!")

    LIST( APPEND python_private
        ${CMAKE_SOURCE_DIR}/private/pybindings/LeptonInjector.cxx
    )

    add_library( pylepton_injector SHARED ${python_private})
    set_target_properties( pylepton_injector PROPERTIES LINKER_LANGUAGE C++ OUTPUT_NAME "leptoninjector" PREFIX "" SUFFIX ".SO")
    target_include_directories( pylepton_injector PRIVATE ${PYTHON_INCLUDE_DIR})
    target_link_libraries(pylepton_injector PRIVATE lepton_injector)
    
    EXECUTE_PROCESS( COMMAND ${PYTHON_EXECUTABLE} -c
        "from distutils.sysconfig import get_python_lib; print(get_python_lib(prefix=''))"
        OUTPUT_VARIABLE PYTHON_MODULE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    INSTALL(TARGETS pylepton_injector LIBRARY DESTINATION ${PYTHON_MODULE_DIR})
ELSE(NOT PYTHON FOUND)
    MESSAGE("-- no python. skipping")
ENDIF(PYTHON_FOUND)

