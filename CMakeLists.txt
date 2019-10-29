cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
cmake_policy(VERSION 3.1.0)


#TODO
# + add support for Apple systems (needs an apple part??)
# + make sure it correctly builds on the cobalt testbed

project(lepton_injector VERSION 1.0.0 DESCRIPTION "lepton injector library")
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_STANDARD 99)
## enable_language(cxx)


find_package(PythonInterp)
find_package(PythonLibs)
set(PYTHON_FOUND PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
if(PYTHON_FOUND)
    MESSAGE("-- Found python at ${PYTHON_INCLUDE_DIR}")
endif()

find_package(HDF5 REQUIRED)
if(HDF5_FOUND)
    MESSAGE("-- Found hdf5 at ${HDF5_INCLUDE_DIR}")
    include_directories(${HDF5_INCLUDE_DIR})
    set(_hdf5_libs hdf5 hdf5_cpp)
endif()

# find boost
set(Boost_USE_STATC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED)
if(Boost_FOUND)
    MESSAGE("-- Found Boost at ${Boost_INCLUDE_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(photospline REQUIRED)

if(photospline_FOUND)
    MESSAGE("-- photospline found")
endif()

include(GNUInstallDirs)

# some apple thing

#core files 
LIST (APPEND core_SOURCES
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
target_include_directories(lepton_injector PUBLIC
    ${CMAKE_SOURCE_DIR}/public/LeptonInjector/
    ${CMAKE_SOURCE_DIR}/public/earthmodel-service/
    ${CMAKE_SOURCE_DIR}/public/phys-services/ 
)
target_link_libraries(lepton_injector ${HDF5_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(lepton_injector cphotospline)

set_target_properties(lepton_injector PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)
###  target_include_directories(lepton_injector PRIVATE .)
install(TARGETS lepton_injector
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# build and install pyhton module
IF(PYTHON_FOUND)
    MESSAGE("-- building python module!")

    LIST( APPEND python_private
        ${CMAKE_SOURCE_DIR}/private/pybindings/LeptonInjector.cxx
    )

    add_library( pylepton_injector SHARED ${python_private})
    set_target_properties( pylepton_injector PROPERTIES LINKER_LANGUAGE CXX OUTPUT_NAME "leptoninjector" PREFIX "" SUFFIX ".SO")
    target_include_directories( pylepton_injector PRIVATE ${PYTHON_INCLUDE_DIR})
    target_link_libraries(pylepton_injector PRIVATE lepton_injector)
    
    EXECUTE_PROCESS( COMMAND ${PYTHON_EXECUTABLE} -c
        "from distutils.sysconfig import get_python_lib; print(get_python_lib(prefix=''))"
        OUTPUT_VARIABLE PYTHON_MODULE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    INSTALL(TARGETS pylepton_injector LIBRARY DESTINATION ${PYTHON_MODULE_DIR})
ELSE(NOT PYTHON FOUND)
    MESSAGE("-- no python. skipping")
ENDIF(PYTHON_FOUND)

