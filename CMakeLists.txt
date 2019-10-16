i3_project(LeptonInjector
          PYTHON_DIR python)


IF(HDF5_FOUND)
    
    message(STATUS "     + hdf5 support")

    SET(LIB_${PROJECT_NAME}_SOURCEFILES private/LeptonInjector/LeptonInjector.cxx private/LeptonInjector/serialization.cxx private/LeptonInjector/h5write.cxx)
    SET(LIB_${PROJECT_NAME}_PROJECTS serialization icetray dataclasses phys-services sim-services earthmodel-service photospline)
    SET(LIB_${PROJECT_NAME}_TOOLS boost python)

    if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)
        LIST(APPEND LIB_${PROJECT_NAME}_PROJECTS tableio)
        LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES private/LeptonInjector/converter/LeptonInjectionConfigurationConverter.cxx)
        ADD_DEFINITIONS(-DUSE_TABLEIO)
        MESSAGE(STATUS "     + tableio support")
    endif(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)


    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

    i3_add_library(${PROJECT_NAME} 
      ${LIB_${PROJECT_NAME}_SOURCEFILES}
      USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS} hdf5
      USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS} 
    )


ELSE(HDF5_FOUND)
    
    colormsg(CYAN "+-- hdf5 library not found, skipping hdf5 support!")

    SET(LIB_${PROJECT_NAME}_SOURCEFILES private/LeptonInjector/LeptonInjector.cxx private/LeptonInjector/serialization.cxx)
    SET(LIB_${PROJECT_NAME}_PROJECTS serialization icetray dataclasses phys-services sim-services earthmodel-service photospline)
    SET(LIB_${PROJECT_NAME}_TOOLS boost python)

    if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)
        LIST(APPEND LIB_${PROJECT_NAME}_PROJECTS tableio)
        LIST(APPEND LIB_${PROJECT_NAME}_SOURCEFILES private/LeptonInjector/converter/LeptonInjectionConfigurationConverter.cxx)
        ADD_DEFINITIONS(-DUSE_TABLEIO)
        MESSAGE(STATUS "     + tableio support")
    endif(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/tableio)

    i3_add_library(${PROJECT_NAME} 
      ${LIB_${PROJECT_NAME}_SOURCEFILES}
      USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS} 
      USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS} 
    )

ENDIF(HDF5_FOUND)

if(BUILD_PYBINDINGS)
  i3_add_pybindings(LeptonInjector
    private/pybindings/LeptonInjector.cxx
    USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS} 
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS} LeptonInjector
    )
endif(BUILD_PYBINDINGS)

i3_test_executable(test 
  private/test/main.cxx
  private/test/RangedInjection.cxx
  private/test/VolumeInjection.cxx
  private/test/MultiInjection.cxx
  private/test/Utils.cxx
  private/test/tools.cxx
  USE_TOOLS gsl
  USE_PROJECTS icetray dataclasses dataio phys-services earthmodel-service LeptonInjector
)

#i3_executable(exercise-final-states 
#  private/exercise-final-states/main.cxx
#  USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
#  USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS} LeptonInjector
#)
#
#SET_SOURCE_FILES_PROPERTIES(
#	private/exercise-final-states/main.cxx
#	PROPERTIES
#	COMPILE_FLAGS -std=c++11
#)
