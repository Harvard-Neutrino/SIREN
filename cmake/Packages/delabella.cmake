find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        if(NOT EXISTS "${PROJECT_SOURCE_DIR}/vendor/delabella/CMakeLists.txt")
            message(STATUS "Submodule update")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init ./
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/delabella
                            RESULT_VARIABLE GIT_SUBMOD_RESULT)
            if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
            endif()
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/vendor/delabella/CMakeLists.txt")
    message(FATAL_ERROR "The delabella submodule was not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/delabella)
#add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/delabella)
#install(TARGETS delabella_shared
#    ARCHIVE DESTINATION lib)
#set(_config_dir share/delabella/cmake)
#install(FILES "${CMAKE_CURRENT_BINARY_DIR}/delabellaConfigVersion.cmake"
#  DESTINATION ${_config_dir}
#)
#install(EXPORT delabellaConfig
#  DESTINATION ${_config_dir}
#)
