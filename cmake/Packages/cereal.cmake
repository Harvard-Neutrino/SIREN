find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        if(NOT EXISTS "${PROJECT_SOURCE_DIR}/vendor/cereal/CMakeLists.txt")
            message(STATUS "Submodule update")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init ./
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/cereal
                            RESULT_VARIABLE GIT_SUBMOD_RESULT)
            if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
            endif()
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/vendor/cereal/CMakeLists.txt")
    message(FATAL_ERROR "The cereal submodule was not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

set(JUST_INSTALL_CEREAL ON CACHE INTERNAL "Cereal just install library.")
set(SKIP_PORTABILITY_TEST ON CACHE INTERNAL "Skip cereal portability tests.")
set(SKIP_PERFORMANCE_COMPARISON ON CACHE INTERNAL "Skip cereal performance comparison.")
add_subdirectory("${PROJECT_SOURCE_DIR}/vendor/cereal" "extern/cereal")
include_directories("${PROJECT_SOURCE_DIR}/vendor/cereal/include")
