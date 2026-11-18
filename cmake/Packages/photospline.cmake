find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        if(NOT EXISTS "${PROJECT_SOURCE_DIR}/vendor/photospline/CMakeLists.txt")
            message(STATUS "Submodule update")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init ./
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/photospline
                            RESULT_VARIABLE GIT_SUBMOD_RESULT)
            if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
            endif()
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/vendor/photospline/CMakeLists.txt")
    message(FATAL_ERROR "The photospline submodule was not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

#add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/photospline EXCLUDE_FROM_ALL)
# Override CMAKE_POLICY_VERSION_MINIMUM before adding subdirectory
set(TEMP_CMAKE_POLICY_VERSION_MINIMUM ${CMAKE_POLICY_VERSION_MINIMUM})
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/photospline)
set(CMAKE_POLICY_VERSION_MINIMUM ${TEMP_CMAKE_POLICY_VERSION_MINIMUM})

# photospline's headers use uint32_t and uint64_t without including <cstdint>.
# The standard headers of GCC 15 no longer bring it in, so give it to
# photospline and to every C++ source compiled against its headers.
if(TARGET photospline AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(photospline PUBLIC
        "$<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:CXX>:SHELL:-include cstdint>>")
endif()
