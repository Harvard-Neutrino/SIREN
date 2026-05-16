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

# Apply local patches to the vendored photospline submodule before including
# it. Upstream photospline @ c6fb3ea (the pin used by SIREN) declares
# bsplvb_simple() and bspline_deriv_nonzero() in the public header but does
# not compile implementations for them -- the inline template code in
# detail/bspline_eval.h calls them, causing undefined-symbol link errors when
# libSIREN is loaded at runtime. Restore the implementations from an earlier
# photospline revision. Each patch is idempotent (--forward skips if already
# applied).
set(_PHOTOSPLINE_PATCH_DIR "${PROJECT_SOURCE_DIR}/cmake/photospline_patches")
if(EXISTS "${_PHOTOSPLINE_PATCH_DIR}")
    file(GLOB _PHOTOSPLINE_PATCHES "${_PHOTOSPLINE_PATCH_DIR}/*.patch")
    list(SORT _PHOTOSPLINE_PATCHES)
    foreach(_patch IN LISTS _PHOTOSPLINE_PATCHES)
        message(STATUS "Applying photospline patch: ${_patch}")
        execute_process(
            COMMAND patch -p1 --forward --silent -i "${_patch}"
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/vendor/photospline"
            RESULT_VARIABLE _patch_result
            OUTPUT_QUIET ERROR_QUIET
        )
        # patch returns 1 if already applied (which is fine); only error on 2+
        if(_patch_result GREATER 1)
            message(FATAL_ERROR "Failed to apply photospline patch ${_patch} (exit ${_patch_result})")
        endif()
    endforeach()
endif()
unset(_PHOTOSPLINE_PATCH_DIR)
unset(_PHOTOSPLINE_PATCHES)

# Override CMAKE_POLICY_VERSION_MINIMUM before adding subdirectory
set(TEMP_CMAKE_POLICY_VERSION_MINIMUM ${CMAKE_POLICY_VERSION_MINIMUM})
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/photospline)
set(CMAKE_POLICY_VERSION_MINIMUM ${TEMP_CMAKE_POLICY_VERSION_MINIMUM})
if(DEFINED SKBUILD)
    if(${CIBUILDWHEEL})
        message(STATUS "Setting photospline install lib dir to: ${CI_INSTALL_PREFIX}/lib")
        message(STATUS "Setting photospline install include dir to: ${CI_INSTALL_PREFIX}/include")
        install(TARGETS photospline
            LIBRARY DESTINATION "${CI_INSTALL_PREFIX}/lib"
            PUBLIC_HEADER DESTINATION "${CI_INSTALL_PREFIX}/include")
    else()
        set_target_properties(photospline PROPERTIES
                BUILD_WITH_INSTALL_RPATH FALSE
                LINK_FLAGS "-Wl,-rpath,$ORIGIN")
        install(TARGETS photospline
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren.libs
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    endif()
    #set_target_properties(spglam PROPERTIES
    #        BUILD_WITH_INSTALL_RPATH FALSE
    #        LINK_FLAGS "-Wl,-rpath,$ORIGIN")
    #install(TARGETS spglam
    #    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren.libs
    #    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
else()
    #install(TARGETS photospline DESTINATION ${CMAKE_INSTALL_LIBDIR})
    #install(TARGETS spglam DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
