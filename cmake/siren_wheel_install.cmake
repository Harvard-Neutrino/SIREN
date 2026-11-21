# Keep CMake's install-time RPATH adjustment, then remove version-name aliases
# before the wheel backend dereferences them into duplicate binary payloads.
set(SIREN_WHEEL_LIBRARY_DIR "siren.libs")

# Directories of the external libraries the wheel links but does not carry.
function(siren_wheel_external_rpath out)
    set(dirs "")
    foreach(target IN ITEMS MARLEY HepMC3::HepMC3)
        if(NOT TARGET ${target})
            continue()
        endif()
        get_target_property(type ${target} TYPE)
        if(type STREQUAL "INTERFACE_LIBRARY")
            continue()
        endif()
        set(properties IMPORTED_LOCATION)
        get_target_property(configs ${target} IMPORTED_CONFIGURATIONS)
        if(configs)
            foreach(config IN LISTS configs)
                string(TOUPPER "${config}" config)
                list(APPEND properties IMPORTED_LOCATION_${config})
            endforeach()
        endif()
        foreach(property IN LISTS properties)
            get_target_property(location ${target} ${property})
            if(location)
                get_filename_component(dir "${location}" DIRECTORY)
                list(APPEND dirs "${dir}")
            endif()
        endforeach()
    endforeach()
    if(dirs)
        list(REMOVE_DUPLICATES dirs)
    endif()
    set(${out} "${dirs}" PARENT_SCOPE)
endfunction()

function(siren_install_wheel_libraries)
    siren_wheel_external_rpath(external_rpath)
    foreach(library IN LISTS ARGN)
        # These may also have native install rules (photospline/spglam). Keep
        # their dependency search paths, especially external @rpath libraries.
        # Wheel repair resolves those dependencies and removes external RPATHs.
        get_target_property(native_rpath ${library} INSTALL_RPATH)
        if(NOT native_rpath)
            set(native_rpath "")
        endif()
        # Prefer the installed sibling over a native-prefix copy during wheel
        # repair too, or the repair tool can graft a second spglam/photospline.
        set(rpath "${SIREN_RPATH_ORIGIN}" ${native_rpath} ${external_rpath})
        set_property(TARGET ${library} PROPERTY INSTALL_RPATH "${rpath}")
        set_target_properties(${library} PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)
    endforeach()
    install(TARGETS ${ARGN}
        LIBRARY DESTINATION ${SIREN_WHEEL_LIBRARY_DIR} COMPONENT PythonWheel EXCLUDE_FROM_ALL
            NAMELINK_SKIP
        RUNTIME DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
    if(UNIX)
        foreach(library IN LISTS ARGN)
            install(CODE "
                set(wheel_directory \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${SIREN_WHEEL_LIBRARY_DIR}\")
                set(real_name \"$<TARGET_FILE_NAME:${library}>\")
                set(soname \"$<TARGET_SONAME_FILE_NAME:${library}>\")
                if(NOT real_name STREQUAL soname)
                    file(REMOVE \"\${wheel_directory}/\${soname}\")
                    file(RENAME \"\${wheel_directory}/\${real_name}\"
                                \"\${wheel_directory}/\${soname}\")
                endif()
            " COMPONENT PythonWheel EXCLUDE_FROM_ALL)
            if(APPLE)
                install(CODE "
                    execute_process(
                        COMMAND \"${Python_EXECUTABLE}\"
                            \"${CMAKE_CURRENT_FUNCTION_LIST_DIR}/wheel_rpath.py\"
                            \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${SIREN_WHEEL_LIBRARY_DIR}/$<TARGET_SONAME_FILE_NAME:${library}>\"
                        RESULT_VARIABLE rpath_result)
                    if(NOT rpath_result EQUAL 0)
                        message(FATAL_ERROR \"Could not order wheel dependency RPATHs\")
                    endif()
                " COMPONENT PythonWheel EXCLUDE_FROM_ALL)
            endif()
        endforeach()
    endif()
endfunction()

function(siren_install_wheel_modules)
    siren_wheel_external_rpath(external_rpath)
    set(rpath "${SIREN_RPATH_ORIGIN}/../${SIREN_WHEEL_LIBRARY_DIR}" ${external_rpath})
    set_target_properties(${ARGN} PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        INSTALL_RPATH "${rpath}"
        INSTALL_RPATH_USE_LINK_PATH FALSE)
    install(TARGETS ${ARGN}
        LIBRARY DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL
        RUNTIME DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
    # Keep the library folder first in the modules' RPATHs.
    if(APPLE AND external_rpath)
        foreach(module IN LISTS ARGN)
            install(CODE "
                execute_process(
                    COMMAND \"${Python_EXECUTABLE}\"
                        \"${CMAKE_CURRENT_FUNCTION_LIST_DIR}/wheel_rpath.py\"
                        \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/siren/$<TARGET_FILE_NAME:${module}>\"
                        \"${SIREN_RPATH_ORIGIN}/../${SIREN_WHEEL_LIBRARY_DIR}\"
                    RESULT_VARIABLE rpath_result)
                if(NOT rpath_result EQUAL 0)
                    message(FATAL_ERROR \"Could not order wheel module RPATHs\")
                endif()
            " COMPONENT PythonWheel EXCLUDE_FROM_ALL)
        endforeach()
    endif()
endfunction()
