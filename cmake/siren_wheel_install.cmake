# Keep CMake's install-time RPATH adjustment, then remove version-name aliases
# before the wheel backend dereferences them into duplicate binary payloads.
set(SIREN_WHEEL_LIBRARY_DIR "siren.libs")

function(siren_install_wheel_libraries)
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
        set_property(TARGET ${library} PROPERTY
            INSTALL_RPATH "${SIREN_RPATH_ORIGIN};${native_rpath}")
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
    set_target_properties(${ARGN} PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        INSTALL_RPATH "${SIREN_RPATH_ORIGIN}/../${SIREN_WHEEL_LIBRARY_DIR}"
        INSTALL_RPATH_USE_LINK_PATH FALSE)
    install(TARGETS ${ARGN}
        LIBRARY DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL
        RUNTIME DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
endfunction()
