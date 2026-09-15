# Keep CMake's install-time RPATH adjustment, then remove version-name aliases
# before the wheel backend dereferences them into duplicate binary payloads.
function(siren_install_wheel_libraries)
    set_target_properties(${ARGN} PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        INSTALL_RPATH "${SIREN_RPATH_ORIGIN}"
        INSTALL_RPATH_USE_LINK_PATH FALSE)
    if(APPLE)
        set_target_properties(${ARGN} PROPERTIES
            INSTALL_NAME_DIR "@rpath"
            BUILD_WITH_INSTALL_NAME_DIR TRUE)
    endif()
    install(TARGETS ${ARGN}
        LIBRARY DESTINATION siren.libs COMPONENT PythonWheel EXCLUDE_FROM_ALL
            NAMELINK_SKIP
        RUNTIME DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
    if(UNIX)
        foreach(library IN LISTS ARGN)
            install(CODE "
                set(wheel_directory \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/siren.libs\")
                set(real_name \"$<TARGET_FILE_NAME:${library}>\")
                set(soname \"$<TARGET_SONAME_FILE_NAME:${library}>\")
                if(NOT real_name STREQUAL soname)
                    file(REMOVE \"\${wheel_directory}/\${soname}\")
                    file(RENAME \"\${wheel_directory}/\${real_name}\"
                                \"\${wheel_directory}/\${soname}\")
                endif()
            " COMPONENT PythonWheel EXCLUDE_FROM_ALL)
        endforeach()
    endif()
endfunction()
