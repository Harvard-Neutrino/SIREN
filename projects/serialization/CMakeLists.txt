add_library(SIREN_serialization INTERFACE)
target_include_directories(SIREN_serialization INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/serialization/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/serialization/>
)

target_link_libraries(SIREN_serialization
    INTERFACE
        siren_compile_options
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/serialization/public/"
    EXPORT ${PROJECT_NAME}Config
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

