
LIST (APPEND serialization_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/serialization/private/serialization.cxx
)

add_library(SIREN_serialization OBJECT ${serialization_SOURCES})
target_include_directories(SIREN_serialization PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/serialization/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/serialization/>
)

apply_siren_compile_options(SIREN_serialization)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/serialization/public/"
    EXPORT ${PROJECT_NAME}Config
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

