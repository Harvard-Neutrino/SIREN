add_library(LI_serialization INTERFACE)
target_include_directories(LI_serialization INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/serialization/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/serialization/>
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/serialization/public/"
    EXPORT ${PROJECT_NAME}Config
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

