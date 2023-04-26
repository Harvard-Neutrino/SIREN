add_library(LI_serialization INTERFACE)
target_include_directories(LI_serialization INTERFACE
    ${PROJECT_SOURCE_DIR}/projects/serialization/public/
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/serialization/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

