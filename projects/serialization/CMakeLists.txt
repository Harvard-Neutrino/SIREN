add_library(LI_serialization INTERFACE)
target_include_directories(LI_serialization INTERFACE
    ${CMAKE_SOURCE_DIR}/projects/serialization/public/
)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/serialization/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

