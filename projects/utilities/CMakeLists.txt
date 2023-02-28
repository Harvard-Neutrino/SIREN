
#core files
LIST (APPEND utilities_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/utilities/private/Interpolator.cxx
    ${CMAKE_SOURCE_DIR}/projects/utilities/private/Random.cxx
)
add_library(LI_utilities SHARED ${utilities_SOURCES})
target_include_directories(LI_utilities PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/utilities/public/
)

install(TARGETS LI_utilities
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/utilities/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

