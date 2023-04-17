execute_process(
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk
    COMMAND touch aclocal.m4 configure
    COMMAND touch Makefile.am Makefile.in aclocal.m4 m4/* rk.pc.in rk/Makefile.am rk/Makefile.in
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    ExternalProject_Add(
      rk
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env CXXFLAGS=-stdlib=libc++\ -fPIC ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk/configure --prefix=${PROJECT_BINARY_DIR}/extern/rk shrext_cmds="${CMAKE_SHARED_LIBRARY_SUFFIX}"
      PREFIX ${PROJECT_BINARY_DIR}/extern/rk
      BUILD_COMMAND make
      BUILD_IN_SOURCE 0
    )
else()
    ExternalProject_Add(
      rk
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env CXXFLAGS=-fPIC ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk/configure --prefix=${PROJECT_BINARY_DIR}/extern/rk shrext_cmds="${CMAKE_SHARED_LIBRARY_SUFFIX}"
      PREFIX ${PROJECT_BINARY_DIR}/extern/rk
      BUILD_COMMAND make
      BUILD_IN_SOURCE 0
    )
endif()

if(NOT SKBUILD)
install(DIRECTORY "${PROJECT_BINARY_DIR}/extern/rk/include/rk"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.hh")

install(DIRECTORY "${PROJECT_BINARY_DIR}/extern/rk/include/rk"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.icc")
endif()
