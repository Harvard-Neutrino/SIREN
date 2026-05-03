add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/crundec3 EXCLUDE_FROM_ALL)
set_property(TARGET crundec_static PROPERTY POSITION_INDEPENDENT_CODE ON)
