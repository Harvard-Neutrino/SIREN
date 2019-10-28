cmake_minimum_required(VERSION 3.9)
project(lepton_injector VERSION 1.0.0 DESCRIPTION "lepton injector library")
include(GNUInstallDirs)

add_library(lepton_injector SHARED
    private/LeptonInjector/Controller.cxx
    private/LeptonInjector/Coordinates.cxx
    private/LeptonInjector/DataWriter.cxx
    private/LeptonInjector/EventProps.cxx
    private/LeptonInjector/LeptonInjector.cxx
    private/LeptonInjector/Particle.cxx
    private/LeptonInjector/Random.cxx
    private/earthmodel-service/EarthModelCalculator.cxx
    private/earthmodel-service/EarthModelService.cxx
    private/phys-services/LICrossSection.cxx
)
configure_file(leptoninjector.pc.in leptoninjector.pc @ONLY)
target_include_directories(lepton_injector PRIVATE .)
install(TARGETS lepton_injector
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION $CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/lepton_injector.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


