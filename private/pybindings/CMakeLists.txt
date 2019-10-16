if(BUILD_PYBINDINGS)
  i3_add_pybindings(shield
    shield.cxx
    USE_TOOLS boost python
    USE_PROJECTS icetray dataclasses shield
    )
endif(BUILD_PYBINDINGS)
