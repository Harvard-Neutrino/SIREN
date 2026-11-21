# Run the documented pip command for the wheel the python_package target built.
# Script mode: invoked by the install_wheel target, never by cmake --install.
file(GLOB wheels "${WHEEL_DIR}/*.whl")
list(LENGTH wheels count)
if(NOT count EQUAL 1)
    message(FATAL_ERROR
        "Expected one wheel in ${WHEEL_DIR} (found ${count}); build the python_package target first")
endif()
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -m pip install --force-reinstall --no-deps ${wheels}
    COMMAND_ECHO STDOUT
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "pip did not install ${wheels} (exit ${result})")
endif()
