cmake_minimum_required(VERSION 3.20)
project(siren_prebuilt LANGUAGES NONE)

# Binaries were installed by the outer build, including their install RPATHs.
# Use the same backend/metadata as a source build to retain native wheel tags.
install(DIRECTORY siren DESTINATION . COMPONENT PythonWheel)
install(DIRECTORY siren.libs DESTINATION . COMPONENT PythonWheel OPTIONAL)
