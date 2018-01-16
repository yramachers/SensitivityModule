# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(SensitivityModule)
find_package(Falaise REQUIRED)

# Build a dynamic library from our sources
add_library(SensitivityModule SHARED SensitivityModule.h SensitivityModule.cpp)

# Link it to the FalaiseModule library
# This ensures the correct compiler flags, include paths
# and linker flags are applied to our dynamic library.
# - Explicit link to ROOT Physics library which is not linked
#   by flreconstruct
target_link_libraries(SensitivityModule
  PUBLIC
    Falaise::FalaiseModule
    ${ROOT_Physics_LIBRARY}
    )

# Configure example pipeline script for use from the build dir
configure_file("SensitivityModuleExample.conf.in" "SensitivityModuleExample.conf" @ONLY)
