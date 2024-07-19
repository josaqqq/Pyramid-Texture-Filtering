cmake_minimum_required(VERSION 3.12)

project(PyramidTextureFiltering)
set(PROJECT_VERSION 1.0.0)

# Add stb_image
include_directories("external")

# Create an executable
add_executable(
    PyramidTextureFiltering
    src/main.cpp
)