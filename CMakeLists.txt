cmake_minimum_required(VERSION 3.10)
project(SeismicInversionTool)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -march=native")

# Find FFTW via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

include_directories(include)

# Add the executable first
add_executable(seismic_inversion
  src/main.cpp
  src/forward_modelling.cpp
)

# Now link libraries to the target
target_link_libraries(seismic_inversion
  PkgConfig::FFTW
)

