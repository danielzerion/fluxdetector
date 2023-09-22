cmake_minimum_required(VERSION 3.25)

project(fluxdetector)
set(EXE fluxdetector)
set(CMAKE_CXX_STANDARD 20)

# ----- Generate compile_commands.json, unless otherwise instructed on CLI --
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
# ----- Ensure that standard C++ headers are found by clangd
if(CMAKE_EXPORT_COMPILE_COMMANDS)
  set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

find_package(Nain4 REQUIRED)

# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)

# Add the executable, and link it to the {Geant,nain}4 libraries
add_executable(${EXE} fluxdetector.cc ${PROJECT_SOURCE_DIR}/materials.cc)
target_link_libraries(${EXE} PRIVATE Nain4 ${Geant4_LIBRARIES})

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${EXE}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
