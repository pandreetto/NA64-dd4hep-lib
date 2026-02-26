cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(NA64-dd4hep-plugins)

find_package(Geant4 REQUIRED)
find_package(DD4hep REQUIRED COMPONENTS DDG4)

include_directories(./Core/include
                    /usr/include/DMG4
                    /usr/include/DMG4/DMParticles
                    /usr/include/DMG4/DMProcesses)

add_subdirectory(Core/src)
add_subdirectory(Geometry)

