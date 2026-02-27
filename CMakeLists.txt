cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(NA64-dd4hep-plugins)

option(INSTALL_COMPACT_FILES "Copy compact files to install area" OFF)

find_package(Geant4 REQUIRED)
find_package(DD4hep REQUIRED COMPONENTS DDG4)

include_directories(./Core/include
                    /usr/include/DMG4
                    /usr/include/DMG4/DMParticles
                    /usr/include/DMG4/DMProcesses)

add_subdirectory(Core/src)
add_subdirectory(Geometry)

if(INSTALL_COMPACT_FILES)
  file(GLOB compat_files "compat/*.xml")
  install(FILES ${compat_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/NA64-dd4hep-geo )
endif()

