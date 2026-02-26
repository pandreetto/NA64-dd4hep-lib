include(${Geant4_USE_FILE})

add_library(na64-dd4hep-plugins Geant4DarkMatterPhysics.cpp)

target_link_libraries(na64-dd4hep-plugins ${Geant4_LIBRARIES}
                                          DD4hep::DDG4)

install(TARGETS na64-dd4hep-plugins DESTINATION lib64)

