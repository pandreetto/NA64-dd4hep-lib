file(GLOB sources *_geo.cpp)

add_dd4hep_plugin(na64geom SHARED ${sources})

target_link_libraries(na64geom DD4hep::DDCore DD4hep::DDRec DD4hep::DDParsers)

install(TARGETS na64geom
        EXPORT ${PROJECT_NAME}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
