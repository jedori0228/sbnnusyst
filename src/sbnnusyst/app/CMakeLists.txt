LIST(APPEND TARGETS_TO_BUILD 
UpdateReweight
ReadCAF
)

foreach(targ ${TARGETS_TO_BUILD})
  add_executable(${targ} ${targ}.cxx)
  target_link_libraries(${targ} 
    nusyst::systproviders
    ${ROOT_LIBRARIES}
    ROOT::Geom
    ROOT::MathMore
    sbnanaobj::StandardRecord
    sbnanaobj::StandardRecordProxy
    sbnanaobj::StandardRecordFlat
    sbnnusyst::commondeps
    sbnnusyst::interface
    sbnnusyst::utility
  )
  target_include_directories(${targ}
    PRIVATE
    $ENV{SRPROXY_INC}
  )

endforeach()

install(TARGETS ${TARGETS_TO_BUILD} DESTINATION bin/)
