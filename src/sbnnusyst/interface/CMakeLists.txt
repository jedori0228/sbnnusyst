SET(IFCE_IMPLFILES
  WeightUpdater.cxx
)

SET(IFCE_HDRFILES
  WeightUpdater.h
)


add_library(sbnnusyst_interface SHARED ${IFCE_IMPLFILES})
add_library(sbnnusyst::interface ALIAS sbnnusyst_interface)

set_target_properties(sbnnusyst_interface PROPERTIES 
  PUBLIC_HEADER "${IFCE_HDRFILES}"
  EXPORT_NAME interface )
target_compile_definitions(sbnnusyst_interface PRIVATE $<$<BOOL:${SBNANAOBJ_KMULTISIM_CAP}>:SBNANAOBJ_KMULTISIM_CAP>)

target_link_libraries(sbnnusyst_interface
PUBLIC
  nusyst::systproviders
  ${ROOT_LIBRARIES}
  ROOT::Geom
  ROOT::MathMore
  sbnanaobj::StandardRecord
  sbnanaobj::StandardRecordProxy
  sbnanaobj::StandardRecordFlat
  sbnnusyst::commondeps
  sbnnusyst::utility
)
target_include_directories(sbnnusyst_interface
  PRIVATE
  $ENV{SRPROXY_INC}
)


#install(TARGETS ${TARGETS_TO_BUILD} DESTINATION interface/)

install(TARGETS sbnnusyst_interface
    EXPORT sbnnusyst-targets
    PUBLIC_HEADER DESTINATION include/sbnnusyst/interface COMPONENT Development)
