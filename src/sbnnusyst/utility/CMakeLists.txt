SET(IFCE_IMPLFILES
  Utilities.cxx
)

SET(IFCE_HDRFILES
  Utilities.h
)


add_library(sbnnusyst_utility SHARED ${IFCE_IMPLFILES})
add_library(sbnnusyst::utility ALIAS sbnnusyst_utility)

set_target_properties(sbnnusyst_utility PROPERTIES 
  PUBLIC_HEADER "${IFCE_HDRFILES}"
  EXPORT_NAME utility )

target_link_libraries(sbnnusyst_utility
PUBLIC
  nusyst::systproviders
  ${ROOT_LIBRARIES}
  ROOT::Geom
  ROOT::MathMore
  sbnanaobj::StandardRecord
  sbnanaobj::StandardRecordProxy
  sbnanaobj::StandardRecordFlat
  sbnnusyst::commondeps
)
target_include_directories(sbnnusyst_utility
  PRIVATE
  $ENV{SRPROXY_INC}
)

install(TARGETS sbnnusyst_utility
    EXPORT sbnnusyst-targets
    PUBLIC_HEADER DESTINATION include/sbnnusyst/utility COMPONENT Development)
