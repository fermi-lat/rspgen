include(UseSWIG)

add_library(
  rspgen STATIC
  src/CircularWindow.cxx
  src/CosineBinner.cxx
  src/GrbResponse.cxx
  src/IResponse.cxx
  src/PointResponse.cxx
  src/RspGenApp.cxx
  src/SpaceCraftCalculator.cxx
)
target_include_directories(
  rspgen
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:>
)

target_link_libraries(
  rspgen
  PRIVATE evtbin irfInterface st_stream tip
  PUBLIC astro st_facilities dataSubselector CLHEP::CLHEP irfLoader
)

add_executable(gtrspgen src/gtrspgen/gtrspgen.cxx)
target_link_libraries(gtrspgen PRIVATE rspgen st_app)
target_include_directories(
  gtrspgen
  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  PUBLIC $<INSTALL_INTERFACE:>
)

###### Tests ######
add_executable(test_rspgen src/test/test_rspgen.cxx)
target_link_libraries(
  test_rspgen
  PRIVATE rspgen astro CLHEP::CLHEP irfLoader dataSubselector evtbin irfInterface
  st_app st_facilities st_stream tip
)
target_include_directories(
  test_rspgen
  PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src> $<INSTALL_INTERFACE:>
)

###############################################################
# Installation
###############################################################
install(DIRECTORY rspgen DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_DATADIR})
install(DIRECTORY pfiles/ DESTINATION ${FERMI_INSTALL_PFILESDIR})

install(
  TARGETS rspgen gtrspgen test_rspgen
  # EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
