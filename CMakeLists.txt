cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

set(CMAKE_CXX_STANDARD 17)

project(sbnnusyst VERSION 1.00.00 LANGUAGES CXX)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE DebWithRelInfo)
endif()

find_package(nusystematics 2.00.01 REQUIRED)
find_package(sbnanaobj 09.20.06.03 REQUIRED)
find_package(GENIE3 REQUIRED)
if(NOT TARGET GENIE3::All)
  message(FATAL_ERROR "Expected find_package(GENIE3 REQUIRED) call to set up target GENIE3::All.")
endif()

## Check if we have fhiclcpp from ups
find_package(fhiclcpp QUIET)
if(NOT fhiclcpp_FOUND)
  find_package(fhiclcppstandalone)
endif()

###### Compiler set up
add_library(sbnnusyst_dependencies INTERFACE)
target_link_libraries(sbnnusyst_dependencies INTERFACE fhiclcpp::fhiclcpp)
target_include_directories(sbnnusyst_dependencies INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include> )
set_target_properties(sbnnusyst_dependencies PROPERTIES EXPORT_NAME commondeps)
add_library(sbnnusyst::commondeps ALIAS sbnnusyst_dependencies)

install(TARGETS sbnnusyst_dependencies EXPORT sbnnusyst-targets)

add_subdirectory(src/sbnnusyst)

add_subdirectory(fcls)

install(EXPORT sbnnusyst-targets
        NAMESPACE sbnnusyst::
        DESTINATION lib/cmake/sbnnusyst)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/setup.sbnnusyst.sh.in 
  ${PROJECT_BINARY_DIR}/setup.sbnnusyst.sh @ONLY)
install(PROGRAMS ${PROJECT_BINARY_DIR}/setup.sbnnusyst.sh DESTINATION bin)
