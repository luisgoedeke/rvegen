cmake_minimum_required(VERSION 3.5)

project(rvegen LANGUAGES CXX)

set(RVEGEN_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/rvegen)
set(RVEGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/rvegen/include)

# Versionning
# ===========

file(STRINGS "${RVEGEN_BASE_DIR}/rvegen_config.h" rvegen_version_defines
     REGEX "#define RVEGEN_VERSION_(MAJOR|MINOR|PATCH)")
foreach(ver ${RVEGEN_version_defines})
    if(ver MATCHES "#define RVEGEN_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(RVEGEN_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()
set(${PROJECT_NAME}_VERSION
    ${RVE_VERSION_MAJOR}.${RVE_VERSION_MINOR}.${RVE_VERSION_PATCH})
message(STATUS "Building revgen v${${PROJECT_NAME}_VERSION}")


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(RVEGEN_HEADER
    ${RVEGEN_INCLUDE_DIR}/shape_base.h
    ${RVEGEN_INCLUDE_DIR}/cylinder.h
    ${RVEGEN_INCLUDE_DIR}/circle.h
    ${RVEGEN_INCLUDE_DIR}/ellipse.h
    ${RVEGEN_INCLUDE_DIR}/rve_shape_input.h
    ${RVEGEN_INCLUDE_DIR}/write_gmsh_geo.h
    ${RVEGEN_INCLUDE_DIR}/check_distance.h
    ${RVEGEN_INCLUDE_DIR}/rve_generator.h
)

# add the library that runs
if(${CMAKE_VERSION} VERSION_LESS "3.19.0")
    add_library(${PROJECT_NAME} INTERFACE)
else()
add_library(${PROJECT_NAME} INTERFACE
    ${RVEGEN_BASE_DIR}/rvegen.h
    ${RVEGEN_BASE_DIR}/rvegen_config.h
    ${RVEGEN_HEADER}
    )
endif()


# add alias so the project can be uses with add_subdirectory
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})


include(GNUInstallDirs)
#include(CTest)


if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(IS_TOPLEVEL_PROJECT TRUE)
else()
  set(IS_TOPLEVEL_PROJECT FALSE)
endif()


OPTION(RVEGEN_INSTALL_LIBRARY
       "Enable installing of rvegen library into default locations"
       ${IS_TOPLEVEL_PROJECT})
OPTION(BUILD_TESTS "rvegen test suite" OFF)
OPTION(BUILD_BENCHMARK "rvegen benchmark" OFF)
OPTION(DOWNLOAD_GTEST "build gtest from downloaded sources" OFF)
OPTION(DOWNLOAD_GBENCHMARK "download google benchmark and build from source" OFF)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

if(BUILD_BENCHMARK)
    add_subdirectory(benchmark)
endif()


# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUIcd build && ctestnstallDirs)
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)


if(RVEGEN_INSTALL_LIBRARY)
  # locations are provided by GNUInstallDirs
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    "${PROJECT_NAME}ConfigVersion.cmake"
    VERSION {${PROJECT_NAME}_VERSION}
    COMPATIBILITY SameMajorVersion)

  configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(
    EXPORT ${PROJECT_NAME}_Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
          DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/rvegen
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()


set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}")

include(CPack)
