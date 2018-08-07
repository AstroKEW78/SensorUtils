cmake_minimum_required(VERSION 3.10)
project(sensorutils VERSION 0.1.1 DESCRIPTION "SensorUtils library")

include(cmake/gtest.cmake)
include(GNUInstallDirs)
include(GoogleTest)

set(CMAKE_CXX_STANDARD 11)

set(COVERAGE OFF CACHE BOOL "Coverage")

add_library(sensorutils SHARED
            src/SensorUtils.cpp)

if(COVERAGE)
    target_compile_options(sensorutils PRIVATE --coverage -O0)
    target_link_libraries(sensorutils PRIVATE --coverage -O0)
endif()

set_target_properties(sensorutils PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER ${PROJECT_SOURCE_DIR}/include/sensorutils/SensorUtils.h
)

find_package(Threads)

# This is the linear algebra library we are using
find_package(Armadillo REQUIRED)

target_include_directories(sensorutils
                           PUBLIC
                           include/sensorutils/
                           include/
                           ${ARMADILLO_INCLUDE_DIRS}
)


target_link_libraries(sensorutils PUBLIC
                      gtest ${CMAKE_THREAD_LIBS_INIT}
                      ${ARMADILLO_LIBRARIES})

install(TARGETS sensorutils
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sensorutils)



include(CTest)

enable_testing()

add_subdirectory(tests)
