project(hartreefock)
cmake_minimum_required(VERSION 2.6)
find_package(Qt4 REQUIRED)
find_package(ROOT REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

include_directories(${QT_INCLUDES}   ${CMAKE_CURRENT_BINARY_DIR} ${ROOT_INCLUDE_DIR})

set(hartreefock_SRCS hartreefock.cpp main.cpp)
qt4_automoc(${hartreefock_SRCS})
add_executable(hartreefock ${hartreefock_SRCS})
target_link_libraries(hartreefock ${QT_QTCORE_LIBRARY} ${ROOT_LIBRARIES})

install(TARGETS hartreefock RUNTIME DESTINATION bin)
