cmake_minimum_required(VERSION 2.0)
project(helloworld)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

LINK_LIBRARIES(boost_python)
INCLUDE_DIRECTORIES("/usr/include/python2.6")
LINK_DIRECTORIES("/usr/lib/python2.6") # Not needed?

ADD_LIBRARY(helloc SHARED hello.cpp)
SET_TARGET_PROPERTIES( helloc PROPERTIES PREFIX "")



