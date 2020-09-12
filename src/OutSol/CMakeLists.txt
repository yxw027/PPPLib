PROJECT(OutLib)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(SET CMP0015 NEW)

file(GLOB SRC_FILES "*.cc")
add_library(OutLib ${SRC_FILES})
target_include_directories(OutLib PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(OutLib CmnLib ${PYTHON_LIBRARIES})