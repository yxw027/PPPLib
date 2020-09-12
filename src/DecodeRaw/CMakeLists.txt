PROJECT(DecodeRawLib)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(SET CMP0015 NEW)

file(GLOB SRC_FILES "*.cc")
add_library(DecodeRawLib ${SRC_FILES})
target_link_libraries(DecodeRawLib CmnLib)