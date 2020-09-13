PROJECT(test_main)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(SET CMP0015 NEW)

add_executable(test_SPP test_SPP.cc)
target_link_libraries(test_SPP PPPLib)

add_executable(test_PPP test_PPP.cc)
target_link_libraries(test_PPP PPPLib)

add_executable(test_PPK test_PPK.cc)
target_link_libraries(test_PPK PPPLib)

add_executable(test_LC_FUSION test_LC_FUSION.cc)
target_link_libraries(test_LC_FUSION PPPLib)

#add_executable(test_TC_FUSION test_TC_FUSION.cc)
#target_link_libraries(test_TC_FUSION PPPLib)

add_executable(FormatFiles FormatFiles.cc)
target_link_libraries(FormatFiles PPPLib)

add_executable(PPPLibMain PPPLibMain.cc)
target_link_libraries(PPPLibMain PPPLib)

