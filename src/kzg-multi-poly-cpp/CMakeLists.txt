cmake_minimum_required(VERSION 3.10)
project(KZG_CPP)

set(CMAKE_CXX_STANDARD 17)  
find_package(mcl REQUIRED)

include_directories(${MCL_INCLUDE_DIRS})

add_executable(kzg main.cpp kzg.cpp polynomial.cpp util.cpp)
target_link_libraries(kzg ${MCL_LIBRARIES})