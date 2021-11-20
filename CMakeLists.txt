cmake_minimum_required(VERSION 3.15)
project(ToyEfficiency)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED COMPONENTS program_options regex)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO Tree)

set(CMAKE_CXX_FLAGS "-Wall")
include(${ROOT_USE_FILE})

include_directories(
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${ROOT_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
)

add_executable(toy_efficiency src/toy_efficiency.cc)
target_link_libraries(toy_efficiency ${Boost_LIBRARIES} ${ROOT_LIBRARIES})