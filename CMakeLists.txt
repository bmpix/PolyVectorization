cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_VERBOSE_MAKEFILE ON)
project(polyvector_thing)
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
add_executable(polyvector_thing)
##### QT #####
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui)
if (WITH_QT)
    target_compile_definitions(polyvector_thing PRIVATE WITH_QT=${WITH_QT})
endif()

if (WITH_GUI)
    target_compile_definitions(polyvector_thing PRIVATE WITH_GUI=${WITH_GUI})
endif()

if ((Qt5_FOUND) AND (WITH_QT OR WITH_GUI))
 message(STATUS "Qt found")
 target_link_libraries(polyvector_thing PUBLIC Qt5::Core Qt5::Widgets Qt5::Gui)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

target_include_directories(polyvector_thing PUBLIC src)
#target_compile_options()

##### OpenMP #####
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(polyvector_thing PUBLIC OpenMP::OpenMP_CXX)
endif()


##### OpenCV #####
find_package(OpenCV REQUIRED)
target_link_libraries(polyvector_thing PUBLIC ${OpenCV_LIBS})

##### Eigen #####
find_package (Eigen3 REQUIRED NO_MODULE)
add_definitions     ( ${EIGEN3_DEFINITIONS} )
include_directories ( ${EIGEN3_INCLUDE_DIRS} )

##### Boost #####
find_package(Boost 1.48 REQUIRED)
message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
target_include_directories(polyvector_thing PUBLIC ${Boost_INCLUDE_DIR})
target_link_libraries(polyvector_thing PUBLIC ${Boost_LIBRARIES})

##### Sources #####
add_subdirectory(src)


