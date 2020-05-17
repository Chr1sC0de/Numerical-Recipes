cmake_minimum_required(VERSION 3.16.1)

project(NumericalRecipes)

file(GLOB test_files "tests/*.cpp")
file(GLOB pybindfiles "source/*.cpp")

find_package(OpenMP)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# the test files to be built

foreach(file ${test_files})

    get_filename_component(filename_component ${file} NAME_WE)

    add_executable(${filename_component} ${file})

    target_include_directories(${filename_component} PUBLIC "include")

    # target_link_libraries(${filename_component} PRIVATE ${VTK_LIBRARIES})

    target_link_libraries(${filename_component} ${Boost_LIBRARIES})

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${filename_component} PUBLIC OpenMP::OpenMP_CXX)
    endif()

endforeach()

# the source files to be built

add_subdirectory(pybind11)

foreach(file ${pybindfiles})

    get_filename_component(filename_component ${file} NAME_WE)

    pybind11_add_module(${filename_component} ${SOURCES} ${file})

    target_include_directories(${filename_component} PUBLIC "include")

    target_link_libraries(${filename_component} ${Boost_LIBRARIES})

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${filename_component} PUBLIC OpenMP::OpenMP_CXX)
    endif()

endforeach()