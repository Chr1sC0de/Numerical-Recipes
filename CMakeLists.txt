cmake_minimum_required(VERSION 3.16.1)

project(NumericalRecipes)

file(GLOB test_files "tests/*.cpp")

find_package(OpenMP)

foreach(file ${test_files})

    get_filename_component(filename_component ${file} NAME_WE)

    add_executable(${filename_component} ${file})

    target_include_directories(${filename_component} PUBLIC "include")

    # target_link_libraries(${filename_component} PRIVATE ${VTK_LIBRARIES})

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${filename_component} PUBLIC OpenMP::OpenMP_CXX)
    endif()

endforeach()