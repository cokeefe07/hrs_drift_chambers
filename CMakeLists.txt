cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(drift_chamber_sims CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)


if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(GSL REQUIRED)
find_package(Garfield REQUIRED)


find_package(OpenMP)


function(add_sim_target name src)
    if(TARGET ${name})
        message(WARNING "${name}: target already defined - ignoring duplicate add_sim_target()")
        return()
    endif()

    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
        message(STATUS "${name}: ${src} not found - skipping this target")
        return()
    endif()

    add_executable(${name} ${src})

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${name} PRIVATE Garfield::Garfield OpenMP::OpenMP_CXX)
        message(STATUS "${name}: OpenMP found - building multithreaded")
    else()
        target_link_libraries(${name} PRIVATE Garfield::Garfield)
        message(STATUS "${name}: OpenMP NOT found - building single-threaded")
    endif()
endfunction()

add_sim_target(strip_cage strip_cage.cc)
add_sim_target(field_cage field_cage.cc)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/run_sweep.sh")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/run_sweep.sh"
         DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
         FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                          GROUP_READ GROUP_EXECUTE
                          WORLD_READ WORLD_EXECUTE)
endif()

