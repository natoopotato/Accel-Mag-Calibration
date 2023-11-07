cmake_minimum_required(VERSION 3.26)
project(accel_calc)
set(CMAKE_C_STANDARD 11)


include_directories("./mathlib")


file(GLOB_RECURSE SOURCES "mathlib/*.*" "Core/*.*")

add_executable(accel_calc ${SOURCES})
if (NOT WIN32)
    target_link_libraries(accel_calc m)
endif()
