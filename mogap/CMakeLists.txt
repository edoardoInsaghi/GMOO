cmake_minimum_required(VERSION 3.10)

project(AirportGAP VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(MAIN_SOURCES
    main.cpp
    mogap.cpp
    chromosome.cpp
)

set(TEST_SOURCES
    tests.cpp
    chromosome.cpp
)

add_executable(main ${MAIN_SOURCES})
add_executable(test ${TEST_SOURCES})

target_compile_options(main PRIVATE -Wall -Wextra -pedantic -funroll-loops -march=native -O3 -mtune=native)
target_compile_options(test PRIVATE -Wall -Wextra -pedantic -funroll-loops -march=native -O3 -mtune=native)
