cmake_minimum_required(VERSION 3.27)
project(Lab2)

set(CMAKE_CXX_STANDARD 17)

add_definitions(-O3 -Os)

add_executable(Lab2 main.cpp
        task1.h
        task2.h)
