cmake_minimum_required(VERSION 3.17)
project(algebra C)

set(CMAKE_C_STANDARD 11)

add_executable(algebra main.c svd.h stuff.h polar.h cholesky.h linear.h)