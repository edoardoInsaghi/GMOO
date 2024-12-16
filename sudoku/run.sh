#!/bin/bash

g++-14 -fopenmp -std=c++17 -o main main.cpp -funroll-loops -O3 -mcpu=apple-m1 -mtune=native
./main
rm main
