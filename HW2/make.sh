#!/usr/bin/env sh

CC=gcc-14 cmake -S . -B build

cmake --build build 

cp build/main .
