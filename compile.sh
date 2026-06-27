#!/bin/bash
python config.py $1 b3s23-a5
echo "Compiling..."
g++ spark.cpp -o spark -O3 -Ofast -flto -march=native -s -Os
echo "Compilation successful."
