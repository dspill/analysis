#!/bin/bash
clear;
make clean
make check && check --offset 3 --step 7 --max 142
