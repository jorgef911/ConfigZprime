#!/bin/env bash
g++ test_root.cc `root-config --cflags` `root-config --libs`  -o test_root
