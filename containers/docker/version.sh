#!/bin/bash

grep "liege VERSION" ../../CMakeLists.txt | sed 's/project (liege VERSION //' | sed 's/ LANGUAGES CXX )//'
