#!/bin/bash

grep "lager VERSION" ../../CMakeLists.txt | sed 's/project (lager VERSION //' | sed 's/ LANGUAGES CXX )//'
