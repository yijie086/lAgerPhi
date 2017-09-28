#!/usr/bin/python

import os, subprocess, sys

current_path = os.path.abspath(os.path.dirname(__file__))
cmake_build = os.path.abspath('@CMAKE_BINARY_DIR@') + '/'
cmake_install_lib = os.path.abspath('@INSTALL_LIB_DIR@') + '/'

## are we intree?
intree = (local_dir.find(cmake_build) is 0)

## get the correct setup_root.C script
setup_root = None
if intree:
    setup_root = cmake_build + "root/setup_root.cc"
else:
    setup_root = cmake_install_lib + "setup_root.cc"

if not os.path.exists(setup_root):
    raise IOError("%s: file not found: %s" % (__file__, setup_root))

if not len(sys.argv) is 2:
    raise IOError("%s: incorrect number of arguments")

subprocess.call("root -e %s -x %s" % (setup_root, argv[1])
