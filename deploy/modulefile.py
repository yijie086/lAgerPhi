#!/usr/bin/env python3

## lAger: General Purpose l/A-event Generator
## Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
## 
## This file is part of lAger.
## 
## lAger is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Shoftware Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## lAger is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with lAger.  If not, see <https://www.gnu.org/licenses/>.
## 

'''
Install modulefile for this container.

Authors:
    - Sylvester Joosten <sjoosten@anl.gov>
'''

import os

## Generic module file
_MODULEFILE='''#%Module1.0#####################################################################
##
## for {name} {version}
##
proc ModulesHelp {{ }} {{
    puts stderr "This module sets up the environment for the {name} container"
}}
module-whatis "{name} {version}"

# For Tcl script use only
set version 4.1.4

prepend-path    PATH    {bindir}
'''

def make_modulefile(project, version, moduledir, bindir):
    '''Configure and install a modulefile for this project.

    Arguments:
        - project: project name
        - version: project version
        - moduledir: root modulefile directory
        - bindir: where executables for this project are located
    '''

    ## create our modulefile
    content = _MODULEFILE.format(name=project, version=version, bindir=bindir)
    fname = '{}/{}'.format(moduledir, version)
    print(' - creating', fname)
    with open(fname, 'w') as file:
        file.write(content)
