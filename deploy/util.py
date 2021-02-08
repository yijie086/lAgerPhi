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
Utility functions for this container

Authors:
    - Sylvester Joosten <sjoosten@anl.gov>
'''

import os

class InvalidArgumentError(Exception):
    pass

def smart_mkdir(dir):
    '''functions as mkdir -p, with a write-check.
    
    Raises an exception if the directory is not writeable.
    '''
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except Exception as e:
            print('ERROR: unable to create directory', dir)
            raise e
    if not os.access(dir, os.W_OK):
        print('ERROR: We do not have the write privileges to', dir)
        raise InvalidArgumentError()

def project_version():
    '''Return the project version based on the current git branch/tag.'''
    ## Shell command to get the current project version
    version_cmd = 'cat VERSION'
    ## Strip will remove the leading \n character
    return os.popen(version_cmd).read().strip()
