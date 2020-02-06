#!/usr/bin/env python3

'''
Deploy the singularity container built by the CI for this version of the software.

The current version is determined from the currently loaded git branch or tag,
unless it is explicitly set by the command line.

Author: Sylvester Joosten <sjoosten@anl.gov>
'''

import os
import argparse
import urllib.request

## Gitlab group and project/program name. 
GROUP_NAME='monte_carlo'
PROJECT_NAME='liege'
PROGRAMS = ['liege', 'root']

## URL for the current container (git tag will be filled in by the script)
CONTAINER_URL = r'https://eicweb.phy.anl.gov/{0}/{1}/-/jobs/artifacts/{2}/raw/build/{1}.sif?job={1}_singularity'

## Singularity bind directive
BIND_DIRECTIVE= '-B {0}:{0}'

## generic launcher bash script to launch the application
LAUNCHER_SCRIPT='''
#!/bin/bash

piped_args=
if [ -p /dev/stdin ]; then
  # If we want to read the input line by line
  while IFS= read line; do
    if [ -z "$piped_args" ]; then
       piped_args="$line"
    else 
       piped_args="$piped_args\n$line"
    fi
  done
fi

if [ $piped_args ]  ; then
  echo -e $piped_args | singularity exec {2} {1} {0} $@
else
  singularity exec {2} {1} {0} $@
fi
'''

class InvalidArgumentError(Exception):
    pass


def project_version():
    ## Shell command to get the current git version
    git_version_cmd = 'git symbolic-ref -q --short HEAD || git describe --tags --exact-match'
    ## Strip will remove the leading \n character
    return os.popen(git_version_cmd).read().strip()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'prefix',
            help='Install prefix. This is where the container will be deployed.')
    parser.add_argument(
            '-v', '--version',
            dest='version',
            default=project_version(),
            help='(opt.) project version. Default: current git branch/tag.')
    parser.add_argument(
            '-f', '--force',
            action='store_true',
            help='Force-overwrite already downloaded container',
            default=False)
    parser.add_argument(
            '-b', '--bind-path',
            dest='bind_paths',
            action='append',
            help='(opt.) extra bind paths for singularity.')

    args = parser.parse_args()

    print('Deploying', PROJECT_NAME, 'version', args.version)

    ## Check if our bind paths are valid
    bind_directive = ''
    if args.bind_paths and len(args.bind_paths):
        print('Singularity bind paths:')
        for path in args.bind_paths:
            print(' -', path)
            if not os.path.exists(path):
                print('ERROR: path', path, 'does not exist.')
                raise InvalidArgumentError()
        bind_directive = ' '.join([BIND_DIRECTIVE.format(path) for path in args.bind_paths])

    ## Create our install prefix if needed and ensure it is writable
    args.prefix = os.path.abspath(args.prefix)
    print('Install prefix:', args.prefix)
    print('Creating install prefix if needed...')
    libdir = '{}/lib'.format(args.prefix)
    bindir = '{}/bin'.format(args.prefix)
    for dir in [libdir, bindir]:
        print(' -', dir)
        if not os.path.exists(dir):
            try:
                os.makedirs(dir)
            except Exception as e:
                print('ERROR: unable to create directory', dir)
                raise e
        if not os.access(dir, os.W_OK):
            print('ERROR: We do not have the write privileges to', dir)
            raise InvalidArgumentError()

    ## At this point we know we can write to our desired prefix and that we have a set of
    ## valid bind paths

    ## Get the container
    ## We want to slightly modify our version specifier: if it leads with a 'v' drop the v
    version = '{}'.format(args.version)
    if version[0] is 'v':
        version = version[1:]
    if args.version[0].isdigit():
        args.version = 'v{}'.format(args.version)
    container = '{}/{}.sif.{}'.format(libdir, PROJECT_NAME, version)
    if not os.path.exists(container) or args.force:
        url = CONTAINER_URL.format(GROUP_NAME, PROJECT_NAME, args.version)
        print('Downloading container from:', url)
        print('Destination:', container)
        urllib.request.urlretrieve(url, container)
    else:
        print('WARNING: Container found at', container, 'run with -f to force a re-download')

    ## configure the application launchers
    print('Configuring applications launchers: ')
    for app in PROGRAMS:
        fname = '{}/{}'.format(bindir, app)
        print(' - creating', fname)
        with open(fname, 'w') as file:
            script = LAUNCHER_SCRIPT.format(app, container, bind_directive)
            file.write(script)
        os.system('chmod +x {}'.format(fname))

    print('Container deployment successful!')
