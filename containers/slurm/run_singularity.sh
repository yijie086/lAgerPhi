#!/bin/bash

module load singularity
singularity exec -B /lcrc:/lcrc /lcrc/project/jlab/images/pcsim_bdw.simg bash "$1" "${@:2}"

