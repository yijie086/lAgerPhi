#!/bin/bash

REG=${CI_PROJECT_NAMESPACE}/${CI_PROJECT_NAME}
echo $REG | tr '[:upper:]' '[:lower:]'
