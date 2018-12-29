#!/bin/bash
set -o nounset
set -o errexit

function print_the_help {
  echo "usage:    run.sh -o <OUTPUT_DIR> "
  echo "options: "
  echo "        -o,--output <OUTPUT_DIR>  output directory"
  echo "        -c,--create               create the output directory if it does not exist"
  echo "        -h,--help                 print this help"
  echo " " 
  exit 
}

function yes_or_no {
  while true; do
    read -p "$* [y/n]: " yn
    case $yn in
      [Yy]*) return 0 ;;
      [Nn]*) echo "No entered" ; return 1 ;;
    esac
  done
}

OUTPUT_DIR=
CREATE_DIR=

if [[ $# -eq 0 ]] ; then
  print_the_help
  exit 
fi


POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    -h|--help)
      shift # past argument
      print_the_help
      ;;
    #-r|--run)
    #  RUNNUM="$2"
    #  shift # past argument
    #  shift # past value
    #  ;;
    -c|--create)
      CREATE_DIR=1
      shift # past argument
      #shift # past value
      ;;
    -o|--output)
      OUTPUT_DIR="$2"
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $1"
      print_the_help
      shift # past argument
      ;;
  esac
done
#set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ ! -d "${OUTPUT_DIR}" ]] ; then
  if [[ -z "${CREATE_DIR}" ]] ; then
    echo "error: output directory, ${OUTPUT_DIR}, not found." 
    exit 73 
  else 
    echo "Creating director ${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}"
  fi 
fi

#yes_or_no "Upload these plots ? " && some_command
# args

tempfile=$(mktemp)
./config_gen.py > $tempfile
pcsim -c $tempfile -o ${OUTPUT_DIR} -r1
rm "$tempfile"


