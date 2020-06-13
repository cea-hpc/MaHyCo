#!/usr/bin/bash
#
# This script returns 0 if the difference with REF_BRANCH has formatted according to STYLE
#
set -euo pipefail 

function usage
{
    echo "Usage: $0 [-h|--help]"
    echo ""
    echo "Checks the formatting of the current branch by comparing it to the reference one"
    echo "Options:"
    echo "-h|--help : print this help and exits"
}

function check_formatting
{
    local readonly REF_BRANCH=$1
    local readonly STYLE=$2
    shift 2
    local readonly FILE_TYPES="$@"

    local OUTPUT=$(git diff -U0 --no-color ${REF_BRANCH} -- ${FILE_TYPES} | python2 $(which clang-format-diff) -p1 -style=${STYLE})
    
    if [[ -z "${OUTPUT}" ]]; then
        echo "Everything correctly formatted!"
        exit 0
    else
        echo "Some files are not correclty formatted!"
        echo "${OUTPUT}"
        exit 1
    fi
}

while getopts ":h:-help" option; do
    case "${option}" in
        h)
            usage
            exit 0
            ;;
        -help)
            usage
            exit 0
            ;;
        *)
            usage
            exit 255
            ;;
    esac
done

if [[ $# -ne 0 ]]; then
    usage
    exit 255
fi

readonly REFERENCE_BRANCH="master"
readonly FILES_TYPES="main.cc eucclhyd_remap/*.cc eucclhyd_remap*.h"
readonly STYLE_EXPECTED="Google"

check_formatting $REFERENCE_BRANCH $STYLE_EXPECTED $FILES_TYPES 