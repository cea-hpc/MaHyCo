#!/usr/bin/bash
#
# This script returns 0 if the difference with REF_BRANCH has formatted according to STYLE
#

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
    local readonly FILE_TYPES=$2
    local readonly STYLE=$3

    local OUTPUT=$(git diff -U0 --no-color ${REF_BRANCH} -- ${FILE_TYPES} | clang-format-diff-8 -p1 -style=${STYLE})
    
    if [[ -z "${OUTPUT}" ]]; then
        echo "Everything correctly formatted!"
        exit 0
    else
        echo "Some files are not correclty formatted!"
        echo "${OUTPUT}"
        exit 1
    fi
}

if [[ $1 = '-h'||$1 = '--help' ]]; then
    usage
    exit 0
elif [[ $# != 0 ]]; then
    usage
    exit 255
fi

readonly REFERENCE_BRANCH="master"
readonly FILES_TYPES="*.cc *.h"
readonly STYLE_EXPECTED="Google"

check_formatting $REFERENCE_BRANCH $FILE_TYPES $STYLE_EXPECTED