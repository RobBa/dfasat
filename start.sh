#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 ini-file input-file"
    exit 1
fi

#echo $param

./cmake-build-debug/dfasat "$1" "$2"
