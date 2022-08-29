#!/bin/bash

tool_name=$1
case_prefix=$2
result_file=$5

echo "Starting biclustering"
eval $3
echo "Evaluating"
eval