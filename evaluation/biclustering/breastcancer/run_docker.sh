#!/bin/bash
docker-compose build
docker-compose up -d
./run_eval.sh
#screen ./run_eval.sh
