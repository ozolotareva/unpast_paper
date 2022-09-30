#!/bin/bash
docker-compose build
docker-compose up -d
#./run_eval.sh
screen -L -Logfile bicluster_BC_eval.log -S bicluster_BC_eval ./run_eval.sh
