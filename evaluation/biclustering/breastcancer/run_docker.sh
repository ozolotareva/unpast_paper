#!/bin/bash
rm bicluster_BC_eval.log
rsync -aHh --info=progress2 ../jbiclustge/ .
docker-compose build
docker-compose up -d
screen -L -Logfile bicluster_BC_eval.log -S bicluster_BC_eval ./run_eval.sh
