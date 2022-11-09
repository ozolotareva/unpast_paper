#!/bin/bash
rm bicluster_eval.log
rsync -aHh --info=progress2 ../jbiclustge/ .
docker-compose build
docker-compose up -d
screen -L -Logfile bicluster_eval.log -S bicluster_eval ./run_eval.sh
