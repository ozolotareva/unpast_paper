#!/bin/bash
# ssh bbb1417@llaima.zbh.uni-hamburg.de
# docker-compose build; docker-compose up -d; screen ./run_eval.sh
screen docker exec -it desmond2_eval_container_clustering2 python run_clusters.py 1
