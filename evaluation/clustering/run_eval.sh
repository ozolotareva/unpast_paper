#!/bin/bash
# ssh bbb1417@llaima.zbh.uni-hamburg.de
# docker-compose build; docker-compose up -d; screen ./run_eval.sh
#docker cp /local/DESMOND2_data_simulated/simulated desmond2_eval_container_clustering2:/root/projects/data/simulated
screen docker exec -it desmond2_eval_container_clustering2 python run_clusters.py 1
# screen docker exec -it desmond2_eval_container_clustering python run_clusters_realdata.py 1
