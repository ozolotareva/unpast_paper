#!/bin/bash
# ssh bbb1417@llaima.zbh.uni-hamburg.de
# docker-compose build; docker-compose up -d; ./run_eval.sh
docker cp /local/DESMOND2_data_simulated/simulated desmond2_eval_container_clustering:/local/DESMOND2_data_simulated/
screen docker exec -it desmond2_eval_container_clustering python run_clusters.py 1

docker cp desmond2_eval_container_clustering:/tmp/ ~/ENCORE_functions/clustering/
