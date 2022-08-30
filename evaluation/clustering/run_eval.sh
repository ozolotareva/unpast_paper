#!/bin/bash
# docker-compose build; docker-compose up -d; ./run_eval.sh
docker cp /local/DESMOND2_data_simulated/simulated desmond2_eval_container_clustering:/local/DESMOND2_data_simulated/
screen docker exec -it desmond2_eval_container_clustering python run_clusters.py 1
docker cp desmond2_eval_container_clustering:/tmp/desmond2_cluster_eval_results/ ~/ENCORE_functions/clustering/summaries/
