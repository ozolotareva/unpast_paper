#!/bin/bash
docker cp /local/DESMOND2_data_simulated/simulated desmond2_eval_container:/local/DESMOND2_data_simulated/
screen docker exec -it desmond2_eval_container python run_biclusters.py 1
docker cp desmond2_eval_container:/tmp/desmond2_bicluster_eval_results/ ~/ENCORE_functions/biclustering/summaries/
