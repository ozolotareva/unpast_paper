#!/bin/bash
docker cp /local/DESMOND2_data/v6/ desmond2_eval_container_bc:/local/DESMOND2_data/
docker exec -it desmond2_eval_container_bc python run_biclusters.py 8
#docker cp desmond2_eval_container:/tmp/desmond2_bicluster_eval_results/ ~/ENCORE_functions/biclustering/summaries/
