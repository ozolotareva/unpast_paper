#!/bin/bash
docker-compose build
docker-compose up -d
screen ./run_eval.sh
