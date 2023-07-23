#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH --time=02:10:00
#SBATCH -A g92-1513   
time /apps/ncbi-blast/2.11.0/bin/blastp -query first100_part0.txt -db zebrafish_db -num_threads 24 > out
time /apps/ncbi-blast/2.11.0/bin/blastp -query first100_part0.txt -db zebrafish_db -num_threads 4 > out
time /apps/ncbi-blast/2.11.0/bin/blastp -query first100_part0.txt -db zebrafish_db -num_threads 1 > out
