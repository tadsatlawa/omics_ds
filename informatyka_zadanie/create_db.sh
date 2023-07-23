#!/bin/bash

/apps/ncbi-blast/2.11.0/bin/makeblastdb -in $1 -dbtype prot -parse_seqids -out $2
