#!/bin/bash

# Function to grep multiple files
grep ">" "$@" | wc -l

