#!/bin/bash -l

# Function to calculate execution time
calculate_execution_time() {
  start_time=$1
  end_time=$2
  echo "Execution Time: $((end_time - start_time)) seconds"
}

# Function to execute command, save output, and calculate execution time
execute_command_and_save_output() {
  local command="$1"
  local output_file="$2"
  echo $command
  start_time=$(date +%s)
  output=$($command)
  end_time=$(date +%s)

  echo "$output" > "$output_file"
  calculate_execution_time "$start_time" "$end_time"
}

execute_command_and_save_output $1 $2 
