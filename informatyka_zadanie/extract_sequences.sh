#!/bin/bash

# Function to split the file on the n-th occurrence of a line beginning
split_file_on_nth_occurrence() {
  local input_file="$1"
  local output_prefix="$2"
  local line_beginning="$3"
  local n="$4"

  local current_occurrence=0
  local output_file="${output_prefix}_part${current_occurrence}.txt"

  # Create an initial output file
  > "$output_file"

  while IFS= read -r line; do
    if [[ "$line" == "$line_beginning"* ]]; then
      ((current_occurrence++))
      if ((current_occurrence > n)); then
        # Start writing to a new output file when the n-th occurrence is reached
	break
        current_occurrence=1
        output_file="${output_prefix}_part${current_occurrence}.txt"
        > "$output_file"
      fi
    fi
    echo "$line" >> "$output_file"
  done < "$input_file"
}

# Check if the required number of arguments is provided
if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <input_file> <output_prefix> <line_beginning> <n>"
  exit 1
fi

input_file="$1"
output_prefix="$2"
line_beginning="$3"
n="$4"

split_file_on_nth_occurrence "$input_file" "$output_prefix" "$line_beginning" "$n"

