#!/bin/bash

# Check for correct number of arguments
if [ $# -ne 2 ]; then
  echo "Usage: $0 /path/to/files max_number"
  exit 1
fi

dir="$1"
max_number="$2"

# Extract numbers after the last underscore before .Rdata
found_numbers=$(find "$dir" -type f -name "*.Rdata" | grep -oP '_\K[0-9]+(?=\.Rdata)' | sort -n | uniq)

expected_numbers=$(seq 1 "$max_number")
# Compare and find missing numbers
missing_numbers=$(comm -23 <(echo "$expected_numbers") <(echo "$found_numbers"))

echo "Missing numbers:"
echo "$missing_numbers"

