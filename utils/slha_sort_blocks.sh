#!/bin/sh

# This script sorts the blocks of the given SLHA input.
#
# Examples:
#
#   slha_sort_blocks.sh input.slha
#
#   cat input.slha | slha_sort_blocks.sh

[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

input=$(cat ${input})

blocks=$(echo "${input}" | awk '{ if (tolower($1) ~ "block") print $2 }' | sort)

for block in ${blocks} ; do
    echo "${input}" | awk -f $(dirname $0)/print_slha_block.awk -v block="${block}" -v omit_comments=1
done
