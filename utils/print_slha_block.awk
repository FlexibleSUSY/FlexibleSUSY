# This script prints a block from a SLHA file.
# Capitalization is ignored.
#
# Parameters:
#
#   block : name of block to print
#   omit_comments : print comment lines (0) (= default) or omit them (1)
#
# Example:
#
#   awk -f print_slha_block.awk -v block=MINPAR input.slha -v omit_comments=0

BEGIN {
   is_block = 0
   if (block == "") {
      print "Error: block name not defined"
      print "   Please define the block name with -v block=<block-name>"
      exit 1
   }
}
{
   # SLHA defines that the first character of a block statement is
   # 'B'.  However, for compatibility we allow for 'b' as well.

   pattern     = "^block[[:blank:]]*" tolower(block) "([^[:graph:]].*)?$"
   not_pattern = "^block[[:blank:]]*.*$"

   # line is a comment
   if (omit_comments == 1 && $0 ~ "^[[:blank:]]*#")
       next

   if (tolower($0) ~ pattern) {
      is_block = 1
   } else if (tolower($0) ~ not_pattern) {
      is_block = 0
   }

   if (is_block)
      print $0
}
