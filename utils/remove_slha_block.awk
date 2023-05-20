# This script removes a block or an entry from a SLHA file.
# Capitalization is ignored.
#
# Examples:
#
#   awk -f config/remove_slha_block.awk -v block=MINPAR input.slha
#
#   awk -f config/remove_slha_block.awk -v block=MINPAR -v entry=1 input.slha
#
#   awk -f config/remove_slha_block.awk -v block=MINPAR -v entry=1:1 input.slha

BEGIN {
   if (block == "") {
      print "Error: block name not defined"
      print "   Please define the block name with -v block=<block-name>"
      exit 1
   }

   is_block = 0
   len = split(entry,k,":");
}
{
   # SLHA defines that the first character of a block statement is
   # 'B'.  However, for compatibility we allow for 'b' as well.


   pattern     = "^block[ \t\n\r\f]*" tolower(block) "([^a-zA-Z0-9_].*)?$";
   not_pattern = "^block[ \t\n\r\f]*.*$";

   if (tolower($0) ~ pattern) {
      is_block = 1
   } else if (tolower($0) ~ not_pattern) {
      is_block = 0
   }

   if (!is_block)
      print $0
   else {
      # check if given indices match the current line
      matches = 1;

      for (i in k) {
         if ($(i) != k[i])
            matches = 0
      }

      if (!matches)
         print $0
   }
}
