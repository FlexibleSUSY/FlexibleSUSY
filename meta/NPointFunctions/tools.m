(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage@"tools`";

unzipRule::usage = "
@brief Expands a set of compact rules into the full one. Example:
       In[1]:= rules = {
                  a -> b,
                  {c, d} -> e
               };
               unzipRule[rules]
       Out[1]= {a -> b, c -> e, d -> e}";

Begin@"`Private`";

unzipRule[rules:{Rule[_|{__}, _]..}] :=
   rules /. Rule[lhs:{__}, rhs_] :> Sequence@@ (Rule[#, rhs]&/@ lhs);
unzipRule // secure;

End[];
Block[{$ContextPath}, EndPackage[]];
