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

Needs["TestSuite`", "TestSuite.m"];
Needs["ThresholdCorrections`", "ThresholdCorrections.m"];

Print["testing InvertRelation[] ..."];

TestEquality[ThresholdCorrections`Private`InvertRelation[A,B,A],
             {A,B}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[A,B + C,A],
             {A,B + C}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[Transpose[A],B,A],
             {A,SARAH`Tp[B]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[SARAH`Tp[A],B,A],
             {A,SARAH`Tp[B]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[ConjugateTranspose[A],B,A],
             {A,SARAH`Adj[B]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[SARAH`Adj[A],B,A],
             {A,SARAH`Adj[B]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[FlexibleSUSY`Diag[A],B,A],
             {A,FlexibleSUSY`Diag[B]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[SARAH`MatMul[A,V],C,A],
             {A,SARAH`MatMul[C,SARAH`Adj[V]]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[SARAH`MatMul[V,A],C,A],
             {A,SARAH`MatMul[SARAH`Adj[V],C]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[SARAH`MatMul[U,A,V],C,A],
             {A,SARAH`MatMul[SARAH`Adj[U],C,SARAH`Adj[V]]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[SARAH`MatMul[SARAH`Adj[U],A,V],C,A],
             {A,SARAH`MatMul[U,C,SARAH`Adj[V]]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[SARAH`MatMul[U,A,SARAH`Adj[V]],C,A],
             {A,SARAH`MatMul[SARAH`Adj[U],C,V]}
            ];

TestEquality[ThresholdCorrections`Private`InvertRelation[SARAH`MatMul[SARAH`Adj[U],A,SARAH`Adj[V]],C,A],
             {A,SARAH`MatMul[U,C,V]}
            ];


Print["testing ExtractSymbols[] ..."];

TestEquality[ThresholdCorrections`Private`ExtractSymbols[a], {a}];
TestEquality[ThresholdCorrections`Private`ExtractSymbols[a[i,j]], {a}];
TestEquality[ThresholdCorrections`Private`ExtractSymbols[a[i,j] + b], {b,a}];
TestEquality[ThresholdCorrections`Private`ExtractSymbols[a[i,j] b], {b,a}];
TestEquality[ThresholdCorrections`Private`ExtractSymbols[SARAH`sum[i,1,3,Y]], {Y}];
TestEquality[ThresholdCorrections`Private`ExtractSymbols[SARAH`sum[i,1,3,Y[i,j]]], {Y}];
TestEquality[ThresholdCorrections`Private`ExtractSymbols[SARAH`sum[i,1,3,Y[i,j] U[i,k]]], {U,Y}];


Print["testing ToMatrixExpression[] ..."];

TestEquality[ThresholdCorrections`Private`ToMatrixExpression[{}],
             Null
            ];

TestEquality[ThresholdCorrections`Private`ToMatrixExpression[{{Y[1,1]}}],
             Y
            ];

TestEquality[ThresholdCorrections`Private`ToMatrixExpression[{{Y[1,1],Y[1,2]},{Y[2,1],Y[2,2]}}],
             Y
            ];

TestEquality[ThresholdCorrections`Private`ToMatrixExpression[{{Y[1,1],Y[2,1]},{Y[1,2],Y[2,2]}}],
             SARAH`Tp[Y]
            ];

TestEquality[ThresholdCorrections`Private`ToMatrixExpression[{{Y[1,1],0},{0,Y[2,2]}}],
             FlexibleSUSY`Diag[Y]
            ];

TestEquality[ThresholdCorrections`Private`ToMatrixExpression[{{SARAH`sum[i,1,2,Y[1,i] U[i,1]],
                                          SARAH`sum[i,1,2,Y[1,i] U[i,2]]},
                                         {SARAH`sum[i,1,2,Y[2,i] U[i,1]],
                                          SARAH`sum[i,1,2,Y[2,i] U[i,2]]}}],
             SARAH`MatMul[Y,U]
            ];

TestEquality[ThresholdCorrections`Private`ToMatrixExpression[{{SARAH`sum[i,1,2,SARAH`sum[k,1,2,V[1,i] Y[i,k] U[k,1]]],
                                          SARAH`sum[i,1,2,SARAH`sum[k,1,2,V[1,i] Y[i,k] U[k,2]]]},
                                         {SARAH`sum[i,1,2,SARAH`sum[k,1,2,V[2,i] Y[i,k] U[k,1]]],
                                          SARAH`sum[i,1,2,SARAH`sum[k,1,2,V[2,i] Y[i,k] U[k,2]]]}}],
             SARAH`MatMul[V,Y,U]
            ];


PrintTestSummary[];
