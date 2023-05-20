(****f* Mathematica/TwoLoopBubbles
* DESCRIPTION
*	* This file defines functions for numerical  calculation of two-loop
* 	vacuum diagrams (bubbles). 
*	* Two-loop self-energy diagrams with different 
*	  masses and the momentum expansion, Nuclear Physics B397 (1993) 123-142  
*
****)

(* Calculation of Two Loop Bubble *)

(* poles are excluded but only poles!! *)

(* Bubble with two equal masses *)

(* dimension -> mass *)

$CalcPrecition = 15;

(* <<MathWorld`SpecialFunctions`; *)

(****u* TwoLoopBubbles/Clausen
* DESCRIPTION
*	* Clausen Function
*	|latex \begin{eqnarray}
*	|latex Cl_n(x) & = & \sum\limits_{k=1}^{\infty}\frac{\sin(k x)}{k^n} = 
*	|latex \frac{i}{2} \left(Li_n(e^{-i x}) - Li_n(e^{i x}) \right), \quad \mathrm{ n~even} \\
*	|latex Cl_n(x) & = & \sum\limits_{k=1}^{\infty}\frac{\cos(k x)}{k^n} =  
*	|latex \frac{1}{2} \left(Li_n(e^{-i x}) - Li_n(e^{i x}) \right), \quad \mathrm{ n~odd} 
*	|latex \end{eqnarray}
*	
****)

ClausenCl[n_Integer?EvenQ,x_] := (PolyLog[n, E^(-I x)] - PolyLog[n, E^(I x)])I/2
ClausenCl[n_Integer?OddQ,x_] :=  (PolyLog[n, E^(-I x)] + PolyLog[n, E^(I x)])/2

(****u* TwoLoopBubbles/Phi
* DESCRIPTION
*	* Phi(z) Aux function
*	|latex \begin{eqnarray}
*	|latex \Phi(z) & = & 4 \ln 4, \quad z = 1\\
*	|latex \Phi(z) & = & 4 \sqrt{\frac{z}{1-z}} Cl_2(2 \arcsin\sqrt{z}), \quad  0\leq z < 1\\
*	|latex \Phi(z) & = & \frac{1}{\lambda}
*	|latex \left( 2 \ln \left(\frac{1-\lambda}{2}\right)^2  
*	|latex 	- 4 Li_2 \left(\frac{1-\lambda}{2}\right) - \ln (4 z )^2  + \frac{\pi^2}{3} \right), 
*	|latex 	\quad\lambda = \sqrt{1 - 1/z}, z > 1
*	|latex \end{eqnarray}
*	
****)

Phi[z_?NumberQ] := N[4 * Sqrt[ z/(1-z) ] * ClausenCl[2, 2 * ArcSin[Sqrt[z]] ],$CalcPrecition] /; (z<1 && z>=0);
Phi[z_?NumberQ] := N[Sqrt[z /(z-1) ] * ( - 4 PolyLog[2, 1/2 * ( 1 - Sqrt[1 - 1/z] ) ]
				       + 2 Log[ 1/2 * (1 - Sqrt[1 - 1/z] ) ]^2  
				       -   Log[ 4 z ]^2 
				       + Pi^2/3 ), $CalcPrecition ] /; z>1;
Phi[1] := N[4 Log [4],$CalcPrecition];

(* mass squared! *)

(* mmu - scale *)
(****u* TwoLoopBubbles/Fin21
* DESCRIPTION
*	* Phi21[mm,mm,MM,mmu] - finite part of two-loop bubble with masses 
*	 two different mm,MM at scale mmu.
*	|latex \begin{eqnarray}
*	|latex \Phi_{21}(m^2,m^2,M^2,\mu^2) & = & 
*	|latex \frac{1}{2} \left((M^2 - 2 m^2) \left( 2 \ln^2 (m^2/\mu^2) - 6 \ln (m^2/\mu^2) + \zeta(2) + 7\right) \\
*	|latex & + & M^2  \ln (M^2/m^2) \left( 6 - 4 \ln (m^2/\mu^2) - \ln (M^2/m^2) \right) \\
*	|latex & + & (4 m^2 - M^2) \Phi(z)\right), \quad z = \frac{M^2}{4 m^2}
*	|latex \end{eqnarray}
*	
****)

Fin21[mm_?NumberQ, mm_?NumberQ, MM_?NumberQ, mmu_?NumberQ] :=
	Module[{z,tmp}, z = MM/(4*mm) ; (* Print[z]; *)
		    tmp = 0.5 * N[-(2 mm + MM ) ( 2 Log[mm/mmu]^2 - 6 Log[mm/mmu]+ Zeta[2] + 7 )
		    	      - Log[MM/mm]^2 MM + 2 Log[MM/mm] * ( 3 - 2 Log[mm/mmu]) * MM 
			      + (4 mm - MM) * Phi[z],$CalcPrecition ];
	  	    tmp
	      ];
(* three different masses *)

LambdaSquared[x_?NumberQ, y_?NumberQ] := (1-x-y)^2-4*x*y;

(* x < 1 and y < 1 *)
Phi[x_?NumberQ, y_?NumberQ] := 
        Module[{lambda,tmp}, lambda = Sqrt[LambdaSquared[x,y]]; 
		tmp := N[ 1/lambda *  
		( 2 Log[ ( 1 + x - y - lambda )/2 ] * Log[ (1 -x +y - lambda )/2 ] - Log[x]*Log[y]
		- 2 PolyLog[2, ( 1 + x - y - lambda )/2 ] - 2 PolyLog[2, (1 - x + y - lambda)/2] + Pi^2/3),$CalcPrecition ];
		tmp
	      ] /; (LambdaSquared[x,y] > 0 );
	      
Phi[x_?NumberQ, y_?NumberQ] := 
	Module[{lambda, tmp}, lambda = Sqrt[-LambdaSquared[x,y]];
		tmp := N[ 2/lambda * ( ClausenCl[2, 2 ArcCos[ (-1 + x + y )/(2 Sqrt[x y] )] ]
		                      +ClausenCl[2, 2 ArcCos[ (1 + x -y )/(2 Sqrt[x] ) ] ]
				      +ClausenCl[2, 2 ArcCos[ (1 - x + y )/(2 Sqrt[y] ) ] ]
				     ), $CalcPrecition];
		tmp
	      ] /; (LambdaSquared[x,y] < 0 );

Fin3[mm1_?NumberQ, mm2_?NumberQ, mm3_?NumberQ, mmu_?NumberQ] :=
	Module[{x,y, tmplist, tmp, mm}, 
		tmplist = Sort[{mm1,mm2,mm3}];
		mm = tmplist[[3]];
		x = tmplist[[1]]/mm;
		y = tmplist[[2]]/mm;
		(* Print[x]; Print[y];Print[ LambdaSquared[x,y]]; *)
		If[LambdaSquared[x,y] == 0,
		     tmp = N[-1/2 * (x * Log[x]^2 
			     + ( (x + y - 1) * Log[y]  + 2 x * ( 2 Log[mm/mmu] - 3) ) * Log[x] 
			     +  y Log[y]^2 
			     + 2 y Log[y] ( 2 Log[mm/mmu] - 3) 
			     + (1 + x + y ) * (2 Log[mm/mmu]^2 - 6 Log[mm/mmu] + Zeta[2] + 7)
			     ) * mm, $CalcPrecition],
		    tmp = N[mm * (
		    	1/2 * ( -2 Log[mm/mmu]^2 + 6 Log[mm/mmu] - LambdaSquared[x,y]*Phi[x,y] 
				- Zeta[2] + Log[x] Log[y] - 7 ) 
		       -1/2 * x * ( Log[x]^2 + ( Log[y] + 4 Log[mm/mmu] - 6 ) Log[x] + 2 Log[mm/mmu]^2 
		       		+ Zeta[2] - 6 Log[mm/mmu] + 7 )
		       -1/2 * y * ( Log[y]^2 + ( Log[x] + 4 Log[mm/mmu] - 6 ) Log[y] + 2 Log[mm/mmu]^2 
		       		+ Zeta[2] - 6 Log[mm/mmu] + 7 )
				), $CalcPrecition]
		  ];
		  tmp
		];
			
(* additional part *)
(* Bubble with 1 mass equal to zero *)

Hmine[mm1_?NumberQ, mm2_?NumberQ] := 
    2 * PolyLog[2, 1-mm1/mm2] + 1/2 * Log[mm1/mm2]^2;

Fin20[mm1_?NumberQ, mm2_?NumberQ, mmu_?NumberQ] :=
	Module[{tmp}, tmp = N[1/2 * ( - (mm1 + mm2) * ( 7 + Zeta[2] )  
		            + 6 * (mm1 * Log[mm1/mmu] + mm2 * Log[mm2/mmu])
			    - 2 * (mm1 * Log[mm1/mmu]^2 + mm2 * Log[mm2/mmu]^2 )
			   +1/2 * (mm1 + mm2) * Log[mm1/mm2]^2 + (mm1-mm2)*Hmine[mm1,mm2] ),$CalcPrecition]; 
		      tmp
	];
