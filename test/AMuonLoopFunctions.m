(* ::Package:: *)

(*One loop functions and contributions for reproducing numerical calculations*)

(*One Loop Functions*)

OneLoopA[x_] := (2 - 9*x + 18*x^2 - 11*x^3 + 6*x^3*Log[x])/(1 - x)^4
OneLoopB[x_] := (2 (1 - 6 x + 3 x^2 + 2 x^3 - 6 x^2 Log[x]))/(1 - x)^4
OneLoopC[x_] := (3 (1 - x^2 + 2 x Log[x]))/(1 - x)^3
OneLoopD[x_] := (16 - 45 x + 36 x^2 - 7 x^3 + 
   6 (2 - 3 x)*Log[x])/(1 - x)^4
OneLoopE[x_] := (2 (2 + 3 x - 6 x^2 + x^3 + 6 x Log[x]))/(1 - x)^4
OneLoopF[x_] := (3 (-3 + 4 x - x^2 - 2 Log[x]))/(2 (1 - x)^3)
OneLoopH[x_] := -(((-1 + x) (41 - 166 x + 113 x^2) - 
    6 x^2 (-15 + 13 x) Log[x]) /(-1 + x)^4)
OneLoopI[x_] := (1 - x + x*Log[x])/(1 - x)^2
OneLoopJ[x_] := (7 - 33 x + 57 x^2 - 31 x^3 + 
   6 x^2 (3 x - 1) Log[x])/(1 - x)^4
OneLoopK[x_] := (1 - 4 x + 3 x^2 - 2 x^2 Log[x])/(1 - x)^3
OneLoopL[x_] := (2 + 27*x - 54*x^2 + 25*x^3 - 
   6*(2 - 9*x + 6*x^2)*Log[x])/(1 - x)^4
OneLoopM[x_] := (4 - 9 x + 5 x^3 + 6 (1 - 2 x) x Log[x])/(1 - x)^4
OneLoopN[x_, y_] := 
 y/(y - 1) (x (x - y)^2 Log[
      x] + (x - 1) ((x - y) (y - 1) - (x - 1) Log[x/y]))/((1 - 
      x)^2 (x - y)^2)
OneLoopNSame[x_] :=  (1 - 4 x + 3 x^2 - 2 x^2 Log[x])/(2 (1 - x)^3)

(* 1 Loop Contributions *)

(* SSF Contributions *)

SSFA1L[x_,SFinleft_,SFinright_,SFoutleft_,SFoutright_] := -(((2 - 9 x + 18 x^2 - 11 x^3 + 
    6 x^3 Log[x]) SFinleft[] SFoutright[])/(18 (-1 + x)^4))

SSFA1R[x_,SFinleft_,SFinright_,SFoutleft_,SFoutright_] := -(((2 - 9 x + 18 x^2 - 11 x^3 + 
    6 x^3 Log[x]) SFinright[] SFoutleft[])/(18 (-1 + x)^4))

SSFA2L[x_,SFinleft_,SFinright_,SFoutleft_,SFoutright_,mF_,mi_,mj_] := (mF (1 - x^2 + 2 x Log[x]) SFinleft[] SFoutleft[])/(
 mj (-1 + x)^3) + ((-1 + 6 x - 3 x^2 - 2 x^3 + 
    6 x^2 Log[x]) SFinright[] SFoutleft[])/(6 (-1 + x)^4) + (
 mi (-1 + 6 x - 3 x^2 - 2 x^3 + 
    6 x^2 Log[x]) SFinleft[] SFoutright[])/(6 mj (-1 + x)^4)

SSFA2R[x_,SFinleft_,SFinright_,SFoutleft_,SFoutright_,mF_,mi_,mj_] := (mi (-1 + 6 x - 3 x^2 - 2 x^3 + 
    6 x^2 Log[x]) SFinright[] SFoutleft[])/(
 6 mj (-1 + x)^4) + ((-1 + 6 x - 3 x^2 - 2 x^3 + 
    6 x^2 Log[x]) SFinleft[] SFoutright[])/(6 (-1 + x)^4) + (
 mF (1 - x^2 + 2 x Log[x]) SFinright[] SFoutright[])/(mj (-1 + x)^3)
 
FFSA1L[x_,SFinleft_,SFinright_,SFoutleft_,SFoutright_] := ((-16 + 45 x - 36 x^2 + 7 x^3 + 
   6 (-2 + 3 x) Log[x]) SFinleft[] SFoutright[])/(18 (-1 + x)^4)
   
FFSA1R[x_,SFinleft_,SFinright_,SFoutleft_,SFoutright_] := ((-16 + 45 x - 36 x^2 + 7 x^3 + 
   6 (-2 + 3 x) Log[x]) SFinright[] SFoutleft[])/(18 (-1 + x)^4)
   
FFSA2L[x_,SFinleft_,SFinright_,SFoutleft_,SFoutright_,mF_,mi_,mj_] := -((mF (3 - 4 x + x^2 + 2 Log[x]) SFinleft[] SFoutleft[])/(
  mj (-1 + x)^3)) - ((2 + 3 x - 6 x^2 + x^3 + 
    6 x Log[x]) SFinright[] SFoutleft[])/(6 (-1 + x)^4) - (
 mi (2 + 3 x - 6 x^2 + x^3 + 6 x Log[x]) SFinleft[] SFoutright[])/(
 6 mj (-1 + x)^4)
 
FFSA2R[x_,SFinleft_,SFinright_,SFoutleft_,SFoutright_,mF_,mi_,mj_] := -((mF (3 - 4 x + x^2 + 2 Log[x]) SFinleft[] SFoutleft[])/(
  mj (-1 + x)^3)) - ((2 + 3 x - 6 x^2 + x^3 + 
    6 x Log[x]) SFinright[] SFoutleft[])/(6 (-1 + x)^4) - (
 mi SFinleft (2 + 3 x - 6 x^2 + x^3 + 6 x Log[x]) SFoutright[])/(
 6 mj (-1 + x)^4)

VVFA1L[x_,VFinleft_,VFinright_,VFoutleft_,VFoutright_] := -((((-1 + x) (41 - 166 x + 113 x^2) - 
    6 x^2 (-15 + 13 x) Log[x]) VFinleft[] VFoutleft[])/(18 (-1 + x)^4))

VVFA1R[x_,VFinleft_,VFinright_,VFoutleft_,VFoutright_] := -((((-1 + x) (41 - 166 x + 113 x^2) + 
    6 (15 - 13 x) x^2 Log[x]) VFinright[] VFoutright[])/(
 18 (-1 + x)^4))

VVFA2L[x_,VFinleft_,VFinright_,VFoutleft_,VFoutright_,mF_,mi_,mj_] := (mi (7 - 33 x + 57 x^2 - 31 x^3 + 
    6 x^2 (-1 + 3 x) Log[x]) VFinleft[] VFoutleft[])/(
 6 mj (-1 + x)^4) + (
 3 mF (-1 + 4 x - 3 x^2 + 2 x^2 Log[x]) VFinleft[] VFoutright[])/(
 mj (-1 + x)^3) + ((7 - 33 x + 57 x^2 - 31 x^3 + 
    6 x^2 (-1 + 3 x) Log[x]) VFinright[] VFoutright[])/(6 (-1 + x)^4)

VVFA2R[x_,VFinleft_,VFinright_,VFoutleft_,VFoutright_,mF_,mi_,mj_] := ((7 - 33 x + 57 x^2 - 31 x^3 + 
    6 x^2 (-1 + 3 x) Log[x]) VFinleft[] VFoutleft[])/(
 6 (-1 + x)^4) + (
 3 mF (-1 + 4 x - 3 x^2 + 2 x^2 Log[x]) VFinright[] VFoutleft[])/(
 mj (-1 + x)^3) + (
 mi (7 - 33 x + 57 x^2 - 31 x^3 + 
    6 x^2 (-1 + 3 x) Log[x]) VFinright[] VFoutright[])/(
 6 mj (-1 + x)^4)

FFVA1L[x_,VFinleft_,VFinright_,VFoutleft_,VFoutright_] := ((2 + 27 x - 54 x^2 + 25 x^3 - 
   6 (2 - 9 x + 6 x^2) Log[x]) VFinleft[] VFoutleft[])/(9 (-1 + x)^4)

FFVA1R[x_,VFinleft_,VFinright_,VFoutleft_,VFoutright_] := ((2 + 27 x - 54 x^2 + 25 x^3 - 
   6 (2 - 9 x + 6 x^2) Log[x]) VFinright[] VFoutright[])/(9 (-1 + x)^4
 )
 
FFVA2L[x_,VFinleft_,VFinright_,VFoutleft_,VFoutright_,mF_,mi_,mj_] := (mi (4 - 9 x + 5 x^3 + 6 (1 - 2 x) x Log[x]) VFinleft[] VFoutleft[])/(
 3 mj (-1 + x)^4) - (
 4 mF (-1 + x^2 - 2 x Log[x]) VFinleft[] VFoutright[])/(
 mj (-1 + x)^3) + ((4 - 9 x + 5 x^3 + 
    6 (1 - 2 x) x Log[x]) VFinright[] VFoutright[])/(3 (-1 + x)^4)

FFVA2R[x_,VFinleft_,VFinright_,VFoutleft_,VFoutright_,mF_,mi_,mj_] := ((4 - 9 x + 5 x^3 + 6 (1 - 2 x) x Log[x]) VFinleft[] VFoutleft[])/(
 3 (-1 + x)^4) - (
 4 mF (-1 + x^2 - 2 x Log[x]) VFinright[] VFoutleft[])/(
 mj (-1 + x)^3) + (
 mi (4 - 9 x + 5 x^3 + 
    6 (1 - 2 x) x Log[x]) VFinright[] VFoutright[])/(3 mj (-1 + x)^4)

VSFA2L[x_,y_,VFinleft_,VFinright_,SFoutleft_,SFoutright_,mF_,mi_,mj_] := (y (x (x - y)^2 Log[
       x] + (-1 + 
        x) ((x - y) (-1 + y) - (-1 + x) x Log[x/
          y])) SFoutleft[] VFinleft[])/(mj (-1 + x)^2 (x - 
     y)^2 (-1 + y))

VSFA2R[x_,y_,VFinleft_,VFinright_,SFoutleft_,SFoutright_,mF_,mi_,mj_] := (y (x (x - y)^2 Log[
       x] + (-1 + 
        x) ((x - y) (-1 + y) - (-1 + x) x Log[x/
          y])) SFoutright[] VFinright[])/(mj (-1 + x)^2 (x - 
     y)^2 (-1 + y))
     
SVFA2L[x_,y_,SFinleft_,SFinright_,VFoutleft_,VFoutright_,mF_,mi_,mj_] := ((x (x - y)^2 y Log[
     x] + (-1 + 
      x) y ((x - y) (-1 + y) - (-1 + x) x Log[x/
        y])) SFinleft[] VFoutright[])/(mj (-1 + x)^2 (x - 
   y)^2 (-1 + y))

SVFA2R[x_,y_,SFinleft_,SFinright_,VFoutleft_,VFoutright_,mF_,mi_,mj_] := ((x (x - y)^2 y Log[
     x] + (-1 + 
      x) y ((x - y) (-1 + y) - (-1 + x) x Log[x/
        y])) SFinright[] VFoutleft[])/(mj (-1 + x)^2 (x - 
   y)^2 (-1 + y))

