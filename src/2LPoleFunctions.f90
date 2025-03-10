! This file is part of SARAH 4.9.3
! Author: Florian Staub, Mark Goodsell

Module Pole2LFunctions

Integer, Parameter :: dp = Selected_real_kind(15,300)
Real(dp), parameter :: epsD = 1E-8_dp, epsM = 1E-10_dp
Real(dp), Parameter :: Pi = 3.1415926535897932384626433832795029_dp           &
   & , Pi2 = Pi**2, Pi3 = Pi2*Pi                                              &
   & , pi2o3 = Pi2 / 3._dp, pi2o6 =0.5_dp * Pi2o3, pi2o12 = 0.5_dp * Pi2o6
Complex(dp), Parameter :: ZeroC=(0._dp,0._dp)

 Interface Smalldp
  Module Procedure SmalldpR, SmalldpC
 End Interface

Contains

!!!!!!!!!!! TADPOLE FUNCTIONS (subroutine wrappers to interface with C/C++)

Subroutine CTfSS(out,x,y,z,Q2)
Implicit None
real(dp), intent(in):: x,y,z,q2
real(dp) :: out
out = TfSS(x,y,z,Q2)
End Subroutine CTfSS

Subroutine CTfSSS(out,x,y,z,Q2)
Implicit None
real(dp), intent(in):: x,y,z,q2
real(dp) :: out
out = TfSSS(x,y,z,Q2)
End Subroutine CTfSSS

Subroutine CTfSSSS(out,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: out
out = TfSSSS(x,y,z,u,Q2)
End Subroutine CTfSSSS

Subroutine CTfSSFF(out,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: out
out = TfSSFF(x,y,z,u,Q2)
End Subroutine CTfSSFF

Subroutine CTfSSFbFb(out,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: out
out = TfSSFbFb(x,y,z,u,Q2)
End Subroutine CTfSSFbFb

Subroutine CTfFFFbS(out,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: out
out = TfFFFbS(x,y,z,u,Q2)
End Subroutine CTfFFFbS

Subroutine CTfFFbFS(out,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: out
out = TfFFbFS(x,y,z,u,Q2)
End Subroutine CTfFFbFS

Subroutine CTfFbFbFbS(out,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: out
out = TfFbFbFbS(x,y,z,u,Q2)
End Subroutine CTfFbFbFbS

Subroutine CTfSV(out,x,Q2)
Implicit None
real(dp), intent(in):: x,q2
real(dp) :: out
out = TfSV(x,Q2)
End Subroutine CTfSV

Subroutine CTfFV(out,x,Q2)
Implicit None
real(dp), intent(in):: x,q2
real(dp) :: out
out = TfFV(x,Q2)
End Subroutine CTfFV

!!!!!!!!!!! TADPOLE FUNCTIONS


real(dp) Function TfSS(x,y,z,Q2)
Implicit None
real(dp), intent(in):: x,y,z,q2
real(dp) :: res

if (smallc(z/q2)) then
   res = 0._dp
else
  res = real(-J0(z,Q2)*BB(x,y,q2),dp)
end if

TfSS=res
end function tfss


real(dp) Function TfSSS(x,y,z,Q2)
Implicit None
real(dp), intent(in):: x,y,z,q2
real(dp) :: res

res=-real(II(x,y,z,q2),dp)

TfSSS=res
end function tfsss

real(dp) Function TfSSSS(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

If(smalldp(x,y)) then
      res=real(-1._dp*IIp(x,z,u,q2),dp)
else
   res=Real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
end if

TfSSSS=res
end function tfssss



real(dp) Function TfSSFF(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

res=real(-(z+u-y)*UU(x,y,z,u,q2)-II(x,z,u,q2)+BB(x,y,q2)*(J0(z,q2)+J0(u,q2)),dp)

TfSSFF=res
end function tfssff

real(dp) Function TfSSFbFb(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

RES=real(UU(x,y,z,u,q2),dp)

TfSSFbFb=res
end function tfssfbfb

real(dp) Function TfFFFbS(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

RES=real(II(x,z,u,q2)-(x+y)/2._dp*UU(x,y,z,u,q2),dp)

TfFFFbS=res
end function TfFFFbS

real(dp) Function TfFFbFS(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

RES=Fffs100(x,y,z,u,q2)

TfFFbFS=res
end function TfFFbFS

real(dp) Function TfFbFbFbS(x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: x,y,z,u,q2
real(dp) :: res

RES=real(UU(x,y,z,u,q2),dp)

TfFbFbFbS=res
end function tffbfbfbs

real(dp) Function TfSV(x,Q2)
Implicit None
real(dp), intent(in):: x,q2
real(dp) :: res

if(smallc(x/q2)) then
   res = 0._dp
else
   res=real(3._dp*x*(-5._dp+4._dp*log(x/q2)-(log(x/q2))**2)-J0(x,q2)+2._dp*x,dp)
end if

TfSV=res
end function tfsv

real(dp) Function TfFV(x,Q2)
Implicit None
real(dp), intent(in):: x,q2
real(dp) :: res


if(smallc(x/q2)) then
res=0._dp
else
   res=real(4._dp*x*(6._dp-7._dp*log(x/q2)+3._dp*(log(x/q2))**2),dp)
end if

TfFV=res
end function tffv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! MASS FUNCTIONS (subroutine wrappers to interface with C/C++)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine CWfSSSS (out,p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: out
out = WfSSSS(p2,x,y,z,u,Q2)
end Subroutine CWfSSSS

Subroutine CXfSSS(out,p2,x,y,z,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,q2
real(dp) :: out
out = XfSSS(p2,x,y,z,Q2)
end Subroutine CXfSSS

Subroutine CYfSSSS(out,p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: out
out = YfSSSS(p2,x,y,z,u,Q2)
end Subroutine CYfSSSS

Subroutine CZfSSSS(out,p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: out
out = ZfSSSS(p2,x,y,z,u,Q2)
end Subroutine CZfSSSS

Subroutine CSfSSS(out,p2,x,y,z,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,q2
real(dp) :: out
out = SfSSS(p2,x,y,z,Q2)
end Subroutine CSfSSS

Subroutine CUfSSSS(out,p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: out
out = UfSSSS(p2,x,y,z,u,Q2)
end Subroutine CUfSSSS

Subroutine CVfSSSSS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfSSSSS(p2,x,y,z,u,v,Q2)
end Subroutine CVfSSSSS

Subroutine CMfSSSSS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfSSSSS(p2,x,y,z,u,v,Q2)
end Subroutine CMfSSSSS

Subroutine CWfSSSV(out,p2,x,Q2)
Implicit None
real(dp), intent(in):: p2,x,q2
real(dp) :: out
out = WfSSSV(p2,x,Q2)
end Subroutine CWfSSSV

Subroutine CMfSSSSV(out,p2,x,y,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,q2
real(dp) :: out
out = MfSSSSV(p2,x,y,Q2)
end Subroutine CMfSSSSV

Subroutine CWfSSFF (out,p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: out
out = WfSSFF(p2,x,y,z,u,Q2)
end Subroutine CWfSSFF

Subroutine CWfSSFbFb (out,p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: out
out = WfSSFbFb(p2,x,y,z,u,Q2)
end Subroutine CWfSSFbFb

Subroutine CMfFbFbFbFbS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfFbFbFbFbS(p2,x,y,z,u,v,Q2)
end Subroutine CMfFbFbFbFbS

Subroutine CMfFFbFbFS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfFFbFbFS(p2,x,y,z,u,v,Q2)
end Subroutine CMfFFbFbFS

Subroutine CMfFFbFFbS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfFFbFFbS(p2,x,y,z,u,v,Q2)
end Subroutine CMfFFbFFbS

Subroutine CMfFFFbFbS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfFFFbFbS(p2,x,y,z,u,v,Q2)
end Subroutine CMfFFFbFbS

Subroutine CMfFFFFS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfFFFFS(p2,x,y,z,u,v,Q2)
end Subroutine CMfFFFFS

Subroutine CMfSFbSFbFb(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfSFbSFbFb(p2,x,y,z,u,v,Q2)
end Subroutine CMfSFbSFbFb

Subroutine CMfSFSFbF(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfSFSFbF(p2,x,y,z,u,v,Q2)
end Subroutine CMfSFSFbF

Subroutine CMfSFSFFb(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = MfSFSFFb(p2,x,y,z,u,v,Q2)
end Subroutine CMfSFSFFb

Subroutine CVfSSSFbFb(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfSSSFbFb(p2,x,y,z,u,v,Q2)
end Subroutine CVfSSSFbFb

Subroutine CVfSSSFF(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfSSSFF(p2,x,y,z,u,v,Q2)
end Subroutine CVfSSSFF

Subroutine CVfFbFbFbFbS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfFbFbFbFbS(p2,x,y,z,u,v,Q2)
end Subroutine CVfFbFbFbFbS

Subroutine CVfFbFFbFS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfFbFFbFS(p2,x,y,z,u,v,Q2)
end Subroutine CVfFbFFbFS

Subroutine CVfFbFFFbS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfFbFFFbS(p2,x,y,z,u,v,Q2)
end Subroutine CVfFbFFFbS

Subroutine CVfFFbFbFS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfFFbFbFS(p2,x,y,z,u,v,Q2)
end Subroutine CVfFFbFbFS

Subroutine CVfFFFbFbS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfFFFbFbS(p2,x,y,z,u,v,Q2)
end Subroutine CVfFFFbFbS

Subroutine CVfFFFFS(out,p2,x,y,z,u,v,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,v,q2
real(dp) :: out
out = VfFFFFS(p2,x,y,z,u,v,Q2)
end Subroutine CVfFFFFS

Subroutine CGfFFV(out,p2,x,y,Q2)
Implicit None
real(dp), intent(in) :: p2,x,y,q2
real(dp) :: out
out = GfFFV(p2,x,y,Q2)
end Subroutine CGfFFV

Subroutine CGfFbFbV(out,p2,x,y,Q2)
Implicit None
real(dp), intent(in) :: p2,x,y,q2
real(dp) :: out
out = GfFbFbV(p2,x,y,Q2)
end Subroutine CGfFbFbV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MASS FUNCTIONS


real(dp) Function WfSSSS (p2,x_in,y_in,z_in,u_in,Q2)
Implicit None
real(dp), intent(in):: p2,x_in,y_in,z_in,u_in,q2
real(dp) :: res,x,y,z,u

x=x_in
y=y_in
z=z_in
u=u_in


If(smalldp(x,y)) then
   res=real(-1._dp*IIp(x,z,u,q2),dp)
else
   res=Real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
end if


WfSSSS=res

end function WfSSSS


real(dp) Function XfSSS(p2,x,y,z,Q2)

Implicit None
real(dp), intent(in):: p2,x,y,z,q2
real(dp) :: res

if (smallc(z/q2)) then
   res = 0._dp
else
  res = real(-J0(z,Q2)*BB(x,y,q2),dp)
end if

XfSSS = res
end function XfSSS


real(dp) Function YfSSSS(p2,x,y,z,u,Q2)
real(dp), intent(in):: p2,x,y,z,u,q2

real(dp) :: res


if(smallc(u/q2)) then
   res = 0._dp
else
   res = real(J0(u,q2)*CC(x,y,z,q2),dp)
end if

YfSSSS=res


end Function YfSSSS


real(dp) Function ZfSSSS(p2,x,y,z,u,Q2)
real(dp), intent(in):: p2,x,y,z,u,q2

real(dp) :: res


res = real(BB(x,y,q2)*BB(z,u,q2),dp)
ZfSSSS=res
end function ZfSSSS

real(dp) Function SfSSS(p2,x,y,z,Q2)
real(dp), intent(in):: p2,x,y,z,q2

real(dp) :: res


res = -real(II(x,y,z,q2),dp)

SfSSS=res
end function SfSSS

real(dp) Function UfSSSS(p2,x,y,z,u,Q2)
real(dp), intent(in):: p2,x,y,z,u,q2

real(dp) :: res

res = UU(x,y,z,u,q2)

UfSSSS=res
end function UfSSSS

real(dp) Function VfSSSSS(p2,x_in,y_in,z_in,u_in,v_in,Q2)
real(dp), intent(in):: p2,x_in,y_in,z_in,u_in,v_in,q2

real(dp) :: res,u,v,x,y,z

if(abs(v_in) .gt. abs(u_in)) then
   v = v_in
   u = u_in
else
   v = u_in
   u = v_in
end if

x=x_in
y=y_in
z=z_in

if(smalldp(x,y)) then
       res=real(-VV0(z,x,u,v,q2),dp)
else
   if(smalldp(y,z)) then
      res=real(-VV0(x,y,u,v,q2),dp)
   else
      res=real((UU(x,y,u,v,q2)-UU(x,z,u,v,q2))/(y-z),dp)
   end if
end if



VfSSSSS=res
end function VfSSSSS

real(dp) Function MfSSSSS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res



res=-real(MM(x,y,z,u,v,q2),dp)

MfSSSSS=res
end function MfSSSSS

!!!!!!!!!!!!!!!!! SCALARS AND VECTORS


real(dp) Function WfSSSV(p2,x,Q2)
real(dp), intent(in):: p2,x,q2

real(dp) :: res

!res=real(3*II(0._dp,x,x,q2)-J0(x,q2)+2._dp*x,dp)
if(smallc(x/q2)) then
   res = 0._dp
else
   res=real(3._dp*x*(-5._dp+4._dp*log(x/q2)-(log(x/q2))**2)-J0(x,q2)+2._dp*x,dp)
end if

WfSSSV=res

end function WfSSSV

real(dp) Function MfSSSSV(p2,x_in,y_in,Q2)
real(dp), intent(in):: p2,x_in,y_in,q2

real(dp) :: res,x,y

if(abs(y_in) .gt. abs(x_in)) then
   y = y_in
   x= x_in
else
   y = x_in
   x =y_in
end if

if(smalldp(x,y)) then
   call saveme(q2,y) !!!! y > x so just maybe ..
   res = real( -1 + 5*log(y/q2) -3 *log(y/q2)**2,dp)
else
   if(smallc(x/q2)) then
     !!! NB the case that they are both small is dealt with above ...
      res=real(12 + (11*(y*log(y/q2))-3*(y*(log(y/q2))**2))/(y),dp)
   else
      res=real(-12 + (11*(x*log(x/q2)-y*log(y/q2))-3*(x*(log(x/q2))**2-y*(log(y/q2))**2))/(x-y),dp)
   end if
end if

MfSSSSV=res

end function MfSSSSV


!!!!!!!!!!!!!!!!! FERMIONS AND SCALARS

real(dp) Function WfSSFF (p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: res

res = real( &
& -(z+u-y)*UU(x,y,z,u,q2)-II(x,z,u,q2)+BB(x,y,q2)*(J0(z,q2)+J0(u,q2)) &
& ,dp)


WfSSFF=res

end function WfSSFF


real(dp) Function WfSSFbFb (p2,x,y,z,u,Q2)
Implicit None
real(dp), intent(in):: p2,x,y,z,u,q2
real(dp) :: res

res=real(-2._dp*WfSSSS(p2,x,y,z,u,Q2),dp)


WfSSFbFb=res

end function WfSSFbFb

real(dp) Function MfFbFbFbFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res



res=2._dp*real(MM(x,y,z,u,v,q2),dp)

MfFbFbFbFbS=res
end function MfFbFbFbFbS


real(dp) Function MfFFbFbFS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=real( &
& (y+z-v)*MM(x,y,z,u,v,q2)-UU(x,z,u,v,q2) - UU(u,y,x,v,q2) +BB(x,z,q2)*BB(y,u,q2) &
& ,dp)

MfFFbFbFS=res
end function MfFFbFbFS

real(dp) Function MfFFbFFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=real( &
& (x+z)*MM(x,y,z,u,v,q2)-UU(y,u,z,v,q2)-UU(u,y,x,v,q2) &
& ,dp)

MfFFbFFbS=res
end function MfFFbFFbS

real(dp) Function MfFFFbFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=real( &
& (x+y-v)*MM(x,y,z,u,v,q2)-UU(x,z,u,v,q2)-UU(y,u,z,v,q2)+BB(x,z,q2)*BB(y,u,q2) &
& ,dp)

MfFFFbFbS=res
end function MfFFFbFbS

real(dp) Function MfFFFFS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res



res=real( &
& (x*u+y*z)*MM(x,y,z,u,v,q2)-x*UU(z,x,y,v,q2)-z*UU(x,z,u,v,q2) &
& -u*UU(y,u,z,v,q2)-y*UU(u,y,x,v,q2)+II(x,u,v,q2)+II(y,z,v,q2) &
& ,dp)

MfFFFFS=res
end function MfFFFFS





real(dp) Function MfSFbSFbFb(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=2._dp*real(MM(x,y,z,u,v,Q2),dp)

MfSFbSFbFb=res
end function MfSFbSFbFb


real(dp) Function MfSFSFbF(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=real((v-x+y)*MM(x,y,z,u,v,q2)+UU(y,u,z,v,q2)-UU(x,z,u,v,q2)-BB(x,z,q2)*BB(y,u,q2),dp)

MfSFSFbF=res
end function MfSFSFbF

real(dp) Function MfSFSFFb(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res=real((y+u)*MM(x,y,z,u,v,q2)-UU(x,z,u,v,q2)-UU(z,x,y,v,q2),dp)

MfSFSFFb=res
end function MfSFSFFb

real(dp) Function VfSSSFbFb(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res = -2._dp*VfSSSSS(p2,x,y,z,u,v,q2)

VfSSSFbFb=res
end function VfSSSFbFb

real(dp) Function VfSSSFF(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res = UU(x,y,u,v,q2) +(z-u-v)*VfSSSSS(p2,x,y,z,u,v,q2)-(J0(u,Q2) +J0(v,Q2))*CC(x,y,z,q2)

VfSSSFF=res
end function VfSSSFF

real(dp) Function VfFbFbFbFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res = real(-2._dp*VfSSSSS(p2,x,y,z,u,v,q2),dp)

VfFbFbFbFbS=res
end function VfFbFbFbFbS

real(dp) Function VfFbFFbFS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res = real(-UU(x,y,u,v,q2) +(v-z-u)*VfSSSSS(p2,x,y,z,u,v,q2)-(J0(v,Q2) -J0(u,Q2))*CC(x,y,z,q2),dp)

VfFbFFbFS=res
end function VfFbFFbFS

real(dp) Function VfFbFFFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res = real(-2._dp*UU(x,y,u,v,q2) -2._dp*z*VfSSSSS(p2,x,y,z,u,v,q2),dp)

VfFbFFFbS=res
end function VfFbFFFbS

real(dp) Function VfFFbFbFS(p2,x_in,y_in,z_in,u_in,v_in,Q2)
real(dp), intent(in):: p2,x_in,y_in,z_in,u_in,v_in,q2

real(dp) :: res,x,y,z,v,u

call atvals(x,y,z,u,v,x_in,y_in,z_in,u_in,v_in)

call saveme(q2,x,y,z,u)


res=Fffs200(x,y,z,u,v,q2)

VfFFbFbFS=res
end function VfFFbFbFS

real(dp) Function VfFFFbFbS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res = real( -UU(x,y,u,v,q2) - UU(y,z,u,v,q2) -(x+z)*VfSSSSS(p2,x,y,z,u,v,q2),dp)

VfFFFbFbS=res
end function VfFFFbFbS

real(dp) Function VfFFFFS(p2,x,y,z,u,v,Q2)
real(dp), intent(in):: p2,x,y,z,u,v,q2

real(dp) :: res

res = Fffs100(y,z,u,v,q2)+x*Fffs200(y,x,z,u,v,q2)


VfFFFFS=res
end function VfFFFFS


!!!!!!!!!!!!!!!!! FERMIONS AND VECTORS


real(dp) Function GfFFV(p2,x_in,y_in,Q2)
real(dp), intent(in) :: p2,x_in,y_in,q2

real(dp) :: res,x,y,z,u,v
call atvals(x,y,z,u,v,x_in,y_in)

if(smalldp(x,y)) then
   if(smallc(x/q2,y/q2)) then
      res = 0._dp 
   else
      res = real(x*(20 - 32*Log(x/Q2) + 24*Log(x/Q2)**2),dp)
   end if
else

res=Real( 2._dp*(  (x+y)*(3._dp*UU(x,y,x,0._dp,q2)     &
     +3._dp*UU(x,y,y,0._dp,q2) - 5._dp*BB(x,y,q2))     &
     -3._dp*II(x,x,0._dp,q2)-3._dp*II(y,y,0._dp,q2)    &
     +5._dp*J0(x,q2)+5._dp*J0(y,q2)-8*(x+y)),dp)

end if



GfFFV = res

end function GfFFV

real(dp) Function GfFbFbV(p2,x_in,y_in,Q2)
real(dp), intent(in) :: p2,x_in,y_in,q2

real(dp) :: res,x,y,z,u,v
call atvals(x,y,z,u,v,x_in,y_in)

if(smalldp(x,y)) then
   if(smallc(x/q2,y/q2)) then
      res = 0._dp !!! since the function will be multiplied by zero anyway
   else
      res = real(-4 - 4*Log(x/Q2) + 12*Log(x/Q2)**2,dp)
   end if
else
   res=real(12._dp*UU(x,y,x,0._dp,q2)+12._dp*UU(x,y,y,0._dp,q2)-20._dp*BB(x,y,q2)-16._dp,dp) 

end if
GfFbFbV = res


end function GfFbFbV




!!!!!!!!!!!!!!!!! FUNCTIONS FOR FIRST DERIVATIVES

real(dp) Function Fffs100(x,y,u,v,Q2)
real(dp), intent(in):: x,y,u,v,q2

real(dp) :: res

res = real(BB(x,y,q2)*(J0(v,q2)-J0(u,q2))+II(x,u,v,q2)-(y+u-v)*UU(x,y,u,v,q2),dp)

Fffs100 =res
end function Fffs100

real(dp) Function Fffs001(x,y,u,v,Q2)
real(dp), intent(in):: x,y,u,v,q2

real(dp) :: res

res = real(-(x+y-u)*UU(x,y,u,v,q2)-II(x,y,v,q2)+BB(u,v,q2)*(J0(x,q2)+J0(y,q2)),dp)

Fffs001 =res
end function Fffs001



!!!!!!!!!!!!!!!!! FUNCTIONS FOR SECOND DERIVATIVES *FROM* FIRST DERIVATIVES

real(dp) Function Fffs200(x,y,z,u,v,Q2)
real(dp), intent(in):: x,y,z,u,v,q2

real(dp) :: res

res = real(CC(x,y,z,q2)*(J0(u,q2)-J0(v,q2))-UU(x,z,u,v,q2)-(y+u-v)*VfSSSSS(0._dp,x,y,z,u,v,q2),dp)

Fffs200 =res
end function Fffs200




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HELPER FUNCTIONS NOW





complex(dp) function BB(m1_in,m2_in,q2)
implicit none
real(dp), intent(in) :: m1_in,m2_in,q2
complex(dp) ::res
real(dp) :: m1,m2

m1=m1_in
m2=m2_in

If(smalldp(m1,m2)) then
   call saveme(q2,m1)
   res=-clog(m1/q2)
else
   res=(J0(m1,q2) - J0(m2,q2))/(m2-m1) 
end if

BB = res

end function BB


complex(dp) function CC(x_in,y_in,z_in,q2)        !!!!! 1/(k+x)(k+y)(k+z)
implicit none
real(dp), intent(in) :: x_in,y_in,z_in,q2
complex(dp) :: res
real(dp) :: x,y,z

x=x_in
y=y_in
z=z_in

If(smalldp(y,z)) then
   If(smalldp(x,y)) then
      call saveme(q2,x)
      res  =1/(2._dp*x)
   else
      res = (-x + y + x*clog(x/y))/(x - y)**2
   end if
else
   res=(BB(x,y,q2) - BB(x,z,q2))/(z-y)
end if

CC = res

end function CC





complex(dp) function BBp(x,y_in,q2)
implicit none
!!! Provides B(x,y')
real(dp), intent(in) :: x,y_in,q2
complex(dp) :: res
real(dp) ::y

y=y_in


If(smalldp(x,y)) then
   call saveme(q2,y)
   res = -1._dp/(2._dp*y)
else
res = (x-y-x*clog(x/y))/((x-y)**2)
end if

BBp=res
end function BBp


complex(dp) function clog(x)
 implicit none
 real(dp), intent(in) :: x

clog=log(dcmplx(x,0._dp))


end function clog





complex(dp) function J0(m1,q2) 
  implicit none
  real(dp), intent(in) ::  m1,q2
  
  if(smallc(m1/q2)) then
     J0=0._dp
  else
     J0 = m1*(clog(m1/q2)-1)
  end if

end function J0


real(dp) function II0(x,q2) 
      !!!! THIS IS I(x,x,0,q2) !!!!!!
      implicit none
      real(dp), intent(in) :: x,q2
      real(dp) :: res
      
      if(abs(x/q2).lt.epsM) then
         res = 0._dp
      else
         res = real(-x* ( (clog(x/Q2))**2 - 4._dp  * cLog(x/q2) + 5._dp ),dp)
      end if

      II0 = res

end function II0



real(dp) Function IIp(x_in,y_in,z_in,Q2)
Implicit None
!!! Calculates the D[II(....),x] 
Real(dp),Intent(in) :: x_in,y_in,z_in,Q2
real(dp) :: x,y,z
Complex(dp) :: res


x=x_in
call saveme(q2,x)


if(abs(z_in) .lt. abs(y_in)) then
   y = z_in
   z = y_in
else
   y = y_in
   z = z_in
end if 

if(smallc(y/q2)) then
   if(smallc(z/q2)) then
      res=real(-0.5_dp - Pi**2/6._dp + cLog(x/Q2) - 0.5_dp*cLog(x/Q2)**2,dp)
   else
      res=real(-CLi2_mod(dcmplx(1 - z/x,0._dp))  - 0.5_dp*(1._dp-cLog(x/Q2))**2,dp)
      
   end if
else
if(smalldp(x,y) .and. smalldp(x,z)) then
res = 0.5_dp*real(-(1._dp- cLog(x/Q2))**2 + phimod(x, x, x),dp)
else
res = real(0.5_dp*(-1._dp + clog(y/Q2)* clog(z/Q2)     &
     - cLog(x/Q2) *(-2._dp + cLog(y/Q2) + cLog(z/Q2))) &
     -0.5_dp*((x-y-z)/z)*phimod(x,y,z),dp)
end if

end if

IIp=res

end function IIp






real(dp) function UU(x,y,z,u,Q2)
implicit none
real(dp), intent(in) :: x,y,z,u,Q2
real(dp) :: res



if(smalldp(x,y)) then
res = -real(IIp(x,z,u,q2),dp)
else
res = real((II(x,z,u,q2) - II(y,z,u,q2))/(y-x),dp)
if(res.ne.res) then
   write(*,*) "Problem in UU: ",x,y,z,u,II(x,z,u,q2),II(y,z,u,q2)
   end if
end if

UU=res
end function UU



real(dp) function VV0(x_in,y_in,u_in,v_in,Q2)
implicit none
real(dp), intent(in) :: x_in,y_in,u_in,v_in,Q2
real(dp) :: res,x,y,z,u,v

integer :: section

call atvals(x,y,u,v,z,x_in,y_in,u_in,v_in)    !!!!! note called with x,y,u,v ...

if(abs(u) .gt. abs(v)) then
   u = v_in
   v = u_in
end if


if(smalldp(x,y)) then
   
   !!!!! 1/2d^2/dx^2 I(x,u,v)
   if(smallc(u/q2)) then
      if(smallc(v/q2)) then 
         !! I(x'',0,0)/2
         call saveme(q2,x) !!!! always singular if x=0
         res = real((1._dp - cLog(x/Q2))/(2._dp*x),dp)
      else
         !! I(x'',0,v)
         if(smalldp(x,v)) then
            !!! I(x'',0,x)/2 
            call saveme(q2,x)
            res = real(-cLog(x/Q2)/(2.*x),dp)
         else
            call saveme(q2,x)
            res = real(-(v*(-v + x)*cLog(v/Q2) + v**2*cLog(v/x) + (v - x)*(-x + (v + x)*cLog(x/Q2)))/(2.*(v - x)*x**2),dp)
         end if
      end if
   else
      !!! u,v are not small 
      call saveme(q2,x) 
      
      res = real((-(deltanm(x,u,v)*(-2 + cLog(u/Q2) + cLog(v/Q2))) &
           + (u + v - x)*((u - v)*cLog(v/u)                        &
           + x*cLog(x**2/(u*v))))/(4.*x*deltanm(x,u,v)),dp)
      res=res+real(u*phimod(x,u,v)/deltanm(x,u,v),dp)
      

   end if
   
else
   res = real((IIp(y,u,v,q2) + UU(x,y,u,v,q2))/(y-x),dp)
   if(res.ne.res) then
         write(*,*) "problem VV0: ",x,y,u,v,IIp(y,u,v,q2),UU(x,y,u,v,q2)
         end if

end if

VV0=res
end function VV0






real(dp) function MM(x_in,y_in,z_in,u_in,v_in,q2)
implicit none
real(dp), intent(in) :: x_in,y_in,z_in,u_in,v_in,Q2
real(dp) :: res,x,y,z,u,v



if(smalldp(x_in,z_in)) then
   ! swap x,z and y,u
  call atvals(y,x,u,z,v,x_in,y_in,z_in,u_in,v_in)
else
   call atvals(x,y,z,u,v,x_in,y_in,z_in,u_in,v_in)
end if

if(smalldp(u,y)) then
   call saveme(q2,x,y)
   if(smalldp(x,z)) then
      ! result is I(x',y',v)
      if(smalldp(x,y)) then      
         !Result is I(x',x',v)
         if(smallc(v/q2)) then
            !Result is I(x',x',0)
            res = real(1/x,dp)
         else
            res=real((v - 2._dp*x)*phimod(x,x,v)/(v*(v-4._dp*x)),dp)
            res=res+real(2._dp*clog(x/v)/(v-4._dp*x),dp)
         end if
      else
         If(smallc(v/q2)) then
            ! result is I(x',y',0) with x =!= y
            res = real(clog(x/y)/(x-y),dp)
         else
            res = real((deltanm(x,y,v)*(v*(cLog(v/x)) + y*phimod(x,y,v)) &
                 - (v - x + y)*(v*((v - x)*cLog(v/x) - y*(cLog(y/v) &
                 + cLog(y/x))) + y*(-v - x + y)*phimod(x,y,v)))/(2.*v*y*deltanm(x,y,v)),dp)
         end if
      end if
      
      
   else
      !! Exploit symmetry of the function here ...
      res=real((IIp(y,x,v,q2)-IIp(y,z,v,q2))/(x-z),dp)
   end if
else

      res=real((UU(x,z,y,v,q2) - UU(x,z,u,v,q2))/(u-y),dp)


end if

MM=res
end function MM


 complex(dp) Function phimod(x,y,z)
 !--------------------------------------------------------------
 ! from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
 ! version by Pietro Slavich
 ! 12.03.02: portation to f90
 !--------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x, y, z
  
  complex(dp) :: xd,yd,zd
  xd=x
  yd=y
  zd=z
  
  If(abs(x).Le.abs(z).And.abs(y).Le.abs(z)) Then
     phimod = myphimod(x,y,z)
  Elseif(abs(z).Le.abs(x).And.abs(y).Le.abs(x)) Then
     phimod = z/x*myphimod(z,y,x)
  Elseif(abs(z).Le.abs(y).And.abs(x).Le.abs(y)) Then
     phimod = z/y*myphimod(z,x,y)
  Endif


  if(phimod.ne.phimod) then
     write(*,*) "wtf is up with phimod? ",x,y,z
     end if

 Contains

  Real(dp) Function myphimod(x,y,z)
   Implicit None
   Real(dp), Intent(in) :: x,y,z
   complex(dp) :: u,v
   Complex(dp) :: clam,cxp,cxm,ccphi

   !     auxiliary variables
   u = x/z
   v = y/z
   If(abs(u)<=epsm) Then
     
     If(abs(v)/=1._dp) Then
       myphimod = real((Log(u)*Log(v)+2._dp*CLi2_mod(1._dp-v))/(1._dp-v),dp)
     Else
       myphimod = real(2._dp-Log(u),dp)
     Endif

     Elseif(abs(v)<=epsM) Then
     If(abs(u)/=1._dp) Then
        myphimod = real((Log(v)*Log(u)+2._dp*CLi2_mod(1._dp-u))/(1._dp-u),dp)
     Else
       myphimod = real(2._dp-Log(v),dp)
     Endif
   Else
     
       clam = Sqrt((1._dp-u-v)**2 - 4._dp*u*v)

    cxp = (1._dp+(u-v)-clam)/2._dp
    cxm = (1._dp-(u-v)-clam)/2._dp
       
     ccphi = (2._dp*Log(cxp)*Log(cxm) - Log(u)*Log(v)  &
& - 2._dp*(CLI2_mod(cxp) + CLI2_mod(cxm)) + Pi2o3)/clam
    myphimod = Real(ccphi,dp)
       
   Endif
   Return

  End Function myphimod
  
 End Function phimod  







real(dp) function deltaNM(x,y,z)
real(dp), intent(in)::  x,y,z

      deltaNM = x**2+y**2+z**2 -2._dp*(x*y+x*z+y*z)
end function deltanm


Complex(dp) Function RFunction(x,y,z)
Implicit None
Complex(dp), Intent(in) :: x,y,z
Complex(dp) :: temp
temp = x**2+y**2+z**2-2._dp*x*y-2._dp*x*z-2._dp*y*z
RFunction = sqrt(temp)

End Function RFunction




Complex(dp) Function II(x_in,y_in,z_in,Q2)
Implicit None
! x_in,y_in,z_in are supposed to be squared masses
! symmetric under interchange of any two of these
Real(dp),Intent(in) :: x_in,y_in,z_in,Q2
Complex(dp) :: res,x,y,z
x = x_in
y = y_in
z = z_in
If (Q2.le.0._dp) Then 
  write(*,*) 'invalid renormalization scale!'
  return
end if
Call SortArguments3C(x,y,z) ! sort x<=y<=z


if(smallc(real(x/q2))) then
   if(smallc(real(y/q2))) then
      if(smallc(real(z/q2))) then
         II = ZeroC ! I(0,0,0)=0
      else  !I(0,0,z)
         II = -0.5_dp*z*(log(z/Q2))**2 + 2*z*log(z/Q2)-2.5_dp*z-pi**2*z/6._dp
      end if
   else !I(0,y,z)
      if(smalldp(y,z)) then ! I(0,z,z)
         II = -z*(log(z/Q2))**2+4._dp*z*log(z/Q2)-5._dp*z
      else !I(0,y,z)
         II = (z-y)*(CLi2_mod(y/z)-log(z/y)*log((z-y)/Q2) &
   & +0.5_dp*(log(z/Q2))**2-pi**2/6._dp)  &
   & -2.5_dp*(z+y)+2*y*log(y/Q2)+2*z*log(z/Q2) & 
   & -z*log(z/Q2)*log(y/Q2)
      end if
   end if

else




res = Cmplx(0._dp,0._dp)
res=res + 0.5_dp*(x-y-z)*log(y/Q2)*log(z/Q2)
res=res + 0.5_dp*(y-x-z)*log(x/Q2)*log(z/Q2)
res=res + 0.5_dp*(z-x-y)*log(x/Q2)*log(y/Q2)
res=res + 2._dp*x*log(x/Q2)
res=res + 2._dp*y*log(y/Q2)
res=res + 2._dp*z*log(z/Q2)
res=res - 2.5_dp*(x+y+z)
res=res - 0.5*XiFunction(x,y,z)
II = res
if(res.ne. res) then
   write(*,*) "Problem in general II: ",x,y,z,xifunction(x,y,z)
end if

End If

End Function II








Complex(dp) Function XiFunction(x0,y0,z0)
Implicit None
Complex(dp), Intent(in) :: x0,y0,z0
Complex(dp) :: x,y,z
Complex(dp) :: res,R
x=x0
y=y0
z=z0
Call SortArguments3C(x,y,z)
R=Rfunction(x,y,z)
res = Cmplx(0.,0.,dp)
res = res + 2._dp*log((z+x-y-R)/(2._dp*z))*log((z+y-x-R)/(2._dp*z))
res = res - log(x/z)*log(y/z)
res = res - 2._dp*CLi2_mod((z+x-y-R)/(2._dp*z))
res = res - 2._dp*CLi2_mod((z+y-x-R)/(2._dp*z))
res = res + (pi**2)/3._dp
res = R*res
XiFunction = res
End Function XiFunction

!! sorts three real arguments a,b,c from small to large
!! and returns them
Subroutine SortArguments3(a,b,c)
Implicit None
Real(dp), Intent(inout) :: a,b,c
Real(dp) :: x,y,z
x=a ! assume a is smallest
y=b
z=c
If (b.le.a) Then
 x=b
 y=a
 z=c
End If
If (c.le.x) Then
 x=c
 y=a
 z=b
End If
! compare last two elements
If (z.le.y) Then
 a=x ! output
 b=z
 c=y
Else
 a=x !output
 b=y
 c=z
End If
End Subroutine SortArguments3




Subroutine SortArguments3C(a,b,c)
Implicit None
Complex(dp), Intent(inout) :: a,b,c
Complex(dp) :: x,y,z
x=a ! assume a is smallest
y=b
z=c
If (Abs(b).le.Abs(a)) Then
 x=b
 y=a
 z=c
End If
If (Abs(c).le.Abs(x)) Then
 x=c
 y=a
 z=b
End If
! compare last two elements
If (Abs(z).le.Abs(y)) Then
 a=x ! output
 b=z
 c=y
Else
 a=x !output
 b=y
 c=z
End If
End Subroutine SortArguments3C


 Complex(DP) Function CLI2_mod(Z)
 !--------------------------------------------------------------
 !     routine for the complex dilogarithm function
 !--------------------------------------------------------------
 Implicit None      
  !-------
  ! input
  !-------
  Complex(DP), Intent(in) :: Z

  !-----------------
  ! local variables
  !-----------------
  Complex(DP) :: sumC, cy, cz, cz2, CLI2
  Real(dp) :: rz, iz, az
  Integer :: i1, jsgn, ipi12
  !---------------------------------------------
  ! Bernulli numbers / (2 n + 1)!
  !---------------------------------------------
  Real(dp), Parameter :: bf(20) = (/ - 1._dp/4._dp, + 1._dp/36._dp        &
       & , - 1._dp/36.e2_dp, + 1._dp/21168.e1_dp,  - 1._dp/108864.e2_dp   &
       & , + 1._dp/52690176.e1_dp,  - 691._dp/16999766784.e3_dp           &
       & , + 1._dp/1120863744.e3_dp,  - 3617._dp/18140058832896.e4_dp     &
       & , + 43867._dp/97072790126247936.e3_dp                            &
       & , - 174611._dp/168600109166641152.e5_dp                          &
       & , + 77683._dp/32432530090601152512.e4_dp                         &
       & , - 236364091._dp/4234560341829359173632.e7_dp                   &
       & , + 657931._dp/5025632054039239458816.e6_dp                      &
       & , - 3392780147._dp/109890470493622010006470656.e7_dp             &
       & , +172.3168255201_dp/2355349904102724211909.3102313472e6_dp      &
       & , -770.9321041217_dp/4428491985594062112714.2791446528e8_dp      &
       & ,  4.1576356446138997196178996207752266734882541595116e-29_dp    &
       & , -9.9621484882846221031940067024558388498548600173945e-31_dp    &
       & ,  2.3940344248961653005211679878937495629342791569329e-32_dp    &
       & /)

  rz = Real(z,dp)
  iz = Aimag(z)
  az = Sqrt(rz**2 + iz**2)
  !------------------
  ! special cases
  !------------------
  If (iz.Eq.0._dp) Then
    If (rz.gt.1._dp) Then
    CLi2_mod = Cmplx( Li2(rz), -pi*log(rz),dp)
    Else
    CLi2_mod = Cmplx( Li2(rz), 0.,dp)
    End If
   Return
  Else If (az.Lt.Epsilon(1._dp)) Then
   CLi2_mod = z
   Return
  End If
  !-----------------------------------------------------
  ! transformation to |z|<1, Re(z)<=0.5
  !-----------------------------------------------------
  If (rz.Le.0.5_dp) Then
   If (az.Gt.1._dp) Then
    cy = -0.5_dp * Log(-z)**2
    cz = - Log(1._dp - 1._dp / z)
    jsgn = -1
    ipi12 = -2
   Else !  (az.le.1._dp)
    cy = 0
    cz = - Log(1._dp - z)
    jsgn = 1
    ipi12 = 0
   End If
  Else ! rz.gt.0.5_dp
   If (az.Le.Sqrt(2*rz)) Then
    cz = - Log(z)
    cy = cz * Log(1._dp - z)
    jsgn = -1
    ipi12 = 2
   Else ! (az.gt.Sqrt(2*rz))
    cy = -0.5_dp * Log(-z)**2
    cz = - Log(1._dp - 1._dp / z)
    jsgn = -1
    ipi12 = -2
   End If
  End If
  !--------------------------------------------
  ! the dilogarithm
  !--------------------------------------------
  cz2 = cz**2
  sumC = cz2 * bf(20)
  Do i1=3,19
   sumC = cz2 * ( sumC + bf(22-i1) )
  End Do
  ! watch the powers of z
  sumC = cz + cz2 * (bf(1) + cz * (bf(2) + sumC))
  If (jsgn.Eq.1) Then
   CLi2 = sumC + cy
  Else
   CLi2 = - sumC + cy
  End If
  If (ipi12.Ne.0) Cli2 = Cli2 + ipi12 * pi2o12
  
  CLI2_mod = CLI2
 End Function CLI2_mod




Logical Function SmalldpC(x,y) 
Implicit None
Complex(dp), Intent(in) :: x,y
If (Max(Abs(x),Abs(y)).lt.epsD) Then
 SmalldpC = .True.
Else If ((Abs(x-y)/Max(Abs(x),Abs(y))).le.epsD) Then
 SmalldpC = .True.
Else
 SmalldpC = .False.
End if
End Function SmalldpC


Logical Function SmalldpR(x,y) 
Implicit None
Real(dp), Intent(in) :: x,y
If (Max(Abs(x),Abs(y)).lt.epsD) Then
 SmalldpR = .True.
Else If ((Abs(x-y)/Max(Abs(x),Abs(y))).le.epsD) Then
 SmalldpR = .True.
Else
 SmalldpR = .False.
End if
End Function SmalldpR






Logical Function smallc(x,y,z,u,v)
Implicit None
real(dp), Intent(in) :: x
real(dp), Intent(in),optional :: y,z,u,v

Logical :: res
res = .False.
smallc=.False.

If(abs(x) .lt. epsM) then
   smallc = .True.
else 
   if(present(y)) then
      if(abs(y) .lt. epsM) then
         smallc = .True.
      end if
   end if
   if(present(z)) then
      if(abs(z) .lt. epsM) then
         smallc = .True.
      end if
   end if
   if(present(u)) then
      if(abs(u) .lt. epsM) then
         smallc = .True.
      end if
   end if
   if(present(v)) then
      if(abs(v) .lt. epsM) then
         smallc = .True.
      end if
   end if
end if





end function smallc


subroutine saveme(q2,x,y,z,u,v)
implicit none

real(dp), intent(in) ::q2
real(dp), Intent(inout) :: x
real(dp), intent(inout), optional :: y,z,u,v

real(dp) :: minm


minm=epsM*Q2

If(abs(x) .lt. minm) then
   x = minm
end if

if(present(y)) then
   if(abs(y) .lt. minM) then
       y = minM
   end if
end if
if(present(z)) then
   if(abs(z) .lt. minM) then
       z = minM
   end if
end if
if(present(u)) then
   if(abs(u) .lt. minM) then
       u = minM
   end if
end if
if(present(v)) then
   if(abs(v) .lt. minM) then
       v = minM
   end if
end if



end subroutine saveme

subroutine atvals(x,y,z,u,v,x_in,y_in,z_in,u_in,v_in)
implicit none

real(dp), intent(in), optional :: x_in,y_in,z_in,u_in,v_in
real(dp), intent(out) :: x,y,z,u,v

x=x_in

If(present(y_in)) then
   y=y_in
else
   y=0._dp
end if
If(present(z_in)) then
   z=z_in
else
   z=0._dp
end if
If(present(u_in)) then
   u=u_in
else
   u=0._dp
end if
If(present(v_in)) then
   v=v_in
else
   v=0._dp
end if

end subroutine atvals

 Real(Dp) Function Li2(X)
 !---------------------------------------------------------------------
 ! Spence-fkt. li2:
 ! local double prec. Version of function dilogarithm, 22. 5. 95
 ! H. Eberl, Convention - see Mathematica-manual side 575
 ! checked
 ! 1.10.2000: change to f90 + implicit none statement
 !---------------------------------------------------------------------
 Implicit None

  Real(Dp), Intent(in) :: x
  Real(Dp), Parameter :: Z1 = 1._dp, HF = Z1/2._dp,              &
        & C(0:19) = (/ 0.42996693560813697_dp,  &
        &              0.40975987533077105_dp,  &
        &             -0.01858843665014592_dp,  &
        &              0.00145751084062268_dp,  &
        &             -0.00014304184442340_dp,  &
        &              0.00001588415541880_dp,  &
        &             -0.00000190784959387_dp,  &
        &              0.00000024195180854_dp,  &
        &             -0.00000003193341274_dp,  &
        &              0.00000000434545063_dp,  &
        &             -0.00000000060578480_dp,  &
        &              0.00000000008612098_dp,  &
        &             -0.00000000001244332_dp,  &
        &              0.00000000000182256_dp,  &
        &             -0.00000000000027007_dp,  &
        &              0.00000000000004042_dp,  &
        &             -0.00000000000000610_dp,  &
        &              0.00000000000000093_dp,  &
        &             -0.00000000000000014_dp,  &
        &              0.00000000000000002_dp  /)
  Real(Dp) :: H,T,Y,S,A,Alfa,B0,B1,B2
  Integer :: i

  If(X .Eq. 1) Then
   H=Pi2o6
  Elseif(X .Eq. -1) Then
   H=-Pi2o12
  Else
   T=-X
   If(T .Le. -2) Then
    Y=-1/(1+T)
    S=1
    A=-Pi2o3+HF*(Log(-T)**2-Log(1+1/T)**2)
   Elseif(T .Lt. -1) Then
    Y=-1-T
    S=-1
    A=Log(-T)
    A=-Pi2o6+A*(A+Log(1+1/T))
   Else If(T .Le. -HF) Then
    Y=-(1+T)/T
    S=1
    A=Log(-T)
    A=-Pi2o6+A*(-HF*A+Log(1+T))
   Else If(T .Lt. 0) Then
    Y=-T/(1+T)
    S=-1
    A=HF*Log(1+T)**2
   Else If(T .Le. 1) Then
    Y=T
    S=1
    A=0
   Else
    Y=1/T
    S=-1
    A=Pi2o6+HF*Log(T)**2
   Endif
   H=Y+Y-1
   ALFA=H+H
   B1=0
   B2=0
   Do I = 19,0,-1
    B0=C(I)+ALFA*B1-B2
    B2=B1
    B1=B0
   Enddo
   H=-(S*(B0-H*B2)+A)
  Endif
  LI2=H

 End Function Li2


End module Pole2LFunctions
