*     Two-loop strong corrections to the Higgs mass matrices in the NMSSM
*
*     code written by Pietro Slavich (slavich@lpthe.jussieu.fr)
*
*     Based on arXiv:0907.4682 by Giuseppe Degrassi and Pietro Slavich
*
*     Latest modifications on 06/07/2010
*
*************************************************************************


      subroutine effpot(lp,mt,mg,T1,T2,st,ct,q2,tanb,vv,l,xx,as,
     $     DMS,DMP) 
     $ bind(C, name="effpot_")

      implicit none

      integer i,j,lp
      double precision mt,mg,T1,T2,st,ct,q2,tanb,vv,l,xx,as,
     $     DMS(0:2,0:2),DMP(0:2,0:2)
      double precision c2t,s2t,At,mu,Xt,ht,sbe,pi,k
      double precision F1t,F2t,F3t,Ft,FA


      pi = 4d0*atan(1.d0)

      if(lp.eq.1) then
         k = 3d0/(16d0*Pi**2)   ! one-loop factor
      elseif(lp.eq.2) then
         k = as/(16d0*Pi**3)    ! two-loop factor
      else
         k = 0d0
      endif

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      mu = l*xx
     
      Xt = (T1-T2)*s2t/2d0/mt    
      At = Xt + mu/tanb         

      sbe = dsin(datan(tanb))

      ht = mt/vv/sbe            ! v ~ 174
 

      if(lp.eq.1) then          !the usual one-loop functions
         
         Ft = T1*(log(T1/q2)-1)-T2*(log(T2/q2)-1)
         
         F1t = log(T1*T2/mt**4)
         
         F2t = log(T1/T2)
         
         F3t = 2-(T1+T2)/(T1-T2)*log(T1/T2)

         FA = Ft
  
      elseif(lp.eq.2) then
         
         call makefuncsNM(mt,mg,T1,T2,s2t,c2t,q2,tanb,At,mu,
     $        F1t,F2t,F3t,Ft,FA)

      endif

      
c     now build up the results


      DMS(0,0) = .5d0*ht**2*mu**2*s2t**2*F3t
     $     + ht**2*tanb*mu*At/(T1-T2)*Ft

      DMS(0,1) = - ht**2*mu*mt*s2t*F2t - .5d0*ht**2*At*mu*s2t**2*F3t
     $     - ht**2*mu*At/(T1-T2)*Ft

      DMS(0,2) = .5d0*ht*l*mu*mt*s2t**2/tanb*F3t
     $     - ht*l*mt*(At - 2*mu/tanb)/(T1-T2)*Ft

      DMS(1,1) = 2d0*ht**2*mt**2*F1t + 2d0*ht**2*At*mt*s2t*F2t
     $     + .5d0*ht**2*At**2*s2t**2*F3t
     $     + ht**2/tanb*mu*At/(T1-T2)*Ft
      
      DMS(1,2) = -.5d0*ht*l*At*mt*s2t**2/tanb*F3t
     $     - ht*l*mt**2*s2t/tanb*F2t - ht*l*mt*At/(T1-T2)/tanb*Ft

      DMS(2,2) = .5d0*l**2*s2t**2*mt**2/tanb**2*F3t
     $     +l**2*mt**2/tanb*At/mu/(T1-T2)*Ft
      
      DMS(2,1) = DMS(1,2)
      DMS(1,0) = DMS(0,1)
      DMS(2,0) = DMS(0,2)

      DMP(0,0) = ht**2*mu*At/(T1-T2)*FA*tanb

      DMP(0,1) = ht**2*mu*At/(T1-T2)*FA

      DMP(0,2) = l*ht*mt*At/(T1-T2)*FA

      DMP(1,1) = ht**2*mu*At/(T1-T2)*FA/tanb
      
      DMP(1,2) = l*ht*mt*At/(T1-T2)*FA/tanb
      
      DMP(2,2) = l**2*mt**2*At/mu/(T1-T2)*FA/tanb

      DMP(1,0) = DMP(0,1)
      DMP(2,0) = DMP(0,2)
      DMP(2,1) = DMP(1,2)


      do i=0,2
         do j=0,2
            DMS(i,j) = k*DMS(i,j)
            DMP(i,j) = k*DMP(i,j)
         enddo
      enddo

      return
      end

*
***********************************************************************
*
      
      subroutine makefuncsNM(mt,mg,T1,T2,s2t,c2t,q,tanb,At,mu,
     $     F1t,F2t,F3t,Ft,FA)

      implicit none

      double precision mt,mg,T1,T2,s2t,c2t,q,tanb,At,mu,
     $     F1t,F2t,F3t,Ft,FA
      
      double precision DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,
     $     DT1t,DT2t,DT1T2,Dtc2t,DT1c2t,DT2c2t,Dcptmptt,
     $     Dtt_1,Dc2t_1,Dc2tc2t_1,Dtc2t_1,Dcptmptt_1,
     $     Dtt_2,Dc2t_2,Dc2tc2t_2,Dtc2t_2,Dcptmptt_2

      call makederivNM(mt,mg,T1,T2,s2t,c2t,q,
     $     DT1,DT1T1,DT1t,DT1c2t,DT1T2,
     $     Dtt_1,Dc2t_1,Dc2tc2t_1,Dtc2t_1,Dcptmptt_1)

      call makederivNM(mt,mg,T2,T1,-s2t,c2t,q,
     $     DT2,DT2T2,DT2t,DT2c2t,DT1T2,
     $     Dtt_2,Dc2t_2,Dc2tc2t_2,Dtc2t_2,Dcptmptt_2)

      Dtt = Dtt_1 + Dtt_2
      Dc2t = Dc2t_1 + Dc2t_2
      Dc2tc2t = Dc2tc2t_1 + Dc2tc2t_2
      Dtc2t = Dtc2t_1 + Dtc2t_2
      Dcptmptt = Dcptmptt_1 + Dcptmptt_2

      F1t = Dtt + DT1T1 + DT2T2 + 2d0*(DT1t + DT2t + DT1T2)

      F2t = DT1T1 - DT2T2 + DT1t - DT2t
     $     -4d0*c2t**2/(T1-T2)*(Dtc2t + DT1c2t + DT2c2t)
      
      F3t = DT1T1 + DT2T2 - 2d0*DT1T2
     $     - 2d0/(T1-T2)*(DT1-DT2)
     $     + 16d0*c2t**2/(T1-T2)**2*(c2t**2*Dc2tc2t + 2d0*Dc2t)
     $     -8d0*c2t**2/(T1-T2)*(DT1c2t-DT2c2t)

      Ft = DT1 - DT2 - 4d0*c2t**2/(T1-T2)*Dc2t      
      
      FA = Ft - 2*mu/tanb/At/(T1-T2)/s2t**2*Dcptmptt
      
      return
      end

*
***********************************************************************
*

      subroutine makederivNM(mt,mg,T1,T2,s2t,c2t,q,
     $     DT1,DT1T1,DT1t,DT1c2t,DT1T2,Dtt,Dc2t,Dc2tc2t,Dtc2t,Dcptmptt)
      
      implicit none

      double precision mt,mg,T1,T2,s2t,c2t,q,
     $     DT1,DT1T1,DT1t,DT1c2t,DT1T2,Dtt,Dc2t,Dc2tc2t,Dtc2t,Dcptmptt
      double precision deltNM,phiNM,II,JJ

      double precision t,g,Logt,Logg,LogT1,LogT2,pphiNM,del,III

      t = mt**2
      g = mg**2
      
      Logt = Log(t/q)
      Logg = Log(g/q)
      LogT1 = Log(T1/q)
      LogT2 = Log(T2/q)
      pphiNM = phiNM(T1,g,t)
      del = deltNM(T1,g,t)
      III = II(q,T1,g,t)
      
      Dc2t = .5d0*JJ(q,T1,T1) - .5d0*JJ(q,T1,T2) + 2*mg*mt/s2t*III

      Dc2tc2t = mg*mt/s2t**3*III
      
      Dcptmptt = -4*mg*mt*s2t*III

      DT1 = -6*T1+2*mg*mt*s2t+4*t*(1-logt+logT1)+4*g*(1-logg+logT1)
     $     +((5-c2t**2)*T1-s2t**2*T2-4*mg*mt*s2t)*logT1
     $     +(-3+c2t**2)*T1*logT1**2+s2t**2*T2*logT1*logT2
     $     -(2*(g+t-T1)-2*mg*mt*s2t)*(logt*(logT1-logg)+logT1*logg)
     $     +(2/t*(del+2*g*t)-2*mg/mt*s2t*(g+t-T1))*pphiNM

      DT1T1= -(1+c2t**2)+4/T1*(g+t-mg*mt*s2t)-s2t**2*T2/T1*(1-logT2)
     $     +(3+c2t**2+8*g*t/del-4*mg*mt*s2t/del*(g+t-T1))*logT1
     $     -4*t/del/T1*(del-g*(g-t-T1)+mg*mt*s2t*(g-t+T1))*logt
     $     -4*g/del/T1*(del+t*(g-t+T1)-mg*mt*s2t*(g-t-T1))*logg
     $     +(-3+c2t**2)*logT1**2+2*(logt*(logT1-logg)+logT1*logg)
     $     -2/t/del*((g+t-T1)*(del-2*g*t)+4*mg**3*mt**3*s2t)*pphiNM

      DT1c2t = (T2*(1-logT2)-T1*(1-logT1))*LogT1
     $     -mg*mt/s2t*(1-2*logT1+logt*(logT1-logg)
     $     +logT1*logg-(g+t-T1)/t*pphiNM)

      DT1t = mg/mt*s2t
     $     +4*g/del*(T1-g-t+2*mg*mt*s2t)*logg
     $     +4/del*(2*g*t-mg*mt*s2t*(g+t-T1))*logt
     $     +2/del*(2*g*(g-t-T1)-mg/mt*s2t*(del-2*t*(t-g-T1)))*logT1
     $     +(-2+mg/mt*s2t)*(logt*(logT1-logg)+logT1*logg)
     $     +1/del/t*(mg/mt*s2t*(del*(T1-g-3*t)+2*t*((t-T1)**2-g**2))
     $     +2*(g-T1)**3+2*t*(del+(2*T1-t)*(g+T1)))*pphiNM

      DT1T2 = s2t**2*logT1*logT2
      
      Dtt = -2-5/2d0*mg/mt**3*s2t*T1+6*logt**2
     $     +4*g/del*(g-t-T1+mt/mg*s2t*(t-g-T1))*logt
     $     -4*g/del*(g-t+T1+mg/mt*s2t*(t-g+T1))*logg
     $     +(8*g*T1/del+2*mg/mt**3*s2t*T1*(1-2*t/del*(t+g-T1)))*logT1
     $     -(2-mg*s2t*T1/2d0/mt**3)*logg*(logt-logT1)
     $     -(2+mg*s2t*T1/2d0/mt**3)*logt*logT1
     $     +mg*s2t/2d0/mt**3*(g+3*t)*logT1*(logt-logg)
     $     -2/del/t*(mg/mt**3*s2t*(del**2/4d0+t*(g-2*t+T1)*del
     $     +t**2*(g-t+T1)**2)-T1*(del+(g+t)*(2*t-T1))-(g-t)**3)*pphiNM

      Dtc2t = -mg/mt/s2t/2d0*(5*T1-4*T1*logT1+(g-3*t)*logT1*(logg-logt)
     $     +T1*(logt*(logT1-logg)+logg*logT1)+(del/t-2*(g-t+T1))*pphiNM)

      return
      end

*
**********************************************************************
*

      double precision function JJ(q,m1,m2) 
      
      implicit none
      double precision q,m1,m2

      JJ = m1*m2*(Log(m1/q)-1)*(Log(m2/q)-1)

      return
      end

*
**********************************************************************
*

      double precision function II(q,m1,m2,m3) 
      
      implicit none
      double precision q,m1,m2,m3,deltNM,phiNM

      II = (m1-m2-m3)/2d0*Log(m2/q)*Log(m3/q)
     $     +(m2-m1-m3)/2d0*Log(m1/q)*Log(m3/q)
     $     +(m3-m1-m2)/2d0*Log(m1/q)*Log(m2/q)
     $     +2*m1*log(m1/q)+2*m2*Log(m2/q)+2*m3*Log(m3/q)
     $     -2.5d0*(m1+m2+m3)-deltNM(m1,m2,m3)/2d0/m3*phiNM(m1,m2,m3)

      return
      end



*
**********************************************************************
*

      function phiNM(x,y,z)

c     from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23

      implicit none
      double precision x,y,z,phiNM,pphiNM,myphiNM
      
      if(x.le.z.and.y.le.z) then
         pphiNM = myphiNM(x,y,z)
      elseif(z.le.x.and.y.le.x) then
         pphiNM = z/x*myphiNM(z,y,x)
      elseif(z.le.y.and.x.le.y) then
         pphiNM = z/y*myphiNM(z,x,y)
      endif

      phiNM = pphiNM
      
      end
      
      function myphiNM(x,y,z)
      
      implicit none

      double precision x,y,z,myphiNM
      double precision u,v
      double precision Pi,SLLi2
      double complex clam,cxp,cxm,SLCLI2,ccphiNM

      parameter (pi = 3.1415926535897932384626433832795029D0)

c     auxiliary variables

      u = x/z
      v = y/z
      
      if(u.le.1d-8) then
         
         if(v.ne.1d0) then
            myphiNM = (log(u)*log(v)+2d0*SLLi2(1d0-v))/(1d0-v)
         else
            myphiNM = 2d0-log(u)
         endif

      elseif(v.le.1d-8) then

         if(u.ne.1d0) then
            myphiNM = (log(v)*log(u)+2d0*SLLi2(1d0-u))/(1d0-u)
         else
            myphiNM = 2d0-log(v)
         endif

      else
         
         if((1d0-u-v)**2.ge.4d0*u*v) then         
            clam = DCMPLX(sqrt((1d0-u-v)**2 - 4d0*u*v),0d0)
         else
            clam = DCMPLX(0d0,sqrt(4d0*u*v - (1d0-u-v)**2))
         endif
         
         cxp = (1d0+(u-v)-clam)/2d0
         cxm = (1d0-(u-v)-clam)/2d0
         
c     phiNM function from eq. (A4)
            
         ccphiNM = (2d0*log(cxp)*log(cxm) - log(u)*log(v) - 
     &        2d0*(SLCLI2(cxp) + SLCLI2(cxm)) + Pi**2/3d0)/clam
         myphiNM = DBLE(ccphiNM)
                     
      endif
      
      return
      end

*
***********************************************************************
*

      function SLLi2(x)

      implicit none

      double complex SLCLI2,z
      double precision x,SLLi2

      z = DCMPLX(x,0d0)
      SLLi2 = DBLE(SLCLI2(z))

      return
      end

*
***********************************************************************
*

      DOUBLE COMPLEX FUNCTION SLCLI2(Z)

c     just call the Dilog routine
      
      DOUBLE COMPLEX Z,Dilog

      SLCLI2 = Dilog(Z)

      return
      end

*
**********************************************************************
*
* Dilog.F
* complex dilogarithm
* this file is part of FeynHiggs
* last modified 20 Oct 05 th


      double complex function Dilog(z)
      implicit none
      double complex z
      
      double complex Dilogsum
      external Dilogsum
      
      double precision absz, abs1z
      double complex t, mlogz
      
      double precision pi, zeta2
      parameter (pi = 3.1415926535897932384626433832795029D0)
      parameter (zeta2 = pi*pi/6D0)
      
      absz = abs(z)
      if( absz .lt. 1D-20 ) then
         Dilog = -log(1 - z)
         return
      endif
      
      abs1z = abs(1 - z)
	if( abs1z .lt. 1D-20 ) then
           Dilog = zeta2
           return
	endif
        
	if( DBLE(z) .gt. .5D0 ) then
           mlogz = -log(z)
           t = zeta2 + mlogz*log(1 - z)
           if( abs1z .gt. 1 ) then
              Dilog = Dilogsum(log(1 - 1/z)) + zeta2 +
     $             .5D0*log(z - 1)**2 + t
           else
	    Dilog = -Dilogsum(mlogz) + t
         endif
      else
         if( absz .gt. 1 ) then
	    Dilog = -Dilogsum(-log(1 - 1/z)) - zeta2 - .5D0*log(-z)**2
         else
	    Dilog = Dilogsum(-log(1 - z))
         endif
      endif
      end
      

************************************************************************

      double complex function Dilogsum(w)
      implicit none
      double complex w
      
      double complex u, t
      integer k
      
      double precision b2, b4, b6, b8, b10, b12, b14
      double precision b16, b18, b20, b22
      parameter (b2 = 1/6D0)
      parameter (b4 = -1/30D0)
      parameter (b6 = 1/42D0)
      parameter (b8 = -1/30D0)
      parameter (b10 = 5/66D0)
      parameter (b12 = -691/2730D0)
      parameter (b14 = 7/6D0)
      parameter (b16 = -3617/510D0)
      parameter (b18 = 43867/798D0)
      parameter (b20 = -174611/330D0)
      parameter (b22 = 854513/138D0)
      
      double precision bernoulliB(11)
      data bernoulliB /b2, b4, b6, b8, b10, b12, b14,
     &     b16, b18, b20, b22/
      
      Dilogsum = w*(1 - .25D0*w)
      if( abs(w) .lt. 1D-10 ) return
      
      u = w
      do k = 1, 11
         u = u*w**2/DBLE(2*k*(2*k + 1))
         t = u*bernoulliB(k)
         Dilogsum = Dilogsum + t
         if( abs(t) .lt. 1D-16*abs(Dilogsum) ) return
      enddo

      end
      
      function deltNM(x,y,z)
      double precision deltNM,x,y,z

      deltNM = x**2+y**2+z**2 -2*(x*y+x*z+y*z)

      return
      end

c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      subroutine swap12NM(M)

      implicit none 
      double precision M(3,3),temp

      temp = M(1,1)
      M(1,1) = M(2,2)
      M(2,2) = temp

      temp = M(1,3)
      M(1,3) = M(2,3)
      M(2,3) = temp

      temp = M(3,1)
      M(3,1) = M(3,2)
      M(3,2) = temp

      return
      end

