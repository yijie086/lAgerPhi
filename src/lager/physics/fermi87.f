!>  ELECTROPRODUCTION YIELDS OF NUCLEONS AND PIONS 
!!  WRITTEN BY J.S. O'CONNELL AND J.W. LIGHTBODY, JR.
!!  NATIONAL BUREAU OF STANDARDS 
!!  SEPTEMBER 1987 
!!
!!  HISTORY:
!!     - transverse scaling region added
!!     - modified by oara to plot like older epc's
!!     - modified slightly by glen warren to compile under linux (g77) sept. 02
!!     - Selected modified subroutines and functions so it can be compiled and linked 
!!         with exisiting codes. Whitney Armstrong 1/10/2014
!!
!! @ingroup EPCV      
!!

*(formerly fermi3) fermi87
      !> Fermi distributions of nucleons 
      !!
      !! @param p     momentum in MeV
      !! @param ia    A=Z+N
      !! @param res   ouput
      !! 
      !! Revision Histroy:
      !!   - Added copy from sgsl_3 function. 1/25/2014 Whitney Armstrong
      !!
      !!  @ingroup EPCV
      !!
      subroutine fermi87(p,ia,res) 
      implicit DOUBLE PRECISION (a-h,o-z) 
c      INTEGER ia
c  p integral over sgsl normalized to 1/4pi 
      if(ia.eq.2)then 
c  begin 2-h
         pp=p/197.3 
         sgs=3.697-7.428*pp-2.257*pp**2 
         sgs=sgs+3.618*pp**3-1.377*pp**4+.221*pp**5-.013*pp**6
         if(sgs.lt.-293.)go to 1
         sgs=exp(sgs) 
         sgs=sgs/.18825/4./3.1416/(197.3)**3
         sgsl_3=sgs/1.
      elseif(ia.eq.3)then 
c  begin 3-he 
         if(-(p/33)**2.lt.-293.)go to 1 
         sgs=2.4101e-6*exp(-p/33) 
         sgs=sgs-1.4461e-6*exp(-(p/33)**2)
         sgs=sgs+1.6871e-10*exp(-(p/493)**2)
         sgsl_3=sgs/2.0d0
      elseif(ia.eq.4)then 
c   begin 4-he
         if(-(p/113.24)**2.lt.-293.)go to 1 
         sgs=1.39066e-6*exp(-(p/113.24)**2) 
         sgs=sgs+3.96476e-9*exp(-(p/390.75)**2) 
         sgsl_3=sgs/2.0d0
         sgsl_3=sgsl_3/2.0d0/3.1416
      elseif(ia.gt.4.and.ia.lt.12)then
         if(-(p/127)**2.lt.-293.)go to 1
         sgs=1.7052e-7*(1.+(p/127)**2)*exp(-(p/127)**2) 
         sgs=sgs+1.7052e-9*exp(-(p/493)**2) 
         sgsl_3=sgs/(float(ia)/2.0d0)
      elseif(ia.eq.12)then
c  begin 12-c 
         if(-(p/127)**2.lt.-293.)go to 1
         sgs=1.7052e-7*(1.+(p/127)**2)*exp(-(p/127)**2) 
         sgs=sgs+1.7052e-9*exp(-(p/493)**2) 
         sgsl_3=sgs/6.
      else
c  begin 16-o 
         if(-(p/120)**2.lt.-293.)go to 1
         sgs=3.0124e-7*(1.+(p/120)**2)*exp(-(p/120)**2) 
         sgs=sgs+1.1296e-9*exp(-(p/493)**2) 
         sgsl_3=sgs/(float(ia)/2.0d0)
      endif 
      res=sgsl_3
      return
    1 sgsl_3=0. 
      res=sgsl_3
      return
      end 
*sgsl_3 
      real*8 function sgsl_3(p) 
      implicit real*8 (a-h,o-z) 
c  p integral over sgsl normalized to 1/4pi 
      common/sp/ia
      common/qd/qdf 
      if(ia.eq.2)then 
c  begin 2-h
         pp=p/197.3 
         sgs=3.697-7.428*pp-2.257*pp**2 
         sgs=sgs+3.618*pp**3-1.377*pp**4+.221*pp**5-.013*pp**6
         if(sgs.lt.-293.)go to 1
         sgs=exp(sgs) 
         sgs=sgs/.18825/4./3.1416/(197.3)**3
!         sgsl=sgs/1.
      elseif(ia.eq.3)then 
c  begin 3-he 
         if(-(p/33)**2.lt.-293.)go to 1 
         sgs=2.4101e-6*exp(-p/33) 
         sgs=sgs-1.4461e-6*exp(-(p/33)**2)
         sgs=sgs+1.6871e-10*exp(-(p/493)**2)
         sgsl_3=sgs/2.0d0
      elseif(ia.eq.4)then 
c   begin 4-he
         if(-(p/113.24)**2.lt.-293.)go to 1 
         sgs=1.39066e-6*exp(-(p/113.24)**2) 
         sgs=sgs+3.96476e-9*exp(-(p/390.75)**2) 
         sgsl_3=sgs/2.0d0
         sgsl_3=sgsl_3/2.0d0/3.1416
      elseif(ia.gt.4.and.ia.lt.12)then
         if(-(p/127)**2.lt.-293.)go to 1
         sgs=1.7052e-7*(1.+(p/127)**2)*exp(-(p/127)**2) 
         sgs=sgs+1.7052e-9*exp(-(p/493)**2) 
         sgsl_3=sgs/(float(ia)/2.0d0)
      elseif(ia.eq.12)then
c  begin 12-c 
         if(-(p/127)**2.lt.-293.)go to 1
         sgs=1.7052e-7*(1.+(p/127)**2)*exp(-(p/127)**2) 
         sgs=sgs+1.7052e-9*exp(-(p/493)**2) 
         sgsl_3=sgs/6.
      else
c  begin 16-o 
         if(-(p/120)**2.lt.-293.)go to 1
         sgs=3.0124e-7*(1.+(p/120)**2)*exp(-(p/120)**2) 
         sgs=sgs+1.1296e-9*exp(-(p/493)**2) 
         sgsl_3=sgs/(float(ia)/2.0d0)
      endif 
      return
    1 sgsl_3=0. 
      return
      end 

