c   July 29, 2005
c
c   Add two real variables:  
c
c     ecosw is the phase difference between secondary and primary eclipse
c
c     tratio  is the ratio of T_2/T_1
c
c
c   August 10, 2004
c
c   This is version 3.  Some improvements:
c
c   o  There is the option to use a Monte Carlo technique to
c      account for fractional pixels.
c
c   o  The user can input and fit for the primary mass, the K-velocity
c      of the primary, the radaius of the primary, and the ratio
c      of the radii.
c
c
c
c
c   UPDATE JUNE 11, 2003
c
c   This version will have most of the two dimensional
c   arrays replaced with 1D arrays.  For example, the
c   code fragment
c
c          do 10 ialf=1,Nalf
c            do 9 ibet=1,ibetlim(ialf)
c               temp=surf(ialf,ibet)
c
c    becomes
c
c          do 10 ialf=1,Nalf
c            do 9 ibet=1,ibetlim(ialf)
c               index=(ialf-i)*ibetlim(ialf)+ibet
c               temp=surf(index)
c          
c
c
c   May 12, 2000
c
c   This is version 2.  Some improvements:
c
c   o  A polar coordinate system is now used.
c
c   o  The star horizons are now more accurately found.
c
c   o  Fractionaly eclipsed pixels are not accounted for.
c
c   o  A more robust reflection effect routine.
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   UPDATE August 13, 2001
c
c   Replace in the *format* statements, the string isw9 (curren....)
c   with 'ielite'
c
c
c   UPDATE NOVEMBER 27, 2006
c
c   Add NRVphase and xRVmod to the list.  This will allow for the
c   velocity curve to have points every dphase when iecheck=5.
c
          subroutine lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $      ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     &      ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &      fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
c   October 6, 1999
c
c   This subroutine will call all of the needed subroutines
c   to set up geometry, etc.
c
c   Adjust the maximum limits as needed.
c
c
c   March 13, 2000
c
c   change the integer flag iswe1 to ionephase.  If ionephase=1, then
c   compute only one phase, given by sw1.
c
c
c   March 20, 2000
c
c   I have speeded up the interpolation of the atmosphere intensities
c   from the table by about a factor of 3.
c
c
c   July 24, 2000
c
c   Add the array obsparm to the argument list.  This contains the 
c   computed physical parameters of stars 1 and 2
c
c   UPDATE September 21, 2008
c
c   Add K_1 and K_2 to the list
c
c   UPDATE October 10, 2008
c
c   Add finc, Q, ecc, arg, Teff1, Teff2 to the list 
c
c   obsparm(1) = mass of star 1      (solar)
c   obsparm(2) = radius of star 1    (solar)
c   obsparm(3) = gravity of star 1   (log cgs)
c   obsparm(4) = V_rot*sin(i) star 1 (km/sec)
c   obsparm(5) = mass of star 2      (solar)
c   obsparm(6) = radius of star 2    (solar)
c   obsparm(7) = gravity of star 2   (log cgs)
c   obsparm(8) = V_rot*sin(i) star 2 (km/sec)
c   obsparm(9) = duration of X-ray eclipse (degrees)
c   obsparm(10) = K_1 (km/sec)
c   obsparm(11) = K_2 (km/sec)
c   obsparm(12) = finc (deg)
c   obsparm(13) = Q
c   obsparm(14) = ecc
c   obsparm(15) = argper (deg)
c   obsparm(16) = Teff1
c   obsparm(17) = Teff2
c

c   UPDATE September 11, 2001
c
c   Add the X-ray eclipse duration in degrees (if any) to obsparm(9).
c
c   **************************************************
c
c   February 5, 2001
c
c   This major revision includes the generalization to eccentric orbits.
c
c
c   NEW BUG August 2, 2001
c
c   Do a global replace and put 'T0' in place of sw4.
c
c   UPDATE June 17, 2002
c
c   Replace all cases of eq.. with eqv..
c   This will make statements with logical orperands like this:
c
c   if(left.eqv..false.)   rather than   if(left.eq..false.)
c
          implicit double precision(a-h,o-z)
c
c   UPDATE MAY 21, 2004
C
C   Add separate ialphmax and ibetmax parameters for star 1 and star 2.
c   The argument lists for copyinty, writepoints, detailref, simpleref
c   initratio, and checkinput need to be modified:  add ialphmax2 and
c   ibetmax2 to the end.
c 
c
          parameter(ialphmax1=1000,ibetmax1=1000)
          parameter(ialphmax2=200,ibetmax2=800)
c
c   Caution:  adjust also the values of tempalf,tempbet in the 
c   subroutine getATMflux!
c
          parameter(Nthetamax=160,Nrmax=160)
c
          parameter (maxlines=1300,maxmu=115)   ! was 1100
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
          parameter(pie=3.14159265358979323d0)
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm to 9
c
c   UPDATE September 20, 2008
c
c   Make the dimension of obsparm 11
c
c   UPDATE October 10, 2008
c
c   Make the dimension of obsparm 17
c
          dimension corr1(8),corr2(8)
          dimension refflux1(8),refflux2(8),rlatflux(8),gimvel(8)
          dimension obsparm(17)
          dimension darkint1(8),darkint2(8)
          dimension mmdx1(ialphmax1,ibetmax1),mmdx2(ialphmax2,ibetmax2)
          dimension dtemp(Nrmax*Nthetamax),dx(Nrmax*Nthetamax),
     &      dy(Nrmax*Nthetamax),dz(Nrmax*Nthetamax),dinty(Nrmax*Nthetamax),
     $      drad(Nrmax),einty(Nthetamax*11),savedinty(Nrmax*Nthetamax),
     &      saveeinty(Nthetamax*11),ibetlim1(ialphmax1),ibetlim2(ialphmax2),
     %      tedge(Nthetamax*11),xedge(Nthetamax*11),yedge(Nthetamax*11),
     &      zedge(Nthetamax*11),dxhoriz(2*Nthetamax),dyhoriz(2*Nthetamax),
     $      dtopx(2*Nthetamax),dtopy(2*Nthetamax),dratio(Nrmax*Nthetamax)
          dimension diskproj(Nrmax*Nthetamax),edgeproj(Nthetamax*11),
     &      dvisib(Nrmax*Nthetamax),evisib(Nthetamax*11),
     &      xskydisk(Nthetamax*Nrmax),yskydisk(Nthetamax*Nrmax),
     &      zskydisk(Nthetamax*Nrmax),
     &      xskyedge(Nthetamax*11),yskyedge(Nthetamax*11)
          dimension x1(ialphmax1*ibetmax1),y1(ialphmax1*ibetmax1),
     $      z1(ialphmax1*ibetmax1),x2(ialphmax2*ibetmax2),
     &      y2(ialphmax2*ibetmax2),
     $      z2(ialphmax2*ibetmax2),visib1(ialphmax1*ibetmax1),
     %      visib2(ialphmax2*ibetmax2),rinty1(ialphmax1*ibetmax1),
     #      rinty2(ialphmax2*ibetmax2),flum1(ialphmax1*ibetmax1),
     %      flum2(ialphmax2*ibetmax2),saveinty1(ialphmax1*ibetmax1),
     $      saveinty2(ialphmax2*ibetmax2),phiar1(ialphmax1*ibetmax1),
     &      phiar2(ialphmax2*ibetmax2),delphi1(ialphmax1*ibetmax1),
     %      delphi2(ialphmax2*ibetmax2),iedgestar1(ialphmax1*ibetmax1),
     #      iedgestar2(ialphmax2*ibetmax2),iedgehor2(ialphmax2*ibetmax2),
     %      iedgehor1(ialphmax1*ibetmax1),delphie1(ialphmax1*ibetmax1),  
     $      delphie2(ialphmax2*ibetmax2)
          dimension phihor1(ialphmax1,2),phihor2(ialphmax2,2)
          dimension surf1(ialphmax1*ibetmax1),surf2(ialphmax2*ibetmax2),
     $      gradx1(ialphmax1*ibetmax1),gradx2(ialphmax2*ibetmax2),
     $      grady1(ialphmax1*ibetmax1),grady2(ialphmax2*ibetmax2),
     $      gradz1(ialphmax1*ibetmax1),gradz2(ialphmax2*ibetmax2)
          dimension rad1(ialphmax1*ibetmax1),rad2(ialphmax2*ibetmax2),
     $      temp1(ialphmax1*ibetmax1),temp2(ialphmax2*ibetmax2),
     $      g1(ialphmax1*ibetmax1),g2(ialphmax2*ibetmax2),Rpol1(ialphmax1),
     #      Rpol2(ialphmax2)
          dimension xhoriz1(4*ibetmax1),yhoriz1(4*ibetmax1),
     $       xhoriz2(4*ibetmax2),yhoriz2(4*ibetmax2),
     %       xsky1(ialphmax1*ibetmax1*4),
     $       ysky1(ialphmax1*ibetmax1*4),xsky2(ialphmax2*ibetmax2*4),
     %       ysky2(ialphmax2*ibetmax2*4),xend1(4),xend2(4),
     $       ratio1(ialphmax1*ibetmax1),ratio2(ialphmax2*ibetmax2),
     $       tempold1(ialphmax1*ibetmax1),tempold2(ialphmax2*ibetmax2),
     $       projarray1(ialphmax1*ibetmax1),projarray2(ialphmax2*ibetmax2),
     %       dumxsky1(ialphmax1*ibetmax1*4),
     $       dumysky1(ialphmax1*ibetmax1*4),dumxsky2(ialphmax2*ibetmax2*4),
     %       dumysky2(ialphmax2*ibetmax2*4),coprat1(ialphmax1*ibetmax1),
     #       coprat2(ialphmax2*ibetmax2)
          dimension xtop2horiz(4*ibetmax2),ytop2horiz(4*ibetmax2)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
          dimension xmod(Nmaxphase),ymodU(Nmaxphase),ymodB(Nmaxphase),
     $      ymodV(Nmaxphase),ymodR(Nmaxphase),ymodI(Nmaxphase),
     $      ymodJ(Nmaxphase),ymodH(Nmaxphase),ymodK(Nmaxphase),
     &      ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      RV1(Nmaxphase),RV2(Nmaxphase),zdcorr(8),dRV1(Nmaxphase),
     %      dRV2(Nmaxphase),third(8),ymods3(Nmaxphase)
          dimension xRVmod(Nmaxphase),yeclipse(720001)  !Nmaxphase
          dimension phistart1(ialphmax1),phistart2(ialphmax2)
c
c
          dimension powercoeff(8,9)
c
c   UPDATE September 10, 2001
c
c   Add these arrays for X-ray eclipse duration computations
c
          dimension xecx(5000),xecy(5000)
c
c   RVG BUG ALERT   May 8, 2001
c
c   Define the variables needed for spots
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
c
c   UPDATE JULY 7, 2004
c
c   Use the sub-random Sobel sequence instead of ran9.  xsob is needed
c   for this.
c
          dimension xsob(2)
c
c   UPDATE DECEMBER 10, 2004
c
c   Add these arrays to compute light curves in time units rather
c   than phase
c 
          dimension timearray(900000)
c
c   UPDATE OCTOBER 10, 2007
c
c   Add this array to keep track of the light curves of the individual
c   components in all band passes.
c
c   UPDATE January 12, 2009
c
c   make fracs fracs1, fracs2, .... fracs8 and put them in the
c   argument of subroutine lightcurve.  In this was one can use
c   Nmaxphase in the dimension statement
c
          dimension compfracs(8,2)   
          dimension fracs1(Nmaxphase,3),fracs2(Nmaxphase,3)
          dimension fracs3(Nmaxphase,3),fracs4(Nmaxphase,3)
          dimension fracs5(Nmaxphase,3),fracs6(Nmaxphase,3)
          dimension fracs7(Nmaxphase,3),fracs8(Nmaxphase,3)

          character*9 extension
c
c   NEW BUG ALERT  July 13, 2001
c
c   Add a new character string and common block for a 'parameter string'
c   This string of parameters will be fed to the genetic code to make it
c   easier to compute uncertainties on the physical quantities like mass
c   and radius.
c
c   UPDATE November 28, 2001
c
c   parmstring was character*199, now should be character*201
c
c   UPDATE January 16, 2002
c
c   parmstring was character*201, now should be character*227
c
c   UPDATE June 7, 2002
c
c   Make the length of parmstring character*237.
c   Also, make a two scratch strings, one of length character*227
c   and one of length 10
c
c   UPDATE October 28, 2002
c
c   Make the length of parmstring character*249, and add 12 to the
c   length of scr1string
c
c   UPDATE October 22, 2008
c
c   Make the length of parmstring character*259, and add 12 to the
c   length of scr1string
c
          character*259 parmstring
c          character*227 scr1string
          character*249 scr1string
          character*10  scr2string
c
          common /stringblock/ parmstring
c
c   UPDATE August 10, 2004
c
c   Add the 8 variables below.
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff
c
c

          common /realblock/ fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     &       alb1,alb2,
     %       rLx,Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,dwavey,
     %       t3,g3,SA3,density,sw1,sw2,sw3,T0,
     %       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     $       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       bigI,bigbeta,sw23,sw24,powercoeff,
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2
c
c   RVG BUG ALERT   May 8, 2001
c
c   Add this common block to all programs
c
          common /spotblock/ spot1parm,spot2parm,spotdparm
c
c   UPDATE August 10, 2004
c
c   Add the 4 variables below.
c
c   UPDATE November 6, 2008
c
c   Add isw25-isw34 below
c
          common /intblock/ Nalph1,Nbet1,Nalph2,Nbet2,
     &       Ntheta,Nradius,Nref,
     $       idraw,iecheck,idint,iatm,ism1,
     $       ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     &       iRVfilt,isw1,isw2,isw3,isw4,
     %       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,idark1,idark2,isw12,isw13,
     %       isw21,isw22,isw23,isw24,
     %       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24 to the block above.
c
c
c   RVG BUG ALERT  June 12, 2001
c
c   Add these common blocks for the model atmosphere variables
c
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,atmint4,
     %       atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin
          common /intatm/  Nlines,Nmu
c
c   October 10, 2007
c
c   Add this common block for the luminosity ratios and disk fractions
c
          common /fracblock/ compfracs
c
c
c   Start here.  First check the input for goofs.
c
c   RVG BUG ALERT:   May 16, 2001
c
c   The recordparm subroutine has been moved to here.  
c
c   UPDATE August 10, 2004
c
c   Add the new variables to the argument list.
c
c
          itide=isw26
          tidephi=sw31
c
c    Update October 15, 2010
c
c    Add Doppler boosting.  beam1=sw33 and beam2=sw34
c

          call recordparm(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     &       Ntheta,Nradius,alb1,alb2,Nref,
     %       rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
     $       idraw,iecheck,idint,iatm,ism1,
     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,
     &       isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey,
     $       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     $       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
     %       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       idark1,idark2,isw12,isw13,
     #       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,
     #       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2,
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
c    July 29, 2005
c
c    If ecosw > 0 and ecc > 0 then set the value of argper
c
          call checkinput(ialphmax1,ibetmax1,Nthetamax,Nrmax,
     &       Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,Ntheta,Nradius,darkbol1,darkbol2,alb1,alb2,
     $       Nref,rLx,Period,fm,separ,
     $       idraw,iecheck,idint,iatm,ism1,
     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     &       ivelout,iXout,ecc,pshift,sw5,sw9,ialphmax2,ibetmax2,ecosw,argper,
     &       temprat)
c
c  UPDATE October 10, 2008
c 
c  Add the routine setdensity which sets the mass ratio given the density
c  and the fill factor of star 1
c

c
c   Record parameters that could be optimized into the ELC.parm file
c
          call recordparm1(fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,
     %       rLx,separ,gamma,t3,g3,SA3,density)
c
c  UPDATE DECEMBER 7, 2009
c
c  Add a very fast analytic mode for transiting planets.  set isw27=1.
c  also needed are frac1,frac2,period,T0,primmass,primK,ecc,argper,
c  omega1,finc,ilaw,dwavex,dwavey,axis_I,axis_beta,bin size
c  for light curves, bin size for velocity curve,ialign,ikeep
c
          if(isw27.ge.1)then
c
c   Determine the masses, radii, and rotational velocities
c
            call analyticscale(Q,finc,period,primmass,primK,ecc,frac1,frac2,
     &        primrad,ratrad,reff1,reff2,separ,vrot1,vrot2,gp1,gp2,omega1,
     $        omega2)
c
c   Determine the reference fluxes.  These will be stored in darkint(k)
c
            call getanalyticint(maxlines,maxmu,Nlines,
     &        atmT,atmg,atmmu,Nmu,
     &        atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &        Tmax,Tmin,gmax,gmin,gp1,darkint1,teff1,
     %        dwavex,dwavey,ilaw,iatm,1,wave,reff1)
c
            call getanalyticint(maxlines,maxmu,Nlines,
     &        atmT,atmg,atmmu,Nmu,
     &        atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &        Tmax,Tmin,gmax,gmin,gp2,darkint2,teff2,
     %        dwavex,dwavey,ilaw,iatm,2,wave,reff2)
c
            call fastanalytic(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $       ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     &       ymodd,RV1,RV2,drv1,drv2,obsparm,NRVphase,xRVmod,
     &       fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,
     $       frac1,frac2,primrad,ratrad,period,T0,primmass,
     #       primK,ecc,argper,omega1,omega2,
     $       finc,ilaw,dwavex,dwavey,bigI,bigbeta,sw29,sw30,isw21,ikeep,
     %       pshift,reff1,reff2,darkint1,darkint2,idark1,
     %       idark2,dphase,iRVfilt,vrot1,vrot2,Q,separ,gamma,isw27,
     &       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     &       timearray,sw9,sw23,sw24,isw7,SA3)
c
            rpole1=1.0d0
            rpole2=1.0d0
            bdist=1.0d0
c
            call parms(1,
     &        Teff2,Q,finc,separ,period,reff1,reff2,
     %        rpole1,rpole2,fill1,fill2,
     $        gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,bdist,ecc)
c
            call parms1(Teff2,Q,finc,separ,period,reff1,reff2,
     %        rpole1,rpole2,fill1,fill2,
     $        gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,
     $        obsparm,bdist,ecc,argper,teff1)
c
            pot1=0.0d0
            pot2=0.0d0
            ave11=0.0d0
            ave12=0.0d0
            ave21=0.0d0
            ave22=0.0d0
            ave1=0.0d0
            ave2=0.0d0
c
            call lineparms1(Teff2,Q,
     %        finc,separ,period,reff1,reff2,
     $        vrot1,vrot2,omega1,omega2,
     $        bdist,ecc,SA3,ave11,ave12,ave21,ave22,ave1,ave2,
     #        parmstring,pot1,pot2)
c
            close(2)
            return
          endif

c   UPDATE JULY 21, 2006
c
c   Add a "fast genetic" mode.  If ifastflag=1, then set 
c     Nalph1=40 
c     Nbet1=14
c     Nalph2=40
c     Nbet2=14
c     dphase=3
c 
c   We need to save the values to reset at the end of the subroutine
c
         iNa1save=Nalph1
         iNa2save=Nalph2
         iNb1save=Nbet1
         iNb2save=Nbet2
         savedphase=dphase
c

         if(ifastflag.ge.1)then
           Nalph1=40
           Nalph2=40
           Nbet1=14
           Nbet2=14
           dphase=10.0d0*savedphase
           if(dphase.gt.3.0d0)dphase=3.0d0
         endif

c
c   Save the value of the flag iecheck so we can reset it at the end of the
c   routine.
c
          iesave=iecheck
c
c  UPDATE May 8, 2006
c
c  Add new flags and variables, including ialign, bigI, bigbeta.
c
c          bigI=sw21
c          bigbeta=sw22
          ialign=isw21
          if(ialign.le.0)then
            bigI=finc
            bigbeta=0.0d0
          endif
c
          ionephase=isw1
          onephase=sw1
          isquare=isw2
          iusepot=isw3
          usepot1=sw2
          usepot2=sw3
          avesep=separ
          bdist=1.0d0
c
          fluxU2=0.0d0
          fluxB2=0.0d0
          fluxV2=0.0d0
          fluxR2=0.0d0
          fluxI2=0.0d0
          fluxJ2=0.0d0
          fluxH2=0.0d0
          fluxK2=0.0d0
c
c   UPDATE May 10, 2006
c
c   If we are in analytic mode (isw12>0), then set iehceck=0 and ism1=0
c
c
          if(isw12.gt.0)then
            ism1=0
            iecheck=1
          endif
c
c   RVG BUG ALERT   May 8, 2001
c
c   Define the 'simpson switch' and the 'gravity exponent switch'
c
          tteff2=0.0d0
          if(sw5.gt.0.0d0)tteff2=teff2
c
c   UPDATE October 13, 2008
c
c   Redefine the isw5 switch.  It will not control how the spot
c   temperature profile is computed.
c
c   ispotprof=0    constant temperature factor
c   ispotprof=1    linear change in temperature profile
c   ispotprof=2    Gaussian change
c
          isimp=0
          ispotprof=isw5
          igrav=isw6
c
c   Count the number of spots on each star.
c
          ispot1=0
          ispot2=0
          ispotd=0
          do 4500 ii=1,2
            if(spot1parm(ii,1).gt.0.0d0)ispot1=ispot1+1          
            if(spot2parm(ii,1).gt.0.0d0)ispot2=ispot2+1          
            if(spotdparm(ii,1).gt.0.0d0)ispotd=ispotd+1          
 4500     continue
c
          if(igrav.eq.1)call gravexp(Teff1,Tgrav1,1)
          if((igrav.eq.2).and.(teff2.gt.0.0d0))call gravexp(Teff2,Tgrav2,2)
          if(igrav.eq.3)then 
            call gravexp(Teff1,Tgrav1,1)           
            if(Teff2.gt.0.0d0)call gravexp(Teff2,Tgrav2,2)
          endif
c
c   RVG BUG ALERT  June 12, 2001
c
c   Load the atmosphere table in the main programs (ELC.for, gridELC.for, ...)
c
cc
cc  If the flag iatm>0, then load the model atmosphere table.
cc
c          if(iatm.ge.1)then
c            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
c     &         atmint,Tmax,Tmin,gmax,gmin)
c             write(2,102)Nlines,Tmax,Tmin
c          endif
c
c   Get the radii, coordinates, gradients, etc.  These do not depend on phase.
c
          ivrt=0
          fillper1=fill1
          fillper2=fill2
          bdist=1.0d0-ecc
          pervol1=1.0d0
          pervol2=1.0d0
c
c   RVG BUG ALERT   April 19, 2001
c
c   If isynch = 1, then set the omega values so that the rotation is
c   synchronous at periastron.
c
          if((ecc.gt.0.0d0).and.(isynch.ge.1))then 
c            sss=dsqrt((1.0d0+ecc)/(1.0d0-ecc**3))
c
c   UPDATE January 13, 2009
c
c   fix the error
c
            sss=dsqrt((1.0d0+ecc)/(1.0d0-ecc)**3)
            omega1=sss
            omega2=sss
          endif
c
c   END BUG
c
c
c   UPDATE JULY 2, 2004
c
c   Initialize the variable jdum here.  It goes in the argument list of
c   get visib.
c
          jdum=-123456
c
c          tstart=1000.0d0
c          tstop=1020.0d0
c          tstep=0.01d0
c
          tstep=sw9
          tstart=sw23
          tstop=sw24

          if(isw7.eq.2)call filltime(ntime,timearray,tstart,tstop,tstep)

cc   UPDATE JULY 4, 2004
c
c   Assign the variable MonteCarlo to isw8.  It will be used
c   in the Monte Carlo routine to compute fractionally eclipsed
c   pixels.
c
          MonteCarlo=isw8
c
c   Initialize the Sobel sequence here
c
          nnn=2          
          call sobseq(nnn,xsob)
          call sobseq(nnn,xsob)
          call sobseq(nnn,xsob)
          nnn=-1
          call sobseq(nnn,xsob)
c
c   UPDATE December 21, 2008
c
c   Here is a new subroutine call.  radfill sets the filling factor based on the 
c   effective radius, rather than the L_1 point distance
c
          radfill1=sw27
          radfill2=sw28
          bdist=1.0d0-ecc
          call setfill(1,Q,fill1,radfill1,omega1,bdist,tidephi,itide)
          call setfill(2,Q,fill2,radfill2,omega2,bdist,tidephi,itide)
c
c   UPDATE September 11, 2001
c
c   Add the iverb flag to the end of setupgeo.
c

          omegatwo=omega2
          iverb=0
          call setupgeo(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $       fill1,omega1,Q,finc,
     $       x1,y1,z1,
     %       surf1,rad1,gradx1,grady1,gradz1,g1,xend1,separ,
     %       Tgrav1,Teff1,reff1,Rl1,Tpole1,Rpol1,Regg1,SA1,pot1,gpole1,
     %       phiar1,isquare,iusepot,usepot1,ivrt,pervol1,fillper1,bdist,
     #       pots1,iverb,mmdx1,primmass,primK,primrad,ratrad,frac1,frac2,
     &       ecc,period,size1,fill2,omegatwo,sw5,tteff2,density,
     $       tidephi,itide,phistart1)
c
          if(teff2.gt.0.0d0)then
            call setupgeo(2,ialphmax2,ibetmax2,Nalph2,Nbet2,ibetlim2,
     $        fill2,omega2,Q,finc,
     $        x2,y2,z2,
     %        surf2,rad2,gradx2,
     &        grady2,gradz2,g2,xend2,separ,
     %        Tgrav2,Teff2,reff2,Rl2,Tpole2,Rpol2,Regg2,SA2,pot2,gpole2,
     $        phiar2,isquare,iusepot,usepot2,ivrt,pervol2,fillper2,bdist,
     $        pots2,iverb,mmdx2,primmass,primK,primrad,ratrad,frac1,frac2,
     &        ecc,period,size1,fill2,omegatwo,sw5,tteff2,density,
     $        tidephi,itide,phistart2)
          else
c
c   UPDATE March 22, 2002
c
c   Add ibetlim2 to the argument list of dummyvalues
c
c
            call dummyvalues(ialphmax2,ibetmax2,Nalph2,Nbet2,
     $        x2,y2,z2,surf2,rad2,gradx2,grady2,gradz2,g2,xend2,darkbol2,
     $        temp2,ibetlim2,mmdx2)
            Regg2=rocheradius(Q)
            reff2=Regg2
            overQ=1.0d0/Q
            call findL1(overQ,omega2,x0,1,bdist,tidephi,itide)
            Rl2=x0
          endif
c
          if(idint.ge.1)call disksetup(Nthetamax,Nrmax,Ntheta,Nradius,
     %       betarim,rinner,router,Regg2,Rl2,separ,
     $       tdisk,xi,dtemp,dx,dy,dz,drad,
     $       tedge,xedge,yedge,zedge,redge,stepr,stepz,bdist,
     #       ivrt,reper,rsper)
c
c
c  If the model atmosphere option is on, then compute the third light.
c
c  UPDATE April 3, 2002
c
c  Add separ to the argument list of thirdlight.
c
          call thirdlight(iatm,t3,g3,SA3,SA1,third,
     &        maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &        atmint1,atmint2,atmint3,atmint4,atmint5,
     &        atmint6,atmint7,atmint8,
     &        Tmax,Tmin,gmax,gmin,
     %        icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,separ,
     &        dwavex,dwavey,ilaw,istar)
c
c   If the orbit is eccentric, then we have to loop over phases here.
c   Otherwise, we loop a bit further down.  
c
c   RVG BUG ALERT   April 19, 2001
c
c   Use the symmetry of the orbit in the mean anomaly M to only compute
c   half the phases for an eccentric orbit.  Thus, do not disable the ism1
c   switch.   Also, define the variable pstep as below (equals dphase by
c   default)
c
c          if(ecc.gt.0.0d0)ism1=0

          dpsave=dphase
          pstart=0.0d0
          pstop=360.0d0-dphase
          pconj=pie
          pconj2=0.0d0
          pstep=dphase
          idcheck=100
          if(ism1.ge.1)pstop=180.0d0

c
c   UPDATE October 18, 2002
c
c   Use the input flags sw7 and sw8 to define a phase range
c   to compute.  Require sw7 > 0 and sw8 > 0  AND  sw7 < sw 8
c
          if((sw7.gt.0.0d0).and.(sw8.gt.0.0d0).and.(sw7.lt.sw8))then
            if((sw8-sw7).lt.dphase)dphase=sw8-sw7
            pstart=sw7
            pstop=sw8-dphase
          endif
          pstartout=onephase
          pstopout=onephase
          if(ecc.gt.0.0d0)then
            pstartout=0.0d0
            pstopout=360.0d0-dphase
            if(ism1.ge.1)pstopout=180.0d0
          endif
c
c   RVG BUG ALERT  April 19, 2001.
c
c   If the ism1 flag is set, then modify the value of pstopout as above.
c
          open(unit=64,file='ELC.phases',status='unknown')
c
c   RVG BUG ALERT  May 3, 2001
c
c   Open a new output file for eccentric orbits.
c
          if(ecc.gt.0.0d0)open(unit=65,file='ELC.eccentric',status='unknown')
c
          itoggle=-1
          icount=0
          togglephase=360.0d0
c
          argrad=argper*pie/180.0d0
c
c  Here is some code adapted from Wilson-Devinney to keep track of
c  phases needed for eccentric orbits.
c
          trc=0.5d0*pie-argrad    
 1139     if(trc.lt.0.d0) trc=trc+2.0d0*pie
          if(trc.lt.0.d0) goto 1139                                           
 1140     if(trc.ge.2.0d0*pie) trc=trc-2.0d0*pie
          if(trc.ge.2.0d0*pie) goto 1140                            
          htrc=0.5d0*trc                                                     
          if(dabs(0.5*pie-htrc).lt.7.d-6) goto 11101              
          if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 11101              
          ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))           
          goto 11103                                                          
11101     ecan=pie
11103     xmc=ecan-ecc*dsin(ecan)                                             
          if(xmc.lt.0.d0) xmc=xmc+2.0d0*pie
          phper=1.d0-xmc/(2.0d0*pie)
          pconj=(xmc+argrad)/(2.0d0*pie)-0.25d0
c
c   UPDATE March 14, 2008
c
c   Make sure the conjunction phase is between 0 and 1
c
          if(pconj.gt.1.0d0)pconj=pconj-1.0d0     
c
c   UPDATE September 10, 2001
c
c   Make this new block to compute the conjunction phase for star 2.
c
          trc=0.5d0*pie-argrad+pie    
 3139     if(trc.lt.0.d0) trc=trc+2.0d0*pie
          if(trc.lt.0.d0) goto 3139                                           
 3140     if(trc.ge.2.0d0*pie) trc=trc-2.0d0*pie
          if(trc.ge.2.0d0*pie) goto 3140                            
          htrc=0.5d0*trc                                                     
          if(dabs(0.5*pie-htrc).lt.7.d-6) goto 31101              
          if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 31101              
          ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))           
          goto 31103                                                          
31101     ecan=pie
31103     xmc=ecan-ecc*dsin(ecan)                                             
          if(xmc.lt.0.d0) xmc=xmc+2.0d0*pie
          phper2=1.d0-xmc/(2.0d0*pie)
          pconj2=(xmc+argrad)/(2.0d0*pie)-0.25d0
c
c   UPDATE March 14, 2008
c
c   Make sure the conjunction phase is between 0 and 1
c
          if(pconj2.gt.1.0d0)pconj2=pconj2-1.0d0     
c
c  UPDATE September 10, 2001
c
c  Add the if-then clauses to account for ikeep=1 for star 1 cna
c  ikeep = 2 for star 2.
c
c          eshift=pconj+phper-0.5
          if(ikeep.eq.1)eshift=phper+pconj-0.5
          if(ikeep.eq.2)eshift=phper2+pconj2
c
          if(ecc.eq.0.0d0)then
            eshift=0.0d0
            pconj=0.0d0
          endif
c    
c  UPDATE May 3, 2006
c
c  Add a "fast transit" mode.  If isw13 > 1, then find the range
c  of latitude rows on the star that are eclipsed.  These values
c  are returned as ialfmin,ialfmax
c
          ialfmin=999999
          ialfmax=-111111
          if(isw13.ge.1)then
            FINCR = (FINC/180.0d0)*pie    !orbital inclination in radians
            pup=360.0d0
            if(ism1.ge.1)pup=180.0d0+dphase
            sf2=dsin(fincr)**2
            thetatan=180.0d0*dacos(dsqrt((1.0d0-(reff1)**2)/sf2))/pie

            do 6666 pp=pup-thetatan+dphase,(0.0d0+dphase),-dphase
              PHASER = (pp/180.0d0)*pie     !orbital phase in radians
c
              delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
              delta=bdist*dsqrt(delta)
              if(delta.gt.reff1+reff2)go to 6666
c
              tt1=dabs(phaser-2.0d0*pie*pconj)
              tt2=dabs(phaser-2.0d0*pie*(pconj+1.0d0))
              if(tt1.le.tt2)then
                diff1=tt1
              else
                diff1=tt2
              endif
              tt1=dabs(phaser-2.0d0*pie*pconj2)
              tt2=dabs(phaser-2.0d0*pie*(pconj2+1.0d0))
              if(tt1.le.tt2)then
                diff2=tt1
              else
                diff2=tt2
              endif              
              if(diff1.lt.diff2)go to 6666

              dummyphase=dmod(pp+180.0d0,360.0d0)
              call gethorizon(2,ialphmax2,ibetmax2,Nalph2,
     %          Nbet2,ibetlim2,dummyphase,finc,Q,
     $          pot2,omega2,x2,y2,z2,rad2,gradx2,grady2,gradz2,xend2,
     $          Nhoriz2,xhoriz2,yhoriz2,phiar2,iedgestar2,delphie2,bdist,
     &          mmdx2,xhmin2,xhmax2,yhmin2,yhmax2,tidephi,itide,phihor2)
              call getalflim(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %          pp,finc,Q,psi0,omega1,
     $          gradx1,grady1,gradz1,x1,y1,z1,
     $          Nhoriz2,xhoriz2,yhoriz2,
     &          separ,bdist,mmdx1,ialfmin,ialfmax,reff1)

6666          continue
6667          format(f9.4,2x,2(i4,1x),2(f10.7,2x))
c
              pp=180.0d0
              dummyphase=dmod(pp+180.0d0,360.0d0)
              call gethorizon(2,ialphmax2,ibetmax2,Nalph2,
     %          Nbet2,ibetlim2,dummyphase,finc,Q,
     $          pot2,omega2,x2,y2,z2,rad2,gradx2,grady2,gradz2,xend2,
     $          Nhoriz2,xhoriz2,yhoriz2,phiar2,iedgestar2,delphie2,bdist,
     &          mmdx2,xhmin,xhmax,yhmin,yhmax,tidephi,itide,phihor2)
              call getalflim(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %          pp,finc,Q,psi0,omega1,
     $          gradx1,grady1,gradz1,x1,y1,z1,
     $          Nhoriz2,xhoriz2,yhoriz2,
     &          separ,bdist,mmdx1,ialfmin,ialfmax,reff1)

              ialfmin=ialfmin-1
              ialfmax=ialfmax+1
c
c   Now we have to find the integrated brightness of the pixels on star
c   1 outside the range ialfmin,ialfmax
c
            if(idark1.le.0)call getrefvisib(1,ialphmax1,ibetmax1,
     &        Nalph1,Nbet1,ibetlim1,
     %        0.0d0,finc,Q,psi0,omega1,
     $        gradx1,grady1,gradz1,visib1,projarray1,separ,bdist,mmdx1)
c
            if(idark1.le.0)call setuptemp(1,
     &         ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     &         g1,Tpole1,Tgrav1,temp1,gpole1,mmdx1)
c
            if(iatm.eq.0)then
              do 1937 jj=1,8
c
                flux1=0.0d0
                if(idark1.le.0)then
                  rpole1=Rpol1(Nalph1/2)
c
                  www=wave(jj)
                  flimbx=dwavex(jj,1)
                  flimby=dwavey(jj,1)
                  flux1=0.0d0
                  call getBBlumcor(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $            www,visib1,projarray1,temp1,surf1,flimbx,flimby,
     &            ilaw,rinty1,flum1,flux1,rldint1,separ,mmdx1,ialfmin,
     $            ialfmax)
                endif
c
                rlatflux(jj)=flux1
 1937          continue
            endif
          endif  !if fast transit mode
c
c
c   UPDATE June 16, 2003
c
c   Add an EBOP mode.  If iecheck=9, then compute 'reference' fluxes
c   for star 1 and star 2.  Then at each phase, check the distance
c   between the two stellar centers.  If the distance is more than
c   the sum of the radii then skip the flux computation.
c
c
c   UPDATE April 24, 2006
c
c   If isw12=1, then compute analytic transits using Mandel and Agol.
c   Compute reference fluxes as in EBOP mode, then compute the ratio
c   of the transit given the separation between centers and the radius
c   ratio
c
c

          if((iecheck.eq.9).or.(isw12.ge.1))then
            if(idark1.le.0)call getrefvisib(1,ialphmax1,ibetmax1,
     &        Nalph1,Nbet1,ibetlim1,
     %        0.0d0,finc,Q,psi0,omega1,
     $        gradx1,grady1,gradz1,visib1,projarray1,separ,bdist,mmdx1)
c
            if(idark2.le.0)call getrefvisib(2,ialphmax2,ibetmax2,
     &        Nalph2,Nbet2,ibetlim2,
     %        180.0d0,finc,Q,psi0,omega2,
     $        gradx2,grady2,gradz2,visib2,projarray2,separ,bdist,mmdx2)
c
            if(idark1.le.0)call setuptemp(1,
     &         ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     &         g1,Tpole1,Tgrav1,temp1,gpole1,mmdx1)
            if((teff2.gt.0.0d0).and.(idark2.le.0))call setuptemp(2,
     #            ialphmax2,ibetmax2,
     %            Nalph2,Nbet2,ibetlim2,g2,Tpole2,Tgrav2,temp2,gpole2,
     $            mmdx2)
c
            if(iatm.eq.0)then
              do 937 jj=1,8
c
c   UPDATE MAY 28, 2010
c
c   Add corr1 and corr2 for analytic mode
c
                corr1(jj)=0.0d0
                corr2(jj)=0.0d0
                flux1=0.0d0
                if(idark1.le.0)then
                  rpole1=Rpol1(Nalph1/2)
c
                  if(jj.eq.1)then
                    call parms(1,
     &                Teff2,Q,finc,separ,period,reff1,reff2,
     %                rpole1,rpole2,fill1,fill2,
     $                gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,
     &                omega2,bdist,ecc)
c
                    call parms1(Teff2,Q,finc,separ,period,reff1,reff2,
     %               rpole1,rpole2,fill1,fill2,
     $               gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,
     $               obsparm,bdist,ecc,argper,teff1)
c
                  endif
                  www=wave(jj)
                  flimbx=dwavex(jj,1)
                  flimby=dwavey(jj,1)
                  flux1=0.0d0
                  call getrefBBflux(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $            www,visib1,projarray1,temp1,surf1,flimbx,flimby,
     &            ilaw,rinty1,flum1,flux1,rldint1,separ,mmdx1)
                endif
c
                flux2=0.0d0
                if(idark2.le.0)then
c
                  rpole2=Rpol2(Nalph2/2)
c
                  if(jj.eq.1)then
                    call parms(1,
     &                Teff2,Q,finc,separ,period,reff1,reff2,
     %                rpole1,rpole2,fill1,fill2,
     $                gp1,gp2,vrot1,vrot2,gscale1,gscale2,
     %                omega1,omega2,bdist,ecc)
c
                    call parms1(Teff2,Q,finc,separ,period,reff1,reff2,
     %                rpole1,rpole2,fill1,fill2,
     $                gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,
     $                obsparm,bdist,ecc,argper,teff1)
c 
                  endif
                  www=wave(jj)
                  flimbx=dwavex(jj,2)
                  flimby=dwavey(jj,2)
                  flux2=0.0d0
                  if(teff2.gt.0.0d0)call getrefBBflux(ialphmax2,
     &              ibetmax2,Nalph2,Nbet2,ibetlim2,
     $              www,visib2,projarray2,temp2,surf2,flimbx,flimby,
     &              ilaw,rinty2,flum2,flux2,rldint2,separ,mmdx2)
                endif
                refflux1(jj)=flux1
                refflux2(jj)=flux2
 937          continue
            endif
            if(iatm.ge.1)then
              fluxU1=0.0d0
              fluxU2=0.0d0
              fluxB1=0.0d0
              fluxB2=0.0d0
              fluxV1=0.0d0
              fluxV2=0.0d0
              fluxR1=0.0d0
              fluxR2=0.0d0
              fluxI1=0.0d0
              fluxI2=0.0d0
              fluxJ1=0.0d0
              fluxJ2=0.0d0
              fluxH1=0.0d0
              fluxH2=0.0d0
              fluxK1=0.0d0
              fluxK2=0.0d0
c
              if(idark1.le.0)then
                rpole1=Rpol1(Nalph1/2)
c
                call parms(1,
     &             Teff2,Q,finc,separ,period,reff1,reff2,
     %             rpole1,rpole2,fill1,fill2,
     $             gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,bdist,ecc)
c
                call parms1(Teff2,Q,finc,separ,period,reff1,reff2,
     %             rpole1,rpole2,fill1,fill2,
     $             gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,
     $             obsparm,bdist,ecc,argper,teff1)
c
                call getATMint(maxlines,maxmu,Nlines,
     &            atmT,atmg,atmmu,Nmu,
     &            atmint1,atmint2,atmint3,atmint4,atmint5,
     &            atmint6,atmint7,atmint8,
     &            Tmax,Tmin,gmax,gmin,gscale1,darkint1,tpole1,gpole1,
     %            dwavex,dwavey,ilaw,iatm,1)
c
                call getrefATMflux(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $             visib1,projarray1,temp1,surf1,g1,rinty1,
     &             flum1,maxlines,maxmu,Nlines,
     &             atmT,atmg,atmmu,Nmu,
     &             atmint1,atmint2,atmint3,atmint4,atmint5,
     &             atmint6,atmint7,atmint8,
     &             Tmax,Tmin,gmax,gmin,gscale1,
     &             fluxU1,fluxB1,fluxV1,fluxR1,fluxI1,fluxJ1,fluxH1,fluxK1,
     &             icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,delphi1,
     %             delphie1,iedgestar1,iedgehor1,darkint1,separ,mmdx1,
     &             dwavex,dwavey,ilaw,iatm,1)
c
c   UPDATE MAY 28, 2010
c
c   Add corr1 and corr2 for analytic mode
c
                corr1(1)=0.0d0
                corr2(1)=0.0d0
                corr1(2)=0.0d0
                corr2(2)=0.0d0
                corr1(3)=0.0d0
                corr2(3)=0.0d0
                corr1(4)=0.0d0
                corr2(4)=0.0d0
                corr1(5)=0.0d0
                corr2(5)=0.0d0
                corr1(6)=0.0d0
                corr2(6)=0.0d0
                corr1(7)=0.0d0
                corr2(7)=0.0d0
                corr1(8)=0.0d0
                corr2(8)=0.0d0

              endif
              if(idark2.le.0)then
                rpole2=Rpol2(Nalph2/2)
c
                call parms(1,
     &             Teff2,Q,finc,separ,period,reff1,reff2,
     %             rpole1,rpole2,fill1,fill2,
     $             gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,bdist,ecc)
c
                call parms1(Teff2,Q,finc,separ,period,reff1,reff2,
     %             rpole1,rpole2,fill1,fill2,
     $             gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,
     $             obsparm,bdist,ecc,argper,teff1)
c
                if(Teff2.gt.0.0d0)call getATMint(maxlines,maxmu,Nlines,
     &            atmT,atmg,atmmu,Nmu,
     &            atmint1,atmint2,atmint3,atmint4,atmint5,
     &            atmint6,atmint7,atmint8,
     &            Tmax,Tmin,gmax,gmin,gscale2,darkint2,tpole2,gpole2,
     &            dwavex,dwavey,ilaw,iatm,2)
c
                if(teff2.gt.0.0d0)call getrefATMflux(ialphmax2,ibetmax2,
     &             Nalph2,Nbet2,ibetlim2,
     $             visib2,projarray2,temp2,surf2,g2,rinty2,
     &             flum2,maxlines,maxmu,Nlines,
     &             atmT,atmg,atmmu,Nmu,
     &             atmint1,atmint2,atmint3,atmint4,atmint5,
     &             atmint6,atmint7,atmint8,
     &             Tmax,Tmin,gmax,gmin,gscale2,
     &             fluxU2,fluxB2,fluxV2,fluxR2,fluxI2,fluxJ2,fluxH2,fluxK2,
     %             icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,delphi2,
     &             delphie2,iedgestar2,iedgehor2,darkint2,separ,mmdx2,
     &             dwavex,dwavey,ilaw,iatm,2)
c
              endif
              refflux1(1)=fluxU1
              refflux2(1)=fluxU2
              refflux1(2)=fluxB1
              refflux2(2)=fluxB2
              refflux1(3)=fluxV1
              refflux2(3)=fluxV2
              refflux1(4)=fluxR1
              refflux2(4)=fluxR2
              refflux1(5)=fluxI1
              refflux2(5)=fluxI2
              refflux1(6)=fluxJ1
              refflux2(6)=fluxJ2
              refflux1(7)=fluxH1
              refflux2(7)=fluxH2
              refflux1(8)=fluxK1
              refflux2(8)=fluxK2
c
            endif
          endif
c
c
c   UPDATE OCTOBER 20, 2005
c
c   Add an "EBOP mode" where only the flux at even phases are computed.
c   Set the refflux=-99 and filter out the negative values later.
c
c
          if(iecheck.eq.5)then
            do 431 jjj=1,8
              refflux1(jjj)=-99.0d0
              refflux2(jjj)=-99.0d0
431         continue
          endif
c 
          ttiny=0.0d0
          if(ecc.le.0.0d0)ttiny=1.0d-6
c
c   UPDATE JULY 22, 2010
c
c   Add this block to get itime=2 working in eccentric mode.
c
          if((isw7.eq.2).and.(ecc.gt.0.0d0))then
            pstartout=360.0d0*(timearray(1)-T0)/period
            pstopout=360.0d0*(timearray(ntime)-T0)/period
            dphase=360.0d0*tstep/period
            ism1=0
          endif

          do 999 phaseout=pstartout,pstopout+ttiny,dphase          

          if(ecc.gt.0.0d0)then
            ivrt=1
c
c   RVG BUG ALERT  April 19, 2001.
c
c   If the ism1 flag is set, then modify the start and stop phase
c   of the inner loop.  Also, define a variable called pstep.
c           
            if(ism1.eq.0)then
              em=phaseout*pie/180.0d0
c
              call getE(em,ecc,bigE)
c
c   RVG BUG ALERT   April 19, 2001
c
c   Add these blocks to ensure that angles are within the set limits
c
11107         if(bigE.lt.0.0)then
                bigE=bigE+2.0d0*pie
                go to 11107
              endif
22207         if(bigE.gt.2.0d0*pie)then
                bigE=bigE-2.0d0*pie
                go to 22207
              endif
              rnu=2.0d0*datan(dsqrt((1.0d0+ecc)/(1.0d0-ecc))*dtan(bigE/2.0d0))
c
c   RVG BUG ALERT   April 19, 2001
c
c   Add these blocks to ensure that angles are within the set limits
c
11105         if(rnu.lt.0.0)then
                rnu=rnu+2.0d0*pie
                go to 11105
              endif
11106         if(rnu.gt.2.0d0*pie)then
                rnu=rnu-2.0d0*pie
                go to 11106
              endif
              bdist=(1.0d0-ecc*dcos(bigE))
c
c   RVG BUG ALERT   April 19, 2001
c
c   Use the dmod function for pstart and pstop.
c
              pstart=dmod(rnu*180.0d0/pie+argper+90.0d0,360.0d0)
              pstop=dmod(rnu*180.0d0/pie+argper+90.0d0,360.0d0)
              pstep=dphase
c
            endif !if isym = 0
c
            if(ism1.ge.1)then
              em=phaseout*pie/180.0d0
              call getE(em,ecc,bigE)
c
c   RVG BUG ALERT   April 19, 2001
c
c   Add these blocks to ensure that angles are within the set limits
c
11108         if(bigE.lt.0.0)then
                bigE=bigE+2.0d0*pie
                go to 11108
              endif
22208         if(bigE.gt.2.0d0*pie)then
                bigE=bigE-2.0d0*pie
                go to 22208
              endif
              rnu=2.0d0*datan(dsqrt((1.0d0+ecc)/(1.0d0-ecc))*dtan(bigE/2.0d0))
c
c   RVG BUG ALERT   April 19, 2001
c
c   Add these blocks to ensure that angles are within the set limits
c
11109         if(rnu.lt.0.0)then
                rnu=rnu+2.0d0*pie
                go to 11109
              endif
11110         if(rnu.gt.2.0d0*pie)then
                rnu=rnu-2.0d0*pie
                go to 11110
              endif
              bdist1=(1.0d0-ecc*dcos(bigE))
c
c   RVG BUG ALERT   April 19, 2001
c
c   Use the dmod function for pstart and pstop.
c
              pstart=dmod(rnu*180.0d0/pie+argper+90.0d0,360.0d0)
c
c   RVG BUG ALERT  April 19, 2001
c
c   Here is the other phase that gives the same binary separation
c
              emnew=-1.0d0*em
c
c   We need to do the periastron and apastron phases once.
c
              if(phaseout.eq.pstartout)emnew=em  
              if(phaseout.eq.pstopout)emnew=em  
              call getE(emnew,ecc,bigEnew)
c
c   RVG BUG ALERT   April 19, 2001
c
c   Add these blocks to ensure that angles are within the set limits
c
11111         if(bigEnew.lt.0.0)then
                bigEnew=bigEnew+2.0d0*pie
                go to 11111
              endif
22211         if(bigEnew.gt.2.0d0*pie)then
                bigEnew=bigEnew-2.0d0*pie
                go to 22211
              endif
              rnunew=2.0d0*datan(dsqrt((1.0d0+ecc)/(1.0d0-ecc))
     %           *dtan(bigEnew/2.0d0))
c
c   RVG BUG ALERT   April 19, 2001
c
c   Add these blocks to ensure that angles are within the set limits
c
11112         if(rnunew.lt.0.0)then
                rnunew=rnunew+2.0d0*pie
                go to 11112
              endif
11113         if(rnunew.gt.2.0d0*pie)then
                rnunew=rnunew-2.0d0*pie
                go to 11113
              endif
              bdist2=(1.0d0-ecc*dcos(bigEnew))
c
c   RVG BUG ALERT   April 19, 2001
c
c   Use the dmod function for pstart and pstop.
c
              pstop=dmod(rnunew*180.0d0/pie+argper+90.0d0,360.0d0)
              pstep=pstop-pstart
c
              if(pstep.eq.0.0d0)pstep=1.0d0
              bdist=bdist1
c
c   We need to do the periastron and apastron phases only once.
c
              if(phaseout.eq.pstartout)pstop=pstart
              if(phaseout.eq.pstopout)pstop=pstart
            endif ! if isym > 1
c
c   UPDATE September 11, 2001
c
c   Add the iverb flag to setupgeo
c
            iverb=0
c
c  
c
            iskip1=0
            fincr=finc*pie/180.0d0
            phaser=pstart*pie/180.0d0
            delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
            delta=bdist*dsqrt(delta)
            if(delta.gt.(reff1*1.02d0+reff2*1.02d0))iskip1=10
c
54325       format(f8.3,2x,f6.4,1x,4(f8.3,1x),i2,2x,i2)
            iskip2=0
            fincr=finc*pie/180.0d0
            phaser=pstop*pie/180.0d0
            delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
            delta=bdist*dsqrt(delta)
c            write(*,*)delta,reff1,reff2,phaser
            if(delta.gt.reff1*1.02d0+reff2*1.02d0)iskip2=10

c       write(*,54325)phase,delta,phaseout,pstart,pstop,phaser,iskip1,iskip2

c   Now, if idark1>0, then we can skip the phase where star 2 is in front
c   (e.g. phase 180).  Likewise, if idark2 > 0, then skip phases near 0.
c
            if(((phase.ge.0.0d0).and.(phase.le.90.0d0)).or.
     $           ((phase.gt.270.0d0).and.(phase.le.360.0d0)))then
              if(idark2.gt.0)skip1=10
            endif
            if((phase.ge.90.0d0).and.(phase.le.270.0d0))then
              if(idark1.gt.0)skip2=10
            endif
c
c
c    Check for analytic mode
c
c            if((ecc.gt.0.0d0).and.(isw12.ge.1))then    !was eq.9
cc              icount=icount+1
cc              phase=pstart
cc              dummyphase=phase+180.0d0
c              phaser=phase*pie/180.0d0
cc              delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
cc              delta=bdist*dsqrt(delta)
cc              zee=delta/reff1
cc              pee=1.0d0/ratrad
c              write(*,*)phase,phaser,delta,'$$'
c              call analyticg(isw12,ilaw,dwavex,dwavey,delta,reff1,ratrad,
c     #              refflux1,refflux2,phaser,pconj,pconj2,gimvel,fincr,vrot1,
c     &              omega1,period,separ,Q,ecc,bigI,bigbeta)
cc
c              icount=icount+1
c              phase=pstart
c              dummyphase=phase+180.0d0
c              if(iatm.eq.0)go to 938
c              if(iatm.ge.1)go to 8938
c            endif
c            
c UPDATE OCTOBER 20, 2005
c
c Add iecheck=5.  In this case, compute only integer phases.
c
            fflag=-99.0d0
            ddum=phase
            if(ecc.gt.0.0d0)ddum=em*180.0d0/2.0d0/pie
            fdiff=(dabs(ddum*0.5d0-dble(dint(ddum*0.5d0))))
            if(fdiff.lt.0.5d0*dphase)fflag=99.0d0
c            write(*,7676)phase,ddum,fdiff,iskip1,iskip2
c
            if((ecc.gt.0.0d0).and.(iecheck.ge.5).and.
     &            ((iskip1.eq.10).or.(iskip2.eq.10)))then    !was eq.9
              if(fflag.ge.90.0d0)then
                iskip1=0
                iskip2=0
                go to 621
              endif
              icount=icount+1
              phase=pstart
              dummyphase=phase+180.0d0
              if(iatm.eq.0)go to 938
              if(iatm.ge.1)go to 8938
            endif
c
c   UPDATE May 10, 2006
c
c   If we are in analytic mode, and the eccentricity is more than 0.0,
c   skip the setupgeo, etc. since it is not needed.
c
            if((isw12.gt.0).and.(ecc.gt.0.0d0))go to 33321
c
621         call setupgeo(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $       fill1,omega1,Q,finc,
     $       x1,y1,z1,
     %       surf1,rad1,gradx1,grady1,gradz1,g1,xend1,separ,
     %       Tgrav1,Teff1,reff1,Rl1,Tpole1,Rpol1,Regg1,SA1,pot1,gpole1,
     %       phiar1,isquare,iusepot,usepot1,ivrt,pervol1,fillper1,bdist,
     #       pots1,iverb,mmdx1,primmass,primK,primrad,ratrad,frac1,frac2,
     &       ecc,period,size1,fill2,omegatwo,sw5,tteff2,density,
     $       tidephi,itide,phistart1)
            if(teff2.gt.0.0d0)then
              call setupgeo(2,ialphmax2,ibetmax2,Nalph2,Nbet2,ibetlim2,
     $          fill2,omega2,Q,finc,
     $          x2,y2,z2,
     %          surf2,rad2,gradx2,
     &          grady2,gradz2,g2,xend2,separ,
     %          Tgrav2,Teff2,reff2,Rl2,Tpole2,Rpol2,Regg2,SA2,pot2,gpole2,
     $          phiar2,isquare,iusepot,usepot2,ivrt,pervol2,fillper2,bdist,
     $          pots2,iverb,mmdx2,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecc,period,size1,fill2,omegatwo,sw5,tteff2,density,
     $          tidephi,itide,phistart2)
            else
c
c   UPDATE June 14, 2002
c
c   Add ibetlim2 to the argument list of dummyvalues
c
              call dummyvalues(ialphmax2,ibetmax2,Nalph2,Nbet2,
     $          x2,y2,z2,surf2,rad2,gradx2,grady2,gradz2,g2,xend2,darkbol2,
     $          temp2,ibetlim2,mmdx2)
              Regg2=rocheradius(Q)
              reff2=Regg2
              overQ=1.0d0/Q
              call findL1(overQ,omega2,x0,1,bdist,tidephi,itide)
              Rl2=x0
            endif
          endif   !end if ecc.gt.0
c
          if(idint.ge.1)call disksetup(Nthetamax,Nrmax,Ntheta,Nradius,
     %       betarim,rinner,router,Regg2,Rl2,separ,
     $       tdisk,xi,dtemp,dx,dy,dz,drad,
     $       tedge,xedge,yedge,zedge,redge,stepr,stepz,bdist,
     #       ivrt,reper,rsper)
c
c   Find the temperatures in the absence of reflection (heating).  The
c   modified temperatures do not depend on phase.
c  
c
          call setuptemp(1,
     &       ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     &       g1,Tpole1,Tgrav1,temp1,gpole1,mmdx1)

          if(teff2.gt.0.0d0)call setuptemp(2,ialphmax2,ibetmax2,
     %            Nalph2,Nbet2,ibetlim2,g2,Tpole2,Tgrav2,temp2,gpole2,
     $            mmdx2)

c
c   RVG BUG ALERT  May 3, 2001 
c
c   Initialize the ratios and copy the temperatures only if Nref>0
c
          if(Nref.ge.-1)then
            call copytemp(ialphmax1,ibetmax1,Nalph1,Nbet1,temp1,
     #          tempold1,mmdx1,ibetlim1)
            call copytemp(ialphmax2,ibetmax2,Nalph2,Nbet2,temp2,
     #          tempold2,mmdx2,ibetlim2)
c
c   Compute the reflection effect and find the modified temperatures. First
c   initialize the ratios...
c
            call initratio(ialphmax1,ibetmax1,
     %           Nalph1,Nbet1,Nalph2,Nbet2,ratio1,ratio2,
     *           coprat1,coprat2,ialphmax2,ibetmax2)
          endif
c
c   UPDATE March 25, 2002
c
c   I have added another routine called simplerefl, which is a much
c   more efficient routine when the stars are nearly point sources.
c   If Nref=0, call this routine and skip to statement 888.  If
c   Nref < 0, skip both.
c
c   If Teff2 < 0, then star 2 is a point source.  Hence we need only 1
c   iteration of the reflection effect.
c
c   Also, if there is no star 2, call simplerefl.
c
          if((Teff2.le.0.0d0).and.(Nref.gt.0))Nref=0
c
          if(Nref.lt.0)go to 888
          if(Nref.eq.0)then
            call simplerefl(ialphmax1,ibetmax1,
     $        Nalph1,ibetlim1,Nalph2,ibetlim2,
     $       x1,y1,z1,gradx1,grady1,gradz1,g1,
     $       x2,y2,z2,gradx2,grady2,gradz2,g2,
     $       temp1,temp2,dbolx,dboly,ilaw,alb1,alb2,teff1,
     $       teff2,Tgrav1,Tgrav2,rLx,idint,redge,betarim,gpole1,gpole2,
     &       Tpole1,Tpole2,bdist,SA1,SA2,rad1,rad2,separ,mmdx1,
     &       mmdx2,ialphmax2,ibetmax2,isw25)
            go to 888
          endif
c
          do 8 jj=1,Nref
c
            call detailrefl(ialphmax1,ibetmax1,
     $        Nalph1,Nbet1,ibetlim1,Nalph2,Nbet2,ibetlim2,ratio1,ratio2,
     $       x1,y1,z1,gradx1,grady1,gradz1,g1,surf1,
     $       x2,y2,z2,gradx2,grady2,gradz2,g2,surf2,
     $       temp1,temp2,tempold1,tempold2,dbolx,dboly,ilaw,alb1,alb2,teff1,
     $       teff2,Tgrav1,Tgrav2,rLx,idint,redge,betarim,gpole1,gpole2,
     &       Tpole1,Tpole2,coprat1,coprat2,bdist,mmdx1,
     &       mmdx2,ialphmax2,ibetmax2)
c
c   UPDATE DECEMBER 17, 2001
c
c   Add if(tdisk.gt.0.0)
c
c
c   UPDATE June 14, 2002
c
c   Comment out this subroutine call for now.
c
c            if(tdisk.gt.0.0d0)then
c              if((idint.ge.1).and.(jj.eq.1))call diskrefl(ialphmax1,ibetmax1,
c     $         Nalph1,Nbet1,ibetlim1,Nalph2,Nbet2,ibetlim2,ratio1,ratio2,
c     $         x1,y1,z1,gradx1,grady1,gradz1,g1,surf1,
c     $         x2,y2,z2,gradx2,grady2,gradz2,g2,surf2,
c     $         temp1,temp2,tempold1,tempold2,dbolx,dboly,ilaw,alb1,alb2,teff1,
c     $         teff2,Tgrav1,Tgrav2,rLx,idint,gpole1,gpole2,
c     %         Tpole1,Tpole2,coprat1,coprat2,Nthetamax,Nrmax,Ntheta,Nradius,
c     %         betarim,rinner,router,reff2,Rl2,separ,
c     $         tdisk,xi,dtemp,dx,dy,dz,drad,
c     $         tedge,xedge,yedge,zedge,redge,stepr,stepz,dratio,
c     %         reper,rsper,ialphmax2,ibetmax2)
c            endif
 8        continue
c
c   RVG BUG ALERT  May 8, 2001   
c
c   Add spots here, if any.
c
c   UPDATE March 25, 2002
c
c   Add a statement label 888 below.
c   
c   UPDATE March 26, 2002
c
c   Get rid of Nbet from the argument list of addstarspot.  Its value
c   is contained within ibetlim.
c
c   UPDATE December 9, 2008
c
c   Set the averages to zero first.
c
888      continue
          ave11=0.0d0
          ave12=0.0d0
          ave21=0.0d0
          ave22=0.0d0
          if(ispot1.gt.0)call addstarspot(1,ialphmax1,ibetmax1,
     $       Nalph1,ibetlim1,temp1,spot1parm,ave11,ave12,omega1,phase,
     %       phiar1,mmdx1,ispotprof)
          if((ispot2.gt.0).and.(teff2.gt.0.0d0))call 
     %       addstarspot(2,ialphmax2,ibetmax2,
     $       Nalph2,ibetlim2,temp2,spot2parm,ave21,ave22,omega2,phase,
     %       phiar2,mmdx2,ispotprof)
c
c   UPDATE DECEMBER 17, 2001
c
c   Add if(tdisk.gt.0.0)
c
c
c   UPDATE March 26, 2002
c
c   Remove betarim, separ, tdisk, xi, dx, dy, dz, drad, 
c   xedge, yedge, zedge, stepr, stepz, bdist, omega, phase
c   from the argument list of adddiskspot.
c
          ave1=0.0d0
          ave2=0.0d0
          if((ispotd.gt.0).and.(idint.gt.0).and.(tdisk.gt.0.0d0))call
     %       adddiskspot(Nthetamax,Nrmax,Ntheta,Nradius,
     %       rinner,router,reff2,Rl2,
     $       dtemp,
     $       tedge,redge,
     #       ivrt,reper,rsper,spotdparm,ave1,ave2)
c
c
c   We can now compute interesting system parameters based in the input
c   numbers.
c
          rpole1=Rpol1(Nalph1/2)
          rpole2=Rpol2(Nalph2/2)
c
          call parms(1,
     &        Teff2,Q,finc,separ,period,reff1,reff2,
     %        rpole1,rpole2,fill1,fill2,
     $        gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,bdist,ecc)
c
          call parms1(Teff2,Q,finc,separ,period,reff1,reff2,
     %        rpole1,rpole2,fill1,fill2,
     $        gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,
     $        obsparm,bdist,ecc,argper,teff1)
c
c   NEW BUG ALERT  July 13, 2001
c
c   Here is a new subroutine call.  The new subroutine is at the end.
c
c   UPDATE January 16, 2001
c
c   add pot1,pot2 to the end of the list
c
c   UPDATE March 26, 2002
c
c   Remove rpole1,rpole2,fill1,fill2 from the argument list.
c
c
c   UPDATE MAy 10, 2006
c
c   Here is the point to skip to if the orbit is eccentric and we
c   are in analytic mode.
c
33321       continue
c
          call lineparms1(Teff2,Q,
     %        finc,separ,period,reff1,reff2,
     $        vrot1,vrot2,omega1,omega2,
     $        bdist,ecc,SA3,ave11,ave12,ave21,ave22,ave1,ave2,
     #        parmstring,pot1,pot2)
c
c
c   Output the necessary information to make contour plots of the temperature
c   and gravity.
c
c
c   RVG BUC ALERT  June 12 2001
c
c   Add the ".and.idraw.ge.1" condition to the if-then block
c
c   UPDATE June 7, 2002
c
c   Add the x, y, and z coordinate arrays to the argument of writetempgrav
c
c   UPDATE March 4, 2010
c
c   Add the separ to the argument list
c
          if((ecc.eq.0.0d0).and.(idraw.ge.1))then
            call writetempgrav(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %           temp1,g1,gscale1,1,x1,y1,z1,mmdx1,phistart1,separ)
            if(teff2.gt.0.0d0)then
              call writetempgrav(ialphmax2,ibetmax2,Nalph2,Nbet2,ibetlim2,
     %           temp2,g2,gscale2,2,x2,y2,z2,mmdx2,phistart2,separ)
            endif
          endif
c
c   We have to make sure we have the correct light curve so we can compute
c   the radial velocity curve if necessary.  If icnRV1.ne.430
c   or if icnRV2.ne.430, then we need to set icn? temporarily to 1.
c
          if(iRVfilt.eq.1)iVsave=icnU
          if(iRVfilt.eq.2)iVsave=icnB
          if(iRVfilt.eq.3)iVsave=icnV
          if(iRVfilt.eq.4)iVsave=icnR
          if(iRVfilt.eq.5)iVsave=icnI
          if(iRVfilt.eq.6)iVsave=icnJ
          if(iRVfilt.eq.7)iVsave=icnH
          if(iRVfilt.eq.8)iVsave=icnK
c
c   UPDATE October 22, 2008
c
c   Always compute at the filter corresponding to iRVfilt
c

c          if((icnRV1.ne.430).or.(icnRV2.ne.430))then
            if(iRVfilt.eq.1)icnU=1
            if(iRVfilt.eq.2)icnB=1
            if(iRVfilt.eq.3)icnV=1
            if(iRVfilt.eq.4)icnR=1
            if(iRVfilt.eq.5)icnI=1
            if(iRVfilt.eq.6)icnJ=1
            if(iRVfilt.eq.7)icnH=1
            if(iRVfilt.eq.8)icnK=1
c          endif
c
c   Initialize the disk correction matrix
c
c   UPDATE September 25, 2008
c
c   Add the if statement to the correction statement.
c
c
          do 88 kk=1,8
            if(idcheck.gt.0)zdcorr(kk)=0.0d0
            darkint1(kk)=1.0d0
            darkint2(kk)=1.0d0
 88       continue
c
c   If we are using the atmosphere table, compute the DINT factors.
c
          if(iatm.ge.1)then
c
c
c
            call getATMint(maxlines,maxmu,Nlines,
     &        atmT,atmg,atmmu,Nmu,
     &        atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &        Tmax,Tmin,gmax,gmin,gscale1,darkint1,tpole1,gpole1,
     $        dwavex,dwavey,ilaw,iatm,1)
            if(Teff2.gt.0.0d0)call getATMint(maxlines,maxmu,Nlines,
     &        atmT,atmg,atmmu,Nmu,
     &        atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &        Tmax,Tmin,gmax,gmin,gscale2,darkint2,tpole2,gpole2,
     &        dwavex,dwavey,ilaw,iatm,2)
           endif
c
c   Loop over phases and turn the binary in space.
c
c   If ism1=1, then do only phases 0.0 to 180.0 and reflect the light curve
c   at the end.
c
          if(ionephase.ge.1)then
            pstart=onephase
            pstop=onephase
          endif
c         
          ipstep=0
          if((isw7.eq.2).and.(ecc.eq.0.0d0))then
            pstart=360.0d0*(timearray(1)-T0)/period
            pstop=360.0d0*(timearray(ntime)-T0)/period
            pstep=360.0d0*tstep/period
            ism1=0
          endif
c          
          do 10 phasein=pstart,pstop+ttiny,pstep
c
            if(ecc.eq.0.0d0)iskip2=10
c
c
c   RVG BUG ALERT  April 19, 2001
c
c   Put in a flag for the case when ecc>0 and ism1>0 to let the program
c   know which phase is being done
c
            ipstep=ipstep+1
c
c   END BUG
c
            phase=dmod(phasein,360.0d0)
            if(phase.lt.0.0d0)phase=phase+360.0d0

            icount=icount+1
c
            if(ecc.eq.0.0d0)then
              if((ispot1.ge.1).and.(isw7.ge.2))then
                call copytemp(ialphmax1,ibetmax1,Nalph1,Nbet1,tempold1,
     #           temp1,mmdx1,ibetlim1)
               call addmovespot(1,ialphmax1,ibetmax1,
     $           Nalph1,ibetlim1,temp1,spot1parm,ave11,ave12,omega1,
     %           phiar1,mmdx1,period,t0,timearray(icount))
              endif
          
              if((ispot2.ge.1).and.(isw7.ge.2).and.(teff2.gt.0.0d0))then
                call copytemp(ialphmax2,ibetmax2,Nalph2,Nbet2,tempold2,
     #            temp2,mmdx2,ibetlim2)
                call addmovespot(2,ialphmax2,ibetmax2,
     $            Nalph2,ibetlim2,temp2,spot2parm,ave21,ave22,omega2,
     %            phiar2,mmdx2,period,t0,timearray(icount))
              endif
            endif

            tphase=phase+360.0d0*(pshift+eshift)
c
c
c   RVG BUG ALERT  April 19, 2001
c
c   Put as the argument to getextension  dmod(extphase,360.0d0), and modify
c   the phase according to the phase shifts.
c
            extphase=phase
            if((ecc.gt.0.0d0).or.(pshift.ne.0.0d0))then
              if((ecc.gt.0.0d0).and.(ism1.eq.0))emphase=180.0d0*em/pie
              if((ecc.gt.0.0d0).and.(ism1.gt.0))then
                if(ipstep.eq.1)emphase=180.0d0*em/pie
                if(ipstep.eq.2)emphase=180.0d0*emnew/pie
              endif
              tshift=pshift+eshift
              extphase=emphase+360.0d0*tshift
c
c   UPDATE September 10, 2001
c
c   Add the if-then clauses
c
              if(ikeep.eq.1)extphase=emphase+360.0d0*(tshift-pconj)
              if(ikeep.eq.2)extphase=emphase+360.0d0*(tshift-pconj2)
            endif
            call getextension(dmod(extphase,360.0d0),extension)
            dummyphase=dmod(phase+180.0d0,360.0d0)  ! this is for star 2
c
c   Check the various toggle switches which tell the code whether to check
c   for eclipses.  For example, if points were eclipsed up until phase
c   15.0, then the code will not check again until phase=180-15.  It
c   will keep checking until phase=180+15 and not check again until
c   phase=360-15.
c
            pdiff1=dabs(phase-(360.0d0-togglephase))
            if(pdiff1.lt.0.00001d0)then
              if(iecheck.ge.1)idcheck=100
c              itoggle=-1
            endif
c
            iskip1=0
            iskip2=0             
            fincr=finc*pie/180.0d0
            phaser=phase*pie/180.0d0
            delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
            delta=bdist*dsqrt(delta)
            if(delta.gt.(reff1*1.07d0+reff2*1.07d0))then    ! was 1.02
              iskip1=10
              iskip2=10
             endif
c            write(*,7777)delta,phaser,iskip1,iskip2
 7777       format('***   ',2(f9.6,1x),i3,1x,i3)
c
c
c   Now, if idark1>0, then we can skip the phase where star 2 is in front
c   (e.g. phase 180).  Likewise, if idark2 > 0, then skip phases near 0.
c
            if(((phase.ge.0.0d0).and.(phase.lt.90.0d0)).or.
     $          ((phase.gt.270.0d0).and.(phase.le.360.0d0)))then
              if(idark2.gt.0)skip1=10
            endif

            if((phase.ge.90.0d0).and.(phase.le.270.0d0))then
              if(idark1.gt.0)skip2=10
            endif
c
c    Check for analytic mode
c
            if(isw12.ge.1)then    !was eq.9
              phaser=phase*pie/180.0d0
              delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
              delta=bdist*dsqrt(delta)
              zee=delta/reff1
              pee=1.0d0/ratrad
c              write(*,4242)phase,delta,reff1
4242          format('phase = ',f7.2,1x,3(f11.7,2x))
              call analyticg(isw12,ilaw,dwavex,dwavey,delta,reff1,ratrad,
     #              refflux1,refflux2,phaser,pconj,pconj2,gimvel,fincr,vrot1,
     #              omega1,period,separ,Q,ecc,bigI,bigbeta,Neclipse1,1,corr1,corr2)
c
c              write(*,*)'Neclipse = ',Neclipse1

              if(iatm.eq.0)go to 938
              if(iatm.ge.1)go to 8938
            endif
c
c   UPDATE OCTOBER 20, 2005
c
c   change to ge.5 (if iecheck=5, then compute the integer phases)
c
            fflag=-99.0d0
            ddum=phase
            if(ecc.gt.0.0d0)ddum=em*180.0d0/2.0d0/pie
            fdiff=(dabs(ddum*0.5d0-dble(dint(ddum*0.5d0))))
            if(fdiff.lt.0.5d0*dphase)fflag=99.0d0
            if((iecheck.ge.5).and.(iskip1.eq.10))then
              if(fflag.ge.90.0d0)then
                iskip1=0
                iskip2=0
                go to 821
              endif
              if(iatm.eq.0)go to 938
              if(iatm.ge.1)go to 8938
            endif
821         call gethorizon(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %        phase,finc,Q,
     $        pot1,omega1,x1,y1,z1,rad1,gradx1,grady1,gradz1,xend1,
     $        Nhoriz1,xhoriz1,yhoriz1,phiar1,iedgestar1,delphie1,bdist,
     &        mmdx1,xhmin1,xhmax1,yhmin1,yhmax1,tidephi,itide,phihor1)
c
            if(teff2.gt.0.0d0)then
              call gethorizon(2,ialphmax2,ibetmax2,Nalph2,
     %            Nbet2,ibetlim2,dummyphase,finc,Q,
     $            pot2,omega2,x2,y2,z2,rad2,gradx2,grady2,gradz2,xend2,
     $            Nhoriz2,xhoriz2,yhoriz2,phiar2,iedgestar2,delphie2,bdist,
     &            mmdx2,xhmin2,xhmax2,yhmin2,yhmax2,tidephi,itide,phihor2)
              if(idint.ge.1)call gettophorizon(2,
     %            ialphmax2,ibetmax2,Nalph2,Nbet2,ibetlim2,dummyphase,finc,Q,
     $            pot2,omega2,x2,y2,z2,rad2,gradx2,grady2,gradz2,xend2,
     $            Ntop2,xtop2horiz,ytop2horiz,phiar2,iedgestar2,delphie2,
     &            bdist,mmdx2)
            else
              call dummyhoriz(ibetmax2,Nbet2,Nhoriz2,xhoriz2,yhoriz2,
     &            Ntop2,xtop2horiz,ytop2horiz)
            endif
c
c    UPDATE May 26, 2004
c
c    If iecheck=9, check to see if the horizon of star 1 overlaps
c    with the horizon of star 2.  If so, then jump to the escape point.
c    ioverlap=999 means there is overlap
c
            ioverlap=-999
            if(iecheck.eq.9)then
              call overlaphoriz(Nhoriz1,xhoriz1,yhoriz1,
     $            Nhoriz2,xhoriz2,yhoriz2,ioverlap)
              if(ioverlap.lt.900)then
                iskip1=10
                iskip2=10
c                write(*,*)phase,'************'
                if(iatm.eq.0)go to 938
                if(iatm.ge.1)go to 8938
              endif
            endif

c
c   UPDATE September 10, 2001
c
c   If Teff < 0 check to see if the center of star 2 is eclipsed.
c
            if(Teff2.le.0.0d0)then
              pppp=dmod(phase,360.0d0)
              if(((pppp.gt.-90.0d0).and.(pppp.lt.90.0d0)).or.
     %           ((pppp.gt.270.0d0).and.(pppp.le.450.0d0)))then
                call getXecl(Nhoriz1,xhoriz1,yhoriz1,ixecl,Q,finc,bdist,phase)
              else
                ixecl=-100
              endif
              xecx(icount)=phase !dmod(phase+360.0d0*(pconj2),360.0d0)
              xecy(icount)=dble(ixecl)
            endif
c
            Ndhoriz=0
            Ndtop=0
            
            if(idint.ge.1)call getdiskhoriz(Nthetamax,Nrmax,Ntheta,Nradius,
     %       Q,phase,finc,xedge,yedge,zedge,Ndhoriz,dxhoriz,dyhoriz,
     %       Ndtop,dtopx,dtopy,bdist)

c
c   Check the visibilities of grid elements.  Note that horizon 2 goes
c   in the argument list for star 1, and vice-versa.
c           
c   UPDATE JULY 2, 2004
c
c   Add the variable jdum to the argument list of get visib.
c
c           
c   UPDATE JULY 4, 2004
c
c   Add the variable MonteCarlo to the argument list of get visib.
c
c  
c
            call getvisib(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     &        phase,finc,Q,pot1,omega1,gradx1,grady1,gradz1,x1,y1,z1,
     %        xend1,visib1,Nhoriz2,xhoriz2,yhoriz2,idint,
     %        Ndhoriz,dxhoriz,dyhoriz,Ndtop,dtopx,dtopy,
     $        Nsky1,xsky1,ysky1,projarray1,iecheck,Neclipse1,phiar1,rad1,
     $        delphi1,separ,iedgehor1,bdist,mmdx1,jdum,MonteCarlo,
     #        isw13,ialfmin,ialfmax,xhmin2,xhmax2,yhmin2,yhmax2,
     $        tidephi,itide,iedgestar1,phihor1,phistart1)
c
c
            if(teff2.gt.0.0d0)call getvisib(2,ialphmax2,ibetmax2,
     $        Nalph2,Nbet2,ibetlim2,dummyphase,finc,Q,pot2,omega2,
     $        gradx2,grady2,gradz2,x2,y2,z2,
     $        xend2,visib2,Nhoriz1,xhoriz1,yhoriz1,idint,
     $        Ndhoriz,dxhoriz,dyhoriz,Ndtop,dtopx,dtopy,
     $        Nsky2,xsky2,ysky2,projarray2,iecheck,Neclipse2,phiar2,rad2,
     $        delphi2,separ,iedgehor2,bdist,mmdx2,jdum,MonteCarlo,
     #        isw13,ialfmin,ialfmax,xhmin1,xhmax1,yhmin1,yhmax1,
     $        tidephi,itide,iedgestar2,phihor2,phistart2)
c
            if(idint.ge.1)call diskvisib(Nrmax,Nthetamax,Nradius,
     %        Ntheta,phase,finc,Q,
     %        betarim,dx,dy,dz,xedge,yedge,zedge,
     &        diskproj,edgeproj,dvisib,evisib,Nskydisk,xskydisk,yskydisk,
     &        zskydisk,
     &        Nskyedge,xskyedge,yskyedge,Ntop2,
     %        xtop2horiz,ytop2horiz,Nhoriz1,xhoriz1,yhoriz1,
     &        Ndtop,dtopx,dtopy,iecheck,Neclipsed,bdist)
c
c   If we are in blackbody mode, we need to
c   loop over filters and find the fluxes.   Otherwise, we can call
c   the atmosphere routines a single time and get the fluxes for the
c   8 filters all at once.
c
            if(iatm.le.0)then
              do 9 jj=1,8
                www=wave(jj)
                flimbx=dwavex(jj,1)
                flimby=dwavey(jj,1)
                flux1=0.0d0
                if(idark1.le.0)then
                  if(isimp.eq.0)then
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.
c
c   UPDATE JULY 4, 2004
c
c   Add MonteCarlo to the argument list
c
                    fluxlat=rlatflux(jj)
                    call getBBflux(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $               www,visib1,projarray1,temp1,surf1,flimbx,flimby,ilaw,
     %               rinty1,flum1,flux1,delphi1,delphie1,
     #               iedgestar1,iedgehor1,rldint1,separ,mmdx1,MonteCarlo,
     $               isw13,ialfmin,ialfmax,fluxlat,1,phiar1,phihor1)
                  else
c
c  UPDATE March 26, 2002
c
c  Get rid of the variables phiar, jj, istar from the argument list
c  of getBBsimp.
c
                    call getBBsimp(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $                www,visib1,projarray1,temp1,surf1,flimbx,flimby,
     &                ilaw,rinty1,flum1,flux1,delphi1,delphie1,
     &                iedgestar1,iedgehor1,rldint1,isimp,separ,mmdx1)
                  endif
                endif
c
c   UPDATE May 24, 2002
c
c   Here is a new subroutine call.
c
c
                if(idraw.eq.1)then
                  if(jj.eq.iRVfilt)call rotkern(ialphmax1,ibetmax1,
     $              Nalph1,ibetlim1,
     &              1,omega1,phase,finc,
     %              Q,flum1,x1,y1,flux1,separ,period,gamma,
     $              rldint1,ecc,argrad,visib1,extension,mmdx1)
                endif
c
                call getvel(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $              1,omega1,phase,finc,
     %              Q,flum1,x1,y1,flux1,separ,period,gamma,vel1,delvel1,
     $              rldint1,ecc,argrad,mmdx1,isw13,ialfmin,ialfmax,fluxlat,
     %              bigI,bigbeta,z1)
c
                flimbx=dwavex(jj,2)
                flimby=dwavey(jj,2)
                flux2=0.0d0
                if(idark2.le.0)then
                  if(teff2.gt.0.0d0)then
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.
c
c   UPDATE JULY 4, 2004
c
c   Add MonteCarlo to the argument list.
c
                    if(isimp.eq.0)then
                      fluxlat=0.0d0
                      call getBBflux(ialphmax2,ibetmax2,Nalph2,Nbet2,
     $                ibetlim2,www,visib2,projarray2,temp2,
     &                surf2,flimbx,flimby,ilaw,
     &                rinty2,flum2,flux2,delphi2,delphie2,
     #                iedgestar2,iedgehor2,rldint2,separ,mmdx2,MonteCarlo,
     &                isw13,ialfmin,ialfmax,fluxlat,2,phiar2,phihor2)
                    else
c
c  UPDATE March 26, 2002
c
c  Get rid of the variables phiar, jj, istar from the argument list
c  of getBBsimp.
c
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.
c
                      call getBBsimp(ialphmax2,ibetmax2,Nalph2,Nbet2,ibetlim2,
     $                  www,visib2,projarray2,temp2,surf2,flimbx,flimby,
     &                  ilaw,rinty2,flum2,flux2,delphi2,delphie2,
     &                  iedgestar2,iedgehor2,rldint2,isimp,separ,mmdx2)
                    endif
                  endif
                endif
c
c  UPDATE August 16, 2001
c
c  Remove the if clause.  If there is no star 2 (as in an X-ray
c  binary, then the sine curve will still be computed.
c
c                if(teff2.gt.0.0d0)
c
                call getvel(ialphmax2,ibetmax2,Nalph2,Nbet2,
     &              ibetlim2,2,omega2,dummyphase,
     %              finc,Q,flum2,x2,y2,flux2,separ,period,gamma,vel2,delvel2,
     %              rldint2,ecc,argrad+pie,mmdx2,
     #              isw13,ialfmin,ialfmax,fluxlat,bigI,bigbeta,z2)
c
c
c   UPDATE May 24, 2002
c
c   Here is a new subroutine call.
c
c
                if(idraw.eq.1)then
                  if(jj.eq.iRVfilt)call rotkern(ialphmax2,ibetmax2,
     $              Nalph2,ibetlim2,
     &              2,omega2,phase,finc,
     %              Q,flum2,x2,y2,flux2,separ,period,gamma,
     $              rldint2,ecc,argrad,visib2,extension,mmdx2)
                endif
c
                dflux=0.0d0
c
c   UPDATE December 17, 2001
c
c   If tdisk is less than 0, and idint is more than 1, then
c   compute the disk geometry and eclipses, but ignore its flux.
c
                if(tdisk.gt.0.0d0)then
                  if(idcheck.ge.1)then
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.
c
                    if(idint.ge.1)call getdiskBBflux(Nrmax,Nthetamax,Nradius,
     %               Ntheta,diskproj,edgeproj,dvisib,evisib,dtemp,tedge,drad,
     $               dinty,einty,stepr,stepz,www,jj,ilaw,dflux,separ)
c
                    zdcorr(jj)=dflux
                  else
                     dflux=zdcorr(jj)
                  endif
                else
                  dflux=0.0d0
                endif
c
c   UPDATE OCTOBER 20, 2005
c
c   Change to .ge.5
c
c
 938            if((iecheck.ge.5).or.(isw12.ge.1))then
                  if((
     &               ((ecc.eq.0.0d0).and.(iskip1.eq.10).and.(iskip2.eq.10))
     &               .or.
     &               ((ecc.gt.0.0d0).and.((iskip1.eq.10).or.(iskip2.eq.10))))
     &               .or.(isw12.ge.1))then 
                     do 9938 ll=1,8
c                       icount=icount+1
c   Add Doppler boosting
c
c
                       call getrefvel(1,omega1,phase,finc,
     %                    Q,separ,period,gamma,vel1,ecc,argrad,isw12,gimvel,
     #                     iRVfilt)
                       call getrefvel(2,omega2,dummyphase,finc,
     %                    Q,separ,period,gamma,vel2,ecc,argrad+pie,isw12,
     #                    gimvel,iRVfilt)
c
                       if(ll.eq.1)ymodU(icount)=refflux1(1)+refflux2(1)+corr1(1)+corr2(1)!+dflux
                       if(ll.eq.2)ymodB(icount)=refflux1(2)+refflux2(2)+corr1(2)+corr2(2)!+dflux
                       if(ll.eq.3)ymodV(icount)=refflux1(3)+refflux2(3)+corr1(3)+corr2(3)!+dflux
                       if(ll.eq.4)ymodR(icount)=refflux1(4)+refflux2(4)+corr1(4)+corr2(4)!+dflux
                       if(ll.eq.5)ymodI(icount)=refflux1(5)+refflux2(5)+corr1(5)+corr2(5)!+dflux
                       if(ll.eq.6)ymodJ(icount)=refflux1(6)+refflux2(6)+corr1(6)+corr2(6)!+dflux
                       if(ll.eq.7)ymodH(icount)=refflux1(7)+refflux2(7)+corr1(7)+corr2(7)!+dflux
                       if(ll.eq.8)ymodK(icount)=refflux1(8)+refflux2(8)+corr1(8)+corr2(8)!+dflux
                       yeclipse(icount)=dble(Neclipse1)
                       RV1(icount)=vel1
                       RV2(icount)=vel2
                       if(ll.eq.iRVfilt)then
                         ymods1(icount)=refflux1(ll)
                         ymods2(icount)=refflux2(ll)
                         drV1(icount)=gimvel(ll)
                       endif
c
c   UPDATE January 12, 2009
c
c   make fracs fracs1, fracs2, fracs3, ... fracs8
c
                       if(ll.eq.1)then
                         fracs1(icount,1)=refflux1(ll)
                         fracs1(icount,2)=refflux2(ll)
                         fracs1(icount,3)=dflux
                       endif
                       if(ll.eq.2)then
                         fracs2(icount,1)=refflux1(ll)
                         fracs2(icount,2)=refflux2(ll)
                         fracs2(icount,3)=dflux
                       endif
                       if(ll.eq.3)then
                         fracs3(icount,1)=refflux1(ll)
                         fracs3(icount,2)=refflux2(ll)
                         fracs3(icount,3)=dflux
                       endif
                       if(ll.eq.4)then
                         fracs4(icount,1)=refflux1(ll)
                         fracs4(icount,2)=refflux2(ll)
                         fracs4(icount,3)=dflux
                       endif
                       if(ll.eq.5)then
                         fracs5(icount,1)=refflux1(ll)
                         fracs5(icount,2)=refflux2(ll)
                         fracs5(icount,3)=dflux
                       endif
                       if(ll.eq.6)then
                         fracs6(icount,1)=refflux1(ll)
                         fracs6(icount,2)=refflux2(ll)
                         fracs6(icount,3)=dflux
                       endif
                       if(ll.eq.7)then
                         fracs7(icount,1)=refflux1(ll)
                         fracs7(icount,2)=refflux2(ll)
                         fracs7(icount,3)=dflux
                       endif
                       if(ll.eq.8)then
                         fracs8(icount,1)=refflux1(ll)
                         fracs8(icount,2)=refflux2(ll)
                         fracs8(icount,3)=dflux
                       endif
c

 9938                continue
                     if(ecc.eq.0.0d0)go to 9939
                     if(ecc.gt.0.0d0)go to 9939
                  endif
                endif
                dop1=(vel1-gamma)/2.99792458d5
                dop2=(vel2-gamma)/2.99792458d5
                if(jj.eq.iRVfilt)then
                  flux1=flux1*(1.0d0-beam1*dop1)
                  flux2=flux2*(1.0d0-beam2*dop2)
                endif

                if(jj.eq.1)ymodU(icount)=flux1+flux2+dflux
                if(jj.eq.2)ymodB(icount)=flux1+flux2+dflux
                if(jj.eq.3)ymodV(icount)=flux1+flux2+dflux
                if(jj.eq.4)ymodR(icount)=flux1+flux2+dflux
                if(jj.eq.5)ymodI(icount)=flux1+flux2+dflux
                if(jj.eq.6)ymodJ(icount)=flux1+flux2+dflux                 
                if(jj.eq.7)ymodH(icount)=flux1+flux2+dflux
                if(jj.eq.8)ymodK(icount)=flux1+flux2+dflux
                yeclipse(icount)=dble(Neclipse1)
c
                if(jj.eq.1)then
                  fracs1(icount,1)=flux1
                  fracs1(icount,2)=flux2
                  fracs1(icount,3)=dflux
                endif
                if(jj.eq.2)then
                  fracs2(icount,1)=flux1
                  fracs2(icount,2)=flux2
                  fracs2(icount,3)=dflux
                endif
                if(jj.eq.3)then
                  fracs3(icount,1)=flux1
                  fracs3(icount,2)=flux2
                  fracs3(icount,3)=dflux
                endif
                if(jj.eq.4)then
                  fracs4(icount,1)=flux1
                  fracs4(icount,2)=flux2
                  fracs4(icount,3)=dflux
                endif
                if(jj.eq.5)then
                  fracs5(icount,1)=flux1
                  fracs5(icount,2)=flux2
                  fracs5(icount,3)=dflux
                endif
                if(jj.eq.6)then
                  fracs6(icount,1)=flux1
                  fracs6(icount,2)=flux2
                  fracs6(icount,3)=dflux
                endif
                if(jj.eq.7)then
                  fracs7(icount,1)=flux1
                  fracs7(icount,2)=flux2
                  fracs7(icount,3)=dflux
                endif
                if(jj.eq.8)then
                  fracs8(icount,1)=flux1
                  fracs8(icount,2)=flux2
                  fracs8(icount,3)=dflux
                endif

                if(jj.eq.iRVfilt)then
c                  write(*,*)rldint1
                  ymods1(icount)=flux1
                  ymods2(icount)=flux2
                  ymodd(icount)=dflux
                  RV1(icount)=vel1
                  RV2(icount)=vel2
                  dRV1(icount)=delvel1
                  dRV2(icount)=delvel2
                  fluxV1=flux1
                  fluxV2=flux2
                  call copyinty(ialphmax1,ibetmax1,Nalph1,Nbet1,
     $              Nalph2,Nbet2,rinty1,saveinty1,rinty2,saveinty2,mmdx1,
     &              mmdx2,ibetlim1,ibetlim2,ialphmax2,ibetmax2)
                  if(idint.ge.1)call copydiskinty(Nrmax,Nthetamax,Nradius,
     %              Ntheta,dinty,savedinty,einty,saveeinty)
c
                endif
 9            continue
            endif
c
c   Do the same block for iatm > 0
c
            if(iatm.ge.1)then
c
c   UPDATE November 7, 2008
c
c   zero out the disk fluxes
c
              dfluxU=0.0d0
              dfluxB=0.0d0
              dfluxV=0.0d0
              dfluxR=0.0d0
              dfluxI=0.0d0
              dfluxJ=0.0d0
              dfluxH=0.0d0
              dfluxK=0.0d0

              jj=1
              www=wave(jj)
              flimbx=dwavex(jj,1)
              flimby=dwavey(jj,1)
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.
c
c   UPDATE JULY 4, 2004
c
c   Add MonteCarlo to the argument list.
c
c 
c
              if(idark1.le.0)then
                call getATMflux(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $             visib1,projarray1,temp1,surf1,g1,rinty1,
     &             flum1,maxlines,maxmu,Nlines,
     &             atmT,atmg,atmmu,Nmu,
     &             atmint1,atmint2,atmint3,atmint4,atmint5,
     $             atmint6,atmint7,atmint8,
     &             Tmax,Tmin,gmax,gmin,gscale1,
     &             fluxU1,fluxB1,fluxV1,fluxR1,fluxI1,fluxJ1,fluxH1,fluxK1,
     &             icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,delphi1,
     %             delphie1,iedgestar1,iedgehor1,darkint1,separ,
     #             mmdx1,MonteCarlo,dwavex,dwavey,ilaw,iatm,1)
c
                if(iRVfilt.eq.1)flux1=fluxU1
                if(iRVfilt.eq.2)flux1=fluxB1
                if(iRVfilt.eq.3)flux1=fluxV1
                if(iRVfilt.eq.4)flux1=fluxR1
                if(iRVfilt.eq.5)flux1=fluxI1
                if(iRVfilt.eq.6)flux1=fluxJ1
                if(iRVfilt.eq.7)flux1=fluxH1
                if(iRVfilt.eq.8)flux1=fluxK1
                call getvel(ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $              1,omega1,phase,finc,
     %              Q,flum1,x1,y1,flux1,separ,period,gamma,vel1,delvel1,
     %              darkint1(iRVfilt),ecc,argrad,mmdx1,
     #              isw13,ialfmin,ialfmax,fluxlat,bigI,bigbeta,z1)
c
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.
c
c
c   UPDATE May 24, 2002
c
c   Here is a new subroutine call.
c
c
              endif
              if(idraw.eq.1)then
                call rotkern(ialphmax1,ibetmax1,
     $              Nalph1,ibetlim1,
     &              1,omega1,phase,finc,
     %              Q,flum1,x1,y1,flux1,separ,period,gamma,
     $              darkint1(iRVfilt),ecc,argrad,visib1,extension,mmdx1)
              endif
c
c   UPDATE JULY 4, 2004
c
c   Add MonteCarlo to the argument list.
c
              flux2=0.0d0
              if(idark2.le.0)then
                if(teff2.gt.0.0d0)call getATMflux(ialphmax2,ibetmax2,
     &             Nalph2,Nbet2,ibetlim2,
     $             visib2,projarray2,temp2,surf2,g2,rinty2,
     &             flum2,maxlines,maxmu,Nlines,
     &             atmT,atmg,atmmu,Nmu,
     &             atmint1,atmint2,atmint3,atmint4,atmint5,
     #             atmint6,atmint7,atmint8,
     #             Tmax,Tmin,gmax,gmin,gscale2,
     &             fluxU2,fluxB2,fluxV2,fluxR2,fluxI2,fluxJ2,fluxH2,fluxK2,
     %             icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,delphi2,
     &             delphie2,iedgestar2,iedgehor2,darkint2,separ,mmdx2,
     #             MonteCarlo,dwavex,dwavey,ilaw,iatm,2)
c
                if(iRVfilt.eq.1)flux2=fluxU2
                if(iRVfilt.eq.2)flux2=fluxB2
                if(iRVfilt.eq.3)flux2=fluxV2
                if(iRVfilt.eq.4)flux2=fluxR2
                if(iRVfilt.eq.5)flux2=fluxI2
                if(iRVfilt.eq.6)flux2=fluxJ2
                if(iRVfilt.eq.7)flux2=fluxH2
                if(iRVfilt.eq.8)flux2=fluxK2
              endif
c
c  UPDATE August 16, 2001
c
c  Remove the if clause.  If there is no star 2 (as in an X-ray
c  binary, then the sine curve will still be computed.
c
c              if(teff2.gt.0.0d0)
c
              call getvel(ialphmax2,ibetmax2,Nalph2,
     &            Nbet2,ibetlim2,2,omega2,dummyphase,
     %            finc,Q,flum2,x2,y2,flux2,separ,period,gamma,vel2,delvel2,
     %            darkint2(iRVfilt),ecc,argrad+pie,mmdx2,
     #            isw13,ialfmin,ialfmax,fluxlat,bigI,bigbeta,z2)
c
c
c   UPDATE May 24, 2002
c
c   Here is a new subroutine call.
c
c
              if(idraw.eq.1)then
                call rotkern(ialphmax2,ibetmax2,
     $              Nalph2,ibetlim2,
     &              2,omega2,phase,finc,
     %              Q,flum2,x2,y2,flux2,separ,period,gamma,
     $              darkint2(iRVfilt),ecc,argrad,visib2,extension,mmdx2)
              endif
c
              dflux=0.0d0
c
c   UPDATE DECEMBER 17, 2001
c
c   If tdisk is less than 0, and if idint is more than 0, then
c   the disk is there for geometrical purposes only (ignore its flux).
c
              if(tdisk.gt.0.0d0)then
                if(idcheck.ge.1)then
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.
c
                  if(idint.ge.1)call getdiskATMflux(Nrmax,Nthetamax,Nradius,
     %              Ntheta,diskproj,edgeproj,dvisib,evisib,dtemp,tedge,drad,
     $              dinty,einty,stepr,stepz,maxlines,maxmu,Nlines,
     &              atmT,atmg,atmmu,Nmu,
     &              atmint1,atmint2,atmint3,atmint4,atmint5,
     #              atmint6,atmint7,atmint8,
     #              Tmax,Tmin,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %              dfluxU,dfluxB,dfluxV,dfluxR,
     &              dfluxI,dfluxJ,dfluxH,dfluxK,iRVfilt,separ,
     &              dwavex,dwavey,ilaw,iatm,1)
c
                  zdcorr(1)=dfluxU
                  zdcorr(2)=dfluxB
                  zdcorr(3)=dfluxV
                  zdcorr(4)=dfluxR
                  zdcorr(5)=dfluxI
                  zdcorr(6)=dfluxJ
                  zdcorr(7)=dfluxH
                  zdcorr(8)=dfluxK
                else
                  dfluxU=zdcorr(1)
                  dfluxB=zdcorr(2)
                  dfluxV=zdcorr(3)
                  dfluxR=zdcorr(4)
                  dfluxI=zdcorr(5)
                  dfluxJ=zdcorr(6)
                  dfluxH=zdcorr(7)
                  dfluxK=zdcorr(8)
                endif
              else
                dfluxU=0.0d0
                dfluxB=0.0d0
                dfluxV=0.0d0
                dfluxR=0.0d0
                dfluxI=0.0d0
                dfluxJ=0.0d0
                dfluxH=0.0d0
                dfluxK=0.0d0
              endif
              if(iRVfilt.eq.1)dflux=dfluxU
              if(iRVfilt.eq.2)dflux=dfluxB
              if(iRVfilt.eq.3)dflux=dfluxV
              if(iRVfilt.eq.4)dflux=dfluxR
              if(iRVfilt.eq.5)dflux=dfluxI
              if(iRVfilt.eq.6)dflux=dfluxJ
              if(iRVfilt.eq.7)dflux=dfluxH
              if(iRVfilt.eq.8)dflux=dfluxK
c
9876          format('@#$ ',8(i4,2x))
c
c   UPDATE OCTOBER 20, 2005
c
c   Change to .ge.5
c
 8938         if((iecheck.ge.5))then
                if(((ecc.eq.0.0d0).and.(iskip1.eq.10).and.(iskip2.eq.10)).or.
     &           ((ecc.gt.0.0d0).and.((iskip1.eq.10).or.(iskip2.eq.10))))then 
                  do 5432 ll=1,8
c                       icount=icount+1
                    if(ll.eq.1)ymodU(icount)=refflux1(1)+refflux2(1)+third(1)+corr1(1)+corr2(1)
                    if(ll.eq.2)ymodB(icount)=refflux1(2)+refflux2(2)+third(2)+corr1(2)+corr2(2)
                    if(ll.eq.3)ymodV(icount)=refflux1(3)+refflux2(3)+third(3)+corr1(3)+corr2(3)
                    if(ll.eq.4)ymodR(icount)=refflux1(4)+refflux2(4)+third(4)+corr1(4)+corr2(4)
                    if(ll.eq.5)ymodI(icount)=refflux1(5)+refflux2(5)+third(5)+corr1(5)+corr2(5)
                    if(ll.eq.6)ymodJ(icount)=refflux1(6)+refflux2(6)+third(6)+corr1(6)+corr2(6)
                    if(ll.eq.7)ymodH(icount)=refflux1(7)+refflux2(7)+third(7)+corr1(7)+corr2(7)
                    if(ll.eq.8)ymodK(icount)=refflux1(8)+refflux2(8)+third(8)+corr1(8)+corr2(8)
                    yeclipse(icount)=dble(Neclipse1)
                    call getrefvel(1,omega1,phase,finc,
     %                  Q,separ,period,gamma,vel1,ecc,argrad,isw12,gimvel,
     #                     iRVfilt)
                    call getrefvel(2,omega2,dummyphase,finc,
     %                  Q,separ,period,gamma,vel2,ecc,argrad+pie,isw12,gimvel,
     #                     iRVfilt)
                    RV1(icount)=vel1
                    RV2(icount)=vel2
                    if(ll.eq.iRVfilt)then
                      ymods1(icount)=refflux1(ll)
                      ymods2(icount)=refflux2(ll)
c                      if(ymodV(icount).eq.0.0d0)write(*,*)phase
                    endif
c
c   UPDATE January 12, 2009
c
c   change fracs to fracs1, fracs2, ...  fracs8
c
                    if(ll.eq.1)then
                      fracs1(icount,1)=refflux1(ll)
                      fracs1(icount,2)=refflux2(ll)
                      fracs1(icount,3)=dflux
                    endif
                    if(ll.eq.2)then
                      fracs2(icount,1)=refflux1(ll)
                      fracs2(icount,2)=refflux2(ll)
                      fracs2(icount,3)=dflux
                    endif
                    if(ll.eq.3)then
                      fracs3(icount,1)=refflux1(ll)
                      fracs3(icount,2)=refflux2(ll)
                      fracs3(icount,3)=dflux
                    endif
                    if(ll.eq.4)then
                      fracs4(icount,1)=refflux1(ll)
                      fracs4(icount,2)=refflux2(ll)
                      fracs4(icount,3)=dflux
                    endif
                    if(ll.eq.5)then
                      fracs5(icount,1)=refflux1(ll)
                      fracs5(icount,2)=refflux2(ll)
                      fracs5(icount,3)=dflux
                    endif
                    if(ll.eq.6)then
                      fracs6(icount,1)=refflux1(ll)
                      fracs6(icount,2)=refflux2(ll)
                      fracs6(icount,3)=dflux
                    endif
                    if(ll.eq.7)then
                      fracs7(icount,1)=refflux1(ll)
                      fracs7(icount,2)=refflux2(ll)
                      fracs7(icount,3)=dflux
                    endif
                    if(ll.eq.8)then
                      fracs8(icount,1)=refflux1(ll)
                      fracs8(icount,2)=refflux2(ll)
                      fracs8(icount,3)=dflux
                    endif


 5432            continue
                  if(ecc.eq.0.0d0)go to 9939
                  if(ecc.gt.0.0d0)go to 9939
                endif
              endif
c
c    Doppler boosting
c
              dop1=(vel1-gamma)/2.99792458d5
              dop2=(vel2-gamma)/2.99792458d5
              if(iRVfilt.eq.1)fluxU1=fluxU1*(1.0d0-beam1*dop1)
              if(iRVfilt.eq.2)fluxB1=fluxB1*(1.0d0-beam1*dop1)
              if(iRVfilt.eq.3)fluxV1=fluxV1*(1.0d0-beam1*dop1)
              if(iRVfilt.eq.4)fluxR1=fluxR1*(1.0d0-beam1*dop1)
              if(iRVfilt.eq.5)fluxI1=fluxI1*(1.0d0-beam1*dop1)
              if(iRVfilt.eq.6)fluxJ1=fluxJ1*(1.0d0-beam1*dop1)
              if(iRVfilt.eq.7)fluxH1=fluxH1*(1.0d0-beam1*dop1)
              if(iRVfilt.eq.8)fluxK1=fluxK1*(1.0d0-beam1*dop1)

              if(iRVfilt.eq.1)fluxU2=fluxU2*(1.0d0-beam2*dop2)
              if(iRVfilt.eq.2)fluxB2=fluxB2*(1.0d0-beam2*dop2)
              if(iRVfilt.eq.3)fluxV2=fluxV2*(1.0d0-beam2*dop2)
              if(iRVfilt.eq.4)fluxR2=fluxR2*(1.0d0-beam2*dop2)
              if(iRVfilt.eq.5)fluxI2=fluxI2*(1.0d0-beam2*dop2)
              if(iRVfilt.eq.6)fluxJ2=fluxJ2*(1.0d0-beam2*dop2)
              if(iRVfilt.eq.7)fluxH2=fluxH2*(1.0d0-beam2*dop2)
              if(iRVfilt.eq.8)fluxK2=fluxK2*(1.0d0-beam2*dop2)

              ymodU(icount)=fluxU1+fluxU2+dfluxU+third(1)
              ymodB(icount)=fluxB1+fluxB2+dfluxB+third(2)
              ymodV(icount)=fluxV1+fluxV2+dfluxV+third(3)
              ymodR(icount)=fluxR1+fluxR2+dfluxR+third(4)
              ymodI(icount)=fluxI1+fluxI2+dfluxI+third(5)
              ymodJ(icount)=fluxJ1+fluxJ2+dfluxJ+third(6)
              ymodH(icount)=fluxH1+fluxH2+dfluxH+third(7)
              ymodK(icount)=fluxK1+fluxK2+dfluxK+third(8)
              yeclipse(icount)=dble(Neclipse1)
c
              fracs1(icount,1)=fluxU1
              fracs2(icount,1)=fluxB1
              fracs3(icount,1)=fluxV1
              fracs4(icount,1)=fluxR1
              fracs5(icount,1)=fluxI1
              fracs6(icount,1)=fluxJ1
              fracs7(icount,1)=fluxH1
              fracs8(icount,1)=fluxK1

              fracs1(icount,2)=fluxU2
              fracs2(icount,2)=fluxB2
              fracs3(icount,2)=fluxV2
              fracs4(icount,2)=fluxR2
              fracs5(icount,2)=fluxI2
              fracs6(icount,2)=fluxJ2
              fracs7(icount,2)=fluxH2
              fracs8(icount,2)=fluxK2

              fracs1(icount,3)=dfluxU
              fracs2(icount,3)=dfluxB
              fracs3(icount,3)=dfluxV
              fracs4(icount,3)=dfluxR
              fracs5(icount,3)=dfluxI
              fracs6(icount,3)=dfluxJ
              fracs7(icount,3)=dfluxH
              fracs8(icount,3)=dfluxK

              ymods1(icount)=flux1
              ymods2(icount)=flux2
              ymods3(icount)=third(iRVfilt)
              ymodd(icount)=dflux
              RV1(icount)=vel1
              dRV1(icount)=delvel1
              RV2(icount)=vel2
              dRV2(icount)=delvel2
c
              call copyinty(ialphmax1,ibetmax1,Nalph1,Nbet1,
     $              Nalph2,Nbet2,rinty1,saveinty1,rinty2,saveinty2,
     %              mmdx1,mmdx2,ibetlim1,ibetlim2,ialphmax2,ibetmax2)
              if(idint.ge.1)call copydiskinty(Nrmax,Nthetamax,Nradius,
     %              Ntheta,dinty,savedinty,einty,saveeinty)

            endif  ! end if iatm.ge.1
C
            if(idraw.eq.1)then

c              write(*,*)'idint = ',idint

              ppp=phase
              dppp=dummyphase
              if((ecc.gt.0.0d0).and.(ism1.eq.0))then  
                 ppp=dmod(em*180.0d0/pie-phase,360.0d0)
                 ppp=pstart
              endif
              if((ecc.gt.0.0d0).and.(ism1.gt.0))then
                if(ipstep.eq.1)ppp=dmod(em*180.0d0/(pie)-phase,360.0d0)
                if(ipstep.eq.2)ppp=dmod(emnew*180.0d0/(pie)-phase,360.0d0)
                ppp=pstart
              endif
              if(ecc.gt.0.0d0)dppp=ppp+180.0d0

c              write(*,*)ppp,phaseout,pstart
              if(idint.ge.1)then
                call hidgrid(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,ppp,
     $           finc,Q,x1,y1,z1,
     $           xend1,visib1,projarray1,g1,gscale1,surf1,
     %           Ndhoriz,dxhoriz,dyhoriz,saveinty1,extension,separ,fluxV1,
     #           reff1,iecheck,temp1,Nhoriz1,xhoriz1,yhoriz1,bdist,mmdx1)
              endif
c
              if(idint.le.0)then
                call hidgrid(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,ppp,
     $           finc,Q,x1,y1,z1,
     $           xend1,visib1,projarray1,g1,gscale1,surf1,
     %           Nhoriz2,xhoriz2,yhoriz2,saveinty1,extension,separ,fluxV1,
     #           reff1,iecheck,temp1,Nhoriz1,xhoriz1,yhoriz1,bdist,mmdx1)
              endif
c
              if(teff2.gt.0.0d0)call hidgrid(2,ialphmax2,ibetmax2,
     &           Nalph2,Nbet2,ibetlim2,dppp,
     %           finc,Q,x2,y2,z2,
     $           xend2,visib2,projarray2,g2,gscale2,surf2,
     %           Nhoriz1,xhoriz1,yhoriz1,saveinty2,extension,separ,fluxV2,
     #           reff2,iecheck,temp2,Nhoriz2,xhoriz2,yhoriz2,bdist,mmdx2)
c
              if(idint.gt.0)call hiddiskgrid(Nrmax,Nthetamax,Nradius,
     %            Ntheta,diskproj,edgeproj,dvisib,evisib,dtemp,tedge,drad,
     $            dx,dy,dz,xedge,yedge,zedge,savedinty,saveeinty,stepr,
     &            stepz,ppp,finc,Q,
     $            Nhoriz1,xhoriz1,yhoriz1,extension,separ,dflux,bdist)
c
              call getcoords(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     &           ppp,finc,Q,
     $           x1,y1,z1,gradx1,grady1,gradz1,temp1,
     %           Ndumsky1,dumxsky1,dumysky1,xend1,
     $           extension,separ,bdist,mmdx1)
c
              if(teff2.gt.0.0d0)call getcoords(2,ialphmax2,ibetmax2,
     &           Nalph2,Nbet2,ibetlim2,dppp,finc,Q,
     $           x2,y2,z2,gradx2,grady2,gradz2,temp2,
     %           Ndumsky2,dumxsky2,dumysky2,xend2,
     $           extension,separ,bdist,mmdx2)
c
              call writepoints(ialphmax1,ibetmax1,Nthetamax,Nrmax,
     $           Nsky1,xsky1,ysky1,Nsky2,xsky2,
     &           ysky2,Nhoriz1,xhoriz1,yhoriz1,Nhoriz2,xhoriz2,yhoriz2,
     &           Nskydisk,xskydisk,yskydisk,zskydisk,
     %           Ntop2,xtop2horiz,ytop2horiz,
     &           Ndtop,dtopx,dtopy,Ndhoriz,dxhoriz,dyhoriz,
     %           Nskyedge,xskyedge,yskyedge,extension,separ,teff2,idint,
     &           ialphmax2,ibetmax2)
            endif
c
c   Next, set various toggle switches based on whether any points were
c   eclipsed.  If iecheck = 1, then the code will not check for eclipses
c   for phases between togglephase and (180-togglephase) and between
c   (180+togglephase) and (360-(togglephase-180)) where togglephase
c   is the phase where the number of eclipsed points is zero for the first
c   time.
c
c   Also, if iecheck was set to 1, the disk integration will be skipped
c   at phases where there is no eclipse.
c
            Ntotal=Neclipse1+Neclipse2+Neclipsed
            if(phase.eq.0.0d0)Ntot0=Ntotal
            if((itoggle.eq.-1).and.(Ntotal.eq.0))then
              itoggle=1
              if(icount.gt.1)then
                if(phase.lt.180.0d0)togglephase=phase        
              else
                if(phase.lt.180.0d0)togglephase=dphase
              endif
              if(iecheck.ge.1)idcheck=-100
            endif
 9939       if(isw7.ge.2)xmod(icount)=timearray(icount)
            if(isw7.le.1)xmod(icount)=phase/360.0d0

c            if(ymodV(icount).le.1.d10)write(*,*)phase,ymodV(icount)
c
c   RVG BUG ALERT
c
c   Change the assignment of the value of xmod for eccentric orbits as
c   below.
c
            if((ecc.gt.0.0d0).and.(ism1.eq.0))xmod(icount)=em/(2.0d0*pie)
            if((ecc.gt.0.0d0).and.(ism1.gt.0))then
              if(ipstep.eq.1)xmod(icount)=em/(2.0d0*pie)
              if(ipstep.eq.2)xmod(icount)=emnew/(2.0d0*pie)
            endif
c
            if((isw7.ge.2).and.(ecc.gt.0.0d0))xmod(icount)=timearray(icount)
            write(64,101)phase,Ntotal,iecheck,idcheck,itoggle,togglephase
            tshift=pshift+eshift
c
c   UPDATE September 10, 2001
c
c   Add the if-then clauses.
c
            if(ikeep.eq.1)tshift=pshift+eshift-pconj
            if(ikeep.eq.2)tshift=pshift+eshift-pconj2
            qqq=dmod(xmod(icount)+tshift,1.0d0)
            overQ=1.0d0/Q
            rrr=pot2/overQ+0.5d0*(overQ-1.0d0)/overQ
c
c   RVG BUG ALERT   May 16, 2001
c
c   Record parameters for eccentric orbits in unit 65.
c
            if(ecc.gt.0.0d0)write(65,1011)qqq,bdist,tpole1,
     #        tpole2,fill1,fill2,rpole1,rpole2,pot1,rrr,SA1,SA2,pots1,pots2,
     #        gscale1,gscale2
c
c   UPDATE OCTOBER 20, 2005
c
c   change to .ge.5
c
            if((iecheck.ge.5).and.(ecc.gt.0.0d0).and.((iskip1.eq.10).
     &          or.(iskip2.eq.10)))go to 999
            if((iecheck.ge.5).and.(ecc.eq.0.0d0).and.((iskip1.eq.10).
     &          and.(iskip2.eq.10)))go to 10
 10       continue
c
 999      continue    ! continue the big loop if eccentric
c
 1011     format(f15.6,1x,f7.5,1x,2(f10.5,1x),1x,2(f7.5,1x),1x,2(f7.5,1x),
     %        2(f9.5,1x),1x,2(f9.6,1x),4(f9.6,1x))

          Nphase=icount
          icountx=icount
c
c
c   If ism1=1, then we have to reflect the light curves to get phases 180
c   to 360-dphase.
c
c   RVG BUG ALERT   April 19, 2001
c
c   If the eccentricity is greater than 0, do not 'finish' the light curve
c   since it is already complete.
c
c   UPDATE October 18, 2002
c
c   Use the input flags sw7 and sw8 to define a phase range
c   to compute.  Require sw7 > 0 and sw8 > 0  AND  sw7 < sw 8
c   If this is true, then don't complete the light curves.   
c
c
          if(isw7.ge.2)go to 1233
          if((sw7.gt.0.0d0).and.(sw8.gt.0.0d0).and.(sw7.lt.sw8))go to 654
c
          if((ionephase.eq.0).and.(ism1.ge.1).and.(ecc.eq.0.0d0))then
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodU,1)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodB,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodV,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodR,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodI,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodJ,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodH,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodK,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymods1,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymods2,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymods3,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,ymodd,0)
            call finishlc(Nmaxphase,icount,dphase,Nphase,xmod,yeclipse,0)
            call finishRV(1,Nmaxphase,icount,dphase,Nphase,xmod,RV1,gamma,
     %          ymods1,dRV1,separ,period,Q,finc)
            call finishRV(2,Nmaxphase,icount,dphase,Nphase,xmod,RV2,gamma,
     %          ymods2,dRV2,separ,period,Q,finc)
          endif            
c
          if(isw24.ge.1)call getfracs(Nmaxphase,icount,
     %         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,
     &         compfracs,dphase,Nphase,
     %         xmod,eshift,pshift,ionephase,ism1,ecc,onephase,sw26)

c   UPDATE October 18, 2002
c
c   Add the statement label 654 below.
c
 654      close(64)
c
c   RVG BUG ALERT   May 16, 2001
c
c   Close the file with the eccentric parameters.
c
          if(ecc.gt.0.0d0)close(65)
c
c   Reset the value of iecheck.
c
          iecheck=iesave
          if(iRVfilt.eq.1)icnU=iVsave
          if(iRVfilt.eq.2)icnB=iVsave
          if(iRVfilt.eq.3)icnV=iVsave
          if(iRVfilt.eq.4)icnR=iVsave
          if(iRVfilt.eq.5)icnI=iVsave
          if(iRVfilt.eq.6)icnJ=iVsave
          if(iRVfilt.eq.7)icnH=iVsave
          if(iRVfilt.eq.8)icnK=iVsave
c
c   Apply the phase shift, if any.
c
c   RVG BUG ALERT   May 4, 2001
c
c   add ionephase.eq.0 to the if() statement
c
c


 1233     if((isw7.ge.2))go to 12345
          if(ionephase.eq.0)then
            if((ecc.gt.0.0d0).or.(pshift.ne.0.0d0))then
              tshift=pshift+eshift
c
c   UPDATE September 10, 2001
c
c   Add the if-then clauses.
c
              if(ikeep.eq.1)tshift=pshift+eshift-pconj
              if(ikeep.eq.2)tshift=pshift+eshift-pconj2
c
c  UPDATE March 26, 2002
c
c  Remove the variables icount and dphase from the argument list of
c  shiftlc (they are not used).
c
              call shiftlc(Nmaxphase,Nphase,xmod,ymodU,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodB,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodV,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodR,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodI,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodJ,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodH,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodK,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymods1,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymods2,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymods3,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodd,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,RV1,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,dRV1,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,RV2,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,yeclipse,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,dRV2,tshift,1)
            endif
          endif
c
c   UPDATE September 10, 2001
c
c   If Teff < 0, then compute the duration of the X-ray eclipse if any.
c   The length of the X-ray eclipse in DEGREES will be stored in the variable
c   obsparm(9).
c
c   UPDATE December 11, 2001
c
c   If ionephase > 0, only one phase was requested (presumably for
c   drawing purposes).  If ionephase > 0, skip to the end (go to 12345)
c

          if(isw23.ge.1)then
            call distorttime(Nmaxphase,Nphase,xmod,ymodV,yeclipse,RV2,gamma,
     %         dphase,pconj)
          endif

          if(ionephase.ge.1)go to 12345
c
c          xecl=0.0d0
          if(Teff2.le.0.0d0)then
c
            if(ecc.eq.0.0d0)then   ! case for circular orbits
              call sort2(icountx,xecx,xecy)
c
c   Check at phase 0.0 to see if there is an eclipse.
c
              bdist=1.0d0
              phase=0.0d0
c
              call gethorizon(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %          phase,finc,Q,
     $          pot1,omega1,x1,y1,z1,rad1,gradx1,grady1,gradz1,xend1,
     $          Nhoriz1,xhoriz1,yhoriz1,phiar1,iedgestar1,delphie1,bdist,
     *          mmdx1,xhmin1,xhmax1,yhmin1,yhmax1,tidephi,itide,phihor1)
c
              call getXecl(Nhoriz1,xhoriz1,yhoriz1,ixecl,Q,finc,bdist,phase)
c
              if(ixecl.le.0)then
                obsparm(9)=0.0d0
                write(2,321)obsparm(9)
                go to 12345  ! no eclipse here
              endif
c
c   Now we have to locate the phase where the X-ray source is not
c   eclipsed and iterate.
c
              phin=0.0
              do 12000 ii=2,icountx
                if(xecy(ii).gt.10.0d0)then
                  phin=xecx(ii)
                  go to 12000
                else
                  phout=xecx(ii)
                  go to 12001
                endif
12000         continue
c
12001         do 12002 ii=1,30
                phase=0.5d0*(phout+phin)
                bdist=1.0d0
c
                call gethorizon(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %            phase,finc,Q,
     $            pot1,omega1,x1,y1,z1,rad1,gradx1,grady1,gradz1,xend1,
     $            Nhoriz1,xhoriz1,yhoriz1,phiar1,iedgestar1,delphie1,bdist,
     &            mmdx1,xhmin1,xhmax1,yhmin1,yhmax1,tidephi,itide,phihor1)
c
                call getXecl(Nhoriz1,xhoriz1,yhoriz1,ixecl,Q,finc,bdist,phase)
c
                if(ixecl.gt.10)phin=phase
                if(ixecl.lt.-10)phout=phase
12002         continue
c
              obsparm(9)=2.0d0*phin
              write(2,321)obsparm(9)
c
            endif       ! endif ecc = 0
c
            if(ecc.gt.0.0d0)then   ! case for eccentric
              call sort2(icountx,xecx,xecy)
c
c   Check at the conjunction phase to see if there is an eclipse.
c
              rnu=-argper-90.0d0
              rnu=rnu*pie/180.0d0
              phase=0.0d0
              bdist=(1.0d0-ecc*ecc)/(1.0d0+ecc*dcos(rnu))
c
              call setupgeo(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $         fill1,omega1,Q,finc,
     $         x1,y1,z1,
     %         surf1,rad1,gradx1,grady1,gradz1,g1,xend1,separ,
     %         Tgrav1,Teff1,reff1,Rl1,Tpole1,Rpol1,Regg1,SA1,pot1,gpole1,
     %         phiar1,isquare,iusepot,usepot1,ivrt,pervol1,
     &         fillper1,bdist,pots1,1,mmdx1,primmass,primK,
     &         primrad,ratrad,frac1,frac2,ecc,period,size1,fill2,
     &         omegatwo,sw5,tteff2,density,tidephi,itide,phistart1)
c 
              call gethorizon(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %          phase,finc,Q,
     $          pot1,omega1,x1,y1,z1,rad1,gradx1,grady1,gradz1,xend1,
     $          Nhoriz1,xhoriz1,yhoriz1,phiar1,iedgestar1,delphie1,bdist,
     &          mmdx1,xhmin1,xhmax1,yhmin1,yhmax1,tidephi,itide,phihor1)
c
              call getXecl(Nhoriz1,xhoriz1,yhoriz1,ixecl,Q,finc,bdist,phase)
c
              if(ixecl.le.0)then
                obsparm(9)=0.0d0
                write(2,321)obsparm(9)
                go to 12345  ! no eclipse here
              endif
c
c   Now we have to locate the phase where the X-ray source is not
c   eclipsed and iterate.
c
              phin=0.0
c
              do 22000 ii=1,icountx
                if(xecy(ii).gt.10.0d0)then
                  phin=xecx(ii)
                  go to 22000
                else
                  phout=xecx(ii)
                  go to 22001
                endif
22000         continue
c
22001         do 22002 ii=1,20
                phase=0.5d0*(phout+phin)
                rnu=pie/180.0d0*dmod((phase-90.0d0-argper),360.0d0)
                if(rnu.lt.0.0d0)rnu=rnu+2.0*pie
                bdist=(1.0d0-ecc*ecc)/(1.0d0+ecc*dcos(rnu))
c
                call setupgeo(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $            fill1,omega1,Q,finc,
     $            x1,y1,z1,
     %            surf1,rad1,gradx1,grady1,gradz1,g1,xend1,separ,
     %            Tgrav1,Teff1,reff1,Rl1,Tpole1,Rpol1,Regg1,SA1,pot1,gpole1,
     %            phiar1,isquare,iusepot,usepot1,ivrt,pervol1,
     &            fillper1,bdist,pots1,1,mmdx1,primmass,primK,
     &            primrad,ratrad,frac1,frac2,ecc,period,size1,fill2,omegatwo,
     &            sw5,tteff2,density,tidephi,itide,phistart1)
c 
                call gethorizon(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %            phase,finc,Q,
     $            pot1,omega1,x1,y1,z1,rad1,gradx1,grady1,gradz1,xend1,
     $            Nhoriz1,xhoriz1,yhoriz1,phiar1,iedgestar1,delphie1,bdist,
     &            mmdx1,xhmin1,xhmax1,yhmin1,yhmax1,tidephi,itide,phihor1)
c
                call getXecl(Nhoriz1,xhoriz1,yhoriz1,ixecl,Q,finc,bdist,phase)
c
                if(ixecl.gt.10)phin=phase
                if(ixecl.lt.-10)phout=phase
22002         continue
c
c   Now figure out the observed phase of the angle corresponding to rnu.
c
              if(rnu.lt.0.0d0)rnu=rnu+2.0d0*pie
              ECAN=2.0d0*datan(dtan(rnu*0.5d0)*
     #           dsqrt((1.0d0-ecc)/(1.0d0+ecc)))
              ECAN=dmod(ECAN,2.0d0*pie)
              if(ECAN.lt.0.0d0)ECAN=ECAN+2.0d0*pie
              obsph1=ECAN-ecc*dsin(ECAN)            ! in radians!
c
c   Now look at the other ingress phase.
c
              do 42000 ii=icountx,1,-1
                if(xecy(ii).gt.10.0d0)then
                  phin=xecx(ii)
                  go to 42000
                else
                  phout=xecx(ii)
                  if(ii.eq.icountx)phin=359.99d0
                  go to 42001
                endif
42000         continue
c
42001         do 42002 ii=1,20
                phase=0.5d0*(phout+phin)
                rnu=pie/180.0d0*dmod((phase-90.0d0-argper),360.0d0)
                if(rnu.lt.0.0d0)rnu=rnu+2.0*pie
                bdist=(1.0d0-ecc*ecc)/(1.0d0+ecc*dcos(rnu))
c
                call setupgeo(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     $            fill1,omega1,Q,finc,
     $            x1,y1,z1,
     %            surf1,rad1,gradx1,grady1,gradz1,g1,xend1,separ,
     %            Tgrav1,Teff1,reff1,Rl1,Tpole1,Rpol1,Regg1,SA1,pot1,gpole1,
     %            phiar1,isquare,iusepot,usepot1,ivrt,pervol1,
     &            fillper1,bdist,pots1,1,mmdx1,primmass,primK,
     $            primrad,ratrad,frac1,frac2,ecc,period,size1,fill2,
     $            omegatwo,sw5,tteff2,density,tidephi,itide,phistart1)
c 
                call gethorizon(1,ialphmax1,ibetmax1,Nalph1,Nbet1,ibetlim1,
     %            phase,finc,Q,
     $            pot1,omega1,x1,y1,z1,rad1,gradx1,grady1,gradz1,xend1,
     $            Nhoriz1,xhoriz1,yhoriz1,phiar1,iedgestar1,delphie1,bdist,
     *            mmdx1,xhmin1,xhmax1,yhmin1,yhmax1,tidephi,itide,phihor1)
c
                call getXecl(Nhoriz1,xhoriz1,yhoriz1,ixecl,Q,finc,bdist,phase)
c
                if(ixecl.gt.10)phin=phase
                if(ixecl.lt.-10)phout=phase
42002         continue
c
c   Figure out the observed phase of the angle rnu.
c
              if(rnu.lt.0.0d0)rnu=rnu+2.0d0*pie
              ECAN=2.0d0*datan(dtan(rnu*0.5d0)*
     #           dsqrt((1.0d0-ecc)/(1.0d0+ecc)))
              ECAN=dmod(ECAN,2.0d0*pie)
              if(ECAN.lt.0.0d0)ECAN=ECAN+2.0d0*pie
              obsph2=ECAN-ecc*dsin(ECAN)    ! in radians
c
              ddd2=dabs(obsph2-obsph1)
              if(ddd2.gt.pie)then
                ddd2=(obsph1+2.0d0*pie)-obsph2
              endif
              obsparm(9)=ddd2*180.0d0/pie
              write(2,321)obsparm(9)
c
            endif       ! endif ecc = 0
c
          endif     ! endif Teff < 0

 100      format(1x,f9.7,3x,f9.5,3x,f8.6)
 101      format(f15.8,10x,i8,2x,i4,1x,i4,2x,i2,1x,f13.8)
 102      format(/'Info:  Model atmospheres will be used',
     &      /'Number of models in the table = ',i4,
     %      /'Maximum tabulated temperature = ',f8.1,2x,
     %      /'Minimum tabulated temperature = ',f8.1)
c
c   UPDATE September 11, 2001
c
c   Add this format statement.
c
 321      format('X-ray eclipse duration = ',f11.6,' degrees')
c
c          call getnorm(Nphase,xmod,ymodV,Vnorm)
c          call lcnorm(Nphase,ymodV,Vnorm)
c          call lcnorm(Nphase,ymods1,Vnorm)
c          call lcnorm(Nphase,ymods2,Vnorm)
c
c   RVG BUG ALERT   May 16, 2001
c
c   Close the unit 2 file here (ELC.out) instead of in the main program.
c
12345     close(2)
c
c
c   UPDATE June 7, 2002
c
c   Add the X-ray eclipse duration to the end of parmstring.
c
c   UPDATE October 28, 2002
c
c   Make the length of parmstring character*249, hence add 12 to the 227 below.
c
c   UPDATE October 22, 2008
c
c   Make the length of parmstring character*259.
c
c
          scr1string=parmstring(1:249)
          write(scr2string,54321)obsparm(9)
          parmstring=scr1string//scr2string
c
54321     format(1x,f9.5)
c
c   UPDATE OCTOBER 24, 2005
c
c   If the iecheck=5 option was used, remove the skipped phases.  These
c   will have negative values in the light curves.
c
c  UPDATE NOVEMBER 22, 2006
c
c  Don't clip the RV curves.  Add the variable NRVphase and xRVmod
c
          NRVphase=Nphase
          do 12346 ii=1,Nphase
            xRVmod(ii)=xmod(ii)
12346     continue
c
          if(iecheck.eq.5)then
            icount=Nphase
            call cliplc(Nmaxphase,icount,dphase,Nphase,xmod,ymodU,
     %         ymodB,ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,
     &         ymods2,ymods3,ymodd,RV1,dRV1,RV2,dRV2)
          endif            
c
c
c   UPDATE JULY 21, 2006
c
c   Add a "fast genetic" mode.  If ifastflag=1, then set 
c     Nalph1=40 
c     Nbet1=14
c     Nalph2=40
c     Nbet2=14
c     dphase=3
c 
c   We need to restore the reset values above.
c
          Nalph1=iNa1save
          Nalph2=iNa2save
          Nbet1=iNb1save
          Nbet2=iNb2save
          dphase=savedphase
c
c   UPDATE JULY 13, 2009
c
c   If requested, bin the light and velocity curves.  sw29 will be
c   the binsize for the photometry, in minutes, and sw30 will be the bin
c   size for the RV curves, in minutes
c 
          if((sw29.gt.0.0d0).and.(icnU.ne.430))call 
     &       binlc(Nmaxphase,Nphase,xmod,ymodU,period,
     #       dphase,sw29)
          if((sw29.gt.0.0d0).and.(icnB.ne.430))call 
     &       binlc(Nmaxphase,Nphase,xmod,ymodB,period,
     #       dphase,sw29)
          if((sw29.gt.0.0d0).and.(icnV.ne.430))call 
     &       binlc(Nmaxphase,Nphase,xmod,ymodV,period,
     #       dphase,sw29)
          if((sw29.gt.0.0d0).and.(icnR.ne.430))call 
     &       binlc(Nmaxphase,Nphase,xmod,ymodR,period,
     #       dphase,sw29)
          if((sw29.gt.0.0d0).and.(icnI.ne.430))call 
     &       binlc(Nmaxphase,Nphase,xmod,ymodI,period,
     #       dphase,sw29)
          if((sw29.gt.0.0d0).and.(icnJ.ne.430))call 
     &       binlc(Nmaxphase,Nphase,xmod,ymodJ,period,
     #       dphase,sw29)
          if((sw29.gt.0.0d0).and.(icnH.ne.430))call 
     &       binlc(Nmaxphase,Nphase,xmod,ymodH,period,
     #       dphase,sw29)
          if((sw29.gt.0.0d0).and.(icnK.ne.430))call 
     &       binlc(Nmaxphase,Nphase,xmod,ymodK,period,
     #       dphase,sw29)
          if(sw29.gt.0.0d0)call binlc(Nmaxphase,Nphase,xmod,ymods1,period,
     #       dphase,sw29)
          if(sw29.gt.0.0d0)call binlc(Nmaxphase,Nphase,xmod,ymods2,period,
     #       dphase,sw29)
          if(sw29.gt.0.0d0)call binlc(Nmaxphase,Nphase,xmod,ymods3,period,
     #       dphase,sw29)
          if((sw29.gt.0.0d0).and.(idint.gt.0))
     &       call binlc(Nmaxphase,Nphase,xmod,ymodd,period,
     #       dphase,sw29)
          if(sw30.gt.0.0d0)call binlc(Nmaxphase,NRVphase,xRVmod,RV1,period,
     #       dphase,sw30)
          if(sw30.gt.0.0d0)call binlc(Nmaxphase,NRVphase,xRVmod,RV2,period,
     #       dphase,sw30)
          if(sw30.gt.0.0d0)call binlc(Nmaxphase,NRVphase,xRVmod,dRV1,period,
     #       dphase,sw30)
          if(sw30.gt.0.0d0)call binlc(Nmaxphase,NRVphase,xRVmod,dRV2,period,
     #       dphase,sw30)
c
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
           subroutine binlc(Nmaxphase,Nphase,xmod,ymod,period,dphase,sw29)
c
c   This subroutine will bin a light curve in phase, using a binsize given
c   by sw29, where the units are minutes.
c
           implicit double precision (a-h,o-z)
c
           dimension xmod(Nmaxphase),ymod(Nmaxphase)
           dimension yinter(900000),y2(900000)
           dimension xpad(900000),ypad(900000)
c
           call addpad(Nphase,xmod,ymod,xpad,ypad)
           Mphase=Nphase*3
c
c   We will interpolate the model, and given each value in xmod, figure
c   out the range of phase needed, and average the interpolated y-values
c
           call spline(xpad,ypad,Mphase,0.0d0,0.0d0,y2)
c
           pstep=sw29/1440.0d0/period
           pwidth=0.5d0*pstep
c
           NN=7
           xxxhigh=1.0d0            !xmod(Nphase)
           do 10 i=1,Mphase
             if(xpad(i).gt.1.0d0)go to 10
             if(xpad(i).lt.0.0d0)go to 10
             kcount=0
             summ=0.0d0        
             jlo=i
c 
             a=xpad(i)-pwidth
             b=xpad(i)+pwidth

             do 9 j=1,NN
               if(j.eq.1)then
                 xxx=a
c                 if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c                 if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
                 call hunt(xpad,Mphase,xxx,jlo)
                 if((jlo.eq.Mphase).or.(jlo.eq.0))then
                   call splint(xpad,ypad,y2,Mphase,xxx,qqqa)
                 else
                   call fastsplint(xpad,ypad,y2,Mphase,xxx,qqqa,jlo,jlo+1)
                 endif
                 xxx=b
c                 if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c                 if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
                 call hunt(xpad,Mphase,xxx,jlo)
                 if((jlo.eq.Mphase).or.(jlo.eq.0))then
                   call splint(xpad,ypad,y2,Mphase,xxx,qqqb)
                 else
                   call fastsplint(xpad,ypad,y2,Mphase,xxx,qqqb,jlo,jlo+1)
                 endif
                 summ=0.5d0*(b-a)*(qqqa+qqqb)
               else
                 it=2**(j-2)
                 tnm=dble(it)
                 del=(b-a)/tnm
                 x=a+0.5d0*del
                 xxx=x
c                 if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c                 if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
                 s=0.0d0
                 do 11 kk=1,it
                   kcount=kcount+1 
                   call hunt(xpad,Mphase,xxx,jlo)
                   if((jlo.eq.Mphase).or.(jlo.eq.0))then
                     call splint(xpad,ypad,y2,Mphase,xxx,qqq)
                   else
                     call fastsplint(xpad,ypad,y2,Mphase,xxx,qqq,jlo,jlo+1)
                   endif
                   s=s+qqq
                   x=x+del
                   xxx=x
c                   if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c                   if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
11               continue
                 summ=0.5d0*(summ+(b-a)*s/tnm)
               endif
c
9           continue

           yinter(i)=summ/pstep

10         continue
c


c           do 9 j=-25,25
c             xxx=dble(j)*0.04d0*pwidth+xmod(i)
c             if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c             if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
c             call hunt(xmod,Nphase,xxx,jlo)
c             if((jlo.eq.Nphase).or.(jlo.eq.0))then
c               call splint(xmod,ymod,y2,Nphase,xxx,qqq)
c             else
c               call fastsplint(xmod,ymod,y2,Nphase,xxx,qqq,jlo,jlo+1)
c             endif
c             summ=summ+qqq
c9          continue
c           yinter(i)=summ/51.0d0
c10         continue
cc

           do 20 i=1,Mphase
             ypad(i)=yinter(i)
20         continue
c
           call removepad(Nphase,xmod,ymod,xpad,ypad)

           return
           end
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine getnorm(N,xx,yy,ynorm)
c
c
c         
          implicit double precision(a-h,o-z)

          dimension yy(N),xx(N)
c
          diffmin=12345.d0
          do 10 i=1,N
            diff=dabs(xx(i)-0.2500d0)
            if(diff.lt.diffmin)then
              ymax=yy(i)
              diffmin=diff
            endif
 10       continue
 15       ynorm=ymax
c
          return
          end
c
c   *******************************************
c
          subroutine lcnorm(N,yy,ynorm)
c
c
c
          implicit double precision(a-h,o-z)
c
          dimension yy(N)
c
          do 10 i=1,N
            yy(i)=yy(i)/ynorm
 10       continue
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine finishlc(Nmaxphase,icount,dphase,Nphase,xxx,yyy,iwrite)
c
c   December 3, 1999
c
c   This routine will use the symmetry of the light curve and fill out
c   the phases 180 to 360-dphase.
c   
          implicit double precision(a-h,o-z)
c
          dimension xxx(Nmaxphase),yyy(Nmaxphase)
c
          i=icount-1
          Nphase=icount
          do 10 phase=180.0d0+dphase,360.0d0-dphase,dphase
            if(iwrite.eq.1)write(64,101)phase
            Nphase=Nphase+1
            xxx(Nphase)=phase/360.0d0
            yyy(Nphase)=yyy(i)
            i=i-1
 10       continue
c
 101      format(f6.2)
          return
          end
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine finishRV(istar,Nmaxphase,icount,dphase,
     &      Nphase,xxx,yyy,gamma,ymods,dRV,separ,period,Q,finc)
c
c   December 3, 1999
c
c   This routine will use the symmetry of the velocity curve and fill out
c   the phases 180 to 360-dphase.
c   
          implicit double precision(a-h,o-z)

          parameter(pie=3.14159265358979323d0)
          dimension xxx(Nmaxphase),yyy(Nmaxphase),ymods(Nmaxphase),
     &                   dRV(Nmaxphase)
c
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
          a=separ*6.9598d5            !separation in km
          p=period*24.00d0*3600.00d0       !period in seconds
          velamp=2.0d0*pie*a/(p)
          fincr=(finc/180.0d0)*pie
          sifinc=dsin(fincr)
c
          i=icount-1
          Nphase=icount
          do 10 phase=180.0d0+dphase,360.0d0-dphase,dphase
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
            siphase=dsin(phaser)
            Nphase=Nphase+1
            xxx(Nphase)=phase/360.0d0
            if(istar.eq.1)then
              dRV(Nphase)=-dRV(i)
            else
              dRV(Nphase)=-dRV(i)       ! was +
            endif
c
c   RVG BUG ALERT  April 27, 2001
c
c   Comment out this block below and replace it with the simple statement
c   below
c
c              vel=overQ/(1.0d0+overQ)*(siphase*sifinc)+dRV(Nphase)
c              vel=velamp*vel
c            if(istar.eq.1)then
c               yyy(Nphase)=vel+gamma
c            else
c              yyy(Nphase)=-vel+gamma
c            endif
c
            yyy(Nphase)=-(yyy(i)-gamma)+gamma     
c
c   END BUG
c
            i=i-1
 10       continue
c
 101      format(f6.2)
          return
          end
c
c  ********************************
c
c  RVG BUG ALERT  May 9, 2001
c
c  Add the spot parameters to the argument list.
c
c
c  UPDATE August 10, 2004
c
c  Add 8 real variables and 4 integer variables to the argument list.
c
c  UPDATE May 8, 2006
c
c  add isw21-isw24, sw21-sw24, powercoeff to the list
c
c  UPDATE November 6, 2008
c
c  Add sw25-sw34 and isw25-isw34 below
c
          subroutine getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     &       Ntheta,Nradius,alb1,alb2,Nref,
     %       rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
     $       idraw,iecheck,idint,iatm,ism1,
     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,
     &       isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey,
     $       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     $       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
     %       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       idark1,idark2,isw12,isw13,
     %       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2,
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
          implicit double precision(a-h,o-z)
c
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
c
c   RVG BUG ALERT  May 9, 2001
c
c   Dimension the spot arrays.
c
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension powercoeff(8,9)
c
c   UPDATE June 22, 2002
c
c   Declare the variable bell to be character*1
c
          character*1 bell
c
          ios=0
          open(unit=1,file='ELC.inp',status='old',err=100,iostat=ios)
c
          read(1,*)Nalph1
          read(1,*)Nbet1
          read(1,*)Nalph2
          read(1,*)Nbet2
          read(1,*)fill1
          read(1,*)fill2
          read(1,*)omega1
          read(1,*)omega2
          read(1,*)dphase
          read(1,*)Q
          read(1,*)finc
          read(1,*)Teff1
          read(1,*)Teff2
          read(1,*)Tgrav1
          read(1,*)Tgrav2
          read(1,*)betarim
          read(1,*)rinner
          read(1,*)router
          read(1,*)tdisk
          read(1,*)xi
          read(1,*)Ntheta
          read(1,*)Nradius
          read(1,*)alb1
          read(1,*)alb2
          read(1,*)Nref
          read(1,*)rLx
          read(1,*)Period
          read(1,*)fm
          read(1,*)separ
          read(1,*)gamma
          read(1,*)t3
          read(1,*)g3
          read(1,*)SA3
          read(1,*)density
          read(1,*)sw1
          read(1,*)sw2
          read(1,*)sw3
          read(1,*)T0
          read(1,*)idraw
          read(1,*)iecheck
          read(1,*)idint
          read(1,*)iatm
          read(1,*)ism1
          read(1,*)icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK
          read(1,*)iRVfilt
          read(1,*)isw1
          read(1,*)isw2
          read(1,*)isw3
          read(1,*)isw4
          read(1,*)ilaw
c
c  Load the limb darkening parameters. 
c         
          do 10 i=1,8
            read(1,*)wave(i),dbolx(i,1),dboly(i,1),dbolx(i,2),dboly(i,2),
     %         dwavex(i,1),dwavey(i,1),dwavex(i,2),dwavey(i,2)
 10       continue
c
c  Look for additional parameters related to the eccentricity.  Use
c  the end=...  option to allow for older input files.
c
          read(1,*,end=99)ecc
          read(1,*)argper
          read(1,*)pshift
          read(1,*)sw5
          read(1,*)sw6
          read(1,*)sw7
          read(1,*)sw8
c
c   UPDATE April 8, 2002
c
c   Correct bug:  change sw5 to sw9 below
c
          read(1,*)sw9
          read(1,*)ikeep
          read(1,*)isynch
          read(1,*)isw5
          read(1,*)isw6
          read(1,*)isw7
          read(1,*)isw8
          read(1,*)isw9
c
c   May 8, 2001
c
c   Now load up parameters describing star (or disk)
c   spots.  The parameters are as follows:
c
c   spot1parm(1,1):  temperature factor for spot 1 on star 1
c   spot1parm(1,2):  latitude of spot 1 on star 1 in DEGREES
c   spot1parm(1,3):  longitude of spot 1 on star 1 in DEGREES
c   spot1parm(1,4):  angular radius of spot 1 on star 1 in DEGREES
c
c   spot1parm(2,1):  temperature factor for spot 2 on star 1
c   spot1parm(2,2):  latitude of spot 2 on star 1 in DEGREES
c   spot1parm(2,3):  longitude of spot 2 on star 1 in DEGREES
c   spot1parm(2,4):  angular radius of spot 2 on star 1 in DEGREES
c
c
c   spot2parm(1,1):  temperature factor for spot 1 on star 2
c   spot2parm(1,2):  latitude of spot 1 on star 2 in DEGREES
c   spot2parm(1,3):  longitude of spot 1 on star 2 in DEGREES
c   spot2parm(1,4):  angular radius of spot 1 on star 2 in DEGREES
c
c   spot2parm(2,1):  temperature factor for spot 2 on star 2
c   spot2parm(2,2):  latitude of spot 2 on star 2 in DEGREES
c   spot2parm(2,3):  longitude of spot 2 on star 2 in DEGREES
c   spot2parm(2,4):  angular radius of spot 2 on star 2 in DEGREES
c
c
c   spotdparm(1,1):  temperature factor for spot 1 on disk  
c   spotdparm(1,2):  azimuth of spot 1 on disk in DEGREES
c   spotdparm(1,3):  radial cutoff for spot 1 on disk (between 0 and 1)
c   spotdparm(1,4):  angular size of spot 1 on disk in DEGREES
c
c   spotdparm(2,1):  temperature factor for spot 2 on disk  
c   spotdparm(2,2):  azimuth of spot 2 on disk in DEGREES
c   spotdparm(2,3):  radial cutoff for spot 2 on disk (between 0 and 1)
c   spotdparm(2,4):  angular size of spot 2 on disk in DEGREES
c
c   If there is an error, abort, setting ispot=0
c
          ios=0
          read(1,*,end=101,err=101)spot1parm(1,1)
          read(1,*,end=101,err=101)spot1parm(1,2)
          read(1,*,end=101,err=101)spot1parm(1,3)
          read(1,*,end=101,err=101)spot1parm(1,4)
c              
          read(1,*,end=101,err=101)spot1parm(2,1)
          read(1,*,end=101,err=101)spot1parm(2,2)
          read(1,*,end=101,err=101)spot1parm(2,3)
          read(1,*,end=101,err=101)spot1parm(2,4)
c              
c              
          read(1,*,end=101,err=101)spot2parm(1,1)
          read(1,*,end=101,err=101)spot2parm(1,2)
          read(1,*,end=101,err=101)spot2parm(1,3)
          read(1,*,end=101,err=101)spot2parm(1,4)
               
          read(1,*,end=101,err=101)spot2parm(2,1)
          read(1,*,end=101,err=101)spot2parm(2,2)
          read(1,*,end=101,err=101)spot2parm(2,3)
          read(1,*,end=101,err=101)spot2parm(2,4)
c              
c              
          read(1,*,end=101,err=101)spotdparm(1,1)
          read(1,*,end=101,err=101)spotdparm(1,2)
          read(1,*,end=101,err=101)spotdparm(1,3)
          read(1,*,end=101,err=101)spotdparm(1,4)
               
          read(1,*,end=101,err=101)spotdparm(2,1)
          read(1,*,end=101,err=101)spotdparm(2,2)
          read(1,*,end=101,err=101)spotdparm(2,3)
          read(1,*,end=101,err=101)spotdparm(2,4)
c
c   UPDATE August 10, 2004
c
c   Add the new variables here.  Abort if there is an error
c
          read(1,*,end=101,err=101)primmass
          read(1,*,end=101,err=101)primK
          read(1,*,end=101,err=101)primrad
          read(1,*,end=101,err=101)ratrad
          read(1,*,end=101,err=101)frac1
          read(1,*,end=101,err=101)frac2
          read(1,*,end=101,err=101)ecosw
          read(1,*,end=101,err=101)temprat
          read(1,*,end=101,err=101)idark1
          read(1,*,end=101,err=101)idark2
          read(1,*,end=101,err=101)isw12
          read(1,*,end=101,err=101)isw13
c
c   UPDATE May 8, 2006
c
c   Add new variables here.  Abort if there is an error.
c
          read(1,*,end=101,err=101)isw21
          read(1,*,end=101,err=101)isw22
          read(1,*,end=101,err=101)isw23
          read(1,*,end=101,err=101)isw24
c
c   Power-series limb darkening coefficients for star 1
c
          read(1,*,end=101,err=101)(powercoeff(1,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(2,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(3,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(4,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(5,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(6,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(7,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(8,k),k=1,9)
c
          read(1,*,end=101,err=101)bigI
          read(1,*,end=101,err=101)bigbeta
          read(1,*,end=101,err=101)sw23
          read(1,*,end=101,err=101)sw24
c
c   UPDATE November 6, 2008
c
c   Add sw25-sw34 and isw25-isw34 below
c
          read(1,*,end=101,err=101)sw25
          read(1,*,end=101,err=101)sw26
          read(1,*,end=101,err=101)sw27
          read(1,*,end=101,err=101)sw28
          read(1,*,end=101,err=101)sw29
          read(1,*,end=101,err=101)sw30
          read(1,*,end=101,err=101)sw31
          read(1,*,end=101,err=101)Tconj
          read(1,*,end=101,err=101)beam1
          read(1,*,end=101,err=101)beam2

          read(1,*,end=101,err=101)isw25
          read(1,*,end=101,err=101)isw26
          read(1,*,end=101,err=101)isw27
          read(1,*,end=101,err=101)isw28
          read(1,*,end=101,err=101)isw29
          read(1,*,end=101,err=101)isw30
          read(1,*,end=101,err=101)isw31
          read(1,*,end=101,err=101)isw32
          read(1,*,end=101,err=101)isw33
          read(1,*,end=101,err=101)isw34

 99       close(1)
c
c   Come here if the input file ELC.inp does not exist.  The subroutine
c   writeinput will make the correct file and set default values.
c
c   RVG BUG ALERT  May 9, 2001
c
c   Update the writeinput routine to include the spot parameters
c
c
c   UPDATE August 10, 2004
c
c   Add 8 real variables and 4 integer variables to the argumemt list
c
 100      if(ios.gt.0)call writeinput(Nalph1,Nbet1,Nalph2,Nbet2,
     $       fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     $       Ntheta,Nradius,
     $       alb1,alb2,Nref,rLx,
     $       Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
     $       idraw,iecheck,idint,iatm,ism1,
     &       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     $       iRVfilt,isw1,isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey,
     $       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     $       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
     $       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       idark1,idark2,isw12,isw13,
     &       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
c   RVG BUG ALERT  May 9, 2001
c
c   Put this if-then block for successful completion.
c
          if(ios.eq.0)then 
            close(1)
            return
          endif

 101      ispot=0               !file ended too soon
          bell=char(7)
          write(*,1002)bell
c
 1000     format(a1,'Error:  File ELC.spot does not exist')
 1001     format(a40)
 1002     format(a1,'Error:  Bad entry in ELC.inp')
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
c   RVG BUG ALERT  May 9, 2001
c
c   Update the writeinput routine to include the spot parameters
c
c   UPDATE August 10, 2004
c
c   Add the 8 real variables and 4 integers to the argument list.
c    
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to argument list
c
c   UPDATE November 6, 2008
c
c   Add sw25-sw34 and isw25-isw34 to the list
c
          subroutine writeinput(Nalph1,Nbet1,Nalph2,Nbet2,
     $       fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     $       Ntheta,Nradius,alb1,alb2,Nref,
     $       rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
     $       idraw,iecheck,idint,iatm,ism1,
     &       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %       iRVfilt,isw1,isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey,
     $       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     $       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
     $       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       idark1,idark2,isw12,isw13,
     %       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
c
c    will write the correctly formatted file ELC.inp and return
c    default parameters
c
          implicit double precision(a-h,o-z)
c
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2),
     %       www(8)
          dimension powercoeff(8,9)
          character*1 bell
c
c   RVG BUG ALERT  May 9, 2001
c
c   Dimension the spot arrays.
c
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
c
          data www/3600.0,4500.0,5550.0,6700.0,8700.0,12000.0,16200.0,22000.0/
c
          bell=char(7)
          write(*,100)bell
c
          Nalph1=40
          Nbet1=14
          Nalph2=40
          Nbet2=14
          fill1=1.00d0
          fill2=0.005d0
          omega1=1.0d0
          omega2=1.0d0
          dphase=3.0d0
          Q=2.0d0
          finc=80.0d0
          Teff1=6500.0d0
          Teff2=6500.0d0
          Tgrav1=0.25d0
          Tgrav2=0.25d0
          betarim=2.0d0
          rinner=0.005d0
          router=0.75d0
          tdisk=30000.0d0
          xi=-0.75d0
          Ntheta=90
          Nradius=60
          alb1=1.0d0
          alb2=1.0d0
          Nref=1
          rLx=0.001d0
          Period=2.62d0
          fm=3.0d0
          separ=5.0d0
          gamma=50.0d0
          idraw=0
          iecheck=1
          idint=1
          iatm=0
          ism1=1
          isw2=0
          ilaw=1
          icnU=1
          icnB=1
          icnV=1
          icnR=1
          icnI=1
          icnJ=1
          icnH=1
          icnK=1
          t3=-5000.0d0
          g3=-5.0d0
          SA3=-0.1d0
          density=0.0d0
          sw1=0.0d0
          sw2=0.0d0
          sw3=0.0d0
          T0=0.0d0
          iRVfilt=3
          isw1=0
          isw2=0
          isw3=0
          isw4=0
          isw5=0
          isw6=0
          isw7=0
          isw8=0
          isw9=0
          ikeep=0
          isynch=0
          ecc=0.0d0
          argper=90.0d0
          pshift=0.0d0
          sw5=0.0d0
          sw6=0.0d0
          sw7=0.0d0
          sw8=0.0d0
          sw9=0.0d0
c
          do 5 i=1,8
            wave(i)=www(i)
            dbolx(i,1)=0.635d0
            dbolx(i,2)=0.635d0
            dboly(i,1)=0.242d0
            dboly(i,2)=0.242d0
 5        continue
c
          do 6 i=1,4
            do 7 j=1,2
              spot1parm(j,i)=-1.0d0
              spot2parm(j,i)=-1.0d0
              spotdparm(j,i)=-1.0d0
 7          continue
 6        continue
c
c   UPDATE August 10, 2004
c 
c   Initialize the new variables to zero.
c
          primmass=0.0d0
          primK=0.0d0
          primrad=0.0d0
          ratrad=0.0d0
c
          frac1=0.0d0
          frac2=0.0d0
          ecosw=0.0d0
          temprat=0.0d0
c
          idark1=0
          idark2=0
          isw12=0
          isw13=0
c
c  UPDATE May 8, 2006
c
c  Add defaults for new variables here.
c
          isw21=0
          isw22=0
          isw23=0
          isw24=0
          bigI=0.0d0
          bigbeta=0.0d0
          sw23=0.0d0
          sw24=0.0d0
c
c    UPDATE November 6, 2008
c
c    Add values for sw25-sw34 and isw25-isw34 below
c
          sw25=0.0d0
          sw26=0.0d0
          sw27=0.0d0
          sw28=0.0d0
          sw29=0.0d0
          sw30=0.0d0
          sw31=0.0d0
          Tconj=0.0d0
          beam1=0.0d0
          beam2=0.0d0

          isw25=0
          isw26=0
          isw27=0
          isw28=0
          isw29=0
          isw30=0
          isw31=0
          isw32=0
          isw33=0
          isw34=0

c
          open(unit=1,file='ELC.inp',status='unknown')
c
          write(1,1000)Nalph1
          write(1,1001)Nbet1
          write(1,1002)Nalph2
          write(1,1003)Nbet2
          write(1,1004)fill1
          write(1,1005)fill2
          write(1,1006)omega1
          write(1,1007)omega2
          write(1,1008)dphase
          write(1,1009)Q
          write(1,1010)finc
          write(1,1011)Teff1
          write(1,1012)Teff2
          write(1,1013)Tgrav1
          write(1,1014)Tgrav2
          write(1,1015)betarim
          write(1,1016)rinner
          write(1,1017)router
          write(1,1018)tdisk
          write(1,2018)xi
          write(1,1019)Ntheta
          write(1,1020)Nradius
          write(1,1042)alb1
          write(1,1043)alb2
          write(1,2043)Nref
          write(1,1021)rLx
          write(1,1023)Period
          write(1,1024)fm
          write(1,1025)separ
          write(1,4025)gamma
          write(1,5000)t3
          write(1,5001)g3
          write(1,5002)SA3
          write(1,5003)density
          write(1,5004)sw1
          write(1,5005)sw2
          write(1,5006)sw3
          write(1,5007)T0
          write(1,1026)idraw
          write(1,1027)iecheck
          write(1,1028)idint
          write(1,4000)iatm
          write(1,4001)ism1
          write(1,4002)icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK
          write(1,5008)iRVfilt
          write(1,5009)isw1
          write(1,5010)isw2
          write(1,5011)isw3
          write(1,5012)isw4
          write(1,3028)ilaw
c
          do 10 i=1,8
            write(1,2000)wave(i),dbolx(i,1),dboly(i,1),dbolx(i,2),dboly(i,2),
     $           dwavex(i,1),dwavey(i,1),dwavex(i,2),dwavey(i,2)
 10       continue
c
          write(1,6000)ecc
          write(1,6001)argper
          write(1,6002)pshift
          write(1,6003)sw5
          write(1,6004)sw6
          write(1,6005)sw7
          write(1,6006)sw8
          write(1,6007)sw9
          write(1,6008)ikeep
          write(1,6009)isynch
          write(1,6010)isw5
          write(1,6011)isw6
          write(1,6012)isw7
          write(1,6013)isw8
          write(1,6014)isw9
c
          write(1,2001)spot1parm(1,1)
          write(1,2002)spot1parm(1,2)
          write(1,2003)spot1parm(1,3)
          write(1,2004)spot1parm(1,4)
c                     
          write(1,2005)spot1parm(2,1)
          write(1,2006)spot1parm(2,2)
          write(1,2007)spot1parm(2,3)
          write(1,2008)spot1parm(2,4)
c                     
c                     
          write(1,2009)spot2parm(1,1)
          write(1,2010)spot2parm(1,2)
          write(1,2011)spot2parm(1,3)
          write(1,2012)spot2parm(1,4)
                      
          write(1,2013)spot2parm(2,1)
          write(1,2014)spot2parm(2,2)
          write(1,2015)spot2parm(2,3)
          write(1,2016)spot2parm(2,4)
c                     
c                     
          write(1,2017)spotdparm(1,1)
          write(1,22018)spotdparm(1,2)
          write(1,2019)spotdparm(1,3)
          write(1,2020)spotdparm(1,4)
                      
          write(1,2021)spotdparm(2,1)
          write(1,2022)spotdparm(2,2)
          write(1,2023)spotdparm(2,3)
          write(1,2024)spotdparm(2,4)
c
c   UPDATE August 10, 2004
c
c   Write out the 8 real and 4 new integer variables here.
c
          write(1,2025)primmass
          write(1,2026)primK
          write(1,2027)primrad
          write(1,2028)ratrat
c
          write(1,2030)frac1
          write(1,2031)frac2
          write(1,2032)ecosw
          write(1,2033)temprat
c
          write(1,2040)idark1
          write(1,2041)idark2
          write(1,2042)isw12
          write(1,3043)isw13
c
          write(1,8001)isw21
          write(1,8002)isw22
          write(1,8003)isw23
          write(1,8004)isw24
c
8001      format(i1,19x,'ialign (0 for rotation aligned with orbit)')
8002      format(i1,19x,'ifastgen (1 for fast genetic mode)')
8003      format(i1,19x,'isw23 (currently inactive)')         
8004      format(i1,19x,'frac switch (>1 to enable ELCratio.???? files)')         
c
          do 85000 kk=1,8
            write(1,85001)(powercoeff(kk,jj),jj=1,9)
85000     continue
c
85001      format(9(f8.5,1x))
c
          write(1,8011)bigI
          write(1,8012)bigbeta
          write(1,8013)sw23
          write(1,8014)sw24
c
8011      format(f11.7,9x,'axis_I (inclination of rotation axis if ialign=1)')
8012      format(f11.7,9x,'axis_beta (angle of rotation axis wrt to orbit if',
     #      ' ialign=1)')          
8013      format(f15.7,5x,'t_start')
8014      format(f15.7,5x,'t_end')         
c
          write(1,8025)sw25
          write(1,8026)sw26
          write(1,8027)sw27
          write(1,8028)sw28
          write(1,8029)sw29
          write(1,8030)sw30
          write(1,8031)sw31
          write(1,8032)Tconj
          write(1,8033)beam1
          write(1,8034)beam2

          write(1,9025)isw25
          write(1,9026)isw26
          write(1,9027)isw27
          write(1,9028)isw28
          write(1,9029)isw29
          write(1,9030)isw30
          write(1,9031)isw31
          write(1,9032)isw32
          write(1,9033)isw33
          write(1,9034)isw34
8025      format(f10.7,19x,'asini error')         
8026      format(f10.8,10x,'reference phase for disk fraction')
8027      format(f10.8,10x,'radfill1 (set to use fill1 in terms of R_eff')
8028      format(f10.8,10x,'radfill2 (set to use fill2 in terms of R_eff')
8029      format(f10.4,10x,'bin size for light curves (minutes)')
8030      format(f10.4,10x,'bin size for RV curves (minutes)')
8031      format(f4.2,16x,'sw31 (currently inactive)')
8032      format(f15.8,5x,'Tconj')
8033      format(f4.2,16x,'beam1 (Doppler boost factor, star 1')
8034      format(f4.2,16x,'beam2 (Doppler boost factor, star 2)')

9025      format(i1,19x,'isw25 (currently inactive)')
9026      format(i1,19x,'isw26 (currently inactive)')
9027      format(i6,14x,'Nterms for fast analytic')
9028      format(i1,19x,'set to 1 to fit for Tconj')
9029      format(i1,19x,'isw29 (currently inactive)')
9030      format(i1,19x,'isw30 (currently inactive)')
9031      format(i1,19x,'isw31 (currently inactive)')
9032      format(i1,19x,'isw32 (currently inactive)')
9033      format(i1,19x,'isw33 (currently inactive)')
9034      format(i1,19x,'isw34 (currently inactive)')

c
 2025     format(f13.9,7x,'primmass (star 1 mass in solar masses)')
 2026     format(f14.9,6x,'primK (K-velocity of star 1 in km/sec)')
 2027     format(f14.9,6x,'primrad (star 1 radius in solar radii)')
 2028     format(f16.9,4x,
     &          'ratrad (ratio of star 1 radius and star 2 radius)')
c
 2030     format(f4.2,16x,'frac1 (fractional radius star 1: R_1/a)')
 2031     format(f4.2,16x,'frac2 (fractional radius star 2: R_2/a)')
 2032     format(f12.9,8x,'ecosw (phase difference between eclipses)')
 2033     format(f10.7,10x,'temprat (T_1/T_2)')
c
 2040     format(i1,19x,'idark1')
 2041     format(i1,19x,'idark2')
 2042     format(i6,14x,'Npoly (0 for numerical)')
 3043     format(i1,19x,'ifasttrans (>0 for fast transit mode)')

          close(1)
c
 100      format(a1,'I can''t find the file ''ELC.inp''!  I''m making',
     $     ' one up and setting default values')
c
 1000     format(i2,18x,'Nalph1')
 1001     format(i2,18x,'Nbet1')
 1002     format(i2,18x,'Nalph2')
 1003     format(i2,18x,'Nbet2')
 1004     format(f6.4,14x,'fill1')
 1005     format(f6.4,14x,'fill2')
 1006     format(f7.2,13x,'omega1')
 1007     format(f7.2,13x,'omega2')
 1008     format(f4.2,16x,'dphase')
 1009     format(f4.2,16x,'Q')
 1010     format(f5.2,15x,'finc')
 1011     format(f6.1,14x,'Teff1')
 1012     format(f6.1,14x,'Teff2')
 1013     format(f4.2,16x,'Tgrav1')
 1014     format(f4.2,16x,'Tgrav2')
 1015     format(f4.2,16x,'betarim')
 1016     format(f5.3,15x,'rinner')
 1017     format(f4.2,16x,'router')
 1018     format(f7.1,13x,'tdisk')
 2018     format(f7.4,13x,'xi')
 1019     format(i3,17x,'Ntheta')
 1020     format(i3,17x,'Nradius')
c
c   UPDATE March 26, 2002
c
c   The variable rLx now means the log10 of the X-ray luminosity.
c
c 1021     format(f10.5,10x,'Lx/Lopt')
 1021     format(f10.5,10x,'log10(Lx)')
 1022     format(f4.2,16x,'W')
 1023     format(f8.6,12x,'Period')
 1024     format(f5.3,15x,'fm')
 1025     format(f5.3,15x,'separ')
 1026     format(i1,19x,'idraw')
 1027     format(i1,19x,'iecheck')
 1028     format(i1,19x,'idint')
 1029     format(i1,19x,'ivelout')
 1030     format(i1,19x,'iXout')
 1040     format(f6.4,14x,'darkbol1')
 1041     format(f6.4,14x,'darkbol2')
 1042     format(f6.4,14x,'alb1')
 1043     format(f6.4,14x,'alb2')
 2043     format(i1,19x,'Nref')
 2000     format(f7.1,3x,8(f7.4,1x))
 3028     format(i2,18x,'ilaw  (1=linear law, 2=logarithmic law,',
     %           ' 3=square root law, 4=quad law, >10 for power series)')
 4025     format(f7.2,13x,'gamma velocity (km/sec)')
 4000     format(i1,19x,'iatm  (0 for BB, 1 for model atmospheres)')
 4001     format(i1,19x,'ism1  (0 for all phases, 1 for 0-180)')
 4002     format(8(i1,1x),4x,'icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK')
 5000     format(f10.2,10x,'t3')
 5001     format(f5.2,15x,'g3')
 5002     format(f12.6,8x,'SA3')
 5003     format(f12.6,8x,'density in g/cc')
 5004     format(f12.6,8x,'onephase')
 5005     format(f12.6,8x,'usepot1')
 5006     format(f12.6,8x,'usepot2')
 5007     format(f12.6,8x,'T0')
 5008     format(i1,19x,'iRVfilt')
 5009     format(i1,19x,'ionephase')
 5010     format(i1,19x,'isquare')
 5011     format(i1,19x,'iusepot')
 5012     format(i1,19x,'ifixgamma')
 6000     format(f10.8,10x,'eccentricity')
 6001     format(f13.8,7x,'argument of peristron in degrees')
 6002     format(f11.8,9x,'pshift')
c
c   UPDATE April 8, 2002
c
c   Change the format statement below
c
 6003     format(f12.6,8x,'asini (projected semimajor axis in seconds)')
c
c   UPDATE May 27, 2002
c
c   Change this format statement.
c
 6004     format(f12.6,8x,'median fit (geneticELC only)')
 6005     format(f12.6,8x,'sw7 (currently inactive)')
 6006     format(f12.6,8x,'sw8 (currently inactive)')
 6007     format(f12.6,8x,'sw9 (currently inactive)')
 6008     format(i1,19x,'ikeep (1 to put eclipse at phase 0.0)')
 6009     format(i1,19x,'isynch (1 to keep rotation synchronous',
     $         ' at periastron)')
c
c   RVG BUG ALERT  May 8, 2001
c
c   Change these two format statements.
c
 6010     format(i1,19x,'ispotprof')
 6011     format(i1,19x,'igrav')
c
c   END BUG
c
 6012     format(i1,19x,'itime')
c
c  UPDATE JULY 4, 2004
c
c  The variable isw8 will be assigned to MonteCarlo, which
c  will be used to determine the number of Monte Carlo
c  integrations done on the pixels.
c
 6013     format(i6,14x,'MonteCarlo (0 for interpolation, >10 ',
     $        'for Monte Carlo)')
 6014     format(i1,19x,'ielite')
c
c  RVG BUG ALERT  May 9, 2001
c
c  Add these format statements.
c
 2001     format(f10.7,10x,'Temperature factor spot 1, star 1')
 2002     format(f11.7, 9x,'Latitude of spot 1, star 1 (degrees)') 
 2003     format(f11.7, 9x,'Longitude of spot 1, star 1 (degrees)') 
 2004     format(f11.7, 9x,'Angular radius of spot 1, star 1 (degrees)') 

 2005     format(f10.7,10x,'Temperature factor spot 2, star 1')
 2006     format(f11.7, 9x,'Latitude of spot 2, star 1 (degrees)') 
 2007     format(f11.7, 9x,'Longitude of spot 2, star 1 (degrees)') 
 2008     format(f11.7, 9x,'Angular radius of spot 2, star 1 (degrees)') 

 2009     format(f10.7,10x,'Temperature factor spot 1, star 2')
 2010     format(f11.7, 9x,'Latitude of spot 1, star 2 (degrees)') 
 2011     format(f11.7, 9x,'Longitude of spot 1, star 2 (degrees)') 
 2012     format(f11.7, 9x,'Angular radius of spot 1, star 2 (degrees)') 

 2013     format(f10.7,10x,'Temperature factor spot 2, star 2')
 2014     format(f11.7, 9x,'Latitude of spot 2, star 2 (degrees)') 
 2015     format(f11.7, 9x,'Longitude of spot 2, star 2 (degrees)') 
 2016     format(f11.7, 9x,'Angular radius of spot 2, star 2 (degrees)') 

 2017     format(f10.7,10x,'Temperature factor spot 1, disk')
22018     format(f11.7, 9x,'Azimuth of spot 1, disk (degrees)') 
 2019     format(f11.7, 9x,'Radial cutoff of spot 1, disk (0 <= r_cut <=1)') 
 2020     format(f11.7, 9x,'Angular size of spot 1, disk (degrees)') 

 2021     format(f10.7,10x,'Temperature factor spot 2, disk')
 2022     format(f11.7, 9x,'Azimuth of spot 2, disk (degrees)') 
 2023     format(f11.7, 9x,'Radial cutoff of spot 2, disk (0 <= r_cut <=1)') 
 2024     format(f11.7, 9x,'Angular size of spot 2, disk (degrees)') 
c

          return
          end
c
c   &&&&&&&&&&
c
          SUBROUTINE POTEN(Q,omega,x,y,z,psi,psix,psixx,psiy,psiz,istar,
     #       bdist,cox,coy,tidephi,itide)
c
c   October 6, 1999
c
c   Here is the subroutine from the original Avni code.  It computes
c   the potential (given Q and x,y,z) and various gradients.
c
c   February 15, 2000
c
c   The iflag is added because the potential computation needs to be
c   modified in some cases, specifically when star 2 is not synchronous.
c
c   UPDATE November 13, 2009
c
c   Add an option to use the potential computed from the equilibrium
c   tide approximation.   The flag itide will be the order of the polynomial
c   expansion.
c
          implicit double precision(a-h,o-z)
c
          dimension pn(0:100),pd(0:100)
          parameter(pie=3.141592653589793d0)

          if(itide.lt.2)then
            RST = DSQRT(X**2 + Y**2 + Z**2)     !dist. from center of 2ndary
            RX = DSQRT((X-bdist)**2 + Y**2 + Z**2)!dist. from center of primary
            A = ((1.0d0+Q)/2.0d0) * OMEGA**2
            RST3 = RST*RST*RST
            RX3 = RX*RX*RX
            PSI = 1.0d0/RST + Q/RX - Q*X/bdist/bdist 
     #       + A*(X**2 + Y**2) !potential (page 45 of  Avni's paper
            PSIY = -Y/RST3 - Q*Y/RX3   + 2.0d0*A*Y        !partial deriv. wrt y
            PSIZ = -Z/RST3 - Q*Z/RX3                     !partial deriv. wrt z
            PSIX = -X/RST3 - Q*(X-bdist)/RX3 -Q/bdist/bdist 
     #           + 2.0d0*A*X !partial deriv wrt x
            if(istar.eq.2)then
               PSIX = -X/RST3 - Q*(X-1.0d0)/RX3 - 2.0d0*A*(1.0d0-X) +1.0d0  
            endif
c
            RST5 = RST3*RST*RST
            RX5 = RX3*RX*RX
            PSIXX = -1.0d0/RST3 + 3.0d0*X**2/RST5 
     $        -Q/RX3 + (3.0d0*Q*(X-bdist)**2)/RX5 +2.0d0*A
 
            RETURN
          endif
c
          tider=pie*tidephi/180.0d0
          psicos=cox*dcos(tider)+coy*dsin(tider)


c
          rr=x*x+y*y+z*z
          w1=(dsqrt(rr))**3
          w1=1.0d0/w1
c
          psi=1.0d0/dsqrt(rr)+Q
          psix=-x*w1
          psiy=-y*w1
          psiz=-z*w1
          psixx=w1*(3.0d0*x*x/rr-1.0d0)
c

          call lpn(itide,psicos,pn,pd)
          do 10 ii=2,itide
            t1=-Q*rr**(0.5d0*dble(ii))*pn(ii)
            psi=psi+t1
            twx=dble(ii)*x*Q*rr**(0.5d0*dble(ii)-1.0d0)*pn(ii)
            twy=dble(ii)*y*Q*rr**(0.5d0*dble(ii)-1.0d0)*pn(ii)
            twz=dble(ii)*z*Q*rr**(0.5d0*dble(ii)-1.0d0)*pn(ii)
            txx=Q*pn(ii)*(rr)**(0.5d0**dble(ii)-1)
            txx=txx*(2.0d0*x*x*(0.5d0*dble(ii)-1.0d0)/rr+dble(ii))
            psix=psix+twx
            psiy=psiy+twy
            psiz=psiz+twz
            psixx=psixx+txx
 10       continue
c
 99       RETURN
          END
c
c  %%%%%%%%%%%%%%%%%%%%%%%
c
c  UPDATE September 11, 2001
c
c  Add the iverb flag to the argument list.  If iverb=1, then don't write
c  to unit 2.
c
c  UPDATE September 21, 2008
c
c  Add sw5 to the argument list
c
          subroutine setupgeo(istar,ialphmax,ibetmax,
     $      Nalph,Nbet,ibetlim,fill,omega,Q,finc,
     $      xarray,yarray,zarray,
     %      surf,radarray,gradx,grady,gradz,garray,
     $      xend,separation,Tgrav,Teff,reff,Rl,Tpole,Rpol,Regg,sarea,pot,
     $      gpole,phiar,isquare,iusepot,usepot,ivrt,pervol,fillper,bdist,
     $      potsum,iverb,mmdx,primmass,primK,primrad,ratrad,frac1,frac2,
     &      ecc,period,size1,fill2,omegatwo,sw5,tteff2,density,
     #      tidephi,itide,phistart)
c
c
c
c   October 6, 1999
c
c   This subroutine will return the radii and surface coordinates
c   for star 1 or 2:
c
c   xarray, yarray, zarray, surf, gradx, grady, gradz, radarray, garray
c
c   The two element array xend will contain the x-coordinates of the
c   two points along the x-axis.  The gradx term is +1 at positive x-axis
c   and -1 at negative x-axis.
c
c   istar is the flag to identify which star is being computed.  If
c   istar=2, then we need to flip the mass ratio.
c
c   February 15, 2000
c
c   I have added 'dummy' arguments for the x-gradient, gravity, and surface
c   elements.  When star 2 has non-synchronous rotation, the form of the
c   x derivative of the potential needs to be modified.  This modified
c   form is only  needed for the reflection effect.  The computation of
c   the light curve is otherwise OK with the standard gradients.
c
c   ********removed May 11, 2000*************
c
c   May 4, 2000
c
c   Change the distribution of grid points. Use a normal polar coordinate
c   scheme.
c
c   ************************
c
c   February 5, 2001
c
c   Generalize to eccentric orbits.  Need extra parameters:
c
c   ivrt=0, normal mode,  
c   ivrt=1, find the filling so that the volume is equal to pervol
c   fillper is the filling factor at periastron   
c
c   UPDATE August 12, 2004
c
c   Add the arguments frac1,frac2,ratrad, etc. to the argument list.
c   These are the fractional radii, and the ratio of the radii.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265359879323d0,twopie=pie+pie)
          dimension xarray(ialphmax*ibetmax),yarray(ialphmax*ibetmax),
     $      zarray(ialphmax*ibetmax),surf(ialphmax*ibetmax),
     &      gradx(ialphmax*ibetmax),xend(4),ibetlim(ialphmax),
     $      grady(ialphmax*ibetmax),gradz(ialphmax*ibetmax),
     %      radarray(ialphmax*ibetmax),garray(ialphmax*ibetmax),
     $      Rpol(ialphmax),phiar(ialphmax*ibetmax)
          dimension mmdx(ialphmax,ibetmax),phistart(ialphmax)
c
c   UPDATE January 11, 2010
c
c   Add a random offset to the phi coordinates in each ialf row.  This
c   will reduce the noise in high precision curves.  We need the sobel
c   sequence for this.
c
          dimension xsob(2)
c
c
c   UPDATE December 21, 2008
c
c   Add this routine if radfill is given.  radfill defines the filling factor
c   in terms of the effective radius.  Adjust fill1 and fill2 accordingly
c
c
c
c
c   UPDATE October 10, 2008
c
c   Add this new call to the routine that sets the mass ratio given
c   the density and fill factor of star 1.
c
          if((istar.eq.1).and.(frac1.le.0.0d0).and.(primmass.le.0.0d0)
     #       .and.(density.gt.0.0d0).and.(primrad.le.0.0d0))then
            call setdensity(fill,omega,bdist,Q,period,density,tidephi,itide)
          endif
c
c
c  UPDATE August 12, 2004
c
c  If istar=1 and IVRT=0, set the separation and Q depending
c  on the values of primmass and primK.
c
          if(istar.eq.1.and.ivrt.eq.0)then
            call setscale(fill1,fill2,omega1,omega2,
     &         Q,finc,Teff1,tteff2,Period,fm,separation,
     &         primmass,primK,primrad,ratrad,sw5,reff1,ecc,bdist)
          endif
c
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
c
c   Set up the grid on the star, define step size and number of steps.
c
          dtheta=pie/dble(Nalph)
          dcostheta=2.0d0/dble(Nalph)
c
c   Initialize the bisection limits.
c
          fbig=fillper
          fsmall=0.0d-8
          fnew=0.5d0*(fbig+fsmall)
c
c   First, find the distance to L1.
c
          call findL1(overQ,omega,x0,1,bdist,tidephi,itide)       
          Rl=x0                           !save the L1 distance
c
c          write(*,6968)bdist,omega,x0,overQ
6968      format(f9.7,1x,f9.6,1x,f9.6,1x,f15.7)
          x=x0
          y=0.0d0
          z=0.0d0
          cox=1.0d0
          coy=0.0d0
          call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,bdist,
     #        cox,coy,tidephi,itide)
          savepsi0=psi     ! value of the potential at x-axis 
          critpsi=psi
c
c
c  If primrad=0, then we can adjust the filling factors to
c  get fractional radii.
c
          if((istar.eq.1).and.(ivrt.eq.0).and.frac1.gt.0.0d0)then
             call findradius(overQ,omega,savepsi0,x0,bdist,reff,tidephi,itide)
             bigfrac=reff
             if(frac1.gt.bigfrac)then
               fill=1.0d0
               size1=reff
               write(2,665)
               go to 666
 665           format('Warning:  frac1 exceeds Roche radius. ',
     $           ' Setting fill1=1.0')
             endif
c
c   Presumably the Roche lobe has a fractional radius large than the
c   requested fractional radius.  Use bisection to find the filling
c   factor fill1.
c
             midfrac=reff
             aa=0.0d0
             bb=1.0d0
             y=0.0d0
             z=0.0d0
             do 664 i=1,60
               fff=(aa+bb)*0.5d0
c               write(*,*)fff,fff*X0
               call POTEN(overQ,omega,x0*fff,y,z,psi,psix,psixx,psiy,
     &             psiz,1,bdist,cox,coy,tidephi,itide)
               call findradius(overQ,omega,psi,x0*fff,bdist,reff,tidephi,itide)
               if(frac1.ge.reff)then
                 aa=fff
               else
                 bb=fff
               endif
c               write(*,*)aa,bb,reff
 664         continue
c             
             fill=fff
             size1=reff
             write(2,663)frac1,fill
 663         format('Info:  frac1 = ',f10.8,', fill1 set to ',f10.8)
             go to 666
          endif
c
          if((istar.eq.2).and.(ivrt.eq.0).and.frac2.gt.0.0d0)then
             call findradius(overQ,omega,savepsi0,x0,bdist,reff,tidephi,itide)
             bigfrac=reff
             if(frac2.gt.bigfrac)then
               fill=1.0d0
               write(2,6665)
               go to 666
 6665          format('Warning:  frac2 exceeds Roche radius. ',
     $           ' Setting fill2=1.0')
             endif
c
c   Presumably the Roche lobe has a fractional radius large than the
c   requested fractional radius.  Use bisection to find the filling
c   factor fill1.
c
             midfrac=reff
             aa=0.0d0
             bb=1.0d0
             y=0.0d0
             z=0.0d0
             do 6664 i=1,60
               fff=(aa+bb)*0.5d0
c               write(*,*)fff,fff*X0
               call POTEN(overQ,omega,x0*fff,y,z,psi,psix,psixx,psiy,
     &             psiz,1,bdist,cox,coy,tidephi,itide)
               call findradius(overQ,omega,psi,x0*fff,bdist,
     &           reff,tidephi,itide)
               if(frac2.ge.reff)then
                 aa=fff
               else
                 bb=fff
               endif
c               write(*,*)aa,bb,reff
 6664        continue
c             
             fill=fff
             write(2,6663)frac2,fill
 6663        format('Info:  frac2 = ',f10.8,', fill2 set to ',f10.8)
             go to 666
          endif
c
c  If primrad>0, then we can adjust the fill1 to
c  get the requested radius, if possible.  However, if frac1>0,
c  then don't set the radius to primrad
c
          if((istar.eq.1).and.(ivrt.eq.0).and.(frac1.le.0.0d0).
     &          and.(primrad.gt.0.0d0))then
c
             call findradius(overQ,omega,savepsi0,x0,bdist,reff,tidephi,itide)
             bigrad=reff*separation
             if(primrad.gt.bigrad)then
               fill=1.0d0
               write(2,4665)
               go to 666
 4665          format('Warning:  primrad*separation exceeds Roche radius. ',
     $           ' Setting fill1=1.0')
             endif
c
c   Presumably the Roche lobe has a  radius larger than the
c   requested radius.  Use bisection to find the filling
c   factor fill1.
c
             midfrac=reff
             aa=0.0d0
             bb=1.0d0
             y=0.0d0
             z=0.0d0
             do 4664 i=1,60
               fff=(aa+bb)*0.5d0
c               write(*,*)fff,fff*X0,istar
               call POTEN(overQ,omega,x0*fff,y,z,psi,psix,psixx,psiy,
     &             psiz,1,bdist,cox,coy,tidephi,itide)
               call findradius(overQ,omega,psi,x0*fff,bdist,reff,
     &           tidephi,itide)
               if(primrad.ge.reff*separation)then
                 aa=fff
               else
                 bb=fff
               endif
c               write(*,*)aa,bb,reff
 4664         continue
c             
             fill=fff
             size1=reff
             write(2,4663)primrad,fill
 4663        format('Info:  primrad = ',f13.8,', fill1 set to ',f10.8)
             go to 666
          endif
c
Cc  If primrad=0 and frac1=0 and ratrad>0, then we can adjust the fill1 to
Cc  get the requested radius, if possible, after checking for the.radius
Cc  of star 2.
Cc
          if((istar.eq.1).and.(ivrt.eq.0).and.(frac1.le.0.0d0).
     &          and.(primrad.le.0.0d0).and.(ratrad.gt.0.0d0).
     $          and.(frac2.le.0.0d0))then
c
             cox=1.0d0
             coy=0.0d0
             y=0.0d0
             z=0.0d0
             call POTEN(overQ,omega,x0*fill,y,z,psi,psix,psixx,psiy,
     &             psiz,1,bdist,cox,coy,tidephi,itide)
             call findradius(overQ,omega,psi,fill*x0,bdist,
     &        reff,tidephi,itide)

             size1=reff
             write(2,3663)ratrad,fill
c             write(*,*)fill,size1
 3663        format('Info:  ratrad = ',f13.8,', fill1 set to ',f10.8)
             go to 666
          endif
c
c  If primrad=0 and frac1=0 and and frac2> and ratrad>0, 
c  then we can adjust the fill1 to
c  get the requested radius, if possible, after checking for the.radius
c  of star 2, which is frac2.
c
          if((istar.eq.1).and.(ivrt.eq.0).and.(frac1.le.0.0d0).
     &          and.(primrad.le.0.0d0).and.(ratrad.gt.0.0d0).
     $          and.(frac2.gt.0.0d0))then
c
c    Now use bisection to adjust fill1 so that reff1/reff2=ratrad.
c
             call findradius(overQ,omega,savepsi0,x0,bdist,
     &         reff,tidephi,itide)
             radneed=ratrad*frac2
             if(reff.lt.radneed)then
               fill=1.0d0
               write(2,2665)ratrad,frac2
               go to 666
 2665          format('Warning:  radrad=',f12.8,' is not possible with'/,
     $                '          frac2=',f8.6,'. ',
     &                 ' Setting fill1=1.0')
             endif
             midfrac=reff
             aa=0.0d0
             bb=1.0d0
             y=0.0d0
             z=0.0d0
             do 2664 i=1,60
               fff=(aa+bb)*0.5d0
c               write(*,*)fff,fff*X0
               call POTEN(overQ,omega,x0*fff,y,z,psi,psix,psixx,psiy,
     &             psiz,1,bdist,cox,coy,tidephi,itide)
               call findradius(overQ,omega,psi,x0*fff,bdist,
     &          reff,tidephi,itide)
               radneed=ratrad*frac2
               if(radneed.ge.reff)then
                 aa=fff
               else
                 bb=fff
               endif
 2664         continue
c             
             fill=fff
             size1=reff
             write(2,2663)ratrad,fill
 2663        format('Info:  ratrad = ',f13.8,', fill1 set to ',f10.8)
             go to 666
          endif
c
c   
c
c  If ratrad>0 and frac1=0
c  then we can adjust the fill2 to
c  get the requested radius, if possible, based on the value of size1
c
          if((istar.eq.2).and.(ivrt.eq.0).and.(ratrad.gt.0.0d0).
     $          and.(frac2.le.0.0d0))then
c
c    Now use bisection to adjust fill2 so that rad2=size1/ratrad.
c
             call findradius(overQ,omega,savepsi0,x0,bdist,
     &         reff,tidephi,itide)
             radneed=size1/ratrad

c             write(*,*)'&&&& ',radneed,size1,reff
             if(reff.lt.radneed)then
               fill=  1.0d0        
               write(2,1665)ratrad
               go to 666
 1665          format('Warning:  radrad=',f12.8,' is not possible.'/,
     $             ' Setting fill2=1.0')
             endif
             midfrac=reff
             aa=0.0d0
             bb=1.0d0
             y=0.0d0
             z=0.0d0
             do 1664 i=1,60
               fff=(aa+bb)*0.5d0
c               write(*,*)fff,fff*X0
               call POTEN(overQ,omega,x0*fff,y,z,psi,psix,psixx,psiy,
     &             psiz,1,bdist,cox,coy,tidephi,itide)
               call findradius(overQ,omega,psi,x0*fff,bdist,
     &           reff,tidephi,itide)
               radneed=size1/ratrad
               if(radneed.ge.reff)then
                 aa=fff
               else
                 bb=fff
               endif
 1664         continue
c             
             fill=fff
             write(2,1663)ratrad,fill
 1663        format('Info:  ratrad = ',f13.8,', fill2 set to ',f10.8)
             go to 666
          endif
c
c  
c   Skip down to here when constraints are set
c
 666      nnn=40
          if(ivrt.eq.0)nnn=1
          do 106 iii=1,nnn   ! here is the loop point when we need to invert
c
            if(ivrt.gt.0)fill=fnew
c
c   If the flag 'iusepot' is 1 or larger, then we must compute the
c   filling factor needed to get the entered value of 'usepot'.
c   If the requested potential is *smaller* than the critical potential,
c   then set the filling factor to 1.0.
c
            if(iusepot.ge.1)then
              if(ivrt.eq.0)then
                call findfill(istar,overQ,omega,critpsi,x0,usepot,psi0,fill,
     %              bdist,tidephi,itide)
                x=fill*x0
              else
                x=fill*x0
                y=0.0d0
                z=0.0d0
                call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,
     #             bdist,cox,coy,tidephi,itide)
                psi0=psi     ! value of the potential at x-axis 
              endif
            else
              x=fill*x0
              y=0.0d0
              z=0.0d0
              call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,
     #           bdist,cox,coy,tidephi,itide)
              psi0=psi     ! value of the potential at x-axis 
            endif
c
            if(ivrt.eq.0)fillper=fill
c
            vol=0.0d0
            sarea=0.0d0
            r=x       
            potsum=0.0d0  
            fincr=finc*pie/180.00000d0
            do 2526 ialf=1,nalph
              theta=-0.5d0*dtheta+dtheta*dble(ialf)
c
c   UPDATE January 11, 2010
c
c   Change to make similar to W-D
c
              ibetlim(ialf)=1+idnint(1.3d0*dble(4*Nbet)*dsin(theta))
c            
              if(mod(ibetlim(ialf),2).eq.1)ibetlim(ialf)=ibetlim(ialf)-1
              if(isquare.ge.1)ibetlim(ialf)=4*Nbet
c
              difflow=(fincr)-140.0d0*dtheta
              diffhigh=(fincr)+140.0d0*dtheta
c
c              if((theta.gt.difflow).and.(theta.lt.diffhigh))then
c                write(*,*)istar,ialf
c                ibetlim(ialf)=4*ibetlim(ialf)
c              endif
c
c              if((ialf.gt.74).and.(ialf.lt.78))ibetlim(ialf)=4*ibetlim(ialf)
c              if(ialf.gt.nalph/2)then
c                if(mod(ialf,2).eq.1)then
c                  ibetlim(ialf)=4*Nbet
c                  ibetlim(nalph-ialf+1)=4*Nbet
c                endif
c              endif
              theta=-0.5d0*dtheta+dtheta*dble(ialf)
 2526       continue
c
c            do 2528 ialf=i,nalph
c              do 2527 ibet=1,ibetlim(ialf)
c                iidx=kount(ialphmax,ialf,ibetlim)+ibet
c                if(ibet.eq.1)write(*,*)ialf,ibet,iidx,ibetlim(ialf)
c 2527         continue
c 2528       continue
c
c   UPDATE January 11, 2010
c
c   Now that ibetlim is set, initalize the random starting phi values
c
          nnn=2
          step=0.25d0/dble(Nalph/2)
          start=0.0d0
          do 765 ialf=1,Nalph
            dphi=twopie/dble(ibetlim(ialf))
            call sobseq(nnn,xsob)
            fac=(xsob(1)-0.5d0)*0.5d0
            fac=0.0d0
            phistart(ialf)=fac*dphi
c            write(*,*)istar,phistart(ialf),(xsob(1)-0.5d0)
c            phistart(ialf)=0.0d0
c            if(dble(ialf/2).eq.(dble(ialf)*0.5d0))phistart(ialf)=fac*dphi 
765       continue
c
c   The following is a quick loop to find the volume.
c
            mcount=0
            do 104 ialf=1,Nalph/2
c
c   UPDATE MARCH 17, 2004
c
c   make the initial value of r smaller here.
c
              r=0.0000001d0

              call rad(overQ,omega,0.0d0,0.0d0,1.0d0,psi0,r,x,y,z,1,bdist,
     $            tidephi,itide)
              theta=-0.5d0*dtheta+dtheta*dble(ialf)
c              ibetlim(ialf)=idnint(dsin(theta)*4*Nbet)
c              if(mod(ibetlim(ialf),2).eq.1)ibetlim(ialf)=ibetlim(ialf)-1
c              if(isquare.ge.1)ibetlim(ialf)=4*Nbet
c              if(ialf.le.nalf/2)then
c                if(mod(ialf,2).eq.0)ibetlim(ialf)=4*Nbet
c              else
c                if(mod(ialf,2).eq.1)ibetlim(ialf)=4*Nbet
c              endif
              dphi=twopie/dble(ibetlim(ialf))
              snth=dsin(theta)
              snth3=snth/3.0d0   
              cnth=dcos(theta)
              DO 105 ibet=1,ibetlim(ialf)/2          !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c                iidx=kount(ialphmax,ialf,ibetlim)+ibet
c 
c                mcount=mcount+1
c                mmdx(ialf,ibet)=mcount
c                iidx=mcount
                iidx=kount(ialphmax,ialf,ibetlim)+ibet
                phi=-0.5d0*dphi+dphi*dble(ibet)
c                phiar(ialf,ibet)=phi
                phiar(iidx)=phi
                cox=dcos(phi)*snth             !*dsin(theta)
                coy=dsin(phi)*snth             !*dsin(theta)
                coz=cnth                       !dcos(theta)
                CALL RAD(overQ,omega,cox,coy,coz,psi0,r,x,y,z,1,bdist,
     $            tidephi,itide)
                VOL = VOL + 4.0d0*R*R*R*dphi*dtheta*snth3
105           CONTINUE    ! continue ibet loop
104         CONTINUE                   ! continue over ialf
c
            diff=vol-pervol          
            if(ivrt.eq.0)then
              pervol=vol !reference volume
              vol=0.0d0
              sarea=0.0d0
              r=x       
              potsum=0.0d0  
              x=fill*x0
              y=0.0d0
              z=0.0d0
              go to 700
            endif
c
c   Here is a little bisection block which compares the volume just found
c   to the reference volume pervol.
c
            if(diff.le.0.0d0)then
              fsmall=fnew
              fnew=0.5d0*(fsmall+fbig)
            else
              fbig=fnew
              fnew=0.5d0*(fsmall+fbig)
            endif
            acc=0.5d0*dabs(fbig-fsmall)/(fbig+fsmall)
            if(acc.lt.1.0d-10)then      
              vol=0.0d0
              sarea=0.0d0
              r=x       
              potsum=0.0d0  
              go to 700
            endif
 106      continue
c
c   We have the correct volume now.  Go and assign the other variables.
c
 700      vol=0.0d0
          sarea=0.0d0
          potsum=0.0d0
          x=fill*x0
          y=0.0d0
          z=0.0d0
c
c   Note the loop below is set for the whole star.  If using one fourth
c   of the star, change the term in the Vol = command to 4.0
c
          initbet=0
          mmcount=0
          if(itide.lt.2)then
          do 1104 ialf=1,Nalph               !Nalph/2
c
c   UPDATE MARCH 17, 2004
c
c   make the initial value of r smaller here.
c
            r=0.0000001d0
            call rad(overQ,omega,0.0d0,0.0d0,1.0d0,psi0,r,x,y,z,1,bdist,
     #         tidephi,itide)
            theta=-0.5d0*dtheta+dtheta*dble(ialf)
c            ibetlim(ialf)=idnint(dsin(theta)*4*Nbet)
c            if(mod(ibetlim(ialf),2).eq.1)ibetlim(ialf)=ibetlim(ialf)-1
c            if(isquare.ge.1)ibetlim(ialf)=4*Nbet
c              if(ialf.le.nalf/2)then
c                if(mod(ialf,2).eq.0)ibetlim(ialf)=4*Nbet
c              else
c                if(mod(ialf,2).eq.1)ibetlim(ialf)=4*Nbet
c              endif
            dphi=twopie/dble(ibetlim(ialf))
            snth=dsin(theta)
            snth3=dsin(theta)/3.0d0  !*0.333333333333333d0
            cnth=dcos(theta)
            DO 1105 ibet=1, ibetlim(ialf)     !ibetlim(ialf)/2          !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c              write(*,*)ialf,ibet,iidx
              mmdx(ialf,ibet)=iidx
              phi=-0.5d0*dphi+dphi*dble(ibet)
              phi=phi+phistart(ialf)
c
c              if(dble(ialf/2).eq.(dble(ialf)*0.5d0))phi=phi+0.25d0*dphi  !DPHI
c
              phiar(iidx)=phi
              cox=dcos(phi)*snth                !*dsin(theta)
              coy=dsin(phi)*snth                !*dsin(theta)
              coz=cnth                          !dcos(theta)
              CALL RAD(overQ,omega,cox,coy,coz,psi0,r,x,y,z,1,bdist,
     $         tidephi,itide)
              call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,
     #          bdist,cox,coy,tidephi,itide)
              radarray(iidx) = R
              garray(iidx) = DSQRT(PSIX**2+PSIY**2+PSIZ**2)
              oneoverg=1.0d0/garray(iidx)
              GRADX(iidx) = -PSIX*oneoverg
              GRADY(iidx) = -PSIY*oneoverg
              GRADZ(iidx) = -PSIZ*oneoverg
              surf(iidx) = COX*GRADX(iidx)+COY*GRADY(iidx)
     $           + COZ*GRADZ(iidx)
c
c   The following check is for large separations.
c
              if(surf(iidx).lt.0.7d0)surf(iidx)=0.7d0
c
              surf(iidx) = R**2 / surf(iidx)
c
c   Add terms to account for the delta phi and delta theta terms.
c
c   Expression with distribution of points linear in theta.
c
              surf(iidx)=surf(iidx)*dphi*dtheta*snth
c
c   Expression with distribution of points linear in cos(theta)
c
c               surf(iidx)=surf(iidx)*dphi*dcostheta
c
c   Keep track of the surface area and volume.
c
              sarea=sarea+surf(iidx)
c
c   Expression with distribution of points linear in theta.
c
              VOL = VOL + 1.0d0*R*R*R*dphi*dtheta*snth3
c
c   Expression with distribution of points linear in cos(theta)
c
c              VOL = VOL + R*R*R*dphi*dcostheta/3.0d0
c
c   Assign x,y,z coordinates of grid point
c
              xarray(iidx)=x  !radius vector times the direction cosine
              yarray(iidx)=y  ! "" ""          ""   ""
              zarray(iidx)=z  

              go to 1105

c
c   Use symmetry to reflect the various quantities.
c 
              I1=nalph-(ialf-1)
              J2=ibetlim(ialf)-(ibet-1)
c              ibetlim(I1)=ibetlim(ialf)
c
              jjdx1=kount(ialphmax,ialf,ibetlim)+j2
              jjdx2=kount(ialphmax,I1,ibetlim)+ibet              
              jjdx3=kount(ialphmax,I1,ibetlim)+j2

              mmdx(ialf,j2)=jjdx1
              mmdx(I1,ibet)=jjdx2
              mmdx(I1,j2)=jjdx3
c
              xarray(jjdx1)=xarray(iidx)
              xarray(jjdx2)=xarray(iidx)
              xarray(jjdx3)=xarray(iidx)
c
              yarray(jjdx1)=-yarray(iidx)
              yarray(jjdx2)=yarray(iidx)
              yarray(jjdx3)=-yarray(iidx)
c
              zarray(jjdx1)=zarray(iidx)
              zarray(jjdx2)=-zarray(iidx)
              zarray(jjdx3)=-zarray(iidx)
c
              surf(jjdx1)=surf(iidx)
              surf(jjdx2)=surf(iidx)
              surf(jjdx3)=surf(iidx)
c
              garray(jjdx1)=garray(iidx)
              garray(jjdx2)=garray(iidx)
              garray(jjdx3)=garray(iidx)
c
              radarray(jjdx1)=radarray(iidx)
              radarray(jjdx2)=radarray(iidx)
              radarray(jjdx3)=radarray(iidx)
c
              gradx(jjdx1)=gradx(iidx)
              gradx(jjdx2)=gradx(iidx)
              gradx(jjdx3)=gradx(iidx)
c
              grady(jjdx1)=-grady(iidx)
              grady(jjdx2)=grady(iidx)
              grady(jjdx3)=-grady(iidx)
c
              gradz(jjdx1)=gradz(iidx)
              gradz(jjdx2)=-gradz(iidx)
              gradz(jjdx3)=-gradz(iidx)

              sarea=sarea+surf(jjdx1)
              sarea=sarea+surf(jjdx2)
              sarea=sarea+surf(jjdx3)
  
              phiar(jjdx1)=twopie-phiar(iidx)
              phiar(jjdx2)=phiar(iidx)
              phiar(jjdx3)=twopie-phiar(iidx)
c
1105        CONTINUE    ! continue ialf loop
c
1104      CONTINUE                   ! continue over ibet
c
          endif !end if itide < 2
c
c  If we are using the tidal approximation, 
c  we cannot use symmetry.  Loop over the full range of 
c  Nalph and Nbet
c
          if(itide.ge.2)then
c
          tider=pie*tidephi/80.0d0

          do 3104 ialf=1,Nalph
c
c   UPDATE MARCH 17, 2004
c
c   make the initial value of r smaller here.
c
            r=0.0000001d0
            call rad(overQ,omega,0.0d0,0.0d0,1.0d0,psi0,r,x,y,z,1,bdist,
     #         tidephi,itide)
            theta=-0.5d0*dtheta+dtheta*dble(ialf)
            dphi=twopie/dble(ibetlim(ialf))
            snth=dsin(theta)
            snth3=dsin(theta)/3.0d0  !*0.333333333333333d0
            cnth=dcos(theta)
            DO 3105 ibet=1,ibetlim(ialf)
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c              write(*,*)ialf,ibet,iidx
              mmdx(ialf,ibet)=iidx
              phi=-0.5d0*dphi+dphi*dble(ibet)
              phi=phi+phistart(ialf)
              phiar(iidx)=phi
              cox=dcos(phi)*snth                !*dsin(theta)
              coy=dsin(phi)*snth                !*dsin(theta)
              coz=cnth                          !dcos(theta)
              CALL RAD(overQ,omega,cox,coy,coz,psi0,r,x,y,z,1,bdist,
     $         tidephi,itide)
              call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,
     #          bdist,cox,coy,tidephi,itide)
c

              radarray(iidx) = R
              garray(iidx) = DSQRT(PSIX**2+PSIY**2+PSIZ**2)
              oneoverg=1.0d0/garray(iidx)
              GRADX(iidx) = -PSIX*oneoverg
              GRADY(iidx) = -PSIY*oneoverg
              GRADZ(iidx) = -PSIZ*oneoverg
              surf(iidx) = COX*GRADX(iidx)+COY*GRADY(iidx)
     $           + COZ*GRADZ(iidx)
c
c   The following check is for large separations.
c
              if(surf(iidx).lt.0.7d0)surf(iidx)=0.7d0
c
              surf(iidx) = R**2 / surf(iidx)
c
c   Add terms to account for the delta phi and delta theta terms.
c
c   Expression with distribution of points linear in theta.
c
              surf(iidx)=surf(iidx)*dphi*dtheta*snth
c
c   Expression with distribution of points linear in cos(theta)
c
c               surf(iidx)=surf(iidx)*dphi*dcostheta
c
c   Keep track of the surface area and volume.
c
              sarea=sarea+surf(iidx)
c
c   Expression with distribution of points linear in theta.
c
              VOL = VOL + 1.0d0*R*R*R*dphi*dtheta*snth3
c
c   Expression with distribution of points linear in cos(theta)
c
c              VOL = VOL + R*R*R*dphi*dcostheta/3.0d0
c
c   Assign x,y,z coordinates of grid point
c
              xarray(iidx)=x  !radius vector times the direction cosine
              yarray(iidx)=y  ! "" ""          ""   ""
              zarray(iidx)=z  
c

3105        CONTINUE    ! continue ialf loop
c
3104      CONTINUE                   ! continue over ibet
c
          endif !end if itide >= 2


 1107     format(4(i2,1x),1x,3(f7.4,1x),1x,g9.4,3x,5(f9.5,1x))
c
          REFF = (0.75d0*VOL/pie) **(1.0d0/3.0d0) 

c          REFF = (0.248732415d0*VOL) **(1.0d0/3.0d0) ! 0.3333333333333333d0
c
          Nalf2=Nalph/2
c
          cox=0.0d0
          coy=0.0d0
          coz=1.0d0
          r=radarray(1)
          CALL RAD(overQ,omega,cox,coy,coz,psi0,r,x,y,z,1,bdist,tidephi,itide)
          call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,
     #       bdist,cox,coy,tidephi,itide)
          div = DSQRT(PSIX**2+PSIY**2+PSIZ**2)
          gpole=div          
          rpol(Nalf2)=r
          zpole=z
c
          potsum=0.0d0
          DO 401 IALF = 1, nalph
            DO 400 IBET = 1, ibetlim(ialf)          !4*NBET
c
c   Keep track of the sums of the normalized gravities to compute the
c   intensity weighted effective temperature.
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              potsum=((garray(iidx)/div)**(4.0d0*Tgrav))
     %               *surf(iidx)+potsum
c
c              write(52,6969)ialf,ibet,surf(iidx),garray(iidx),gradx(iidx)
 400        continue
401       CONTINUE
 6969     format(2(i5,2x),3(f15.10,1x))
c
c   keep track of the point on the star on x-axis
c
          cox=1.0d0
          coy=0.0d0
          coz=0.0d0
          r=fill*x0
          CALL RAD(overQ,omega,cox,coy,coz,psi0,r,x,y,z,1,bdist,tidephi,itide)
          xend(3)=r
c         
          cox=-1.0d0
          coy=0.0d0
          coz=0.0d0
          CALL RAD(overQ,omega,cox,coy,coz,psi0,r,x,y,z,1,bdist,tidephi,itide)
          xend(4)=r
c
          Tpole=Teff*dsqrt(dsqrt(sarea/potsum))
c
c  UPDATE September 11, 2001
c 
c  Put the if-then clauses around the write statements.
c
          pot=psi0
          if(istar.eq.1)then
            if(iverb.eq.0)then
              write(2,999)Q,1.0d0/Q
              write(2,1000)Rl,Rl*fill,reff,fill*rocheradius(1.0d0/Q)
              write(2,1002)reff*separation
              write(2,1005)rocheradius(1.0d0/Q)
              write(2,1006)Rpol(Nalf2)
              write(2,1004)xend(3),xend(4)
              write(2,1007)savepsi0,psi0
              write(3,2000)reff
              write(3,2001)Rpol(Nalf2)
            endif
            Regg=rocheradius(1.0d0/Q)
c
c   NEW BUG ALERT  July 13, 2001
c
c   Change the indices of xend to 3,4
c
            if(iverb.eq.0)then
              write(3,2002)xend(3)
              write(3,2003)xend(4)
              write(3,2004)psi0
              write(3,2005)sarea
              write(3,2006)vol
              write(3,2007)Tpole
            endif
          endif
          if(istar.eq.2)then
c
            psi0=psi0/overQ+0.5d0*(overQ-1.0d0)/overQ
            savepsi0=savepsi0/overQ+0.5d0*(overQ-1.0d0)/overQ
            if(iverb.eq.0)then
              write(2,1001)Rl,Rl*fill,reff,fill*rocheradius(Q)
              write(2,1002)reff*separation
              write(2,1005)rocheradius(Q)
              write(2,1006)Rpol(Nalf2)
              write(2,1004)xend(3),xend(4)
              write(2,1007)savepsi0,psi0
              write(3,3000)reff
              write(3,3001)Rpol(Nalf2)
            endif
            Regg=rocheradius(Q)
c
c   NEW BUG ALERT  July 13, 2001
c
c   Change the indices of xend to 3,4
c
            if(iverb.eq.0)then
              write(3,3002)xend(3)
              write(3,3003)xend(4)
              write(3,3004)psi0
              write(3,3005)sarea
              write(3,3006)vol
              write(3,3007)Tpole
            endif
          endif
c
c   Make the xend array contain the poles of the star for compatability
c   with the plotting subroutine
c
          xend(1)=zpole
          xend(2)=-zpole
c         
c   UPDATE September 11, 2001
c
c   Put the if-then clause around the write statement.
c
          if(iverb.eq.0)write(2,1003)sarea,vol,Teff,Tpole
c
 999      format(/'Q = ',f8.5,'  1/Q = ',f15.7)
 1000     format(/'star 1:'/'L1 = ',f6.4,',  fill*L1 = ',f6.4,
     &     ', r_eff = ',f9.7,
     $     ', r_eff (Eggleton) = ',f9.7)
 1001     format(/'star 2:'/'L1 = ',f6.4,',  fill*L1 = ',f6.4,
     &     ', r_eff = ',f9.7,
     $     ', r_eff (Eggleton) = ',f9.7)
 1002     format('effective radius in solar units = ',f9.5)
c
 1003     format('surface area = ',e16.9,2x,'volume = ',e16.9,2x,/
     %     'effective temperature = ',f9.3,2x,'polar temperature = ',f9.3)
 1004     format('x(point) = ',f10.8,2x,'x(end) = ',f10.8)
 1005     format('Roche lobe effective dimensionless radius ',
     &      '(Eggleton formula) = ',f7.5)
 1006     format('polar radius = ',f10.8)
 1007     format('potential at L1 = ',f12.6,2x,'potential at fill*L1 = ',
     %      f12.6)
c
 2000     format( f14.11,11x,'r_eff (star 1)')
 3000     format( f14.11,11x,'r_eff (star 2)')
 2001     format( f14.11,11x,'r_pole (star 1)')
 3001     format( f14.11,11x,'r_pole (star 2)')
 2002     format( f14.11,11x,'x(point) (star 1)')
 3002     format( f14.11,11x,'x(point) (star 2)')
 2003     format( f14.11,11x,'x(end) (star 1)')
 3003     format( f14.11,11x,'x(end) (star 2)')
 2004     format(f17.10, 8x,'potential at fill1*L1')     
 3004     format(f17.10, 8x,'potential at fill2*L1')
 2005     format(e16.9, 9x,'surface area (star 1)')
 3005     format(e16.9, 9x,'surface area (star 2)')
 2006     format(e16.9, 9x,'volume (star 1)')
 3006     format(e16.9, 9x,'volume (star 2)')
 2007     format(f15.9,10x,'polar temperature (star 1)')
 3007     format(f15.9,10x,'polar temperature (star 2)')
c

          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          double precision function rocheradius(q)
c
c      Finds the radius of the Roche lobe around mass 1 as given in 
c      Pringle and Wade where q = M1/M2
c
          implicit double precision(a-h,o-z)
c
          q13=q**(1.0d0/3.0d0)
          q23=q13*q13
          t1=dlog(1.0d0+q13)
          rocheradius=0.49d0*q23/(0.6d0*q23+t1)
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine findL1(Q,omega,x0,iflag,bdist,tidephi,itide)
c
c   October 6, 1999
c
c   This routine returns the distance to L1 when given Q.  This is more
c   or less stolen from the Avni code.
c
c   February 15, 2000
c
c   The iflag is added because the potential computation needs to be
c   modified in some cases, specifically when star 2 is not synchronous.
c
          implicit double precision(a-h,o-z)

          x=0.2d0*bdist 
          x=0.99d0*bdist 
c          if(Q.lt.1.0d-3)x=0.95d0*bdist
          if(bdist.gt.1.7d0)x=0.450d0*bdist
c          write(*,*)'x initial = ',x,Q
          y=0.0d0
          z=0.0d0
c
c   UPDATE November 14, 2009
c
c   Make the modifications to use the tidal potential option.
c   cox=1.0 and coy=0.0
c
          cox=1.0d0
          coy=0.0d0
c
c   UPDATE MARCH 5, 2008
C
c   Make the number of loops 40, to ensure convergence for extreme values
c   of the mass ratio or omega.
c
          do 10 i=1,60
c            write (*,*)i,x,dx   
            call POTEN(Q,omega,x,y,z,psi,psix,psixx,psiy,psiz,iflag,
     #         bdist,cox,coy,tidephi,1)
            dx=-psix/psixx
            if(dabs(dx).lt.1.0d-16)go to 15
            x=x+dx
            dx=dabs(dx)
 10       continue
 15       x0=x           ! here is the distance to L1
c
c          write(*,*)bdist,x0
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          SUBROUTINE RAD(Q,omega,cox,coy,coz,psi0,r,x,y,z,iflag,bdist,
     $         tidephi,itide)
c
c   October 6, 1999
c
c   This routine finds the radius of a given point on the Roche lobe.
c   From Avni's original code.
c
c   February 15, 2000
c
c   The iflag is added because the potential computation needs to be
c   modified in some cases, specifically when star 2 is not synchronous.
c
c
          implicit double precision(a-h,o-z)
c
          do 10 i=1,190
            X=R*COX
            Y=R*COY
            Z=R*COZ
c
c            CALL POTEN(Q,omega,x,y,z,psi,psix,psixx,psiy,psiz,iflag,bdist)
c            PSIR=PSIX*COX + PSIY*COY + PSIZ*COZ
c            DR=-(PSI-PSI0)/PSIR
c            IF(dabs(dr).lt.1.0e-19) GO TO 15
c            R=R+DR
c            DR=DABS(DR)
c
c
c  UPDATE March 22, 2002
c
c  Remove the variable coy from the argument list of spherepot.
c

             call spherepot(Q,omega,cox,coz,r,psi,dpsidr,bdist,tidephi,itide)
             rnew=r-(psi-psi0)/dpsidr
             dr=dabs(rnew-r)
c             IF(dabs(dr).lt.1.0d-19) GO TO 15
             IF(dabs(dr).lt.1.0d-16) GO TO 15
             r=rnew
 10       continue
c
c   if we made to this point, the radius did not converge.  Try bisection
c   instead.
c  
          rsave=r
          call findL1(Q,omega,x0,iflag,bdist,tidephi,itide)
c
c   The radius needed is between 0.0 and x0      
c
          aa=1.0d-15
          bb=x0
c
          do 30 i=1,90
            call spherepot(Q,omega,cox,coz,aa,psi,dpsidr,bdist,tidephi,itide)
            aapsi=psi-psi0     !should be positive
            call spherepot(Q,omega,cox,coz,bb,psi,dpsidr,bdist,tidephi,itide)
            bbpsi=psi-psi0     !should be negative
            rmid=0.50*(aa+bb)
            call spherepot(Q,omega,cox,coz,rmid,psi,dpsidr,bdist,tidephi,itide)
            rmidpsi=psi-psi0
            if(aapsi.eq.0.0d0)then
              r=aa
              go to 31
            endif
            if(bbpsi.eq.0.0d0)then
              r=bb
              go to 31
            endif            
            if(aapsi*rmidpsi.gt.0.0d0)then
              aa=rmid
            else
              bb=rmid
            endif
c            write(*,*)aapsi,bbpsi,rmid
30        continue
          r=rmid
31        X=R*COX
          Y=R*COY
          Z=R*COZ
c
c          write(*,*)'%%%%%%%% ',r-rsave

 15       CALL POTEN(Q,omega,x,y,z,psi,psix,psixx,psiy,psiz,iflag,bdist,
     #         cox,coy,tidephi,itide)


          RETURN
          END
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  RVG BUG ALERT  May 9, 2001
c
c  Add the spot parameters to this subroutine
c
c  UPDATE August 10, 2004
c
c  Add the 8 new real and 4 new integer variables to the list.
c
c  UPDATE May 8, 2006
c
c  Add isw21-isw24, sw21-sw24, powercoeff to list.
c
          subroutine recordparm(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     &       Ntheta,Nradius,alb1,alb2,Nref,
     %       rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
     $       idraw,iecheck,idint,iatm,ism1,
     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,
     &       isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey,
     $       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     $       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
     $       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       idark1,idark2,isw12,isw13,
     #       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,
     #       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2,
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)

c
c    Will record the parameters used in the file ELC.out, and will also
c    record interesting computed parameters.
c
c
          implicit double precision(a-h,o-z)
c
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
          dimension powercoeff(8,9)
c
c   RVG BUG ALERT  May 9, 2001
c
c   Dimension the spot arrays.
c
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
c
          open(unit=2,file='ELC.out',status='unknown')
c
          write(2,*)'This is ELC version 3.0 (August 13, 2004)'
          write(2,1000)Nalph1
          write(2,1001)Nbet1
          write(2,1002)Nalph2
          write(2,1003)Nbet2
          write(2,1004)fill1
          write(2,1005)fill2
          write(2,1006)omega1
          write(2,1007)omega2
          write(2,1008)dphase
          write(2,1009)Q
          write(2,1010)finc
          write(2,1011)Teff1
          write(2,1012)Teff2
          write(2,1013)Tgrav1
          write(2,1014)Tgrav2
          write(2,1015)betarim
          write(2,1016)rinner
          write(2,1017)router
          write(2,1018)tdisk
          write(2,2018)xi
          write(2,1019)Ntheta
          write(2,1020)Nradius
          write(2,1042)alb1
          write(2,1043)alb2
          write(2,2043)Nref
          write(2,1021)rLx
          write(2,1023)Period
          if(fm.gt.1.0d-4)then
            write(2,1024)fm
          else
            write(2,9024)fm
          endif
 9024     format(1pe16.9,4x,'fm')
          write(2,1025)separ
          write(2,4025)gamma
          write(2,5000)t3
          write(2,5001)g3
          write(2,5002)SA3
          write(2,5003)density
          write(2,5004)sw1
          write(2,5005)sw2
          write(2,5006)sw3
          write(2,5007)T0
          write(2,1026)idraw
          write(2,1027)iecheck
          write(2,1028)idint
          write(2,4000)iatm
          write(2,4001)ism1
c
c  RVG BUG ALERT   May 16, 2001
c
c  Modify the icn? flags so that they are either 0 or 1, which conforms
c  to the output format statement.
c
c   Set the icn? control numbers back to zeros and 1s.
c
          if(icnU.eq.430)then
            iU=0
          else
            iU=1
          endif
          if(icnB.eq.430)then
            iB=0
          else
            iB=1
          endif
          if(icnV.eq.430)then
            iV=0
          else
            iV=1
          endif
          if(icnR.eq.430)then
            iR=0
          else
            iR=1
          endif
          if(icnI.eq.430)then
            iI=0
          else
            iI=1
          endif
          if(icnJ.eq.430)then
            iJ=0
          else
            iJ=1
          endif
          if(icnH.eq.430)then
            iH=0
          else
            iH=1
          endif
          if(icnK.eq.430)then
            iK=0
          else
            iK=1
          endif
c
          write(2,4002)iU,iB,iV,iR,iI,iJ,iH,iK
c
c  END BUG
c
          write(2,5008)iRVfilt
          write(2,5009)isw1
          write(2,5010)isw2
          write(2,5011)isw3
          write(2,5012)isw4
          write(2,3028)ilaw
c
          do 10 i=1,8
            write(2,2000)wave(i),dbolx(i,1),dboly(i,1),dbolx(i,2),dboly(i,2),
     $           dwavex(i,1),dwavey(i,1),dwavex(i,2),dwavey(i,2)
 10       continue
c
          write(2,6000)ecc
          write(2,6001)argper
          write(2,6002)pshift
          write(2,6003)sw5
          write(2,6004)sw6
          write(2,6005)sw7
          write(2,6006)sw8
          write(2,6007)sw9
          write(2,6008)ikeep
          write(2,6009)isynch
          write(2,6010)isw5
          write(2,6011)isw6
          write(2,6012)isw7
          write(2,6013)isw8
          write(2,6014)isw9
c
          write(2,2001)spot1parm(1,1)
          write(2,2002)spot1parm(1,2)
          write(2,2003)spot1parm(1,3)
          write(2,2004)spot1parm(1,4)
c                     
          write(2,2005)spot1parm(2,1)
          write(2,2006)spot1parm(2,2)
          write(2,2007)spot1parm(2,3)
          write(2,2008)spot1parm(2,4)
c                     
c                     
          write(2,2009)spot2parm(1,1)
          write(2,2010)spot2parm(1,2)
          write(2,2011)spot2parm(1,3)
          write(2,2012)spot2parm(1,4)
                      
          write(2,2013)spot2parm(2,1)
          write(2,2014)spot2parm(2,2)
          write(2,2015)spot2parm(2,3)
          write(2,2016)spot2parm(2,4)
c                     
c                     
          write(2,2017)spotdparm(1,1)
          write(2,22018)spotdparm(1,2)
          write(2,2019)spotdparm(1,3)
          write(2,2020)spotdparm(1,4)
                      
          write(2,2021)spotdparm(2,1)
          write(2,2022)spotdparm(2,2)
          write(2,2023)spotdparm(2,3)
          write(2,2024)spotdparm(2,4)
c
c   UPDATE August 10, 2004
c
c   Write out the 8 real and 4 new integer variables here.
c
          write(2,2025)primmass
          write(2,2026)primK
          write(2,2027)primrad
          write(2,2028)ratrad
c
          write(2,2030)frac1
          write(2,2031)frac2
          write(2,2032)ecosw
          write(2,2033)temprat
c
          write(2,2040)idark1
          write(2,2041)idark2
          write(2,2042)isw12
          write(2,3043)isw13
c
          write(2,8001)isw21
          write(2,8002)isw22
          write(2,8003)isw23
          write(2,8004)isw24
c
8001      format(i1,19x,'ialign (0 for rotation aligned with orbit)')
8002      format(i1,19x,'ifastgen (1 for fast genetic mode)')
8003      format(i1,19x,'isw23 (currently inactive)')         
8004      format(i1,19x,'frac switch (>1 to enable ELCratio.???? files)')         
c
          do 85000 kk=1,8
              write(2,85001)powercoeff(kk,1),powercoeff(kk,2),
     $       powercoeff(kk,3),powercoeff(kk,4),powercoeff(kk,5),
     $       powercoeff(kk,6),powercoeff(kk,7),powercoeff(kk,8),
     #       powercoeff(kk,9)
85000     continue
c
85001      format(9(f7.4,1x))
c
          write(2,8011)bigI
          write(2,8012)bigbeta
          write(2,8013)sw23
          write(2,8014)sw24
c
8011      format(f11.7,9x,'axis_I (inclination of rotation axis if ialign=1)')
8012      format(f11.7,9x,'axis_beta (angle of rotation axis wrt to orbit if',
     #      ' ialign=1)')          
8013      format(f15.7,5x,'t_start')
8014      format(f15.7,5x,'t_end')         
c
c
c  UPDATE November 6, 2008
c
c  write the new variables sw25-sw34 and isw25-isw34 here
c
          write(2,8025)sw25
          write(2,8026)sw26
          write(2,8027)sw27
          write(2,8028)sw28
          write(2,8029)sw29
          write(2,8030)sw30
          write(2,8031)sw31
          write(2,8032)Tconj
          write(2,8033)beam1
          write(2,8034)beam2

          write(2,9025)isw25
          write(2,9026)isw26
          write(2,9027)isw27
          write(2,9028)isw28
          write(2,9029)isw29
          write(2,9030)isw30
          write(2,9031)isw31
          write(2,9032)isw32
          write(2,9033)isw33
          write(2,9034)isw34
8025      format(f10.7,19x,'asini error')
8026      format(f10.8,10x,'reference phase for disk fraction')
8027      format(f10.8,10x,'radfill1 (set to use fill1 in terms of R_eff')
8028      format(f10.8,10x,'radfill2 (set to use fill2 in terms of R_eff')
8029      format(f10.4,10x,'bin size for light curves (minutes)')
8030      format(f10.4,10x,'bin size for RV curves (minutes)')
8031      format(f4.2,16x,'sw31 (currently inactive)')
8032      format(f15.8,5x,'Tconj')
8033      format(f4.2,16x,'beam1 (Doppler boost factor, star 1)')
8034      format(f4.2,16x,'beam2 (Doppler boost factor, star 2)')

9025      format(i1,19x,'isw25 (currently inactive)')
9026      format(i1,19x,'isw26 (currently inactive)')
9027      format(i6,14x,'Nterms for fast analytic')
9028      format(i1,19x,'set to 1 to fit for Tconj')
9029      format(i1,19x,'isw29 (currently inactive)')
9030      format(i1,19x,'isw30 (currently inactive)')
9031      format(i1,19x,'isw31 (currently inactive)')
9032      format(i1,19x,'isw32 (currently inactive)')
9033      format(i1,19x,'isw33 (currently inactive)')
9034      format(i1,19x,'isw34 (currently inactive)')


c
 2025     format(f13.9,7x,'primmass (star 1 mass in solar masses)')
 2026     format(f14.9,6x,'primK (K-velocity of star 1 in km/sec)')
 2027     format(f14.9,6x,'primrad (star 1 radius in solar radii)')
 2028     format(f16.9,4x,
     &          'ratrad (ratio of star 1 radius and star 2 radius)')
c
 2030     format(f4.2,16x,'frac1 (fractional radius star 1: R_1/a)')
 2031     format(f4.2,16x,'frac2 (fractional radius star 2: R_2/a)')
 2032     format(f12.9,8x,'ecosw (phase difference between eclipses)')
 2033     format(f10.7,10x,'temprat (T_2/T_1)')
c
 2040     format(i1,19x,'idark1')
 2041     format(i1,19x,'idark2')
 2042     format(i6,14x,'Npoly (0 for numerical)')
 3043     format(i1,19x,'ifasttrans (>0 for fast transit mode)')


 1000     format(i4,16x,'Nalph1')
 1001     format(i3,17x,'Nbet1')
 1002     format(i4,16x,'Nalph2')
 1003     format(i3,17x,'Nbet2')
 1004     format(f11.9,9x,'fill1')
 1005     format(f11.9,9x,'fill2')
 1006     format(f10.6,10x,'omega1')
 1007     format(f10.6,10x,'omega2')
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the format statement for line 1008.
c
 1008     format(f11.6,9x,'dphase')
 1009     format(f15.10,5x,'Q')
 1010     format(f8.5,12x,'finc')
 1011     format(f9.2,11x,'Teff1')
 1012     format(f9.2,11x,'Teff2')
c
c  RVG BUG ALERT   May 16, 2001
c
c  Change the format statements for Tgrav:
c
 1013     format(f8.6,12x,'Tgrav1')
 1014     format(f8.6,12x,'Tgrav2')
c
 1015     format(f8.5,12x,'betarim')
 1016     format(f8.6,12x,'rinner')
 1017     format(f8.6,12x,'router')
 1018     format(f7.1,13x,'tdisk')
 2018     format(f7.4,13x,'xi')
 1019     format(i3,17x,'Ntheta')
 1020     format(i3,17x,'Nradius')
c
c   UPDATE March 26, 2002
c
c   The meaning of rLx is now log10(Lx)
c
c 1021     format(f10.5,10x,'Lx/Lopt')
 1021     format(f10.5,10x,'log10(Lx)')
 1022     format(f4.2,16x,'W')
 1023     format(f16.10,4x,'Period')
 1024     format(f8.5,12x,'fm')
 1025     format(f14.5,6x,'separ')
 1026     format(i1,19x,'idraw')
 1027     format(i1,19x,'iecheck')
 1028     format(i1,19x,'idint')
 1029     format(i1,19x,'ivelout')
 1030     format(i1,19x,'iXout')
 1040     format(f6.4,14x,'darkbol1')
 1041     format(f6.4,14x,'darkbol2')
 1042     format(f6.4,14x,'alb1')
 1043     format(f6.4,14x,'alb2')
c
c   UPDATE April 15, 2002
c
c   Change the format statement below to i2,18x
c
 2043     format(i2,18x,'Nref')
 2000     format(f7.1,3x,8(f7.4,1x))
 3028     format(i2,18x,'ilaw  (1=linear law, 2=logarithmic law,',
     %           ' 3=square root law, 4=quad law, >10 for power series)')
 4000     format(i1,19x,'iatm')
 4001     format(i1,19x,'ism1')
 4002     format(8(i1,1x),4x,'icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK')
 4025     format(f11.6,9x,'gamma velocity')
 5000     format(f10.2,10x,'t3')
 5001     format(f10.5,10x,'g3')
 5002     format(f12.6,8x,'SA3')
5003      format(f12.6,8x,'density in g/cc')
 5004     format(f12.6,8x,'onephase')
 5005     format(f12.6,8x,'usepot1')
 5006     format(f12.6,8x,'usepot2')
 5007     format(f16.9,4x,'T0 ')
 5008     format(i1,19x,'iRVfilt')
 5009     format(i1,19x,'ionephase')
 5010     format(i1,19x,'isquare')
 5011     format(i1,19x,'iusepot')
 5012     format(i1,19x,'ifixgamma (currently inactive)')
 6000     format(f10.8,10x,'eccentricity')
 6001     format(f13.8,7x,'argument of peristron in degrees')
 6002     format(f11.8,9x,'pshift')
c
c   UPDATE April 8, 2002
c
c   Change format statement below.
c
 6003     format(f12.6,8x,'asini (projected semimajor axis in seconds)')
c
c   UPDATE May 27, 2002
c
c   Change this format statement
c
 6004     format(f12.6,8x,'median fit (geneticELC only)')
 6005     format(f12.6,8x,'sw7 (currently inactive)')
 6006     format(f12.6,8x,'sw8 (currently inactive)')
 6007     format(f12.6,8x,'sw9 (currently inactive)')
 6008     format(i1,19x,'ikeep (1 to put eclipse at phase 0.0)')
 6009     format(i1,19x,'isynch (1 to keep rotation synchronous',
     $         ' at periastron)')
 6010     format(i1,19x,'ispotprof')
 6011     format(i1,19x,'igrav')
 6012     format(i1,19x,'itime')
c
c  UPDATE JULY 4, 2004
c
c  The variable isw8 will be assigned to MonteCarlo, which
c  will be used to determine the number of Monte Carlo
c  integrations done on the pixels.
c
 6013     format(i6,14x,'MonteCarlo (0 for interpolation, >10 ',
     $        'for Monte Carlo)')
c
 6014     format(i1,19x,'ielite')
c
c  RVG BUG ALERT  May 9, 2001
c
c  Add these format statements.
c
 2001     format(f10.7,10x,'Temperature factor spot 1, star 1')
 2002     format(f11.7, 9x,'Latitude of spot 1, star 1 (degrees)') 
 2003     format(f11.7, 9x,'Longitude of spot 1, star 1 (degrees)') 
 2004     format(f11.7, 9x,'Angular radius of spot 1, star 1 (degrees)') 

 2005     format(f10.7,10x,'Temperature factor spot 2, star 1')
 2006     format(f11.7, 9x,'Latitude of spot 2, star 1 (degrees)') 
 2007     format(f11.7, 9x,'Longitude of spot 2, star 1 (degrees)') 
 2008     format(f11.7, 9x,'Angular radius of spot 2, star 1 (degrees)') 

 2009     format(f10.7,10x,'Temperature factor spot 1, star 2')
 2010     format(f11.7, 9x,'Latitude of spot 1, star 2 (degrees)') 
 2011     format(f11.7, 9x,'Longitude of spot 1, star 2 (degrees)') 
 2012     format(f11.7, 9x,'Angular radius of spot 1, star 2 (degrees)') 

 2013     format(f10.7,10x,'Temperature factor spot 2, star 2')
 2014     format(f11.7, 9x,'Latitude of spot 2, star 2 (degrees)') 
 2015     format(f11.7, 9x,'Longitude of spot 2, star 2 (degrees)') 
 2016     format(f11.7, 9x,'Angular radius of spot 2, star 2 (degrees)') 

 2017     format(f10.7,10x,'Temperature factor spot 1, disk')
22018     format(f11.7, 9x,'Azimuth of spot 1, disk (degrees)') 
 2019     format(f11.7, 9x,'Radial cutoff of spot 1, disk (0 <= r_cut <=1)') 
 2020     format(f11.7, 9x,'Angular size of spot 1, disk (degrees)') 

 2021     format(f10.7,10x,'Temperature factor spot 2, disk')
 2022     format(f11.7, 9x,'Azimuth of spot 2, disk (degrees)') 
 2023     format(f11.7, 9x,'Radial cutoff of spot 2, disk (0 <= r_cut <=1)') 
 2024     format(f11.7, 9x,'Angular size of spot 2, disk (degrees)') 
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine setuptemp(istar,ialphmax,ibetmax,
     $       Nalf,Nbet,ibetlim,gmatrix,Tpole,Tgrav,tmatrix,gpole,mmdx)
c
c   This routine will assign the temperatures of the grid points
c   based on the polar temperature and the gravity darkening law.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265359879323d0)
          dimension gmatrix(ialphmax*ibetmax),tmatrix(ialphmax*ibetmax),
     %       ibetlim(ialphmax),mmdx(ialphmax,ibetmax)
c
          tmax=-123456.0d0
          tmin=123456.0d0
          gmin=123456.0d0
          gmax=-123456.0d0
          NALF2 = NALF/2
          DCOAL = 2.0d0/NALF              ! step size in alpha (latitude)
          DBETA = (pie/2.0d0)/NBET        ! step size in longitude
          DIV = gpole                  ! gravity at the pole
c
          DO 10 IALF = 1, NALF
            DO 9 IBET = 1, ibetlim(ialf)    !4*NBET
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              G_div = gmatrix(iidx)/DIV
              tmatrix(iidx)=((G_div**Tgrav*Tpole))
c
c
c
c
c
c  
              if(tmatrix(iidx).lt.tmin)then
                tmin=tmatrix(iidx)
                iamin=ialf
                ibmin=ibet
              endif
              if(tmatrix(iidx).gt.tmax)then
                tmax=tmatrix(iidx)
                iamax=ialf
                ibmax=ibet
              endif
              if(gmatrix(iidx).lt.gmin)then
                gmin=gmatrix(iidx)
                iagmin=ialf
                ibgmin=ibet
              endif
              if(gmatrix(iidx).gt.gmax)then
                gmax=gmatrix(iidx)
                iagmax=ialf
                ibgmax=ibet
              endif
 9          continue
 10       continue
c
          write(2,100)istar,tmin,iamin,ibmin,tmax,iamax,ibmax
          write(2,101)istar,gmin,iagmin,ibgmin,gmax,iagmax,ibgmax
c
 100      format(/'star ',i1,':  min temp = ',f11.1,'     (at ialf = ',i4,
     %       ' ibet = ',i4,')',/'         max temp = ',f11.1,
     *       '     (at ialf = ',i4,
     %       ' ibet = ',i4,')')
 101      format('star ',i1,':  min grav = ',f11.1,'     (at ialf = ',i4,
     %       ' ibet = ',i4,', program units)',
     &                     /'         max grav = ',f11.1,'     (at ialf = ',i4,
     %       ' ibet = ',i4,', program units)')
c
          return
          end
c
c  ************************
c
          subroutine getcoords(istar,ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     $      phase,finc,Q,
     $      xarray,yarray,zarray,gradx,grady,gradz,temp,
     $      Ncoords,xcoords,ycoords,xend,extension,separation,bdist,
     &      mmdx)
c
c  October 9, 1999
c
c  This subroutine will return the sky coordinates of the star in the array
c  xcoords(1:Ncoords),ycoords(1:Ncoords).  Only the points that are visible
c  to the observer are included, and eclipses are not accounted for.
c
c  Set istar=1 to do star 1, istar=2 to do star2
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xarray(ialphmax*ibetmax),yarray(ialphmax*ibetmax),
     $        zarray(ialphmax*ibetmax),xcoords(ialphmax*ibetmax*4),
     $        gradx(ialphmax*ibetmax),grady(ialphmax*ibetmax),
     $        gradz(ialphmax*ibetmax),ycoords(ialphmax*ibetmax*4),xend(4),
     $        temp(ialphmax*ibetmax),ibetlim(ialphmax),
     #        mmdx(ialphmax,ibetmax)
c
          character*9 extension
c
          if(istar.eq.1)open(unit=39,file='star1temp.'//extension,
     %               status='unknown')
          if(istar.eq.2)open(unit=39,file='star2temp.'//extension,
     %               status='unknown')

c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c          
          NBET4 = NBET*4
          DBETA = (pie/2.0d0)/NBET        ! step size in longitude
          AZ = DCOS(FINCR)
          IF (AZ.LT.0.0d0) AZ = 0.0d0
          AX = -DSIN(FINCR)*DCOS(PHASER)    ! l in Wilson & Sofia
          AY = DSIN(FINCR)*DSIN(PHASER)     ! m in Wilson & Sofia
          A2 = DACOS(AX)
          A3 = DSIN(A2)
          IF (A3.LT.0.0d0) A3=0.0d0
          IF (A3.EQ.0.0d0) GO TO 508
          B1=AZ/DSIN(A2)
          IF(B1.GT.1.0d0) B1=1.0d0
          BETA = DASIN(B1)              !beta is the angle between the surface
          GO TO 509                    !normal and the radius vector.
508       BETA = 0.0d0
509       KBETA = BETA/DBETA + 0.50d0
          KBETA1 = KBETA + 1
          KBETAN = KBETA + 2*NBET
c
          Ncoords=0
c
          if(istar.eq.1)then
            xx=0.0d0
            yy=0.0d0
            zz=0.0d0
          endif
          if(istar.eq.2)then
            xx=1.0d0
            yy=0.0d0
            zz=0.0d0
          endif
          xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)    ! projected coords
          yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
          Ncoords=Ncoords+1
          xcoords(Ncoords)=xp*separation
          ycoords(Ncoords)=yp*separation


c   Check the visibility of the nose and end
c
          do 10 i=1,2
            if(i.eq.1)proj=AX
            if(i.eq.2)proj=-AX
            if(proj.gt.0.0)then
              xp=xtran(xend(i+2),0.0d0,0.0d0,phase,fincr,Q,istar,bdist) 
              yp=ytran(xend(i+2),0.0d0,0.0d0,phase,fincr,Q,istar,bdist)    
              Ncoords=Ncoords+1
              xcoords(Ncoords)=xp*separation
              ycoords(Ncoords)=yp*separation
            endif
 10       continue
c
          DO 501 IALF = 1, NALF
            DO 502 IBET = 1,ibetlim(ialf)
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              PROJ = AX * GRADX(iidx) + AY*GRADY(iidx) + 
     1	        AZ*GRADZ(iidx)
              IF (PROJ.LT.0.0d0) GO TO 502    ! is the surface element visible?
              xx=xarray(iidx)
              yy=yarray(iidx)
              zz=zarray(iidx)
              xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) ! projected coords
              yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)    
              Ncoords=Ncoords+1
              xcoords(Ncoords)=xp*separation
              ycoords(Ncoords)=yp*separation
c
c   Record the x,y,z coordinates of the nearby points.  These points
c   will be used for area filling.
c
              xx1=xp*separation
              yy1=yp*separation
              if(ibet.gt.1)then
c                iidx=kount(ialphmax,ialf,ibetlim)+(ibet-1)
                iidx=mmdx(ialf,ibet-1)
                xx=xarray(iidx)
                yy=yarray(iidx)
                zz=zarray(iidx)
              else
c                iidx=kount(ialphmax,ialf,ibetlim)+ibetlim(ialf)
                iidx=mmdx(ialf,ibetlim(ialf))
                xx=xarray(iidx)
                yy=yarray(iidx)
                zz=zarray(iidx)
              endif
              xx2=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
              yy2=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
              if(ialf.gt.1)then
                if(ibet.gt.1)then
                  izz=ialf-1
                  jzz=ibet-1
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xx=xarray(iidx)
                  yy=yarray(iidx)
                  zz=zarray(iidx)
                else
                  izz=ialf-1
                  jzz=ibetlim(ialf)
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xx=xarray(iidx)
                  yy=yarray(iidx)
                  zz=zarray(iidx)
                endif
              else
                xx=xend(1)
                yy=0.0d0
                zz=0.0d0
              endif
              xx3=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
              yy3=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
              if(ialf.gt.1)then
               izz=ialf-1
               jzz=ibet
c               iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                xx=xarray(iidx)
                yy=yarray(iidx)
                zz=zarray(iidx)
              else
                xx=xend(1)
                yy=0.0d0
                zz=0.0d0
              endif
              xx4=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
              yy4=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
c
              izz=ialf
              jzz=ibet
c              iidx=kount(ialphmax,izz,ibetlim)+jzz
              iidx=mmdx(izz,jzz)
              write(39,69)temp(iidx),xx1,yy1,xx2,yy2,xx3,yy3,xx4,yy4
              if(ialf.eq.Nalf)then
                if(ibet.gt.1)then
                  izz=ialf
                  jzz=ibet-1
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xx=xarray(iidx)
                  yy=yarray(iidx)
                  zz=zarray(iidx)
                else
                  izz=ialf
                  jzz=ibetlim(ialf)
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xx=xarray(iidx)
                  yy=yarray(iidx)
                  zz=zarray(iidx)
                endif
                xx2=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
                yy2=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
c
                xx=xend(2)
                yy=0.0d0
                zz=0.0d0
                xx3=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
                yy3=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation    
c
                xx=xend(2)
                yy=0.0d0
                zz=0.0d0
                xx4=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
                yy4=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)*separation
c
                izz=ialf
                jzz=ibet
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(39,69)temp(iidx),xx1,yy1,xx2,yy2,xx3,yy3,xx4,yy4
              endif
c
502         CONTINUE
c
 501      continue         ! continue the alpha loop
c
 69       format(f7.1,1x,8(f7.4,2x))
 70       format(e16.9,3x,1x,8(f7.4,1x))
c
          close(39)
          close(40)
c
          return
          end
c
c  %%%%%%%%%%%%%%%%%%%%%
c
          double precision function xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c   will return the coordinate of a point (xx,yy,zz) projected on the sky
c
c   (xx,yy,zz) refers to the coordinates in the rotating system
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q

c          if(istar.eq.2)xx=xx+bdist
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif

          qphaser=((phaser))
c
c          xtran=(-xx*dsin(qphaser)-yy*dcos(qphaser))
c     $      -bdist*dsin(qphaser)*(overQ/(1.0d0+overQ))

          xtran=-(xx*dsin(qphaser)+yy*dcos(qphaser))+
     $      bdist*(overQ/(1.0d0+overQ))*dsin(qphaser)

c
c
c   RVG BUG ALERT   May 2, 2001
c
c   Comment out this if-then statement.
c
c          if(phase.gt.180.0d0)xtran=-xtran
c
c 
c   Added February 8, 2001
c
c          xtran=xtran*bdist
c
          return
          end
c
c  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          double precision function ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c   will return the coordinate of a point (xx,yy,zz) projected on the sky
c
c   (xx,yy,zz) refers to the coordinates in the rotating system
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q

c          if(istar.eq.2)xx=xx+bdist
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians).
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif

          qphaser=((phaser))

c          ytran=(xx*dcos(fincr)*dcos(qphaser)-yy*dcos(fincr)*dsin(qphaser)
c     $      +zz*dsin(fincr))
c     $      -bdist*dcos(fincr)*dcos(qphaser)*(overQ/(1.0d0+overQ))

          ytran=-(-xx*dcos(fincr)*dcos(qphaser)+yy*dcos(fincr)*dsin(qphaser)
     $      -zz*dsin(fincr))+
     $     bdist*(-(overQ/(1.0d0+overQ))*dcos(fincr)*dcos(qphaser))

c
c   RVG BUG ALERT   May 2, 2001
c
c   Comment out this statement.
c
c          if(phase.gt.180.0d0)ytran=-ytran
c
c   Added February 8, 2001
c
c          ytran=ytran*bdist
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine sortcircle(N,xcir,ycir)
c
c   October 7, 1999
c
c   This routine will arrange the (x,y) points of the polygon in order by
c   sorting by theta in polar coordinates.  Written by Orosz circa 1996.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xcir(N),ycir(N),theta(10000)
c
          call getmaxvalues(N,xcir,ycir,xmin,xmax,ymin,ymax)  !get the extreme
c                                                            !values
          xcenter=(xmax+xmin)/2.0d0
          ycenter=(ymax+ymin)/2.0d0
c
c   convert to polar coordinates         
c
          do 10 i=1,N
            xxx=xcir(i)-xcenter
            yyy=ycir(i)-ycenter
            if(xxx.eq.0.0)then
              if(yyy.lt.0.0)then
                theta(i)=3.0d0*pie/2.0d0
                go to 10
              endif
              if(yyy.ge.0.0d0)then
                theta(i)=pie/2.0d0
                go to 10
              endif
            endif
            if((yyy.ge.0.0d0).and.(xxx.ge.0.0d0))then
              theta(i)=(atan2(yyy,xxx))
            endif
            if((yyy.ge.0.0d0).and.(xxx.lt.0.0d0))then
              theta(i)=(atan2(yyy,xxx))
            endif
            if((yyy.lt.0.0d0).and.(xxx.lt.0.0d0))then
              theta(i)=2.0d0*pie+(atan2(yyy,xxx))
            endif
            if((yyy.lt.0.0d0).and.(xxx.ge.0.0d0))then
              theta(i)=2.0d0*pie+atan2(yyy,xxx)
            endif
 10       continue
c
          call sort3(N,theta,xcir,ycir)  !sort by theta and swap x and y also
c
          return
          end
c
c ==========================================================================
c
          subroutine getmaxvalues(N,xdata,ydata,xmin,xmax,ymin,ymax)
c
          implicit double precision(a-h,o-z)
c
          dimension xdata(N),ydata(N)
c
          xmin=100000000000.0d0
          ymin=100000000000.0d0
          xmax=-10000000.0d0
          ymax=-10000000.0d0
c
          do 10 i=1,N
            t1=xdata(i)
            q1=ydata(i)
            if(t1.gt.xmax)xmax=t1
            if(t1.lt.xmin)xmin=t1
            if(q1.gt.ymax)ymax=q1
            if(q1.lt.ymin)ymin=q1
 10       continue
          return
          end
c
c============================================================
c
          SUBROUTINE SORT3(N,RA,RB,rc)
c
c  October 7, 1999
c
c  This routine will sort an array RA, and rearrange the arrays rb and rc.
c  Taken from Numerical Recipes.
c
          implicit double precision(a-h,o-z)

          DIMENSION RA(N),RB(N),rc(N)
c
          L=N/2+1
          IR=N
10        CONTINUE
          IF(L.GT.1)THEN
            L=L-1
            RRA=RA(L)
            RRB=RB(L)
            rrc=rc(L)
          ELSE
            RRA=RA(IR)
            RRB=RB(IR)
            rrc=rc(IR)
            RA(IR)=RA(1)
            RB(IR)=RB(1)
            rc(IR)=RC(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
              RA(1)=RRA
              RB(1)=RRB
              rc(1)=rrc
              RETURN
            ENDIF
          ENDIF
          I=L
          J=L+L
20        IF(J.LE.IR)THEN
            IF(J.LT.IR)THEN
              IF(RA(J).LT.RA(J+1))J=J+1
            ENDIF
            IF(RRA.LT.RA(J))THEN
              RA(I)=RA(J)
              RB(I)=RB(J)
              rc(i)=rc(j)
              I=J
              J=J+J
            ELSE
              J=IR+1
            ENDIF
            GO TO 20
          ENDIF
          RA(I)=RRA
          RB(I)=RRB
          rc(i)=rrc
          GO TO 10
          END
c
c   ################################################
c
          subroutine dump(iunit,N,x,y)
c
          implicit double precision(a-h,o-z)
c
          dimension x(N),y(N)
          character*2 ext
c
          write(ext,100)iunit
          open(unit=iunit,file='fort.'//ext,status='unknown')
          do 10 i=1,N
            write(iunit,200)x(i),y(i)
 10       continue
c
          close(iunit)
 100      format(i2)
 200      format(f7.5,3x,f10.8)
          return
          end
c
c   ################################################
c 
          subroutine gethorizon(istar,ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     &      phase,finc,Q,
     $      psi0,omega,xarray,yarray,zarray,radarray,gradx,grady,gradz,xend,
     $      Nhoriz,xhoriz,yhoriz,phiar,iedgestar,delphiedge,bdist,mmdx,
     #      xhmin,xhmax,yhmin,yhmax,tidephi,itide,phihor)
c
c  October 9, 1999
c
c  This subroutine will return the sky coordinates of the horizon of the star
c  xhoriz(1:Nhoriz),yhoriz(1:Nhoriz).  For phases near conjunction, go along
c  the alpha direction and record when the "projection factor" turns negative.
c  For phases near quadrature, so a similar loop in the beta direction.
c
c  Set istar=1 to do star 1, istar=2 to do star2
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
c   UPDATE JULY 2, 2004
c
c   Add this parameter for the temparary arrays.
c
          parameter(itemp=10000)
          dimension xarray(ialphmax*ibetmax),yarray(ialphmax*ibetmax),
     $        zarray(ialphmax*ibetmax),xhoriz(4*ibetmax),yhoriz(4*ibetmax),
     $        gradx(ialphmax*ibetmax),grady(ialphmax*ibetmax),xend(4),
     $        gradz(ialphmax*ibetmax),arrproj(itemp),
     &        xdummy(itemp),ydummy(itemp),
     $        xfirstring(itemp),yfirstring(itemp),phiar(ialphmax*ibetmax),
     $        xlastring(itemp),ylastring(itemp),radarray(ialphmax*ibetmax),
     $        ibetlim(ialphmax),iedgestar(ialphmax*ibetmax),
     #        delphiedge(ialphmax*ibetmax),savex(2),savey(2)
          dimension mmdx(ialphmax,ibetmax),phihor(ialphmax,2)
 
c
c  initialize
c
          do 1 i=1,itemp
            xdummy(i)=0.0d0
            ydummy(i)=0.0d0
 1        continue
c
          do 3 ialf=1,nalf
            phihor(ialf,1)=-99.00d0
            phihor(ialf,2)=-99.00d0
            do 2 ibet=1,ibetlim(ialf)   !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              iedgestar(iidx)=0
              delphiedge(iidx)=-99.99d0
 2          continue
 3        continue
c
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c          
          NBET4 = NBET*4
          DBETA = (pie/2.0d0)/NBET        ! step size in longitude
          AZ = DCOS(FINCR)
          IF (AZ.LT.0.0d0) AZ = 0.0d0
          AX = -DSIN(FINCR)*DCOS(PHASER)    ! l in Wilson & Sofia
          AY = DSIN(FINCR)*DSIN(PHASER)     ! m in Wilson & Sofia
cc
          Ndummy=0
          Nfirst=0
          Nlast=0
c
c   Loop over alpha first.  
c
          do 950 ialf=1,nalf       
            do 949 ibet=1,ibetlim(ialf)     !nbet4   
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              arrproj(ibet)=AX*GRADX(iidx) + AY*GRADY(iidx) + 
     &	            AZ*GRADZ(iidx)               
 949        continue
c
c   We have the array of projection factors along a given direction.
c   Now find out where the sign change is.  This is where the line of sight
c   has moved over the horizon.
c
            ihorcount=1
            do 940 ibet=1,ibetlim(ialf)    !nbet4-1
              if(ibet.lt.ibetlim(ialf))then
                index=ibet+1
              else
                index=1
              endif
              rsign=arrproj(ibet)*arrproj(index)
              if(rsign.le.0.0d0)then     ! crossed over
                if(arrproj(ibet).gt.0.0d0)then
                  izz=ialf
                  jzz=ibet
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xvis=xarray(iidx)
                  yvis=yarray(iidx)
                  zvis=zarray(iidx)
                  rvis=radarray(iidx)
                  phivis=phiar(iidx)
                  iedgestar(iidx)=10
                  izz=ialf
                  jzz=index
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xhid=xarray(iidx)
                  yhid=yarray(iidx)
                  zhid=zarray(iidx)
                  rhid=radarray(iidx)
                  phihid=phiar(iidx)
                  iedgestar(iidx)=-10
                  index=ibet
                else
                  izz=ialf
                  jzz=index
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xvis=xarray(iidx)
                  yvis=yarray(iidx)
                  zvis=zarray(iidx)
                  rvis=radarray(iidx)
                  phivis=phiar(iidx)
                  iedgestar(iidx)=20     ! just appeared at limb
c
                  izz=ialf
                  jzz=ibet
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xhid=xarray(iidx)
                  yhid=yarray(iidx)
                  zhid=zarray(iidx)
                  rhid=radarray(iidx)
                  phihid=phiar(iidx)
                  iedgestar(iidx)=-10     ! just went behind limb
                endif
c
                call acchor(overQ,psi0,omega,
     $              xvis,yvis,zvis,rvis,phivis,xhid,yhid,zhid,rhid,
     $              phihid,ax,ay,az,xacc,yacc,zacc,bdist,tidephi,itide)
c
                xx=xacc  
                yy=yacc  
                zz=zacc  
c
                phihor(ialf,ihorcount)=phivis
                ihorcount=ihorcount+1
                izz=ialf
                jzz=ibet
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                diff=(phiar(iidx)-phivis)
                if(dabs(diff).lt.2.0d0)then
                  izz=ialf
                  jzz=index
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  delphiedge(iidx)=dabs(phiar(iidx)-phivis)
                endif
                if(dabs(diff).ge.2.0d0)then
                  izz=ialf
                  jzz=index
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  if(phiar(iidx).gt.4.0)then
                    phitemp=2.0*pie-phiar(iidx)
                    delphiedge(iidx)=(phitemp+phivis)
                  endif
                  if(phivis.gt.4.0d0)then
                    phitemp=2.0*pie-phivis
                    izz=ialf
                    jzz=index
c                    iidx=kount(ialphmax,izz,ibetlim)+jzz
                    iidx=mmdx(izz,jzz)
                    delphiedge(iidx)=(phitemp+phiar(iidx))
                  endif
                endif
c
                xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)  
                yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)    
                Ndummy=Ndummy+1
                xdummy(Ndummy)=xp
                ydummy(Ndummy)=yp
              endif
 940        continue
 950      continue
c
c   Now put beta on the outside loop and look for points when the line
c   of sight passes over the horizon.
c
          kk=1
          do 1950 ibet=1,ibetlim(kk)   !nbet4
            do 1949 ialf=kk,nalf
              izz=ialf
              jzz=ibet
c              iidx=kount(ialphmax,izz,ibetlim)+jzz  ! was (kk)
               iidx=mmdx(izz,jzz)
              arrproj(ialf)=AX*GRADX(iidx) + AY*GRADY(iidx) + 
     1	        AZ*GRADZ(iidx)               
 1949       continue
c
c   Now find out where the sign change is.  This is where the line of sight
c   has moved over the horizon.
c
            do 1940 ialf=kk,nalf-1
              if(ialf.lt.nalf)then
                index=ialf+1
              else
                index=kk
              endif
              rsign=arrproj(ialf)*arrproj(index)
              if(rsign.le.0.0d0)then     ! crossed over
                if(arrproj(ialf).gt.0.0)then
                  izz=ialf
                  jzz=ibet
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xvis=xarray(iidx)
                  yvis=yarray(iidx)
                  zvis=zarray(iidx)
                  rvis=radarray(iidx)
                  phivis=phiar(iidx)
                  iedgestar(iidx)=30     ! just appeared at limb
c
                  izz=index
                  jzz=ibet
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xhid=xarray(iidx)
                  yhid=yarray(iidx)
                  zhid=zarray(iidx)
                  rhid=radarray(iidx)
                  phihid=phiar(iidx)
                else
                  izz=index
                  jzz=ibet
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xvis=xarray(iidx)
                  yvis=yarray(iidx)
                  zvis=zarray(iidx)
                  rvis=radarray(iidx)
                  phivis=phiar(iidx)
                  iedgestar(iidx)=30     ! just appeared at limb
c
                  izz=ialf
                  jzz=ibet
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xhid=xarray(iidx)
                  yhid=yarray(iidx)
                  zhid=zarray(iidx)
                  rhid=radarray(iidx)
                  phihid=phiar(iidx)
                endif
c
                call acchor(overQ,psi0,omega,
     $              xvis,yvis,zvis,rvis,phivis,xhid,yhid,zhid,rhid,
     $              phihid,ax,ay,az,xacc,yacc,zacc,bdist,tidephi,itide)
c
                xx=xacc  
                yy=yacc  
                zz=zacc  
c
                xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)  
                yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)    
                Ndummy=Ndummy+1
                xdummy(Ndummy)=xp
                ydummy(Ndummy)=yp
              endif
 1940       continue
 1950     continue
 1960     continue
c
c   Check the visibility of the nose and see if it should be in
c   the horzon
c
          if((phase.eq.90.0d0).or.(phase.eq.270.0d0))then   ! both ends visible
            xp=xtran(xend(3),0.0d0,0.0d0,phase,fincr,Q,istar,bdist) 
            yp=ytran(xend(3),0.0d0,0.0d0,phase,fincr,Q,istar,bdist)   
            Ndummy=Ndummy+1      ! outside the first ring.  Include the
            xdummy(Ndummy)=xp    ! point in the horizon
            ydummy(Ndummy)=yp
            xp=xtran(xend(4),0.0d0,0.0d0,phase,fincr,Q,istar,bdist) 
            yp=ytran(xend(4),0.0d0,0.0d0,phase,fincr,Q,istar,bdist)    
c            Ndummy=Ndummy+1      ! outside the first ring.  Include the
c            xdummy(Ndummy)=xp    ! point in the horizon
c            ydummy(Ndummy)=yp
c            go to 69
          endif
c
c   Now 'sort' the horizon points so that a regular polygon is made.  If
c   we are near the quadrature phases, then use the old algorithm.  For
c   other phases, use the new algorithm which returns horizon points for
c   every degree in polar coordinates.
c
c          if(((phase.gt.80.0d0).and.(phase.lt.100.0d0)).or.
c     %       ((phase.gt.260.0d0).and.(phase.lt.280.0d0)).and.idraw.ge.1)then
c 69         call sortcircle(Ndummy,xdummy,ydummy)
cc
cc   The above array may have repeated points.  Remove them.
cc
c            call uniquepoint(Ndummy,xdummy,ydummy,Nhoriz,xhoriz,yhoriz)
c          else
c
            call newsortcircle(Ndummy,xdummy,ydummy,Nhoriz,
     $         xhoriz,yhoriz,ibetmax,xhmin,xhmax,yhmin,yhmax)

c            write(*,*)istar,phase,Ndummy,Nhoriz

c          endif
c
          return
          end
c
c  %%%%%%%%%%%%%%%%%%%%%
c
          subroutine uniquepoint(Nbig,xbig,ybig,Ntrim,xtrim,ytrim)
c
c   October 8, 1999
c
c   This routine will take the large arrays xbig and ybig and remove
c   possible duplicate points.  It is assumed that the arrays xbig
c   and ybig have been passed to 'sortcircle' so that duplicate points
c   are in adjacent array positions.
c
          implicit double precision(a-h,o-z)
c
          dimension xbig(Nbig),ybig(Nbig),xtrim(Nbig),ytrim(Nbig)
c
c   UPDATE JULY 2, 2004
c
c   Change the threshold to 1.0d-5
c
          Ntrim=0
          do 10 i=2,Nbig
            t1=xbig(i)-xbig(i-1)
            t2=ybig(i)-ybig(i-1)
            distance=dsqrt(t1*t1+t2*t2)
            if(distance.gt.1.0d-5)then
              Ntrim=Ntrim+1
              xtrim(Ntrim)=xbig(i-1)
              ytrim(Ntrim)=ybig(i-1)
            endif
 10       continue
c
c   Are the two points at the end the same?  If not, then write the last one.
c
          t1=xbig(Nbig)-xbig(Nbig-1)
          t2=ybig(Nbig)-ybig(Nbig-1)
          distance=dsqrt(t1*t1+t2*t2)
          if(distance.gt.1.0d-5)then
            Ntrim=Ntrim+1
            xtrim(Ntrim)=xbig(Nbig)
            ytrim(Ntrim)=ybig(Nbig)
          endif
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine insidecircle(N,xcir,ycir,xp,yp,iyes,icut)
c
c
c   UPDATE May 8, 2006
c
c   Get rid of the call to checkbetween since that takes up lots of CPU time.
c   Use a simple if-then statement.
c
c    given a polygon with N points (they need to be in order), this
c    routine will check to see if the point (xp,yp) is inside the polygon
c
c    iyes=0 means that the point is outside
c    iyes=100 means that the point is inside
c
          implicit double precision(a-h,o-z)
c
          dimension xcir(N),ycir(N)
c
          iyes=0
c
c          call getmaxvalues(N,xcir,ycir,xmin,xmax,ymin,ymax)  !get the extreme
cc                                                             !values
c          if(xp.le.xmin)then         ! check to see if the
c             iyes=0
c             icut=-100               ! point is anywhere
c             return                  ! near the projected polygon
c          endif
cc
c          if(xp.ge.xmax)then
c             iyes=0
c             icut=-100
c             return
c          endif
cc
c          if(yp.le.ymin)then
c             iyes=0
c             icut=2
c             return
c          endif
cc
c          if(yp.ge.ymax)then
c             iyes=0
c             icut=0
c             return
c          endif
c
c    Now draw a line parallel to the y-axis to see how many sides are 
c    intersected -- if the number is even, then the point is outside.
c
          icut=0
          do 10 i=1,N-1
c
c            xx1=xcir(i)
c            xx2=xcir(i+1)
c            call checkbetween(xp,xx1,xx2,between)
c
c            if(between)then
            if(((xcir(i).le.xp).and.(xp.le.xcir(i+1))).or.
     &          ((xcir(i+1).le.xp).and.(xp.le.xcir(i))))then  
              call getline(yline,xp,yp,xcir(i),xcir(i+1),ycir(i),
     $            ycir(i+1))
              if((yp.le.yline))icut=icut+1
            endif
 10       continue
c
c          xx1=xcir(1)
c          xx2=xcir(N)
c
c          call checkbetween(xp,xx1,xx2,between)
c          if(between)then
c
          if(((xcir(1).le.xp).and.(xp.le.xcir(N))).or.
     &          ((xcir(N).le.xp).and.(xp.le.xcir(1))))then  
            call getline(yline,xp,yp,xcir(1),xcir(N),ycir(1),
     $            ycir(N))
            if((yp.le.yline))icut=icut+1
          endif
c
          if((icut.eq.1).or.(icut.eq.3).or.(icut.eq.5))then
            iyes=100    ! inside
          endif
          if((icut.eq.0).or.(icut.eq.2).or.(icut.eq.4))then
            iyes=0      ! outside
          endif
c
c         if((yp.lt.ymax).and.(iyes.eq.0))icut=2
c
          return
          end
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   UPDATE May 8, 2006
c
c   retire this subroutine, since it takes up too much CPU time.
c
c          subroutine checkbetween(xp,x1,x2,between)
cc
cc    October 8, 1999
cc
cc    This routine will check to see if the point xp is bewteen x1 and
cc    x2 on the x-axis.  between = .true. if this is the case.
cc    This is from 1995-1996 changes to the Avni code
cc    by Orosz.
cc
c          implicit double precision(a-h,o-z)
c
c          logical between
cc
c          between=.false.
cc
c          if((xp.ge.x1).and.(xp.lt.x2))between=.true.
c          if((xp.gt.x2).and.(xp.le.x1))between=.true.
cc          
cc    the points x1 and x2 may come in reverse order, so two statements
cc    are needed
cc
c          return
c          end
c
c ==========================================================================
c
          subroutine getline(yline,xp,yp,x1,x2,y1,y2)
c
c    October 8, 1999
c
c    Check to see if the point (xp,yp) is above the line connecting
c    (x1,y1) and (x2,y2).  This is from 1995-1996 changes to the Avni code
c    by Orosz.
c
c
          implicit double precision(a-h,o-z)
c
         if(x1-x2.eq.0.0d0)then
           yline=(y1)
c           write(*,*)'error in getline: x1=x2'
           return
         endif
c
         slope=(y1-y2)/(x1-x2)
c
         yline=slope*(xp-x1)+y1
c
         return
         end
c
c  *************************************************************************
c
          subroutine getvisib(istar,ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     %      phase,finc,Q,psi0,omega,
     $      gradx,grady,gradz,xarray,yarray,zarray,
     $      xend,visib,Nhoriz,xhoriz,yhoriz,idint,Ndhoriz,dxhoriz,dyhoriz,
     %      Ndtop,dtopx,dtopy,
     $      Ncoords,xcoords,ycoords,projarray,iecheck,Neclipse,phiar,
     %      radarray,delphi,separ,iedgehor,bdist,mmdx,jdum,MonteCarlo,
     #      isw13,ialfmin,ialfmax,xhmin,xhmax,yhmin,yhmax,tidephi,itide,
     #      iedgestar,phihor,phistart)

c
c  October 9, 1999
c
c  This routine will compute the 'projection' factor of each grid element on
c  the star (istar=1 to do star 1, istar=2 to do star2), check for eclipses
c  (the horizon of the other body is in xhoriz,yhoriz), and return the
c  sky coordinates of the visible points.  Set iecheck = -1 to skip the
c  check for 
c  eclipses.
c
c  UPDATE JULY 4, 2004
c
c  Add jdum and MonteCarlo to the argument list.  If MonteCarlo > 10,
c  then use Monte Carlo integration to determine fractionally
c  eclipsed pixels.  If MonteCarlo < 10, then proceed as before
c  and use interpolation in getBBflux and getATMflux.
c
c   UPDATE May 3, 2006
c
c   Add a "fast transit" mode.  If isw13 > 0, then consider ialf ranges
c   only between ialfmin and ialfmax
c
c   UPDATE May 8, 2006
c
c   Add the minimum and maximum x and y-values of the horizon.  Do
c   an if-then check to determine of insidecircle should be called.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          parameter(twopie=2.0d0*pie)
          dimension visib(ialphmax*ibetmax),xcoords(ialphmax*ibetmax*4),
     $        gradx(ialphmax*ibetmax),grady(ialphmax*ibetmax),
     $        gradz(ialphmax*ibetmax),ycoords(ialphmax*ibetmax*4),
     $        xhoriz(Nhoriz),yhoriz(Nhoriz),xend(4),xarray(ialphmax*ibetmax),
     $        yarray(ialphmax*ibetmax),zarray(ialphmax*ibetmax),
     $        projarray(ialphmax*ibetmax),dxhoriz(Ndhoriz),dyhoriz(Ndhoriz),
     %        dtopx(Ndtop),dtopy(Ndtop),ibetlim(ialphmax),
     $        phiar(ialphmax*ibetmax),ihid(200000),  !ialphmax*ibetmax
     %        radarray(ialphmax*ibetmax),delphi(ialphmax*ibetmax),
     #        iedgehor(ialphmax*ibetmax),mmdx(ialphmax,ibetmax),
     $        iedgestar(ialphmax*ibetmax)
          dimension phihor(ialphmax,2),phistart(ialphmax)
c
c   UPDATE JULY 7, 2004
c
c   Use the sub-random Sobel sequence instead of ran9.  xsob is needed
c   for this.
c
          dimension xsob(2)
c
c
c  initialize the visibities!
c
          do 1 i=1,nalf
            do 2 j=1,ibetlim(i)  !4*nbet
c              iidx=(i-1)*ibetlim(i)+j
              iidx=mmdx(i,j)
              visib(iidx)=0.0d0
              projarray(iidx)=-1.0d0
              ihid(iidx)=-1
              delphi(iidx)=-999.9d0
              iedgehor(iidx)=-999
 2          continue
 1        continue
c
c    UPDATE May 3, 2006
c
c    Keep track of the smallest and largest values of ialf on star 1 that
c    are eclipsed.
c
          iasmall=999999
          iabig=-111111
c
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c          
          NBET4 = NBET*4
          DBETA = (pie/2.0d0)/NBET        ! step size in longitude
c
          AZ = DCOS(FINCR)
          IF (AZ.LT.0.0d0) AZ = 0.0d0
          AX = -DSIN(FINCR)*DCOS(PHASER)    ! l in Wilson & Sofia
          AY = DSIN(FINCR)*DSIN(PHASER)     ! m in Wilson & Sofia
c
c   Check to see of the star in question is in front.  If so, then simply
c   find the projection factors.
c
          infront=0
          Ncoords=0
          Neclipse=0
          if((istar.eq.1).and.((phase.ge.0.0d0).
     #           and.(phase.lt.90.0d0)))infront=1
          if((istar.eq.1).and.((phase.ge.270.0d0).
     #           and.(phase.le.360.0d0)))infront=1
          if((istar.eq.2).and.((phase.ge.0.0d0).
     #           and.(phase.lt.90.0d0)))infront=1
          if((istar.eq.2).and.((phase.ge.270.0d0).
     $           and.(phase.le.360.0d0)))infront=1
c
          do 10 i=1,2
            if(i.eq.1)proj=AX
            if(i.eq.2)proj=-AX
            if(proj.gt.0.0d0)then
              xp=xtran(xend(i+2),0.0d0,0.0d0,phase,fincr,Q,istar,bdist) 
              yp=ytran(xend(i+2),0.0d0,0.0d0,phase,fincr,Q,istar,bdist)
c
c   If we are looking at star 2 and there is a disk, then check to see if
c   the points are *inside* the top horizon of the disk.
c
              if((idint.ge.1).and.(istar.eq.2))then
                if(iecheck.ge.0)then
                  iyes=-100
                  call insidecircle(Ndtop,dtopx,dtopy,xp,yp,iyes,icut)
                  if(iyes.ne.100)go to 10   ! the point is outside the
                endif                       ! top horizon, so it is beneath the
              endif                         ! rim.

              if(infront.eq.1)then        
                Ncoords=Ncoords+1    
                xcoords(Ncoords)=xp
                ycoords(Ncoords)=yp
              else
                if(iecheck.ge.0)then
                  iyes=-100
                  call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                  if(iyes.eq.100)then
                    go to 10       ! eclipsed by star 2
                  endif
                endif
c
c   UPDATE JUNE 17, 2005
c
c   Add if.idint.gt.0 to if-then
c
                if((istar.eq.1).and.(idint.gt.0))then 
                  if(iecheck.ge.0)then 
                    call insidecircle(Ndhoriz,dxhoriz,dyhoriz,xp,yp,iyes,icut)
                    if(iyes.eq.100)then
                      go to 10     ! eclipsed by disk
                    endif
                  endif
                endif
              
                Ncoords=Ncoords+1
                xcoords(Ncoords)=xp
                ycoords(Ncoords)=yp
              endif
            endif
 10       continue
c

          iimin=123456
          iimax=-12345
          DO 501 IALF = 1, NALF
            if((isw13.gt.0).and.(istar.eq.1))then
              if((ialf.lt.ialfmin).or.(ialf.gt.ialfmax))go to 501
            endif
            DO 502 IBET = 1,ibetlim(ialf)      !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              PROJ = AX * GRADX(iidx) + AY*GRADY(iidx) + 
     1	        AZ*GRADZ(iidx)
              projarray(iidx)=proj
              IF (PROJ.LT.0.) GO TO 502    ! is the surface element visible?
              xx=xarray(iidx)
              yy=yarray(iidx)
              zz=zarray(iidx)
              xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) ! projected coords
              yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c   Check to see of the point in question is eclipsed by the other star, or
c   in the case of star 1, eclipsed by the disk, or
c   in the case of a point on the bottom half of star 2, eclipsed by the disk. 
c   icut=2 for points outside and below the horizon.
c 
              if(infront.eq.1)then 
                if((idint.ge.1).and.(istar.eq.2).and.(zz.lt.0.0d0))then
                  iyes=-100
                  call insidecircle(Ndhoriz,dxhoriz,dyhoriz,xp,yp,iyes,icut)
                  if(iyes.eq.100)then
                    ihid(iidx)=2     !eclipsed by disk
                    go to 502   
                  endif
                endif
                if((idint.ge.1).and.(istar.eq.2).and.(zz.ge.0.0d0))then
                  icut=-100
                  iyes=-100
                  call insidecircle(Ndtop,dtopx,dtopy,xp,yp,iyes,icut)
                  if((iyes.eq.0).and.(icut.eq.2))then
                    ihid(iidx)=3
                    go to 502
                  endif
                endif
                Ncoords=Ncoords+1  
                xcoords(Ncoords)=xp
                ycoords(Ncoords)=yp
                visib(iidx)=proj
              else
                iyes=-100
                if(iecheck.ge.0)then
c
c   UPDATE May 8, 2006
c
c   Add these if-then statements to see if insidecircle should be called.
c              
                  if(((xhmin.lt.xp).and.(xp.lt.xhmax)).and.
     #                ((yhmin.lt.yp).and.(yp.lt.yhmax)))then
                    call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                    if(iyes.eq.100)then
                      Neclipse=Neclipse+1
                      ihid(iidx)=1
c
c   UPDATE JULY 1, 2004
c
c   Add this to determine the latitude rows above and below the eclipsing
c   horizon.
c
                      if(ialf.lt.iimin)iimin=ialf
                      if(ialf.gt.iimax)iimax=ialf
                      go to 502
                    endif
                  endif
                endif
                if((idint.ge.1).and.(istar.eq.2).and.(zz.lt.0.0d0))then
                  iyes=-100
                  call insidecircle(Ndhoriz,dxhoriz,dyhoriz,xp,yp,iyes,icut)
                  if(iyes.eq.100)then
                    ihid(iidx)=2
                    go to 502
                  endif   
                endif
                if((idint.ge.1).and.(istar.eq.1).and.(iecheck.ge.0))then
                  iyes=-100
                  call insidecircle(Ndhoriz,dxhoriz,dyhoriz,xp,yp,iyes,icut)
                  if(iyes.eq.100)then
                    Neclipse=Neclipse+1
                    ihid(iidx)=2
                    go to 502   
                  endif
                endif
c
c   Finally, if the star in question is star 2 and there is a disk, we
c   must check to see if the point is inside the top horizon of the disk.
c   Points outside this horizon and below it are beneath the rim.  Points
c   outside this horizon and above it are still visible.  icut=2
c   for points outside and below.
c
                if((idint.ge.1).and.(istar.eq.2).and.(zz.ge.0.0d0))then
                  icut=-100
                  iyes=-100
                  call insidecircle(Ndtop,dtopx,dtopy,xp,yp,iyes,icut)
                  if((iyes.eq.0).and.(icut.eq.2))then
                    ihid(iidx)=3
                    go to 502
                  endif
                endif
                Ncoords=Ncoords+1
                xcoords(Ncoords)=xp
                ycoords(Ncoords)=yp
                visib(iidx)=proj
              endif                    ! end if in front
c
502         CONTINUE
 501      continue         ! continue the alpha loop
c
c          if(Neclipse.eq.0)return
c
c   Now we have to go along the beta direction and find out which visible
c   point is nearest to the eclipsing horizon.
c
          do 602 ialf=1,Nalf 
            if((isw13.gt.0).and.(istar.eq.1))then
              if((ialf.lt.ialfmin).or.(ialf.gt.ialfmax))go to 602
            endif
            do 601 ibet=1,ibetlim(ialf)      !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              if(projarray(iidx).le.0.0d0)go to 601
              if(ibet.lt.ibetlim(ialf))then
                index=ibet+1
              else
                index=1
              endif
c
              izz=ialf
              jzz=index
c              jjdx=kount(ialphmax,izz,ibetlim)+jzz
              jjdx=mmdx(izz,jzz)
              if(projarray(jjdx).le.0.0d0)go to 601
              izz=ialf
              jzz=ibet
c              iidx=kount(ialphmax,izz,ibetlim)+jzz
              iidx=mmdx(izz,jzz)
              isign=ihid(iidx)*ihid(jjdx)
c
              if(isign.le.-1)then  ! crossed through the horizon
                if(ihid(iidx).eq.-1)then
                  izz=ialf
                  jzz=ibet
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xvis=xarray(iidx)
                  yvis=yarray(iidx)
                  zvis=zarray(iidx)
                  rvis=radarray(iidx)
                  phivis=phiar(iidx)
                  iedgehor(iidx)=10
c
                  izz=ialf
                  jzz=index
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xhid=xarray(iidx)
                  yhid=yarray(iidx)
                  zhid=zarray(iidx)
                  rhid=radarray(iidx)
                  phihid=phiar(iidx)
                  iedgehor(iidx)=-10
                  index=ibet
                else
                  izz=ialf
                  jzz=index
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xvis=xarray(iidx)
                  yvis=yarray(iidx)
                  zvis=zarray(iidx)
                  rvis=radarray(iidx)
                  phivis=phiar(iidx)
                  iedgehor(iidx)=20
c
                  izz=ialf
                  jzz=ibet
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xhid=xarray(iidx)
                  yhid=yarray(iidx)
                  zhid=zarray(iidx)
                  rhid=radarray(iidx)
                  phihid=phiar(iidx)
                  iedgehor(iidx)=-10
                endif                
                xpv=xtran(xvis,yvis,zvis,phase,fincr,Q,istar,bdist)  
                ypv=ytran(xvis,yvis,zvis,phase,fincr,Q,istar,bdist)    
                xph=xtran(xhid,yhid,zhid,phase,fincr,Q,istar,bdist)  
                yph=ytran(xhid,yhid,zhid,phase,fincr,Q,istar,bdist)    
c
c                if(isign.eq.-1)call clip(Nhoriz,xhoriz,yhoriz,
c     $             xpv,ypv,xph,yph,xedge,yedge)
c                if(isign.eq.-2)call clip(Ndhoriz,dxhoriz,dyhoriz,
c     $             xpv,ypv,xph,yph,xedge,yedge)
c                if(isign.eq.-3)call clip(Ndtop,dtopx,dtopy,
c     $             xpv,ypv,xph,yph,xedge,yedge)
c
c   UPDATE March 26, 2002
c
c   Remove the variable separ from the argument list of accphi.
c
c
c   UPDATE JULY 4, 2004
c
c   If MonteCarlo<10, then we need these calls to use interpolation.
c
                if(MonteCarlo.lt.10)then
                  if(isign.eq.-1)call accphi(Q,psi0,omega,phase,fincr,istar,
     &               xvis,yvis,zvis,rvis,phivis,
     $               xhid,yhid,zhid,rhid,phihid,Nhoriz,xhoriz,yhoriz,
     %               phiacc,bdist,tidephi,itide)
c
                  if(isign.eq.-2)call accphi(Q,psi0,omega,phase,fincr,istar,
     &               xvis,yvis,zvis,rvis,phivis,
     $               xhid,yhid,zhid,rhid,phihid,Ndhoriz,dxhoriz,dyhoriz,
     %               phiacc,bdist,tidephi,itide)
c
                  if(isign.eq.-3)call accphi(Q,psi0,omega,phase,fincr,istar,
     &               xvis,yvis,zvis,rvis,phivis,
     $               xhid,yhid,zhid,rhid,phihid,Ndtop,dtopx,dtopy,
     %               phiacc,bdist,tidephi,itide)
c
c                 iidx=kount(ialphmax,ialf,ibetlim)+index
                  iidx=mmdx(ialf,index)
                  diff=(phiar(iidx)-phiacc)
c
                  if(dabs(diff).lt.2.0d0)then
                    delphi(iidx)=(phiar(iidx)-phiacc)
                  endif
                  if(dabs(diff).ge.2.0d0)then
                    if(phiar(iidx).gt.4.0d0)then
                      phitemp=2.0d0*pie-phiar(iidx)
                      delphi(iidx)=(phitemp+phiacc)
                    endif
                    if(phiacc.gt.4.0d0)then
                      phitemp=2.0d0*pie-phiacc
                      delphi(iidx)=(phitemp+phiar(iidx))
                    endif
                  endif
                endif      !    if MonteCarlo < 10
              endif
c
c
601        continue
602      continue
c

c   UPDATE JULY 4, 2004
c
c   if -10 < MonteCarlo < 10, then we are done.
c
         if((MonteCarlo.lt.10).and.(MonteCarlo.gt.-10))return

c   UPDATE JULY 1, 2004
c
c   Add a loop which will use Monte Carlo integration to
c   determine what fraction of a pixel is partially eclipsed
c   by the horizon of the star in front.  The variable MonteCarlo
c   sets the number of integrations.
c
          dtheta=pie/dble(Nalf)
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
c
          iimin=iimin-1
          if(MonteCarlo.gt.10)then
            if(iimin.lt.1)iimin=1
          endif
          iimax=iimax+1
          if(MonteCarlo.gt.10)then
            if(iimax.gt.Nalf)iimax=Nalf
          endif
c
          do 702 ialf=1,Nalf
            if((isw13.gt.0).and.(istar.eq.1))then
              if((ialf.lt.ialfmin).or.(ialf.gt.ialfmax))go to 7701
            endif
            r=0.0000001d0
            theta=-0.5d0*dtheta+dtheta*dble(ialf)
            do 701 ibet=1,ibetlim(ialf)      !4*Nbet
              iidx=mmdx(ialf,ibet)
c
              if(projarray(iidx).lt.0.0d0)go to 7701
c
c   UPDATE July 14, 2004
c
c   We have two options here.  If MonteCarlo > 10, then check
c   all pixels near the eclipsing horizon and get fractional
c   areas.  If MonteCarlo < -10, then check only the latitude row just
c   above the eclipsing horizon and the latitude row just below the
c   eclipsing horizon.
c
              if(MonteCarlo.gt.10)then
                if((iedgehor(iidx).eq.-10).or.(iedgehor(iidx).gt.5)
     %            .or.(ialf.eq.iimin).or.(ialf.eq.iimax)
     %            .or.(ialf.eq.iimin+1).or.(ialf.eq.iimax-1))then
                  dphi=twopie/dble(ibetlim(ialf))
                  phi=-0.5d0*dphi+dphi*dble(ibet)
                  phi=phi+phistart(ialf)
c
c             if(dble(ialf/2).eq.(dble(ialf)*0.5d0))phi=phi+0.25d0*dphi    !DPHI
c
c
c   This loop will find random theta,phi locations on the pixel in
c   question, compute x,y,z coordinates on the star, and determine
c   if that location is eclipsed.  It will then determine the fraction
c   of the pixel that is eclipsed.
c
                  Nloop=MonteCarlo
                  NNhid=0
c
c   Do a check of the pixel corners and edges to see if they are eclipsed.
c   If all corners and all edges are eclipsed, then assume that the pixel
c   is completely hidden.  Likewise, if none of the corners or edges
c   are eclipsed, assume the pixel is completely visible
c
                  do 901 kk=-3,3
                    do 900 ll=-3,3
                      phinew=(dble(ll)/6.0d0)*dphi+phi
                      thetanew=(dble(kk)/6.0d0)*dtheta+theta          
                      snth=dsin(thetanew)
                      snth3=dsin(thetanew)/3.0d0  !*0.333333333333333d0
                      cnth=dcos(thetanew)
                      cox=dcos(phinew)*snth                !*dsin(theta)
                      coy=dsin(phinew)*snth                !*dsin(theta)
                      coz=cnth                          !dcos(theta)
                      CALL RAD(overQ,omega,cox,coy,coz,psi0,r,xx,yy,zz,1,
     $                 bdist,tidephi,itide)
                      call POTEN(overQ,omega,xx,yy,zz,psi,
     $                   psix,psixx,psiy,psiz,1,bdist,
     #                   cox,coy,tidephi,itide)
                      xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) 
                      yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
                      iyes=-1000
                      call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                      if(iyes.eq.100)then
                        NNhid=NNhid+1
                      endif
 900                continue
 901              continue
c
c                  write(*,*)istar,ialf,ibet,NNhid
                  Nkeep=NNhid
                  if(NNhid.eq.0)go to 801   ! completely visible
                  if(NNhid.eq.49)then
                    NNhid=Nloop
                    go to 801        ! completely hidden
                  endif
                  NNhid=0
c
c   UPDATE July 7, 2004
c
c   Use the Sobel sequence instead of ran9.
c
                  nnn=2
                  do 800 jj=1,Nloop
c  
c                    phinew=(ran9(jdum)-0.5d0)*dphi+phi
c                    thetanew=(ran9(jdum)-0.5d0)*dtheta+theta          
c  
                    call sobseq(nnn,xsob)
                    phinew=(xsob(1)-0.5d0)*dphi+phi
                    thetanew=(xsob(2)-0.5d0)*dtheta+theta          
c  
                    snth=dsin(thetanew)
                    snth3=dsin(thetanew)/3.0d0  !*0.333333333333333d0
                    cnth=dcos(thetanew)
                    cox=dcos(phinew)*snth                !*dsin(theta)
                    coy=dsin(phinew)*snth                !*dsin(theta)
                    coz=cnth                          !dcos(theta)
                    CALL RAD(overQ,omega,cox,coy,coz,psi0,r,xx,yy,zz,1,
     $                bdist,tidephi,itide)
                    call POTEN(overQ,omega,xx,yy,zz,psi,
     $                   psix,psixx,psiy,psiz,1,bdist,cox,coy,tidephi,itide)
                    xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) 
                    yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c                    write(*,1001)xp,yp,ialf,ibet
                    iyes=-1000
                    call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                    if(iyes.eq.100)then
                      NNhid=NNhid+1
                    endif
 800              continue
 801              nvis=Nloop-NNhid
                  visib(iidx)=projarray(iidx)*(dble(nvis)/dble(Nloop))
                endif
              endif  ! end if MonteCarlo > 10
c
c  Check only the latitude row just above and just below.
c
              if(MonteCarlo.lt.-10)then
                if((ialf.eq.iimin).or.(ialf.eq.iimax))then
                  dphi=twopie/dble(ibetlim(ialf))
                  phi=-0.5d0*dphi+dphi*dble(ibet)
                  phi=phi+phistart(ialf)
c
c              if(dble(ialf/2).eq.(dble(ialf)*0.5d0))phi=phi+0.25d0*dphi  !DPHI
c
c
c   This loop will find random theta,phi locations on the pixel in
c   question, compute x,y,z coordinates on the star, and determine
c   if that location is eclipsed.  It will then determine the fraction
c   of the pixel that is eclipsed.
c
                  Nloop=(abs(MonteCarlo))
                  NNhid=0
c
c   Do a check of the pixel corners and edges to see if they are eclipsed.
c   If all corners and all edges are eclipsed, then assume that the pixel
c   is completely hidden.  Likewise, if none of the corners or edges
c   are eclipsed, assume the pixel is completely visible
c
                  do 5901 kk=-3,3
                    do 5900 ll=-3,3
                      phinew=(dble(ll)/6.0d0)*dphi+phi
                      thetanew=(dble(kk)/6.0d0)*dtheta+theta          
                      snth=dsin(thetanew)
                      snth3=dsin(thetanew)/3.0d0  !*0.333333333333333d0
                      cnth=dcos(thetanew)
                      cox=dcos(phinew)*snth                !*dsin(theta)
                      coy=dsin(phinew)*snth                !*dsin(theta)
                      coz=cnth                          !dcos(theta)
                      CALL RAD(overQ,omega,cox,coy,coz,psi0,r,xx,yy,zz,1,
     $                   bdist,tidephi,itide)
                      call POTEN(overQ,omega,xx,yy,zz,psi,
     $                   psix,psixx,psiy,psiz,1,bdist,cox,coy,tidephi,itide)
                      xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) 
                      yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
                      iyes=-1000
                      call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                      if(iyes.eq.100)then
                        NNhid=NNhid+1
                      endif
 5900               continue
 5901             continue
c
                  Nkeep=NNhid
                  if(NNhid.eq.0)go to 5801   ! completely visible
                  if(NNhid.eq.49)then
                    NNhid=Nloop
                    go to 5801        ! completely hidden
                  endif
                  NNhid=0
c
c   UPDATE July 7, 2004
c
c   Use the Sobel sequence instead of ran9.
c
                  nnn=2
                  do 5800 jj=1,Nloop
c
c                  phinew=(ran9(jdum)-0.5d0)*dphi+phi
c                  thetanew=(ran9(jdum)-0.5d0)*dtheta+theta          
c
                    call sobseq(nnn,xsob)
                    phinew=(xsob(1)-0.5d0)*dphi+phi
                    thetanew=(xsob(2)-0.5d0)*dtheta+theta          
c  
                    snth=dsin(thetanew)
                    snth3=dsin(thetanew)/3.0d0  !*0.333333333333333d0
                    cnth=dcos(thetanew)
                    cox=dcos(phinew)*snth                !*dsin(theta)
                    coy=dsin(phinew)*snth                !*dsin(theta)
                    coz=cnth                          !dcos(theta)
                    CALL RAD(overQ,omega,cox,coy,coz,psi0,r,xx,yy,zz,1,
     $               bdist,tidephi,itide)
                    call POTEN(overQ,omega,xx,yy,zz,psi,
     $                   psix,psixx,psiy,psiz,1,bdist,cox,coy,tidephi,itide)
                    xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) 
                    yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c                    write(*,1001)xp,yp,ialf,ibet
                    iyes=-1000
                    call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                    if(iyes.eq.100)then
                      NNhid=NNhid+1
                    endif
 5800             continue
 5801             nvis=Nloop-NNhid
c                  write(*,987)dble(nvis)/dble(Nloop),ialf,ibet,istar
                  visib(iidx)=projarray(iidx)*(dble(nvis)/dble(Nloop))
                endif
              endif ! end if MonteCarlo < -10
c
c  UPDATE January 9, 2010
c
c  Check points near the limb and find the fractions of points
c  visible.
c
7701         go to 701    !this part really does not work, so skip it.

                if(MonteCarlo.gt.10)then
                if((iedgestar(iidx).eq.-10).or.(iedgestar(iidx).gt.5))then
                  dphi=twopie/dble(ibetlim(ialf))
                  phi=-0.5d0*dphi+dphi*dble(ibet)
                  phi=phi+phistart(ialf)
c
c              if(dble(ialf/2).eq.(dble(ialf)*0.5d0))phi=phi+0.25d0*dphi   !DPHI
c
c
c   This loop will find random theta,phi locations on the pixel in
c   question, compute x,y,z coordinates on the star, and determine
c   if that location is visible.  It will then determine the fraction
c   of the pixel that is visible.
c
                  Nloop=40
                  NNhid=0
c
c   Do a check of the pixel corners and edges to see if they are eclipsed.
c   If all corners and all edges are eclipsed, then assume that the pixel
c   is completely hidden.  Likewise, if none of the corners or edges
c   are eclipsed, assume the pixel is completely visible
c
c                  do 9901 kk=-3,3
c                    do 9900 ll=-3,3
c                      phinew=(dble(ll)/6.0d0)*dphi+phi
c                      thetanew=(dble(kk)/6.0d0)*dtheta+theta          
c                      snth=dsin(thetanew)
c                      snth3=dsin(thetanew)/3.0d0  !*0.333333333333333d0
c                      cnth=dcos(thetanew)
c                      cox=dcos(phinew)*snth                !*dsin(theta)
c                      coy=dsin(phinew)*snth                !*dsin(theta)
c                      coz=cnth                          !dcos(theta)
c                      CALL RAD(overQ,omega,cox,coy,coz,psi0,r,xx,yy,zz,1,
c     $                 bdist,tidephi,itide)
c                      call POTEN(overQ,omega,xx,yy,zz,psi,
c     $                   psix,psixx,psiy,psiz,1,bdist,
c     #                   cox,coy,tidephi,itide)
c                      xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) 
c                      yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c                      iyes=-1000
c                      call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
c                      if(iyes.eq.100)then
c                        NNhid=NNhid+1
c                      endif
c 9900               continue
c 9901             continue
c
c                  write(*,*)istar,ialf,ibet,NNhid
c                  Nkeep=NNhid
c
c                  if(NNhid.eq.0)go to 8801   ! completely visible
c                  if(NNhid.eq.49)then
c                    NNhid=Nloop
c                    go to 8801        ! completely hidden
c                  endif
c
                  NNhid=0
                  psum=0.0d0
c
c   Use the Sobel sequence instead of ran9.
c
                  nnn=2
                  do 8800 jj=1,Nloop
c  
                    call sobseq(nnn,xsob)
                    phinew=(xsob(1)-0.5d0)*dphi+phi
                    thetanew=(xsob(2)-0.5d0)*dtheta+theta          
c  
                    snth=dsin(thetanew)
                    snth3=dsin(thetanew)/3.0d0  !*0.333333333333333d0
                    cnth=dcos(thetanew)
                    cox=dcos(phinew)*snth                !*dsin(theta)
                    coy=dsin(phinew)*snth                !*dsin(theta)
                    coz=cnth                          !dcos(theta)
                    CALL RAD(overQ,omega,cox,coy,coz,psi0,r,xx,yy,zz,1,
     $                bdist,tidephi,itide)
                    call POTEN(overQ,omega,xx,yy,zz,psi,
     $                   psix,psixx,psiy,psiz,1,bdist,cox,coy,tidephi,itide)
                    gravity=dsqrt(PSIX**2+PSIY**2+PSIZ**2)
                    GX = -PSIX/gravity
                    GY = -PSIY/gravity
                    GZ = -PSIZ/gravity
                    proj=ax*gx+ay*gy+az*gz
                    write(*,*)proj,visib(iidx)
                    if(proj.lt.0.0d0)then
                      NNhid=NNhid+1
c                      go to 8800
                    endif
c
c   Check to see if the pixel is eclipsed
c
                    if((iedgehor(iidx).eq.-10).or.(iedgehor(iidx).gt.5)
     %                  .or.(ialf.eq.iimin).or.(ialf.eq.iimax)
     %                 .or.(ialf.eq.iimin+1).or.(ialf.eq.iimax-1))then
c
                      xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) 
                      yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c                     write(*,1001)xp,yp,ialf,ibet
                      iyes=-1000
                      call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                      if(iyes.eq.100)then
                        NNhid=NNhid+1
c                        go to 8800
                        proj=0.0d0
                      endif
                    endif
c                    if(proj.lt.0.0d0)proj=0.0d0
                    psum=psum+proj
 8800             continue
 8801             nvis=Nloop-NNhid
                  write(*,*)Nloop,NNhid,visib(iidx),psum/dble(Nloop)
                  if(psum.lt.0.0d0)psum=0.0d0
                  visib(iidx)=psum !projarray(iidx)*(dble(nvis)/dble(Nloop))
                endif
              endif  ! end if MonteCarlo > 10
c
 701        continue
 702      continue
c
c
c
 987      format(f8.6,2x,3(i4,1x))
c
 1001     format(2(f11.8,1x),'ialf=',i3,1x,'ibet=',i3)
          return
          end
c
c  *************************************************************************
c
          subroutine getBBflux(ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     $      wave,visib,projarray,temp,surf,flimbx,flimby,ilaw,rinty,
     &      flum,flux,delphi,delphiedge,iedgestar,iedgehor,rldint,
     &      separ,mmdx,MonteCarlo,isw13,ialfmin,ialfmax,fluxlat,istar,phiarr,
     #      phihor)
c
c  October 11, 1999
c
c  This routine will compute the intensities of each element, given the
c  temperatures (temp(ialf,ibet)) and the input wavelength.  It will then
c  integrate the flux given the visibilities (visib) and surface elements
c  (surf).   
c
c  The projarray contains the cosine mu terms for each element.  The visib
c  array contains the cosine mu terms for each element, except if the point
c  is eclipsed in which case the visib=0.
c 
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.  Then scale the fluxes
c   by (separ*solarrad)**2
c
c   
c   UPDATE JULY 4, 2004
c
c   Add the variable MonteCarlo to the argument list.  If MonteCarlo < 10,
c   then proceed as before.  If Monte Carlo > 10, then the fractional
c   pixels were computed in getvisib via Monte Carlo integration.  In
c   that case, we can skip some steps below.
c
          implicit double precision(a-h,o-z)
c
c
c   Set these to the value of ialphmax,ibetmax
c
          integer tempalf,tempbet
          parameter(tempalf=3000,tempbet=3000,itab=tempalf*tempbet)

          parameter(pie=3.14159265358979323d0)
          dimension visib(ialphmax*ibetmax),delphi(ialphmax*ibetmax),
     $        surf(ialphmax*ibetmax),ibetlim(ialphmax),
     $        temp(ialphmax*ibetmax),flum(ialphmax*ibetmax),
     $        rinty(ialphmax*ibetmax),projarray(ialphmax*ibetmax),
     %        iedgehor(ialphmax*ibetmax),iedgestar(ialphmax*ibetmax),
     $        delphiedge(ialphmax*ibetmax),saveflum(itab),
     #        iflag(itab),mmdx(ialphmax,ibetmax),phiarr(ialphmax*ibetmax)
          dimension xrow(2000),yrow(2000),phihor(ialphmax,2)
c
          if(tempalf.lt.ialphmax)then
            write(*,*)'dimension error in getBBflux ',ialphmax,ibetmax
            stop
          endif
          if(tempbet.lt.ibetmax)then
            write(*,*)'dimension error in getBBflux ',ialphmax,ibetmax
            stop
          endif
c
          dint=pie*(1.0d0-flimbx/3.0d0)
          if(ilaw.eq.2)then
            dint=pie*(1.0d0-flimbx/3.0d0+2.0d0*flimby/9.0d0)
          endif
          if(ilaw.eq.3)then
            dint=pie*(1.0d0-flimbx/3.0d0-flimby/5.0d0)
          endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic law, ilaw=4
c
          if(ilaw.eq.4)then
            dint=pie*(1.0d0-flimbx/3.0d0-flimby/6.0d0)           
          endif

          flux=0.0d0
          C2 = 1.4384d8          ! 1.4384 * 10.**8      ! hc/(k*1e-8)
          C1 = 1.191044d35       ! 2hc^2/((1e-8)**5)
c
c   Initialize the flum matrix.
c
c
          do 2 ialf=1,nalf
            do 1 ibet=1,ibetlim(ialf)        !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              flum(iidx)=0.0d0
              rinty(iidx)=0.0d0
              saveflum(iidx)=0.0d0
              iflag(iidx)=-99
 1          continue
 2        continue
c
c   Compute the intensity values of pixels just behind the edge or the
c   eclipsing horizon for use in fractional eclipse corrections.
c
          wavemu=wave/10000.0d0
          c1=3.74185
          c2=14.3883
c
c   UPDATE JULY 4, 2004
c
c   This loop is no longer needed when using the Monte Carlo
c   routine to compute fractionally eclipsed pixels.
c
          if(MonteCarlo.lt.10)then
            do 4 ialf=1,nalf
              if((istar.eq.1).and.(isw13.gt.0))then
                if((ialf.lt.ialfmin).or.(ialf.gt.ialfmax))go to 4
              endif
              do 3 ibet=1,ibetlim(ialf)
c
cc   UPDATE June 11, 2003
cc
cc   change the 2D arrays into 1D
cc
c                iidx=kount(ialphmax,ialf,ibetlim)+ibet
c   
                iidx=mmdx(ialf,ibet)
                if((iedgehor(iidx).eq.-10).or.(iedgehor(iidx).gt.
     #               5).or.(iedgestar(iidx).eq.-10).or.(iedgestar(iidx)
     $               .gt.5).or.(delphi(iidx).gt.-10.0d0))then
c                  C3 = C2/(WAVE*TEMP(iidx))
c
                  tkkelv=temp(iidx)/1000.0d0
                  C3 = C2/(wavemu*tkkelv)
                  saveflum(iidx)=C1/(dexp(c3)-1.0d0)/wavemu**5
                  dark=(1.0d0-flimbx+flimbx*dabs(projarray(iidx)))
                  if(ilaw.eq.2)dark=dark-flimby*dabs(projarray(iidx))*
     %                  dlog(dabs(projarray(iidx)))
                  if(ilaw.eq.3)dark=dark-flimby*(1.0d0-
     &              dsqrt(dabs(projarray(iidx))))
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic law, ilaw=4
c
                  if(ilaw.eq.4)dark=dark-flimby*(1.0d0-
     &              dabs(projarray(iidx)))**2
                  saveflum(iidx)=saveflum(iidx)*dark
                  saveflum(iidx)=surf(iidx)*saveflum(iidx)*
     $              dabs(projarray(iidx))
                endif
 3            continue
 4          continue
          endif
cc
          sumcor1=0.0d0
          sumcor2=0.0d0
          dtheta=0.5d0*pie/dble(Nalf)
c
          corr1max=-12345.
          corr2max=-12345.
          fmax=-12345.
          rowflux=0.0d0
          DO 10 ialf=1,nalf
            if((istar.eq.1).and.(isw13.gt.0))then
              if((ialf.lt.ialfmin).or.(ialf.gt.ialfmax))go to 10
            endif
            theta=-dtheta+2.0d0*dtheta*dble(ialf)
            sitheta=dsin(theta)
            dphi=pie/dble(ibetlim(ialf))
            DO 9 ibet = 1,ibetlim(ialf)               !4*nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c 
              iidx=mmdx(ialf,ibet)
              corr1=0.0d0
              corr2=0.0d0
              if((projarray(iidx).le.0.0d0))then
                xrow(ibet)=phiarr(iidx)
                yrow(ibet)=0.0d0
                go to 9
              endif
              tkkelv=temp(iidx)/1000.0d0
              C3 = C2/(wavemu*tkkelv)
              flum(iidx)=C1/(dexp(c3)-1.0d0)/wavemu**5
              dark=(1.0d0-flimbx+flimbx*projarray(iidx))
              if(ilaw.eq.2)dark=dark-flimby*projarray(iidx)*
     %               dlog(projarray(iidx))
              if(ilaw.eq.3)dark=dark-flimby*(1.0-dsqrt(projarray(iidx)))
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic law, ilaw=4
c
               if(ilaw.eq.4)dark=dark-flimby*(1.0d0
     #               -dabs(projarray(iidx)))**2
c
              flum(iidx)=flum(iidx)*dark
              rinty(iidx)=flum(iidx) ! save intensities for plotting
              saveflum(iidx)=surf(iidx)*flum(iidx)*
     $           projarray(iidx)
              flum(iidx)=surf(iidx)*flum(iidx)*visib(iidx)
c
c
c   Correct for fractional pixels near the limb.
c
c   UPDATE JULY 4, 2004
c
c   These corrections are no longer needed when using the Monte
c   Carlo routine to compute fractionally eclipsed pixels.
c
 8            iidx=mmdx(ialf,ibet)
              if(MonteCarlo.gt.10)then
                corr1=0.0d0
                corr2=0.0d0
                go to 867
              endif
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
              if(iedgestar(iidx).eq.10)then
                frac=0.5d0*(dabs(delphiedge(iidx))-dphi)/dphi
                if(frac.lt.0.0d0)then
                  corr2=frac*saveflum(iidx)
                  saveflum(iidx)=saveflum(iidx)+corr2
                else
                  if(ibet.lt.ibetlim(ialf))then
c                    iidx=kount(ialphmax,ialf,ibetlim)+ibet+1
                    iidx=mmdx(ialf,ibet+1)
                    corr2=frac*saveflum(iidx)
                    saveflum(iidx)=saveflum(iidx)+corr2
                  else
c                    iidx=kount(ialphmax,ialf,ibetlim)+1
                    iidx=mmdx(ialf,1)
                    corr2=frac*saveflum(iidx)
                    saveflum(iidx)=saveflum(iidx)+corr2
                  endif
                endif
              endif
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
              iidx=mmdx(ialf,ibet)
              if(iedgestar(iidx).eq.20)then
                frac=0.5d0*(dabs(delphiedge(iidx))-dphi)/dphi
                if(frac.lt.0.0d0)then
                  corr2=frac*saveflum(iidx)
                  saveflum(iidx)=saveflum(iidx)+corr2
                else
                  if(ibet.gt.1)then
c                    iidx=kount(ialphmax,ialf,ibetlim)+(ibet-1)
                    iidx=mmdx(ialf,ibet-1)
                    corr2=frac*saveflum(iidx)
                    saveflum(iidx)=saveflum(iidx)+corr2
                  else
c                    iidx=kount(ialphmax,ialf,ibetlim)+ibetlim(ialf)
                    iidx=mmdx(ialf,ibetlim(ialf))
                    corr2=frac*saveflum(iidx)
                    saveflum(iidx)=saveflum(iidx)+corr2
                  endif
                endif
              endif
c
c
c   Check for fractional pixels near the horizon of the star in front,
c   in the beta direction (along constant latitude rows).
c
8677              iidx=mmdx(ialf,ibet)
              if(delphi(iidx).ge.-10.0d0)then
                frac=0.5d0*(dabs(delphi(iidx))-dphi)/dphi
                if(frac.lt.0.0d0)then
                  corr1=frac*saveflum(iidx)
c                  if(ibet.lt.ibetlim(ialf))then 
cc                    jjdx=kount(ialphmax,ialf,ibetlim)+ibet+1   
c                    jjdx=mmdx(ialf,ibet+1)
c                    corr1=frac*saveflum(jjdx)
c                  else
cc                    jjdx=kount(ialphmax,ialf,ibetlim)+1   
c                    jjdx=mmdx(ialf,1)
c                    corr1=frac*saveflum(jjdx)
c                  endif
                else
                  if(iedgehor(iidx).eq.10)then
                    if(ibet.lt.ibetlim(ialf))then 
c                      jjdx=kount(ialphmax,ialf,ibetlim)+ibet+1   
                      jjdx=mmdx(ialf,ibet+1)
                      corr1=frac*saveflum(jjdx)
                    else
c                      jjdx=kount(ialphmax,ialf,ibetlim)+1   
                      jjdx=mmdx(ialf,1)
                      corr1=frac*saveflum(jjdx)
                    endif
                  endif
c                  iidx=kount(ialphmax,ialf,ibetlim)+ibet
                  iidx=mmdx(ialf,ibet)
                  if(iedgehor(iidx).eq.20)then
                    if(ibet.gt.1)then    
c                      jjdx=kount(ialphmax,ialf,ibetlim)+(ibet-1)/dphi   
                      jjdx=mmdx(ialf,ibet-1)
                      corr1=frac*saveflum(jjdx)
                    else
c                      jjdx=kount(ialphmax,ialf,ibetlim)+ibetlim(ialf)   
                      jjdx=mmdx(ialf,ibetlim(ialf))
                      corr1=frac*saveflum(jjdx) 
                    endif
                  endif
                endif
              endif
c
              iidx=mmdx(ialf,ibet)
867           flux=flux+flum(iidx)  +corr1 

              xrow(ibet)=phiarr(iidx)
              yrow(ibet)=flum(iidx)
              flum(iidx)=flum(iidx) +corr1 
c
 9          continue
c            if(phihor(ialf,1).ge.0.0d0)then
c              xrow(ibetlim(ialf)+1)=phihor(ialf,1)
c              yrow(ibetlim(ialf)+1)=0.0d0
c            endif
c            if(phihor(ialf,2).ge.0.0d0)then
c              xrow(ibetlim(ialf)+1)=phihor(ialf,2)
c              yrow(ibetlim(ialf)+1)=0.0d0
c            endif
c
            fred1=phihor(ialf,1)
            fred2=phihor(ialf,2)
            Irow=ibetlim(ialf)
            call edgecor(Irow,xrow,yrow,
     #           fred1,fred2,Nbet,corr1)
            flux=flux+corr1
c
c            if((istar.eq.1).and.(ialf.eq.32))then
c              open(unit=33,file='row.dat',status='unknown')
c              do 888 jjj=1,ibetlim(ialf)
c                write(33,*)xrow(jjj),yrow(jjj)
c888           continue
c                fred=0.0d0
c                if(phihor(ialf,1).ge.0.0d0)write(33,*)phihor(ialf,1),fred
c                if(phihor(ialf,2).ge.0.0d0)write(33,*)phihor(ialf,2),fred
c              close(33)
c            endif
c
c   UPDATE January 10, 2010
c
c   Use a routine to sample the phi vs. flum 
c   curve and integrate using the trapizoid
c   rule
c
c           call binphi(ibetlim(ialf),xrow,yrow,rowsum,Nbet)
c           rowflux=rowflux+rowsum
c
 10       continue
c
c
c   Scale the light curve by the integral of the limb darkening law
c   for compatibility with Wilson-Devinney.
c        
c          if(montecarlo.gt.10)flux=rowflux  

          if((istar.eq.1).and.(isw13.gt.0))flux=flux+fluxlat
          flux=pie*flux/dint
c
c
c   UPDATE April 3, 2002
c
c   Scale the fluxes.
c
          solarrad=6.9598d10
          flux=flux*(separ*solarrad)**2
c
          rldint=dint
          return
          end
c
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine detailrefl(ialphmax1,ibetmax1,
     $      Nalph1,Nbet1,ibetlim1,Nalph2,Nbet2,ibetlim2,ratio1,ratio2,
     $      xarray1,yarray1,zarray1,gradx1,grady1,gradz1,garray1,surf1,
     $      xarray2,yarray2,zarray2,gradx2,grady2,gradz2,garray2,surf2,
     $      temp1,temp2,tempold1,tempold2,dbolx,dboly,ilaw,alb1,alb2,teff1,
     $      teff2,Tgrav1,Tgrav2,rLx,idint,redge,betarim,gpole1,gpole2,
     %      Tpole1,Tpole2,coprat1,coprat2,bdist,mmdx1,mmdx2,
     *      ialphmax2,ibetmax2)
c
c   October 15, 1999
c
c    This routine will alter the temperatures on star istar by
c    means of 'detailed reflection' (R. E. Wilson 1990, ApJ, 356, 613). 
c    The arrays ratio1 and ratio2 store the 'reflection' factors for
c    each surface element.  These can be iterated to any accuracy desired.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xarray1(ialphmax1*ibetmax1),yarray1(ialphmax1*ibetmax1),
     $      zarray1(ialphmax1*ibetmax1),surf1(ialphmax1*ibetmax1),
     &      gradx1(ialphmax1*ibetmax1),temp1(ialphmax1*ibetmax1),
     $      grady1(ialphmax1*ibetmax1),gradz1(ialphmax1*ibetmax1),
     %      garray1(ialphmax1*ibetmax1),ratio1(ialphmax1*ibetmax1),
     $      tempold1(ialphmax1*ibetmax1),ibetlim1(ialphmax1),
     $      ibetlim2(ialphmax2),coprat1(ialphmax1*ibetmax1),
     $      coprat2(ialphmax2*ibetmax2)
          dimension xarray2(ialphmax2*ibetmax2),yarray2(ialphmax2*ibetmax2),
     $      zarray2(ialphmax2*ibetmax2),surf2(ialphmax2*ibetmax2),
     &      gradx2(ialphmax2*ibetmax2),temp2(ialphmax2*ibetmax2),
     $      grady2(ialphmax2*ibetmax2),gradz2(ialphmax2*ibetmax2),
     %      garray2(ialphmax2*ibetmax2),ratio2(ialphmax2*ibetmax2),
     $      tempold2(ialphmax2*ibetmax2),mmdx2(ialphmax2,ibetmax2)
          dimension dbolx(8,2),dboly(8,2),mmdx1(ialphmax1,ibetmax1)
c
c   Start with star 1 and compute the flux from star 2.
c
          darkbolx1=dbolx(1,1)
          darkboly1=dboly(1,1)
          darkbolx2=dbolx(1,2)
          darkboly2=dboly(1,2)
c
          dtheta1=pie/(1.0d0*nalph1)
          dtheta2=pie/(1.0d0*nalph2)
c
c   Define the integrated  limb darkening coefficients.  The equation is
c
c   dint=2*pi*int_0^1{mu*(1-x*(1-mu))d(mu)}  for the linear law, etc.
c
          dint1=pie*(1.0d0-darkbolx1/3.0d0)
          dint2=pie*(1.0d0-darkbolx2/3.0d0)
          if(ilaw.eq.2)then
            dint1=pie*(1.0d0-darkbolx1/3.0d0+2.0d0*darkboly1/9.0d0)
            dint2=pie*(1.0d0-darkbolx2/3.0d0+2.0d0*darkboly2/9.0d0)
          endif
          if(ilaw.eq.3)then
            dint1=pie*(1.0d0-darkbolx1/3.0d0-darkboly1/5.0d0)
            dint2=pie*(1.0d0-darkbolx2/3.0d0-darkboly2/5.0d0)
          endif
c
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic law, ilaw=4
c
          if(ilaw.eq.4)then
            dint1=pie*(1.0d0-darkbolx1/3.0d0-darkboly1/6.0d0)
            dint2=pie*(1.0d0-darkbolx2/3.0d0-darkboly2/6.0d0)
          endif
c
          if(teff2.gt.0.0d0)then
            C1=(Tpole2/Tpole1)**(4)*(dint1/dint2)
            C1=C1*alb1/dint1
          endif
          nalf12=nalph1/2
          nalf22=nalph2/2
          DIV1 = gpole1    ! gravity at the pole
          DIV2 = gpole2    ! gravity at the pole
c
c   If teff2 < 0, then star 2 does not exist (X-ray binary mode usually).
c   In this case, we simply set C1=rLx, which is the ratio of the X-ray
c   luminosity to that of the star 1's optical luminosity.  The coordinates,
c   gradients, etc of star 2 have been adjusted correctly above.  
c
c   For the mode with 2 normal stars where star 2 is reasonably spherical,
c   you can input rLx = [(Teff2/Teff1)**4]*[(surface_area2/surface_area1)],
c   set Teff2=-Teff2, and get basically the same light curve.  This
c   is the "inverse square law" discussed by Wilson.
c
          if(teff2.le.0.0d0)then
            C1=rLx*alb1
            C1=C1/dint1
          endif
          diff1max=-1.0d0
          diff2max=-1.0d0
c
          T4g2=4.0d0*Tgrav2
          T4g1=4.0d0*Tgrav1
          do 10 ialf=1,Nalph1/2  !Nalph1/2,1,-1
            do 9 ibet=1,ibetlim1(ialf)/2  !4*Nbet1
              summ=0.0
              kcount=0
              do 8 i=1,Nalph2       !1,Nalph2
                do 7 j=1,ibetlim2(i)/2
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c                  iidx=(ialf-1)*ibetlim1(ialf)+ibet
c                  iidx=kount(ialphmax,ialf,ibetlim1)+ibet
                  iidx=mmdx1(ialf,ibet)
c                  jjdx=(i-1)*ibetlim2(i)+j
c                  jjdx=kount(ialphmax,i,ibetlim2)+j
                  jjdx=mmdx2(i,j)
c
                  xflip2=bdist-xarray2(jjdx)     ! move the center of the
                  yflip2=-yarray2(jjdx)          ! second star
                  zflip2=zarray2(jjdx)
                  dist1=(xarray1(iidx)-xflip2)**2 +
     %                  (yarray1(iidx)-yflip2)**2 +
     $                  (zarray1(iidx)-zflip2)**2
                  dist1=dsqrt(dist1)
                  term1=(xflip2-xarray1(iidx))*gradx1(iidx)+
     $                  (yflip2-yarray1(iidx))*grady1(iidx)+
     #                  (zflip2-zarray1(iidx))*gradz1(iidx)
                  foreshort1=(term1/(dist1)) 
c

c   UPDATE March 22, 2002
c
c   If Teff < 0, then bail out here.
c
                  if(teff2.lt.0.0d0)go to 99
c
                  if(foreshort1.le.0.0d0)go to 7
c
                  dist2=dist1            
c
                  xflip1=bdist-xarray1(iidx)
                  yflip1=-yarray1(iidx)
                  zflip1=zarray1(iidx)
                  term2=(xflip1-xarray2(jjdx))*gradx2(jjdx)+
     $                  (yflip1-yarray2(jjdx))*grady2(jjdx)+
     #                  (zflip1-zarray2(jjdx))*gradz2(jjdx)
                  foreshort2=(term2/(dist2))  
c 
                  if(foreshort2.le.-0.05d0)go to 8   ! leave this beta loop
                  if(foreshort2.le.0.0d0)go to 7
                  if(teff2.le.0.0d0)go to 7
c
c                  term3=(garray2(i,j)/div2)**(T4g2)
c                  term3=term3*surf2(i,j)*coprat2(i,j)
c                  term3=term3*foreshort1*foreshort2/(dist2*dist2)
c
                  term3=(garray2(jjdx)/div2)**(T4g2)*
     $     surf2(jjdx)*coprat2(jjdx)*foreshort1*foreshort2/(dist2*dist2)
c
                  if(ilaw.eq.1)then
                    term3=term3*(1.0d0-darkbolx2+darkbolx2*foreshort2)
                  endif
                  if(ilaw.eq.2)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly2*foreshort2*dlog(foreshort2)
                    else
                      ttt=0.0d0
                    endif
                    term3=term3*(1.0d0-darkbolx2*(1.0d0-foreshort2)-ttt)
                  endif
                  if(ilaw.eq.3)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly2*(1.0-dsqrt(foreshort2))
                    else
                      ttt=darkboly2
                    endif
                    term3=term3*(1.0d0-darkbolx2*(1.0d0-foreshort2)-ttt)
                  endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
                  if(ilaw.eq.4)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly2*(1.0-(foreshort2))**2
                    else
                      ttt=darkboly2
                    endif
                    term3=term3*(1.0d0-darkbolx2*(1.0d0-foreshort2)-ttt)
                  endif

c
c   If there is a disk, we need to check which lines of sight are blocked
c
                  if(idint.ge.1)then
                    zrim=redge*dtan(betarim*0.017453293d0)
                    xA=xarray1(iidx)
                    yA=yarray1(iidx)
                    zA=zarray1(iidx)
                    xB=bdist-xarray2(jjdx)     
                    yB=-yarray2(jjdx)          
                    zB=zarray2(jjdx)           
                    call zheight(xA,yA,zA,xB,yB,zB,redge,zcross)
                    if((zB.ge.0.0d0).and.(zcross.lt.zrim))term3=0.0d0
                    if((zB.lt.0.0d0).and.(zcross.gt.-zrim))term3=0.0d0
                  endif
                  summ=summ+term3
 7              continue
 8            continue
c
 887          do 888 i=1,Nalph2
                do 777 j=ibetlim2(i),ibetlim2(i)/2+1,-1
c                  iidx=(ialf-1)*ibetlim1(ialf)+ibet
c                  jjdx=(i-1)*ibetlim2(i)+j
c                  iidx=kount(ialphmax,ialf,ibetlim1)+ibet
c                  jjdx=kount(ialphmax,i,ibetlim2)+j
                  iidx=mmdx1(ialf,ibet)
                  jjdx=mmdx2(i,j)
                  xflip2=bdist-xarray2(jjdx)     ! move the center of the
                  yflip2=-yarray2(jjdx)          ! second star
                  zflip2=zarray2(jjdx)
                  dist1=(xarray1(iidx)-xflip2)**2 +
     %                  (yarray1(iidx)-yflip2)**2 +
     $                  (zarray1(iidx)-zflip2)**2
                  dist1=dsqrt(dist1)
                  term1=(xflip2-xarray1(iidx))*gradx1(iidx)+
     $                  (yflip2-yarray1(iidx))*grady1(iidx)+
     #                  (zflip2-zarray1(iidx))*gradz1(iidx)
                  foreshort1=(term1/(dist1)) 
c
                  if(foreshort1.le.0.0d0)go to 777
c
                  dist2=dist1            
c
                  xflip1=bdist-xarray1(iidx)
                  yflip1=-yarray1(iidx)
                  zflip1=zarray1(iidx)
                  term2=(xflip1-xarray2(jjdx))*gradx2(jjdx)+
     $                  (yflip1-yarray2(jjdx))*grady2(jjdx)+
     #                  (zflip1-zarray2(jjdx))*gradz2(jjdx)
                  foreshort2=(term2/(dist2))  
c
                  if(foreshort2.le.0.0d0)go to 888
                  if(teff2.le.0.0d0)go to 777
c
c                  term3=(garray2(i,j)/div2)**(T4g2)
c                  term3=term3*surf2(i,j)*coprat2(i,j)
c                  term3=term3*foreshort1*foreshort2/(dist2*dist2)
c
                  term3=(garray2(jjdx)/div2)**(T4g2)*
     %    surf2(jjdx)*coprat2(jjdx)*foreshort1*foreshort2/(dist2*dist2)
c
                  if(ilaw.eq.1)then
                    term3=term3*(1.0d0-darkbolx2+darkbolx2*foreshort2)
                  endif
                  if(ilaw.eq.2)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly2*foreshort2*dlog(foreshort2)
                    else
                      ttt=0.0d0
                    endif
                    term3=term3*(1.0d0-darkbolx2*(1.0d0-foreshort2)-ttt)
                  endif
                  if(ilaw.eq.3)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly2*(1.0-dsqrt(foreshort2))
                    else
                      ttt=darkboly2
                    endif
                    term3=term3*(1.0d0-darkbolx2*(1.0d0-foreshort2)-ttt)
                  endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
                  if(ilaw.eq.4)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly2*(1.0-(foreshort2))**2
                    else
                      ttt=darkboly2
                    endif
                    term3=term3*(1.0d0-darkbolx2*(1.0d0-foreshort2)-ttt)
                  endif

c
c   If there is a disk, we need to check which lines of sight are blocked
c
                  if(idint.ge.1)then
c                    iidx=(ialf-1)*ibetlim1(ialf)+ibet
c                    jjdx=(i-1)*ibetlim2(i)+j
                    zrim=redge*dtan(betarim*0.017453293d0)
                    xA=xarray1(iidx)
                    yA=yarray1(iidx)
                    zA=zarray1(iidx)
                    xB=bdist-xarray2(jjdx)     
                    yB=-yarray2(jjdx)          
                    zB=zarray2(jjdx)           
                    call zheight(xA,yA,zA,xB,yB,zB,redge,zcross)
                    if((zB.ge.0.0d0).and.(zcross.lt.zrim))term3=0.0d0
                    if((zB.lt.0.0d0).and.(zcross.gt.-zrim))term3=0.0d0
                  endif
                  summ=summ+term3
 777            continue
 888          continue

c
c   Check to see if the point on star 1 can actually see any point on
c   star 2.  If not, then the summ will be zero.
c
c   UPDATE March 22, 2002
c
c   Delete this statement
c
c              if((teff2.le.0.0d0).and.(idint.lt.1))summ=foreshort1 ! no disk
c
              if(summ.eq.0.0d0)go to 10  
c
c   If Teff2 < 0, we assume star 2 is a point source.  Thus the integration
c   over the surface of star 2 reduces to the foreshorting angle of
c   the specific ialf-ibet element on star 1.  If there is a disk we have
c   to compute which points on star A are shielded by the disk rim.
c
c
c   UPDATE March 22, 2002
c
c   Comment this out, and move statement 99 one line down
c
c 99           if((teff2.le.0.0d0).and.(idint.lt.1))summ=foreshort1 ! no disk
 99           if((teff2.le.0.0d0).and.(idint.ge.1))then
                zrim=redge*dtan(betarim*0.017453293d0)
                xA=xarray1(iidx)
                yA=yarray1(iidx)
                zA=zarray1(iidx)
                xB=bdist
                yB=0.0d0
                zB=0.0d0
                call zheight(xA,yA,zA,xB,yB,zB,redge,zcross)
c
c   UPDATE March 22, 2002
c
c   Comment this out.
c
c                summ=foreshort1
                if(summ.lt.0.0d0)summ=0.0d0
                if((zA.ge.0.0d0).and.(zcross.lt.zrim))summ=0.0d0
                if((zA.lt.0.0d0).and.(zcross.gt.-zrim))summ=0.0d0
              endif
c                    
              FpBoverFA=C1/((garray1(iidx)/div1)**(T4g1))*summ
              ratio1(iidx)=1.0d0+FpBoverFA
              tnew=tempold1(iidx)*ratio1(iidx)**0.25d0
              diff=tnew-tempold1(iidx)
              if(diff.gt.diff1max)diff1max=diff
              temp1(iidx)=tnew
              if(summ.gt.0.0d0)jcount=jcount+1
 9          continue
 10       continue
c
c   Now go to star 2 and compute the irradation from star 1
c
 11       if(teff2.le.0.0d0)go to 50
          C2=(Tpole1/Tpole2)**(4)*(dint2/dint1)
          C2=C2*alb2/dint2
c
          do 20 ialf=1,Nalph2/2  !Nalph2/2,1,-1
            do 19 ibet=1,ibetlim2(ialf)/2   !Nbet2
              summ=0.0d0
              do 18 i=1,Nalph1
                do 17 j=1,ibetlim1(i)/2
c                  iidx=(ialf-1)*ibetlim2(ialf)+ibet
c                  jjdx=(i-1)*ibetlim1(i)+j
c                  iidx=kount(ialphmax,ialf,ibetlim2)+ibet
c                  jjdx=kount(ialphmax,i,ibetlim1)+j
c
                  iidx=mmdx2(ialf,ibet)
                  jjdx=mmdx1(i,j)
                  xflip1=bdist-xarray1(jjdx)     ! move the center of the
                  yflip1=-yarray1(jjdx)        ! first star
                  zflip1=zarray1(jjdx)
                  dist1=(xflip1-xarray2(iidx))**2 +
     %                  (yflip1-yarray2(iidx))**2 +
     $                  (zflip1-zarray2(iidx))**2
                  dist1=dsqrt(dist1)
                  term1=(xflip1-xarray2(iidx))*gradx2(iidx)+
     $                  (yflip1-yarray2(iidx))*grady2(iidx)+
     #                  (zflip1-zarray2(iidx))*gradz2(iidx)
                  foreshort1=(term1/(dist1))    
                  if(foreshort1.lt.0.0d0)go to 17
c
                  dist2=dist1  
                  
                  xflip2=bdist-xarray2(iidx)
                  yflip2=-yarray2(iidx)
                  zflip2=zarray2(iidx)
c
                  term2=(xflip2-xarray1(jjdx))*gradx1(jjdx)+
     $                  (yflip2-yarray1(jjdx))*grady1(jjdx)+
     #                  (zflip2-zarray1(jjdx))*gradz1(jjdx)
                  foreshort2=(term2/(dist2)) 
                  if(foreshort2.lt.0.0d0)go to 18
c
c                  term3=(garray1(i,j)/div1)**(T4g1)
c                  term3=term3*surf1(i,j)*coprat1(i,j)
c                  term3=term3*foreshort1*foreshort2/(dist2*dist2)
                  term3=(garray1(jjdx)/div1)**(T4g1)*
     &     surf1(jjdx)*coprat1(jjdx)*foreshort1*foreshort2/(dist2*dist2)
c
                  if(ilaw.eq.1)then
                    term3=term3*(1.0d0-darkbolx1+darkbolx1*foreshort2)
                  endif
                  if(ilaw.eq.2)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly1*foreshort2*dlog(foreshort2)
                    else
                      ttt=0.0d0
                    endif
                    term3=term3*(1.0d0-darkbolx1*(1.0d0-foreshort2)-ttt)
                  endif
                  if(ilaw.eq.3)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly1*(1.0d0-dsqrt(foreshort2))
                    else
                      ttt=darkboly1
                    endif
                    term3=term3*(1.0d0-darkbolx1*(1.0d0-foreshort2)-ttt)
                  endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
                  if(ilaw.eq.4)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly1*(1.0d0-(foreshort2))**2
                    else
                      ttt=darkboly1
                    endif
                    term3=term3*(1.0d0-darkbolx1*(1.0d0-foreshort2)-ttt)
                  endif
c
c   If there is a disk, we need to check which lines of sight are blocked.
c
                  if(idint.ge.1)then
                    zrim=redge*dtan(betarim*0.017453293d0)
                    xA=xarray1(jjdx)
                    yA=yarray1(jjdx)
                    zA=zarray1(jjdx)
                    xB=bdist-xarray2(iidx) ! move the center of the
                    yB=-yarray2(iidx)    ! second star to be consistent
                    zB=zarray2(iidx)     ! with the disk coordinates
                    call zheight(xA,yA,zA,xB,yB,zB,redge,zcross)
                    if((zB.ge.0.0).and.(zcross.lt.zrim))term3=0.0
                    if((zB.lt.0.0).and.(zcross.gt.-zrim))term3=0.0
                  endif
                  summ=summ+term3
 17             continue
 18           continue
c
              do 1888 i=1,Nalph1
                do 1777 j=ibetlim1(i),ibetlim1(i)/2+1,-1  
c                  iidx=(ialf-1)*ibetlim2(ialf)+ibet
c                  jjdx=(i-1)*ibetlim1(i)+j
c                  iidx=kount(ialphmax,ialf,ibetlim2)+ibet
c                  jjdx=kount(ialphmax,i,ibetlim1)+j
                  iidx=mmdx2(ialf,ibet)
                  jjdx=mmdx1(i,j)
                  xflip1=bdist-xarray1(jjdx)     ! move the center of the
                  yflip1=-yarray1(jjdx)        ! first star
                  zflip1=zarray1(jjdx)
                  dist1=(xflip1-xarray2(iidx))**2 +
     %                  (yflip1-yarray2(iidx))**2 +
     $                  (zflip1-zarray2(iidx))**2
                  dist1=dsqrt(dist1)
                  term1=(xflip1-xarray2(iidx))*gradx2(iidx)+
     $                  (yflip1-yarray2(iidx))*grady2(iidx)+
     #                  (zflip1-zarray2(iidx))*gradz2(iidx)
                  foreshort1=(term1/(dist1))    
                  if(foreshort1.lt.0.0d0)go to 1777
c
                  dist2=dist1  
                  
                  xflip2=bdist-xarray2(iidx)
                  yflip2=-yarray2(iidx)
                  zflip2=zarray2(iidx)
c
                  term2=(xflip2-xarray1(jjdx))*gradx1(jjdx)+
     $                  (yflip2-yarray1(jjdx))*grady1(jjdx)+
     #                  (zflip2-zarray1(jjdx))*gradz1(jjdx)
                  foreshort2=(term2/(dist2)) 
                  if(foreshort2.lt.0.0d0)go to 1888
c
c                  term3=(garray1(jjdx)/div1)**(T4g1)
c                  term3=term3*surf1(jjdx)*coprat1(jjdx)
c                  term3=term3*foreshort1*foreshort2/(dist2*dist2)
                  term3=(garray1(jjdx)/div1)**(T4g1)*
     &     surf1(jjdx)*coprat1(jjdx)*foreshort1*foreshort2/(dist2*dist2)
c
                  if(ilaw.eq.1)then
                    term3=term3*(1.0d0-darkbolx1+darkbolx1*foreshort2)
                  endif
                  if(ilaw.eq.2)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly1*foreshort2*dlog(foreshort2)
                    else
                      ttt=0.0d0
                    endif
                    term3=term3*(1.0d0-darkbolx1*(1.0d0-foreshort2)-ttt)
                  endif
                  if(ilaw.eq.3)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly1*(1.0d0-dsqrt(foreshort2))
                    else
                      ttt=darkboly1
                    endif
                    term3=term3*(1.0d0-darkbolx1*(1.0d0-foreshort2)-ttt)
                  endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
                  if(ilaw.eq.4)then
                    if(foreshort2.gt.0.0d0)then
                      ttt=darkboly1*(1.0d0-(foreshort2))**2
                    else
                      ttt=darkboly1
                    endif
                    term3=term3*(1.0d0-darkbolx1*(1.0d0-foreshort2)-ttt)
                  endif
c
c   If there is a disk, we need to check which lines of sight are blocked.
c
                  if(idint.ge.1)then
                    zrim=redge*dtan(betarim*0.017453293d0)
                    xA=xarray1(jjdx)
                    yA=yarray1(jjdx)
                    zA=zarray1(jjdx)
                    xB=bdist-xarray2(iidx) ! move the center of the
                    yB=-yarray2(iidx)    ! second star to be consistent
                    zB=zarray2(iidx)     ! with the disk coordinates
                    call zheight(xA,yA,zA,xB,yB,zB,redge,zcross)
                    if((zB.ge.0.0).and.(zcross.lt.zrim))term3=0.0
                    if((zB.lt.0.0).and.(zcross.gt.-zrim))term3=0.0
                  endif
                  summ=summ+term3
 1777           continue
 1888         continue
c
c   Check to see if the point on star 2 can actually see any point on
c   star 1.  If not, then the summ will be zero.
c
              if(summ.eq.0.0d0)go to 20  
c
              FpAoverFB=C2/((garray2(iidx)/div2)**(T4g2))*summ
              ratio2(iidx)=1.0d0+FpAoverFB
              tnew=tempold2(iidx)*ratio2(iidx)**0.25d0
              diff=tnew-tempold2(iidx)
              if(diff.gt.diff2max)diff2max=diff
              temp2(iidx)=tnew
 19         continue
 20       continue
c
c   Now use symmetry to fill in the other quadrants on the star.  
c
 50       continue
c
          DO 401 IALF = 1, nalph1/2
            DO 400 IBET = 1, ibetlim1(ialf)/2
              I1=nalph1-(ialf-1)
              J2=ibetlim1(ialf)-(ibet-1)
c              iidx=(ialf-1)*ibetlim1(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim1)+ibet
               iidx=mmdx1(ialf,ibet)
c
              izz=ialf
              jzz=j2
c              jjdx=(izz-1)*ibetlim1(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim1)+jzz
              jjdx=mmdx1(izz,jzz)
              ratio1(jjdx)=ratio1(iidx)
              temp1(jjdx)=temp1(iidx)

              izz=I1
              jzz=ibet
c              jjdx=(izz-1)*ibetlim1(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim1)+jzz
              jjdx=mmdx1(izz,jzz)
              ratio1(jjdx)=ratio1(iidx)
              temp1(jjdx)=temp1(iidx)
c
              izz=I1
              jzz=j2
c              jjdx=(izz-1)*ibetlim1(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim1)+jzz
              jjdx=mmdx1(izz,jzz)
              ratio1(jjdx)=ratio1(iidx)
              temp1(jjdx)=temp1(iidx)
c
 400        continue
401       CONTINUE
c
          if(teff2.le.0.0)go to 999
          DO 501 IALF = 1, nalph2/2
            DO 500 IBET = 1, ibetlim2(ialf)/2
              I1=nalph1-(ialf-1)
              J2=ibetlim2(ialf)-(ibet-1)
c              iidx=(ialf-1)*ibetlim2(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim2)+ibet
              iidx=mmdx2(ialf,ibet)
              I1=nalph2-(ialf-1)
              J2=ibetlim2(ialf)-(ibet-1)

              izz=ialf
              jzz=j2
c              jjdx=(izz-1)*ibetlim2(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim2)+jzz
              jjdx=mmdx2(izz,jzz)
              ratio2(jjdx)=ratio2(iidx)
              temp2(jjdx)=temp2(iidx)
c
              izz=I1
              jzz=ibet
c              jjdx=(izz-1)*ibetlim2(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim2)+jzz
              jjdx=mmdx2(izz,jzz)
              ratio2(jjdx)=ratio2(iidx)
              temp2(jjdx)=temp2(iidx)

              izz=I1
              jzz=j2
c              jjdx=(izz-1)*ibetlim2(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim2)+jzz
              jjdx=mmdx2(izz,jzz)
              ratio2(jjdx)=ratio2(iidx)
              temp2(jjdx)=temp2(iidx)
c
 500        continue
501       CONTINUE
c
c   Update the ratio arrays.
c
          call copytemp(ialphmax1,ibetmax1,Nalph1,Nbet1,ratio1,coprat1,mmdx1,
     $          ibetlim1)
          call copytemp(ialphmax2,ibetmax2,Nalph2,Nbet2,ratio2,coprat2,mmdx2,
     $          ibetlim2)
c
          if(diff2max.lt.0.0d0)diff2max=0.0
          if(diff1max.lt.0.0d0)diff1max=0.0

          write(2,72)diff2max
 999      write(2,71)diff1max
c
 71       format(/'maximum temperature change for star 1 = ',f15.9)
 72       format('maximum temperature change for star 2 = ',f15.9)

          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine initratio(ialphmax1,ibetmax1,
     #          Nalph1,Nbet1,Nalph2,Nbet2,ratio1,ratio2,
     &          coprat1,coprat2,ialphmax2,ibetmax2)
c
c   October 15, 1999
c
c   Initialize the ratio arrays to 1.0.
c
          implicit double precision(a-h,o-z)
c
          dimension ratio1(ialphmax1*ibetmax1),ratio2(ialphmax2*ibetmax2)
          dimension coprat1(ialphmax1*ibetmax1),coprat2(ialphmax2*ibetmax2)
c
          do 10 ialf=1,Nalph1
            do 9 ibet=1,4*Nbet1
              iidx=(ialf-1)*4*Nbet1+ibet
              ratio1(iidx)=1.0d0
              coprat1(iidx)=1.0d0
 9          continue
 10       continue
c
          do 20 ialf=1,Nalph2
            do 19 ibet=1,4*Nbet2
              iidx=(ialf-1)*4*Nbet2+ibet
              ratio2(iidx)=1.0d0
              coprat2(iidx)=1.0d0
 19         continue
 20       continue
c
          return
          end
c
c   
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     
          subroutine copytemp(ialphmax,ibetmax,Nalph,Nbet,tempold,
     &      tempnew,mmdx,ibetlim)
c
c   October 18, 1999
c
c   Copy the temperature array into a new array.  We need to keep
c   the original temperatures for the reflection effect routine
c
          implicit double precision(a-h,o-z)
c
          dimension tempold(ialphmax*ibetmax),tempnew(ialphmax*ibetmax)
          dimension mmdx(ialphmax,ibetmax),ibetlim(ialphmax)
c
          do 9 ialf=1,nalph
            do 10 ibet=1,ibetlim(ialf)
c              iidx=(ialf-1)*4*Nbet+ibet
              iidx=mmdx(ialf,ibet)
              tempnew(iidx)=tempold(iidx)
 10         continue
 9        continue
c
          return
          end
c
c   **********************************
c
          subroutine clip(Nhoriz,xhoriz,yhoriz,xvis,yvis,xhid,yhid,
     #       xedge,yedge)
c
c   October 21, 1999
c
c   This routine assumes that the point xvis,yvis is outside the polygon
c   (xhoriz,yhoriz) and the point xhid,yhid is inside the horizon.  It
c   then computes the coordinate of the point along the line joining the
c   two points that intersects the horizon.
c
          implicit double precision(a-h,o-z)
c
          dimension xhoriz(Nhoriz),yhoriz(Nhoriz)
c
c   Find the slope
c
          xvsave=xvis
          yvsave=yvis
          xhsave=xhid
          yhsave=yhid
          rise=yvis-yhid
          run=xvis-xhid
          if(run.eq.0)then       !special case
            xedge=xvis
            do 10 i=1,20
              iyes=-100
              yedge=(yvis+yhid)/2.0d0
              call insidecircle(Nhoriz,xhoriz,yhoriz,xedge,yedge,iyes,icut)
              if(iyes.eq.100)then       ! inside horizon
                yhid=yedge
              else                      ! outside horizon
                yvis=yedge
              endif
 10         continue
c
          go to 99
          endif
c
          slope=rise/run
          do 20 i=1,20
            iyes=-100
            xedge=(xvis+xhid)/2.0d0
            yedge=slope*(xedge-xhsave)+yhsave
            call insidecircle(Nhoriz,xhoriz,yhoriz,xedge,yedge,iyes,icut)
            if(iyes.eq.100)then       ! inside horizon
              xhid=xedge
            else                      ! outside horizon
              xvis=xedge
            endif
 20       continue
c
 99       yvis=yvsave
          xvis=xvsave
          xhid=xhsave
          yhid=yhsave
c
 100      format(6(f9.6,1x))
          return
          end
c
c  *************************************************************************
c
          subroutine hidgrid(istar,
     %      ialphmax,ibetmax,Nalf,Nbet,ibetlim,phase,finc,Q,
     $      xarray,yarray,zarray,
     $      xend,visib,projarray,garray,gscale,surf,
     $      Nhoriz,xhoriz,yhoriz,rinty,extension,
     $      separation,flux,reff,iecheck,tarray,Nh,xh,yh,bdist,mmdx)
c
c   January 14, 2000
c
c   This routine will output files used for various external plotting
c   packages.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension visib(ialphmax*ibetmax),ibetlim(ialphmax),
     $        xhoriz(Nhoriz),yhoriz(Nhoriz),xend(4),xarray(ialphmax*ibetmax),
     $        yarray(ialphmax*ibetmax),zarray(ialphmax*ibetmax),
     $        rinty(ialphmax*ibetmax),garray(ialphmax*ibetmax),
     $        projarray(ialphmax*ibetmax),surf(ialphmax*ibetmax),
     &        tarray(ialphmax*ibetmax),savex(20000),savey(20000),
     &        xh(Nh),yh(Nh),conex(5),coney(5),mmdx(ialphmax,ibetmax)
c
          character*9 extension
c
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c          
          if(istar.eq.1)open(unit=40,
     %         file='star1inty.'//extension,status='unknown')
          if(istar.eq.2)open(unit=40,
     %         file='star2inty.'//extension,status='unknown')

c
c   Check to see of the star in question is in front.  
c
          infront=0
          if((istar.eq.1).and.((phase.ge.0.0d0).
     #        and.(phase.lt.90.0d0)))infront=1
          if((istar.eq.1).and.((phase.ge.270.0d0).
     %        and.(phase.le.360.0d0)))infront=1
          if((istar.eq.2).and.((phase.ge.0.0d0).
     #        and.(phase.lt.90.0d0)))infront=1
          if((istar.eq.2).and.((phase.ge.270.0d0).
     $        and.(phase.le.360.0d0)))infront=1
c
c   Find the sky coordinates of the center of mass of the star.  This
c   will be recorded as the first line of the star?inty.???.?? file.
c
          if(istar.eq.1)then
            xx=0.0d0
            yy=0.0d0
            zz=0.0d0
            ddphase=phase
          endif
          if(istar.eq.2)then
            xx=0.0d0                     !FIX JULY 21, 2004 (was 1.0)
            yy=0.0d0
            zz=0.0d0
            ddphase=dmod(phase+180.0d0,360.0d0)  ! this is for star 2
          endif

          xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)    ! projected coords
          yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)

          write(40,68)xp*separation,yp*separation,flux,reff*separation,
     %       separation

          isave=0
          DO 501 IALF = 1, NALF
            DO 502 IBET = 1,ibetlim(ialf)        !4*NBET
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
              iidx=mmdx(ialf,ibet)
c
c
c   UPDATE March 4, 2010
c
c   Record the radius of each point.
c
c
              RRRR=sqrt(xarray(iidx)**2+yarray(iidx)**2+zarray(iidx)**2)
              iv1=1
              IF (projarray(iidx).le.0.0d0) go to 502 ! is the surface 
              xx=xarray(iidx)                   ! element visible?
              yy=yarray(iidx)
              zz=zarray(iidx)
              xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)   ! projected coords
              yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              if(infront.ne.1)then    
                iyes=-100
                iv1=1
                if(iecheck.ge.0) 
     $             call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                if((iyes.eq.100))iv1=0    ! point could be visible
              endif                       ! but is eclipsed
c       
c   Record the x,y,z coordinates of the nearby points.  These points
c   will be used for area filling
c
              xx1=xp
              yy1=yp
c
              isave=isave+1
              savex(isave)=xp*separation
              savey(isave)=yp*separation
c
              if(ibet.gt.1)then
                izz=ialf
                jzz=ibet-1
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                xx=xarray(iidx)
                yy=yarray(iidx)
                zz=zarray(iidx)
              else
                izz=ialf
                jzz=ibetlim(ialf)
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                xx=xarray(iidx)
                yy=yarray(iidx)
                zz=zarray(iidx)
              endif
              xx2=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy2=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              if(ialf.gt.1)then
                if(ibet.gt.1)then
                  izz=ialf-1
                  jzz=min(ibet-1,ibetlim(ialf-1))
c                  iidx=(izz-1)*ibetlim(ialf)+jzz
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xx=xarray(iidx)
                  yy=yarray(iidx)
                  zz=zarray(iidx)
                else
                  izz=ialf-1
                  jzz=ibetlim(ialf-1)
c                  iidx=(izz-1)*ibetlim(ialf)+jzz
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xx=xarray(iidx)
                  yy=yarray(iidx)
                  zz=zarray(iidx)
                endif
              else
                xx=0.0d0
                yy=0.0d0
                zz=xend(1)
              endif
              xx3=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy3=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              if(ialf.gt.1)then
                izz=ialf-1
                jzz=min(ibet,ibetlim(ialf-1))
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                xx=xarray(iidx)
                yy=yarray(iidx)
                zz=zarray(iidx)
              else
                xx=0.0d0
                yy=0.0d0
                zz=xend(1)
              endif
              xx4=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy4=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c   Check the corners for eclipsed points.
c
              iv2=1
              iv3=1
              iv4=1
              if(infront.ne.1)then
                iyes=-100
                iv2=1
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx2,yy2,iyes,icut)
                if((iyes.eq.100))iv2=0      
                iv3=1
                iyes=-100
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx3,yy3,iyes,icut)
                if(iyes.eq.100)iv3=0      
                iv4=1
                iyes=-100
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx4,yy4,iyes,icut)
                if(iyes.eq.100)iv4=0      
              endif
              if(iecheck.le.-1)then
                iv1=1
                iv2=1
                iv3=1
                iv4=1
              endif
c
c   There are 13 possibilities for which corners were hidden.  Do each
c   case separately.
c
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c              
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     %            xx4new*separation,yy4new*separation,  
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3*separation,yy3*separation
                go to 99
              endif
c              
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=5
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.1))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
c
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=4
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c 
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                ncorner=5
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     $            xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx3*separation,yy3*separation,
     $            xx4*separation,yy4*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=5
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                ncorner=5
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4new*separation,yy4new*separation,
     &            xx1new*separation,yy1new*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                ncorner=4
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     %             xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=4
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=4
                izz=ialf
                jzz=ibet
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     %             xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
cc   Add an extra step for the back end.
cc
 99           continue
c              endif   ! if Nalf = 1
502         CONTINUE
 501      CONTINUE               
c
c   Add the L1 point for plotting if phase=90.0 or 180.0
c
          if(phase.eq.90.0)then
            do 600 i=1,Nh
              isave=isave+1
              savex(isave)=xh(i)*separation
              savey(isave)=yh(i)*separation
 600        continue
c
            call sort3(isave,savex,savey,savey)
c
c  We want the 5 points with the smallest x-coordinates
c
            do 601 i=1,5
              conex(i)=savex(i)
              coney(i)=savey(i)
 601        continue
c
            call sortcircle(5,conex,coney)
c
            xx1=conex(1)
            xx2=conex(2)
            xx3=conex(3)
            xx4=conex(4)
            xx5=conex(5)
c
            yy1=coney(1)
            yy2=coney(2)
            yy3=coney(3)
            yy4=coney(4)
            yy5=coney(5)
c
            ialf=Nalf/2
            ibet=1
            ncorner=5
c            iidx=(ialf-1)*ibetlim(ialf)+ibet
c            iidx=kount(ialphmax,ialf,ibetlim)+ibet
            iidx=mmdx(ialf,ibet)
            write(40,69)rinty(iidx),ncorner,projarray(iidx),
     %            surf(iidx),dlog10(gscale*garray(iidx)),
     %            ialf,ibet,tarray(iidx),RRRR,
     %            xx1,yy1,
     %            xx2,yy2,
     %            xx3,yy3,
     %            xx4,yy4,
     %            xx5,yy5
          endif
          
          
 68       format(2(f14.8,3x),e16.9,f10.4,2x,f10.4)
 69       format(e16.9,1x,i3,1x,f9.7,1x,e12.6,1x,f8.5,1x,
     %       2(i3,1x),1x,f10.3,2x,f11.6/,10(f11.6,1x))
c
          close(40)
c
          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine getextension(phase,extension)
c
c   October 25, 1999
c
c   extension = '001.00'  for phase = 1.0
c   extension = '011.22'  for phase = 11.22
c   extension = '123.45'  for phase = 123.45  etc.
c
          implicit double precision(a-h,o-z)
c
          character*9 extension
c
          if((phase.ge.0.0d0).and.(phase.lt.10.0d0))write(extension,100)phase
          if((phase.ge.10.0d0).and.(phase.lt.100.0d0))write(extension,101)phase
          if((phase.ge.100.0d0).and.(phase.lt.1000.0d0))
     $         write(extension,102)phase

 100      format('00',f7.5)
 101      format('0',f8.5)
 102      format(f9.5)
c
          return
          end
c
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine checkinput(ialphmax1,ibetmax1,Nthetamax,Nrmax,
     &       Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,Ntheta,Nradius,darkbol1,darkbol2,alb1,alb2,
     $       Nref,rLx,Period,fm,separ,
     $       idraw,iecheck,idint,iatm,ism1,icnU,icnB,icnV,icnR,icnI,
     #       icnJ,icnH,icnK,ivelout,iXout,ecc,pshift,sw5,
     %       sw9,ialphmax2,ibetmax2,ecosw,argper,temprat)
c
c   November 1, 1999
c
c   This routine will check the input parameter values and flag any illegal
c   entries (Nalf should be positive, etc.).
c
c   UPDATE April 8, 2002
c
c   Add the variable sw5 to the argument list.  If teff2 < 0 and
c   isw5 > 0, then the separation will be set as follows:
c
c   rkns = 2*pie*sw5*c/Period   ! K-velocity of pulsar, if sw5 is
c                                  projected semimajor axis in seconds
c
c          velamp=(Q+1.0d0)*rkns
c          sifinc=dsin(finc*3.141592653589793d0/180.0d0)
c          p=period*24.0d0*3600.0d0
c          a=velamp/(2.0d0*3.141592653588783d0)*p/sifinc
c          separ=a/6.959d5
c
c
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)

           character*1 bell
c
           bell=char(7)
           if(Nalph1.lt.1)then
             write(*,96)bell
             stop
           endif
c
           if(Nalph2.lt.1)then
             write(*,97)bell
             stop
           endif
c
           if(Nbet1.lt.1)then
             write(*,98)bell
             stop
           endif
c
           if(Nbet2.lt.1)then
             write(*,99)bell
             stop
           endif
c
           if(Nalph1.lt.10)then
             write(*,92)bell
           endif
c
           if(Nalph2.lt.10)then
             write(*,93)bell
           endif
c
           if(Nbet1.lt.5)then
             write(*,94)bell
           endif
c
           if(Nbet2.lt.5)then
             write(*,95)bell
           endif
c
           bell=char(7)
           if(Nalph1.gt.ialphmax1)then
             write(*,100)bell,ialphmax1
             stop
           endif
c
           if(Nalph2.gt.ialphmax2)then
             write(*,101)bell,ialphmax2
             stop
           endif
c
           if(4*Nbet1.gt.ibetmax1)then
             write(*,102)bell,ibetmax1
             stop
           endif
c
           if(4*Nbet2.gt.ibetmax2)then
             write(*,103)bell,ibetmax2
             stop
           endif
c
           if(fill1.gt.1.0d0)then
             write(*,104)bell
             stop
           endif
c
           if(fill2.gt.1.0d0)then
             write(*,104)bell
             stop
           endif
c
           if(omega1.le.0.0d0)then
             write(*,106)bell
             omega1=dabs(omega1) !stop
           endif           
c
           if(omega2.le.0.0d0)then
             write(*,107)bell
             stop
           endif    
c
           if((dphase.le.0.0d0).or.(dphase.gt.360.0d0))then
             write(*,108)bell
             stop
           endif    
c
           if(Q.le.0.0)then
             write(*,109)bell
             stop
           endif
c
           if((finc.lt.0.0d0).or.(finc.gt.90.0d0))then
             write(*,110)bell
             stop
           endif
c
           if(Teff1.lt.0.0d0)then
             write(*,111)bell
             stop
           endif
c
           if(Teff2.lt.0.0d0)then
             write(2,112)
           endif
c
           if(Tgrav1.le.0.0d0)then
             write(*,113)bell
           endif
c
           if(Tgrav2.le.0.0d0)then
             write(*,114)bell
           endif
c
           if(betarim.lt.0.0d0)then
             write(*,115)bell
             stop
           endif
c
           if(rinner.lt.0.0d0)then
             write(*,116)bell
             stop
           endif
c
           if(router.lt.0.0d0)then
             write(*,117)bell
             stop
           endif
c
           if(temprat.gt.0.0d0)then
             Teff2=Teff1*temprat
             write(2,2118)Teff2
           endif
c
           if((rinner.lt.fill2).and.(Teff2.gt.0.0d0).and.(idint.ge.1))then
             write(*,118)bell,fill2
             write(2,118)bell,fill2
             rinner=fill2
           endif
c
           if(rinner.gt.1.0d0)then
             write(*,119)bell
             stop
           endif
c
           if(router.gt.1.0d0)then
             write(*,120)bell
             stop
           endif
c
c           if((idint.ge.1).and.(tdisk.lt.0.0d0))then
c             write(*,121)bell
c             stop
c           endif
c
           if((idint.ge.1).and.(Ntheta.lt.1))then
             write(*,122)bell
             stop
           endif
c
           if((idint.ge.1).and.(Ntheta.gt.Nthetamax))then
             write(*,123)bell,Nthetamax
             stop
           endif
c
           if((idint.ge.1).and.(Nradius.lt.1))then
             write(*,124)bell
             stop
           endif
c
           if((idint.ge.1).and.(Nradius.gt.Nrmax))then
             write(*,125)bell,Nrmax
             stop
           endif
c
           if(alb1.lt.0.0d0)then
             write(*,126)bell
             stop
           endif
c
           if(alb2.lt.0.0d0)then
             write(*,127)bell
             stop
           endif
c
           if(alb1.gt.1.0d0)then
             write(*,128)bell
             stop
           endif
c
           if(alb2.gt.1.0d0)then
             write(*,129)bell
             stop
           endif
c
           if((separ.lt.0.0d0).and.(Period.le.0.0d0))then
             write(*,130)bell
             stop
           endif
c
           if((separ.lt.0.0d0).and.(fm.le.0.0d0))then
             write(*,131)bell
             stop
           endif
c
c   RVG BUG ALERT  April 19, 2001
c
c   Add the variable ecc to the argument list of getradius.
c
           if(separ.lt.0.0d0)then
             write(2,132)
             call getradius(Q,finc,rad_in_cm,fm,period,ecc)
             separ=rad_in_cm/6.9598d10
             write(2,133)separ
           endif
c
           if(idraw.gt.0)then
             write(2,134)
           endif
c
           if(iecheck.lt.0)then
             write(2,135)
           endif
c
           if(idint.le.0)then
             write(2,136)
           endif
c
c           if(Nref.lt.0)then
c             write(*,137)bell
c             Nref=0
c           endif
c
           if(Nref.gt.5)then
             write(*,138)bell
           endif
c
           if((ecc.lt.0.0d0).or.(ecc.ge.1.0d0))then
             write(*,141)bell
             ecc=0.0d0
           endif
c
           if((pshift.lt.-1.0d0))then
             write(*,142)bell
             pshift=-1.0d0
           endif
           if((pshift.gt.1.0d0))then
             write(*,142)bell
             pshift=1.0d0
           endif

           icount=0
           if(icnU.eq.0)then
             icnU=430
           else
             icount=icount+1
           endif
           if(icnB.eq.0)then
             icnB=430
           else
             icount=icount+1
           endif
           if(icnV.eq.0)then
             icnV=430
           else
             icount=icount+1
           endif
           if(icnR.eq.0)then
             icnR=430
           else
             icount=icount+1
           endif
           if(icnI.eq.0)then
             icnI=430
           else
             icount=icount+1
           endif
           if(icnJ.eq.0)then
             icnJ=430
           else
             icount=icount+1
           endif
           if(icnH.eq.0)then
             icnH=430
           else
             icount=icount+1
           endif
           if(icnK.eq.0)then
             icnK=430
           else
             icount=icount+1
           endif
c
c        
c    UPDATE June 16, 2003
c
c    If in EBOP mode, set ism1=0
c
           if((iecheck.eq.9).and.(ecc.gt.0.0d0))then
             ism1=0
             write(2,9292)
           endif

           if(icount.eq.0)write(*,140)bell
c
c    July 29, 2005
c
c    If ecosw > 0 and ecc > 0 then set the value of argper
c
          if((ecc.gt.0.0d0).and.(ecosw.gt.0.0d0))then
            call getom(ecc,ecosw,argper)
            write(2,88)argper
          endif 
c
88         format('Info:  The value of argper has been set to ',
     %           f8.4,' degrees')
 92        format(a1,'Warning:  Nalph1 is probably too small (try something',
     %                ' greater than 10)')
 93        format(a1,'Warning:  Nalph2 is probably too small (try something',
     %                ' greater than 10)')
 94        format(a1,'Warning:  Nbet1 is probably too small (try something',
     %                ' greater than 6')
 95        format(a1,'Warning:  Nbet2 is probably too small (try something',
     %                ' greater than 6')
 96        format(a1,'Error:  Nalph1 is less than 1')
 97        format(a1,'Error:  Nalph2 is less than 1')
 98        format(a1,'Error:  Nbet1 is less than 1')
 99        format(a1,'Error:  Nbet2 is less than 1')
 100       format(a1,'Error:  Nalph1 exceeds the maximum limit (currently ',
     %                  i4,')')
 101       format(a1,'Error:  Nalph2 exceeds the maximum limit (currently ',
     %                  i4,')')
 102       format(a1,'Error:  4*Nbet1 exceeds the maximum limit (currently ',
     %                  i4,')')
 103       format(a1,'Error:  4*Nbet2 exceeds the maximum limit (currently ',
     &                  i4,')')
 104       format(a1,'Error:  fill1 exceeds 1.0')
 105       format(a1,'Error:  fill2 exceeds 1.0')
 106       format(a1,'Error:  omega1 is less than 0.0')
 107       format(a1,'Error:  omega1 is less than 0.0')
 108       format(a1,'Error:  dphase is out of bounds (0.0 < dphase < 360.0)')
 109       format(a1,'Error:  Q is negative')
 110       format(a1,'Error:  finc is out of bounds (0.0 <= finc <= 90.0)')
 111       format(a1,'Error:  Teff1 is negative')
 112       format(/'Info:  Teff2 is negative---star 2 will be invisible')
 2118      format(/'Info:  Teff2 is set to ',f10.4,' K')
 113       format(a1,'Warning:  Tgrav1 is negative')
 114       format(a1,'Warning:  Tgrav2 is negative')
 115       format(a1,'Error:  betarim is negative')
 116       format(a1,'Error:  rinner is negative')
 117       format(a1,'Error:  router is negative')
 118       format(a1,'Warning:  rinner is in inside star 2---setting rinner ',
     %         'to fill2 (currently ',f5.3,')')
 119       format(a1,'Error:  rinner exceeds 1.0')
 120       format(a1,'Error:  router exceeds 1.0')
 121       format(a1,'Error:  tdisk is negative')
 122       format(a1,'Error:  Ntheta is less than 1')
 123       format(a1,'Error:  Ntheta exceeds the maximum (currently ',i3,')')
 124       format(a1,'Error:  Nradius is less than 1')
 125       format(a1,'Error:  Nradius exceeds the maximum (currently ',i3,')')
 126       format(a1,'Error:  alb1 is less than 0')
 127       format(a1,'Error:  alb2 is less than 0')
 128       format(a1,'Error:  alb1 exceeds 1.0')
 129       format(a1,'Error:  alb2 exceeds 1.0')
 130       format(a1,'Error:  Period is less than 0.0')
 131       format(a1,'Error:  fm is less than 0.0')
 132       format('Info:  separation is less than 0.0---computing ',
     %             'the separation from Period and ',
     %            '         f(M)')
 133       format('Info:  the separation has been set to ',f9.5,' solar radii')
 134       format('Info:  output files for drawing codes will be written')
 135       format('Info:  eclipse checking is turned off')
 136       format(/'Info:  there will be no accretion disk')
 137       format(a1,'Error:  Nref is less than zero---setting Nref=0')
 138       format(a1,'Warning:  Nref is rather large---long execution time')
 140       format(a1,'Error:  No filters are specified')
 141       format(a1,'Error:  eccentricity is out of range')
 142       format(a1,'Error:  pshift is out of range')
c
c   UPDATE April 8, 2002
c
1199       format('Info:  M_1 fixed at ',f6.3,' solar masses.  The  
     %        separation has been set to ',f13.7,
     $         ' solar radii')
c
 9292      format('Info:  Use ism1=0 for iececk=9 and ecc > 0.0')
           return
           end
c
c &&&&&&&&&&&&&&&&
c
           subroutine getradius(Q,finc,rad_in_cm,fm,period,ecc)
c
c   November 1, 1999
c 
c   This routine will return the radius of the orbit in cm, 
c   given the mass function fm in solar masses, the period in days, the
c   inclination finc in degrees, and mass ratio Q.  It is assumed that
c   the mass function of star 2 was measured from the motion of star 1.
c
c   RVG BUG ALERT  April 19, 2001
c
c   Add the variable ecc to the argument list of getradius.  Define fmnew
c   as down below.
c
c
          implicit double precision(a-h,o-z)
c
           eccfac=dsqrt(1.0d0-ecc*ecc)**3
           fmnew=fm/eccfac
c
c    UPDATE Oct31, 2002
c
c    Define fmnew as fm
c
           fmnew=fm
           fincr=finc*0.017453292d0        ! radians
           ppp=period*24.0d0               ! period in hours
           ovq=1.0d0/Q
           sinecubed=(dsin(fincr)**3)
           x_mass=fmnew*(1.0d0+ovq)*(1.0d0+ovq)/sinecubed   ! mass of star 2
           sec_mass=x_mass*ovq                       ! mass of star 2
           total_mass=x_mass+sec_mass   !total mass in solar masses
           coef=3.518847d10
           rad_in_cm=(coef)*(ppp*ppp*total_mass)**(0.33333333333333d0)
c 
           return
           end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine disksetup(Nthetamax,Nrmax,Ntheta,Nradius,
     %       betarim,rinner,router,reff2,Rl2,separ,
     $       tdisk,xi,dtemp,dx,dy,dz,drad,
     $       tedge,xedge,yedge,zedge,redge,stepr,stepz,bdist,
     #       ivrt,reper,rsper)
c
c   November 1, 1999
c
c   This routine will return the x, y, and z-coordinates of the grid
c   points on the disk.  Also, the temperatures of each point are returned.
c
c   November 29, 1999
c
c   Use as the radial coordinate zeta=2*dsqrt(r).  This will give more
c   points near the inner edge where the temperature is changing rapidly.
c
c
c   February 12, 2001
c
c   Put in the general case for eccentric orbits.  The flat ivrt=0 means
c   to setup the disk at periastron.  The values of rsmall and redge
c   are saved and used for other phases.
c
c
          implicit double precision(a-h,o-z)
c
          parameter (pie=3.14159265358979323d0)
          dimension dtemp(Nrmax*Nthetamax),dx(Nrmax*Nthetamax),
     &      dy(Nrmax*Nthetamax),dz(Nrmax*Nthetamax),
     &      drad(Nrmax),xedge(Nthetamax*11),yedge(Nthetamax*11),
     &      zedge(Nthetamax*11),tedge(Nthetamax*11)
c
c   Start with the disk face.  It is assumed that the lower face is exactly
c   the same as the upper face, but with a negative z-value.  In practice,
c   however, we never see the bottom face.
c          
          radcon=pie/180.0d0
          if(ivrt.eq.0)then
            redge=router*reff2              ! radius of outer edge in x units
            rsmall=rinner*Rl2               ! radius of inner edge in x units
            reper=redge
            rsper=rsmall
          else
            redge=reper
            rsmall=rsper
          endif
          betarad=betarim*radcon       ! radians
          steptheta=360.0d0/dble(ntheta)
          stepr=(redge-rsmall)/dble(Nradius-1)
c
c   Transform r into zeta.
c
          zetain=2.0d0*dsqrt(rsmall)
          zetaout=2.0d0*dsqrt(redge)
          stepzeta=(zetaout-zetain)/dble(Nradius-1)
c
c  UPDATE DECEMBEE 17, 2001
c
c  If tdisk is negative, we should ignore the disk flux.
c  If this is the case, then assign a dummy disk temperature.
c
          tdummy=10000.0
          if(tdisk.gt.0.0d0)tdummy=tdisk
c
          theta=0.0d0
          zeta=zetain
          do 10 ir=1,Nradius
            zeta=zetain+dble(ir-1)*stepzeta
            drad(ir)=zeta
            r=0.25d0*zeta*zeta
            do 9 ithet=1,Ntheta              ! theta goes from zero to 360-step
              theta=dble(ithet)*steptheta -0.5*steptheta  ! degrees
              thetar=theta*radcon            ! radians
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(ir-1)*Ntheta+ithet
c
              dx(iidx)=bdist-r*dcos(thetar)
              dy(iidx)=r*dsin(thetar)
              dz(iidx)=r*dtan(betarad)
c
c  UPDATE DECEMBEE 17, 2001
c
c  If tdisk is negative, we should ignore the disk flux.
c  If this is the case, then assign a dummy disk temperature.
c
              dtemp(iidx)=tdummy*(r/rsmall)**(xi)
c              dtemp(ir,ithet)=tdisk*(r/rsmall)**(xi)
              tlast=dtemp(iidx)
 9          continue
 10       continue
c
          zrim=redge*dtan(betarad)      ! z coordinate of outer edge
          stepz=zrim*0.2d0             ! 5 steps above and 5 steps below
c                                    ! the plane
          do 20 iz=-5,5
            z=dble(iz)*stepz
            do 19 ithet=1,Ntheta             ! theta goes from zero to 360-step
              theta=dble(ithet)*steptheta-0.5*steptheta   ! degrees
              thetar=theta*radcon         ! radians
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(iz+6-1)*Ntheta+ithet
c
              xedge(iidx)=bdist-redge*dcos(thetar)
              yedge(iidx)=redge*dsin(thetar)
              zedge(iidx)=z
              tedge(iidx)=tlast
 19         continue
 20       continue
c
          write(2,100)rsmall*separ,redge*separ
c
c   UPDATE DECEMBER 17, 2001
c
c   Add an if-then clause in case the disk temperature is negative.
c
          if(tdisk.gt.0.0d0)then
            iidx=(6-1)*Ntheta+1
            write(2,101)tlast,xedge(iidx),yedge(iidx),zedge(iidx)
          else
            write(2,102)
          endif
c
c   Scale the step size in zeta by 1/cos^2(betarim) for use in the
c   flux summing routine.
c
          stepr=stepzeta/(dcos(betarad)*dcos(betarad))

 100      format(//'inner disk radius = ',f11.6,' solar radii',
     &      5x,/'outer disk radius = ',f11.6,' solar radii')
 101      format('temperature at outer edge = ',f7.1,3x,3(f6.4,1x))
 102      format('disk temperature negative, will ignore disk flux')
c
          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine dummyvalues(ialphmax2,ibetmax2,Nalph2,Nbet2,
     $        x2,y2,z2,surf2,rad2,gradx2,grady2,gradz2,g2,xend2,darkbol2,
     #        temp2,ibetlim2,mmdx2)
c
c  November 1, 1999
c
c  This routine is called when Teff2 < 0.  In this case, star 2 is invisible
c  (usually X-ray binary mode).  The coordinates are set to 0.0, the
c  gradients are set to 0.0, the gravities are set to 1.0, and the
c  limbdarkening coefficients are set to 0.0
c
c  UPDATE March 22, 2002
c
c  Add ibetlim2 to the argument list.  Dimension it below:
c
          implicit double precision(a-h,o-z)

          dimension x2(ialphmax2*ibetmax2),y2(ialphmax2*ibetmax2),
     %      z2(ialphmax2*ibetmax2),surf2(ialphmax2*ibetmax2),
     &      rad2(ialphmax2*ibetmax2),gradx2(ialphmax2*ibetmax2),
     &      grady2(ialphmax2*ibetmax2),gradz2(ialphmax2*ibetmax2),
     $      g2(ialphmax2*ibetmax2),xend2(4),temp2(ialphmax2*ibetmax2),
     %      ibetlim2(ialphmax2),mmdx2(ialphmax2,ibetmax2)
c
          mmcount=0
          do 10 ialf=1,Nalph2
c
c   UPDATE March 22, 2002
c
c   Define the values of ibetlim2, and assign gradients values of 1
c
            ibetlim2(ialf)=4*Nbet2
            do 9 ibet=1,4*Nbet2
c              iidx=(ialf-1)*4*Nbet2+ibet
              mmcount=mmcount+1
              mmdx2(ialf,ibet)=mmcount
              iidx=mmcount
              x2(iidx)=0.0d0          
              y2(iidx)=0.0d0          
              z2(iidx)=0.0d0          
              gradx2(iidx)=1.0d0          
              grady2(iidx)=1.0d0          
              gradz2(iidx)=1.0d0
              g2(iidx)=1.0d0      
              surf2(iidx)=1.762429d-2          
              temp2(iidx)=1.0d0          
 9          continue
 10       continue
c
          xend2(1)=9.9999999d0
          xend2(2)=9.9999999d0
          reff=9.9999999d0
          psi0=-9.999999d0
          Tpole=-9.99999d0
c
          darkbol2=0.0d0
c
          write(3,3000)reff
          write(3,3001)reff
          write(3,3002)xend2(1)
          write(3,3003)xend2(2)
          write(3,3004)psi0
          write(3,3005)reff
          write(3,3006)reff
          write(3,3007)Tpole
c
 3000     format( f9.6,11x,'R_eff (star 2)')
 3001     format( f9.6,11x,'R_pole (star 2)')
 3002     format( f9.6,11x,'x(point) (star 2)')
 3003     format( f9.6,11x,'x(end) (star 2)')
 3004     format(f12.6, 8x,'potential at fill2*L1')
 3005     format(e16.9, 4x,'surface area (star 2)')
 3006     format(e16.9, 4x,'volume (star 2)')
 3007     format(f10.4,10x,'polar temperature (star 2)')

          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine dummyhoriz(ibetmax2,Nbet2,Nhoriz2,
     %       xhoriz2,yhoriz2,Ntop2,xtop2horiz,ytop2horiz)
c
c  November 1, 1999
c
c  This routine is called if Teff2 < 0, in which case star 2 is invisible
c  (X-ray binary mode usually).  The horizon of star 2 is moved off to
c  a suitably large distance so that star 1 is never eclipsed.
c
c
          implicit double precision(a-h,o-z)
c
          dimension xhoriz2(ibetmax2),yhoriz2(ibetmax2),
     %       xtop2horiz(ibetmax2),ytop2horiz(ibetmax2)
c
          Nhoriz2=Nbet2
          Ntop2=Nbet2
c
          steptheta=360.0d0/dble(Nbet2)
          r=10.0
          do 10 i=1,Nbet2
            theta=dble(i)*steptheta
            thetar=theta*3.14159265358979d0/180.0d0
            x=r*dcos(thetar)
            y=r*dsin(thetar)
            xhoriz2(i)=x+9.99d2
            yhoriz2(i)=y+9.99d2
            xtop2horiz(i)=x+9.99d2
            ytop2horiz(i)=y+9.99d2
 10       continue
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine zheight(xA,yA,zA,xB,yB,zB,rout,zcross)
c
c   This program assumes there are two stars in the normal coordinate
c   frame, and that there is a disk of radius r_outer centered on
c   star 2.  Given a point on star 1 (xA,yA,zA), and on star 2 (xB,yB,zB),
c   this routine determines the z-value of where the line joining the two
c   points crosses the disk edge.
c
c   Use Newton's method to find the t such that dziff=0
c
c
          implicit double precision(a-h,o-z)
c
          t=0.5d0
          do 10 i=1,6  
            diff=dabs(t-tnew)          
            tnew=t-zdiff(xA,yA,zA,xB,yB,zB,rout,t)
     &             /zprimediff(xA,yA,zA,xB,yB,zB,rout,t)
            t=tnew
 10       continue
c
          x=xA+(xB-xA)*t
          y=yA+(yB-yA)*t
          z=zA+(zB-zA)*t
c
          zcross=z          
c
          return
          end
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          double precision function zdiff(xA,yA,zA,xB,yB,zB,rout,ttt)
c
c   November 2, 1999
c
c   When given two grid points xA,yA,zA and xB,yB,zB, parametric equations
c   in ttt can be written that give the x,y,z coordinates of points along
c   the line joining the two grid points.
c   This function computes the difference between the radius of the
c   x,y coordinates
c   of a point given by ttt
c   on a line joining grid points xA,yA,zA and xB,yB,zB and
c   the radius rout.
c
c
          implicit double precision(a-h,o-z)
c
          t1=(xA+(xB-xA)*ttt-1.0)**2
          t2=(yA+(yB-yA)*ttt)**2
          t3=0.0d0  
          zdiff=dsqrt(t1+t2+t3)-rout
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          double precision function zprimediff(xA,yA,zA,xB,yB,zB,rout,ttt)
c
c   November 2, 1999
c
c   The derivitive wrt to ttt of zdiff given above.   
c
c
          implicit double precision(a-h,o-z)
c
          t1=(xA+(xB-xA)*ttt-1)**2
          t2=(yA+(yB-yA)*ttt)**2
          t3=0.0d0           !(zA+(zB-zA)*ttt)**2
c
          t4=2.0d0*(xA+(xB-xA)*ttt-1.0d0)*(xB-xA)
          t5=2.0d0*(yA+(yB-yA)*ttt)*(yB-yA)
          t6=0.0d0  
c
          zprimediff=0.5d0*(t4+t5+t6)/dsqrt(t1+t2+t3)
c
          return
          end
c
c   ################################################
c 
          subroutine gettophorizon(istar,
     %      ialphmax,ibetmax,Nalf,Nbet,ibetlim,phase,finc,Q,psi0,omega,
     $      xarray,yarray,zarray,radarray,gradx,grady,gradz,xend,
     $      Nhoriz,xhoriz,yhoriz,phiar,iedgestar,delphiedge,bdist,
     *      mmdx)
c
c  November 3, 1999
c
c  This routine will consider only positive z-values of the star and
c  find the apparent horizon.  This horizon for star 2 is needed if there
c  is a disk---points on the disk that project (in sky coordinates)
c  inside this horizon are behind star 2.
c
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xarray(ialphmax*ibetmax),yarray(ialphmax*ibetmax),
     $        zarray(ialphmax*ibetmax),xhoriz(4*ibetmax),yhoriz(4*ibetmax),
     $        gradx(ialphmax*ibetmax),grady(ialphmax*ibetmax),xend(4),
     $        gradz(ialphmax*ibetmax),arrproj(10000),xdummy(10000),
     &        ydummy(10000),
     $        xfirstring(10000),yfirstring(10000),radarray(ialphmax*ibetmax),
     $        xlastring(10000),ylastring(10000),ibetlim(ialphmax),
     $        phiar(ialphmax*ibetmax),iedgestar(ialphmax*ibetmax),
     $        delphiedge(ialphmax*ibetmax),mmdx(ialphmax,ibetmax)
c
c  Initialize.
c
          do 1 i=1,5000
            xdummy(i)=0.0d0
            ydummy(i)=0.0d0
 1        continue
c
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c          
          NBET4 = NBET*4
          DBETA = (pie/2.0d0)/NBET        ! step size in longitude
          AZ = DCOS(FINCR)
          IF (AZ.LT.0.0d0) AZ = 0.0d0
          AX = -DSIN(FINCR)*DCOS(PHASER)    ! l in Wilson & Sofia
          AY = DSIN(FINCR)*DSIN(PHASER)     ! m in Wilson & Sofia
          A2 = DACOS(AX)
          A3 = DSIN(A2)
          IF (A3.LT.0.0d0) A3=0.0d0
          IF (A3.EQ.0.0d0) GO TO 508
          B1=AZ/DSIN(A2)
          IF(B1.GT.1.0d0) B1=1.d0
          BETA = DASIN(B1)               !beta is the angle between the surface
          GO TO 509                      !normal and the radius vector.
508       BETA = 0.0d0
509       KBETA = BETA/DBETA + 0.5d0
          KBETA1 = KBETA + 1
          KBETAN = KBETA + 2*NBET
c
          Ndummy=0
          Nfirst=0
          Nlast=0
c
c   Loop over alpha first.  
c
          do 950 ialf=1,nalf/2       
            do 949 ibet=1,ibetlim(ialf)    !2*Nbet   
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c               iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              arrproj(ibet)=AX*GRADX(iidx) + AY*GRADY(iidx) + 
     1	            AZ*GRADZ(iidx)               
c
c   Check the last latitude row before the xy plane and include that
c   in the top horizon
c
            if(ialf.eq.nalf/2)then
              if(arrproj(ibet).ge.0.0d0)then
                izz=ialf
                jzz=1
c                iidx=(izz-1)*ibetlim(ialf)+jzz
c                iidx=kount(ialphmax,izz,ibetlim)+jzz
                iidx=mmdx(izz,jzz)
                xx=xarray(iidx)
                yy=yarray(iidx)
                zz=zarray(iidx)
                xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)  
                yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)    
                Ndummy=Ndummy+1
                xdummy(Ndummy)=xp
                ydummy(Ndummy)=yp
              endif
            endif
 949        continue
c
c   We have the array of projection factors along a given direction.
c   Now find out where the sign change is.  This is where the line of sight
c   has moved over the horizon.
c
            do 940 ibet=1,ibetlim(ialf)  !  -1    !nbet4-1
              if(ibet.lt.ibetlim(ialf))then
                index=ibet+1
              else
                index=1
              endif
              rsign=arrproj(ibet)*arrproj(index)
              if(rsign.le.0.0d0)then     ! crossed over
                if(arrproj(ibet).gt.0.0d0)then
                  izz=ialf
                  jzz=ibet
c                  iidx=(izz-1)*ibetlim(ialf)+jzz
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xvis=xarray(iidx)
                  yvis=yarray(iidx)
                  zvis=zarray(iidx)
                  rvis=radarray(iidx)
                  phivis=phiar(iidx)
c
                  izz=ialf
                  jzz=index
c                  iidx=(izz-1)*ibetlim(ialf)+jzz
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xhid=xarray(iidx)
                  yhid=yarray(iidx)
                  zhid=zarray(iidx)
                  rhid=radarray(iidx)
                  phihid=phiar(iidx)
                else
                  izz=ialf
                  jzz=index
c                  iidx=(izz-1)*ibetlim(ialf)+jzz
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xvis=xarray(iidx)
                  yvis=yarray(iidx)
                  zvis=zarray(iidx)
                  rvis=radarray(iidx)
                  phivis=phiar(iidx)
c
                  izz=ialf
                  jzz=ibet
c                  iidx=(izz-1)*ibetlim(ialf)+jzz
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  xhid=xarray(iidx)
                  yhid=yarray(iidx)
                  zhid=zarray(iidx)
                  rhid=radarray(iidx)
                  phihid=phiar(iidx)
                endif
c
                xx=xvis
                yy=yvis
                zz=zvis
c
c   RVG BUG ALERT  April 20, 2001
c
c   Leave out the call to acchor here, since it is not really needed.
c
c                call acchor(overQ,psi0,omega,
c     $              xvis,yvis,zvis,rvis,phivis,xhid,yhid,zhid,rhid,
c     $              phihid,ax,ay,az,xacc,yacc,zacc,bdist)
cc
c                xx=xacc  
c                yy=yacc  
c                zz=zacc  
c
                xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist)  
                yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)    
                Ndummy=Ndummy+1
                xdummy(Ndummy)=xp
                ydummy(Ndummy)=yp
              endif
 940        continue
 950      continue
c
c   RVG BUG ALERT   April 20, 2001
c
c   Comment out everything until statement 69 below.
c
c   Check the visibility of the nose and see if it should be in
c   the horzon
c
c          if((phase.eq.90.0d0).or.(phase.eq.270.0d0))then ! both ends visible
c            xp=xtran(xend(1),0.0d0,0.0d0,phase,fincr,Q,istar,bdist) 
c            yp=ytran(xend(1),0.0d0,0.0d0,phase,fincr,Q,istar,bdist) 
c            Ndummy=Ndummy+1      ! outside the first ring.  Include the
c            xdummy(Ndummy)=xp    ! point in the horizon
c            ydummy(Ndummy)=yp
c            xp=xtran(xend(2),0.0d0,0.0d0,phase,fincr,Q,istar,bdist) 
c            yp=ytran(xend(2),0.0d0,0.0d0,phase,fincr,Q,istar,bdist)    
c            Ndummy=Ndummy+1      ! outside the first ring.  Include the
c            xdummy(Ndummy)=xp    ! point in the horizon
c            ydummy(Ndummy)=yp
c            go to 69
c          endif
c
c   Now 'sort' the horizon points so that a regular polygon is made.
c
 69       call sortcircle(Ndummy,xdummy,ydummy)
c
c   The above array may have repeated points.  Remove them.
c
          call uniquepoint(Ndummy,xdummy,ydummy,Nhoriz,xhoriz,yhoriz)
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine getdiskhoriz(Nthetamax,Nrmax,Ntheta,Nradius,
     %       Q,phase,finc,xedge,yedge,zedge,Ndhoriz,dxhoriz,dyhoriz,
     %       Ndtop,dtopx,dtopy,bdist)
c
c   November 3, 1999
c
c   This routine will find the horizon of the disk.  The horizon is simply
c   the lower rim for the front part (proj>0) and the upper rim for the
c   back part (proj<0).
c
c
          implicit double precision(a-h,o-z)
c
          parameter (pie=3.14159265358979323d0)
          dimension xedge(Nthetamax*11),yedge(Nthetamax*11),
     &      zedge(Nthetamax*11),dxhoriz(2*Nthetamax),dyhoriz(2*Nthetamax),
     %      dtopx(2*Nthetamax),dtopy(2*Nthetamax)       
          dimension index(2)
c
          radcon=pie/180.0d0
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          fincr=finc*radcon
          sifinc=dsin(fincr)
c
          Ndhoriz=0
          Ndtop=0
          kkk=0
          steptheta=360.0d0/dble(ntheta)
          do 19 ithet=1,Ntheta              ! theta goes from zero to 360-step
            theta=dble(ithet)*steptheta-0.5*steptheta   ! degrees
            thetar=theta*radcon             ! radians
            angdiff=(theta-phase)*radcon
            proj=sifinc*dcos(angdiff)  
            if(proj.gt.0.0)then
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(1-1)*Ntheta+ithet
c
              xx=xedge(iidx)
              yy=yedge(iidx)
              zz=zedge(iidx)
              xp=diskxtran(xx,yy,zz,phase,fincr,Q,1,bdist)  
              yp=diskytran(xx,yy,zz,phase,fincr,Q,1,bdist)  
              Ndhoriz=Ndhoriz+1
              dxhoriz(Ndhoriz)=xp
              dyhoriz(Ndhoriz)=yp
            else
              iidx=(11-1)*Ntheta+ithet
              xx=xedge(iidx)
              yy=yedge(iidx)
              zz=zedge(iidx)
              xp=diskxtran(xx,yy,zz,phase,fincr,Q,1,bdist)  
              yp=diskytran(xx,yy,zz,phase,fincr,Q,1,bdist)  
              Ndhoriz=Ndhoriz+1
              dxhoriz(Ndhoriz)=xp
              dyhoriz(Ndhoriz)=yp
            endif
            iidx=(11-1)*Ntheta+ithet
            xx=xedge(iidx)
            yy=yedge(iidx)
            zz=zedge(iidx)
            xp=diskxtran(xx,yy,zz,phase,fincr,Q,1,bdist)  
            yp=diskytran(xx,yy,zz,phase,fincr,Q,1,bdist)  
            Ndtop=Ndtop+1
            dtopx(Ndtop)=xp
            dtopy(Ndtop)=yp
c
            if(ithet.gt.1)then
              rsign=proj*rlastproj
              if(rsign.lt.0.0d0)then
                kkk=kkk+1
                index(kkk)=ithet-1
              endif
            endif
              rlastproj=proj
 19       continue
c
c          do 50 iz=1,11
c            do 49 ii=1,2
c              xx=xedge(index(ii),iz)
c              yy=yedge(index(ii),iz)
c              zz=zedge(index(ii),iz)
c              xp=xtran(xx,yy,zz,phase,fincr,Q,1,bdist)  
c              yp=ytran(xx,yy,zz,phase,fincr,Q,1,bdist)  
c              Ndhoriz=Ndhoriz+1
c              dxhoriz(Ndhoriz)=xp
c              dyhoriz(Ndhoriz)=yp
c 49         continue
c 50       continue
c           
          call sortcircle(Ndhoriz,dxhoriz,dyhoriz)
          call sortcircle(Ndtop,dtopx,dtopy)
c
          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine diskvisib(Nrmax,Nthetamax,Nradius,Ntheta,phase,finc,Q,
     %      betarim,dx,dy,dz,xedge,yedge,zedge,
     &      diskproj,edgeproj,dvisib,evisib,Nskydisk,xskydisk,yskydisk,
     #      zskydisk,
     &      Nskyedge,xskyedge,yskyedge,Ntop2,
     %      xtop2horiz,ytop2horiz,Nhoriz1,xhoriz1,yhoriz1,
     &      Ndtop,dtopx,dtopy,iecheck,Neclipse,bdist)
c
c   November 8, 1999
c
c   This routine will return an array of projection angles for the
c   disk (diskproj) and for the edge (edgeproj).  It will also check
c   for eclipses by star 1, the top of star 2, and the rim of the disk.
c   If any part of the disk is blocked, then the visibility of that
c   element is set to zero (dvisib, evisib).  For elements that are not
c   blocked, the visibilities are set to the projection.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension diskproj(Nrmax*Nthetamax),edgeproj(Nthetamax*11),
     &      dvisib(Nrmax*Nthetamax),evisib(Nthetamax*11),
     &      xskydisk(Nthetamax*Nrmax),yskydisk(Nthetamax*Nrmax),
     &      zskydisk(Nthetamax*Nrmax),
     &      xskyedge(Nthetamax*11),yskyedge(Nthetamax*11),
     $      xtop2horiz(Ntop2),ytop2horiz(Ntop2),
     @      xhoriz1(Nhoriz1),yhoriz1(Nhoriz1),dx(Nrmax*Nthetamax),
     &      dy(Nrmax*Nthetamax),dz(Nrmax*Nthetamax),xedge(Nthetamax*11),
     &      yedge(Nthetamax*11),zedge(Nthetamax*11),
     &      dtopx(2*Nthetamax),dtopy(2*Nthetamax)
c
c  initialize the visibilities
c
          radcon=pie/180.0d0
          do 2 i=1,Nrmax
            do 1 j=1,Nthetamax
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(i-1)*Ntheta+j
c
              dvisib(iidx)=0.0d0
              diskproj(iidx)=0.0d0
 1          continue
 2        continue
c
          do 4 i=1,Nthetamax
            do 3 j=1,11
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(j-1)*Ntheta+i
c
              evisib(iidx)=0.0d0
              edgeproj(iidx)=0.0d0
 3          continue
 4        continue
c
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c          
          Nskydisk=0
          Nskyedge=0
          Neclipse=0
c
          infront=0
          if((phase.ge.90.0d0).and.(phase.le.270.0d0))infront=1
c
c   Start with the disk face.  Loop over theta and r, compute the projection
c   factors from the standard coordinates, check for eclipses, etc.
c
          steptheta=360.0d0/dble(ntheta)
          betarad=betarim*radcon                  ! radians
c
c   Put the trig calls outside of the loop.
c
          cophase=dcos(phaser)
          siphase=dsin(phaser)
          cofinc=dcos(fincr)
          sifinc=dsin(fincr)
          cobet=dcos(betarad)
          sibet=dsin(betarad)
          do 20 ir=1,Nradius
            do 19 ithet=1,Ntheta             ! theta goes from zero to 360-step
              theta=dble(ithet)*steptheta-0.5*steptheta  ! degrees
              thetar=theta*radcon            ! radians
c
               iidx=(ir-1)*Ntheta+ithet
c               
c              t1=-dcos(phaser)*dsin(fincr)*dsin(betarad)*dcos(thetar)
c              t2=-dsin(phaser)*dsin(fincr)*dsin(betarad)*dsin(thetar)
c              t3=dcos(fincr)*dcos(betarad)
c
              t1=-cophase*sifinc*sibet*dcos(thetar)
              t2=-siphase*sifinc*sibet*dsin(thetar)
              t3=cofinc*cobet
              proj=t1+t2+t3
c
c   This projection factor accounts for points below the rim (in cases of
c   large beta_rim and high inclination).
c
c   RVG BUG ALERT
c
c   old statement:
c
c            front=dcos(phaser-thetar)
c 
c   The correct new statement:
c
              angdiff=(theta-phase)*radcon
              front=sifinc*dcos(angdiff)  !sifinc*dcos(phaser-thetar)
c
c   END BUG
c
              if((ir.eq.Nradius).and.(front.gt.0.0).and.(proj.le.0.0))then
                proj=sifinc*dcos(angdiff)   !dcos(phaser-thetar)
              endif
              diskproj(iidx)=proj
              if(proj.le.0.0)go to 19
              xx=dx(iidx)
              yy=dy(iidx)
              zz=dz(iidx)

c              xp=xtran(xx,yy,zz,phase,fincr,Q,1,bdist)
c              yp=ytran(xx,yy,zz,phase,fincr,Q,1,bdist)
c
              xp=diskxtran(xx,yy,zz,phase,fincr,Q,1,bdist)
              yp=diskytran(xx,yy,zz,phase,fincr,Q,1,bdist)
c
              if(infront.ne.1)then
                if(iecheck.ge.0)then
                  iyes=-100
                  call insidecircle(Nhoriz1,xhoriz1,yhoriz1,xp,yp,iyes,icut)
                  if(iyes.eq.100)then
                    Neclipse=Neclipse+1
                    go to 19   !eclipsed by star 1
                  endif
                endif
              endif
c
c   Check to see of the top of the star 2 blocks the point in question.
c   Front < 0 for points that are further away than star 2.  These points
c   could be blocked by star 2.  Points with front > 0 in the inner ring
c   may be inside the top horizon of star 2 due to roundoff error.  The
c   if-then block should make sure that points in the inner ring that
c   are in front of star 2 are included.
c
              if(front.lt.0.0)then
                iyes=-100
                call insidecircle(Ntop2,xtop2horiz,ytop2horiz,
     %                  xp,yp,iyes,icut)
                if(iyes.eq.100)go to 19
              endif
c
c   Finally, check to see of the point in question is beneath the disk
c   rim.  If so, then it is invisible to the observer.  To check this
c   we simply see of the point in question is *inside* the top horizon of
c   the disk.
c
              if(ir.ne.Nradius)then
                iyes=-100
                call insidecircle(Ndtop,dtopx,dtopy,xp,yp,iyes,icut)
                if(iyes.ne.100)go to 19 
              endif
              Nskydisk=Nskydisk+1
              xskydisk(Nskydisk)=xp
              yskydisk(Nskydisk)=yp
              zskydisk(Nskydisk)=proj
              dvisib(iidx)=proj
 19         continue
 20       continue
c
          do 30 ithet=1,Ntheta             ! theta goes from zero to 360-step
            theta=dble(ithet)*steptheta-0.5*steptheta  ! degrees
            thetar=theta*radcon            ! radians
c
            do 29 iz=1,11
c
c   RVG BUG ALERT
c
c   Old statement:
c
c              proj=dsin(fincr)*dcos(phaser-thetar)            
c
c   The correct new statements:
c
              angdiff=(theta-phase)*radcon
              proj=sifinc*dcos(angdiff)  !sifinc*dcos(phaser-thetar)
c
c   END BUG
c
              iidx=(iz-1)*Ntheta+ithet

              edgeproj(iidx)=proj
              if(proj.le.0.0d0)go to 29
              xx=xedge(iidx)
              yy=yedge(iidx)
              zz=zedge(iidx)
c
c              xp=xtran(xx,yy,zz,phase,fincr,Q,1,bdist)
c              yp=ytran(xx,yy,zz,phase,fincr,Q,1,bdist)
c
              xp=diskxtran(xx,yy,zz,phase,fincr,Q,1,bdist)
              yp=diskytran(xx,yy,zz,phase,fincr,Q,1,bdist)
c
              if(infront.ne.1)then
                if(iecheck.ge.0)then
                  iyes=-100
                  call insidecircle(Nhoriz1,xhoriz1,yhoriz1,xp,yp,iyes,icut)
                  if(iyes.eq.100)then
                    Neclipse=Neclipse+1
                    go to 29   !eclipsed by star 1
                  endif
                endif
              endif
              Nskyedge=Nskyedge+1
              xskyedge(Nskyedge)=xp
              yskyedge(Nskyedge)=yp
              evisib(iidx)=proj
 29         continue
 30       continue
c
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine getdiskBBflux(Nrmax,Nthetamax,Nradius,
     %        Ntheta,diskproj,edgeproj,dvisib,evisib,dtemp,tedge,drad,
     $        dinty,einty,stepr,stepz,www,iwww,ilaw,dflux,separ)
c
c   November 8, 1999
c
c   This routine will integrate the flux from the visible parts of
c   the disk and return the result as dflux.
c
c   November 16, 1999
c
c   This routine uses passband specific limb darkening coefficients taken
c   from Van Hamme (1993, AJ, 106, 2096).  I assume log(g)=5.
c
c   November 29, 1999
c
c   Use as the radial coordinate zeta=2.0*dsqrt(r).  This will
c   put more grid points at smaller radii where the temperature
c   might be changing rapidly.  r*d(r) = (1/8)*zeta^3*d(zeta).
c   The variable 'stepr' is actually the stepsize in zeta, and the
c   array 'drad' actually contains zeta values.
c
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.  Then scale the fluxes
c   by (separ*solarrad)**2.
c
          implicit double precision(a-h,o-z)
c
          dimension diskproj(Nrmax*Nthetamax),edgeproj(Nthetamax*11),
     &      dvisib(Nrmax*Nthetamax),evisib(Nthetamax*11),
     &      dinty(Nrmax*Nthetamax),einty(Nthetamax*11),
     &      dtemp(Nrmax*Nthetamax),drad(Nrmax),
     &      tedge(Nthetamax*11)
c
          parameter(pie=3.14159265358979323d0)
c
          dflux=0.0
          C2 = 1.4384d8             ! hc/(k*1e-8)
          C1 = 1.191044d35          ! 2hc^2/((1e-8)**5)
c
          steptheta=360.0d0/dble(ntheta)
          deg2rad=pie/180.0d0
c
          c1=3.74185
          c2=14.3883
          wavemu=www/10000.0d0
          DO 10 ir=1,Nradius
            iidx=(ir-1)*Ntheta+1
            if(iwww.eq.1)call flcU50(dtemp(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.2)call flcB50(dtemp(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.3)call flcV50(dtemp(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.4)call flcR50(dtemp(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.5)call flcI50(dtemp(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.6)call flcJ50(dtemp(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.7)call flcH50(dtemp(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.8)call flcK50(dtemp(iidx),ilaw,flimbx,flimby)
c
            ddint=pie*(1.0d0-flimbx/3.0d0)
            if(ilaw.eq.2)then
              ddint=pie*(1.0d0-flimbx/3.0d0+2.0d0*flimby/9.0d0)
            endif
            if(ilaw.eq.3)then
              ddint=pie*(1.0d0-flimbx/3.0d0-flimby/5.0d0)
            endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
            if(ilaw.eq.4)then
              ddint=pie*(1.0d0-flimbx/3.0d0-flimby/6.0d0)
            endif
c
            ringflux=0.0d0
            DO 9 ithet=1,Ntheta
              iidx=(ir-1)*Ntheta+ithet
              dinty(iidx)=0.0d0
              if(diskproj(iidx).le.0.0)go to 9
c              C3 = C2/(www*dtemp(iidx))
              c3=c2/(wavemu*dtemp(iidx)/1000.0d0)
              flum=C1/(dexp(c3)-1.0d0)/wavemu**5
              if(ilaw.le.1)dark=(1.0d0-flimbx+flimbx*diskproj(iidx))
              if(ilaw.eq.2)then
                dark=1.0d0-flimbx*(1.0d0-diskproj(iidx))
                dark=dark-flimby*diskproj(iidx)*dlog(diskproj(iidx))
              endif
              if(ilaw.eq.3)then
                dark=1.0d0-flimbx*(1.0d0-diskproj(iidx))
                dark=dark-flimby*(1.0d0-dsqrt(diskproj(iidx)))
              endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
              if(ilaw.eq.4)then
                dark=1.0d0-flimbx*(1.0d0-diskproj(iidx))
                dark=dark-flimby*(1.0d0-(diskproj(iidx)))**2
              endif
c
              flum=flum*dark
              dinty(iidx)=flum ! save intensities for plotting
c
c   Here is the old term when the steps were linear in r:
c
c              flum=drad(ir)*flum*dvisib(ir,ithet)*stepr*steptheta
c
c   Here is the expression for steps in the zeta coordinate:
c
c   April 19, 2001:   RVG BUG ALERT
c
c   Convert the steptheta into radians!
c
              flum=((1.0d0/8.0d0)*drad(ir)**3)*flum*
     &                 dvisib(iidx)*stepr*steptheta*deg2rad
              ringflux=ringflux+flum
 9          continue
            dflux=ringflux+dflux
 10       continue
c
          do 20 ithet=1,Ntheta
            iidx=(1-1)*Ntheta+ithet
            if(iwww.eq.1)call flcU50(tedge(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.2)call flcB50(tedge(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.3)call flcV50(tedge(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.4)call flcR50(tedge(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.5)call flcI50(tedge(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.6)call flcJ50(tedge(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.7)call flcH50(tedge(iidx),ilaw,flimbx,flimby)
            if(iwww.eq.8)call flcK50(tedge(iidx),ilaw,flimbx,flimby)
c
            ringflux=0.0d0
            ddint=pie*(1.0d0-flimbx/3.0d0)
            if(ilaw.eq.2)then
              ddint=pie*(1.0d0-flimbx/3.0d0+2.0d0*flimby/9.0d0)
            endif
            if(ilaw.eq.3)then
              ddint=pie*(1.0d0-flimbx/3.0d0-flimby/5.0d0)
            endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
            if(ilaw.eq.4)then
              ddint=pie*(1.0d0-flimbx/3.0d0-flimby/6.0d0)
            endif
c
            do 19 iz=1,11
              iidx=(iz-1)*Ntheta+ithet
c              jjdx=(ir-1)*Ntheta+ithet
              einty(iidx)=0.0d0
              if(edgeproj(iidx).lt.0.0)go to 19
c              C3 = C2/(www*tedge(ithet,iz))
              c3=c2/(wavemu*tedge(iidx)/1000.0)
              flum=C1/(dexp(c3)-1.0d0)/wavemu**5
              if(ilaw.le.1)dark=(1.0d0-flimbx+flimbx*edgeproj(iidx))
              if(ilaw.eq.2)then
                dark=1.0d0-flimbx*(1.0d0-edgeproj(jjdx))
                dark=dark-flimby*edgeproj(iidx)*dlog(edgeproj(iidx))
              endif
              if(ilaw.eq.3)then
                dark=1.0d0-flimbx*(1.0d0-edgeproj(iidx))
                dark=dark-flimby*(1.0d0-dsqrt(edgeproj(iidx)))
              endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
              if(ilaw.eq.4)then
                dark=1.0d0-flimbx*(1.0d0-edgeproj(iidx))
                dark=dark-flimby*(1.0d0-(edgeproj(iidx)))**2
              endif
c
              flum=flum*dark
              einty(iidx)=flum ! save intensities for plotting
c
c   April 19, 2001:   RVG BUG ALERT
c
c   Convert the steptheta into radians!
c
              flum=flum*evisib(iidx)*stepz*steptheta*deg2rad
c
c   The upper rim points are in the face integration.  If iz=11, then don't
c   add to dflux.
c
              if(iz.lt.11)ringflux=ringflux+flum          
 19         continue
c
c   April 17, 2001
c
c   Add this scaling to make the flux consistent with the stellar flux
c   integration.
c            
c            dflux=pie*ringflux/ddint+dflux
c
            dflux=ringflux+dflux
c
 20       continue
c
c   UPDATE April 3, 2002
c
c   Scale the flux
c
          solarrad=6.9598d10
          dflux=dflux*(separ*solarrad)**2
c
          return
          end
c
c  %%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine wlinmod(icount,xmod,ymod,fileout,isw7)
c
c  November 10, 1999
c
c  This routine will open the output file and record the given model
c  with flux in linear units.  The phases will be scaled to go from 0.0 to
c  1.0, and an extra phase will be added.
c
c
c
          implicit double precision(a-h,o-z)

          dimension xmod(icount),ymod(icount)
          character*(*) fileout
c
          open(unit=20,file=fileout,status='unknown')
c
          do 10 i=1,icount
            write(20,100)xmod(i),ymod(i)
 10       continue
c
          if(isw7.ge.2)return
          do 20 i=1,icount
            write(20,100)xmod(i)+1.0d0,ymod(i)
 20       continue
c
          close(20)
c
 100      format(f23.15,3x,1pe21.14)
c
          return
          end
c
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine getvel(ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     &      istar,omega,phase,finc,
     %      Q,flum,xcoords,ycoords,flux,separ,period,gamma,vel,delvel,
     $      rldint,ecc,argrad,mmdx,isw13,ialfmin,ialfmax,fluxlat,
     $      bigI,bigbeta,zcoords)

c
c    November 12, 1999
c
c    This routine will compute the flux-weighted radial velocity of the
c    star in question.
c
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)    
          dimension flum(ialphmax*ibetmax),xcoords(ialphmax*ibetmax),
     &      ycoords(ialphmax*ibetmax),ibetlim(ialphmax),
     &      mmdx(ialphmax,ibetmax),zcoords(ialphmax*ibetmax)
c
c   RVG BUG ALERT   April 23, 2001
c
c   Move the definitions of argfac, efact, and dint to the top here
c
          argfac=ecc*dcos(argrad)
          efact=1.0d0/dsqrt(1.0d0-ecc*ecc)
          dint=rldint/pie
          vel=0.0d0
          delvel=0.0d0
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
c
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c
          siphase=dsin(phaser)
          cophase=dcos(phaser)
          sifinc=dsin(fincr)
          cofinc=dcos(fincr)
c
c   Compute the expected K velocity, which is the circular velocity times 
c   dsin(finc).
c
c    
          a=separ*6.9598d5            !separation in km
          p=period*24.00d0*3600.00d0       !period in seconds
          velamp=2.0d0*pie*a/p*efact
c
c    Check to see if the star has any flux before going further!
c        
          if(flux.le.0.0d0)then
c
c   RVG BUG ALERT   April 23, 2001
c
c   Move this if-then block a bit further down.
c   Set delvel=0.0d0 in the case when the flux is 0.  Also, set
c   the velocity to the velocity expected from orbital motion so the
c   output velocity curve is smooth during the eclipse
c
            ppp = (PHASE/180.0d0)*pie    
            siphase=dsin(ppp)
            vel=overQ/(1.0d0+overQ)*(siphase+argfac)*sifinc
            vel=vel*velamp+gamma
            delvel=0.0d0
            return
          endif
c
c   UPDATE October 21, 2002
c
c   Scale the fluxes by (separ*solarrad)**2
c
          sI=dsin(pie*bigI/180.0d0)
          bigbetar=pie*bigbeta/180.0d0
          cB=dcos(bigbetar)
          sB=dsin(bigbetar)

          solarrad=6.9598d10
          sscale=(separ*solarrad)**2
          do 10 ialf=1,Nalf
            if((isw13.gt.0).and.(istar.eq.0))then
              if((ialf.lt.ialfmin).or.(ialf.gt.ialfmax)) go to 10
            endif
            do 9 ibet=1,ibetlim(ialf)        !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              xxx=xcoords(iidx)
              yyy=ycoords(iidx)
              zzz=zcoords(iidx)
              xp=-xxx*siphase-yyy*cophase
              yp=xxx*cofinc*cophase-yyy*cofinc*siphase+zzz*sifinc
              v1=omega*sI*(xp*cB+yp*sB)
              v1=v1*flum(iidx)*sscale
              delvel=delvel+v1
 9          continue
 10       continue
c
          PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c
c   RVG BUG FIX   May 2, 2001
c
c   The following if-then statement is no longer needed.
c
c          if(phase.gt.180.0d0)delvel=-1.0d0*delvel

          siphase=dsin(phaser)
c
c   Bug fix, February 28, 2000
c
c          vel=velamp*(siphase*sifinc+delvel/flux)+gamma
c
c   RVG BUG ALERT   APRIL 23, 2001
c
c   Move this correction of delvel to the top, and change the last
c   part of the first vel equation.
c
          delvel=delvel/flux/dint
          vel=overQ/(1.0d0+overQ)*(siphase+argfac)*sifinc+delvel
          vel=velamp*vel+gamma
c
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine writepoints(ialphmax1,ibetmax1,Nthetamax,Nrmax,
     %           Nsky1,xsky1,ysky1,Nsky2,xsky2,
     &           ysky2,Nhoriz1,xhoriz1,yhoriz1,Nhoriz2,xhoriz2,yhoriz2,
     &           Nskydisk,xskydisk,yskydisk,zskydisk,
     &           Ntop2,xtop2horiz,ytop2horiz,
     &           Ndtop,dtopx,dtopy,Ndhoriz,dxhoriz,dyhoriz,
     %           Nskyedge,xskyedge,yskyedge,extension,separ,teff2,
     &           idint,ialphmax2,ibetmax2)
c
c   November 18, 1999
c
c   This routine will simply write the x,y sky coordinates of all of the
c   components for use in quick plotting routines.  The horizons will contain
c   potentially eclipsed points, but all other coordinates should contain
c   only visible points.
c
c
          implicit double precision(a-h,o-z)
c
          dimension xsky1(ialphmax1*ibetmax1*4),ysky1(ialphmax1*ibetmax1*4),
     &      xsky2(ialphmax2*ibetmax2*4),ysky2(ialphmax2*ibetmax2*4),
     %      xhoriz1(Nhoriz1),yhoriz1(Nhoriz1),xhoriz2(Nhoriz2),
     &      yhoriz2(Nhoriz2),xskydisk(Nskydisk),yskydisk(Nskydisk),
     &      zskydisk(Nskydisk),
     &      xtop2horiz(Ntop2),ytop2horiz(Ntop2),
     &      dtopx(2*Nthetamax),dtopy(2*Nthetamax),
     &      dxhoriz(Ndhoriz),dyhoriz(Ndhoriz),xskyedge(Nskyedge),
     &      yskyedge(Nskyedge)
c
          character*9 extension
c
          if(Nsky1.gt.0)then
            open(unit=40,file='star1coo.'//extension,status='unknown')
c
            do 10 i=1,Nsky1
              write(40,100)separ*xsky1(i),separ*ysky1(i)
 10         continue
            close(40)
c
            open(unit=40,file='star1horiz.'//extension,status='unknown')
c
            do 15 i=1,Nhoriz1
              write(40,100)separ*xhoriz1(i),separ*yhoriz1(i)
 15         continue
            close(40)
          endif             ! endif Nsky1
c
          if(Nsky2.gt.0)then
            if(teff2.gt.0.0d0)then
              open(unit=40,file='star2coo.'//extension,status='unknown')
c
              do 20 i=1,Nsky2
                write(40,100)separ*xsky2(i),separ*ysky2(i)
 20           continue
              close(40)
c
              open(unit=40,file='star2horiz.'//extension,status='unknown')
c
              do 30 i=1,Nhoriz2
                write(40,100)separ*xhoriz2(i),separ*yhoriz2(i)
 30           continue
              close(40)
c
              if(idint.ge.1)then
                open(unit=40,file='star2tophor.'//extension,status='unknown')
c
                do 40 i=1,Ntop2
                  write(40,100)separ*xtop2horiz(i),separ*ytop2horiz(i)
 40             continue
                close(40)
              endif         ! endif idint > 0
            endif         ! endif teff2 > 0
          endif         ! endif Nsky2 > 0
c
c          write(*,*)'idint writepoints = ',idint,Nskydisk

          if(idint.gt.0)then
            if(Nskydisk.gt.0)then
              open(unit=40,file='diskcoo.'//extension,status='unknown')
c
              do 50 i=1,Nskydisk
                write(40,101)separ*xskydisk(i),separ*yskydisk(i),zskydisk(i)
 50           continue
              close(40)
            endif
c
            if(Nskyedge.gt.0)then
              open(unit=40,file='diskedge.'//extension,status='unknown')
c
              do 60 i=1,Nskyedge
                write(40,100)separ*xskyedge(i),separ*yskyedge(i)
 60           continue
              close(40)
            endif
c
            if(Ndtop.gt.0)then
              open(unit=40,file='disktophor.'//extension,status='unknown')
c
              do 70 i=1,Ndtop
                write(40,100)separ*dtopx(i),separ*dtopy(i)
 70           continue
              close(40)
            endif
c
            if(Ndhoriz.gt.0)then
              open(unit=40,file='diskhoriz.'//extension,status='unknown')
c
              do 80 i=1,Ndhoriz
                write(40,100)separ*dxhoriz(i),separ*dyhoriz(i)
 80           continue
              close(40)
            endif
c
          endif
c
 100      format(2(f15.9,5x))
 101      format(3(f15.9,5x))
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&
c
c   November 19, 1999
c
c   This routine was taken from Numerical Recipes, second edition.
c
c
      SUBROUTINE LOCATE(XX,N,X,J)
c
c   Taken from Numerical Recipes.
c
          implicit double precision(a-h,o-z)

      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END
c
c  &&&&&&&&&&&&&&
c
          subroutine computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &       atmT,atmg,atmmu,Nmu,
     &       atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &       Tmax,Tmin,gmax,gmin,outinty,
     $       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     $       dwavex,dwavey,ilaw,iatm,istar)
c
c   This routine will return values of the specific intensity for the
c   8 filters based on the input values of Tin, gin, and rmuin.
c
c   UPDATE SEPTEMBER 11, 2009
c
c   Modify this routine to use a parameterized limb darkening law when
c   iatm=2.  When iatm=2, compute the intensity at mu=1, then use
c   the limb darkening law to find I(mu).
c
          implicit double precision(a-h,o-z)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines),outinty(8),tempgh(400),
     &       tempgl(400),ghinty(400,8),yscratch(400),tscratch(2),glinty(400,8),
     &       tinty(2,8),y2scratch(2),ingh(400),ingl(400),gnew(2)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
          dimension dwavex(8,2),dwavey(8,2)
c
c          call locate(atmT,Nlines,Tin,indexT)
c
c    UPDATE October 23, 2009
c
c     If iatm=2, then we are using a combination of model atmospheres
c     and a limb darkening law.  The model atmosphere intensity for 
c     mu=1 will be found, and then a limb darkening law will be
c     used to find the intensity at other values of mu.
c
c     First, we need to save the input value of mu
c
          if(iatm.ge.2)rsavemu=rmuin 
c
c
          call hunt(atmT,Nlines,Tin,itguess)
c
          indexT=itguess
          tscratch(1)=atmT(indexT+1)
c
c   Now search the Tvalues equal to atmT(indexT+1) and find the array
c   of gravities.
c
          Nhigh=0
          ingh(1)=1
          ingl(1)=1
          do 1 i=0,Nlines-(indexT+1)
            if(atmT(indexT+1+i).eq.atmT(indexT+1))then
              Nhigh=Nhigh+1
              tempgh(Nhigh)=atmg(indexT+1+i)
              ingh(Nhigh)=indexT+1+i
            else
              go to 5
            endif
 1       continue
c
5       call locate(tempgh,Nhigh,gin,indexgh)
c
c   Here is the case when the input log(g) value is in between
c   two tabulated log(g) values:
c
         if((indexgh.lt.Nhigh).and.(indexgh.gt.0))then
           gnew(1)=tempgh(indexgh)
           gnew(2)=tempgh(indexgh+1)
         endif 
c
c   If the input log(g) is larger than the largest tabulated log(g), then
c   wing it:  the polint routine will attempt to extrapolate.  This
c   is not too bad if the input log(g) is not too much larger than the
c   largest entry.
c
         if((indexgh.eq.Nhigh).and.(Nhigh.gt.1))then
           gnew(1)=tempgh(Nhigh-1)
           gnew(2)=tempgh(Nhigh)
         endif 
c
c   Same, but for an input log(g) value smaller than the smallest table
c   entry for that temperature:
c
         if(indexgh.eq.0)then
           gnew(1)=tempgh(1)
           gnew(2)=tempgh(2)
           indexgh=1               !UPDATE June 14, 2003
         endif 
c
c   Get the intensities for the gravities corresponding to indexgh
c   and indexgh+1.  The index numbers for these entries in the table
c   are ingh(indexgh) and ingh(indexgh+1).
c
c         write(*,*)Tin,gin,ingh(indexgh),rmuin,imuguess

         if(iatm.ge.2)rmuin=1.0d0
         call indexinty(ingh(indexgh),maxlines,maxmu,atmmu,
     &     atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &     Nmu,rmuin,
     %     outinty,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,imuguess)
         do 9 j=1,8
           ghinty(1,j)=outinty(j)
 9       continue
c
c          write(*,69)Tin,gin,rmuin,outinty(2)
c
         if(iatm.ge.2)rmuin=1.0d0
         call indexinty(ingh(indexgh+1),maxlines,maxmu,atmmu,
     &     atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &     Nmu,rmuin,
     %     outinty,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,imuguess)
         do 10 j=1,8
           ghinty(2,j)=outinty(j)
 10      continue
c
c          write(*,69)Tin,gin,rmuin,outinty(2)
c
          m=2
          do 20 i=1,8
            if ((i.eq.1).and.(icnU.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if ((i.eq.2).and.(icnB.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if ((i.eq.3).and.(icnV.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if ((i.eq.4).and.(icnR.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if ((i.eq.5).and.(icnI.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if ((i.eq.6).and.(icnJ.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if ((i.eq.7).and.(icnH.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if ((i.eq.8).and.(icnK.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            do 19 j=1,2
              yscratch(j)=ghinty(j,i)
 19         continue
            diff=dabs(gnew(1)-gnew(2))
            if(diff.gt.1.0d-5)then
c
c              call polint(gnew,yscratch,m,gin,qqq,dy)
c
              qqq=((gnew(2)-gin)*yscratch(1)+(gin-gnew(1))*yscratch(2))
     #             /(gnew(2)-gnew(1))

              outinty(i)=qqq
              tinty(1,i)=qqq
            else
              outinty(i)=yscratch(1)
              tinty(1,i)=yscratch(1)
            endif
 20       continue
c
c   Now search the Tvalues equal to atmT(indexT+1) and find the gravities.
c
          tscratch(2)=atmT(indexT)
c
          Nlow=0
          do 1000 i=0,indexT-1
            if(atmT(indexT-i).eq.atmT(indexT))then
              Nlow=Nlow+1
              tempgl(Nlow)=atmg(indexT-i)
              ingl(Nlow)=indexT-i
            else
              go to 1500
            endif
 1000     continue
 1500     call locate(tempgl,Nlow,gin,indexgl)
c
c
c   Here is the case when the input log(g) value is in between
c   two tabulated log(g) values:
c
         if((indexgl.lt.Nlow).and.(indexgl.gt.0))then
           gnew(1)=tempgl(indexgl)
           gnew(2)=tempgl(indexgl+1)
         endif 
c
c   If the input log(g) is larger than the largest tabulated log(g), then
c   wing it:  the polint routine will attempt to extrapolate.  This
c   is not too bad if the input log(g) is not too much larger than the
c   largest entry.
c
         if(indexgl.eq.Nlow)then
           gnew(1)=tempgl(Nlow-1)
           gnew(2)=tempgl(Nlow)
         endif 
c
c   Same, but for an input log(g) value smaller than the smallest table
c   entry for that temperature:
c
         if(indexgl.eq.0)then
           gnew(1)=tempgl(1)
           gnew(2)=tempgl(2)
           indexgl=1
         endif 
c
         if(iatm.ge.2)rmuin=1.0d0
         call indexinty(ingl(indexgl),maxlines,maxmu,atmmu,
     &     atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &     Nmu,rmuin,
     %     outinty,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,imuguess)
          do 90 j=1,8
            glinty(1,j)=outinty(j)
 90       continue
c
c          write(*,69)Tin,gin,rmuin,outinty(2)
c
         if(iatm.ge.2)rmuin=1.0d0
         call indexinty(ingl(indexgl+1),maxlines,maxmu,atmmu,
     &     atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &     Nmu,rmuin,
     %     outinty,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,imuguess)
          do 91 j=1,8
            glinty(2,j)=outinty(j)
 91       continue
c
c          write(*,69)Tin,gin,rmuin,outinty(2)
c
69        format(f8.2,2x,f6.3,2x,f8.6,3x,1pe13.5)
          m=2
c
          do 200 i=1,8
            if ((i.eq.1).and.(icnU.eq.430))then
              outinty(i)=1.0d0
              go to 200
            endif
            if ((i.eq.2).and.(icnB.eq.430))then
              outinty(i)=1.0d0
              go to 200
            endif
            if ((i.eq.3).and.(icnV.eq.430))then
              outinty(i)=1.0d0
              go to 200
            endif
            if ((i.eq.4).and.(icnR.eq.430))then
              outinty(i)=1.0d0
              go to 200
            endif
            if ((i.eq.5).and.(icnI.eq.430))then
              outinty(i)=1.0d0
              go to 200
            endif
            if ((i.eq.6).and.(icnJ.eq.430))then
              outinty(i)=1.0d0
              go to 200
            endif
            if ((i.eq.7).and.(icnH.eq.430))then
              outinty(i)=1.0d0
              go to 200
            endif
            if ((i.eq.8).and.(icnK.eq.430))then
              outinty(i)=1.0d0
              go to 200
            endif
            do 190 j=1,2
              yscratch(j)=glinty(j,i)
 190      continue
            diff=dabs(gnew(1)-gnew(2))
            if(diff.gt.1.0d-5)then
c
c              call polint(gnew,yscratch,m,gin,qqq,dy)
c
              qqq=((gnew(2)-gin)*yscratch(1)+(gin-gnew(1))*yscratch(2))
     #             /(gnew(2)-gnew(1))
              outinty(i)=qqq
              tinty(2,i)=qqq
            else
              outinty(i)=yscratch(1)
              tinty(2,i)=yscratch(1)
            endif
 200      continue
c
c   Finally, take the final pass and interpolate between T
c
          do 300 i=1,8
            if ((i.eq.1).and.(icnU.eq.430))then
              outinty(i)=1.0d0
              go to 300
            endif
            if ((i.eq.2).and.(icnB.eq.430))then
              outinty(i)=1.0d0
              go to 300
            endif
            if ((i.eq.3).and.(icnV.eq.430))then
              outinty(i)=1.0d0
              go to 300
            endif
            if ((i.eq.4).and.(icnR.eq.430))then
              outinty(i)=1.0d0
              go to 300
            endif
            if ((i.eq.5).and.(icnI.eq.430))then
              outinty(i)=1.0d0
              go to 300
            endif
            if ((i.eq.6).and.(icnJ.eq.430))then
              outinty(i)=1.0d0
              go to 300
            endif
            if ((i.eq.7).and.(icnH.eq.430))then
              outinty(i)=1.0d0
              go to 300
            endif
            if ((i.eq.8).and.(icnK.eq.430))then
              outinty(i)=1.0d0
              go to 300
            endif
            do 290 j=1,2
              y2scratch(j)=tinty(j,i)
 290        continue
            if(tscratch(1).ne.tscratch(2))then
c
c              call polint(tscratch,y2scratch,m,Tin,qqq,dy)
c
              qqq=((tscratch(2)-Tin)*y2scratch(1)+(Tin-tscratch(1))
     #             *y2scratch(2))/(tscratch(2)-tscratch(1))

              outinty(i)=qqq
            else
              outinty(i)=y2scratch(1)
            endif
 300      continue
c
c   UPDATE October 23, 2009
c
c   If iatm=2, we are using a limb darkening law to find the specific
c   intensity at values of mu different than 1.0.  Apply the limb
c   darkening correction here, based on ilaw and the user supplied
c   coefficients.
c
          if(iatm.ge.2)then
            do 400 i=1,8
              fx=dwavex(i,istar)
              fy=dwavey(i,istar)
              dark=(1.0d0-fx+fx*rsavemu)
              if(ilaw.eq.2)dark=dark-fy*dabs(rsavemu)*dlog(dabs(rsavemu))
              if(ilaw.eq.3)dark=dark-fy*(1.0d0-dsqrt(dabs(rsavemu)))
              if(ilaw.eq.4)dark=dark-fy*(1.0d0-dabs(rsavemu))**2
              outinty(i)=outinty(i)*dark
400         continue
            rmuin=rsavemu
          endif
          return
          end
c
c   &&&&&&&&&&&&&&&&&
c
          subroutine indexinty(index,maxlines,maxmu,atmmu,
     &        atmint1,atmint2,atmint3,atmint4,atmint5,
     #        atmint6,atmint7,atmint8,Nmu,rmuin,
     %        outinty,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,imuguess)
c
c   November 25, 1999
c
c   This subroutine will read the return values of the intensity at the
c   angle rmuin for model index in the table.
c
c
          implicit double precision(a-h,o-z)

          dimension atmmu(maxlines,maxmu),Nmu(maxlines)
          dimension xmu(400),ymu(400),outinty(8)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
c   RVG BUG ALERT  June 13, 2001
c
c   This if-then clause seems to be needed...
c
          if((index.le.0))then
            do 1 i=1,8
              outinty(i)=0.0d0
 1          continue
            return
          endif
c
c          do 2 i=1,128
c            xmu(i)=0.0d0
c2         continue
c
c   Copy the mu values to a one dimensional array
c
          k=8
c          iset=0
          do 10 i=1,Nmu(index)
            xmu(i)=atmmu(index,i)
c            if(iset.eq.1)go to 10
c            if(rmuin.le.xmu(i))then
c              iset=1
c              imuguess=i-1
c             endif
c            qqq=atmint(index,i,k)
c            write(66,6666)xmu(i),qqq,index,i
 10       continue
6666      format(f7.5,3x,1pe14.5,2x,i4,2x,i4)
c             
          N=Nmu(index)
c
          call locate(xmu,N,rmuin,muindex)
c
c          call hunt(xmu,N,rmuin,imuguess)
c
c          muindex=imuguess
c          if(muindex.lt.1)muindex=1
c          if(muindex.gt.N)muindex=N
          m=2
          k=min(max(muindex-(m-1)/2,1),N+1-m)
          if(k.ge.N)k=N-1
          if(k.lt.1)k=1
          do 20 i=1,8
            if((i.eq.1).and.(icnU.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if((i.eq.2).and.(icnB.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if((i.eq.3).and.(icnV.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if((i.eq.4).and.(icnR.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if((i.eq.5).and.(icnI.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if((i.eq.6).and.(icnJ.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if((i.eq.7).and.(icnH.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif
            if((i.eq.8).and.(icnK.eq.430))then
              outinty(i)=1.0d0
              go to 20
            endif

            if(i.eq.1)then
              qqq=((xmu(k+1)-rmuin)*atmint1(index,k)+
     &               (rmuin-xmu(k))*atmint1(index,k+1))/
     $               (xmu(k+1)-xmu(k))
              outinty(i)=dabs(qqq)
            endif
            if(i.eq.2)then
              qqq=((xmu(k+1)-rmuin)*atmint2(index,k)+
     &               (rmuin-xmu(k))*atmint2(index,k+1))/
     $               (xmu(k+1)-xmu(k))
              outinty(i)=dabs(qqq)
            endif
            if(i.eq.3)then
              qqq=((xmu(k+1)-rmuin)*atmint3(index,k)+
     &               (rmuin-xmu(k))*atmint3(index,k+1))/
     $               (xmu(k+1)-xmu(k))
              outinty(i)=dabs(qqq)
            endif
            if(i.eq.4)then
              qqq=((xmu(k+1)-rmuin)*atmint4(index,k)+
     &               (rmuin-xmu(k))*atmint4(index,k+1))/
     $               (xmu(k+1)-xmu(k))
              outinty(i)=dabs(qqq)
            endif
            if(i.eq.5)then
              qqq=((xmu(k+1)-rmuin)*atmint5(index,k)+
     &               (rmuin-xmu(k))*atmint5(index,k+1))/
     $               (xmu(k+1)-xmu(k))
              outinty(i)=dabs(qqq)
            endif
            if(i.eq.6)then
              qqq=((xmu(k+1)-rmuin)*atmint6(index,k)+
     &               (rmuin-xmu(k))*atmint6(index,k+1))/
     $               (xmu(k+1)-xmu(k))
              outinty(i)=dabs(qqq)
            endif
            if(i.eq.7)then
              qqq=((xmu(k+1)-rmuin)*atmint7(index,k)+
     &               (rmuin-xmu(k))*atmint7(index,k+1))/
     $               (xmu(k+1)-xmu(k))
              outinty(i)=dabs(qqq)
            endif
            if(i.eq.8)then
              qqq=((xmu(k+1)-rmuin)*atmint8(index,k)+
     &               (rmuin-xmu(k))*atmint8(index,k+1))/
     $               (xmu(k+1)-xmu(k))
              outinty(i)=dabs(qqq)
            endif

c           if(outinty(i).lt.0.0d0)outinty(i)=0.0d0
 20       continue
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&
c
          subroutine loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &       atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &       Tmax,Tmin,gmax,gmin)
c
c   November 25, 1999
c
c   This routine will read the file with the model atmosphere data.  The
c   name is assumed to be 'ELC.atm', the the form is assumed to be the
c   following:
c
c   Teff   g
c   N_mu
c   mu1    intyU intyB intyV intyR intyI intyJ intyH intyK
c   mu2    intyU intyB intyV intyR intyI intyJ intyH intyK
c   ...
c
c
c   UPDATE January 9, 2009
c
c   make the variable atmint two dimensional, and have 8 copies called
c   atmint1, atmint2, ...
c
c
          implicit double precision(a-h,o-z)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)

          character*40 blank
c
c   UPDATE June 14, 2002
c
c   Declare the variable bell to be character*1
c
          character*1 bell
c
          ios=0
          open(unit=19,file='ELC.atm',status='old',err=999,iostat=ios)
c
c   Read the header(s).
c
          nskip=0
          do 5 i=1,100
            read(19,6)blank
            if(blank(1:2).eq.'#'.or.blank(1:2).eq.'!'.or.blank(1:2).eq.'%')
     %        nskip=nskip+1
 5        continue
 6        format(a40)
c
          rewind(19)
c
          if(nskip.le.0)go to 8
          do 7 i=1,Nskip
            read(19,6)blank
 7        continue
c
c   Attempt to read in the intensity values.
c
 8        Tmax=-1000.0d0
          Tmin=10000000.0d00
          gmax=-10000.0d0
          gmin=11111.0d0
          do 10 i=1,maxlines
            read(19,*,end=15)atmT(i),atmg(i)
            read(19,*,end=15)Nmu(i)
            if(atmT(i).gt.Tmax)Tmax=atmT(i)
            if(atmT(i).lt.Tmin)Tmin=atmT(i)
            if(atmg(i).gt.gmax)gmax=atmg(i)
            if(atmg(i).lt.gmin)gmin=atmg(i)
            do 9 j=1,Nmu(i)
c              read(19,*,end=15)atmmu(i,j),(atmint(i,j,k),k=1,8)
              read(19,*,end=15)atmmu(i,j),atmint1(i,j),atmint2(i,j),
     $          atmint3(i,j),atmint4(i,j),atmint5(i,j),atmint6(i,j),
     &          atmint7(i,j),atmint8(i,j)

c
c   Add this line to check for repeated mu values.  If there 
c   is a repeat, then add 0.0001 (February 9, 2000).
c
              if((j.gt.1).and.(atmmu(i,j).eq.atmmu(i,j-1)))then
c                atmmu(i,j)=atmmu(i,j)+0.0001d0
                 atmmu(i,j-1)=atmmu(i,j-1)-0.00005
              endif
 9          continue
 10       continue
 15       close(19)
c
          Nlines=i-1
c
 999      if(ios.ne.0)then
            write(*,100)bell
          endif
c
 100      format(a1,'Error:  I can''t find the file ''ELC.atm''!')
c

          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
c
c   Taken from Numerical Recipes.
c
          implicit double precision(a-h,o-z)
c
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=DABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.0d0)then
             write(*,*)'pause in polint ',den,x,xa(i),ya(i)
          endif
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
c
c  &&&&&&&&&&&&&&&&&
c
          subroutine getATMflux(ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     $      visib,projarray,temp,surf,garray,rinty,
     &      flum,maxlines,maxmu,Nlines,
     &      atmT,atmg,atmmu,Nmu,
     &      atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &      Tmax,Tmin,gmax,gmin,gscale,
     &      fluxU,fluxB,fluxV,fluxR,fluxI,fluxJ,fluxH,fluxK,
     $      icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,delphi,
     %      delphiedge,iedgestar,iedgehor,darkint,separ,mmdx,MonteCarlo,
     &      dwavex,dwavey,ilaw,iatm,istar)
c
c  November 30, 1999
c
c  This routine will return fluxes computed from model atmospheres.
c  The flum  and rinty arrays will contain the intensities for the V band.
c  The rinty array contains the specific intensities, while the
c  flum array contains the intensities weighted by the area and the surf
c  vectors.
c
c  The projarray contains the cosine mu terms for each element.  The visib
c  array contains the cosine mu terms for each element, except if the point
c  is eclipsed in which case the visib=0.
c 
c  The parameter gscale is used to convert the gravities in program units
c  into cgs units.  This number is G*M/(a*a).
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.  Then scale the fluxes
c   by (separ*solarrad)**2
c
c   
c   UPDATE JULY 4, 2004
c
c   Add the variable MonteCarlo to the argument list.  If MonteCarlo < 10,
c   then proceed as before.  If Monte Carlo > 10, then the fractional
c   pixels were computed in getvisib via Monte Carlo integration.  In
c   that case, we can skip some steps below.
c
c   UPDATE SEPTEMBER 11, 2009
c
c   Add the ability to use a parameterized limb darkening law to this
c   routine.  computeinty is modified so that when iatm=2, the flux
c   at mu=1 is found, then I(mu)=I_0*ld_law
c

          implicit double precision(a-h,o-z)

          parameter(pie=3.141592653589793d0)
c
c   Set these to the value of ialphmax,ibetmax
c
          integer tempalf,tempbet
          parameter(tempalf=3000,tempbet=3000,itab=tempalf*tempbet)

          dimension visib(ialphmax*ibetmax),ibetlim(ialphmax),
     $        surf(ialphmax*ibetmax),garray(ialphmax*ibetmax),
     $        temp(ialphmax*ibetmax),flum(ialphmax*ibetmax),
     $        rinty(ialphmax*ibetmax),projarray(ialphmax*ibetmax),
     $        delphi(ialphmax*ibetmax),delphiedge(ialphmax*ibetmax),
     $        iedgestar(ialphmax*ibetmax),iedgehor(ialphmax*ibetmax),
     &        mmdx(ialphmax,ibetmax)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines),outinty(8),
     #       corr1(8),corr2(8),saveflum(itab,8),darkint(8)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
          dimension dwavex(8,2),dwavey(8,2)
c
c  
c
          if(tempalf.lt.ialphmax)then
            write(*,*)'dimension error in getATMflux'
            write(*,*)'tempalf = ',tempalf,'  ialphmax = ',ialphmax
            stop
          endif
          if(tempbet.lt.ibetmax)then
            write(*,*)'dimension error in getATMflux'
            write(*,*)'tempbet = ',tempbet,'  ibetmax = ',ibetmax
            stop
          endif
c
          fluxU=0.0d0
          fluxB=0.0d0
          fluxV=0.0d0
          fluxR=0.0d0
          fluxI=0.0d0
          fluxJ=0.0d0
          fluxH=0.0d0
          fluxK=0.0d0
          corr1(1)=0.0d0
          corr1(2)=0.0d0
          corr1(3)=0.0d0
          corr1(4)=0.0d0
          corr1(5)=0.0d0
          corr1(6)=0.0d0
          corr1(7)=0.0d0
          corr1(8)=0.0d0
          corr2(1)=0.0d0
          corr2(2)=0.0d0
          corr2(3)=0.0d0
          corr2(4)=0.0d0
          corr2(5)=0.0d0
          corr2(6)=0.0d0
          corr2(7)=0.0d0
          corr2(8)=0.0d0
c
c   Initialize the flum matrix.
c
          do 2 ialf=1,nalf
            do 1 ibet=1,ibetlim(ialf)       !4*Nbet
c              iidx=(ialf-1)*4*Nbet+ibet
              iidx=mmdx(ialf,ibet)
              flum(iidx)=0.0d0
              rinty(iidx)=0.0d0
              saveflum(iidx,1)=0.0d0
              saveflum(iidx,2)=0.0d0
              saveflum(iidx,3)=0.0d0
              saveflum(iidx,4)=0.0d0
              saveflum(iidx,5)=0.0d0
              saveflum(iidx,6)=0.0d0
              saveflum(iidx,7)=0.0d0
              saveflum(iidx,8)=0.0d0
 1          continue
 2        continue
c
          Nalf2=Nalf/2
c
c   Find the rough place in the atmosphere table.
c
          Tin=temp(1)
          call locate(atmT,Nlines,Tin,indexT)
          itguess=indexT
          imuguess=1
c
c   Loop for fractional pixels near the edge.
c
c
c   UPDATE JULY 4, 2004
c
c   if MonteCarlo > 10, we can skip this loop.
c
          if(MonteCarlo.lt.10)then
            do 4 ialf=1,nalf
              do 3 ibet=1,ibetlim(ialf)
c                iidx=(ialf-1)*ibetlim(ialf)+ibet
c                iidx=kount(ialphmax,ialf,ibetlim)+ibet
                iidx=mmdx(ialf,ibet)
                if((iedgehor(iidx).eq.-10).or.(iedgehor(iidx).gt.
     #               5).or.(iedgestar(iidx).eq.-10).
     $               or.(iedgestar(iidx).gt.5).or.(delphi(iidx).
     $               gt.-10.0d0))then
c
                  Tin=temp(iidx)
                  gin=dlog10(gscale*garray(iidx))
                  rmuin=dabs(projarray(iidx))
c                  write(*,*)gin,istar
c  
                  call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &             atmT,atmg,atmmu,Nmu,
     &             atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &             atmint7,atmint8,Tmax,Tmin,gmax,gmin,outinty,
     %             icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     &             dwavex,dwavey,ilaw,iatm,istar)
c
                  do 88 k=1,8
                    saveflum(iidx,k)=outinty(k)*surf(iidx)*
     #                 dabs(projarray(iidx))
 88               continue
                endif
 3            continue
 4          continue
          endif
c
          DO 10 ialf=1,nalf
            DO 9 ibet = 1,ibetlim(ialf)      !4*nbet
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
              iidx=mmdx(ialf,ibet)
              corr1(1)=0.0d0
              corr1(2)=0.0d0
              corr1(3)=0.0d0
              corr1(4)=0.0d0
              corr1(5)=0.0d0
              corr1(6)=0.0d0
              corr1(7)=0.0d0
              corr1(8)=0.0d0
              corr2(1)=0.0d0
              corr2(2)=0.0d0
              corr2(3)=0.0d0
              corr2(4)=0.0d0
              corr2(5)=0.0d0
              corr2(6)=0.0d0
              corr2(7)=0.0d0
              corr2(8)=0.0d0
              dphi=pie/dble(ibetlim(ialf))
              if(projarray(iidx).le.0.0d0)go to 9
              Tin=temp(iidx)
              gin=dlog10(gscale*garray(iidx))
              rmuin=projarray(iidx)
c
c                write(*,*)ialf,ibet,Tin,gin
c
              call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,
     &           atmint5,atmint6,atmint7,atmint8,
     &           Tmax,Tmin,gmax,gmin,outinty,
     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     #           dwavex,dwavey,ilaw,iatm,istar)
c
              rinty(iidx)=outinty(iRVfilt) ! save intensities for plotting
c
              do 8 k=1,8
                saveflum(iidx,k)=
     #             outinty(k)*surf(iidx)*projarray(iidx)
                outinty(k)=outinty(k)*surf(iidx)*visib(iidx)
 8            continue
c
c   Check for fractional pixels near the horizon of the star in front,
c   in the beta direction (along constant latitude rows).
c
c
c   UPDATE JULY 4, 2004
c
c   If MonteCarlo > 10, we can skip this step since the fractional
c   pixel corrections were done in getvisib.
c
              if(MonteCarlo.gt.10)then
                corr1(1)=0.0d0
                corr1(2)=0.0d0
                corr1(3)=0.0d0
                corr1(4)=0.0d0
                corr1(5)=0.0d0
                corr1(6)=0.0d0
                corr1(7)=0.0d0
                corr1(8)=0.0d0
                go to 867
              endif
c
              if(delphi(iidx).ge.-10.0d0)then
                frac=0.5d0*(dabs(delphi(iidx))-dphi)/dphi
                if(frac.lt.0.0d0)then
                  corr1(1)=frac*saveflum(iidx,1)
                  corr1(2)=frac*saveflum(iidx,2)
                  corr1(3)=frac*saveflum(iidx,3)
                  corr1(4)=frac*saveflum(iidx,4)
                  corr1(5)=frac*saveflum(iidx,5)
                  corr1(6)=frac*saveflum(iidx,6)
                  corr1(7)=frac*saveflum(iidx,7)
                  corr1(8)=frac*saveflum(iidx,8)
                else
                  if(iedgehor(iidx).eq.10)then
                    if(ibet.lt.ibetlim(ialf))then  
                      izz=ialf
                      jzz=ibet+1
c                      iidx=(izz-1)*ibetlim(ialf)+jzz
c                      iidx=kount(ialphmax,izz,ibetlim)+jzz
                      iidx=mmdx(izz,jzz)
                      corr1(1)=frac*saveflum(iidx,1)
                      corr1(2)=frac*saveflum(iidx,2)
                      corr1(3)=frac*saveflum(iidx,3)
                      corr1(4)=frac*saveflum(iidx,4)
                      corr1(5)=frac*saveflum(iidx,5)
                      corr1(6)=frac*saveflum(iidx,6)
                      corr1(7)=frac*saveflum(iidx,7)
                      corr1(8)=frac*saveflum(iidx,8)
                    else
                      izz=ialf
                      jzz=1
c                      iidx=(izz-1)*ibetlim(ialf)+jzz
c                      iidx=kount(ialphmax,izz,ibetlim)+jzz
                      iidx=mmdx(izz,jzz)
                      corr1(1)=frac*saveflum(iidx,1)
                      corr1(2)=frac*saveflum(iidx,2)
                      corr1(3)=frac*saveflum(iidx,3)
                      corr1(4)=frac*saveflum(iidx,4)
                      corr1(5)=frac*saveflum(iidx,5)
                      corr1(6)=frac*saveflum(iidx,6)
                      corr1(7)=frac*saveflum(iidx,7)
                      corr1(8)=frac*saveflum(iidx,8)
                    endif
                  endif

                  izz=ialf
                  jzz=ibet
c                  iidx=(izz-1)*ibetlim(ialf)+jzz
c                  iidx=kount(ialphmax,izz,ibetlim)+jzz
                  iidx=mmdx(izz,jzz)
                  if(iedgehor(iidx).eq.20)then
                    if(ibet.gt.1)then    
                      izz=ialf
                      jzz=ibet-1
c                      iidx=(izz-1)*ibetlim(ialf)+jzz
c                      iidx=kount(ialphmax,izz,ibetlim)+jzz
                      iidx=mmdx(izz,jzz)
                      corr1(1)=frac*saveflum(iidx,1)
                      corr1(2)=frac*saveflum(iidx,2)
                      corr1(3)=frac*saveflum(iidx,3)
                      corr1(4)=frac*saveflum(iidx,4)
                      corr1(5)=frac*saveflum(iidx,5)
                      corr1(6)=frac*saveflum(iidx,6)
                      corr1(7)=frac*saveflum(iidx,7)
                      corr1(8)=frac*saveflum(iidx,8)
                    else
                      izz=ialf
                      jzz=ibetlim(ialf)
c                      iidx=(izz-1)*ibetlim(ialf)+jzz
c                      iidx=kount(ialphmax,izz,ibetlim)+jzz
                      iidx=mmdx(izz,jzz)
                      corr1(1)=frac*saveflum(iidx,1) 
                      corr1(2)=frac*saveflum(iidx,2) 
                      corr1(3)=frac*saveflum(iidx,3) 
                      corr1(4)=frac*saveflum(iidx,4) 
                      corr1(5)=frac*saveflum(iidx,5) 
                      corr1(6)=frac*saveflum(iidx,6) 
                      corr1(7)=frac*saveflum(iidx,7) 
                      corr1(8)=frac*saveflum(iidx,8) 
                    endif
                  endif
                endif
              endif
c
 6969         format(a10,e16.7,2x,2(i2,2x),f9.6)
c
              izz=ialf
              jzz=ibet
c              iidx=(izz-1)*ibetlim(ialf)+jzz
c              iidx=kount(ialphmax,izz,ibetlim)+jzz
               iidx=mmdx(izz,jzz)
              if(iedgestar(iidx).eq.10)then
                frac=0.5d0*(dabs(delphiedge(iidx))-dphi)/dphi
                if(frac.lt.0.0d0)then
                  izz=ialf
                  jzz=ibet
c                  iidx=(izz-1)*ibetlim(ialf)+jzz
                  iidx=mmdx(izz,jzz)
                  corr2(1)=frac*saveflum(iidx,1)
                  corr2(2)=frac*saveflum(iidx,2)
                  corr2(3)=frac*saveflum(iidx,3)
                  corr2(4)=frac*saveflum(iidx,4)
                  corr2(5)=frac*saveflum(iidx,5)
                  corr2(6)=frac*saveflum(iidx,6)
                  corr2(7)=frac*saveflum(iidx,7)
                  corr2(8)=frac*saveflum(iidx,8)
                else
                  if(ibet.lt.ibetlim(ialf))then
                    izz=ialf
                    jzz=ibet+1
c                    iidx=(izz-1)*ibetlim(ialf)+jzz
c                    iidx=kount(ialphmax,izz,ibetlim)+jzz
                    iidx=mmdx(izz,jzz)
                    corr2(1)=frac*saveflum(iidx,1)
                    corr2(2)=frac*saveflum(iidx,2)
                    corr2(3)=frac*saveflum(iidx,3)
                    corr2(4)=frac*saveflum(iidx,4)
                    corr2(5)=frac*saveflum(iidx,5)
                    corr2(6)=frac*saveflum(iidx,6)
                    corr2(7)=frac*saveflum(iidx,7)
                    corr2(8)=frac*saveflum(iidx,8)
                  else
                    izz=ialf
                    jzz=1
c                    iidx=(izz-1)*ibetlim(ialf)+jzz
c                    iidx=kount(ialphmax,izz,ibetlim)+jzz
                    iidx=mmdx(izz,jzz)
                    corr2(1)=frac*saveflum(iidx,1)
                    corr2(2)=frac*saveflum(iidx,2)
                    corr2(3)=frac*saveflum(iidx,3)
                    corr2(4)=frac*saveflum(iidx,4)
                    corr2(5)=frac*saveflum(iidx,5)
                    corr2(6)=frac*saveflum(iidx,6)
                    corr2(7)=frac*saveflum(iidx,7)
                    corr2(8)=frac*saveflum(iidx,8)
                  endif
                endif
              endif
c
              izz=ialf
              jzz=ibet
c              iidx=(izz-1)*ibetlim(ialf)+jzz
c              iidx=kount(ialphmax,izz,ibetlim)+jzz
              iidx=mmdx(izz,jzz)
              if(iedgestar(iidx).eq.20)then
                frac=0.5d0*(dabs(delphiedge(iidx))-dphi)/dphi
                if(frac.lt.0.0d0)then
                  corr2(1)=frac*saveflum(iidx,1)
                  corr2(2)=frac*saveflum(iidx,2)
                  corr2(3)=frac*saveflum(iidx,3)
                  corr2(4)=frac*saveflum(iidx,4)
                  corr2(5)=frac*saveflum(iidx,5)
                  corr2(6)=frac*saveflum(iidx,6)
                  corr2(7)=frac*saveflum(iidx,7)
                  corr2(8)=frac*saveflum(iidx,8)
                else
                  if(ibet.gt.1)then
                    izz=ialf
                    jzz=ibet-1
c                    iidx=(izz-1)*ibetlim(ialf)+jzz
c                    iidx=kount(ialphmax,izz,ibetlim)+jzz
                    iidx=mmdx(izz,jzz)
                    corr2(1)=frac*saveflum(iidx,1)
                    corr2(2)=frac*saveflum(iidx,2)
                    corr2(3)=frac*saveflum(iidx,3)
                    corr2(4)=frac*saveflum(iidx,4)
                    corr2(5)=frac*saveflum(iidx,5)
                    corr2(6)=frac*saveflum(iidx,6)
                    corr2(7)=frac*saveflum(iidx,7)
                    corr2(8)=frac*saveflum(iidx,8)
                  else
                    izz=ialf
                    jzz=ibetlim(ialf)
c                    iidx=(izz-1)*ibetlim(ialf)+jzz
c                    iidx=kount(ialphmax,izz,ibetlim)+jzz
                    iidx=mmdx(izz,jzz)
                    corr2(1)=frac*saveflum(iidx,1)
                    corr2(2)=frac*saveflum(iidx,2)
                    corr2(3)=frac*saveflum(iidx,3)
                    corr2(4)=frac*saveflum(iidx,4)
                    corr2(5)=frac*saveflum(iidx,5)
                    corr2(6)=frac*saveflum(iidx,6)
                    corr2(7)=frac*saveflum(iidx,7)
                    corr2(8)=frac*saveflum(iidx,8)
                  endif
                endif
              endif
c
 867          fluxU=fluxU+outinty(1)+corr1(1) !+corr2(1)
              fluxB=fluxB+outinty(2)+corr1(2)!+corr2(2)
              fluxV=fluxV+outinty(3)+corr1(3)!+corr2(3)
              fluxR=fluxR+outinty(4)+corr1(4)!+corr2(4)
              fluxI=fluxI+outinty(5)+corr1(5)!+corr2(5)
              fluxJ=fluxJ+outinty(6)+corr1(6)!+corr2(6)
              fluxH=fluxH+outinty(7)+corr1(7)!+corr2(7)
              fluxK=fluxK+outinty(8)+corr1(8)!+corr2(8)
              izz=ialf
              jzz=ibet
c              iidx=(izz-1)*ibetlim(ialf)+jzz
c              iidx=kount(ialphmax,izz,ibetlim)+jzz
               iidx=mmdx(izz,jzz)
              flum(iidx)=outinty(iRVfilt)+corr1(iRVfilt) !+corr2(iRVfilt) 
 9          continue
 10       continue
c
c          if(darkint(1).ne.0.0d0)fluxU=pie*fluxU/darkint(1)
c          if(darkint(2).ne.0.0d0)fluxB=pie*fluxB/darkint(2)
c          if(darkint(3).ne.0.0d0)fluxV=pie*fluxV/darkint(3)
c          if(darkint(4).ne.0.0d0)fluxR=pie*fluxR/darkint(4)
c          if(darkint(5).ne.0.0d0)fluxI=pie*fluxI/darkint(5)
c          if(darkint(6).ne.0.0d0)fluxJ=pie*fluxJ/darkint(6)
c          if(darkint(7).ne.0.0d0)fluxH=pie*fluxH/darkint(7)
c          if(darkint(8).ne.0.0d0)fluxK=pie*fluxK/darkint(8)
c
c  UPDATE April 3, 2002
c
c  Scale the fluxes.
c
          solarrad=6.9598d10
          fluxU=fluxU*(separ*solarrad)**2
          fluxB=fluxB*(separ*solarrad)**2
          fluxV=fluxV*(separ*solarrad)**2
          fluxR=fluxR*(separ*solarrad)**2
          fluxI=fluxI*(separ*solarrad)**2
          fluxJ=fluxJ*(separ*solarrad)**2
          fluxH=fluxH*(separ*solarrad)**2
          fluxK=fluxK*(separ*solarrad)**2
c
          return
          end
c
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine parms(iunit,
     $        Teff2,Q,
     %        finc,separ,period,reff1,reff2,rpole1,rpole2,fill1,fill2,
     $        gp1,gp2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,bdist,
     $        ecc)
c
c   November 30, 1999
c
c   This routine will compute the component masses, radii, etc. based
c   on the mass ratio, inclination, orbital separation, and orbital period.
c   Set iunit=1 to print to the output file, or 0 to print to screen.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
           
           fincr=finc*pie/180.0d0        ! radians
           ppp=period*24.0d0                            ! period in hours
           coef=3.518847d10
           coef3=4.35713636d31
           sifinc=dsin(fincr)
c
c   Use the formula separ = coef*(perid*period*total_mass)**(1/3) to
c   solve for the total mass in solar masses.  The separation is
c   entered in solar masses, so (R_sun/coef)**3=7.737294491.
c
           total_mass=(separ)**(3)*7.737294491d0/(ppp*ppp)
c
           rM1=total_mass/(1.0d0+Q)
           rM2=Q*rM1
c
c   Use the value of GM_sun found from the solar system.
c
           gmsun=1.32712440018d20  !mks units
c
           solarrad=6.9598d8
           p=period*86400.0d0
           rM1=(separ*solarrad)**(3)*4.0d0*pie*pie
           rM1=rM1/(gmsun*p*p*(1.0d0+Q))
           rM2=Q*rM1

           R1=reff1*separ
           if(Teff2.gt.0.0d0)then
             R2=reff2*separ
           else
             R2=0.00d0
           endif
c
           gsun=2.739910d4
           gpole1=gsun*rM1/(R1*R1)  
           if(Teff2.gt.0.0d0)then
             gpole2=gsun*rM2/(R2*R2)  
           else
             gpole2=1.0d0
           endif
c
           gscale1=27397.726d0*rM1/(separ*separ*bdist*bdist)
           gscale2=27397.726d0*rM2/(separ*separ*bdist*bdist)

           gscale1=27397.726d0*rM1/(separ*separ)
           gscale2=27397.726d0*rM2/(separ*separ)
c
           fact=solarrad*2.0d0*pie/86400.0d0/1.0d3  !50.613093d0
c
           vrot1=fact*R1/period*sifinc
           vrot2=fact*R2/period*sifinc
c
c   UPDATE MARCH 4, 2005
c
c   Was (separ/(1.0d0+Q))*bdist)
c
           a2=(separ/(1.0d0+Q))
c
c   UPDATE September 12, 2001
c 
c   Bug fix, change a1=(separ-a2)*bdist  to  a1=separ(1.0d0-bdist/(1.0d0+Q))
c

           aa1=separ*(1.0d0-bdist/(1.0d0+Q))
           a1=separ-a2
c
c   RVG BUG ALERT  April 20, 2001
c
c   Change efact below
c
c           efact=1.0d0/dsqrt(1.0-ecc*ecc)
c
           efact=1.0d0/dsqrt(1.0-ecc*ecc)
           velK1=fact*a1/period*sifinc*efact
           velK2=fact*a2/period*sifinc*efact
c
c
c   NEW BUG August 10, 2001
c
c   Add a correction factor to the rotational velocities in the
c   case of eccentric orbits.
c
           hutfac=(1.0d0+7.5d0*ecc*ecc+5.625d0*ecc**4+
     #          0.3125d0*ecc**6)/((1.0d0+3.0d0*ecc*ecc+
     $          3.0d0/8.0d0*ecc**4)*dsqrt((1.0d0-ecc*ecc)**3))


           if(iunit.ge.1)write(2,100)rM1,R1,dlog10(gpole1),rM2,R2,
     $        dlog10(gpole2),
     $        period,a1,a2,separ,velK1,velK2,vrot1,vrot2,
     &        omega1*vrot1*hutfac,omega2*vrot2*hutfac
           if(iunit.le.0)write(*,100)rM1,R1,dlog10(gpole1),rM2,R2,
     $        dlog10(gpole2),
     $        period,a1,a2,separ,velK1,velK2,vrot1,vrot2,
     &        omega1*vrot1*hutfac,omega2*vrot2*hutfac

c
 100       format(/'M1 = ',f6.3,' M_sun;',1x,'R1 = ',f7.3,' R_sun;',1x,
     %          'log(g1) = ',f5.3,' cgs', 
     &        /'M2 = ',f6.3,' M_sun;',1x,'R2 = ',f7.3,' R_sun;',1x,
     %          'log(g2) = ',f5.3,' cgs', 
     %        /'P = ',f11.6,' d;',1x,'a1 = ',f7.3,' R_sun;',1x,
     $        'a2 = ',f8.3,' R_sun;',1x,'a = ',f8.3,' R_sun',
     %        /'K1 = ',f8.3,' km/sec;',1x,'K2 = ',f8.3,' km/sec;',1x,
     &        /'V1_rot*sin(i) = ',f8.3,' km/s;',1x,'V2_rot*sin(i) = ',
     %        f8.3,' km/s',/'V1_rot*sin(i) = ',f8.3,' km/s;',
     %        1x,'V2_rot*sin(i) = ',
     %        f8.3,' km/s (scaled by omegas)')
c
           return
           end
c
c
c  %%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine wmagmod(icount,xmod,ymod,fileout,isw7)
c
c  November 10, 1999
c
c  This routine will open the output file and record the given model
c  with flux in linear units.  The phases will be scaled to go from 0.0 to
c  1.0, and an extra phase will be added.
c
          implicit double precision(a-h,o-z)
c
          dimension xmod(icount),ymod(icount)
c
c   UPDATE June 17, 2002
c
c   Change the declaration of fileout to character*(*)
c
          character*(*) fileout
c
          open(unit=20,file=fileout,status='unknown')
c
c   UPDATE April 3, 2002
c
c   Add a variable called zeropoint, instead of 40.0d0
c
          zeropoint=75.0d0
          err=0.005d0
          do 10 i=1,icount
            write(20,100)xmod(i),-2.5d0*dlog10(ymod(i))+zeropoint,err
 10       continue
c
          if(isw7.ge.2)return
          do 20 i=1,icount
            write(20,100)xmod(i)+1.0d0,-2.5d0*dlog10(ymod(i))+zeropoint,err
 20       continue
c
          close(20)
c
 100      format(f16.10,3x,f13.9,3x,f9.6)
c
          return
          end
c
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine  getdiskATMflux(Nrmax,Nthetamax,Nradius,
     %            Ntheta,diskproj,edgeproj,dvisib,evisib,dtemp,tedge,drad,
     $            dinty,einty,stepr,stepz,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,atmint5,
     #           atmint6,atmint7,atmint8,
     #           Tmax,Tmin,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           dfluxU,dfluxB,dfluxV,dfluxR,
     &           dfluxI,dfluxJ,dfluxH,dfluxK,iRVfilt,separ,
     &           dwavex,dwavey,ilaw,iatm,istar)
c
c
c  November 30, 1999
c
c  This routine will return fluxes computed from model atmospheres.
c  The flum  and rinty arrays will contain the intensities for the V band.
c  The rinty array contains the specific intensities, while the
c  flum array contains the intensities weighted by the area and the surf
c  vectors.
c
c  The projarray contains the cosine mu terms for each element.  The visib
c  array contains the cosine mu terms for each element, except if the point
c  is eclipsed in which case the visib=0.
c 
c
c   UPDATE April 3, 2002
c
c   Add separ to the argument list of getBBflux, getATMflux,
c   getdiskBBflux, getdiskATMflux, and getBBsimp.  Then scale the fluxes
c   by (separ*solarrad)**2
c
c   UPDATE SEPTEMBER 11, 2009
c
c   If iatm=2, then use a parameterized limb darkening law to 
c   compute the flux at angles less than mu=1.
c
          implicit double precision(a-h,o-z)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines),outinty(8)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)

          dimension diskproj(Nrmax*Nthetamax),edgeproj(Nthetamax*11),
     &      dvisib(Nrmax*Nthetamax),evisib(Nthetamax*11),
     &      dinty(Nrmax*Nthetamax),einty(Nthetamax*11),
     &      dtemp(Nrmax*Nthetamax),drad(Nrmax),
     &      tedge(Nthetamax*11)
c
          dimension dwavex(8,2),dwavey(8,2)

          parameter(pie=3.14159265358979323d0)
c
          isaveatm=iatm
          if(iatm.ge.2)iatm=1
          dfluxU=0.0d0
          dfluxB=0.0d0
          dfluxV=0.0d0
          dfluxR=0.0d0
          dfluxI=0.0d0
          dfluxJ=0.0d0
          dfluxH=0.0d0
          dfluxK=0.0d0
c
          steptheta=360.0d0/dble(ntheta)
          deg2rad=pie/180.0d0
c
c   Find the rough place in the table
c
          Tin=dtemp(1)
          call locate(atmT,Nlines,Tin,indexT)
          itguess=indexT
          imuguess=1          

          DO 10 ir=1,Nradius
            DO 9 ithet=1,Ntheta
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(ir-1)*Ntheta+ithet
c
              if(diskproj(iidx).le.0.0d0)go to 9
              Tin=dtemp(iidx)
              if(Tin.lt.Tmin)Tin=Tmin+1.0
              if((Tin.lt.6800.0d0))then
                gin=3.9d0   !UPDATE AUG-04-2008  was 4.9
              else
                gin=3.9d0   !UPDATE AUG-04-2008  was 4.9
              endif
              rmuin=diskproj(iidx)
c
              call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &           atmint7,atmint8,Tmax,Tmin,gmax,gmin,outinty,
     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     &           dwavex,dwavey,ilaw,iatm,1)
c
              dinty(iidx)=outinty(iRVfilt)  ! save intensities for plotting
c
              if(dtemp(iidx).lt.Tmin)then
                scale=(dtemp(iidx)/Tmin)**4
              else
                scale=1.0d0
              endif
c
c   April 18, 2001  RGV BUG ALERT
c
c   Convert the steptheta into radians!
c
              do 8 k=1,8
                outinty(k)=outinty(k)*((1.0d0/8.0d0)*drad(ir)**3)*
     &                 dvisib(iidx)*stepr*steptheta*scale*deg2rad
 8            continue
c
              dfluxU=dfluxU+outinty(1)
              dfluxB=dfluxB+outinty(2)
              dfluxV=dfluxV+outinty(3)
              dfluxR=dfluxR+outinty(4)
              dfluxI=dfluxI+outinty(5)
              dfluxJ=dfluxJ+outinty(6)
              dfluxH=dfluxH+outinty(7)
              dfluxK=dfluxK+outinty(8)
 9          continue
 10       continue
c
          do 20 ithet=1,Ntheta
            do 19 iz=1,11
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(iz-1)*Ntheta+ithet
c
              einty(iidx)=0.0d0
              if(edgeproj(iidx).lt.0.0)go to 19
              Tin=tedge(iidx)
              if(Tin.lt.Tmin)Tin=Tmin+1.0
              if((Tin.lt.6800.0d0))then
                gin=3.9d0               !UPDATE AUG-04-2008  was 4.9
              else
                gin=3.9d0               !UPDATE AUG-04-2008  was 4.9
              endif
              rmuin=edgeproj(iidx)
c
              call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &           atmint7,atmint8,
     &           Tmax,Tmin,gmax,gmin,outinty,
     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     &           dwavex,dwavey,ilaw,iatm,1)
c
              einty(iidx)=outinty(iRVfilt)  ! save intensities for plotting
c
              if(tedge(iidx).lt.Tmin)then
                scale=(tedge(iidx)/Tmin)**4
              else
                scale=1.0d0
              endif
c
c   April 18, 2001  RGV BUG ALERT
c
c   Convert the steptheta into radians!
c
              do 18 k=1,8
                outinty(k)=outinty(k)*
     %             evisib(iidx)*stepz*steptheta*scale*deg2rad
 18           continue
c
c   The upper rim points are in the face integration.  If iz=11, then don't
c   add to dflux.
c
              if(iz.lt.11)then          
                dfluxU=dfluxU+outinty(1)
                dfluxB=dfluxB+outinty(2)
                dfluxV=dfluxV+outinty(3)
                dfluxR=dfluxR+outinty(4)
                dfluxI=dfluxI+outinty(5)
                dfluxJ=dfluxJ+outinty(6)
                dfluxH=dfluxH+outinty(7)
                dfluxK=dfluxK+outinty(8)
              endif
 19         continue
 20       continue
c
c   UPDATE April 3, 2002
c
c   Scale the fluxes.
c
          solarrad=6.9598d10
          dfluxU=dfluxU*(separ*solarrad)**2
          dfluxB=dfluxB*(separ*solarrad)**2
          dfluxV=dfluxV*(separ*solarrad)**2
          dfluxR=dfluxR*(separ*solarrad)**2
          dfluxI=dfluxI*(separ*solarrad)**2
          dfluxJ=dfluxJ*(separ*solarrad)**2
          dfluxH=dfluxH*(separ*solarrad)**2
          dfluxK=dfluxK*(separ*solarrad)**2
c
          if(isaveatm.ge.2)iatm=isaveatm
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine thirdlight(iatm,t3,g3,SA3,SA1,third,
     %       maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &       atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &       Tmax,Tmin,gmax,gmin,
     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,separ,
     %       dwavex,dwavey,ilaw,istar)
c
c  January 7, 2000
c 
c  This subroutine will compute a 'third light' flux using the
c  model atmosphere table.  The user specifies the temperature,
c  gravity, and fractional surface area of the third star.  The
c  fractional surface area is given in terms of the ratio of the
c  surface area of star 3 to that of star 1.  Finally, the third light
c  that is added to the output light curves is given by
c
c  third(ifilt) = SA3*inty(ifilt)*2*pi*pi*SA1)
c
c  
c  December 19, 2000
c
c  Bug fix:  third(ifilt) = SA1*SA3*inty(ifilt,mu=1)/4.0
c
c  UPDATE April 2, 2002
c
c  Add separ to the argument list of thirdlight.  Then scale the fluxes
c  by (separ*solarrad)**2
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
          dimension third(8)
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %      Nmu(maxlines),outinty(8)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
          dimension dwavex(8,2),dwavey(8,2)

c
c   Initialize
c
          do 10 i=1,8
            third(i)=0.0d0
 10       continue
c
c   Check for negative values of the input parameters.  If
c   any are negative, then the third light option is turned off.
c
          if(t3.le.0.0d0)then
            write(2,100)
            return
          endif
c
          if(g3.le.0.0d0)then
            write(2,100)
            return
          endif
c
          if(SA3.le.0.0d0)then
            write(2,100)
            return
          endif
c
c   If the atmosphere option is not on, then set the third light
c   to zero.
c 
          if(iatm.le.0d0)then
            write(2,100)
            return
          endif
c
c   Now compute the third light.  Define a series of mu values
c   from 0.01 to 1.0 and integrate.
c
          Tin=t3
          gin=g3
          itguess=1
          imuguess=1
c
c          do 20 i=1,100
c            rmuin=dble(i)/100.0d0
c            call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
c     &           atmT,atmg,atmmu,Nmu,
c     &           atmint,Tmax,Tmin,gmax,gmin,outinty,
c     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess)
cc       
c            do 19 jj=1,8
c              third(jj)=third(jj)+outinty(jj)*rmuin*0.01d0
c 19         continue
c 20       continue
c
          rmuin=1.0d0
          call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &           atmint7,atmint8,
     &           Tmax,Tmin,gmax,gmin,outinty,
     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     &           dwavex,dwavey,ilaw,iatm,istar)
c
c   UPDATE April 3, 2002
c
c   Scale the fluxes by (separ*solarrad)**2
c
          solarrad=6.9598d10
          DO 30 i=1,8
            third(i)=outinty(i)*SA1*SA3/4.0d0*(separ*solarrad)**2
c            third(i)=third(i)*SA1*SA3/4.0d0*(separ*solarrad)**2
 30       continue
c
 100      format(/'Info:  There is no third light')
c
          return
          end
c
c
c
c  *************************************************************************
c
          subroutine hiddiskgrid(Nrmax,Nthetamax,Nradius,
     %            Ntheta,diskproj,edgeproj,dvisib,evisib,dtemp,tedge,drad,
     $     dx,dy,dz,xxedge,yyedge,zzedge,dinty,einty,stepr,stepz,phase,finc,Q,
     $      Nhoriz,xhoriz,yhoriz,extension,separation,flux,bdist)
c
c   January 14, 2000
c
c   This routine will output files used for various external plotting
c   packages.  This routine is for the disk.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xhoriz(Nhoriz),yhoriz(Nhoriz)
          dimension diskproj(Nrmax*Nthetamax),edgeproj(Nthetamax*11),
     &      dvisib(Nrmax*Nthetamax),evisib(Nthetamax*11),
     &      dinty(Nrmax*Nthetamax),einty(Nthetamax*11),
     &      dtemp(Nrmax*Nthetamax),drad(Nrmax),dy(Nrmax*Nthetamax),
     &      tedge(Nthetamax*11),dx(Nrmax*Nthetamax),dz(Nrmax*Nthetamax),
     %      xxedge(Nthetamax*11),yyedge(Nthetamax*11),
     &      zzedge(Nthetamax*11)
c
c
          character*9 extension
c
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c          
          open(unit=40,file='diskinty.'//extension,status='unknown')
c
c   Check to see if star 1 is in front.  
c
          istar=1
          infront=0
          if((istar.eq.1).and.((phase.ge.0.0d0).
     #        and.(phase.lt.90.0d0)))infront=1
          if((istar.eq.1).and.((phase.ge.270.0d0).
     #        and.(phase.le.360.0d0)))infront=1
          if((istar.eq.2).and.((phase.ge.0.0).
     #        and.(phase.lt.90.0d0)))infront=1
          if((istar.eq.2).and.((phase.ge.270.0d0).
     $        and.(phase.le.360.0d0)))infront=1
c
c   Find the sky coordinates of the center of mass of the disk.  This
c   will be recorded as the first line of the diskinty.???.?? file.
c
          xx=1.0d0
          yy=0.0d0
          zz=0.0d0
          xp=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist)    
          yp=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c          if(phase.gt.180.0d0)then
c            yp=-yp
c          endif
c
          write(40,68)xp*separation,yp*separation,flux
c
          DO 501 ir=1,Nradius-1
            DO 502 ithet=1,Ntheta
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(ir-1)*Ntheta+ithet
c
              iv1=1
              IF (diskproj(iidx).le.0.0d0) go to 502 ! is the surface 
              xx=dx(iidx)                   ! element visible?
              yy=dy(iidx)
              zz=dz(iidx)
              xp=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist) 
              yp=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c              if(phase.gt.180.0d0)then
c                yp=-yp
c              endif
c
              if(infront.eq.1)then    
                iyes=-100
                iv1=1
                call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                if((iyes.eq.100))iv1=0    ! point could be visible
              endif                       ! but is eclipsed
c       
c   Record the x,y,z coordinates of the nearby points.  These points
c   will be used for area filling
c
              xx1=xp
              yy1=yp
c
              if(ithet.gt.1)then
                izz=ir
                jzz=ithet-1
                iidx=(izz-1)*Ntheta+jzz
                xx=dx(iidx)
                yy=dy(iidx)
                zz=dz(iidx)
              else
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                xx=dx(iidx)
                yy=dy(iidx)
                zz=dz(iidx)
              endif
              xx2=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy2=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c              if(phase.gt.180.0d0)then
c                yy2=-yy2
c              endif
c
              if(ithet.gt.1)then
                izz=ir+1
                jzz=ithet-1
                iidx=(izz-1)*Ntheta+jzz
                xx=dx(iidx)
                yy=dy(iidx)
                zz=dz(iidx)
              else
                izz=ir+1
                jzz=Ntheta
                iidx=(izz-1)*Ntheta+jzz
                xx=dx(iidx)
                yy=dy(iidx)
                zz=dz(iidx)
              endif
              xx3=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy3=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c              if(phase.gt.180.0d0)then
c                yy3=-yy3
c              endif
c
              izz=ir+1
              jzz=ithet
              iidx=(izz-1)*Ntheta+jzz
              xx=dx(iidx)                   
              yy=dy(iidx)
              zz=dz(iidx)
              xx4=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy4=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c              if(phase.gt.180.0d0)then
c                yy4=-yy4
c              endif
c
c   Check the corners for eclipsed points.
c
              iv2=1
              iv3=1
              iv4=1
              if(infront.eq.1)then
                iyes=-100
                iv2=1
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx2,yy2,iyes,icut)
                if((iyes.eq.100))iv2=0      
                iv3=1
                iyes=-100
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx3,yy3,iyes,icut)
                if(iyes.eq.100)iv3=0      
                iv4=1
                iyes=-100
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx4,yy4,iyes,icut)
                if(iyes.eq.100)iv4=0      
              endif
c
c   There are 13 possibilities for which corners were hidden.  Do each
c   case separately.
c
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c              
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     %            xx4new*separation,yy4new*separation,  
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3*separation,yy3*separation
                go to 99
              endif
c              
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=5
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.1))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=4
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c 
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                ncorner=5
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     $            xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx3*separation,yy3*separation,
     $            xx4*separation,yy4*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=5
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                ncorner=5
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4new*separation,yy4new*separation,
     &            xx1new*separation,yy1new*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                ncorner=4
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4*separation,yy4*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=4
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
c
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=4
                izz=ir
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)dinty(iidx),ncorner,diskproj(iidx),
     %            dtemp(iidx),
     %            ir,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 99
              endif
 99           continue
502         CONTINUE
501       CONTINUE               
c
          DO 5001 iz=1,10
            DO 5002 ithet=1,Ntheta
              iidx=(iz-1)*Ntheta+ithet
              iv1=1
              IF (edgeproj(iidx).le.0.0d0) go to 5002 ! is the surface 
              xx=xxedge(iidx)                   ! element visible?
              yy=yyedge(iidx)
              zz=zzedge(iidx)
              xp=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist)   
              yp=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c              
c              if(phase.gt.180.0d0)then
c                yp=-yp
c              endif
c
              if(infront.eq.1)then    
                iyes=-100
                iv1=1
                call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                if((iyes.eq.100))iv1=0    ! point could be visible
              endif                       ! but is eclipsed
c       
c   Record the x,y,z coordinates of the nearby points.  These points
c   will be used for area filling
c
              xx1=xp
              yy1=yp
c
              if(ithet.gt.1)then
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                xx=xxedge(iidx)
                yy=yyedge(iidx)
                zz=zzedge(iidx)
              else
                izz=iz
                jzz=Ntheta
                iidx=(izz-1)*Ntheta+jzz
                xx=xxedge(iidx)
                yy=yyedge(iidx)
                zz=zzedge(iidx)
              endif
              xx2=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy2=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c              if(phase.gt.180.0d0)then
c                yy2=-yy2
c              endif
c
              if(ithet.gt.1)then
                izz=iz+1
                jzz=ithet-1
                iidx=(izz-1)*Ntheta+jzz
                xx=xxedge(iidx)
                yy=yyedge(iidx)
                zz=zzedge(iidx)
              else
                izz=iz+1
                jzz=Ntheta
                iidx=(izz-1)*Ntheta+jzz
                xx=xxedge(iidx)
                yy=yyedge(iidx)
                zz=zzedge(iidx)
              endif
              xx3=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy3=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c              if(phase.gt.180.0d0)then
c                yy3=-yy3
c              endif
c
              izz=iz+1
              jzz=ithet
              iidx=(izz-1)*Ntheta+jzz
              xx=xxedge(iidx)                   
              yy=yyedge(iidx)
              zz=zzedge(iidx)
              xx4=diskxtran(xx,yy,zz,phase,fincr,Q,istar,bdist)
              yy4=diskytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c              if(phase.gt.180.0d0)then
c                yy4=-yy4
c              endif
c
c   Check the corners for eclipsed points.
c
              iv2=1
              iv3=1
              iv4=1
              if(infront.eq.1)then
                iyes=-100
                iv2=1
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx2,yy2,iyes,icut)
                if((iyes.eq.100))iv2=0      
                iv3=1
                iyes=-100
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx3,yy3,iyes,icut)
                if(iyes.eq.100)iv3=0      
                iv4=1
                iyes=-100
                call insidecircle(Nhoriz,xhoriz,yhoriz,xx4,yy4,iyes,icut)
                if(iyes.eq.100)iv4=0      
              endif
c
c   There are 13 possibilities for which corners were hidden.  Do each
c   case separately.
c
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation
                go to 999
              endif
c              
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     %            xx4new*separation,yy4new*separation,  
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3*separation,yy3*separation
                go to 999
              endif
c              
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=5
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation,
     %            xx4new*separation,yy4new*separation
                go to 999
              endif
c
              if((iv1.eq.0).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.1))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4*separation,yy4*separation
                go to 999
              endif
c
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4new*separation,yy4new*separation
                go to 999
              endif
c
              if((iv1.eq.0).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.0))then
                ncorner=4
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx1,yy1,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     $            xx1new*separation,yy1new*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 999
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.1))then
                ncorner=4
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4*separation,yy4*separation
                go to 999
              endif
c 
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.1).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx2,yy2,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                ncorner=5
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     $            xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx3*separation,yy3*separation,
     $            xx4*separation,yy4*separation
                go to 999
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=5
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation,
     %            xx4*separation,yy4*separation
                go to 999
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.1).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx3,yy3,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx1new=xedge
                yy1new=yedge
                ncorner=5
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3*separation,yy3*separation,
     %            xx4new*separation,yy4new*separation,
     &            xx1new*separation,yy1new*separation
                go to 999
              endif
c
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.1))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx4,yy4,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                ncorner=4
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4*separation,yy4*separation
                go to 999
              endif
c
              if((iv1.eq.1).and.(iv2.eq.1).and.(iv3.eq.0).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx2,yy2,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=4
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2*separation,yy2*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 999
              endif
c
              if((iv1.eq.1).and.(iv2.eq.0).and.(iv3.eq.0).and.(iv4.eq.0))then  
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx2,yy2,
     #              xedge,yedge)
                xx2new=xedge
                yy2new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx3,yy3,
     #              xedge,yedge)
                xx3new=xedge
                yy3new=yedge
                call clip(Nhoriz,xhoriz,yhoriz,xx1,yy1,xx4,yy4,
     #              xedge,yedge)
                xx4new=xedge
                yy4new=yedge
                ncorner=4
                izz=iz
                jzz=ithet
                iidx=(izz-1)*Ntheta+jzz
                write(40,69)einty(iidx),ncorner,edgeproj(iidx),
     %            tedge(iidx),
     %            iz,ithet,iv1,iv2,iv3,iv4,
     %             xx1*separation,yy1*separation,
     %            xx2new*separation,yy2new*separation,
     %            xx3new*separation,yy3new*separation,
     %            xx4new*separation,yy4new*separation
                go to 999
              endif
c
 999          continue
 5002       continue
 5001     continue
 68       format(2(f9.4,3x),e16.9)
 69       format(e16.9,3x,i3,4x,f6.4,2x,e12.6,1x,1x,
     %       2(i3,1x),6x,4(i1,1x)/,10(f9.4,1x))
c
          close(40)
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine hunt(xx,n,x,jlo)
c
c   This routine was taken from Numerical Recipes, second edition.
c
          implicit double precision(a-h,o-z)

          integer jlo,n
c          real*8  x,xx(n)
          dimension xx(n)
          integer inc,jhi,jm
          logical*2 ascnd
          ascnd=xx(n).gt.xx(1)
          if(jlo.le.0.or.jlo.gt.n)then
            jlo=0
            jhi=n+1
            go to 3
          endif
          inc=1
          if((x.ge.xx(jlo)).eqv.ascnd)then
 1          jhi=jlo+inc
            if(jhi.gt.n)then
              jhi=n+1
            else if((x.ge.xx(jhi)).eqv.ascnd)then
              jlo=jhi
              inc=inc+inc
              go to 1
            endif
          else
            jhi=jlo
 2          jlo=jhi-inc
            if(jlo.lt.1)then
              jlo=0
            else if((x.lt.xx(jlo)).eqv.ascnd)then
              jhi=jlo
              inc=inc+inc
              go to 2
            endif
          endif
 3        if(jhi-jlo.eq.1)return
          jm=(jhi+jlo)/2
          if((x.gt.xx(jm)).eqv.ascnd)then
            jlo=jm
          else
            jhi=jm
          endif
          go to 3
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine recordparm1(fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,
     %       rLx,separ,gamma,t3,g3,SA3,density)
c
c    Will record the parameters used in the file ELC.parm, and will also
c    record interesting computed parameters.
c
          implicit double precision(a-h,o-z)
c
          open(unit=3,file='ELC.parm',status='unknown')
c
          write(3,1004)fill1
          write(3,1005)fill2
          write(3,1006)omega1
          write(3,1007)omega2
          write(3,1009)Q
          write(3,1010)finc
          write(3,1011)Teff1
          write(3,1012)Teff2
          write(3,1015)betarim
          write(3,1016)rinner
          write(3,1017)router
          write(3,1018)tdisk
          write(3,2018)xi
          write(3,1021)rLx
          write(3,1025)separ
          write(3,4025)gamma
          write(3,5000)t3
          write(3,5001)g3
          write(3,5002)SA3
          write(3,5003)density
c
 1000     format(i5,20x,'Nalph1')
 1001     format(i5,20x,'Nbet1')
 1002     format(i5,20x,'Nalph2')
 1003     format(i5,20x,'Nbet2')
 1004     format(f15.13,10x,'fill1')
 1005     format(f15.13,10x,'fill2')
 1006     format(f15.9,10x,'omega1')
 1007     format(f15.9,10x,'omega2')
 1008     format(f15.12,10x,'dphase')
 1009     format(f18.15,7x,'Q')
 1010     format(f15.12,10x,'finc')
 1011     format(f11.4,14x,'Teff1')
 1012     format(f11.4,14x,'Teff2')
 1013     format(f9.7,16x,'Tgrav1')
 1014     format(f9.7,16x,'Tgrav2')
 1015     format(f13.10,12x,'betarim')
 1016     format(f13.10,12x,'rinner')
 1017     format(f13.10,12x,'router')
 1018     format(f12.3,13x,'tdisk')
 2018     format(f12.5,13x,'xi')
 1019     format(i8,17x,'Ntheta')
 1020     format(i8,17x,'Nradius')
c
c  UPDATE March 26, 2002
c
c  The meaning of rLx is now log10(Lx)
c
 1021     format(f15.10,10x,'log10(Lx)')
 1022     format(f9.7,16x,'W')
 1023     format(f20.12,5x,'Period')
 1024     format(f20.13,5x,'fm')
 1025     format(f23.10,2x,'separ')
 1026     format(i2,23x,'idraw')
 1027     format(i2,23x,'iecheck')
 1028     format(i2,23x,'idint')
 1029     format(i2,23x,'ivelout')
 1030     format(i2,23x,'iXout')
 1040     format(f11.9,14x,'darkbol1')
 1041     format(f11.9,14x,'darkbol2')
 1042     format(f11.9,14x,'alb1')
 1043     format(f11.9,14x,'alb2')
 2043     format(i2,23x,'Nref')
 2000     format(f7.1,3x,8(f6.3,2x))
 3028     format(i2,23x,'ilaw  (1=linear law, 2=logarithmic law,',
     %           ' 3=square root law)')
 4000     format(i2,23x,'iatm')
 4001     format(i2,23x,'ism1')
 4002     format(8(i1,1x),4x,'icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK')
 4025     format(f12.6,13x,'gamma velocity')
 5000     format(f15.4,10x,'t3')
 5001     format(f15.5,10x,'g3')
 5002     format(f17.8,8x,'SA3')
 5003     format(f17.6,8x,'density in g/cc')
 5004     format(f12.6,8x,'onephase')
 5005     format(f12.6,8x,'sw2 (currently inactive)')
 5006     format(f12.6,8x,'sw3 (currently inactive)')
 5007     format(f12.6,8x,'T0')
 5008     format(i1,19x,'iRVfilt')
 5009     format(i1,19x,'ionephase')
 5010     format(i1,19x,'isw2 (currently inactive)')
 5011     format(i1,19x,'isw3 (currently inactive)')
 5012     format(i1,19x,'isw4 (currently inactive)')
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine parms1(Teff2,Q,
     %        finc,separ,period,reff1,reff2,rpole1,rpole2,fill1,fill2,
     $        gpole1,gpole2,vrot1,vrot2,gscale1,gscale2,omega1,omega2,
     $        obsparm,bdist,ecc,argper,teff1)
c
c   April 26, 2000
c
c   This routine will compute the component masses, radii, etc. based
c   on the mass ratio, inclination, orbital separation, and orbital period.
c   These values will be written to ELC.parm (unit=3)
c
c
c   July 24, 2000
c
c   Add the array obsparm to the argument list.  This contains the 
c   computed physical parameters of stars 1 and 2
c
c   UPDATE September 21, 2008
c
c   Add K_1 and K_2 to the list
c
c   obsparm(1) = mass of star 1      (solar)
c   obsparm(2) = radius of star 1    (solar)
c   obsparm(3) = gravity of star 1   (log cgs)
c   obsparm(4) = V_rot*sin(i) star 1 (km/sec)
c   obsparm(5) = mass of star 2      (solar)
c   obsparm(6) = radius of star 2    (solar)
c   obsparm(7) = gravity of star 2   (log cgs)
c   obsparm(8) = V_rot*sin(i) star 2 (km/sec)
c   obsparm(9) = duration of X-ray eclipse (degrees)
c   obsparm(10) = K_1 (km/sec)
c   obsparm(11) = K_2 (km/sec)
c
c   UPDATE September 11, 2001
c
c   Change the dimemsion of obsparm to 9.  The X-ray eclipse duration
c   (in degrees) will be stored there
c
c   UPDATE September 21, 2008
c
c   Make the dimension of obsparm 11
c
c   UPDATE October 10, 2008
c
c   Add finc, Q, ecc, argper, teff1, teff2 to the list in obsparm
c 
          implicit double precision(a-h,o-z)
          dimension obsparm(17)
c
          parameter(pie=3.14159265358979323d0)
c
           fpsq=4.0d0*pie*pie
           fincr=finc*3.141592653589879d0/180.0d0        ! radians
           ppp=period*24.0d0               ! period in hours
           coef=3.518847d10
           coef3=4.35713636d31
           sifinc=dsin(fincr)
c
c   Use the formula separ = coef*(perid*period*total_mass)**(1/3) to
c   solve for the total mass in solar masses.  The separation is
c   entered in solar masses, so (R_sun/coef)**3=7.737294491.
c
           total_mass=(separ)**(3)*7.737294491d0/(ppp*ppp)
c
           gmsun=1.32712440018d20   !m^3/sec^2
           smet=separ*6.9598d8
           total_mass=smet*smet*smet*fpsq/(period*86400.0d0)**2/gmsun
           solarrad=6.9598d10
c
           rM1=total_mass/(1.0d0+Q)
           rM2=Q*rM1           
c
           R1=reff1*separ
           if(Teff2.gt.0.0d0)then
             R2=reff2*separ
           else
             R2=0.00d0
           endif
c
           gsun=2.739910d4
           gpole1=gsun*rM1/(R1*R1)  
           if(Teff2.gt.0.0d0)then
             gpole2=gsun*rM2/(R2*R2)  
           else
             gpole2=1.0d0
           endif
c
           gscale1=27397.726d0*rM1/(separ*separ*bdist*bdist)
           gscale2=27397.726d0*rM2/(separ*separ*bdist*bdist)
           gscale1=27397.726d0*rM1/(separ*separ)
           gscale2=27397.726d0*rM2/(separ*separ)
c
           fact=solarrad*2.0d0*pie/86400.0d0/1.0d5  !50.613093d0
           vrot1=fact*R1/period*sifinc
           vrot2=fact*R2/period*sifinc
c
c   UPDATE MARCH 4, 2005
c
c   change a2 and a1 below
c
           a2=separ/(1.0d0+Q)*bdist
c
c   UPDATE September 12, 2001
c 
c   Bug fix, change a1=(separ-a2)*bdist  to  a1=separ(1.0d0-bdist/(1.0d0+Q))
c

           a1=separ*(1.0d0-bdist/(1.0d0+Q))
c
           a2=separ/(1.0d0+Q)
           a1=separ-a2
           efact=1.0d0/dsqrt(1.0-ecc*ecc)
           velK1=fact*a1/period*sifinc*efact
           velK2=fact*a2/period*sifinc*efact
c
c           if(iunit.ge.1)write(2,100)rM1,R1,dlog10(gpole1),rM2,R2,
c     $        dlog10(gpole2),
c     $        period,a1,a2,separ,velK1,velK2,vrot1,vrot2,
c     &        omega1*vrot1,omega2*vrot2
c           if(iunit.le.0)write(*,100)rM1,R1,dlog10(gpole1),rM2,R2,
c     $        dlog10(gpole2),
c     $        period,a1,a2,separ,velK1,velK2,vrot1,vrot2,
c     &        omega1*vrot1,omega2*vrot2
c

c
c   NEW BUG August 10, 2001
c
c   Add a correction factor to the rotational velocities in the
c   case of eccentric orbits.
c
           hutfac=(1.0d0+7.5d0*ecc*ecc+5.625d0*ecc**4+
     #          0.3125d0*ecc**6)/((1.0d0+3.0d0*ecc*ecc+
     $          3.0d0/8.0d0*ecc**4)*dsqrt((1.0d0-ecc*ecc)**3))

           obsparm(1)=rM1
           obsparm(5)=rM2
           obsparm(2)=R1
           obsparm(6)=R2
           obsparm(3)=dlog10(gpole1)
           obsparm(7)=dlog10(gpole2)
           obsparm(4)=vrot1*omega1*hutfac
           obsparm(8)=vrot2*omega2*hutfac
c
c   UPDATE September 21, 2008
c
c   Add the K-velocities to obsparm
c
           obsparm(10)=velK1
           obsparm(11)=velK2
c
           obsparm(12)=finc
           obsparm(13)=Q
           obsparm(14)=ecc
           obsparm(15)=argper
           obsparm(16)=teff1
           obsparm(17)=teff2

c
c
c   UPDATE September 11, 2001
c
c   initialize obsparm(9)
c
c           obsparm(9)=0.0d0
c
           write(3,1000)rM1
           write(3,2000)rM2
           write(3,1001)R1
           write(3,2001)R2
           write(3,1002)dlog10(gpole1)
           write(3,2002)dlog10(gpole2)
           write(3,1003)a1
           write(3,1004)a2
           write(3,1005)separ
           write(3,1006)velK1
           write(3,2006)velK2
           write(3,1007)vrot1
           write(3,2007)vrot2
           write(3,1008)vrot1*omega1*hutfac
           write(3,2008)vrot2*omega2*hutfac
c
 1000      format( f15.10,10x,'mass of star 1 in solar masses')
 2000      format( f15.10,10x,'mass of star 2 in solar masses')
 1001      format(f15.10,10x,'radius of star 1 in solar radii')
 2001      format(f15.10,10x,'radius of star 2 in solar radii')
 1002      format( f9.6,16x,'log(g) of star 1, cgs')
 2002      format( f9.6,16x,'log(g) of star 2, cgs')
 1003      format(f11.5, 14x,'a1 (solar radii)')
 1004      format(f11.5, 14x,'a2 (solar radii)')
 1005      format(f15.10, 10x,'a (solar radii)')
 1006      format( f15.10,10x,'K_1 (km/sec)')
 2006      format( f15.10,10x,'K_2 (km/sec)')
 1007      format( f9.4,16x,'V1_rot*sin(i) (km/sec)')
 2007      format( f9.4,16x,'V2_rot*sin(i) (km/sec)')
 1008      format( f9.4,16x,'V1_rot*sin(i), scaled by omega1 (km/sec)')
 2008      format( f9.4,16x,'V2_rot*sin(i), scaled by omega2 (km/sec)')
c
           close(3)
c
           return
           end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine copyinty(ialphmax1,ibetmax1,Nalph1,Nbet1,
     $              Nalph2,Nbet2,rinty1,saveinty1,rinty2,saveinty2,
     &            mmdx1,mmdx2,ibetlim1,ibetlim2,ialphmax2,ibetmax2)
c
c   May 1, 2000
c
c   This routine is needed to fix a bug in the output routine---the
c   incorrect intensities for the BB mode were being written.
c
          implicit double precision(a-h,o-z)
c
          dimension rinty1(ialphmax1*ibetmax1),
     $        rinty2(ialphmax2*ibetmax2),
     $        saveinty1(ialphmax1*ibetmax1),saveinty2(ialphmax2*ibetmax2)
          dimension mmdx1(ialphmax1,ibetmax1),mmdx2(ialphmax2,ibetmax2)
          dimension ibetlim1(ialphmax1),ibetlim2(ialphmax2)
c
          do 10 ialf=1,Nalph1
            do 9 ibet=1,ibetlim1(ialf)

c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*4*Nbet1+ibet
c
              iidx=mmdx1(ialf,ibet)
              saveinty1(iidx)=rinty1(iidx)
 9          continue
 10       continue
c
          do 20 ialf=1,Nalph2
            do 19 ibet=1,ibetlim2(ialf)
c              iidx=(ialf-1)*4*Nbet2+ibet
              iidx=mmdx2(ialf,ibet)
              saveinty2(iidx)=rinty2(iidx)
 19         continue
 20       continue
c
          return
          end
c
c    &&&&&&&******************%%%%%%%%%%%%%
c
          subroutine copydiskinty(Nrmax,Nthetamax,Nradius,
     %              Ntheta,dinty,savedinty,einty,saveeinty)
c
c
c
          implicit double precision(a-h,o-z)
c
          dimension dinty(Nrmax*Nthetamax),einty(Nthetamax*11),
     &      savedinty(Nrmax*Nthetamax),saveeinty(Nthetamax*11)
c
          do 10 ir=1,Nradius
            do 9 ithet=1,Ntheta
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(ir-1)*Ntheta+ithet
c
              savedinty(iidx)=dinty(iidx)
 9          continue
 10       continue
c
          do 20 ithet=1,Ntheta
            do 19 iz=1,11
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(iz-1)*Ntheta+ithet
c
              saveeinty(iidx)=einty(iidx)
 19         continue
 20       continue
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  UPDATE March 22, 2002
c
c  Remove the variable coy from the argument list of spherepot.
c
c
          subroutine spherepot(Q,omega,cox,coz,r,psi,dpsidr,bdist,
     #        tidephi,itide)
c
c   May 3, 2000
c
c   This routine will return the value of the potential and its
c   derivative at r in the direction (cox,coy,coz)
c
          implicit double precision(a-h,o-z)
c
          dimension pn(0:100),pd(0:100)
          parameter(pie=3.141592653589793d0)
c
          if(itide.lt.2)then
            t1=(bdist*bdist-2.0d0*cox*r*bdist+r*r)
            t2=0.5d0*omega*omega*(1.0d0+Q)*(1.0d0-coz*coz)
            psi=1.0d0/r+Q*(1.0d0/dsqrt(t1)-cox*r/(bdist*bdist))+r*r*t2
            dpsidr=-1.0d0/(r*r)+Q*((dsqrt(t1)**3)*(cox*bdist-r)
     %         -cox/(bdist*bdist))
     $         +t2*2.0d0*r
            return
          endif
c
          coy=coz
          tider=pie*tidephi/180.0d0
          psicos=cox*dcos(tider)+coy*dsin(tider)


          call lpn(itide,psicos,pn,pd)
c
          psi=1.0d0/(r)+Q
          dpsidr=-1.0d0/(r*r)
c

          do 10 ii=2,itide
            psi=psi-Q*(r**(ii))*pn(ii)
            dpsidr=dpsidr+dble(ii)*Q*(rr**(ii-1))*pn(ii)
 10       continue
c
 99       return
          end
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine acchor(Q,psi0,omega,xvis,yvis,zvis,rvis,phivis,
     $              xhid,yhid,zhid,rhid,phihid,
     %              ax,ay,az,xacc,yacc,zacc,bdist,tidephi,itide)
c
c  May 5, 2000
c
c  This subroutine will take two points along a latitude row, where the
c  first is visible and the second is hidden, and iterate to find a
c  more accurate horizon.  We can find coz, which is kept constant.  We
c  bisect on phi, and check the visibility at each iteration.  The
c  x,y,z coordinates of the horizon are returned as xacc,yacc,zacc.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
          if(itide.lt.2)then
            tider=0.0d0
          else
            tider=pie*tidephi/180.0d0
          endif
          do 10 i=1,25
            cozvis=zvis/rvis
            cozhid=zhid/rhid
            coz=0.5d0*(cozvis+cozhid)
            diff=100.0d0
            theta=dacos(coz)
            sithet=dsin(theta)
c
c   Find a new phi, which is the average of the two phi values.
c   Then given theta and phi, compute cox,coy,coz, etc.
c
            diffphi=dabs(phivis-phihid)
            if(diffphi.gt.1.0d0)then
              if(phivis.gt.2.0d0)phivis=phivis-2.0d0*pie
              if(phihid.gt.2.0d0)phihid=phihid-2.0d0*pie
            endif
            diff=phivis-phihid
            phinew=0.5d0*(phivis+phihid)
            if(phinew.lt.0.0)phinew=phinew+2.0d0*pie
            rnew=0.5d0*(rvis+rhid)
            cox=dcos(phinew)*sithet
            coy=dsin(phinew)*sithet
c
            call rad(Q,omega,cox,coy,coz,psi0,rnew,x,y,z,1,bdist,
     $        tidephi,itide)
c
c   UPDATE March 26, 2002
c
c   Remove psixx and istar (the 1) from the argument list of fastPOT.
c
            call fastPOT(Q,omega,x,y,z,psi,psix,psiy,psiz,bdist,
     $         cox,coy,tidephi,itide)
            gravity=DSQRT(PSIX**2+PSIY**2+PSIZ**2)
            GX = -PSIX/gravity
            GY = -PSIY/gravity
            GZ = -PSIZ/gravity
c
            proj=ax*gx+ay*gy+az*gz
            if(proj.gt.0.0d0)then
              xvis=x
              yvis=y
              zvis=z
              rvis=rnew
              phivis=phinew
            else
              xhid=x
              yhid=y
              zhid=z
              rhid=rnew
              phihid=phinew
            endif
            diff=phivis-phihid
 10       continue
c          write(*,*)diff,proj
c
 15       xacc=xvis
          yacc=yvis
          zacc=zvis
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   UPDATE March 26, 2002
c
c   Remove psixx and istar from the argument list of fastPOT.
c
          SUBROUTINE fastPOT(Q,omega,x,y,z,psi,psix,psiy,psiz,bdist,
     #      cox,coy,tidephi,itide)
c
c   May 6, 2000
c
c   This routine computes the potential and its derivatives.
c
          implicit double precision(a-h,o-z)
          dimension pn(0:100),pd(0:100)
          parameter(pie=3.141592653589793d0)
          
          if(itide.lt.2)then
            RST = DSQRT( X**2 + Y**2 + Z**2 )     
            RX = DSQRT((X-bdist)**2 + Y**2 + Z**2 ) 
            A = ((1.0d0+Q)/2.0d0) * OMEGA**2
            RST3 = RST*RST*RST
            RX3 = RX*RX*RX
            PSI = 1.0d0/RST + Q/RX - Q*X/bdist/bdist
     #        + A*(X**2 + Y**2) !potential page 45 of
                                                       !Avni's paper

            PSIY = -Y/RST3 - Q*Y/RX3      +2.*A*Y        !partial deriv. wrt y
            PSIZ = -Z/RST3 - Q*Z/RX3                     !partial deriv. wrt z
            PSIX = -X/RST3 - Q*(X-bdist)/RX3 -Q/bdist/bdist 
     $        + 2.0d0*A*X !partial deriv. wrt x
c
            RETURN
          endif
c
c
c
          tider=pie*tidephi/180.0d0
          psicos=cox*dcos(tider)+coy*dsin(tider)

c
          rr=x*x+y*y+z*z
          w1=(dsqrt(rr))**3
          w1=1.0d0/w1
c
          psi=1.0d0/dsqrt(rr)+Q
          psix=-x*w1
          psiy=-y*w1
          psiz=-z*w1

          call lpn(itide,psicos,pn,pd)
          do 10 ii=2,itide
            t1=-Q*rr**(0.5d0*dble(ii))*pn(ii)
            psi=psi+t1
            twx=dble(ii)*x*Q*rr**(0.5d0*dble(ii)-1.0d0)*pn(ii)
            twy=dble(ii)*y*Q*rr**(0.5d0*dble(ii)-1.0d0)*pn(ii)
            twz=dble(ii)*z*Q*rr**(0.5d0*dble(ii)-1.0d0)*pn(ii)
c
            psix=psix+twx
            psiy=psiy+twy
            psiz=psiz+twz
 10       continue
c
 99       return
          END
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c   UPDATE March 26, 2002
c
c   Remove the variable separ from the argument list of accphi.
c
          subroutine accphi(Q,psi0,omega,phase,fincr,istar,
     %              xvis,yvis,zvis,rvis,phivis,
     $              xhid,yhid,zhid,rhid,phihid,Nhoriz,xhoriz,yhoriz,
     %              phiacc,bdist,tidephi,itide)
c
c  May 5, 2000
c
c  This subroutine will take two points along a latitude row, where the
c  first is visible and the second is hidden, and iterate to find the
c  phi value of the horizon crossing.
c
          implicit double precision(a-h,o-z)
c
          dimension xhoriz(Nhoriz),yhoriz(Nhoriz)
c
          parameter(pie=3.14159265358979323d0)
c
c  Take care of the case where phi1 is slightly larger than zero and
c  phi2 is near 2*pi.
c

          if(itide.lt.2)then
            tider=0.0d0
          else
            tider=pie*tidephi/180.0d0
          endif
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
          cozvis=zvis/rvis
          cozhid=zhid/rhid
          coz=0.5d0*(cozvis+cozhid)
          theta=dacos(coz)
          sithet=dsin(theta)
c
          do 10 i=1,20
c
c   Find a new phi, which is the average of the two phi values.
c   Then given theta and phi, compute cox,coy,coz, etc.
c

            diffphi=dabs(phivis-phihid)
            if(diffphi.gt.2.0d0)then
              if(phivis.gt.2.0d0)phivis=phivis-2.0d0*pie
              if(phihid.gt.2.0d0)phihid=phihid-2.0d0*pie
            endif
            phinew=0.5d0*(phivis+phihid)
            if(phinew.lt.0.0d0)phinew=phinew+2.0d0*pie
            rnew=0.5d0*(rvis+rhid)
            cox=dcos(phinew)*sithet
            coy=dsin(phinew)*sithet
c
            call rad(overQ,omega,cox,coy,coz,psi0,rnew,x,y,z,1,bdist,
     $       tidephi,itide)
c
            xpnew=xtran(x,y,z,phase,fincr,Q,istar,bdist)
            ypnew=ytran(x,y,z,phase,fincr,Q,istar,bdist)
            iyes=-100
            call insidecircle(Nhoriz,xhoriz,yhoriz,xpnew,ypnew,iyes,icut)
c
            if(iyes.eq.100)then
              xhid=x
              yhid=y
              zhid=z
              rhid=rnew
              phihid=phinew
            else
              xvis=x
              yvis=y
              zvis=z
              rvis=rnew
              phivis=phinew
            endif
            diff=phivis-phihid
 10       continue
c
 15       phiacc=phivis
c
          return
          end
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c   UPDATE March 26, 2002
c
c   Comment out this routine, since it is not being used.
c
c          subroutine acctheta(Q,psi0,omega,phase,fincr,istar,
c     %              xvis,yvis,zvis,rvis,thetavis,
c     $              xhid,yhid,zhid,rhid,thetahid,Nhoriz,xhoriz,yhoriz,
c     %              thetaacc,separ,phivis,phihid,bdist)
cc
cc  May 12, 2000
cc
cc  This subroutine will take two points along a longitude row (constant ialf),
cc   where the
cc  first is visible and the second is hidden, and iterate to find the
cc  theta value of the horizon crossing.
cc
c          implicit double precision(a-h,o-z)
cc
c          dimension xhoriz(Nhoriz),yhoriz(Nhoriz)
cc
c          pie=3.141592653589793d0
cc
cc  Take care of the case where phi1 is slightly larger than zero and
cc  phi2 is near 2*pi.
cc
c          overQ=Q
c          if(istar.eq.2)overQ=1.0d0/Q
c          phinew=0.5d0*(phivis+phihid)
c          do 10 i=1,14
cc
cc   Find a new phi, which is the average of the two phi values.
cc   Then given theta and phi, compute cox,coy,coz, etc.
cc
c
c            thetanew=0.5d0*(thetavis+thetahid)
c            rnew=0.5d0*(rvis+rhid)
c            cox=dcos(phinew)*dsin(thetanew)
c            coy=dsin(phinew)*dsin(thetanew)
c            coz=dcos(thetanew)
cc
c            call rad(overQ,omega,cox,coy,coz,psi0,rnew,x,y,z,1,bdist)
cc
c            xpnew=xtran(x,y,z,phase,fincr,Q,istar,bdist)
c            ypnew=ytran(x,y,z,phase,fincr,Q,istar,bdist)
c            iyes=-100
c            call insidecircle(Nhoriz,xhoriz,yhoriz,xpnew,ypnew,iyes,icut)
cc
c            if(iyes.eq.100)then
c              xhid=x
c              yhid=y
c              zhid=z
c              rhid=rnew
c              thetahid=thetanew
c            else
c              xvis=x
c              yvis=y
c              zvis=z
c              rvis=rnew
c              thetavis=thetanew
c            endif
c            diff=thetavis-thetahid
c 222        format(i2,2x,3(f9.6,2x))
c 10       continue
cc
c 15       thetaacc=thetavis
cc
cc          write(46,*)separ*xpnew,separ*ypnew
cc
c          return
c          end
c
c  %%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine findfill(istar,overQ,omega,critpsi,x0,usepot,psi0,fill,
     %           bdist,tidephi,itide)
c
c  May 16, 2000.
c
c  This routine will find the filling factor needed to achieve the
c  surface potential usepot.
c
          implicit double precision (a-h,o-z)
c
          crit=critpsi
          if(istar.eq.2)then
            crit=critpsi/overQ+0.5d0*(overQ-1.0d0)/overQ
          endif
          if(usepot.le.crit)then
            fill=1.0
            write(2,100)istar,usepot,crit,istar
            psi0=critpsi
            return
          endif
c
          xsmall=0.5d0*x0
          xbig=x0
          y=0.0d0
          z=0.0d0
c
c   UPDATE November 14, 2009
c
c   Make the modifications needed to have the tidal approximation.  
c   Initialize cox=1.0 and coy=0.0
c
          cox=1.0d0
          coy=0.0d0
c
          psi0=critpsi
          psismall=critpsi
          do 10 i=1,45
c            write(*,200)xsmall,xbig,psismall,usepot,critpsi
            call POTEN(overQ,omega,xsmall,y,z,psismall,
     $          psix,psixx,psiy,psiz,1,bdist,cox,coy,tidephi,itide)
            psmall=psismall
            if(istar.eq.2)then
              psmall=psismall/overQ+0.5d0*(overQ-1.0d0)/overQ
            endif
            if(psmall.gt.usepot)then
              xsmall=0.5*(xbig-xsmall)+xsmall
            else
              delta=0.5*(xbig-xsmall)
              xbig=xsmall
              xsmall=xbig-delta
            endif
 10       continue
c
          psi0=psismall
          fill=xsmall/x0
          write(2,201)istar,fill

 100      format(/'Info:  The value of usepot',i1,' = ',f11.5,' is less ',
     %           /'than the critical potential = ',f11.5,
     $           '. Setting fill',i1,'=1.0')
c
 200      format(2(f12.9,2x),3(f13.9,2x))
 201      format(/'Info:  fill',i1,' has been set to ',f9.7)
          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine getATMint(maxlines,maxmu,Nlines,
     &      atmT,atmg,atmmu,Nmu,
     &      atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &      Tmax,Tmin,gmax,gmin,gscale,darkint,tpole,gpole,
     %      dwavex,dwavey,ilaw,iatm,istar)
c
c  May 16, 2000
c
c  This subroutine will evaluate the integral:
c
c  dint = int^1_0 (I(T,g,mu)*mu*du)
c
c  The light curves are scaled by DINT for compatability with
c  Wilson and Devinney.
c 
c  The parameter gscale is used to convert the gravities in program units
c  into cgs units.  This number is G*M/(a*a).
c
          implicit double precision(a-h,o-z)

          parameter(pie=3.14159265358979323d0)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines),outinty(8),
     #       darkint(8),summ(8),rnorm(8)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
          dimension dwavex(8,2),dwavey(8,2)

c
c   Find the rough place in the atmosphere table.
c
          Tin=tpole
          gin=dlog10(gscale*gpole)
          call locate(atmT,Nlines,Tin,indexT)
          itguess=indexT
          imuguess=1
c
          do 1 i=1,8
            summ(i)=0.0d0
 1        continue
c
          do 4 i=100,1,-1
            rmuin=dble(i-1)/100.0d0+0.005d0
            call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &           atmint7,atmint8,
     &           Tmax,Tmin,gmax,gmin,outinty,
     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     *           dwavex,dwavey,ilaw,iatm,istar)
c
            do 88 k=1,8
              if(i.eq.100)rnorm(k)=outinty(k)
              summ(k)=summ(k)+2.0d0*0.01*outinty(k)*rmuin
 88         continue
 4        continue
c
          DO 10 i=1,8
            rmuin=1.0d0
            call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &           atmint7,atmint8,
     &           Tmax,Tmin,gmax,gmin,outinty,
     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     &           dwavex,dwavey,ilaw,iatm,istar)
c
            darkint(i)=pie*summ(i)/outinty(i)
 10       continue
c
c          write(*,*)darkint(1),darkint(2),darkint(3)
          return
          end

c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine newsortcircle(N,xcir,ycir,Nhoriz,xhoriz,yhoriz,ibetmax,
     #       xhmin,xhmax,yhmin,yhmax)
c
c   May 18, 2000
c
c   This routine will take the x,y points of the star's horizon, sort them
c   in polar coordinates, resample by interpolation, and return new arrays
c   xhoriz,yhoriz.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xcir(N),ycir(N),theta(45000),rr(45000),
     $       dumt1(45000),dumr1(45000),xhoriz(4*ibetmax),yhoriz(4*ibetmax),
     &       dumt2(45000),dumr2(45000)
c
          call getmaxvalues(N,xcir,ycir,xmin,xmax,ymin,ymax)  !get the extreme
c                                                            !values
          xcenter=(xmax+xmin)/2.0d0
          ycenter=(ymax+ymin)/2.0d0
c
c   convert to polar coordinates         
c
          do 10 i=1,N
            xxx=xcir(i)-xcenter
            yyy=ycir(i)-ycenter
            if(xxx.eq.0.0)then
              if(yyy.lt.0.0)then
                theta(i)=3.0d0*pie/2.0d0
                rr(i)=dsqrt(xxx*xxx+yyy*yyy)
                go to 10
              endif
              if(yyy.ge.0.0d0)then
                theta(i)=pie/2.0d0
                rr(i)=dsqrt(xxx*xxx+yyy*yyy)
                go to 10
              endif
            endif
            if((yyy.ge.0.0d0).and.(xxx.ge.0.0d0))then
              theta(i)=(atan2(yyy,xxx))
              rr(i)=dsqrt(xxx*xxx+yyy*yyy)
            endif
            if((yyy.ge.0.0d0).and.(xxx.lt.0.0d0))then
              theta(i)=(atan2(yyy,xxx))
              rr(i)=dsqrt(xxx*xxx+yyy*yyy)
            endif
            if((yyy.lt.0.0d0).and.(xxx.lt.0.0d0))then
              rr(i)=dsqrt(xxx*xxx+yyy*yyy)
              theta(i)=2.0d0*pie+(atan2(yyy,xxx))
            endif
            if((yyy.lt.0.0d0).and.(xxx.ge.0.0d0))then
              rr(i)=dsqrt(xxx*xxx+yyy*yyy)
              theta(i)=2.0d0*pie+atan2(yyy,xxx)
            endif
 10       continue
c
          call sort3(N,theta,rr,ycir)  !sort by theta and swap x and y also
c
c   Now make a dummy array with theta going from -360 to 720 degrees.
c   This will ensure the interpolation is smooth near the boundaries.
c
          icount=0
          do 20 i=1,N
            icount=icount+1
            dumt1(icount)=theta(i)-2.0d0*pie
            dumr1(icount)=rr(i)
 20       continue
c
          do 30 i=1,N
            icount=icount+1
            dumt1(icount)=theta(i)
            dumr1(icount)=rr(i)
 30       continue
c
          do 40 i=1,N
            icount=icount+1
            dumt1(icount)=theta(i)+2.0d0*pie
            dumr1(icount)=rr(i)
 40       continue
c
          Ndum=icount
c
c   Remove possible repeated points.
c
          icount=0
          do 45 i=2,Ndum
            diff=dabs(dumt1(i)-dumt1(i-1))
            if(diff.gt.1.0d-10)then
              icount=icount+1
              dumt2(icount)=dumt1(i-1)
              dumr2(icount)=dumr1(i-1)
            endif
 45       continue
c
          Ndum=icount

          radcon=pie/180.0d0
c
          index=N
          m=3
          Nhoriz=360
c
c   Update May 8, 2006
c
c   Find the minimum and maximum x and y values
c
          xhmin=1.d20
          xhmax=-1.d20
          yhmin=1.d20
          yhmax=-1.d20
c
          do 50 i=1,360
            angle=dble(i)*radcon
            call hunt(dumt2,Ndum,angle,index)
            k=min(max(index-(m-1)/2,1),Ndum+1-m)
            call polint(dumt2(k),dumr2(k),m,angle,qqq,dy)
c
c   Now assign x,y coordinates based on the qqq (radius) and angle.
c
            xhoriz(i)=qqq*dcos(angle)+xcenter
            yhoriz(i)=qqq*dsin(angle)+ycenter
            if(xhoriz(i).lt.xhmin)xhmin=xhoriz(i)
            if(yhoriz(i).lt.yhmin)yhmin=yhoriz(i)
            if(xhoriz(i).gt.xhmax)xhmax=xhoriz(i)
            if(yhoriz(i).gt.yhmax)yhmax=yhoriz(i)
 50       continue

          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   UPDATE June 14, 2002
c
c   Comment out this subroutine for now.
c
c
c          subroutine diskrefl(ialphmax1,ibetmax1,
c     $      Nalph1,Nbet1,ibetlim1,Nalph2,Nbet2,ibetlim2,ratio1,ratio2,
c     $      xarray1,yarray1,zarray1,gradx1,grady1,gradz1,garray1,surf1,
c     $      xarray2,yarray2,zarray2,gradx2,grady2,gradz2,garray2,surf2,
c     $      temp1,temp2,tempold1,tempold2,dbolx,dboly,ilaw,alb1,alb2,teff1,
c     $      teff2,Tgrav1,Tgrav2,rLx,idint,gpole1,gpole2,
c     %      Tpole1,Tpole2,coprat1,coprat2,Nthetamax,Nrmax,Ntheta,Nradius,
c     %      betarim,rinner,router,reff2,Rl2,separ,
c     $      tdisk,xi,dtemp,dx,dy,dz,drad,
c     $      tedge,xedge,yedge,zedge,redge,stepr,stepz,dratio,reper,
c     &      rsper,ialphmax2,ibetmax2)
cc
cc   July 27, 2000
cc
cc    This routine will alter the temperatures on the disk by
cc    means of a modified version of the 'detailed reflection' 
cc    (R. E. Wilson 1990, ApJ, 356, 613). 
cc
c          implicit double precision(a-h,o-z)
cc
c          parameter(pie=3.14159265358979323d0)
c          dimension xarray1(ialphmax1*ibetmax1),yarray1(ialphmax1*ibetmax1),
c     $      zarray1(ialphmax1*ibetmax1),surf1(ialphmax1*ibetmax1),
c     &      gradx1(ialphmax1*ibetmax1),temp1(ialphmax1*ibetmax1),
c     $      grady1(ialphmax1*ibetmax1),gradz1(ialphmax1*ibetmax1),
c     %      garray1(ialphmax1*ibetmax1),ratio1(ialphmax1*ibetmax1),
c     $      tempold1(ialphmax1*ibetmax1),ibetlim1(ialphmax1),
c     $      ibetlim2(ialphmax2),coprat1(ialphmax1*ibetmax1),
c     $      coprat2(ialphmax2*ibetmax2)
c          dimension xarray2(ialphmax2*ibetmax2),yarray2(ialphmax2*ibetmax2),
c     $      zarray2(ialphmax2*ibetmax2),surf2(ialphmax2*ibetmax2),
c     &      gradx2(ialphmax2*ibetmax2),temp2(ialphmax2*ibetmax2),
c     $      grady2(ialphmax2*ibetmax2),gradz2(ialphmax2*ibetmax2),
c     %      garray2(ialphmax2*ibetmax2),ratio2(ialphmax2*ibetmax2),
c     $      tempold2(ialphmax2*ibetmax2)
c          dimension dbolx(8,2),dboly(8,2)
c          dimension dtemp(Nrmax,Nthetamax),dx(Nrmax,Nthetamax),
c     &      dy(Nrmax,Nthetamax),dz(Nrmax,Nthetamax),dratio(Nrmax,Nthetamax),
c     &      drad(Nrmax),xedge(Nthetamax,11),yedge(Nthetamax,11),
c     &      zedge(Nthetamax,11),tedge(Nthetamax,11)
cc
cc   Start with the disk face.  It is assumed that the lower face is exactly
cc   the same as the upper face, but with a negative z-value.  In practice,
cc   however, we never see the bottom face.
cc          
c          radcon=pie/180.0d0
c          diff2max=-1234.
c          iii=1
c          if(iii.eq.0)then
c            redge=router*reff2              ! radius of outer edge in x units
c            rsmall=rinner*Rl2               ! radius of inner edge in x units
c            reper=redge
c            rsper=rsmall
c          else
c            redge=reper
c            rsmall=rsper
c          endif
cc
c          betarad=betarim*radcon       ! radians
c          steptheta=360.0d0/dble(ntheta)
c          stepr=(redge-rsmall)/dble(Nradius-1)
cc
cc   Transform r into zeta.
cc
c          zetain=2.0d0*dsqrt(rsmall)
c          zetaout=2.0d0*dsqrt(redge)
c          stepzeta=(zetaout-zetain)/dble(Nradius-1)
cc
c          theta=0.0d0
c          zeta=zetain
cc
cc   Start with star 1 and compute the flux from star 2.
cc
c          darkbolx1=dbolx(1,1)
c          darkboly1=dboly(1,1)
c          darkbolx2=dbolx(1,2)
c          darkboly2=dboly(1,2)
cc
c          dtheta1=pie/(1.0d0*nalph1)
c          dtheta2=pie/(1.0d0*nalph2)
cc
cc   Define the integrated  limb darkening coefficients.  The equation is
cc
cc   dint=2*pi*int_0^1{mu*(1-x*(1-mu))d(mu)}  for the linear law, etc.
cc
c          dint1=pie*(1.0d0-darkbolx1/3.0d0)
c          dint2=pie*(1.0d0-darkbolx2/3.0d0)
c          if(ilaw.eq.2)then
c            dint1=pie*(1.0d0-darkbolx1/3.0d0+2.0d0*darkboly1/9.0d0)
c            dint2=pie*(1.0d0-darkbolx2/3.0d0+2.0d0*darkboly2/9.0d0)
c          endif
c          if(ilaw.eq.3)then
c            dint1=pie*(1.0d0-darkbolx1/3.0d0-darkboly1/5.0d0)
c            dint2=pie*(1.0d0-darkbolx2/3.0d0-darkboly2/5.0d0)
c          endif
c          if(teff2.gt.0.0d0)then
c            C1=(Tpole2/Tpole1)**(4)*(dint1/dint2)
c            C1=C1*alb1/dint1
c          endif
c          nalf12=nalph1/2
c          nalf22=nalph2/2
c          DIV1 = gpole1    ! gravity at the pole
c          DIV2 = gpole2    ! gravity at the pole
ccc
cc
c 11       continue
cc
c          do 2 ir=1,Nradius
c            do 1 ithet=1,Ntheta
c              dratio(ir,ithet)=1.0d0
c 1          continue
c 2        continue
cc
c          dint2=1.0
c          T4g1=4.0d0*Tgrav1
c          T4g2=4.0d0*Tgrav2
c          sbet=dsin(betarad)
c          cbet=dcos(betarad)
c          do 20 ir=1,Nradius
c            zeta=zetain+dble(ir-1)*stepzeta
c            r=0.25d0*zeta*zeta
c            do 19 ithet=1,Ntheta ! theta goes from zero to 360-step
c              theta=dble(ithet)*steptheta-0.5*steptheta  ! degrees
c              thetar=theta*radcon            ! radians
c              cthet=dcos(thetar)
c              sthet=dsin(thetar)
c              summ=0.0d0
c              C2=(Tpole1/dtemp(ir,ithet))**(4)*(dint2/dint1)
c              C2=C2*alb2/dint2
cc
c              do 18 i=1,Nalph1/2
c                do 17 j=1,ibetlim1(i)*4
c                  xflip1=xarray1(i,j)     
c                  yflip1=yarray1(i,j)     
c                  zflip1=zarray1(i,j)
c                  dist1=(xflip1-dx(ir,ithet))**2 +
c     %                  (yflip1-dy(ir,ithet))**2 +
c     $                  (zflip1-dz(ir,ithet))**2
c                  dist1=dsqrt(dist1)
c                  term1=(xflip1-dx(ir,ithet))*sbet*cthet+
c     $                  (yflip1-dy(ir,ithet))*(-sthet*sbet)+
c     #                  (zflip1-dz(ir,ithet))*cbet
c                  foreshort1=(term1/(dist1))    
cc
c                  if(foreshort1.le.0.0d0)go to 17
cc
c                  dist2=dist1  
cc                  
c                  xflip2=dx(ir,ithet)
c                  yflip2=dy(ir,ithet)
c                  zflip2=dz(ir,ithet)
cc
c                  term2=(-xflip2+xarray1(i,j))*gradx1(i,j)+
c     $                  (-yflip2+yarray1(i,j))*grady1(i,j)+
c     #                  (-zflip2+zarray1(i,j))*gradz1(i,j)
c                  foreshort2=(term2/(dist2)) 
c                  if(foreshort2.le.0.0d0)go to 17
cc
cc                  term3=(garray1(i,j)/div1)**(T4g1)
cc                  term3=term3*surf1(i,j)*coprat1(i,j)
cc                  term3=term3*foreshort1*foreshort2/(dist2*dist2)
c                  term3=(garray1(i,j)/div1)**(T4g1)*
c     &      surf1(i,j)*coprat1(i,j)*foreshort1*foreshort2/(dist2*dist2)
c                  if(ilaw.eq.1)then
c                    term3=term3*(1.0d0-darkbolx1+darkbolx1*foreshort2)
c                  endif
c                  if(ilaw.eq.2)then
c                    if(foreshort2.gt.0.0d0)then
c                      ttt=darkboly1*foreshort2*dlog(foreshort2)
c                    else
c                      ttt=0.0d0
c                    endif
c                    term3=term3*(1.0d0-darkbolx1*(1.0d0-foreshort2)-ttt)
c                  endif
c                  if(ilaw.eq.3)then
c                    if(foreshort2.gt.0.0d0)then
c                      ttt=darkboly1*(1.0d0-dsqrt(foreshort2))
c                    else
c                      ttt=darkboly1
c                    endif
c                    term3=term3*(1.0d0-darkbolx1*(1.0d0-foreshort2)-ttt)
c                  endif
c                  summ=summ+term3
c 17             continue
c 18           continue
cc
cc   Check to see if the point on star 2 can actually see any point on
cc   the disk.  If not, then the summ will be zero.
cc
c              if(summ.eq.0.0d0)go to 19
cc
c              FpAoverFB=C2*summ
c              dratio(ir,ithet)=1.0d0+FpAoverFB
c              tnew=dtemp(ir,ithet)*dratio(ir,ithet)**0.25d0
c              diff=tnew-dtemp(ir,ithet)
c              if(diff.gt.diff2max)diff2max=diff
c              dtemp(ir,ithet)=tnew
cc              write(*,6969)ir,ithet,theta,dx(ir,ithet),
cc     &           dy(ir,ithet),dz(ir,ithet),dtemp(ir,ithet)
c 19         continue
c 20       continue
c 6969     format(2(i3,1x),f7.3,3x,3(f8.5,2x),f7.1)
cc
cc   Now use symmetry to fill in the other quadrants on the star.  
cc
c 50       continue
cc
c          do 2000 iz=-5,5
c            z=dble(iz)*stepz
c            do 1900 ithet=1,Ntheta        !theta goes from zero to 360-step
c              theta=dble(ithet)*steptheta-0.5*steptheta   ! degrees
c              thetar=theta*radcon         ! radians
cc              xedge(ithet,iz+6)=bdist-redge*dcos(thetar)
cc              yedge(ithet,iz+6)=redge*dsin(thetar)
cc              zedge(ithet,iz+6)=z
c              tedge(ithet,iz+6)=dtemp(Nradius,ithet)
cc
c              if(theta.gt.300.0d0.and.theta.lt.320.0d0)tedge(ithet,iz+6)=
c     %          30.0d0*diff2max+tedge(ithet,iz+6)
cc
c 1900       continue
c 2000     continue
cc
c          if(diff2max.lt.0.0d0)diff2max=0.0
c          if(diff1max.lt.0.0d0)diff1max=0.0
cc
c          write(2,72)diff2max
c 999      continue
cc          write(2,71)diff1max
cc
c 71       format(/'maximum temperature change for star 1 = ',f13.5)
c 72       format('maximum temperature change for the disk = ',f13.5)
c
c          return
c          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine writetempgrav(ialphmax,ibetmax,Nalph,Nbet,ibetlim,
     %        tarray,garray,gscale,istar,xx,yy,zz,mmdx,phistart,separ)
c
c   October 6, 2000
c
c   This subroutine will output the values of the temperature and
c   gravity in physical units [i.e. log(g) in cgs] for each array element.
c   It will also output the angles for use in simple contouring programs.
c
c
c   UPDATE June 7, 2002
c
c   Add the x, y, and z coordinate arrays to the argument of writetempgrab
c
          implicit double precision (a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension garray(ialphmax*ibetmax),tarray(ialphmax*ibetmax),
     %       ibetlim(ialphmax),xx(ialphmax*ibetmax),
     &       yy(ialphmax*ibetmax),zz(ialphmax*ibetmax),
     #       mmdx(ialphmax,ibetmax),phistart(ialphmax)
c
          dtheta=pie/dble(Nalph)
          dcostheta=2.0d0/dble(Nalph)
c
          if(istar.eq.1)then
            open(unit=27,file='star1tempgrav.dat',status='unknown')
          else
            open(unit=27,file='star2tempgrav.dat',status='unknown')
          endif
c
c          isquare=1
          do 10 ialf=1,Nalph
            theta=-0.5d0*dtheta+dtheta*dble(ialf)
c            ibetlim(ialf)=idnint(dsin(theta)*4*Nbet)
c            if(isquare.ge.1)ibetlim(ialf)=4*Nbet
            dphi=2.0d0*pie/dble(ibetlim(ialf))
            DO 9 ibet=1,ibetlim(ialf)          !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              phi=-0.5d0*dphi+dphi*dble(ibet)
              phi=phi+phistart(ialf)
c
c            if(dble(ialf/2).eq.(dble(ialf)*0.5d0))phi=phi+0.25d0*dphi   !DPHI
c
              grav=dlog10(gscale*garray(iidx))
              conphi=phi
              if(phi.gt.pie)conphi=phi-2.0d0*pie
              contheta=theta-0.5d0*pie
              if(ibet.le.2*Nbet)then
                iconbet=ibet+2*Nbet
              else
                iconbet=ibet-2*Nbet
              endif
c
c   UPDATE June 7, 2002
c
c   Add the x, y, and z coordinates to the write statement.  Modify
c   format statement 100 below (add 3(f7.4,1x) to the end).
c
              RRR=dsqrt(xx(iidx)**2+yy(iidx)**2+zz(iidx)**2)
              RRR=RRR*separ
              write(27,100)tarray(iidx),grav,theta,phi,ialf,ibet,
     %         contheta,conphi,ialf,iconbet,
     %         xx(iidx),yy(iidx),zz(iidx),RRR
 9          continue
 10       continue
          close(27)
c
 100      format(f12.5,1x,f8.5,1x,2(f8.6,1x),1x,
     #          2(i3,1x),1x,f8.5,2x,f8.5,2x,2(i3,1x),3(f7.4,1x),1x,f12.6)
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  UPDATE March 26, 2002
c
c  Remove the variables icount and dphase from the argument list of
c  shiftlc (they are not used).
c
          subroutine shiftlc(Nmaxphase,Nphase,xxx,yyy,pshift,
     %       ilast)
c
c   February 5, 2001
c
c   This routine will apply a phase shift to a light or velocity curve

c   
          implicit double precision(a-h,o-z)
c
          parameter(NNdum=720001)
          dimension xxx(Nmaxphase),yyy(Nmaxphase),ydum(NNdum),xdum(NNdum)
c
c          write(*,*)Nphase
          do 5 i=1,Nphase
            xdum(i)=xxx(i)
 5        continue
c
          tiny=1.0d-7
          do 10 i=1,Nphase
            xdum(i)=xdum(i)+pshift
c            if(xdum(i).ge.(1.0d0-tiny))xdum(i)=xdum(i)-1.0d0
c            if(xdum(i).lt.0.0d0)xdum(i)=xdum(i)+1.0d0
            xdum(i)=dmod(xdum(i),1.0d0)
            ydum(i)=yyy(i)
 10       continue
c
c   RVG BUG ALERT  April 19, 2001
c
c   Add the loop below to ensure that the final x-values are always between
c   0.0 and 1.0.  In some cases where the eccentricity is large, the 
c   first pass above might not be enough.
c
          xmax=-999.d0
          xmin=999.d0
 99       do 11 i=1,Nphase
            if(xdum(i).gt.xmax)xmax=xdum(i)
            if(xdum(i).lt.xmin)xmin=xdum(i)
            if(xdum(i).ge.1.0d0)xdum(i)=xdum(i)-1.0d0
            if(xdum(i).lt.0.0d0)xdum(i)=xdum(i)+1.0d0
 11       continue
c
          if(xmax.ge.2.0d0)go to 99
          if(xmin.lt.-1.0d0)go to 99
c
c   END BUG
c
          call sort3(Nphase,xdum,yyy,ydum)
c
          if(ilast.eq.1)then
            do 15 i=1,Nphase
              xxx(i)=xdum(i)
 15         continue
          endif
c
 101      format(f6.2)
          return
          end
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine getE(em,ecc,bigE)
c
c    Solve for bigE in  em = bigE - ecc*sin(bigE)  (i.e. Kepler's Eq.)
c
          implicit double precision (a-h,o-z)
          Eold=em
c
          do 10 i=1,20
            top=Eold-ecc*dsin(Eold)-em
            bottom=1.0d0-ecc*dcos(Eold)
            Enew=Eold-top/bottom
            diff=dabs(Enew-Eold)
            if(diff.lt.1.0d-15)go to 15
            Eold=Enew
 10       continue
c
 15       bigE=Eold
c
          return
c
          end
c
c  *************************************************************************
c
c  UPDATE March 26, 2002
c
c  Get rid of the variables phiar, jj, istar from the argument list
c  of getBBsimp.
c
c  UPDATE June 22, 2002
c
c  Add separ to the argument list.
c
c
          subroutine getBBsimp(ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     $      wave,visib,projarray,temp,surf,flimbx,flimby,ilaw,rinty,
     &      flum,flux,delphi,delphiedge,iedgestar,iedgehor,rldint,
     %      isimp,separ,mmdx)
c
c  October 11, 1999
c
c  This routine will compute the intensities of each element, given the
c  temperatures (temp(ialf,ibet)) and the input wavelength.  It will then
c  integrate the flux given the visibilities (visib) and surface elements
c  (surf).   
c
c  The projarray contains the cosine mu terms for each element.  The visib
c  array contains the cosine mu terms for each element, except if the point
c  is eclipsed in which case the visib=0
c
c  April 2, 2001
c
c  This version will integrate the flux along latitude rows using a
c  Simpson's rule summation. 
c 
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension visib(ialphmax*ibetmax),delphi(ialphmax*ibetmax),
     $        surf(ialphmax*ibetmax),ibetlim(ialphmax),
     $        temp(ialphmax*ibetmax),flum(ialphmax*ibetmax),
     $        rinty(ialphmax*ibetmax),projarray(ialphmax*ibetmax),
     %        iedgehor(ialphmax*ibetmax),iedgestar(ialphmax*ibetmax),
     $        delphiedge(ialphmax*ibetmax),saveflum(1000*1000),
     #        iflag(1000*1000),mmdx(ialphmax,ibetmax)
          dimension yflux(1000),xflux(1000),icase(1000),mark(2,1000),
     #        yend(2),xend(2),yvalue(2)
c
c   UPDATE March 26, 2002.
c
c   Remove the variable flagx1 from the list of the declaration statement.
c
          logical x1assign,x2assign
c
c   UPDATE June 22, 2002
c
c   Declare flag as a logical variable.
c
          logical flag
c
c
c  The array icase(ialf) contains information about the arrangement of
c  the horizons:
c
c  icase=1 : all phi points for a given ialf are visible 
c
c  icase=2 : phi(1)<phi(2)  and all points inbetween are visible
c
c  icase=3 : phi(1)<phi(2)  points between phi(2),2*pi and 0,phi(1)
c            are visible
c
c  The variables yend(1) and yend(2) contain the corrections for the end 
c  points:
c
c      icase=1  yend(1)=0.5*yflux(1)     yend(2)=0.5*y(N)
c      etc.
c
c   UPDATE June 17, 2002
c
c   Add this dummy assignment for separ to supress compiler warning
c   about unused variable.
c
C
C   UPDATE June 11, 2003
c
c   Comment out the routine from here on.
C
          stop
c
c          ddd1=separ
cc
c          do 111 i=1,1000
c            xflux(i)=0.0d0
c            yflux(i)=0.0d0
c            icase(i)=0
c            mark(1,i)=0
c            mark(2,i)=0
c 111      continue
cc
c          dint=pie*(1.0d0-flimbx/3.0d0)
c          if(ilaw.eq.2)then
c            dint=pie*(1.0d0-flimbx/3.0d0+2.0d0*flimby/9.0d0)
c          endif
c          if(ilaw.eq.3)then
c            dint=pie*(1.0d0-flimbx/3.0d0-flimby/5.0d0)
c          endif
cc
c          itotal=0
c          igood=0
c          flux=0.0d0
c          C2 = 1.4384d8          ! 1.4384 * 10.**8      ! hc/(k*1e-8)
c          C1 = 1.191044d35       ! 2hc^2/((1e-8)**5)
cc
cc   Initialize the flum matrix.
cc
cc
c          do 2 ialf=1,nalf
c            do 1 ibet=1,4*Nbet
c              iidx=(ialf-1)*4*Nbet+ibet
c              flum(iidx)=0.0d0
c              rinty(iidx)=0.0d0
c              saveflum(iidx)=0.0d0
c              iflag(iidx)=-99
c 1          continue
c 2        continue
cc
cc   Compute the intensity values of pixels just behind the edge or the
cc   eclipsing horizon for use in fractional eclipse corrections.
cc
c          wavemu=wave/10000.0d0
c          c1=3.74185
c          c2=14.3883
c          do 4 ialf=1,nalf
c            flag=.true.
c            imark=0
c            do 3 ibet=1,ibetlim(ialf)
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c              if(projarray(iidx).lt.0.0d0)flag=.false.
c              if(iedgestar(iidx).ge.10)then
c                imark=imark+1
c                mark(imark,ialf)=ibet
c              endif
c              if((iedgehor(ialf,ibet).eq.-10).or.(iedgehor(ialf,ibet).gt.
c     #             5).or.(iedgestar(ialf,ibet).eq.-10).or.(iedgestar(ialf,ibet)
c     $             .gt.5).or.(delphi(ialf,ibet).gt.-10.0d0))then
c                C3 = C2/(WAVE*TEMP(IALF,IBET))
c
c                tkkelv=temp(ialf,ibet)/1000.0d0
c                C3 = C2/(wavemu*tkkelv)
c                saveflum(ialf,ibet)=C1/(dexp(c3)-1.0d0)/wavemu**5
c                dark=(1.0d0-flimbx+flimbx*dabs(projarray(ialf,ibet)))
c                if(ilaw.eq.2)dark=dark-flimby*dabs(projarray(ialf,ibet))*
c     %                dlog(dabs(projarray(ialf,ibet)))
c                if(ilaw.eq.3)dark=dark-flimby*(1.0d0-
c     &            dsqrt(dabs(projarray(ialf,ibet))))
c                saveflum(ialf,ibet)=saveflum(ialf,ibet)*dark
c                saveflum(ialf,ibet)=surf(ialf,ibet)*saveflum(ialf,ibet)*
c     $            dabs(projarray(ialf,ibet))
c              endif
c 3          continue
c            if(flag.eqv..true.)icase(ialf)=1
c 4        continue
cc
c          do 44 ialf=1,nalf
c            if(icase(ialf).eq.1)go to 44
c            flag=.true.
c            do 33 ibet=mark(1,ialf),mark(2,ialf)
c              if(projarray(ialf,ibet).lt.0.0d0)flag=.false.              
c 33         continue
c            if(flag.eqv..true.)icase(ialf)=2
c 44       continue
cc
c          do 444 ialf=1,nalf
c            if(icase(ialf).eq.1)go to 444
c            if(icase(ialf).eq.2)go to 444
c            flag=.true.
c            do 333 ibet=mark(2,ialf),ibetlim(ialf)
c              if(projarray(ialf,ibet).lt.0.0d0)flag=.false.              
c 333         continue
c            do 334 ibet=1,mark(1,ialf)
c              if(projarray(ialf,ibet).lt.0.0d0)flag=.false.              
c 334        continue
c            if(flag.eqv..true.)icase(ialf)=3
c 444      continue
cc
c          do 555 ialf=1,nalf
c            if(icase(ialf).eq.0)write(*,*)'icase=0 ',icase(ialf)
c 555      continue
cc
c          sumcor1=0.0d0
c          sumcor2=0.0d0
c          totalsum=0.0d0
c          dtheta=0.5d0*pie/dble(Nalf)
cc
c          flag=.false.
c          icount=0
c          phisum=0.0d0
c          DO 10 ialf=1,nalf
c            theta=-dtheta+2.0d0*dtheta*dble(ialf)
c            sitheta=dsin(theta)
cc
cc   Caution!  The dphi defined here is exactly half the dphi in the
cc   subroutine setupgeo!
cc
c            dphi=pie/dble(ibetlim(ialf))
c            if(flag.eqv..true.)then
c              if(isimp.eq.2)then
c                call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                totalsum=totalsum+rowsum  
c              else
c                call simpson(icount,yflux,rowsum,dphi)
c                totalsum=totalsum+rowsum+yend(1)+yend(2)
c              endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
c
cc              if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 1',ialf,ibet
cc              if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 1',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c              itotal=itotal+icount
c              icount=0
c              rowsum=0.0d0
c              yend(1)=0.0d0
c              yend(2)=0.0d0
c              flag=.false.
c            endif
cc
cc   Initialize the flag which indicates that we are going along a
cc   visible row.
cc
c            flag=.false.
c            icount=0
c            rowsum=0.0d0
c            yend(1)=0.0d0
c            yend(2)=0.0d0
c            phisum=0.0d0
c            x1assign=.false.
c            x2assign=.false.
c            if(icase(ialf).eq.1)then
c              DO 9 ibet = 1,ibetlim(ialf)              !4*nbet
c                phiabs=dabs(delphi(ialf,ibet))
c                edgeabs=dabs(delphiedge(ialf,ibet))
c                if((projarray(ialf,ibet).le.0.0d0))then
cc
cc   We have gone over the horizon.  If previous ibet values were visible,
cc   then call the Simpson routine.  Otherwise, go on to the next ibet.
cc
c                  if(flag.eqv..false.)go to 9
c                  if(isimp.eq.2)then
c                    call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                    totalsum=totalsum+rowsum  
c                  else
c                    call simpson(icount,yflux,rowsum,dphi)
c                    totalsum=totalsum+rowsum+yend(1)+yend(2)
c                  endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
c
c                  if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 2',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 2',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c                  itotal=itotal+icount
c                  icount=0
c                  flag=.false.
c                  phisum=0.0d0
c                  go to 9
c                endif
c                tkkelv=temp(ialf,ibet)/1000.0d0
c                C3 = C2/(wavemu*tkkelv)
c                flum(ialf,ibet)=C1/(dexp(c3)-1.0d0)/wavemu**5
c                dark=(1.0d0-flimbx+flimbx*projarray(ialf,ibet))
c                if(ilaw.eq.2)dark=dark-flimby*projarray(ialf,ibet)*
c     %               dlog(projarray(ialf,ibet))
c                if(ilaw.eq.3)dark=dark-flimby*(1.0-dsqrt(projarray(ialf,ibet)))
c                flum(ialf,ibet)=flum(ialf,ibet)*dark
c                rinty(ialf,ibet)=flum(ialf,ibet) ! save intys for plotting
c                saveflum(ialf,ibet)=surf(ialf,ibet)*flum(ialf,ibet)*
c     $            projarray(ialf,ibet)
c                flum(ialf,ibet)=surf(ialf,ibet)
c     #             *flum(ialf,ibet)*visib(ialf,ibet)
cc
c                if((visib(ialf,ibet).le.0.0d0))then
cc
cc   We have a point that would have been visible (proj>0), but it is eclipsed.
cc   Call the Simpson routine using the current string of points.
cc
c                  if(flag.eqv..false.)go to 9
c                  if(isimp.eq.2)then
c                    call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                    totalsum=totalsum+rowsum  
c                  else
c                    call simpson(icount,yflux,rowsum,dphi)
c                    totalsum=totalsum+rowsum+yend(1)+yend(2)
c                  endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
c
c                  if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 3',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 3',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c                  itotal=itotal+icount
c                  icount=0
c                  flag=.false.
c                  rowsum=0.0d0
c                  yend(1)=0.0d0
c                  yend(2)=0.0d0
c                  phisum=0.0d0
c                endif
c                if((visib(ialf,ibet).gt.0.0d0))then
cc
cc   We have a point that is visible.  Add it to the current string of
cc   points.
cc
c                  flag=.true.
c                  icount=icount+1
c                  yflux(icount)=flum(ialf,ibet)
c                  xflux(icount)=phisum
c                  phisum=phisum+2.0d0*dphi
cc
cc   In case 1 xend(1) might not necessairly be assigned since it may not
cc   be next to an eclipsing horizon or the edge horizon.
cc
c                  xend(1)=-2.0d0*dphi              
c                  yend(1)=0.5d0*yflux(1)/dabs(xend(1))
c                  x1assign=.true.
c                  if(icount.eq.1)then
c                    if(iedgestar(ialf,ibet).ge.10)then 
cc
cc   If this point is the first in the string and at the horizon, then
cc   the connecting end piece is a triangle.  The area is 0.5 times the
cc   y-value times the phi distance.  
cc
c                      yvalue(1)=0.0d0
c                      yend(1)=0.25d0*yflux(1)/dphi*edgeabs 
c                      xend(1)=xflux(1)-edgeabs
c                      x1assign=.true.
c                    else 
c                      if(iedgehor(ialf,ibet).eq.10)then
cc
cc   If the first point is next to an eclipsing horizon, then the end
cc   piece will be a trapizoid.  Use interpolation to find the y-value
cc   at the horizon.
cc
c                        if(ibet.lt.ibetlim(ialf))then    
c                          y1=saveflum(ialf,ibet+1)
c                          y0=yflux(icount)
c                          slope=0.5d0*phiabs/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y0+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        else
c                          y1=saveflum(ialf,1)
c                          y0=yflux(icount)
c                          slope=0.5d0*phiabs/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y0+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        endif
c                      endif
c                      if(iedgehor(ialf,ibet).eq.20)then
cc
cc   If the first point is next to an eclipsing horizon, then the end
cc   piece will be a trapizoid.  Use interpolation to find the y-value
cc   at the horizon.
cc
c                        if(ibet.gt.1)then    
c                          y0=saveflum(ialf,ibet-1)
c                          y1=yflux(icount)
c                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y1+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        else
c                          y0=saveflum(ialf,ibetlim(ialf))
c                          y1=yflux(icount)
c                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y1+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        endif
c                      endif
c                    endif  ! end if iedgehor > 10
c                  endif ! end if icount=1
c                  if(iedgestar(ialf,ibet).ge.10)then 
cc
cc   The point in the string (not the first) is at the horizon.  The
cc   second end piece is a triangle.
c
c                    yvalue(2)=0.0d0
c                    yend(2)=0.25d0*yflux(icount)/dphi*edgeabs 
c                    xend(2)=xflux(icount)+edgeabs
c                    x2assign=.true.
c                  else
cc
cc   Otherwise, check to see if it near an eclipsing horizon.
cc
c                    if(iedgehor(ialf,ibet).eq.10)then
c                      if(ibet.lt.ibetlim(ialf))then    
c                        y1=saveflum(ialf,ibet+1)
c                        y0=yflux(icount)
c                        slope=0.5d0*phiabs/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y0+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      else
c                        y1=saveflum(ialf,1)
c                        y0=yflux(icount)
c                        slope=0.5d0*phiabs/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yend(2)=0.5d0*(y0+yinterp)*slope
c                        yvalue(2)=yinterp
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      endif
c                    endif
c                    if(iedgehor(ialf,ibet).eq.20)then
c                      if(ibet.gt.1)then    
c                        y0=saveflum(ialf,ibet-1)
c                        y1=yflux(icount)
c                        slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y1+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      else
c                        y0=saveflum(ialf,ibetlim(ialf))
c                        y1=yflux(icount)
c                        slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yend(2)=0.5d0*(y1+yinterp)*slope
c                        yvalue(2)=yinterp
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      endif  ! end if ibet > 1
c                    endif    ! end if iedgehor > 20 
c                  endif      ! end if iedgestar > 10
c                endif        ! end if visib > 0
cc
c                flux=flux+flum(ialf,ibet)   
c                flum(ialf,ibet)=flum(ialf,ibet)
c                igood=igood+1
c 9            continue
c            endif            ! end if icase=1
cc
c            if(icase(ialf).eq.2)then
c              DO 79 ibet = mark(1,ialf),mark(2,ialf)              !4*nbet
c                phiabs=dabs(delphi(ialf,ibet))
c                edgeabs=dabs(delphiedge(ialf,ibet))
c                if((projarray(ialf,ibet).le.0.0d0))then
c                  if(flag.eqv..false.)go to 79
c                  if(isimp.eq.2)then
c                    call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                    totalsum=totalsum+rowsum  
c                  else
c                    call simpson(icount,yflux,rowsum,dphi)
c                    totalsum=totalsum+rowsum+yend(1)+yend(2)
c                  endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
c
c                  if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 4',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 4',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c                  itotal=itotal+icount
c                  icount=0
c                  flag=.false.
c                  phisum=0.0d0
c                  go to 79
c                endif
c                tkkelv=temp(ialf,ibet)/1000.0d0
c                C3 = C2/(wavemu*tkkelv)
c                flum(ialf,ibet)=C1/(dexp(c3)-1.0d0)/wavemu**5
c                dark=(1.0d0-flimbx+flimbx*projarray(ialf,ibet))
c                if(ilaw.eq.2)dark=dark-flimby*projarray(ialf,ibet)*
c     %               dlog(projarray(ialf,ibet))
c                if(ilaw.eq.3)dark=dark-flimby*(1.0-dsqrt(projarray(ialf,ibet)))
c                flum(ialf,ibet)=flum(ialf,ibet)*dark
c                rinty(ialf,ibet)=flum(ialf,ibet) ! save intys for plotting
c                saveflum(ialf,ibet)=surf(ialf,ibet)*flum(ialf,ibet)*
c     $            projarray(ialf,ibet)
c                flum(ialf,ibet)=surf(ialf,ibet)
c     #             *flum(ialf,ibet)*visib(ialf,ibet)
cc
c                if((visib(ialf,ibet).le.0.0d0))then
c                  if(flag.eqv..false.)go to 79
c                  if(isimp.eq.2)then
c                    call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                    totalsum=totalsum+rowsum  
c                  else
c                    call simpson(icount,yflux,rowsum,dphi)
c                    totalsum=totalsum+rowsum+yend(1)+yend(2)
c                  endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
cc
c
c                  if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 4',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 4',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c                  itotal=itotal+icount
c                  icount=0
c                  flag=.false.
c                  rowsum=0.0d0
c                  yend(1)=0.0d0
c                  yend(2)=0.0d0
c                  phisum=0.0d0
c                endif
c                if((visib(ialf,ibet).gt.0.0d0))then
c                  flag=.true.
c                  icount=icount+1
c                  yflux(icount)=flum(ialf,ibet)
c                  xflux(icount)=phisum
c                  phisum=phisum+2.0d0*dphi
c                  if(icount.eq.1)then
c                    if(iedgestar(ialf,ibet).ge.10)then 
cc
cc   If this point is the first in the string and at the horizon, then
cc   the connecting end piece is a triangle.  The area is 0.5 times the
cc   y-value times the phi distance.  
cc
c                      yvalue(1)=0.0d0
c                      yend(1)=0.25d0*yflux(1)/dphi*edgeabs
c                      xend(1)=xflux(1)-edgeabs 
c                      x1assign=.true.
c                    else
c                      if(iedgehor(ialf,ibet).eq.10)then
cc
cc   If the first point is next to an eclipsing horizon, then the end
cc   piece will be a trapizoid.  Use interpolation to find the y-value
cc   at the horizon.
cc
c                        if(ibet.lt.ibetlim(ialf))then    
c                          y1=saveflum(ialf,ibet+1)
c                          y0=yflux(icount)
c                          slope=0.5d0*phiabs/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y0+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs 
c                      x1assign=.true.
c                        else
c                          y1=saveflum(ialf,1)
c                          y0=yflux(icount)
c                          slope=0.5d0*phiabs/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y0+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        endif
c                      endif
c                      if(iedgehor(ialf,ibet).eq.20)then
cc
cc   If the first point is next to an eclipsing horizon, then the end
cc   piece will be a trapizoid.  Use interpolation to find the y-value
cc   at the horizon.
cc
c                        if(ibet.gt.1)then    
c                          y0=saveflum(ialf,ibet-1)
c                          y1=yflux(icount)
c                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y1+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        else
c                          y0=saveflum(ialf,ibetlim(ialf))
c                          y1=yflux(icount)
c                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y1+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        endif
c                      endif
c                    endif
c                  endif ! end if icount=1
c                  if(iedgestar(ialf,ibet).ge.10)then 
cc
cc   The point in the string (not the first) is at the horizon.  The
cc   second end piece is a triangle.
cc
c                    yvalue(2)=0.0d0
c                    yend(2)=0.25d0*yflux(icount)/dphi*edgeabs
c                    xend(2)=xflux(icount)+edgeabs 
c                    x2assign=.true.
c                  else
cc
cc   Otherwise, check to see if it near an eclipsing horizon.
cc
c                    if(iedgehor(ialf,ibet).eq.10)then
c                      if(ibet.lt.ibetlim(ialf))then    
c                        y1=saveflum(ialf,ibet+1)
c                        y0=yflux(icount)
c                        slope=0.5d0*phiabs/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y0+yinterp)*slope
c                        xend(2)=xflux(icount)+delphi(ialf,ibet) 
c                    x2assign=.true.
c                      else
c                        y1=saveflum(ialf,1)
c                        y0=yflux(icount)
c                        slope=0.5d0*phiabs/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y0+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      endif
c                    endif
c                    if(iedgehor(ialf,ibet).eq.20)then
c                      if(ibet.gt.1)then    
c                        y0=saveflum(ialf,ibet-1)
c                        y1=yflux(icount)
c                        slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y1+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      else
c                        y0=saveflum(ialf,ibetlim(ialf))
c                        y1=yflux(icount)
c                        slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y1+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      endif
c                    endif
c                  endif   
c                endif            ! end if visib > 0
cc
c                flux=flux+flum(ialf,ibet)    
c                flum(ialf,ibet)=flum(ialf,ibet)
c                igood=igood+1
c 79           continue
c            endif
cc
c            if(icase(ialf).eq.3)then
c              DO 89 ibet = mark(2,ialf),ibetlim(ialf)
c                phiabs=dabs(delphi(ialf,ibet))
c                edgeabs=dabs(delphiedge(ialf,ibet))
c                if((projarray(ialf,ibet).le.0.0d0))then
c                  if(flag.eqv..false.)go to 89
c                  if(isimp.eq.2)then
c                    call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                    totalsum=totalsum+rowsum  
c                  else
c                    call simpson(icount,yflux,rowsum,dphi)
c                    totalsum=totalsum+rowsum+yend(1)+yend(2)
c                  endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
cc
c                  if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 5',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 5',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c                  itotal=itotal+icount
c                  icount=0
c                  flag=.false.
c                  phisum=0.0d0
c                  go to 89
c                endif
c                tkkelv=temp(ialf,ibet)/1000.0d0
c                C3 = C2/(wavemu*tkkelv)
c                flum(ialf,ibet)=C1/(dexp(c3)-1.0d0)/wavemu**5
c                dark=(1.0d0-flimbx+flimbx*projarray(ialf,ibet))
c                if(ilaw.eq.2)dark=dark-flimby*projarray(ialf,ibet)*
c     %               dlog(projarray(ialf,ibet))
c                if(ilaw.eq.3)dark=dark-flimby*(1.0-dsqrt(projarray(ialf,ibet)))
c                flum(ialf,ibet)=flum(ialf,ibet)*dark
c                rinty(ialf,ibet)=flum(ialf,ibet) ! save intys for plotting
c                saveflum(ialf,ibet)=surf(ialf,ibet)*flum(ialf,ibet)*
c     $            projarray(ialf,ibet)
c                flum(ialf,ibet)=surf(ialf,ibet)
c     #             *flum(ialf,ibet)*visib(ialf,ibet)
cc
c                if((visib(ialf,ibet).le.0.0d0))then
c                  if(flag.eqv..false.)go to 89
c                  if(isimp.eq.2)then
c                    call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                    totalsum=totalsum+rowsum  
c                  else
c                    call simpson(icount,yflux,rowsum,dphi)
c                    totalsum=totalsum+rowsum+yend(1)+yend(2)
c                  endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                rite(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
cc
c                  if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 6',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 6',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c
c                  itotal=itotal+icount
c                  icount=0
c                  flag=.false.
c                  rowsum=0.0d0
c                  yend(1)=0.0d0
c                  yend(2)=0.0d0
c                  phisum=0.0d0
c                endif
c                if((visib(ialf,ibet).gt.0.0d0))then
c                  flag=.true.
c                  icount=icount+1
c                  yflux(icount)=flum(ialf,ibet)
c                  xflux(icount)=phisum
c                  phisum=phisum+2.0d0*dphi
c                  if(icount.eq.1)then
c                    if(iedgestar(ialf,ibet).ge.10)then 
cc
cc   If this point is the first in the string and at the horizon, then
cc   the connecting end piece is a triangle.  The area is 0.5 times the
cc   y-value times the phi distance.  
cc
c                      yvalue(1)=0.0d0
c                      yend(1)=0.25d0*yflux(1)/dphi*edgeabs
c                      xend(1)=xflux(1)-edgeabs 
c                      x1assign=.true.
c                    else
c                      if(iedgehor(ialf,ibet).eq.10)then
cc
cc   If the first point is next to an eclipsing horizon, then the end
cc   piece will be a trapizoid.  Use interpolation to find the y-value
cc   at the horizon.
cc
c                        if(ibet.lt.ibetlim(ialf))then    
c                          y1=saveflum(ialf,ibet+1)
c                          y0=yflux(icount)
c                          slope=0.5d0*phiabs/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y0+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        else
c                          y1=saveflum(ialf,1)
c                          y0=yflux(icount)
c                          slope=0.5d0*phiabs/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                          yend(1)=0.5d0*(y0+yinterp)*slope
c                        endif
c                      endif
c                      if(iedgehor(ialf,ibet).eq.20)then
cc
cc   If the first point is next to an eclipsing horizon, then the end
cc   piece will be a trapizoid.  Use interpolation to find the y-value
cc   at the horizon.
cc
c                        if(ibet.gt.1)then    
c                          y0=saveflum(ialf,ibet-1)
c                          y1=yflux(icount)
c                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y1+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        else
c                          y0=saveflum(ialf,ibetlim(ialf))
c                          y1=yflux(icount)
c                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                          yend(1)=0.5d0*(y1+yinterp)*slope
c                        endif
c                      endif
c                    endif   ! end if iedgehor > 10
c                  endif ! end if icount=1
c                  if(iedgestar(ialf,ibet).ge.10)then 
cc
cc   The point in the string (not the first) is at the horizon.  The
cc   second end piece is a triangle.
cc
c                    yvalue(2)=0.0d0
c                    yend(2)=0.25d0*yflux(icount)/dphi*edgeabs
c                    xend(2)=xflux(icount)+edgeabs
c                    x2assign=.true.
c                  else
cc
cc   Otherwise, check to see if it near an eclipsing horizon.
cc
c                    if(iedgehor(ialf,ibet).eq.10)then
c                      if(ibet.lt.ibetlim(ialf))then    
c                        y1=saveflum(ialf,ibet+1)
c                        y0=yflux(icount)
c                        slope=0.5d0*phiabs/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y0+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      else
c                        y1=saveflum(ialf,1)
c                        y0=yflux(icount)
c                        slope=0.5d0*phiabs/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y0+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      endif
c                    endif
c                    if(iedgehor(ialf,ibet).eq.20)then
c                      if(ibet.gt.1)then    
c                        y0=saveflum(ialf,ibet-1)
c                        y1=yflux(icount)
c                        slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y1+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      else
c                        y0=saveflum(ialf,ibetlim(ialf))
c                        y1=yflux(icount)
c                        slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y1+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      endif
c                    endif    ! end if iedgehor > 20
c                  endif      ! end if iedgestar > 10
c                endif        ! end if visib > 1
cc
c                flux=flux+flum(ialf,ibet)    
c                flum(ialf,ibet)=flum(ialf,ibet)
c                igood=igood+1
c 89           continue
cc
cc   Continue with the string accross the break in the coordinates.
cc
c              DO 90 ibet = 1,mark(1,ialf)
c                phiabs=dabs(delphi(ialf,ibet))
c                edgeabs=dabs(delphiedge(ialf,ibet))
c                if((projarray(ialf,ibet).le.0.0d0))then
c                  if(flag.eqv..false.)go to 90
c                  if(isimp.eq.2)then
c                    call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                    totalsum=totalsum+rowsum  
c                  else
c                    call simpson(icount,yflux,rowsum,dphi)
c                    totalsum=totalsum+rowsum+yend(1)+yend(2)
c                  endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
c
c                  if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 7',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 7',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c                  itotal=itotal+icount
c                  icount=0
c                  flag=.false.
c                  phisum=0.0d0
c                  go to 90
c                endif
c                tkkelv=temp(ialf,ibet)/1000.0d0
c                C3 = C2/(wavemu*tkkelv)
c                flum(ialf,ibet)=C1/(dexp(c3)-1.0d0)/wavemu**5
c                dark=(1.0d0-flimbx+flimbx*projarray(ialf,ibet))
c                if(ilaw.eq.2)dark=dark-flimby*projarray(ialf,ibet)*
c     %               dlog(projarray(ialf,ibet))
c                if(ilaw.eq.3)dark=dark-flimby*(1.0-dsqrt(projarray(ialf,ibet)))
c                flum(ialf,ibet)=flum(ialf,ibet)*dark
c                rinty(ialf,ibet)=flum(ialf,ibet) ! save intys for plotting
c                saveflum(ialf,ibet)=surf(ialf,ibet)*flum(ialf,ibet)*
c     $            projarray(ialf,ibet)
c                flum(ialf,ibet)=surf(ialf,ibet)
c     #             *flum(ialf,ibet)*visib(ialf,ibet)
cc
c                if((visib(ialf,ibet).le.0.0d0))then
c                  if(flag.eqv..false.)go to 90
c                  if(isimp.eq.2)then
c                    call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                    totalsum=totalsum+rowsum  
c                  else
c                    call simpson(icount,yflux,rowsum,dphi)
c                    totalsum=totalsum+rowsum+yend(1)+yend(2)
c                  endif
cc
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
cc
c
c                  if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 8',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 8',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c                  itotal=itotal+icount
c                  icount=0
c                  flag=.false.
c                  rowsum=0.0d0
c                  yend(1)=0.0d0
c                  yend(2)=0.0d0
c                  phisum=0.0d0
c                endif
c                if((visib(ialf,ibet).gt.0.0d0))then
c                  flag=.true.
c                  icount=icount+1
c                  yflux(icount)=flum(ialf,ibet)
c                  xflux(icount)=phisum
c                  phisum=phisum+2.0d0*dphi
c                  if(icount.eq.1)then
c                    if(iedgestar(ialf,ibet).ge.10)then 
cc
cc   If this point is the first in the string and at the horizon, then
cc   the connecting end piece is a triangle.  The area is 0.5 times the
cc   y-value times the phi distance.  
cc
c                      yvalue(1)=0.0d0
c                      yend(1)=0.25d0*yflux(1)/dphi*edgeabs
c                      xend(1)=xflux(1)-edgeabs
c                      x1assign=.true.
c                    else
c                      if(iedgehor(ialf,ibet).eq.10)then
cc
cc   If the first point is next to an eclipsing horizon, then the end
cc   piece will be a trapizoid.  Use interpolation to find the y-value
cc   at the horizon.
cc
c                        if(ibet.lt.ibetlim(ialf))then    
c                          y1=saveflum(ialf,ibet+1)
c                          y0=yflux(icount)
c                          slope=0.5d0*phiabs/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y0+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        else
c                          y1=saveflum(ialf,1)
c                          y0=yflux(icount)
c                          slope=0.5d0*phiabs/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y0+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        endif
c                      endif
c                      if(iedgehor(ialf,ibet).eq.20)then
cc
cc   If the first point is next to an eclipsing horizon, then the end
cc   piece will be a trapizoid.  Use interpolation to find the y-value
cc   at the horizon.
cc
c                        if(ibet.gt.1)then    
c                          y0=saveflum(ialf,ibet-1)
c                          y1=yflux(icount)
c                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y1+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        else
c                          y0=saveflum(ialf,ibetlim(ialf))
c                          y1=yflux(icount)
cc                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                          yinterp=y0+slope*(y1-y0)
c                          yvalue(1)=yinterp
c                          yend(1)=0.5d0*(y1+yinterp)*slope
c                          xend(1)=xflux(icount)-phiabs
c                      x1assign=.true.
c                        endif
c                      endif
c                    endif
c                  endif ! end if icount=1
c                  if(iedgestar(ialf,ibet).ge.10)then 
cc
cc   The point in the string (not the first) is at the horizon.  The
cc   second end piece is a triangle.
cc
c                    yvalue(2)=0.0d0
c                    yend(2)=0.25d0*yflux(icount)/dphi*edgeabs
c                    xend(2)=xflux(icount)+edgeabs
c                    x2assign=.true.
c                  else
cc
cc   Otherwise, check to see if it near an eclipsing horizon.
cc
c                    if(iedgehor(ialf,ibet).eq.10)then
c                      if(ibet.lt.ibetlim(ialf))then    
c                        y1=saveflum(ialf,ibet+1)
c                        y0=yflux(icount)
c                        slope=0.5d0*phiabs/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y0+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      else
c                        y1=saveflum(ialf,1)
c                        y0=yflux(icount)
c                        slope=0.5d0*phiabs/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y0+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      endif
c                    endif
c                    if(iedgehor(ialf,ibet).eq.20)then
c                      if(ibet.gt.1)then    
c                        y0=saveflum(ialf,ibet-1)
c                        y1=yflux(icount)
c                        slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y1+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      else
c                        y0=saveflum(ialf,ibetlim(ialf))
c                        y1=yflux(icount)
c                        slope=0.5d0*dabs(phiabs-2.0d0*dphi)/dphi
c                        yinterp=y0+slope*(y1-y0)
c                        yvalue(2)=yinterp
c                        yend(2)=0.5d0*(y1+yinterp)*slope
c                        xend(2)=xflux(icount)+phiabs
c                    x2assign=.true.
c                      endif
c                    endif
c                  endif
c                endif    ! endif visib > 0
cc
c                flux=flux+flum(ialf,ibet)   
c                flum(ialf,ibet)=flum(ialf,ibet)
c                igood=igood+1
c 90           continue
c            endif
c 10       continue
cc
c            if(flag.eqv..true.)then
c              if(isimp.eq.2)then
c                call simpson1(icount,xflux,yflux,xend,yvalue,rowsum,dphi) 
c                totalsum=totalsum+rowsum  
c              else
c                call simpson(icount,yflux,rowsum,dphi)
c                totalsum=totalsum+rowsum+yend(1)+yend(2)
c              endif
cc              if(jj.eq.3.and.ialf.eq.5)then
cc                call dump(33,icount,xflux,yflux)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc                call dump(34,2,xend,yvalue)
ccc                write(*,*)rowsum,dphi,rowsum+yend(1)+yend(2)
cc              endif
c 
c                 if(x1assign.eqv..false.)write(*,*)'xend(1) not assigned spot 9',ialf,ibet
c                  if(x2assign.eqv..false.)write(*,*)'xend(2) not assigned spot 9',ialf,ibet
c                  x1assign=.false.
c                  x2assign=.false.
c
c
c
c              itotal=itotal+icount
c              icount=0
c              rowsum=0.0d0
c              yend(1)=0.0d0
c              yend(2)=0.0d0
c              flag=.false.
c              phisum=0.0d0
c            endif
cc
cc   Scale the light curve by the integral of the limb darkening law
cc   for compatibility with Wilson-Devinney.
cc        
c          totalsum=pie*totalsum/dint
c          flux=totalsum
cc
c          rldint=dint
cc
cc
c 666      format(a6,2x,i2,2x,f9.7,2x,2(i3,1x))
c 667      format(3(i2,2x),2x,f8.6,2x,f11.6,3x,i4)
c
          return
          end
c
c   ************************
c
          subroutine simpson(N,y,total,dphi)
c
c   April 2, 2001
c
c   This routine will sum the y-values in y(1:N) using Simpson's rule.
c
          implicit double precision (a-h,o-z)
          dimension y(N)
c
c   UPDATE March 26, 2002
c
c   Assign a variable fake the value of dphi to supress compiler
c   warning about an unusued variable.
c
          fake=dphi
c
          over3=1.0d0/3.0d0
          over8=3.0d0/8.0d0
          over6=7.0d0/6.0d0
          over24=23.0d0/24.0d0
c
          total=0.0d0
c
          if(N.eq.1)then
            total=y(1)
c     $          +0.5*y(1)
            return
          endif
c
          if(N.eq.2)then
            total=0.5d0*(y(1)+y(2))
c     $          +0.5d0*(y(1)+y(2))
            return
          endif
c
          if(N.eq.3)then
            total=(over3*y(1)+4.0d0*over3*y(2)+over3*y(3))
c     &           +0.5d0*(y(1)+y(3))
            return
          endif
c
          if(N.eq.4)then
            total=(over8*y(1)+3.0d0*over8*y(2)+3.0d0*over8*y(3)
     $            +over8*y(4))
c     $            +0.5d0*(y(1)+y(4))
            return
          endif
c
          if(N.eq.5)then
            total=(over3*y(1)+4.0d0*over3*y(2)+2.0*over3*y(3)
     $            +4.0d0*over3*y(4)+over3*y(5))
c     %            +0.5d0*(y(1)+y(5))
            return
          endif
c
          if(N.gt.5)then
            do 30 i=1,N
              weight=1.0d0
              if(i.eq.1)weight=over8
              if(i.eq.N)weight=over8
              if(i.eq.2)weight=over6
              if(i.eq.N-1)weight=over6
              if(i.eq.3)weight=over24
              if(i.eq.N-2)weight=over24
              total=total+weight*y(i)
 30         continue
c              total=total+0.25d0*(y(1)+y(N))
            return
          endif
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&
c
c
c  RVG BUG ALERT   May 2, 2001
c
c  Move spline and splint from optimizesubs.for to here:
c
         subroutine spline(x,y,n,yp1,ypn,y2)
c
c   November 12, 1999
c
c   This is a spline interpolation routine taken from NUMERICAL RECIPES.
c
         implicit double precision (a-h,o-z)
c
         integer n,NMAX
c         REAL*8 yp1,ypn,x(n),y(n),y2(n)
         dimension x(n),y(n),y2(n)
         parameter(NMAX=900000)
         integer i,k
c         REAL*8 p,qn,sig,un,u(NMAX)
         dimension u(NMAX)
         if(yp1.gt.0.99d30)then
           y2(1)=0.0
           u(1)=0.0
         else
           y2(1)=-0.5
           u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
         endif
         do 11 i=2,n-1
           sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
           p=sig*y2(i-1)+2.0
           y2(i)=(sig-1.0)/p
           u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     #       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
 11      continue
         if(ypn.gt..99e30)then
           qn=0.0
           un=0.0
         else
           qn=0.5
           un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
         endif
         y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
         do 12 k=n-1,1,-1
           y2(k)=y2(k)*y2(k+1)+u(k)
 12      continue
         return
         end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
         subroutine splint(xa,ya,y2a,n,x,y)
c
c   November 12, 1999
c
c   This is a spline interpolation routine taken from Numerical Recipes.
c
         implicit double precision (a-h,o-z)
c
         integer n
c         real*8 x,y,xa(n),y2a(n),ya(n)
         dimension xa(n),y2a(n),ya(n)
         integer k,khi,klo
c         real*8 a,b,h
         klo=1
         khi=n
 1       if(khi-klo.gt.1)then
           k=(khi+klo)/2
           if(xa(k).gt.x)then
             khi=k
           else
             klo=k
           endif
           go to 1
         endif
         h=xa(khi)-xa(klo)
         if(h.eq.0.0)pause 'bad xa input in splint'
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         y=a*ya(klo)+b*ya(khi)+
     $     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0
         return
         end
c
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
         subroutine fastsplint(xa,ya,y2a,n,x,y,klo,khi)
c
c   November 12, 1999
c
c   This is a spline interpolation routine taken from Numerical Recipes.
c
         implicit double precision (a-h,o-z)
c
         integer n
c         real*8 x,y,xa(n),y2a(n),ya(n)
         dimension xa(n),y2a(n),ya(n)
         integer k,khi,klo
c         real*8 a,b,h
c         klo=1
c         khi=n
c         k=(khi-klo)/2
c         if((xa(k).gt.x).and.(khi-klo.gt.1))then
c           go to 2
c         else
c           khi=n
c         endif
 1       if(khi-klo.gt.1)then
           k=(khi+klo)/2
           if(xa(k).gt.x)then
             khi=k
           else
             klo=k
           endif
           go to 1
         endif
 2       h=xa(khi)-xa(klo)
         if(h.eq.0.0d0)then
           write(*,*)klo,khi,n
           pause 'bad xa input in fastsplint'
         endif
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         y=a*ya(klo)+b*ya(khi)+
     $     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
         return
         end
c
c   &&&&&&&&&&&&&&&&&
c
          subroutine simpson1(N,x,y,xend,yend,total,dphi)
c
c    April 2, 2001
c
c   This routine will sum the y-values in y(1:N) using splines.
c
          implicit double precision (a-h,o-z)
          dimension y(N),x(N),xend(2),yend(2),y2(1000)
          dimension xdummy(1000),ydummy(1000),w(8),wx(8)
c
c
          xdummy(1)=xend(1)
          ydummy(1)=yend(1)
c
          do 10 i=1,N
            xdummy(i+1)=x(i)
            ydummy(i+1)=y(i)
 10       continue
c
          xdummy(N+2)=xend(2)
          ydummy(N+2)=yend(2)
c
          newN=N+2
c
          yp1=1.0d31

          call spline(xdummy,ydummy,newN,yp1,yp1,y2)
c
          m=8*newN+1
          h=(xdummy(newN)-xdummy(1))/dble(m-1)
c
          over3=1.0d0/3.0d0
          fover3=4.0d0/3.0d0
          tover3=2.0d0/3.0d0
c
          total=0.0d0
          xstart=xdummy(1)
          klo=1
          khi=m

          wx(1)=0.9602898565d0
          wx(2)=-0.9602898565d0
          w(1)=0.1012285363d0
          w(2)=0.1012285363d0
c
          wx(3)=0.7966664774d0
          wx(4)=-0.7966664774d0
          w(3)=0.2223810345d0
          w(4)=0.2223810345d0
c
          wx(5)=0.5255324099d0
          wx(6)=-0.5255324099d0
          w(5)=0.3137066459d0
          w(6)=0.3137066459d0
c
          wx(7)=0.1834346425d0
          wx(8)=-0.1834346425d0
          w(7)=0.3626837834
          w(8)=0.3626837834
c
          a=xdummy(1)
          b=xdummy(newN)
c
          do 30 i=1,8
            xxx=0.5d0*(b+a+wx(i)*(b-a))
            call splint(xdummy,ydummy,y2,newN,xxx,yout)
            total=total+w(i)*yout
 30       continue
          total=total*0.5d0*(b-a)
          total=0.5d0*total/dphi
c
c          do 20 i=1,m
c            xxx=xstart+dble(i-1)*h
c            weight=1.0d0
c            if(mod(i,2).eq.0)weight=fover3
c            if(mod(i,2).eq.1)weight=tover3
c            if(i.eq.1)weight=over3
c            if(i.eq.m)weight=over3
cc            call fastsplint(xdummy,ydummy,y2,newN,xxx,yout,klo,m)
c            call splint(xdummy,ydummy,y2,newN,xxx,yout)
c            total=total+weight*yout
c 20       continue
c
c          total=0.5d0*dabs(h)*total/dphi
          return
          end
c
c &&&&&&&&&&&&&&&&&&&&
c
c  RVG BUG ALERT   May 2, 2001
c
c  These are two new functions:
c
          double precision function diskxtran(xx,yy,zz,phase,
     &                fincr,Q,istar,bdist)
c
c   will return the coordinate of a point (xx,yy,zz) projected on the sky
c
c   (xx,yy,zz) refers to the coordinates in the rotating system
c
c   UPDATE March 26, 2002
c
c   The variables zz and fincr are not used in diskxtran (they are in
c   diskytran).  To suppress compiler warnings, define two fake variables
c   and assign their values to fincr and zz.  This is simpler than
c   changing the argument list of calls to diskxtran.
c
c   Also, get rid of historical text (comment out code).
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
          fake1=zz
          fake2=fincr
c
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
c
          PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians

          qphaser=phaser
c
          diskxtran=-(xx*dsin(qphaser)+yy*dcos(qphaser))+
     $      bdist*(overQ/(1.0d0+overQ))*dsin(qphaser)
c
          return
          end
c
c  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          double precision function diskytran(xx,yy,
     &        zz,phase,fincr,Q,istar,bdist)
c
c   will return the coordinate of a point (xx,yy,zz) projected on the sky
c
c   (xx,yy,zz) refers to the coordinates in the rotating system
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
          PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians

          qphaser=phaser
c
          diskytran=-(-xx*dcos(fincr)*dcos(qphaser)+
     #      yy*dcos(fincr)*dsin(qphaser)
     $      -zz*dsin(fincr))+
     $     bdist*(-(overQ/(1.0d0+overQ))*dcos(fincr)*dcos(qphaser))
c
c          if(phaser.gt.180.0d0)ytran=-ytran
c
c   Added February 8, 2001
c
c          ytran=ytran*bdist
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   UPDATE March 26, 2002
c
c   Get rid of this unused routine.
c
c          subroutine loadspot(Ns1,Ns2,Nd,spot1parm,spot2parm,
c     %       spotdparm,ispot)
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine addstarspot(istar,ialphmax,ibetmax,
     $       Nalph,ibetlim,tmatrix,spotparm,ave1,ave2,omega,phase,
     %       phiar,mmdx,ispotprof)
c
c   This routine will assign the temperatures of the grid points
c   of the stars that are covered by spots.  The underlying temperatures
c   are simply scaled by the temperature spot factor.  
c   
c   UPDATE March 26, 2002
c
c   Get rid of Nbet from the argument list of addstarspot.  Its value
c   is contained within ibetlim.
c
c    UPDATE June 17, 2002
c
c    Add these dummy assignments for phase and omega to supress
c    compiler warnings about unused variables.
c
c    UPDATE October 13, 2008
c
c    Add ispotprof flag:
c
c    ispotprof=0    temperature factor same as before
c    ispotprof=1    linear profile for the temperature factor
c    ispotprof=2    Gaussian profile for temperature factor
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265359879323d0)
          dimension tmatrix(ialphmax*ibetmax),phiar(ialphmax*ibetmax),
     %       ibetlim(ialphmax),spotparm(2,4),mmdx(ialphmax,ibetmax)
c
c
c    UPDATE June 17, 2002
c
c    Add these dummy assignments for phase and omega to supress
c    compiler warnings about unused variables.
c
          ddd1=phase
          ddd2=omega
c
          radcon=pie/180.0d0
          halfpie=0.5*pie
c
          fac1=spotparm(1,1)
          fac2=spotparm(2,1)
          rlat1=radcon*spotparm(1,2)-halfpie
          rlat2=radcon*spotparm(2,2)-halfpie
          rlong1=radcon*spotparm(1,3)
          rlong2=radcon*spotparm(2,3)
          if(rlong1.gt.pie)rlong1=rlong1-2.0d0*pie
          if(rlong2.gt.pie)rlong2=rlong2-2.0d0*pie
          rad1=radcon*spotparm(1,4)
          rad2=radcon*spotparm(2,4)
c
          if((fac1.lt.0.0d0).and.(fac2.lt.0.0d0))return !no valid factors
c
          dtheta=pie/dble(Nalph)
          summ1=0.0d0
          summ2=0.0d0
          icount1=0
          icount2=0
          icount3=0
          DO 10 IALF = 1, nalph
            theta=-0.5d0*dtheta+dtheta*dble(ialf)
            rlat=theta-halfpie
            DO 9 IBET = 1, ibetlim(ialf)    !4*NBET
c
              icount3=icount3+1
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              rlong=phiar(iidx)
              if(rlong.gt.pie)rlong=rlong-2.0d0*pie
c
              sepang1=dsin(rlat)*dsin(rlat1)+dcos(rlat)*dcos(rlat1)*
     %           dcos(rlong-rlong1)
              sepang1=dabs(dacos(sepang1))
              if((sepang1.le.rad1).and.(fac1.gt.0.0d0))then
                if(ispotprof.eq.0)then
                  icount1=icount1+1
                  tmatrix(iidx)=tmatrix(iidx)*fac1
                  summ1=summ1+tmatrix(iidx)
                endif
                if(ispotprof.eq.1)then
                  if(rad1.le.0.0d0)then
                    go to 9
                  else
                    icount1=icount1+1
                    slope=(1.0d0-fac1)/rad1
                    fff=slope*(sepang1-rad1)+1.0d0
                    tmatrix(iidx)=tmatrix(iidx)*fff
                    summ1=summ1+tmatrix(iidx)
                  endif
                endif
                if(ispotprof.eq.2)then
                  if(rad1.le.0.0d0)then
                    go to 9
                  else
                    icount1=icount1+1
                    fff=(fac1-1.0d0)*dexp(-4.5d0*(sepang1/rad1)**2)+1.0d0
                    tmatrix(iidx)=tmatrix(iidx)*fff
                    summ1=summ1+tmatrix(iidx)
                  endif
                endif
              endif
c
              sepang2=dsin(rlat)*dsin(rlat2)+dcos(rlat)*dcos(rlat2)*
     %           dcos(rlong-rlong2)
              sepang2=dabs(dacos(sepang2))
              if((sepang2.le.rad2).and.(fac2.gt.0.0d0))then
                if(ispotprof.eq.0)then
                  icount2=icount2+1
                  tmatrix(iidx)=tmatrix(iidx)*fac2
                  summ2=summ2+tmatrix(iidx)
                endif
                if(ispotprof.eq.1)then
                  if(rad2.le.0.0d0)then
                    go to 9
                  else
                    icount2=icount2+1
                    slope=(sepang2-fac2)/rad2
                    fff=slope*(1.0d0-rad2)+1.0d0
                    tmatrix(iidx)=tmatrix(iidx)*fff
                    summ2=summ2+tmatrix(iidx)
                  endif
                endif
                if(ispotprof.eq.2)then
                  if(rad2.le.0.0d0)then
                    go to 9
                  else
                    icount2=icount2+1
                    fff=(fac2-1.0d0)*dexp(-4.5d0*(sepang2/rad2)**2)+1.0d0
                    tmatrix(iidx)=tmatrix(iidx)*fff
                    summ2=summ2+tmatrix(iidx)
                  endif
                endif
              endif

 9          continue
 10       continue
c
          ave1=0.0d0
          ave2=0.0d0
          if(icount1.gt.0)ave1=summ1/dble(icount1)
          if(icount2.gt.0)ave2=summ2/dble(icount2)

          if(ave1.gt.0.0d0)write(2,100)istar,ave1,icount1,icount3
          if(ave2.gt.0.0d0)write(2,101)istar,ave2,icount2,icount3

c
 100      format(/'star ',i1,', spot 1:    average temperature ',
     &       f9.3,',',/19x,'number of grid points = ',i6,' out of ',i6)
 101      format(/'star ',i1,', spot 2:    average temperature ',
     &       f9.3,',',/19x,'number of grid points = ',i6,' out of ',i6)
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   UPDATE March 26, 2002
c
c   Remove betarim, separ, tdisk, xi, dx, dy, dz, drad, 
c   xedge, yedge, zedge, stepr, stepz, bdist, omega, phase
c   from the argument list, and from the appropriate dimension
c   statements below.
c
          subroutine adddiskspot(Nthetamax,Nrmax,Ntheta,Nradius,
     %       rinner,router,reff2,Rl2,
     $       dtemp,
     $       tedge,redge,
     #       ivrt,reper,rsper,spotparm,ave1,ave2)
c
c   May 8, 2001
c
c   This routine will modify the temperatures on the disk that are
c   within spots.
c
          implicit double precision(a-h,o-z)
c
          parameter (pie=3.14159265358979323d0)
          dimension dtemp(Nrmax*Nthetamax),tedge(Nthetamax*11)
          dimension spotparm(2,4)
c
c   Start with the disk face.  It is assumed that the lower face is exactly
c   the same as the upper face, but with a negative z-value.  In practice,
c   however, we never see the bottom face.
c          
          radcon=pie/180.0d0
c
          if(ivrt.eq.0)then
            redge=router*reff2              ! radius of outer edge in x units
            rsmall=rinner*Rl2               ! radius of inner edge in x units
            reper=redge
            rsper=rsmall
          else
            redge=reper
            rsmall=rsper
          endif
c
          fac1=spotparm(1,1)
          az1=radcon*spotparm(1,2)
          cut1=spotparm(1,3)
          width1=radcon*spotparm(1,4)
c
          fac2=spotparm(2,1)
          az2=radcon*spotparm(2,2)
          cut2=spotparm(2,3)
          width2=radcon*spotparm(2,4)
c

          steptheta=360.0d0/dble(ntheta)
c
c   Transform r into zeta.
c
          zetain=2.0d0*dsqrt(rsmall)
          zetaout=2.0d0*dsqrt(redge)
          stepzeta=(zetaout-zetain)/dble(Nradius-1)
c
          theta=0.0d0
          zeta=zetain
          icount1=0
          icount2=0
          summ1=0.0d0
          summ2=0.0d0
c
          if((cut1.ge.1.0d0).and.(cut2.ge.1.0d0))go to 99
          if((cut1.ge.1.0d0).and.(fac2.le.0.0d0))go to 99
          if((cut2.ge.1.0d0).and.(fac1.le.0.0d0))go to 99
c
          do 10 ir=1,Nradius
            zeta=zetain+dble(ir-1)*stepzeta
            r=0.25d0*zeta*zeta
            ratt=r/redge
            do 9 ithet=1,Ntheta              ! theta goes from zero to 360-step
              theta=dble(ithet)*steptheta -0.5*steptheta  ! degrees
              thetar=theta*radcon            ! radians
              angdiff1=dacos(dcos(thetar-az1))
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(ir-1)*Ntheta+ithet
c

              if((angdiff1.le.width1).and.(fac1.gt.0.0d0).and.
     #            (ratt.ge.cut1))then
                dtemp(iidx)=dtemp(iidx)*fac1
                icount1=icount1+1
                summ1=summ1+dtemp(iidx)
              endif
              angdiff2=dacos(dcos(thetar-az2))
              if((angdiff2.le.width2).and.(fac2.gt.0.0d0).and.
     #            (ratt.ge.cut2))then
                dtemp(iidx)=dtemp(iidx)*fac2
                icount2=icount2+1
                summ2=summ2+dtemp(iidx)
              endif
 9          continue
 10       continue
c
 99       do 20 iz=-5,5
            do 19 ithet=1,Ntheta             ! theta goes from zero to 360-step
              theta=dble(ithet)*steptheta-0.5*steptheta   ! degrees
              thetar=theta*radcon         ! radians
              angdiff1=dacos(dcos(thetar-az1))
c
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
              iidx=(iz+6-1)*Ntheta+ithet
c
              if((angdiff1.le.width1).and.(fac1.gt.0.0d0))then
                tedge(iidx)=tedge(iidx)*fac1
                icount1=icount1+1
                summ1=summ1+tedge(iidx)
              endif
              angdiff2=dacos(dcos(thetar-az2))
              if((angdiff2.le.width2).and.(fac2.gt.0.0d0))then
                tedge(iidx)=tedge(iidx)*fac2
                icount2=icount2+1
                summ2=summ2+tedge(iidx)
              endif
 19         continue
 20       continue
c
          ave1=0.0d0
          ave2=0.0d0
          if(icount1.gt.0)ave1=summ1/dble(icount1)
          if(icount2.gt.0)ave2=summ2/dble(icount2)

          if(ave1.gt.0.0d0)write(2,100)ave1,icount1
          if(ave2.gt.0.0d0)write(2,101)ave2,icount2
c
 100      format(/'disk spot 1:    average temperature ',
     &       f9.3,',  number of grid points = ',i4)
 101      format(/'disk spot 2:    average temperature ',
     &       f9.3,',  number of grid points = ',i4)
c
          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine gravexp(teff,tgrav,istar)
c 
c   May 16, 2001
c
c   This routine will return a value of the gravity darkening exponent
c   Tgrav based on the input effective temperature teff.
c
c   Reference:  A. Claret, 2000, A&A, 359, 289
c
          implicit double precision (a-h,o-z)
c
          dimension xarray(45),yarray(45),y2(45)
c
c   xarray = log(t),  yarray = beta_1, where Tgrav=beta_1/4.0
c

          data xarray/3.3010,3.3640,3.4217,3.4490,3.4849,3.5097,3.5180,3.5361,
     $     3.5470,3.5640,3.5920,3.6280,3.6633,3.7031,3.7424,3.7557,
     $     3.7707,3.7844,3.7964,3.8116,3.8251,3.8458,3.8553,3.8648,
     $     3.8736,3.8802,3.8868,3.8912,3.8956,3.9022,3.9088,3.9491,
     $     4.0249,4.0960,4.1630,4.2285,4.2891,4.3468,4.4008,4.4551,
     $     4.4966,4.5376,4.5749,4.6072,4.6350/
c
          data yarray/0.2150,0.2200,0.2150,0.2000,0.1900,0.1800,0.1700,0.1800,
     $     0.1800,0.1900,0.2300,0.3400,0.4100,0.4345,0.4048,0.3888,
     $     0.3675,0.3463,0.3239,0.3004,0.2753,0.2499,0.2047,0.1599,
     $     0.1731,0.3150,0.5862,0.7292,0.8160,0.9350,0.9435,0.9857,
     $     0.9985,0.9962,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,
     $     1.0000,1.0000,1.0000,1.0000,1.0000/
c
          tlog=dlog10(teff)
c
          if(tlog.lt.3.301d0)then
            tgrav=0.25d0*0.2150
            write(2,100)istar,tgrav
            return
          endif
c
          if(tlog.gt.4.2d0)then
            tgrav=0.25d0
            write(2,100)istar,tgrav
            return
          endif
c
c   Set up the splines.
c
          yp1=1.0d31
          call spline(xarray,yarray,45,yp1,yp1,y2)
          call splint(xarray,yarray,y2,45,tlog,yout)
c
          tgrav=0.25d0*yout
          write(2,100)istar,tgrav
c
 100      format('Info:  The gravity darkening exponent for star',i1,
     $      ' has been set to ',f9.7)
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
c   UPDATE March 26, 2002
c
c   Remove rpole1,rpole2,fill1,fill2 from the argument list.
c
          subroutine lineparms1(Teff2,Q,
     %        finc,separ,period,reff1,reff2,
     $        vrot1,vrot2,omega1,omega2,
     $        bdist,ecc,SA3,ave11,ave12,ave21,ave22,ave1,ave2,parmstring,
     $        pot1,pot2)
c
c   July 13, 2001
c
c   This subroutine is similar to parms1, except that it writes the
c   computed physical quantities on a single line.
c
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm to 9.
c
c   UPDATE January 16, 2001
c
c   add pot1,pot2 to the end of the list
c
          implicit double precision(a-h,o-z)
c
c   UPDATE March 26, 2002
c
c   obsparm is not used here, so comment out its declaration.
c
c          dimension obsparm(9)
c
c   UPDATE November 28, 2001
c
c   Change parmstring from character*199 to character*201
c
c
c   UPDATE January 16, 2002
c
c   parmstring was character*201, now should be character*227
c
c   UPDATE June 7, 2002
c
c   Make the length of parmstring 237
c
c
c   UPDATE October 28, 2002
c
c   Make the length of parmstring character*249
c
c   UPDATE October 22, 2008
c
c   Make the length of parmstring character*259
c
          character*259 parmstring
c
          parameter(pie=3.14159265358979323d0)
c
           fincr=finc*3.141592653589879d0/180.0d0        ! radians
           ppp=period*24.0d0               ! period in hours
           coef=3.518847d10
           coef3=4.35713636d31
           sifinc=dsin(fincr)
c
c   Use the formula separ = coef*(perid*period*total_mass)**(1/3) to
c   solve for the total mass in solar masses.  The separation is
c   entered in solar masses, so (R_sun/coef)**3=7.737294491.
c
           total_mass=(separ)**(3)*7.737294491d0/(ppp*ppp)

c   Use the value of GM_sun found from the solar system.
c
           gmsun=1.32712440018d20  !mks units
           solarrad=6.9598d8
           p=period*86400.0d0
           rM1=(separ*solarrad)**(3)*4.0d0*pie*pie
           rM1=rM1/(gmsun*p*p*(1.0d0+Q))
           rM2=Q*rM1
c 
           R1=reff1*separ
           if(Teff2.gt.0.0d0)then
             R2=reff2*separ
           else
             R2=0.00d0
           endif
c
           gsun=2.739910d4
           gpole1=gsun*rM1/(R1*R1)  
           if(Teff2.gt.0.0d0)then
             gpole2=gsun*rM2/(R2*R2)  
           else
             gpole2=1.0d0
           endif
c
c           gscale1=27397.726d0*rM1/(separ*separ*bdist*bdist)
c           gscale2=27397.726d0*rM2/(separ*separ*bdist*bdist)
c
           fact=solarrad*2.0d0*pie/86400.0d0/1.0d3  !50.613093d0
c
           vrot1=fact*R1/period*sifinc
           vrot2=fact*R2/period*sifinc
c
           a2=separ/(1.0d0+Q)*bdist
c
c   UPDATE September 12, 2001
c 
c   Bug fix, change a1=(separ-a2)*bdist  to  a1=separ(1.0d0-bdist/(1.0d0+Q))
c

           a1=separ*(1.0d0-bdist/(1.0d0+Q))
c
c   UPDATE MARCH 4, 2005
c
c   Change a2 and a1 below.
c
           a2=separ/(1.0d0+Q)
           a1=separ-a2
           efact=1.0d0/dsqrt(1.0-ecc*ecc)
           velK1=fact*a1/period*sifinc*efact
           velK2=fact*a2/period*sifinc*efact
c
           ttt=SA3
           if(ttt.lt.0.0d0)ttt=0.0d0
           R3=R1*dsqrt(ttt)

c
c   NEW BUG August 10, 2001
c
c   Add a correction factor to the rotational velocities in the
c   case of eccentric orbits.
c
           hutfac=(1.0d0+7.5d0*ecc*ecc+5.625d0*ecc**4+
     #          0.3125d0*ecc**6)/((1.0d0+3.0d0*ecc*ecc+
     $          3.0d0/8.0d0*ecc**4)*dsqrt((1.0d0-ecc*ecc)**3))
c
c   UPDATE January 16, 2001
c
c   add pot1,pot2 to the end of the list
c
c   UPDATE OCTOBER 21, 2005
C
C   Modify pot2 for consistency with ELC.out.
c
c
           overQ=1.0d0/Q
           pppp=pot2/overQ+0.5d0*(overQ-1.0d0)/overQ
           write(parmstring,5000)rM1,rM2,R1,R2,R3,
     $        dlog10(gpole1),dlog10(gpole2),
     #        a1,a2,separ,velk1,velk2,vrot1,vrot2,
     $        vrot1*omega1*hutfac,vrot2*omega2*hutfac,
     $        ave11,ave12,ave21,ave22,ave1,ave2,pot1,pppp
c
c  UPDATE November 28, 2001
c
c  Change the format statement below.  The first field was 2(f7.4,1x)
c  and now should be 2(f8.4,1x)
c
c  UPDATE January 16, 2001
c
c  Add 2(f12.6,1x) to the end
c
c  UPDATE October 28, 2002
c 
c  Change the first two format statements to f11.8
c  Also change 2(f8.3,1x) at the end (velK) to 2(f10.5,1x)
c
 5000      format(2(f14.8,1x),3(f12.6,1x),2(f7.3,1x),3(f11.6,1x),2(f10.5,1x),
     &        4(f7.3,1x),6(f7.1,1x),2(f12.6,1x))
c
c           write(3,1000)rM1
c           write(3,2000)rM2
c           write(3,1001)R1
c           write(3,2001)R2
c           write(3,1002)dlog10(gpole1)
c           write(3,2002)dlog10(gpole2)
c           write(3,1003)a1
c           write(3,1004)a2
c           write(3,1005)separ
c           write(3,1006)velK1
c           write(3,2006)velK2
c           write(3,1007)vrot1
c           write(3,2007)vrot2
c           write(3,1008)vrot1*omega1
c           write(3,2008)vrot2*omega2
c
 1000      format( f9.5,11x,'mass of star 1 in solar masses')
 2000      format( f9.5,11x,'mass of star 2 in solar masses')
 1001      format(f10.5,10x,'radius of star 1 in solar radii')
 2001      format(f10.5,10x,'radius of star 2 in solar radii')
 1002      format( f9.6,11x,'log(g) of star 1, cgs')
 2002      format( f9.6,11x,'log(g) of star 2, cgs')
 1003      format(f11.5, 9x,'a1 (solar radii)')
 1004      format(f11.5, 9x,'a2 (solar radii)')
 1005      format(f11.5, 9x,'a (solar radii)')
 1006      format( f9.4,11x,'K_1 (km/sec)')
 2006      format( f9.4,11x,'K_2 (km/sec)')
 1007      format( f9.4,11x,'V1_rot*sin(i) (km/sec)')
 2007      format( f9.4,11x,'V2_rot*sin(i) (km/sec)')
 1008      format( f9.4,11x,'V1_rot*sin(i), scaled by omega1 (km/sec)')
 2008      format( f9.4,11x,'V2_rot*sin(i), scaled by omega2 (km/sec)')
c
c
           return
           end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine getXecl(Nhoriz,xhoriz,yhoriz,ixecl,Q,finc,bdist,phase)
c
c   UPDATE September 10, 2001
c
c   This routine will check to see if the center of star 2 is eclipsed
c   by star 1.  ixecl=100 if eclipsed, zero otherwise.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xhoriz(Nhoriz),yhoriz(Nhoriz)
c
          PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians

          xp=xtran(bdist,0.0d0,0.0d0,phase,fincr,Q,istar,bdist) 
          yp=ytran(bdist,0.0d0,0.0d0,phase,fincr,Q,istar,bdist)
c
c   If we are looking at star 2 and there is a disk, then check to see if
c   the points are *inside* the top horizon of the disk.
c
          ixecl=-100
          iyes=-100
          call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
          if(iyes.eq.100)ixecl=100
c
          return
          end
c
c
c============================================================
c
          SUBROUTINE SORT2(N,RA,RB)
c
c   UPDATE September 10, 2001
c
c   This is a new subroutine, similar to sort3.
c
c   Taken from Numerical Recipes.
c
          implicit double precision(a-h,o-z)

          DIMENSION RA(N),RB(N)
c
          L=N/2+1
          IR=N
10        CONTINUE
          IF(L.GT.1)THEN
            L=L-1
            RRA=RA(L)
            RRB=RB(L)
c            rrc=rc(L)
          ELSE
            RRA=RA(IR)
            RRB=RB(IR)
c            rrc=rc(IR)
            RA(IR)=RA(1)
            RB(IR)=RB(1)
c            rc(IR)=RC(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
              RA(1)=RRA
              RB(1)=RRB
c              rc(1)=rrc
              RETURN
            ENDIF
          ENDIF
          I=L
          J=L+L
20        IF(J.LE.IR)THEN
            IF(J.LT.IR)THEN
              IF(RA(J).LT.RA(J+1))J=J+1
            ENDIF
            IF(RRA.LT.RA(J))THEN
              RA(I)=RA(J)
              RB(I)=RB(J)
c              rc(i)=rc(j)
              I=J
              J=J+J
            ELSE
              J=IR+1
            ENDIF
            GO TO 20
          ENDIF
          RA(I)=RRA
          RB(I)=RRB
c          rc(i)=rrc
          GO TO 10
          END
c
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine simplerefl(ialphmax1,ibetmax1,
     $      Nalph1,ibetlim1,Nalph2,ibetlim2,
     $      xarray1,yarray1,zarray1,gradx1,grady1,gradz1,garray1,
     $      xarray2,yarray2,zarray2,gradx2,grady2,gradz2,garray2,
     $      temp1,temp2,dbolx,dboly,ilaw,alb1,alb2,teff1,
     $      teff2,Tgrav1,Tgrav2,rLx,idint,redge,betarim,gpole1,gpole2,
     %      Tpole1,Tpole2,bdist,SA1,SA2,rad1,rad2,separ,mmdx1,mmdx2,
     %      ialphmax2,ibetmax2,isw25)
c
c    UPDATE March 22, 2002
c
c    This routine is a simplified and faster version of the detailed
c    reflection, valid for nearly spherical stars.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xarray1(ialphmax1*ibetmax1),yarray1(ialphmax1*ibetmax1),
     $      zarray1(ialphmax1*ibetmax1),
     &      gradx1(ialphmax1*ibetmax1),temp1(ialphmax1*ibetmax1),
     $      grady1(ialphmax1*ibetmax1),gradz1(ialphmax1*ibetmax1),
     %      garray1(ialphmax1*ibetmax1),
     $      ibetlim1(ialphmax1),rad2(ialphmax2*ibetmax2),
     $      ibetlim2(ialphmax2),rad1(ialphmax1*ibetmax1)
          dimension xarray2(ialphmax2*ibetmax2),yarray2(ialphmax2*ibetmax2),
     $      zarray2(ialphmax2*ibetmax2),
     &      gradx2(ialphmax2*ibetmax2),temp2(ialphmax2*ibetmax2),
     $      grady2(ialphmax2*ibetmax2),gradz2(ialphmax2*ibetmax2),
     %      garray2(ialphmax2*ibetmax2),mmdx1(ialphmax1,ibetmax1)
          dimension dbolx(8,2),dboly(8,2),mmdx2(ialphmax2,ibetmax2)
c
c   Start with star 1 and compute the flux from star 2.
c
          darkbolx1=dbolx(1,1)
          darkboly1=dboly(1,1)
          darkbolx2=dbolx(1,2)
          darkboly2=dboly(1,2)
c
          dtheta1=pie/(1.0d0*nalph1)
          dtheta2=pie/(1.0d0*nalph2)
c
c   Define the integrated  limb darkening coefficients.  The equation is
c
c   dint=2*pi*int_0^1{mu*(1-x*(1-mu))d(mu)}  for the linear law, etc.
c
          dint1=pie*(1.0d0-darkbolx1/3.0d0)
          dint2=pie*(1.0d0-darkbolx2/3.0d0)
          if(ilaw.eq.2)then
            dint1=pie*(1.0d0-darkbolx1/3.0d0+2.0d0*darkboly1/9.0d0)
            dint2=pie*(1.0d0-darkbolx2/3.0d0+2.0d0*darkboly2/9.0d0)
          endif
          if(ilaw.eq.3)then
            dint1=pie*(1.0d0-darkbolx1/3.0d0-darkboly1/5.0d0)
            dint2=pie*(1.0d0-darkbolx2/3.0d0-darkboly2/5.0d0)
          endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
          if(ilaw.eq.4)then
            dint1=pie*(1.0d0-darkbolx1/3.0d0-darkboly1/6.0d0)
            dint2=pie*(1.0d0-darkbolx2/3.0d0-darkboly2/6.0d0)
          endif
c
          if(teff2.gt.0.0d0)then
            C1=(Tpole2/Tpole1)**(4)*(dint1/dint2)
            C1=C1*alb1/dint1
          endif
          nalf12=nalph1/2
          nalf22=nalph2/2
          DIV1 = gpole1    ! gravity at the pole
          DIV2 = gpole2    ! gravity at the pole
c
c   If teff2 < 0, we are in X-ray binary mode.  The parameter rLx is the
c   log10 of the X-ray luminosity.  Compute the bolometric luminosity
c   of star 1, and compute what surface area it would need to have a
c   luminosity equal to that of Lx.   Put this fake surface area as
c   the area of star 2
c
          if(teff2.le.0.0d0)then
            sigma=5.675d-5
            rbol=SA1*(separ*6.9598d10)**2*sigma*teff1**4
            SA2=SA1*10.0d0**(rLx)/rbol
c            SA2=SA1*rLx
            C1=alb1/dint2
          endif
          diff1max=-1.0d0
          diff2max=-1.0d0
c
          T4g2=4.0d0*Tgrav2
          T4g1=4.0d0*Tgrav1

          do 10 ialf=1,Nalph1/2  !Nalph1/2,1,-1
            do 9 ibet=1,ibetlim1(ialf)/2  !4*Nbet1
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim1(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim1)+ibet
              iidx=mmdx1(ialf,ibet)
c
              summ=0.0d0
              term1=(bdist-xarray1(iidx))*gradx1(iidx)
     $                   -(yarray1(iidx))*grady1(iidx)
     #                   -(zarray1(iidx))*gradz1(iidx)
              if(term1.le.0.0d0)go to 10
c
              dist1=dsqrt(rad1(iidx)**2+
     $              bdist*(bdist-2.0d0*xarray1(iidx)))
              foreshort1=(term1/(dist1)) 
c
              summ=foreshort1*(0.25d0*SA2)*(dint2/pie)/(dist1*dist1)
c
c   UPDATE August 7, 2008
c
c   If the flag isw25 is 0, then
c   Add a foreshortening correction for the X-ray heating (the disk
c   is assumed to be a thin disk in the plane, not a point source).
c
c   If isw25 = 1, then assume a point source.
c
              if((teff2.lt.0.0d0).and.(isw25.eq.0))then
                xA=xarray1(iidx)
                yA=yarray1(iidx)
                zA=zarray1(iidx)
                xB=bdist
                yB=0.0d0
                zB=0.0d0
                dist=dsqrt(xA*xA+yA*yA+(zA-bdist)*(zA-bdist))
                xshort=dabs(zA)/dist
c                write(*,6969)zA,xshort*foreshort1
                summ=summ*xshort
              endif
 6969         format(2(f9.6,3x))

              if((idint.ge.1))then
                zrim=redge*dtan(betarim*0.017453293d0)
                xA=xarray1(iidx)
                yA=yarray1(iidx)
                zA=zarray1(iidx)
                xB=bdist
                yB=0.0d0
                zB=0.0d0
                call zheight(xA,yA,zA,xB,yB,zB,redge,zcross)
                if(summ.lt.0.0d0)summ=0.0d0
                if((zA.ge.0.0d0).and.(zcross.lt.zrim))summ=0.0d0
                if((zA.lt.0.0d0).and.(zcross.gt.-zrim))summ=0.0d0
              endif
c                           
              FpBoverFA=C1/((garray1(iidx)/div1)**(T4g1))*summ
              tnew=temp1(iidx)*(1.0d0+FpBoverFA)**0.25d0
              diff=tnew-temp1(iidx)
              if(diff.gt.diff1max)diff1max=diff
              temp1(iidx)=tnew
 9          continue
 10       continue
c
c   Now go to star 2 and compute the irradation from star 1
c
 11       if(teff2.le.0.0d0)go to 50
          C2=(Tpole1/Tpole2)**(4)*(dint2/dint1)
          C2=C2*alb2/dint2
c
          do 20 ialf=1,Nalph2/2  !Nalph2/2,1,-1
            do 19 ibet=1,ibetlim2(ialf)/2   !Nbet2
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim2(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim2)+ibet
c
              iidx=mmdx2(ialf,ibet)
              summ=0.0d0
              term1=(bdist-xarray2(iidx))*gradx2(iidx)
     $                   -(yarray2(iidx))*grady2(iidx)
     #                   -(zarray2(iidx))*gradz2(iidx)
              if(term1.le.0.0d0)go to 20
c
              dist1=dsqrt(rad2(iidx)**2+
     $          bdist*(bdist-2.0d0*xarray2(iidx)))

              foreshort1=(term1/(dist1))    

              summ=foreshort1*(0.25d0*SA1)*(dint1/pie)/(dist1*dist1)
c
              if((idint.ge.1))then
                zrim=redge*dtan(betarim*0.017453293d0)
                xA=xarray2(iidx)
                yA=yarray2(iidx)
                zA=zarray2(iidx)
                xB=bdist
                yB=0.0d0
                zB=0.0d0
                call zheight(xA,yA,zA,xB,yB,zB,redge,zcross)
                if(summ.lt.0.0d0)summ=0.0d0
                if((zA.ge.0.0d0).and.(zcross.lt.zrim))summ=0.0d0
                if((zA.lt.0.0d0).and.(zcross.gt.-zrim))summ=0.0d0
              endif
c                           
              FpAoverFB=C2/((garray2(iidx)/div2)**(T4g2))*summ
              tnew=temp2(iidx)*(1.0d0+FpAoverFB)**0.25d0
              diff=tnew-temp2(iidx)
              if(diff.gt.diff2max)diff2max=diff
              temp2(iidx)=tnew
c
 19         continue
 20       continue
c
c   Now use symmetry to fill in the other quadrants on the star.  
c
 50       continue
c
          DO 401 IALF = 1, nalph1/2
            DO 400 IBET = 1, ibetlim1(ialf)/2
              I1=nalph1-(ialf-1)
              J2=ibetlim1(ialf)-(ibet-1)
              izz=ialf
              jzz=ibet
c              iidx=(izz-1)*ibetlim1(ialf)+jzz
c              iidx=kount(ialphmax,izz,ibetlim1)+jzz
              iidx=mmdx1(izz,jzz)
c
              izz=ialf
              jzz=j2
c              jjdx=(izz-1)*ibetlim1(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim1)+jzz
              jjdx=mmdx1(izz,jzz)
              temp1(jjdx)=temp1(iidx)

              izz=I1
              jzz=ibet
c              jjdx=(izz-1)*ibetlim1(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim1)+jzz
              jjdx=mmdx1(izz,jzz)
              temp1(jjdx)=temp1(iidx)
c
              izz=I1
              jzz=j2
c              jjdx=(izz-1)*ibetlim1(ialf)+jzz
c              iidx=kount(ialphmax,izz,ibetlim1)+jzz
              jjdx=mmdx1(izz,jzz)
              temp1(jjdx)=temp1(iidx)
 400        continue
401       CONTINUE
c
          if(teff2.le.0.0)go to 999
          DO 501 IALF = 1, nalph2/2
            DO 500 IBET = 1, ibetlim2(ialf)/2
              I1=nalph2-(ialf-1)
              J2=ibetlim2(ialf)-(ibet-1)
              izz=ialf
              jzz=ibet
c              iidx=(izz-1)*ibetlim2(ialf)+jzz
c              iidx=kount(ialphmax,izz,ibetlim2)+jzz
              iidx=mmdx2(izz,jzz)

              izz=ialf
              jzz=j2
c              jjdx=(izz-1)*ibetlim2(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim2)+jzz
              jjdx=mmdx2(izz,jzz)
              temp2(jjdx)=temp2(iidx)

              izz=I1
              jzz=ibet
c              jjdx=(izz-1)*ibetlim2(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim2)+jzz
              jjdx=mmdx2(izz,jzz)
              temp2(jjdx)=temp2(iidx)

              izz=I1
              jzz=j2
c              jjdx=(izz-1)*ibetlim2(ialf)+jzz
c              jjdx=kount(ialphmax,izz,ibetlim2)+jzz
              jjdx=mmdx2(izz,jzz)
              temp2(jjdx)=temp2(iidx)

 500        continue
501       CONTINUE
c
          if(diff2max.lt.0.0d0)diff2max=0.0
          if(diff1max.lt.0.0d0)diff1max=0.0

          write(2,72)diff2max
 999      write(2,71)diff1max
c
 71       format(/'maximum temperature change for star 1 = ',f13.5)
 72       format('maximum temperature change for star 2 = ',f13.5)

          return
          end
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine rotkern(ialphmax,ibetmax,Nalf,ibetlim,
     &      istar,omega,phase,finc,
     %      Q,flum,xcoords,ycoords,flux,separ,period,gamma,
     $      rldint,ecc,argrad,visib,extension,mmdx)
c
c    UPDATE May 22, 2002
c
c    This routine will compute a rotational broadening kernel for
c    star istar.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)    
          dimension flum(ialphmax*ibetmax),xcoords(ialphmax*ibetmax),
     &      ycoords(ialphmax*ibetmax),ibetlim(ialphmax),
     &      visib(ialphmax*ibetmax),mmdx(ialphmax,ibetmax)
c
          dimension binx(1001),biny(1001),xscratch(10000),yscratch(10000)
c
          character*9 extension
c
c
c    Check to see if the star has any flux before going further!
c        
          if(flux.le.0.0d0)return
c
          if(istar.eq.1)open(unit=98,file='star1rotkern.'//extension,
     %               status='unknown')
          if(istar.eq.2)open(unit=98,file='star2rotkern.'//extension,
     %               status='unknown')

          argfac=ecc*dcos(argrad)
          efact=1.0d0/dsqrt(1.0d0-ecc*ecc)
          dint=rldint/pie
          vel=0.0d0
          delvel=0.0d0
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
          PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
          siphase=dsin(phaser)
          cophase=dcos(phaser)
          sifinc=dsin(fincr)
c
c   Compute the expected K velocity, which is the circular velocity times 
c   dsin(finc).
c
c    
          a=separ*6.9598d5            !separation in km
          p=period*24.00d0*3600.00d0       !period in seconds
          velamp=2.0d0*pie*a/p*efact
c
c    First, find the limits of the minimum and maximum rotational
c    velocity on the star.
c
          rmin=123456.
          rmax=-123456.
          do 60 ialf=1,Nalf
            do 59 ibet=1,ibetlim(ialf)        !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
              iidx=mmdx(ialf,ibet)
c
              v1=omega*(-xcoords(iidx)*siphase*sifinc-
     %            ycoords(iidx)*cophase*sifinc)
              if(v1*velamp.gt.rmax)rmax=v1*velamp
              if(v1*velamp.lt.rmin)rmin=v1*velamp
 59          continue
 60       continue
c
c   Next, define equal bins in velocity.
c
          nbin=501
          binsize=(rmax-rmin)/dble(nbin-1)
c
          do 1120 i=1,nbin
            binx(i)=dble(i-1)*binsize+rmin
            biny(i)=0.0d0
1120      continue
c
c   Finally, loop over the bins.
c   In the inner loop, loop over Nalph (go along latitude
c   rows) and save up curves of velocity and intensity.  In
c   general, these curves will either intersect the velocity
c   of the bin twice or not at all.  If there are intersections,
c   interpolate the intensities to the velocity of the bin velocity
c   and accumulate the sums.
c
            dtheta=pie/dble(Nalf)
            do 199 ialf=1,Nalf 
              theta=-0.5d0*dtheta+dtheta*dble(ialf)
              do 190 ibet=1,ibetlim(ialf)
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c                iidx=(ialf-1)*ibetlim(ialf)+ibet
c                iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
                iidx=mmdx(ialf,ibet)
                v1=omega*(-xcoords(iidx)*siphase*sifinc-
     %             ycoords(iidx)*cophase*sifinc)
                xscratch(ibet)=velamp*v1
                if(visib(iidx).ne.0.0d0)then
                  yscratch(ibet)=flum(iidx)/visib(iidx)
                else
                  yscratch(ibet)=0.0d0
                endif
 190          continue
c
              xscratch(ibetlim(ialf)+1)=xscratch(1)
              yscratch(ibetlim(ialf)+1)=yscratch(1)
c
              icount=0

             do 200 i=1,nbin

              do 191 ibet=1,ibetlim(ialf)
                if(((binx(i).gt.xscratch(ibet)).and.
     $              (binx(i).le.xscratch(ibet+1))).or.
     &              ((binx(i).le.xscratch(ibet)).and.
     $              (binx(i).gt.xscratch(ibet+1))))then
                  x0=xscratch(ibet)
                  x1=xscratch(ibet+1)
                  y0=yscratch(ibet)
                  y1=yscratch(ibet+1)
                  xx=binx(i)
c
                  yy=((x1-xx)*y0+(xx-x0)*y1)/(x1-x0)
                  weight=0.666666666666666666
                  if(mod(ialf,2).eq.0)weight=1.3333333333333333
                  if((ialf.eq.1).or.(ialf.eq.Nalf))weight=0.333333333333
                  biny(i)=biny(i)+weight*yy
                  icount=icount+1
                  jcount=jcount+1
                endif
 191          continue

 200          continue
 199        continue
c
          area=0.0d0
          do 300 i=1,nbin
            area=area+binsize*biny(i)
 300      continue
c
          solarrad=6.9598d10
          vel=overQ/(1.0d0+overQ)*(siphase+argfac)*sifinc
          vel=velamp*vel+gamma
          do 301 i=1,nbin
            biny(i)=biny(i)/flux*(separ*solarrad)**2
            write(98,500)binx(i),biny(i),vel
 301      continue
c
          close(98)
c
 500      format(f12.6,2x,1pe16.9,0pf12.6)
          return
          end
c
c  ##################################
c
c
          integer function kount(ialphmax,ialf,ibetlim)
c
c   June 11, 2003
c
c   this little function will keep track of the
c   index for the new 1D arrays.
c
          dimension ibetlim(ialphmax)
c
          kount=0
          if(ialf.le.1)then
            kount=0
            return
          endif
c
          kount=0
          do 10 i=2,ialf
            kount=kount+ibetlim(i-1)
 10       continue
c          
          return
          end
c
c

          subroutine getrefvisib(istar,ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     %      phase,finc,Q,psi0,omega,
     $      gradx,grady,gradz,visib,projarray,separ,bdist,mmdx)
c
c   June 16, 2003
c
c   This routine will return the 'reference visibilities' for use
c   in the EBOP mode.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension visib(ialphmax*ibetmax),
     $        gradx(ialphmax*ibetmax),grady(ialphmax*ibetmax),
     $        gradz(ialphmax*ibetmax),
     $        projarray(ialphmax*ibetmax),ibetlim(ialphmax),
     $        mmdx(ialphmax,ibetmax)
c
c  initialize the visibities!
c
          do 1 i=1,nalf
            do 2 j=1,ibetlim(i)  !4*nbet
c              iidx=(i-1)*ibetlim(i)+j
              iidx=mmdx(i,j)
              visib(iidx)=0.0d0
              projarray(iidx)=-1.0d0
 2          continue
 1        continue
c
c
c   RVG BUG ALERT   May 2, 2001
c
c   Change the definition of phaser to the simplified form below (i.e.
c   phaser is simply the phase in radians.
c
c          if(phase.gt.180.0d0)then
c            phaser=-(phase)*pie/180.0d0
c          else
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
c          endif
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c          
          NBET4 = NBET*4
          DBETA = (pie/2.0d0)/NBET        ! step size in longitude
c
          AZ = DCOS(FINCR)
          IF (AZ.LT.0.0d0) AZ = 0.0d0
          AX = -DSIN(FINCR)*DCOS(PHASER)    ! l in Wilson & Sofia
          AY = DSIN(FINCR)*DSIN(PHASER)     ! m in Wilson & Sofia
c
c   Check to see of the star in question is in front.  If so, then simply
c   find the projection factors.
c
c
          DO 501 IALF = 1, NALF
            DO 502 IBET = 1,ibetlim(ialf)      !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              PROJ = AX * GRADX(iidx) + AY*GRADY(iidx) + 
     1	        AZ*GRADZ(iidx)
              projarray(iidx)=proj
              visib(iidx)=proj
c              write(*,*)istar,ialf,ibet,proj
c
 502        continue
 501      continue
          return
          end
c
c  *************************************************************************
c
c  *************************************************************************
c
          subroutine getrefBBflux(ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     $      wave,visib,projarray,temp,surf,flimbx,flimby,ilaw,rinty,
     &      flum,flux,rldint,
     &      separ,mmdx)
c
c     June 16, 2003
c 
c     This routine will return the 'reference BB flux' for use in
c     the EBOP mode.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension visib(ialphmax*ibetmax),
     $        surf(ialphmax*ibetmax),ibetlim(ialphmax),
     $        temp(ialphmax*ibetmax),flum(ialphmax*ibetmax),
     $        rinty(ialphmax*ibetmax),projarray(ialphmax*ibetmax),
     %        mmdx(ialphmax,ibetmax)
c
          dint=pie*(1.0d0-flimbx/3.0d0)
          if(ilaw.eq.2)then
            dint=pie*(1.0d0-flimbx/3.0d0+2.0d0*flimby/9.0d0)
          endif
          if(ilaw.eq.3)then
            dint=pie*(1.0d0-flimbx/3.0d0-flimby/5.0d0)
          endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
          if(ilaw.eq.4)then
            dint=pie*(1.0d0-flimbx/3.0d0-flimby/6.0d0)
          endif

          wavemu=wave/10000.0d0
          C2 = 1.4384d8          ! 1.4384 * 10.**8      ! hc/(k*1e-8)
          C1 = 1.191044d35       ! 2hc^2/((1e-8)**5)
c
c   Initialize the flum matrix.
c
c
          c1=3.74185
          c2=14.3883

          flux=0.0d0
          do 2 ialf=1,nalf
            do 1 ibet=1,ibetlim(ialf)        !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              flum(iidx)=0.0d0
              rinty(iidx)=0.0d0
 1          continue
 2        continue
          sumcor1=0.0d0
          sumcor2=0.0d0
          dtheta=0.5d0*pie/dble(Nalf)
c
          DO 10 ialf=1,nalf
            theta=-dtheta+2.0d0*dtheta*dble(ialf)
            sitheta=dsin(theta)
            dphi=pie/dble(ibetlim(ialf))
            DO 9 ibet = 1,ibetlim(ialf)               !4*nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c 
              iidx=mmdx(ialf,ibet)
              corr1=0.0d0
              corr2=0.0d0
              if((projarray(iidx).le.0.0d0))go to 9
c              C3 = C2/(WAVE*TEMP(iidx))
c
              tkkelv=temp(iidx)/1000.0d0
              C3 = C2/(wavemu*tkkelv)
              flum(iidx)=C1/(dexp(c3)-1.0d0)/wavemu**5
              dark=(1.0d0-flimbx+flimbx*projarray(iidx))
              if(ilaw.eq.2)dark=dark-flimby*projarray(iidx)*
     %               dlog(projarray(iidx))
              if(ilaw.eq.3)dark=dark-flimby*(1.0-dsqrt(projarray(iidx)))
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
              if(ilaw.eq.4)dark=dark-flimby*(1.0-(projarray(iidx)))**2
c
              flum(iidx)=flum(iidx)*dark
              rinty(iidx)=flum(iidx) ! save intensities for plotting
              flum(iidx)=surf(iidx)*flum(iidx)*visib(iidx)
              flux=flux+flum(iidx)
 9          continue
 10       continue
c
c   Scale the light curve by the integral of the limb darkening law
c   for compatibility with Wilson-Devinney.
c        
c   Scale the light curve by the integral of the limb darkening law
c   for compatibility with Wilson-Devinney.
c        
          flux=pie*flux/dint
c
c
c   UPDATE April 3, 2002
c
c   Scale the fluxes.
c
          solarrad=6.9598d10
          flux=flux*(separ*solarrad)**2
c
          rldint=dint
          return
          end
c
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
          subroutine getrefvel(istar,omega,phase,finc,
     %      Q,separ,period,gamma,vel,ecc,argrad,isw12,gimvel,
     #                     iRVfilt)
c
c    November 12, 1999
c
c    This routine will compute the flux-weighted radial velocity of the
c    star in question.
c
c
          implicit double precision(a-h,o-z)
c
          dimension gimvel(8)

          parameter(pie=3.14159265358979323d0)    
c
          argfac=ecc*dcos(argrad)
          efact=1.0d0/dsqrt(1.0d0-ecc*ecc)
          dint=rldint/pie
          vel=0.0d0
          delvel=0.0d0
          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q
c
          PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
          FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c
          siphase=dsin(phaser)
          cophase=dcos(phaser)
          sifinc=dsin(fincr)
c
c   Compute the expected K velocity, which is the circular velocity times 
c   dsin(finc).
c
c    
          a=separ*6.9598d5            !separation in km
          p=period*24.00d0*3600.00d0       !period in seconds
          velamp=2.0d0*pie*a/p*efact
          ppp = (PHASE/180.0d0)*pie    
          siphase=dsin(ppp)
          vel=overQ/(1.0d0+overQ)*(siphase+argfac)*sifinc
          if(isw12.gt.0)then
            vel=vel+gimvel(iRVfilt)
          endif
          vel=vel*velamp+gamma
          delvel=0.0d0
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&
c  &&&&&&&&&&&&&&&&&
c
          subroutine getrefATMflux(ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     $      visib,projarray,temp,surf,garray,rinty,
     &      flum,maxlines,maxmu,Nlines,
     &      atmT,atmg,atmmu,Nmu,
     &      atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &      Tmax,Tmin,gmax,gmin,gscale,
     &      fluxU,fluxB,fluxV,fluxR,fluxI,fluxJ,fluxH,fluxK,
     $      icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,delphi,
     %      delphiedge,iedgestar,iedgehor,darkint,separ,mmdx,
     &      dwavex,dwavey,ilaw,iatm,istar)
c
c
c     June 16, 2003
c 
c     This routine will return the 'reference ATM flux' for use in
c     the EBOP mode.
c
          implicit double precision(a-h,o-z)

          parameter(pie=3.141592653589793d0)
c
c   Set these to the value of ialphmax,ibetmax
c
          integer tempalf,tempbet
          parameter(tempalf=3000,tempbet=3000,itab=tempalf*tempbet)

          dimension visib(ialphmax*ibetmax),ibetlim(ialphmax),
     $        surf(ialphmax*ibetmax),garray(ialphmax*ibetmax),
     $        temp(ialphmax*ibetmax),flum(ialphmax*ibetmax),
     $        rinty(ialphmax*ibetmax),projarray(ialphmax*ibetmax),
     $        delphi(ialphmax*ibetmax),delphiedge(ialphmax*ibetmax),
     $        iedgestar(ialphmax*ibetmax),iedgehor(ialphmax*ibetmax),
     &        mmdx(ialphmax,ibetmax)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines),outinty(8),
     #       corr1(8),corr2(8),saveflum(itab,8),darkint(8)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
          dimension dwavex(8,2),dwavey(8,2)
c  
c
          if(tempalf.lt.ialphmax)then
            write(*,*)'dimension error in getATMflux'
            stop
          endif
          if(tempbet.lt.ibetmax)then
            write(*,*)'dimension error in getATMflux'
            stop
          endif
c
          fluxU=0.0d0
          fluxB=0.0d0
          fluxV=0.0d0
          fluxR=0.0d0
          fluxI=0.0d0
          fluxJ=0.0d0
          fluxH=0.0d0
          fluxK=0.0d0
          corr1(1)=0.0d0
          corr1(2)=0.0d0
          corr1(3)=0.0d0
          corr1(4)=0.0d0
          corr1(5)=0.0d0
          corr1(6)=0.0d0
          corr1(7)=0.0d0
          corr1(8)=0.0d0
          corr2(1)=0.0d0
          corr2(2)=0.0d0
          corr2(3)=0.0d0
          corr2(4)=0.0d0
          corr2(5)=0.0d0
          corr2(6)=0.0d0
          corr2(7)=0.0d0
          corr2(8)=0.0d0
c
c   Initialize the flum matrix.
c
          do 2 ialf=1,nalf
            do 1 ibet=1,ibetlim(ialf)       !4*Nbet
c              iidx=(ialf-1)*4*Nbet+ibet
              iidx=mmdx(ialf,ibet)
              flum(iidx)=0.0d0
              rinty(iidx)=0.0d0
              saveflum(iidx,1)=0.0d0
              saveflum(iidx,2)=0.0d0
              saveflum(iidx,3)=0.0d0
              saveflum(iidx,4)=0.0d0
              saveflum(iidx,5)=0.0d0
              saveflum(iidx,6)=0.0d0
              saveflum(iidx,7)=0.0d0
              saveflum(iidx,8)=0.0d0
 1          continue
 2        continue
c
          Nalf2=Nalf/2
c
c   Find the rough place in the atmosphere table.
c
          Tin=temp(1)
          call locate(atmT,Nlines,Tin,indexT)
          itguess=indexT
          imuguess=1
c
c
          DO 10 ialf=1,nalf
            DO 9 ibet = 1,ibetlim(ialf)      !4*nbet
              iidx=mmdx(ialf,ibet)
              dphi=pie/dble(ibetlim(ialf))
              if(projarray(iidx).le.0.0d0)go to 9
              Tin=temp(iidx)
              gin=dlog10(gscale*garray(iidx))
              rmuin=projarray(iidx)
              call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &           atmint7,atmint8,
     &           Tmax,Tmin,gmax,gmin,outinty,
     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     #           dwavex,dwavey,ilaw,iatm,istar)
c
              rinty(iidx)=outinty(iRVfilt) ! save intensities for plotting
c
              do 8 k=1,8
                saveflum(iidx,k)=
     #             outinty(k)*surf(iidx)*projarray(iidx)
                outinty(k)=outinty(k)*surf(iidx)*visib(iidx)
 8            continue
c
              fluxU=fluxU+outinty(1)
              fluxB=fluxB+outinty(2)
              fluxV=fluxV+outinty(3)
              fluxR=fluxR+outinty(4)
              fluxI=fluxI+outinty(5)
              fluxJ=fluxJ+outinty(6)
              fluxH=fluxH+outinty(7)
              fluxK=fluxK+outinty(8)
              izz=ialf
              jzz=ibet
              iidx=mmdx(izz,jzz)
              flum(iidx)=outinty(iRVfilt)+corr1(iRVfilt)+corr2(iRVfilt) 
 9          continue
 10       continue
c
c          if(darkint(1).ne.0.0d0)fluxU=pie*fluxU!/darkint(1)
c          if(darkint(2).ne.0.0d0)fluxB=pie*fluxB!/darkint(2)
c          if(darkint(3).ne.0.0d0)fluxV=pie*fluxV!/darkint(3)
c          if(darkint(4).ne.0.0d0)fluxR=pie*fluxR!/darkint(4)
c          if(darkint(5).ne.0.0d0)fluxI=pie*fluxI!/darkint(5)
c          if(darkint(6).ne.0.0d0)fluxJ=pie*fluxJ!/darkint(6)
c          if(darkint(7).ne.0.0d0)fluxH=pie*fluxH!/darkint(7)
c          if(darkint(8).ne.0.0d0)fluxK=pie*fluxK!/darkint(8)
c
c  UPDATE April 3, 2002
c
c  Scale the fluxes.
c
          solarrad=6.9598d10
          fluxU=fluxU*(separ*solarrad)**2
          fluxB=fluxB*(separ*solarrad)**2
          fluxV=fluxV*(separ*solarrad)**2
          fluxR=fluxR*(separ*solarrad)**2
          fluxI=fluxI*(separ*solarrad)**2
          fluxJ=fluxJ*(separ*solarrad)**2
          fluxH=fluxH*(separ*solarrad)**2
          fluxK=fluxK*(separ*solarrad)**2
c
          return
          end
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine overlaphoriz(Nhoriz1,xhoriz1,yhoriz1,
     $            Nhoriz2,xhoriz2,yhoriz2,ioverlap)
c
c   UPDATE May 26, 2004
c
c   This routine will take the horizons for star 1 and star 2 and
c   see if there is an overlap on the sky.  If so, then ioverlap=999
c         
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension xhoriz1(Nhoriz1),yhoriz1(Nhoriz1)
          dimension xhoriz2(Nhoriz2),yhoriz2(Nhoriz2)
c
          do 10 i=1,Nhoriz1
            iyes=-1
            xp=xhoriz1(i)
            yp=yhoriz1(i)
            call insidecircle(Nhoriz2,xhoriz2,yhoriz2,xp,yp,iyes,icut)
            if(iyes.eq.100)ioverlap=999
 10       continue
c
          do 20 i=1,Nhoriz2
            iyes=-1
            xp=xhoriz2(i)
            yp=yhoriz2(i)
            call insidecircle(Nhoriz1,xhoriz1,yhoriz1,xp,yp,iyes,icut)
            if(iyes.eq.100)ioverlap=999
 20       continue
c        
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
c
      FUNCTION ran9(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      double precision ran9,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     $     NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.e-16,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran9=min(AM*iy,RNMX)
      return
      END
c
c  *****************************************
c
          subroutine sobseq(n,x)
          integer n,MAXBIT,MAXDIM
          double precision x
          dimension x(*)
          parameter(MAXBIT=30,MAXDIM=6)
          integer i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT),
     $       iv(MAXBIT*MAXDIM),ix(MAXDIM),mdeg(MAXDIM)
          double precision fac
          SAVE ip,mdeg,ix,iv,in,fac
          equivalence (iv,iu)
          data ip /0,1,1,2,1,4/, mdeg/1,2,3,3,4,4/, ix/6*0/
          data iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
c
c
          if(n.eq.-1)then
            in=0
c
            do 55 ll=1,MAXDIM
              do 54 mm=1,MAXBIT
                iu(ll,mm)=0
 54                       continue
 55                                continue
c
            do 99 kkk=1,6
              ix(kkk)=0
              iv(kkk)=1
 99                    continue

           iv(7)=3
           iv(8)=1
           iv(9)=3
           iv(10)=3
           iv(11)=1
           iv(12)=1
           iv(13)=5
           iv(14)=7
           iv(15)=7
           iv(16)=3
           iv(17)=3
           iv(18)=5
           iv(19)=15
           iv(20)=11
           iv(21)=5
           iv(22)=15
           iv(23)=13
           iv(24)=9
c
           do 66 kkk=0,155
             iv(kkk+25)=0
 66                   continue
c
          endif
c
          if(n.lt.0)then
            do 14  k=1,MAXDIM
              do 11 j=1,mdeg(k) 
                iu(k,j)=iu(k,j)*2**(MAXBIT-j)
 11           continue
              do 13 j=mdeg(k)+1,MAXBIT
                ipp=ip(k) 
                i=iu(k,j-mdeg(k))
                i=ieor(i,i/2**mdeg(k))
                do 12 l=mdeg(k)-1,1,-1
                  if(iand(ipp,1).ne.0)i=ieor(i,iu(k,j-l))
                  ipp=ipp/2
 12             continue
                iu(k,j)=i
 13           continue
 14         continue
            fac=1.0d0/2.0d0**MAXBIT
            in=0
          else
            im=in
            do 15 j=1,MAXBIT 
              if(iand(im,1).eq.0)go to 1
              im=im/2
 15         continue
            pause 'MAXBIT too small in sobseq'
 1          im=(j-1)*MAXDIM
            do 16 k=1,min(n,MAXDIM)
              ix(k)=ieor(ix(k),iv(im+k))
              x(k)=ix(k)*fac
 16         continue
            in=in+1
          endif
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine setscale(fill1,fill2,omega1,omega2,
     &       Q,finc,Teff1,Teff2,Period,fm,separ,
     &       primmass,primK,primrad,ratrad,sw5,reff1,ecc,bdist)
c
c   August 10, 2004
c
c   This a new subroutine to set constraints specified by the
c   new variables primmass, primK, etc.
c
          implicit double precision (a-h,o-z)  

c   We can set most of the constraints specified by primmas, primK,
c   primrad here.  ratrad needs to be dealt with in getinput.
c
c   There are various combinations of parameters that one can adjust,
c   and depending on which constraints are specified, different
c   ways are used:
c
c   if sw5 > 0, then bail out (MSP mode).
c
c   primmass > 0 only:  set the separation.
c   primK > 0 only:     compute f(M), use getradius() to set separation.
c   primrad > 0 only:   scale the separation
c
c
c
          parameter(pie=3.14159265358979323d0)
c
           solarmass=1.9889d33
           smmks=1.9889d30    !solar mass in kg
           solarrad=6.9598d10
           Gcgs=6.67259d-8     !G in cgs
           Gmks=6.67259d-11   !G in mks
           p=period*24.0d0*3600.0d0
           gmsun=1.32712440018d20  !mks units
c
c 
c   MSP mode:
c
c   Add the variable sw5 to the argument list.  If teff2 < 0 and
c   sw5 > 0, then the separation will be set as follows:
c
c   rkns = 2*pie*sw5*c/Period   ! K-velocity of pulsar, if sw5 is
c                                  projected semimajor axis in seconds
c
           if((teff2.le.0.0d0).and.(sw5.gt.0.0d0))then
             speedlight=2.997924548d5
             rkns=2.0d0*pie*sw5*speedlight/Period
             efact=dsqrt(1.0d0-ecc*ecc)
c
c   UPDATE November 14, 2008
c
c   If primK >0, set the mass ratio
c
             if(primK.gt.0.0d0)Q=primK*86400.0d0/(rkns*efact)
             vkcgs=primK*100000.0d0*efact
             fincr=finc*pie/180.0d0
c             write(*,*)p,q
             separ=vkcgs*p*(1.0d0+Q)/(2.0*pie*dsin(fincr)*Q)/solarrad
             separ=(Q+1.0d0)*sw5*speedlight/(dsin(fincr))/solarrad
             separ=separ*1.0d5
             write(2,199)separ
             if(primK.gt.0.0d0)write(2,1999)primK,Q
c             write(*,*)separ,vkcgs,finc,fincr
             return
           endif
c
  199      format('Info:  MSP mode:  separation has been set to ',
     %         f13.7,7x,' solar radii')
 1999      format('Info:  MSP mode:  primK = ',f13.7,' the mass ratio has'/,
     #         6x,' been set to',f13.7)
c
c   primmass > 0 only
c
           if((primmass.gt.0.0d0).and.(primK.le.0.0d0).and.
     $       (primrad.le.0.0d0))then
             rmass=primmass*solarmass
             separ=(Gcgs*p*p*(1.0d0+Q)*rmass
     $          /(4.0d0*pie*pie))**(1.d0/3.0d0)/solarrad
c
             separ=(gmsun*p*p*primmass*(1.0d0+Q)/(4.0d0*pie*pie))**(1.0d0/3.0d0)
             separ=separ/solarrad*100.0d0
             write(2,1199)primmass,separ
             return
           endif
c
c
c
c   primmass > 0 and primrad > 0
c
           if((primmass.gt.0.0d0).and.(primK.le.0.0d0).and.
     $       (primrad.gt.0.0d0))then
             rmass=primmass*solarmass
             separ=(Gcgs*p*p*(1.0d0+Q)*rmass
     $          /(4.0d0*pie*pie))**(1.d0/3.0d0)/solarrad
c
             separ=(gmsun*p*p*primmass*(1.0d0+Q)/(4.0d0*pie*pie))**(1.0d0/3.0d0)
             separ=separ/solarrad*100.0d0
             write(2,1199)primmass,separ
             return
           endif
c
c   primK > 0 only
c
           if((primmass.le.0.0d0).and.(primK.gt.0.0d0).and.
     $       (primrad.le.0.0d0))then
c
c   UPDATE NOVEMBER 1
c
c   add the efact  =  dsqrt(1-e*e)
c
             efact=dsqrt(1.0d0-ecc*ecc)
             vkcgs=primK*100000.0d0*efact
             fincr=finc*pie/180.0d0
             separ=vkcgs*p*(1.0d0+Q)/(2.0*pie*dsin(fincr)*Q)/solarrad
c             fm=p*(vkmks**3.0)/(2.0d0*pie*Gmks)/smmks
c             call getradius(Q,finc,rad_in_cm,fm,period,ecc)
c             separ=rad_in_cm/solarrad
             write(2,1188)primK,separ
             return
           endif
c
c
c   primK > 0 and primrad > 0
c
           if((primmass.le.0.0d0).and.(primK.gt.0.0d0).and.
     $       (primrad.gt.0.0d0))then
c
c   UPDATE NOVEMBER 1
c
c   add the efact  =  dsqrt(1-e*e)
c
             efact=dsqrt(1.0d0-ecc*ecc)
             vkcgs=primK*100000.0d0*efact
             fincr=finc*pie/180.0d0
             separ=vkcgs*p*(1.0d0+Q)/(2.0*pie*dsin(fincr)*Q)/solarrad
c             fm=p*(vkmks**3.0)/(2.0d0*pie*Gmks)/smmks
c             call getradius(Q,finc,rad_in_cm,fm,period,ecc)
c             separ=rad_in_cm/solarrad
             write(2,1188)primK,separ
             return
           endif
c

c   primrad > 0 only
c
c           if((primmass.le.0.0d0).and.(primK.le.0.0d0).and.
c     $       (primrad.gt.0.0d0))then
c               separ=primrad/reff1
c               return
c           endif
c
c   primmass > 0 and primK > 0.  Solve for Q and separ
c
           if((primmass.gt.0.0d0).and.(primK.gt.0.0d0))then
             rmass=primmass*solarmass
c
c
c   UPDATE NOVEMBER 1
c
c   add the efact  =  dsqrt(1-e*e)
c
             efact=dsqrt(1.0d0-ecc*ecc)
c
             dqhi=7.
             dqlo=-7.
             do 555 kk=1,35
               Qhigh=10.0d0**dqhi
               Qlow=10.0d0**dqlo
               Qmid=10.0**((dqhi+dqlo)*0.5d0)
               aa=vfcn(Qlow,period,finc,primmass,primK,ecc)
               bb=vfcn(Qhigh,period,finc,primmass,primK,ecc)
               cc=vfcn(Qmid,period,finc,primmass,primK,ecc)
               if(aa*cc.lt.0.0d0)then
                 dqhi=(dqhi+dqlo)*0.5d0
               else
                 dqlo=(dqhi+dqlo)*0.50
               endif
 555           continue
c
             do 556 kk=1,25
               Qmid=(Qhigh+Qlow)*0.5d0
               aa=vfcn(Qlow,period,finc,primmass,primK,ecc)
               bb=vfcn(Qhigh,period,finc,primmass,primK,ecc)
               cc=vfcn(Qmid,period,finc,primmass,primK,ecc)
               if(aa*cc.lt.0.0d0)then
                 Qhigh=(Qhigh+Qlow)*0.5d0
               else
                 Qlow=(Qhigh+Qlow)*0.50
               endif
 556          continue

              Q=Qmid
c              separ=(Gcgs*p*p*(1.0d0+Q)*rmass/(4.0d0*pie*pie))
c              separ1=separ**(1.0d0/3.0d0)/solarrad
              separ=(gmsun*p*p*primmass*
     &              (1.0d0+Q)/(4.0d0*pie*pie))**(1.0d0/3.0d0)
              separ=separ/solarrad*100.0d0              
c
c
c   Use the formula separ = coef*(perid*period*total_mass)**(1/3) to
c   solve for the total mass in solar masses.  The separation is
c   entered in solar masses, so (R_sun/coef)**3=7.737294491.
c
              total_mass=primmass*(1.0d0+Q)
c
              ppp=period*24.0d0
              separ1=(total_mass*ppp*ppp/7.737294491d0)**(1.0d0/3.0d0)
c
c           total_mass=(separ)**(3)*7.737294491d0/(ppp*ppp)
c
             write(2,1177)primmass,primK,Q,separ
             

             return
           endif
c
 1177      format(/'Info:  M_1 fixed at ',f9.6,' solar masses and ',
     &        'K_1 fixed at ',f11.6,' km/sec.',/
     &        'Q is set to ',f12.7,' and the separation is ',f13.7,
     &        ' solar radii')
1199       format(/'Info:  M_1 fixed at ',f9.6,' solar masses.  The',  
     %        ' separation',/ '       has been set to ',f13.7,
     $         ' solar radii')
 1188      format(/'Info:  K_1 fixed at ',f14.8,' km/sec.  The',  
     %        ' separation',/ '       has been set to ',f13.7,
     $         ' solar radii')  
c
            return
            end
c
c   *******
c
            function vfcn(Q,period,finc,primmass,primK,ecc)
c
c
c
c
c   UPDATE NOVEMBER 1
c
c   add the efact  =  dsqrt(1-e*e)
c
            implicit double precision (a-h,o-z)
c
           parameter(pie=3.14159265358979323d0)
c
             efact=dsqrt(1.0d0-ecc*ecc)
c
           solarmass=1.9889d33
           smmks=1.9889d30    !solar mass in kg
           solarrad=6.9598d10
           Gcgs=6.67259d-8     !G in cgs
           Gmks=6.67259d-11   !G in mks
           p=period*24.0d0*3600.0d0
           fincr=finc*pie/180.0d0
c
           gmsun=1.32712440018d20  !mks units
c
           tt1=2.0d0*pie*dsin(fincr)/p*(Q/(Q+1.0d0))
           tt2=gmsun*p*p*(1.0d0+Q)*primmass/(4.0d0*pie*pie)
c
           total_mass=primmass*(1.0d0+Q)
c
           ppp=period*24.0d0
           separ=(total_mass*ppp*ppp/7.737294491d0)**(1.0d0/3.0d0)*solarrad
           tt2=tt2**(1.0d0/3.0d0)*100.0d0
c
           vfcn=tt1*tt2-primK*100000.0d0*efact
c
           return
           end
c
c
c
           subroutine findradius(overQ,omega,psi0,x0,bdist,reff,
     #        tidephi,itide)
c
c   This will find the effective radius.
c
           implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
c
           twopie=2.0d0*pie
           Na=120
           Nb=40
           dtheta=pie/dble(Na)
           dcostheta=2.0d0/dble(Na)
 700       vol=0.0d0
           sarea=0.0d0
           potsum=0.0d0
c
           if(itide.lt.2)then
           do 104 ialf=1,Na/2
c
c   UPDATE MARCH 17, 2004
c
c   make the initial value of r smaller here.
c
             r=1.0d-15
             x=x0
             y=0.0d0
             z=0.0d0
             call rad(overQ,omega,0.0d0,0.0d0,1.0d0,psi0,r,x,y,z,1,bdist,
     $          tidephi,itide)
             theta=-0.5d0*dtheta+dtheta*dble(ialf)
             dphi=twopie/dble(4*Nb)
             snth=dsin(theta)
             snth3=snth*(1.0d0/3.0d0)
             cnth=dcos(theta)
             DO 105 ibet=1,Nb*2         !4*Nbet
               phi=-0.5d0*dphi+dphi*dble(ibet)
c
c              if(dble(ialf/2).eq.(dble(ialf)*0.5d0))phi=phi+0.25d0*dphi   !DPHI
c
               cox=dcos(phi)*snth             !*dsin(theta)
               coy=dsin(phi)*snth             !*dsin(theta)
               coz=cnth                       !dcos(theta)
               CALL RAD(overQ,omega,cox,coy,coz,psi0,r,x,y,z,1,bdist,
     $            tidephi,itide)
               VOL = VOL + 4.0d0*R*R*R*dphi*dtheta*snth3
c               write(*,*)'ialf = ',ialf,' ibet = ',ibet,'  ',r
105          CONTINUE    ! continue ibet loop
104        CONTINUE                   ! continue over ialf
c
           endif   !end if itide < 2
c
           if(itide.ge.2)then
           do 304 ialf=1,Na
c
c   UPDATE MARCH 17, 2004
c
c   make the initial value of r smaller here.
c
             r=1.0d-15
             x=x0
             y=0.0d0
             z=0.0d0
             call rad(overQ,omega,0.0d0,0.0d0,1.0d0,psi0,r,x,y,z,1,bdist,
     $          tidephi,itide)
             theta=-0.5d0*dtheta+dtheta*dble(ialf)
             dphi=twopie/dble(4*Nb)
             snth=dsin(theta)
             snth3=snth*(1.0d0/3.0d0)
             cnth=dcos(theta)
             DO 305 ibet=1,Nb*4         !4*Nbet
               phi=-0.5d0*dphi+dphi*dble(ibet)
c
c              if(dble(ialf/2).eq.(dble(ialf)*0.5d0))phi=phi+0.25d0*dphi  !DPHI
c
               cox=dcos(phi)*snth             !*dsin(theta)
               coy=dsin(phi)*snth             !*dsin(theta)
               coz=cnth                       !dcos(theta)
               CALL RAD(overQ,omega,cox,coy,coz,psi0,r,x,y,z,1,bdist,
     $            tidephi,itide)
               VOL = VOL + 1.0d0*R*R*R*dphi*dtheta*snth3
c               write(*,*)'ialf = ',ialf,' ibet = ',ibet,'  ',r
305          CONTINUE    ! continue ibet loop
304        CONTINUE                   ! continue over ialf
c
           endif   !end if itide >= 2

c
c
c           REFF = (0.248732415d0*VOL) **(1.0d0/3.0d0) 
c
           REFF = (0.75d0*VOL/pie) **(1.0d0/3.0d0) 
c
           return
           end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
           subroutine filltime(ntime,timearray,tstart,tstop,tstep)
c
           implicit double precision (a-h,o-z)
c
           dimension timearray(900000)
c
           ntime=0
           do 10 i=1,900000
             tt=tstart+dble(i-1)*tstep
             if(tt.gt.tstop) go to 15
             ntime=ntime+1
             timearray(ntime)=tt
 10        continue
 15        return
c
           end
c
c  &&&&&&&&&&&&&&&&
c
c
          subroutine addmovespot(istar,ialphmax,ibetmax,
     $       Nalph,ibetlim,tmatrix,spotparm,ave1,ave2,omega,
     %       phiar,mmdx,period,t0,ttime)
c
c   This routine will assign the temperatures of the grid points
c   of the stars that are covered by spots.  The underlying temperatures
c   are simply scaled by the temperature spot factor.  
c   
c   UPDATE March 26, 2002
c
c   Get rid of Nbet from the argument list of addstarspot.  Its value
c   is contained within ibetlim.
c
c    UPDATE June 17, 2002
c
c    Add these dummy assignments for phase and omega to supress
c    compiler warnings about unused variables.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265359879323d0)
          dimension tmatrix(ialphmax*ibetmax),phiar(ialphmax*ibetmax),
     %       ibetlim(ialphmax),spotparm(2,4),mmdx(ialphmax,ibetmax)
c
c
c    UPDATE June 17, 2002
c
c    Add these dummy assignments for phase and omega to supress
c    compiler warnings about unused variables.
c
          ddd1=phase
          ddd2=omega
c
          radcon=pie/180.0d0
          halfpie=0.5*pie
c
          poff=(ttime-T0)/period*(omega-1.0d0)*2.0d0*pie
          fac1=spotparm(1,1)
          fac2=spotparm(2,1)
          rlat1=radcon*spotparm(1,2)-halfpie
          rlat2=radcon*spotparm(2,2)-halfpie
          rlong1=dmod(radcon*spotparm(1,3)+poff,2.0d0*pie)
          rlong2=dmod(radcon*spotparm(2,3)+poff,2.0d0*pie)
          if(rlong1.gt.pie)rlong1=rlong1-2.0d0*pie
          if(rlong2.gt.pie)rlong2=rlong2-2.0d0*pie
          rad1=radcon*spotparm(1,4)
          rad2=radcon*spotparm(2,4)
c
          if((fac1.lt.0.0d0).and.(fac2.lt.0.0d0))return !no valid factors
c
          dtheta=pie/dble(Nalph)
          summ1=0.0d0
          summ2=0.0d0
          icount1=0
          icount2=0
          DO 10 IALF = 1, nalph
            theta=-0.5d0*dtheta+dtheta*dble(ialf)
            rlat=theta-halfpie
            DO 9 IBET = 1, ibetlim(ialf)    !4*NBET
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=(ialf-1)*ibetlim(ialf)+ibet
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              rlong=phiar(iidx)
              if(rlong.gt.pie)rlong=rlong-2.0d0*pie
c
              sepang1=dsin(rlat)*dsin(rlat1)+dcos(rlat)*dcos(rlat1)*
     %           dcos(rlong-rlong1)
              sepang1=dabs(dacos(sepang1))
              if((sepang1.le.rad1).and.(fac1.gt.0.0d0))then
                icount1=icount1+1
                tmatrix(iidx)=tmatrix(iidx)*fac1
                summ1=summ1+tmatrix(iidx)
              endif
c
              sepang2=dsin(rlat)*dsin(rlat2)+dcos(rlat)*dcos(rlat2)*
     %           dcos(rlong-rlong2)
              sepang2=dabs(dacos(sepang2))
              if((sepang2.le.rad2).and.(fac2.gt.0.0d0))then
                icount2=icount2+1
                tmatrix(iidx)=tmatrix(iidx)*fac2
                summ2=summ2+tmatrix(iidx)
              endif

 9          continue
 10       continue
c
          ave1=0.0d0
          ave2=0.0d0
          if(icount1.gt.0)ave1=summ1/dble(icount1)
          if(icount2.gt.0)ave2=summ2/dble(icount2)

c          if(ave1.gt.0.0d0)write(2,100)istar,ave1,icount1
c          if(ave2.gt.0.0d0)write(2,101)istar,ave2,icount2

c
 100      format(/'star ',i1,', spot 1:    average temperature ',
     &       f9.3,',  number of grid points = ',i4)
 101      format(/'star ',i1,', spot 2:    average temperature ',
     &       f9.3,',  number of grid points = ',i4)
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine getom(ecc,ecosw,argper)
c
c   July 29, 2005
c
c   This routine will return the value of omega when given the eccentricity
c   and the phase difference between secondary and primary eclipse.
c
c   It is assumed ecosw is in phase units, so multiply by 2pi to get radians
c   The value of omega will be returned in degrees.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265359879323d0)
c
          if(ecosw.eq.0.0d0)return
          if(ecc.le.0.0d0)return
          if(ecc.ge.1.0d0)return
c
          ppp=ecosw*2.0d0*pie
          call getX(ppp,ecc,xxx)
c
          ttt=dtan((xxx-pie)*0.5d0)
          www=ttt*dsqrt(1.0d0-ecc*ecc)/ecc
          if(www.lt.-1.0d0)www=-1.0d0
          if(www.gt.1.0d0)www=1.0d0
          argrad=dacos(www)
          argper=argrad*180.0d0/pie
          return
          end
c
c
c
          subroutine getX(em,ecc,bigE)
c
c    Solve for bigE in  em = bigE - sin(bigE)
c
          implicit double precision (a-h,o-z)
          Eold=em
c
          do 10 i=1,20
            top=Eold-dsin(Eold)-em
            bottom=1.0d0-dcos(Eold)
            Enew=Eold-top/bottom
            diff=dabs(Enew-Eold)
            if(diff.lt.1.0d-15)go to 15
            Eold=Enew
 10       continue
c
 15       bigE=Eold
c
          return
c
          end
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine cliplc(Nmaxphase,icount,dphase,Nphase,xmod,ymodU,
     %         ymodB,ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,
     &         ymods2,ymods3,ymodd,RV1,dRV1,RV2,dRV2)

c
          implicit double precision(a-h,o-z)
c
          dimension xmod(Nmaxphase),ymodU(Nmaxphase)
          dimension ymodB(Nmaxphase),ymodV(Nmaxphase)
          dimension ymodR(Nmaxphase),ymodI(Nmaxphase)
          dimension ymodJ(Nmaxphase),ymodH(Nmaxphase)
          dimension ymodK(Nmaxphase),ymods1(Nmaxphase)
          dimension ymods2(Nmaxphase),ymods3(Nmaxphase)
          dimension ymodd(Nmaxphase),RV1(Nmaxphase)
          dimension dRV1(Nmaxphase),RV2(Nmaxphase),dRV2(Nmaxphase)
c
          jcount=0
          do 10 i=1,icount
            if(ymods1(i).ge.0.0d0)then
              jcount=jcount+1
              xmod(jcount)=xmod(i)
              ymodU(jcount)=ymodU(i)
              ymodB(jcount)=ymodB(i)
              ymodV(jcount)=ymodV(i)
              ymodR(jcount)=ymodR(i)
              ymodI(jcount)=ymodI(i)
              ymodJ(jcount)=ymodJ(i)
              ymodH(jcount)=ymodH(i)
              ymodK(jcount)=ymodK(i)
              ymods1(jcount)=ymods1(i)
              ymods2(jcount)=ymods2(i)
              ymods3(jcount)=ymods3(i)
              ymodd(jcount)=ymodd(i)
c
c   UPDATE NOVEMBER 22, 2006
c 
c   Don't clip the velocity curves
c
c              RV1(jcount)=RV1(i)
c              RV2(jcount)=RV2(i)
c              dRV1(jcount)=dRV1(i)
c              dRV2(jcount)=dRV2(i)
            endif
10        continue
c
          Nphase=jcount
c
          return
          end
c 
c   Add subroutines to interpolate limb darkening coefficients for the disk.
c
          double precision function rlambdae(pee,zee)
c
c     April 24, 2006
c
c     Gives light curve of the transit of a uniform disk across a uniform
c     star.  pee is the size ratio (radius of planet divided by radius of
c     star). zee is the normalized separation of centers (d/radius of star)
c

          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265359879323d0)
c
          rlambdae=0.0d0
c
          if((1.0d0+pee).lt.zee)then
            rlambdae=0.0d0
c            write(*,99)rlambdae,pee,zee
99          format('block 1',2x,2(f9.5,2x))
            return
          endif
c
          z1=dabs(1.0d0-pee)
          z2=1.0d0+pee
          if((z1.lt.zee).and.(zee.le.z2))then
            rappa1=dacos(min((1.0-pee*pee+zee*zee)/(2.0d0*zee),1.0d0))
            rappa0=dacos(min((pee*pee+zee*zee-1.0d0)/(2.0d0*zee*pee),1.0d0))
            t1=0.5d0*dsqrt((4.0d0*zee*zee-(1.0d0+zee*zee-pee*pee)**2))
            rlambdae=(pee*pee*rappa0+rappa1-t1)/pie
c            write(*,199)rlambdae,pee,zee
199         format('block 2',3x,3(f9.4,3x))
            return
          endif
c
          if(zee.le.(1.0d0-pee))then
            rlambdae=pee*pee
            return
          endif
c
          if(zee.le.(pee-1.0d0))then
            rlambdae=1.0d0
            return
          endif
c
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine analytic(isw12,ilaw,dwavex,dwavey,pee,zee,
     #       refflux1,refflux2)
c
c   April 24, 2006
c
c   Will compute analytic transits (Mandel & Agol)
c
c 
          implicit double precision(a-h,o-z)
c 
          dimension dwavex(8,2),dwavey(8,2),refflux1(8),refflux2(8)
c
          rlam=rlambdae(pee,zee)
c          write(*,100)rlam
100       format('rlam = ',5x,3(f15.6,2x),1pe16.8)
c
          do 10 i=1,8
            refflux2(i)=-refflux1(i)*(rlam)
10        continue
c
          return
          end
c
c
c
          subroutine analyticg(isw12,ilaw,dwavex,dwavey,delta,reff1,ratrad,
     #       refflux1,refflux2,phaser,pconj,pconj2,gimvel,fincr,vrot,omega,
     &       period,separ,Q,ecc,bigI,bigbeta,Neclipse,istar,corr1,corr2)
c
c   April 24, 2006
c
c   Will compute analytic transits according to Gimenez.  The number isw12 is
c   the number of terms in the summation.
c
c 
          implicit double precision(a-h,o-z)
c 
          dimension dwavex(8,2),dwavey(8,2),refflux1(8),refflux2(8),cx1(10000)
          dimension corr1(8),corr2(8)
          dimension cx2(10000),c(10),alpha(10),cv1(10000),cv2(10000),alfone(10)
          dimension gimvel(8)
c
          parameter(pie=3.14159265359879323d0)
c
c
c

           Neclipse=99
           solarrad=6.9598d8
           p=period*86400.0d0
           fact=solarrad*2.0d0*pie/86400.0d0/1.0d3  !50.613093d0
           a2=(separ/(1.0d0+Q))
           a1=separ-a2
           efact=1.0d0/dsqrt(1.0-ecc*ecc)
           sifinc=dsin(fincr)
           velK2=fact*a2/period*sifinc*efact
c
c   add rotational velocity corrections
c
           
           sI=dsin(bigI*pie/180.0d0)
           sB=dsin(bigbeta*pie/180.0d0)
           cB=dcos(bigbeta*pie/180.0d0)
           vstar=-vrot*sI*(sB*dcos(fincr)*dcos(phaser)-cB*dsin(phaser))*omega
c
c   Do the case for linear first (ilaw=1)
c
          alpha(1)=0.0d0
          alpha(2)=0.0d0
          alpha(3)=0.0d0
          if(istar.eq.1)then
            rstar=reff1
            rplanet=rstar/ratrad
          endif
          if(istar.eq.2)then
            rplanet=reff1
            rstar=rplanet/ratrad

            rstar=reff1
            rplanet=rstar/ratrad
          endif
          bee=rplanet/(rstar+rplanet)
          cee=delta/(rplanet+rstar)
          if(delta.gt.(rstar+rplanet))then
            do 99 i=1,8
c              refflux2(i)=0.0d0
              corr1(i)=0.0d0
              corr2(i)=0.0d0
              gimvel(i)=0.0d0
99          continue
            Neclipse=0
            return
          endif

c   October 17, 2010
c
c   Update to make sure the nearest conjunction phase is found
c
          tt1=dabs(phaser-2.0d0*pie*pconj)
          tt2=dabs(phaser-2.0d0*pie*(pconj+1.0d0))
          tt3=dabs(phaser-2.0d0*pie*(pconj-1.0d0))
c
          small=123456789.0d0
          if(tt1.le.small)small=tt1
          if(tt2.le.small)small=tt2
          if(tt3.le.small)small=tt3
          diff1=small
c
c          tt1=dabs(phaser-2.0d0*pie*pconj)
c          tt2=dabs(phaser-2.0d0*pie*(pconj+1.0d0))
c          if(tt1.le.tt2)then
c            diff1=tt1
c          else
c            diff1=tt2
c          endif

          tt1=dabs(phaser-2.0d0*pie*pconj2)
          tt2=dabs(phaser-2.0d0*pie*(pconj2+1.0d0))
          tt3=dabs(phaser-2.0d0*pie*(pconj2-1.0d0))
c
          small=123456789.0d0
          if(tt1.le.small)small=tt1
          if(tt2.le.small)small=tt2
          if(tt3.le.small)small=tt3
          diff2=small

c          if(tt1.le.tt2)then
c            diff2=tt1
c          else
c            diff2=tt2
c          endif              

          if(istar.eq.1)then
            if(diff1.lt.diff2)then
              do 999 i=1,8
                corr1(i)=0.0d0   !was refflux2
                gimvel(i)=0.0d0
999            continue
              return
            endif
          endif

          if(istar.eq.2)then
            if(diff2.lt.diff1)then
              do 9999 i=1,8
c                refflux2(i)=0.0d0    !was commented out
                corr2(i)=0.0d0    
                gimvel(i)=0.0d0
9999           continue
              return
            endif
          endif

          if(ilaw.eq.1)then
c
            do 2 n=0,1
              rnu=dble(n+2)*0.5d0
c
              pee=rnu+2.0d0
              cue=rnu+1.0d0
              alf=pee-cue
              beta=cue-1.0d0
              xx=1.0d0-(2.0d0*(1.0d0-bee))
c
              call jacobi_poly(isw12,beta,alf,xx,cx1)
              cue=1.0d0
              alf=pee-cue
              beta=cue-1.0d0
              xx=1.0d0-2.0d0*cee*cee
c
              call jacobi_poly(isw12,beta,alf,xx,cx2)
c
              pee=rnu+3.0d0
              cue=2.0d0
              alf=pee-cue
              beta1=cue-1.0d0
              xx=1.0d0-(2.0d0*(cee*cee))
              call jacobi_poly(isw12,beta1,alf,xx,cv2)

              pee=rnu+3.0d0
              cue=rnu+2.0d0
              alf=pee-cue
              beta=cue-1.0d0
              xx=1.0d0-(2.0d0*(1.0d0-bee))
              call jacobi_poly(isw12,beta,alf,xx,cv1)

              t1=bee*bee*((1.0d0-cee*cee)**(rnu+1.0d0))
     &             /(rnu*dexp(gamma_log(rnu+1.0d0)))
c
              v1=dexp(gamma_log(rnu)-2.0d0*gamma_log(rnu+2.0d0))
              v1=v1*cee*bee*bee*(1.0d0-bee)*(1.0d0-cee*cee)**(rnu+1.0d0)

              summ=0.0d0
              summv=0.0d0
              do 3 j=0,isw12
                dj=dble(j)

                t3=(gamma_log(dj+1.0d0)+gamma_log(rnu+1.0d0)-
     %              gamma_log(dj+rnu+1.0d0))
                t3=exp(t3)

                t2=dexp(gamma_log(rnu+dble(j+1))-gamma_log(dble(j+2)))
                summ=summ+(-1.0d0)**j*(2.0d0*dble(j)+rnu+2.0d0)*t2*
     #              cx1(j+1)**2*cx2(j+1)*t3*t3
c
                v3=dexp(gamma_log(rnu+dj+3.0d0)-gamma_log(dj+1.0d0))
                v4=dexp(gamma_log(dj+1.0d0)+gamma_log(beta+1.0d0)-
     #              gamma_log(dj+1.0d0+beta))
                v5=dexp(gamma_log(dj+1.0d0)+gamma_log(beta1+1.0d0)-
     #              gamma_log(dj+1.0d0+beta1))

                v4=v4*v4*v5

                summv=summv+
     $             v3*v4*cv1(j+1)*cv1(j+1)*cv2(j+1)*
     #             (-1.0d0)**j*(2.0d0*dj+rnu+3.0d0)
c

 3            continue
              alpha(n+1)=t1*summ
              alfone(n+1)=v1*summv
 2          continue

            do 1 i=1,8
              you1=dwavex(i,istar)   
              c(1)=(1.0d0-you1)/(1.0d0-you1/3.0D0)
              c(2)=you1/(1.0d0-you1/3.0D0)      
              atot=alpha(1)*c(1)+alpha(2)*c(2)

c              if(istar.eq.1)refflux2(i)=-refflux1(i)*(atot)
c              if(istar.eq.2)refflux1(i)=-refflux2(i)*(atot)

              if(istar.eq.1)corr1(i)=-refflux1(i)*(atot)
              if(istar.eq.2)corr2(i)=-refflux2(i)*(atot)
c
              vtot=alfone(1)*c(1)+alfone(2)*c(2)
              if(delta.ne.0.0d0)then
                delvel=vstar/delta*(vtot/(1.0d0-atot))
              else
                delvel=0.0d0
              endif
              if(velK2.ne.0.0d0)then
                gimvel(i)=delvel/velK2
              else
                gimvel(i)=0.0d0
              endif
1           continue
          endif  ! end if ilaw=1
c

          if(ilaw.eq.4)then   !quadratic
c
            do 200 n=0,2
              rnu=dble(n+2)*0.5d0
c
              pee=rnu+2.0d0
              cue=rnu+1.0d0
              alf=pee-cue
              beta=cue-1.0d0
              xx=1.0d0-(2.0d0*(1.0d0-bee))
              call jacobi_poly(isw12,beta,alf,xx,cx1)
              cue=1.0d0
              alf=pee-cue
              beta=cue-1.0d0
              xx=1.0d0-2.0d0*cee*cee
              call jacobi_poly(isw12,beta,alf,xx,cx2)
              t1=bee*bee*((1.0d0-cee*cee)**(rnu+1.0d0))
     &             /(rnu*dexp(gamma_log(rnu+1.0d0)))
c
              pee=rnu+3.0d0
              cue=2.0d0
              alf=pee-cue
              beta1=cue-1.0d0
              xx=1.0d0-(2.0d0*(cee*cee))
              call jacobi_poly(isw12,beta1,alf,xx,cv2)

              pee=rnu+3.0d0
              cue=rnu+2.0d0
              alf=pee-cue
              beta=cue-1.0d0
              xx=1.0d0-(2.0d0*(1.0d0-bee))
              call jacobi_poly(isw12,beta,alf,xx,cv1)
c
              v1=dexp(gamma_log(rnu)-2.0d0*gamma_log(rnu+2.0d0))
              v1=v1*cee*bee*bee*(1.0d0-bee)*(1.0d0-cee*cee)**(rnu+1.0d0)

              summ=0.0d0
              summv=0.0d0
              summ=0.0d0
              do 300 j=0,isw12
                dj=dble(j)

c                  t3=(gamma_log(dj+1.0d0)+gamma_log(dj+rnu+2.0d0)-
c     %              gamma_log(2.0d0*dj+rnu+2.0d0))
c                  t3=exp(3.0d0*t3)

                t3=(gamma_log(dj+1.0d0)+gamma_log(rnu+1.0d0)-
     %              gamma_log(dj+rnu+1.0d0))
                t3=exp(t3)

                t2=dexp(gamma_log(rnu+dble(j+1))-gamma_log(dble(j+2)))
                summ=summ+(-1.0d0)**j*(2.0d0*dble(j)+rnu+2.0d0)*t2*
     #              cx1(j+1)**2*cx2(j+1)*t3*t3
c
                v3=dexp(gamma_log(rnu+dj+3.0d0)-gamma_log(dj+1.0d0))
                v4=dexp(gamma_log(dj+1.0d0)+gamma_log(beta+1.0d0)-
     #              gamma_log(dj+1.0d0+beta))
                v5=dexp(gamma_log(dj+1.0d0)+gamma_log(beta1+1.0d0)-
     #              gamma_log(dj+1.0d0+beta1))

                v4=v4*v4*v5

                summv=summv+
     $             v3*v4*cv1(j+1)*cv1(j+1)*cv2(j+1)*
     #             (-1.0d0)**j*(2.0d0*dj+rnu+3.0d0)
c
300           continue
              alpha(n+1)=t1*summ
              alfone(n+1)=v1*summv
200         continue
c
            do 100 i=1,8
              you1=dwavex(i,1)+2.0d0*dwavey(i,1)
              you2=-1.0d0*dwavey(i,1)   
              c(1)=(1.0d0-you1-you2)/(1.0d0-you1/3.0d0-you2/2.0d0)
              c(2)=you1/(1.0d0-you1/3.0d0-you2/2.0d0)      
              c(3)=you2/(1.0d0-you1/3.0d0-you2/2.0d0)      
              atot=alpha(1)*c(1)+alpha(2)*c(2)+alpha(3)*c(3)

c              refflux2(i)=-refflux1(i)*(atot)

              if(istar.eq.1)corr1(i)=-refflux1(i)*(atot)
              if(istar.eq.2)corr2(i)=-refflux2(i)*(atot)
c
              vtot=alfone(1)*c(1)+alfone(2)*c(2)+alfone(3)*c(3)
              if(delta.ne.0.0d0)then
                delvel=vstar/delta*(vtot/(1.0d0-atot))
              else
                delvel=0.0d0
              endif
              if(velK2.ne.0.0d0)then
                gimvel(i)=delvel/velK2
              else
                gimvel(i)=0.0d0
              endif
100           continue
          endif  ! end if ilaw=4
c
          return
          end
c
c   %$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$
c

c
          double precision function gamma_log ( x )
c
c  April 24, 2006
c
c  Routine to compute the natural log of the gamma function.  Received by way
c  of A. Gimenez.  Converted to FORTRAN 77 by Orosz.
c
c*******************************************************************************
c
c! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
c
c  Discussion:
c
c    Computation is based on an algorithm outlined in references 1 and 2.
c    The program uses rational functions that theoretically approximate
c    log ( GAMMA(X) ) to at least 18 significant decimal digits.  The
c    approximation for 12 < X is from reference 3, while approximations
c    for X < 12.0 are similar to those in reference 1, but are unpublished.
c    The accuracy achieved depends on the arithmetic system, the compiler,
c    intrinsic functions, and proper selection of the machine-dependent
c    constants.
c
c  Modified:
c
c    16 June 1999
c
c  Authors:
c
c    W. J. Cody and L. Stoltz
c    Argonne National Laboratory
c
c  Reference:
c
c    # 1)
c    W. J. Cody and K. E. Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
c    Mathematics of Computation,
c    Volume 21, 1967, pages 198-203.
c
c    # 2)
c    K. E. Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
c
c    # 3)
c    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, real ( kind = 8 ) X, the argument of the Gamma function. 
c    X must be positive.
c
c    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma 
c    function of X.  If X <= 0.0, or if overflow would occur, the
c    program returns the value HUGE().
c
c  Machine-dependent constants:
c
c       - radix for the floating-point representation.
c
c    MAXEXP - the smallest positive power of BETA that overflows.
c
c    XBIG   - largest argument for which LN(GAMMA(X)) is representable
c             in the machine, i.e., the solution to the equation
c             LN(GAMMA(XBIG)) = BETA**MAXEXP.
c
c    XINF   - largest machine representable floating-point number;
c             approximately BETA**MAXEXP.
c
c    FRTBIG - Rough estimate of the fourth root of XBIG
c
c
c    Approximate values for some important machines are:
c
c                              BETA      MAXEXP         XBIG
c
c    CRAY-1        (S.P.)        2        8191       9.62D+2461
c    Cyber 180/855
c      under NOS   (S.P.)        2        1070       1.72D+319
c    IEEE (IBM/XT,
c      SUN, etc.)  (S.P.)        2         128       4.08D+36
c    IEEE (IBM/XT,
c    SUN, etc.)  (D.P.)        2        1024       2.55D+305
c    IBM 3033      (D.P.)       16          63       4.29D+73
c    VAX D-Format  (D.P.)        2         127       2.05D+36
c    VAX G-Format  (D.P.)        2        1023       1.28D+305
c
c
c                            FRTBIG
c
c    CRAY-1        (S.P.)   3.13D+615
c    Cyber 180/855
c      under NOS   (S.P.)   6.44D+79
c    IEEE (IBM/XT,
c      SUN, etc.)  (S.P.)   1.42D+9
c    IEEE (IBM/XT,
c      SUN, etc.)  (D.P.)   2.25D+76
c    IBM 3033      (D.P.)   2.56D+18
c    VAX D-Format  (D.P.)   1.20D+9
c    VAX G-Format  (D.P.)   1.89D+76
c
         implicit none
         double precision c,corr,d1,d2,d4,eps,frtbig,p1,p2,p4
         double precision pnt68,q1,q2,q4,res,sqrtpi,x,xbig,xden,xm1,xm2,xm4
         double precision xnum,xsq
c
         integer i
c
         dimension c(7),p1(8),p2(8),p4(8),q1(8),q2(8),q4(8)
c         
         parameter(d1=-5.772156649015328605195174D-01)
         parameter(d2=4.227843350984671393993777D-01,pnt68=0.6796875d+00)
         parameter(d4=1.791759469228055000094023D+00,frtbig=1.42d+09)
         parameter(sqrtpi = 0.9189385332046727417803297D+00,xbig=4.08d+36)
c
         data c/-1.910444077728D-03, 8.4171387781295D-04, -5.952379913043012D-04, 
     #            7.93650793500350248D-04, -2.777777777777681622553D-03, 
     #            8.333333333333333331554247D-02, 5.7083835261D-03/
c          
         data p1/4.945235359296727046734888D+00,2.018112620856775083915565D+02, 
     %          2.290838373831346393026739D+03,1.131967205903380828685045D+04, 
     #          2.855724635671635335736389D+04,3.848496228443793359990269D+04,
     #          2.637748787624195437963534D+04,7.225813979700288197698961D+03/
c
         data p2/4.974607845568932035012064D+00,5.424138599891070494101986D+02, 
     #           1.550693864978364947665077D+04,1.847932904445632425417223D+05,
     #           1.088204769468828767498470D+06,3.338152967987029735917223D+06, 
     #           5.106661678927352456275255D+06,3.074109054850539556250927D+06/
c
         data p4/1.474502166059939948905062D+04,2.426813369486704502836312D+06, 
     #           1.214755574045093227939592D+08,2.663432449630976949898078D+09,
     #           2.940378956634553899906876D+10,1.702665737765398868392998D+11,
     #           4.926125793377430887588120D+11,5.606251856223951465078242D+11/
c
         data q1/6.748212550303777196073036D+01,1.113332393857199323513008D+03, 
     #           7.738757056935398733233834D+03,2.763987074403340708898585D+04,
     #           5.499310206226157329794414D+04,6.161122180066002127833352D+04,
     #           3.635127591501940507276287D+04,8.785536302431013170870835D+03/
c
         data q2/1.830328399370592604055942D+02,7.765049321445005871323047D+03,
     #           1.331903827966074194402448D+05,1.136705821321969608938755D+06,
     #           5.267964117437946917577538D+06,1.346701454311101692290052D+07,
     #           1.782736530353274213975932D+07,9.533095591844353613395747D+06/
c
         data q4/2.690530175870899333379843D+03,6.393885654300092398984238D+05, 
     #           4.135599930241388052042842D+07,1.120872109616147941376570D+09, 
     #           1.488613728678813811542398D+10,1.016803586272438228077304D+11,
     #           3.417476345507377132798597D+11,4.463158187419713286462081D+11 /
c
c  Return immediately if the argument is out of range.
c
         if ((x.le.0.0D+00).or.(xbig.lt.x)) then
           gamma_log = 2.55d+305
           return
         end if
         eps = 1.d-36  !epsilon ( eps )
         if (x.le.eps) then
           res = -dlog ( x )
         else if (x.le.1.5D+00) then
           if (x.lt.pnt68) then
             corr = - dlog(x)
             xm1 = x
           else
            corr = 0.0D+00
            xm1 = (x - 0.5D+00) - 0.5D+00
           end if
           if ((x.le.0.5D+00).or.(pnt68.le.x)) then
             xden = 1.0D+00
             xnum = 0.0D+00
             do i = 1, 8
               xnum = xnum * xm1 + p1(i)
               xden = xden * xm1 + q1(i)
             end do
             res = corr+(xm1*(d1+xm1*(xnum/xden)))
           else
             xm2=(x-0.5D+00) - 0.5D+00
             xden = 1.0D+00
             xnum = 0.0D+00
             do i = 1, 8
               xnum = xnum * xm2 + p2(i)
               xden = xden * xm2 + q2(i)
             end do
             res=corr+xm2*(d2+xm2*(xnum/xden))
           end if
         else if (x.le.4.0D+00) then
           xm2 = x - 2.0D+00
           xden = 1.0D+00
           xnum = 0.0D+00
           do i = 1, 8
             xnum = xnum * xm2 + p2(i)
             xden = xden * xm2 + q2(i)
           end do
           res=xm2*(d2+xm2*(xnum/xden))
         else if (x.le.12.0D+00) then
           xm4 = x - 4.0D+00
           xden = - 1.0D+00
           xnum = 0.0D+00
           do i = 1, 8
             xnum = xnum * xm4 + p4(i)
             xden = xden * xm4 + q4(i)
           end do
           res=d4+xm4*(xnum/xden)
         else
           res = 0.0D+00
           if(x.le.frtbig) then
             res = c(7)
             xsq = x * x
             do i = 1, 6
               res=res/xsq + c(i)
             end do
           end if
           res=res/x
           corr=dlog( x)
           res = res + sqrtpi - 0.5D+00 * corr
           res = res + x * ( corr - 1.0D+00 )
         end if
         gamma_log = res
         return
         end
c
c
c
         subroutine jacobi_poly(n,alpha,beta,x,cx)
c
c   April 24, 2006
c
c   Subroutine to evaluate the Jacobi polynomials at x.  Routine received
c   from A. Gimenez and converted to FORTRAN 77 by Orosz.
c
c*******************************************************************************
c
c! JACOBI_POLY evaluates the Jacobi polynomials at X.
c
c  Differential equation:
c
c    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
c
c  Recursion:
c
c    P(0,ALPHA,BETA,X) = 1,
c
c    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
c
c    P(N,ALPHA,BETA,X)  = 
c      ( 
c        (2*N+ALPHA+BETA-1) 
c        * ((ALPHA**2-BETA**2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X) 
c        * P(N-1,ALPHA,BETA,X)
c        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
c      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
c
c  Restrictions:
c
c    -1 < ALPHA
c    -1 < BETA
c
c  Norm:
c
c    Integral ( -1 <= X <= 1 ) ( 1 - X )**ALPHA * ( 1 + X )**BETA 
c      * P(N,ALPHA,BETA,X)**2 dX 
c    = 2**(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
c      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
c
c  Special values:
c
c    P(N,ALPHA,BETA)(1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
c
c  Modified:
c
c    01 October 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz and Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964.
c
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.  Note
c    that polynomials 0 through N will be computed.
c
c    Input, real ( kind = 8 ) ALPHA, one of the parameters defining the Jacobi
c    polynomials, ALPHA must be greater than -1.
c
c    Input, real ( kind = 8 ) BETA, the second parameter defining the Jacobi
c    polynomials, BETA must be greater than -1.
c
c    Input, real ( kind = 8 ) X, the point at which the polynomials are 
c    to be evaluated.
c
c    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 Jacobi
c    polynomials at the point X.
c
          implicit none
          integer n,i
          double precision alpha,beta,cx,c1,c2,c3,c4,r_i,x
          dimension cx(n+1)
c
          if(alpha.le.-1.0D+00) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
            write ( *, '(a,g14.6)' ) '  Illegal input value of ALPHA = ', alpha
            write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
            stop
          end if
          if(beta.le.-1.0D+00) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
            write ( *, '(a,g14.6)' ) '  Illegal input value of BETA = ', beta
            write ( *, '(a)' ) '  But BETA must be greater than -1.'
            stop
          end if
         if(n.lt.0) then
           return
         end if
         cx(1) = 1.0D+00
         if(n.eq.0) then
           return
         end if
         cx(1+1) = (1.0D+00+0.5D+00*(alpha+beta))*x+0.5D+00*(alpha-beta)
         do i = 2, n
           r_i = dble(i) 
           c1=2.0D+00*r_i*(r_i+alpha+beta)*(2.0D+00*r_i-2.0D+00+alpha+beta)
           c2=(2.0D+00*r_i-1.0D+00+alpha+beta)*(2.0D+00*r_i+alpha+beta) 
     #      *(2.0D+00*r_i-2.0D+00+alpha+beta)
           c3=(2.0D+00*r_i-1.0D+00+alpha+beta)*(alpha+beta)*(alpha-beta)
           c4=-2.0D+00*(r_i-1.0D+00+alpha)*(r_i-1.0D+00+beta)  
     #      *(2.0D+00*r_i+alpha+beta)
           cx(i+1)=((c3+c2*x)*cx(i-1+1)+c4*cx(i-2+1))/c1
         end do
         return
         end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine getalflim(istar,ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     %      phase,finc,Q,psi0,omega,
     $      gradx,grady,gradz,xarray,yarray,zarray,
     $      Nhoriz,xhoriz,yhoriz,
     &      separ,bdist,mmdx,ialfmin,ialfmax,reff1)
c
c  October 9, 1999
c
c  This routine will compute the 'projection' factor of each grid element on
c  the star (istar=1 to do star 1, istar=2 to do star2), check for eclipses
c  (the horizon of the other body is in xhoriz,yhoriz), and return the
c  sky coordinates of the visible points.  Set iecheck = -1 to skip the
c  check for 
c  eclipses.
c
c  UPDATE JULY 4, 2004
c
c  Add jdum and MonteCarlo to the argument list.  If MonteCarlo > 10,
c  then use Monte Carlo integration to determine fractionally
c  eclipsed pixels.  If MonteCarlo < 10, then proceed as before
c  and use interpolation in getBBflux and getATMflux.
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          parameter(twopie=2.0d0*pie)
          dimension gradx(ialphmax*ibetmax),grady(ialphmax*ibetmax),
     $        gradz(ialphmax*ibetmax),
     $        xhoriz(Nhoriz),yhoriz(Nhoriz),xarray(ialphmax*ibetmax),
     $        yarray(ialphmax*ibetmax),zarray(ialphmax*ibetmax),
     %        ibetlim(ialphmax),
     #        mmdx(ialphmax,ibetmax)
c
c
c    Keep track of the smallest and largest values of ialf on star 1 that
c    are eclipsed.
c
c
            PHASER = (PHASE/180.0d0)*pie     !orbital phase in radians
            FINCR = (FINC/180.0d0)*pie       !orbital inclination in radians
c
c            delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
c            delta=bdist*dsqrt(delta)
c
c            if(delta.gt.reff1+reff2)go to 999
c          
            NBET4 = NBET*4
            DBETA = (pie/2.0d0)/NBET        ! step size in longitude
c
            AZ = DCOS(FINCR)
            IF (AZ.LT.0.0d0) AZ = 0.0d0
            AX = -DSIN(FINCR)*DCOS(PHASER)    ! l in Wilson & Sofia
            AY = DSIN(FINCR)*DSIN(PHASER)     ! m in Wilson & Sofia
c
c   Check to see of the star in question is in front.  If so, then simply
c   find the projection factors.
c
c            diff1=dabs(phaser-pconj)
c            diff2=dabs(phaser-pconj2)
c            if(diff1.lt.diff2)go to 999

            iimax=-12345
            iimin=12345
            DO 501 IALF = 1, NALF
              DO 502 IBET = 1,ibetlim(ialf)      !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
                iidx=mmdx(ialf,ibet)
                PROJ = AX * GRADX(iidx) + AY*GRADY(iidx) + 
     1	          AZ*GRADZ(iidx)
                IF (PROJ.LT.0.) GO TO 502    ! is the surface element visible?
                xx=xarray(iidx)
                yy=yarray(iidx)
                zz=zarray(iidx)
                xp=xtran(xx,yy,zz,phase,fincr,Q,istar,bdist) ! projected coords
                yp=ytran(xx,yy,zz,phase,fincr,Q,istar,bdist)
c
c   Check to see of the point in question is eclipsed by the other star, or
c   in the case of star 1, eclipsed by the disk, or
c   in the case of a point on the bottom half of star 2, eclipsed by the disk. 
c   icut=2 for points outside and below the horizon.
c 
                iyes=-100
                call insidecircle(Nhoriz,xhoriz,yhoriz,xp,yp,iyes,icut)
                if(iyes.eq.100)then
                  if(ialf.lt.ialfmin)ialfmin=ialf
                  if(ialf.lt.iimin)iimin=ialf
                  if(ialf.gt.ialfmax)ialfmax=ialf
                  if(ialf.gt.iimax)iimax=ialf
                  go to 502
                endif
c
502           CONTINUE
501         continue         ! continue the alpha loop
c
999         continue
c
          dsi=dsin(fincr)*dsin(fincr)
          ttop=1.0d0-reff1*reff1
c
c
 987      format(f8.6,2x,3(i4,1x))
c
 1001     format(2(f11.8,1x),'ialf=',i3,1x,'ibet=',i3)
          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine getBBlumcor(ialphmax,ibetmax,Nalf,Nbet,ibetlim,
     $      wave,visib,projarray,temp,surf,flimbx,flimby,ilaw,rinty,
     &      flum,flux,rldint,
     &      separ,mmdx,ialfmin,ialfmax)
c
c   May 3, 2006
c
c   This routine will return the integrate flux on star 1 outside the
c   ialpha range of ialfmin,ialfmax
c
c
c
          implicit double precision(a-h,o-z)
c
          parameter(pie=3.14159265358979323d0)
          dimension visib(ialphmax*ibetmax),
     $        surf(ialphmax*ibetmax),ibetlim(ialphmax),
     $        temp(ialphmax*ibetmax),flum(ialphmax*ibetmax),
     $        rinty(ialphmax*ibetmax),projarray(ialphmax*ibetmax),
     %        mmdx(ialphmax,ibetmax)
c
          dint=pie*(1.0d0-flimbx/3.0d0)
          if(ilaw.eq.2)then
            dint=pie*(1.0d0-flimbx/3.0d0+2.0d0*flimby/9.0d0)
          endif
          if(ilaw.eq.3)then
            dint=pie*(1.0d0-flimbx/3.0d0-flimby/5.0d0)
          endif
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
          if(ilaw.eq.4)then
            dint=pie*(1.0d0-flimbx/3.0d0-flimby/6.0d0)
          endif

          wavemu=wave/10000.0d0
          C2 = 1.4384d8          ! 1.4384 * 10.**8      ! hc/(k*1e-8)
          C1 = 1.191044d35       ! 2hc^2/((1e-8)**5)
c
c   Initialize the flum matrix.
c
c
          c1=3.74185
          c2=14.3883

          flux=0.0d0
          do 2 ialf=1,nalf
            do 1 ibet=1,ibetlim(ialf)        !4*Nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c
              iidx=mmdx(ialf,ibet)
              flum(iidx)=0.0d0
              rinty(iidx)=0.0d0
 1          continue
 2        continue
          sumcor1=0.0d0
          sumcor2=0.0d0
          dtheta=0.5d0*pie/dble(Nalf)
c
          DO 10 ialf=1,nalf
            if((ialf.ge.ialfmin).and.(ialf.le.ialfmax))go to 10
            theta=-dtheta+2.0d0*dtheta*dble(ialf)
            sitheta=dsin(theta)
            dphi=pie/dble(ibetlim(ialf))
            DO 9 ibet = 1,ibetlim(ialf)               !4*nbet
c
c   UPDATE June 11, 2003
c
c   change the 2D arrays into 1D
c
c              iidx=kount(ialphmax,ialf,ibetlim)+ibet
c 
              iidx=mmdx(ialf,ibet)
              corr1=0.0d0
              corr2=0.0d0
              if((projarray(iidx).le.0.0d0))go to 9
c              C3 = C2/(WAVE*TEMP(iidx))
c
              tkkelv=temp(iidx)/1000.0d0
              C3 = C2/(wavemu*tkkelv)
              flum(iidx)=C1/(dexp(c3)-1.0d0)/wavemu**5
              dark=(1.0d0-flimbx+flimbx*projarray(iidx))
              if(ilaw.eq.2)dark=dark-flimby*projarray(iidx)*
     %               dlog(projarray(iidx))
              if(ilaw.eq.3)dark=dark-flimby*(1.0-dsqrt(projarray(iidx)))
c
c   UPDATE JULY 21, 2004
c
c   Add a quadratic limb darkening law, ilaw=4
c
c
              if(ilaw.eq.4)dark=dark-flimby*(1.0-(projarray(iidx)))**2
c
              flum(iidx)=flum(iidx)*dark
              rinty(iidx)=flum(iidx) ! save intensities for plotting
              flum(iidx)=surf(iidx)*flum(iidx)*visib(iidx)
              flux=flux+flum(iidx)
 9          continue
 10       continue
c
c   Scale the light curve by the integral of the limb darkening law
c   for compatibility with Wilson-Devinney.
c        
c   Scale the light curve by the integral of the limb darkening law
c   for compatibility with Wilson-Devinney.
c        
c          flux=pie*flux/dint
c
c
c   UPDATE April 3, 2002
c
c   Scale the fluxes.
c
c          solarrad=6.9598d10
c          flux=flux*(separ*solarrad)**2
cc
c          rldint=dint
          return
          end
c
c
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine newcomputeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &       atmT,atmg,atmmu,Nmu,
     &       atmint,Tmax,Tmin,gmax,gmin,outinty,
     #       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess)
c
          implicit double precision (a-h,o-z)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       atmint(maxlines,maxmu,8),Nmu(maxlines),outinty(8),tempgh(200),
     &       tempgl(200),ghinty(200,8),yscratch(200),tscratch(2),glinty(200,8),
     &       tinty(2,8),y2scratch(2)
c

          dimension scalex1(200,20),scaley1(200,20)
          dimension scalex2(200,20),scaley2(200,20)
          dimension scalex3(200,20),scaley3(200,20)
          dimension scalex4(200,20),scaley4(200,20)
          dimension scalex5(200,20),scaley5(200,20)
          dimension scalex6(200,20),scaley6(200,20)
          dimension scalex7(200,20),scaley7(200,20)
          dimension scalex8(200,20),scaley8(200,20)
c          dimension stdscalex(100),stdscaley(100)
c          dimension stdscalex(50),stdscaley(50)
          dimension stdscalex(36),stdscaley(36)
          dimension Ngmu(20)
          dimension gscratch1(2),gscratch2(2),gscratch3(2),gscratch4(2)
          dimension gscratch5(2),gscratch6(2),gscratch7(2),gscratch8(2)
          dimension gx1(20),gy1(20)
          dimension gx2(20),gy2(20)
          dimension gx3(20),gy3(20)
          dimension gx4(20),gy4(20)
          dimension gx5(20),gy5(20)
          dimension gx6(20),gy6(20)
          dimension gx7(20),gy7(20)
          dimension gx8(20),gy8(20)
          dimension finalx1(200),finaly1(200)
          dimension finalx2(200),finaly2(200)
          dimension finalx3(200),finaly3(200)
          dimension finalx4(200),finaly4(200)
          dimension finalx5(200),finaly5(200)
          dimension finalx6(200),finaly6(200)
          dimension finalx7(200),finaly7(200)
          dimension finalx8(200),finaly8(200)
          dimension scrx(200),scry(200)
          dimension peakx(20)
          dimension peaky1(20),peaky2(20),peaky3(20),peaky4(20)
          dimension peaky5(20),peaky6(20),peaky7(20),peaky8(20)
          data stdscalex/0.003,0.01d0,0.03d0,0.06d0,0.09d0,0.12d0,
     %                   0.15d0,0.18d0,0.21d0,0.24d0,0.27d0,
     #                   0.30d0,0.33d0,0.36d0,0.39d0,
     %                   0.42d0,0.45d0,0.48d0,0.51d0,
     #                   0.54d0,0.57d0,0.60d0,
     %                   0.63d0,0.66d0,0.69d0,0.72d0,0.75d0,
     #                   0.78d0,0.81d0,0.84d0,0.87d0,
     %                   0.90d0,0.93d0,0.96d0,0.99,1.00d0/

c          data stdscalex/0.02d0,0.04d0,
c     #                   0.06d0,0.08d0,0.10d0,
c     %                   0.12d0,0.14d0,0.16d0,
c     #                   0.18d0,0.20d0,
c     %                   0.22d0,0.24d0,0.26d0,
c     #                   0.28d0,0.30d0,
c     %                   0.32d0,0.34d0,0.36d0,
c     #                   0.38d0,0.40d0,
c     %                   0.42d0,0.44d0,0.46d0,
c     #                   0.48d0,0.50d0,
c     %                   0.52d0,0.54d0,0.56d0,
c     #                   0.58d0,0.60d0,
c     %                   0.62d0,0.64d0,0.66d0,
c     #                   0.68d0,0.70d0,
c     %                   0.72d0,0.74d0,0.76d0,
c     #                   0.78d0,0.80d0,
c     %                   0.82d0,0.84d0,0.86d0,
c     #                   0.88d0,0.90d0,
c     %                   0.92d0,0.94d0,0.96d0,
c     #                   0.98d0,1.00d0/

C          data stdscalex/0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,
C     #                   0.06d0,0.07d0,0.08d0,0.09d0,0.10d0,
C     %                   0.11d0,0.12d0,0.13d0,0.14d0,0.15d0,
C     #                   0.16d0,0.17d0,0.18d0,0.19d0,0.20d0,
C     %                   0.21d0,0.22d0,0.23d0,0.24d0,0.25d0,
C     #                   0.26d0,0.27d0,0.28d0,0.29d0,0.30d0,
C     %                   0.31d0,0.32d0,0.33d0,0.34d0,0.35d0,
C     #                   0.36d0,0.37d0,0.38d0,0.39d0,0.40d0,
C     %                   0.41d0,0.42d0,0.43d0,0.44d0,0.45d0,
C     #                   0.46d0,0.47d0,0.48d0,0.49d0,0.50d0,
C     %                   0.51d0,0.52d0,0.53d0,0.54d0,0.55d0,
C     #                   0.56d0,0.57d0,0.58d0,0.59d0,0.60d0,
C     %                   0.61d0,0.62d0,0.63d0,0.64d0,0.65d0,
C     #                   0.66d0,0.67d0,0.68d0,0.69d0,0.70d0,
C     %                   0.71d0,0.72d0,0.73d0,0.74d0,0.75d0,
C     #                   0.76d0,0.77d0,0.78d0,0.79d0,0.80d0,
C     %                   0.81d0,0.82d0,0.83d0,0.84d0,0.85d0,
C     #                   0.86d0,0.87d0,0.88d0,0.89d0,0.90d0,
C     %                   0.91d0,0.92d0,0.93d0,0.94d0,0.95d0,
C     #                   0.96d0,0.97d0,0.98d0,0.99d0,1.00d0/
c
c
c         if(rmuin.gt.0.69d0)then
c            call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
c     &        atmT,atmg,atmmu,Nmu,
c     &        atmint,Tmax,Tmin,gmax,gmin,outinty)
c            return
c          endif
c
c          write(*,*)icnU,icnB,icnV,icnR

          call locate(atmT,Nlines,Tin,indexT)
c
c          call hunt(atmT,Nlines,Tin,itguess)
c
c          indexT=itguess
c
c   indexT is index where the temperature is less than Tin.
c
c   For all temperature values equal to atmT(indexT), produce arrays
c   of (inty,mu), and interpolate mu(inty=stdval) vs. g curves to
c   produce the limb darkening curve for log(g)=gin.  Then given this,
c   produce the intensity at mu=rmuin.
c
          tscratch(1)=atmT(indexT+1)
c
c   Now search the Tvalues equal to atmT(indexT+1) and find intensities
c   for all of the g values.
c
          igcount=0
          do 600 i=0,Nlines-(indexT+1)
            if(atmT(indexT+1+i).eq.atmT(indexT+1))then
              igcount=igcount+1
              tempgl(igcount)=atmg(indexT+1+i)
            endif
 600      continue
          call locate(tempgl,igcount,gin,index1)
          if(index1.le.2)then
            gl1=tempgl(1)
            gl2=tempgl(2)
            gl3=tempgl(3)
          endif
          if(index1.ge.igcount)then
            gl1=tempgl(igcount-2)
            gl2=tempgl(igcount-1)
            gl3=tempgl(igcount)
          endif
          if((index1.gt.2).and.(index1.lt.igcount))then
            gl1=tempgl(index1-1)
            gl2=tempgl(index1)
            gl3=tempgl(index1+1)
          endif
c
          N=36
          Nhigh=0
          igcount=0
          do 10 i=0,Nlines-(indexT+1)
            if(atmT(indexT+1+i).eq.atmT(indexT+1))then
              fredg=atmg(indexT+1+i)
              if((fredg.eq.gl1).or.(fredg.eq.gl2).or.(fredg.eq.gl3))then
                igcount=igcount+1
                tempgh(igcount)=atmg(indexT+1+i)
c                Ngmu(igcount)=Nmu(indexT+1+i)
                Ngmu(igcount)=0
                do 1 j=1,Nmu(indexT+1+i)
                  tempmu=atmmu(indexT+1+i,j)
                  dl=rmuin  -1.0d0
                  dh=rmuin  +1.0d0
                  if((j.eq.1).or.(j.eq.Nmu(indexT+1+i)).or.
     $              ((tempmu.ge.dl).and.(tempmu.le.dh)))then
c
                    Ngmu(igcount)=Ngmu(igcount)+1
                    jjj=Ngmu(igcount)
c
                    if(icnU.ne.430)then
                      scalex1(jjj,igcount)=atmint(indexT+1+i,j,1)
                      scaley1(jjj,igcount)=tempmu
                    endif
                    if(icnB.ne.430)then
                      scalex2(jjj,igcount)=atmint(indexT+1+i,j,2)
                      scaley2(jjj,igcount)=tempmu
                    endif
                    if(icnV.ne.430)then
                      scalex3(jjj,igcount)=atmint(indexT+1+i,j,3)
                      scaley3(jjj,igcount)=tempmu
                    endif
                    if(icnR.ne.430)then
                      scalex4(jjj,igcount)=atmint(indexT+1+i,j,4)
                      scaley4(jjj,igcount)=tempmu
                    endif
                    if(icnI.ne.430)then
                      scalex5(jjj,igcount)=atmint(indexT+1+i,j,5)
                      scaley5(jjj,igcount)=tempmu
                    endif
                    if(icnJ.ne.430)then
                      scalex6(jjj,igcount)=atmint(indexT+1+i,j,6)
                      scaley6(jjj,igcount)=tempmu
                    endif
                    if(icnH.ne.430)then
                      scalex7(jjj,igcount)=atmint(indexT+1+i,j,7)
                      scaley7(jjj,igcount)=tempmu
                    endif
                    if(icnK.ne.430)then
                      scalex8(jjj,igcount)=atmint(indexT+1+i,j,8)
                      scaley8(jjj,igcount)=tempmu
                    endif
                  endif
1               continue
c                write(*,*)jjj,igcount
              endif
            else
              go to 15
            endif
10        continue
c
c   Establish the intensity at mu=1 for log(g)=gin
c
15         do 500 ig=1,igcount
             if(icnU.ne.430)peaky1(ig)=scalex1(Ngmu(ig),ig)
             if(icnB.ne.430)peaky2(ig)=scalex2(Ngmu(ig),ig)
             if(icnV.ne.430)peaky3(ig)=scalex3(Ngmu(ig),ig)
             if(icnR.ne.430)peaky4(ig)=scalex4(Ngmu(ig),ig)
             if(icnI.ne.430)peaky5(ig)=scalex5(Ngmu(ig),ig)
             if(icnJ.ne.430)peaky6(ig)=scalex6(Ngmu(ig),ig)
             if(icnH.ne.430)peaky7(ig)=scalex7(Ngmu(ig),ig)
             if(icnK.ne.430)peaky8(ig)=scalex8(Ngmu(ig),ig)
500        continue
c
           call locate(tempgh,igcount,gin,index1)
           m=2
           k=min(max(index1-(m-1)/2,1),igcount+1-m)
           if(icnU.ne.430)call polint(tempgh(k),peaky1(k),m,gin,qqq1,dy)
           if(icnB.ne.430)call polint(tempgh(k),peaky2(k),m,gin,qqq2,dy)
           if(icnV.ne.430)call polint(tempgh(k),peaky3(k),m,gin,qqq3,dy)
           if(icnR.ne.430)call polint(tempgh(k),peaky4(k),m,gin,qqq4,dy)
           if(icnI.ne.430)call polint(tempgh(k),peaky5(k),m,gin,qqq5,dy)
           if(icnJ.ne.430)call polint(tempgh(k),peaky6(k),m,gin,qqq6,dy)
           if(icnH.ne.430)call polint(tempgh(k),peaky7(k),m,gin,qqq7,dy)
           if(icnK.ne.430)call polint(tempgh(k),peaky8(k),m,gin,qqq8,dy)
c
c   Now loop over the scalex1 array and find the mu where the intensity
c   is equal to the peak, 0.98*peak, 0.96*peak, etc.
c
           do 20 i=1,N
             do 11 ig=1,igcount
               if(icnU.ne.430)scale1=scalex1(Ngmu(ig),ig)*stdscalex(i)
               if(icnB.ne.430)scale2=scalex2(Ngmu(ig),ig)*stdscalex(i)
               if(icnV.ne.430)scale3=scalex3(Ngmu(ig),ig)*stdscalex(i)
               if(icnR.ne.430)scale4=scalex4(Ngmu(ig),ig)*stdscalex(i)
               if(icnI.ne.430)scale5=scalex5(Ngmu(ig),ig)*stdscalex(i)
               if(icnJ.ne.430)scale6=scalex6(Ngmu(ig),ig)*stdscalex(i)
               if(icnH.ne.430)scale7=scalex7(Ngmu(ig),ig)*stdscalex(i)
               if(icnK.ne.430)scale8=scalex8(Ngmu(ig),ig)*stdscalex(i)
c
               ifred=0
               rlowmu=rmuin-0.1d0
               rhighmu=rmuin+0.1d0

               if(icnU.ne.430)then
                 do 101 j=1,Ngmu(ig)
                   scrx(j)=scalex1(j,ig)
                   scry(j)=scaley1(j,ig)
101              continue
                 call locate(scrx,Ngmu(ig),scale1,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale1,qqq,dy)
                 gx1(ig)=tempgh(ig)
                 gy1(ig)=qqq
               endif
c
               if(icnB.ne.430)then
                 do 102 j=1,Ngmu(ig)
                   scrx(j)=scalex2(j,ig)
                   scry(j)=scaley2(j,ig)
102              continue
                 call locate(scrx,Ngmu(ig),scale2,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale2,qqq,dy)
                 gx2(ig)=tempgh(ig)
                 gy2(ig)=qqq
               endif
c
               if(icnV.ne.430)then
                 do 103 j=1,Ngmu(ig)
                   scrx(j)=scalex3(j,ig)
                   scry(j)=scaley3(j,ig)
103              continue
                 call locate(scrx,Ngmu(ig),scale3,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale3,qqq,dy)
                 gx3(ig)=tempgh(ig)
                 gy3(ig)=qqq
               endif
c
               if(icnR.ne.430)then
                 do 104 j=1,Ngmu(ig)
                   scrx(j)=scalex4(j,ig)
                   scry(j)=scaley4(j,ig)
104              continue
                 call locate(scrx,Ngmu(ig),scale4,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale4,qqq,dy)
                 gx4(ig)=tempgh(ig)
                 gy4(ig)=qqq
               endif
c
               if(icnI.ne.430)then
                 do 105 j=1,Ngmu(ig)
                   scrx(j)=scalex5(j,ig)
                   scry(j)=scaley5(j,ig)
105              continue
                 call locate(scrx,Ngmu(ig),scale5,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale5,qqq,dy)
                 gx5(ig)=tempgh(ig)
                 gy5(ig)=qqq
               endif
c
               if(icnJ.ne.430)then
                 do 106 j=1,Ngmu(ig)
                   scrx(j)=scalex6(j,ig)
                   scry(j)=scaley6(j,ig)
106              continue
                 call locate(scrx,Ngmu(ig),scale6,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale6,qqq,dy)
                 gx6(ig)=tempgh(ig)
                 gy6(ig)=qqq
               endif
c
               if(icnH.ne.430)then
                 do 107 j=1,Ngmu(ig)
                   scrx(j)=scalex7(j,ig)
                   scry(j)=scaley7(j,ig)
107              continue
                 call locate(scrx,Ngmu(ig),scale7,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale7,qqq,dy)
                 gx7(ig)=tempgh(ig)
                 gy7(ig)=qqq
               endif
c
               if(icnK.ne.430)then
                 do 108 j=1,Ngmu(ig)
                   scrx(j)=scalex8(j,ig)
                   scry(j)=scaley8(j,ig)
108              continue
                 call locate(scrx,Ngmu(ig),scale8,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale8,qqq,dy)
                 gx8(ig)=tempgh(ig)
                 gy8(ig)=qqq
               endif
c
11           continue
c
             if(icnU.ne.430)then
               call polint(gx1,gy1,igcount,gin,qqq,dy)
               finalx1(i)=qqq
               finaly1(i)=stdscalex(i)*qqq1
             endif
c
             if(icnB.ne.430)then
               call polint(gx2,gy2,igcount,gin,qqq,dy)
               finalx2(i)=qqq
               finaly2(i)=stdscalex(i)*qqq2
             endif
c
             if(icnV.ne.430)then
               call polint(gx3,gy3,igcount,gin,qqq,dy)
               finalx3(i)=qqq
               finaly3(i)=stdscalex(i)*qqq3
             endif
c
             if(icnR.ne.430)then
               call polint(gx4,gy4,igcount,gin,qqq,dy)
               finalx4(i)=qqq
               finaly4(i)=stdscalex(i)*qqq4
             endif
c
             if(icnI.ne.430)then
               call polint(gx5,gy5,igcount,gin,qqq,dy)
               finalx5(i)=qqq
               finaly5(i)=stdscalex(i)*qqq5
             endif
c
             if(icnJ.ne.430)then
               call polint(gx6,gy6,igcount,gin,qqq,dy)
               finalx6(i)=qqq
               finaly6(i)=stdscalex(i)*qqq6
             endif
c
             if(icnH.ne.430)then
               call polint(gx7,gy7,igcount,gin,qqq,dy)
               finalx7(i)=qqq
               finaly7(i)=stdscalex(i)*qqq7
             endif
c
             if(icnK.ne.430)then
               call polint(gx8,gy8,igcount,gin,qqq,dy)
               finalx8(i)=qqq
               finaly8(i)=stdscalex(i)*qqq8
             endif
c
20         continue
c
           if(icnU.ne.430)then
             call locate(finalx1,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx1(k),finaly1(k),m,rmuin,qqq,dy)
               gscratch1(1)=qqq
             else
               gscratch1(1)=finaly1(1)
             endif
           endif
c
           if(icnB.ne.430)then
             call locate(finalx2,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx2(k),finaly2(k),m,rmuin,qqq,dy)
               gscratch2(1)=qqq
             else
               gscratch2(1)=finaly2(1)
             endif
           endif
c
           if(icnV.ne.430)then
             call locate(finalx3,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx3(k),finaly3(k),m,rmuin,qqq,dy)
               gscratch3(1)=qqq
             else
               gscratch3(1)=finaly3(1)
             endif
           endif
c
           if(icnR.ne.430)then
             call locate(finalx4,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx4(k),finaly4(k),m,rmuin,qqq,dy)
               gscratch4(1)=qqq
             else
               gscratch4(1)=finaly4(1)
             endif
           endif
c
           if(icnI.ne.430)then
             call locate(finalx5,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx5(k),finaly5(k),m,rmuin,qqq,dy)
               gscratch5(1)=qqq
             else
               gscratch5(1)=finaly5(1)
             endif
           endif
c
           if(icnJ.ne.430)then
             call locate(finalx6,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx6(k),finaly6(k),m,rmuin,qqq,dy)
               gscratch6(1)=qqq
             else
               gscratch6(1)=finaly6(1)
             endif
           endif
c
           if(icnH.ne.430)then
             call locate(finalx7,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx7(k),finaly7(k),m,rmuin,qqq,dy)
               gscratch7(1)=qqq
             else
               gscratch7(1)=finaly7(1)
             endif
           endif
c
           if(icnK.ne.430)then
             call locate(finalx8,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx8(k),finaly8(k),m,rmuin,qqq,dy)
               gscratch8(1)=qqq
             else
               gscratch8(1)=finaly8(1)
             endif
           endif
c
c   Now search the Tvalues equal to atmT(indexT+1) and find intensities
c   for all of the g values.
c
          tscratch(2)=atmT(indexT)
c
c   Now search the Tvalues equal to atmT(indexT+1) and find intensities
c   for all of the g values.
c
          igcount=0
          do 601 i=0,Nlines-(indexT+1)
            if(atmT(indexT+1+i).eq.atmT(indexT+1))then
              igcount=igcount+1
              tempgl(igcount)=atmg(indexT+1+i)
            endif
 601     continue
          call locate(tempgl,igcount,gin,index1)
          if(index1.le.2)then
            gl1=tempgl(1)
            gl2=tempgl(2)
            gl3=tempgl(3)
          endif
          if(index1.ge.igcount)then
            gl1=tempgl(igcount-2)
            gl2=tempgl(igcount-1)
            gl3=tempgl(igcount)
          endif
          if((index1.gt.2).and.(index1.lt.igcount))then
            gl1=tempgl(index1-1)
            gl2=tempgl(index1)
            gl3=tempgl(index1+1)
          endif
c
          Nhigh=0
          igcount=0
          do 30 i=0,indexT-1
            if(atmT(indexT-i).eq.atmT(indexT))then
              fredg=atmg(indexT-i)
              if((fredg.eq.gl1).or.(fredg.eq.gl2).or.(fredg.eq.gl3))then
                igcount=igcount+1
                tempgh(igcount)=atmg(indexT-i)
c                Ngmu(igcount)=Nmu(indexT-i)
                Ngmu(igcount)=0
                do 31 j=1,Nmu(indexT-i)
                  dl=rmuin -1.0d0
                  dh=rmuin +1.0d0
                  if((j.eq.1).or.(j.eq.Nmu(indexT+1+i)).or.
     $              ((tempmu.ge.dl).and.(tempmu.le.dh)))then
                    Ngmu(igcount)=Ngmu(igcount)+1
                    jjj=Ngmu(igcount)
                    if(icnU.ne.430)then
                      scalex1(jjj,igcount)=atmint(indexT-i,j,1)
                      scaley1(jjj,igcount)=atmmu(indexT-i,j)
                    endif
                    if(icnB.ne.430)then
                      scalex2(jjj,igcount)=atmint(indexT-i,j,2)
                      scaley2(jjj,igcount)=atmmu(indexT-i,j)
                    endif
                    if(icnV.ne.430)then
                      scalex3(jjj,igcount)=atmint(indexT-i,j,3)
                      scaley3(jjj,igcount)=atmmu(indexT-i,j)
                    endif
                    if(icnR.ne.430)then
                      scalex4(jjj,igcount)=atmint(indexT-i,j,4)
                      scaley4(jjj,igcount)=atmmu(indexT-i,j)
                    endif
                    if(icnI.ne.430)then
                      scalex5(jjj,igcount)=atmint(indexT-i,j,5)
                      scaley5(jjj,igcount)=atmmu(indexT-i,j)
                    endif
                    if(icnJ.ne.430)then
                      scalex6(jjj,igcount)=atmint(indexT-i,j,6)
                      scaley6(jjj,igcount)=atmmu(indexT-i,j)
                    endif
                    if(icnH.ne.430)then
                      scalex7(jjj,igcount)=atmint(indexT-i,j,7)
                      scaley7(jjj,igcount)=atmmu(indexT-i,j)
                    endif
                    if(icnK.ne.430)then
                      scalex8(j,igcount)=atmint(indexT-i,j,8)
                      scaley8(j,igcount)=atmmu(indexT-i,j)
                    endif
                  endif
31              continue
              endif
            else
              go to 25
            endif
30        continue
25        do 501 ig=1,igcount
            if(icnU.ne.430)peaky1(ig)=scalex1(Ngmu(ig),ig)
            if(icnB.ne.430)peaky2(ig)=scalex2(Ngmu(ig),ig)
            if(icnV.ne.430)peaky3(ig)=scalex3(Ngmu(ig),ig)
            if(icnR.ne.430)peaky4(ig)=scalex4(Ngmu(ig),ig)
            if(icnI.ne.430)peaky5(ig)=scalex5(Ngmu(ig),ig)
            if(icnJ.ne.430)peaky6(ig)=scalex6(Ngmu(ig),ig)
            if(icnH.ne.430)peaky7(ig)=scalex7(Ngmu(ig),ig)
            if(icnK.ne.430)peaky8(ig)=scalex8(Ngmu(ig),ig)
501       continue
c
          call locate(tempgh,igcount,gin,index1)
          m=2
          k=min(max(index1-(m-1)/2,1),igcount+1-m)
          if(icnU.ne.430)call polint(tempgh(k),peaky1(k),m,gin,qqq1,dy)
          if(icnB.ne.430)call polint(tempgh(k),peaky2(k),m,gin,qqq2,dy)
          if(icnV.ne.430)call polint(tempgh(k),peaky3(k),m,gin,qqq3,dy)
          if(icnR.ne.430)call polint(tempgh(k),peaky4(k),m,gin,qqq4,dy)
          if(icnI.ne.430)call polint(tempgh(k),peaky5(k),m,gin,qqq5,dy)
          if(icnJ.ne.430)call polint(tempgh(k),peaky6(k),m,gin,qqq6,dy)
          if(icnH.ne.430)call polint(tempgh(k),peaky7(k),m,gin,qqq7,dy)
          if(icnK.ne.430)call polint(tempgh(k),peaky8(k),m,gin,qqq8,dy)
c
c   Now loop over the scalex1 array and find the mu where the intensity
c   is equal to the peak, 0.98*peak, 0.96*peak, etc.
c
           do 40 i=1,N
             do 21 ig=1,igcount
               if(icnU.ne.430)scale1=scalex1(Ngmu(ig),ig)*stdscalex(i)
               if(icnB.ne.430)scale2=scalex2(Ngmu(ig),ig)*stdscalex(i)
               if(icnV.ne.430)scale3=scalex3(Ngmu(ig),ig)*stdscalex(i)
               if(icnR.ne.430)scale4=scalex4(Ngmu(ig),ig)*stdscalex(i)
               if(icnI.ne.430)scale5=scalex5(Ngmu(ig),ig)*stdscalex(i)
               if(icnJ.ne.430)scale6=scalex6(Ngmu(ig),ig)*stdscalex(i)
               if(icnH.ne.430)scale7=scalex7(Ngmu(ig),ig)*stdscalex(i)
               if(icnK.ne.430)scale8=scalex8(Ngmu(ig),ig)*stdscalex(i)
c
               if(icnU.ne.430)then
                 do 201 j=1,Ngmu(ig)
                   scrx(j)=scalex1(j,ig)
                   scry(j)=scaley1(j,ig)
201              continue
                 call locate(scrx,Ngmu(ig),scale1,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale1,qqq,dy)
                 gx1(ig)=tempgh(ig)
                 gy1(ig)=qqq
               endif
c
               if(icnB.ne.430)then
                 do 202 j=1,Ngmu(ig)
                   scrx(j)=scalex2(j,ig)
                   scry(j)=scaley2(j,ig)
202              continue
                 call locate(scrx,Ngmu(ig),scale2,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale2,qqq,dy)
                 gx2(ig)=tempgh(ig)
                 gy2(ig)=qqq
               endif
c
               if(icnV.ne.430)then
                 do 203 j=1,Ngmu(ig)
                   scrx(j)=scalex3(j,ig)
                   scry(j)=scaley3(j,ig)
203              continue
                 call locate(scrx,Ngmu(ig),scale3,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale3,qqq,dy)
                 gx3(ig)=tempgh(ig)
                 gy3(ig)=qqq
               endif
c
               if(icnR.ne.430)then
                 do 204 j=1,Ngmu(ig)
                   scrx(j)=scalex4(j,ig)
                   scry(j)=scaley4(j,ig)
204              continue
                 call locate(scrx,Ngmu(ig),scale4,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale4,qqq,dy)
                 gx4(ig)=tempgh(ig)
                 gy4(ig)=qqq
               endif
c
               if(icnI.ne.430)then
                 do 205 j=1,Ngmu(ig)
                   scrx(j)=scalex5(j,ig)
                   scry(j)=scaley5(j,ig)
205              continue
                 call locate(scrx,Ngmu(ig),scale5,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale5,qqq,dy)
                 gx5(ig)=tempgh(ig)
                 gy5(ig)=qqq
               endif
c
               if(icnJ.ne.430)then
                 do 206 j=1,Ngmu(ig)
                   scrx(j)=scalex6(j,ig)
                   scry(j)=scaley6(j,ig)
206              continue
                 call locate(scrx,Ngmu(ig),scale6,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale6,qqq,dy)
                 gx6(ig)=tempgh(ig)
                 gy6(ig)=qqq
               endif
c
               if(icnH.ne.430)then
                 do 207 j=1,Ngmu(ig)
                   scrx(j)=scalex7(j,ig)
                   scry(j)=scaley7(j,ig)
207              continue
                 call locate(scrx,Ngmu(ig),scale7,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale7,qqq,dy)
                 gx7(ig)=tempgh(ig)
                 gy7(ig)=qqq
               endif
c
               if(icnK.ne.430)then
                 do 208 j=1,Ngmu(ig)
                   scrx(j)=scalex8(j,ig)
                   scry(j)=scaley8(j,ig)
208              continue
                 call locate(scrx,Ngmu(ig),scale8,index1)
                 m=2
                 k=min(max(index1-(m-1)/2,1),Ngmu(ig)+1-m)
                 call polint(scrx(k),scry(k),m,scale8,qqq,dy)
                 gx8(ig)=tempgh(ig)
                 gy8(ig)=qqq
               endif
c
21           continue
c
             if(icnU.ne.430)then
               call polint(gx1,gy1,igcount,gin,qqq,dy)
               finalx1(i)=qqq
               finaly1(i)=stdscalex(i)*qqq1
             endif
c
             if(icnB.ne.430)then
               call polint(gx2,gy2,igcount,gin,qqq,dy)
               finalx2(i)=qqq
               finaly2(i)=stdscalex(i)*qqq2
             endif
c
             if(icnV.ne.430)then
               call polint(gx3,gy3,igcount,gin,qqq,dy)
               finalx3(i)=qqq
               finaly3(i)=stdscalex(i)*qqq3
             endif
c
             if(icnR.ne.430)then
               call polint(gx4,gy4,igcount,gin,qqq,dy)
               finalx4(i)=qqq
               finaly4(i)=stdscalex(i)*qqq4
             endif
c
             if(icnI.ne.430)then
               call polint(gx5,gy5,igcount,gin,qqq,dy)
               finalx5(i)=qqq
               finaly5(i)=stdscalex(i)*qqq5
             endif
c
             if(icnJ.ne.430)then
               call polint(gx6,gy6,igcount,gin,qqq,dy)
               finalx6(i)=qqq
               finaly6(i)=stdscalex(i)*qqq6
             endif
c
             if(icnH.ne.430)then
               call polint(gx7,gy7,igcount,gin,qqq,dy)
               finalx7(i)=qqq
               finaly7(i)=stdscalex(i)*qqq7
             endif
c
             if(icnK.ne.430)then
               call polint(gx8,gy8,igcount,gin,qqq,dy)
               finalx8(i)=qqq
               finaly8(i)=stdscalex(i)*qqq8
             endif
c
40         continue
c
           if(icnU.ne.430)then
             call locate(finalx1,N,rmuin,muindex)
             if(muindex.gt.0)then  !(rmuin.ge.finalx1(1))then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx1(k),finaly1(k),m,rmuin,qqq,dy)
               gscratch1(2)=qqq
             else
               gscratch1(2)=finaly1(1)
             endif
           endif
c
           if(icnB.ne.430)then
             call locate(finalx2,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx2(k),finaly2(k),m,rmuin,qqq,dy)
               gscratch2(2)=qqq
             else
               gscratch2(2)=finaly2(1)
             endif
           endif
c
           if(icnV.ne.430)then
             call locate(finalx3,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx3(k),finaly3(k),m,rmuin,qqq,dy)
               gscratch3(2)=qqq
             else
               gscratch3(2)=finaly3(1)
             endif
           endif
c
           if(icnR.ne.430)then
             call locate(finalx4,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx4(k),finaly4(k),m,rmuin,qqq,dy)
               gscratch4(2)=qqq
             else
               gscratch4(2)=finaly4(1)
             endif
           endif
c
           if(icnI.ne.430)then
             call locate(finalx5,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx5(k),finaly5(k),m,rmuin,qqq,dy)
               gscratch5(2)=qqq
             else
               gscratch5(2)=finaly5(1)
             endif
           endif
c
           if(icnJ.ne.430)then
             call locate(finalx6,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx6(k),finaly6(k),m,rmuin,qqq,dy)
               gscratch6(2)=qqq
             else
               gscratch6(2)=finaly6(1)
             endif
           endif
c
           if(icnH.ne.430)then
             call locate(finalx7,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx7(k),finaly7(k),m,rmuin,qqq,dy)
               gscratch7(2)=qqq
             else
               gscratch7(2)=finaly7(1)
             endif
           endif
c
           if(icnK.ne.430)then
             call locate(finalx8,N,rmuin,muindex)
             if(muindex.gt.0)then
               m=2
               k=min(max(muindex-(m-1)/2,1),N+1-m)
               call polint(finalx8(k),finaly8(k),m,rmuin,qqq,dy)
               gscratch8(2)=qqq
             else
               gscratch8(2)=finaly8(1)
             endif
           endif
c
c   Finally, take the final pass and interpolate between T
c
           if(icnU.ne.430)then
             call polint(tscratch,gscratch1,m,Tin,qqq,dy)
             outinty(1)=qqq
           endif
           if(icnB.ne.430)then
             call polint(tscratch,gscratch2,m,Tin,qqq,dy)
             outinty(2)=qqq
           endif
           if(icnV.ne.430)then
             call polint(tscratch,gscratch3,m,Tin,qqq,dy)
             outinty(3)=qqq
           endif
           if(icnR.ne.430)then
             call polint(tscratch,gscratch4,m,Tin,qqq,dy)
             outinty(4)=qqq
           endif
           if(icnI.ne.430)then
             call polint(tscratch,gscratch5,m,Tin,qqq,dy)
             outinty(5)=qqq
           endif
           if(icnJ.ne.430)then
             call polint(tscratch,gscratch6,m,Tin,qqq,dy)
             outinty(6)=qqq
           endif
           if(icnH.ne.430)then
             call polint(tscratch,gscratch7,m,Tin,qqq,dy)
             outinty(7)=qqq
           endif
           if(icnK.ne.430)then
             call polint(tscratch,gscratch8,m,Tin,qqq,dy)
             outinty(8)=qqq
           endif
c
 300      continue
c

c
          return
          end
c
c
c
          subroutine distorttime(Nmaxphase,Nphase,xmod,ymod,yeclipse,RV,gamma,
     $        dphase,pconj)
c
c   February 5, 2001
c
c   This routine will apply a phase shift to a light or velocity curve

c   
          implicit double precision(a-h,o-z)
c
          parameter(NNdum=720001,speed=2.99792458d5)
          dimension xmod(Nmaxphase),ymod(Nmaxphase),yeclipse(Nmaxphase)
          dimension RV(Nmaxphase)
c

          do 10 i=1,Nphase
c            write(*,*)xmod(i),yeclipse(i)
             if(yeclipse(i).gt.0.0d0)then
               if(xmod(i).lt.0.5d0)then
                 ttt=(RV(i)-gamma)/speed*(xmod(i)-pconj)
                 write(*,*)ttt
                 xmod(i)=xmod(i)+ttt
                else
                 ttt=(RV(i)-gamma)/speed*(xmod(i)-(pconj+1.0d0))
                 write(*,*)ttt
                 xmod(i)=xmod(i)+ttt
                endif
             endif
 10       continue

c
 101      format(f6.2)
          return
          end
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine getfracs(Nmaxphase,icount,
     %      fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,
     &         compfracs,dphase,Nphase,
     %         xmod,eshift,pshift,ionephase,ism1,ecc,onephase,sw26)
c
c   October 10, 2007
c
c   This is a new routine that will return the luminosity ratios in
c   each bandpass.  compfracs(i,1) will have L_2/L_1, and
c   compfracs(i,2) will have L_disk/L_tot.
c
          implicit double precision(a-h,o-z)
c
          dimension compfracs(8,2),xmod(Nmaxphase)
          dimension fracs1(Nmaxphase,3),fracs2(Nmaxphase,3)
          dimension fracs3(Nmaxphase,3),fracs4(Nmaxphase,3)
          dimension fracs5(Nmaxphase,3),fracs6(Nmaxphase,3)
          dimension fracs7(Nmaxphase,3),fracs8(Nmaxphase,3)
c
          parameter(kkk=720001)  !Nmaxphase
          dimension scratch1(kkk),scratch2(kkk),xscratch(kkk)
c

          do 10 i=1,8     !loop over filters
c
c   UPDATE October 27, 2008
c
c   First, 'finish' the fraction curves
c
            iii=icount-1
            NNN=icount
            if((ionephase.eq.0).and.(ism1.ge.1).and.(ecc.eq.0.0d0))then
              do 1 phase=180.0d0+dphase,360.0d0-dphase,dphase
                NNN=NNN+1
                if(i.eq.1)then
                  fracs1(NNN,1)=fracs1(iii,1)
                  fracs1(NNN,2)=fracs1(iii,2)
                  fracs1(NNN,3)=fracs1(iii,3)
                endif
                if(i.eq.2)then
                  fracs2(NNN,1)=fracs2(iii,1)
                  fracs2(NNN,2)=fracs2(iii,2)
                  fracs2(NNN,3)=fracs2(iii,3)
                endif
                if(i.eq.3)then
                  fracs3(NNN,1)=fracs3(iii,1)
                  fracs3(NNN,2)=fracs3(iii,2)
                  fracs3(NNN,3)=fracs3(iii,3)
                endif
                if(i.eq.4)then
                  fracs4(NNN,1)=fracs4(iii,1)
                  fracs4(NNN,2)=fracs4(iii,2)
                  fracs4(NNN,3)=fracs4(iii,3)
                endif
                if(i.eq.5)then
                  fracs5(NNN,1)=fracs5(iii,1)
                  fracs5(NNN,2)=fracs5(iii,2)
                  fracs5(NNN,3)=fracs5(iii,3)
                endif
                if(i.eq.6)then
                  fracs6(NNN,1)=fracs6(iii,1)
                  fracs6(NNN,2)=fracs6(iii,2)
                  fracs6(NNN,3)=fracs6(iii,3)
                endif
                if(i.eq.7)then
                  fracs7(NNN,1)=fracs7(iii,1)
                  fracs7(NNN,2)=fracs7(iii,2)
                  fracs7(NNN,3)=fracs7(iii,3)
                endif
                if(i.eq.8)then
                  fracs8(NNN,1)=fracs8(iii,1)
                  fracs8(NNN,2)=fracs8(iii,2)
                  fracs8(NNN,3)=fracs8(iii,3)
                endif
                xscratch(NNN)=phase/360.0d0
                xscratch(iii)=xscratch(NNN)
                iii=iii-1
 1            continue
            endif
c
c  UPDATE October 27, 2008
c
c  If the variable sw26 is positive, then reference the disk fraction at 
c  sw26.   Otherwise, take the median over the whole orbit.
c   
c
            if(sw26.le.0.0d0)then
              do 8 jj=1,NNN
                if(i.eq.1)then
                  aaa=fracs1(jj,1)
                  bbb=fracs1(jj,2)
                  ccc=fracs1(jj,3)
                endif
                if(i.eq.2)then
                  aaa=fracs2(jj,1)
                  bbb=fracs2(jj,2)
                  ccc=fracs2(jj,3)
                endif
                if(i.eq.3)then
                  aaa=fracs3(jj,1)
                  bbb=fracs3(jj,2)
                  ccc=fracs3(jj,3)
                endif
                if(i.eq.4)then
                  aaa=fracs4(jj,1)
                  bbb=fracs4(jj,2)
                  ccc=fracs4(jj,3)
                endif
                if(i.eq.5)then
                  aaa=fracs5(jj,1)
                  bbb=fracs5(jj,2)
                  ccc=fracs5(jj,3)
                endif
                if(i.eq.6)then
                  aaa=fracs6(jj,1)
                  bbb=fracs6(jj,2)
                  ccc=fracs6(jj,3)
                endif
                if(i.eq.7)then
                  aaa=fracs7(jj,1)
                  bbb=fracs7(jj,2)
                  ccc=fracs7(jj,3)
                endif
                if(i.eq.8)then
                  aaa=fracs8(jj,1)
                  bbb=fracs8(jj,2)
                  ccc=fracs8(jj,3)
                endif
                ddd=aaa+bbb+ccc
                if(aaa.eq.0.0d0)then
                  scratch1(jj)=0.0d0
                else      
                  scratch1(jj)=bbb/aaa
                endif
                if(ddd.eq.0.0d0)then
                  scratch2(jj)=0.0d0
                else      
                  scratch2(jj)=ccc/ddd
                endif
8             continue
c  
              if(icount.le.2)then
                compfracs(i,1)=0.0d0
                compfracs(i,2)=0.0d0
                go to 10
              endif

              call sort2(NNN,scratch1,scratch2)
c
              q1=dble(NNN/2)
              q2=dble(NNN)/2.0
              if(q1.eq.q2)then
                rmedian1=(scratch1(NNN/2)+scratch1(NNN/2+1))/2.0
              else
                rmedian1=scratch1(NNN/2+1)
              endif
c
              call sort2(NNN,scratch2,scratch1)
c
              q1=dble(NNN/2)
              q2=dble(NNN)/2.0
              if(q1.eq.q2)then
                 rmedian2=(scratch2(NNN/2)+scratch2(NNN/2+1))/2.0
              else
                rmedian2=scratch2(NNN/2+1)
              endif
c
              if(rmedian1.lt.0.0d0)rmedian1=0.0d0
              if(rmedian2.lt.0.0d0)rmedian2=0.0d0

              compfracs(i,1)=rmedian1
              compfracs(i,2)=rmedian2
            else
              diffmin=123456789.0d0
              do 9 jj=1,NNN
                fred=dmod(xscratch(jj)+eshift+pshift,1.0d0)
                if(fred.lt.0.0d0)fred=fred+1.0d0
                if(fred.gt.1.0d0)fred=fred-1.0d0
                diff=dabs(fred-sw26)
                if(diff.le.diffmin)then
                  diffmin=diff

                 if(i.eq.1)then
                   aaa=fracs1(jj,1)
                   bbb=fracs1(jj,2)
                   ccc=fracs1(jj,3)
                 endif
                 if(i.eq.2)then
                   aaa=fracs2(jj,1)
                   bbb=fracs2(jj,2)
                   ccc=fracs2(jj,3)
                 endif
                 if(i.eq.3)then
                   aaa=fracs3(jj,1)
                   bbb=fracs3(jj,2)
                   ccc=fracs3(jj,3)
                 endif
                 if(i.eq.4)then
                   aaa=fracs4(jj,1)
                   bbb=fracs4(jj,2)
                   ccc=fracs4(jj,3)
                 endif
                 if(i.eq.5)then
                   aaa=fracs5(jj,1)
                   bbb=fracs5(jj,2)
                   ccc=fracs5(jj,3)
                 endif
                 if(i.eq.6)then
                   aaa=fracs6(jj,1)
                   bbb=fracs6(jj,2)
                   ccc=fracs6(jj,3)
                 endif
                 if(i.eq.7)then
                   aaa=fracs7(jj,1)
                   bbb=fracs7(jj,2)
                   ccc=fracs7(jj,3)
                 endif
                 if(i.eq.8)then
                   aaa=fracs8(jj,1)
                   bbb=fracs8(jj,2)
                   ccc=fracs8(jj,3)
                 endif
c
                  ddd=aaa+bbb+ccc
                  compfracs(i,1)=bbb/aaa
                  if(ddd.eq.0.0d0)then
                    compfracs(i,2)=0.0d0
                  else
                    compfracs(i,2)=ccc/ddd
                  endif  
                endif
9             continue
            endif
c
10        continue
c
          return
          end
c 
c
c
c
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
          subroutine setdensity(fill1,omega,bdist,Q,period,density,tidephi,itide)

c     &       Q,finc,Teff1,Teff2,Period,fm,separ,
c     &       primmass,primK,primrad,ratrad,sw5,reff1,ecc,bdist,frac1,
c     $       psi0,x0,density,tidephi,itide)
c
c   October 10, 2008
c
c   This a new subroutine to set the mass ratio given the fill factor
c   and density.
c
          implicit double precision (a-h,o-z)  
c
          parameter(pie=3.14159265358979323d0)
c
           solarmass=1.9889d33
           smmks=1.9889d30    !solar mass in kg
           solarrad=6.9598d10
           Gcgs=6.67259d-8     !G in cgs
           Gmks=6.67259d-11   !G in mks
           p=period*24.0d0*3600.0d0
           gmsun=1.32712440018d20  !mks units
           ddd=density*1000.0d0
c
           istar=1
c
c   Have an ititial guess at Q
c
           Qhigh=0.1d0
           overQh=Qhigh
           Qlow=1.0d-20
           overQl=Qlow
           icount=1

999       Qmid=0.5d0*(Qhigh+Qlow)
          overQh=Qhigh
          overQm=Qmid
          overQl=Qlow

          call findL1(overQh,omega,x0h,1,bdist,tidephi,itide) 

          Rl=x0h                           !save the L1 distance
          x=fill1*x0h
          y=0.0d0
          z=0.0d0

          call POTEN(overQh,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,bdist,
     #      cox,coy,tidephi,itide)
          savepsi0=psi     ! value of the potential at x-axis 
          critpsi=psi
          call findradius(overQh,omega,savepsi0,x0h*fill1,bdist,
     &      reffhigh,tidephi,itide)


          call findL1(overQl,omega,x0l,1,bdist,tidephi,itide) 

          Rl=x0l                           !save the L1 distance
          x=fill1*x0l
          y=0.0d0
          z=0.0d0
          call POTEN(overQl,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,
     #     bdist,cox,coy,tidephi,itide)
          savepsi0=psi     ! value of the potential at x-axis 
          critpsi=psi
          call findradius(overQl,omega,savepsi0,x0l*fill1,bdist,
     &         refflow,tidephi,itide)

          call findL1(overQm,omega,x0m,1,bdist,tidephi,itide) 

          Rl=x0m                           !save the L1 distance
          x=fill1*x0m
          y=0.0d0
          z=0.0d0
          call POTEN(overQm,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,
     #     bdist,cox,coy,tidephi,itide)
          savepsi0=psi     ! value of the potential at x-axis 
          critpsi=psi
          call findradius(overQm,omega,savepsi0,x0m*fill1,bdist,
     &       reffmid,tidephi,itide)

c          write(*,455)reffhigh,reffmid,refflow
c          write(*,*)gmsun,p

          denhigh=3.0d0*pie*smmks/(reffhigh**3*gmsun*
     %         p*p*(1.0d0+Qhigh))
          denmid=3.0d0*pie*smmks/(reffmid**3*gmsun*
     %         p*p*(1.0d0+Qmid))
          denlow=3.0d0*pie*smmks/(refflow**3*gmsun*
     %         p*p*(1.0d0+Qlow))

c
          flow=denlow-ddd
          fmid=denmid-ddd
          fhigh=denhigh-ddd
c
          if(fhigh*fmid.gt.0.0d0)then
            Qhigh=Qmid
          else
            Qlow=Qmid
          endif
c
          icount=icount+1
          Q=Qmid
c
          if((dabs(fmid).gt.1.0d-12).and.(icount.le.200))go to 999
c
          if(fmid.le.1.0d-12)then
            write(2,1000)density,fill1,Q
          else
            write(2,1001)density,fill1,Q
          endif

1000      format('Info:  density = ',f7.5,' g/cc and fill1 = ',f7.5/,
     %      '       The mass ratio has been set to ',1pe16.9)
1001      format('Warning:  a density of ',f7.5,' g/cc and fill1 = ',f7.5/,
     %      '         is not possible.  The mass ratio has been set to ',
     %        1pe16.9)
          return
          end
c
c   *******
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&c
c
         function gammq(a,x)
        
         implicit double precision(a-h,o-z)
         double precision a,gammq,x
         double precision gammcf,gamser,gln

         if((x.lt.0.0d0).or.(a.le.0.0d0))pause 'bad argument in gammq'
c
         if(x.lt.a+1.0d0)then
           call gser(gamser,a,x,gln)
           gammq=1.0d0-gamser
         else
           call gcf(gammcf,a,x,gln)
           gammq=gammcf
         endif
         return
         end
c
c  ********c
c
      SUBROUTINE gcf(gammcf,a,x,gln)
      implicit double precision(a-h,o-z)
      INTEGER ITMAX
      double precision a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
      INTEGER i
      double precision an,b,c,d,del,h,gammln
      gln=gamma_log(a)

      b=x+1.0d0-a
      c=1.0d0/FPMIN
      d=1.0d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.0d0
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.0d0/d
        del=d*c
        h=h*del
        if(dabs(del-1.0d0).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=dexp(-x+a*dlog(x)-gln)*h
      return
      END
c
c
c
      SUBROUTINE gser(gamser,a,x,gln)
      implicit double precision(a-h,o-z)
      INTEGER ITMAX
      double precision a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
      INTEGER n
      double precision ap,del,sum,gamm_log
      gln=gamma_log(a)

      if(x.le.0.0d0)then
        if(x.lt.0.0d0)pause 'x < 0 in gser'
        gamser=0.0d0
        return
      endif
      ap=a
      sum=1.0d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.0d0
        del=del*x/ap
        sum=sum+del
        if(dabs(del).lt.dabs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*dexp(-x+a*dlog(x)-gln)
      return
      END
c
c
c   &&&&&&&&&&&&&&&&&&&&&&
c
          subroutine setfill(istar,Q,fill,radfill,omega,bdist,tidephi,itide)
c
c
c 
          implicit double precision(a-h,o-z)
c
          if(radfill.ge.1.0d0)return
          if(radfill.le.0.0d0)return

          overQ=Q
          if(istar.eq.2)overQ=1.0d0/Q

          fbig=1.0d0
          fsmall=0.0d-8
c
c   First, find the distance to L1.
c
          call findL1(overQ,omega,x0,1,bdist,tidephi,itide)       
          Rl=x0                           !save the L1 distance
          x=x0
          y=0.0d0
          z=0.0d0
c
c   UPDATE November 14, 2009
c
c   initialize cox and coy for the routine that computes the tidal
c   equilibrium option.
c
          cox=1.0d0
          coy=0.0d0
          call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,bdist,
     $      cox,coy,tidephi,itide)
          savepsi0=psi     ! value of the potential at x-axis 
          call findradius(overQ,omega,savepsi0,x0,bdist,reff,tidephi,itide)
          reffA=reff
c
          call findL1(overQ,omega,x0,1,bdist,tidephi,itide)       
          Rl=x0                           !save the L1 distance
          x=x0*fsmall
          y=0.0d0
          z=0.0d0
          call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,bdist,
     $       cox,coy,tidephi,itide)
          savepsi0=psi     ! value of the potential at x-axis 
          call findradius(overQ,omega,savepsi0,x0,bdist,reff,tidephi,itide)
          reffB=reff
c
          do 10 i=1,20
            fnew=0.5d0*(fbig+fsmall)
            call findL1(overQ,omega,x0,1,bdist,tidephi,itide)       
            Rl=x0                           !save the L1 distance
            x=x0*fnew
            y=0.0d0
            z=0.0d0
            call POTEN(overQ,omega,x,y,z,psi,psix,psixx,psiy,psiz,1,
     #          bdist,cox,coy,tidephi,itide)
            savepsi0=psi     ! value of the potential at x-axis 
            call findradius(overQ,omega,savepsi0,x0,bdist,reff,tidephi,itide)
            reffnew=reff
            if((reffA*radfill-reffnew).gt.0.0d0)then
              fsmall=fnew
            else
              fbig=fnew
            endif
10        continue
c
          fill=fnew
          if(istar.eq.1)write(2,665)radfill,fill
          if(istar.eq.2)write(2,666)radfill,fill

665       format('Info:  radfill1 = ',f10.8,', fill1 set to ',f10.8)
666       format('Info:  radfill2 = ',f10.8,', fill2 set to ',f10.8)
          return
          end

c
c  &&&&&&&&&&&&&&&&&&&&
c

        SUBROUTINE LPN(N,X,PN,PD)
C
C       ===============================================
C       Purpose: Compute Legendre polynomials Pn(x)
C                and their derivatives Pn'(x)
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,...)
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C       ===============================================
C
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PN(0:N),PD(0:N)
        PN(0)=1.0D0
        PN(1)=X
        PD(0)=0.0D0
        PD(1)=1.0D0
        P0=1.0D0
        P1=X
        DO 10 K=2,N
           PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
           PN(K)=PF
           IF (DABS(X).EQ.1.0D0) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
           ENDIF
           P0=P1
 10        P1=PF
        RETURN
        END
c
c
c   ##########################
c
c
          subroutine analyticscale(Q,finc,period,primmass,primK,
     %        ecc,frac1,frac2,
     &        primrad,ratrad,reff1,reff2,separ,vrot1,vrot2,gpole1,gpole2,   
     &        omega1,omega2)
c
c   December 9, 2009
c
c   This routine will return the radii and ratio of radii, which will
c   be needed for the fastanalytic routine.
c
          implicit double precision(a-h,o-z)        
c
          parameter(pie=3.14159265358979323d0)
c
          bdist=1.0d0
          solarmass=1.9889d33
          smmks=1.9889d30    !solar mass in kg
          solarrad=6.9598d10
          Gcgs=6.67259d-8     !G in cgs
          Gmks=6.67259d-11   !G in mks
          p=period*24.0d0*3600.0d0
          gmsun=1.32712440018d20  !mks units
          fincr=finc*pie/180.0d0
          gsun=2.739910d4

c
c   If frac1 and primrad are zero, then there is not enough information.
c   Write an error message and stop.
c
c   Likewise, if ratrad and frac2 are both zero, write an error message
c   and stop.
c
          if((frac1.le.0.0d0).and.(primrad.le.0.0d0))then
            write(*,*)'There is not enough information to determine ',
     &       'the primary radius.  Set frac1 > 0 or primrad > 0'        
            stop
          endif
c
          if((frac2.le.0.0d0).and.(ratrad.le.0.0d0))then
            write(*,*)'There is not enough information to determine ',
     &       'the secondary radius.  Set frac2 > 0 or ratrad > 0'        
            stop
          endif
c
c   primmass > 0 only, set the separation
c
          if((primmass.gt.0.0d0).and.(primK.le.0.0d0))then
            rmass=primmass*solarmass
            separ=(Gcgs*p*p*(1.0d0+Q)*rmass
     $          /(4.0d0*pie*pie))**(1.d0/3.0d0)/solarrad
            separ=(gmsun*p*p*primmass*(1.0d0+Q)/(4.0d0*pie*pie))**(1.0d0/3.0d0)
            separ=separ/solarrad*100.0d0
            write(2,1199)primmass,separ
            go to 99
          endif
c
c
c   primK > 0 only, set the separation
c
          if((primmass.le.0.0d0).and.(primK.gt.0.0d0))then
            efact=dsqrt(1.0d0-ecc*ecc)
            vkcgs=primK*100000.0d0*efact
            separ=vkcgs*p*(1.0d0+Q)/(2.0*pie*dsin(fincr)*Q)/solarrad
            write(2,1188)primK,separ
            go to 99
          endif
c
c   primmass > 0 and primK > 0.  Solve for Q and separ
c
          if((primmass.gt.0.0d0).and.(primK.gt.0.0d0))then
            rmass=primmass*solarmass
            efact=dsqrt(1.0d0-ecc*ecc)
            dqhi=7.
            dqlo=-7.
            do 555 kk=1,35
              Qhigh=10.0d0**dqhi
              Qlow=10.0d0**dqlo
              Qmid=10.0**((dqhi+dqlo)*0.5d0)
              aa=vfcn(Qlow,period,finc,primmass,primK,ecc)
              bb=vfcn(Qhigh,period,finc,primmass,primK,ecc)
              cc=vfcn(Qmid,period,finc,primmass,primK,ecc)
              if(aa*cc.lt.0.0d0)then
                dqhi=(dqhi+dqlo)*0.5d0
              else
                dqlo=(dqhi+dqlo)*0.50
              endif
 555        continue
c
            do 556 kk=1,25
              Qmid=(Qhigh+Qlow)*0.5d0
              aa=vfcn(Qlow,period,finc,primmass,primK,ecc)
              bb=vfcn(Qhigh,period,finc,primmass,primK,ecc)
              cc=vfcn(Qmid,period,finc,primmass,primK,ecc)
              if(aa*cc.lt.0.0d0)then
                Qhigh=(Qhigh+Qlow)*0.5d0
              else
                Qlow=(Qhigh+Qlow)*0.50
              endif
 556         continue
c
             Q=Qmid
             separ=(gmsun*p*p*primmass*
     &              (1.0d0+Q)/(4.0d0*pie*pie))**(1.0d0/3.0d0)
             separ=separ/solarrad*100.0d0              
c
c   Use the formula separ = coef*(perid*period*total_mass)**(1/3) to
c   solve for the total mass in solar masses.  The separation is
c   entered in solar masses, so (R_sun/coef)**3=7.737294491.
c
             total_mass=primmass*(1.0d0+Q)
c
             ppp=period*24.0d0
             separ1=(total_mass*ppp*ppp/7.737294491d0)**(1.0d0/3.0d0)
             write(2,1177)primmass,primK,Q,separ
             go to 99
           endif
c
1177      format(/'fast analytic mode:',
     &        /'Info:  M_1 fixed at ',f9.6,' solar masses and ',
     &        'K_1 fixed at ',f11.6,' km/sec.',/
     &        'Q is set to ',f12.7,' and the separation is ',f13.7,
     &        ' solar radii')
1199      format(/'fast analytic mode:',
     %            /'Info:  M_1 fixed at ',f9.6,' solar masses.  The',  
     %        ' separation',/ '       has been set to ',f13.7,
     $         ' solar radii')  
1188      format(/'fast analytic mode:',
     #         /'Info:  K_1 fixed at ',f14.8,' km/sec.  The',  
     %        ' separation',/ '       has been set to ',f13.7,
     $         ' solar radii')  
c

99        continue
c
c   We should now have the separation in solar radii.  Figure out if
c   there is enough information to get the radii from the input values
c   of frac1, frac2, primrad, and ratrad.  If frac1 is specified, then
c   find the primary radius from that.  If frac2 is specified, then
c   find the secondary radius from that.
c
          if(frac1.gt.0.0d0)then
            reff1=frac1
            if(frac2.gt.0.0d0)then
              reff2=frac2
              go to 199
            endif
            if(ratrad.gt.0.0d0)then
              reff2=reff1/ratrad
              go to 199
            endif
          endif
c
          if(primrad.gt.0.0d0)then
            reff1=primrad/separ
            if(frac2.gt.0.0d0)then
              reff2=frac2
              go to 199
            endif
            if(ratrad.gt.0.0d0)then
              reff2=reff1/ratrad
              go to 199
            endif
          endif
c
199       continue
c
         if(ratrad.eq.0.0d0)ratrad=reff1/reff2
c
c   We have to figure out the gravities and the rotational velocities.
c
          fpsq=4.0d0*pie*pie
          ppp=period*24.0d0               ! period in hours
          coef=3.518847d10
          coef3=4.35713636d31
          sifinc=dsin(fincr)
          smet=separ*6.9598d8
          total_mass=smet*smet*smet*fpsq/(period*86400.0d0)**2/gmsun
          rM1=total_mass/(1.0d0+Q)
          rM2=Q*rM1           
          R1=reff1*separ
          R2=reff2*separ
c
          gpole1=gsun*rM1/(R1*R1)  
          gpole2=gsun*rM2/(R2*R2)  
c
          gscale1=27397.726d0*rM1/(separ*separ*bdist*bdist)
          gscale2=27397.726d0*rM2/(separ*separ*bdist*bdist)
c
          fact=solarrad*2.0d0*pie/86400.0d0/1.0d5  !50.613093d0
          vrot1=fact*R1/period*sifinc
          vrot2=fact*R2/period*sifinc

          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine getanalyticint(maxlines,maxmu,Nlines,
     &      atmT,atmg,atmmu,Nmu,
     &      atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &      Tmax,Tmin,gmax,gmin,gpole,darkint,tpole,
     %      dwavex,dwavey,ilaw,iatm,istar,wave,reff)
c
c  December 9, 2009
c
c  This subroutine will evaluate the integral:
c
c  dint = int^1_0 (I(T,g,mu)*mu*du)
c
c  The dint values will go into the reference fluxes for the fast
c  analytic mode.  
c
c  if iatm=0, then use black bodies
c
c
          implicit double precision(a-h,o-z)

          parameter(pie=3.14159265358979323d0)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines),outinty(8),
     #       darkint(8),summ(8),rnorm(8)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
          dimension dwavex(8,2),dwavey(8,2),wave(8)
c
          if(iatm.gt.0)then
c
c   Find the rough place in the atmosphere table.
c
            Tin=tpole
            gin=dlog10(gpole)
            call locate(atmT,Nlines,Tin,indexT)
            itguess=indexT
            imuguess=1
c
            do 1 i=1,8
              summ(i)=0.0d0
 1          continue
c
            do 4 i=100,1,-1
              rmuin=dble(i-1)/100.0d0+0.005d0
              call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &           atmT,atmg,atmmu,Nmu,
     &           atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &           atmint7,atmint8,
     &           Tmax,Tmin,gmax,gmin,outinty,
     %           icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,itguess,imuguess,
     *           dwavex,dwavey,ilaw,iatm,istar)
c
              do 88 k=1,8
                if(i.eq.100)rnorm(k)=outinty(k)
                summ(k)=summ(k)+2.0d0*pie*0.01*outinty(k)*rmuin
                darkint(k)=summ(k)*reff*reff
 88           continue
 4          continue
            return
          endif
c
          if(iatm.le.0)then 
            do 200 jj=1,8
              flimbx=dwavex(jj,istar)
              flimby=dwavey(jj,istar)
              ddint=pie*(1.0d0-flimbx/3.0d0)
              if(ilaw.eq.2)then
                ddint=pie*(1.0d0-flimbx/3.0d0+2.0d0*flimby/9.0d0)
              endif
              if(ilaw.eq.3)then
                ddint=pie*(1.0d0-flimbx/3.0d0-flimby/5.0d0)
              endif
              if(ilaw.eq.4)then
                ddint=pie*(1.0d0-flimbx/3.0d0-flimby/6.0d0)           
              endif
c
              wavemu=wave(jj)/10000.0d0
              c1=3.74185
              c2=14.3883
              tkkelv=tpole/1000.0d0
              C3 = C2/(wavemu*tkkelv)
              darkint(jj)=ddint*C1/(dexp(c3)-1.0d0)/wavemu**5*reff*reff
200         continue
          endif
c
          return
c
          end
c
c  &&&&&&&&&&&&&
c
          subroutine fastanalytic(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $       ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     &       ymodd,RV1,RV2,drv1,drv2,obsparm,NRVphase,xRVmod,
     &       fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,
     $       frac1,frac2,primrad,ratrad,period,T0,primmass,
     #       primK,ecc,argper,omega1,omega2,
     $       finc,ilaw,dwavex,dwavey,bigI,bigbeta,sw29,sw30,ialign,ikeep,
     &       pshift,reff1,reff2,darkint1,darkint2,
     #       idark1,idark2,dphase,iRVfilt,vrot1,vrot2,Q,separ,gamma,isw27,
     &       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     &       timearray,sw9,sw23,sw24,isw7,SA3)
c
c   UPDATE December 9, 2009
c
c   This is a new routine that will compute light curves using the 
c   analytic expressions of Giminez.  This routine tries to do as little
c   as possible when computing the curves.
c
c
          implicit double precision(a-h,o-z)
c
          dimension dwavex(8,2),dwavey(8,2)
          dimension xmod(Nmaxphase),ymodU(Nmaxphase),ymodB(Nmaxphase),
     $      ymodV(Nmaxphase),ymodR(Nmaxphase),ymodI(Nmaxphase),
     $      ymodJ(Nmaxphase),ymodH(Nmaxphase),ymodK(Nmaxphase),
     &      ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      RV1(Nmaxphase),RV2(Nmaxphase),dRV1(Nmaxphase),
     %      dRV2(Nmaxphase),ymods3(Nmaxphase)
          dimension xRVmod(Nmaxphase)
          dimension fracs1(Nmaxphase,3),fracs2(Nmaxphase,3)
          dimension fracs3(Nmaxphase,3),fracs4(Nmaxphase,3)
          dimension fracs5(Nmaxphase,3),fracs6(Nmaxphase,3)
          dimension fracs7(Nmaxphase,3),fracs8(Nmaxphase,3)
          dimension refflux1(8),refflux2(8),third(8)
          dimension corr1(8),corr2(8)
          dimension darkint1(8),darkint2(8),gimvel(8)
          dimension obsparm(17)

          dimension timearray(900000)

          parameter(pie=3.141592653589793d0)
c
          tstep=sw9
          tstart=sw23
          tstop=sw24



          if(isw7.eq.2)call filltime(ntime,timearray,tstart,tstop,tstep)

          if(ialign.le.0)then
            bigI=finc
            bigbeta=0.0d0
          endif
          do 9 kk=1,8
            refflux1(kk)=darkint1(kk)
            refflux2(kk)=darkint2(kk)
            corr1(kk)=0.0d0
            corr2(kk)=0.0d0
            if(idark2.ge.1)darkint2(kk)=0.0d0
            third(kk)=0.0d0
9         continue
c
          if(ecc.le.0.0d0)argper=90.0d0
          bdist=1.0d0
          dpsave=dphase
          pstart=0.0d0
          pstop=360.0d0-dphase
          pconj=pie
          pconj2=0.0d0
          pstep=dphase
          pstartout=0.0d0
          pstopout=0.0d0
          if(ecc.gt.0.0d0)then
            pstartout=0.0d0
            pstopout=360.0d0-dphase
          endif
          argrad=argper*pie/180.0d0
          fincr=finc*pie/180.0d0
          dflux=0.0d0
c
          icount=0
c
c  Here is some code adapted from Wilson-Devinney to keep track of
c  phases needed for eccentric orbits.
c
          trc=0.5d0*pie-argrad
 1139     if(trc.lt.0.d0) trc=trc+2.0d0*pie
          if(trc.lt.0.d0) goto 1139
 1140     if(trc.ge.2.0d0*pie) trc=trc-2.0d0*pie
          if(trc.ge.2.0d0*pie) goto 1140
          htrc=0.5d0*trc
          if(dabs(0.5*pie-htrc).lt.7.d-6) goto 11101
          if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 11101
          ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
          goto 11103
11101     ecan=pie
11103     xmc=ecan-ecc*dsin(ecan)
          if(xmc.lt.0.d0) xmc=xmc+2.0d0*pie
          phper=1.d0-xmc/(2.0d0*pie)
          pconj=(xmc+argrad)/(2.0d0*pie)-0.25d0
          if(pconj.gt.1.0d0)pconj=pconj-1.0d0
c
c   Make this new block to compute the conjunction phase for star 2.
c
          trc=0.5d0*pie-argrad+pie
 3139     if(trc.lt.0.d0) trc=trc+2.0d0*pie
          if(trc.lt.0.d0) goto 3139
 3140     if(trc.ge.2.0d0*pie) trc=trc-2.0d0*pie
          if(trc.ge.2.0d0*pie) goto 3140
          htrc=0.5d0*trc
          if(dabs(0.5*pie-htrc).lt.7.d-6) goto 31101
          if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 31101
          ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
          goto 31103
31101     ecan=pie
31103     xmc=ecan-ecc*dsin(ecan)
          if(xmc.lt.0.d0) xmc=xmc+2.0d0*pie
          phper2=1.d0-xmc/(2.0d0*pie)
          pconj2=(xmc+argrad)/(2.0d0*pie)-0.25d0
          if(pconj2.gt.1.0d0)pconj2=pconj2-1.0d0
c
          eshift=0.0d0
          if(ikeep.eq.1)eshift=phper+pconj-0.5d0
          if(ikeep.eq.2)eshift=phper2+pconj2
          if(ecc.eq.0.0d0)then
            eshift=0.0d0
            pconj=0.0d0
          endif
c
          ttiny=0.0d0
          if(ecc.le.0.0d0)ttiny=1.0d-6
c
          if((isw7.eq.2).and.(ecc.gt.0.0d0))then
            pstartout=360.0d0*(timearray(1)-T0)/period
            pstopout=360.0d0*(timearray(ntime)-T0)/period
            dphase=360.0d0*tstep/period
            ism1=0
          endif

          do 999 phaseout=pstartout,pstopout+ttiny,dphase
c
            if(ecc.gt.0.0d0)then
              em=phaseout*pie/180.0d0
              call getE(em,ecc,bigE)
11107         if(bigE.lt.0.0)then
                bigE=bigE+2.0d0*pie
                go to 11107
              endif
22207         if(bigE.gt.2.0d0*pie)then
                bigE=bigE-2.0d0*pie
                go to 22207
              endif
              rnu=2.0d0*datan(dsqrt((1.0d0+ecc)/(1.0d0-ecc))*dtan(bigE/2.0d0))
11105         if(rnu.lt.0.0)then
                rnu=rnu+2.0d0*pie
                go to 11105
              endif
11106         if(rnu.gt.2.0d0*pie)then
                rnu=rnu-2.0d0*pie
                go to 11106
              endif
              bdist=(1.0d0-ecc*dcos(bigE))
              pstart=dmod(rnu*180.0d0/pie+argper+90.0d0,360.0d0)
              pstop=dmod(rnu*180.0d0/pie+argper+90.0d0,360.0d0)
            endif
            pstep=dphase
c
c            phaser=pstart*pie/180.0d0
c            delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
c            delta=bdist*dsqrt(delta)
c

          ipstep=0
          if((isw7.eq.2).and.(ecc.eq.0.0d0))then
            pstart=360.0d0*(timearray(1)-T0)/period
            pstop=360.0d0*(timearray(ntime)-T0)/period
            pstep=360.0d0*tstep/period
            ism1=0
          endif

            do 10 phasein=pstart,pstop+ttiny,pstep
              icount=icount+1
              phase=dmod(phasein,360.0d0)
              if(phase.lt.0.0d0)phase=phase+360.0d0
              tphase=phase+360.0d0*(pshift+eshift)
              extphase=phase
              if((ecc.gt.0.0d0).or.(pshift.ne.0.0d0))then
                if((ecc.gt.0.0d0).and.(ism1.eq.0))emphase=180.0d0*em/pie
                if((ecc.gt.0.0d0).and.(ism1.gt.0))then
                  if(ipstep.eq.1)emphase=180.0d0*em/pie
                  if(ipstep.eq.2)emphase=180.0d0*emnew/pie
                endif
                tshift=pshift+eshift
                extphase=emphase+360.0d0*tshift
                if(ikeep.eq.1)extphase=emphase+360.0d0*(tshift-pconj)
                if(ikeep.eq.2)extphase=emphase+360.0d0*(tshift-pconj2)
              endif
              dummyphase=dmod(phase+180.0d0,360.0d0)  ! this is for star 2
c
              phaser=phase*pie/180.0d0
              delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
              delta=bdist*dsqrt(delta)
c
              call analyticg(isw27,ilaw,dwavex,dwavey,delta,reff1,ratrad,
     #              refflux1,refflux2,phaser,pconj,pconj2,gimvel,fincr,vrot1,
     #              omega1,period,separ,Q,ecc,bigI,bigbeta,Neclipse1,1,
     &              corr1,corr2)

              phaser=phase*pie/180.0d0
              delta=(cos(fincr)**2+(sin(fincr)*sin(phaser))**2)
              delta=bdist*dsqrt(delta)

              dummyrat= 1.0d0/ratrad
              call analyticg(isw27,ilaw,dwavex,dwavey,delta,reff2,dummyrat,
     #              refflux1,refflux2,phaser,pconj,pconj2,gimvel,fincr,vrot2,
     #              omega1,period,separ,Q,ecc,bigI,bigbeta,Neclipse1,2,
     &              corr1,corr2)

              do 5432 ll=1,8
                if(ll.eq.1)ymodU(icount)=third(1)+darkint1(1)+darkint2(1)+corr1(1)+corr2(1)
                if(ll.eq.2)ymodB(icount)=third(2)+darkint1(2)+darkint2(2)+corr1(2)+corr2(2)
                if(ll.eq.3)ymodV(icount)=third(3)+darkint1(3)+darkint2(3)+corr1(3)+corr2(3)
                if(ll.eq.4)ymodR(icount)=third(4)+darkint1(4)+darkint2(4)+corr1(4)+corr2(4)
                if(ll.eq.5)ymodI(icount)=third(5)+darkint1(5)+darkint2(5)+corr1(5)+corr2(5)
                if(ll.eq.6)ymodJ(icount)=third(6)+darkint1(6)+darkint2(6)+corr1(6)+corr2(6)
                if(ll.eq.7)ymodH(icount)=third(7)+darkint1(7)+darkint2(7)+corr1(7)+corr2(7)
                if(ll.eq.8)ymodK(icount)=third(8)+darkint1(8)+darkint2(8)+corr1(8)+corr2(8)

c                if(ll.eq.1)ymodU(icount)=ymodU(icount)+darkint1(1)+darkint2(1)
c                if(ll.eq.2)ymodB(icount)=ymodB(icount)+darkint1(2)+darkint2(2)
c                if(ll.eq.3)ymodV(icount)=ymodV(icount)+darkint1(3)+darkint2(3)
c                if(ll.eq.4)ymodR(icount)=ymodR(icount)+darkint1(4)+darkint2(4)
c                if(ll.eq.5)ymodI(icount)=ymodI(icount)+darkint1(5)+darkint2(5)
c                if(ll.eq.6)ymodJ(icount)=ymodK(icount)+darkint1(6)+darkint2(6)
c                if(ll.eq.7)ymodH(icount)=ymodH(icount)+darkint1(7)+darkint2(7)
c                if(ll.eq.8)ymodK(icount)=ymodK(icount)+darkint1(8)+darkint2(8)

c
                call getrefvel(1,omega1,phase,finc,
     %                  Q,separ,period,gamma,vel1,ecc,argrad,isw27,gimvel,
     #                     iRVfilt)
                call getrefvel(2,omega2,dummyphase,finc,
     %                  Q,separ,period,gamma,vel2,ecc,argrad+pie,isw27,gimvel,
     #                     iRVfilt)
                RV1(icount)=vel1
                RV2(icount)=vel2
                if(ll.eq.iRVfilt)then
                  ymods1(icount)=darkint1(ll)+corr1(ll)   !refflux1(ll)
                  ymods2(icount)=darkint2(ll)+corr2(ll)   !refflux2(ll)
                endif
c
                if(ll.eq.1)then
                  fracs1(icount,1)=refflux1(ll)
                  fracs1(icount,2)=refflux2(ll)
                  fracs1(icount,3)=dflux
                endif
                if(ll.eq.2)then
                  fracs2(icount,1)=refflux1(ll)
                  fracs2(icount,2)=refflux2(ll)
                  fracs2(icount,3)=dflux
                endif
                if(ll.eq.3)then
                  fracs3(icount,1)=refflux1(ll)
                  fracs3(icount,2)=refflux2(ll)
                  fracs3(icount,3)=dflux
                endif
                if(ll.eq.4)then
                  fracs4(icount,1)=refflux1(ll)
                  fracs4(icount,2)=refflux2(ll)
                  fracs4(icount,3)=dflux
                endif
                if(ll.eq.5)then
                  fracs5(icount,1)=refflux1(ll)
                  fracs5(icount,2)=refflux2(ll)
                  fracs5(icount,3)=dflux
                endif
                if(ll.eq.6)then
                  fracs6(icount,1)=refflux1(ll)
                  fracs6(icount,2)=refflux2(ll)
                  fracs6(icount,3)=dflux
                endif
                if(ll.eq.7)then
                  fracs7(icount,1)=refflux1(ll)
                  fracs7(icount,2)=refflux2(ll)
                  fracs7(icount,3)=dflux
                endif
                if(ll.eq.8)then
                  fracs8(icount,1)=refflux1(ll)
                  fracs8(icount,2)=refflux2(ll)
                  fracs8(icount,3)=dflux
                endif

 5432         continue
c
              if(isw7.le.1)xmod(icount)=phase/360.0d0   
              if(ecc.gt.0.0d0)xmod(icount)=em/(2.0d0*pie)
              xRVmod(icount)=xmod(icount)          
              if(isw7.ge.2)xmod(icount)=timearray(icount)
c
  10        continue
c
 999      continue    ! continue the big loop if eccentric
c
          Nphase=icount
          NRVphase=icount
c
          if(isw7.le.1)then
            if((ecc.gt.0.0d0).or.(pshift.ne.0.0d0))then
              tshift=pshift+eshift
c
c   UPDATE September 10, 2001
c
c   Add the if-then clauses.
c  
              if((ecc.gt.0.0d0).and.(ikeep.eq.1))tshift=pshift+eshift-pconj
              if((ecc.gt.0.0d0).and.(ikeep.eq.2))tshift=pshift+eshift-pconj2
              if(icnU.ne.430)call shiftlc(Nmaxphase,Nphase,xmod,ymodU,tshift,0)
              if(icnB.ne.430)call shiftlc(Nmaxphase,Nphase,xmod,ymodB,tshift,0)
              if(icnV.ne.430)call shiftlc(Nmaxphase,Nphase,xmod,ymodV,tshift,0)
              if(icnR.ne.430)call shiftlc(Nmaxphase,Nphase,xmod,ymodR,tshift,0)
              if(icnI.ne.430)call shiftlc(Nmaxphase,Nphase,xmod,ymodI,tshift,0)
              if(icnJ.ne.430)call shiftlc(Nmaxphase,Nphase,xmod,ymodJ,tshift,0)
              if(icnH.ne.430)call shiftlc(Nmaxphase,Nphase,xmod,ymodH,tshift,0)
              if(icnK.ne.430)call shiftlc(Nmaxphase,Nphase,xmod,ymodK,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymods1,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymods2,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymods3,tshift,0)
              call shiftlc(Nmaxphase,Nphase,xmod,ymodd,tshift,1)
              if(icnRV1.ne.430)call shiftlc(Nmaxphase,Nphase,xRVmod,
     &             RV1,tshift,0)
              if(icnRV2.ne.430)call shiftlc(Nmaxphase,Nphase,xRVmod,
     $             dRV1,tshift,0)
              if(icnRV1.ne.430)call shiftlc(Nmaxphase,Nphase,xRVmod,
     &             RV2,tshift,0)
              if(icnRV2.ne.430)call shiftlc(Nmaxphase,Nphase,xRVmod,
     &            dRV2,tshift,1)
            endif
          endif
c
c   If requested, bin the light and velocity curves.  sw29 will be
c   the binsize for the photometry, in minutes, and sw30 will be the bin
c   size for the RV curves, in minutes
c
          if(sw29.gt.0.0d0)then
            if(isw7.le.1)then
              if(icnU.ne.430)call binlc(Nmaxphase,Nphase,xmod,ymodU,period,
     #        dphase,sw29)
            else
              if(icnU.ne.430)call binlctime(Nmaxphase,Nphase,xmod,ymodU,period,
     #        tstep,sw29)
            endif
          endif
          if(sw29.gt.0.0d0)then
            if(isw7.le.1)then
              if(icnB.ne.430)call binlc(Nmaxphase,Nphase,xmod,ymodB,period,
     #        dphase,sw29)
            else
              if(icnB.ne.430)call binlctime(Nmaxphase,Nphase,xmod,ymodB,period,
     #        tstep,sw29)
            endif
          endif
          if(sw29.gt.0.0d0)then
            if(isw7.le.1)then
              if(icnV.ne.430)call binlc(Nmaxphase,Nphase,xmod,ymodV,period,
     #        dphase,sw29)
            else
              if(icnV.ne.430)call binlctime(Nmaxphase,Nphase,xmod,ymodV,period,
     #        tstep,sw29)
            endif
          endif
          if(sw29.gt.0.0d0)then
            if(isw7.le.1)then
              if(icnR.ne.430)call binlc(Nmaxphase,Nphase,xmod,ymodR,period,
     #        dphase,sw29)
            else
              if(icnR.ne.430)call binlctime(Nmaxphase,Nphase,xmod,ymodR,period,
     #        tstep,sw29)
            endif
          endif
          if(sw29.gt.0.0d0)then
            if(isw7.le.1)then
              if(icnI.ne.430)call binlc(Nmaxphase,Nphase,xmod,ymodI,period,
     #        dphase,sw29)
            else
              if(icnI.ne.430)call binlctime(Nmaxphase,Nphase,xmod,ymodI,period,
     #        tstep,sw29)
            endif
          endif
          if(sw29.gt.0.0d0)then
            if(isw7.le.1)then
              if(icnJ.ne.430)call binlc(Nmaxphase,Nphase,xmod,ymodJ,period,
     #        dphase,sw29)
            else
              if(icnJ.ne.430)call binlctime(Nmaxphase,Nphase,xmod,ymodJ,period,
     #        tstep,sw29)
            endif
          endif
          if(sw29.gt.0.0d0)then
            if(isw7.le.1)then
              if(icnH.ne.430)call binlc(Nmaxphase,Nphase,xmod,ymodH,period,
     #        dphase,sw29)
            else
              if(icnH.ne.430)call binlctime(Nmaxphase,Nphase,xmod,ymodH,period,
     #        tstep,sw29)
            endif
          endif
          if(sw29.gt.0.0d0)then
            if(isw7.le.1)then
              if(icnK.ne.430)call binlc(Nmaxphase,Nphase,xmod,ymodK,period,
     #        dphase,sw29)
            else
              if(icnK.ne.430)call binlctime(Nmaxphase,Nphase,xmod,ymodK,period,
     #        tstep,sw29)
            endif
          endif
          if(sw29.gt.0.0d0)then
             if(isw7.le.1)then
               call binlc(Nmaxphase,Nphase,xmod,ymods1,period,
     #         dphase,sw29)
             else
               call binlc(Nmaxphase,Nphase,xmod,ymods1,period,
     #         dphase,sw29)
             endif
          endif
          if(sw29.gt.0.0d0)then
             if(isw7.le.1)then
               call binlc(Nmaxphase,Nphase,xmod,ymods2,period,
     #         dphase,sw29)
             else
               call binlctime(Nmaxphase,Nphase,xmod,ymods2,period,
     #         tstep,sw29)
             endif
          endif
          if(sw29.gt.0.0d0)then
             if(isw7.le.1)then
               call binlc(Nmaxphase,Nphase,xmod,ymods3,period,
     #         dphase,sw29)
             else
               call binlctime(Nmaxphase,Nphase,xmod,ymods3,period,
     #         tstep,sw29)
             endif
          endif
          if(sw30.gt.0.0d0)then
            if(icnRV1.ne.430)call binlc(Nmaxphase,NRVphase,xRVmod,RV1,period,
     #       dphase,sw30)
          endif
          if(sw30.gt.0.0d0)then
            if(icnRV2.ne.430)call binlc(Nmaxphase,NRVphase,xRVmod,RV2,period,
     #       dphase,sw30)
          endif
          if(sw30.gt.0.0d0)then
            if(icnRV1.ne.430)call binlc(Nmaxphase,NRVphase,xRVmod,dRV1,period,
     #       dphase,sw30)
          endif
          if(sw30.gt.0.0d0)then
            if(icnRV2.ne.430)call binlc(Nmaxphase,NRVphase,xRVmod,dRV2,period,
     #       dphase,sw30)
          endif
c
          rmedian=refflux1(1)+refflux2(1)+third(1)
          rmedian=SA3*rmedian/(1.0d0-SA3)
          do 1234 jj=1,Nphase
            ymodU(jj)=ymodU(jj)+rmedian
            ymods3(jj)=rmedian            
1234      continue
c          
          return
          end
c
c   &&%&%&%&%&%&%&%&%@@#@$@#@$#$#@@@
c
          subroutine addpad(Nphase,xmod,ymod,xpad,ypad)
c
c   This routine will return a padded pair of arrays with phases going
c   from -1 to 2.
c
           implicit double precision(a-h,o-z)

           dimension xmod(Nphase),ymod(Nphase),xpad(Nphase*3),ypad(Nphase*3)
c
           icount=0
           do 10 i=1,Nphase
             icount=icount+1
             xpad(icount)=xmod(i)-1.0d0
             ypad(icount)=ymod(i)
10         continue

           do 20 i=1,Nphase
             icount=icount+1
             xpad(icount)=xmod(i)
             ypad(icount)=ymod(i)
20         continue

           do 30 i=1,Nphase
             icount=icount+1
             xpad(icount)=xmod(i)+1.0d0
             ypad(icount)=ymod(i)
30         continue
c
           return
           end     
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
          subroutine removepad(Nphase,xmod,ymod,xpad,ypad)
c
c   This routine will remove the pad.
c
           implicit double precision(a-h,o-z)

           dimension xmod(Nphase),ymod(Nphase),xpad(Nphase*3),ypad(Nphase*3)
c
           icount=0
           do 10 i=Nphase+1,2*Nphase
             icount=icount+1
             ymod(icount)=ypad(i)
10         continue
c
           return
           end     
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
           subroutine binphi(Nphi,xmod,ymod,rowsum,Nbet)
c
c   This subroutine will bin a light curve in phase, using a binsize given
c   by sw29, where the units are minutes.
c
           implicit double precision (a-h,o-z)
c
           dimension xmod(Nphi),ymod(Nphi)
           dimension yinter(900000),y2(900000)
           dimension xpad(900000),ypad(900000)
           parameter(pie=3.1415926535897932384d0,twopie=pie+pie)
c
           call sort2(Nphi,xmod,ymod)
           call addphipad(Nphi,xmod,ymod,xpad,ypad)
           Mphase=Nphi*3
c
c   We will interpolate the model, and given each value in xmod, figure
c   out the range of phase needed, and average the interpolated y-values
c
           m=3
c
           iover=8
           NN=4*Nbet
           NN=NN*iover
           dphi=twopie/dble(NN)
           summ=0.0d0
c
           do 99 i=1,NN
             phi=-0.5d0*dphi+dphi*dble(i)
             xxx=phi
             call hunt(xpad,Mphase,xxx,index)
             k=min(max(index-(m-1)/2,1),Mphase+1-m)
c
c               call polint(xpad(k),ypad(k),m,xxx,qqqa,dy)
             call intep(xxx,qqqa,xpad(k),ypad(k),m,ier)

c
c             qqqa=((xpad(k+1)-xxx)*ypad(k)+(xxx-xpad(k))*ypad(k+1))
c     #           /(xpad(k+1)-xpad(k))
c
             summ=summ+qqqa
99         continue
c
c           kcount=0
c           summ=0.0d0
c           jlo=1
cc
c           a=0.0d0
c           b=2.0d0*pie
c           m=2
c           do 9 j=1,NN
c             if(j.eq.1)then
c               xxx=a
c               call hunt(xpad,Mphase,xxx,index)
c               k=min(max(index-(m-1)/2,1),Mphase+1-m)
cc
cc               call polint(xpad(k),ypad(k),m,xxx,qqqa,dy)
cc
c               qqqa=((xpad(k+1)-xxx)*ypad(k)+(xxx-xpad(k))*ypad(k+1))
c     #           /(xpad(k+1)-xpad(k))
cc
c               xxx=b
c               call hunt(xpad,Mphase,xxx,index)
c               k=min(max(index-(m-1)/2,1),Mphase+1-m)
cc
cc               call polint(xpad(k),ypad(k),m,xxx,qqqb,dy)
cc
c               qqqb=((xpad(k+1)-xxx)*ypad(k)+(xxx-xpad(k))*ypad(k+1))
c     #           /(xpad(k+1)-xpad(k))
cc
c               summ=0.5d0*(b-a)*(qqqa+qqqb)
c             else
c               it=2**(j-2)
c               tnm=dble(it)
c               del=(b-a)/tnm
c               x=a+0.5d0*del
c               xxx=x
c               s=0.0d0
c               do 11 kk=1,it
c                 kcount=kcount+1
c                 call hunt(xpad,Mphase,xxx,index)
c                 k=min(max(index-(m-1)/2,1),Mphase+1-m)
cc
cc                 call polint(xpad(k),ypad(k),m,xxx,qqq,dy)
cc
c                  qqq=((xpad(k+1)-xxx)*ypad(k)+(xxx-xpad(k))*ypad(k+1))
c     #              /(xpad(k+1)-xpad(k))
cc
c                 s=s+qqq
c                 x=x+del
c                 xxx=x
c11             continue
c               summ=0.5d0*(summ+(b-a)*s/tnm)
cc             endif
cc
c9           continue

           rowsum=summ/dble(iover)

10         continue
c
           return
           end
c
c   &&%&%&%&%&%&%&%&%@@#@$@#@$#$#@@@
c
          subroutine addphipad(Nphase,xmod,ymod,xpad,ypad)
c
c   This routine will return a padded pair of arrays with phases going
c   from -1 to 2.
c
           implicit double precision(a-h,o-z)
           parameter(pie=3.14159265358979323d0)

c           dimension xmod(Nphase),ymod(Nphase),xpad(Nphase*3),ypad(Nphase*3)
           dimension xmod(2000),ymod(2000),xpad(9999),ypad(9999)
c
           icount=0
           do 10 i=1,Nphase
             icount=icount+1
             xpad(icount)=xmod(i)-2.0d0*pie
             ypad(icount)=ymod(i)
10         continue

           do 20 i=1,Nphase
             icount=icount+1
             xpad(icount)=xmod(i)
             ypad(icount)=ymod(i)
20         continue

           do 30 i=1,Nphase
             icount=icount+1
             xpad(icount)=xmod(i)+2.0d0*pie
             ypad(icount)=ymod(i)
30         continue
c
           return
           end
c
c   &&&&&&&&&&&&&&&&&&&&&&&@#$%$$#$@@&@&@&@@@@@@@&&&&&&&&
c
           subroutine edgecor(Nrow,xrow,yrow,
     #           phor1,phor2,Nbet,corr)
c
           implicit double precision (a-h,o-z)
c
           dimension xrow(2000),yrow(2000),xpad(9999),ypad(9999)
c
           parameter(pie=3.1415926535897932384d0)
c
           corr=0.0d0
           corr1=0.0d0
           corr2=0.0d0
c
c    First, sort the array and pad it (make the angle go from
c    -2pi to 4pi to avoid edge effects
c
           call sort2(Nrow,xrow,yrow)
           call addphipad(Nrow,xrow,yrow,xpad,ypad)
           Npad=Nrow*3
c
c   Find the phi value of the visible pixel nearest to the
c   first limb (phor1)
c
           dphi=pie/dble(Nrow)
           halfdphi=0.5d0*dphi
           diffsmall=123456789.0d0
           index=1
c
c           write(*,*)'phihor = ',phor1
           if(phor1.gt.0.0d0)then
             do 10 i=1,Npad
               if(ypad(i).le.0.0d0)go to 10
               diff=dabs(phor1-xpad(i))
               if(diff.lt.diffsmall)then
                 index=i
                 diffsmall=diff
               endif
10           continue
c
c             write(*,101)index,ypad(index),phor1
101          format('&&&&&&& ',i4,3x,1pe17.7,2x,(0pf15.12,2x))
c
c   Make sure the difference is less than dphi
c
             if(diffsmall.gt.dphi)go to 99
c
c   if the difference is greater than 0.5dphi, then add a small 
c   correction
c
             if(diffsmall.lt.halfphi)then
               corr1=ypad(index)*(diffsmall-dphi)/dphi
c             else
               t1=diffsmall-halfphi
               t2=(t1/halfphi)**2
               corr1=t2*0.5d0*diffsmall*ypad(index)
             endif  
           endif
c
99         continue
c
c   Check the other limb
c
           index=1
           diffsmall=123456789.0d0
           if(phor2.gt.0.0d0)then
             do 20 i=1,Npad
               if(ypad(i).le.0.0d0)go to 20
               diff=dabs(phor2-xpad(i))
               if(diff.lt.diffsmall)then
                 index=i
                 diffsmall=diff
               endif
20           continue
c
c   Make sure the difference is less than dphi
c
             if(diffsmall.gt.dphi)go to 999
c
c   if the difference is greater than 0.5dphi, then add a small 
c   correction
c
             if(diffsmall.lt.halfphi)then
               corr2=ypad(index)*(diffsmall-dphi)/dphi
c             else
               t1=diffsmall-halfphi
               t2=(t1/halfphi)**2
               corr2=t2*0.5d0*diffsmall*ypad(index)
             endif  
           endif
c
999        corr=corr1+corr2
c
c           write(*,*)corr1,corr2
c           if(corr1.gt.0.0d0)write(*,*)'corr1'
c           if(corr2.gt.0.0d0)write(*,*)'corr2'
c
           corr=0.0d0

           return
           end
c
c
c   ##################################
c
           subroutine intep(xp,p,x,f,n,ier)
c
           implicit double precision(a-h,o-z)
           double precision lp1,lp2,l1,l2
           dimension f(N),x(N)
           ier=1
           i0=1
           iup=0
           if(x(2).lt.x(1))iup=1
           N1=N-1

           if((xp.ge.x(N).and.iup.eq.0).or.(xp.le.x(N).and.iup.eq.1))then
5            p=f(N)
             go to 6
           else if((xp.le.x(1).and.iup.eq.0).or.(xp.ge.x(1).and.iup.eq.1))then
             p=f(1)
6            ier=2
             return
           endif
c
c           entry eintep(xp,p,x,f,n,ier)
8          do 1 i=i0,N
             if(xp.lt.x(i).and.iup.eq.0)go to 2
             if(xp.gt.x(i).and.iup.eq.1)go to 2
1          continue
           go to 5
2          i=i-1
           if(i.eq.i0-1)go to 4
           i0=i+1
           lp1=1.0d0/(x(i)-x(i+1))
           lp2=1.0d0/(x(i+1)-x(i))
           if(i.eq.1)fp1=(f(2)-f(1))/(x(2)-x(1))
           if(i.eq.1)go to 3
           fp1=(f(i+1)-f(i-1))/(x(i+1)-x(i-1))
3          if(i.ge.n1)fp2=(f(n)-f(n-1))/(x(n)-x(n-1))         
           if(i.ge.n1)go to 4
           fp2=(f(i+2)-f(i))/(x(i+2)-x(i))
4          xpi1=xp-x(i+1)
           xpi=xp-x(i)
           l1=xpi1*lp1
           l2=xpi*lp2
           p=f(i)*(1.0d0-2.0d0*lp1*xpi)*l1*l1+f(i+1)*(1.0d0-2.0d0*lp2*xpi1)
     $         *l2*l2+fp2*xpi1*l2*l2+fp1*xpi*l1*l1
           if(p.lt.0.0d0)p=0.0d0
           return
           end
c
c  &&&&&&&&&&&&&&&&&&&
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
c
c
          subroutine getcontimes(finc,period,ecc,argper,T0,tconj1,tconj2)
c
c   November 5, 2008
c
c   This routine will compute the conjunction times given the input parameters.
c   It basically finds the value of the true anamoly that minimizes 
c   the projected
c   separation of centers.
c
          implicit double precision (a-h,o-z)
c
          pie=3.141592653589793d0
c
          fincr=finc*pie/180.0d0
          omegar=argper*pie/180.0d0
          ft=1.5d0*pie-omegar
3139      if(ft.lt.0.0d0)ft=ft+2.0d0*pie
          if(ft.lt.0.0d0)go to 3139
3140      if(ft.ge.2.0d0*pie)ft=ft-2.0d0*pie
          if(ft.ge.2.0d0*pie)go to 3140

          guess=ft
          do 22 jj=1,20
            guess=guess-deltap(fincr,omegar,ecc,guess)/
     &            deltapp(fincr,omegar,ecc,guess)
22        continue

          ft=guess
          htrc=0.5d0*ft
          if(dabs(0.5d0*pie-htrc).lt.7.0d-6)go to 31101
          if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 31101
          ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
          goto 31103
31101     ecan=pie
31103     xmc=ecan-ecc*dsin(ecan)
c          if(xmc.lt.0.0d0)xmc=xmc+2.0d0*pie
          deltaT=(xmc*period)/(2.0d0*pie)
          ttt=deltaT/period
          tconj1=t0+deltaT
c
          ft=0.5d0*pie-omegar
4139      if(ft.lt.0.0d0)ft=ft+2.0d0*pie
          if(ft.lt.0.0d0)go to 4139
4140      if(ft.ge.2.0d0*pie)ft=ft-2.0d0*pie
          if(ft.ge.2.0d0*pie)go to 4140

          guess=ft
          do 422 jj=1,20
            guess1=guess
            guess=guess-deltap(fincr,omegar,ecc,guess)/
     &            deltapp(fincr,omegar,ecc,guess)
422       continue

          ft=guess
          htrc=0.5d0*ft
          if(dabs(0.5d0*pie-htrc).lt.7.0d-6)go to 41101
          if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 41101
          ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
          goto 41103
41101     ecan=pie
41103     xmc=ecan-ecc*dsin(ecan)
c          if(xmc.lt.0.0d0)xmc=xmc+2.0d0*pie
          deltaT=(xmc*period)/(2.0d0*pie)
          ttt=deltaT/period
          tconj2=t0+deltaT
c
          return
          end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
          function delta(fincr,omegar,ecc,theta)
c
          implicit double precision (a-h,o-z)
c
          t1=1.0d0-ecc*ecc
          t2=1.0d0+ecc*dcos(theta)
          t3=1.0d0-(dsin(fincr)*dsin(theta+omegar))**2
          delta=t1/t2*dsqrt(t3)
c
          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          function deltap(fincr,omegar,ecc,theta)
c
          implicit double precision (a-h,o-z)
c
          tiny=1.0d-4
c
          t1=delta(fincr,omegar,ecc,theta+tiny)
          t2=delta(fincr,omegar,ecc,theta-tiny)
          deltap=(t1-t2)/(2.0d0*tiny)

          return
          end
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          function deltapp(fincr,omegar,ecc,theta)
c
          implicit double precision (a-h,o-z)
c
          tiny=1.0d-4
c
          t1=deltap(fincr,omegar,ecc,theta+tiny)
          t2=deltap(fincr,omegar,ecc,theta-tiny)
          deltapp=(t1-t2)/(2.0d0*tiny)
          return
          end
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
c
c
          subroutine getT0(finc,period,ecc,argper,T0,Tconj)
c
c   January 30, 2010.
c
c   This routine is the inverse of getcontimes.  Given a time of transit,
c   it will figure out the T0 needed.
c
          implicit double precision (a-h,o-z)
c
          pie=3.141592653589793d0
c
          fincr=finc*pie/180.0d0
          omegar=argper*pie/180.0d0
          ft=1.5d0*pie-omegar
3139      if(ft.lt.0.0d0)ft=ft+2.0d0*pie
          if(ft.lt.0.0d0)go to 3139
3140      if(ft.ge.2.0d0*pie)ft=ft-2.0d0*pie
          if(ft.ge.2.0d0*pie)go to 3140

          guess=ft
          do 22 jj=1,20
            guess=guess-deltap(fincr,omegar,ecc,guess)/
     &            deltapp(fincr,omegar,ecc,guess)
22        continue

          ft=guess
          htrc=0.5d0*ft
          if(dabs(0.5d0*pie-htrc).lt.7.0d-6)go to 31101
          if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 31101
          ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
          goto 31103
31101     ecan=pie
31103     xmc=ecan-ecc*dsin(ecan)
c          if(xmc.lt.0.0d0)xmc=xmc+2.0d0*pie
          deltaT=(xmc*period)/(2.0d0*pie)
          ttt=deltaT/period

c
          ft=0.5d0*pie-omegar
4139      if(ft.lt.0.0d0)ft=ft+2.0d0*pie
          if(ft.lt.0.0d0)go to 4139
4140      if(ft.ge.2.0d0*pie)ft=ft-2.0d0*pie
          if(ft.ge.2.0d0*pie)go to 4140

          guess=ft
          do 422 jj=1,20
            guess1=guess
            guess=guess-deltap(fincr,omegar,ecc,guess)/
     &            deltapp(fincr,omegar,ecc,guess)
422       continue

          ft=guess
          htrc=0.5d0*ft
          if(dabs(0.5d0*pie-htrc).lt.7.0d-6)go to 41101
          if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 41101
          ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
          goto 41103
41101     ecan=pie
41103     xmc=ecan-ecc*dsin(ecan)
c          if(xmc.lt.0.0d0)xmc=xmc+2.0d0*pie
          deltaT=(xmc*period)/(2.0d0*pie)
          ttt=deltaT/period

          t0=tconj-deltaT
c
          return
          end
c
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
           subroutine binlctime(Nmaxphase,Nphase,xmod,ymod,period,tstep,sw29)
c
c   This subroutine will bin a light curve in time, using a binsize given
c   by sw29, where the units are minutes.
c
           implicit double precision (a-h,o-z)
c
           dimension xmod(Nmaxphase),ymod(Nmaxphase)
           dimension yinter(900000),y2(900000)
           dimension xpad(900000),ypad(900000)
c
           NNN=1
           do kk=1,Nphase
              xpad(kk)=xmod(kk)
              ypad(kk)=ymod(kk)
           enddo
           Mphase=Nphase
c
c   We will interpolate the model, and given each value in xmod, figure
c   out the range of phase needed, and average the interpolated y-values
c
           call spline(xpad,ypad,Mphase,0.0d0,0.0d0,y2)
c
           pstep=sw29/1440.0d0/period
           pstep=sw29/60.0d0/24.0d0
           pwidth=0.5d0*pstep
c
           NN=7
           xxxhigh=1.0d0            !xmod(Nphase)
           do 10 i=1,Mphase
c             if(xpad(i).gt.1.0d0)go to 10
c             if(xpad(i).lt.0.0d0)go to 10
             kcount=0
             summ=0.0d0        
             jlo=i
c 
             a=xpad(i)-pwidth
             b=xpad(i)+pwidth

             do 9 j=1,NN
               if(j.eq.1)then
                 xxx=a
c                 if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c                 if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
                 call hunt(xpad,Mphase,xxx,jlo)
                 if((jlo.eq.Mphase).or.(jlo.eq.0))then
                   call splint(xpad,ypad,y2,Mphase,xxx,qqqa)
                 else
                   call fastsplint(xpad,ypad,y2,Mphase,xxx,qqqa,jlo,jlo+1)
                 endif
                 xxx=b
c                 if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c                 if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
                 call hunt(xpad,Mphase,xxx,jlo)
                 if((jlo.eq.Mphase).or.(jlo.eq.0))then
                   call splint(xpad,ypad,y2,Mphase,xxx,qqqb)
                 else
                   call fastsplint(xpad,ypad,y2,Mphase,xxx,qqqb,jlo,jlo+1)
                 endif
                 summ=0.5d0*(b-a)*(qqqa+qqqb)
               else
                 it=2**(j-2)
                 tnm=dble(it)
                 del=(b-a)/tnm
                 x=a+0.5d0*del
                 xxx=x
c                 if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c                 if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
                 s=0.0d0
                 do 11 kk=1,it
                   kcount=kcount+1 
                   call hunt(xpad,Mphase,xxx,jlo)
                   if((jlo.eq.Mphase).or.(jlo.eq.0))then
                     call splint(xpad,ypad,y2,Mphase,xxx,qqq)
                   else
                     call fastsplint(xpad,ypad,y2,Mphase,xxx,qqq,jlo,jlo+1)
                   endif
                   s=s+qqq
                   x=x+del
                   xxx=x
c                   if(xxx.lt.0.0d0)xxx=1.0d0+xxx   !dabs(xxx)
c                   if(xxx.gt.xxxhigh)xxx=xxx-1.0d0
11               continue
                 summ=0.5d0*(summ+(b-a)*s/tnm)
               endif
c
9           continue

           yinter(i)=summ/pstep

10         continue
c
           do 20 i=1,Mphase
             ymod(i)=yinter(i)
20         continue
c
c           call removepad(Nphase,xmod,ymod,xpad,ypad)

           return
           end
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  &&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&
c
      include 'limbdarksubs.for'
c
