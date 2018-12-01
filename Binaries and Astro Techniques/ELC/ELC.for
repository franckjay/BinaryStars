         program ELC
c
c   October 6, 1999
c
c   This is a total rewrite of the Avni code.  
c
          implicit double precision(a-h,o-z)
c
          parameter(Nmaxphase=720001)    ! was 1500000
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm to 9.
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm to 11
c
c
c   UPDATE October 10, 2008
c
c   change the dimensions of obsparm, obv, sobv, and eobv to 11
c
          dimension obsparm(17)
          dimension xmod(Nmaxphase),ymodU(Nmaxphase),ymodB(Nmaxphase),
     $      ymodV(Nmaxphase),ymodR(Nmaxphase),ymodI(Nmaxphase),
     $      ymodJ(Nmaxphase),ymodH(Nmaxphase),ymodK(Nmaxphase),
     &      ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      RV1(Nmaxphase),RV2(Nmaxphase),ymods3(Nmaxphase),
     $      drv1(Nmaxphase),drv2(Nmaxphase)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
          dimension xRVmod(Nmaxphase)
c
c   RVG BUG ALERT  June 12, 2001
c
c   Dimension the variables needed for the model atmospheres here.
c
          parameter (maxlines=1300,maxmu=115)  ! was 960
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
c
c   RVG BUG ALERT   May 8, 2001
c
c   Define the variables needed for spots
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
c
c
c   UPDATE JULY 7, 2004
c
c   This array is needed for the sub-random Sobel sequence, used in
c   place of ran9
c
          dimension xsob(2)
c
c   UPDATE May 8, 2006
c
c   Add this array for power-law limb darkening coefficients.  
c   8=filter index, 9=coefficient c1,c2,c3 ...
c
          dimension powercoeff(8,9)
c
c
c   UPDATE January 12, 2009
c
c   make fracs fracs1, fracs2, .... fracs8 and put them in the
c   argument of subroutine lightcurve.  In this was one can use
c   Nmaxphase in the dimension statement
c
          dimension fracs1(Nmaxphase,3),fracs2(Nmaxphase,3)
          dimension fracs3(Nmaxphase,3),fracs4(Nmaxphase,3)
          dimension fracs5(Nmaxphase,3),fracs6(Nmaxphase,3)
          dimension fracs7(Nmaxphase,3),fracs8(Nmaxphase,3)
c
c   NEW BUG ALERT  July 13, 2001
c
c   Add a new character string and common block for a 'parameter string'
c   This string of parameters will be fed to the genetic code to make it
c   easier to compute uncertainties on the physical quantities like mass
c   and radius.
c
c   UPDATE June 14, 2002
c
c   Make the length of parmstring character*237.
c
c   UPDATE October 31, 2002
c
c   Make the length of parmstring 249
c
c   UPDATE October 22, 2008
c
c   Make the length of parmstring 259
c
          character*259 parmstring
c
          common /stringblock/ parmstring
c
c   NEW BUG August 2, 2001
c
c   Change 'sw4' to 'T0'
c
c   UPDATE August 10, 2004
c
c   Add the 8 variables below.
c
c   UPDATE May 8, 2006
c
c   Add sw21-sw24, powercoeff below.
c  
c   UPDATE November 6, 2008
c
c   Add sw25-sw34
c
          common /realblock/ fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     &       alb1,alb2,
     %       rLx,Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,dwavey,
     %       t3,g3,SA3,density,sw1,sw2,sw3,T0,
     %       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     &         primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       bigI,bigbeta,sw23,sw24,powercoeff,
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2
c
c   RVG BUG ALERT   May 8, 2001
c
c   Add this common block to all programs
c
          common /spotblock/ spot1parm,spot2parm,spotdparm
c
c
c   UPDATE August 10, 2004
c
c   Add the 4 switches below.
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
c
c   RVG BUG ALERT  June 12, 2001
c
c   Add these common blocks for the model atmosphere variables
c
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,atmint4,
     %       atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin
          common /intatm/  Nlines,Nmu
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
c
        ifastflag=0   !disable fast genetic mode
c
c   RVG BUG ALERT  May 9, 2001
c
c   Add the spot parameters to the arguments of getinput and recordparm
c
c   UPDATE August 10, 2004
c
c   Add the 8 real variables and 4 integer variables to the list.
c
          call getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
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
     %       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
c
c   RVG BUG ALERT   May 16, 2001
c
c   Move the recordparm subroutine to inside the lcsubs.for file.
c
c          call recordparm(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
c     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
c     $       rinner,router,tdisk,xi,
c     &       Ntheta,Nradius,alb1,alb2,Nref,
c     %       rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
c     $       idraw,iecheck,idint,iatm,ism1,
c     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,
c     &       isw2,isw3,isw4,
c     &       ilaw,wave,dbolx,dboly,dwavex,dwavey,
c     $       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
c     $       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
c     %       spotdparm)
c
c   RVG BUG ALERT  June 12, 2001
c
c   Load the atmosphere table here, rather than in the subroutine
c
c  If the flag iatm>0, then load the model atmosphere table.
c
          if(iatm.ge.1)then
            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &         atmint1,atmint2,atmint3,atmint4,atmint5,
     @         atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin)
          endif
          write(*,*)Nlines
c
c
c   UPDATE JULY 7, 2004
c
c   Initialize the Sobel sequence here.
c
          nnn=-11
          call sobseq(nnn,xsob)
c
c   UPDATE January 30, 2010
c
c   if isw28 is greater than 0, then set the value of T0 based
c   on the time of transit Tconj.
c
          if(isw28.gt.0)then
            call getT0(finc,period,ecc,argper,T0,Tconj)
          endif
c
          call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $      ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %      ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c


          call wlinmod(Nphase,xmod,ymodU,'modelU.linear',isw7)
          call wlinmod(Nphase,xmod,ymodB,'modelB.linear',isw7)
          call wlinmod(Nphase,xmod,ymodV,'modelV.linear',isw7)
          call wlinmod(Nphase,xmod,ymodR,'modelR.linear',isw7)
          call wlinmod(Nphase,xmod,ymodI,'modelI.linear',isw7)
          call wlinmod(Nphase,xmod,ymodJ,'modelJ.linear',isw7)
          call wlinmod(Nphase,xmod,ymodH,'modelH.linear',isw7)
          call wlinmod(Nphase,xmod,ymodK,'modelK.linear',isw7)
c
          call wmagmod(Nphase,xmod,ymodU,'modelU.mag',isw7)
          call wmagmod(Nphase,xmod,ymodB,'modelB.mag',isw7)
          call wmagmod(Nphase,xmod,ymodV,'modelV.mag',isw7)
          call wmagmod(Nphase,xmod,ymodR,'modelR.mag',isw7)
          call wmagmod(Nphase,xmod,ymodI,'modelI.mag',isw7)
          call wmagmod(Nphase,xmod,ymodJ,'modelJ.mag',isw7)
          call wmagmod(Nphase,xmod,ymodH,'modelH.mag',isw7)
          call wmagmod(Nphase,xmod,ymodK,'modelK.mag',isw7)

          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(idint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',isw7)
c
          call wlinmod(NRVphase,xRVmod,RV1,'star1.RV',isw7)
          call wlinmod(NRVphase,xRVmod,RV2,'star2.RV',isw7)

          call wlinmod(NRVphase,xRVmod,dRV1,'star1.delRV',isw7)
          call wlinmod(NRVphase,xRVmod,dRV2,'star2.delRV',isw7)
c

c          
c   RVG BUG ALERT   May 16, 2001
c
c   Close the input file within the lightcurve subroutine, not here.
c
c          close(2)   ! close the output file
c
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'



