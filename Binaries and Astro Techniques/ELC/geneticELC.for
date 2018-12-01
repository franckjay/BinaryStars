          program geneticELC
c
c   May 11, 2001
c
c   This program will read in the light curves and parameters specified
c   in the 'gridloop.opt' file and optimize the fits using a genetic
c   algorithm (P. Charbonneau, 1995, ApJS, 101, 309).
c
c   
c   July 10, 2001
c
c   Program modified to include "black sheep" (A. Bobinger, 2000,
c   A&A 357, 1170).  After each population replacement and ranking,
c   replace the worst members with copies of the best member where
c   the genes are varied by small random amounts.  Specifically, the
c   number of black sheep = int(Npop*[1-pcross]/2.0), and 
c   parray(j,i)=parray(j,i)+0.05*gasdev(idum)
c
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
          implicit double precision (a-h,o-z)
c
c   UPDATE August 13, 2001
c
c   Change the number of maximum variables to 25 (was 19).  Change
c   the 19 to 25 in the dimension of parmarray,dummy,fakevar and svar
c
          parameter(Nmaxphase=50000,Ndatamax=200000)
          parameter(Nvmax=25)
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm,obv,eobv,sobv to 9.
c
c
c   UPDATE March 19, 2002
c
c   Add the parameter 'ifixgamma' to the end of the argument list
c   of checklcfit.  It is missing from some calls, but not others.
c
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv  to 11
c
c
c   UPDATE October 10, 2008
c
c   change the dimensions of obsparm, obv, sobv, and eobv to 17
c
          dimension obsparm(17),obv(17),eobv(17)
          dimension xmod(Nmaxphase),ymodU(Nmaxphase),ymodB(Nmaxphase),
     $      ymodV(Nmaxphase),ymodR(Nmaxphase),ymodI(Nmaxphase),
     $      ymodJ(Nmaxphase),ymodH(Nmaxphase),ymodK(Nmaxphase),
     &      ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      RV1(Nmaxphase),RV2(Nmaxphase),ymods3(Nmaxphase)
          dimension xdataU(Ndatamax),ydataU(Ndatamax),errU(Ndatamax),
     &      ydataB(Ndatamax),errB(Ndatamax),ydataV(Ndatamax),errV(Ndatamax),
     &      ydataR(Ndatamax),errR(Ndatamax),ydataI(Ndatamax),errI(Ndatamax),
     &      ydataJ(Ndatamax),errJ(Ndatamax),ydataH(Ndatamax),errH(Ndatamax),
     $      ydataK(Ndatamax),errK(Ndatamax),yRV1(Ndatamax),errRV1(Ndatamax),
     %      yRV2(Ndatamax),errRV2(Ndatamax),xdataB(Ndatamax),xdataV(Ndatamax),
     &      xdataR(Ndatamax),xdataI(Ndatamax),xdataJ(Ndatamax),
     &      xdataH(Ndatamax),xdataK(Ndatamax),xRV1(Ndatamax),xRV2(Ndatamax),
     $      drv1(Nmaxphase),drv2(Nmaxphase)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %      sigvar(Nvmax),stepsave(Nvmax)
          character*40 Udatafile,svar(nvmax),Hdatafile,Kdatafile,
     &             RV1file,RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*40 fakevar(Nvmax),sobv(17),command
          dimension parmarray(Nvmax,10000),chiarr(10000),indxchi(10000)
          dimension dummy(Nvmax)
          dimension xRVmod(Nmaxphase)
c
c   UPDATE August 13, 2001
c
c   Make the variable 'line' 200 characters
c
c   UPDATE October 28, 2002
c
c   Make the length of line 300 (was 200)
c
c
          character*300 line        ! was 132
          character*7 extension
c
c   RVG BUG ALERT  June 12, 2001
c
c   Dimension the variables needed for the model atmospheres here.
c
          parameter (maxlines=1300,maxmu=115)   !was 1100
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)

c
c   RVG BUG ALERT   May 8, 2001
c
c   Define variables for the spots
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
c
c   Here are some variables needed for the genetic part.
c
          parameter(nd=5,ipmax=1024,nmax=32,idmax=8,
     %       pmutmn=0.0005d0,pmutmx=0.25d0)

          dimension ign1(nmax*idmax),ign2(nmax*idmax),ifit(ipmax),jfit(ipmax)
          dimension fitns(ipmax),ph(nmax,2),oldph(nmax,ipmax),
     &        rnewph(nmax,ipmax)
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
c   NEW BUG ALERT  July 13, 2001
c
c   Add a new character string and common block for a 'parameter string'
c   This string of parameters will be fed to the genetic code to make it
c   easier to compute uncertainties on the physical quantities like mass
c   and radius.
c
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
c   Make the length of parmstring character*237.
c
c   UPDATE October 28, 2002
c
c   Make the length of parmstring character*249.
c
c   UPDATE October 22, 20028
c
c   Make the length of parmstring character*259.
c

          character*259 parmstring
c
          common /stringblock/ parmstring
c
c   NEW BUG August 2, 2001
c
c   New variables are needed to allow for the fitting of period and T0.
c
          dimension savxdataU(Ndatamax),savydataU(Ndatamax),saverrU(Ndatamax),
     &      savydataB(Ndatamax),saverrB(Ndatamax),savydataV(Ndatamax),
     #      saverrV(Ndatamax),
     &      savydataR(Ndatamax),saverrR(Ndatamax),savydataI(Ndatamax),
     $      saverrI(Ndatamax),
     &      savydataJ(Ndatamax),saverrJ(Ndatamax),savydataH(Ndatamax),
     #      saverrH(Ndatamax),
     $      savydataK(Ndatamax),saverrK(Ndatamax),savyRV1(Ndatamax),
     $      saverrRV1(Ndatamax),
     %      savyRV2(Ndatamax),saverrRV2(Ndatamax),savxdataB(Ndatamax),
     $      savxdataV(Ndatamax),
     &      savxdataR(Ndatamax),savxdataI(Ndatamax),savxdataJ(Ndatamax),
     &      savxdataH(Ndatamax),savxdataK(Ndatamax),savxRV1(Ndatamax),
     #      savxRV2(Ndatamax)
c
c
c   UPDATE OCTOBER 10, 2007
c
c   Add this array to keep track of the light curves of the individual
c   components in all band passes.
c
          dimension compfracs(8,2),ochidisk(8)
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
          common /fracblock/ compfracs          
c            
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list, and change 'sw4'
c   to 'T0'
c
c   UPDATE August 10, 2004
c
c   Add the 8 variables below.
c
c   UPDATE May 8, 2006
c
c   Add sw21-sw24 and powercoeff below.
c
c   UPDATE November 6, 2008
c
c   Add sw25-sw34 below
c
          common /realblock/ fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     &       alb1,alb2,
     %       rLx,Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,dwavey,
     %       t3,g3,SA3,density,sw1,sw2,sw3,T0,
     %       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     &       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
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
          common /ranblock/ idum
c
c   RVG BUG ALERT  June 12, 2001
c
c   Add these common blocks for the model atmosphere variables
c
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,
     &      atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin
          common /intatm/  Nlines,Nmu
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
c
c   RVG BUG ALERT  May 9, 2001
c
c   Add the spot parameters to getinput and recordparm.
c
c   UPDATE August 10, 2004
c
c   Add the 8 realand 4 integer variables to the list.
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
     %       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,
     &       temprat,idark1,idark2,isw12,isw13,
     %       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff, 
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2,
     %       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
c
c   UPDATE October 28, 2002
c
c   Define a new variable called sepsave, which is equal to the
c   value of separ given in ELC.inp
c
           savesep=separ

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
c
c   RVG BUG ALERT  June 12, 2001
c
c   Load the atmosphere table here, rather than in the subroutine
c
c  If the flag iatm>0, then load the model atmosphere table.
c
          gmin=0.0d0
          gmax=0.0d0
          Tmin=0.0d0
          Tmax=0.0d0
          if(iatm.ge.1)then
            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &         atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,
     &         atmint8,Tmax,Tmin,gmax,gmin)
          endif

c
c   NEW BUG August 2, 2001
c
c   Use these flags when fitting for the period and/or T0
c
          itime=isw7
c
c   UPDATE JULY 30, 2004
c
c   set Nbin = 0
c
          Nbin=0
c
c   UPDATE May 27, 2002
c
c   Add a new variable to enable median fitting.
c   if sw6 > 1.0, then use median fitting
c
          rmed=sw6 
c
c   Disable the draw option
c
          idraw=0
c
c   Do we force the gamma to be the input value?
c
          ifixgamma=isw4
c
c   UPDATE April 15, 2002
c
c   Define the values of gamma1 and gamma2
c
          gamma1=gamma
          gamma2=gamma
c
c   Fix the inner disk radius at fill2.
c
          if(teff2.gt.0.0)rinner=fill2
c
c   UPDATE August 13, 2001
c
c   If the flag isw9 (=ielite) is more than 0, insert the parameters
c   set in the ELC.inp file into the population.
c
c
c   UPDATE JULY 7, 2004
c
c   Initialize the Sobel sequence here.
c
          mmm=-11
          call sobseq(mmm,xsob)
c
c   Initialize the steps.
c
          do 1 i=1,nvmax
            Nstep(i)=1
            var(i)=0.0d0
            svar(i)='none'
 1        continue
c
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
          call initNdata(NdataU,NdataB,NdataV,NdataR,
     &        NdataI,NdataJ,NdataH,NdataK,
     %        NRV1,NRV2,Nobv)
c
c   Get the parameters for the grid search.
c
          call getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,Nvmax,Nvar,
     $      svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 5, 2005
c
c   Add also the ability to include the luminosity ratio in the total chi^2
c
          ifrac=-99
          ilum=-99
          do 111 i=1,Nobv
            icn=icnvrt(sobv(i)(1:2))
            if((icn.eq.104).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.112).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.113).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.114).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.115).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.116).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.117).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.118).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.119).and.(idint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.369))ilum=99 ! luminosity ratio = yes
 111      continue
c
c   Figure out which data files were specified.
c
          icnU=icnvrt(Udatafile(1:2))
          icnB=icnvrt(Bdatafile(1:2))
          icnV=icnvrt(Vdatafile(1:2))
          icnR=icnvrt(Rdatafile(1:2))
          icnI=icnvrt(Idatafile(1:2))
          icnJ=icnvrt(Jdatafile(1:2))
          icnH=icnvrt(Hdatafile(1:2))
          icnK=icnvrt(Kdatafile(1:2))
          icnRV1=icnvrt(RV1file(1:2))
          icnRV2=icnvrt(RV2file(1:2))
c
c   icn? = 430 indicates that no file was specified.  Load the files listed.
c
          if(icnU.ne.430)call loaddata(Ndatamax,Udatafile,
     &        NdataU,xdataU,ydataU,errU)
          if(icnB.ne.430)call loaddata(Ndatamax,Bdatafile,
     &        NdataB,xdataB,ydataB,errB)
          if(icnV.ne.430)call loaddata(Ndatamax,Vdatafile,
     &        NdataV,xdataV,ydataV,errV)
          if(icnR.ne.430)call loaddata(Ndatamax,Rdatafile,
     &        NdataR,xdataR,ydataR,errR)
          if(icnI.ne.430)call loaddata(Ndatamax,Idatafile,
     &        NdataI,xdataI,ydataI,errI)
          if(icnJ.ne.430)call loaddata(Ndatamax,Jdatafile,
     &        NdataJ,xdataJ,ydataJ,errJ)
          if(icnH.ne.430)call loaddata(Ndatamax,Hdatafile,
     &        NdataH,xdataH,ydataH,errH)
          if(icnK.ne.430)call loaddata(Ndatamax,Kdatafile,
     &        NdataK,xdataK,ydataK,errK)
          if(icnRV1.ne.430)call loaddata(Ndatamax,RV1file,
     &        NRV1,xRV1,yRV1,errRV1)
          if(icnRV2.ne.430)call loaddata(Ndatamax,RV2file,
     &        NRV2,xRV2,yRV2,errRV2)
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or T0, we have to make copies
c   of the arrays.
c
          if(icnU.ne.430)call kopydata(Ndatamax,NdataU,
     $       xdataU,ydataU,errU,isavNU,
     $       savxdataU,savydataU,saverrU)
          if(icnB.ne.430)call kopydata(Ndatamax,NdataB,
     $       xdataB,ydataB,errB,isavNB,
     $       savxdataB,savydataB,saverrB)
          if(icnV.ne.430)call kopydata(Ndatamax,NdataV,
     $       xdataV,ydataV,errV,isavNV,
     $       savxdataV,savydataV,saverrV)
          if(icnI.ne.430)call kopydata(Ndatamax,NdataI,
     $       xdataI,ydataI,errI,isavNI,
     $       savxdataI,savydataI,saverrI)
          if(icnJ.ne.430)call kopydata(Ndatamax,NdataJ,
     $       xdataJ,ydataJ,errJ,isavNJ,
     $       savxdataJ,savydataJ,saverrJ)
          if(icnH.ne.430)call kopydata(Ndatamax,NdataH,
     $       xdataH,ydataH,errH,isavNH,
     $       savxdataH,savydataH,saverrH)
          if(icnK.ne.430)call kopydata(Ndatamax,NdataK,
     $       xdataK,ydataK,errK,isavNK,
     $       savxdataK,savydataK,saverrK)
          if(icnR.ne.430)call kopydata(Ndatamax,NdataR,
     $       xdataR,ydataR,errR,isavNR,
     $       savxdataR,savydataR,saverrR)
          if(icnRV1.ne.430)call kopydata(Ndatamax,NRV1,
     $       xRV1,yRV1,errRV1,isavRV1,
     $       savxRV1,savyRV1,saverrRV1)
          if(icnRV2.ne.430)call kopydata(Ndatamax,NRV2,
     $       xRV2,yRV2,errRV2,isavRV2,
     $       savxRV2,savyRV2,saverrRV2)
c
c   Sort the data files by phase just to be on the safe side.
c
c   UPDATE April 15, 2002
c
c   Put an if-then clause to cover the case when Ndata=1
c
          if((icnU.ne.430).and.(NdataU.gt.1))
     #         call sort3(NdataU,xdataU,ydataU,errU)
          if((icnB.ne.430).and.(NdataB.gt.1))
     #         call sort3(NdataB,xdataB,ydataB,errB)
          if((icnV.ne.430).and.(NdataV.gt.1))
     #         call sort3(NdataV,xdataV,ydataV,errV)
          if((icnR.ne.430).and.(NdataR.gt.1))
     #         call sort3(NdataR,xdataR,ydataR,errR)
          if((icnI.ne.430).and.(NdataI.gt.1))
     #         call sort3(NdataI,xdataI,ydataI,errI)
          if((icnJ.ne.430).and.(NdataJ.gt.1))
     #         call sort3(NdataJ,xdataJ,ydataJ,errJ)
          if((icnH.ne.430).and.(NdataH.gt.1))
     #         call sort3(NdataH,xdataH,ydataH,errH)
          if((icnK.ne.430).and.(NdataK.gt.1))
     #         call sort3(NdataK,xdataK,ydataK,errK)
          if((icnRV1.ne.430).and.(NRV1.gt.1))
     #         call sort3(NRV1,xRV1,yRV1,errRV1)
          if((icnRV2.ne.430).and.(NRV2.gt.1))
     #         call sort3(NRV2,xRV2,yRV2,errRV2)
c
          isvel1=0
          isvel2=0
          if(icnRV1.ne.430)isvel1=299
          if(icnRV2.ne.430)isvel2=299
c
          idum=-1234567
c
c   UPDATE November 6, 2008
c 
c   if sw25 > 0, use that as an error for asini.
c
          saveasini=sw5
c
c   Start the looping here.  
c
          chi1=0.0d0
          chi2=0.0d0
          chi3=0.0d0
          small=9.0d30
          Ndattot=0
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+NdataK+
     %        NRV1+NRV2
c
c   UPDATE February 5, 2003
c  
c   Check the value of ielete.  If it is negative, then we need to
c   go into recover mode.
c

c
c   Open an output file for the fitting statistics.
c
          open(unit=38,file='geneticELC.out',status='unknown')
          open(unit=45,file='generation.1000000',status='unknown')
c
c   October 10, 2007
c
c   Add this additional output file if the flag is set.
c                      
          if(isw24.ge.1)open(unit=47,file='ELCratio.1000000',status='unknown')

c
c   UPDATE October 12, 2001
c
c   Open another output file for chi^2 values
c
          open(unit=55,file='chi.1000000',status='unknown')
c
c   NEW BUG ALERT July 13, 2001
c
c   Open a new output file for parameters.
c
          open(unit=46,file='ELCparm.1000000',status='unknown')

          Nterms=0
          do 699 i7=1,Nvar
            var(i7)=vstart(i7)     !initialize the variables
            kkk=icnvrt(svar(i7)(1:2))
            if(kkk.eq.430)then
              go to 749
            else
              Nterms=Nterms+1
            endif
 699      continue
c
          
 749      continue
c
          write(*,*)'Nterms, Nvar, Ndattot',Nterms,Nvar,Ndattot

c
c   Here are some variables needed for the genetic part.
c
          pcross=0.85d0 
          pmut=0.005d0
          fdif=1.0d0
c   
c
c   Define the random parameter sets.  Nstep(1) is the number of random
c   sets to define.
c

          np=Nstep(1)     ! number of members in the population 
          ngen=Nstep(2)   ! number of generations

          do 99 j=1,Nterms
            vx1=vstart(j)
            vx2=vstep(j)
            vxmult=(vx2-vx1)
            vxadd=vx1
            do 98 i=1,np                 
              urand=ran1(idum)
              varx=urand*vxmult+vxadd
              parmarray(j,i)=varx
              oldph(j,i)=urand
 98         continue
 99       continue

c
c   Evaluate the fitness of the initial population.
c
          do 750 i16=1,np
c
c   UPDATE JULY 21, 2006
c
c   Add a fast mode where the light curve is first computed using a very
c   coarse grid and a large dphase.  Then compare the chi^2 to the current
c   lowest.  If the chi^2 is very large, then use the chi^2 from the coarse
c   grid.  If the chi^2 is close to the smallest value, repeat the light
c   curve computation with the normal grid.
c
c   We need to always compute the first model, so set ifastflag=0 when i16=1
c
            if(isw22.ge.1)ifastflag=1
            if(i16.eq.1)ifastflag=0
            do 760 j=1,Nterms
              var(j)=parmarray(j,i16)
 760        continue
c
c   UPDATE August 16, 2001
c
c   If the isw9 flag is 1 or greater, then alter the first population member
c   and put the parameters inside ELC.inp into the var array
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
            if((isw9.eq.1).and.(i16.eq.1))then
c
c   UPDATE November 6, 2002
c
c   Add the limb darkening coefficients dwavex and dwavey
c   to the argument list of varassign.
c
              call varassign(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     %          pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %          ecosw,temprat,bigI,bigbeta,powercoeff,density,
     %          Tconj,beam1,beam2,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
              do 799 j=1,Nterms
                vx1=vstart(j)
                vx2=vstep(j)
                vxmult=(vx2-vx1)
                vxadd=vx1
                varx=var(j)
                urand=(varx-vxadd)/vxmult
                if(urand.gt.1.0d0)urand=1.0d0
                if(urand.lt.0.0d0)urand=0.0d0
                parmarray(j,1)=varx
                oldph(j,1)=urand
 799         continue
c
            endif
c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
c   UPDATE October 28, 2002
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              if(savesep.lt.0.0d0)separ=savesep 
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
              call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $         rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     %         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &         dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %         ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
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
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
c  February 12, 2000
c
c  Save the filling factors in case there is an eccentric orbit
c
              remf1=fill1
              remf2=fill2
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or T0 (itime>0), then we
c   have to fold the saved arrays.
c
              if(itime.gt.0)then
                if(icnU.ne.430)call phasefold(Ndatamax,isavNU,
     $            savxdataU,savydataU,saverrU,
     #            Nbin,NdataU,xdataU,ydataU,errU,period,T0,ikeep,ecc,argper,isw7)
                if(icnB.ne.430)call phasefold(Ndatamax,isavNB,
     $            savxdataB,savydataB,saverrB,
     #            Nbin,NdataB,xdataB,ydataB,errB,period,T0,ikeep,ecc,argper,isw7)
                if(icnV.ne.430)call phasefold(Ndatamax,isavNV,
     $            savxdataV,savydataV,saverrV,
     #            Nbin,NdataV,xdataV,ydataV,errV,period,T0,ikeep,ecc,argper,isw7)
                if(icnR.ne.430)call phasefold(Ndatamax,isavNR,
     $            savxdataR,savydataR,saverrR,
     #            Nbin,NdataR,xdataR,ydataR,errR,period,T0,ikeep,ecc,argper,isw7)
                if(icnI.ne.430)call phasefold(Ndatamax,isavNI,
     $            savxdataI,savydataI,saverrI,
     #            Nbin,NdataI,xdataI,ydataI,errI,period,T0,ikeep,ecc,argper,isw7)
                if(icnJ.ne.430)call phasefold(Ndatamax,isavNJ,
     $            savxdataJ,savydataJ,saverrJ,
     #            Nbin,NdataJ,xdataJ,ydataJ,errJ,period,T0,ikeep,ecc,argper,isw7)
                if(icnH.ne.430)call phasefold(Ndatamax,isavNH,
     $            savxdataH,savydataH,saverrH,
     #            Nbin,NdataH,xdataH,ydataH,errH,period,T0,ikeep,ecc,argper,isw7)
                if(icnK.ne.430)call phasefold(Ndatamax,isavNK,
     $            savxdataK,savydataK,saverrK,
     #            Nbin,NdataK,xdataK,ydataK,errK,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV1.ne.430)call phasefold(Ndatamax,isavRV1,
     $            savxRV1,savyRV1,saverrRV1,
     #            Nbin,NRV1,xRV1,yRV1,errRV1,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV2.ne.430)call phasefold(Ndatamax,isavRV2,
     $            savxRV2,savyRV2,saverrRV2,
     #            Nbin,NRV2,xRV2,yRV2,errRV2,period,T0,ikeep,ecc,argper,isw7)
              endif
c
c  Get the models.
c
c  UPDATE November 6, 2008
c
c  If sw25>0, use that as an error for asini (sw5)
c
               if(sw25.gt.0.0d0)sw5=saveasini+gasdev(idum)*sw25
c
1122           call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $         ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
c
              fill1=remf1
              fill2=remf2
c
c
c   Get the chi^2 values
c
              call initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &            chisqH,chisqK,chisqRV1,chisqRV2)
c
              if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %         ydataU,errU,chisqU,zeroU,ifixgamma,isw7)
              if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %         ydataB,errB,chisqB,zeroB,ifixgamma,isw7)
              if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %         ydataV,errV,chisqV,zeroV,ifixgamma,isw7)
              if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %         ydataR,errR,chisqR,zeroR,ifixgamma,isw7)
              if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %         ydataI,errI,chisqI,zeroI,ifixgamma,isw7)
              if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %         ydataJ,errJ,chisqJ,zeroJ,ifixgamma,isw7)
              if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %         ydataH,errH,chisqH,zeroH,ifixgamma,isw7)
              if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %         ydataK,errK,chisqK,zeroK,ifixgamma,isw7)

              if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
                call checkRVfit(islc,NRVphase,xRVmod,
     %            RV1,NRVphase,xRVmod,RV2,
     #            NRV1,xRV1,yRV1,errRV1,
     %            NRV2,xRV2,yRV2,errRV2,
     &            chisqRV1,chisqRV2,gam,ifixgamma)
                gamma1=gam
                gamma2=gam
              else
                if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,RV1,
     &            NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,ggamma,ifixgamma,isw7)
                 gamma1=ggamma
                 if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,
     %             RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,ggamma,ifixgamma,isw7)
                 gamma2=gamma
              endif

c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
              ochi=0.0d0
              if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c    UPDATE FEBRUARY 5, 2005
c  
c    Add also the luminosity ratio.
c
              if(ifrac.gt.10)call diskchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)

c
              if(ilum.gt.10)call lrchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
              chi1=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
c
              if(ifastflag.ge.1)then
                if(chi1.lt.3.0d0*small*dabs(dble(Ndattot-Nterms)))then
                  ifastflag=0
                  go to 1122
                endif
              endif
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
              call writevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
 5555         format(a)        ! UPDATE August 13, 2001: was format(a132)

c
c   UPDATE May 27, 2002
c
c   If rmed > 1, then we are doing median fitting.  In that case,
c   change the label of the 'chi1' variable.
c
              if(rmed.ge.1.0d0)then
                call printmed(chi1)
              else
                call printchi(chi1)
              endif
c
              chiarr(i16)=chi1
              fitns(i16)=1.0d0/chi1    ! the fitness is  prop.to. chi^2
c
c              write(*,*)chiarr(i16),fitns(i16)
c
              call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
              if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 5, 2005
c
c   Add the luminosity ratio also
c
              if(ifrac.gt.10)call wdiskobschi(ochidisk,Nobv,sobv)
              if(ilum.gt.10)call wlrobschi(ochilr)
c
c   UPDATE October 12, 2001
c
c   Write the chi^2 values to the new output file (unit 55).
c
             call chiline(i16,chi1,ochidisk,
     #        icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %        icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &        chisqJ,chisqH,chisqK,
     &        chisqRV1,chisqRV2,Nobv,sobv,obv,eobv,obsparm,ochilr)
c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of newwritevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients dwavex and dwavey to the
c   argument list of newwritevar.
c
              call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,chi1,i16,line,alb1,alb2,
     %          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
c   NEW BUG ALERT July 13, 2001
c
c   Write to the new output file.
c
c   UPDATE October 28, 2002
c
c   Add the string 0.0d0 to the end of 9999999.
c
c
              cccc=chi1
              if(cccc.gt.9999999.0d0)cccc=9999999.0d0
              write(46,4646)i16,cccc,parmstring
c
c   UPDATE January 16, 2002
c
c   parmstring was character*201, now should be character*227. Thus
c   change a199 to a227
c
c   UPDATE June 7, 2002
c
c   Make the length of parmstring character*237, and hence the format
c   field to a237
c
c   UPDATE October 22, 2008
c
c   Make the length of parmstring character*259, and hence the format
c   field to a237.  Also, add digits to chi^2
c
 4646         format(i4,1x,f16.7,1x,a259)
c
              chi1=chi1/dabs(dble(Ndattot-Nterms))
c
c   UPDATE August 13, 2001
c
c   Truncate the string line to the last blank character.
c           
              lll=lnblnk(line)
              write(45,5555)line(1:lll)
c
c   October 10, 2007
c
c   write the ratios if the flag is set
c
              if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ijk,1),ijk=1,8),
     #           (compfracs(ijk,2),ijk=1,8)
c
5566          format(i5,1x,f14.4,1x,16(1pe13.6,1x))
c
c   UPDATE June 17, 2002
c
c   define nnn=0.
c
              nnn=0
              write(*,444)i16,nnn
c
c              write(*,*)chisqU,chisqB,chisqV
c              write(*,*)chisqR,chisqI,chisqJ
c              write(*,*)chisqH,chisqK,chisqRV1
c              write(*,*)chi1,Nterms,Ndattot
c
c              chi11=dsqrt(chi1)*dsqrt(chi1)
c              write(*,*)chi1,chi11

6543          format(2(f19.4,2x),2x,I7,3x,i7)
              if(chi1.lt.small)then
                small=chi1
                s1=var(1)
                s2=var(2)
                s3=var(3)
                s4=var(4)
                s5=var(5)
                s6=var(6)
                s7=var(7)
                s8=var(8)
                s9=var(9)
                s10=var(10)
                s11=var(11)
                s12=var(12)
                s13=var(13)
                s14=var(14)
                s15=var(15)
                s16=var(16)
c
c  UPDATE August 13, 2001
c
c  Pad out the 's' variables to 25
c
                s17=var(17)
                s18=var(18)
                s19=var(19)
                s20=var(20)
                s21=var(21)
                s22=var(22)
                s23=var(23)
                s24=var(24)
                s25=var(25)
              endif
 750        continue
            close(45)
c
c   UPDATE October 12, 2001
c
c   Close the new output file (unit 55).
c
            close(55)
c 
c   Sort the chiarr and print out parameters for the best Nstep(2) sets.
c
          call rnkpop(np,fitns,ifit,jfit)
          call indexx(Np,chiarr,indxchi)
          write(38,8850)nnn,chiarr(indxchi(1)),chiarr(indxchi(np/2)),pmut
cc
c
c   reset the variables at their best values and print the chi^2
c
          var(1)=s1
          var(2)=s2
          var(3)=s3
          var(4)=s4
          var(5)=s5
          var(6)=s6
          var(7)=s7
          var(8)=s8
          var(9)=s9
          var(10)=s10
          var(11)=s11
          var(12)=s12
          var(13)=s13
          var(14)=s14
          var(15)=s15
          var(16)=s16
c
c  UPDATE August 13, 2001
c
c  Pad out the 's' variables to 25
c
          var(17)=s17
          var(18)=s18
          var(19)=s19
          var(20)=s20
          var(21)=s21
          var(22)=s22
          var(23)=s23
          var(24)=s24
          var(25)=s25

c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
c   UPDATE October 28, 2002
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
          if(savesep.lt.0.0d0)separ=savesep
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
          call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $      omega2,Q,finc,Teff1,Teff2,betarim,
     $      rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
     &       pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &       dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
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
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
c  February 12, 2000
c
c  Save the filling factors in case there is an eccentric orbit
c
          remf1=fill1
          remf2=fill2
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or T0 (itime>0), then we
c   have to fold the saved arrays.
c
              if(itime.gt.0)then
                if(icnU.ne.430)call phasefold(Ndatamax,isavNU,
     $            savxdataU,savydataU,saverrU,
     #            Nbin,NdataU,xdataU,ydataU,errU,period,T0,ikeep,ecc,argper,isw7)
                if(icnB.ne.430)call phasefold(Ndatamax,isavNB,
     $            savxdataB,savydataB,saverrB,
     #            Nbin,NdataB,xdataB,ydataB,errB,period,T0,ikeep,ecc,argper,isw7)
                if(icnV.ne.430)call phasefold(Ndatamax,isavNV,
     $            savxdataV,savydataV,saverrV,
     #            Nbin,NdataV,xdataV,ydataV,errV,period,T0,ikeep,ecc,argper,isw7)
                if(icnR.ne.430)call phasefold(Ndatamax,isavNR,
     $            savxdataR,savydataR,saverrR,
     #            Nbin,NdataR,xdataR,ydataR,errR,period,T0,ikeep,ecc,argper,isw7)
                if(icnI.ne.430)call phasefold(Ndatamax,isavNI,
     $            savxdataI,savydataI,saverrI,
     #            Nbin,NdataI,xdataI,ydataI,errI,period,T0,ikeep,ecc,argper,isw7)
                if(icnJ.ne.430)call phasefold(Ndatamax,isavNJ,
     $            savxdataJ,savydataJ,saverrJ,
     #            Nbin,NdataJ,xdataJ,ydataJ,errJ,period,T0,ikeep,ecc,argper,isw7)
                if(icnH.ne.430)call phasefold(Ndatamax,isavNH,
     $            savxdataH,savydataH,saverrH,
     #            Nbin,NdataH,xdataH,ydataH,errH,period,T0,ikeep,ecc,argper,isw7)
                if(icnK.ne.430)call phasefold(Ndatamax,isavNK,
     $            savxdataK,savydataK,saverrK,
     #            Nbin,NdataK,xdataK,ydataK,errK,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV1.ne.430)call phasefold(Ndatamax,isavRV1,
     $            savxRV1,savyRV1,saverrRV1,
     #            Nbin,NRV1,xRV1,yRV1,errRV1,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV2.ne.430)call phasefold(Ndatamax,isavRV2,
     $            savxRV2,savyRV2,saverrRV2,
     #            Nbin,NRV2,xRV2,yRV2,errRV2,period,T0,ikeep,ecc,argper,isw7)
              endif
c
cc   UPDATE JULY 21, 2006
c
c   Add a fast mode where the light curve is first computed using a very
c   coarse grid and a large dphase.  Then compare the chi^2 to the current
c   lowest.  If the chi^2 is very large, then use the chi^2 from the coarse
c   grid.  If the chi^2 is close to the smallest value, repeat the light
c   curve computation with the normal grid.
c
c   Here we want the coarse grid.
c
            ifastflag=0

c  Get the models.
c
c
c  UPDATE November 6, 2008
c
c  If sw25>0, use that as an error for asini (sw5)
c
          if(sw25.gt.0.0d0)sw5=saveasini+gasdev(idum)*sw25
c
          call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $         ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
          fill1=remf1
          fill2=remf2
c
c   Get the chi^2 values
c
c   NOTE March 19, 2002
c
c   The first argument to checklcfit is a 9, indicating that a light
c   curve is being fit.  The last argument is a zero, and is a place
c   holder for 'ifixgamma'.
c
          if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %      ydataU,errU,chisqU,zeroU,0,isw7)
          if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %      ydataB,errB,chisqB,zeroB,0,isw7)
          if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %      ydataV,errV,chisqV,zeroV,0,isw7)
          if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %      ydataR,errR,chisqR,zeroR,0,isw7)
          if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %      ydataI,errI,chisqI,zeroI,0,isw7)
          if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %      ydataJ,errJ,chisqJ,zeroJ,0,isw7)
          if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %      ydataH,errH,chisqH,zeroH,0,isw7)
          if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %      ydataK,errK,chisqK,zeroK,0,isw7)
          if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
            call checkRVfit(islc,NRVphase,xRVmod,RV1,NRVphase,xRVmod,RV2,
     #             NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     &             chisqRV1,chisqRV2,gam,ifixgamma)
            gamma1=gam
            gamma2=gam
          else
            if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,RV1,NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,ggamma,ifixgamma,isw7)
            gamma1=ggamma
            if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,ggamma,ifixgamma,isw7)
            gamma2=ggamma
          endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
          call writevar(Nvmax,svar,var,fill1,fill2,omega1,
     $        omega2,Q,finc,Teff1,Teff2,betarim,
     $        rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &        t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     #        period,T0,alb1,alb2,
     #        dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
          ochi=0.0d0
          if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 5, 2005
c
c   Add the luminosity ratio.
c
              if(ifrac.gt.10)call diskchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)
c
              if(ilum.gt.10)call lrchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
          chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
          call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
c
          if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 10, 2005
c
c   Add in the luminosity ratio also.
c
              if(ifrac.gt.10)call wdiskobschi(ochidisk,Nobv,sobv)
              if(ilum.gt.10)call wlrobschi(ochilr)
c
          write(*,*)' '
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
          call writevar(Nvmax,fakevar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercroeff to list above.
c
          write(*,*)' '
          write(*,1002)chiall
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
          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(idint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',isw7)
c
          call wlinmod(NRVphase,xRVmod,RV1,'star1.RV',isw7)
          call wlinmod(NRVphase,xRVmod,RV2,'star2.RV',isw7)
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or the T0, write the current
c   folded light curves
c
          if(itime.gt.0)then
            if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,ydataU,
     %              errU,'ELCdataU.fold')
            if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,ydataB,
     %              errB,'ELCdataB.fold')
            if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,ydataV,
     %              errV,'ELCdataV.fold')
            if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,ydataR,
     %              errR,'ELCdataR.fold')
            if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,ydataI,
     %              errI,'ELCdataI.fold')
            if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,ydataJ,
     %              errJ,'ELCdataJ.fold')
            if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,ydataH,
     %              errH,'ELCdataH.fold')
            if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,ydataK,
     %              errK,'ELCdataK.fold')
            if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,yRV1,
     %              errRV1,'ELCdataRV1.fold')
            if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,yRV2,
     %              errRV2,'ELCdataRV2.fold')
          endif
c
c   Main generation  loop
c
c   UPDATE February 5, 2003
c
c   Add the ability to recover from a crash.  At the end of each generation,
c   write the value of np, ig, iseed.  Then write the values of
c   each gene to a file.
c
c


          do 10000 ig=1,ngen
c
c   Main population loop
c
            write(extension,3333)ig+1000000
 3333       format(i7)
            open(unit=45,file='generation.'//extension,status='unknown')
c
            if(isw24.ge.1)open(unit=47,file='ELCratio.'//extension,status='unknown')
c
c   UPDATE October 12, 2001
c
c   Open another output file for chi^2 values
c
          open(unit=55,file='chi.'//extension,status='unknown')
c
c   NEW BUG ALERT July 13, 2001
c
c   Open a new output file for parameters.
c
          open(unit=46,file='ELCparm.'//extension,status='unknown')
c
            do 20000 ip=1,np/2
c
              n=Nterms
c
c   1.  Pick two parents.
c
              call select(np,jfit,fdif,ip1)
21000         call select(np,jfit,fdif,ip2)
              if(ip1.eq.ip2)go to 21000
c
c   2.  Encode parent phenotypes.
c
              do 19000 kkk=1,Nterms
                dummy(kkk)=oldph(kkk,ip1)
19000         continue
              call encode(n,nd,dummy,ign1)

              do 19001 kkk=1,Nterms
                dummy(kkk)=oldph(kkk,ip2)
19001         continue
              call encode(n,nd,dummy,ign2)
c
c   3.  Breed.
c
              call cross(n,nd,pcross,ign1,ign2)
              call mutate(n,nd,pmut,ign1)
              call mutate(n,nd,pmut,ign2)
c
c   4.  Decode offspring genotypes.
c
              call decode(n,nd,ign1,dummy)
              do 19002 kkk=1,Nterms
                ph(kkk,1)=dummy(kkk)
19002        continue

              call decode(n,nd,ign2,dummy)
              do 19003 kkk=1,Nterms
                ph(kkk,2)=dummy(kkk)
19003        continue

c             
c   5.  Insert into population.
c
              call genrep(nmax,n,np,ip,ph,rnewph)
c
20000       continue
c
c   We now must replace the old population by the new, and evaluate
c   and rank the fitness.  First, put the fittest member of the old
c   population into the new one.
c
             do 44000 k=1,Nterms
               rnewph(k,1)=oldph(k,ifit(np))
44000        continue
c
             do 22000 i=1,np
               do 33000 k=1,Nterms
                 oldph(k,i)=rnewph(k,i)
33000          continue
22000        continue
c
            do 999 j=1,Nterms
              vx1=vstart(j)
              vx2=vstep(j)
              vxmult=(vx2-vx1)
              vxadd=vx1
              do 998 i=1,np                 
                urand=oldph(j,i)
                varx=urand*vxmult+vxadd
                parmarray(j,i)=varx
 998          continue
 999        continue
c
c   Get the new fitness.
c
          do 7750 i16=1,np
c
c
c   UPDATE JULY 21, 2006
c
c   Add a fast mode where the light curve is first computed using a very
c   coarse grid and a large dphase.  Then compare the chi^2 to the current
c   lowest.  If the chi^2 is very large, then use the chi^2 from the coarse
c   grid.  If the chi^2 is close to the smallest value, repeat the light
c   curve computation with the normal grid.
c
c   We need to always compute the first model, so set ifastflag=0 when i16=1
c
            if(isw22.ge.1)ifastflag=1
            if(i16.eq.1)ifastflag=0
c
            do 7760 j=1,Nterms
              var(j)=parmarray(j,i16)
 7760       continue

c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
c   UPDATE October 28, 2002
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              if(savesep.lt.0.0d0)separ=savesep
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
              call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $         omega2,Q,finc,Teff1,Teff2,betarim,
     $         rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     %         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &         dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %         ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
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
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
c  February 12, 2000
c
c  Save the filling factors in case there is an eccentric orbit
c
              remf1=fill1
              remf2=fill2
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or T0 (itime>0), then we
c   have to fold the saved arrays.
c
              if(itime.gt.0)then
                if(icnU.ne.430)call phasefold(Ndatamax,isavNU,
     $            savxdataU,savydataU,saverrU,
     #            Nbin,NdataU,xdataU,ydataU,errU,period,T0,ikeep,ecc,argper,isw7)
                if(icnB.ne.430)call phasefold(Ndatamax,isavNB,
     $            savxdataB,savydataB,saverrB,
     #            Nbin,NdataB,xdataB,ydataB,errB,period,T0,ikeep,ecc,argper,isw7)
                if(icnV.ne.430)call phasefold(Ndatamax,isavNV,
     $            savxdataV,savydataV,saverrV,
     #            Nbin,NdataV,xdataV,ydataV,errV,period,T0,ikeep,ecc,argper,isw7)
                if(icnR.ne.430)call phasefold(Ndatamax,isavNR,
     $            savxdataR,savydataR,saverrR,
     #            Nbin,NdataR,xdataR,ydataR,errR,period,T0,ikeep,ecc,argper,isw7)
                if(icnI.ne.430)call phasefold(Ndatamax,isavNI,
     $            savxdataI,savydataI,saverrI,
     #            Nbin,NdataI,xdataI,ydataI,errI,period,T0,ikeep,ecc,argper,isw7)
                if(icnJ.ne.430)call phasefold(Ndatamax,isavNJ,
     $            savxdataJ,savydataJ,saverrJ,
     #            Nbin,NdataJ,xdataJ,ydataJ,errJ,period,T0,ikeep,ecc,argper,isw7)
                if(icnH.ne.430)call phasefold(Ndatamax,isavNH,
     $            savxdataH,savydataH,saverrH,
     #            Nbin,NdataH,xdataH,ydataH,errH,period,T0,ikeep,ecc,argper,isw7)
                if(icnK.ne.430)call phasefold(Ndatamax,isavNK,
     $            savxdataK,savydataK,saverrK,
     #            Nbin,NdataK,xdataK,ydataK,errK,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV1.ne.430)call phasefold(Ndatamax,isavRV1,
     $            savxRV1,savyRV1,saverrRV1,
     #            Nbin,NRV1,xRV1,yRV1,errRV1,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV2.ne.430)call phasefold(Ndatamax,isavRV2,
     $            savxRV2,savyRV2,saverrRV2,
     #            Nbin,NRV2,xRV2,yRV2,errRV2,period,T0,ikeep,ecc,argper,isw7)
              endif
c
c  Get the models.
c
c
c  UPDATE November 6, 2008
c
c  If sw25>0, use that as an error for asini (sw5)
c
               if(sw25.gt.0.0d0)sw5=saveasini+gasdev(idum)*sw25
c

1133           call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $         ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
              fill1=remf1
              fill2=remf2
c
c   Get the chi^2 values
c
              call initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &            chisqH,chisqK,chisqRV1,chisqRV2)
c
c   NOTE March 19, 2002
c
c   Add 'ifixgamma' to the end of the argument list for checklcfit.  For
c   light curves (first argument '9'), this parameter is ignored.
c
              if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %         ydataU,errU,chisqU,zeroU,ifixgamma,isw7)
              if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %         ydataB,errB,chisqB,zeroB,ifixgamma,isw7)
              if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %         ydataV,errV,chisqV,zeroV,ifixgamma,isw7)
              if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %         ydataR,errR,chisqR,zeroR,ifixgamma,isw7)
              if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %         ydataI,errI,chisqI,zeroI,ifixgamma,isw7)
              if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %         ydataJ,errJ,chisqJ,zeroJ,ifixgamma,isw7)
              if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %         ydataH,errH,chisqH,zeroH,ifixgamma,isw7)
              if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %         ydataK,errK,chisqK,zeroK,ifixgamma,isw7)
              if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
                call checkRVfit(islc,NRVphase,xRVmod,
     %            RV1,NRVphase,xRVmod,RV2,
     #            NRV1,xRV1,yRV1,errRV1,
     %            NRV2,xRV2,yRV2,errRV2,
     &            chisqRV1,chisqRV2,gam,ifixgamma)
                gamma1=gam
                gamma2=gam
              else
                if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,
     &            RV1,NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,ggamma,ifixgamma,isw7)
                 gamma1=ggamma
                 if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,
     &             RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,ggamma,ifixgamma,isw7)
                 gamma2=ggamma
              endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
              ochi=0.0d0
              if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio.
c   
              if(ifrac.gt.10)call diskchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)
c
              if(ilum.gt.10)call lrchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
              chi1=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
c
              if(ifastflag.ge.1)then
                if(chi1.lt.3.0d0*small*dabs(dble(Ndattot-Nterms)))then
                  ifastflag=0
                  go to 1133
                endif
              endif
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
              call writevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
c
c   UPDATE May 27, 2002
c
c   If rmed > 1, then we are doing median fitting.  In that case,
c   change the label of the 'chi1' variable.
c
              if(rmed.ge.1.0d0)then
                call printmed(chi1)
              else
                call printchi(chi1)
              endif
c
              chiarr(i16)=chi1
              fitns(i16)=1.0d0/chi1    ! the fitness is prop.to. chi^2
c
              call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
              if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio also.
c
              if(ifrac.gt.10)call wdiskobschi(ochidisk,Nobv,sobv)
              if(ilum.gt.10)call wlrobschi(ochilr)
c
c   UPDATE October 12, 2001
c
c   Write the chi^2 values to the new output file (unit 55).
c
              call chiline(i16,chi1,ochidisk,
     #          icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %          icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &        chisqJ,chisqH,chisqK,
     &          chisqRV1,chisqRV2,Nobv,sobv,obv,eobv,obsparm,ochilr)
c
c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of newwritevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients dwavex and dwavey to the
c   argument list of newwritevar.
c
              call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,chi1,i16,line,alb1,alb2,
     %          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
c
c   NEW BUG ALERT July 13, 2001
c
c   Write to the new output file.
c
c   UPDATE October 28, 2002
c
c   Add 0.0d0 to the string 9999999 below
c
              cccc=chi1
              if(cccc.gt.9999999.0d0)cccc=9999999.0d0
              write(46,4646)i16,cccc,parmstring
c
              chi1=chi1/dabs(dble(Ndattot-Nterms))
c
c
c   UPDATE August 13, 2001
c
c   Truncate the string line to the last blank character.
c           
              lll=lnblnk(line)
              write(45,5555)line(1:lll)
c
c   October 10, 2007
c
c   Write the fractions if requested
c
              if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ijk,1),ijk=1,8),
     #           (compfracs(ijk,2),ijk=1,8)
c
              write(*,444)i16,ig
c              if((i16.eq.1).and.(ig.eq.2))stop
 444          format('iteration number = ',i5,'  generation  number = ',i7)
              if(chi1.lt.small)then
                small=chi1
                s1=var(1)
                s2=var(2)
                s3=var(3)
                s4=var(4)
                s5=var(5)
                s6=var(6)
                s7=var(7)
                s8=var(8)
                s9=var(9)
                s10=var(10)
                s11=var(11)
                s12=var(12)
                s13=var(13)
                s14=var(14)
                s15=var(15)
                s16=var(16)
c
c  UPDATE August 13, 2001
c
c  Pad out the 's' variables to 25
c
                s17=var(17)
                s18=var(18)
                s19=var(19)
                s20=var(20)
                s21=var(21)
                s22=var(22)
                s23=var(23)
                s24=var(24)
                s25=var(25)
              endif
 7750       continue
c 
c   Sort the chiarr and print out parameters for the best Nstep(2) sets.
c
          call rnkpop(np,fitns,ifit,jfit)
          call indexx(Np,chiarr,indxchi)
c
c   NEW BUG ALERT July 10, 2001
c
c   Add the "black sheep stuff here.  Take the worst individuals and
c   replace them with mutated copies of the best guy
c
          Nsheep=dint(0.6d0*dble(np)*(1.0d0-pcross))+2
c
c   Increase the number of black sheep once the number of generations 
c   is larger than 75% of the population number.
c
          lll=dint(0.75d0*dble(np))
          if(ig.gt.lll)Nsheep=dint(0.7d0*dble(np)*(1.0d0-pcross))+2
          if(ig.gt.2*lll)Nsheep=dint(0.8d0*dble(np)*(1.0d0-pcross))+2
          if(ig.gt.3*lll)Nsheep=dint(1.0d0*dble(np)*(1.0d0-pcross))+2
c
          i16=i16-1               ! adjust the counter
c
          do 9950 isheep=np,np-Nsheep,-1
            do 7999 j=1,Nterms
              vx1=vstart(j)
              vx2=vstep(j)
              vxmult=(vx2-vx1)
              vxadd=vx1
              urand=oldph(j,indxchi(1))+0.05d0*gasdev(idum)
              if(urand.gt.1.0d0)urand=1.0d0
              if(urand.lt.0.0d0)urand=0.0d0
              varx=urand*vxmult+vxadd
              parmarray(j,indxchi(isheep))=varx
              var(j)=varx
              oldph(j,indxchi(isheep))=urand
 7999        continue
c
c   We now have created a mutated copy of the fittest individual for the
c   member at the bottom.  Evaluate the fitness, etc.
c
              if(isw22.ge.1)ifastflag=1
c
              i16=i16+1
              iii=indxchi(isheep)
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
c   UPDATE October 28, 2002
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              if(savesep.lt.0.0d0)separ=savesep
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
              call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $         rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     %         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &         dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %         ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
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
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
c  February 12, 2000
c
c  Save the filling factors in case there is an eccentric orbit
c
              remf1=fill1
              remf2=fill2
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or T0 (itime>0), then we
c   have to fold the saved arrays.
c
              if(itime.gt.0)then
                if(icnU.ne.430)call phasefold(Ndatamax,isavNU,
     $            savxdataU,savydataU,saverrU,
     #            Nbin,NdataU,xdataU,ydataU,errU,period,T0,ikeep,ecc,argper,isw7)
                if(icnB.ne.430)call phasefold(Ndatamax,isavNB,
     $            savxdataB,savydataB,saverrB,
     #            Nbin,NdataB,xdataB,ydataB,errB,period,T0,ikeep,ecc,argper,isw7)
                if(icnV.ne.430)call phasefold(Ndatamax,isavNV,
     $            savxdataV,savydataV,saverrV,
     #            Nbin,NdataV,xdataV,ydataV,errV,period,T0,ikeep,ecc,argper,isw7)
                if(icnR.ne.430)call phasefold(Ndatamax,isavNR,
     $            savxdataR,savydataR,saverrR,
     #            Nbin,NdataR,xdataR,ydataR,errR,period,T0,ikeep,ecc,argper,isw7)
                if(icnI.ne.430)call phasefold(Ndatamax,isavNI,
     $            savxdataI,savydataI,saverrI,
     #            Nbin,NdataI,xdataI,ydataI,errI,period,T0,ikeep,ecc,argper,isw7)
                if(icnJ.ne.430)call phasefold(Ndatamax,isavNJ,
     $            savxdataJ,savydataJ,saverrJ,
     #            Nbin,NdataJ,xdataJ,ydataJ,errJ,period,T0,ikeep,ecc,argper,isw7)
                if(icnH.ne.430)call phasefold(Ndatamax,isavNH,
     $            savxdataH,savydataH,saverrH,
     #            Nbin,NdataH,xdataH,ydataH,errH,period,T0,ikeep,ecc,argper,isw7)
                if(icnK.ne.430)call phasefold(Ndatamax,isavNK,
     $            savxdataK,savydataK,saverrK,
     #            Nbin,NdataK,xdataK,ydataK,errK,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV1.ne.430)call phasefold(Ndatamax,isavRV1,
     $            savxRV1,savyRV1,saverrRV1,
     #            Nbin,NRV1,xRV1,yRV1,errRV1,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV2.ne.430)call phasefold(Ndatamax,isavRV2,
     $            savxRV2,savyRV2,saverrRV2,
     #            Nbin,NRV2,xRV2,yRV2,errRV2,period,T0,ikeep,ecc,argper,isw7)
              endif
c
c  Get the models.
c
c
c  UPDATE November 6, 2008
c
c  If sw25>0, use that as an error for asini (sw5)
c
               if(sw25.gt.0.0d0)sw5=saveasini+gasdev(idum)*sw25
c
1144           call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $         ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
c
              fill1=remf1
              fill2=remf2
c
c
c   Get the chi^2 values
c
              call initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &            chisqH,chisqK,chisqRV1,chisqRV2)
c
c
c   NOTE March 19, 2002
c
c   Add 'ifixgamma' to the end of the argument list for checklcfit.  For
c   light curves (first argument '9'), this parameter is ignored.
c
              if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %         ydataU,errU,chisqU,zeroU,ifixgamma,isw7)
              if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %         ydataB,errB,chisqB,zeroB,ifixgamma,isw7)
              if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %         ydataV,errV,chisqV,zeroV,ifixgamma,isw7)
              if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %         ydataR,errR,chisqR,zeroR,ifixgamma,isw7)
              if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %         ydataI,errI,chisqI,zeroI,ifixgamma,isw7)
              if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %         ydataJ,errJ,chisqJ,zeroJ,ifixgamma,isw7)
              if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %         ydataH,errH,chisqH,zeroH,ifixgamma,isw7)
              if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %         ydataK,errK,chisqK,zeroK,ifixgamma,isw7)
              if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
                call checkRVfit(islc,NRVphase,xRVmod,
     %            RV1,NRVphase,xRVmod,RV2,
     #            NRV1,xRV1,yRV1,errRV1,
     %            NRV2,xRV2,yRV2,errRV2,
     &            chisqRV1,chisqRV2,gam,ifixgamma)
                gamma1=gam
                gamma2=gam
              else
                if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,
     &            RV1,NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,ggamma,ifixgamma,isw7)
                 gamma1=ggamma
                 if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,
     &             RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,ggamma,ifixgamma,isw7)
                 gamma2=ggamma
              endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
              ochi=0.0d0
              if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio.
c
              if(ifrac.gt.10)call diskchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)
c
              if(ilum.gt.10)call lrchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
              chi1=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
c
              if(ifastflag.ge.1)then
                if(chi1.lt.3.0d0*small*dabs(dble(Ndattot-Nterms)))then
                  ifastflag=0
                  go to 1144
                endif
              endif
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
              call writevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c   UPDATE May 27, 2002
c
c   If rmed > 1, then we are doing median fitting.  In that case,
c   change the label of the 'chi1' variable.
c
              if(rmed.ge.1.0d0)then
                call printmed(chi1)
              else
                call printchi(chi1)
              endif
c
              chiarr(iii)=chi1
              fitns(iii)=1.0d0/chi1    ! the fitness is prop.to. chi^2
c
              call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
              if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio.
c
              if(ifrac.gt.10)call wdiskobschi(ochidisk,Nobv,sobv)
              if(ilum.gt.10)call wlrobschi(ochilr)
c
c   UPDATE October 12, 2001
c
c   Write the chi^2 values to the new output file (unit 55).
c
             call chiline(i16,chi1,ochidisk,
     #        icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %        icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &        chisqJ,chisqH,chisqK,
     &        chisqRV1,chisqRV2,Nobv,sobv,obv,eobv,obsparm,ochilr)
c
c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of newwritevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients dwavex and dwavey to the
c   argument list of newwritevar.
c
              call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,chi1,i16,line,alb1,alb2,
     %          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
c
c
c   NEW BUG ALERT July 13, 2001
c
c   Write to the new output file.
c
c   UPDATE Oct 28, 2002
c
c   Add 0.0d0 to the 9999999.9 string below
c
              cccc=chi1
              if(cccc.gt.9999999.0d0)cccc=9999999.0d0
              write(46,4646)i16,cccc,parmstring
c
              chi1=chi1/dabs(dble(Ndattot-Nterms))
c
c   UPDATE August 13, 2001
c
c   Truncate the string line to the last blank character.
c           
              lll=lnblnk(line)
              write(45,5555)line(1:lll)
c
c   October 10, 2007
c
c   Write the fractions if requested.
c
              if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ijk,1),ijk=1,8),
     #           (compfracs(ijk,2),ijk=1,8)

              write(*,444)i16,ig

              if(chi1.lt.small)then
                small=chi1
                s1=var(1)
                s2=var(2)
                s3=var(3)
                s4=var(4)
                s5=var(5)
                s6=var(6)
                s7=var(7)
                s8=var(8)
                s9=var(9)
                s10=var(10)
                s11=var(11)
                s12=var(12)
                s13=var(13)
                s14=var(14)
                s15=var(15)
                s16=var(16)
c
c  UPDATE August 13, 2001
c
c  Pad out the 's' variables to 25
c
                s17=var(17)
                s18=var(18)
                s19=var(19)
                s20=var(20)
                s21=var(21)
                s22=var(22)
                s23=var(23)
                s24=var(24)
                s25=var(25)
              endif
c
 9950       continue     ! do 9950 isheep=np,np-Nsheep,-1
c
c   Now rerank the population, and continue with the normal genetic 
c   procedure.
c
c   END NEW BUG
c
          call rnkpop(np,fitns,ifit,jfit)
          call indexx(Np,chiarr,indxchi)
c
c   UPDATE November 1, 2002
c
c   Add a new loop to fine tune the variables.  If we have had more than
c   20 generations then
c   step through them ONE
c   AT A TIME and tweak by 0.005*gasdev(idum).
c
          if(ig.lt.20)go to 9952 

          Nsheep=Nterms
          ikount=0
c
          do 9951 isheep=np,np-Nsheep,-1
            ikount=ikount+1
            if(isw22.ge.1)ifastflag=1
            do 7998 j=1,Nterms
              vx1=vstart(j)
              vx2=vstep(j)
              vxmult=(vx2-vx1)
              vxadd=vx1
              urand=oldph(j,indxchi(1))
              if(j.eq.ikount)urand=oldph(j,indxchi(1))+0.005d0*gasdev(idum)
              if(urand.gt.1.0d0)urand=1.0d0
              if(urand.lt.0.0d0)urand=0.0d0
              varx=urand*vxmult+vxadd
              parmarray(j,indxchi(isheep))=varx
              var(j)=varx
              oldph(j,indxchi(isheep))=urand
 7998       continue
c
c   We now have created a mutated copy of the fittest individual for the
c   member at the bottom.  Evaluate the fitness, etc.
c
            i16=i16+1
            iii=indxchi(isheep)
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
c   UPDATE October 28, 2002
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
            if(savesep.lt.0.0d0)separ=savesep
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
            call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $         rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     %         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &         dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %         ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
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
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
c  February 12, 2000
c
c  Save the filling factors in case there is an eccentric orbit
c
            remf1=fill1
            remf2=fill2
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or T0 (itime>0), then we
c   have to fold the saved arrays.
c
            if(itime.gt.0)then
              if(icnU.ne.430)call phasefold(Ndatamax,isavNU,
     $            savxdataU,savydataU,saverrU,
     #            Nbin,NdataU,xdataU,ydataU,errU,period,T0,ikeep,ecc,argper,isw7)
              if(icnB.ne.430)call phasefold(Ndatamax,isavNB,
     $            savxdataB,savydataB,saverrB,
     #            Nbin,NdataB,xdataB,ydataB,errB,period,T0,ikeep,ecc,argper,isw7)
              if(icnV.ne.430)call phasefold(Ndatamax,isavNV,
     $            savxdataV,savydataV,saverrV,
     #            Nbin,NdataV,xdataV,ydataV,errV,period,T0,ikeep,ecc,argper,isw7)
              if(icnR.ne.430)call phasefold(Ndatamax,isavNR,
     $            savxdataR,savydataR,saverrR,
     #            Nbin,NdataR,xdataR,ydataR,errR,period,T0,ikeep,ecc,argper,isw7)
              if(icnI.ne.430)call phasefold(Ndatamax,isavNI,
     $            savxdataI,savydataI,saverrI,
     #            Nbin,NdataI,xdataI,ydataI,errI,period,T0,ikeep,ecc,argper,isw7)
              if(icnJ.ne.430)call phasefold(Ndatamax,isavNJ,
     $            savxdataJ,savydataJ,saverrJ,
     #            Nbin,NdataJ,xdataJ,ydataJ,errJ,period,T0,ikeep,ecc,argper,isw7)
              if(icnH.ne.430)call phasefold(Ndatamax,isavNH,
     $            savxdataH,savydataH,saverrH,
     #            Nbin,NdataH,xdataH,ydataH,errH,period,T0,ikeep,ecc,argper,isw7)
              if(icnK.ne.430)call phasefold(Ndatamax,isavNK,
     $            savxdataK,savydataK,saverrK,
     #            Nbin,NdataK,xdataK,ydataK,errK,period,T0,ikeep,ecc,argper,isw7)
              if(icnRV1.ne.430)call phasefold(Ndatamax,isavRV1,
     $            savxRV1,savyRV1,saverrRV1,
     #            Nbin,NRV1,xRV1,yRV1,errRV1,period,T0,ikeep,ecc,argper,isw7)
              if(icnRV2.ne.430)call phasefold(Ndatamax,isavRV2,
     $            savxRV2,savyRV2,saverrRV2,
     #            Nbin,NRV2,xRV2,yRV2,errRV2,period,T0,ikeep,ecc,argper,isw7)
            endif
c
c  Get the models.
c
c
c  UPDATE November 6, 2008
c
c  If sw25>0, use that as an error for asini (sw5)
c
               if(sw25.gt.0.0d0)sw5=saveasini+gasdev(idum)*sw25
c
1155           call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $         ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
c
            fill1=remf1
            fill2=remf2
c
c
c   Get the chi^2 values
c
            call initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &            chisqH,chisqK,chisqRV1,chisqRV2)
c
c
c   NOTE March 19, 2002
c
c   Add 'ifixgamma' to the end of the argument list for checklcfit.  For
c   light curves (first argument '9'), this parameter is ignored.
c
            if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %         ydataU,errU,chisqU,zeroU,ifixgamma,isw7)
            if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %         ydataB,errB,chisqB,zeroB,ifixgamma,isw7)
            if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %         ydataV,errV,chisqV,zeroV,ifixgamma,isw7)
            if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %         ydataR,errR,chisqR,zeroR,ifixgamma,isw7)
            if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %         ydataI,errI,chisqI,zeroI,ifixgamma,isw7)
            if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %         ydataJ,errJ,chisqJ,zeroJ,ifixgamma,isw7)
            if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %         ydataH,errH,chisqH,zeroH,ifixgamma,isw7)
            if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %         ydataK,errK,chisqK,zeroK,ifixgamma,isw7)
            if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
              call checkRVfit(islc,NRVphase,xRVmod,
     %            RV1,NRVphase,xRVmod,RV2,
     #            NRV1,xRV1,yRV1,errRV1,
     %            NRV2,xRV2,yRV2,errRV2,
     &            chisqRV1,chisqRV2,gam,ifixgamma)
              gamma1=gam
              gamma2=gam
            else
              if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,RV1,NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,ggamma,ifixgamma,isw7)
              gamma1=ggamma
              if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,ggamma,ifixgamma,isw7)
              gamma2=ggamma
            endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
            ochi=0.0d0
            if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio also.
c
            if(ifrac.gt.10)call diskchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)
c
            if(ilum.gt.10)call lrchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
            chi1=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
c
              if(ifastflag.ge.1)then
                if(chi1.lt.3.0d0*small*dabs(dble(Ndattot-Nterms)))then
                  ifastflag=0
                  go to 1155
                endif
              endif
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
            call writevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
c   UPDATE May 27, 2002
c
c   If rmed > 1, then we are doing median fitting.  In that case,
c   change the label of the 'chi1' variable.
c
            if(rmed.ge.1.0d0)then
              call printmed(chi1)
            else
              call printchi(chi1)
            endif
c
            chiarr(iii)=chi1
            fitns(iii)=1.0d0/chi1    ! the fitness is prop.to. chi^2
c
            call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
            if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio.
c
            if(ifrac.gt.10)call wdiskobschi(ochidisk,Nobv,sobv)
            if(ilum.gt.10)call wlrobschi(ochilr)
c
c   UPDATE October 12, 2001
c
c   Write the chi^2 values to the new output file (unit 55).
c
           call chiline(i16,chi1,ochidisk,
     #        icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %        icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &        chisqJ,chisqH,chisqK,
     &        chisqRV1,chisqRV2,Nobv,sobv,obv,eobv,obsparm,ochilr)
c
c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of newwritevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients dwavex and dwavey to the
c   argument list of newwritevar.
c
           call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,chi1,i16,line,alb1,alb2,
     %          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
c
c
c   NEW BUG ALERT July 13, 2001
c
c   Write to the new output file.
c
c   UPDATE Oct 28, 2002
c
c   Add 0.0d0 to the 9999999.9 string below
c
            cccc=chi1
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(46,4646)i16,cccc,parmstring
c
            chi1=chi1/dabs(dble(Ndattot-Nterms))
c
c   UPDATE August 13, 2001
c
c   Truncate the string line to the last blank character.
c           
            lll=lnblnk(line)
            write(45,5555)line(1:lll)
c
c   October 10, 2007
c
c   Write the fractions if requested
c
            if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ijk,1),ijk=1,8),
     #           (compfracs(ijk,2),ijk=1,8)

            write(*,444)i16,ig

            if(chi1.lt.small)then
              small=chi1
              s1=var(1)
              s2=var(2)
              s3=var(3)
              s4=var(4)
              s5=var(5)
              s6=var(6)
              s7=var(7)
              s8=var(8)
              s9=var(9)
              s10=var(10)
              s11=var(11)
              s12=var(12)
              s13=var(13)
              s14=var(14)
              s15=var(15)
              s16=var(16)
c
c  UPDATE August 13, 2001
c
c  Pad out the 's' variables to 25
c
              s17=var(17)
              s18=var(18)
              s19=var(19)
              s20=var(20)
              s21=var(21)
              s22=var(22)
              s23=var(23)
              s24=var(24)
              s25=var(25)
            endif
c
 9951     continue             ! do 9950 isheep=np,np-Nterms,-1
c
c   Now rerank the population, and continue with the normal genetic 
c   procedure.
c
c   END NEW BUG
c
          call rnkpop(np,fitns,ifit,jfit)
          call indexx(Np,chiarr,indxchi)
c
c   UPDATE November 1, 2002
c 
c   Add the statement label 9952 below.
c
 9952     write(38,8850)ig,chiarr(indxchi(1)),chiarr(indxchi(np/2)),pmut
c
 8850     format(i4,2x,2(f13.7,2x),2x,f6.4)
          call adjmut(np,fitns,ifit,pmutmn,pmutmx,pmut)
c
c   reset the variables at their best values and print the chi^2
c
          var(1)=s1
          var(2)=s2
          var(3)=s3
          var(4)=s4
          var(5)=s5
          var(6)=s6
          var(7)=s7
          var(8)=s8
          var(9)=s9
          var(10)=s10
          var(11)=s11
          var(12)=s12
          var(13)=s13
          var(14)=s14
          var(15)=s15
          var(16)=s16
c
c  UPDATE August 13, 2001
c
c  Pad out the 's' variables to 25
c
          var(17)=s17
          var(18)=s18
          var(19)=s19
          var(20)=s20
          var(21)=s21
          var(22)=s22
          var(23)=s23
          var(24)=s24
          var(25)=s25
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
c   UPDATE October 28, 2002
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
          if(savesep.lt.0.0d0)separ=savesep
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
          call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $      omega2,Q,finc,Teff1,Teff2,betarim,
     $      rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
     &       pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &       dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
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
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
c  February 12, 2000
c
c  Save the filling factors in case there is an eccentric orbit
c
          remf1=fill1
          remf2=fill2
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or T0 (itime>0), then we
c   have to fold the saved arrays.
c
              if(itime.gt.0)then
                if(icnU.ne.430)call phasefold(Ndatamax,isavNU,
     $            savxdataU,savydataU,saverrU,
     #            Nbin,NdataU,xdataU,ydataU,errU,period,T0,ikeep,ecc,argper,isw7)
                if(icnB.ne.430)call phasefold(Ndatamax,isavNB,
     $            savxdataB,savydataB,saverrB,
     #            Nbin,NdataB,xdataB,ydataB,errB,period,T0,ikeep,ecc,argper,isw7)
                if(icnV.ne.430)call phasefold(Ndatamax,isavNV,
     $            savxdataV,savydataV,saverrV,
     #            Nbin,NdataV,xdataV,ydataV,errV,period,T0,ikeep,ecc,argper,isw7)
                if(icnR.ne.430)call phasefold(Ndatamax,isavNR,
     $            savxdataR,savydataR,saverrR,
     #            Nbin,NdataR,xdataR,ydataR,errR,period,T0,ikeep,ecc,argper,isw7)
                if(icnI.ne.430)call phasefold(Ndatamax,isavNI,
     $            savxdataI,savydataI,saverrI,
     #            Nbin,NdataI,xdataI,ydataI,errI,period,T0,ikeep,ecc,argper,isw7)
                if(icnJ.ne.430)call phasefold(Ndatamax,isavNJ,
     $            savxdataJ,savydataJ,saverrJ,
     #            Nbin,NdataJ,xdataJ,ydataJ,errJ,period,T0,ikeep,ecc,argper,isw7)
                if(icnH.ne.430)call phasefold(Ndatamax,isavNH,
     $            savxdataH,savydataH,saverrH,
     #            Nbin,NdataH,xdataH,ydataH,errH,period,T0,ikeep,ecc,argper,isw7)
                if(icnK.ne.430)call phasefold(Ndatamax,isavNK,
     $            savxdataK,savydataK,saverrK,
     #            Nbin,NdataK,xdataK,ydataK,errK,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV1.ne.430)call phasefold(Ndatamax,isavRV1,
     $            savxRV1,savyRV1,saverrRV1,
     #            Nbin,NRV1,xRV1,yRV1,errRV1,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV2.ne.430)call phasefold(Ndatamax,isavRV2,
     $            savxRV2,savyRV2,saverrRV2,
     #            Nbin,NRV2,xRV2,yRV2,errRV2,period,T0,ikeep,ecc,argper,isw7)
              endif
c
c
c  Get the models.
c
          ifastflag=0  !need regular grid here
c
c
c  UPDATE November 6, 2008
c
c  If sw25>0, use that as an error for asini (sw5)
c
          if(sw25.gt.0.0d0)sw5=saveasini+gasdev(idum)*sw25
c
          call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $         ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
          fill1=remf1
          fill2=remf2
c
c   NOTE March 19, 2002
c
c   These calls to checklcfit have a zero at the end, which is a place
c   holder for
c   'ifixgamma'. For
c   light curves (first argument '9'), this parameter is ignored.
c
c
c   Get the chi^2 values
c
          if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %      ydataU,errU,chisqU,zeroU,0,isw7)
          if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %      ydataB,errB,chisqB,zeroB,0,isw7)
          if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %      ydataV,errV,chisqV,zeroV,0,isw7)
          if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %      ydataR,errR,chisqR,zeroR,0,isw7)
          if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %      ydataI,errI,chisqI,zeroI,0,isw7)
          if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %      ydataJ,errJ,chisqJ,zeroJ,0,isw7)
          if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %      ydataH,errH,chisqH,zeroH,0,isw7)
          if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %      ydataK,errK,chisqK,zeroK,0,isw7)
          if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
            call checkRVfit(islc,NRVphase,xRVmod,RV1,NRVphase,xRVmod,RV2,
     #             NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     &             chisqRV1,chisqRV2,gam,ifixgamma)
            gamma1=gam
            gamma2=gam
          else
            if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,RV1,NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,ggamma,ifixgamma,isw7)
            gamma1=ggamma
            if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,ggamma,ifixgamma,isw7)
            gamma2=ggamma
          endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
          call writevar(Nvmax,svar,var,fill1,fill2,omega1,
     $        omega2,Q,finc,Teff1,Teff2,betarim,
     $        rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &        t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $        period,T0,alb1,alb2,
     #        dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
          ochi=0.0d0
          if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio.
c
              if(ifrac.gt.10)call diskchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)
c
              if(ilum.gt.10)call lrchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
          chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
          call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
c
          if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio.
c
              if(ifrac.gt.10)call wdiskobschi(ochidisk,Nobv,sobv)
              if(ilum.gt.10)call wlrobschi(ochilr)
c
          write(*,*)' '
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
          call writevar(Nvmax,fakevar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
          write(*,*)' '
          write(*,1002)chiall
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
          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(idint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',isw7)
c
          call wlinmod(NRVphase,xRVmod,RV1,'star1.RV',isw7)
          call wlinmod(NRVphase,xRVmod,RV2,'star2.RV',isw7)
c
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or the T0, write the current
c   folded light curves
c
          if(itime.gt.0)then
            if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,ydataU,
     %              errU,'ELCdataU.fold')
            if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,ydataB,
     %              errB,'ELCdataB.fold')
            if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,ydataV,
     %              errV,'ELCdataV.fold')
            if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,ydataR,
     %              errR,'ELCdataR.fold')
            if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,ydataI,
     %              errI,'ELCdataI.fold')
            if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,ydataJ,
     %              errJ,'ELCdataJ.fold')
            if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,ydataH,
     %              errH,'ELCdataH.fold')
            if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,ydataK,
     %              errK,'ELCdataK.fold')
            if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,yRV1,
     %              errRV1,'ELCdataRV1.fold')
            if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,yRV2,
     %              errRV2,'ELCdataRV2.fold')
          endif
c
c   UPDATE August 10, 2004
c
c   Add 8 real and 4 integer variables to the list.
c 
c   UPDATE November 6, 2008
c
c   Add sw25-sw34 and isw25-isw34 to the list
c
          call writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     &       omega1,
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
     #       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,
     $       temprat,idark1,idark2,isw12,isw13,
     #       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2,
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
            call recordloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $        Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,
     %        RV2file,Nvmax,Nvar,
     $        svar,var,vstart,stepsave,Nstep,Nobv,sobv,obv,eobv)
c
            write(command,101)ig+1000000
            call system(command)
            write(command,102)ig+1000000
            call system(command)
c
            close(45)
            close(47)
c
c  UPDATE October 12, 2001
c
c  Close the new output file (unit 55).
c
            close(55)
c
10000     continue
c
          call indexx(Np,chiarr,indxchi)
c
c   reset the variables at their best values and print the chi^2
c
          var(1)=s1
          var(2)=s2
          var(3)=s3
          var(4)=s4
          var(5)=s5
          var(6)=s6
          var(7)=s7
          var(8)=s8
          var(9)=s9
          var(10)=s10
          var(11)=s11
          var(12)=s12
          var(13)=s13
          var(14)=s14
          var(15)=s15
          var(16)=s16
c
c  UPDATE August 13, 2001
c
c  Pad out the 's' variables to 25
c
          var(17)=s17
          var(18)=s18
          var(19)=s19
          var(20)=s20
          var(21)=s21
          var(22)=s22
          var(23)=s23
          var(24)=s24
          var(25)=s25
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c
c   UPDATE October 28, 2002
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
          if(savesep.lt.0.0d0)separ=savesep
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
          ifastflag=0
c
          call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $      omega2,Q,finc,Teff1,Teff2,betarim,
     $      rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
     &       pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &       dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
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
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c
c  February 12, 2000
c
c  Save the filling factors in case there is an eccentric orbit
c
          remf1=fill1
          remf2=fill2
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or T0 (itime>0), then we
c   have to fold the saved arrays.
c
              if(itime.gt.0)then
                if(icnU.ne.430)call phasefold(Ndatamax,isavNU,
     $            savxdataU,savydataU,saverrU,
     #            Nbin,NdataU,xdataU,ydataU,errU,period,T0,ikeep,ecc,argper,isw7)
                if(icnB.ne.430)call phasefold(Ndatamax,isavNB,
     $            savxdataB,savydataB,saverrB,
     #            Nbin,NdataB,xdataB,ydataB,errB,period,T0,ikeep,ecc,argper,isw7)
                if(icnV.ne.430)call phasefold(Ndatamax,isavNV,
     $            savxdataV,savydataV,saverrV,
     #            Nbin,NdataV,xdataV,ydataV,errV,period,T0,ikeep,ecc,argper,isw7)
                if(icnR.ne.430)call phasefold(Ndatamax,isavNR,
     $            savxdataR,savydataR,saverrR,
     #            Nbin,NdataR,xdataR,ydataR,errR,period,T0,ikeep,ecc,argper,isw7)
                if(icnI.ne.430)call phasefold(Ndatamax,isavNI,
     $            savxdataI,savydataI,saverrI,
     #            Nbin,NdataI,xdataI,ydataI,errI,period,T0,ikeep,ecc,argper,isw7)
                if(icnJ.ne.430)call phasefold(Ndatamax,isavNJ,
     $            savxdataJ,savydataJ,saverrJ,
     #            Nbin,NdataJ,xdataJ,ydataJ,errJ,period,T0,ikeep,ecc,argper,isw7)
                if(icnH.ne.430)call phasefold(Ndatamax,isavNH,
     $            savxdataH,savydataH,saverrH,
     #            Nbin,NdataH,xdataH,ydataH,errH,period,T0,ikeep,ecc,argper,isw7)
                if(icnK.ne.430)call phasefold(Ndatamax,isavNK,
     $            savxdataK,savydataK,saverrK,
     #            Nbin,NdataK,xdataK,ydataK,errK,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV1.ne.430)call phasefold(Ndatamax,isavRV1,
     $            savxRV1,savyRV1,saverrRV1,
     #            Nbin,NRV1,xRV1,yRV1,errRV1,period,T0,ikeep,ecc,argper,isw7)
                if(icnRV2.ne.430)call phasefold(Ndatamax,isavRV2,
     $            savxRV2,savyRV2,saverrRV2,
     #            Nbin,NRV2,xRV2,yRV2,errRV2,period,T0,ikeep,ecc,argper,isw7)
              endif
c
c
c  Get the models.
c
c
c  UPDATE November 6, 2008
c
c  If sw25>0, use that as an error for asini (sw5)
c
          if(sw25.gt.0.0d0)sw5=saveasini+gasdev(idum)*sw25
c
          call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $         ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
          fill1=remf1
          fill2=remf2
c
c   NOTE March 19, 2002
c
c   These calls to checklcfit have a zero at the end of the argument, 
c   which is a place holder for
c   'ifixgamma'. For
c   light curves (first argument '9'), this parameter is ignored.
c
c   Get the chi^2 values
c
          if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %      ydataU,errU,chisqU,zeroU,0,isw7)
          if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %      ydataB,errB,chisqB,zeroB,0,isw7)
          if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %      ydataV,errV,chisqV,zeroV,0,isw7)
          if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %      ydataR,errR,chisqR,zeroR,0,isw7)
          if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %      ydataI,errI,chisqI,zeroI,0,isw7)
          if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %      ydataJ,errJ,chisqJ,zeroJ,0,isw7)
          if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %      ydataH,errH,chisqH,zeroH,0,isw7)
          if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %      ydataK,errK,chisqK,zeroK,0,isw7)
          if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
            call checkRVfit(islc,NRVphase,xRVmod,RV1,NRVphase,xRVmod,RV2,
     #             NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     &             chisqRV1,chisqRV2,gam,ifixgamma)
            gamma1=gam
            gamma2=gam
          else
            if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,RV1,NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,ggamma,ifixgamma,isw7)
            gamma1=ggamma
            if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,ggamma,ifixgamma,isw7)
            gamma2=ggamma
          endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
          call writevar(Nvmax,svar,var,fill1,fill2,omega1,
     $        omega2,Q,finc,Teff1,Teff2,betarim,
     $        rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &        t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $        period,T0,alb1,alb2,
     #        dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
          ochi=0.0d0
          if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio.
c
              if(ifrac.gt.10)call diskchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)
c
              if(ilum.gt.10)call lrchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
          chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
          call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
c
          if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c   NEW BUG  July 24, 2001
c
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 4, 2005
c
c   Add in the luminosity ratio
c
              if(ifrac.gt.10)call wdiskobschi(ochidisk,Nobv,sobv)
              if(ilum.gt.10)call wlrobschi(ochilr)
c
          write(*,*)' '
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add dwavex and dwavey (limb darkening coefficients) to
c   the argument list of writevar.
c
          call writevar(Nvmax,fakevar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
          write(*,*)' '
          write(*,1002)chiall
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
          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(idint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',isw7)
c
          call wlinmod(NRVphase,xRVmod,RV1,'star1.RV',isw7)
          call wlinmod(NRVphase,xRVmod,RV2,'star2.RV',isw7)
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or the T0, write the current
c   folded light curves
c
          if(itime.gt.0)then
            if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,ydataU,
     %              errU,'ELCdataU.fold')
            if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,ydataB,
     %              errB,'ELCdataB.fold')
            if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,ydataV,
     %              errV,'ELCdataV.fold')
            if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,ydataR,
     %              errR,'ELCdataR.fold')
            if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,ydataI,
     %              errI,'ELCdataI.fold')
            if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,ydataJ,
     %              errJ,'ELCdataJ.fold')
            if(icnH.ne.430) call wELCdata(Ndatamax,NdataH,xdataH,ydataH,
     %              errH,'ELCdataH.fold')
            if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,ydataK,
     %              errK,'ELCdataK.fold')
            if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,yRV1,
     %              errRV1,'ELCdataRV1.fold')
            if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,yRV2,
     %              errRV2,'ELCdataRV2.fold')
          endif
c
c   The output file is closed within the light curve subroutine now.
c
c          close(2)
          close(38)
c
 1001     format(1x,2(f15.7,2x),a15)
 1002     format(1x,'S = ',f15.6)
 1003     format(1x,8(f11.6,2x))

 101      format('cp gridELC.opt gridELC.',i7)
 102      format('cp gridELC.inp ELC.',i7)

          stop
          end
c
c &^&%$*$*$*$^#@^#^$^&
c
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     $     NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.e-16,RNMX=1.-EPS)
      INTEGER j,k,iivv(NTAB),iiyy
      SAVE iivv,iiyy
      DATA iivv /NTAB*0/, iiyy /0/
      if (idum.le.0.or.iiyy.eq.0) then
         idum=max(-idum,1)
         do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iivv(j)=idum
         enddo
         iiyy=iivv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iiyy/NDIV
      iiyy=iivv(j)
      iivv(j)=idum
      ran1=min(AM*iiyy,RNMX)
      return
      END
c
c     &&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine indexx(n,arr,indx)
c
          integer n,indx(n),m,nstack
          real*8 arr(n)
          parameter(m=7,nstack=50)
          integer i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
          real*8 a
c
          do 11 j=1,n
            indx(j)=j
 11       continue
c
          jstack=0
          l=1
          ir=n
 1        if(ir-l.lt.m)then
            do 13 j=l+1,ir
              indxt=indx(j)
              a=arr(indxt)
              do 12 i=j-1,1,-1
                if(arr(indx(i)).le.a)go to 2
                indx(i+1)=indx(i)
 12           continue
              i=0
 2            indx(i+1)=indxt
 13         continue
            if(jstack.eq.0)return
            ir=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
          else
            k=(l+ir)/2
            itemp=indx(k)
            indx(k)=indx(l+1)
            indx(l+1)=itemp
            if(arr(indx(l+1)).gt.arr(indx(ir)))then
              itemp=indx(l+1)
              indx(l+1)=indx(ir)
              indx(ir)=itemp
            endif
            if(arr(indx(l)).gt.arr(indx(ir)))then
              itemp=indx(l)
              indx(l)=indx(ir)
              indx(ir)=itemp
            endif
            if(arr(indx(l+1)).gt.arr(indx(l)))then
              itemp=indx(l+1)
              indx(l+1)=indx(l)
              indx(l)=itemp
            endif
            i=l+1
            j=ir
            indxt=indx(l)
            a=arr(indxt)
 3          continue
            i=i+1
            if(arr(indx(i)).lt.a)go to 3
 4          continue
            j=j-1
            if(arr(indx(j)).gt.a)go to 4
            if(j.lt.i) go to 5
            itemp=indx(i)
            indx(i)=indx(j)
            indx(j)=itemp
            go to 3
 5          indx(l)=indx(j)
            indx(j)=indxt
            jstack=jstack+2
            if(jstack.gt.nstack)pause 'nstack too small in indexx'
            if(ir-i+1.ge.j-l)then
              istack(jstack)=ir
              istack(jstack-1)=i
              ir=j-1
            else
              istack(jstack)=j-1
              istack(jstack-1)=l
              l=i
            endif
          endif
          go to 1
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine rnkpop(n,arr,indx,irank)
c
          integer n,indx(n),m,nstack,irank(n)
          real*8 arr(n)
          parameter(m=7,nstack=50)
          integer i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
          real*8 a
c
          call indexx(n,arr,indx)
c
          do 1 i=1,n
            irank(indx(i))=n-i+1
 1        continue
          return
          end
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine select(np,jfit,fdif,idad)
c
          implicit double precision(a-h,o-z)
c
          dimension jfit(np)
c
          common /ranblock/ idum
c
          np1=np+1
          urand=ran1(idum)
          dice=urand*dble(np*np1)
          rtfit=0.0d0
          do 1 i=1,np
            rtfit=rtfit+dble(np1)+fdif*dble(np1-2*jfit(i))
            if(rtfit.ge.dice)then
              idad=i
              go to 2
            endif
 1        continue
 2        return
          end
c
c  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine encode(n,nd,ph,ign)
c
          double precision ph(n),z
          dimension ign(nd*n)
c
          z=10.0d0**nd
          ii=0
c
          do 1 i=1,n
            ip=dint(ph(i)*z)
            do 2 j=nd,1,-1
              ign(ii+j)=mod(ip,10)
              ip=ip/10
 2          continue
            ii=ii+nd
 1        continue
          return
          end
c
c   *************************************
c
          subroutine cross(n,nd,pcross,ign1,ign2)
c
c          double precision pcross,urand
          implicit double precision (a-h,o-z)

          dimension ign1(nd*n),ign2(nd*n)
c
          common /ranblock/ idum
c
          urand=ran1(idum)
          if(urand.lt.pcross)then
            urand=ran1(idum)
            ispl=dint(urand*dble(n*nd))+1
            do 10 i=ispl,n*nd
              it=ign2(i)
              ign2(i)=ign1(i)
              ign1(i)=it
 10         continue
          endif
          return
          end
c
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine mutate(n,nd,pmut,ign)
c
c          double precision pmut,urand
          implicit double precision(a-h,o-z)
          dimension ign(n*nd)
c
          common /ranblock/ idum

          do 10 i=1,n*nd
            urand=ran1(idum)
            if(urand.lt.pmut)then
              urand=ran1(idum)
              ign(i)=int(dint(urand*10.0d0))
            endif
 10       continue
          return
          end
c
c   ##################################
c
          subroutine decode(n,nd,ign,ph)
c
          dimension ign(n*nd)
          double precision z,ph(n)
c
          z=10.0d0**(-nd)
          ii=0
          do 1 i=1,n
            ip=0
            do 2 j=1,nd
              ip=10*ip+ign(ii+j)
 2          continue
            ph(i)=dble(ip)*z
            ii=ii+nd
 1        continue
          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine genrep(ndim,n,np,ip,ph,rnewph)
c
          double precision ph(ndim,2),rnewph(ndim,np)
c
          i1=2*ip-1
          i2=i1+1
          do 1 k=1,n
            rnewph(k,i1)=ph(k,1)
            rnewph(k,i2)=ph(k,2)
 1        continue
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&
c 
          subroutine adjmut(np,fitns,ifit,pmutmn,pmutmx,pmut)
c
          implicit double precision (a-h,o-z)
          double precision pmutmn,pmutmx,pmut,fitns,rat,rprod
c
          dimension fitns(np),ifit(np)
c
          rdiflo=0.05d0
          rdifhi=0.25d0
          delta=1.5d0
c
          rdif=dabs(fitns(ifit(np))-fitns(ifit(np/2)))/
     $      dabs(fitns(ifit(np))+fitns(ifit(np/2)))
c
          if(rdif.le.rdiflo)then
            rprod=pmut*delta
c            pmut=min(pmutmx,prod)
            if(pmutmx.lt.rprod)pmut=pmutmx
            if(rprod.le.pmutmx)pmut=rprod 
          endif
          if(rdif.ge.rdifhi)then
            rat=pmut/delta 
c            pmut=max(pmutmn,rat)
            if(pmutmn.gt.rat)pmut=pmutmn
            if(rat.ge.pmutmn)pmut=rat
          endif
c
 100      format(6(f10.6,2x))
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   NEW BUG ALERT July 10, 2001
c
c   Add this function 'gasdev'
c
          function gasdev(idum)
c
          implicit double precision(a-h,o-z)

          save iset,gset
          data iset/0/
c
          if(iset.eq.0)then
 1          v1=2.0d0*ran1(idum)-1.0d0
            v2=2.0d0*ran1(idum)-1.0d0
            rsq=v1**2+v2**2
            if(rsq.ge.1.0d0.or.rsq.eq.0.0d0)go to 1
            fac=dsqrt(-2.0d0*dlog(rsq)/rsq)
            gset=v1*fac
            gasdev=v2*fac
            iset=1
          else
            gasdev=gset
            iset=0
          endif
          return
          end
c
c &&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
