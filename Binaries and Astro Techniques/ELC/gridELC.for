          program gridELC
c
c   December 24, 1999
c
c   This is an optimizer code based on the 'grid search' technique
c   discussed in Bevington (1969).  This program will call the
c   subroutines contained within 'lcsubs.for' and 'optimizesubs.for'.
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
          implicit double precision (a-h,o-z)
c
          parameter(Nmaxphase=50000,Ndatamax=200000)
          parameter(Nvmax=25)
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm,obv,eobv,sobv to 9.
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
     &      xdataH(Ndatamax),xdataK(Ndatamax),xRV1(Ndatamax),xRV2(Ndatamax)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %      sigvar(Nvmax),stepsave(Nvmax),
     $      drv1(Nmaxphase),drv2(Nmaxphase)
          dimension xRVmod(Nmaxphase)
c
          character*40 Udatafile,svar(Nvmax),Hdatafile,Kdatafile,RV1file,
     #                 RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*40 fakevar(Nvmax),sobv(17)
c
          character*300 line        ! was 132
c
          character*259 parmstring
c
          common /stringblock/ parmstring
c
c   RVG BUG ALERT  June 12, 2001
c
c   Dimension the variables needed for the model atmospheres here.
c
          parameter (maxlines=1300,maxmu=115)
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
c   UPDATE OCTOBER 10, 2007
c
c   Add this array to keep track of the light curves of the individual
c   components in all band passes.
c
          dimension compfracs(8,2),ochidisk(8)
          common /fracblock/ compfracs          
c
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
c   Add sw21-sw24, powercoeff below
c
c   UPDATE November 6, 2008
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
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   Disable median fitting (if selected by sw6 > 0)
c
         rmed=0.0d0
c
        ifastflag=0   !disable fast genetic mode
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
c
c   RVG BUG ALERT  May 9, 2001
c
c   Add the spot parameters to the arguments of getinput and recordparm.
c
c   UPDATE August 10, 2004
c
c   Add the 8 real and 4 integer variables to the list below.
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
c   May 16, 2001
c
c   This subroutine is now called within the lightcurve subroutine
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
c
c   RVG BUG ALERT  June 12, 2001
c
c   Load the atmosphere table here, rather than in the subroutine
c
c  If the flag iatm>0, then load the model atmosphere table.
c
          if(iatm.ge.1)then
            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &      atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &      Tmax,Tmin,gmax,gmin)
          endif

c
c   NEW BUG August 2, 2001
c
c   Use these flags when fitting for the period and/or T0
c
          itime=isw7
          Nbin=isw8
c
c   Do we force the gamma to be the input value?
c
          ifixgamma=isw4
c
c
c   UPDATE JULY 30, 2004
c
c   set Nbin = 0
c
          Nbin=0
c
c   UPDATE JULY 7, 2004
c
c   Initialize the Sobel sequence here.
c
          nnn=-11
          call sobseq(nnn,xsob)
c
c
c   UPDATE April 15, 2002
c
c   Define the values of gamma1 and gamma2
c
          gamma1=gamma
          gamma2=gamma
c
c
c   Disable the draw option
c
          idraw=0
c
c   Fix the inner disk radius at fill2.
c
          if(teff2.gt.0.0d0)rinner=fill2
c
c   UPDATE AUGUST 4, 2004
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              savesep=separ
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
c
c   Get the parameters for the grid search.
c
          call getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,Nvmax,Nvar,
     $      svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)

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
          gamma1=gamma
          gamma2=gamma
c
c   Define the variables whose values will be printed on the screen at the
c   end of the iteration.
c
          fakevar(1)='mass'
          fakevar(2)='incl'
          fakevar(3)='t1'
          if(teff2.gt.1.0)then
            fakevar(4)='t2'
            fakevar(6)='f2'
            fakevar(12)='o2'
          else
            fakevar(4)='nothing'
            fakevar(6)='nothing'
            fakevar(12)='lx'
          endif
          fakevar(5)='f1'
          if(idint.ge.1)then
            fakevar(7)='rinner'
            fakevar(8)='router'
            fakevar(9)='tdisk'
            fakevar(10)='beta'
            fakevar(11)='xi'
          else
            fakevar(7)='nothing'
          endif
          fakevar(13)='o1'
          fakevar(14)='separ'
c
c
          open(unit=45,file='generation.0000',status='unknown')
          open(unit=46,file='ELCparm.0000',status='unknown')
c
c
c   October 10, 2007
c
c   Add this additional output file if the flag is set.
c                      
          if(isw24.ge.1)open(unit=47,file='ELCratio.0000',status='unknown')
c
c
c   Start the grid search here.  This code is more or less from Bevington
c   (1969).
c
          chi1=0.0
          chi2=0.0
          chi3=0.0
          small=1000000.0
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+NdataK+
     %        NRV1+NRV2+Nobv
c
          Nterms=0
c
c   UPDATE May 30, 2002
c
c   Change the upper loop limit to max(1,Nvar)
c
          do 699 i7=1,max(1,Nvar)
            var(i7)=vstart(i7)     !initialize the variables
            kkk=icnvrt(svar(i7)(1:2))
            stepsave(i7)=vstep(i7)
            if(kkk.eq.430)then
              go to 749
            else
              Nterms=Nterms+1
            endif
 699      continue
c
c   UPDATE May 30, 2002
c
c   Change the if-then statement to Nvar.eq.0
c
 749      if(Nvar.eq.0)go to 69   !no valid variables specified
c
          write(*,*)'Nterms, Nvar, Ndattot',Nterms,Nvar,Ndattot
          i16=0
          do 750 ijk=1,Nstep(1)
            do 700 i7=1,Nterms    
c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
c   UPDATE AUGUST 4, 2004
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              if(savesep.lt.0.0d0)separ=savesep
c
              call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $         rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
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
c  Find the chi square at first two grid search points.
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
              call initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &            chisqH,chisqK,chisqRV1,chisqRV2)
c
              if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %         ydataU,errU,chisqU,zeroU,0,isw7)
              if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %         ydataB,errB,chisqB,zeroB,0,isw7)
              if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %         ydataV,errV,chisqV,zeroV,0,isw7)
              if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %         ydataR,errR,chisqR,zeroR,0,isw7)
              if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %         ydataI,errI,chisqI,zeroI,0,isw7)
              if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %         ydataJ,errJ,chisqJ,zeroJ,0,isw7)
              if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %         ydataH,errH,chisqH,zeroH,0,isw7)
              if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %         ydataK,errK,chisqK,zeroK,0,isw7)
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
     %            yRV1,errRV1,chisqRV1,gamma,ifixgamma,isw7)
                gamma1=gamma
                if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,
     %             NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,gamma,ifixgamma,isw7)
                gamma2=gamma
              endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
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
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     &          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
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
              write(*,*)'chi1 = ',chi1
              chi1=chi1/dabs(dble(Ndattot-Nterms))
              call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
              if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
                s17=var(17)
                s18=var(18)
                s19=var(19)
                s20=var(20)
                s21=var(21)
                s22=var(22)
                s23=var(23)
                s24=var(24)
                s25=var(25)
c                write(*,*)small,s1,s2,s3,s4,s5,s6,s7
              endif
c 
c   UPDATE September 3, 2004
c
c   Add the output like geneticELC has (e.g. ELCparm.0000
c   and generation.0000
c
              i16=i16+1
              chikkk=chi1*dabs(dble(Ndattot-Nterms))
              call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,chikkk,i16,line,alb1,alb2,
     %          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
              cccc=chikkk
              if(cccc.gt.9999999.0d0)cccc=9999999.0d0
              write(46,4646)i16,cccc,parmstring
              lll=lnblnk(line)
              write(45,5555)line(1:lll)
c
c   October 10, 2007
c
c   write the ratios if the flag is set
c
              if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ljk,1),ljk=1,8),
     #           (compfracs(ljk,2),ljk=1,8)
c
5566          format(i5,1x,f14.4,1x,16(1pe13.6,1x))
c
c   UPDATE October 22, 2008
c
c   Updated format 4646 below
c
 4646         format(i4,1x,f16.7,1x,a259)
 5555         format(a)        ! UPDATE August 13, 2001: was format(a132)

              fn=0.0
              delta=vstep(i7) 
              var(i7)=var(i7)+delta
c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
c   UPDATE AUGUST 4, 2004
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              if(savesep.lt.0.0d0)separ=savesep
c
              call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $         omega2,Q,finc,Teff1,Teff2,betarim,
     $         rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
     &         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
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
              if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %         ydataU,errU,chisqU,zeroU,0,isw7)
              if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %         ydataB,errB,chisqB,zeroB,0,isw7)
              if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %         ydataV,errV,chisqV,zeroV,0,isw7)
              if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %         ydataR,errR,chisqR,zeroR,0,isw7)
              if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %         ydataI,errI,chisqI,zeroI,0,isw7)
              if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %         ydataJ,errJ,chisqJ,zeroJ,0,isw7)
              if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %         ydataH,errH,chisqH,zeroH,0,isw7)
              if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %         ydataK,errK,chisqK,zeroK,0,isw7)
              if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
                call checkRVfit(islc,NRVphase,xRVmod,RV1,NRVphase,xRVmod,RV2,
     #             NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     &             chisqRV1,chisqRV2,gam,ifixgamma)
                gamma1=gam
                gamma2=gam
              else
                if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,
     &            RV1,NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,gamma,ifixgamma,isw7)
                gamma1=gamma
                if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,
     %             RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,gamma,ifixgamma,isw7)
                gamma2=gamma
              endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
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
              ochi=0.0d0
              if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
              chi2=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
              write(*,*)'chi2 = ',chi2
              chi2=chi2/dabs(dble(Ndattot-Nterms))
              call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
              if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)

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
              if(chi2.lt.small)then
                small=chi2
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
              i16=i16+1
              chikkk=chi2*dabs(dble(Ndattot-Nterms))
              call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,chikkk,i16,line,alb1,alb2,
     %          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
              cccc=chikkk
              if(cccc.gt.9999999.0d0)cccc=9999999.0d0
              write(46,4646)i16,cccc,parmstring
              lll=lnblnk(line)
              write(45,5555)line(1:lll)
c
c
c   October 10, 2007
c
c   Write the fractions if requested
c
              if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ljk,1),ljk=1,8),
     #           (compfracs(ljk,2),ljk=1,8)
c
              diff=chi1-chi2
              if(diff.gt.0.0)go to 5061
c
              delta=-delta
              var(i7)=var(i7)+delta
c  
              csave=chi1
              chi1=chi2
              chi2=csave
c 
c   incriment the variable var(i7) until chi square increases
c
 5061         fn=fn+1.0
c 
              var(i7)=var(i7)+delta
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
c   UPDATE AUGUST 4, 2004
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              if(savesep.lt.0.0d0)separ=savesep
c
              call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
     &          pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
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
              if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %         ydataU,errU,chisqU,zeroU,0,isw7)
              if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %         ydataB,errB,chisqB,zeroB,0,isw7)
              if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %         ydataV,errV,chisqV,zeroV,0,isw7)
              if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %         ydataR,errR,chisqR,zeroR,0,isw7)
              if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %         ydataI,errI,chisqI,zeroI,0,isw7)
              if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %         ydataJ,errJ,chisqJ,zeroJ,0,isw7)
              if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %         ydataH,errH,chisqH,zeroH,0,isw7)
              if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %         ydataK,errK,chisqK,zeroK,0,isw7)
              if((icnRV1.ne.430).and.(icnRV2.ne.430).and.(ifixgamma.ge.2))then
                call checkRVfit(islc,NRVphase,xRVmod,RV1,NRVphase,xRVmod,RV2,
     #             NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     &             chisqRV1,chisqRV2,gam,ifixgamma)
                gamma1=gam
                gamma2=gam
              else
                if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,
     &            RV1,NRV1,xRV1,
     %            yRV1,errRV1,chisqRV1,gamma,ifixgamma,isw7)
                gamma1=gamma
                if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,
     $             RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,gamma,ifixgamma,isw7)
                gamma2=gamma
              endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
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
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
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
              chi3=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
              write(*,*)'chi3 = ',chi3
              chi3=chi3/dabs(dble(Ndattot-Nterms))
              call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
              if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
              if(chi3.lt.small)then
                small=chi3
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
              i16=i16+1
              chikkk=chi3*dabs(dble(Ndattot-Nterms))
              call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,chikkk,i16,line,alb1,alb2,
     %          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
              cccc=chikkk
              if(cccc.gt.9999999.0d0)cccc=9999999.0d0
              write(46,4646)i16,cccc,parmstring
              lll=lnblnk(line)
              write(45,5555)line(1:lll)
c
c   October 10, 2007
c
c   Write the fractions if requested
c
              if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ljk,1),ljk=1,8),
     #           (compfracs(ljk,2),ljk=1,8)
c
              diff23=chi3-chi2
              if(diff23.lt.0)then
                chi1=chi2
                chi2=chi3
                go to 5061
              endif          
c
c  Find the minimum of parabola defined by the last three points.
c
 5081         if(chi3-chi2.eq.0.0)go to 700
              delta=delta*(1.0/(1.0+(chi1-chi2)/(chi3-chi2))+0.5)
              var(i7)=var(i7)-delta
              free=dabs(dble(Ndattot-Nterms))
              sigvar(i7)=vstep(i7)*sqrt(2.0/(free*(chi3-2.0*chi2+chi1)))
              vstep(i7)=vstep(i7)*fn/3.0
c
              if((Nterms.eq.1).or.(i7.eq.Nterms))then
c
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list
c
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
c   UPDATE AUGUST 4, 2004
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              if(savesep.lt.0.0d0)separ=savesep
c
                call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     $            omega2,Q,finc,Teff1,Teff2,betarim,
     $            rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
     &            pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &            dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %            ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
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
                call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $           ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     &           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
                fill1=remf1
                fill2=remf2
c
c   Get the chi^2 values
c
                if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,
     &            NdataU,xdataU,ydataU,errU,chisqU,zeroU,0,isw7)
                if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,
     &            NdataB,xdataB,ydataB,errB,chisqB,zeroB,0,isw7)
                if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,
     &            NdataV,xdataV,ydataV,errV,chisqV,zeroV,0,isw7)
                if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,
     &            NdataR,xdataR,ydataR,errR,chisqR,zeroR,0,isw7)
                if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,
     %            NdataI,xdataI,ydataI,errI,chisqI,zeroI,0,isw7)
                if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,
     &            NdataJ,xdataJ,ydataJ,errJ,chisqJ,zeroJ,0,isw7)
                if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,
     &            NdataH,xdataH,ydataH,errH,chisqH,zeroH,0,isw7)
                if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,
     &            NdataK,xdataK,ydataK,errK,chisqK,zeroK,0,isw7)
                if((icnRV1.ne.430).and.(icnRV2.ne.430).
     #               and.(ifixgamma.ge.2))then
                  call checkRVfit(islc,NRVphase,xRVmod,RV1,NRVphase,
     &               xRVmod,RV2,
     #               NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     &               chisqRV1,chisqRV2,gam,ifixgamma)
                  gamma1=gam
                  gamma2=gam
                else
                  if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,RV1,
     &              NRV1,xRV1,
     %              yRV1,errRV1,chisqRV1,gamma,ifixgamma,isw7)
                  gamma1=gamma
                  if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,
     &               NRV2,xRV2,yRV2,errRV2,chisqRV2,gamma,ifixgamma,isw7)
                  gamma2=gamma
                endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
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
     $            omega2,Q,finc,Teff1,Teff2,betarim,
     $            rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,
     &            gamma2,t3,g3,SA3,ecc,argper,pshift,spot1parm,
     %            spot2parm,spotdparm,period,T0,alb1,alb2,
     #            dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
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
                chi4=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
                write(*,*)'chi4 = ',chi4
                chi4=chi4/dabs(dble(Ndattot-Nterms))
                call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %            icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &            chisqJ,chisqH,chisqK,
     &            chisqRV1,chisqRV2)              
                if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
                if(chi4.lt.small)then
                  small=chi4
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
                i16=i16+1
                chikkk=chi4*dabs(dble(Ndattot-Nterms))
                call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $           omega2,Q,finc,Teff1,Teff2,betarim,
     $           rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %           t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $           period,T0,chikkk,i16,line,alb1,alb2,
     %           dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
                cccc=chikkk
                if(cccc.gt.9999999.0d0)cccc=9999999.0d0
                write(46,4646)i16,cccc,parmstring
                lll=lnblnk(line)
                write(45,5555)line(1:lll)
c
c   October 10, 2007
c
c   Write the fractions if requested
c
               if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ljk,1),ljk=1,8),
     #            (compfracs(ljk,2),ljk=1,8)
c
              endif

 700        continue  !loop over the 7 terms
c
            write(*,*)' '
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
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
     #          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
            write(*,*)' '
c
            do 7749 iq=1,Nterms
              write(*,1001)var(iq),sigvar(iq),svar(iq)(1:15)
 7749       continue
c
            write(*,1002)small
            call wlinmod(Nphase,xmod,ymodU,'modelU.linear',isw7)
            call wlinmod(Nphase,xmod,ymodB,'modelB.linear',isw7)
            call wlinmod(Nphase,xmod,ymodV,'modelV.linear',isw7)
            call wlinmod(Nphase,xmod,ymodR,'modelR.linear',isw7)
            call wlinmod(Nphase,xmod,ymodI,'modelI.linear',isw7)
            call wlinmod(Nphase,xmod,ymodJ,'modelJ.linear',isw7)
            call wlinmod(Nphase,xmod,ymodH,'modelH.linear',isw7)
            call wlinmod(Nphase,xmod,ymodK,'modelK.linear',isw7)
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
c
c   Reset the variables at their best values.
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
 750      continue  ! loop over grid searches
c
 751      continue                  ! come here when delta chi is small
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
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
c   UPDATE November 6, 2002
c
c   Add limb darkening coefficients to the argument list of
c   assignvar and varassign, i.e. dwavey and dwavey
c
c   UPDATE AUGUST 4, 2004
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              if(savesep.lt.0.0d0)separ=savesep
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
c  Get the models.
c
 69       call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $         ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
          fill1=remf1
          fill2=remf2
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
     %            yRV1,errRV1,chisqRV1,gamma,ifixgamma,isw7)
            gamma1=gamma
            if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,NRV2,xRV2,
     %             yRV2,errRV2,chisqRV2,gamma,ifixgamma,isw7)
            gamma2=gamma
          endif
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
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


          chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
          call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
c
          if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
          write(*,*)' '
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the argument list of writevar
c
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
     #          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
          write(*,*)' '
c
c          call parms(0,
c     &        Teff2,Q,finc,separ,period,reff1,reff2,
c     %        rpole1,rpole2,fill1,fill2,
c     $        gpole1,gpole2,vrot1,vrot2,gscale1,gscale2,omega1,omega2)
c
          do 3749 iq=1,Nterms
            write(*,1001)var(iq),sigvar(iq),svar(iq)(1:15)
 3749     continue
c
          write(*,1002)chiall
c
              i16=i16+1
              chikkk=chi4*dabs(dble(Ndattot-Nterms))
              call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,chiall,i16,line,alb1,alb2,
     %          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
              cccc=chiall
              if(cccc.gt.9999999.0d0)cccc=9999999.0d0
              write(46,4646)i16,cccc,parmstring
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
c          close(2)   ! close the output file
c
c   Finally, make a file similar to ELC.inp with the current parameters.
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
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
          call recordloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,Nvmax,Nvar,
     $      svar,var,vstart,stepsave,Nstep,Nobv,sobv,obv,eobv)
c
 1001     format(1x,2(f15.7,2x),a15)
 1002     format(1x,'S = ',f15.6)
 1003     format(1x,8(f11.6,2x))

          close(45)
          close(46)
          close(47)
          stop
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
