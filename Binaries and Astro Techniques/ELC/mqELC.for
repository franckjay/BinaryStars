         program mqELC
c
c   December 24, 1999
c
c   This is an optimizer code based on the Levenberg-Marquardt method.
c   The code is from Bevington (1969).
c
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
          implicit double precision (a-h,o-z)
          parameter(Nmaxphase=50000,Ndatamax=100000)
          parameter(Nvmax=25)
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm, obv, sobv and eobv to 9.
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
c
          dimension xmod(Nmaxphase),ymodU(Nmaxphase),ymodB(Nmaxphase),
     $      ymodV(Nmaxphase),ymodR(Nmaxphase),ymodI(Nmaxphase),
     $      ymodJ(Nmaxphase),ymodH(Nmaxphase),ymodK(Nmaxphase),
     &      ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      RV1(Nmaxphase),RV2(Nmaxphase),ymods3(Nmaxphase)
          dimension ymodU1(Nmaxphase),ymodB1(Nmaxphase),
     $      ymodV1(Nmaxphase),ymodR1(Nmaxphase),ymodI1(Nmaxphase),
     $      ymodJ1(Nmaxphase),ymodH1(Nmaxphase),ymodK1(Nmaxphase),
     &      RV11(Nmaxphase),RV21(Nmaxphase)
          dimension ymodU2(Nmaxphase),ymodB2(Nmaxphase),
     $      ymodV2(Nmaxphase),ymodR2(Nmaxphase),ymodI2(Nmaxphase),
     $      ymodJ2(Nmaxphase),ymodH2(Nmaxphase),ymodK2(Nmaxphase),
     &      RV12(Nmaxphase),RV22(Nmaxphase)
          dimension xdataU(Ndatamax),ydataU(Ndatamax),errU(Ndatamax),
     &      ydataB(Ndatamax),errB(Ndatamax),ydataV(Ndatamax),errV(Ndatamax),
     &      ydataR(Ndatamax),errR(Ndatamax),ydataI(Ndatamax),errI(Ndatamax),
     &      ydataJ(Ndatamax),errJ(Ndatamax),ydataH(Ndatamax),errH(Ndatamax),
     $      ydataK(Ndatamax),errK(Ndatamax),yRV1(Ndatamax),errRV1(Ndatamax),
     %      yRV2(Ndatamax),errRV2(Ndatamax),xdataB(Ndatamax),xdataV(Ndatamax),
     &      xdataR(Ndatamax),xdataI(Ndatamax),xdataJ(Ndatamax),
     &      xdataH(Ndatamax),xdataK(Ndatamax),xRV1(Ndatamax),xRV2(Ndatamax)
         dimension riden(Nvmax,Nvmax),beta(Nvmax),alpha(Nvmax,Nvmax),
     $      b(Nvmax),sigmaa(Nvmax),array(Nvmax,Nvmax),deriv(10000,Nvmax),
     $      drv1(Nmaxphase),drv2(Nmaxphase)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %      sigvar(Nvmax)
          dimension xRVmod(Nmaxphase)
          character*40 Udatafile,svar(19),Hdatafile,Kdatafile,RV1file,RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*40 fakevar(19),sobv(17)  
c
c   RVG BUG ALERT  June 12, 2001
c
c   Dimension the variables needed for the model atmospheres here.
c
          parameter (maxlines=1100,maxmu=99)
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
          dimension compfracs(8,2),ochidisk(8)   
          dimension fracs1(Nmaxphase,3),fracs2(Nmaxphase,3)
          dimension fracs3(Nmaxphase,3),fracs4(Nmaxphase,3)
          dimension fracs5(Nmaxphase,3),fracs6(Nmaxphase,3)
          dimension fracs7(Nmaxphase,3),fracs8(Nmaxphase,3)
c
          common /fracblock/ compfracs

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
c   Add sw21-sw24, powercoeff below.
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
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,sw33,sw34

c
c   RVG BUG ALERT   May 8, 2001
c
c   Add this common block to all programs
c
          common /spotblock/ spot1parm,spot2parm,spotdparm
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
c   RVG BUG ALERT  June 12, 2001
c
c   Add these common blocks for the model atmosphere variables
c
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,atmint4,
     &      atmint5,atmint6,atmint7,atmint8,
     &      Tmax,Tmin,gmax,gmin
          common /intatm/  Nlines,Nmu
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
        ifastflag=0   !disable fast genetic mode
c
c   RVG BUG ALERT  May 9, 2001
c
c   Add the spot parameters to getinput and recordparm
c
c   UPDATE August 10, 2004
c
c   Add the 8 real and 4 integer variables below.
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
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,sw33,sw34,
     %       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
c
c   May 16, 2001
c
c   This routine is called within the lightcurve subroutine now.
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
     &         atmin1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,
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
c   Disable median fitting.
c
          rmed=0.0d0
c
c   Disable the draw option
c
          idraw=0
c
c   NEW BUG ALERT July 13, 2001
c
c   Define the ifixgamma flag
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
c   UPDATE JULY 7, 2004
c
c   Initialize the Sobel sequence here.
c
          nnn=-11
          call sobseq(nnn,xsob)
c
c   Get the parameters for the grid search.
c
          call getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,Nvmax,Nvar,
     $      svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)
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
          if(t3.gt.0.0d0.and.g3.gt.0.0d0.and.SA3.gt.0.0d0.and.iatm.gt.0)then
            fakevar(13)='t3'
            fakevar(14)='g3'
            fakevar(15)='SA3'
          endif
c
c   Start the grid search here.  This code is more or less from Bevington
c   (1969).
c
          chi1=0.0d0
          chi2=0.0d0
          chi3=0.0d0
          flamda=0.001d0
          small=1000000.0d0
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+NdataK+
     %        NRV1+NRV2
c
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
 749      if(Nterms.eq.0)go to 69   !no valid variables specified
c
          write(*,*)Nterms,Nvar,Ndattot
          do 750 ijk=1,Nstep(1)
c
 31         do 34 j=1,nterms
              beta(j)=0.0d0
              do 34 k=1,j
 34         alpha(j,k)=0.0d0
c
c
c  Find the chi square at first point.
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
     %       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
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
            call initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &            chisqH,chisqK,chisqRV1,chisqRV2)
c
c   UPDATE March 19, 2002
c
c   Add 'ifixgamma' to the argument list for checklcfit.  For light
c   curves (first argument '9'), this parameter is ignored.
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
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     &          t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $          period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
            ochi=0.0d0
            if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
            endif
c 
            do 50 j=1,nterms
              delta=vstep(j) 
              oldvar=var(j)
              var(j)=var(j)+delta
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
     %       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
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
              call lightcurve(Nphase,Nmaxphase,xmod,ymodU1,ymodB1,
     $          ymodV1,ymodR1,ymodI1,ymodJ1,ymodH1,ymodK1,ymods1,ymods2,
     %          ymods3,ymodd,RV11,RV21,drv1,drv2,obsparm,ifastflag,
     &          NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
              fill1=remf1
              fill2=remf2
c
c   Get the chi^2 values
c
c   UPDATE March 19, 2002
c
c   Add 'ifixgamma' to the argument list of checklcfit.  For light curves
c   (first argument '9'), this parameter is ignored.
c
              if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU1,NdataU,
     %          xdataU,ydataU,errU,chisqU,zeroU1,ifixgamma,isw7)
              if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB1,NdataB,
     &          xdataB,ydataB,errB,chisqB,zeroB1,ifixgamma,isw7)
              if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV1,NdataV,
     %          xdataV,ydataV,errV,chisqV,zeroV1,ifixgamma,isw7)
              if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR1,NdataR,
     &          xdataR,ydataR,errR,chisqR,zeroR1,ifixgamma,isw7)
              if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI1,NdataI,
     &          xdataI,ydataI,errI,chisqI,zeroI1,ifixgamma,isw7)
              if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ1,NdataJ,
     &          xdataJ,ydataJ,errJ,chisqJ,zeroJ1,ifixgamma,isw7)
              if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH1,NdataH,
     %          xdataH,ydataH,errH,chisqH,zeroH1,ifixgamma,isw7)
              if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK1,NdataK,
     &          xdataK,ydataK,errK,chisqK,zeroK1,ifixgamma,isw7)
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
                 if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,NRV2,
     &             xRV2,
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
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma11,isvel2,
     %          gamma21,t3,g3,SA3,ecc,argper,pshift,spot1parm,
     &          spot2parm,spotdparm,period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
              ochi=0.0d0
              if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
              endif
c
              var(j)=oldvar
              var(j)=var(j)-delta
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
     $          rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
     &          pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
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
c
c  February 12, 2000
c
c  Save the filling factors in case there is an eccentric orbit
c
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
              call lightcurve(Nphase,Nmaxphase,xmod,ymodU2,ymodB2,
     $          ymodV2,ymodR2,ymodI2,ymodJ2,ymodH2,ymodK2,ymods1,ymods2,
     %          ymods3,
     %          ymodd,RV12,RV22,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8)
c
              fill1=remf1
              fill2=remf2
c
c   Get the chi^2 values
c   
c
c   UPDATE March 19, 2002
c
c   Add 'ifixgamma' to the argument list of checklcfit.  For light curves
c   (first argument '9'), this parameter is ignored.
c
              if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU2,NdataU,
     $           xdataU,ydataU,errU,chisqU,zeroU2,ifixgamma,isw7)
              if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB2,NdataB,
     $           xdataB,ydataB,errB,chisqB,zeroB2,ifixgamma,isw7)
              if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV2,NdataV,
     $           xdataV,ydataV,errV,chisqV,zeroV2,ifixgamma,isw7)
              if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR2,NdataR,
     $           xdataR,ydataR,errR,chisqR,zeroR2,ifixgamma,isw7)
              if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI2,NdataI,
     $           xdataI,ydataI,errI,chisqI,zeroI2,ifixgamma,isw7)
              if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ2,NdataJ,
     $           xdataJ,ydataJ,errJ,chisqJ,zeroJ2,ifixgamma,isw7)
              if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH2,NdataH,
     $           xdataH,ydataH,errH,chisqH,zeroH2,ifixgamma,isw7)
              if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK2,NdataK,
     $           xdataK,ydataK,errK,chisqK,zeroK2,ifixgamma,isw7)
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
                 if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,RV2,NRV2,
     %             xRV2,yRV2,errRV2,chisqRV2,gamma,ifixgamma,isw7)
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
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma12,isvel2,
     &          gamma22,t3,g3,SA3,ecc,argper,pshift,spot1parm,
     ^          spot2parm,spotdparm,period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
              ochi=0.0d0
              if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
              endif
c
c   Get the derivatives.
c
              jcount=0
              if(icnU.ne.430)then 
                do 6677 m=1,NdataU
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodU,
     $            NdataU,xdataU,ydataU,errU,m,yout,zeroU)
                  call evalfit(9,Nphase,xmod,ymodU1,
     $            NdataU,xdataU,ydataU,errU,m,yout1,zeroU1)
                  call evalfit(9,Nphase,xmod,ymodU2,
     $            NdataU,xdataU,ydataU,errU,m,yout2,zeroU2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     %             (ydataU(m)-yout)*deriv(jcount,j)/(errU(m)*errU(m))
                  do 6646 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errU(m)*errU(m))
 6646             continue
 6677           continue
              endif
c
              if(icnB.ne.430)then 
                do 6678 m=1,NdataB
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodB,
     $            NdataB,xdataB,ydataB,errB,m,yout,zeroB)
                  call evalfit(9,Nphase,xmod,ymodB1,
     $            NdataB,xdataB,ydataB,errB,m,yout1,zeroB1)
                  call evalfit(9,Nphase,xmod,ymodB2,
     $            NdataB,xdataB,ydataB,errB,m,yout2,zeroB2)
c
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataB(m)-yout)*deriv(jcount,j)/(errB(m)*errB(m))
                  do 6647 k=1,j
                    alpha(j,k)=alpha(j,k)+
     $                deriv(jcount,j)*deriv(jcount,k)/(errB(m)*errB(m))
 6647             continue
 6678           continue
              endif
c
              if(icnV.ne.430)then 
                do 6679 m=1,NdataV
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodV,
     $            NdataV,xdataV,ydataV,errV,m,yout,zeroV)
                  call evalfit(9,Nphase,xmod,ymodV1,
     $            NdataV,xdataV,ydataV,errV,m,yout1,zeroV1)
                  call evalfit(9,Nphase,xmod,ymodV2,
     $            NdataV,xdataV,ydataV,errV,m,yout2,zeroV2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     %            (ydataV(m)-yout)*deriv(jcount,j)/(errV(m)*errV(m))
                  do 6648 k=1,j
                    alpha(j,k)=alpha(j,k)+
     $              deriv(jcount,j)*deriv(jcount,k)/(errV(m)*errV(m))
 6648             continue
 6679           continue
              endif
c
              if(icnR.ne.430)then 
                do 6680 m=1,NdataR
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodR,
     $            NdataR,xdataR,ydataR,errR,m,yout,zeroR)
                  call evalfit(9,Nphase,xmod,ymodR1,
     $            NdataR,xdataR,ydataR,errR,m,yout1,zeroR1)
                  call evalfit(9,Nphase,xmod,ymodR2,
     $            NdataR,xdataR,ydataR,errR,m,yout2,zeroR2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataR(m)-yout)*deriv(jcount,j)/(errR(m)*errR(m))
                  do 6649 k=1,j
                    alpha(j,k)=alpha(j,k)+
     $               deriv(jcount,j)*deriv(jcount,k)/(errR(m)*errR(m))
 6649             continue
 6680           continue
              endif
c
              if(icnI.ne.430)then 
                do 6681 m=1,NdataI
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodI,
     $            NdataI,xdataI,ydataI,errI,m,yout,zeroI)
                  call evalfit(9,Nphase,xmod,ymodI1,
     $            NdataI,xdataI,ydataI,errI,m,yout1,zeroI1)
                  call evalfit(9,Nphase,xmod,ymodI2,
     $            NdataI,xdataI,ydataI,errI,m,yout2,zeroI2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataI(m)-yout)*deriv(jcount,j)/(errI(m)*errI(m))
                  do 6650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errI(m)*errI(m))
6650              continue
6681            continue
              endif
c
              if(icnJ.ne.430)then 
                do 16681 m=1,NdataJ
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodJ,
     $            NdataJ,xdataJ,ydataJ,errJ,m,yout,zeroJ)
                  call evalfit(9,Nphase,xmod,ymodJ1,
     $            NdataJ,xdataJ,ydataJ,errJ,m,yout1,zeroJ1)
                  call evalfit(9,Nphase,xmod,ymodJ2,
     $            NdataJ,xdataJ,ydataJ,errJ,m,yout2,zeroJ2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataJ(m)-yout)*deriv(jcount,j)/(errJ(m)*errJ(m))
                  do 16650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errJ(m)*errJ(m))
16650           continue
16681         continue
              endif
c
              if(icnK.ne.430)then 
                do 26681 m=1,NdataK
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodK,
     $            NdataK,xdataK,ydataK,errK,m,yout,zeroK)
                  call evalfit(9,Nphase,xmod,ymodK1,
     $            NdataK,xdataK,ydataK,errK,m,yout1,zeroK1)
                  call evalfit(9,Nphase,xmod,ymodK2,
     $            NdataK,xdataK,ydataK,errK,m,yout2,zeroK2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataK(m)-yout)*deriv(jcount,j)/(errK(m)*errK(m))
                  do 26650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errK(m)*errK(m))
26650           continue
26681         continue
              endif
c
              if(icnH.ne.430)then 
                do 36681 m=1,NdataH
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodH,
     $               NdataH,xdataH,ydataH,errH,m,yout,zeroH)
                  call evalfit(9,Nphase,xmod,ymodH1,
     $               NdataH,xdataH,ydataH,errH,m,yout1,zeroH1)
                  call evalfit(9,Nphase,xmod,ymodH2,
     $               NdataH,xdataH,ydataH,errH,m,yout2,zeroH2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataH(m)-yout)*deriv(jcount,j)/(errH(m)*errH(m))
                  do 36650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errH(m)*errH(m))
36650           continue
36681         continue
              endif
c 
              if(icnRV1.ne.430)then 
                do 46681 m=1,NRV1
                  jcount=jcount+1
                  call evalfit(0,Nphase,xmod,RV1,
     $            NRV1,xRV1,yRV1,errRV1,m,yout,gamma1)
                  call evalfit(0,Nphase,xmod,RV11,
     $            NRV1,xRV1,yRV1,errRV1,m,yout1,gamma11)
                  call evalfit(0,Nphase,xmod,RV12,
     $            NRV1,xRV1,yRV1,errRV1,m,yout2,gamma12)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $              (xRV1(m)-yout)*deriv(jcount,j)/(errRV1(m)*
     #              errRV1(m))
                  do 46650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %                 deriv(jcount,j)*deriv(jcount,k)/(errRV1(m)*
     %                 errRV1(m))
46650             continue
46681           continue
              endif
c
              if(icnRV2.ne.430)then 
                do 56681 m=1,NRV2
                  jcount=jcount+1
                  call evalfit(0,Nphase,xmod,RV2,
     %            NRV2,xRV2,yRV2,errRV2,m,yout,gamma2)
                  call evalfit(0,Nphase,xmod,RV21,
     $            NRV2,xRV2,yRV2,RV2,m,yout1,gamma21)
                  call evalfit(0,Nphase,xmod,RV22,
     $            NRV2,xRV2,yRV2,errRV2,m,yout2,gamma22)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $               (xRV2(m)-yout)*deriv(jcount,j)/(errRV2(m)*
     %               errRV2(m))
                  do 56650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %                 deriv(jcount,j)*deriv(jcount,k)/(errRV2(m)*
     &                 errRV2(m))
56650             continue
56681           continue
              endif
c
              var(j)=oldvar
 50         continue               !continue the loop over Nterms
c
c   we have filled up the beta and alpha matrices
c
 51         do 53 j=1,nterms
            do 53 k=1,j
 53         alpha(k,j)=alpha(j,k)
c
c   invert the matrices
c
 71         do 74 j=1,nterms
            do 73 k=1,nterms
            riden(j,k)=0.0d0
 73         array(j,k)=alpha(j,k)/sqrt(alpha(j,j)*alpha(k,k))
 74         array(j,j)=(1.0d0+flamda) 
            do 652 j=1,nterms
              riden(j,j)=1.0d0
 652        continue
            call gaussj(array,Nterms,Nvmax,riden,Nterms,Nvmax)
c
 81         do 84 j=1,nterms
            b(j)=var(j)
            do 84 k=1,nterms
 84         b(j)=b(j)+beta(k)*array(j,k)/sqrt(alpha(j,j)*alpha(k,k))
c
c   Evaluate the function and chi^2.  Note the array 'b' is used instead of
c   'var'.
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
            call assignvar(Nvmax,svar,b,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,SA3,ecc,argper,
     &         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &         dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %         ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
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
c
c   UPDATE March 19, 2002
c
c   Add 'ifixgamma' to the argument list of checklcfit.  For light curves
c   (first argument '9'), this parameter is ignored.
c
            if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %         ydataU,errU,chisqU,zeroU2,ifixgamma,isw7)
            if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %         ydataB,errB,chisqB,zeroB2,ifixgamma,isw7)
            if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %         ydataV,errV,chisqV,zeroV2,ifixgamma,isw7)
            if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %         ydataR,errR,chisqR,zeroR2,ifixgamma,isw7)
            if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %         ydataI,errI,chisqI,zeroI2,ifixgamma,isw7)
            if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %         ydataJ,errJ,chisqJ,zeroJ2,ifixgamma,isw7)
            if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %         ydataH,errH,chisqH,zeroH2,ifixgamma,isw7)
            if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %         ydataK,errK,chisqK,zeroK2,ifixgamma,isw7)
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
            call writevar(Nvmax,svar,b,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma12,isvel2,
     &          gamma22,t3,g3,SA3,ecc,argper,pshift,spot1parm,
     %          spot2parm,spotdparm,period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
            ochi=0.0d0
            if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
            chi4=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
            write(*,*)'chi4 = ',chi4
            chi4=chi4/dabs(dble(Ndattot-Nterms))
            call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
            if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
            if(chi4.lt.small)then
              small=chi4
              s1=b(1)
              s2=b(2)
              s3=b(3)
              s4=b(4)
              s5=b(5)
              s6=b(6)
              s7=b(7)
              s8=b(8)
              s9=b(9)
              s10=b(10)
              s11=b(11)
              s12=b(12)
              s13=b(13)
              s14=b(14)
              s15=b(15)
              s16=b(16)
            endif
            diff=chi1-chi4
c
            if(diff) 95, 101, 101
 95         flamda=10.0d0*flamda
            write(*,*)'flamda = ',flamda
            if(flamda.gt.1.0d12)stop
            go to 71
c
 101       write(*,*)' '
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
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,
     &          gamma2,t3,g3,SA3,ecc,argper,pshift,spot1parm,
     &          spot2parm,spotdparm,period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
           write(*,*)' '
c
            do 103 j=1,nterms
              var(j)=b(j)
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
              vstep(j)=0.9*vstep(j)
              sigmaa(j)=sqrt(array(j,j)/alpha(j,j))
              write(*,1001)var(j),sigmaa(j),svar(j)(1:15)
 103        continue
c
c
            flamda=flamda/10.0
            write(*,*)'flamda = ',flamda
c
            write(*,2002)small
            call wlinmod(Nphase,xmod,ymodU,'modelU.linear')
            call wlinmod(Nphase,xmod,ymodB,'modelB.linear')
            call wlinmod(Nphase,xmod,ymodV,'modelV.linear')
            call wlinmod(Nphase,xmod,ymodR,'modelR.linear')
            call wlinmod(Nphase,xmod,ymodI,'modelI.linear')
            call wlinmod(Nphase,xmod,ymodJ,'modelJ.linear')
            call wlinmod(Nphase,xmod,ymodH,'modelH.linear')
            call wlinmod(Nphase,xmod,ymodK,'modelK.linear')
            call wlinmod(NRVphase,xRVmod,RV1,'star1.RV')
            call wlinmod(NRVphase,xRVmod,RV2,'star2.RV')
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
 750      continue  ! loop over iterations
c
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
     &         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &       dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
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
c   RVG BUG ALERT   April 26, 2001
c
c   Set fill1 to remf1.
c
          fill1=remf1
          fill2=remf2
c
c   Get the chi^2 values
c
c
c   UPDATE March 19, 2002
c
c   Add 'ifixgamma' to the argument list of checklcfit.  For light curves
c   (first argument '9'), this parameter is ignored.
c
          if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     %      ydataU,errU,chisqU,zeroU,ifixgamma,isw7)
          if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     %      ydataB,errB,chisqB,zeroB,ifixgamma,isw7)
          if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     %      ydataV,errV,chisqV,zeroV,ifixgamma,isw7)
          if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     %      ydataR,errR,chisqR,zeroR,ifixgamma,isw7)
          if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     %      ydataI,errI,chisqI,zeroI,ifixgamma,isw7)
          if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     %      ydataJ,errJ,chisqJ,zeroJ,ifixgamma,isw7)
          if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     %      ydataH,errH,chisqH,zeroH,ifixgamma,isw7)
          if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     %      ydataK,errK,chisqK,zeroK,ifixgamma,isw7)
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
     $        period,T0,alb1,alb2,
     #        dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
          ochi=0.0d0
          if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
          chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
          call writechi(icnU,icnB,icnV,icnR,icnI,icnK,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2)              
c
          if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
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
     $          rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,
     &          gamma2,t3,g3,SA3,ecc,argper,pshift,spot1parm,
     &          spot2parm,spotdparm,period,T0,alb1,alb2,
     #          dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &          ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj)
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
c     $        gpole1,gpole2,vrot1,vrot2,gscale1,gscale2)
c
          do 3749 iq=1,Nterms
            write(*,1001)var(iq),sigmaa(iq),svar(iq)(1:15)
 3749     continue
c

          write(*,1002)chiall
c
          call wlinmod(Nphase,xmod,ymodU,'modelU.linear')
          call wlinmod(Nphase,xmod,ymodB,'modelB.linear')
          call wlinmod(Nphase,xmod,ymodV,'modelV.linear')
          call wlinmod(Nphase,xmod,ymodR,'modelR.linear')
          call wlinmod(Nphase,xmod,ymodI,'modelI.linear')
          call wlinmod(Nphase,xmod,ymodJ,'modelJ.linear')
          call wlinmod(Nphase,xmod,ymodH,'modelH.linear')
          call wlinmod(Nphase,xmod,ymodK,'modelK.linear')
c
          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear')
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear')
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear')
          if(idint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear')
c
          call wlinmod(NRVphase,xRVmod,RV1,'star1.RV')
          call wlinmod(NRVphase,xRVmod,RV2,'star2.RV')
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
c   UPDATE June 17, 2002
c
c   Add a comma after dwavey!
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
     #       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       idark1,idark2,isw12,isw13,
     #       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,sw33,sw34,
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
          call recordloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,Nvmax,Nvar,
     $      svar,var,vstart,stepsave,Nstep,Nobv,sobv,obv,eobv)
c
 1001     format(1x,2(f15.7,2x),a15)
 1002     format(1x,'S = ',f15.6)
 2002     format(1x,'smallest reduced chi^2 = ',f15.6)
 1003     format(1x,8(f11.6,2x))

 69       stop
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
