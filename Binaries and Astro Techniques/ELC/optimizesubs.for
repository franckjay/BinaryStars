c
c    November 12, 1999
c
c    These subroutines will perform the necessary tasks for the optimizer
c    routines, such as variable assignment, chi^2 checking, etc.
c
          subroutine getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,Nvmax,Nvar,
     $      svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)
c
c    November 12, 1999
c
c    This routine will read the input file to set up the optimizer (either
c    the 'grid search' or the Levenburg-Marquardt routine).
c
          implicit double precision (a-h,o-z)
c
c   UPDATE September 11, 2001
c
c   Change the dimensions of obs,eobv,sobv to 9
c
c   UPDATE September 21, 2008
c
c   make the dimension of obv, eobv, and sobv 11
c
c   UPDATE November 17, 2008
c
c   make the dimension of obv, epob, and sobv 17

          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %        obv(17),eobv(17)
          character*40 Udatafile,svar(Nvmax),Hdatafile,Kdatafile,
     %        RV1file,RV2file,sobv(17)
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*1 bell
c
c   UPDATE September 11, 2001
c
c   Change the limit to 9
c
c   UPDATE September 21, 2008
c
c   make the limit 11
c
          do 1 i=1,11
            sobv(i)=' '
            obv(i)=-99.
            eobv(i)=-99.
 1        continue
          Nobv=0
c
          do 2 i=1,Nvmax
            svar(i)='zz'
 2        continue
c
          bell=char(7)
          iso=0
          open(unit=11,file='gridloop.opt',status='old',err=1200,iostat=ios)
c
          read(11,100)Udatafile
          read(11,100)Bdatafile
          read(11,100)Vdatafile
          read(11,100)Rdatafile
          read(11,100)Idatafile
          read(11,100)Jdatafile
          read(11,100)Hdatafile
          read(11,100)Kdatafile
          read(11,100)RV1file
          read(11,100)RV2file
c
          read(11,*)Nvar
          if(Nvar.eq.0)go to 999
          if(Nvar.gt.Nvmax)then
            write(*,200)bell
            stop
          endif
c
          do 10 i=1,Nvar
            read(11,100)svar(i)
 10       continue
c
          do 20 i=1,Nvar
            read(11,*)vstart(i),vstep(i),Nstep(i)
            var(i)=vstart(i)
 20       continue
c
c  Look here for possible observed parameters
c
 999      read(11,*,end=30,err=30)Nobv
c
          do 25 i=1,Nobv
            read(11,100)sobv(i)
 25       continue
          do 26 i=1,Nobv
            read(11,*)obv(i),eobv(i)
 26       continue

 30       close(11)
c
          istop=0
 1200     if(ios.ne.0)then    ! error in opening the file
            call makeloopopt(Nvmax)
          endif
c
 100      format(a40)
 200      format(a1,'Error:  too many variables in ''gridloop.opt''')
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine makeloopopt(Nvmax)
c
c   This routine will write a correctly formatted 'gridloop.opt' file
c   in the event one does not exist.
c     
          implicit double precision (a-h,o-z)

          character*1 bell
c
          bell=char(7)
          write(*,100)bell
c 
          open(unit=11,file='gridloop.opt',status='new')
c
          write(11,1000)
          write(11,1001)
          write(11,1002)
          write(11,1003)
          write(11,1004)
          write(11,1005)
          write(11,1006)
          write(11,1007)
          write(11,1008)
          write(11,1009)
c
          Nvar=Nvmax
          write(11,1010)Nvar
c
          write(11,1011)
          write(11,1012)
          write(11,1013)
          write(11,1014)
          write(11,1015)
          write(11,1016)
          write(11,1017)
          write(11,1018)
          write(11,1019)
          write(11,1020)
          write(11,1021)
          write(11,1022)
          write(11,1023)
          write(11,1024)
          write(11,1025)
          write(11,1026)
          write(11,1030)
          write(11,1031)
          write(11,1032)
          write(11,1033)
          write(11,1034)
          write(11,1035)
          write(11,1036)
          write(11,1037)
          write(11,1038)
          write(11,1039)
          write(11,1040)
          write(11,1041)
          write(11,1042)
          write(11,1043)
          write(11,1044)
          write(11,1045)
c
          write(11,2000)
          write(11,2001)
          write(11,2002)
          write(11,2003)
          write(11,3004)
          write(11,2004)
          write(11,2005)
          write(11,2006)
          write(11,2007)
          write(11,2008)
          write(11,2011)
          write(11,2012)
          write(11,2014)
          write(11,2015)
c
          idummy=10
          if(idummy.eq.10)then
            write(*,2020)
            stop
          endif
c
 100      format(a1,'Error:  file ''gridloop.opt'' not found!  ',
     %      'I''m making one up!')
c
 1000     format('put_your_file_for_U_data_here')
 1001     format('put_your_file_for_B_data_here')
 1002     format('put_your_file_for_V_data_here')
 1003     format('put_your_file_for_R_data_here')
 1004     format('put_your_file_for_I_data_here')
 1005     format('put_your_file_for_J_data_here')
 1006     format('put_your_file_for_H_data_here')
 1007     format('put_your_file_for_K_data_here')
 1008     format('put_your_file_for_RV1_data_here')
 1009     format('put_your_file_for_RV2_data_here')
 1010     format(i2,18x,'Number of variables to adjust')
 1011     format('inclination')
 1012     format('mass ratio')
 1013     format('f1 (fill 1)')
 1014     format('f2 (fill 2)')
 1015     format('o1 (omega1)')
 1016     format('o2 (omega2)')
 1017     format('rinner [inner disk radius (same units as fill2)]')
 1018     format('router [outer disk radius','
     %              (in units of Rl volume of star 2)]')
 1019     format('Tdisk (inner disk temperature)')
 1020     format('betarim (half-opening angle of disk rim in degrees)')
 1021     format('T1 (T_eff of star 1)')
 1022     format('T2 (T_eff of star 2)')
 1023     format('xi')
 1024     format('Lx (L_x/L_opt)')
 1025     format('separation (solar radii)')
 1026     format('gamma (gamma velocity in km/sec)')
 1030     format('70.0     1.0     3')
 1031     format('2.0      0.25    3')
 1032     format('0.70     0.1     3')
 1033     format('0.70     0.1     3')
 1034     format('1.0      0.1     3')
 1035     format('1.0      0.1     3')
 1036     format('0.01     0.001   3')
 1037     format('0.75     0.01    3')
 1038     format('20000.0  100.0   3')
 1039     format('4.0      0.5     3')
 1040     format('6500.0   100.0   3')
 1041     format('7500.0   100.0   3')
 1042     format('-0.75    0.05    3')
 1043     format('10.0     1.0     3')
 1044     format('5.0      0.5     3')
 1045     format('-69.0    10.0    3')
c
 2000   format('##############################################')
 2001   format('#    The first 8 entry MUST be a string indicating ',
     $      'the name of the FOLDED')
 2002   format('#    data files.  Put ''none'' if there is no file.')
 2003   format('#')
 3004   format('#    The next entry MUST be a positive integer (Nvar)')
 2004   format('#    The next Nvar entries MUST be strings.  Put the name ',
     $      'of the variable you want ')
 2005   format('#    in the innermost loop first, the name of the ', 
     $      'variable in the second') 
 2006   format('#    loop second, etc.',  
     $      ' If you do not want to loop over all variables, ')
 2007   format('#    put the string ''NOTHING'' in the unwanted positions.')
 2008   format('#    The names of the allowed variables that can be ',
     $      'put into loops are listed')
 2011   format('#    The next Nvar lines contain the starting values',
     $      ' of the variables') 
 2012   format('#    the step size, and the number of loops.')
 2014   format('#')
 2015   format('##############################################')
c
 2020   format('Edit the file ''gridloop.opt'' and restart the program.')
c
          return
          end
c
c   **************************************
c
      INTEGER  FUNCTION  ICNVRT (ASTRING)
c
c    November 12, 1999
c
c    This function was stolen directly out of Peter Stetson's DAOPHOT IIe
c    source code.  The following are his comments:
c
C
C=======================================================================
C
C This little function is supposed to take two ASCII characters and
C express them as an integer in the range 0-(32**2-1) without 
C distinguishing upper and lower case:
C
C AA = Aa = aA = aa = 0, AB = Ab = aB = ab = 1, BA = Ba = bA = ba = 32,
C etc.
C
C Argument
C
C ASTRING is a character string containing two ASCII characters.
C
C=======================================================================
C
      IMPLICIT NONE
      CHARACTER*2 ASTRING
C
C-----------------------------------------------------------------------
C
      ICNVRT=32*MOD(ICHAR(ASTRING(1:1))-1,32)+
     .     MOD(ICHAR(ASTRING(2:2))-1,32)
      RETURN
C
      END!
C
c    End of stolen code.
c
c
c    *************************************************************
c
          subroutine loaddata(Ndatamax,filein,Ndata,xdata,ydata,err)
c
c    November 12, 1999
c
c    This general routine will load a data file.  It assumes that there
c    are three columns:   phase    flux or RV     error
c
          implicit double precision (a-h,o-z)
c
          dimension xdata(Ndatamax),ydata(Ndatamax),err(Ndatamax)
c
          character*40 filein
          character*1 bell
c
          bell=char(7)
          ios=0
c
          Ndata=0
          open(unit=20,file=filein,status='old',err=999,iostat=ios)
c
          do 10 i=1,Ndatamax
            read(20,*,end=15)xdata(i),ydata(i),err(i)
 10       continue
 15       close(20)
c
          Ndata=i-1
c
 999      if(ios.ne.0)then
            write(*,10000)bell,filein
            stop
          endif
c
10000     format(a1,'Error:  data file not found (',a40,')')
c  
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  RVG BUG ALERT   May 8, 2001
c
c  This subroutine has been updated to include the spot parameters
c
c  NEW BUG August 2, 2001
c
c  This routine has been updated to include the period and phase zeropoint.
c
c
c  UPDATE January 15, 2002
c
c  Update the routine to include the albedos of the two stars (alb1, alb2)
c
c  UPDATE November 6, 2002
c
c  Add limb darkening coefficients dwavex and dwavey to the
c  argument list of assignvar and varassign.
c
c  UPDATE August 13, 2004
c
c  Add primmass,primK,primrad,ratrad,frac1,frac2
c
c  UPDATE JULY 29, 2005
c
c  Add ecosw and temprat to the list.
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c  UPDATE October 10, 2008
c
c  Add the density.
c
          subroutine assignvar(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     &       pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,   
     $       dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c   November 12, 1999
c
c   This routine will determine which variables need to be changed
c   based on the string codes in svar(1:Nvarmax)
c
          implicit double precision (a-h,o-z)
c
          character*40 svar(Nvarmax)
          dimension var(Nvarmax)
c
c   RVG BUG ALERT   May 8, 2001
c
c   Dimension the spot parameter arrays.
c       
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
c
c
c  UPDATE November 6, 2002
c
c  dimension the limb darkening variables
c
          dimension dwavex(8,2),dwavey(8,2)
          dimension powercoeff(8,9)

          write(*,*)' '
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
            if(kk.eq.144)then   !beam1, use string e1
              beam1=var(i)
            endif
c
            if(kk.eq.145)then   !beam2, use string e2
              beam2=var(i)
            endif
c
            if(kk.eq.610)then   !Tconj, use string tc
              Tconj=var(i)
            endif
c
            if(kk.eq.100)then   !density, use string de
              density=var(i)
            endif
c
            if(kk.eq.8)then  !bigI, use string ai
              bigI=var(i)
            endif
c
            if(kk.eq.1)then  !bigbeta, use string ab
              bigbeta=var(i)
            endif
c
            if(kk.eq.111)then  !ecosw, use string dphi
              ecosw=var(i)
            endif
c
            if(kk.eq.612)then  !temprat
              temprat=var(i)
            endif
c
            if(kk.eq.492)then  !primmass, use string pm
              primmass=var(i)
              if(primmass.lt.0.0d0)primmass=0.0d0
            endif
c
            if(kk.eq.497)then  !primrad, use string pr
              primrad=var(i)
              if(primrad.lt.0.0d0)primrad=0.0d0
            endif
c
            if(kk.eq.490)then  !primK, use string pk
              primK=var(i)
              if(primK.lt.0.0d0)primK=0.0d0
            endif
c
            if(kk.eq.544)then  !ratrad, use string ra
              ratrad=var(i)
              if(ratrad.lt.0.0d0)ratrad=0.0d0
            endif
c
            if(kk.eq.528)then  !frac1, use string q1
              frac1=var(i)
              if(frac1.lt.0.0d0)frac1=0.0d0
            endif
c
            if(kk.eq.529)then  !frac2, use string q2
              frac2=var(i)
              if(frac2.lt.0.0d0)frac2=0.0d0
            endif
c
c
c   UPDATE January 15, 2002
c
c   Here are the assignments for the albedos
c
            if(kk.eq.368)then  !alb1, use string l1
              alb1=var(i)
              if(alb1.lt.0.0d0)alb1=0.0d0
            endif
c
            if(kk.eq.369)then  !alb2, use string l2
              alb2=var(i)
              if(alb2.lt.0.0d0)alb2=0.0d0
            endif
c
c   NEW BUG August 2, 2001
c
c   Here are the assignments for the period and T0
c
            if(kk.eq.484)then  !period
              period=var(i)
            endif
c
            if(kk.eq.623)then ! T0
              T0=var(i)
            endif
c
c   RVG BUG ALERT   May 8, 2001
c
c   Here is the new block for the assignments of spot parameters.
c
            if(kk.eq.112)then   ! temperature factor spot 1 on disk
              spotdparm(1,1)=var(i)
              if(spotdparm(1,1).le.-10.0d0)spotdparm(1,1)=-1.0d0
            endif
c
            if(kk.eq.116)then   ! temperature factor spot 2 on disk
              spotdparm(2,1)=var(i)
              if(spotdparm(2,1).le.-10.0d0)spotdparm(2,1)=-1.0d0
            endif
c
            if(kk.eq.113)then  ! azimuth spot 1 on disk
              spotdparm(1,2)=dmod(var(i),360.0d0)
            endif   
c
            if(kk.eq.117)then  ! azimuth spot 2 on disk
              spotdparm(2,2)=dmod(var(i),360.0d0)
            endif   
c
            if(kk.eq.114)then    ! cutoff radius for spot 1 on disk
              spotdparm(1,3)=var(i)
              if(spotdparm(1,3).gt.1.0d0)spotdparm(1,3)=1.0d0
              if(spotdparm(1,3).lt.0.0d0)spotdparm(1,3)=0.0d0
            endif
c
            if(kk.eq.118)then    ! cutoff radius for spot 2 on disk
              spotdparm(2,3)=var(i)
              if(spotdparm(2,3).gt.1.0d0)spotdparm(2,3)=1.0d0
              if(spotdparm(2,3).lt.0.0d0)spotdparm(2,3)=0.0d0
            endif
c
            if(kk.eq.115)then    ! angular width of spot 1 on disk
              spotdparm(1,4)=var(i)
              if(spotdparm(1,4).lt.0.0d0)spotdparm(1,4)=0.0d0
            endif
c
            if(kk.eq.119)then    ! angular width of spot 2 on disk
              spotdparm(2,4)=var(i)
              if(spotdparm(2,4).lt.0.0d0)spotdparm(2,4)=0.0d0
            endif
c
c
            if(kk.eq.80)then         ! temperature factor spot 1, star 2
              spot2parm(1,1)=var(i)
            endif
c
            if(kk.eq.84)then         ! temperature factor spot 2, star 2
              spot2parm(2,1)=var(i)
            endif
c
            if(kk.eq.48)then         ! temperature factor spot 1, star 1
              spot1parm(1,1)=var(i)
            endif
c
            if(kk.eq.52)then         ! temperature factor spot 2, star 1
              spot1parm(2,1)=var(i)
            endif
c
            if(kk.eq.81)then         ! latitude spot 1, star 2
              spot2parm(1,2)=dmod(var(i),180.0d0)
            endif
c
            if(kk.eq.85)then         ! latitude spot 2, star 2
              spot2parm(2,2)=dmod(var(i),180.0d0)
            endif
c
            if(kk.eq.49)then         ! latitude spot 1, star 1
              spot1parm(1,2)=dmod(var(i),180.0d0)
            endif
c
            if(kk.eq.53)then         ! latitude spot 2, star 1
              spot1parm(2,2)=dmod(var(i),180.0d0)
            endif
c
            if(kk.eq.82)then         ! longitude spot 1, star 2
              spot2parm(1,3)=dmod(var(i),360.0d0)
            endif
c
            if(kk.eq.86)then         ! longitude spot 2, star 2
              spot2parm(2,3)=dmod(var(i),360.0d0)
            endif
c
            if(kk.eq.50)then         ! longitude spot 1, star 1
              spot1parm(1,3)=dmod(var(i),360.0d0)
            endif
c
            if(kk.eq.54)then         ! longitude spot 2, star 1
              spot1parm(2,3)=dmod(var(i),360.0d0)
            endif
c
            if(kk.eq.83)then         ! radius spot 1, star 2
              spot2parm(1,4)=var(i)
            endif
c
            if(kk.eq.87)then         ! radius spot 2, star 2
              spot2parm(2,4)=var(i)
            endif
c
            if(kk.eq.51)then         ! radius spot 1, star 1
              spot1parm(1,4)=var(i)
            endif
c
            if(kk.eq.55)then         ! radius spot 2, star 1
              spot1parm(2,4)=var(i)
            endif

            if(kk.eq.498)then
              pshift=var(i)
              if(pshift.gt.1.0d0)pshift=1.0d0
              if(pshift.lt.-1.0d0)pshift=-1.0d0
            endif
c
            if(kk.eq.269)then
              finc=var(i)
              if(finc.gt.90.0)finc=90.0
              if(finc.lt.0.0)finc=0.0
            endif
c
            if(kk.eq.384)then
              Q=var(i)
              if(Q.le.0.0)Q=0.001
            endif
c
            if(kk.eq.130)then
              ecc=var(i)
              if(ecc.le.0.0d0)ecc=0.000d0
              if(ecc.ge.1.0d0)ecc=0.9999d0
            endif
c
            if(kk.eq.17)then
              argper=var(i)
c              if(argper.le.0.0d0)argper=0.000d0
c              if(argper.gt.360.0d0)argper=360.0d0
            endif
c
            if(kk.eq.176)then
              fill1=var(i)
              if(fill1.gt.1.0)fill1=1.0
            endif
c
            if(kk.eq.552)then
              rinner=var(i)
              if((rinner.lt.fill2).and.(Teff2.gt.0.0))rinner=fill2
            endif
c
            if(kk.eq.177)then
              fill2=var(i)
              if(fill2.gt.1.0)fill2=1.0
              if(teff2.gt.0.0)rinner=fill2
              if(fill2.lt.0.0)then
                fill2=0.000001
                rinner=fill2
              endif
            endif
c
            if(kk.eq.464)then
              omega1=var(i)
            endif
c
            if(kk.eq.465)then
              omega2=var(i)
            endif
c
            if(kk.eq.558)then
              router=var(i)
              if(router.gt.1.0)router=1.0
            endif
c
            if(kk.eq.611)then
              Tdisk=var(i)
              if(Tdisk.lt.100.0)Tdisk=100.
            endif
c
            if(kk.eq.36)then
              betarim=var(i)
              if(betarim.lt.0.0)betarim=0.0
            endif
c
            if(kk.eq.624)then
              Teff1=var(i)
              if(Teff1.lt.100.0)Teff1=100.
            endif
c
            if(kk.eq.625)then
              Teff2=var(i)
              if(Teff2.lt.100.0)Teff2=100.
            endif
c
            if(kk.eq.744)then
              xi=var(i)
            endif
c
            if(kk.eq.375)then
              rLx=var(i)
              if(rLx.lt.0.0)rLx=0.0
            endif
c
            if(kk.eq.580)then
              separ=var(i)
              if(separ.lt.0.01)separ=0.01
            endif
c
            if(kk.eq.192)then
              gamma=var(i)
            endif
c
            if(kk.eq.626)then
              t3=var(i)
              if(t3.lt.0.01)t3=0.01
            endif
c
            if(kk.eq.210)then
              g3=var(i)
              if(g3.lt.0.01)t3=0.01
            endif
c
            if(kk.eq.576)then
              SA3=var(i)
c
c   UPDATE January 16, 2001
c
c   comment out this if-then statement
c
c              if(SA3.lt.0.01)t3=0.01
            endif
c
c   UPDATE November 6, 2002
c
c   Here are the assignments for the limb darkening parameters.
c   use strings x1, x2, to x8 for the x-coefficient for star 1
c   and y1, y2, to y8 for the y-coefficient for star 1.
c
c   Use z1, z2, to z8 for the x-coefficient for star 2 and
c   use w1, w2, to w8 for the y-coefficient for star 2.
c
            if(kk.eq.752)dwavex(1,1)=var(i)
            if(kk.eq.753)dwavex(2,1)=var(i)
            if(kk.eq.754)dwavex(3,1)=var(i)
            if(kk.eq.755)dwavex(4,1)=var(i)
            if(kk.eq.756)dwavex(5,1)=var(i)
            if(kk.eq.757)dwavex(6,1)=var(i)
            if(kk.eq.758)dwavex(7,1)=var(i)
            if(kk.eq.759)dwavex(8,1)=var(i)
c
            if(kk.eq.784)dwavey(1,1)=var(i)
            if(kk.eq.785)dwavey(2,1)=var(i)
            if(kk.eq.786)dwavey(3,1)=var(i)
            if(kk.eq.787)dwavey(4,1)=var(i)
            if(kk.eq.788)dwavey(5,1)=var(i)
            if(kk.eq.789)dwavey(6,1)=var(i)
            if(kk.eq.790)dwavey(7,1)=var(i)
            if(kk.eq.791)dwavey(8,1)=var(i)
c
            if(kk.eq.816)dwavex(1,2)=var(i)
            if(kk.eq.817)dwavex(2,2)=var(i)
            if(kk.eq.818)dwavex(3,2)=var(i)
            if(kk.eq.819)dwavex(4,2)=var(i)
            if(kk.eq.820)dwavex(5,2)=var(i)
            if(kk.eq.821)dwavex(6,2)=var(i)
            if(kk.eq.822)dwavex(7,2)=var(i)
            if(kk.eq.823)dwavex(8,2)=var(i)
c
            if(kk.eq.720)dwavey(1,2)=var(i)
            if(kk.eq.721)dwavey(2,2)=var(i)
            if(kk.eq.722)dwavey(3,2)=var(i)
            if(kk.eq.723)dwavey(4,2)=var(i)
            if(kk.eq.724)dwavey(5,2)=var(i)
            if(kk.eq.725)dwavey(6,2)=var(i)
            if(kk.eq.726)dwavey(7,2)=var(i)
            if(kk.eq.727)dwavey(8,2)=var(i)
c
 10       continue
c
 999      format(7x,i4,2x,a40,3x,f15.8)
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&
c
          subroutine checklcfit(islc,Nmodel,xmodel,
     %            yinput,Ndata,xdata,ydata,err,chisq,zero,
     #            ifixgamma,itime)
c
c   November 12, 1999
c
c   This routine will convert the linear model to magnitudes and adjust
c   the zeropoint to give the best chi^2.  The best chi^2 value is
c   returned, as well as the zeropoint.  Set islc=1 for light curves
c   (linear models are converted to magnitudes) and islc=0 for RV curves.
c
c
c   January 24, 2001
c
c   Added the parameter ifixgamma.  If this is 1 or more and we are
c   fitting a velocity curve then the zero point is forced to be the input
c   gamma.
c
c

          implicit double precision (a-h,o-z)
c
          dimension xmodel(Nmodel),yinput(Nmodel),xdata(Ndata),ydata(Ndata),
     %      err(Ndata),ymodel(900000)
          dimension xinter(900000),yinter(900000),y2(900000),ydummy(900000)
          dimension xpad(900000),ypad(900000)
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   First, convert to magnitudes.
c
c
           if(itime.le.1)then
             call addpad(Nmodel,xmodel,yinput,xpad,ypad)
             Mmodel=Nmodel*3
           else
             do jj=1,Nmodel
               xpad(jj)=xmodel(jj)
               ypad(jj)=yinput(jj)
             enddo
             Mmodel=Nmodel
           endif
 
c            open(unit=88,file='testmodel',status='unknown')
            do 10 i=1,Mmodel
              if(islc.gt.0)then
c                write(88,*)xpad(i),ypad(i)
c
c  NEW_BUG  July 5, 2001
c
c  Check to see that the argument of log10 is positive.
c
                if(ypad(i).gt.0.0d0)then
                  ymodel(i)=-2.5d0*dlog10(ypad(i))
                else
                  ymodel(i)=-99.9d0
                endif
              else
                ymodel(i)=ypad(i)
              endif
 10         continue
c            close(88)
c
c   Find the maximum and minimum y-values of the data
c
          ymin=1000.
          ymax=-1000.
          do 20 i=1,Ndata
            if(ydata(i).lt.ymin)ymin=ydata(i)
            if(ydata(i).gt.ymax)ymax=ydata(i)
 20       continue
c

c
c   Interpolate the model so that we have y-values at all observed phases.
c
          call spline(xpad,ymodel,Mmodel,0.0d0,0.0d0,y2)
c
c          open(unit=88,file='testdata',status='unknown')
          do 30 i=1,Ndata
c            write(88,*)xdata(i),ydata(i),err(i)
            call splint(xpad,ymodel,y2,Mmodel,xdata(i),qqq)
            xinter(i)=xdata(i)
            yinter(i)=qqq
 30       continue
c          close(88)
c
c   Now find the optimal zero point that will give the lowest chi^2.
c
          call getmean(Ndata,xdata,ydata,dataave)
          call getmean(Ndata,xinter,yinter,rmodelave)

          savezero=zero
          zero=(dataave-rmodelave)
          step=abs(rmodelave-(dataave))/1000.0d0

          chi1=0.0d0
          chi2=0.0d0
          chi3=0.0d0
          small=1.0d35
c
          do 750 i=1,40
            zero1=zero
            call offset(Ndata,xinter,yinter,ydummy,zero1)
            call getchi(Ndata,ydata,err,ydummy,chi1)
            if(chi1.lt.small)then
              small=chi1
              zerosmall=zero1
            endif
            fn=0.0d0
            zero2=zero+step
            call offset(Ndata,xinter,yinter,ydummy,zero2)
            call getchi(Ndata,ydata,err,ydummy,chi2)
            if(chi2.lt.small)then
              small=chi2
              zerosmall=zero2
            endif
            diff=chi1-chi2
            if(diff.gt.0.0d0)go to 5061
            step=-step
c  
            csave=chi1
            chi1=chi2
            chi2=csave
            zsave=zero1
            zero1=zero2
            zero2=zsave
c
 5061       fn=fn+1.0d0
c
            zero3=zero2+step
            call offset(Ndata,xinter,yinter,ydummy,zero3)
            call getchi(Ndata,ydata,err,ydummy,chi3)
            if(chi3.lt.small)then
              small=chi3
              zerosmall=zero3
            endif
            diff23=chi3-chi2
            if(diff23.lt.0.0d0)then
              chi1=chi2
              chi2=chi3
              zero1=zero2
              zero2=zero3
              zero=zero3
              go to 5061
            endif          
c
c         find the minimum of parabola defined by the last three points
c
 5081       if((chi3-chi2).eq.0.0d0)go to 999
            step=step*(1.0d0/(1.0d0+(chi1-chi2)/(chi3-chi2))+0.5d0)
            zero=zero2-step
            step=step*fn/3.0d0
            call offset(Ndata,xinter,yinter,ydummy,zero)
            call getchi(Ndata,ydata,err,ydummy,chi4)
            if(chi4.lt.small)then
              small=chi4
              zerosmall=zero
            endif
c
 750      continue  ! loop over grid searches
c
 751      continue                  ! come here when delta chi is small
c
 999      zero=zerosmall
c
          if((islc.le.0).and.(ifixgamma.ge.1))zero=savezero
          call offset(Ndata,xinter,yinter,ydummy,zero)
          call getchi(Ndata,ydata,err,ydummy,chi3)
c
c          write(*,*)zero
          chisq=chi3
c
          return
c
          end
c
c  -----------------------------------------------------
c
          subroutine getchi(N,y,err,ydum,chisq)
c
          implicit double precision (a-h,o-z)
c
          dimension y(N),ydum(N),err(N)
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 

          chisq=0.0d0
c
c   If the variable rmed > 1, then define chisq to be the
c   absolute deviation
c
          do 10 i=1,N
            if(rmed.ge.1.0d0)then
              chisq=chisq+dabs((y(i)-ydum(i))/err(i))
            else
              chisq=chisq+(y(i)-ydum(i))*(y(i)-ydum(i))/(err(i)*err(i))
            endif
 10       continue
          return
          end
c
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine offset(N,x,y,yoff,off)
c
c   Will offset the y values by offset and return new array yoff.
c
          implicit double precision (a-h,o-z)
c
          dimension x(N),y(N),yoff(N)
c
          do 10 i=1,N
            yoff(i)=y(i)+off
 10       continue
c
          return
          end
c
c  *********************************************************
c
          subroutine getmean(N,x,y,average)
c
          implicit double precision (a-h,o-z)
c
          dimension x(N),y(N)
c
          summ=0.0
          average=0.0
c
          do 10 i=1,N
            summ=summ+y(i)
c            write(*,*)x(i),y(i),summ
 10       continue
c 
          average=summ/float(N)
c
c          write(*,*)'***'
c          write(*,*)average,N
 
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &            chisqH,chisqK,chisqRV1,chisqRV2)
c
          implicit double precision (a-h,o-z)
c
          chisqU=0.d0
          chisqB=0.d0
          chisqV=0.d0
          chisqR=0.d0
          chisqI=0.d0
          chisqJ=0.d0
          chisqH=0.d0
          chisqK=0.d0
          chisqRV1=0.d0
          chisqRV2=0.d0
c
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %       icnRV1,icnRV2,chiU,chiB,chiV,chiR,chiI,chiJ,chiH,chiK,
     &       chiRV1,chiRV2)
c
c   November 15, 1999
c
c   This routine will write chi^2 values to the screen in a compact way.
c
          implicit double precision (a-h,o-z)
c
          character*17 section(10)
          character*80 line
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   UPDATE May 27, 2002.
c
c   If the variable rmed > 1, then we are doing median fitting.
c   In that case, use different format statements.
c
          iline=0
          iwrite=0
c
          do 10 i=1,10
            section(i)='1234567890123456'
            section(i)='                '
 10       continue
          line='  -->'
c
          if(icnU.ne.430)then
            iline=iline+1
            if(rmed.ge.1.0d0)then
              write(section(iline),2000)chiU
            else
              write(section(iline),100)chiU
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnB.ne.430)then
            iline=iline+1
            if(rmed.ge.1.0d0)then
              write(section(iline),201)chiB
            else
              write(section(iline),101)chiB
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnV.ne.430)then
            iline=iline+1
            if(rmed.ge.1.0d0)then
              write(section(iline),202)chiV
            else
              write(section(iline),102)chiV
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnR.ne.430)then
            iline=iline+1
            if(rmed.ge.1.0d0)then
              write(section(iline),203)chiR
            else
              write(section(iline),103)chiR
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
            if(iline.eq.4)then
              write(*,200)line
              iwrite=100
              iline=0
              line='  -->'
            endif
          endif
          if(icnI.ne.430)then
            iline=iline+1
            if(rmed.ge.1.0d0)then
              write(section(iline),204)chiI
            else
              write(section(iline),104)chiI
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
            if(iline.eq.4)then
              write(*,200)line
              iwrite=100
              iline=0
              line='  -->'
            endif
          endif
          if(icnJ.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              write(section(iline),205)chiJ
            else
              write(section(iline),105)chiJ
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
            if(iline.eq.4)then
              write(*,200)line
              iwrite=100
              line='  -->'
              iline=0
            endif
          endif
          if(icnH.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              write(section(iline),206)chiH
            else
              write(section(iline),106)chiH
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
            if(iline.eq.4)then
              write(*,200)line
              iwrite=100
              iline=0
              line='  -->'
            endif
          endif
          if(icnK.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              write(section(iline),207)chiK
            else
              write(section(iline),107)chiK
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
            if(iline.eq.4)then
              write(*,200)line
              iwrite=100
              iline=0
              line='  -->'
            endif
          endif
          if(icnRV1.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              write(section(iline),208)chiRV1
            else
              write(section(iline),108)chiRV1
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
            if(iline.eq.4)then
              write(*,200)line
              iwrite=100
              iline=0
              line='  -->'
            endif
          endif
          if(icnRV2.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              write(section(iline),209)chiRV2
            else
              write(section(iline),109)chiRV2
            endif
            K=lnblnk(line)
            line=line(1:K)//section(iline)
            if(iline.eq.4)then
              write(*,200)line
              iwrite=100
              iline=0
              line='  -->'
            endif
          endif
c
          if(iline.gt.0.and.iwrite.eq.0)write(*,200)line
c
c   UPDATE May 27, 2002
c
c   Define the new format statements below (200-209).
c
 100      format(' chiU=',f11.4)
 101      format(' chiB=',f11.4)
 102      format(' chiV=',f11.4)
 103      format(' chiR=',f11.4)
 104      format(' chiI=',f11.4)
 105      format(' chiJ=',f11.4)
 106      format(' chiH=',f11.4)
 107      format(' chiK=',f11.4)
 108      format(' chiRV1=',f9.2)
 109      format(' chiRV2=',f9.2)
 2000     format(' medU=',f11.4)
 201      format(' medB=',f11.4)
 202      format(' medV=',f11.4)
 203      format(' medR=',f11.4)
 204      format(' medI=',f11.4)
 205      format(' medJ=',f11.4)
 206      format(' medH=',f11.4)
 207      format(' medK=',f11.4)
 208      format(' medRV1=',f9.2)
 209      format(' medRV2=',f9.2)
 200      format(a80)
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  RVG BUG ALERT  May 9, 2001
c
c  Add the spot parameters to the argument list.
c
c  NEW BUG ALERT August 2, 2001
c
c  Replace 'sw4' with T0
c
c  UPDATE August 10, 2004
c
c  Add 8 real and 4 integer variables to the argument list.
c
c  May 8, 2006
c
c  Add isw21-isw24, sw21-sw24, powercoeff to list.
c
c  UPDATE November 6, 2008
c
c  Add sw25-sw34 and isw25-isw34 below
c
          subroutine writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
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
     &       temprat,idark1,idark2,isw12,isw13,
     #       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,Tconj,beam1,beam2,
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
c
c   November 15, 1999
c
c   This routine will write a file called 'gridELC.inp' which will be
c   of the same format as ELC.inp.  The parameters written will be those
c   found by gridELC.
c
          implicit double precision (a-h,o-z)
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
          open(unit=1,file='gridELC.inp',status='unknown')
c
c   Set the icn? control numbers back to zeros and 1s.
c
          if(icnU.eq.430)then
            icnU=0
          else
            icnU=1
          endif
          if(icnB.eq.430)then
            icnB=0
          else
            icnB=1
          endif
          if(icnV.eq.430)then
            icnV=0
          else
            icnV=1
          endif
          if(icnR.eq.430)then
            icnR=0
          else
            icnR=1
          endif
          if(icnI.eq.430)then
            icnI=0
          else
            icnI=1
          endif
          if(icnJ.eq.430)then
            icnJ=0
          else
            icnJ=1
          endif
          if(icnH.eq.430)then
            icnH=0
          else
            icnH=1
          endif
          if(icnK.eq.430)then
            icnK=0
          else
            icnK=1
          endif
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
          if(fm.gt.1.0d-4)then
            write(1,1024)fm
          else
            write(1,9024)fm
          endif
 9024     format(1pe16.9,4x,'fm')
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
c   RVG BUG ALERT  May 9, 2001
c
c   Add the spot parameters here.
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
          write(1,2028)ratrad
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
              write(1,85001)powercoeff(kk,1),powercoeff(kk,2),
     $       powercoeff(kk,3),powercoeff(kk,4),powercoeff(kk,5),
     $       powercoeff(kk,6),powercoeff(kk,7),powercoeff(kk,8),
     #       powercoeff(kk,9)
85000     continue
c
85001      format(9(f7.4,1x))
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
c  UPDATE November 6, 2008
c
c  write the new variables sw25-sw34 and isw25-isw34 here
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
8029      format(f10.4,16x,'bin size for light curves (minutes)')
8030      format(f10.4,16x,'bin size for RV curves (minutes)')
8031      format(f4.2,16x,'sw31 (currently inactive)')
8032      format(f15.8,5x,'Tconj')
8033      format(f4.2,16x,'beam1 (Doppler boosting factor, star 1)')
8034      format(f4.2,16x,'beam2 (Doppler boosting factor, star 2)')

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
 2030     format(f14.12,6x,'frac1 (fractional radius star 1: R_1/a)')
 2031     format(f14.12,6x,'frac2 (fractional radius star 2: R_2/a)')
 2032     format(f12.9,8x,'ecosw (phase difference between eclipses)')
 2033     format(f10.7,10x,'temprat (T_2/T_1)')
c
 2040     format(i1,19x,'idark1')
 2041     format(i1,19x,'idark2')
 2042     format(i6,14x,'Npoly (0 for numerical integration)')
 3043     format(i1,19x,'ifasttrans (>0 for fast transit mode)')

          close(1)
c
 100      format(a1,'I can''t find the file ''ELC.inp''!  I''m making',
     $     ' one up and setting default values')
c
 1000     format(i4,16x,'Nalph1')
 1001     format(i3,17x,'Nbet1')
 1002     format(i4,16x,'Nalph2')
 1003     format(i3,17x,'Nbet2')
 1004     format(f11.9,9x,'fill1')
 1005     format(f11.9,9x,'fill2')
 1006     format(f10.6,10x,'omega1')
 1007     format(f10.6,10x,'omega2')
 1008     format(f11.6,9x,'dphase')
 1009     format(f15.10,5x,'Q')
 1010     format(f8.5,12x,'finc')
 1011     format(f9.2,11x,'Teff1')
 1012     format(f9.2,11x,'Teff2')
 1013     format(f7.5,13x,'Tgrav1')
 1014     format(f7.5,13x,'Tgrav2')
 1015     format(f8.5,12x,'betarim')
 1016     format(f8.6,12x,'rinner')
 1017     format(f8.6,12x,'router')
 1018     format(f7.1,13x,'tdisk')
 2018     format(f7.4,13x,'xi')
 1019     format(i3,17x,'Ntheta')
 1020     format(i3,17x,'Nradius')
 1021     format(f10.5,10x,'Lx/Lopt')
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
 5003     format(f12.6,8x,'density in g/cc')
 5004     format(f12.6,8x,'onephase')
 5005     format(f12.6,8x,'usepot1')
 5006     format(f12.6,8x,'usepot2')
 5007     format(f16.9,4x,'T0')
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
c   Change the format statement below:
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
c  will be used in the Monte Carlo routine to determine
c  fractionally eclipsed pixels.
c
 6013     format(i6,14x,'MonteCarlo (0 for interpolation, >10 for ',
     $     'Monte Carlo integration)')
 6014     format(i1,19x,'ielite')
c
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
c
c  NEW BUG ALERT  August 2, 2001
c
c  Reset the values of icnU, etc. in case they were 430
c
c
c   Set the icn? control numbers back to zeros and 1s.
c
          if(icnU.eq.0)then
            icnU=430
          endif
          if(icnB.eq.0)then
            icnB=430
          endif
          if(icnV.eq.0)then
            icnV=430
          endif
          if(icnR.eq.0)then
            icnR=430
          endif
          if(icnI.eq.0)then
            icnI=430
          endif
          if(icnJ.eq.0)then
            icnJ=430
          endif
          if(icnH.eq.0)then
            icnH=430
          endif
          if(icnK.eq.0)then
            icnK=430
          endif

c
          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine recordloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,Nvmax,Nvar,
     $      svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)
c
c    November 15, 1999
c
c    This routine will write an input file to set up the optimizer (either
c    the 'grid search' or the Levenburg-Marquardt routine), based on
c    the final parameters derived by the optimizer.
c
c     UPDATE September 21, 2008
c
c     make the dimension of sobv, eobv, and obv 11
c
c
          implicit double precision (a-h,o-z)
c
          dimension obv(11),eobv(11)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax)
          character*40 Udatafile,svar(Nvmax),Hdatafile,Kdatafile,
     %        RV1file,RV2file,sobv(11)
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*1 bell
c
          bell=char(7)
          iso=0
          open(unit=11,file='gridELC.opt',status='unknown',iostat=ios)
c
          write(11,100)Udatafile
          write(11,100)Bdatafile
          write(11,100)Vdatafile
          write(11,100)Rdatafile
          write(11,100)Idatafile
          write(11,100)Jdatafile
          write(11,100)Hdatafile
          write(11,100)Kdatafile
          write(11,100)RV1file
          write(11,100)RV2file
c
          write(11,500)Nvar
c
          do 10 i=1,Nvar
            write(11,100)svar(i)
            if((svar(i)(1:2).eq.'f1').or.(svar(i)(1:2).eq.'F1'))then
              if(var(i).gt.1.0d0)var(i)=1.0d0
            endif
            if((svar(i)(1:2).eq.'f2').or.(svar(i)(1:2).eq.'F2'))then
              if(var(i).gt.1.0d0)var(i)=1.0d0
            endif
 10       continue
c
          do 20 i=1,Nvar
            write(11,501)var(i),vstep(i),Nstep(i)
 20       continue
c
          if(Nobv.lt.1)go to 30
          write(11,500)Nobv
          do 25 i=1,Nobv
            write(11,100)sobv(i)
 25       continue
c
          do 26 i=1,Nobv
            write(11,502)obv(i),eobv(i)
 26       continue
c
 30       close(11)
c
 100      format(a40)
 200      format(a1,'Error:  too many variables in ''gridloop.opt''')
 500      format(i2)
 501      format(f17.9,3x,f17.9,3x,i7)
 502      format(f16.9,3x,f15.10)
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  NEW BUG August 2, 2001
c
c  Add the period and T0 to the argument list.
c
c  UPDATE January 15, 2002
c
c  Add alb1 and alb2 to the argument list
c
c  UPDATE November 6, 2002
c
c  Add the limb darkening coefficients dwavex and dwavey to the
c  argument list of writevar.  Dimension them as dwavex(8,2)
c  and dwavey(8,2).
c
          subroutine writevar(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %       t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $       period,T0,alb1,alb2,dwavex,dwavey,
     &       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     &       bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c  UPDATE October 10, 2008
c
c  Add the density to the list above
c
c   November 15, 1999
c
c   This routine will determine which variables need to be printed on the
c   screen based on the string codes in svar(1:Nvarmax).  The output
c   printed on the screen will be compact.
c
          implicit double precision (a-h,o-z)
c
          character*80 line
          character*40 svar(Nvarmax)
          character*16 string16
          character*10 string10
          character*11 string11
          character*12 string12
          character*13 string13
          character*14 string14
          character*15 string15
c
c   UPDATE August 2, 2004
c
c   Add string 17 for period and T0.
c
          character*17 string17
          character*18 string18
          character*19 string19
          character*22 string22
          character*23 string23

          dimension var(Nvarmax)
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension dwavex(8,2),dwavey(8,2)
c
          line='**'
          iline=0
          iwrite=0
          ilength=2
          iargper=0              !flag for ecosw being assigned.
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
            if(kk.eq.144)then          !beam1, tag e1
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,781)beam1
              K=lnblnk(line)
              line=line(1:K)//string14
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.145)then
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,782)beam2
              K=lnblnk(line)
              line=line(1:K)//string14
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.492)then
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,777)primmass
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.497)then
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,778)primrad
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.612)then    ! temprat
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,1778)temprat
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.490)then
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,779)primK
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.544)then
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,776)ratrad
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.528)then
              iline=iline+1
              iwrite=0
              ilength=ilength+16
              write(string16,775)frac1
              K=lnblnk(line)
              line=line(1:K)//string16
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.100)then
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,780)density
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.529)then
              iline=iline+1
              iwrite=0
              ilength=ilength+16
              write(string16,774)frac2
              K=lnblnk(line)
              line=line(1:K)//string16
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
 775        format(' frac1=',f9.6)
 774        format(' frac2=',f9.6)
 776        format(' ratrad=',f11.6)
 777        format(' primmass=',f9.5)
 778        format(' primrad=',f10.5)
 779        format(' primK=',f12.8)
780         format(' density=',f10.7)
781         format(' beam1=',f7.3)
782         format(' beam2=',f7.3)
c
1778        format(' temprat=',f10.7)
 
c   UPDATE August 13, 2001
c
c   Do a global replace and change if(ilength.gt.60) to if(ilength.gt.51)
c
c
c   UPDATE January 15, 2002
c
c   Here are the blocks for alb1 and alb2
c
            if(kk.eq.368)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2202)alb1
              K=lnblnk(line)
              line=line(1:K)//string12
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.369)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2203)alb2
              K=lnblnk(line)
              line=line(1:K)//string12
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.111)then   !ecosw, use string dphi
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3202)ecosw
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
              iargper=100
            endif
c

            if(kk.eq.610)then           !Tconj string
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,7701)Tconj
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.623)then           !T0 string
              iline=iline+1
              iwrite=0
              ilength=ilength+16
              write(string16,701)T0
              K=lnblnk(line)
              line=line(1:K)//string16
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.484)then      !period string
              iline=iline+1
              iwrite=0
              ilength=ilength+17
              write(string17,700)period
              K=lnblnk(line)
              line=line(1:K)//string17
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.269)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,200)finc
              K=lnblnk(line)
              line=line(1:K)//string10
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   UPDATE May 8, 2006
c
c   Add bigI and bigbeta
c
            if(kk.eq.8)then          !bigI
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,88200)bigI
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1)then          !bigbeta
              iline=iline+1
              iwrite=0
              ilength=ilength+17
              write(string17,88201)bigbeta
              K=lnblnk(line)
              line=line(1:K)//string17
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif

c
c
c   UPDATE JULY @7, 2004
c
c   Change Q to string length 14 (f11.8).
c
            if(kk.eq.384)then
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,201)Q
              K=lnblnk(line)
              line=line(1:K)//string14
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   UPDATE AUGUST 2, 2004
c
c   Make fill1 and fill2 string length 14
c
            if(kk.eq.176)then
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,202)fill1
              K=lnblnk(line)
              line=line(1:K)//string14
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.177)then
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,203)fill2
              K=lnblnk(line)
              line=line(1:K)//string14
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.552)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,204)rinner
              K=lnblnk(line)
              line=line(1:K)//string12
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.464)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,205)omega1
              K=lnblnk(line)
              line=line(1:K)//string12
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.465)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,206)omega2
              K=lnblnk(line)
              line=line(1:K)//string12
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.558)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,207)router
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.626)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,217)t3  
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.210)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,218)g3  
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.36)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,211)betarim
              K=lnblnk(line)
              line=line(1:K)//string12
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.624)then
              iline=iline+1
              iwrite=0
              ilength=ilength+16
              write(string16,209)Teff1
              K=lnblnk(line)
              line=line(1:K)//string16
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.625)then
              iline=iline+1
              iwrite=0
              ilength=ilength+16
              write(string16,210)Teff2
              K=lnblnk(line)
              line=line(1:K)//string16
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.744)then
              iline=iline+1
              iwrite=0
              ilength=ilength+11
              write(string11,212)xi
              K=lnblnk(line)
              line=line(1:K)//string11
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.611)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,208)tdisk 
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.375)then
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,213)rLx
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.580)then
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,214)separ
              K=lnblnk(line)
              line=line(1:K)//string14
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.576)then
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,219)SA3
              K=lnblnk(line)
              line=line(1:K)//string14
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.498)then
              iline=iline+1
              iwrite=0
              ilength=ilength+15
              write(string15,220)pshift
              K=lnblnk(line)
              line=line(1:K)//string15
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.130)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,221)ecc
              K=lnblnk(line)
              line=line(1:K)//string10
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.17)then
              iline=iline+1
              iwrite=0
              ilength=ilength+15
              write(string15,222)argper
              K=lnblnk(line)
              line=line(1:K)//string15
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.48)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,401)spot1parm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.49)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,402)spot1parm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.50)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,403)spot1parm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.51)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,404)spot1parm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.52)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,405)spot1parm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.53)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,406)spot1parm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.54)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,407)spot1parm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.55)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,408)spot1parm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.80)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,501)spot2parm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.81)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,502)spot2parm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.82)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,503)spot2parm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.83)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,504)spot2parm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.84)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,505)spot2parm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.85)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,506)spot2parm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.86)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,507)spot2parm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.87)then
              iline=iline+1
              iwrite=0
              ilength=ilength+23
              write(string23,508)spot2parm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string23
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.112)then
              iline=iline+1
              iwrite=0
              ilength=ilength+22
              write(string22,601)spotdparm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string22
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.113)then
              iline=iline+1
              iwrite=0
              ilength=ilength+22
              write(string22,602)spotdparm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string22
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.114)then
              iline=iline+1
              iwrite=0
              ilength=ilength+22
              write(string22,603)spotdparm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string22
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.115)then
              iline=iline+1
              iwrite=0
              ilength=ilength+22
              write(string22,604)spotdparm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string22
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.116)then
              iline=iline+1
              iwrite=0
              ilength=ilength+22
              write(string22,605)spotdparm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string22
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.117)then
              iline=iline+1
              iwrite=0
              ilength=ilength+22
              write(string22,606)spotdparm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string22
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.118)then
              iline=iline+1
              iwrite=0
              ilength=ilength+22
              write(string22,607)spotdparm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string22
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.119)then
              iline=iline+1
              iwrite=0
              ilength=ilength+22
              write(string22,608)spotdparm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string22
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   UPDATE November 6, 2002
c
c   Add the limb darkening coefficients here.
c
c   x-coefficient, star 1
c
            if(kk.eq.752)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,801)dwavex(1,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.753)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,802)dwavex(2,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.754)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,803)dwavex(3,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.755)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,804)dwavex(4,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.756)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,805)dwavex(5,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.757)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,806)dwavex(6,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.758)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,807)dwavex(7,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.759)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,808)dwavex(8,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   x-coefficient, star 2
c
            if(kk.eq.816)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,2801)dwavex(1,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.817)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,2802)dwavex(2,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.818)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,2803)dwavex(3,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.819)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,2804)dwavex(4,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.820)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,2805)dwavex(5,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.821)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,2806)dwavex(6,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.822)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,2807)dwavex(7,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.823)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,2808)dwavex(8,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   y-coefficient, star 1
c
            if(kk.eq.784)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,1801)dwavey(1,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.785)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,1802)dwavey(2,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.786)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,1803)dwavey(3,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.787)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,1804)dwavey(4,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.788)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,1805)dwavey(5,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.789)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,1806)dwavey(6,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.790)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,1807)dwavey(7,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.791)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,1808)dwavey(8,1)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   y-coefficient, star 2
c
            if(kk.eq.720)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3801)dwavey(1,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.721)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3802)dwavey(2,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.722)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3803)dwavey(3,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.723)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3804)dwavey(4,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.724)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3805)dwavey(5,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.725)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3806)dwavey(6,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.726)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3807)dwavey(7,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.727)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,3808)dwavey(8,2)
              K=lnblnk(line)
              line=line(1:K)//string13
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
 10       continue
c
c
c   Update July 29, 2005
c
c   If the variable ecosw was assigned, then also print out the argument
c   of periastron.
c
            if(iargper.gt.0)then
              call getom(ecc,ecosw,argper)
              iline=iline+1
              iwrite=0
              ilength=ilength+15
              write(string15,222)argper
              K=lnblnk(line)
              line=line(1:K)//string15
              if(ilength.gt.51)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif


c
c UPDATE June 10, 2003
c
c Add two more digits to gamma, use string 15
c
          if(isvel1.gt.0)then
            iline=iline+1
            iwrite=0
            ilength=ilength+15
            write(string15,215)gamma1
            K=lnblnk(line)
            line=line(1:K)//string15
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=2
                line=' *'
            endif
          endif
c
          if(isvel2.gt.0)then
            iline=iline+1
            iwrite=0
            ilength=ilength+15
            write(string15,216)gamma2
            K=lnblnk(line)
            line=line(1:K)//string15
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=2
                line=' *'
            endif
          endif
c
          if(ilength.gt.8.and.iwrite.eq.0)write(*,500)line
c
200       format(' i=',f7.4)
88200     format(' bigI=',f7.4)
88201     format(' bigbeta=',f8.4)
c
c   UPDATE JULY 27, 2004
c
c   Change this format statement to allow for more decimal
c   places to be displayed for Q.
c
 201      format(' Q=',f11.8)
 202      format(' fi1=',f9.7)
 2202     format(' al1=',f7.5)
 203      format(' fi2=',f9.7)
 2203     format(' al2=',f7.5)
 204      format(' rin=',f7.5)
 205      format(' om1=',f7.4)
 206      format(' om2=',f7.4)
 207      format(' rout=',f7.5)
 208      format(' Td=',f9.2)
 209      format(' Teff1=',f9.2)
 210      format(' Teff2=',f9.2)
 211      format(' beta=',f6.3)
 212      format(' xi=',f7.4)
 213      format(' Lx/Lopt=',1pe10.4)
 214      format(' sep=',f9.4)
c
3202     format(' dphi=',f7.5)
 
c   UPDATE June 10, 2003
c
c   Add two more digits to gamma1 and gamma2 (f9.4)
c
 215      format(' gam1=',f9.4)
 216      format(' gam2=',f9.4)
 217      format(' t3=',f9.2)
 218      format(' g3=',f9.2)
 219      format(' SA3=',f9.4)
 220      format(' pshift=',f7.4)
 221      format(' e=',f7.5)
 222      format(' argper=',f7.3)
c
c   RVG BUG ALERT   May 8, 2001
c
c   Add these statements for spots
c
 401      format(' TF_spot1_star1=',f7.4)
 402      format(' lat_spot1_star1=',f6.2)
 403      format(' lon_spot1_star1=',f6.2)
 404      format(' rad_spot1_star1=',f6.3)
 405      format(' TF_spot2_star1=',f7.4)
 406      format(' lat_spot2_star1=',f6.2)
 407      format(' lon_spot2_star1=',f6.2)
 408      format(' rad_spot2_star1=',f6.3)
 501      format(' TF_spot1_star2=',f7.4)
 502      format(' lat_spot1_star2=',f6.2)
 503      format(' lon_spot1_star2=',f6.2)
 504      format(' rad_spot1_star2=',f6.3)
 505      format(' TF_spot2_star2=',f7.4)
 506      format(' lat_spot2_star2=',f6.2)
 507      format(' lon_spot2_star2=',f6.2)
 508      format(' rad_spot2_star2=',f6.3)
 601      format(' TF_spot1_disk=',f7.4)
 602      format(' azi_spot1_disk=',f6.2)
 603      format(' cut_spot1_disk=',f6.4)
c
c  UPDATE November 28, 2001
c
c  Make format statement 604 f6.2 (was f6.3)
c
c
 604      format(' wid_spot1_disk=',f6.2)
 605      format(' TF_spot2_disk=',f7.4)
 606      format(' azi_spot2_disk=',f6.2)
 607      format(' cut_spot2_disk=',f6.4)
 608      format(' wid_spot2_disk=',f6.3)
c
c   NEW BUG August 2, 2001
c  
c   Here are the format statements for the period and T0
c
c   UPDATE August 2, 2004
c
c   Add two digits to the period and T0
c
 700      format(' P=',f14.9)
 701      format(' T0=',f12.6)
 7701      format(' Tconj=',f12.6)
c
c
c  UPDATE November 6, 2002
c
c  Here are the format statements for the limb darkening coefficients.
c
c
 801      format(' x1(U)=',f6.3)
 802      format(' x1(B)=',f6.3)
 803      format(' x1(V)=',f6.3)
 804      format(' x1(R)=',f6.3)
 805      format(' x1(I)=',f6.3)
 806      format(' x1(J)=',f6.3)
 807      format(' x1(H)=',f6.3)
 808      format(' x1(K)=',f6.3)
c
 1801     format(' y1(U)=',f6.3)
 1802     format(' y1(B)=',f6.3)
 1803     format(' y1(V)=',f6.3)
 1804     format(' y1(R)=',f6.3)
 1805     format(' y1(I)=',f6.3)
 1806     format(' y1(J)=',f6.3)
 1807     format(' y1(H)=',f6.3)
 1808     format(' y1(K)=',f6.3)
c
 2801      format(' x2(U)=',f6.3)
 2802      format(' x2(B)=',f6.3)
 2803      format(' x2(V)=',f6.3)
 2804      format(' x2(R)=',f6.3)
 2805      format(' x2(I)=',f6.3)
 2806      format(' x2(J)=',f6.3)
 2807      format(' x2(H)=',f6.3)
 2808      format(' x2(K)=',f6.3)
c
 3801     format(' y2(U)=',f6.3)
 3802     format(' y2(B)=',f6.3)
 3803     format(' y2(V)=',f6.3)
 3804     format(' y2(R)=',f6.3)
 3805     format(' y2(I)=',f6.3)
 3806     format(' y2(J)=',f6.3)
 3807     format(' y2(H)=',f6.3)
 3808     format(' y2(K)=',f6.3)
c
 500      format(a80)
 999      format(7x,i4,2x,a40,3x,f15.8)
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&
c
          subroutine evalfit(islc,Nmodel,xmodel,yinput,
     $      Ndata,xdata,ydata,errdata,index,yout,bestoff)
c
c   This routine will return the value of the model at a specific
c   datapoint
c
          implicit double precision (a-h,o-z)
c
          dimension xmodel(Nmodel),yinput(Nmodel),xdata(Ndata),ydata(Ndata)
          dimension xinter(100000),yinter(100000),errdata(Ndata)
          dimension ydummy(100000),y2(100000),ymodel(100000)
c
          do 10 i=1,Nmodel
            if(islc.gt.0)then
c
c  NEW_BUG  July 5, 2001
c
c  Check to see that the argument of log10 is positive.
c
              if(yinput(i).gt.0.0d0)then
                ymodel(i)=-2.5d0*dlog10(yinput(i))
              else
                ymodel(i)=-99.9d0
              endif
            else
              ymodel(i)=yinput(i)
            endif
 10       continue
c
          ymin=1000.0d0
          ymax=-1000.0d0
c
          call spline(xmodel,ymodel,Nmodel,0.0d0,0.0d0,y2)
c
          do 30 i=1,Ndata
            call splint(xmodel,ymodel,y2,Nmodel,xdata(i),qqq)
            xinter(i)=xdata(i)
            yinter(i)=qqq
 30       continue
c
          call offset(Ndata,xinter,yinter,ydummy,bestoff)
c
          yout=ydummy(index)
c
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
       SUBROUTINE MATINV(ARRAY,NORDER,DET)
       implicit double precision(a-h)
       implicit double precision(o-z)
       double precision array(20,20)
       DIMENSION IK(20),JK(20)
10     DET=1.d0
11     DO 100 K=1,NORDER
C        FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
       AMAX=0.d0
21     DO 30 I=K,NORDER
       DO 30 J=K,NORDER
23     IF(DABS(AMAX)-DABS(ARRAY(I,J))) 24,24,30
24     AMAX=ARRAY(I,J)
       IK(K)=I
       JK(K)=J
30     CONTINUE
C        INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY(K,K)
31     IF(AMAX) 41,32,41
32     DET=0.d0
       GOTO 140
41     I=IK(K)
       IF(I-K) 21,51,43
43     DO 50 J=1,NORDER
       SAVE=ARRAY(K,J)
       ARRAY(K,J)=ARRAY(I,J)
50     ARRAY(I,J)=-SAVE
51     J=JK(K)
       IF(J-K) 21,61,53
53     DO 60 I=1,NORDER
       SAVE=ARRAY(I,K)
       ARRAY(I,K)=ARRAY(I,J)
60     ARRAY(I,J)=-SAVE
C        ACCUMULATE ELEMENTS OF INVERSE MATRIX
61     DO 70 I=1,NORDER
       IF(I-K) 63,70,63
63     ARRAY(I,K)=-ARRAY(I,K)/AMAX
70     CONTINUE
71     DO 80 I=1,NORDER
       DO 80 J=1,NORDER
       IF(I-K) 74,80,74
74     IF(J-K) 75,80,75
75     ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J)
80     CONTINUE
81     DO 90 J=1,NORDER
       IF(J-K) 83,90,83
83     ARRAY(K,J)=ARRAY(K,J)/AMAX
90     CONTINUE
       ARRAY(K,K)=1.d0/AMAX
100    DET=DET*AMAX
C        RESTORE ORDERING OF MATRIX
101    DO 130 L=1,NORDER
       K=NORDER-L+1
       J=IK(K)
       IF(J-K) 111,111,105
105    DO 110 I=1,NORDER
       SAVE=ARRAY(I,K)
       ARRAY(I,K)=-ARRAY(I,J)
110    ARRAY(I,J)=SAVE
111    I=JK(K)
       IF(I-K) 130,130,113
113    DO 120 J=1,NORDER
       SAVE=ARRAY(K,J)
       ARRAY(K,J)=-ARRAY(I,J)
120    ARRAY(I,J)=SAVE
130    CONTINUE
140    RETURN
       END
c
c  #########
c
c
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      implicit double precision (a-h,o-z)
      integer m,mp,n,np,nmax
      real*8 a(np,np),b(np,np)
      PARAMETER (NMAX=50)
c      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      integer i,icol,irow,j,k,l,ll,indxc(nmax),indxr(nmax),ipiv(nmax)
      real*8 big,dum,pivinv
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix in gaussj'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix in gaussj.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE 
        ENDIF
24    CONTINUE
      RETURN
      END
c
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
          subroutine obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c  July 24, 2000
c
c  This routine will take the observed variables specified in the
c  gridloop.opt file and compute the additional chi^2
c
c  UPDATE September 11, 2001
c
c  Change the dimensions of obsparm, obv, sobv, and eobv to 9.  The duration
c  of the X-ray eclipse in degrees is in obsparm(9).
c
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv  to 11
c
c   UPDATE October 10, 2008
c
c   Make the dimensions of obv, obsparm, sobv, eobv to 17
c
          implicit double precision (a-h,o-z)
c
          dimension obv(17),eobv(17),obsparm(17)
          character*40 sobv(17)
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
          ochi=0.0d0
          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
            if(icn.eq.400)index=1
            if(icn.eq.401)index=5
            if(icn.eq.560)index=2
            if(icn.eq.561)index=6
            if(icn.eq.208)index=3
c
c   UPDATE September 21, 2008
c
c   Here are the strings for K1 and K2
c
            if(icn.eq.336)index=10
            if(icn.eq.337)index=11
c
c   UPDATE October 10, 2008
c
c   Here are the strings for incl, mass, ecc, arg, t1, t2
c
            if(icn.eq.269)index=12
            if(icn.eq.384)index=13
            if(icn.eq.130)index=14
            if(icn.eq.17)index=15
            if(icn.eq.624)index=16
            if(icn.eq.625)index=17

c
c   UPDATE February 21, 2002
c
c   Change icn.eq.208  to  icn.eq.209
c
            if(icn.eq.209)index=7
            if(icn.eq.688)index=4
            if(icn.eq.689)index=8
c
c  UPDATE September 11, 2001
c
c  Add this if block.
c
            if(icn.eq.740)index=9
c
c   NEW BUG ALERT  July 24, 2001
c
c   Add this if-then block
c
            if(icn.eq.104)index=-1
            if(icn.eq.112)index=-1
            if(icn.eq.113)index=-1
            if(icn.eq.114)index=-1
            if(icn.eq.115)index=-1
            if(icn.eq.116)index=-1
            if(icn.eq.117)index=-1
            if(icn.eq.118)index=-1
            if(icn.eq.119)index=-1
c
c   UPDATE FEBRUARY 4, 2005
c
c
c   Add this if-then block
c
            if(icn.eq.369)index=-1
c
c   UPDATE May 27, 2002
c
c   If the variable rmed > 1, then define chisq as the absolute deviation,
c   rather than the normal chi^2.
c
c   UPDATE June 3, 2002
c
c   Add this if-then statement.   See below for details.
c
            if((index.eq.9).and.(obv(i).lt.0.0d0))go to 1234
            if(index.ge.1)then
              if(rmed.ge.1.0d0)then
                chisq=dabs((obv(i)-obsparm(index))/(eobv(i)))
              else
                chisq=(obv(i)-obsparm(index))*(obv(i)-obsparm(index))
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+chisq
            endif
c
c   UPDATE June 3, 2002
c
c   If the observed duration of the X-ray eclipse is negative, then
c   the number is meant as an upper limit.  That is, if obv = -10.0,
c   then the eclipse duration is less than 10 degrees.
c
 1234       if((index.eq.9).and.(obv(i).lt.0.0d0))then
c
c   Here is the case where the upper limit is less than the computed
c   eclipse duration.
c
              if(dabs(obv(i)).gt.obsparm(index))go to 1001
c
              if(rmed.ge.1.0d0)THEN
                chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
              else
                chisq=(dabs(obv(i))-obsparm(index))*
     #            (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
              endif
              ochi=ochi+chisq
            endif
c
c   END BUG
c
 1001     continue
c
          return
          end
c
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
          subroutine wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c  July 24, 2000
c
c  This routine will take the observed variables specified in the
c  gridloop.opt file and compute the additional chi^2
c
c
c  UPDATE September 11, 2001
c
c  Change the dimensions of obsparm, obv, sobv, and eobv to 9.  The duration
c  of the X-ray eclipse in degrees is in obsparm(9).
c
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv  to 11
c
c   UPDATE October 10, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv to 17
c
c
          implicit double precision (a-h,o-z)
c
          dimension obv(17),eobv(17),obsparm(17)
          character*40 sobv(17)
c
          character*21 section(17)
          character*80 line
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
          iline=0
          iwrite=0
c
          do 10 i=1,10
            section(i)='1234567890123456'
            section(i)='                '
 10       continue
          line='  -->'
c
          ochi=0.0d0
          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
c
c   NEW BUG  July 24, 2001
c
c   Escape if icn=104 (disk fraction requested).
c
            if(icn.eq.104)go to 1001
c
c   UPDATE FEBRUARY 4, 2005
c
c   Escape if icn=369 (luminosity ratio requested).
c
            if(icn.eq.369)go to 1001
c
c   END BUG
c
c   NEW BUG August 2, 2001
c 
c   Put an escape if one of the variables is not valid.
c
            index=-1
            if(icn.eq.400)index=1
            if(icn.eq.401)index=5
            if(icn.eq.560)index=2
            if(icn.eq.561)index=6
            if(icn.eq.208)index=3
c
c   UPDATE September 21, 2008
c
c   Add K1 and K2
c
            if(icn.eq.336)index=10
            if(icn.eq.337)index=11
c
c
c   UPDATE October 10, 2008
c
c   Here are the strings for incl, mass, ecc, arg, t1, t2
c
            if(icn.eq.269)index=12
            if(icn.eq.384)index=13
            if(icn.eq.130)index=14
            if(icn.eq.17)index=15
            if(icn.eq.624)index=16
            if(icn.eq.625)index=17
c
c   UPDATE February 21, 2002
c
c   Change icn.eq.208  to  icn.eq.209
c
            if(icn.eq.209)index=7
            if(icn.eq.688)index=4
            if(icn.eq.689)index=8
c
c   UPDATE September 11, 2001
c
c   Add this if block
c
            if(icn.eq.740)index=9
c
c   Here is the escape.
c
            if(index.lt.0)go to 1001
c
c   UPDATE May 27, 2002
c
c   If the variable rmed > 1, then define chisq as the absolute
c   deviation, rather than the normal chi^2.
c
c
c   UPDATE June 3, 2002
c
c   If the observed duration of the X-ray eclipse is negative, then
c   the number is meant as an upper limit.  That is, if obv = -10.0,
c   then the eclipse duration is less than 10 degrees.
c
            if(index.eq.9)then
              if(obv(i).lt.0.0d0)then
c
c   Here is the case where the upper limit is less than the computed
c   eclipse duration.
c
                if(dabs(obv(i)).gt.obsparm(index))then
                  chisq=0.0d0
                  iline=iline+1
                  go to 1234
                else
                  if(rmed.ge.1.0d0)THEN
                    chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
                  else
                    chisq=(dabs(obv(i))-obsparm(index))*
     #                (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
                  endif
                  iline=iline+1
                  go to 1234
                endif  
              else
                if(rmed.ge.1.0d0)THEN
                  chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
                else
                  chisq=(dabs(obv(i))-obsparm(index))*
     #                (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
                endif
                iline=iline+1
                go to 1234
              endif
            endif
c
            if(rmed.ge.1.0d0)then
              chisq=dabs((obv(i)-obsparm(index))/(eobv(i)))
            else
              chisq=(obv(i)-obsparm(index))*(obv(i)-obsparm(index))
     $          /(eobv(i)*eobv(i))
            endif
            iline=iline+1
c
c   UPDATE May 27, 2002
c
c   If the variable rmed > 1, then define chisq as the absolute
c   deviation, rather than the normal chi^2.  In that case, write
c   using different format statements.
c
 1234       if(rmed.ge.1.0d0)then
              if(index.eq.1)then
                write(section(iline),91)chisq
                iwrite=0
              endif
              if(index.eq.2)then
                write(section(iline),92)chisq
                iwrite=0
              endif
              if(index.eq.3)then
                write(section(iline),93)chisq
                iwrite=0
              endif
              if(index.eq.4)then
                write(section(iline),94)chisq
                iwrite=0
              endif
              if(index.eq.5)then
                write(section(iline),95)chisq
                iwrite=0
              endif
              if(index.eq.6)then
                write(section(iline),96)chisq
                iwrite=0
              endif
              if(index.eq.7)then
                write(section(iline),97)chisq
                iwrite=0
              endif
              if(index.eq.8)then
                write(section(iline),98)chisq
                iwrite=0
              endif
              if(index.eq.9)then
                write(section(iline),99)chisq
                iwrite=0
              endif
              if(index.eq.10)then
                write(section(iline),302)chisq
                iwrite=0
              endif
              if(index.eq.11)then
                write(section(iline),303)chisq
                iwrite=0
              endif
              if(index.eq.12)then
                write(section(iline),992)chisq
                iwrite=0
              endif
              if(index.eq.13)then
                write(section(iline),993)chisq
                iwrite=0
              endif
              if(index.eq.14)then
                write(section(iline),994)chisq
                iwrite=0
              endif
              if(index.eq.15)then
                write(section(iline),995)chisq
                iwrite=0
              endif
              if(index.eq.16)then
                write(section(iline),996)chisq
                iwrite=0
              endif
              if(index.eq.17)then
                write(section(iline),997)chisq
                iwrite=0
              endif
            else
              if(index.eq.1)then
                write(section(iline),1)chisq
                iwrite=0
              endif
              if(index.eq.2)then
                write(section(iline),2)chisq
                iwrite=0
              endif
              if(index.eq.3)then
                write(section(iline),3)chisq
                iwrite=0
              endif
              if(index.eq.4)then
                write(section(iline),4)chisq
                iwrite=0
              endif
              if(index.eq.5)then
                write(section(iline),5)chisq
                iwrite=0
              endif
              if(index.eq.6)then
                write(section(iline),6)chisq
                iwrite=0
              endif
              if(index.eq.7)then
                write(section(iline),7)chisq
                iwrite=0
              endif
              if(index.eq.8)then
                write(section(iline),8)chisq
                iwrite=0
              endif
              if(index.eq.9)then
                write(section(iline),9)chisq
                iwrite=0
              endif
              if(index.eq.10)then
                write(section(iline),300)chisq
                iwrite=0
              endif
              if(index.eq.11)then
                write(section(iline),301)chisq
                iwrite=0
              endif
              if(index.eq.12)then
                write(section(iline),882)chisq
                iwrite=0
              endif
              if(index.eq.13)then
                write(section(iline),883)chisq
                iwrite=0
              endif
              if(index.eq.14)then
                write(section(iline),884)chisq
                iwrite=0
              endif
              if(index.eq.15)then
                write(section(iline),885)chisq
                iwrite=0
              endif
              if(index.eq.16)then
                write(section(iline),886)chisq
                iwrite=0
              endif
              if(index.eq.17)then
                write(section(iline),887)chisq
                iwrite=0
              endif

            endif
c
c   UPDATE September 11, 2001
c
c   Add this if block
c
c   UPDATE May 27, 2002
c
c   Move this if-statement to the above if-then-else block.
c
c            if(index.eq.9)write(section(iline),9)chisq

            K=lnblnk(line)
            line=line(1:K)//section(iline)
            ochi=ochi+chisq
            if(iline.eq.3)then
              write(*,200)line
              iwrite=100
              iline=0
              line='  -->'
            endif
 1001     continue
c
          if(iline.gt.0.and.iwrite.eq.0)write(*,200)line
c
 1        format(' chi(m1)=',f12.6)
 2        format(' chi(r1)=',f12.6)
 3        format(' chi(g1)=',f12.6)
 4        format(' chi(v1)=',f12.6)
 5        format(' chi(m2)=',f12.6)
882       format(' chi(inc)=',f11.6)
883       format(' chi(Q)=',f13.6)
884       format(' chi(ecc)=',f10.6)
885       format(' chi(arg)=',f10.6)
886       format(' chi(t1)=',f11.6)
887       format(' chi(t2)=',f11.6)

c
c   NEW BUG August 8, 2001
c
c   Correct the format statement below (was chi(r3), is not chi(r2))
c
 6        format(' chi(r2)=',f11.6)
 7        format(' chi(g2)=',f11.6)
 8        format(' chi(v2)=',f11.6)
c
c   UPDATE September 21, 2008
c
c   Here are the formats for K1 and K2
c
 300      format(' chi(k1)=',f11.6)
 301      format(' chi(k2)=',f11.6)
c
c   UPDATE September 11, 2001
c
c   Add this format statement
c
 9        format(' chi(duras)=',f8.3)
c
c   UPDATE May 27, 2002
c 
c   Add these new format statements for median fitting
c
 91        format(' med(m1)=',f12.6)
 92        format(' med(r1)=',f12.6)
 93        format(' med(g1)=',f12.6)
 94        format(' med(v1)=',f12.6)
 95        format(' med(m2)=',f12.6)
 96        format(' med(r2)=',f12.6)
 97        format(' med(g2)=',f12.6)
 98        format(' med(v2)=',f12.6)
 99        format(' med(duras)=',f9.3)
992        format(' med(inc)=',f11.6)
993        format(' med(Q)=',f13.6)
994        format(' med(ecc)=',f11.6)
995        format(' med(arg)=',f11.6)
996        format(' med(t1)=',f12.6)
997        format(' med(t2)=',f12.6)
c
 302       format(' med(k1)=',f11.6)
 303       format(' med(k2)=',f11.6)

 200      format(a80)
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&
c
          subroutine checkRVfit(islc,Nmodel1,xmodel1,
     %            ymodel1,Nmodel2,xmodel2,ymodel2,
     #            Ndata1,xdata1,ydata1,err1,
     %            Ndata2,xdata2,ydata2,err2,
     &            chisq1,chisq2,zero,ifixgamma)
c
c   April 11, 2001
c
c   This routine will fit two velocity curves simultaneously, using
c   a common gamma as a free parameter.  The flag ifixgamma must be 2
c   or larger
c
c

          implicit double precision (a-h,o-z)
c
           dimension xmodel1(Nmodel1),ymodel1(Nmodel1),
     #      xdata1(Ndata1),ydata1(Ndata1),
     %      err1(Ndata1),ymodel2(Nmodel2),xmodel2(Nmodel2),
     #      xdata2(Ndata2),ydata2(Ndata2),err2(Ndata2)
          dimension xinter1(100000),yinter1(100000),y21(100000),ydummy1(100000)
          dimension xinter2(100000),yinter2(100000),y22(100000),ydummy2(100000)
c
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   Find the maximum and minimum y-values of the data
c
          ymin1=10000.
          ymax1=-10000.
          do 20 i=1,Ndata1
            if(ydata1(i).lt.ymin1)ymin1=ydata1(i)
            if(ydata1(i).gt.ymax1)ymax1=ydata1(i)
 20       continue
c
          ymin2=10000.
          ymax2=-10000.
          do 21 i=1,Ndata2
            if(ydata2(i).lt.ymin2)ymin2=ydata2(i)
            if(ydata2(i).gt.ymax2)ymax2=ydata2(i)
 21       continue
c
c   Interpolate the model so that we have y-values at all observed phases.
c
          call spline(xmodel1,ymodel1,Nmodel1,0.0d0,0.0d0,y21)
          call spline(xmodel2,ymodel2,Nmodel2,0.0d0,0.0d0,y22)
c
          do 30 i=1,Ndata1
            call splint(xmodel1,ymodel1,y21,Nmodel1,xdata1(i),qqq)
            xinter1(i)=xdata1(i)
            yinter1(i)=qqq
c            write(67,*)xinter1(i),yinter1(i)
c            write(*,*)xinter1(i),yinter1(i)
 30       continue
c
          do 31 i=1,Ndata2
            call splint(xmodel2,ymodel2,y22,Nmodel2,xdata2(i),qqq)
            xinter2(i)=xdata2(i)
            yinter2(i)=qqq
c            write(68,*)xinter2(i),yinter2(i)
c            write(*,*)'****',xinter1(i),yinter1(i)
 31       continue
c
c   Now find the optimal zero point that will give the lowest chi^2.
c
          call getmean(Ndata1,xdata1,ydata1,dataave1)
          call getmean(Ndata1,xinter1,yinter1,rmodelave1)
          call getmean(Ndata2,xdata2,ydata2,dataave2)
          call getmean(Ndata2,xinter2,yinter2,rmodelave2)

          savezero=zero
          zero=(dataave1-rmodelave1)
          step=abs(rmodelave1-dataave1)/1000.0d0

          chi1=0.0d0
          chi2=0.0d0
          chi3=0.0d0
          chisq1=0.0d0
          chisq2=0.0d0
          small=1.0d35
c
          do 750 i=1,20
            zero1=zero
            call offset(Ndata1,xinter1,yinter1,ydummy1,zero1)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi11)
            call offset(Ndata2,xinter2,yinter2,ydummy2,zero1)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi12)
            chi1=chi11+chi12
            if(chi1.lt.small)then
              small=chi1
              zerosmall=zero1
            endif
            fn=0.0d0
            zero2=zero+step
            call offset(Ndata1,xinter1,yinter1,ydummy1,zero2)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi21)
            call offset(Ndata2,xinter2,yinter2,ydummy2,zero2)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi22)
            chi2=chi21+chi22
            if(chi2.lt.small)then
              small=chi2
              zerosmall=zero2
            endif
            diff=chi1-chi2
            if(diff.gt.0.0d0)go to 5061
            step=-step
c  
            csave=chi1
            chi1=chi2
            chi2=csave
            zsave=zero1
            zero1=zero2
            zero2=zsave
c
 5061       fn=fn+1.0d0
c
            zero3=zero2+step
            call offset(Ndata1,xinter1,yinter1,ydummy1,zero3)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi31)
            call offset(Ndata2,xinter2,yinter2,ydummy2,zero3)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi32)
            chi3=chi31+chi32
            if(chi3.lt.small)then
              small=chi3
              zerosmall=zero3
            endif
            diff23=chi3-chi2
            if(diff23.lt.0.0d0)then
              chi1=chi2
              chi2=chi3
              zero1=zero2
              zero2=zero3
              zero=zero3
              go to 5061
            endif          
c
c         find the minimum of parabola defined by the last three points
c
 5081       if((chi3-chi2).eq.0.0d0)go to 999
            step=step*(1.0d0/(1.0d0+(chi1-chi2)/(chi3-chi2))+0.5d0)
            zero=zero2-step
            step=step*fn/3.0d0
            call offset(Ndata1,xinter1,yinter1,ydummy1,zero)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi41)
            call offset(Ndata2,xinter2,yinter2,ydummy2,zero)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi42)
            chi4=chi41+chi42
            if(chi4.lt.small)then
              small=chi4
              zerosmall=zero
            endif
c
 750      continue  ! loop over grid searches
c
 751      continue                  ! come here when delta chi is small
c
 999      zero=zerosmall
c
          call offset(Ndata1,xinter1,yinter1,ydummy1,zero)
          call getchi(Ndata1,ydata1,err1,ydummy1,chisq1)
          call offset(Ndata2,xinter2,yinter2,ydummy2,zero)
          call getchi(Ndata2,ydata2,err2,ydummy2,chisq2)
c
          return
c
          end
c
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  NEW BUG August 2, 2001
c
c  Add the period and T0 to the argument list
c
c  UPDATE August 13, 2001
c
c  Make the variable line 200 characters.
c
c
c  UPDATE January 15, 2002
c
c  Add alb1 and alb2 to the list of variables.
c
c  UPDATE November 6, 2002
c
c  Add dwavex and dwavey (limb darkening coefficients) to the 
c  argument list of newwritevar.
c  Dimension them as dwavex(8,2), dwavey(8,2).
c
c  July 29, 2005
c
c  Add ecosw and temprat to the list.
c
          subroutine newwritevar(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     %       t3,g3,SA3,ecc,argper,pshift,spot1parm,spot2parm,spotdparm,
     $       period,T0,chisq,nnn,line,alb1,alb2,dwavex,dwavey,
     @       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     %       bigI,bigbeta,powercoeff,density,sw5,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
c   November 15, 1999
c
c   This routine will determine which variables need to be printed on the
c   screen based on the string codes in svar(1:Nvarmax).  The output
c   printed on the screen will be compact.
c
          implicit double precision (a-h,o-z)
c
c   UPDATE October 28, 2002
c
c   Make the length of line 300 (was 200)
c
c
          character*300 line        ! was 132
          character*40 svar(Nvarmax)
          character*10 string10
          character*11 string11
          character*12 string12
          character*13 string13
          character*14 string14
          character*15 string15
          character*16 string16
          character*17 string17
          character*18 string18
          character*19 string19
          character*22 string22
          character*23 string23
          character*7 string7
          character*8 string8
          character*5 string5
          character*6 string6
          character*9 string9

          dimension var(Nvarmax)
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension dwavex(8,2),dwavey(8,2)
c
          line=' '
          iline=0
          iwrite=0
          ilength=1
c
          ilength=ilength+7
          write(string7,2000)nnn
          K=lnblnk(line)
          line=line(1:K)//string7
c
          ilength=ilength+15
c
c   NEW_BUG July 5, 2001
c
c   Check that the value of chi^ does not exceed 99,999,999
c
          tempchi=chisq
          if(chisq.gt.99999999.0d0)tempchi=99999999.9d0
          write(string17,2001)tempchi
c
c          write(string12,2001)chisq
c
          K=lnblnk(line)
          line=line(1:K)//string17
c

 2000     format(' ',i6)
 2001     format(' ',f16.7)

          iargper=0
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
            if(kk.eq.144)then   !use string e1
              iline=iline+1
              iwrite=0
              ilength=ilength+9
              write(string9,781)beam1
              K=lnblnk(line)
              line=line(1:K)//string9
            endif
c
            if(kk.eq.145)then   !use string e2
              iline=iline+1
              iwrite=0
              ilength=ilength+9
              write(string9,781)beam2
              K=lnblnk(line)
              line=line(1:K)//string9
            endif
c
            if(kk.eq.492)then   !use string pm
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,777)primmass
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.497)then   !use string pr
              iline=iline+1
              iwrite=0
              ilength=ilength+11
              write(string11,778)primrad
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.490)then   !use string pk
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,779)primk
              K=lnblnk(line)
              line=line(1:K)//string13
            endif
c
            if(kk.eq.544)then   !use string ra
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,776)ratrad
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.528)then   !use string q1
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,203)frac1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.529)then   !use string q2
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,203)frac2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.612)then   !temprat
              iline=iline+1
              iwrite=0
              ilength=ilength+11
              write(string11,778)temprat
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.100)then   !density
              iline=iline+1
              iwrite=0
              ilength=ilength+11
              write(string11,780)density
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c

 776        format(' ',f11.6)
 777        format(' ',f9.5)
 778        format(' ',f10.5)
 779        format(' ',f12.8)
780         format(' ',f10.7)
781         format(' ',f8.4)
c
c   UPDATE January 15, 2002
c
c   Here are the blocks for alb1 and alb2
c
c   UPDATE MARCH 20, 2004
c
c   USE string 10 in the two blocks below
c
            if(kk.eq.368)then   !use string l1
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,202)alb1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.369)then   !use string l2
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,202)alb2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.484)then
              iline=iline+1
              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 15 (add two decimals in format 700 below).
c
              ilength=ilength+15
              write(string15,700)period
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.623)then
              iline=iline+1
              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 14 (add two decimals in format 701 below).
c
              ilength=ilength+14
c
c  UPDATE May 21, 2002
c
c  Change (string11,200) to (string12,701)
c
              write(string14,701)T0
              K=lnblnk(line)
c
c   UPDATE May 27, 2002
c
c   Change string11 to string12.
c
              line=line(1:K)//string14
            endif

            if(kk.eq.610)then          !Tconj
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,701)Tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.269)then
              iline=iline+1
              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 9 (add one decimal in format 200 below).
c
c   UPDATE September 20, 2007
c
c   Use length 12 (add three decimals in format 20000 below).
c
              ilength=ilength+12
              write(string12,20000)finc
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
c   UPDATE May 8, 2006
c
c   Add bigI and bigbeta
c
            if(kk.eq.8)then     !bigI, string aI
              iline=iline+1
              iwrite=0
              ilength=ilength+9
              write(string9,200)bigI
              K=lnblnk(line)
              line=line(1:K)//string9
            endif
c
            if(kk.eq.1)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,88201)bigbeta
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.384)then
              iline=iline+1
              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 12 (add four decimals in format 201 below).
c
              ilength=ilength+12
              write(string12,201)Q
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.176)then
              iline=iline+1
              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 10 (add twp decimals in format 202 below).
c
              ilength=ilength+10
              write(string10,202)fill1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.177)then
              iline=iline+1
              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 10 (add twp decimals in format 203 below).
c
              ilength=ilength+10
              write(string10,203)fill2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.552)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,204)rinner
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.464)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,205)omega1
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.465)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,206)omega2
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.558)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,207)router
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.626)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,217)t3  
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.210)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,218)g3  
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.36)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,211)betarim
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.624)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,209)Teff1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.625)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,210)Teff2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.744)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,212)xi
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.611)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,208)tdisk 
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.375)then
              iline=iline+1
              iwrite=0
              ilength=ilength+11
              write(string11,213)rLx
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.580)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,214)separ
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.576)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,219)SA3
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
c  UPDATE May 24, 2002
c
c  Make the length 9, as in ilength=ilength+9, and
c  write(string9,223), and //string9
c
            if(kk.eq.498)then
              iline=iline+1
              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 13 (add four decimals in format 223 below).
c
              ilength=ilength+13
              write(string13,223)pshift
              K=lnblnk(line)
              line=line(1:K)//string13
            endif
c
            if(kk.eq.130)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,221)ecc
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.17)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,222)argper
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.48)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,401)spot1parm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.49)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,402)spot1parm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.50)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,403)spot1parm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.51)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,404)spot1parm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.52)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,405)spot1parm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.53)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,406)spot1parm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.54)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,407)spot1parm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.55)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,408)spot1parm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.80)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,501)spot2parm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.81)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,502)spot2parm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.82)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,503)spot2parm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.83)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,504)spot2parm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.84)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,505)spot2parm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.85)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,506)spot2parm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.86)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,507)spot2parm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.87)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,508)spot2parm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.112)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,601)spotdparm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.113)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,602)spotdparm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.114)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,603)spotdparm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.115)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,604)spotdparm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.116)then
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,605)spotdparm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.117)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,606)spotdparm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.118)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,607)spotdparm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.119)then
              iline=iline+1
              iwrite=0
              ilength=ilength+7
              write(string7,608)spotdparm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
c
c   UPDATE November 6, 2002
c
c   Add the limb darkening coefficients here.
c
c   x-coefficient, star 1
c
            if(kk.eq.752)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,801)dwavex(1,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.753)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,802)dwavex(2,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.754)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,803)dwavex(3,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.755)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,804)dwavex(4,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.756)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,805)dwavex(5,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.757)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,806)dwavex(6,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.758)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,807)dwavex(7,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.759)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,808)dwavex(8,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
c   x-coefficient, star 2
c
            if(kk.eq.816)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2801)dwavex(1,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.817)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2802)dwavex(2,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.818)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2803)dwavex(3,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.819)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2804)dwavex(4,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.820)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2805)dwavex(5,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.821)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2806)dwavex(6,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.822)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2807)dwavex(7,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.823)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,2808)dwavex(8,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
c   y-coefficient, star 1
c
            if(kk.eq.784)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,1801)dwavey(1,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.785)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,1802)dwavey(2,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.786)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,1803)dwavey(3,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.787)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,1804)dwavey(4,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.788)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,1805)dwavey(5,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.789)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,1806)dwavey(6,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.790)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,1807)dwavey(7,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.791)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,1808)dwavey(8,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
c   y-coefficient, star 2
c
            if(kk.eq.720)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,3801)dwavey(1,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.721)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,3802)dwavey(2,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.722)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,3803)dwavey(3,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.723)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,3804)dwavey(4,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.724)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,3805)dwavey(5,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.725)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,3806)dwavey(6,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.726)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,3807)dwavey(7,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.727)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,3808)dwavey(8,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.111)then   ! ecosw, string dphi
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,202)ecosw
              K=lnblnk(line)
              line=line(1:K)//string10
              iargper=100
            endif
c
 10       continue
c
c  July 29, 2005
c
c  If the variable ecosw is assigned, then also print out the argument
c  of periastron
c 
c
            if(iargper.gt.0)then
              call getom(ecc,ecosw,argper)
              iline=iline+1
              iwrite=0
              ilength=ilength+8
              write(string8,222)argper
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
c   UPDATE AUGUST 10, 2005
c
c   Make the gammas string10
c
          if(isvel1.gt.0)then
            iline=iline+1
            iwrite=0
            ilength=ilength+10
            write(string10,215)gamma1
            K=lnblnk(line)
            line=line(1:K)//string10
          endif
c
          if(isvel2.gt.0)then
            iline=iline+1
            iwrite=0
            ilength=ilength+10
            write(string10,215)gamma2
            K=lnblnk(line)
            line=line(1:K)//string10
          endif
c
c   UPDATE November 6, 2008
c
c   if asini (sw5) > 0, record it here
c
          if(sw5.gt.0.0d0)then
            iline=iline+1
            iwrite=0
            ilength=ilength+12
            write(string12,20000)sw5
            K=lnblnk(line)
            line=line(1:K)//string12
          endif
c
3215      format(' ',f11.8)
c
c    UPDATE November 5, 2008
c   
c    If ecc > 0, record the times of inferior and superior
c    conjunction.
c
          if(ecc.gt.0.0d0)then
            call getcontimes(finc,period,ecc,argper,T0,tconj1,tconj2)
            iline=iline+1
            iwrite=0
            ilength=ilength+14
            write(string14,701)tconj1
            K=lnblnk(line)
            line=line(1:K)//string14
            iline=iline+1
            iwrite=0
            ilength=ilength+14
            write(string14,701)tconj2
            K=lnblnk(line)
            line=line(1:K)//string14
          endif
c
c
c
c   UPDATE October 20, 2002
c
c   Use length 9 (add one decimal in format 200 below).
c
 200      format(' ',f8.5)
20000     format(' ',f11.7)
c
c   UPDATE October 20, 2002
c
c   Use length 12 (add four decimals in format 201 below).
c
 201      format(' ',f11.8)
c
c   UPDATE October 20, 2002
c
c   Use length 10 (add twp decimals in format 202 below).
c
 202      format(' ',f9.7)
c
c   UPDATE October 20, 2002
c
c   Use length 10 (add twp decimals in format 203 below).
c
 203      format(' ',f9.7)
 204      format(' ',f7.5)
 205      format(' ',f7.4)
 206      format(' ',f7.4)
 207      format(' ',f7.5)
 208      format(' ',f9.2)
 209      format(' ',f9.2)
 210      format(' ',f9.2)
 211      format(' ',f6.3)
 212      format(' ',f7.4)
 213      format(' ',1pe10.4)
 214      format(' ',f9.4)
 215      format(' ',f9.4)
 216      format(' ',f7.2)
 217      format(' ',f9.2)
 218      format(' ',f9.2)
 219      format(' ',f9.4)
 220      format(' ',f7.4)
 221      format(' ',f7.5)
 222      format(' ',f7.3)
c
c  UPDATE May 24, 2002
c
c  Add this format statement
c
c
c   UPDATE October 20, 2002
c
c   Use length 13 (add four decimals in format 223 below).
c
 223      format(' ',f12.9)
c
c   RVG BUG ALERT   May 8, 2001
c
c   Add these statements for spots
c
 401      format(' ',f7.4)
 402      format(' ',f6.2)
 403      format(' ',f6.2)
 404      format(' ',f6.3)
 405      format(' ',f7.4)
 406      format(' ',f6.2)
 407      format(' ',f6.2)
 408      format(' ',f6.3)
 501      format(' ',f7.4)
 502      format(' ',f6.2)
 503      format(' ',f6.2)
 504      format(' ',f6.3)
 505      format(' ',f7.4)
 506      format(' ',f6.2)
 507      format(' ',f6.2)
 508      format(' ',f6.3)
 601      format(' ',f7.4)
 602      format(' ',f6.2)
 603      format(' ',f6.4)
 604      format(' ',f6.3)
 605      format(' ',f7.4)
 606      format(' ',f6.2)
 607      format(' ',f6.4)
 608      format(' ',f6.3)
c
c   NEW BUG August 2, 2001
c
c   Add these format statements for period and T0
c
c
c   UPDATE October 20, 2002
c
c   Use length 15 (add two decimals in format 700 below).
c

 700      format(' ',f14.9)
c
c   UPDATE May 27, 2002
c
c   Change the format from f10.4 to f11.5.
c
c   UPDATE October 20, 2002
c
c   Use length 14 (add two decimals in format 701 below).
c
 701      format(' ',f13.7)

 500      format(a132)
 999      format(7x,i4,2x,a40,3x,f15.8)
c
c
c  UPDATE November 6, 2002
c
c  Here are the format statements for the limb darkening coefficients.
c
c
 801      format(' ',f11.8)
 802      format(' ',f11.8)
 803      format(' ',f11.8)
 804      format(' ',f11.8)
 805      format(' ',f11.8)
 806      format(' ',f11.8)
 807      format(' ',f11.8)
 808      format(' ',f11.8)
c
 1801     format(' ',f11.8)
 1802     format(' ',f11.8)
 1803     format(' ',f11.8)
 1804     format(' ',f11.8)
 1805     format(' ',f11.8)
 1806     format(' ',f11.8)
 1807     format(' ',f11.8)
 1808     format(' ',f11.8)
c
 2801     format(' ',f11.8)
 2802     format(' ',f11.8)
 2803     format(' ',f11.8)
 2804     format(' ',f11.8)
 2805     format(' ',f11.8)
 2806     format(' ',f11.8)
 2807     format(' ',f11.8)
 2808     format(' ',f11.8)
c
 3801     format(' ',f11.8)
 3802     format(' ',f11.8)
 3803     format(' ',f11.8)
 3804     format(' ',f11.8)
 3805     format(' ',f11.8)
 3806     format(' ',f11.8)
 3807     format(' ',f11.8)
 3808     format(' ',f11.8)
c
88201     format(' ',f9.5)

          return
          end
c
c
c  NEW BUG  July 24, 2001
c
c  Here is a new subroutine that will compute the chi^2 of the disk fraction.
c
         subroutine diskchi(Nphase,Nmaxphase,ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)
c
c
c
          implicit double precision (a-h,o-z)
c
c  UPDATE September 21, 2008
c
c  make the dimension of obv, eobv, and sobv 11
c
          dimension ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      ymods3(Nmaxphase)
          dimension obv(17),eobv(17)
          dimension xr(200000),rat(200000)
c
c   UPDATE October 27, 2008
c
c   Make ochidisk have dimension 8.  This will allow the user to fit
c   for the disk fraction in more than 1 band at a time.  Also, put
c   the array compfracs into the argument list above
c
          dimension ochidisk(8),compfracs(8,2) 
          character*40 sobv(11)
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
         do 1 i=1,8
            ochidisk(i)=0.0d0
1        continue

c
c          iii=Nphase/8
c          frac=ymodd(iii)/(ymods1(iii)+ymods2(iii)+ymods3(iii)+ymodd(iii))
c
c          icount=0
c          do 10 i=1,Nphase           
c            icount=icount+1
c            rat(icount)=ymodd(i)/(ymods1(i)+ymods2(i)+ymods3(i)+ymodd(i))
c            xr(icount)=ymods1(i)
c 10       continue
cc
c          if(icount.le.2)then
c            ochilr=9999999999.0d0
c            ochi=ochi+ochilr
c            return
c          endif
c
c          call sort2(icount,rat,xr)
cc
c          q1=dble(icount/2)
c          q2=dble(icount)/2.0d0
c          if(q1.eq.q2)then
c            rmedian=(rat(icount/2)+rat(icount/2+1))/2.0d0
c          else
c            rmedian=rat(icount/2+1)
c          endif
c          frac=rmedian
cc  

          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
c
c  UPDATE May 27, 2002
c
c  If the variable rmed > 1, then define ochidisk to be the absolute
c  deviation, rather than the normal chi^2.
c
c  UPDATE October 27, 2008
c
c  use the tags d1 to d8 to fit for the disk fraction in filter 1, filter 2,
c  etc.  String d1=112, d2=113, ... d9=119
c
            if(icn.eq.112)then
              frac=compfracs(1,2)
              if(rmed.ge.1.0d0)then
                ochidisk(1)=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochidisk(1)=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+ochidisk(1)
            endif
c
            if(icn.eq.113)then
              frac=compfracs(2,2)
              if(rmed.ge.1.0d0)then
                ochidisk(2)=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochidisk(2)=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+ochidisk(2)
            endif
c
            if(icn.eq.114)then
              frac=compfracs(3,2)
              if(rmed.ge.1.0d0)then
                ochidisk(3)=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochidisk(3)=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+ochidisk(3)
            endif
c
            if(icn.eq.115)then
              frac=compfracs(4,2)
              if(rmed.ge.1.0d0)then
                ochidisk(4)=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochidisk(4)=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+ochidisk(4)
            endif
c
            if(icn.eq.116)then
              frac=compfracs(5,2)
              if(rmed.ge.1.0d0)then
                ochidisk(5)=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochidisk(5)=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+ochidisk(5)
            endif
c
            if(icn.eq.117)then
              frac=compfracs(6,2)
              if(rmed.ge.1.0d0)then
                ochidisk(6)=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochidisk(6)=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+ochidisk(6)
            endif
c
            if(icn.eq.118)then
              frac=compfracs(7,2)
              if(rmed.ge.1.0d0)then
                ochidisk(7)=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochidisk(7)=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+ochidisk(7)
            endif
c
            if(icn.eq.119)then
              frac=compfracs(8,2)
              if(rmed.ge.1.0d0)then
                ochidisk(8)=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochidisk(8)=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+ochidisk(8)
            endif
c
 1001     continue
c

          return
          end
c
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c
c  NEW BUG  July 24, 2001
c
c  Here is a new subroutine that will compute the chi^2 of the disk fraction.
c
          subroutine wdiskobschi(ochidisk,Nobv,sobv)
c
c  July 24, 2000
c
c  This routine will take the observed variables specified in the
c  gridloop.opt file and compute the additional chi^2
c
          implicit double precision (a-h,o-z)
c
c  UPDATE September 21, 2008
c
c  Make the dimension of obv, eobv, sobv to 11
c
          dimension obv(17),eobv(17),ochidisk(8)
          character*40 sobv(Nobv)
c
          character*22 section
          character*80 line
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
          iline=0
          iwrite=0
c
          section='1234567890123456'
          section='                '
          line='  -->'
c
c          ochi=0.0d0 
c
c   UPDATE May 27, 2002
c
c   If rmed > 1, then we are doing median fitting.  In that case, write
c   to a different format statement.
c
c
c   UPDATE October 27, 2008
c
c   Update this routine so that one can specify which filter the disk
c   fraction is, and also have the ability to specify two or more disk
c   fractions.
c
           do 99 ii=1,Nobv
             icn=icnvrt(sobv(ii)(1:2))

             if(icn.eq.112)then   ! string d1
               if(rmed.ge.1.0d0)then
                 write(section,21)ochidisk(1)
               else
                 write(section,11)ochidisk(1)
               endif
               K=lnblnk(line)
               line=line(1:K)//section
c
 11            format(' chi_Ufrac=',f11.6)
 21            format(' med_Ufrac=',f11.6)
             endif
c
             if(icn.eq.113)then   ! string d2
               if(rmed.ge.1.0d0)then
                 write(section,22)ochidisk(2)
               else
                 write(section,12)ochidisk(2)
               endif
               K=lnblnk(line)
               line=line(1:K)//section

 12            format(' chi_Bfrac=',f11.6)
 22            format(' med_Bfrac=',f11.6)
             endif
c
             if(icn.eq.114)then   ! string d3
               if(rmed.ge.1.0d0)then
                 write(section,23)ochidisk(3)
               else
                 write(section,13)ochidisk(3)
               endif
               K=lnblnk(line)
               line=line(1:K)//section

 13            format(' chi_Vfrac=',f11.6)
 23            format(' med_Vfrac=',f11.6)
             endif
c
             if(icn.eq.115)then   ! string d4
               if(rmed.ge.1.0d0)then
                 write(section,24)ochidisk(4)
               else
                 write(section,14)ochidisk(4)
               endif
               K=lnblnk(line)
               line=line(1:K)//section

 14            format(' chi_Rfrac=',f11.6)
 24            format(' med_Rfrac=',f11.6)
             endif
c
             if(icn.eq.116)then   ! string d5
               if(rmed.ge.1.0d0)then
                 write(section,25)ochidisk(5)
               else
                 write(section,15)ochidisk(5)
               endif
               K=lnblnk(line)
               line=line(1:K)//section

 15            format(' chi_Ifrac=',f11.6)
 25            format(' med_Ifrac=',f11.6)
             endif
c
             if(icn.eq.117)then   ! string d6
               if(rmed.ge.1.0d0)then
                 write(section,26)ochidisk(6)
               else
                 write(section,16)ochidisk(6)
               endif
               K=lnblnk(line)
               line=line(1:K)//section

 16            format(' chi_Jfrac=',f11.6)
 26            format(' med_Jfrac=',f11.6)
             endif
c
             if(icn.eq.118)then   ! string d7
               if(rmed.ge.1.0d0)then
                 write(section,27)ochidisk(7)
               else
                 write(section,17)ochidisk(7)
               endif
               K=lnblnk(line)
               line=line(1:K)//section

 17            format(' chi_Hfrac=',f11.6)
 27            format(' med_Hfrac=',f11.6)
             endif
c
             if(icn.eq.119)then   ! string d8
               if(rmed.ge.1.0d0)then
                 write(section,28)ochidisk(8)
               else
                 write(section,18)ochidisk(8)
               endif
               K=lnblnk(line)
               line=line(1:K)//section

 18            format(' chi_Kfrac=',f11.6)
 28            format(' med_Kfrac=',f11.6)
             endif
c

99        continue

          write(*,200)line

 1        format(' chi(frac)=',f11.6)
c
c   UPDATE May 27, 2002
c
c   Add this new format statement for median fitting.
c
 2        format(' med(frac)=',f11.6)
 200      format(a80)

          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  UPDATE JUNE 19, 2008
c
c  Add ikeep, ecc, and argper as arguments.
c
c   NEW BUG August 2, 2001
c
c   These new routines will allow for the fitting of period and T0
c
          subroutine phasefold(Ndatamax,isavN,savx,savy,saverr,
     #      Nbin,Ndata,xdata,ydata,errdata,period,T0,ikeep,ecc,argper,
     #      itime)
c
c   This routine will take data arrays in savex,savey,saverr 
c   (time, mag, magerr), and fold them into the arrays
c   xdata,ydata,errdata.  The latter arrays are the ones normally
c   used in the optimization programs.
c 
          implicit double precision(a-h,o-z)
c
          dimension savx(Ndatamax),savy(Ndatamax),saverr(Ndatamax)
          dimension xdata(Ndatamax),ydata(Ndatamax),errdata(Ndatamax)
c
          parameter(pie=3.14159265358979323d0)
c
c   UPDATE JUNE 19, 2008
c
c   Figure out the conjunction phases if the orbit is eccentric.
c
          if(itime.eq.2)then
            if(Ndata.gt.1)call sort3(Ndata,xdata,ydata,errdata)
            return
          endif

          eshift=0.0d0
          pconj=pie
          pconj2=0.0d0

          if((ecc.gt.0.0d0).and.(ikeep.gt.0))then
            argrad=argper*pie/180.0d0
c
c  Here is some code adapted from Wilson-Devinney to keep track of
c  phases needed for eccentric orbits.
c
            trc=0.5d0*pie-argrad
 1139       if(trc.lt.0.d0) trc=trc+2.0d0*pie
            if(trc.lt.0.d0) goto 1139
 1140       if(trc.ge.2.0d0*pie) trc=trc-2.0d0*pie
            if(trc.ge.2.0d0*pie) goto 1140
            htrc=0.5d0*trc
            if(dabs(0.5*pie-htrc).lt.7.d-6) goto 11101
            if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 11101
            ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
            goto 11103
11101       ecan=pie
11103       xmc=ecan-ecc*dsin(ecan)
            if(xmc.lt.0.d0) xmc=xmc+2.0d0*pie
            phper=1.d0-xmc/(2.0d0*pie)
            pconj=(xmc+argrad)/(2.0d0*pie)-0.25d0
c
c   Make sure the conjunction phase is between 0 and 1
c
            if(pconj.gt.1.0d0)pconj=pconj-1.0d0

            trc=0.5d0*pie-argrad+pie
 3139       if(trc.lt.0.d0) trc=trc+2.0d0*pie
            if(trc.lt.0.d0) goto 3139
 3140       if(trc.ge.2.0d0*pie) trc=trc-2.0d0*pie
            if(trc.ge.2.0d0*pie) goto 3140
            htrc=0.5d0*trc
            if(dabs(0.5*pie-htrc).lt.7.d-6) goto 31101
            if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 31101
            ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
            goto 31103
31101       ecan=pie
31103       xmc=ecan-ecc*dsin(ecan)
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
            if(ikeep.eq.1)eshift=phper+pconj-0.5
            if(ikeep.eq.2)eshift=phper2+pconj2
          endif

c   Here is the case where no binning is requested.
c
          pshift=0.0d0
          tshift=0.0d0
          if(ikeep.eq.1)tshift=pshift+eshift-pconj
          if(ikeep.eq.2)tshift=pshift+eshift-pconj2

          if(nbin.eq.0)then
            do 10 i=1,isavN
              t1=(savx(i)-T0)/period+tshift
              xdata(i)=dmod(t1,1.0d0)
              ydata(i)=savy(i)
              errdata(i)=saverr(i)
 10         continue
            Ndata=isavN
c
           do 15 i=1,Ndata
             if(xdata(i).lt.0.0d0)xdata(i)=xdata(i)+1.0d0
             if(xdata(i).gt.1.0d0)xdata(i)=xdata(i)-1.0d0
             if(xdata(i).lt.0.0d0)xdata(i)=xdata(i)+1.0d0
             if(xdata(i).gt.1.0d0)xdata(i)=xdata(i)-1.0d0
             if(xdata(i).lt.0.0d0)xdata(i)=xdata(i)+1.0d0
             if(xdata(i).gt.1.0d0)xdata(i)=xdata(i)-1.0d0
 15        continue
c
c   Sort the data by phase.
c
c   UPDATE April 15, 2002
c
c   Put an if-then clause to cover the case when Ndata=1
c
c
            if(Ndata.gt.1)call sort3(Ndata,xdata,ydata,errdata)
c
c   We are done, so leave.
c
            return
          endif
c

          return
          end
c
c
c
          subroutine kopydata(Ndatamax,Ndata,xdata,ydata,errdata,isavN,
     $       savx,savy,saverr)
c
c   August 2, 2001
c
c   When the period and/or T0 is being adjusted, the data is entered
c   in time units.  We have to save the original arrays because the
c   data can be multiply folded during a single run of a program.
c
          implicit double precision(a-h,o-z)
c
          dimension savx(Ndatamax),savy(Ndatamax),saverr(Ndatamax)
          dimension xdata(Ndatamax),ydata(Ndatamax),errdata(Ndatamax)
c
          isavN=Ndata
c
          do 10 i=1,Ndata
            savx(i)=xdata(i)
            savy(i)=ydata(i)
            saverr(i)=errdata(i)
 10       continue
c
          return
c
          end
c
c  ##########################################################         
c
          subroutine wELCdata(Ndatamax,Ndata,xdata,ydata,err,fileout)
c
c  NEW BUG August 2, 2001
c
c  This routine is new.  When the period and/or T0 is being fitted for,
c  the folded data files will be written.
c
          implicit double precision(a-h,o-z)
c
          dimension xdata(Ndatamax),ydata(Ndatamax),err(Ndatamax)
c
c   UPDATE November 18, 2002
c
c   Change to character*(*) below.
c
          character*(*) fileout

          open(unit=51,file=fileout,status='unknown')
c
          do 10 i=1,Ndata
            write(51,100)xdata(i),ydata(i),err(i)
 10       continue
c
          close(51)
c
 100      format(f14.11,3x,f18.13,3x,f17.13)
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   UPDATE August 13, 2001
c
c   This is a new subroutine, which is the inverse of assignvar.  That it,
c   the routine has values for the parameters (fill1, fill2, omega1, ...)
c   and it will assign them to the var string based on the variables
c   selected in gridloop.opt file.
c
c
c  UPDATE January 15, 2002
c
c  Update the routine to include the albedos of the two stars (alb1, alb2)
c
c
c  UPDATE November 6, 2002
c
c  Add limb darkening coefficients dwavex and dwavey to the
c  argument list of assignvar and varassign.
c
c  UPDATE August 13, 2004
c  
c  Add primmass,primK,primrad,ratrad,frac1,frac2
c
          subroutine varassign(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     &       pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     #       dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c  UPDATE October 10, 2008
c
c  Add the density to the list
c

          implicit double precision (a-h,o-z)
c
          character*40 svar(Nvarmax)
          dimension var(Nvarmax)
c       
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension powercoeff(8,9)
c
c
c  UPDATE November 6, 2002
c
c  Dimension the limb darkening variables.
c
          dimension dwavex(8,2),dwavey(8,2)
c
          write(*,*)' '
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
            if(kk.eq.144)then  !beam1, use string e1
              var(i)=beam1
            endif
c
            if(kk.eq.145)then  !Tconj, use string e2
              var(i)=beam2
            endif
c
            if(kk.eq.610)then  !Tconj, use string tc
              var(i)=Tconj
            endif
c
            if(kk.eq.100)then  !density, use string de
              var(i)=density
            endif
c
            if(kk.eq.8)then  !bigI, use string ai
              var(i)=bigI
            endif
c
            if(kk.eq.1)then  !bigbeta, use string ab
              var(i)=bigbeta
            endif
c
            if(kk.eq.111)then  !ecosw, use string dphi
              var(i)=ecosw
            endif
c
            if(kk.eq.612)then  !temprat
              var(i)=temprat
            endif
c
            if(kk.eq.492)then  !primmass, use string pm
              var(i)=primmass
            endif
c
            if(kk.eq.497)then  !primrad, use string pr
              var(i)=primrad
            endif
c
            if(kk.eq.490)then  !primK, use string pk
              var(i)=primK
            endif
c
            if(kk.eq.544)then  !ratrad, use string ra
              var(i)=ratrad
            endif
c
            if(kk.eq.528)then  !frac1, use string q1
              var(i)=frac1
            endif
c
            if(kk.eq.529)then  !frac2, use string q2
              var(i)=frac2
            endif
c
c
c  UPDATE January 15, 2002
c
c  Here are the assignments for alb1 and alb2
c
            if(kk.eq.368)then  !alb1
              var(i)=alb1
            endif
c
            if(kk.eq.369)then  !alb2
              var(i)=alb2
            endif
c
            if(kk.eq.484)then  !period
              var(i)=period
            endif
c
            if(kk.eq.623)then ! T0
              var(i)=T0
            endif
c
            if(kk.eq.112)then   ! temperature factor spot 1 on disk
              var(i)=spotdparm(1,1)
            endif
c
            if(kk.eq.116)then   ! temperature factor spot 2 on disk
              var(i)=spotdparm(2,1)
            endif
c
            if(kk.eq.113)then  ! azimuth spot 1 on disk
              var(i)=spotdparm(1,2)
            endif   
c
            if(kk.eq.117)then  ! azimuth spot 2 on disk
              var(i)=spotdparm(2,2)
            endif   
c
            if(kk.eq.114)then    ! cutoff radius for spot 1 on disk
              var(i)=spotdparm(1,3)
            endif
c
            if(kk.eq.118)then    ! cutoff radius for spot 2 on disk
              var(i)=spotdparm(2,3)
            endif
c
            if(kk.eq.115)then    ! angular width of spot 1 on disk
              var(i)=spotdparm(1,4)
            endif
c
            if(kk.eq.119)then    ! angular width of spot 2 on disk
              var(i)=spotdparm(2,4)
            endif
c
c
            if(kk.eq.80)then         ! temperature factor spot 1, star 2
              var(i)=spot2parm(1,1)
            endif
c
            if(kk.eq.84)then         ! temperature factor spot 2, star 2
              var(i)=spot2parm(2,1)
            endif
c
            if(kk.eq.48)then         ! temperature factor spot 1, star 1
              var(i)=spot1parm(1,1)
            endif
c
            if(kk.eq.52)then         ! temperature factor spot 2, star 1
              var(i)=spot1parm(2,1)
            endif
c
            if(kk.eq.81)then         ! latitude spot 1, star 2
              var(i)=spot2parm(1,2)
            endif
c
            if(kk.eq.85)then         ! latitude spot 2, star 2
              var(i)=spot2parm(2,2)
            endif
c
            if(kk.eq.49)then         ! latitude spot 1, star 1
              var(i)=spot1parm(1,2)
            endif
c
            if(kk.eq.53)then         ! latitude spot 2, star 1
              var(i)=spot1parm(2,2)
            endif
c
            if(kk.eq.82)then         ! longitude spot 1, star 2
              var(i)=spot2parm(1,3)
            endif
c
            if(kk.eq.86)then         ! longitude spot 2, star 2
              var(i)=spot2parm(2,3)
            endif
c
            if(kk.eq.50)then         ! longitude spot 1, star 1
              var(i)=spot1parm(1,3)
            endif
c
            if(kk.eq.54)then         ! longitude spot 2, star 1
              var(i)=spot1parm(2,3)
            endif
c
            if(kk.eq.83)then         ! radius spot 1, star 2
              var(i)=spot2parm(1,4)
            endif
c
            if(kk.eq.87)then         ! radius spot 2, star 2
              var(i)=spot2parm(2,4)
            endif
c
            if(kk.eq.51)then         ! radius spot 1, star 1
              var(i)=spot1parm(1,4)
            endif
c
            if(kk.eq.55)then         ! radius spot 2, star 1
              var(i)=spot1parm(2,4)
            endif

            if(kk.eq.498)then
              var(i)=pshift
            endif
c
            if(kk.eq.269)then
              var(i)=finc
            endif
c
            if(kk.eq.384)then
              var(i)=Q
            endif
c
            if(kk.eq.130)then
              var(i)=ecc
            endif
c
            if(kk.eq.17)then
              var(i)=argper
            endif
c
            if(kk.eq.176)then
              var(i)=fill1
            endif
c
            if(kk.eq.552)then
              var(i)=rinner
            endif
c
            if(kk.eq.177)then
              var(i)=fill2
            endif
c
            if(kk.eq.464)then
              var(i)=omega1
            endif
c
            if(kk.eq.465)then
              var(i)=omega2
            endif
c
            if(kk.eq.558)then
              var(i)=router
            endif
c
            if(kk.eq.611)then
              var(i)=Tdisk
            endif
c
            if(kk.eq.36)then
              var(i)=betarim
            endif
c
            if(kk.eq.624)then
              var(i)=Teff1
            endif
c
            if(kk.eq.625)then
              var(i)=Teff2
            endif
c
            if(kk.eq.744)then
              var(i)=xi
            endif
c
            if(kk.eq.375)then
              var(i)=rLx
            endif
c
            if(kk.eq.580)then
              var(i)=separ
            endif
c
            if(kk.eq.192)then
              var(i)=gamma
            endif
c
            if(kk.eq.626)then
              var(i)=t3
            endif
c
            if(kk.eq.210)then
              var(i)=g3
            endif
c
            if(kk.eq.576)then
              var(i)=SA3
            endif
c
            if(kk.eq.752)var(i)=dwavex(1,1)
            if(kk.eq.753)var(i)=dwavex(2,1)
            if(kk.eq.754)var(i)=dwavex(3,1)
            if(kk.eq.755)var(i)=dwavex(4,1)
            if(kk.eq.756)var(i)=dwavex(5,1)
            if(kk.eq.757)var(i)=dwavex(6,1)
            if(kk.eq.758)var(i)=dwavex(7,1)
            if(kk.eq.759)var(i)=dwavex(8,1)
c
            if(kk.eq.784)var(i)=dwavey(1,1)
            if(kk.eq.785)var(i)=dwavey(2,1)
            if(kk.eq.786)var(i)=dwavey(3,1)
            if(kk.eq.787)var(i)=dwavey(4,1)
            if(kk.eq.788)var(i)=dwavey(5,1)
            if(kk.eq.789)var(i)=dwavey(6,1)
            if(kk.eq.790)var(i)=dwavey(7,1)
            if(kk.eq.791)var(i)=dwavey(8,1)
c
            if(kk.eq.816)var(i)=dwavex(1,2)
            if(kk.eq.817)var(i)=dwavex(2,2)
            if(kk.eq.818)var(i)=dwavex(3,2)
            if(kk.eq.819)var(i)=dwavex(4,2)
            if(kk.eq.820)var(i)=dwavex(5,2)
            if(kk.eq.821)var(i)=dwavex(6,2)
            if(kk.eq.822)var(i)=dwavex(7,2)
            if(kk.eq.823)var(i)=dwavex(8,2)
c
            if(kk.eq.720)var(i)=dwavey(1,2)
            if(kk.eq.721)var(i)=dwavey(2,2)
            if(kk.eq.722)var(i)=dwavey(3,2)
            if(kk.eq.723)var(i)=dwavey(4,2)
            if(kk.eq.724)var(i)=dwavey(5,2)
            if(kk.eq.725)var(i)=dwavey(6,2)
            if(kk.eq.726)var(i)=dwavey(7,2)
            if(kk.eq.727)var(i)=dwavey(8,2)
c
 10       continue
c
 999      format(7x,i4,2x,a40,3x,f15.8)
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&
c
c  UPDATE October 12, 2001
c
c  Here is a new subroutine that will write the chi^2 values in a
c  single line
c
c  UPDATE FEBRUARY 4, 2005
c
c  Add ochilr to the end of the list.
c

          subroutine chiline(iter,chiall,ochidisk,
     #       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %       icnRV1,icnRV2,chiU,chiB,chiV,chiR,chiI,chiJ,chiH,chiK,
     &       chiRV1,chiRV2,Nobv,sobv,obv,eobv,obsparm,ochilr)
c
c   November 15, 1999
c
c   This routine will write chi^2 values to the screen in a compact way.
c
          implicit double precision (a-h,o-z)
c
c   UPDATE September 21, 2008
c
c   make the dimension of obsparm, obv, eobv, sobv 11
c
c   UPDATE October 10, 2008
c
c   Make the dimension of obsparm, obv, eobv, sobv 17
c
c   UPDATE October 27, 2008
c
c   ochidisk is now an array with dimension 8.
c
c
          dimension ochidisk(8)
          dimension obv(17),eobv(17),obsparm(17)
          character*40 sobv(17)
          character*13 section(20)
          character*80 line
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
          iline=1
          iwrite=0
c
          do 10 i=1,10
            section(i)='1234567890123456'
            section(i)='                '
 10       continue
          line='  -->'
c
c    UPDATE Oct 28, 2002
c
c    Add d0 to the end of the 999999 strings below.
c
c
          cccc=chiall
          if(cccc.gt.9999999.0d0)cccc=9999999.0d0
          write(section(1),100)cccc
          
          if(icnU.ne.430)then
            iline=iline+1
            cccc=chiU
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),100)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnB.ne.430)then
            iline=iline+1
            cccc=chiB
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),101)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnV.ne.430)then
            iline=iline+1
            cccc=chiV
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),102)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnR.ne.430)then
            iline=iline+1
            cccc=chiR
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),103)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnI.ne.430)then
            iline=iline+1
            cccc=chiI
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),104)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnJ.ne.430)then
            iline=iline+1
            iwrite=0
            cccc=chiJ
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),105)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnH.ne.430)then
            iline=iline+1
            iwrite=0
            cccc=chiH
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),106)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnK.ne.430)then
            iline=iline+1
            iwrite=0
            cccc=chiK
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),107)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnRV1.ne.430)then
            iline=iline+1
            iwrite=0
            cccc=chiRV1
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),108)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnRV2.ne.430)then
            iline=iline+1
            iwrite=0
            cccc=chiRV2
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),109)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
c
          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
c
c   NEW BUG  July 24, 2001
c
c   Escape if icn=104 (disk fraction requested).
c
c            if(icn.eq.104)then
c              iline=iline+1
c              cccc=ochidisk(1)
c              if(cccc.gt.9999999.d0)cccc=9999999.0d0
c              write(section(iline),100)cccc
c              go to 1001
c            endif
c
            if(icn.eq.112)then
              iline=iline+1
              cccc=ochidisk(1)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.113)then
              iline=iline+1
              cccc=ochidisk(2)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.114)then
              iline=iline+1
              cccc=ochidisk(3)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.115)then
              iline=iline+1
              cccc=ochidisk(4)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.116)then
              iline=iline+1
              cccc=ochidisk(5)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.117)then
              iline=iline+1
              cccc=ochidisk(6)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.118)then
              iline=iline+1
              cccc=ochidisk(7)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.119)then
              iline=iline+1
              cccc=ochidisk(8)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif

c
c   Escape if icn=369 (luminosity ratio requested).
c
            if(icn.eq.369)then
              iline=iline+1
              cccc=ochilr
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
c   END BUG
c
c   NEW BUG August 2, 2001
c 
c   Put an escape if one of the variables is not valid.
c
            index=-1
            if(icn.eq.400)index=1
            if(icn.eq.401)index=5
            if(icn.eq.560)index=2
            if(icn.eq.561)index=6
            if(icn.eq.208)index=3
c  
            if(icn.eq.336)index=10
            if(icn.eq.337)index=11
c
c
c   UPDATE October 10, 2008
c
c   Here are the strings for incl, mass, ecc, arg, t1, t2
c
            if(icn.eq.269)index=12
            if(icn.eq.384)index=13
            if(icn.eq.130)index=14
            if(icn.eq.17)index=15
            if(icn.eq.624)index=16
            if(icn.eq.625)index=17
c
c   UPDATE May 22, 2002
c
c   Change the 208 to 209 in the argument of if.
c
            if(icn.eq.209)index=7
            if(icn.eq.688)index=4
            if(icn.eq.689)index=8
c
c   UPDATE September 11, 2001
c
c   Add this if block
c
            if(icn.eq.740)index=9
c
c   Here is the escape.
c
            if(index.lt.0)go to 1001
c
c   UPDATE June 3, 2002
c
c   There are two changes.  First, if the variable rmed > 1,
c   then we are doing median fitting.  chisq becomes the absolute
c   deviation. 
c   Second, if the observed duration of the X-ray eclipse is negative, then
c   the number is meant as an upper limit.  That is, if obv = -10.0,
c   then the eclipse duration is less than 10 degrees.
c
            if(index.eq.9)then
              if(obv(i).lt.0.0d0)then
c
c   Here is the case where the upper limit is less than the computed
c   eclipse duration.
c
                if(dabs(obv(i)).gt.obsparm(index))then
                  chisq=0.0d0
                  iline=iline+1
                  go to 1234
                else
                  if(rmed.ge.1.0d0)THEN
                    chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
                  else
                    chisq=(dabs(obv(i))-obsparm(index))*
     #                (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
                  endif
                  iline=iline+1
                  go to 1234
                endif  
              else
                if(rmed.ge.1.0d0)THEN
                  chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
                else
                  chisq=(dabs(obv(i))-obsparm(index))*
     #                (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
                endif
                iline=iline+1
                go to 1234
              endif
            endif
c
            if(rmed.ge.1.0d0)then
              chisq=dabs((obv(i)-obsparm(index))/(eobv(i)))
            else
              chisq=(obv(i)-obsparm(index))*(obv(i)-obsparm(index))
     $          /(eobv(i)*eobv(i))
            endif
            iline=iline+1
c
c   UPDATE October 28, 2002
c
c   Add 0.0d0 to the end of the 9999999 numbers
c
 1234       if(chisq.gt.9999999.0d0)chisq=9999999.0d0
            if(index.eq.1)write(section(iline),100)chisq
            if(index.eq.2)write(section(iline),100)chisq
            if(index.eq.3)write(section(iline),100)chisq
            if(index.eq.4)write(section(iline),100)chisq
            if(index.eq.5)write(section(iline),100)chisq
            if(index.eq.6)write(section(iline),100)chisq
            if(index.eq.7)write(section(iline),100)chisq
            if(index.eq.8)write(section(iline),100)chisq
c
            if(index.eq.10)write(section(iline),100)chisq
            if(index.eq.11)write(section(iline),100)chisq
            if(index.eq.12)write(section(iline),100)chisq
            if(index.eq.13)write(section(iline),100)chisq
            if(index.eq.14)write(section(iline),100)chisq
            if(index.eq.15)write(section(iline),100)chisq
            if(index.eq.16)write(section(iline),100)chisq
            if(index.eq.17)write(section(iline),100)chisq
c
c   UPDATE September 11, 2001
c
c   Add this if block
c
            if(index.eq.9)write(section(iline),100)chisq

            K=lnblnk(line)
            line=line(1:K)//section(iline)
 1001     continue
c          
          write(55,300)iter,(section(j),j=1,iline)
c
c          if(iline.gt.0.and.iwrite.eq.0)write(*,200)line
c
 99       format(i4,' ')
 100      format(' ',f12.4)
 101      format(' ',f12.4)
 102      format(' ',f12.4)
 103      format(' ',f12.4)
 104      format(' ',f12.4)
 105      format(' ',f12.4)
 106      format(' ',f12.4)
 107      format(' ',f12.4)
 108      format(' ',f12.4)
 109      format(' ',f12.4)
 200      format(a80)
 300      format(i6,1x,20(a13))
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
c  UPDATE FEBRUARY 4, 2005
c
c  Here is a new subroutine that will compute the chi^2 of the luminosity
c  fraction.
c
         subroutine lrchi(Nphase,Nmaxphase,ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
c
          implicit double precision (a-h,o-z)
c
c   UPDATE September 21, 2008
c
c   Change the dimensions of obv, eobv, and sobv to 11
c
          dimension ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      ymods3(Nmaxphase)
          dimension obv(11),eobv(11)
          dimension xr(200000),rat(200000)
          character*40 sobv(11)
c
         common /medblock/ rmed 
c
          icount=0
          do 10 i=1,Nphase           
            if(ymods1(i).gt.0.0d0)then
              icount=icount+1
              rat(icount)=ymods2(i)/ymods1(i)
              xr(icount)=ymods3(i)
            endif
 10       continue
c
          if(icount.le.2)then
            ochilr=9999999999.0d0
            ochi=ochi+ochilr
            return
          endif

          call sort2(icount,rat,xr)
c
          q1=dble(icount/2)
          q2=dble(icount)/2.0
          if(q1.eq.q2)then
            rmedian=(rat(icount/2)+rat(icount/2+1))/2.0
          else
            rmedian=rat(icount/2+1)
          endif
          frac=rmedian
c  
          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
c
c  UPDATE May 27, 2002
c
c  If the variable rmed > 1, then define ochidisk to be the absolute
c  deviation, rather than the normal chi^2.
c
            if(icn.eq.369)then        ! string lr = 369
              if(rmed.ge.1.0d0)then
                ochilr=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochilr=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
            endif
 1001     continue
c
          ochi=ochi+ochilr

          return
          end
c
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c  UPDATE FEBRUARY 4, 2005
c
c  Here is a new subroutine that will compute the chi^2 of the luminosity
c  fraction.
c
          subroutine wlrobschi(ochilr)
c
c  July 24, 2000
c
c  This routine will take the observed variables specified in the
c  gridloop.opt file and compute the additional chi^2
c
          implicit double precision (a-h,o-z)
c
c  UPDATE September 21, 2008
c
c  change the dimension of obv, eobv, and sobv to 11
c
          dimension obv(11),eobv(11)
          character*40 sobv(11)
c
          character*26 section
          character*80 line
c
         common /medblock/ rmed 
c
          iline=0
          iwrite=0
c
          section='1234567890123456'
          section='                '
          line='  -->'
c
c   If rmed > 1, then we are doing median fitting.  In that case, write
c   to a different format statement.
c
           if(rmed.ge.1.0d0)then
             write(section,2)ochilr
           else
             write(section,1)ochilr
           endif
           K=lnblnk(line)
           line=line(1:K)//section
           write(*,200)line
c
 1        format(' chi(lum_rat)=',f11.6)
c
c   UPDATE May 27, 2002
c
c   Add this new format statement for median fitting.
c
 2        format(' med(lum_rat)=',f11.6)
 200      format(a80)

          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine initNdata(NdataU,NdataB,NdataV,NdataR,
     &        NdataI,NdataJ,NdataH,NdataK,
     %        NRV1,NRV2,Nobv)
c
          NdataU=0
          NdataB=0
          NdataV=0
          NdataR=0
          NdataI=0
          NdataJ=0
          NdataH=0
          NdataK=0
          NRV1=0
          NRV2=0
          Nobv=0
c
          return
          end
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printchi(chi1)
c
          implicit double precision(a-h,o-z)
c
          if(chi1.gt.9.999999999999999D20)then
            write(*,100)
            return
          endif
          if(chi1.gt.9.999999999999999D19)then
            write(*,101)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D18)then
            write(*,102)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D17)then
            write(*,103)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D16)then
            write(*,104)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D15)then
            write(*,105)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D14)then
            write(*,106)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D13)then
            write(*,107)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D12)then
            write(*,108)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D11)then
            write(*,109)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D10)then
            write(*,110)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D9)then
            write(*,111)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D8)then
            write(*,112)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D7)then
            write(*,113)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D6)then
            write(*,114)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D5)then
            write(*,115)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D4)then
            write(*,116)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D3)then
            write(*,117)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D2)then
            write(*,118)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D1)then
            write(*,119)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D0)then
            write(*,120)chi1
            return
          endif
          write(*,121)chi1

100       format('chi1 = 999999999999999999999.999999')
101       format('chi1 = ',f28.6)
102       format('chi1 = ',f27.6)
103       format('chi1 = ',f26.6)
104       format('chi1 = ',f25.6)
105       format('chi1 = ',f24.6)
106       format('chi1 = ',f23.6)
107       format('chi1 = ',f22.6)
108       format('chi1 = ',f21.6)
109       format('chi1 = ',f20.6)
110       format('chi1 = ',f19.6)
111       format('chi1 = ',f18.6)
112       format('chi1 = ',f17.6)
113       format('chi1 = ',f16.6)
114       format('chi1 = ',f15.6)
115       format('chi1 = ',f14.6)
116       format('chi1 = ',f13.6)
117       format('chi1 = ',f12.6)
118       format('chi1 = ',f11.6)
119       format('chi1 = ',f10.6)
120       format('chi1 = ',f9.6)
121       format('chi1 = ',f8.6)

c
          return
          end
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printmed(chi1)
c
          implicit double precision(a-h,o-z)
c
          if(chi1.gt.9.999999999999999D20)then
            write(*,100)
            return
          endif
          if(chi1.gt.9.999999999999999D19)then
            write(*,101)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D18)then
            write(*,102)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D17)then
            write(*,103)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D16)then
            write(*,104)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D15)then
            write(*,105)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D14)then
            write(*,106)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D13)then
            write(*,107)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D12)then
            write(*,108)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D11)then
            write(*,109)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D10)then
            write(*,110)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D9)then
            write(*,111)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D8)then
            write(*,112)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D7)then
            write(*,113)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D6)then
            write(*,114)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D5)then
            write(*,115)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D4)then
            write(*,116)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D3)then
            write(*,117)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D2)then
            write(*,118)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D1)then
            write(*,119)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D0)then
            write(*,120)chi1
            return
          endif
          write(*,121)chi1

100       format('med1 = 999999999999999999999.999999')
101       format('med1 = ',f28.6)
102       format('med1 = ',f27.6)
103       format('med1 = ',f26.6)
104       format('med1 = ',f25.6)
105       format('med1 = ',f24.6)
106       format('med1 = ',f23.6)
107       format('med1 = ',f22.6)
108       format('med1 = ',f21.6)
109       format('med1 = ',f20.6)
110       format('med1 = ',f19.6)
111       format('med1 = ',f18.6)
112       format('med1 = ',f17.6)
113       format('med1 = ',f16.6)
114       format('med1 = ',f15.6)
115       format('med1 = ',f14.6)
116       format('med1 = ',f13.6)
117       format('med1 = ',f12.6)
118       format('med1 = ',f11.6)
119       format('med1 = ',f10.6)
120       format('med1 = ',f9.6)
121       format('med1 = ',f8.6)

c
          return
          end
c
c  ##########
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printsmall(chi1)
c
          implicit double precision(a-h,o-z)
c
          if(chi1.gt.9.999999999999999D20)then
            write(*,100)
            return
          endif
          if(chi1.gt.9.999999999999999D19)then
            write(*,101)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D18)then
            write(*,102)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D17)then
            write(*,103)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D16)then
            write(*,104)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D15)then
            write(*,105)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D14)then
            write(*,106)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D13)then
            write(*,107)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D12)then
            write(*,108)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D11)then
            write(*,109)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D10)then
            write(*,110)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D9)then
            write(*,111)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D8)then
            write(*,112)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D7)then
            write(*,113)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D6)then
            write(*,114)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D5)then
            write(*,115)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D4)then
            write(*,116)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D3)then
            write(*,117)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D2)then
            write(*,118)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D1)then
            write(*,119)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D0)then
            write(*,120)chi1
            return
          endif
          write(*,121)chi1

100       format('chi_small = 99999999999999999999.999999')
101       format('chi_small = ',f28.6)
102       format('chi_small = ',f27.6)
103       format('chi_small = ',f26.6)
104       format('chi_small = ',f25.6)
105       format('chi_small = ',f24.6)
106       format('chi_small = ',f23.6)
107       format('chi_small = ',f22.6)
108       format('chi_small = ',f21.6)
109       format('chi_small = ',f20.6)
110       format('chi_small = ',f19.6)
111       format('chi_small = ',f18.6)
112       format('chi_small = ',f17.6)
113       format('chi_small = ',f16.6)
114       format('chi_small = ',f15.6)
115       format('chi_small = ',f14.6)
116       format('chi_small = ',f13.6)
117       format('chi_small = ',f12.6)
118       format('chi_small = ',f11.6)
119       format('chi_small = ',f10.6)
120       format('chi_small = ',f9.6)
121       format('chi_small = ',f8.6)

c
          return
          end
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printsmallmed(chi1)
c
          implicit double precision(a-h,o-z)
c
          if(chi1.gt.9.999999999999999D20)then
            write(*,100)
            return
          endif
          if(chi1.gt.9.999999999999999D19)then
            write(*,101)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D18)then
            write(*,102)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D17)then
            write(*,103)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D16)then
            write(*,104)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D15)then
            write(*,105)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D14)then
            write(*,106)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D13)then
            write(*,107)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D12)then
            write(*,108)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D11)then
            write(*,109)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D10)then
            write(*,110)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D9)then
            write(*,111)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D8)then
            write(*,112)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D7)then
            write(*,113)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D6)then
            write(*,114)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D5)then
            write(*,115)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D4)then
            write(*,116)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D3)then
            write(*,117)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D2)then
            write(*,118)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D1)then
            write(*,119)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D0)then
            write(*,120)chi1
            return
          endif
          write(*,121)chi1

100       format('med_small = 99999999999999999999.999999')
101       format('med_small = ',f28.6)
102       format('med_small = ',f27.6)
103       format('med_small = ',f26.6)
104       format('med_small = ',f25.6)
105       format('med_small = ',f24.6)
106       format('med_small = ',f23.6)
107       format('med_small = ',f22.6)
108       format('med_small = ',f21.6)
109       format('med_small = ',f20.6)
110       format('med_small = ',f19.6)
111       format('med_small = ',f18.6)
112       format('med_small = ',f17.6)
113       format('med_small = ',f16.6)
114       format('med_small = ',f15.6)
115       format('med_small = ',f14.6)
116       format('med_small = ',f13.6)
117       format('med_small = ',f12.6)
118       format('med_small = ',f11.6)
119       format('med_small = ',f10.6)
120       format('med_small = ',f9.6)
121       format('med_small = ',f8.6)

c
          return
          end




