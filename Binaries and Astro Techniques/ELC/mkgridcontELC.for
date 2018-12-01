          program makecontourELC
c
c    November 16, 1999
c
c    This program will drive gridELC in order to make a contour
c    plot of chi^2 in some two-parameter plane (for example Q vs i).
c    The user specifies a list of other variables to be optimized at
c    each point---the values of var1,var2 are fixed at the values
c    corresponding to the grid location.
c
          implicit double precision (a-h,o-z)

          parameter(Maxval=500)
          dimension var1val(Maxval),var2val(Maxval),idone(Maxval,Maxval),
     %             var(2)
          character*40 command,donefile,svar(2)
c
c    Open the input file.  The first line is the name of the parameter for
c    x-axis variable.  The second line is the name of the parameter for the
c    y-axis variable.  The names of the variables are the same as the names
c    in the gridloop.opt file. The third line has three entries: starting value
c    for x, the stepsize, and the Number of steps.  The fourth line is
c    similar to the third, but for the y-axis.  The fifth line has two integer
c    entries which are the index of the starting point for the x array and
c    the index for the starting point of the y-array.  The sixth line is a 
c    flag:  negative values mean we should start from scratch while positive
c    values mean there are grid points already filled in.  In the latter case
c    we have to give the name of the file on line 7.  This file should have
c    two columns of integers.  These numbers are the i,j index values of the
c    previously finished points.
c 
c    For example, suppose
c    var1=Q and var2=finc.  The solution was optimized at Q=5 and i=60.  We
c    want to make a contour map with 1 <= Q <= 10 in steps of 1 and
c    50 <= i <= 70 in steps of 1 degree.  This is the first run.
c    The input file will look like this:
c
c    mass ratio
c    inclination
c    1.0  1.0  10
c    50.0  1.0  21
c    5  11
c    -1
c
c    
c    At each point on the grid, the axis variables will be set and fixed.
c    gridELC will be called to optimized a list of other variables.  We
c    start the iteration at the central value, using the list of variables
c    as the initial guess.  We then work our way outwards, using the variables
c    from the nearest best solution as the initial guess.    At each point,
c    the output from gridELC is saved in a file called chi.III.JJJ   Here
c    III starts at 500 and goes to 500 + Nx.  Similarly, JJJ starts at 500
c    and goes to 500 + Ny.  
c
c    THE USER MUST RUN gridELC AT THE CENTRAL STARTING POINT FIRST, SAVING
c    THE SCREEN OUTPUT TO A FILE.  Then
c    rename the screen output as chi.AAA.BBB, where AAA=ixmin+500 and 
c    BBB=iymin+500.  In the above example, AAA=505 and BBB=511.  Also,
c    rename the gridELC.opt file to gridELC.AAA.BBB.
c
c    When everything is done, we extract the chi^2 values from the chi.III.JJJ
c    files and make the contour map.
c
          open(unit=30,file='cont.inp',status='old')
          read(30,200)svar(1)
          read(30,200)svar(2)
          read(30,*)start1,step1,N1
          read(30,*)start2,step2,N2
          read(30,*)imin,jmin
          iflag=-100
          read(30,*)iflag
c
c   If iflag is greater than zero, we have to read the name of the file
c   containing the index values of the finished points.  Open that file
c   and get the list of index values.
c
          if(iflag.gt.0)then
            read(30,200)donefile
            close(30)
            open(unit=30,file=donefile,status='old')
            do 5 i=1,100000
              read(30,*,end=6)j,k
              idone(j-500,k-500)=1
 5          continue
 6          close(30)
          else
            close(30)
          endif
c
c   It is always assumed that the central point is finished.
c
          idone(imin,jmin)=1
c 
c   Fill the arrays with the parameter values.
c 
          call setval(start1,step1,N1,var1val)
          call setval(start2,step2,N2,var2val)
c
          imin=imin+500
          jmin=jmin+500
          if(iflag.gt.0)go to 11    !skip this if things have been done before
          do 10 i=imin-1,imin+1
            do 9 j=jmin-1,jmin+1
              if(i.eq.imin.and.j.eq.jmin)then
                go to 9
              else
                call setparm(i,j,idone,var1val,var2val,N1,N2,svar)
                write(command,100)i,j
                write(*,200)command
                call system(command)
                idone(i-500,j-500)=1
                write(command,101)i,j
                write(*,200)command
                call system(command)
                write(command,102)i,j
                write(*,200)command
                call system(command)
                write(command,103)i,j
                write(*,200)command
                call system(command)
                write(command,104)i,j
                write(*,200)command
                call system(command)
                write(command,105)i,j
                write(*,200)command
                call system(command)
                write(command,106)i,j
                write(*,200)command
                call system(command)
              endif
 9          continue
 10       continue
c
 11       continue
c
c   We are done with the first pass.  Now increase the grow radius and
c   loop again.  However, we must adjust 'gridloop.opt' based on the
c   best fit of a nearby grid point.
c
          do 1000 igrow=2,max(N1,N2)
c 
            do 999 i=max(501,imin-igrow),min(500+N1,igrow+imin)
              do 998 j=max(501,jmin-igrow),min(500+N2,igrow+jmin)
                if(idone(i-500,j-500).eq.1)go to 998
                call setparm(i,j,idone,var1val,var2val,N1,N2,svar)
                write(command,100)i,j                
                write(*,200)command
                call system(command)
                write(command,101)i,j
                write(*,200)command
                call system(command)
                write(command,102)i,j
                write(*,200)command
                call system(command)
                write(command,103)i,j
                write(*,200)command
                call system(command)
                write(command,104)i,j
                write(*,200)command
                call system(command)
                write(command,105)i,j
                write(*,200)command
                call system(command)
                write(command,106)i,j
                write(*,200)command
                call system(command)
 998          continue
 999        continue
 1000     continue
c
          open(unit=32,file='makecontour.x',status='unknown')
          do 500 i=1,N1
            write(32,300)var1val(i)
 500      continue
          close(32)
c
          open(unit=32,file='makecontour.y',status='unknown')
          do 501 i=1,N2
            write(32,300)var2val(i)
 501      continue
          close(32)
c


 100      format('gridELC >> chi.',i3,'.',i3)
 101      format('cp gridELC.opt gridELC.',i3,'.',i3)
 102      format('cp gridELC.inp ELC.',i3,'.',i3)
 103      format('cp ELC.phases phase.',i3,'.',i3)
 104      format('cp ELC.parm parm.',i3,'.',i3)
 105      format('cp ELCparm.0000 ELCparm.',i3,'.',i3)
 106      format('cp generation.0000 generation.',i3,'.',i3)
 200      format(a40)
 300      format(f11.6)
c
          end
c
c   **************************************************************
c          
          subroutine setparm(i,j,idone,var1val,var2val,N1,N2,svar)
c
c   given the k,l index of the grid point, check all of the neighbor points
c   that have been done and get the values of the parameters that has the
c   lowest chi^2
c
          implicit double precision (a-h,o-z)

          dimension idone(500,500),var1val(N1),var2val(N2)
          dimension var(2)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
c
c   Dimension the spot arrays.
c
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          character*40 test(12),command,svar(2)
          character*38 line
          character*80 kommand
c
c   First, assign the variables and write a new ELC.inp file.
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
     %       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,sw12,sw13,
     &       idark1,idark2,isw12,isw13,
     %       isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff)
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
c
          var(1)=var1val(i-500)
          var(2)=var2val(j-500)
          write(*,*)var(1),var(2)
c
             call assignvar(2,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $         rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     %         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &         dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     &         bigI,bigbeta,powercoeff,density)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
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
     #       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,sw12,sw13,
     &       idark1,idark2,isw12,isw13)
c
c  Now we have to find the nearby points that have been done and look for
c  the point with the lowest value of S
c
          Smin=1000222.0d0
          do 10 k=max(501,i-1),min(N1+500,i+1)
            do 9 l=max(501,j-1),min(N2+500,j+1)
              if(k.eq.i.and.l.eq.j)go to 9
              if(idone(k-500,l-500).eq.1)then
                command='/bin/rm -f Sfile'
                call system(command)
                write(command,100)k,l
                write(*,200)command
                call system(command)
                open(unit=30,file='Sfile',status='unknown')
                read(30,101)S
                close(30)
                if(S.lt.Smin)then
                  Smin=S
                  isave=k
                  jsave=l
                endif
              endif
 9          continue
 10       continue
c
c   Now we have the chi^2 file we want with index isave,jsave.  Copy
c   the file gridELC.isave.jsave (e.g. gridELC.510.515) to gridloop.opt.
c
          write(command,102)isave,jsave
                write(*,200)command
          call system(command)
c
          idone(i-500,j-500)=1
c
          call writevar(2,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     $       line)
c
          write(kommand,300)line,i,j
          call system(kommand)
          
 100      format('tail -2 chi.',i3,'.',i3,'  | grep S > Sfile')
 101      format(1x,4x,f15.6)
 102      format('cp gridELC.',i3,'.',i3,' gridloop.opt')
 200      format(a40)
 300      format('echo "',a38,'" > chi.',i3,'.',i3)
          return
          end          
c
c *******************************************************
c
          subroutine setval(start,step,N,varval)
c
          implicit double precision (a-h,o-z)

          dimension varval(N)
c
          varval(1)=start
          x=start
          do 10 i=2,N
            x=x+step
            varval(i)=x
            write(*,*)i,varval(i)
 10       continue
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
          subroutine assignvar(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     &       pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,   
     $       dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     *       bigI,bigbeta,powercoeff)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
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

          write(*,*)' '
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
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
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
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
c    *************************************************************
c
          subroutine writevar(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,
     $       rinner,router,tdisk,xi,rLx,separ,isvel1,gamma1,isvel2,gamma2,
     $       line)
c
c   November 15, 1999
c
c   This routine will determine which variables need to be printed on the
c   screen based on the string codes in svar(1:Nvarmax).  The output
c   printed on the screen will be compact.
c
          implicit double precision (a-h,o-z)

          character*38 line
          character*40 svar(Nvarmax)
          character*16 string16
          character*10 string10
          character*11 string11
          character*12 string12
          character*13 string13
          character*14 string14
          character*19 string19
          dimension var(Nvarmax)
c
          line='**'
          iline=0
          iwrite=0
          ilength=2
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
            if(kk.eq.269)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,200)finc
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.384)then
              iline=iline+1
              iwrite=0
              ilength=ilength+10
              write(string10,201)Q
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.176)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,202)fill1
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.177)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,203)fill2
              K=lnblnk(line)
              line=line(1:K)//string13
            endif
c
            if(kk.eq.552)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,204)rinner
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.464)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,205)omega1
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.465)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,206)omega2
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.558)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,207)router
              K=lnblnk(line)
              line=line(1:K)//string13
            endif
c
            if(kk.eq.611)then
              iline=iline+1
              iwrite=0
              ilength=ilength+13
              write(string13,208)Tdisk
              K=lnblnk(line)
              line=line(1:K)//string13
            endif
c
            if(kk.eq.36)then
              iline=iline+1
              iwrite=0
              ilength=ilength+12
              write(string12,211)betarim
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.624)then
              iline=iline+1
              iwrite=0
              ilength=ilength+16
              write(string16,209)Teff1
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.625)then
              iline=iline+1
              iwrite=0
              ilength=ilength+16
              write(string16,210)Teff2
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.744)then
              iline=iline+1
              iwrite=0
              ilength=ilength+11
              write(string11,212)xi
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.377)then
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,213)rLx
              K=lnblnk(line)
              line=line(1:K)//string19
            endif
c
            if(kk.eq.580)then
              iline=iline+1
              iwrite=0
              ilength=ilength+14
              write(string14,214)separ
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
 10       continue
c
c
 200      format(' i=',f7.4)
 201      format(' Q=',f7.4)
 202      format(' fi1=',f7.5)
 203      format(' fi2=',f7.5)
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
 215      format(' gam1=',f7.2)
 216      format(' gam2=',f7.2)
 500      format(a80)
 999      format(7x,i4,2x,a40,3x,f15.8)
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&
c
c
c  Add the spot parameters to the argument list.
c
c  NEW BUG ALERT August 2, 2001
c
c  Replace 'sw4' with T0
c
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
     #       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,sw12,sw13,
     &       idark1,idark2,isw12,isw13)
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
c
c   RVG BUG ALERT  May 9, 2001
c
c   Dimension the spot arrays.
c
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
c 
          open(unit=1,file='ELC.inp',status='unknown')
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
          write(1,2032)sw12
          write(1,2033)sw13
c
          write(1,2040)idark1
          write(1,2041)idark2
          write(1,2042)isw12
          write(1,22043)isw13
c
 2025     format(f13.9,7x,'primmass (star 1 mass in solar masses)')
 2026     format(f14.9,6x,'primK (K-velocity of star 1 in km/sec)')
 2027     format(f14.9,6x,'primrad (star 1 radius in solar radii)')
 2028     format(f16.9,4x,
     &          'ratrad (ratio of star 1 radius and star 2 radius)')
c
 2030     format(f14.12,6x,'frac1 (fractional radius star 1: R_1/a)')
 2031     format(f14.12,6x,'frac2 (fractional radius star 2: R_2/a)')
 2032     format(f4.2,16x,'sw12 (currently unactive)')
 2033     format(f4.2,16x,'sw13 (currently unactive)')
c
 2040     format(i1,19x,'idark1')
 2041     format(i1,19x,'idark2')
 2042     format(i1,19x,'isw12 (currently unactive)')
22043     format(i1,19x,'isw13 (currently unactive)')
c
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
 1008     format(f4.2,16x,'dphase')
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
 3028     format(i1,19x,'ilaw  (1=linear law, 2=logarithmic law,',
     %           ' 3=square root law)')
 4000     format(i1,19x,'iatm')
 4001     format(i1,19x,'ism1')
 4002     format(8(i1,1x),4x,'icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK')
 4025     format(f11.6,9x,'gamma velocity')
 5000     format(f10.2,10x,'t3')
 5001     format(f10.5,10x,'g3')
 5002     format(f12.6,8x,'SA3')
 5003     format(f12.6,8x,'density')
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
 6010     format(i1,19x,'isimp')
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
     %       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,sw12,sw13,
     &       idark1,idark2,isw12,isw13,
     %       isw21,usw22,isw23,isw24,sw21,sw22,sw23,sw24,powercoeff)
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
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
          read(1,*,end=101,err=101)sw12
          read(1,*,end=101,err=101)sw13
          read(1,*,end=101,err=101)idark1
          read(1,*,end=101,err=101)idark2
          read(1,*,end=101,err=101)isw12
          read(1,*,end=101,err=101)isw13
c
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
     $       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,sw12,sw13,
     &       idark1,idark2,isw12,isw13)
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
     $       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,sw12,sw13,
     &       idark1,idark2,isw12,isw13)
c
c    will write the correctly formatted file ELC.inp and return
c    default parameters
c
          implicit double precision(a-h,o-z)
c
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2),
     %       www(8)
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
          density=-0.01d0
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
          sw12=0.0d0
          sw13=0.0d0
c
          idark1=0
          idark2=0
          isw12=0
          isw13=0
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
          write(1,2028)ratrad
c
          write(1,2030)frac1
          write(1,2031)frac2
          write(1,2032)sw12
          write(1,2033)sw13
c
          write(1,2040)idark1
          write(1,2041)idark2
          write(1,2042)isw12
          write(1,3043)isw13

c
 2025     format(f13.9,7x,'primmass (star 1 mass in solar masses)')
 2026     format(f14.9,6x,'primK (K-velocity of star 1 in km/sec)')
 2027     format(f14.9,6x,'primrad (star 1 radius in solar radii)')
 2028     format(f16.9,4x,
     &          'ratrad (ratio of star 1 radius and star 2 radius)')
c
 2030     format(f4.2,16x,'frac1 (fractional radius star 1: R_1/a)')
 2031     format(f4.2,16x,'frac2 (fractional radius star 2: R_2/a)')
 2032     format(f4.2,16x,'sw12 (currently unactive)')
 2033     format(f4.2,16x,'sw13 (currently unactive)')
c
 2040     format(i1,19x,'idark1')
 2041     format(i1,19x,'idark2')
 2042     format(i1,19x,'isw12 (currently unactive)')
 3043     format(i1,19x,'isw13 (currently unactive)')

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
 3028     format(i1,19x,'ilaw  (1=linear law, 2=logarithmic law,',
     %           ' 3=square root law)')
 4025     format(f7.2,13x,'gamma velocity (km/sec)')
 4000     format(i1,19x,'iatm  (0 for BB, 1 for model atmospheres)')
 4001     format(i1,19x,'ism1  (0 for all phases, 1 for 0-180)')
 4002     format(8(i1,1x),4x,'icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK')
 5000     format(f10.2,10x,'t3')
 5001     format(f5.2,15x,'g3')
 5002     format(f12.6,8x,'SA3')
 5003     format(f12.6,8x,'density')
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
 6010     format(i1,19x,'isimp')
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


