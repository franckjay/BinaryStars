         program findgridminELC
c
c    April 24, 2000
c
c    This program is similar to mkgridcontELC, except that the program
c    searches for the minimum chi^2 on an edge and fills in the next
c    point.  Thus this program will tend to find the minimum of the grid
c    much faster.
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
c
c   NEW BUG August 2, 2001
c
c   Add the period and T0 to the assignvar argument list, and change 'sw4' to
c   'T0'
c
          open(unit=30,file='cont.inp',status='unknown')
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
            open(unit=30,file=donefile,status='unknown')
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
c
          do 1000 igrow=1,N1*N2
                call newsetparm(i,j,idone,var1val,var2val,N1,N2,svar)
                write(command,100)i,j                
                write(*,*)command
                if(i.gt.500+N1)go to 1000
                if(j.gt.500+N2)go to 1000
                call system(command)
                write(command,101)i,j
                write(*,*)command
                call system(command)
                write(command,102)i,j
                write(*,*)command
                call system(command)
                write(command,103)i,j
                write(*,*)command
                call system(command)
                write(command,104)i,j
                write(*,*)command
                call system(command)
                write(command,105)i,j
                write(*,200)command
                call system(command) 
                write(command,106)i,j
                write(*,200)command
                call system(command)
                write(command,107)i,j
                write(*,200)command
                call system(command)
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
 107      format('cp ELCratio.0000 ELCratio.',i3,'.',i3)
 200      format(a40)
 300      format(f11.6)
c
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
 10       continue
          return
          end
c
c   **************************************************************
c          
          subroutine newsetparm(i,j,idone,var1val,var2val,N1,N2,svar)
c
c   given the k,l index of the grid point, check all of the neighbor points
c   that have been done and get the values of the parameters that has the
c   lowest chi^2
c
          implicit double precision (a-h,o-z)
          dimension idone(500,500),var1val(N1),var2val(N2)
          dimension var(2),powercoeff(8,9)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          character*40 test(12),command,svar(2)
          character*38 line
          character*80 kommand
          logical edge,up,down,left,right
c
c   Loop over the grid and look for edge points
c

          Smin=1000222.0d0
          do 2 k=501,N1+500
            do 1 l=501,N2+500
              edge=.false.
              left=.false.
              right=.false.
              up=.false.
              down=.false.
              if(idone(k-500,l-500).eq.1)then      !done point
                if(k.eq.501)left=.true.
                if(k.eq.N1+500)right=.true.
                if(l.eq.501)down=.true.
                if(l.eq.N2+500)up=.true.
                if(idone(k-500-1,l-500).eq.1)left=.true.
                if(idone(k-500+1,l-500).eq.1)right=.true.
                if(idone(k-500,l-500-1).eq.1)down=.true.
                if(idone(k-500,l-500+1).eq.1)up=.true.
                if(up.eqv..false.)edge=.true.
                if(down.eqv..false.)edge=.true.
                if(left.eqv..false.)edge=.true.
                if(right.eqv..false.)edge=.true.
                if(edge.eqv..true.)then
                  command='rm -f Sfile'
                  call system(command)
                  write(command,100)k,l
                  write(*,*)command
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
              endif
 1          continue
 2        continue
c
          write(command,102)isave,jsave
          write(*,*)command
          call system(command)
c
c   Now set the variables to the grid point nearest the minimum edge
c
          if((isave.gt.501).and.(idone(isave-500-1,jsave-500).eq.0))then
            i=isave-1
            j=jsave
            go to 3
          endif
c
          if(idone(isave-500+1,jsave-500).eq.0)then
            i=isave+1
            j=jsave
            go to 3
          endif
c
          if((jsave.gt.501).and.(idone(isave-500,jsave-500-1).eq.0))then
            i=isave
            j=jsave-1
            go to 3
          endif
c
          if(idone(isave-500,jsave-500+1).eq.0)then
            i=isave
            j=jsave+1
            go to 3
          endif
c
c   Assign the variables and write a new ELC.inp file.
c
           
           command='cat ELC.inp'
           write(*,*)command
           call system(command)
c
3          command='cat ELC.inp'
           write(*,*)command
           call system(command)

           call getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $        omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
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
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,sw32,sw33,sw34,
     %       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c
         write(*,*)bigI,bigbeta
         write(*,*)'test  ',isw25,isw26,isw27,isw28
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
c
          var(1)=var1val(i-500)
          var(2)=var2val(j-500)
c
          write(*,*)'var(1), var(2) ',var(1),var(2),i,j
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
           call assignvar(2,svar,var,fill1,fill2,omega1,
     $          omega2,Q,finc,Teff1,Teff2,betarim,
     $         rinner,router,tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     %         pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     &         dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     %         ecosw,temprat,bigI,bigbeta,powercoeff,density)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
          write(*,*)bigI,bigbeta
          write(*,*)'isw25=',isw25,isw26,isw27,isw28,isw29
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
     %       sw25,sw26,sw27,sw28,sw29,sw30,sw31,sw32,sw33,sw34,
     #       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34)
c


          command='cp gridELC.inp ELC.inp'
          write(*,*)command
          call system(command)
c
          idone(i-500,j-500)=1
c
          if(i.gt.500+N1)return
          if(j.gt.500+N2)return
c
          line='  '
          write(kommand,300)line,i,j
          write(*,*)kommand
          call system(kommand)
          
 100      format('tail -1 chi.',i3,'.',i3,'  | grep S > Sfile')
 101      format(1x,4x,f15.6)
 102      format('cp gridELC.',i3,'.',i3,' gridloop.opt')
 200      format(a40)
 300      format('echo "',a38,'" > chi.',i3,'.',i3)
          return
          end          
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
c
