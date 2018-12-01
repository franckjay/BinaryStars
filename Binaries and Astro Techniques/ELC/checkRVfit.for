          program checkfit
c
c   This will use the command line arguments.
c
c
c
          dimension xmodel(9500),ymodel(9500),xdata(9500),ydata(9500)
          dimension xinter(9500),yinter(9500),err(9500)
          dimension errplot(9500),y2(9500)
          dimension xdummy(9500),ydummy(9500),xplot(9500),yplot(9500)
          character*40 filein,fileout
          character*130 txt
c
          M=iargc()    ! find out how many arguments are on the command line
c
          call getvariables(M,filein,fileout)
c
          open(unit=20,file=filein)
c
          ymin=1000.
          ymax=-1000.
          do 10 i=1,9500                            ! read in the model here
            read(20,*,end=15)xmodel(i),ymodel(i)
c            ymodel(i)=-2.5*log10(ymodel(i))
 10       continue
 15       Nmodel=i-1
          close(20)
c
          call sort2(Nmodel,xmodel,ymodel)         ! sort by phase
c
          open(unit=23,file=fileout)      ! read in the data
c
          ymin=100000.0
          ymax=-100000.0
          do 20 i=1,9500
            read(23,*,end=25)xdata(i),ydata(i),err(i)
            if(ydata(i).lt.ymin)ymin=ydata(i)-err(i)
            if(ydata(i).gt.ymax)ymax=ydata(i)+err(i)
            xplot(i)=xdata(i)
            yplot(i)=ydata(i)
            errplot(i)=err(i)
 20       continue
 25       Ndata=i-1
          close(23)
c
          call sort3(Ndata,xdata,ydata,err)  !sort by phase
c
c   now we have to find the value of the model at each of the observed
c   phase values
c
          call spline(xmodel,ymodel,Nmodel,0.0,0.0,y2)
c
          do 30 i=1,Ndata
            call splint(xmodel,ymodel,y2,Nmodel,xdata(i),qqq)
            xinter(i)=xdata(i)
            yinter(i)=qqq
 30       continue
c
c    Now we have an array xinter and yinter containing the model resampled
c    at the same x-values as the data.  Perform the chi^2 fit.
c
          yhi=ymax
          ylow=ymin
c
c          call pgbegin(0,'?',1,1)
c          call pgscf(2)
c          call pgsch(1.2)
c          call pgenv(0.0,2.0,yhi+0.02,ylow-0.02,0,0)
c
c          do 2000 i=1,Ndata
c            call pgpoint(1,xplot(i),yplot(i),17)
c            call pgpoint(1,xplot(i)+1.0,yplot(i),17)
c            et=yplot(i)+errplot(i)
c            eb=yplot(i)-errplot(i)
c            call pgerry(1,xplot(i),et,eb,1.0)
c            call pgerry(1,xplot(i)+1.0,et,eb,1.0)
c 2000     continue
c
          call getmean(Ndata,xdata,ydata,dataave)
          call getmean(Ndata,xinter,yinter,rmodelave)
c

          zero=(dataave-rmodelave)
          step=abs(rmodelave-dataave)/1000.0
c
          chi1=0.0
          chi2=0.0
          chi3=0.0
          small=1000000.0
c
          do 750 i=1,20
            zero1=zero
c            write(*,*)chi1,chi2,chi3
            call offset(Ndata,xinter,yinter,ydummy,zero1)
            call getchi(Ndata,ydata,err,ydummy,chi1)
c            write(*,*)'chi1 = ',chi1,zero1
            if(chi1.lt.small)then
              small=chi1
              zerosmall=zero1
            endif
            fn=0.0
            zero2=zero+step
            call offset(Ndata,xinter,yinter,ydummy,zero2)
            call getchi(Ndata,ydata,err,ydummy,chi2)
c            write(*,*)'chi2 = ',chi2,zero2
            if(chi2.lt.small)then
              small=chi2
              zerosmall=zero2
            endif
            diff=chi1-chi2
            if(diff.gt.0.0)go to 5061
            step=-step
c  
            csave=chi1
            chi1=chi2
            chi2=csave
            zsave=zero1
            zero1=zero2
            zero2=zsave
            
c
 5061       fn=fn+1.0
c
            zero3=zero2+step
            call offset(Ndata,xinter,yinter,ydummy,zero3)
            call getchi(Ndata,ydata,err,ydummy,chi3)
c            write(*,*)'chi3 = ',chi3,zero3
            if(chi3.lt.small)then
              small=chi3
              zerosmall=zero3
            endif
            diff23=chi3-chi2
            if(diff23.lt.0)then
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
 5081       if((chi3-chi2).eq.0.0)go to 999
            step=step*(1.0/(1.0+(chi1-chi2)/(chi3-chi2))+0.5)
            zero=zero2-step
            step=step*fn/3.0
            call offset(Ndata,xinter,yinter,ydummy,zero)
            call getchi(Ndata,ydata,err,ydummy,chi4)
c            write(*,*)'chi4 = ',chi4,zero
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
          call offset(Ndata,xinter,yinter,ydummy,zero)
          call getchi(Ndata,ydata,err,ydummy,chi3)
c
          chisq=chi3
c
          open(unit=37,file='residuals')
          do 200 i=1,Ndata
            qqq=xdata(i)
            write(37,300)qqq,(ydata(i)-ydummy(i)),err(i)
 200      continue

          call offset(Nmodel,xmodel,ymodel,ydummy,zero)
c          call pgline(Nmodel,xmodel,ydummy)
c
          open(unit=33,file='model.out')
          do 3333 i=1,Nmodel
            write(33,1003)xmodel(i),ydummy(i)
 3333     continue
          close(33)
c
c
c          call pglabel('Phase','Relative Magnitude',' ')
c
c
 1000     format(1x,'i=',f4.1,1x,'Q=',f4.1,2x,'beta=',f4.1,1x,
     $      'T_d=',f7.1,1x,'L_x=',1pe10.3,1x,'W='0pf4.2,1x,
     $      'r_d=',f5.3,1x,'xi=',f5.2,1x,'chi^2=',f12.4,2x)
c
 1001     format('i=',f4.1,1x,'Q=',f4.1,1x,'\gb=',f4.1,1x,'T\dd\u=',
     $      f7.1,1x,'L\dx\u=',1pe10.3,1x,'W=',0pf4.2,1x,'r\dd\u=',
     $      f5.3,1x,'\gc=',f5.2,2x,'\gx\u2\d=',f13.4)
 1002     format(a130)
 1003     format(1x,2(f15.6,2x))
 300      format(1x,2(f15.9,3x),2x,f10.7)
 1004     format(1x,'chi^2 = ',f15.6,3x,'zero = ',f9.5)
c
          write(*,1004)chisq,zero

c          call pgend

          end
c
c    ===================================================
c
          subroutine getvariables(M,filein,fileout)
c
          character*40 filein,fileout
          integer M
c
          if(M.eq.0)then   !case of no arguments
c
            write(*,*)'enter the file name of the file with the model ',
     $       'in LINEAR units'
            read(*,100)filein
c
            write(*,*)'enter the file name of the file with the FOLDED ',
     $        'data in MAGNITUDES'
            read(*,100)fileout
c
          endif
c
          if(M.eq.1)then        !case of only one argument
c
            call getarg(1,filein)     ! assign filein the value of the entered
                                      ! variable 

            write(*,*)'enter the file name of the file with the FOLDED ',
     $        'data in MAGNITUDES'
            read(*,100)fileout
c
          endif
c
          if(M.eq.2)then        ! case of both variables specified
c
            call getarg(1,filein)
            call getarg(2,fileout)
          endif
c
 100      format(a40)
          return
          end

c
c
c &&&&&&&&&&&&&&&&&&&&
c
      SUBROUTINE SORT2(N,RA,RB)
      DIMENSION RA(N),RB(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END
c
c
c
c
      SUBROUTINE SORT3(N,RA,RB,rc)
      DIMENSION RA(N),RB(N),rc(N)
      L=N/2+1
      IR=N
10    CONTINUE
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
20      IF(J.LE.IR)THEN
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
c  -----------------------------------------------------
c
          subroutine getchi(N,y,err,ydum,chisq)
c
          dimension y(N),ydum(N),err(N)
          chisq=0.0
c
          do 10 i=1,N
            chisq=chisq+(y(i)-ydum(i))*(y(i)-ydum(i))/(err(i)*err(i))
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
c  &&&&&&&&&&&&&&&&&
c
         subroutine spline(x,y,n,yp1,ypn,y2)
c
c   November 12, 1999
c
c   This is a spline interpolation routine taken from NUMERICAL RECIPES.
c
         integer n,NMAX
         REAL yp1,ypn,x(n),y(n),y2(n)
         parameter(NMAX=5000)
         integer i,k
         REAL p,qn,sig,un,u(NMAX)
         if(yp1.gt..99e30)then
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
         integer n
         real x,y,xa(n),y2a(n),ya(n)
         integer k,khi,klo
         real a,b,h
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
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
