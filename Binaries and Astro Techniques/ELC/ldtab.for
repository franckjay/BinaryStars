          program ldtab
c
c   November 25, 1999
c
c   This program will return an intensity from a table
c   of compiled values.  One gives a temperature, gravity, and a mu.
c   The code first locates nearby values of log(g) and T in the table
c   and interpolates to the input mu.  Then a small matrix of intensities
c   in the T-log(g) plane is generated and the final intenisity value
c   is arrived at using a two dimensional interpolation scheme.
c
          parameter (maxlines=700,maxmu=64)
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       atmint(maxlines,maxmu,8),Nmu(maxlines),outinty(8)
          dimension xout(100),yout(100)
c
          call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &       atmint,Tmax,Tmin,gmax,gmin)
c
          write(*,*)'enter T, logg, and the filter index'
          read(*,*)Tin,gin,idx
c
          do 10 i=100,1,-1
            rmuin=float(i)/100.0
            xout(i)=rmuin
c            write(*,*)rmuin
            call computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &       atmT,atmg,atmmu,Nmu,
     &       atmint,Tmax,Tmin,gmax,gmin,outinty)
c
            yout(i)=outinty(idx)
 10       continue
c
          open(unit=28,file='ldtab.out',status='unknown')
          do 21 i=1,100
            write(28,102)xout(i),yout(i)
 21       continue
c
          close(28)
c
          open(unit=28,file='ldtab.normout',status='unknown')
          do 20 i=1,100
            yout(i)=yout(i)/yout(100)
            write(28,101)xout(i),yout(i)
 20       continue
c
          close(33)
c
 100      format(5(1pe13.6,1x))
 101      format(f4.2,3x,f9.6)
 102      format(f4.2,3x,1pe16.9)
          end
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine computeinty(Tin,gin,rmuin,maxlines,maxmu,Nlines,
     &       atmT,atmg,atmmu,Nmu,
     &       atmint,Tmax,Tmin,gmax,gmin,outinty)
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       atmint(maxlines,maxmu,8),Nmu(maxlines),outinty(8),tempgh(200),
     &       tempgl(200),ghinty(200,8),yscratch(200),tscratch(2),glinty(200,8),
     &       tinty(2,8),y2scratch(2)


          call locate(atmT,Nlines,Tin,indexT)
c          write(*,*)indexT
          tscratch(1)=atmT(indexT+1)
c
c   Now search the Tvalues equal to atmT(indexT+1) and find intensities
c   for all of the g values.
c
          Nhigh=0
          do 10 i=0,Nlines-(indexT+1)
            if(atmT(indexT+1+i).eq.atmT(indexT+1))then
c              write(*,*)atmg(indexT+1+i),atmT(indexT+1+i) 
              call indexinty(indexT+1+i,maxlines,maxmu,atmmu,atmint,Nmu,rmuin,
     %          outinty)
c              write(*,*)(outinty(j),j=1,8)
              write(33,*)atmg(indexT+1+i),outinty(5),atmT(indexT+1+i)
              Nhigh=Nhigh+1
              tempgh(Nhigh)=atmg(indexT+1+i)
              do 9 j=1,8
                ghinty(Nhigh,j)=outinty(j)
 9            continue
            else
              go to 15
            endif
 10       continue
 15       call locate(tempgh,Nhigh,gin,indexgh)
c
c   The input value of gin is between atmg(indexT+Nhigh) and 
c   atmg(indexT+nhigh+1).
c
          m=2
          k=min(max(indexgh-(m-1)/2,1),Nhigh+1-m)
          do 20 i=1,8
            do 19 j=1,Nhigh
              yscratch(j)=ghinty(j,i)
 19         continue
            call polint(tempgh(k),yscratch(k),m,gin,qqq,dy)
            outinty(i)=qqq
            tinty(1,i)=qqq
 20       continue
c          write(*,*)(outinty(j),j=1,8)
          write(34,*)gin,outinty(5)
c
c   Now search the Tvalues equal to atmT(indexT+1) and find intensities
c   for all of the g values.
c
          tscratch(2)=atmT(indexT)
          Nlow=0
          do 100 i=0,indexT-1
            if(atmT(indexT-i).eq.atmT(indexT))then
c              write(*,*)atmg(indexT-i),atmT(indexT-i)
              call indexinty(indexT-i,maxlines,maxmu,atmmu,atmint,Nmu,rmuin,
     %          outinty)
              write(33,*)atmg(indexT-i),outinty(5),atmT(indexT-i)
              Nlow=Nlow+1
              tempgl(Nlow)=atmg(indexT-i)
              do 90 j=1,8
                glinty(Nlow,j)=outinty(j)
 90           continue
            else
              go to 150
            endif
 100      continue
 150      call locate(tempgl,Nlow,gin,indexgl)
c
c   The input value of gin is between atmg(indexT-Nlow) and atmg(indexT-Nlow-1).
c
          m=2
          k=min(max(indexgl-(m-1)/2,1),Nlow+1-m)
          do 200 i=1,8
            do 190 j=1,Nlow
              yscratch(j)=glinty(j,i)
 190      continue
            call polint(tempgl(k),yscratch(k),m,gin,qqq,dy)
            outinty(i)=qqq
            tinty(2,i)=qqq
 200      continue
c          write(*,*)(outinty(j),j=1,8)
          write(34,*)gin,outinty(5)
c
c   Finally, take the final pass and interpolate between T
c
          do 300 i=1,8
            do 290 j=1,2
              y2scratch(j)=tinty(j,i)
 290      continue
            call polint(tscratch,y2scratch,m,Tin,qqq,dy)
            outinty(i)=qqq
 300      continue
c
          write(35,*)gin,outinty(1),outinty(2),outinty(3),outinty(4),outinty(5)
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&
c
          subroutine indexinty(index,maxlines,maxmu,atmmu,atmint,Nmu,rmuin,
     %        outinty)
c
c   November 25, 1999
c
c   This subroutine will read the return values of the intensity at the
c   angle rmuin for model index in the table.
c
          dimension atmmu(maxlines,maxmu),
     %       atmint(maxlines,maxmu,8),Nmu(maxlines)
          dimension xmu(128),ymu(128),outinty(8)
c
c   Copy the mu values to a one dimensional array
c
          do 10 i=1,Nmu(index)
            xmu(i)=atmmu(index,i)
 10       continue
c             
          N=Nmu(index)
          call locate(xmu,N,rmuin,muindex)
          m=2
          k=min(max(muindex-(m-1)/2,1),N+1-m)
          do 20 i=1,8
            do 19 j=1,Nmu(index)
              ymu(j)=atmint(index,j,i)
 19         continue
            call polint(xmu(k),ymu(k),m,rmuin,qqq,dy)
            if(muindex.eq.0)qqq=atmint(index,1,i)
            outinty(i)=qqq
 20       continue
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&
c
      SUBROUTINE LOCATE(XX,N,X,J)
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
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &       atmint,Tmax,Tmin,gmax,gmin)
c
c   November 25, 1999
c
c   This routine will read the file with the model atmosphere data.  The
c   name is assumed to be 'Avni.atm', the the form is assumed to be the
c   following:
c
c   Teff   g
c   N_mu
c   mu1    intyU intyB intyV intyR intyI intyJ intyH intyK
c   mu2    intyU intyB intyV intyR intyI intyJ intyH intyK
c   ...
c
c
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       atmint(maxlines,maxmu,8),Nmu(maxlines)
c
          character*40 blank
          ios=0
          open(unit=19,file='ELC.atm',status='old',err=999,iostat=ios)
c
c   Attempt to read in the values
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
          Tmax=-1000.0
          Tmin=10000000.0
          gmax=-10000.0
          gmin=11111.0
 8        do 10 i=1,maxlines
            read(19,*,end=15)atmT(i),atmg(i)
            read(19,*,end=15)Nmu(i)
            if(atmT(i).gt.Tmax)Tmax=atmT(i)
            if(atmT(i).lt.Tmin)Tmin=atmT(i)
            if(atmg(i).gt.gmax)gmax=atmg(i)
            if(atmg(i).lt.gmin)gmin=atmg(i)
            do 9 j=1,Nmu(i)
              read(19,*,end=15)atmmu(i,j),(atmint(i,j,k),k=1,8)
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
 100      format(a1,'Error:  I can''t find the file ''Avni.atm''!')
c
          return
          end
c
c
c
c
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
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
          IF(DEN.EQ.0.)then
            write(*,*)ho,hp
            PAUSE
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

