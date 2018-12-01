          program plotsquare
c
c
c
          double precision fase
          dimension phase(500000),xmodel(500000),ymodel(500000)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
          character*9 extension
          character*25 fileout
c
          igrid=0
          write(*,*)'enter 1 to plot the grid'
          read(*,*)igrid
          if(igrid.lt.1)then
            imag=0
            write(*,*)'enter 1 to plot intensities as magnitudes'
            read(*,*)imag
            write(*,*)'enter the low and high value'
            read(*,*)tlow,thigh
            xxxmin=1.e37
            xxxmax=-1.0e30
          endif
c
          imod=0
          write(*,*)'enter the filter of the plot light curve'
          write(*,*)'to plot, or 0 to skip'
          read(*,*)imod

          if(imod.eq.1)open(unit=20,file='modelU.mag',status='old')
          if(imod.eq.2)open(unit=20,file='modelB.mag',status='old')
          if(imod.eq.3)open(unit=20,file='modelV.mag',status='old')
          if(imod.eq.4)open(unit=20,file='modelR.mag',status='old')
          if(imod.eq.5)open(unit=20,file='modelI.mag',status='old')
          if(imod.eq.6)open(unit=20,file='modelJ.mag',status='old')
          if(imod.eq.7)open(unit=20,file='modelH.mag',status='old')
          if(imod.eq.8)open(unit=20,file='modelK.mag',status='old')

          ymin=10000000.0
          ymax=-100093.0
          N=0
          if(imod.gt.0)then
            do 7 i=1,500000
              read(20,*,end=8)xmodel(i),ymodel(i)
              if(ymodel(i).gt.ymax)ymax=ymodel(i)
              if(ymodel(i).lt.ymin)ymin=ymodel(i)
 7          continue
 8          close(20)
            N=i-1
          endif
c
          call getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     &       Ntheta,Nradius,alb1,alb2,Nref,
     %       rLx,Period,fm,separ,gamma,t3,g3,SA3,xecl,sw1,sw2,sw3,sw4,
     $       idraw,iecheck,idint,iatm,ism1,
     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,
     &       isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey)
c
          write(*,*)'enter the phase in degrees'
          read(*,*)fase
          call getextension(fase,extension)
          call pgbegin(0,'?',1,1)
          call pgscf(2)
          call pgsch(1.0)
          rightcorner=0.1+0.80*8.5/11.0
          call pgvport(0.1,rightcorner,0.12,0.92)
          xlow=-(separ+0.19*separ)
          ylow=-(separ+0.19*separ)
          xhigh=(separ+0.19*separ)
          yhigh=(separ+0.19*separ)
          call pgwindow(xlow,xhigh,ylow,yhigh)
          call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
          call pglabel('x (R\d\(2281)\u)','y (R\d\(2281)\u)',' ') 
          call plotstar(fase,Teff2,idint,imag,igrid,tlow,thigh,xxxmin,xxxmax)
          call pgsci(1)
          call pgsch(1.0)
          if(imod.gt.0)then
            call pgvport(rightcorner+0.02,0.95,0.6,0.92)
            ylow=ymax+(ymax-ymin)/10.0
            yhigh=ymin-(ymax-ymin)/10.0
            call pgwindow(0.0,1.0,ylow,yhigh)
            call pgbox('bcnst',0.0,0,'bcmst',0.0,0)
            diffmin=100000000.
            do 10 i=1,N
              diff=abs(fase/360.0-xmodel(i))
              if(diff.lt.diffmin)then
                diffmin=diff
                index=i
              endif
 10         continue
            call pgpoint(1,xmodel(index),ymodel(index),17)
            call pgline(index,xmodel,ymodel)
          endif
          call pgend
c
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine plotstar(phase,Teff2,idint,imag,igrid,
     $       tlow,thigh,xxxmin,xxxmax)
c
          double precision phase
          dimension x(6),y(6)
          character*9 extension
c
          if(igrid.lt.1)then
            write(*,*)'enter the color table:'
            write(*,*)'0  =   b-w_linear'
            write(*,*)'3  =   red_temperature'
            write(*,*)'5  =   std_gamma-II'
            read(*,*)itab
c
            if(itab.eq.3)then
              open(unit=20,file='red_temperature.tab',status='old')
              go to 6
            endif
c
            if(itab.eq.5)then
              open(unit=20,file='std_gamma-II.tab',status='old')
              go to 6
            endif

            open(unit=20,file='b-w_linear.tab',status='old')
          
c
c          open(unit=20,file='redtemp.0.5')
c
 6          do 7 i=16,255+15      !,1,-1
              read(20,*)icr,icg,icb
              cr=float(icr)/255.0
              cg=float(icg)/255.0
              cb=float(icb)/255.0
              call pgscr(i,cr,cg,cb)  
 7          continue
c
            close(20)
          endif
c
          call getextension(phase,extension)
c
          open(unit=39,file='star1inty.'//extension,status='unknown')
          read(39,*)xcent,ycent,flux
          do 10 i=1,1000000
            read(39,*,end=15)temp,ncorner,zz1,zz2,zz3,ialf,ibet
            if(imag.ge.1)temp=-2.5*log10(temp)+40.0
            if(ncorner.eq.4)then
              read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4)
              x(5)=x(1)
              y(5)=y(1)
            endif
            if(ncorner.eq.5)then
              read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5)
              x(6)=x(1)
              y(6)=y(1)
            endif
            if(igrid.ge.1)then
              call pgline(ncorner+1,x,y)
            else
              icolor=nint(((temp-tlow)/(thigh-tlow))*255)
c            write(*,*)icolor,temp
              if(icolor.gt.255)icolor=255
              if(icolor.lt.0)icolor=0
              call pgsci(icolor+16)
              call pgslw(3)
              call pgline(ncorner+1,x,y)
              call pgslw(1)
              call pgpoly(ncorner+1,x,y)
              if(temp.le.xxxmin)then
                xxxmin=temp
                iamin=ialf
                ibmin=ibet
              endif
              if(temp.ge.xxxmax)then
                xxxmax=temp
                iamax=ialf
                ibmax=ibet
              endif
            endif
 10       continue
 15       close(39)
c
          write(*,*)'star 1 min, max =',xxxmin,xxxmax
          write(*,*)iamin,ibmin,iamax,ibmax
c
          xxxmax=-12345.
          if(teff2.gt.0.0)then
            open(unit=39,file='star2inty.'//extension,status='unknown')
            read(39,*)xcent,ycent,flux
            do 20 i=1,1000000
              read(39,*,end=25)temp,ncorner
              if(imag.ge.1)temp=-2.5*log10(temp)+40.0
              if(ncorner.eq.4)then
                read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4)
                x(5)=x(1)
                y(5)=y(1)
              endif
              if(ncorner.eq.5)then
                read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5)
                x(6)=x(1)
                y(6)=y(1)
              endif
              if(igrid.ge.1)then
                call pgline(ncorner+1,x,y)
              else
                icolor=nint(((temp-tlow)/(thigh-tlow))*255 )
c                write(*,*)icolor,temp
                if(icolor.gt.255)icolor=255
                if(icolor.lt.0)icolor=0
                call pgsci(icolor+16)
                call pgslw(3)
                call pgline(ncorner+1,x,y)
                call pgslw(1)
                call pgpoly(ncorner+1,x,y)
                if(temp.le.xxxmin)xxxmin=temp
                if(temp.ge.xxxmax)xxxmax=temp
              endif
 20         continue
 25         close(39)
            write(*,*)'star 2 min, max = ',xxxmin,xxxmax
          endif
c 
          xxxmax=-12345.
          if(idint.gt.0)then
            open(unit=39,file='diskinty.'//extension,status='unknown')
            read(39,*)xcent,ycent,flux
            do 30 i=1,1000000
              read(39,*,end=35)temp,ncorner
              if(imag.ge.1)then
                if(temp.gt.0.0)then
                  temp=-2.5*log10(temp)+40.0
                else
                  temp=0.0
                endif
              endif
              if(ncorner.eq.4)then
                read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4)
                x(5)=x(1)
                y(5)=y(1)
              endif
              if(ncorner.eq.5)then
                read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5)
                x(6)=x(1)
                y(6)=y(1)
              endif
              if(igrid.ge.1)then
                call pgline(ncorner+1,x,y)
              else
                icolor=nint(((temp-tlow)/(thigh-tlow))*255 )
c                write(*,*)icolor,temp
                if(icolor.gt.255)icolor=255
                if(icolor.lt.0)icolor=0
                call pgsci(icolor+16)
                call pgslw(3)
                call pgline(ncorner+1,x,y)
                call pgslw(1)
                call pgpoly(ncorner+1,x,y)
                if(temp.le.xxxmin)xxxmin=temp
                if(temp.ge.xxxmax)xxxmax=temp
              endif
 30         continue
 35         close(39)
            write(*,*)'disk min, max = ',xxxmin,xxxmax
          endif
c 
          return
c
          end

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
          double precision phase
          character*9 extension
c
          if((phase.ge.0.0).and.(phase.lt.10.0))write(extension,100)phase
          if((phase.ge.10.0).and.(phase.lt.100.0))write(extension,101)phase
          if((phase.ge.100.0).and.(phase.lt.1000.0))write(extension,102)phase

 100      format('00',f7.5)
 101      format('0',f8.5)
 102      format(f9.5)
c
          return
          end
c
          subroutine erasestar(phase,tlow,thigh,xxxmin,xxxmax)
c
          double precision phase
          dimension x(6),y(6)
          character*9 extension
c
          call getextension(phase,extension)
c
          open(unit=39,file='star1inty.'//extension,status='unknown')
          do 10 i=1,10000
            read(39,*,end=15)temp,ncorner
            if(ncorner.eq.4)then
              read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4)
              x(5)=x(1)
              y(5)=y(1)
            endif
            if(ncorner.eq.5)then
              read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5)
              x(6)=x(1)
              y(6)=y(1)
            endif
c            icolor=int((temp-tlow)/(thigh-tlow))*255
c            write(*,*)icolor,temp
c            if(icolor.gt.255)icolor=255
c            if(icolor.lt.0)icolor=0
            call pgsci(0)
            call pgline(ncorner+1,x,y)
            call pgpoly(ncorner+1,x,y)
            if(temp.le.xxxmin)xxxmin=temp
            if(temp.ge.xxxmax)xxxmax=temp
 10       continue
 15       close(39)
c
          write(*,*)xxxmin,xxxmax
c
          open(unit=39,file='star2inty.'//extension,status='unknown')
          do 20 i=1,10000
            read(39,*,end=25)temp,ncorner
            if(ncorner.eq.4)then
              read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4)
              x(5)=x(1)
              y(5)=y(1)
            endif
            if(ncorner.eq.5)then
              read(39,*)x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4),x(5),y(5)
              x(6)=x(1)
              y(6)=y(1)
            endif
c            icolor=int((temp-tlow)/(thigh-tlow))*255 
c            write(*,*)icolor,temp
c            if(icolor.gt.255)icolor=255
c            if(icolor.lt.0)icolor=0
            call pgsci(0)
            call pgline(ncorner+1,x,y)
            call pgpoly(ncorner+1,x,y)
            if(temp.le.xxxmin)xxxmin=temp
            if(temp.ge.xxxmax)xxxmax=temp
 20       continue
 25       close(39)
c 
          write(*,*)xxxmin,xxxmax
          return
c
          end
c
c
c &&&&&&&&&&&&&&&&&&&&&&&
c
c
          subroutine getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     &       Ntheta,Nradius,alb1,alb2,Nref,
     %       rLx,Period,fm,separ,gamma,t3,g3,SA3,xecl,sw1,sw2,sw3,sw4,
     $       idraw,iecheck,idint,iatm,ism1,
     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,
     &       isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey)
c
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2)
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
          read(1,*)xecl
          read(1,*)sw1
          read(1,*)sw2
          read(1,*)sw3
          read(1,*)sw4
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
          close(1)
c
c   Come here if the input file ELC.inp does not exist.  The subroutine
c   writeinput will make the correct file and set default values.
c
 100      if(ios.gt.0)call writeinput(Nalph1,Nbet1,Nalph2,Nbet2,
     $       fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     $       Ntheta,Nradius,
     $       alb1,alb2,Nref,rLx,
     $       Period,fm,separ,gamma,t3,g3,SA3,xecl,sw1,sw2,sw3,sw4,
     $       idraw,iecheck,idint,iatm,ism1,
     &       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     $       iRVfilt,isw1,isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey)
c
          if(ios.eq.0)close(1)
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine writeinput(Nalph1,Nbet1,Nalph2,Nbet2,
     $       fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,
     $       Ntheta,Nradius,alb1,alb2,Nref,
     $       rLx,Period,fm,separ,gamma,t3,g3,SA3,xecl,sw1,sw2,sw3,sw4,
     $       idraw,iecheck,idint,iatm,ism1,
     &       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %       iRVfilt,isw1,isw2,isw3,isw4,
     &       ilaw,wave,dbolx,dboly,dwavex,dwavey)
c
c    will write the correctly formatted file ELC.inp and return
c    default parameters
c
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,2),dwavey(8,2),
     %       www(8)
          character*1 bell
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
          fill1=1.00
          fill2=0.005
          omega1=1.0
          omega2=1.0
          dphase=3.0
          Q=2.0
          finc=80.0
          Teff1=6500.0
          Teff2=6500.0
          Tgrav1=0.25
          Tgrav2=0.25
          betarim=2.0
          rinner=0.005
          router=0.75
          tdisk=30000.0
          xi=-0.75
          Ntheta=90
          Nradius=60
          alb1=1.0
          alb2=1.0
          Nref=1
          rLx=0.001
          Period=2.62
          fm=3.0
          separ=5.0
          gamma=50.0
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
          t3=-5000.0
          g3=-5.0
          SA3=-0.1
          xecl=-0.01
          sw1=0.0
          sw2=0.0
          sw3=0.0
          sw4=0.0
          iRVfilt=3
          isw1=0
          isw2=0
          isw3=0
          isw4=0
c
          do 5 i=1,8
            wave(i)=www(i)
            dbolx(i,1)=0.635
            dbolx(i,2)=0.635
            dboly(i,1)=0.242
            dboly(i,2)=0.242
 5        continue
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
          write(1,5003)xecl
          write(1,5004)sw1
          write(1,5005)sw2
          write(1,5006)sw3
          write(1,5007)sw4
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
 1006     format(f4.2,16x,'omega1')
 1007     format(f4.2,16x,'omega2')
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
 1021     format(f10.5,10x,'Lx/Lopt')
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
 2000     format(f7.1,3x,8(f6.3,2x))
 3028     format(i1,19x,'ilaw  (1=linear law, 2=logarithmic law,',
     %           ' 3=square root law)')
 4025     format(f7.2,13x,'gamma velocity (km/sec)')
 4000     format(i1,19x,'iatm  (0 for BB, 1 for model atmospheres)')
 4001     format(i1,19x,'ism1  (0 for all phases, 1 for 0-180)')
 4002     format(8(i1,1x),4x,'icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK')
 5000     format(f10.2,10x,'t3')
 5001     format(f5.2,15x,'g3')
 5002     format(f12.6,8x,'SA3')
 5003     format(f12.6,8x,'xecl')
 5004     format(f12.6,8x,'sw1 (currently inactive)')
 5005     format(f12.6,8x,'sw2 (currently inactive)')
 5006     format(f12.6,8x,'sw3 (currently inactive)')
 5007     format(f12.6,8x,'sw4 (currently inactive)')
 5008     format(i1,19x,'iRVfilt')
 5009     format(i1,19x,'isw1 (currently inactive)')
 5010     format(i1,19x,'isw2 (currently inactive)')
 5011     format(i1,19x,'isw3 (currently inactive)')
 5012     format(i1,19x,'isw4 (currently inactive)')
c
          return
          end


