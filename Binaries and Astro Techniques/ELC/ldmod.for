          program ldmod
c
c   May 18, 2000
c
c   This routine will read an output intensity file from ELC
c   and write a file containing mu,intensity pairs for the visible
c   pixels of the chosen star at the chosen phase.
c
          dimension x(6),y(6)
          character*6 extension
c
          open(unit=20,file='ldmod.out',status='unknown')
c
          write(*,*)'enter the star number'
          read(*,*)istar
          if(istar.ne.2)istar=1
          write(*,*)'enter the phase in degrees'
          read(*,*)fase
c
          call getextension(fase,extension)
c
          if(istar.eq.1)then
            open(unit=39,file='star1inty.'//extension,status='unknown')
          else
            open(unit=39,file='star2inty.'//extension,status='unknown')
          endif
c
          read(39,*,end=15)xcent,ycent,flux
          do 10 i=1,1000000
            read(39,*,end=15)temp,ncorner,rmu
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
            write(20,100)rmu,temp
 10       continue
 15       close(39)
c
 100      format(f6.4,3x,e16.6)
c
          end
c
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
          character*6 extension
c
          if((phase.ge.0.0).and.(phase.lt.10.0))write(extension,100)phase
          if((phase.ge.10.0).and.(phase.lt.100.0))write(extension,101)phase
          if((phase.ge.100.0).and.(phase.lt.1000.0))write(extension,102)phase

 100      format('00',f4.2)
 101      format('0',f5.2)
 102      format(f6.2)
c
          return
          end
c
