!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE tuball

      integer, parameter :: maxstep=500
      integer :: inr,ix 
      double precision :: std,am,avg,xv,xydl(24),y(150),dt
      double precision :: rootbp(4100)
      


CONTAINS

	subroutine initsys

	use affall

	implicit NONE

        double precision :: x,yy,z,w,xx
	integer :: lmn,i

        character*9 :: na1='prtstart1',na2 ='timestop1',na3='prtstart2'
        character*9 :: na4='timestop2',na11='prtstart3',na12='timestop3'
        character*2 :: na5 = 'G0', na10 = 'BP'
        character*4 :: na6 = 'gnef'
        character*9 :: na9 = 'xgi'
        character*5 :: na7 ='jobno'

!  read system wide parameters


	open (unit=22,file='inputsys.dat')

         read (22,*) prtstart(1),timestop(1),prtstart(2),timestop(2),&
              prtstart(3),timestop(3), gnef,xgi,g0,inr,jobno 

	write (*,'(//)')
        write (*,'(x,3(A9,2x,A9,2x),3x,A4,7x,A2,9x,A3,5x,A5)') na1, &
             na2,na3,na4 ,na11,na12 
        write (*,'(D9.2,6(2x,D9.2),2x,I5)') (prtstart(i),timestop(i),i=1,3)
	write (*,'(//)')
        write (*,'(2x,A4,12x,A2,9x,A3,5x,A5)') na6,na9,na5,na7
        write (*,'((2x,D9.2,2(3x,D9.2),2x,I5))') gnef,xgi,g0,jobno
	write (*,'(//)')

        read (22,*) std,am,ix
	read (22,*) dt, numberOfNephrons, xintv
        read (22,*) x
        read (22,*) (artRad(i),i=1,9)
        read (22,*) (artLen(i),i=1,9)
        close (22)
        phin=x

        return
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine initmm
        use affall


!  reads initial values for  myogenic model

        implicit NONE

        double precision :: deltat, xint, yint, betax, zyx,xj
        double precision :: rea,dum,ydl,beta(15) 

        integer :: i,j,k


!  read parameter values

        open(unit=21,file='inputash.dat')
        read (21,*) v1, v2, v4, v5
        read (21,*) v6, ca3,ca4,gpd
        read (21,*) vel,vk,vca,cond
        read (21,*) gel, gk, gca, kd
        read (21,*) bt, alpha, betax,kca
        read (21,*) qm, caim, psim, xkpsi
        read (21,*) x1,x2,x3
        read (21,*) x4,x5,x6,x7
        read (21,*) x8,x9,u1,u2
        read (21,*) u3,deltat,dum,y1
        read (21,*) y2,y3,y4,xnuref
        read (21,*) cairef, ax,bx,cx
        read (21,*) dx, xpi, Area, we
        read (21,*) S, sig0, wm, tau
        read (21,*) res0
        gk=gk/cond
        gel=gel/cond
        gca=gca/cond
        gpd=gpd/cond
        alpha=alpha*cond
	cond=1.0D0

        close(21)

        psiref=(cairef**qm)/((caim**qm)+(cairef**qm))
        xomref=psiref/(psim+psiref)

        open(unit=21,file='input.dat')

        read (21,*) (y(i), i=1,4)
        read (21,*) y(5),y(6),y(10)
	read (21,*) ydl, y(11),y(12),v7
        close(21)

        do i=1,4;y(i+6)=y(i);end do
           
        do i=1,10 !numberOfNephrons
           j=(i-1)*12
            do k=1,12;y(j+k)=y(k);enddo                                                            
            enddo

            rootbp=0.00D0

            open (unit=23, file='rand115.dat')

            do i=1,4100
               read(23,*) zyx
               rootbp(i)=zyx
            enddo
            close(23)

!            do i=1,4100
!               do j=1,3
!                  xj = real(j)
!                  rootbp(i+j) = rootbp(i) + (rootbp(i+4)-rootbp(i))*xj/4.0D0
!               enddo
!            enddo        

            
          end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gauss
        use affall

       implicit NONE


       real*8 :: yrand(10),xrand(9),sumart,sumartrand,xart,zyx
       integer :: i,j,n
       integer, allocatable :: seed(:)

  call random_seed(size=n)
  allocate(seed(n))
  
  seed=ix

  call random_seed(put=seed)

           call random_number(yrand)
           call random_number(xrand)

           do i =1,numberOfNephrons
              zyx = yrand(i)/2.0D0

              if (zyx .ge. .375D0) then
                 xydl(i) = 3.0D-02
              else if ((zyx .ge. .25D0) .and. (zyx .lt. .375D0)) then
                 xydl(i) = 2.0D-02
              else if ((zyx .ge. .125D0) .and. (zyx .lt. .25D0)) then
                 xydl(i) = 1.0D-02
              else
                 xydl(i) = 0.0D-02
              endif
              
         enddo
         do i=1,numberOfNephrons; write (*,'(2x,I2,3x,F5.2)') i,xydl(i);enddo

            sumart =0.0D0
            sumartrand=0.0D0

            do i = 1,9;sumart=sumart+ artLen(i);enddo

            do i = 1,9
               xart= xrand(i)-0.5D0
               artLen(i)=(1.0D0+0.5D0*xart)*artLen(i)
            enddo

          end subroutine gauss

        end MODULE tuball
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
