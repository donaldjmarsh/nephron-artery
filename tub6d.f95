!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	MODULE tub6d_F90_m



!        local scalars

         real*8, private :: xt1
         real*8, private :: xt2
         real*8, private :: ydl6
         real*8, private :: res1201
         real*8, private :: gpress
         real*8, private :: g
         real*8, private :: h
         real*8, private :: ras
         real*8, private :: ras0
         real*8, private :: ras1
         real*8, private :: ras2
         real*8, private :: xrad
         real*8, private :: pt6
         real*8, private :: euknorm6
         real*8, private :: amm6

         real*8, private, parameter :: rhopde=1.0D0
         real*8, private, parameter :: tgfmax=0.455D0
         real*8, private, parameter :: dtgf=0.91D0
         real*8, private, parameter :: tpnacl=0.044D0
         real*8, private, parameter :: xgain=150.D0
         real*8, private, parameter :: kf=2.5D0
         real*8, private, parameter :: ca=57.D0
         real*8, private, parameter :: ha=0.5D0
         real*8, private, parameter ::  ks=3.3D0
         real*8, private, parameter ::  res=2.0D-02
         real*8, private, parameter :: rv=7.02D-02
         real*8, private, parameter :: Km=0.02D0
         real*8, private, parameter ::  K1=5.6D-7
         real*8, private, parameter ::  K2=1.3D0
         real*8, private, parameter ::  ph2o=2.D-5
         real*8, private, parameter :: bmm6=-0.1631D0
         real*8, private, parameter :: cmm6=-.00294D0
         real*8, private, parameter :: eta=7.2D-3
         real*8, private, parameter :: bf=-0.1631D0
         real*8, private, parameter :: cf=-0.00294D0

         integer, private :: irlen
         integer, private :: imd
         integer, private :: ilendl
         integer, private :: ilal
         integer, private :: ioffset
         integer, private :: istart
         integer, private :: jump
         integer, private :: jvmp
         integer, parameter, private :: npt = 400
         integer, parameter, private :: nefid = 6
         integer, parameter, private :: r0 = 195
         
!        local arrays

         real*8, private ::                                              &
          p(400,2)     ,q(0:400,2)     ,nacl(0:400,2)     ,vm(400)     , &
          cint(400)    ,comp(400)      ,reab(400)         ,pna(400)    , &
          radt(400)    ,xtub6(15)      ,v(400)            , &
          a(400,4)     ,b(400,4)       ,c(400,4)          ,d(400,2)    , &
          u(0:400)

         integer, private :: ilen(7)


!         public subroutines and functions

!          PUBLIC ::                                              &




!         private subroutines and functions

          PRIVATE ::                                              &

          bitri6    ,BOUND6    ,dfdpg6    ,feedb6   ,     &
          glom6     ,LUBKSB6   ,ludcmp6   ,f6       ,     &
          MNEWT6    ,MMNEWT6   ,nephron6  ,prflcl6  ,     &
          RADIUS6   ,REABGRAD6 ,tubgeom6  ,PDE6     ,     &
          chlorid6  ,stm16     ,stm26     ,stm36    ,     &
          stm46     ,stm56

!! Beginning of tub6 interface.                                                                 
!! ..                                                                                                
!! .. Generic Interface Blocks ..                                                                    
      INTERFACE tub6_F90

!        initall6 downloads all tubular parameters and initial conditions

         MODULE PROCEDURE initall6 

!        calc6 runs the calculations for the tubule

         MODULE PROCEDURE calc6
         
         END INTERFACE

CONTAINS

	subroutine prflcl6
	use affall
	use tuball

! subroutine for writing initial values for p,q, and nacl vectors

	implicit NONE

	integer :: i 
	double precision :: xqalh,del,xpfix
	character*12 initpq, initchl

!
!  Data files for output and the initial values 
!
      initpq='tgfpq.dat'
      initchl='tgfchl.dat'
!
!   Get initial values for the flow and pressure in the tubular system 
!   from file.
!

      open(40, FILE=initpq)

      if (irlen .le. r0) then
        read(40,'(11F7.2)')(p(i,1),i=1,irlen+1)
        read(40,'(11F7.2)')(q(i,1),i=0,irlen)
      else
        read(40,'(11F7.2)')(p(i,1),i=1,r0+1)
        read(40,'(11F7.2)')(q(i,1),i=0,r0)
      endif
      close(40)

      xqalh=q(r0,1)
      
      open(30,FILE=initchl)
      if (irlen .le. r0) then
        read(30,'(13F6.3)')(nacl(i,1),i=nint(1.D0/dz),irlen)
      else
        read(30,'(13F6.3)')(nacl(i,1),i=nint(1.D0/dz),r0)
      endif
      close(30)


      ! check for length of loop of Henle

	if (irlen .gt. r0) then

    ! set NaCl in tubular fluid to cint, and initialize p and q 

	  jump= 2*ilen(7)

	  do i=irlen+1,ilen(5)+1,-1
	    nacl(i-1,1)=nacl(i-jump -1,1)
	    p(i,1)=p(i-jump,1)
	  end do

	  do i=ilen(5),ilen(2)
	    nacl(i,1)=cint(i)
	  enddo

          jvmp=ilen(6)-jump
         do i = ilen(6),ilen(2)+1,-1
          nacl(i,1)=cint(jvmp)
          jvmp = jvmp+1
         enddo

          do i=ilen(2)+1,ilen(4)
            q(i,1) = xqalh
          enddo

          ! fixup p and q values 

	do i=ilen(1), irlen+1
	  del = p(i,1)-p(i-1,1)
	  if (del .gt. 0.0) then
	  istart = i
	  exit
	  endif
	enddo
	xpfix = (p(irlen+1,1)-p(istart-1,1))/float(irlen-istart+2)

	do i=istart, irlen+1
	  p(i,1)=p(i-1,1)+xpfix
	enddo

       ! set vm in talh

	do i=ilen(2)+1, ilen(6)
	  vm(i)=0.0D0
	enddo
	endif

      do i=1,irlen+1
        q(i-1,1)=q(i-1,1)*conv2
        q(i-1,2)=q(i-1,1)
        p(i,1)=p(i,1)/conv1
        p(i,2)=p(i,1)
      enddo

      do i=ilen(1),irlen
        nacl(i,2)=nacl(i,1)
      enddo

      ! talh

	ioffset = 0

	if (ilen(7) .gt. 0) then
	   ioffset=ilen(2)-ilen(1)
	  do i=ilen(1)+1, ilen(2)
	   pna(i)=2.96D-6
	   pna(i+ioffset)=2.9D-7
	  end do
	endif

       ! Distal segment

       do i=ilen(3)+nint(0.2D0/dz)+1,irlen
         vm(i)=6.52D-8
	enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine tubgeom6 
	use affall
	use tuball

       ! subroutine to calculate tubule length and interstitial NaCl


	implicit NONE

        integer :: i

	ilen(1)=nint(xp/dz)
	ilen(2)=nint(xdl/dz)
	ilen(3)=nint(xmd/dz)
	ilen(4)=nint(xd/dz)
	ilen(5)=nint(xal/dz)
        ilendl=nint(xdl/dz)

	ydl6=ydl6+xydl(6)

       ilen(7)=nint(ydl6/dz)

       irlen=ilen(4)+2*ilen(7)

       ilen(4)=ilen(4)+2*ilen(7) ! distance to distal segment end
       ilen(3)=ilen(3)+2*ilen(7) !distance to MD
       ilen(2)=ilen(2)+ilen(7) ! distance to DLH end
       ilen(5)=ilen(2) -ilen(7) !distance to short DLH end
       ilen(6)=ilen(2)+ilen(7) !distance to tALH end and TALH start

       !  calculate medullary interstitial NaCl values

       ilal=nint(xal/dz)+2*ilen(7)
       ioffset = ilal-ilen(1)-1

	do i=ilen(1)+1,ilen(2)
         cint(i)=2.93D-01+2.D-01*tanh(2.88D-02*(real((i-ilendl),8)))
	 cint(i+ioffset)=cint(i)
	 ioffset=ioffset-2
	 vm(i)=0.D0
	enddo

       return
       end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine initall6 
	use affall
	use tuball

        !  The subroutine to initialize the tubule model leaving 
        !  initialization of the myogenic model to the main program.

	implicit NONE

	integer :: r0,ntrial,ilp 
        integer :: i, itstep, nfstep
        double precision :: res123, rea, tolx,tolf


	character*12 initpq,initchl

	xtub6(1)=35.1D0;xtub6(2)=106.9D0;xtub6(3)=55.D0
	xtub6(4)=0.209D0;xtub6(5)=.59D0;xtub6(6)=81.5D0
	xtub6(7:15)=0.0D0

        ras0=0.045D0

	xt1=0.0D0;xt2=0.0D0

	do i=1,npt
	 pna(i)=3.4D-07
	 vm(i)=1.0D-07
	 cint(i)=0.15D0
	 enddo

       ras1=ras0
       ras2=ras0

        gpress=xtub6(3)/conv1

	call tubgeom6 
	call prflcl6
	imd = ilen(3)

	return
	end subroutine initall6
!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc6(tubTime,jcount)

        use mpi
        use affall
        use tuball

        !  subroutine to calculate tubule dynamics

        implicit NONE

        integer :: i,j,ntrial,jcount,r,itstep,ierr
!        ireq(X) X= 2*6
        integer :: status(MPI_STATUS_SIZE),ireq(12)
        real*8 :: sum(3),tubTime 
        real*8 :: gfr,sbf,spf,signal,xmacd,ytime,tpress
        real*8 :: intv 
        real(kind=8):: bufSend6(3)
        
        !   Start iteration loop. Itstep counts iteration steps. 
        !   Exit when euclidian norms of pressure and 
        !   flow, and NaCl concentration changes < 1.D-3, and 1.D-5 
        !   respectively, and relative change in afferent resistance 
        !   < 1.D-6, or number of iterations exceed maxstep.

      itstep=0;euknorm6=1.D0;ras=1.0

      !  Store previous value for afferent resistance.

      ras=ras2

      !  Solve the glomerular and tubular equations

      call nephron6

      gpress=xtub6(3)/conv1

      !  Calculate the afferent resistance.


      xt2=feedb6(nacl(imd,1),nacl(imd,2))


      ras2=rasAff(nefid)


      !  Increase iteration counter

      itstep=itstep+1

      !  Send message and stop if iteration counter exceeds maxstep.

      if (itstep.gt.maxstep) then
         stop 'itstep exceeds maxstep in MAIN'
      endif

      ! Update all values for next time step.
       
      ras1=ras2; xt1=xt2
      
      do i=1,irlen+1;q(i-1,1)=q(i-1,2);p(i,1)=p(i,2);enddo

      do j=nint(1.D0/dz),irlen; nacl(j,1)=nacl(j,2); enddo

         ! Transmit data to main, check time for output

	sbf=xtub6(2);signal=xt2

	bufSend6(1)=sbf;bufSend6(2)=signal;bufSend6(3)=gpress


!  ireq(X) X related to send sequence

        call MPI_ISEND(bufSend6,3,MPI_DOUBLE_PRECISION, &
          0,6,MPI_COMM_WORLD,ireq(11),ierr)

 	if (jflag .eq. 0) then
	 if (jcount.EQ.iprint) then

                
            if (tubTime.GE.startprt) then   

	     ytime=tubTime-startprt
            gfr=xtub6(1);spf=xtub6(2)
            tpress=p(ilen(1),1)*conv1;xmacd=nacl(imd,1)*1.0d03

            write(lout(nefid),'(f8.3,7(2x,g12.4))') ytime ,tpress,xmacd,2.0D0*spf,  &
                 pressNod(nodeNumber(6))/dyn,(nacl(imd,1)+nacl(imd,2))/2.0D0, &
                 q(ilen(1),1)/conv2,gfr
            

          endif 
         endif       
	endif
!  ireq(X) X related to send sequence

	call MPI_WAIT(ireq(11),status,ierr)


      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine nephron6
      use affall
      use tuball

      !  Calculate pressure and flow as a two-point   
      !  non-linear boundary value problem

       implicit NONE


      integer :: i,j,k,r, ntrial, n, itstep 
      real*8 xnacl(400) 
      real*8 pa 
      real*8 sum(3), tolx, tolf 

      a=0D0;b=0D0;c=0D0;d=0D0;u=0D0;v=0D0

      do i=1,401;u(i-1)=0.D0;enddo
      do i=1,400;v(i)=0.D0;enddo


!  Parameters:  
!   maxstep=maximum number of iterations for estimation
!     of non-linear coefficients.
!      maxstep=50

!   n=number of equations for glomerulus
!   ntrial=number of max iterations for glomerulus
!   tolx,tolf=required precision for glomerulus solution

      n=6;ntrial=30;tolx=5.D-7;tolf=5.D-7

      !   Initialize intermediate variables for PDE

      u(irlen)=q(irlen,1);v(irlen+1)=p(irlen+1,1)

!   Get initial tubular radius

      call radius6

! Start iteration loop. 

      itstep=0;euknorm6=1.D0;sum(3)=1.D0

      do while (((euknorm6.GT.1D-3).OR.(sum(3).GT.1.D-5)) &
               .AND.(itstep.LE.maxstep))
!
!! Calculate initial guess for the GFR, q(0,2), using the 
!! Bowmans space pressure, p(1,1) as initial guess for 
!! pressure at next time step, p(1,2). Guess is placed in 
!! q(0,2) as first boundary condition. Equations solved by 
!! multidimensional Newton-Raphson method in subroutine
!! mnewt. *
!
      pt6=p(1,2)*conv1

      call MMNEWT6(ntrial,xtub6,n,tolx,tolf)

      q(0,2)=xtub6(1)*conv2


!  Get tubular reabsorption rates

      call reabgrad6 

! PDE gets coefficient values for finite difference 
! scheme for pressure and flow.

      call pde6


!  Bitri solves bitridiagonal system of equations using 
!    algorithm given in von Rosenbergs book.

      call bitri6 


!  Use the new value of Bowmans capsule pressure for
!  new guess for GFR. System of glomerular equations
!  in subroutine glom is solved by subroutine mnewt.

      pt6=v(1)*conv1

      call MMNEWT6(ntrial,xtub6,n,tolx,tolf)

      u(0)=xtub6(1)*conv2


!  Calculate euclidian norm of change in solution
!  vector, then replace the old guess with the new guess.

      sum(1)=0.0D0
      sum(2)=0.0D0

      do i=1,irlen+1
         sum(1)=sum(1)+dabs((q(i-1,2)-u(i-1))/q(i-1,2))
         sum(2)=sum(2)+dabs((p(i,2)-v(i))/p(i,2))
      enddo

      euknorm6=sum(1)+sum(2)

      do i=1,irlen+1;q(i-1,2)=u(i-1);p(i,2)=v(i);enddo

!  Calculate new value for tubular radius,using new pressure values.

      call radius6


!  Calculate tubular chloride concentration and afferent 
!  resistance. Then replace the values at old time step with
!  new ones and start next time step. If convergence fails  
!  send message and stop.

!  First store old value then calculate new guess for the tubular 
!  chloride concentration

      do i=ilen(1),irlen;xnacl(i)=nacl(i,2);enddo
        
      call chlorid6

!  Calculate the euclidian norm for the NaCl concentration.

      sum(3)=0.D0

      do i=ilen(1),irlen
         sum(3)=sum(3)+dabs((nacl(i,2)-xnacl(i))/xnacl(i))
      enddo


!  Increase iteration counter

      itstep=itstep+1

!  Start next iteration step for the loop.

     enddo

!  Send message and stop if iteration counter exceeds maxstep.

      if (itstep.gt.maxstep) then

       do i=ilen(1),irlen+1
	 write (*,*) i, nacl(i,1), q(i-1,1)
	enddo
         stop 'itstep exceeds maxstep in nephron6'
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine bitri6 

!  Solve bi-tridiagonal set of equations resulting from finite
!  difference scheme. a, b, c, d, and e are coefficients for 
!  equations, u and v are solution vectors and r the number 
!  of spatial steps.



      implicit NONE

      integer,parameter :: np=400


      integer r, ntrial, n,i 
      real*8 lam(np,4), gam(np,2), bet(4), my, delt(2)
      real*8 x(15), tolx, tolf,bound6

      data lam /1600*0.0D0/


! Calculate coefficients needed in the algorithm

      bet(1)=b(1,1)
      bet(2)=b(1,2)
      bet(3)=b(1,3)
      bet(4)=b(1,4)


      my=bet(1)*bet(4)-bet(2)*bet(3)

      delt(1)=d(1,1)
      delt(2)=d(1,2)

      do i=1,irlen

        gam(i,1)=(bet(4)*delt(1)-bet(2)*delt(2))/my
        gam(i,2)=(bet(1)*delt(2)-bet(3)*delt(1))/my

	  if( i /= irlen ) then 

        lam(i,2)=(bet(4)*c(i,2)-bet(2)*c(i,4))/my
        lam(i,4)=(bet(1)*c(i,4)-bet(3)*c(i,2))/my
        
        bet(1)=b(i+1,1)
        bet(2)=b(i+1,2)-a(i+1,1)*lam(i,2)
        bet(3)=b(i+1,3)
        bet(4)=b(i+1,4)-a(i+1,3)*lam(i,2)

        my=bet(1)*bet(4)-bet(2)*bet(3)

        if (i.eq.irlen-1) then

          g=(bet(4)*b(irlen,2)-bet(2)*(-b(irlen,4)))/my
           h=(bet(4)*(d(irlen,1)-a(irlen,1)*gam(irlen-1,1))- bet(2)&
                *(d(irlen,2)-a(irlen,3)*gam(irlen-1,1)))/my

!  Solve equations for second boundary condition. Equations 
!  contained in subroutine bound1 are solved by subroutine 
!  newt using the Newton-Raphson method. 

          ntrial=40;tolx=1D-12;tolf=1D-12;n=2
    
          x(1)=u(irlen);x(2)=v(irlen+1)

          call mnewt6(ntrial,x,n,tolx,tolf)
          
          v(irlen+1)=x(2)

          d(i+1,1)=d(irlen,1)+b(irlen,2)*v(irlen+1)
          d(i+1,2)=d(irlen,2)-b(irlen,4)*v(irlen+1)
        
        endif
 
        delt(1)=d(i+1,1)-a(i+1,1)*gam(i,1)
        delt(2)=d(i+1,2)-a(i+1,3)*gam(i,1)

        endif
      enddo


      !  Substitution back for dependent variables u and v.

      u(irlen)=gam(irlen,1);v(irlen)=gam(irlen,2)

      do i=irlen-1,1,-1
        u(i)=gam(i,1)-lam(i,2)*v(i+1)
        v(i)=gam(i,2)-lam(i,4)*v(i+1)
      enddo

      return
      end    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine RADIUS6
      use affall

!  Calculates tubular radius in different segments, given 
!  the pressure and the spatial index i, and total number 
!  of spatial indexes, r.

      implicit NONE

      integer,parameter :: np=400
      real*8 pcomp, pxint, a1, a2, a3, a4, pavg ,pint,x
      integer i,r 

!  Values for slope and intercept of proximal tubule
!  compliance function, comp and xint.

      pcomp=9.9758D-9;pxint=9.9825D-4

!  Values for the hyperbolic tangent function for loop compliance

      a1=1.1571D-3;a2=1.05D-3;a3=2.154D-4;a4=1.4D4

!  Interstitial pressure = 5 mmHg in in dyn/cm^2.

      pint=5*1.33322368D3

!  Tubular radius and compliance.

      do i=1,irlen

!  Tubular pressure at ith step

      pavg=(p(i,1)+p(i,2)+p(i+1,1)+p(i+1,2))/4.D0

!  Calculate position in tubule for radius.

      if ((real(i,8)-0.5D0).le. ilen(1)) then

!  Proximal tubule
 
        radt(i)=pcomp*(pavg-pint)+pxint
        comp(i)=pcomp

      else if(i .gt. ilen(1)) then 

!  Loop of Henle

!  Interceptet er en anelse lavere end brugt i artiklen, hvor
!  int=8.4D-4, Formaalet er at faa en anelse hoejere Proksimalt
!  tryk, saaledes at middelvaaerdien for gain saenkes.


        radt(i)=1.D-8*pavg+7.9D-4
        comp(i)=1.D-8

      endif

      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine BOUND6(x,alpha,beta)


        !  Solves equations for the boundary conditions

      implicit NONE

      real*8 x(15), alpha(15,15), beta(15), slope, inter

      slope=2.33D-8
      inter=1.738D-3

      beta(1)=-(x(1)-h-g*x(2))
      beta(2)=-(-x(1)+x(2)*(slope*x(2)+inter)**4)

      alpha(1,1)=1.D0
      alpha(1,2)=-g
      alpha(2,1)=-1.D0
      alpha(2,2)=(slope*x(2)+inter)**4+4*slope*x(2)*(slope*x(2)+inter)**3
 
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine REABGRAD6

        use affall
        use tuball
        !
        !  Calculates tubular reabsorption rate in different nephron
        !  segments.

        implicit NONE

        integer :: np=400
        integer :: i,r

        real(8) :: xdel,xK1,xdelp

        ! K1 and K2 are parameters of the exponential function
        ! determining the rate of proximal reabsorption.
        ! Phso is the water permeability.


        xdel=real((ilen(2)-ilen(1)),8)

        do  i=1,irlen

           if (i .le. ilen(1)) then

              
!  Proximal tubule

     if (pressNod(nodeNumber(6))/dyn .gt. 85.0D0) then
        xdelp = (pressNod(nodeNumber(6))/dyn - 85.D0)
        xK1 = K1*(1.0D0 - varphi*xdelp)
     else
        xK1=K1
     endif

     reab(i)=xK1*dexp(-K2*(dfloat(i)-0.5D0)*dz)

     
!           reab(i)=K1*dexp(-K2*(dfloat(i)-0.5D0)*dz)

	   else if (i .le. ilen(2)) then

!  Descending loop of Henle

       reab(i)=ph2o*2.D0*(cint(i)+(0.15D0/(2.D0*xdel) &
               -(nacl(i,1)+nacl(i,2))/2.D0))

        else
 
!  Ascending loop of Henle
         
           reab(i)=0.0D0
        endif
      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MNEWT6(ntrial,x,n,tolx,tolf)

!  Given an initial guess X for a root in N dimension, take 
!  NTRIAL Newton-Raphson steps to improve the root. Stop if 
!  the root converges in either the summed variable increments 
!  TOLX or summed function values TOLX. USRFUN is a user 
!  supplied subroutine, which returns the matrix of partial 
!  derivatives, and the negative of the function values.
!
        implicit NONE

      integer,parameter :: np=15
      integer :: i,k,n,ntrial 
      integer :: indx(np)
      double precision ::  x(np),alpha(np,np),beta(np)
      double precision :: errf,errx,tolf,tolx 

      first: do k=1,ntrial

             call bound6(x,alpha,beta)

      
          errf=0.D0
          second: do  i=1,n
              errf=errf+abs(beta(i)) 
	   enddo second
          if (errf.le.tolf) return

          !--- Solve linear equations using LU decomposition.

          call ludcmp6(alpha,n,np,indx) 
          call lubksb6(alpha,n,np,indx,beta)
          errx=0.D0
          third: do i=1,n
              errx=errx+abs(beta(i)) 
              x(i)=x(i)+beta(i)
	   enddo third
          if (errx.le.tolx) return
	enddo first
        stop 'MNEWT6 exceeding max number of iterations'
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine MMNEWT6(ntrial,x,n,tolx,tolf)

      implicit NONE

      integer,parameter :: np=15
      integer :: i,k,ntrial,n
      integer :: indx(np)
      real(kind=8) ::  errf, errx, tolx,tolf,xp
      real(kind=8) ::  x(np),alpha(np,np),beta(np)!,indx(np)

      first: do k=1,ntrial

          call glom6(x,alpha,beta)

	  errf=0.D0

          second: do i=1,n
              errf=errf+abs(beta(i))
	  enddo second
          if (errf.le.tolf) return

          !--- Solve linear equations using LU decomposition.
          
          call ludcmp6(alpha,n,np,indx) 
          call lubksb6(alpha,n,np,indx,beta)

          errx=0.D0

          third: do i=1,n
              errx=errx+abs(beta(i))
              x(i)=x(i)+beta(i)
          enddo third

          if (errx.le.tolx) return

          amm6=x(3)-pt6

          xp=amm6+bmm6*x(6)+cmm6*x(6)*x(6) 

          do while (xp.LE.0.0D0)
          fourth: do i=1,n
                 beta(i)=0.5D0*beta(i)
                 x(i)=x(i)-beta(i)
              enddo fourth
              amm6=x(3)-pt6
              xp=amm6+bmm6*x(6)+cmm6*x(6)*x(6)
          enddo 

  	  enddo first

      write (6,*) xp
      stop 'MMNEWT6 exceeding max number of iterations'

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ludcmp6(a,n,np,indx) 

!--- Given an NxN matrix A, with physical dimension NP, this 
! routine replaces it by the LU decomposition of a rowwise 
! permutation of itself. A and N are input. A is output, 
! arranged as in equation (2.3.14) above ( in NUMERICAL RECIPES )
! INDX is an output vector which records the row permutation 
! ffected by the partial pivoting; D is output as +-1 depending 
! on whether the number of row interchanges was even or odd, 
! respectively. This routine is used in combination with LUBKSB 
! to solve linear equations or invert a matrix.


        implicit NONE

!--- Largest expected N, and a small number.

      integer, parameter :: nmax=100
      integer :: i,j,k,n,np,imax,d
      integer :: indx(np)
      real(kind=8) :: aamax, dum, sum 
      real(kind=8) :: tiny=5.0D-15
  
!--- Dimension vv states the implicit scaling of each row.

      real(kind=8) :: a(np,np),vv(nmax)

     first: do i=1,n
          aamax=0.D0
          second: do  j=1,n
              if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))

	  enddo second
          if (aamax.eq.0.D0) then
             
            write(*,*)(a(i,j),j=1,n)
            write(*,*) i,j, n  
	stop 'Singular matrix - 6 '
	call exit
          endif
          vv(i)=1./aamax
	enddo first
!--- This is the loop over columns of Crout's method

      third: do j=1,n
          if (j.gt.1) then
              fourth :do i=1,j-1
                  sum=a(i,j)
                  if (i.gt.1) then
                      fifth: do k=1,i-1
                          sum=sum-a(i,k)*a(k,j)
		      enddo fifth
                      a(i,j)=sum
                  endif
	     enddo fourth

          endif

          aamax=0.

	  sixth: do i=j,n
              sum=a(i,j)
              if (j.gt.1) then
                  seventh: do k=1,j-1
                      sum=sum-a(i,k)*a(k,j)
		  enddo seventh
                  a(i,j)=sum
              endif
              dum=vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax=i
                  aamax=dum
              endif
	  enddo sixth

          if (j.ne.imax) then
              eighth: do k=1,n
                  dum=a(imax,k)
                  a(imax,k)=a(j,k)
                  a(j,k)=dum
	      enddo eighth

              vv(imax)=vv(j)
          endif
          indx(j)=imax
          if (j.ne.n) then
              if(a(j,j).eq.0.) a(j,j)=tiny
              dum=1./a(j,j)
              ninth: do i=j+1,n
                  a(i,j)=a(i,j)*dum
	      enddo ninth
          endif

	  enddo third

      if (a(n,n).eq.0.D0) a(n,n)=tiny
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine LUBKSB6(a,n,np,indx,b)
      

!--- Solves the set of N linear equations A*X=B. Here A is input,
!    not as the matrix A but as its LU decomposition, deter-
!    mined by the routine LUDCMP. INDX is input as the permuta-
!    tion vector returned by LUDCMP. B is input as the right-hand
!    side vector B, and returns with the solution vector X. 
!    A, N, NP and INDX are not modified by this routine and can 
!    be left in for succesive calls with different right-hand 
!    sides B. This routine takes into account the possibility 
!    that B will begin with many zero elements, so it is 
!    efficient for use in matrix inversion.

        implicit NONE

      integer :: n,np,i,ii,j,ll
      integer :: indx(np)
      real(kind=8) :: sum
      real(kind=8) :: a(np,np),b(n)

      ii=0
      first: do i=1,n
          ll=indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if (ii.ne.0) then
             second:  do j=ii,i-1
                  sum=sum-a(i,j)*b(j)
	     enddo second
          else if (sum.ne.0.D0) then
              ii=i
          endif
          b(i)=sum
	  enddo first
      third: do i=n,1,-1
          sum=b(i)
          if(i.lt.n) then
              fourth: do j=i+1,n
                  sum=sum-a(i,j)*b(j)
          enddo fourth
          endif
          b(i)=sum/a(i,i)
	  enddo third
      return
 
    end subroutine LUBKSB6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine PDE6 
      use affall
      use tuball

!   This subroutine calculates the coefficient of the finite
!   difference scheme needed to solve the PDE's for the flow
!   and pressure.

      implicit NONE
      
      integer,parameter :: np=400
      integer :: i,r
!      real*8 :: eta=7.2D-3


!            eta: viscosity of tubular fluid, (g/cm*s)


!   Calculate values for the non-linear coefficients of the PDE.
!   Comp(i) and reab(i) are tubular compliance and reabsorption

	
      do i=1,irlen
         b(i,1)=1.D0/dt + 4.D0*eta/(rhopde*radt(i)*radt(i))
         b(i,2)=-pi*radt(i)*radt(i)/(rhopde*dz)
         b(i,3)=1.D0/(2.D0*pi*radt(i)*comp(i)*dz)
         b(i,4)=1.D0/dt

         a(i,1)=b(i,1)
         a(i,3)=-b(i,3)

         c(i,2)=-b(i,2)
         c(i,4)=b(i,4)

         d(i,1)=b(i,2)*(p(i+1,1)-p(i,1))+(2.D0*b(i,4)-b(i,1)) &
              *(q(i,1)+q(i-1,1))

         d(i,2)=-b(i,3)*(q(i,1)-q(i-1,1)) + b(i,4)*(p(i+1,1)+p(i,1)) &
              - reab(i)/(pi*radt(i)*comp(i))

      enddo

      d(1,1)=d(1,1)-b(1,1)*q(0,2)
      d(1,2)=d(1,2)+b(1,3)*q(0,2)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine chlorid6
      use affall
      use tuball

! Calculates the tubular concentration of chloride

      implicit NONE


      integer :: r, flag, i
      integer,parameter :: np=400

      real*8 b1(np), b2(np), b3(np) 
      real*8 b4(np), b5(np), a1,a2,a3,a4,a5

 
! Find the tubular chloride concentrations.

      do i=ilen(1)+1,irlen

! Calculate the coefficients of the finite difference scheme.

      a1=(q(i,1)+q(i-1,1)+q(i,2)+q(i-1,2))/(8.D0*dz)
      a2=pi*radt(i)*radt(i)/(2.D0*dt)
      a3=reab(i)/4.D0 
      a4=Pna(i)/4.D0

      b1(i)=a1-a2+a3-a4
      b2(i)=a1-a2-a3+a4
      b3(i)=a1+a2+a3-a4
      b4(i)=a1+a2-a3+a4
      b5(i)=pna(i)*cint(i)
      enddo

      ! Iterate

      do i=ilen(1)+1,irlen


        a5=b5(i)-vm(i)*(1.D0/(1.D0+Km/ &
        ((nacl(i-1,1)+nacl(i,1)+nacl(i-1,2)+nacl(i,2))/4.D0)))

        nacl(i,2)=(b1(i)*nacl(i-1,2)-b2(i)*nacl(i,1) &
         +b3(i)*nacl(i-1,1)+a5)/b4(i)
        

      enddo

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 FUNCTION feedb6(nacl1,nacl2)

!     TGF response of macula densa and afferent arteriole

      real*8 nacl1,nacl2, nacl 

      nacl=(nacl1+nacl2)/2.D0
      feedb6=tgfmax-dtgf/(1.D0+dexp(xgain*(nacl-tpnacl)))

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 FUNCTION f6(ce,gpf)

!   This subroutine gives the value of the expression 
!   integ(1./ce^2*(pg-pt-pi(c))0-kf/gpf*ca

      implicit NONE
      real*8 :: ce, gpf,af, qr
      real*8 :: integ1, integ2 

      af=gpress-pt6
!      bf=-0.1631D0
!      cf=-0.00294D0

      qr=4.D0*af*cf-bf*bf 
  

      integ1=(bf/(2.D0*af*af))*log((af+bf*ce+cf*ce*ce)/(ce*ce))- &
      1.D0/(af*ce)+(bf*bf/(2.D0*af*af)-cf/af)*stm16(ce,af,bf,cf)

      integ2=(bf/(2.D0*af*af))*log((af+bf*ca+cf*ca*ca)/(ca*ca))- &
      1.D0/(af*ca)+(bf*bf/(2.D0*af*af)-cf/af)*stm16(ca,af,bf,cf)

      f6=integ1-integ2-kf/(gpf*ca)


      return
      end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 FUNCTION dfdpg6(CA,CE)  


!  Calculates gradient of dCE/dP, givet ved den implicitte
!  funktion.

      real*8 ca,ce ,a
      real*8 integ1,integ2,integ 

!  Calculation of integral of 1/C^2*(PG-PT-PI(C))^2 in [CA,CE]

      a=gpress-pt6
      integ1=stm56(ce,a,bf,cf)
      integ2=stm56(ca,a,bf,cf)
      integ=integ1-integ2

!  Calculation af dF/dPg i punktet CE

      dfdpg6=-integ
!
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine glom6(xtub,alfa,beta)
      use affall

!    This subroutine gives values of the algebraic equations and the
!    related first derivatives of the glomerular model. In this version the efferent arteriole
!    is dependent on the hydrostatic pressure, as it undergoes passive dilatation (compliance set
!    to 0 in this version).

      implicit NONE
      integer :: i,j,k
      real*8 :: alfa(15,15), beta(15),xtub(15)
      real*8 :: sngfr, gpf,  re, he, ce ,xres,yres,zres
      real*8 :: bcomp, acomp, pg0, aa, bb
      real*8 :: eps=1.D-05

      zres = pressNod(nodeNumber(6))/dyn
      if (zres .gt. 85.0D0) then
         yres = (zres - 85.D0)
         xres=res*(1.0D0-vartheta*yres)
      else
         xres =res
      endif


      do i=1,15
         beta(i)=0.0D0
         do j=1,15
           alfa(i,j)=0.0D0
         enddo
      enddo
      
      sngfr=xtub(1)
      gpf=xtub(2)
      gpress=xtub(3)

      re=xtub(4)
      he=xtub(5)
      ce=xtub(6)

      pg=gpress

      if (ce .le. ca) ce = ce + eps
      aa=0.1631D0
      bb=0.00294D0
      pg0=55.00
      bcomp=0.0
      acomp=2.D0/3.D0*bcomp

      beta(1)=sngfr-(1.-ca/ce)*gpf
      beta(2)=gpf-(1.-ha)*(pressNod(nodeNumber(6))/dyn-gpress)/(ras*exp(ks*ha))
      beta(3)=gpress-re*(gpf/(1-ha)-sngfr)-gpf*rv/(1-ha)
      beta(4)=re-(xres+acomp*(gpress-pg0))*exp(ks*he)
      beta(5)=he-1./(1+(ca/ce)*((1-ha)/ha))
      beta(6)=f6(ce,gpf)

      alfa(1,1)=-1.D0
      alfa(1,2)=1.D0-ca/ce
      alfa(1,6)=gpf*(ca/(ce*ce))
      alfa(2,2)=-1.D0
      alfa(2,3)=-(1.D0-ha)/(ras*exp(ks*ha))
      alfa(3,1)=-re
      alfa(3,2)=(re+rv)/(1-ha)
      alfa(3,3)=-1.D0
      alfa(3,4)=gpf/(1.D0-ha)-sngfr
      alfa(4,3)=acomp*exp(ks*he)
      alfa(4,4)=-1.D0
      alfa(4,5)=(xres+acomp*(gpress-pg0))*ks*exp(ks*he)
      alfa(5,5)=-1.D0
      alfa(5,6)=(ca*((1.D0-ha)/ha)/(ce+ca*((1.D0-ha)/ha))**2)
      alfa(6,2)=-kf/(ca*gpf*gpf)
      alfa(6,3)=-dfdpg6(ca,ce)
      alfa(6,6)=-1.D0/(ce*ce*(gpress-pt6-aa*ce-bb*ce*ce))

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Real*8 FUNCTION stm16(x,a,b,c)

!   Denne funktion beregner stamfunktionen til funktionen 1/X^2,
!   hvor X = ( a + bx + cx^2 ). q = 4ac-b^2. 
  
      real*8 a,b,c,x,q
  
      q=4.D0*a*c-b*b
  
      stm16=(-1.D0/sqrt(-q))*log((2.D0*c*x+b+sqrt(-q))/ &
       (sqrt(-q)-2.D0*c*x-b))
 
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 FUNCTION stm26(x,a,b,c)

      real*8 x,a,b,c 

      stm26=(1.D0/(2.D0*a))*log(x*x/(a+b*x+c*x*x))- &
      (b/(2.D0*a))*stm16(x,a,b,c)
  
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 FUNCTION stm36(x,a,b,c)
      real*8 x,a,b,c,q
  
      q=4.D0*a*c-b*b
      stm36=(2.D0*c*x+b)/(q*(a+b*x+c*x*x)) + &
      (2.D0*c/q)*stm16(x,a,b,c)
  
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 FUNCTION stm46(x,a,b,c)

      real*8 x,a,b,c
  
      stm46=1.D0/(2.D0*a*(a+b*x+c*x*x))-(b/(2.D0*a))* &
      stm36(x,a,b,c)+stm26(x,a,b,c)/a
  
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 FUNCTION stm56(x,a,b,c)

      double precision :: x,a,b,c 
  
      stm56=-1.D0/(a*x*(a+b*x+c*x*x))-(2.D0*b/a)*  &
      stm46(x,a,b,c)-(3.D0*c/a)*stm36(x,a,b,c)
  
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end module

