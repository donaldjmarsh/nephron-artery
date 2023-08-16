!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               program opsplit

   ! This version uses tubular length [1.95,1.97,1.99,2.01,2.03] cm
                 
               use mpi
	       use dvode_f90_m
	       use affall
	       use tuball
               use tub1_F90_m
               use tub2_F90_m
               use tub3_F90_m
               use tub4_F90_m
               use tub5_F90_m
               use tub6_F90_m
               use tub7_F90_m
               use tub8_F90_m
               use tub9_F90_m
               use tub10_F90_m

!  This program calculates pressure and flow in a tubule as a 
!  two-point   non-linear boundary value problem, with two 
!  afferent arteriolar segments solved with dvode

       implicit NONE

       external diffun

       integer :: i,j,k,klm,lmn,mno,jcount,kcount,lcount,mcount
       integer :: r,ntrial, itstep, neq 
       integer :: myid, ierr, numprocs
       integer(kind=4), dimension(MPI_STATUS_SIZE,1) :: stat_recv
!      ireq dimensioned to 2*numberOfNephrons-1, i.e. the number of IRECV requests
       integer(kind=4) :: ireq(19)
       integer :: iopar, iopt, iout, istate, itask, &
         leniw, lenrw, liw, lrw, mband, meth,miter, &
          ml, mu, nerr, nfe, nfea, nje, nout, nqu, nst, &
           iwork(44)

  
       double precision :: rad,bifuradadd,bifuraddiff,xnodeelec
       real(kind=8) :: artgRad1,artgRad2,artgRad3,artgRad4,artgRad5, &
            artgRad6,artgRad7,artgRad8,artgRad9

!  dimensions of rtol and atol must be changed with numberOfNephrons

       double precision :: rwork(814), rtol(120),atol(120),sysTime,time 
       double precision :: hu,xtime

       double precision :: t,  dtout, ras0=0.045D0,tout,tout2
       real(kind=8) :: avrad(50), flowsum, flow1,flow23,flow4,flow5,flow67
       real(kind=8) :: flow8,flow910,xdist
       real(kind=8) :: ydist,pydist,nodediff(9),nodeadd(9),nodemeandiff(9)
       real(kind=8),parameter :: poisLenPi=(8.0D0*4.6D-02)/3.14159265D0
       real(kind=8) ::flowcgs=1.0D-06/60.0D0
       real(kind=8) :: ha=0.50D0,rha=2.0D0
       
       real(kind=8) :: bufRecv1(3),bufRecv2(3)
       real(kind=8) :: bufRecv3(3),bufRecv4(3)
       real(kind=8) :: bufRecv5(3),bufRecv6(3)
       real(kind=8) :: bufRecv7(3),bufRecv8(3)
       real(kind=8) :: bufRecv9(3),bufRecv10(3)
       
      character*3 mr,nr,irun
      character*2 filno,nyid

      TYPE(VODE_OPTS) :: OPTIONS

      artRad = 0.D0;artLen =0.D0      
 
      y=0.0D0

      call initmm 

      myid = 0
      call MPI_INIT ( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr)

      if (myid .eq. 0) then

      call initsys 

      neq=12*numberOfNephrons

         do i=1,numberOfNephrons
            k=12*(i-1)
            do j=1,12
              atol(k+j)=1.D-8
              rtol(k+j)=1.D-8
            enddo
         enddo

        OPTIONS = SET_OPTS(dense_J=.true.,RELERR_VECTOR=RTOL, &
        ABSERR_VECTOR=ATOL)

	 do i=1,25;lout(i) =i+20;potDiff(i)=-4.545D01;enddo
	 do i=1,24;pressGlom(i)=55.D0;enddo

!      assign matrix connections for node PD calculations

     call inverse


!       vary nephron length

        call gauss

       nerr=0
       rho = 1.0D0
       avrad = 0.0D0 

       dtout=dt

       iprint =(xintv/dt)


      do i=1,numberOfNephrons
      pressArt(i)=dfloat(inr)
      do j=1,2;affRad(i,j)=1.5D-03;enddo
         enddo

!      Assign node numbers for nephron PD calculations

        nodeNumber=0

        do i=1,numberOfNephrons
           nodeselect: select case (i)
           case (1)
              nodeNumber(i) = 3
           case (2:3)
              nodeNumber(i) = 4
           case (4)
              nodeNumber(i) = 5
           case (5)
              nodeNumber(i) = 6
           case(6:7)
              nodeNumber(i) = 8
           case(8)
              nodeNumber(i) = 9
           case(9:10)
              nodeNumber(i) = 10
           end select nodeselect
        enddo

        do j=1,numberOfNephrons
           PDnodeselect: select case (j)
           case (1)
              PDnodeNumber(j) = 2
           case (2:3)
              PDnodeNumber(j) = 3
           case (4)
              PDnodeNumber(j) = 4
           case (5)
              PDnodeNumber(j) = 5
           case(6:7)
              PDnodeNumber(j) = 7
           case(8)
              PDnodeNumber(j) = 8
           case(9:10)
              PDnodeNumber(j) = 9
           end select PDnodeselect
        enddo
        

       rootpress= dyn*dble(inr)
       write (*,*) inr
       
       do i =1,10;pressNod(i)=rootpress;enddo
          
    endif

       call MPI_BCAST(startprt,1,MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)

       call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)

       call MPI_BCAST(dtout,1,MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)

       call MPI_BCAST(stoptime,1,MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)


       call MPI_BCAST(inr,1,MPI_CHARACTER,0, &
            MPI_COMM_WORLD, ierr)

       call MPI_BCAST(jobno,1,MPI_INTEGER,0, &
            MPI_COMM_WORLD, ierr)


       call MPI_BCAST(iprint,1,MPI_INTEGER,0, &
            MPI_COMM_WORLD, ierr)

       call MPI_BCAST(lout,25,MPI_INTEGER,0, &
            MPI_COMM_WORLD, ierr)

       call MPI_BCAST(nodeNumber,15,MPI_INTEGER,0, &
            MPI_COMM_WORLD, ierr)

       call MPI_BCAST(xydl,24,MPI_DOUBLE_PRECISION,0, &
            MPI_COMM_WORLD, ierr)


       call MPI_BCAST(gnef, 1, MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr )

       call MPI_BCAST(pressNod,10,MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr )


              call MPI_BARRIER( MPI_COMM_WORLD, ierr)


!       open output files and initialize nephron

              if (myid == 1) call initall1
              if (myid == 2) call initall2
              if (myid == 3) call initall3
              if (myid == 4) call initall4
              if (myid == 5) call initall5
              if (myid == 6) call initall6
              if (myid == 7) call initall7
              if (myid == 8) call initall8
              if (myid == 9) call initall9
              if (myid == 10) call initall10
              
              call MPI_BARRIER(MPI_COMM_WORLD, ierr )
              
!  Initialize time and counters.
       
       nodemeandiff(1)=4.55D0
       nodemeandiff(2)=.063D0
       nodemeandiff(3)=.86D0
       nodemeandiff(4)=.61D0
       nodemeandiff(5)=.610D0
       nodemeandiff(6)=.630D0
       nodemeandiff(7)=.430D0
       nodemeandiff(8)=1.27D0
       nodemeandiff(9)=.245D0
       
              xdist=1.0D0;lcount=1;ydist=1.0D0;ncount=0;pydist =.02D0
              lflag=0
              iterateradii: do while (lcount .lt. 3)

       t=0.0d0
       sysTime=0.0;jflag=0;jcount=0;kcount=0;mcount=0
       istate=1; itask=1
       nodeadd = 0.0D0

       startprt = prtstart(lcount)
       stoptime = timestop(lcount)

       call MPI_BCAST(startprt,1,MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)

       call MPI_BCAST(stoptime,1,MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr)

              call MPI_BARRIER(MPI_COMM_WORLD, ierr )


       ! assign names and numbers to output files

       if (myid .eq. 0) then
          call outname(irun,nr)
          open(unit=20,FILE='ntw0'//trim(nr)//irun)
       endif


        if (myid .eq. 1) then
        call outname(irun,nr)
         open(lout(1),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
         endif

         if (myid .eq. 2) then
            call outname(irun,nr)
            open(lout(2),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
         endif

         if (myid .eq. 3) then
            call outname(irun,nr)
            open(lout(3),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
         endif

         if (myid .eq. 4) then
            call outname(irun,nr)
            open(lout(4),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
         endif

         if (myid .eq. 5) then
            call outname(irun,nr)
            open(lout(5),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
         endif

          if (myid .eq. 6) then
            call outname(irun,nr)
            open(lout(6),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
            endif
         
          if (myid .eq. 7) then
            call outname(irun,nr)
            open(lout(7),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
           endif

          if (myid .eq. 8) then
            call outname(irun,nr)
            open(lout(8),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
           endif


          if (myid .eq. 9) then
            call outname(irun,nr)
            open(lout(9),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
           endif

          if (myid .eq. 10) then
            call outname(irun,nr)
            open(lout(10),FILE='ntw'//trim(nyid(myid))//trim(nr)//irun)
           endif

         ! initialize run

              
       tout=sysTime+dt
       tout2=tout/2.D0

       ! solve afferent arteriolar ODEs

       if (myid == 0) then       

        call dvode_F90(diffun,neq,y,t,tout2,itask, &
          istate,options)
        
        tout=tout2+dt

        hu = rwork(11)
        nqu = iwork(14)


        if (istate .lt. 0) go to 175

     endif
     
     ! calculate neprhon blood flows for network vascular node values

     flow1=0.0D0;flow23=0.0D0;flow4=0.0D0;flow5=0.0D0;flow67=0.0D0
     flow8=0.0D0;flow910=0.0D0;flowsum=0.D0


     j=12

         firstRadii: do i=1,numberOfNephrons
	 do k=1,2
	 if (k .eq. 1) then	 	 
	  affRad(i,k)=rad(y(((i-1)*j)+5))

	 else 
	  affRad(i,k)=rad(y(((i-1)*j)+11))

	  endif
         enddo
         enddo firstRadii

        ! rasAff calculated from segres in dvode

         CALL MPI_BCAST(rasAff,24,MPI_DOUBLE_PRECISION,0, &
              MPI_COMM_WORLD, ierr)


      therun: do while(sysTime.LE.stoptime)
	jcount=jcount+1

         CALL MPI_BCAST(jcount,1,MPI_INTEGER,0, &
           MPI_COMM_WORLD, ierr)

        call MPI_BARRIER(MPI_COMM_WORLD, ierr )

        if (myid == 1) then
           call calc1(sysTime,jcount)
        endif

        if (myid ==2) then
           call calc2(sysTime,jcount)
        endif

        if (myid == 3) then
           call calc3(sysTime,jcount)
        endif

        if (myid == 4) then
           call calc4(sysTime,jcount)
        endif


        if (myid == 5) then
           call calc5(sysTime,jcount)
        endif

        if (myid == 6) then
           call calc6(sysTime,jcount)
        endif

           if (myid == 7) then
           call calc7(sysTime,jcount)
        endif

        if (myid == 8) then
           call calc8(sysTime,jcount)
        endif

        if (myid == 9) then
           call calc9(sysTime,jcount)
        endif

        if (myid == 10) then
           call calc10(sysTime,jcount)
        endif

        call MPI_BARRIER ( MPI_COMM_WORLD,  ierr)

        if (lcount .gt. 1) then

        endif
        if (myid == 0) then
           call MPI_IRECV(bufRecv1,3,MPI_DOUBLE_PRECISION, &
                1,1,MPI_COMM_WORLD,ireq(1),ierr)

           call MPI_IRECV(bufRecv2,3,MPI_DOUBLE_PRECISION, &
                2,2,MPI_COMM_WORLD,ireq(2),ierr)

           call MPI_IRECV(bufRecv3,3,MPI_DOUBLE_PRECISION, &
                3,3,MPI_COMM_WORLD,ireq(3),ierr)

           call MPI_IRECV(bufRecv4,3,MPI_DOUBLE_PRECISION, &
                4,4,MPI_COMM_WORLD,ireq(4),ierr)

           call MPI_IRECV(bufRecv5,3,MPI_DOUBLE_PRECISION, &
                5,5,MPI_COMM_WORLD,ireq(5),ierr)

           call MPI_IRECV(bufRecv6,3,MPI_DOUBLE_PRECISION, &
                6,6,MPI_COMM_WORLD,ireq(6),ierr)

           call MPI_IRECV(bufRecv7,3,MPI_DOUBLE_PRECISION, &
                7,7,MPI_COMM_WORLD,ireq(7),ierr)

           call MPI_IRECV(bufRecv8,3,MPI_DOUBLE_PRECISION, &
                8,8,MPI_COMM_WORLD,ireq(8),ierr)

           call MPI_IRECV(bufRecv9,3,MPI_DOUBLE_PRECISION, &
                9,9,MPI_COMM_WORLD,ireq(9),ierr)

           call MPI_IRECV(bufRecv10,3,MPI_DOUBLE_PRECISION, &
                10,10,MPI_COMM_WORLD,ireq(10),ierr)

           call MPI_WAITALL(10,ireq,stat_recv,ierr)

           
           rbf(1)=rha*bufRecv1(1)
           rbf(2)=rha*bufRecv2(1)
           rbf(3)=rha*bufRecv3(1)
           rbf(4)=rha*bufRecv4(1)
           rbf(5)=rha*bufRecv5(1)
           rbf(6)=rha*bufRecv6(1)
           rbf(7)=rha*bufRecv7(1)
           rbf(8)=rha*bufRecv8(1)
           rbf(9)=rha*bufRecv9(1)
           rbf(10)=rha*bufRecv10(1)

           yt(1)=bufRecv1(2)
           yt(2)=bufRecv2(2)
           yt(3)=bufRecv3(2)
           yt(4)=bufRecv4(2)
           yt(5)=bufRecv5(2)
           yt(6)=bufRecv6(2)
           yt(7)=bufRecv7(2)
           yt(8)=bufRecv8(2)
           yt(9)=bufRecv9(2)
           yt(10)=bufRecv10(2)

           pressGlom(1)=bufRecv1(3)
           pressGlom(2)=bufRecv2(3)
           pressGlom(3)=bufRecv3(3)
           pressGlom(4)=bufRecv4(3)
           pressGlom(5)=bufRecv5(3)
           pressGlom(6)=bufRecv6(3)
           pressGlom(7)=bufRecv7(3)
           pressGlom(8)=bufRecv8(3)
           pressGlom(9)=bufRecv9(3)
           pressGlom(10)=bufRecv10(3)

           flow1=0.0D0;flow23=0.0D0;flow4=0.0D0;flow5=0.0D0;flow67=0.0D0
           flow8=0.0D0;flow910=0.0D0;flowsum=0.0D0


           do mno=1,numberOfNephrons
              flowsum=flowsum+rbf(mno)
              if (mno == 1) then
                 flow1=flow1+rbf(mno)
              else if (mno == 2) then
                 flow23=flow23+rbf(mno)
              else if (mno == 3) then
                 flow23=flow23+rbf(mno)
              else if (mno == 4) then
                 flow4=flow4+rbf(mno)
              else if (mno == 5) then
                 flow5=flow5+rbf(mno)
              else if (mno == 6) then
                 flow67=flow67+rbf(mno)
              else if (mno == 7) then
                 flow67=flow67+rbf(mno)
              else if (mno == 8) then
                 flow8=flow8+rbf(mno)
              else if (mno == 9) then
                 flow910=flow910+rbf(mno)
              else
                 flow910=flow910+rbf(mno)
              endif
           enddo


           flowsum=flowcgs*flowsum
           flow1 =flowcgs*flow1
           flow23=flowcgs*flow23
           flow4 =flowcgs*flow4
           flow5=flowcgs*flow5
           flow67 =flowcgs*flow67
           flow8 =flowcgs*flow8
           flow910=flowcgs*flow910
           
           pressNod(1)=rootpress
           pressNod(2)=pressNod(1)-(flowsum*poisLenPi*artLen(1)/(artRad(1)**4.0D0))
           pressNod(3)=pressNod(2)-((flow1+flow23)*poisLenPi*artLen(2)/(artRad(2)**4.0D0))
           pressNod(4)=pressNod(3)-(flow23*poisLenPi*artLen(3)/&
                (artRad(3)**4.0D0))
           pressNod(5)=pressNod(2)-((flowsum-flow1-flow23)*poisLenPi*artLen(4)/ &
                (artRad(4)**4.0D0))
           pressNod(6)=pressNod(5)-((flowsum-flow1-flow23-flow4)*poisLenPi*artLen(5)/ &
                (artRad(5)**4.0D0))
           pressNod(7)=pressNod(6)-((flowsum -flow1-flow23-flow4-flow5))*poisLenPi*artLen(6)/ &
             (artRad(6)**4.00)
           pressNod(8)=pressNod(7)-(flow67*poisLenPi*artLen(7)/(artRad(7)**4.0D0))
           pressNod(9)=pressNod(7)-((flow8+flow910)*poisLenPi*artLen(8)/(artRad(8)**4.0D0))
           pressNod(10)=pressNod(9)-(flow910*poisLenPi*artLen(9)/(artRad(9)**4.0D0))


            ! calculate pressure differences between contiguous vascular nodes

            nodediff(1)=(rootpress-pressNod(2))/dyn
            nodediff(2)=(pressNod(2)-pressNod(3))/dyn
            nodediff(3)=(pressNod(3)-pressNod(4))/dyn
            nodediff(4)=(pressNod(2)-pressNod(5))/dyn
            nodediff(5)=(pressNod(5)-pressNod(6))/dyn          
            nodediff(6)=(pressNod(6)-pressNod(7))/dyn
            nodediff(7)=(pressNod(7)-pressNod(8))/dyn
            nodediff(8)=(pressNod(7)-pressNod(9))/dyn
            nodediff(9)=(pressNod(9)-pressNod(10))/dyn 

            
           

           call dvode_F90(diffun,neq,y,t,tout,itask, &
             istate,options)

        hu = rwork(11)
        nqu = iwork(14)

        if (istate .lt. 0) then
           write (*,*) sysTime,istate
           go to 175
        endif

	j=12

         laterRadii: do i=1,numberOfNephrons
	 do k=1,2
	 if (k .eq. 1) then	 	 
	  affRad(i,k)=rad(y(((i-1)*j)+5))
	 else 
	  affRad(i,k)=rad(y(((i-1)*j)+11))
	  endif
         enddo
         enddo laterRadii

         endif

         call MPI_BCAST(rasAff,24,MPI_DOUBLE_PRECISION,0, &
           MPI_COMM_WORLD, ierr)


         call MPI_BCAST(potNode,15,MPI_DOUBLE_PRECISION,0, &
           MPI_COMM_WORLD, ierr)

         call MPI_BCAST(pressArt,15,MPI_DOUBLE_PRECISION,0, &
              MPI_COMM_WORLD, ierr)

         call MPI_BCAST(pressNod,10,MPI_DOUBLE_PRECISION, 0, &
            MPI_COMM_WORLD, ierr )

         

         jflag=0

        sysTime=sysTime+dt
	tout=tout+dtout

         do i = 1,9
          nodeadd(i)=nodeadd(i)+nodediff(i)
         enddo
        mcount=mcount+1
        if (mod(jcount,iprint) == 0) then
          jcount=0

         endif
 
       if (jcount .eq. 0 .and. sysTime .ge. startprt) then
          write (20,'(4x,F8.3,20(2x,g15.6))') sysTime-startprt, & 
               (potNode(i),i=1,9),(pressNod(j)/dyn,j=1,10)
       endif
 
         kcount=kcount+1


          do lmn=1,numberOfNephrons
          mno=12*(lmn-1)
          avrad(lmn)=avrad(lmn)+rad(y(mno+5))
       enddo

      enddo therun




        call MPI_BARRIER ( MPI_COMM_WORLD,  ierr)


!     avrad(i) is the radius of the ith afferent arteriole

      if (myid == 0) then
      avrad(1)= avrad(1)/dble(kcount)
      avrad(2)= avrad(2)/dble(kcount)
      avrad(3)= avrad(3)/dble(kcount)
      avrad(4)= avrad(4)/dble(kcount)
      avrad(5)= avrad(5)/dble(kcount)
      avrad(6)= avrad(6)/dble(kcount)
      avrad(7)= avrad(7)/dble(kcount)
      avrad(8)= avrad(8)/dble(kcount)
      avrad(9)= avrad(9)/dble(kcount)
      avrad(10)= avrad(10)/dble(kcount)

      if (lcount == 1) then

      artgRad1=bifuradadd(avrad(1),artRad(2))
      artgRad2=bifuradadd(artRad(3),artRad(4))
      artgRad3=bifuradadd(avrad(2),avrad(3))
      artgRad4=bifuradadd(avrad(4),artRad(5))
      artgRad5=bifuradadd(artRad(6),artRad(9))
      artgRad6=bifuradadd(artRad(7),avrad(7))
      artgRad7=bifuradadd(avrad(5),avrad(6))
      artgRad8=bifuradadd(avrad(9),avrad(10))
      artgRad9=bifuradadd(avrad(8),artRad(8))

        xdist = abs((artRad(1)-artgRad1))/artRad(1) + &
                abs((artRad(2)-artgRad2))/artRad(2) + &
                abs((artRad(3)-artgRad3))/artRad(3) + &
                abs((artRad(4)-artgRad4))/artRad(4) + &
                abs((artRad(5)-artgRad5))/artRad(5) + &
                abs((artRad(6)-artgRad6))/artRad(6) + &
                abs((artRad(7)-artgRad7))/artRad(7) + &
                abs((artRad(8)-artgRad8))/artRad(8) + &
                abs((artRad(9)-artgRad9))/artRad(9)

      artRad(1)=artgRad1
      artRad(2)=artgRad2
      artRad(3)=artgRad3
      artRad(4)=artgRad4
      artRad(5)=artgRad5
      artRad(6)=artgRad6
      artRad(7)=artgRad7
      artRad(8)=artgRad8
      artRad(9)=artgRad9

   endif
   
      do i =1,9; nodeadd(i)=nodeadd(i)/dble(mcount); enddo


      
         ydist = abs((nodeadd(1)-nodemeandiff(1))/nodemeandiff(1)) + &
                 abs((nodeadd(2)-nodemeandiff(2))/nodemeandiff(2)) + &
                 abs((nodeadd(3)-nodemeandiff(3))/nodemeandiff(3)) + &
                 abs((nodeadd(4)-nodemeandiff(4))/nodemeandiff(4)) + &
                 abs((nodeadd(5)-nodemeandiff(5))/nodemeandiff(5)) + &
                 abs((nodeadd(6)-nodemeandiff(6))/nodemeandiff(6)) + &
                 abs((nodeadd(7)-nodemeandiff(7))/nodemeandiff(7)) + &
                 abs((nodeadd(8)-nodemeandiff(8))/nodemeandiff(8)) + &
                 abs((nodeadd(9)-nodemeandiff(9))/nodemeandiff(9)) 

            write (*,*) "ydist : ",ydist

            do i = 1,9
               nodemeandiff(i)=nodeadd(i)
            enddo

         call inverse

      endif


      if (lcount == 2 ) then
         lcount = 3
      else if (ydist .le. pydist) then
         lcount = 2
         ncount = 1
         lflag = 1
      endif

      if (ncount == 1) write (*,*) sysTime,jcount,lflag

175   continue
      if (istate .lt. 0) then
      nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (miter .eq. 2) nfea = nfe - neq*nje
      if (miter .eq. 3) nfea = nfe - nje
       write (*,*) nerr,nst,nfe
       write (*,*) nje,lenrw,leniw

  write (*,*) 'ABNORMAL EXIT: ISTATE = ', istate,' Processor 1'
      endif

      if (myid == 0) close(20)
      if (myid == 1) close(21)
      if (myid == 2) close(22)
      if (myid == 3) close(23)
      if (myid == 4) close(24)
      if (myid == 5) close(25)
      if (myid == 6) close(26)
      if (myid == 7) close(27)
      if (myid == 8) close(28)
      if (myid == 9) close(29)
      if (myid == 10) close(30)
      
      end do iterateradii

      open(unit=31,File='runlog.dat',STATUS="old",POSITION="APPEND")
      write(31,*) 'Job number :',jobno
      write(31,*) 'Root Pressure : ',inr
      write(31,*) '      Neph. No.     xydl'
        do i=1,numberOfNephrons;write(31,'(2x,I2,3x,F5.2)') i,xydl(i);enddo

      write(31,*) '      Art. No.     artLen       artRad'
          
         do i =1,numberOfNodes
            write (31,*) i,artLen(i),artRad(i)
         enddo
       close(31)         
     end program opsplit



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function nyid(myid)
	integer :: myid
	character*2 oyid,nyid

       if (myid .lt. 10) then
	  write (oyid,'(i1)') myid
         else
	  write(oyid,'(i2)') myid
	endif

	nyid=oyid

	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine outname(irun,nr)
        use affall
	use tuball

        implicit NONE

	character*3 nr,irun

        if (jobno .lt. 10) then
           write (irun,'(i1)') jobno

          else if (jobno .lt. 100) then
           write (irun,'(i2)') jobno
	   
          elseif (jobno .lt. 1000) then
           write (irun,'(i3)') jobno

         endif 

       if (inr .lt. 100) then
	  write (nr,'(i2)') inr
         else
	  write(nr,'(i3)')inr
	endif

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double precision function rad(x)

        use affall

        implicit NONE
!
        double precision ::  x

            ! calculate radius from y(5) and y(11)

        if (x .le. 6.88D-02) STOP 'Needs cell area adjustment - main'

        rad=0.5D0*(x/xpi - Area/x)

        end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
        function bifuradadd(x1,x2)

          implicit NONE
          real(kind=8) :: bifuradadd
          real(kind=8) :: x1,x2

          bifuradadd=(x1**3+x2**3)**(1.D0/3.D0)

          end function bifuradadd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
        function bifuraddiff(x1,x2)

          implicit NONE
          real(kind=8) :: bifuraddiff
          real(kind=8) :: x1,x2

          bifuraddiff=(x1**3-x2**3)**(1.D0/3.D0)

          end function bifuraddiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          function xnodeelec(i)

            use affall
            
            implicit NONE
            real(kind=8) :: xnodeelec
            integer :: i
            xnodeelec = (potNode(i)-potNode(i+1))*xgi*2.0D0*xpi*artRad(i)/artLen(i)


         end function xnodeelec
