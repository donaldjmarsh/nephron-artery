	MODULE	affall
	

	real(kind=8), parameter :: dyn=1333.22D0
	real(kind=8), parameter :: xp=1.0D0
	real(kind=8), parameter :: xdl=1.3D0
	real(kind=8), parameter :: xal=1.6D0
	real(kind=8), parameter :: xmd=1.80D0
	real(kind=8), parameter :: xd=1.95D0


	real(kind=8), parameter :: x0=1.5D-01
	real(kind=8), parameter :: yx=9.0D-01
	real(kind=8), parameter :: pi=3.141592653589D0
	real(kind=8), parameter :: conv1=7.50061682704D-4
	real(kind=8), parameter :: conv2=1.66666666667D-8
        real(kind=8), parameter :: dz=1.0D-02

	integer :: lout(25),ntime,iprint,ncount
	integer :: jobno,numberOfNephrons,numberOfNodes,jflag,lflag
        integer :: nodeNumber(15),PDnodeNumber(15)

        real(kind=8) :: v1,v2,v4,v5,v6,ca3,ca4,vel,vk,vca
	real(kind=8) :: gel,gk,cond,gca ,gpd,gnef ,g0, xgi
	real(kind=8) :: kd,bt,alpha,kca ,rho 
        real(kind=8) :: qm, caim, psim, xkpsi,cairef
	real(kind=8) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
	real(kind=8) :: u1,u2,u3
	real(kind=8) :: psiref,xomref
	real(kind=8) :: y1,y2,y3,y4 
	real(kind=8) :: ax,bx,cx,dx
        real(kind=8) :: phin(50)
	real(kind=8) :: xpi,Area,we,S,sig0,wm
	real(kind=8) :: tau,sigyref,xnuref,res12
	real(kind=8) :: ap,pg,segres(50,0:3),apress(50,2),yt(50),v7
	real(kind=8) :: startprt,stoptime,prtstart(2),timestop(2)
	real(kind=8) :: res0,affRad(50,2),pressGlom(50),rasAff(10)
	real(kind=8) :: pressArt(15),artPress(15),xpress(15)
	real(kind=8) :: potDiff(25),potNode(15),upn(15),rbf(50)
	real(kind=8) :: vPN1(3),vPN2(3),vPN3(3),vPN4(3),vPN5(3),vPN6(3)
	real(kind=8) :: vPN7(3),vPN8(3),vPN9(3),vPN10(3),vPN11(3),vPN12(3)
	real(kind=8) :: vPN13(3),vPN14(3),vPN15(3),vPN16(3),vPN17(3),vPN18(3)
	real(kind=8) :: vPN19(3),vPN20(3),vPN21(3),vPN22(3),vPN23(3),vPN24(3)
        real(kind=8) :: DM(9,9),CM(9,11),Matrix3(9,11)
        real(kind=8) :: artLen(9),artRad(9)
        real(kind=8) :: g1,g2,g3,g4,g5,g6,g7,g8,g9,gend
        real(kind=8) :: Vec1(11),Vec2(9),Vec3(11)

	real(kind=8) :: zsum,psum 
        real(kind=8) :: odediv(15)
        real(kind=8) :: xintv,suma,sumb

        real(kind=8) :: rootpress,pressNod(25)
        real*8 kr(2)/1.0D-13,0.50D-13 /
        
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE INVERSE


        IMPLICIT NONE

        INTEGER, PARAMETER :: n=9,np=9
        
        REAL(kind=8) :: artArea(np)
        REAL(kind=8) :: A(9,9),d,B(n),w(np,np)
        INTEGER :: i, j,INDX(np), ERRORFLAG

      real(kind=8), parameter :: xpi=3.14159265D0

      real*8 :: master(9,20) 
      
      real*8 :: c1112,c1312,c1512
      real*8 :: c113,c1213,c1413
      real*8 :: c214,c314,c1314
      real*8 :: c415,c1215,c1615
      real*8 :: c516,c1516,c1716
      real*8 :: c1617,c1817
      real*8 :: c618,c718,c1718
      real*8 :: c819,c1719,c1917
      real*8 :: c920,c1020,c1920,c2019

         numberOfNodes = 9


      do i = 1,numberOfNodes
         artArea(i)=xgi*2.0D0*xpi*artRad(i)/artLen(i)
      enddo

         g1 = artArea(1)
         g2 = artArea(2)
         g3 = artArea(3)
         g4 = artArea(4)
         g5 = artArea(5)
         g6 = artArea(6)
         g7 = artArea(7)
         g8 = artArea(8)
         g9 = artArea(9)

         gend = 1.0D0/(((1.0D0/g1) + (1.0D0/g0)))

         c1112=gend;c1312=g5;c1512=g2
         c113=gnef;c1213=g5;c1413=g6
         c214=gnef;c314=gnef;c1314=g6
         c415=gnef;c1215=g2;c1615=g3
         c516=gnef;c1516=g3;c1716=g4
         c1617=g4;c1817=g7;c1917=g8
         c618=gnef;c718=gnef;c1718=g7
         c819=gnef;c1719=g8;c2019=g9
         c920=gnef;c1020=gnef;c1920=g9
                
         
         master=0.0D0

         master(1,11)=c1112;master(1,13)=c1312;master(1,15)=c1512
         master(2,1)=c113;master(2,12)=c1213;master(2,14)=c1413
         master(3,2)=c214;master(3,3)=c314;master(3,13)=c1314
         master(4,4)=c415;master(4,11)=c1215;master(4,16)=c1615
         master(5,5)=c516;master(5,15)=c1516;master(5,17)=c1716
         master(6,16)=c1617;master(6,18)=c1817;master(6,19)=c1917
         master(7,6)=c618;master(7,7)=c718;master(7,17)=c1718
         master(8,8)=c819;master(8,17)=c1719;master(8,20)=c2019
         master(9,9)=c920;master(9,10)=c1020;master(9,19)=c1920
         master(1,12)=-(c1112+c1312+c1512)
         master(2,13)=-(c113+c1213+c1413)
         master(3,14)=-(c214+c314+c1314)
         master(4,15)=-(c415+c1215+c1615)
         master(5,16)=-(c516+c1516+c1716)
         master(6,17)=-(c1617+c1817+c1917)
         master(7,18)=-(c618+c718+c1718)
         master(8,19)=-(c819+c1719+c2019)
         master(9,20)=-(c920+c1020+c1920)
         
         DM=0.0D0

       
       do i =1,numberOfNodes
          do j =1,numberOfNodes
             DM(i,j)=master(i,j+11)
           enddo
       enddo
       
      CM=0.0D0

      do i =1,numberOfNodes
         do j = 1,numberOfNephrons+1
           CM(i,j)=master(i,j)
         enddo
      enddo


        do i =1,numberOfNodes;do j =1,numberOfNodes;A(i,j)=DM(i,j);enddo;enddo

         Matrix3 =0.0D0

        first: do i=1,n
         second: do j = 1,n
          w(i,j)=0.0D0
          enddo second
         w(i,i)=1.D0
         enddo first

       call lucduc(A,N, NP,INDX,D)

        do j=1,n
         call ludkud(A,N,NP,INDX,w(1,J))
         enddo


        Matrix3=MATMUL(w,CM)


        END SUBROUTINE INVERSE

!!!!!!!!!!!!!!!!!!!!                                                                                                    
      SUBROUTINE LUCDUC(A,N,NP,INDX,D)

      implicit NONE

      integer, parameter :: nmax=100, tiny=1.0D-20
      integer :: i,n,np,indx(n),imax,j,k

      real(kind=8) :: A(np,np),vv(nmax), aamax,d, dum,sum

      d=1
      do i=1,n
          aamax=0.D0
          do  j=1,n
              if (abs(A(i,j)).gt.aamax) aamax=abs(A(i,j))
          enddo
          if (aamax.eq.0.D0) then
            write(6,*)(A(i,j),j=1,n)
            write(6,*) i, n
            stop 'Singular matrix - affall '
          endif
          vv(i)=1./aamax
         enddo
       first: do j=1,n
          if (j.gt.1) then
              do i=1,j-1
                  sum=A(i,j)
                  if (i.gt.1) then
                      do k=1,i-1
                          sum=sum-A(i,k)*A(k,j)
                      enddo
                      A(i,j)=sum
                  endif
             enddo
          endif
          aamax=0
          second: do i=j,n
              sum=A(i,j)
              if (j.gt.1) then
                  third: do  k=1,j-1
                      sum=sum-A(i,k)*A(k,j)
                  enddo third
                  A(i,j)=sum
              endif
              dum=vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax=i
                  aamax=dum
             endif
              enddo second
          if (j.ne.imax) then
              fourth: do  k=1,n
                  dum=A(imax,k)
                  A(imax,k)=A(j,k)
                  A(j,k)=dum
              enddo fourth
              d=-d
              vv(imax)=vv(j)
          endif
          indx(j)=imax
          if (j.ne.n) then
              if(A(j,j).eq.0.) A(j,j)=tiny
              dum=1./A(j,j)
              fifth: do i=j+1,n
                  A(i,j)=A(i,j)*dum
              enddo fifth
          endif
          enddo first
      if (A(n,n).eq.0.D0) A(n,n)=tiny
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE LUDKUD(A,N,NP,INDX,B)


      implicit NONE
      
      integer :: ii,ll,i,j,n,np
      integer :: indx(n)
      real(kind=8) :: sum
      real(kind=8) :: A(np,np),b(n)

      ii=0
      first: do i=1,n
          ll=indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if (ii.ne.0) then
              second: do  j=ii,i-1
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

      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end module
