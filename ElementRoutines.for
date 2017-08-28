c
c  The following are all utility routines used
c  specifically in element codes
c
c  Shawn Chester, June 2009
c  Shawn Chester, modified Sept. 2010
c
************************************************************************

      subroutine xint1D1pt(xi,w,nIntPt)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 1D elements
      !  using 1 gauss point for integration as shown in the picture
      !
      !              |---------> xi (=xi_1)
      !        o-----1-----o
      !                    
      !      xi=-1       xi=1
      !
      !  xi(nIntPt,1): xi coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,1), w(1)


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 2.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0

      return
      end subroutine xint1D1pt

************************************************************************

      subroutine xint1D2pt(xi,w,nIntPt)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 1D elements
      !  using 2 gauss points for integration as shown in the picture
      !
      !              |---------> xi (=xi_1)
      !        o--1--|--2--o
      !                     
      !      xi=-1       xi=1
      !
      !  xi(nIntPt,1): xi coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(2,1), w(2)


      ! Number of Gauss points
      !
      nIntPt = 2


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(2,1) =  dsqrt(1.d0/3.d0)

      return
      end subroutine xint1D2pt

************************************************************************

      subroutine xint1D3pt(xi,w,nIntPt)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 1D elements
      !  using 3 gauss points for integration as shown in the picture
      !
      !              |---------> xi (=xi_1)
      !        o--1--2--3--o
      !                     
      !      xi=-1       xi=1
      !
      !  xi(nIntPt,1): xi coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(3,1),w(3)


      ! Number of Gauss points
      !
      nIntPt = 3


      ! Gauss weights
      !
      w(1) = 5.d0/9.d0
      w(2) = 8.d0/9.d0
      w(3) = 5.d0/9.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = -dsqrt(3.d0/5.d0)
      xi(2,1) = 0.d0
      xi(3,1) =  dsqrt(3.d0/5.d0) 

      return
      end subroutine xint1D3pt

************************************************************************

      subroutine xintTri1pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration as shown in the picture
      !
      !              A eta (=xi_2)
      !              |
      !              o
      !              |\
      !              | \
      !              |  \
      !              |   \
      !              |    \
      !              | 1   \
      !              |      \
      !              |       \
      !              1--------o--> xi (=xi_1)
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,2), w(1)


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 1.d0/2.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 1.d0/3.d0
      xi(1,2) = 1.d0/3.d0

      return
      end subroutine xintTri1pt

************************************************************************

      subroutine xintTri3pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 3 gauss point for integration as shown in the picture
      !
      !              A eta (=xi_2)
      !              |
      !              o
      !              |\
      !              | \
      !              |  \
      !              |   \
      !              | 3  \
      !              |     \
      !              | 1  2 \
      !              |       \
      !              1--------o--> xi (=xi_1)
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(3,2), w(3)


      ! Number of Gauss points
      !
      nIntPt = 3


      ! Gauss weights
      !
      w(1) = 1.d0/6.d0
      w(2) = 1.d0/6.d0
      w(3) = 1.d0/6.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 1.d0/6.d0
      xi(1,2) = 1.d0/6.d0
      xi(2,1) = 2.d0/3.d0
      xi(2,2) = 1.d0/6.d0
      xi(3,1) = 1.d0/6.d0
      xi(3,2) = 2.d0/3.d0

      return
      end subroutine xintTri3pt

************************************************************************

      subroutine xintTri3pt_b(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 3 gauss point for integration as shown in the picture
      !
      !              A eta (=xi_2)
      !              |
      !              o
      !              |\
      !              | \
      !              |  \
      !              |   \
      !              |    1
      !              2     \
      !              |      \
      !              |       \
      !              1---3----o--> xi (=xi_1)
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(3,2), w(3)


      ! Number of Gauss points
      !
      nIntPt = 3


      ! Gauss weights
      !
      w(1) = 1.d0/6.d0
      w(2) = 1.d0/6.d0
      w(3) = 1.d0/6.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 1.d0/2.d0
      xi(1,2) = 1.d0/2.d0
      xi(2,1) = 0.d0
      xi(2,2) = 1.d0/6.d0
      xi(3,1) = 1.d0/6.d0
      xi(3,2) = 2.d0/3.d0

      return
      end subroutine xintTri3pt_b

************************************************************************

      subroutine xintTri4pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 3 gauss point for integration as shown in the picture
      !
      !              A eta (=xi_2)
      !              |
      !              o
      !              |\
      !              | \
      !              |  \
      !              |4  \
      !              |    \
      !              |  1  \
      !              |      \
      !              | 2  3  \
      !              1--------o--> xi (=xi_1)
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(4,2), w(4)


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = -27.d0/96.d0
      w(2) = 25.d0/96.d0
      w(3) = 25.d0/96.d0
      w(4) = 25.d0/96.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 1.d0/3.d0
      xi(1,2) = 1.d0/3.d0
      xi(2,1) = 1.d0/5.d0
      xi(2,2) = 1.d0/5.d0
      xi(3,1) = 3.d0/5.d0
      xi(3,2) = 1.d0/5.d0
      xi(4,1) = 1.d0/5.d0
      xi(4,2) = 3.d0/5.d0

      return
      end subroutine xintTri4pt

************************************************************************

      subroutine xintTri6pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 3 gauss point for integration as shown in the picture
      !
      !              A eta (=xi_2)
      !              |
      !              o
      !              |\
      !              | \
      !              |  \
      !              |   \
      !              |    \
      !              |     \
      !              |      \
      !              |       \
      !              1--------o--> xi (=xi_1)
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(6,2), w(6)


      ! Number of Gauss points
      !
      nIntPt = 6


      ! Gauss weights
      !
      w(1) = 0.109951743655322d0
      w(2) = 0.109951743655322d0
      w(3) = 0.109951743655322d0
      w(4) = 0.223381589678011d0
      w(5) = 0.223381589678011d0
      w(6) = 0.223381589678011d0


      ! Gauss pt location in master element
      !
      xi(1,1) = 0.091576213509771d0
      xi(1,2) = 0.091576213509771d0
      xi(2,1) = 0.816847572980459d0
      xi(2,2) = 0.091576213509771d0
      xi(3,1) = 0.091576213509771d0
      xi(3,2) = 0.816847572980459d0
      xi(4,1) = 0.445948490915965d0
      xi(4,2) = 0.108103018168070d0
      xi(5,1) = 0.445948490915965d0
      xi(5,2) = 0.445948490915965d0
      xi(6,1) = 0.108103018168070d0
      xi(6,2) = 0.445948490915965d0



      return
      end subroutine xintTri6pt

************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration as shown in the picture
      !
      !              A eta (=xi_2)
      !              |
      !              |
      !        o-----------o
      !        |     |     |
      !        |     |     |
      !        |     1-----|---> xi (=xi_1)
      !        |           |
      !        |           |
      !        o-----------o
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,2), w(1)


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0

      return
      end subroutine xint2D1pt
      
************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration as shown
      !
      !              A eta (=xi_2)
      !              |
      !              |
      !        o-----------o
      !        |     |     |
      !        |  3  |  4  |
      !        |     ------|---> xi (=xi_1)
      !        |           |
      !        |  1     2  |
      !        o-----------o
      !
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(4,2), w(4)


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)

      return
      end subroutine xint2D4pt

************************************************************************

      subroutine xint2D4ptNodal(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration as shown
      !
      !              A eta (=xi_2)
      !              |
      !              |
      !        3-----------4
      !        |     |     |
      !        |     |     |
      !        |     ------|---> xi (=xi_1)
      !        |           |
      !        |           |
      !        1-----------2
      !
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(4,2), w(4)


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -1.d0
      xi(1,2) = -1.d0
      xi(2,1) = 1.d0
      xi(2,2) = -1.d0
      xi(3,1) = -1.d0
      xi(3,2) = 1.d0
      xi(4,1) = 1.d0
      xi(4,2) = 1.d0

      return
      end subroutine xint2D4ptNodal

************************************************************************

      subroutine xint2D9pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 9 gauss ponits for integration as shown
      !
      !              A eta (=xi_2)
      !              |
      !              |
      !        o-----o-----o
      !        |  7  8  9  |
      !        |     |     |
      !        o  4  5--6--o---> xi (=xi_1)
      !        |           |
      !        |  1  2  3  |
      !        o-----o-----o
      !
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(9,2), w(9)


      ! Number of Gauss points
      !
      nIntPt = 9


      ! Gauss weights
      !
      w(1) = (5.d0/9.d0)*(5.d0/9.d0)
      w(2) = (5.d0/9.d0)*(8.d0/9.d0)
      w(3) = (5.d0/9.d0)*(5.d0/9.d0)
      w(4) = (5.d0/9.d0)*(8.d0/9.d0)
      w(5) = (8.d0/9.d0)*(8.d0/9.d0)
      w(6) = (5.d0/9.d0)*(8.d0/9.d0)
      w(7) = (5.d0/9.d0)*(5.d0/9.d0)
      w(8) = (5.d0/9.d0)*(8.d0/9.d0)
      w(9) = (5.d0/9.d0)*(5.d0/9.d0)
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(3.d0/5.d0)
      xi(1,2) = -dsqrt(3.d0/5.d0)
      xi(2,1) = 0.d0
      xi(2,2) = -dsqrt(3.d0/5.d0)
      xi(3,1) = dsqrt(3.d0/5.d0)
      xi(3,2) = -dsqrt(3.d0/5.d0)
      xi(4,1) = -dsqrt(3.d0/5.d0)
      xi(4,2) = 0.d0
      xi(5,1) = 0.d0
      xi(5,2) = 0.d0
      xi(6,1) = dsqrt(3.d0/5.d0)
      xi(6,2) = 0.d0
      xi(7,1) = -dsqrt(3.d0/5.d0)
      xi(7,2) = dsqrt(3.d0/5.d0)
      xi(8,1) = 0.d0
      xi(8,2) = dsqrt(3.d0/5.d0)
      xi(9,1) = dsqrt(3.d0/5.d0)
      xi(9,2) = dsqrt(3.d0/5.d0)

      return
      end subroutine xint2D9pt
      
************************************************************************

      subroutine xintTet1pt(xi,w,nIntPt)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a single gauss point for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer nIntPt,nDim

      real*8 xi(1,3),w(1)


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 1.d0/6.d0
     

      ! Gauss pt locations in master element
      !
      xi(1,1) = 1.d0/4.d0
      xi(1,2) = 1.d0/4.d0
      xi(1,3) = 1.d0/4.d0

      return
      end subroutine xintTet1pt

************************************************************************

      subroutine xintTet4pt(xi,w,nIntPt)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 4 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer nIntPt,nDim

      real*8 xi(4,3),w(4)


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0/24.d0
      w(2) = 1.d0/24.d0
      w(3) = 1.d0/24.d0
      w(4) = 1.d0/24.d0


      ! Gauss pt locations in master element
      !
      xi(1,1) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(1,2) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(1,3) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(2,1) = (5.d0 + 3.d0*dsqrt(5.d0))/20.d0
      xi(2,2) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(2,3) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(3,1) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(3,2) = (5.d0 + 3.d0*dsqrt(5.d0))/20.d0
      xi(3,3) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(4,1) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(4,2) = (5.d0 - dsqrt(5.d0))/20.d0
      xi(4,3) = (5.d0 + 3.d0*dsqrt(5.d0))/20.d0

      return
      end subroutine xintTet4pt

************************************************************************

      subroutine xintTet5pt(xi,w,nIntPt)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 5 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer nIntPt,nDim

      real*8 xi(5,3),w(5)

      ! Number of Gauss points
      !
      nIntPt = 5


      ! Gauss weights
      !
      w(1) = -2.d0/15.d0
      w(2) = 3.d0/40.d0
      w(3) = 3.d0/40.d0
      w(4) = 3.d0/40.d0
      w(5) = 3.d0/40.d0
     

      ! Gauss pt locations in master element
      !
      xi(1,1) = 1.d0/4.d0
      xi(1,2) = 1.d0/4.d0
      xi(1,3) = 1.d0/4.d0
      xi(2,1) = 1.d0/6.d0
      xi(2,2) = 1.d0/6.d0
      xi(2,3) = 1.d0/6.d0
      xi(3,1) = 1.d0/6.d0
      xi(3,2) = 1.d0/6.d0
      xi(3,3) = 1.d0/2.d0
      xi(4,1) = 1.d0/6.d0
      xi(4,2) = 1.d0/2.d0
      xi(4,3) = 1.d0/6.d0
      xi(5,1) = 1.d0/2.d0
      xi(5,2) = 1.d0/6.d0
      xi(5,3) = 1.d0/6.d0


      return
      end subroutine xintTet5pt

************************************************************************

      subroutine xint3D1pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(1,3),w(1)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0

      return
      end subroutine xint3D1pt

************************************************************************

      subroutine xint3D2pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(2,3),w(2)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 2


      ! Gauss weights
      !
      w(1) = 4.d0
      w(2) = 4.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = 0.d0
      xi(2,2) = 0.d0
      xi(2,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D2pt

************************************************************************

      subroutine xint3D4pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 4 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(4,3),w(4)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 2.d0
      w(2) = 2.d0
      w(3) = 2.d0
      w(4) = 2.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) =  0.d0
      xi(1,2) =  dsqrt(2.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) =  0.d0
      xi(2,2) =  dsqrt(2.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) =  dsqrt(2.d0/3.d0)
      xi(3,2) =  0.d0
      xi(3,3) =  dsqrt(1.d0/3.d0)
      xi(4,1) = -dsqrt(2.d0/3.d0)
      xi(4,2) =  0.d0
      xi(4,3) =  dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D4pt
      
************************************************************************

      subroutine xint3D6pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 6 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(6,3),w(6)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 6


      ! Gauss weights
      !
      w(1) = 4.d0/3.d0
      w(2) = 4.d0/3.d0
      w(3) = 4.d0/3.d0
      w(4) = 4.d0/3.d0
      w(5) = 4.d0/3.d0
      w(6) = 4.d0/3.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 1.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0
      xi(2,1) = -1.d0
      xi(2,2) = 0.d0
      xi(2,3) = 0.d0
      xi(3,1) = 0.d0
      xi(3,2) = 1.d0
      xi(3,3) = 0.d0
      xi(4,1) = 0.d0
      xi(4,2) = -1.d0
      xi(4,3) = 0.d0
      xi(5,1) = 0.d0
      xi(5,2) = 0.d0
      xi(5,3) = 1.d0
      xi(6,1) = 0.d0
      xi(6,2) = 0.d0
      xi(6,3) = -1.d0


      return
      end subroutine xint3D6pt
      
************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(8,3),w(8)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

************************************************************************

      subroutine xint3D12pt(xi,w,nIntPt)


      !
      ! THIS SUBROUTINE IS NOT 100% GOOD YET
      !
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 12 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(12,3),w(12)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 12


      ! Gauss weights
      !
      w(1) = 2.d0/3.d0
      w(2) = 2.d0/3.d0
      w(3) = 2.d0/3.d0
      w(4) = 2.d0/3.d0
      w(5) = 2.d0/3.d0
      w(6) = 2.d0/3.d0
      w(7) = 2.d0/3.d0
      w(8) = 2.d0/3.d0
      w(9) = 2.d0/3.d0
      w(10) = 2.d0/3.d0
      w(11) = 2.d0/3.d0
      w(12) = 2.d0/3.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = -dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = 0.d0
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = 0.d0
      xi(7,1) = dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = 0.d0
      xi(8,1) = -dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = 0.d0
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = dsqrt(1.d0/3.d0)
      xi(3,1) = dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = dsqrt(1.d0/3.d0)
      xi(4,1) = -dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = dsqrt(1.d0/3.d0)



      return
      end subroutine xint3D12pt

************************************************************************

      subroutine xint3D14pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 14 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim,fac

      real*8 xi(14,3),w(14)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 14


      ! Gauss weights
      !
      w(1)  = 320.d0/361.d0
      w(2)  = 320.d0/361.d0
      w(3)  = 320.d0/361.d0
      w(4)  = 320.d0/361.d0
      w(5)  = 320.d0/361.d0
      w(6)  = 320.d0/361.d0
      w(7)  = 121.d0/361.d0
      w(8)  = 121.d0/361.d0
      w(9)  = 121.d0/361.d0
      w(10) = 121.d0/361.d0
      w(11) = 121.d0/361.d0
      w(12) = 121.d0/361.d0
      w(13) = 121.d0/361.d0
      w(14) = 121.d0/361.d0


      ! Gauss pt locations in master element
      !
      xi(1,1)  =  dsqrt(19.d0/30.d0)
      xi(1,2)  =  0.d0
      xi(1,3)  =  0.d0
      xi(2,1)  = -dsqrt(19.d0/30.d0)
      xi(2,2)  =  0.d0
      xi(2,3)  =  0.d0
      xi(3,1)  =  0.d0
      xi(3,2)  =  dsqrt(19.d0/30.d0)
      xi(3,3)  =  0.d0
      xi(4,1)  =  0.d0
      xi(4,2)  = -dsqrt(19.d0/30.d0)
      xi(4,3)  =  0.d0
      xi(5,1)  =  0.d0
      xi(5,2)  =  0.d0
      xi(5,3)  =  dsqrt(19.d0/30.d0)
      xi(6,1)  =  0.d0
      xi(6,2)  =  0.d0
      xi(6,3)  = -dsqrt(19.d0/30.d0)
      xi(7,1)  =  dsqrt(19.d0/33.d0)
      xi(7,2)  =  dsqrt(19.d0/30.d0)
      xi(7,3)  =  dsqrt(19.d0/30.d0)
      xi(8,1)  = -dsqrt(19.d0/30.d0)
      xi(8,2)  =  dsqrt(19.d0/30.d0)
      xi(8,3)  =  dsqrt(19.d0/30.d0)
      xi(9,1)  =  dsqrt(19.d0/30.d0)
      xi(9,2)  = -dsqrt(19.d0/30.d0)
      xi(9,3)  =  dsqrt(19.d0/30.d0)
      xi(10,1) = -dsqrt(19.d0/30.d0)
      xi(10,2) = -dsqrt(19.d0/30.d0)
      xi(10,3) =  dsqrt(19.d0/30.d0)
      xi(11,1) =  dsqrt(19.d0/30.d0)
      xi(11,2) =  dsqrt(19.d0/30.d0)
      xi(11,3) = -dsqrt(19.d0/30.d0)
      xi(12,1) =  dsqrt(19.d0/30.d0)
      xi(12,2) = -dsqrt(19.d0/30.d0)
      xi(12,3) = -dsqrt(19.d0/30.d0)
      xi(13,1) = -dsqrt(19.d0/30.d0)
      xi(13,2) = -dsqrt(19.d0/30.d0)
      xi(13,3) = -dsqrt(19.d0/30.d0)
      xi(14,1) = -dsqrt(19.d0/30.d0)
      xi(14,2) =  dsqrt(19.d0/30.d0)
      xi(14,3) = -dsqrt(19.d0/30.d0)


      return
      end subroutine xint3D14pt
      
************************************************************************

      subroutine xint3D27pt(xi,w,nIntPt)
      
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 27 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      
      implicit none

      integer nIntPt,nDim

      real*8 xi(27,3), w(27)


      ! Number of Gauss points
      !
      nIntPt = 27


      ! Gauss weights
      !
      w(1)  = (5.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(2)  = (8.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(3)  = (5.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(4)  = (5.d0/9.0)*(8.d0/9.0)*(5.d0/9.d0)
      w(5)  = (8.d0/9.0)*(8.d0/9.0)*(5.d0/9.d0)
      w(6)  = (5.d0/9.0)*(8.d0/9.0)*(5.d0/9.d0)
      w(7)  = (5.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(8)  = (8.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(9)  = (5.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(10) = (5.d0/9.0)*(5.d0/9.0)*(8.d0/9.d0)
      w(11) = (8.d0/9.0)*(5.d0/9.0)*(8.d0/9.d0)
      w(12) = (5.d0/9.0)*(5.d0/9.0)*(8.d0/9.d0)
      w(13) = (5.d0/9.0)*(8.d0/9.0)*(8.d0/9.d0)
      w(14) = (8.d0/9.0)*(8.d0/9.0)*(8.d0/9.d0)
      w(15) = (5.d0/9.0)*(8.d0/9.0)*(8.d0/9.d0)
      w(16) = (5.d0/9.0)*(5.d0/9.0)*(8.d0/9.d0)
      w(17) = (8.d0/9.0)*(5.d0/9.0)*(8.d0/9.d0)
      w(18) = (5.d0/9.0)*(5.d0/9.0)*(8.d0/9.d0)
      w(19) = (5.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(20) = (8.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(21) = (5.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(22) = (5.d0/9.0)*(8.d0/9.0)*(5.d0/9.d0)
      w(23) = (8.d0/9.0)*(8.d0/9.0)*(5.d0/9.d0)
      w(24) = (5.d0/9.0)*(8.d0/9.0)*(5.d0/9.d0)
      w(25) = (5.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(26) = (8.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      w(27) = (5.d0/9.0)*(5.d0/9.0)*(5.d0/9.d0)
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(3.d0/5.d0)
      xi(1,2) = -dsqrt(3.d0/5.d0)
      xi(1,3) = -dsqrt(3.d0/5.d0)
      xi(2,1) = 0.d0
      xi(2,2) = -dsqrt(3.d0/5.d0)
      xi(2,3) = -dsqrt(3.d0/5.d0)
      xi(3,1) = dsqrt(3.d0/5.d0)
      xi(3,2) = -dsqrt(3.d0/5.d0)
      xi(3,3) = -dsqrt(3.d0/5.d0)
      xi(4,1) = -dsqrt(3.d0/5.d0)
      xi(4,2) = 0.d0
      xi(4,3) = -dsqrt(3.d0/5.d0)
      xi(5,1) = 0.d0
      xi(5,2) = 0.d0
      xi(5,3) = -dsqrt(3.d0/5.d0)
      xi(6,1) = dsqrt(3.d0/5.d0)
      xi(6,2) = 0.d0
      xi(6,3) = -dsqrt(3.d0/5.d0)
      xi(7,1) = -dsqrt(3.d0/5.d0)
      xi(7,2) = dsqrt(3.d0/5.d0)
      xi(7,3) = -dsqrt(3.d0/5.d0)
      xi(8,1) = 0.d0
      xi(8,2) = dsqrt(3.d0/5.d0)
      xi(8,3) = -dsqrt(3.d0/5.d0)
      xi(9,1) = dsqrt(3.d0/5.d0)
      xi(9,2) = dsqrt(3.d0/5.d0)
      xi(9,3) = -dsqrt(3.d0/5.d0)
      xi(10,1) = -dsqrt(3.d0/5.d0)
      xi(10,2) = -dsqrt(3.d0/5.d0)
      xi(10,3) = 0.d0
      xi(11,1) = 0.d0
      xi(11,2) = -dsqrt(3.d0/5.d0)
      xi(11,3) = 0.d0
      xi(12,1) = dsqrt(3.d0/5.d0)
      xi(12,2) = -dsqrt(3.d0/5.d0)
      xi(12,3) = 0.d0
      xi(13,1) = -dsqrt(3.d0/5.d0)
      xi(13,2) = 0.d0
      xi(13,3) = 0.d0
      xi(14,1) = 0.d0
      xi(14,2) = 0.d0
      xi(14,3) = 0.d0
      xi(15,1) = dsqrt(3.d0/5.d0)
      xi(15,2) = 0.d0
      xi(15,3) = 0.d0
      xi(16,1) = -dsqrt(3.d0/5.d0)
      xi(16,2) = dsqrt(3.d0/5.d0)
      xi(16,3) = 0.d0
      xi(17,1) = 0.d0
      xi(17,2) = dsqrt(3.d0/5.d0)
      xi(17,3) = 0.d0
      xi(18,1) = dsqrt(3.d0/5.d0)
      xi(18,2) = dsqrt(3.d0/5.d0)
      xi(18,3) = 0.d0
      xi(19,1) = -dsqrt(3.d0/5.d0)
      xi(19,2) = -dsqrt(3.d0/5.d0)
      xi(19,3) = dsqrt(3.d0/5.d0)
      xi(20,1) = 0.d0
      xi(20,2) = -dsqrt(3.d0/5.d0)
      xi(20,3) = dsqrt(3.d0/5.d0)
      xi(21,1) = dsqrt(3.d0/5.d0)
      xi(21,2) = -dsqrt(3.d0/5.d0)
      xi(21,3) = dsqrt(3.d0/5.d0)
      xi(22,1) = -dsqrt(3.d0/5.d0)
      xi(22,2) = 0.d0
      xi(22,3) = dsqrt(3.d0/5.d0)
      xi(23,1) = 0.d0
      xi(23,2) = 0.d0
      xi(23,3) = dsqrt(3.d0/5.d0)
      xi(24,1) = dsqrt(3.d0/5.d0)
      xi(24,2) = 0.d0
      xi(24,3) = dsqrt(3.d0/5.d0)
      xi(25,1) = -dsqrt(3.d0/5.d0)
      xi(25,2) = dsqrt(3.d0/5.d0)
      xi(25,3) = dsqrt(3.d0/5.d0)
      xi(26,1) = 0.d0
      xi(26,2) = dsqrt(3.d0/5.d0)
      xi(26,3) = dsqrt(3.d0/5.d0)
      xi(27,1) = dsqrt(3.d0/5.d0)
      xi(27,2) = dsqrt(3.d0/5.d0)
      xi(27,3) = dsqrt(3.d0/5.d0)


      return
      end subroutine xint3D27pt
      
************************************************************************

      subroutine calcShape1DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 2-node linear 1D element as shown
      !
      !         |---------> xi
      !   1-----------2
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,1),sh(2),dshxi(2,1),xi

      real*8 zero,one,half
      parameter(zero=0.d0,one=1.d0,half=1.d0/2.d0)
      

      ! Calculate the shape functions and derivatives:
      !
      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      !
      ! The shape functions
      !
      sh(1) = half*(one - xi)
      sh(2) = half*(one + xi)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -half
      dshxi(2,1) =  half


      return
      end subroutine calcShape1DLinear

************************************************************************

      subroutine calcShape1DQuad(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 2-node linear 1D element as shown
      !
      !         |---------> xi
      !   1-----3-----2
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,1),sh(3),dshxi(3,1),xi

      real*8 zero,one,half,two
      parameter(zero=0.d0,one=1.d0,half=1.d0/2.d0,two=2.d0)
      

      ! Calculate the shape functions and derivatives:
      !
      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      !
      ! The shape functions
      !
      sh(1) = half*xi*(one - xi)
      sh(2) = half*xi*(one + xi)
      sh(3) = one - xi*xi
      !
      ! The first derivatives
      !
      dshxi(1,1) = xi - half
      dshxi(2,1) = xi + half
      dshxi(3,1) = -two*xi


      return
      end subroutine calcShape1DQuad

************************************************************************

      subroutine calcShapeTriLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 3-node linear 2D element as shown
      !
      !              A eta (=xi_2)
      !              |
      !              3
      !              |\
      !              | \
      !              |  \
      !              |   \
      !              |    \
      !              |     \
      !              |      \
      !              |       \
      !              1--------2--> xi (=xi_1)
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,2),sh(3),dshxi(3,2),d2shxi(3,2,2),xi,eta

      real*8 zero,one,four,two
      parameter(zero=0.d0,one=1.d0,two=2.d0,four=4.d0)
      

      ! Calculate the shape functions and derivatives:
      !
      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      !
      ! The shape functions
      !
      sh(1) = one - xi - eta
      sh(2) = xi
      sh(3) = eta
      !
      ! The first derivatives
      !
      dshxi(1,1) = -one
      dshxi(1,2) = -one
      dshxi(2,1) = one
      dshxi(2,2) = zero
      dshxi(3,1) = zero
      dshxi(3,2) = one
      !
      ! The second derivatives
      !
      d2shxi = zero


      return
      end subroutine calcShapeTriLinear

************************************************************************

      subroutine calcShapeTriQuad(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 6-node quadratic 2D element as shown
      !
      !              A eta (=xi_2)
      !              |
      !              3
      !              |\
      !              | \
      !              |  \
      !              |   5
      !              6    \
      !              |     \
      !              |      \
      !              |       \
      !              1---4----2--> xi (=xi_1)
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,2),sh(6),dshxi(6,2),d2shxi(6,2,2),xi,eta
      real*8 lam

      real*8 zero,one,four,two
      parameter(zero=0.d0,one=1.d0,two=2.d0,four=4.d0)
      

      ! Calculate the shape functions and derivatives:
      !
      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      !
      ! lam is just a useful factor
      !
      lam = one - xi - eta
      !
      ! The shape functions
      !
      sh(1) = lam*(two*lam - one)
      sh(2) = xi*(two*xi - one)
      sh(3) = eta*(two*eta - one)
      sh(4) = four*xi*lam
      sh(5) = four*xi*eta
      sh(6) = four*eta*lam
      !
      ! The first derivatives
      !
      dshxi(1,1) = one - four*lam
      dshxi(1,2) = one - four*lam
      dshxi(2,1) = four*xi - one
      dshxi(2,2) = zero
      dshxi(3,1) = zero
      dshxi(3,2) = four*eta - one
      dshxi(4,1) = four*(lam - xi)
      dshxi(4,2) = -four*xi
      dshxi(5,1) = four*eta
      dshxi(5,2) = four*xi
      dshxi(6,1) = -four*eta
      dshxi(6,2) = four*(lam - eta)
      !
      ! The second derivatives
      !
      d2shxi = zero


      return
      end subroutine calcShapeTriQuad

************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 4-node linear 2D element as shown
      !
      !         A eta
      !         |
      !         |
      !   4-----------3
      !   |     |     |
      !   |     |     |
      !   |     ------|---> xi
      !   |           |
      !   |           |
      !   1-----------2
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),d2shxi(4,2,2),xi,eta

      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Calculate the shape functions and derivatives:
      !
      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      !
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = fourth
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(2,1,2) = -fourth
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,2) = fourth
      d2shxi(3,2,1) = d2shxi(3,1,2)
      d2shxi(4,1,2) = -fourth
      d2shxi(4,2,1) = d2shxi(4,1,2)

      return
      end subroutine calcShape2DLinear

************************************************************************

      subroutine calcShape2DQuad(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node quadratic 2D element as shown
      !
      !         A eta
      !         |
      !         |
      !   4-----7-----3
      !   |     |     |
      !   |     |     |
      !   8     ------6---> xi
      !   |           |
      !   |           |
      !   1-----5-----2
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt,i

      real*8 xi_int(nIntPt,2),sh(8),dshxi(8,2)
      real*8 d2shxi(8,2,2),xi,eta

      real*8 zero,one,two,half,fourth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      !
      ! The shape functions
      !
      sh(1) = -fourth*(one - xi)*(one - eta)*(one + xi + eta)
      sh(2) = -fourth*(one + xi)*(one - eta)*(one - xi + eta)
      sh(3) = -fourth*(one + xi)*(one + eta)*(one - xi - eta)
      sh(4) = -fourth*(one - xi)*(one + eta)*(one + xi - eta)
      sh(5) = half*(one - xi)*(one + xi)*(one - eta)
      sh(6) = half*(one - eta)*(one + eta)*(one + xi)
      sh(7) = half*(one - xi)*(one + xi)*(one + eta)
      sh(8) = half*(one - eta)*(one + eta)*(one - xi)
      !
      ! The first derivatives
      !
      dshxi(1,1) = fourth*(one - eta)*(two*xi + eta)
      dshxi(1,2) = fourth*(one - xi)*(two*eta + xi)
      dshxi(2,1) = fourth*(one - eta)*(two*xi - eta)
      dshxi(2,2) = fourth*(one + xi)*(two*eta - xi)
      dshxi(3,1) = fourth*(one + eta)*(two*xi + eta)
      dshxi(3,2) = fourth*(one + xi)*(two*eta + xi)
      dshxi(4,1) = fourth*(one + eta)*(two*xi - eta)
      dshxi(4,2) = fourth*(one - xi)*(two*eta - xi)
      dshxi(5,1) = -xi*(one - eta)
      dshxi(5,2) = -half*(one - xi)*(one + xi)
      dshxi(6,1) = half*(one - eta)*(one + eta)
      dshxi(6,2) = -eta*(one + xi)
      dshxi(7,1) = -xi*(one + eta)
      dshxi(7,2) = half*(one - xi)*(one + xi)
      dshxi(8,1) = -half*(one - eta)*(one + eta)
      dshxi(8,2) = -eta*(one - xi)
      !
      ! The second derivatives
      !
      d2shxi(1,1,1) = half*(one - eta)
      d2shxi(1,2,2) = half*(one - xi)
      d2shxi(1,1,2) = fourth*(one - two*xi - two*eta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(2,1,1) = half*(one - eta)
      d2shxi(2,2,2) = half*(one + xi)
      d2shxi(2,1,2) = fourth*(-one - two*xi + two*eta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,1) = half*(one + eta)
      d2shxi(3,2,2) = half*(one + xi)
      d2shxi(3,1,2) = fourth*(one + two*xi + two*eta)
      d2shxi(3,2,1) = d2shxi(3,1,2)
      d2shxi(4,1,1) = half*(one + eta)
      d2shxi(4,2,2) = half*(one - xi)
      d2shxi(4,1,2) = fourth*(-one + two*xi - two*eta)
      d2shxi(4,2,1) = d2shxi(4,1,2)
      d2shxi(5,1,1) = eta - one
      d2shxi(5,2,2) = zero
      d2shxi(5,1,2) = xi
      d2shxi(5,2,1) = d2shxi(5,1,2)
      d2shxi(6,1,1) = zero
      d2shxi(6,2,2) = -(one + xi)
      d2shxi(6,1,2) = -eta
      d2shxi(6,2,1) = d2shxi(6,1,2)
      d2shxi(7,1,1) = -(one + eta)
      d2shxi(7,2,2) = zero
      d2shxi(7,1,2) = -xi
      d2shxi(7,2,1) = d2shxi(7,1,2)
      d2shxi(8,1,1) = zero
      d2shxi(8,2,2) = xi - one
      d2shxi(8,1,2) = eta
      d2shxi(8,2,1) = d2shxi(8,1,2)

      return
      end subroutine calcShape2DQuad

************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |       
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt,i,j

      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3)
      real*8 d2shxi(8,3,3),xi,eta,zeta

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)
      
      return
      end subroutine calcShape3DLinear
            
*************************************************************************

      subroutine calcShape3DQuad(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 20-node quadratic 3D element as shown
      !
      !
      !
      !         8------15-------7
      !        /|              /|
      !       / |             / |
      !      16 |            14 |
      !     /   20          /   19   zeta
      !    /    |          /    |    
      !   5-------13------6     |    |     eta
      !   |     |         |     |    |     /
      !   |     4------11-|-----3    |    /
      !   |    /          |    /     |   /
      !   17  /           18  /      |  /
      !   |  12           |  10      | /
      !   | /             | /        |/
      !   |/              |/         O--------- xi
      !   1-------9-------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,3),sh(20),dshxi(20,3)
      real*8 d2shxi(20,3,3),xi,eta,zeta

      real*8 zero,one,two,half,fourth,eighth,three
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0,three=3.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1)  = -eighth*(one-xi)*(one-eta)*(one-zeta)*(two+xi+eta+zeta)
      sh(2)  = -eighth*(one+xi)*(one-eta)*(one-zeta)*(two-xi+eta+zeta)
      sh(3)  = -eighth*(one+xi)*(one+eta)*(one-zeta)*(two-xi-eta+zeta)
      sh(4)  = -eighth*(one-xi)*(one+eta)*(one-zeta)*(two+xi-eta+zeta)
      sh(5)  = -eighth*(one-xi)*(one-eta)*(one+zeta)*(two+xi+eta-zeta)
      sh(6)  = -eighth*(one+xi)*(one-eta)*(one+zeta)*(two-xi+eta-zeta)
      sh(7)  = -eighth*(one+xi)*(one+eta)*(one+zeta)*(two-xi-eta-zeta)
      sh(8)  = -eighth*(one-xi)*(one+eta)*(one+zeta)*(two+xi-eta-zeta)
      sh(9)  = fourth*(one-xi)*(one+xi)*(one-eta)*(one-zeta)
      sh(10) = fourth*(one+xi)*(one+eta)*(one-eta)*(one-zeta)
      sh(11) = fourth*(one+xi)*(one-xi)*(one+eta)*(one-zeta)
      sh(12) = fourth*(one-xi)*(one+eta)*(one-eta)*(one-zeta)
      sh(13) = fourth*(one-xi)*(one+xi)*(one-eta)*(one+zeta)
      sh(14) = fourth*(one+xi)*(one-eta)*(one+eta)*(one+zeta)
      sh(15) = fourth*(one+xi)*(one-xi)*(one+eta)*(one+zeta)
      sh(16) = fourth*(one-xi)*(one+eta)*(one-eta)*(one+zeta)
      sh(17) = fourth*(one-xi)*(one-eta)*(one-zeta)*(one+zeta)
      sh(18) = fourth*(one+xi)*(one-eta)*(one-zeta)*(one+zeta)
      sh(19) = fourth*(one+xi)*(one+eta)*(one-zeta)*(one+zeta)
      sh(20) = fourth*(one-xi)*(one+eta)*(one-zeta)*(one+zeta)
      !
      ! The first derivatives
      !
      dshxi(1,1)=eighth*(zeta-one)*(eta-one)*(eta+zeta+two*xi+one)
      dshxi(1,2)=eighth*(zeta-one)*(xi-one)*(xi+two*eta+zeta+one)
      dshxi(1,3)=eighth*(eta-one)*(xi-one)*(xi+eta+two*zeta+one)
      dshxi(2,1)=-eighth*(zeta-one)*(eta-one)*(eta+zeta-two*xi+one)
      dshxi(2,2)=eighth*(zeta-one)*(one+xi)*(xi-zeta-two*eta-one)
      dshxi(2,3)=eighth*(eta-one)*(one+xi)*(xi-eta-two*zeta-one)
      dshxi(3,1)=-eighth*(zeta-one)*(eta+one)*(eta-zeta+two*xi-one)
      dshxi(3,2)=-eighth*(zeta-one)*(one+xi)*(xi-zeta+two*eta-one)
      dshxi(3,3)=-eighth*(eta+one)*(one+xi)*(xi+eta-two*zeta-one)
      dshxi(4,1)=eighth*(zeta-one)*(eta+one)*(eta-zeta-two*xi-one)
      dshxi(4,2)=-eighth*(zeta-one)*(xi-one)*(xi+zeta-two*eta+one)
      dshxi(4,3)=-eighth*(eta+one)*(xi-one)*(xi-eta+two*zeta+one)
      dshxi(5,1)=-eighth*(zeta+one)*(eta-one)*(eta-zeta+two*xi+one)
      dshxi(5,2)=-eighth*(zeta+one)*(xi-one)*(xi-zeta+two*eta+one)
      dshxi(5,3)=-eighth*(eta-one)*(xi-one)*(xi+eta-two*zeta+one)
      dshxi(6,1)=eighth*(zeta+one)*(eta-one)*(eta-zeta-two*xi+one)
      dshxi(6,2)=-eighth*(zeta+one)*(xi+one)*(xi+zeta-two*eta-one)
      dshxi(6,3)=-eighth*(eta-one)*(xi+one)*(xi-eta+two*zeta-one)
      dshxi(7,1)=eighth*(zeta+one)*(eta+one)*(eta+zeta+two*xi-one)
      dshxi(7,2)=eighth*(zeta+one)*(xi+one)*(xi+zeta+two*eta-one)
      dshxi(7,3)=eighth*(eta+one)*(xi+one)*(xi+eta+two*zeta-one)
      dshxi(8,1)=-eighth*(zeta+one)*(eta+one)*(eta+zeta-two*xi-one)
      dshxi(8,2)=eighth*(zeta+one)*(xi-one)*(xi-zeta-two*eta+one)
      dshxi(8,3)=eighth*(eta+one)*(xi-one)*(xi-eta-two*zeta+one)
      dshxi(9,1)=-half*xi*(zeta-one)*(eta-one)
      dshxi(9,2)=-fourth*(xi-one)*(xi+one)*(zeta-one)
      dshxi(9,3)=-fourth*(xi-one)*(xi+one)*(eta-one)
      dshxi(10,1)=fourth*(eta-one)*(eta+one)*(zeta-one)
      dshxi(10,2)=half*eta*(xi+one)*(zeta-one)
      dshxi(10,3)=fourth*(eta-one)*(eta+one)*(xi+one)
      dshxi(11,1)=half*xi*(zeta-one)*(eta+one)
      dshxi(11,2)=fourth*(xi-one)*(xi+one)*(zeta-one)
      dshxi(11,3)=fourth*(xi-one)*(xi+one)*(eta+one)
      dshxi(12,1)=-fourth*(eta -one)*(eta+one)*(zeta-one)
      dshxi(12,2)=-half*eta*(zeta-one)*(xi-one)
      dshxi(12,3)=-fourth*(eta-one)*(eta+one)*(xi-one)
      dshxi(13,1)=half*xi*(zeta+one)*(eta-one)
      dshxi(13,2)=fourth*(xi-one)*(xi+one)*(zeta+one)
      dshxi(13,3)=fourth*(xi-one)*(xi+one)*(eta-one)
      dshxi(14,1)=-fourth*(eta-one)*(eta+one)*(zeta+one)
      dshxi(14,2)=-half*eta*(zeta+one)*(xi+one)
      dshxi(14,3)=-fourth*(eta-one)*(eta+one)*(xi+one)
      dshxi(15,1)=-half*xi*(eta+one)*(zeta+one)
      dshxi(15,2)=-fourth*(xi-one)*(xi+one)*(zeta+one)
      dshxi(15,3)=-fourth*(xi-one)*(xi+one)*(eta+one)
      dshxi(16,1)=fourth*(eta-one)*(eta+one)*(zeta+one)
      dshxi(16,2)=half*eta*(zeta+one)*(xi-one)
      dshxi(16,3)=fourth*(eta-one)*(eta+one)*(xi-one)
      dshxi(17,1)=-fourth*(zeta-one)*(zeta+one)*(eta-one)
      dshxi(17,2)=-fourth*(zeta-one)*(zeta+one)*(xi-one)
      dshxi(17,3)=-half*zeta*(eta-one)*(xi-one)
      dshxi(18,1)=fourth*(zeta-one)*(zeta+one)*(eta-one)
      dshxi(18,2)=fourth*(zeta-one)*(zeta+one)*(xi+one)
      dshxi(18,3)=half*zeta*(eta-one)*(xi+one)
      dshxi(19,1)=-fourth*(zeta-one)*(zeta+one)*(eta+one)
      dshxi(19,2)=-fourth*(zeta-one)*(zeta+one)*(xi+one)
      dshxi(19,3)=-half*zeta*(eta+one)*(xi+one)
      dshxi(20,1)=fourth*(zeta-one)*(zeta+one)*(eta+one)
      dshxi(20,2)=fourth*(zeta-one)*(zeta+one)*(xi-one)
      dshxi(20,3)=half*zeta*(eta+one)*(xi-one)
      !
      ! The second derivatives
      !
      !  if you need them, you can compute them, good luck
      !
      d2shxi = zero
      
      return
      end subroutine calcShape3DQuad


************************************************************************

      subroutine calcShapeTetLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 4-node linear 3D element as shown
      !
      !             4.              zeta
      !           / |  .             |
      !          /  |    .           |
      !         /   |      .         |
      !        /    |        .       |
      !       /    .1---------3      o----------- eta
      !      /   .        . `       .
      !     /  .      . `         . (origin at node 1)
      !    / .    . `           .
      !    2 .  `             xi
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,3),sh(4),dshxi(4,3),d2shxi(4,3,3),xi,eta,zeta

      real*8 zero,one,four,two
      parameter(zero=0.d0,one=1.d0,two=2.d0,four=4.d0)
      

      ! Calculate the shape functions and derivatives:
      !
      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = one - xi - eta - zeta
      sh(2) = xi
      sh(3) = eta
      sh(4) = zeta
      !
      ! The first derivatives
      !
      dshxi(1,1) = -one
      dshxi(1,2) = -one
      dshxi(1,3) = -one
      dshxi(2,1) = one
      dshxi(2,2) = zero
      dshxi(2,3) = zero
      dshxi(3,1) = zero
      dshxi(3,2) = one
      dshxi(3,3) = zero
      dshxi(4,1) = zero
      dshxi(4,2) = zero
      dshxi(4,3) = one
      !
      ! The second derivatives
      !
      d2shxi = zero


      return
      end subroutine calcShapeTetLinear

************************************************************************

      subroutine calcShapeTetQuad(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 10-node quadratic 3D element as shown
      !
      !             4.              zeta
      !           / |  .             |
      !          /  |    .           |
      !         /   8     10         |
      !        9    |        .       |
      !       /    .1----7----3      o----------- eta
      !      /   .        . `       .
      !     /  5      6 `         . (origin at node 1)
      !    / .    . `           .
      !    2 .  `             xi
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,3),sh(10),dshxi(10,3),d2shxi(10,3,3),xi,eta
      real*8 lam,zeta

      real*8 zero,one,four,two,three
      parameter(zero=0.d0,one=1.d0,two=2.d0,four=4.d0,three=3.d0)
      

      ! Calculate the shape functions and derivatives:
      !
      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! lam is just a useful factor
      !
      lam = one - xi - eta - zeta
      !
      ! The shape functions
      !
      sh(1)  = lam*(two*lam - one)
      sh(2)  = xi*(two*xi - one)
      sh(3)  = eta*(two*eta - one)
      sh(4)  = zeta*(two*zeta - one)
      sh(5)  = four*xi*lam
      sh(6)  = four*xi*eta
      sh(7)  = four*eta*lam
      sh(8)  = four*zeta*lam
      sh(9)  = four*xi*zeta
      sh(10) = four*eta*zeta
      !
      ! The first derivatives
      !
      dshxi(1,1) = four*(xi + eta + zeta) - three
      dshxi(1,2) = four*(xi + eta + zeta) - three
      dshxi(1,3) = four*(xi + eta + zeta) - three
      dshxi(2,1) = four*xi - one
      dshxi(2,2) = zero
      dshxi(2,3) = zero
      dshxi(3,1) = zero
      dshxi(3,2) = four*eta - one
      dshxi(3,3) = zero
      dshxi(4,1) = zero
      dshxi(4,2) = zero
      dshxi(4,3) = four*zeta - one
      dshxi(5,1) = -four*(-one + two*xi + eta + zeta)
      dshxi(5,2) = -four*xi
      dshxi(5,3) = -four*xi
      dshxi(6,1) = four*eta
      dshxi(6,2) = four*xi
      dshxi(6,3) = zero
      dshxi(7,1) = -four*eta
      dshxi(7,2) = -four*(-one + xi + two*eta + zeta)
      dshxi(7,3) = -four*eta
      dshxi(8,1) = -four*zeta
      dshxi(8,2) = -four*zeta
      dshxi(8,3) = -four*(-one + xi + eta + two*zeta)
      dshxi(9,1) = four*zeta
      dshxi(9,2) = zero
      dshxi(9,3) = four*xi
      dshxi(10,1) = zero
      dshxi(10,2) = four*zeta
      dshxi(10,3) = four*eta
      !
      ! The second derivatives
      !
      d2shxi = zero


      return
      end subroutine calcShapeTetQuad

************************************************************************

      subroutine calcShape2DQuadFull(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 9-node quadratic 2D element as shown
      !
      !         A eta
      !         |
      !         |
      !   4-----7-----3
      !   |     |     |
      !   |     |     |
      !   8     9-----6---> xi
      !   |           |
      !   |           |
      !   1-----5-----2
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,2),sh(9),dshxi(9,2)
      real*8 xi,eta,xi2,eta2

      real*8 zero,one,two,half,fourth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      xi2 = xi**2
      eta2 = eta**2
      !
      ! The shape functions
      !
      sh(1) = fourth*(xi - one)*(eta - one)*xi*eta
      sh(2) = fourth*(xi + one)*(eta - one)*xi*eta
      sh(3) = fourth*(xi + one)*(eta + one)*xi*eta
      sh(4) = fourth*(xi - one)*(eta + one)*xi*eta
      sh(5) = half*eta*(one - xi2)*(eta - one)
      sh(6) = half*xi*(one - eta2)*(xi + one)
      sh(7) = half*eta*(one - xi2)*(eta + one)
      sh(8) = half*xi*(one - eta2)*(xi - one)
      sh(9) = (one - xi2)*(one - eta2)
      !
      ! The first derivatives
      !
      dshxi(1,1) = fourth*(two*xi - one)*eta*(eta - one)
      dshxi(1,2) = fourth*xi*(xi - one)*(two*eta - one)
      dshxi(2,1) = fourth*(two*xi + one)*eta*(eta - one)
      dshxi(2,2) = fourth*xi*(xi + one)*(two*eta - one)
      dshxi(3,1) = fourth*(two*xi + one)*eta*(eta + one)
      dshxi(3,2) = fourth*xi*(xi + one)*(two*eta + one)
      dshxi(4,1) = fourth*(two*xi - one)*eta*(eta + one)
      dshxi(4,2) = fourth*xi*(xi - one)*(two*eta + one)
      dshxi(5,1) = -xi*eta*(eta - one)
      dshxi(5,2) = half*(one - xi2)*(two*eta - one)
      dshxi(6,1) = half*(two*xi + one)*(one - eta2)
      dshxi(6,2) = -eta*xi*(xi + one)
      dshxi(7,1) = -xi*eta*(eta + one)
      dshxi(7,2) = half*(one - xi2)*(two*eta + one)
      dshxi(8,1) = half*(two*xi - one)*(one - eta2)
      dshxi(8,2) = -eta*xi*(xi - one)
      dshxi(9,1) = -two*xi*(one - eta2)
      dshxi(9,2) = -two*(one - xi2)*eta
      
      return
      end subroutine calcShape2DQuadFull

************************************************************************

      subroutine calcShape3DQuadFull(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses the 27-node quadratic 2D element 
      ! The first 20 nodes are in the same convention as the 20 node
      ! brick in Abaqus (see calcShape3DQuad() above). The next
      ! 6 nodes are on the faces, in the order of the face
      ! numbering convention in Abaqus (see the manual, vol 4, ch 28).
      ! The last node is at the centroid.
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt

      real*8 xi_int(nIntPt,3),sh(27),dshxi(27,3)
      real*8 xi,eta,zeta
      real*8 Nx(3),Ny(3),Nz(3)
      real*8 dNx(3),dNy(3),dNz(3)

      real*8 one,two,half
      parameter(one=1.d0,two=2.d0,half=0.5d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      
      Nx(1) = half*xi*(xi - one)
      Nx(2) = half*xi*(xi + one)
      Nx(3) = one - xi**2
      
      dNx(1) = xi - half
      dNx(2) = xi + half
      dNx(3) = -two*xi

      Ny(1) = half*eta*(eta - one)
      Ny(2) = half*eta*(eta + one)
      Ny(3) = one - eta**2

      dNy(1) = eta - half
      dNy(2) = eta + half
      dNy(3) = -two*eta

      Nz(1) = half*zeta*(zeta - one)
      Nz(2) = half*zeta*(zeta + one)
      Nz(3) = one - zeta**2

      dNz(1) = zeta - half
      dNz(2) = zeta + half
      dNz(3) = -two*zeta

      !
      ! The shape functions
      !
      sh( 1) = Nx(1)*Ny(1)*Nz(1)
      sh( 2) = Nx(2)*Ny(1)*Nz(1)
      sh( 3) = Nx(2)*Ny(2)*Nz(1)
      sh( 4) = Nx(1)*Ny(2)*Nz(1)
      sh( 5) = Nx(1)*Ny(1)*Nz(2)
      sh( 6) = Nx(2)*Ny(1)*Nz(2)
      sh( 7) = Nx(2)*Ny(2)*Nz(2)
      sh( 8) = Nx(1)*Ny(2)*Nz(2)
      sh( 9) = Nx(3)*Ny(1)*Nz(1)
      sh(10) = Nx(2)*Ny(3)*Nz(1)
      sh(11) = Nx(3)*Ny(2)*Nz(1)
      sh(12) = Nx(1)*Ny(3)*Nz(1)
      sh(13) = Nx(3)*Ny(1)*Nz(2)
      sh(14) = Nx(2)*Ny(3)*Nz(2)
      sh(15) = Nx(3)*Ny(2)*Nz(2)
      sh(16) = Nx(1)*Ny(3)*Nz(2)
      sh(17) = Nx(1)*Ny(1)*Nz(3)
      sh(18) = Nx(2)*Ny(1)*Nz(3)
      sh(19) = Nx(2)*Ny(2)*Nz(3)
      sh(20) = Nx(1)*Ny(2)*Nz(3)
      sh(21) = Nx(3)*Ny(3)*Nz(3)
      sh(22) = Nx(3)*Ny(3)*Nz(1)
      sh(23) = Nx(3)*Ny(3)*Nz(2)
      sh(24) = Nx(3)*Ny(1)*Nz(3)
      sh(25) = Nx(2)*Ny(3)*Nz(3)
      sh(26) = Nx(3)*Ny(2)*Nz(3)
      sh(27) = Nx(1)*Ny(3)*Nz(3)
      !
      ! The first derivatives
      !
      dshxi( 1,1) = dNx(1)*Ny(1)*Nz(1)
      dshxi( 1,2) = Nx(1)*dNy(1)*Nz(1)
      dshxi( 1,3) = Nx(1)*Ny(1)*dNz(1)
      dshxi( 2,1) = dNx(2)*Ny(1)*Nz(1)
      dshxi( 2,2) = Nx(2)*dNy(1)*Nz(1)
      dshxi( 2,3) = Nx(2)*Ny(1)*dNz(1)
      dshxi( 3,1) = dNx(2)*Ny(2)*Nz(1)
      dshxi( 3,2) = Nx(2)*dNy(2)*Nz(1)
      dshxi( 3,3) = Nx(2)*Ny(2)*dNz(1)
      dshxi( 4,1) = dNx(1)*Ny(2)*Nz(1)
      dshxi( 4,2) = Nx(1)*dNy(2)*Nz(1)
      dshxi( 4,3) = Nx(1)*Ny(2)*dNz(1)
      dshxi( 5,1) = dNx(1)*Ny(1)*Nz(2)
      dshxi( 5,2) = Nx(1)*dNy(1)*Nz(2)
      dshxi( 5,3) = Nx(1)*Ny(1)*dNz(2)
      dshxi( 6,1) = dNx(2)*Ny(1)*Nz(2)
      dshxi( 6,2) = Nx(2)*dNy(1)*Nz(2)
      dshxi( 6,3) = Nx(2)*Ny(1)*dNz(2)
      dshxi( 7,1) = dNx(2)*Ny(2)*Nz(2)
      dshxi( 7,2) = Nx(2)*dNy(2)*Nz(2)
      dshxi( 7,3) = Nx(2)*Ny(2)*dNz(2)
      dshxi( 8,1) = dNx(1)*Ny(2)*Nz(2)
      dshxi( 8,2) = Nx(1)*dNy(2)*Nz(2)
      dshxi( 8,3) = Nx(1)*Ny(2)*dNz(2)
      dshxi( 9,1) = dNx(3)*Ny(1)*Nz(1)
      dshxi( 9,2) = Nx(3)*dNy(1)*Nz(1)
      dshxi( 9,3) = Nx(3)*Ny(1)*dNz(1)
      dshxi(10,1) = dNx(2)*Ny(3)*Nz(1)
      dshxi(10,2) = Nx(2)*dNy(3)*Nz(1)
      dshxi(10,3) = Nx(2)*Ny(3)*dNz(1)
      dshxi(11,1) = dNx(3)*Ny(2)*Nz(1)
      dshxi(11,2) = Nx(3)*dNy(2)*Nz(1)
      dshxi(11,3) = Nx(3)*Ny(2)*dNz(1)
      dshxi(12,1) = dNx(1)*Ny(3)*Nz(1)
      dshxi(12,2) = Nx(1)*dNy(3)*Nz(1)
      dshxi(12,3) = Nx(1)*Ny(3)*dNz(1)
      dshxi(13,1) = dNx(3)*Ny(1)*Nz(2)
      dshxi(13,2) = Nx(3)*dNy(1)*Nz(2)
      dshxi(13,3) = Nx(3)*Ny(1)*dNz(2)
      dshxi(14,1) = dNx(2)*Ny(3)*Nz(2)
      dshxi(14,2) = Nx(2)*dNy(3)*Nz(2)
      dshxi(14,3) = Nx(2)*Ny(3)*dNz(2)
      dshxi(15,1) = dNx(3)*Ny(2)*Nz(2)
      dshxi(15,2) = Nx(3)*dNy(2)*Nz(2)
      dshxi(15,3) = Nx(3)*Ny(2)*dNz(2)
      dshxi(16,1) = dNx(1)*Ny(3)*Nz(2)
      dshxi(16,2) = Nx(1)*dNy(3)*Nz(2)
      dshxi(16,3) = Nx(1)*Ny(3)*dNz(2)
      dshxi(17,1) = dNx(1)*Ny(1)*Nz(3)
      dshxi(17,2) = Nx(1)*dNy(1)*Nz(3)
      dshxi(17,3) = Nx(1)*Ny(1)*dNz(3)
      dshxi(18,1) = dNx(2)*Ny(1)*Nz(3)
      dshxi(18,2) = Nx(2)*dNy(1)*Nz(3)
      dshxi(18,3) = Nx(2)*Ny(1)*dNz(3)
      dshxi(19,1) = dNx(2)*Ny(2)*Nz(3)
      dshxi(19,2) = Nx(2)*dNy(2)*Nz(3)
      dshxi(19,3) = Nx(2)*Ny(2)*dNz(3)
      dshxi(20,1) = dNx(1)*Ny(2)*Nz(3)
      dshxi(20,2) = Nx(1)*dNy(2)*Nz(3)
      dshxi(20,3) = Nx(1)*Ny(2)*dNz(3)
      dshxi(21,1) = dNx(3)*Ny(3)*Nz(3)
      dshxi(21,2) = Nx(3)*dNy(3)*Nz(3)
      dshxi(21,3) = Nx(3)*Ny(3)*dNz(3)
      dshxi(22,1) = dNx(3)*Ny(3)*Nz(1)
      dshxi(22,2) = Nx(3)*dNy(3)*Nz(1)
      dshxi(22,3) = Nx(3)*Ny(3)*dNz(1)
      dshxi(23,1) = dNx(3)*Ny(3)*Nz(2)
      dshxi(23,2) = Nx(3)*dNy(3)*Nz(2)
      dshxi(23,3) = Nx(3)*Ny(3)*dNz(2)
      dshxi(24,1) = dNx(3)*Ny(1)*Nz(3)
      dshxi(24,2) = Nx(3)*dNy(1)*Nz(3)
      dshxi(24,3) = Nx(3)*Ny(1)*dNz(3)
      dshxi(25,1) = dNx(2)*Ny(3)*Nz(3)
      dshxi(25,2) = Nx(2)*dNy(3)*Nz(3)
      dshxi(25,3) = Nx(2)*Ny(3)*dNz(3)
      dshxi(26,1) = dNx(3)*Ny(2)*Nz(3)
      dshxi(26,2) = Nx(3)*dNy(2)*Nz(3)
      dshxi(26,3) = Nx(3)*Ny(2)*dNz(3)
      dshxi(27,1) = dNx(1)*Ny(3)*Nz(3)
      dshxi(27,2) = Nx(1)*dNy(3)*Nz(3)
      dshxi(27,3) = Nx(1)*Ny(3)*dNz(3)
      
      return
      end subroutine calcShape3DQuadFull

*************************************************************************

      subroutine mapShape1D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for the 2-node
      !  linear 1D element.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,1),dsh(nNode,1),coords(1,nNode),Le,detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian
      !
      Le = coords(1,2) - coords(1,1)
      detmapJ = half*Le


      ! Calculate first derivatives wrt x, y, z
      !
      do i=1,nNode
         dsh(i,1) = dshxi(i,1)/detmapJ
      enddo

      return
      end subroutine mapShape1D

*************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 4-node
      !  linear and 8-node quadratic quadrillateral 2D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode)
      real*8 mapJ(2,2),mapJ_inv(2,2),detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do

c$$$      write(*,*) 'dshxi='
c$$$      do i=1,nNode
c$$$         write(*,*) dshxi(i,:)
c$$$      enddo
c$$$
c$$$      write(*,*) 'coords='
c$$$      do i=1,nNode
c$$$         write(*,*) coords(:,i)
c$$$      enddo
c$$$
c$$$      write(*,*) 'mapJ='
c$$$      write(*,*) mapJ(1,1),mapj(1,2)
c$$$      write(*,*) mapJ(2,1),mapJ(2,2)


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape2D

*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 4-node
      !  linear and 8-node quadratic 2D elements.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a ``heat transfer'' and 
      !  ``static''step MCRD=2, but for ``coupled-temperature-displacement''
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode)
      real*8 mapJ(2,2),mapJ_inv(2,2),detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


c$$$      write(*,*) 'dshxi='
c$$$      do i=1,nNode
c$$$         write(*,*) dshxi(i,:)
c$$$      enddo
c$$$
c$$$      write(*,*) 'coords='
c$$$      do i=1,nNode
c$$$         write(*,*) coords(:,i)
c$$$      enddo
c$$$
c$$$      write(*,*) 'mapJ='
c$$$      write(*,*) mapJ(1,1),mapj(1,2)
c$$$      write(*,*) mapJ(2,1),mapJ(2,2)


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape2Da

*************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real*8 mapJ(3,3),mapJ_inv(3,3),detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)

c$$$      if(stat.eq.0) then
c$$$         write(*,*) 'dshxi='
c$$$         do i=1,nNode
c$$$            write(*,*) dshxi(i,:)
c$$$         enddo
c$$$         
c$$$         write(*,*) 'coords='
c$$$         do i=1,nNode
c$$$            write(*,*) coords(:,i)
c$$$         enddo
c$$$         
c$$$         write(*,*) 'mapJ='
c$$$         write(*,*) mapJ(1,1),mapj(1,2),mapj(1,3)
c$$$         write(*,*) mapJ(2,1),mapJ(2,2),mapj(2,3)
c$$$         write(*,*) mapJ(3,1),mapJ(3,2),mapj(3,3)
c$$$
c$$$         call xit
c$$$      endif


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape3D

************************************************************************

        subroutine matInv3D(A,A_INV,DET_A,istat)
        !
        ! Returns A_INV, the inverse and DET_A, the determinent
        ! Note that the det if of the original matrix, not the
        ! inverse <-- Prof Anand's subroutine has been modified.
        !
        IMPLICIT REAL*8 (A-H,O-Z)
      
        DIMENSION A(3,3), A_INV(3,3)
      
        PARAMETER(ZERO=0., ONE=1.)

        istat = 1
      
        DET_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     &          A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     &          A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
        
        IF (DET_A .LE. ZERO) THEN
          WRITE(*,*) 'WARNING: SUBROUTINE matInv3:'
          WRITE(*,*) 'WARNING: DET of MAT=',DET_A
          istat = 0
          return
        ENDIF
            
        DET_A_INV = ONE/DET_A
          
        A_INV(1,1) = DET_A_INV*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
        A_INV(1,2) = DET_A_INV*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
        A_INV(1,3) = DET_A_INV*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
        A_INV(2,1) = DET_A_INV*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
        A_INV(2,2) = DET_A_INV*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
        A_INV(2,3) = DET_A_INV*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
        A_INV(3,1) = DET_A_INV*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
        A_INV(3,2) = DET_A_INV*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
        A_INV(3,3) = DET_A_INV*(A(1,1)*A(2,2)-A(2,1)*A(1,2))

        RETURN
        END subroutine matInv3D

****************************************************************************

        subroutine matInv2D(A,A_INV,DET_A,istat)
        !
        ! Returns A_INV, the inverse and DET_A
        !
        IMPLICIT REAL*8 (A-H,O-Z)
      
        DIMENSION A(2,2), A_INV(2,2)
      
        PARAMETER(ZERO=0., ONE=1.)

        istat = 1
      
        DET_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
        IF (DET_A .LE. ZERO) THEN
          WRITE(*,*) 'WARNING: SUBROUTINE M2INV:'
          WRITE(*,*) 'WARNING: DET of MAT=',DET_A
          istat = 0
          return
        ENDIF
            
        DET_A_INV = ONE/DET_A
          
        A_INV(1,1) = DET_A_INV*A(2,2)
        A_INV(1,2) = -DET_A_INV*A(1,2)
        A_INV(2,1) = -DET_A_INV*A(2,1)
        A_INV(2,2) = DET_A_INV*A(1,1)

        return
        end subroutine matInv2D

****************************************************************************
