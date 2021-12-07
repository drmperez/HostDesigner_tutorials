c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Adapted from TINKER code center.f (https://dasher.wustl.edu) with 
c     permission from Jay Ponder
c
c     'center' -- moves weighted centroid of each molecule to origin
c
c     passed variables:
c     integer    n1         number of atoms in molecule 1
c     real       x1,y1,z1   coords of molecule 1
c     integer    n2         number of atoms in molecule 2
c     real       x2,y2,z2   coords of molecule 2
c     real       xmid-zmid  weighted centroid of first molecule
c     integer    nfit       number of atom pairs to overlap
c     integer    ifit(1,i)  list of atoms to overlap in molecule 1
c     integer    ifit(2,i)  list of atoms to overlap in molecule 2
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine center (n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid,
     &                   nfit,ifit)
      implicit none
      integer i,k,n1,n2
      real xmid,ymid,zmid,norm
      real x1(*),y1(*),z1(*)
      real x2(*),y2(*),z2(*) 

      integer nfit,ifit(2,*)

c-----------------------------------------------------------------------
c     Find weighted centroid of second structure, translate to origin
c-----------------------------------------------------------------------

      xmid = 0.0e0
      ymid = 0.0e0
      zmid = 0.0e0
      norm = real(nfit)
      do i = 1, nfit
         k = ifit(2,i)
         xmid = xmid + x2(k) 
         ymid = ymid + y2(k) 
         zmid = zmid + z2(k) 
      end do

      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n2
         x2(i) = x2(i) - xmid
         y2(i) = y2(i) - ymid
         z2(i) = z2(i) - zmid
      end do

c-----------------------------------------------------------------------
c     Find weighted centroid of first structure, translate to origin
c-----------------------------------------------------------------------

      xmid = 0.0e0
      ymid = 0.0e0
      zmid = 0.0e0
      norm = real(nfit)
      do i = 1, nfit
         k = ifit(1,i)
         xmid = xmid + x1(k) 
         ymid = ymid + y1(k) 
         zmid = zmid + z1(k) 
      end do

      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n1
         x1(i) = x1(i) - xmid
         y1(i) = y1(i) - ymid
         z1(i) = z1(i) - zmid
      end do

      return
      end
