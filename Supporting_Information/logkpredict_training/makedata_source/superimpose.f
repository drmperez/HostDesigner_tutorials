c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Adapted from TINKER code impose.f (https://dasher.wustl.edu) with 
c     permission from Jay Ponder
c
c     'superimpose' -- overlays two sets of atoms, returns root mean
c     squared displacement of best superposition       
c
c     passed variables:
c     integer    n1         number of atoms in molecule 1
c     real       x1,y1,z1   coords of molecule 1
c     integer    n2         number of atoms in molecule 2
c     real       x2,y2,z2   coords of molecule 2
c     integer    nfit       number of atom pairs to overlap
c     integer    ifit(1,i)  list of atoms to overlap in molecule 1
c     integer    ifit(2,i)  list of atoms to overlap in molecule 2
c     real       rmsd       average displacement of superimposed atoms
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine superimpose(n1,x1,y1,z1,n2,x2,y2,z2,nfit,ifit,rmsd)

      implicit none
      
      integer i,n1,n2
      integer nfit,ifit(2,*)
      
      real x1(*),y1(*),z1(*)
      real x2(*),y2(*),z2(*)
      real xmid,ymid,zmid
      real dx,dy,dz,dist,rmsd

c-----------------------------------------------------------------------
c     Position two groups of atoms at their centers of coords or mass
c-----------------------------------------------------------------------

      call center(n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid,nfit,ifit)

c-----------------------------------------------------------------------
c     Rotate to get best overlap
c-----------------------------------------------------------------------

      call quatfit(x1,y1,z1,n2,x2,y2,z2,nfit,ifit)

c-----------------------------------------------------------------------
c     Translate both molecules back to where the first one was
c-----------------------------------------------------------------------

      do i = 1, n1
         x1(i) = x1(i) + xmid
         y1(i) = y1(i) + ymid
         z1(i) = z1(i) + zmid
      end do

      do i = 1, n2
         x2(i) = x2(i) + xmid
         y2(i) = y2(i) + ymid
         z2(i) = z2(i) + zmid
      end do
      
c-----------------------------------------------------------------------
c     Calculate rmsd 
c-----------------------------------------------------------------------
      
      rmsd = 0.0e0
      
      do i = 1, nfit
         dx = x1(ifit(1,i)) - x2(ifit(2,i))
         dy = y1(ifit(1,i)) - y2(ifit(2,i))
         dz = z1(ifit(1,i)) - z2(ifit(2,i))
         dist = dx*dx+dy*dy+dz*dz
         rmsd = rmsd + dist
      end do
      
      rmsd = sqrt(rmsd/nfit)
            
      return      
      end
