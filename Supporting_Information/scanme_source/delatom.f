c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'delatom' -- removes atom at serial number spot from molecule 
c   
c     passed variables:
c     integer spot      serial number of atom to be removed
c     integer n         number of atoms
c     integer itype     atom type
c     integer n12,i12   connectivity arrays
c     integer bo12      bond orders
c     real x,y,z        coordinates
c     char pichar       pi string
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine delatom(spot,n,itype,n12,i12,bo12,x,y,z,pichar)
      
      implicit none
      include 'params.i'

      integer i,j,spot
      integer n,itype(*),n12(*),i12(MAXVAL,*),bo12(MAXVAL,*)
      integer tn12,ti12(MAXVAL),tbo12(MAXVAL)

      real x(*),y(*),z(*)

      character*1 pichar(*)

      n = n - 1

      do i = spot, n
         itype(i) = itype(i+1)
         n12(i) = n12(i+1)
         do j = 1, n12(i)
            i12(j,i) = i12(j,i+1)
            bo12(j,i) = bo12(j,i+1)
         end do
         x(i) = x(i+1)
         y(i) = y(i+1)
         z(i) = z(i+1)
         pichar(i) = pichar(i+1)
      end do

      do i = 1, n
         tn12 = 0
         do j = 1, n12(i)
            if(i12(j,i).lt.spot) then
               tn12 = tn12 + 1
               ti12(tn12) = i12(j,i)
               tbo12(tn12) = bo12(j,i)
            else if (i12(j,i).gt.spot) then
               tn12 = tn12 + 1
               ti12(tn12) = i12(j,i) - 1
               tbo12(tn12) = bo12(j,i)
            end if
         end do
         n12(i) = tn12
         do j = 1, n12(i)
            i12(j,i) = ti12(j)
            bo12(j,i) = tbo12(j)
         end do
      end do
         
      return
      end
