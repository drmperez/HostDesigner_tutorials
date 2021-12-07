c***********************************************************************
c***********************************************************************
c     'delatom' -- removes atom at serial number spot from molecule 
c   
c     passed variables:
c     integer spot      serial number of atom to be removed
c     integer n         number of atoms
c     integer itype     atom type
c     integer n12,i12   connectivity arrays
c     integer bo12      bond orders
c     real x,y,z        coordinates
c     char atlab        atom labels
c     char pichar       pi string
c***********************************************************************
c***********************************************************************

      subroutine delatom(spot,n,itype,n12,i12,bo12,x,y,z,atlab,pichar)
      
      implicit none
      include 'sizes.i'

      integer i,j,spot
      integer n,itype(*),n12(*),i12(MAXVAL,*),bo12(MAXVAL,*)
      integer tn12,ti12(MAXVAL),tbo12(MAXVAL)

      real x(*),y(*),z(*)

      character*1 pichar(*)
      character*2 atlab(*)

c***********************************************************************
c     delete atom at spot
c***********************************************************************

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
         atlab(i) = atlab(i+1)
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

c***********************************************************************
c     fix atom types
c***********************************************************************
c     do i = 1, n
c        if(itype(i).eq.147) itype(i) = 47
c        if(itype(i).eq.205) itype(i) = 6
c        if(itype(i).eq.208) itype(i) = 8
c        if(itype(i).eq.237) itype(i) = 37
c        if(itype(i).eq.238) itype(i) = 37
c        if(itype(i).eq.239) itype(i) = 37
c        if(itype(i).eq.241) itype(i) = 41
c        if(itype(i).eq.247) itype(i) = 47
c        if(itype(i).eq.267) itype(i) = 167
c        if(itype(i).eq.268) itype(i) = 168
c        if(itype(i).eq.269) itype(i) = 69
c        if(itype(i).eq.270) itype(i) = 165
c        if(itype(i).eq.271) itype(i) = 166
c        if(itype(i).eq.272) itype(i) = 168
c        if(itype(i).eq.278) itype(i) = 176
c     end do
         
      return
      end
