c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Adapted from TINKER code molecule.f (https://dasher.wustl.edu) 
c     with permission from Jay Ponder
c 
c     'nummol' -- returns number of molecules in structure
c
c     passed variables:
c     integer  n          number of atoms in structure
c     integer  n12,i12    connectivity arrays
c     integer  npart      number of molecules in structure
c     integer  part       which molecule atom belongs to
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine nummol(n,n12,i12,npart,part)

      implicit none
      include 'params.i'

      integer n,n12(*),i12(MAXVAL,*),npart
      integer i,j,k,iatt
      integer ip,jp,kp,part(MAXATOM)

      npart = 0
      do i = 1, MAXATOM
         part(i) = 0
      end do

      do i = 1, n
         if(part(i).eq.0) then
            npart = npart + 1
            part(i) = npart
         end if
         ip = part(i)
         do iatt = 1, n12(i)
            j = i12(iatt,i)
            jp = part(j)
            if(jp.eq.0) then
               part(j) = ip
            else if(ip.lt.jp) then
               npart = npart - 1
               do k = 1, n
                  kp = part(k)
                  if(kp.eq.jp) then
                     part(k) = ip
                  else if(kp.gt.jp) then
                     part(k) = kp - 1
                  end if
               end do
            else if(ip.gt.jp) then
               npart = npart - 1
               do k = 1, n
                  kp = part(k)
                  if(kp.eq.ip) then
                     part(k) = jp
                  else if(kp.gt.ip) then
                     part(k) = kp - 1
                  end if
               end do
               ip = jp
            end if
         end do
      end do
 
      return
      end
