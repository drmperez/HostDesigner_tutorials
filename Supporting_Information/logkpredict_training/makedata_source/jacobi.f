c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Adapted from TINKER code jacobi.f (https://dasher.wustl.edu) with
c     permission from Jay Ponder
c
c     'jacobi' - performs matrix diagonalization of real symmetric
c     matrix using Jacobi rotations
c
c     passed variables:
c     integer  n    logical dimension of the matrix to be diagonalized
c     integer  np   physical dimension of the matrix storage area
c     real     a    input with matrix to be diagonalized; only upper
c                   upper triangle and diagonal are required
c     real     d    returned with the eigenvalues in ascending order
c     real     v    returned with the eigenvectors of the matrix
c     real     b    temporary work vector
c     real     z    temporary work vector
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine jacobi (n,np,a,d,v,b,z)

      implicit none

      integer i,j,k,ip,iq,n,np,nrot,maxrot

      real sm,tresh,s,c,t,theta,tau,h,g,p
      real a(np,np),d(np),v(np,np),b(np),z(np)

c-----------------------------------------------------------------------
c     Initialize
c-----------------------------------------------------------------------
 
      maxrot = 1000
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0e0
         end do
         v(ip,ip) = 1.0e0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0e0
      end do
 
c-----------------------------------------------------------------------
c     Perform jacobi rotations
c-----------------------------------------------------------------------
 
      do i = 1, maxrot
         sm = 0.0e0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. (0.0e0))  goto 10
         if (i .lt. 4) then
            tresh = 0.2e0*sm / n*n
         else
            tresh = 0.0e0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0e0 * abs(a(ip,iq))
               if ((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     &                    .and.(abs(d(iq))+g.eq.abs(d(iq)))) then
                  a(ip,iq) = 0.0e0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if ((abs(h)+g).eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5e0*h / a(ip,iq)
                     t = 1.0e0 / (abs(theta)+sqrt(1.0e0+theta*theta))
                     if (theta .lt. 0.0e0)  t = -t
                  end if
                  c = 1.0e0 / sqrt(1.0e0+t*t)
                  s = t * c
                  tau = s / (1.0e0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0e0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0e0
         end do
      end do

   10 continue
      if (nrot.eq.maxrot) then
         write(6,'(a)') ' '
         write(6,'(a)') 'Warning: subroutine jacobi:'
         write(6,'(a)') 'matrix diagonalization not converged'
         write(6,'(a)') ' '
      end if

      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do
      return
      end
