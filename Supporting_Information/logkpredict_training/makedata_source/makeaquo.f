c***********************************************************************
c***********************************************************************
      subroutine makeaquo(mtype,mcoord,filename)
c***********************************************************************
c     makeaquo: creates pcmodel file for metal aquo ion geometry  
c
c     passed variables:
c     integer    mtype     mm3 atom type for metal ion
c     integer    mcoord    aquo ion coordination number
c     char*30    filename  name of the pcmodel file
c
c     output:    <filename>.pcm file
c***********************************************************************
c***********************************************************************

      implicit none
      include 'sizes.i'

      integer i,j,n
      integer itype(MAXAT),n12(MAXAT),i12(MAXVAL,MAXAT)
      integer bo12(MAXVAL,MAXAT)
      integer mtype,mcoord
      
      real x(MAXAT),y(MAXAT),z(MAXAT)

      character*1 pichar(MAXAT)
      character*40 filename
      character*50 title

c***********************************************************************
c     initialize variables
c***********************************************************************

      do i = 1, MAXAT
         n12(i) = 0
         pichar(i) = ' '
         do j = 1, MAXVAL
            i12(j,i) = 0
            bo12(j,i) = 0
         end do
      end do

c***********************************************************************
c     if-then block of aquo ion types
c***********************************************************************

      if(mcoord.eq.4) then
         if((mtype.eq.338).or.(mtype.eq.369)) then       ! square planar
            n = 7
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            itype(6) = 201
            itype(7) = 201
            n12(1) = 6
            do i = 1, 6
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            do i = 2, 7
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) =  1.442e0
            y(2) = -1.442e0
            z(2) =  0.000e0
            x(3) =  1.442e0
            y(3) =  1.442e0
            z(3) =  0.000e0
            x(4) = -1.442e0
            y(4) =  1.442e0
            z(4) =  0.000e0
            x(5) = -1.442e0
            y(5) = -1.442e0
            z(5) =  0.000e0
            x(6) =  1.442e0
            y(6) = -1.442e0
            z(6) =  0.000e0
            x(7) =  1.442e0
            y(7) =  1.442e0
            z(7) =  0.000e0
            x(6) =  0.000e0
            y(6) =  0.000e0
            z(6) = -1.000e0
            x(7) =  0.000e0
            y(7) =  0.000e0
            z(7) =  1.000e0
         else                                              ! tetrahedral
            n = 5
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            n12(1) = 4
            do i = 1, 4
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            do i = 2, 5
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) = -1.120e0
            y(2) =  0.000e0
            z(2) = -1.584e0
            x(3) = -1.120e0
            y(3) =  0.000e0
            z(3) =  1.584e0
            x(4) =  1.120e0
            y(4) = -1.584e0
            z(4) =  0.000e0
            x(5) =  1.120e0
            y(5) =  1.584e0
            z(5) =  0.000e0
         end if
      elseif(mcoord.eq.5) then                ! bi-capped trigonal prism
         n = 6
         itype(1) = mtype
         itype(2) = 206
         itype(3) = 206
         itype(4) = 206
         itype(5) = 206
         itype(6) = 206
         n12(1) = 5
         do i = 1, 5
            i12(i,1) = i + 1
            bo12(i,1) = 1
         end do
         do i = 2, 6
            n12(i) = 1
            i12(1,i) = 1
            bo12(1,i) = 1
         end do
         x(1) =  0.000e0
         y(1) =  0.000e0
         z(1) =  0.000e0
         x(2) = -0.707e0
         y(2) =  0.707e0
         z(2) = -1.732e0
         x(3) = -0.707e0
         y(3) =  0.707e0
         z(3) =  1.732e0
         x(4) = -1.414e0
         y(4) = -1.414e0
         z(4) =  0.000e0
         x(5) =  1.414e0
         y(5) = -1.414e0
         z(5) =  0.000e0
         x(6) =  1.414e0
         y(6) =  1.414e0
         z(6) =  0.000e0
      elseif(mcoord.eq.6) then               
         if((mtype.eq.322).or.(mtype.eq.323)) then
            n = 7                             ! J-T distorted octahedron
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            itype(6) = 205
            itype(7) = 205
            n12(1) = 6
            do i = 1, 6
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            do i = 2, 7
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) =  1.442e0
            y(2) = -1.442e0
            z(2) =  0.000e0
            x(3) =  1.442e0
            y(3) =  1.442e0
            z(3) =  0.000e0
            x(4) = -1.442e0
            y(4) =  1.442e0
            z(4) =  0.000e0
            x(5) = -1.442e0
            y(5) = -1.442e0
            z(5) =  0.000e0
            x(6) =  0.000e0
            y(6) =  0.000e0
            z(6) = -2.300e0
            x(7) =  0.000e0
            y(7) =  0.000e0
            z(7) =  2.300e0
         elseif((mtype.eq.312).or.(mtype.eq.331).or.(mtype.eq.364)) then
            n = 7                                  ! mono-oxo octahedron
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            itype(6) = 206
            itype(7) = 207
            n12(1) = 6
            do i = 1, 6
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            bo12(6,1) = 2
            do i = 2, 7
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            bo12(1,7) = 2
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) =  1.515e0
            y(2) =  1.508e0
            z(2) =  0.000e0
            x(3) =  0.218e0
            y(3) =  0.000e0
            z(3) = -2.128e0
            x(4) = -1.484e0
            y(4) =  1.508e0
            z(4) = -0.307e0
            x(5) = -1.484e0
            y(5) = -1.508e0
            z(5) = -0.307e0
            x(6) =  1.516e0
            y(6) = -1.508e0
            z(6) =  0.000e0
            x(7) = -0.176e0
            y(7) =  0.000e0
            z(7) =  1.717e0
         elseif((mtype.eq.332).or.(mtype.eq.365)) then       ! cis-dioxo
            n = 7
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            itype(6) = 207
            itype(7) = 207
            n12(1) = 6
            do i = 1, 6
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            bo12(5,1) = 2
            bo12(6,1) = 2
            do i = 2, 7
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            bo12(1,6) = 2
            bo12(1,7) = 2
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) =  1.588e0
            y(2) =  1.482e0
            z(2) =  0.000e0
            x(3) =  1.588e0
            y(3) = -1.482e0
            z(3) =  0.000e0
            x(4) =  0.367e0
            y(4) =  0.169e0
            z(4) = -2.134e0
            x(5) = -1.352e0
            y(5) =  1.674e0
            z(5) = -0.257e0
            x(6) = -1.224e0
            y(6) = -1.173e0
            z(6) = -0.377e0
            x(7) = -0.017e0
            y(7) =  0.125e0
            z(7) =  1.732e0
         elseif(mtype.eq.367) then                         ! trans-dioxo
            n = 7
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            itype(6) = 207
            itype(7) = 207
            n12(1) = 6
            do i = 1, 6
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            bo12(5,1) = 2
            bo12(6,1) = 2
            do i = 2, 7
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            bo12(1,6) = 2
            bo12(1,7) = 2
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) =  1.527e0
            y(2) =  1.527e0
            z(2) =  0.000e0
            x(3) =  1.527e0
            y(3) = -1.527e0
            z(3) =  0.000e0
            x(4) = -1.527e0
            y(4) = -1.527e0
            z(4) =  0.000e0
            x(5) = -1.527e0
            y(5) =  1.527e0
            z(5) =  0.000e0
            x(6) =  0.000e0
            y(6) =  0.000e0
            z(6) = -1.730e0
            x(7) =  0.000e0
            y(7) =  0.000e0
            z(7) =  1.730e0
         else                                               ! octahedron
            n = 7
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            itype(6) = 206
            itype(7) = 206
            n12(1) = 6
            do i = 1, 6
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            do i = 2, 7
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) =  0.000e0
            y(2) =  0.000e0
            z(2) = -2.300e0
            x(3) = -1.442e0
            y(3) =  1.442e0
            z(3) =  0.000e0
            x(4) =  0.000e0
            y(4) =  0.000e0
            z(4) =  2.300e0
            x(5) = -1.442e0
            y(5) = -1.442e0
            z(5) =  0.000e0
            x(6) =  1.442e0
            y(6) = -1.442e0
            z(6) =  0.000e0
            x(7) =  1.442e0
            y(7) =  1.442e0
            z(7) =  0.000e0
         end if
      elseif(mcoord.eq.7) then 
         if((mtype.eq.382).or.(mtype.eq.383).or.(mtype.eq.384).or.
     &   (mtype.eq.386).or.(mtype.eq.387)) then            ! trans-dioxo
            n = 8
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            itype(6) = 206
            itype(7) = 207
            itype(8) = 207
            n12(1) = 7
            do i = 1, 7
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            bo12(6,1) = 2
            bo12(7,1) = 2
            do i = 2, 8
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            bo12(1,7) = 2
            bo12(1,8) = 2
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) =  2.008e0
            y(2) =  1.459e0
            z(2) =  0.000e0
            x(3) =  2.008e0
            y(3) = -1.459e0
            z(3) =  0.000e0
            x(4) = -0.767e0
            y(4) =  2.361e0
            z(4) =  0.000e0
            x(5) = -2.482e0
            y(5) =  0.000e0
            z(5) =  0.000e0
            x(6) = -0.767e0
            y(6) = -2.361e0
            z(6) =  0.000e0
            x(7) =  0.000e0
            y(7) =  0.000e0
            z(7) =  1.833e0
            x(8) =  0.000e0
            y(8) =  0.000e0
            z(8) = -1.833e0
         else                               ! mono-capped trigonal prism
            n = 8
            itype(1) = mtype
            itype(2) = 206
            itype(3) = 206
            itype(4) = 206
            itype(5) = 206
            itype(6) = 206
            itype(7) = 206
            itype(8) = 206
            n12(1) = 7
            do i = 1, 7
               i12(i,1) = i + 1
               bo12(i,1) = 1
            end do
            do i = 2, 8
               n12(i) = 1
               i12(1,i) = 1
               bo12(1,i) = 1
            end do
            x(1) =  0.000e0
            y(1) =  0.000e0
            z(1) =  0.000e0
            x(2) =  0.471e0
            y(2) = -0.092e0
            z(2) = -2.100e0
            x(3) = -1.055e0
            y(3) =  1.770e0
            z(3) = -0.629e0
            x(4) = -0.284e0
            y(4) =  1.012e0
            z(4) =  1.879e0
            x(5) = -1.808e0
            y(5) = -0.849e0
            z(5) = -0.802e0
            x(6) = -0.554e0
            y(6) = -1.561e0
            z(6) =  1.399e0
            x(7) =  1.666e0
            y(7) = -1.364e0
            z(7) =  0.000e0
            x(8) =  1.666e0
            y(8) =  1.364e0
            z(8) =  0.000e0
         end if
      elseif(mcoord.eq.8) then                        ! square antiprism
         n = 9
         itype(1) = mtype
         itype(2) = 206
         itype(3) = 206
         itype(4) = 206
         itype(5) = 206
         itype(6) = 206
         itype(7) = 206
         itype(8) = 206
         itype(9) = 206
         n12(1) = 8
         do i = 1, 8
            i12(i,1) = i + 1
            bo12(i,1) = 1
         end do
         do i = 2, 9
            n12(i) = 1
            i12(1,i) = 1
            bo12(1,i) = 1
         end do
         x(1) =  0.000e0
         y(1) =  0.000e0
         z(1) =  0.000e0
         x(2) =  0.666e0
         y(2) =  0.000e0
         z(2) =  2.321e0
         x(3) = -0.858e0
         y(3) =  2.042e0
         z(3) =  0.962e0
         x(4) = -0.219e0
         y(4) =  1.444e0
         z(4) = -1.923e0
         x(5) = -0.219e0
         y(5) = -1.444e0
         z(5) = -1.923e0
         x(6) = -2.382e0
         y(6) =  0.000e0
         z(6) = -0.398e0
         x(7) = -0.858e0
         y(7) = -2.042e0
         z(7) =  0.961e0
         x(8) =  1.935e0
         y(8) =  1.444e0
         z(8) =  0.000e0
         x(9) =  1.935e0
         y(9) = -1.444e0
         z(9) =  0.000e0
      elseif(mcoord.eq.9) then            ! tricapped trigonal antiprism
         n = 10
         itype(1) = mtype
         itype(2) = 206
         itype(3) = 206
         itype(4) = 206
         itype(5) = 206
         itype(6) = 206
         itype(7) = 206
         itype(8) = 206
         itype(9) = 206
         itype(10) = 206
         n12(1) = 9
         do i = 1, 9
            i12(i,1) = i + 1
            bo12(i,1) = 1
         end do
         do i = 2, 10
            n12(i) = 1
            i12(1,i) = 1
            bo12(1,i) = 1
         end do
         x(1) =  0.000e0
         y(1) =  0.000e0
         z(1) =  0.000e0
         x(2) =  0.497e0
         y(2) =  0.843e0
         z(2) = -2.358e0
         x(3) = -0.643e0
         y(3) =  2.470e0
         z(3) =  0.000e0
         x(4) =  0.962e0
         y(4) =  0.178e0
         z(4) =  2.358e0
         x(5) = -1.871e0
         y(5) =  0.437e0
         z(5) =  1.698e0
         x(6) = -2.260e0
         y(6) =  0.125e0
         z(6) = -1.179e0
         x(7) = -0.656e0
         y(7) = -2.167e0
         z(7) =  1.179e0
         x(8) = -0.230e0
         y(8) = -1.908e0
         z(8) = -1.698e0
         x(9) =  2.101e0
         y(9) = -1.460e0
         z(9) =  0.000e0
         x(10) =  2.101e0
         y(10) =  1.460e0
         z(10) =  0.000e0
      end if

c***********************************************************************
c     write pcmodel file
c***********************************************************************

      title = ' '//filename

      call writepcm(filename,title,n,x,y,z,itype,n12,i12,bo12,pichar)

      return
      end
