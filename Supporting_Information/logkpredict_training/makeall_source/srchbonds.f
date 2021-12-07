c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'srchbonds' -- returns searchable bonds and flag for whether
c     ring system is searchable
c
c     passed variables:
c     integer  n          number of atoms in structure
c     integer  n12,i12    connectivity arrays
c     integer  bo12       bond order arrays
c     char*2   atlab      atom label
c     integer  nsb        number searchable bonds
c     integer  isb        ith atom of searchable bond
c     integer  jsb        jth atom of searchable bond
c     logical  doall      true only if called from nrotbond
c     logical  srchring   .true. if ring should be searched
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine srchbonds(n,n12,i12,bo12,atlab,nsb,isb,jsb,doall,
     & srchring)

      implicit none
      include 'sizes.i'

      integer n,n12(*),i12(MAXVAL,*),bo12(MAXVAL,*)
      integer i,j,k,cycle,nt,ith,jth
      integer ii,jj,start
      integer nx,n12x(MAXAT),i12x(MAXVAL,MAXAT)
      integer tn12,ti12(MAXVAL)
      integer nsb,isb(MAXSB),jsb(MAXSB)
      integer nrb,irb(MAXAT),jrb(MAXAT)
      integer ip(MAXAT),ones(MAXAT),nsingle,nbridge
      integer nmol,n12r(MAXAT),i12r(MAXVAL,MAXAT)
      integer nparts,parts(MAXAT),nsp3

      character*2 atlab(*)

      logical srchring,doall

c-----------------------------------------------------------------------
c     Initialize
c-----------------------------------------------------------------------

      srchring = .false.
      nsb = 0
      nmol = 0
      do i = 1, MAXSB
         isb(i) = 0
         jsb(i) = 0
      end do
      do i = 1, MAXAT
         ip(i) = 0
         n12x(i) = 0
         do j = 1, MAXVAL
            i12x(j,i) = 0
         end do
      end do

      do i = 1, n
         ip(i) = i
         n12x(i) = n12(i)
         do j = 1, n12x(i)
            i12x(j,i) = i12(j,i)
         end do
      end do

c-----------------------------------------------------------------------
c     Start removing terminal atoms
c-----------------------------------------------------------------------

      nx = n
      cycle = 0

   50 continue

c-----------------------------------------------------------------------
c     Find all terminal atoms
c-----------------------------------------------------------------------

      nt = 0
      do i = 1, nx
         if(n12x(i).eq.1) then
            nt = nt + 1
            ones(nt) = i
         end if
      end do

      if(nt.eq.0) goto 150

      nt = nt + 1
      cycle = cycle + 1

c-----------------------------------------------------------------------
c     Loop to remove all terminal atoms
c-----------------------------------------------------------------------

  100 continue

      nt = nt - 1
      k = ones(nt)

      ith = k                                            ! terminal atom
      jth = i12x(1,ith)                ! atom connected to terminal atom

      ith = ip(ith)
      jth = ip(jth)

      call ratconn(k,nx,ip,n12x,i12x)

c-----------------------------------------------------------------------
c     Test if terminal atom is connected by a single bond
c-----------------------------------------------------------------------

      if(nx.ge.1) then
         if(doall) then
            if((atlab(ith).ne.'H ').and.(atlab(jth).ne.'H ')) then
               call addbond(ith,jth,n12,i12,bo12,nsb,isb,jsb)
            end if
         else
            if(cycle.eq.2) then
               if(atlab(ith).ne.'C ') then         ! not CH3, CH2, or CH
                  call addbond(ith,jth,n12,i12,bo12,nsb,isb,jsb)
               end if
            else if(cycle.gt.2) then
               call addbond(ith,jth,n12,i12,bo12,nsb,isb,jsb)
            end if
         end if
      end if

      if(nt.gt.1) goto 100

c-----------------------------------------------------------------------
c     Go again if any terminal atoms remain
c-----------------------------------------------------------------------

      if(nx.gt.1) then
         do i = 1, nx
            if(n12x(i).eq.1) goto 50        ! if any terminal atoms left
         end do
      end if

  150 continue

c-----------------------------------------------------------------------
c     If nx > 2, then there is a ring system, check to see if there are
c     any rotatable bonds between individual rings
c-----------------------------------------------------------------------

      if(nx.gt.2) then                            ! there are ring atoms

c-----------------------------------------------------------------------
c     Check ring bonds to see if in a ring or connecting two rings
c-----------------------------------------------------------------------

         call nummol(nx,n12x,i12x,nmol,parts)
         call nringbonds(nx,n12x,i12x,nrb,irb,jrb)
         start = 1

  200    continue
         do ii = start, nrb
            ith = irb(ii)
            jth = jrb(ii)

c-----------------------------------------------------------------------
c     Remove the ith-jth bond, temp store result in n12r,i12r
c-----------------------------------------------------------------------

            do i = 1, nx
               if(i.eq.ith) then
                  tn12 = 0
                  do j = 1, n12x(i)
                     if(i12x(j,i).ne.jth) then
                        tn12 = tn12 + 1
                        ti12(tn12) = i12x(j,i)
                     end if
                  end do
                  n12r(i) = tn12
                  do j = 1, n12r(i)
                     i12r(j,i) = ti12(j)
                  end do
               else if(i.eq.jth) then
                  tn12 = 0
                  do j = 1, n12x(i)
                     if(i12x(j,i).ne.ith) then
                        tn12 = tn12 + 1
                        ti12(tn12) = i12x(j,i)
                     end if
                  end do
                  n12r(i) = tn12
                  do j = 1, n12r(i)
                     i12r(j,i) = ti12(j)
                  end do
               else
                  n12r(i) = n12x(i)
                  do j = 1, n12r(i)
                     i12r(j,i) = i12x(j,i)
                  end do
               end if
            end do

c-----------------------------------------------------------------------
c     Number of molecules?
c-----------------------------------------------------------------------

            call nummol(nx,n12r,i12r,nparts,parts)

c-----------------------------------------------------------------------
c     If number molecules one bigger after removing ith-jth bond
c-----------------------------------------------------------------------

            if(nparts.eq.nmol+1) then            

c-----------------------------------------------------------------------
c     Add lost bond to rotatable bond list?
c-----------------------------------------------------------------------

               ith = ip(ith)
               jth = ip(jth)
               call addbond(ith,jth,n12,i12,bo12,nsb,isb,jsb)

c-----------------------------------------------------------------------
c     Update the n12x,i12x arrays
c-----------------------------------------------------------------------

               do i = 1, nx
                  n12x(i) = n12r(i)
                  do j = 1, n12x(i)
                     i12x(j,i) = i12r(j,i)
                  end do
               end do

               nmol = nmol + 1

c-----------------------------------------------------------------------
c     Remove any terminal atoms from the parts (see prior loop above)
c-----------------------------------------------------------------------

  250          continue

               nt = 0
               do i = 1, nx
                  if(n12x(i).eq.1) then
                     nt = nt + 1
                     k = i
                     goto 300
                  end if
               end do
  300          continue

               if(nt.gt.0) then

                  ith = ip(k)
                  jth = ip(i12x(1,k))

                  call ratconn(k,nx,ip,n12x,i12x)  
                  call addbond(ith,jth,n12,i12,bo12,nsb,isb,jsb)

                  if(nt.gt.1) goto 300

                  do i = 1, nx
                     if(n12x(i).eq.1) goto 250 
                  end do

               end if

               goto 400

            end if

         end do

         goto 500                                          ! normal exit

  400    continue
         call nringbonds(nx,n12x,i12x,nrb,irb,jrb)
         start = ii
         goto 200

  500    continue

      end if

c-----------------------------------------------------------------------
c     Sort rotatable bond list by ith bond
c-----------------------------------------------------------------------

      start = 1
      do j = start, nsb
         ii = isb(j)
         jj = jsb(j)
         do i = j-1, start,-1
            if(isb(i).le.ii) goto 600
            isb(i+1) = isb(i)
            jsb(i+1) = jsb(i)
         end do
         i = start - 1
  600    isb(i+1) = ii
         jsb(i+1) = jj
      end do

c-----------------------------------------------------------------------
c     Get final list of ring bonds
c-----------------------------------------------------------------------

      if(nrb.gt.0) then

         call nringbonds(nx,n12x,i12x,nrb,irb,jrb)

c-----------------------------------------------------------------------
c     Number of ring single bonds
c-----------------------------------------------------------------------

         nsingle = 0

         do i = 1, nrb
            ii = ip(irb(i))
            jj = ip(jrb(i))
            if((n12(ii).eq.4).or.(n12(jj).eq.4)) then
                nsingle = nsingle + 1
            end if
         end do

         if(nsingle.ge.4) srchring = .true.

c-----------------------------------------------------------------------
c     Number of sp2 C atoms in the ring system
c-----------------------------------------------------------------------

         nsp3 = 0
         do i = 1, nx
            ii = ip(i)
            if((atlab(ii).eq.'C ').and.(n12(ii).eq.4)) then
               nsp3 = nsp3 + 1
            end if
         end do

         if(nsp3.lt.4) srchring = .false.

c-----------------------------------------------------------------------
c     Number of bridging bonds if single ring molecule
c-----------------------------------------------------------------------

         if(nmol.eq.1) then
            nbridge = 0
            do i = 1, nx
               if(n12x(i).ge.3) nbridge = nbridge + 1
            end do

            if( (nx.lt.8).and.(nbridge.ge.2).or.
     &          (nx.lt.11).and.(nbridge.ge.4).or.
     &          (nx.lt.14).and.(nbridge.ge.6).or.
     &          (nx.lt.17).and.(nbridge.ge.8) ) then
               srchring = .false.
            end if
         end if

      end if

      return
      end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subroutine addbond - adds rotatable bond to list
c
c     Passed:
c     integer ith,jth   bond atoms
c     integer n12,i12   original connectivity lists
c     integer bo12      bond orders
c     integer nsb       number single bonds
c     integer isb       ith atom of bond
c     integer jsb       jth atom of bond
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine addbond(ith,jth,n12,i12,bo12,nsb,isb,jsb)

      implicit none
      include 'sizes.i'
      integer n12(*),i12(MAXVAL,*),bo12(MAXVAL,*)
      integer i,ith,jth,iswap
      integer nsb,isb(*),jsb(*)

      if(jth.lt.ith) then
         iswap = ith
         ith = jth
         jth = iswap
      end if

      do i = 1, nsb
         if((isb(i).eq.ith).and.(jsb(i).eq.jth)) return
      end do

      do i = 1, n12(ith)
         if(i12(i,ith).eq.jth) then
            if(bo12(i,ith).eq.1) then
               nsb = nsb + 1
               if(nsb.gt.MAXSB) then
            write(6,'(a)') 'Maximum number of searchable bonds has been'
            write(6,'(a)') 'exceeded.  If you really want to search '
            write(6,'(a,i3,a)') 'hits with > ',MAXSB,' rotatable bonds'
            write(6,'(a)') 'increase MAXSP in params.i and remake code'
                  stop
               end if
               isb(nsb) = ith
               jsb(nsb) = jth
            end if
         end if
      end do
      
      return
      end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subroutine nringbonds - returns number of ring bonds
c
c     Passed:
c     integer n         number of atoms
c     integer n12,i12   connectivity arrays
c     integer nrb       number of bonds
c     integer irb       ith atom of bond
c     integer jrb       jth atom of bond
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine nringbonds(n,n12,i12,nrb,irb,jrb)
      implicit none
      include 'sizes.i'
      integer i,j,k,n,n12(*),i12(MAXVAL,*)
      integer nrb,irb(*),jrb(*),ith,jth

      nrb = 0
      do i = 1, MAXAT
         irb(i) = 0
         jrb(i) = 0
      end do

      do i = 1, n
         do j = 1, n12(i)
            ith = i
            jth = i12(j,i)
            if(ith.gt.jth) then
               jth = i
               ith = i12(j,i)
            end if
            do k = 1, nrb
               if((irb(k).eq.ith).and.(jrb(k).eq.jth)) goto 10
            end do
            nrb = nrb + 1
            irb(nrb) = ith
            jrb(nrb) = jth
   10    continue
         end do
      end do

      return
      end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subroutine ratconn - rm atom k from connectivity list
c
c     Passed:
c     integer k         atom to be removed
c     integer n         number of atoms
c     integer ip        pointer list
c     integer n12,i12   connectivity arrays
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine ratconn(k,n,ip,n12,i12)
      implicit none
      include 'sizes.i'
      integer k,n,ip(*),n12(*),i12(MAXVAL,*)
      integer i,j,tn12,ti12(MAXVAL)

      n = n - 1
      do i = k, n
         ip(i) = ip(i+1)
         n12(i) = n12(i+1)
         do j = 1, n12(i+1)
            i12(j,i) = i12(j,i+1)
         end do
      end do
      do i = 1, n
         tn12 = 0
         do j = 1, n12(i)
            if(i12(j,i).lt.k) then
               tn12 = tn12 + 1
               ti12(tn12) = i12(j,i)
            else if(i12(j,i).gt.k) then
               tn12 = tn12 + 1
               ti12(tn12) = i12(j,i) - 1
            end if
         end do
         n12(i) = tn12
         do j = 1, n12(i)
            i12(j,i) = ti12(j)
         end do
      end do
      return
      end
