c***********************************************************************
c***********************************************************************
c     'confanal' -- this routine sets up and performs conformational
c     analysis on passed in input structure
c***********************************************************************
c***********************************************************************
c     Passed variables:
c     char*30    infile     pcmodel input file
c     logical    dosrch     search if true or only optimize if false
c     integer    n          number of atoms
c     real       x,y,z      coordinates
c     integer    n12,i12    connectivity arrays
c     integer    bo12       bond orders
c     char*2     atlab      atom labels
c     char*1     pichar     pi atom array
c     real       elow       energy of global minimum
c***********************************************************************
c***********************************************************************

      subroutine confanal(infile,dosrch,n,x,y,z,type,n12,i12,bo12,atlab,
     &  pichar,elow)
      
      implicit none
      include 'sizes.i'

      integer i,j
      integer n,n12(MAXAT),i12(MAXVAL,MAXAT),type(MAXAT)
      integer bo12(MAXVAL,MAXAT)
      integer nsb,isb(MAXSB),jsb(MAXSB),sincr(MAXSB),nheavy
      integer maxsbond,starttime,endtime,time,icstate
      
      real x(MAXAT),y(MAXAT),z(MAXAT)
      real e1,e2,elow

      character*1 pichar(MAXAT)
      character*2 atlab(MAXAT)
      character*40 infile,outfile                  ! name of pcmod files
      character*120 cmd

      logical flag,srchring,doall,dosrch

c***********************************************************************
c     Get molecular info
c***********************************************************************

      doall = .false.
      call srchbonds(n,n12,i12,bo12,atlab,nsb,isb,jsb,doall,srchring)

      do i = 1, nsb
         if(n12(isb(i)).eq.4.and.n12(jsb(i)).eq.4) then
            sincr(i) = 120.0e0
         else
            sincr(i) = 60.0e0
         end if
      end do

      nheavy = 0
      do i = 1, n
         if(atlab(i).ne.'H ') then
            nheavy = nheavy + 1
         end if
      end do

      maxsbond = 8
      if(nsb.gt.maxsbond) then
         icstate = 1                             ! too many search bonds
         goto 200                                         ! write to log
      end if

      if(.not.dosrch) then
         nsb = 0
         srchring = .false.
      else
         srchring = .true.
      end if

      if(infile(1:19).eq.'salicylaldehyde.pcm') then       ! special case
         nsb = 1
         isb(1) = 1
         jsb(1) = 7
         sincr(1) = 180.0e0
         srchring = .false.
      end if

      if( (nsb.eq.0).and.(.not.srchring) ) then
         icstate = 2                                 ! nothing to search
         outfile = 'outfile.pcm'
         open(unit=33,file='conpcm')
            write(33,'(a)') 'mode opt'
            write(33,'(a,a)') 'infile ',infile
            write(33,'(a,a)') 'outfile ',outfile
            write(33,'(a)') 'forcefield mm3'
            write(33,'(a)') 'dielec 4.0'
         close(33)
         call system('mengine')
         goto 200
      end if

c-----------------------------------------------------------------------
c     Write engine control file
c-----------------------------------------------------------------------

      outfile = 'outfile.pcm'
      call writesrchcon(infile,outfile,nsb,isb,jsb,sincr,srchring)

c-----------------------------------------------------------------------
c     Run mengine with preset time limit in seconds
c-----------------------------------------------------------------------

      starttime = time()
      open(unit=33,file='fixme',status='unknown')
      write(33,'(a)') 'limit cputime 1000'
      write(33,'(a)') 'mengine conpcm'
      close(33)
      cmd = 'csh < fixme'
      call system(cmd)
      endtime = time()

c-----------------------------------------------------------------------
c     Check whether output file is produced
c-----------------------------------------------------------------------

      inquire(file=outfile,exist=flag)
      if(.not.flag) then
         icstate = 3                                     ! exceeded time
      else
         icstate = 4                                        ! successful
      end if

c-----------------------------------------------------------------------
c     Write results
c-----------------------------------------------------------------------

  200 continue

      if(icstate.eq.1) then                      ! too many search bonds

         write(6,'(a)') ' '
         write(6,'(a)') 'SEARCH NOT ATTEMPTED'
         write(6,'(a)') 'exceeded maximum of searchable bonds'
         write(6,'(a,i3)') 'control variable maxsbond = ',maxsbond
         write(6,'(a)') ' '

      else if(icstate.eq.2) then                     ! nothing to search
         write(6,'(a)') ' '
         write(6,'(a)') 'SEARCH NOT ATTEMPTED'
         write(6,'(a)') 'No rotatable bonds and no searchable rings'
         write(6,'(a)') ' '
         call readpcm(outfile,n,type,x,y,z,e1,e2,n12,i12,bo12,i,j,
     &    pichar)
         elow = e2
         call system('rm outfile.pcm')

      else if(icstate.eq.3) then                   ! exceeded time limit

         write(6,'(a)') ' '
         write(6,'(a)') 'SEARCH ATTEMPTED BUT FAILED'
         write(6,'(a)') 'No output produced within time limit'
         write(6,'(a,i7)') 'runtime(sec) = ',endtime-starttime
         write(6,'(a)') ' '

      else if(icstate.eq.4) then                      ! search completed

         call getlow(outfile,x,y,z,elow) 
         write(6,'(a)') ' '
         write(6,'(a)') 'SEARCH COMPLETED'
         write(6,'(a,i7)') 'runtime(sec) = ',endtime-starttime
         write(6,'(a)') ' '

      endif

c-----------------------------------------------------------------------
c     Bottom of routine
c-----------------------------------------------------------------------

      cmd = 'rm fixme*'
      call system(cmd)
      if((icstate.eq.3).or.(icstate.eq.4)) then
         cmd = 'rm pcmod*'
         call system(cmd)
         cmd = 'rm conpcm*'
         call system(cmd)
         cmd = 'rm scratch'
         call system(cmd)
         cmd = 'rm outfile.pcm'
         call system(cmd)
      end if
 
c-----------------------------------------------------------------------
c     Bottom of routine
c-----------------------------------------------------------------------

      return
      end
