c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     searchme -- reads in a pcmodel formatted structure and creates
c     the conpcm file needed to run a conformer search
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      program searchme
      
      implicit none
      
      include 'params.i'

      integer i,j
      integer n,n12(MAXATOM),i12(MAXVAL,MAXATOM),itype(MAXATOM)
      integer bo12(MAXVAL,MAXATOM)
      integer nsb,isb(MAXSB),jsb(MAXSB),sincr(MAXSB)
      integer ifield
      
      real x(MAXATOM),y(MAXATOM),z(MAXATOM)

      character*1 pichar(MAXATOM),addpar
      character*80 filename,outfile,mmpar

      logical flag,srchring,norun

c-----------------------------------------------------------------------
c     Initialize 
c-----------------------------------------------------------------------

      flag = .false.
      srchring = .false.
      norun = .false.
      ifield = 0
      n = 0
      do i = 1, MAXATOM 
         x(i) = 0.00e0
         y(i) = 0.00e0
         z(i) = 0.00e0
         itype(i) = 0
         pichar(i) = ' '
         n12(i) = 0
         do j = 1, MAXVAL
            i12(j,i) = 0
            bo12(j,i) = 0
         end do
      end do

c-----------------------------------------------------------------------
c     Get input pcmodel file name
c-----------------------------------------------------------------------

      write(6,'(a)') 'Enter name of PCModel structure file'
      read(5,'(a80)') filename
      inquire(file=filename,exist=flag)
      if(.not.flag) then
         write(6,'(a)') 'PCModel structure file not found'
         write(6,'(a)') ' '
         stop
      end if

c-----------------------------------------------------------------------
c     Read pcmodel file
c-----------------------------------------------------------------------

      call readpcm(filename,ifield,n,x,y,z,itype,n12,i12,bo12,pichar)

c-----------------------------------------------------------------------
c     Check ifield value
c-----------------------------------------------------------------------

      if((ifield.eq.3).or.(ifield.eq.7)) then
         continue
      else
         write(6,'(a,i1,a)') 'ATOMTYPES ',ifield,' not supported'
         write(6,'(a)') 'Allowed values are:'
         write(6,'(a)') '  ATOMTYPES 3      (MM3)'
         write(6,'(a)') '  ATOMTYPES 7     (MMFF94)'
         write(6,'(a)') ' '
         stop
      end if
         
c-----------------------------------------------------------------------
c     Get search info
c-----------------------------------------------------------------------

      call srchbonds(ifield,n,n12,i12,bo12,nsb,isb,jsb,itype,srchring)

      do i = 1, nsb
         if(n12(isb(i)).eq.4.and.n12(jsb(i)).eq.4) then
            sincr(i) = 120
         else
            sincr(i) = 60
         end if
      end do

c-----------------------------------------------------------------------
c     Warn user if there are too many search bonds
c-----------------------------------------------------------------------

      if(nsb.ge.20) then
         write(6,'(a)') ' '
         write(6,'(a)') 'With 20 or more rotatable bonds, this '//
     &   'structure is too flexible to search'
         write(6,'(a)') 'Exiting code without writing conpcm file'
         write(6,'(a)') ' '
         stop
      end if

c-----------------------------------------------------------------------
c     Tell user molecule is rigid, nothing to search
c-----------------------------------------------------------------------

      if((nsb.eq.0).and.(.not.srchring)) then
         write(6,'(a)') ' '
         write(6,'(a)') 'This molecule is rigid and conformational '//
     &   'search is not necessary'
         write(6,'(a)') 'Exiting code without writing conpcm file'
         write(6,'(a)') ' '
         stop
      end if

c-----------------------------------------------------------------------
c     Ask for output file name, ask if there is an added parameter file
c-----------------------------------------------------------------------

      write(6,'(a)') ' '
      write(6,'(a)') 'Enter name of output file:'
      read(5,*) outfile
      write(6,'(a)') ' '

      addpar = ' '
  333 write(6,'(a)') 'Use an added parameter file:(y/n)'
      read(5,*) addpar
      if(addpar.eq.'y') then
         write(6,'(a)') 'Enter name of added parameter file:'
         read(5,*) mmpar
      elseif(addpar.eq.'n') then
         mmpar = '_nogo'
      else
         goto 333
      end if         

c-----------------------------------------------------------------------
c     Write conpcm file
c-----------------------------------------------------------------------

      call writesrchcon(filename,outfile,mmpar,ifield,nsb,isb,jsb,
     &  sincr,srchring)

c-----------------------------------------------------------------------
c     Run mengine
c-----------------------------------------------------------------------

      if(nsb.ge.9) then
         write(6,'(a)') ' '
         write(6,'(a)') 'WARNING: FLEXIBLE STRUCTURE WITH 9 OR MORE'//
     &   ' ROTATABLE BONDS'
         write(6,'(a)') 'Although conpcm file has been written, '//
     &   'conformer search was not started'
         write(6,'(a)') 'If you choose to run this search, then expect'//
     &   ' run time to exceed 30 min'
         write(6,'(a)') ' '
         norun = .true.
         write(6,'(a)') ' '
      elseif(nsb.ge.12) then
         write(6,'(a)') ' '
         write(6,'(a)') 'WARNING: FLEXIBLE STRUCTURE WITH 12 or MORE'//
     &   ' ROTATABLE BONDS'
         write(6,'(a)') 'Although conpcm file has been written, '//
     &   'conformer search was not started'
         write(6,'(a)') 'If you choose to run this search, then expect'//
     &   ' run time to exceed 4 hr'
         write(6,'(a)') ' '
         norun = .true.
      end if
     
      if(.not.norun) then
         call system('mengine')
         call system('rm scratch*')
      end if

      end
