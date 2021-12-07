c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Computes ligand strain in lanthanide complexes
c     Input - single ligand bound to lanthanide metal ion
c     Output - global_min.pcm contains results of ligand conf search
c            - table with entries of metal, Shannon ionic radius at 
c              CN=8, and ligand strain where ligand strain =
c              E(bound) - E(minimum) in kcal/mol
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      program scanme

      implicit none
      include 'params.i'

      integer i,j,ii,spot
      integer n,n12(MAXATOM),i12(MAXVAL,MAXATOM),itype(MAXATOM)
      integer nu,n12u(MAXATOM),i12u(MAXVAL,MAXATOM),itypeu(MAXATOM)
      integer bo12(MAXVAL,MAXATOM),bo12u(MAXVAL,MAXATOM),ifield
      integer met_type(MAXMET),lenpar

      real x(MAXATOM),y(MAXATOM),z(MAXATOM)
      real xu(MAXATOM),yu(MAXATOM),zu(MAXATOM),IE,FE,Efree
      real met_radius(MAXMET),strain(MAXMET),Ebound(MAXMET)

      character*1 yesno,pichar(MAXATOM),picharu(MAXATOM)
      character*2 met_label(MAXMET)
      character*50 title
      character*80 filename,infile,outfile,mmpar

      logical flag,dosearch

c----------------------------------------------------------------------
c     initialize  
c-----------------------------------------------------------------------

      flag = .false.
      dosearch = .true.
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
      nu = 0
      do i = 1, MAXATOM
         xu(i) = 0.00e0
         yu(i) = 0.00e0
         zu(i) = 0.00e0
         itypeu(i) = 0
         picharu(i) = ' '
         n12u(i) = 0
         do j = 1, MAXVAL
            i12u(j,i) = 0
            bo12u(j,i) = 0
         end do
      end do

      data met_label/  'La',    'Ce',    'Pr',    'Nd',    'Pm',
     &                 'Sm',    'Eu',    'Gd',    'Tb',    'Dy',
     &                 'Ho',    'Er',    'Tm',    'Yb',    'Lu',
     &                 'Be',    'Mg',    'Ca',    'Sr',    'Ba',
     &                 'Li',    'Na',    'K ',    'Rb',    'Cs'/
      data met_radius/ 1.160,   1.143,   1.126,   1.109,   1.093,
     &                 1.079,   1.066,   1.053,   1.040,   1.027,
     &                 1.015,   1.004,   0.994,   0.985,   0.977,
     &                 0.450,   0.720,   1.000,   1.180,   1.350,
     &                 0.760,   1.020,   1.380,   1.520,   1.670/
      data met_type/   346,     347,     349,     350,     351,
     &                 352,     354,     355,     356,     357,
     &                 358,     359,     360,     361,     362,
     &                 302,     304,     307,     328,     345,  
     &                 301,     303,     306,     327,     344/

c----------------------------------------------------------------------
c     get name of input file and added parameter file if used 
c-----------------------------------------------------------------------

      write(6,'(a)') 'Enter name of PCModel input structure file'
      read(5,'(a80)') filename
      inquire(file=filename,exist=flag)
      if(.not.flag) then
         write(6,'(a)') 'PCModel input structure file not found'
         write(6,'(a)') ' '
         stop
      end if

      yesno = ' '
    5 write(6,'(a)') 'Conformer search the free ligand? (y/n)'
      read(5,'(a1)') yesno
      if(yesno.eq.'y') then
         dosearch = .true. 
      elseif(yesno.eq.'n') then
         dosearch = .false. 
      else
        goto 5
      end if

      mmpar = '_nogo'
      lenpar = 80
c     yesno = ' '
c     lenpar = 0
c  10 write(6,'(a)') 'Read in added parameter file:(y/n)'
c     read(5,'(a1)') yesno
c     if(yesno.eq.'y') then
c        write(6,'(a)') 'Enter name of added parameter file:'
c        read(5,*) mmpar
c        lenpar = 80
c        do i = 80,1,-1
c           if(mmpar(i:i).eq.' ') lenpar = lenpar - 1
c        end do
c     elseif(yesno.eq.'n') then
c        mmpar = '_nogo'
c     else
c        goto 10
c     end if

c----------------------------------------------------------------------
c     read PCModel input file
c     original molecule description goes in normal arrays
c     modified molecule (loss of metal) goes in suffix u arrays
c-----------------------------------------------------------------------

      call readpcm(filename,ifield,n,x,y,z,itype,n12,i12,bo12,pichar,
     & IE,FE)

      spot = 0
      nu = n
      do i = 1, nu
         xu(i) = x(i)
         yu(i) = y(i)
         zu(i) = z(i)
         itypeu(i) = itype(i)
         if(itypeu(i).ge.300) spot = i
         picharu(i) = pichar(i)
         n12u(i) = n12(i)
         do j = 1, n12u(i)
            i12u(j,i) = i12(j,i)
            bo12u(j,i) = bo12(j,i)
         end do
      end do

      if(spot.eq.0) then
         write(6,'(a)') ' '
         write(6,'(a)') 'Input structure does not contain metal ion'
         stop
      end if

c----------------------------------------------------------------------
c     remove metal, write infile, run searchme to get E(free)
c-----------------------------------------------------------------------

      if(dosearch) then

         call delatom(spot,nu,itypeu,n12u,i12u,bo12u,xu,yu,zu,picharu)
         infile = 'srch_input.pcm'
         title = ' ligand without metal ion'
         call writepcm(infile,title,nu,xu,yu,zu,itypeu,n12u,i12u,bo12u,
     &    picharu)

         write(6,'(a)') ' '
         write(6,'(a)') 'Starting conformer search . . .'
         write(6,'(a)') ' '

         call searchme(infile,mmpar,ifield,nu,xu,yu,zu,itypeu,n12u,i12u,
     &    bo12u,picharu,Efree)
         call system('rm srch_input.pcm')

      else

         call delatom(spot,nu,itypeu,n12u,i12u,bo12u,xu,yu,zu,picharu)
         infile = 'input.pcm'
         title = ' ligand binding conformer'
         call writepcm(infile,title,nu,xu,yu,zu,itypeu,n12u,i12u,bo12u,
     &    picharu)

         open(unit=20,file='conpcm')
            write(20,'(a)') 'mode opt'
            write(20,'(a)') 'infile input.pcm'
            write(20,'(a)') 'outfile binding_form.pcm'
            write(20,'(a)') 'forcefield mm3'
            if(mmpar(1:5).ne.'_nogo') then
               write(10,'(2a)') 'addpar ', mmpar(1:lenpar)
            end if
            write(20,'(a)') 'dielec 4.0'
         close(20)

         call system('mengine')
         outfile = 'binding_form.pcm'
         call readpcm(outfile,ifield,nu,xu,yu,zu,itypeu,n12u,i12u,
     &    bo12u,picharu,IE,Efree)
         call system('rm input.pcm')

c debug
      print*,'binding form energy = ',Efree
 
      end if

c----------------------------------------------------------------------
c     loop over lanthanides, calc E(bound)
c-----------------------------------------------------------------------

      do ii = 1, 15

c----------------------------------------------------------------------
c     load original input molecule into working arrays
c-----------------------------------------------------------------------

         nu = n
         do i = 1, nu
            xu(i) = x(i)
            yu(i) = y(i)
            zu(i) = z(i)
            itypeu(i) = itype(i)
            picharu(i) = pichar(i)
            n12u(i) = n12(i)
            do j = 1, n12u(i)
               i12u(j,i) = i12(j,i)
               bo12u(j,i) = bo12(j,i)
            end do
         end do
 
c----------------------------------------------------------------------
c     update metal atom type and input file
c-----------------------------------------------------------------------

         infile = 'swapmetal.pcm'
         outfile = 'output.pcm'

         itypeu(spot) = met_type(ii)
         call writepcm(infile,title,nu,xu,yu,zu,itypeu,n12u,i12u,bo12u,
     &    picharu)

c----------------------------------------------------------------------
c     write conpcm file for optimization
c-----------------------------------------------------------------------

         open(unit=20,file='conpcm')
            write(20,'(a)') 'mode opt'
            write(20,'(a)') 'infile swapmetal.pcm'
            write(20,'(a)') 'outfile output.pcm'
            write(20,'(a)') 'forcefield mm3'
            if(mmpar(1:5).ne.'_nogo') then
               write(10,'(2a)') 'addpar ', mmpar(1:lenpar)
            end if
            write(20,'(a)') 'dielec 4.0'
         close(20)

c----------------------------------------------------------------------
c     optimize and read in updated coordinates
c-----------------------------------------------------------------------

         call system('mengine')
         call readpcm(outfile,ifield,nu,xu,yu,zu,itypeu,n12u,i12u,
     &    bo12u,picharu,IE,FE)

c----------------------------------------------------------------------
c     remove metal ion and update input file
c-----------------------------------------------------------------------

         infile = 'sansmetal.pcm'
         outfile = 'spe.pcm'

         call delatom(spot,nu,itypeu,n12u,i12u,bo12u,xu,yu,zu,picharu)
         call writepcm(infile,title,nu,xu,yu,zu,itypeu,n12u,i12u,bo12u,
     &    picharu)

c----------------------------------------------------------------------
c     write conpcm file for single point energy
c-----------------------------------------------------------------------

         open(unit=20,file='conpcm')
            write(20,'(a)') 'mode single'
            write(20,'(a)') 'infile sansmetal.pcm'
            write(20,'(a)') 'outfile spe.pcm'
            write(20,'(a)') 'forcefield mm3'
            if(mmpar(1:5).ne.'_nogo') then
               write(10,'(2a)') 'addpar ', mmpar(1:lenpar)
            end if
            write(20,'(a)') 'dielec 4.0'
         close(20)

         call system('mengine')
         call readpcm(outfile,ifield,nu,xu,yu,zu,itypeu,n12u,i12u,bo12u,
     &    picharu,Ebound(ii),FE)
         strain(ii) = Ebound(ii) - Efree

c debug
      print*,'Ebound = ',Ebound(ii)
      print*,'Ebound - Efree = ',strain(ii)
      print*,' '

         write(6,'(a,a2,a,F10.3)') 'Metal ion ',met_label(ii),
     &    ' ligand strain = ',strain(ii)
         write(6,'(a)') ' '

c----------------------------------------------------------------------
c     bottom of loop over lanthanides
c-----------------------------------------------------------------------

      end do

c----------------------------------------------------------------------
c     print results table
c-----------------------------------------------------------------------
 
      open(unit=30,file='table')
         write(30,'(a)') 'Metal     Radius(CN=8)     Ligand Strain '
         do ii = 1, 15
            write(30,'(x,a2,6x,F10.3,6x,F10.3,6x,F10.3)') met_label(ii),
     &       met_radius(ii),strain(ii)
         end do
      close(30)

c----------------------------------------------------------------------
c     tidy directory 
c-----------------------------------------------------------------------

      call system('rm pcmod*')
      call system('rm spe.pcm')
      call system('rm swapmetal.pcm')
      call system('rm sansmetal.pcm')
      call system('rm output.pcm')

      end
