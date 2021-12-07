c***********************************************************************
c***********************************************************************
      subroutine writesrchcon(infile,outfile,nsb,isb,jsb,sincr,srchring)
c***********************************************************************
c***********************************************************************

      implicit none

      integer i
      integer nsb,isb(*),jsb(*),sincr(*)

      character*40 infile,outfile

      logical srchring

c***********************************************************************
c     open control file
c***********************************************************************

      open(unit=20, file='conpcm')
     
c***********************************************************************
c     set mode 
c***********************************************************************

      if((nsb.eq.0).and.(.not.srchring)) then
         write(20,'(a)') 'mode opt'
      else
         write(20,'(a)') 'mode search'
      end if

c***********************************************************************
c     name of input and output files for first cycle
c***********************************************************************

      write(20,'(2a)') 'infile ',infile

      if((nsb.eq.0).and.(.not.srchring)) then
         write(20,'(2a)') 'outfile ', outfile
      else
         write(20,'(a)') 'outfile scratch'
      end if

c***********************************************************************
c     level of printing
c***********************************************************************

      write(20,'(a)') 'iprint 0'
      write(20,'(a)') 'jprint 0'

c***********************************************************************
c     force field
c***********************************************************************

      write(20,'(a)') 'forcefield mm3'
      write(20,'(a)') 'dielec 4.0'

c***********************************************************************
c     search control
c***********************************************************************

      if((nsb.gt.0).or.(srchring)) then

c***********************************************************************
c     method: conformer search method, 1=search by bonds only,
c     2=search by cartesians only, 3=search by both methods(default)
c     maxpush: maximum number of conformers minimized, default = 100000
c     dopi: must  be there to get engine to run pi calculations
c     chig: maintain chirality
c***********************************************************************

         if(srchring) then
            if(nsb.eq.0) then
               write(20,'(a)') 'method 2'
            else if(nsb.gt.0) then
               write(20,'(a)') 'method 3'
            end if
         else
            if(nsb.eq.0) then
               write(20,'(a)') 'method 2'
            else if(nsb.gt.0) then
               write(20,'(a)') 'method 1'
            end if
         end if

c***********************************************************************
c     Maxpush: maximum number of conformers minimized, default = 100000
c***********************************************************************

         write(20,'(a,i6)') 'max ',500

c***********************************************************************
c     ewindow1: energy cutoff for first pass
c***********************************************************************

         write(20,'(a)') 'ewindow 5.00'

c***********************************************************************
c     which non-ring bonds are we rotating
c***********************************************************************

         write(20,'(a,i3)') 'NumBonds ', nsb

         do i = 1, nsb
            write(20,'(a, x, i4, x, i4, x, i4)') 'Bond', isb(i),
     &      jsb(i), sincr(i)
         enddo

c***********************************************************************
c     variables affecting second cycle
c***********************************************************************

         write(20,'(a)') 'Second_cycle'

c***********************************************************************
c     name of input and output files for second cycle
c***********************************************************************

         write(20,'(a)') 'infile scratch '
         write(20,'(2a)') 'outfile ',outfile

c***********************************************************************
c     ewindow2: energy cutoff for second pass
c***********************************************************************

         write(20,'(a)') 'ewindow 5.00'

      end if

      close(20)
      return
      end
