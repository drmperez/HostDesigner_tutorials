c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'writesrchcon' -- writes control file for mengine conf search
c
c     passed variables:
c     char     input,output    names of input and output file
c     char     mmpar           name of added parameter file
c     integer  ifield          mm3=3, mmff94=7
c     integer  nsb             number of bonds to search
c     integer  isb             serial number for ith atom in bond
c     integer  jsb             serial number for jth atom in bond
c     integer  sincr           increment to drive the bond
c     logical  srchring       
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine writesrchcon(infile,outfile,mmpar,ifield,nsb,isb,jsb,
     &   sincr,srchring)

      implicit none
      include 'params.i'
      
      integer i,nsb,isb(*),jsb(*),sincr(*),ifield
      integer lenin,lenout,lenpar
      character*80 infile,outfile,mmpar
      logical srchring

      open(unit=10,file='conpcm',status='unknown')
            
c-----------------------------------------------------------------------
c     Get lengths of input and output strings
c-----------------------------------------------------------------------

      lenin = 80
      do i = 80,1,-1
         if(infile(i:i).eq.' ') lenin = lenin - 1      
      end do

      lenout = 80
      do i = 80,1,-1
         if(outfile(i:i).eq.' ') lenout = lenout - 1      
      end do

      lenpar = 80
      do i = 80,1,-1
         if(mmpar(i:i).eq.' ') lenpar = lenpar - 1      
      end do

c-----------------------------------------------------------------------
c     Set mode to search
c-----------------------------------------------------------------------
    
      write(10,'(a)') 'mode search'
      
c-----------------------------------------------------------------------
c     Input and output file names for first cycle
c-----------------------------------------------------------------------
      
      write(10,'(2a)') 'infile ',infile(1:lenin)
      write(10,'(2a)') 'outfile scratch'
      
c-----------------------------------------------------------------------
c     Force field
c-----------------------------------------------------------------------
      
      if(ifield.eq.3) then
         write(10,'(a)') 'forcefield mm3'
      else if(ifield.eq.7) then
         write(10,'(a)') 'forcefield mmff94'
      end if
      
c-----------------------------------------------------------------------
c     Name of added parameter file
c-----------------------------------------------------------------------

      if(mmpar(1:5).ne.'_nogo') then
         write(10,'(2a)') 'addpar ', mmpar(1:lenpar)
      end if

c-----------------------------------------------------------------------
c     Dielec: resets dielectric constant from default of 1.5 for 
c     MMX and MM3 or 1.0 for other force fields.
c-----------------------------------------------------------------------

      if(ifield.eq.3) then
         write(10,'(a)') 'dielec 4.0'
      elseif(ifield.eq.7) then
         write(10,'(a)') 'dielec 1.0'
      end if
      
c-----------------------------------------------------------------------
c     Set maximum number of trial structures used in search 
c-----------------------------------------------------------------------

         write(10,'(a)') 'max 5000'

c-----------------------------------------------------------------------
c     Method: conformer search method, 1=search by bonds only,
c     2=search by cartesians only, 3=search by both methods(default)
c-----------------------------------------------------------------------

      if(srchring) then
         if(nsb.eq.0) then
            write(10,'(a)') 'method 2'
         else if(nsb.gt.0) then
            write(10,'(a)') 'method 3'
         end if
      else
         if(nsb.eq.0) then
            write(10,'(a)') 'method 2'
         else if(nsb.gt.0) then
            write(10,'(a)') 'method 1'
         end if
      end if
      
c-----------------------------------------------------------------------
c     Ewindow: energy cutoff for first pass, default = 3.5 kcal/mol
c-----------------------------------------------------------------------

      write(10,'(a,f5.2)') 'ewindow 3.5'

c-----------------------------------------------------------------------
c     Searchable bonds
c-----------------------------------------------------------------------

      write(10,'(a,i3)') 'NumBonds ',nsb
      do i = 1, nsb
         write(10,'(a,x,i4,x,i4,x,i4)') 'Bond',isb(i),jsb(i),sincr(i)
      end do

c-----------------------------------------------------------------------
c     Variables affecting second cycle
c-----------------------------------------------------------------------

      write(10,'(a)') 'Second_cycle'

c-----------------------------------------------------------------------
c     Name of input and output files for second cycle
c-----------------------------------------------------------------------
      
      write(10,'(2a)') 'infile scratch '
      write(10,'(2a)') 'outfile ',outfile(1:lenout)
      
c-----------------------------------------------------------------------
c     Ewindow: energy cutoff for second pass, default = 3.5 kcal/mol
c-----------------------------------------------------------------------

      write(10,'(a,f5.2)') 'ewindow 3.5'
      
      close(10)

      return
      end
