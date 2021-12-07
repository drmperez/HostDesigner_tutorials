c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     searchme -- creates conpcm file and runs conformer search
c     returns energy for global minimum and writes pcmodel file named
c     global_min.pcm containing structure
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine searchme(infile,mmpar,ifield,n,x,y,z,itype,n12,i12,
     & bo12,pichar,Efree)
      
      implicit none
      
      include 'params.i'

      integer i
      integer n,n12(MAXATOM),i12(MAXVAL,MAXATOM),itype(MAXATOM)
      integer bo12(MAXVAL,MAXATOM)
      integer nsb,isb(MAXSB),jsb(MAXSB),sincr(MAXSB)
      integer ifield,lenpar,lenin
      
      real x(MAXATOM),y(MAXATOM),z(MAXATOM),IE,Efree

      character*1 pichar(MAXATOM)
      character*80 infile,outfile,mmpar

      logical flag,srchring,norun

c-----------------------------------------------------------------------
c     Initialize 
c-----------------------------------------------------------------------

      flag = .false.
      srchring = .false.
      norun = .false.

      lenin = 80
      do i = 80,1,-1
         if(infile(i:i).eq.' ') lenin = lenin - 1
      end do

      lenpar = 80
      do i = 80,1,-1
         if(mmpar(i:i).eq.' ') lenpar = lenpar - 1
      end do

c-----------------------------------------------------------------------
c     Check ifield value
c-----------------------------------------------------------------------

      if(ifield.eq.3) then
         continue
      else
         write(6,'(a,i1,a)') 'ATOMTYPES ',ifield,' not supported'
         write(6,'(a)') 'Allowed values are:'
         write(6,'(a)') '  ATOMTYPES 3      (MM3)'
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
c     Stop if there are too many search bonds
c-----------------------------------------------------------------------

      if(nsb.ge.9) then
         write(6,'(a)') ' '
         write(6,'(a)') 'With 9 or more rotatable bonds, this '//
     &   'structure is too flexible to search'
         write(6,'(a)') 'Exiting code without writing conpcm file'
         write(6,'(a)') ' '
         stop
      end if

c-----------------------------------------------------------------------
c     If molecule is rigid, nothing to search, optimize and get energy
c-----------------------------------------------------------------------

      if((nsb.eq.0).and.(.not.srchring)) then
         open(unit=10,file='conpcm')
            write(10,'(a)') 'mode opt'
            write(10,'(2a)') 'infile ',infile(1:lenin)
            write(10,'(a)') 'outfile global_min.pcm'
            write(10,'(a)') 'forcefield mm3'
            if(mmpar(1:5).ne.'_nogo') then
               write(10,'(2a)') 'addpar ', mmpar(1:lenpar)
            end if
            write(10,'(a)') 'dielec 4.0'
         close(10)
         call system('mengine')
         outfile = 'global_min.pcm'
         call readpcm(outfile,ifield,n,x,y,z,itype,n12,i12,bo12,pichar,
     &    IE,Efree)
c-----------------------------------------------------------------------
c     else run search and get energy
c-----------------------------------------------------------------------
      else
         outfile = 'global_min.pcm'
         call writesrchcon(infile,outfile,mmpar,ifield,nsb,isb,jsb,
     &     sincr,srchring)
         open(unit=11,file='fixme',status='unknown')
            write(11,'(a)') 'limit cputime 10000'
            write(11,'(a)') 'mengine conpcm'
         close(11)
         call system('csh < fixme')
         call system('rm scratch*')
         call system('rm fixme')
         call readpcm(outfile,ifield,n,x,y,z,itype,n12,i12,bo12,pichar,
     &    IE,Efree)
      end if

      return
      end
