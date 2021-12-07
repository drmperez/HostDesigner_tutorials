c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      params.i  
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      integer MAXVAL     ! maximum valence, number of bonds per atom
      integer MAXATOM    ! maximum number of atoms in link
      integer MAXSB      ! maximum number of searchable bonds
      integer MAXMET     ! number of metal ions in data strutures

      parameter (MAXVAL = 10)
      parameter (MAXATOM = 200)
      parameter (MAXSB = 100)
      parameter (MAXMET = 25)
