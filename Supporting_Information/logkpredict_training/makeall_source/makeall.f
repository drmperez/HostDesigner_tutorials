c***********************************************************************
c***********************************************************************
c     Program makeall 
c
c     automates process used to create training data files for all
c     ligands.  input consists of a series of directories, each 
c     containing three files:
c         input                makedata control file
c         free_ligand.pcm      global minimum ligand structure
c         <lig_name>.pcm       metal-ligand complex 
c
c***********************************************************************
c***********************************************************************

      program makeall

      implicit none
      include 'sizes.i'
      include 'metals.i'

      integer i,j,ii,jj,length,start
      integer nlogk,qm(MAXMET),ql,nrot
      integer natomin,itypein(MAXAT)
      integer n12in(MAXAT),i12in(MAXVAL,MAXAT),bo12in(MAXVAL,MAXAT)
      integer n,itype(MAXAT)
      integer n12(MAXAT),i12(MAXVAL,MAXAT),bo12(MAXVAL,MAXAT)
      integer naq,itypeaq(MAXAT)
      integer n12aq(MAXAT),i12aq(MAXVAL,MAXAT),bo12aq(MAXVAL,MAXAT)
      integer nm,itypem(MAXAT)
      integer n12m(MAXAT),i12m(MAXVAL,MAXAT),bo12m(MAXVAL,MAXAT)
      integer miss,fail,spot,mtype,mcoord
      integer nfit,ifit(2,3),ndonor
      integer nbond,ibond(2*MAXAT),jbond(2*MAXAT),bo(2*MAXAT)
      integer nsb,isb(MAXAT),jsb(MAXAT),ibo(MAXAT),idt1,idt2
      integer nchg,zatom(10),zchg(10)
      integer nsd,nn

      real logk(MAXMET),logkin(MAXMET),e1,e2
      real istrength(MAXMET),mrdhe,mrdhc,msolv
      real xin(MAXAT),yin(MAXAT),zin(MAXAT)
      real x(MAXAT),y(MAXAT),z(MAXAT)
      real xaq(MAXAT),yaq(MAXAT),zaq(MAXAT)
      real xm(MAXAT),ym(MAXAT),zm(MAXAT)
      real dele,eligand,eaquo,msir,rmsd
      real ebound,ligstrain,estart,efinal,innerstrain

      character*1 picharin(MAXAT),pichar(MAXAT),picharaq(MAXAT)
      character*1 picharm(MAXAT)
      character*1 char1
      character*2 matlab(MAXMET)
      character*2 atlabin(MAXAT),atlab(MAXAT),atlabaq(MAXAT)
      character*2 atlabm(MAXAT)
      character*40 ligname
      character*40 pcmodname,filename,subname,subdir(MAXDIR)
      character*50 title
      character*80 cmd

c***********************************************************************
c     load metal data
c***********************************************************************

      call loadmet

c***********************************************************************
c     get names for all ligand subdirectories in working directory
c***********************************************************************


      do i = 1, MAXDIR
         subdir(i) = '                                        '
      end do

      call system('ls > junk')

      nsd = 0
      open(unit=20,file='junk')
    5    continue
         read(20,'(a40)',end=10,err=10) subname
         if(subname(1:4).ne.'junk') then
            nsd = nsd + 1
            subdir(nsd) = subname
         end if
         goto 5
   10 continue
      close(20)
      call system('rm junk')

c***********************************************************************
c     loop over ligand subdirectories
c***********************************************************************

      do nn = 1, nsd 

         write(6,'(2a)') 'working on ',subdir(nn)

         length = 40
         do i = 40, 1, -1
            if(subdir(nn)(i:i).eq.' ') length = length - 1
         end do

         cmd = 'mv '//subdir(nn)(1:length)//'/'//
     &    subdir(nn)(1:length)//'+M.pcm .'
         call system(cmd)
         cmd = 'mv '//subdir(nn)(1:length)//'/free_ligand.pcm .'
         call system(cmd)
         cmd = 'mv '//subdir(nn)(1:length)//'/input .'
         call system(cmd)

c***********************************************************************
c     load metal data and print out table
c***********************************************************************

c***********************************************************************
c***********************************************************************
c     read input file:
c     pcmodname         name of pcmod input file    40 char string
c     ligname           name of ligand              40 char string
c     ql                charge on ligand            integer
c     nrot              number restricted rots      integer
c     idt1,idt2         donor types                 integer
c     nlogk             # of log K values           integer
c     matlab()          metal label                 2 char string
c     qm()              metal ion charge            integer
c     logk(i)           1:1 formation constant      real
c     istrength(i)      ionic strength              real
c***********************************************************************
c***********************************************************************

      do i = 1, MAXMET
         matlab(i) = '  '
         logk(i) = 0.0e0
         logkin(i) = 0.0e0
         istrength(i) = 0.0e0
         qm(i) = 0
      end do
      pcmodname = '                                        '
      ligname = '                                        '
      nlogk = 0
      ql = 0
      nrot = 0
      idt1 = 0
      idt2 = 0

      open(unit=10,file='input')
         read(10,*,end=999,err=999) pcmodname
         read(10,*,end=999,err=999) ligname,ql,nrot
         read(10,*) nlogk
         if(nlogk.gt.MAXMET) then
            write(6,'(a)') 'Number of logK values exceeds maximum'
            goto 1000
         end if
         do i = 1, nlogk
            read(10,*,end=999,err=999) matlab(i),qm(i),logkin(i),
     &       istrength(i)
         end do
      close(10)

c***********************************************************************
c     apply Davies correction: exp log K –> adjusted log K at I = 0.0
c***********************************************************************

      do i = 1, nlogk
         call davies(ql,qm(i),istrength(i),logkin(i),0.0,logk(i))
      end do

c***********************************************************************
c     read input geometry for an isolated chelate ring
c***********************************************************************

      call readpcm(pcmodname,natomin,itypein,xin,yin,zin,e1,e2,n12in,
     & i12in,bo12in,miss,fail,picharin)
      call atomtyper(natomin,itypein,atlabin)

c***********************************************************************
c     optimize global minimum for ligand and get energy 
c     must already be a file named free_ligand.pcm containing input
c***********************************************************************

      open(unit=11,file='conpcm')
         write(11,'(a)') 'mode opt'
         write(11,'(a)') 'infile free_ligand.pcm'
         write(11,'(a)') 'outfile junk.pcm'
         write(11,'(a)') 'forcefield mm3'
         write(11,'(a)') 'dielec 4.0'
      close(11)
      call system('mengine')
      filename = 'junk.pcm'
      call readpcm(filename,n,itype,x,y,z,e1,eligand,n12,i12,bo12,
     & miss,fail,pichar)
      call system('rm junk.pcm')

c***********************************************************************
c     write global minimum ligand geometry to ligand.xyz file 
c***********************************************************************

      call atomtyper(n,itype,atlab)

      length = 40
      do i = 40, 1, -1
         if(ligname(i:i).eq.' ') length = length - 1
      end do
      filename = ligname(1:length)//'.xyz'
      open(unit=15,file=filename)
         write(15,'(i5)') n
         write(15,'(a)') ligname(1:length)
         do i = 1, n
            write(15,'(a2,3F10.4)') atlab(i),x(i),y(i),z(i)
         end do
      close(15)

c***********************************************************************
c***********************************************************************
c     loop over metal ion list
c***********************************************************************
c***********************************************************************

      do ii = 1, nlogk

c***********************************************************************
c     properties of the iith metal
c***********************************************************************

         write(char1,'(i1)') qm(ii)
         mtype = 0
         do i = 1, MAXMTYPE
            if((matlab(ii).eq.metatlab(i)).and.
     &      (qm(ii).eq.metcharge(i))) then
               mtype = mettype(i)
               mcoord = metcn(i)
               msir = sir(i)
               mrdhe = rdhe(i)
               mrdhc = rdhc(i)
               msolv = metsolv(i)
            end if
         end do
         if(mtype.eq.0) then
            write(6,'(a)') 'Unable to define atom type for metal'
            write(6,'(3a,i1)') 'with label ',matlab(ii),'and charge ',
     &       qm(ii)
            goto 1000
         end if
         if((mrdhe.gt.9000).or.(mrdhc.gt.9000).or.(msolv.gt.9000)) then
            write(6,'(a,a2,i3)') 'Skipping metal ',matlab(ii),qm(ii)
            write(6,'(a)') 'requisite parameters not available'
            goto 500          
         end if
            
c***********************************************************************
c     create metal aquo ion pcmod input file
c***********************************************************************

         filename = 'infile.pcm'
         call makeaquo(mtype,mcoord,filename)

c***********************************************************************
c     optimize aquo ion
c***********************************************************************
         
         if(matlab(ii)(2:2).eq.' ') then
            filename = matlab(ii)(1:1)//'__'//char1//'+_aquo_ion.pcm'
         else
            filename = matlab(ii)//'_'//char1//'+_aquo_ion.pcm'
         end if
         open(unit=11,file='conpcm')
            write(11,'(a)') 'mode opt'
            write(11,'(a)') 'infile infile.pcm'
            write(11,'(2a)') 'outfile ',filename 
            write(11,'(a)') 'forcefield mm3'
            write(11,'(a)') 'dielec 4.0'
         close(11)
         call system('mengine')

         call readpcm(filename,naq,itypeaq,xaq,yaq,zaq,e1,eaquo,n12aq,
     &    i12aq,bo12aq,miss,fail,picharaq)

         do i = 1, naq
            if(itypeaq(i).ge.300) then
               atlabaq(i) = matlab(ii)
            elseif(itypeaq(i).eq.201) then
               atlabaq(i) = 'Lp'
            else
               atlabaq(i) = 'O '
            end if
         end do

c***********************************************************************
c     create hydrated  metal-ligand complex
c     algorithm assumes that input structure consists of single
c     metal ion + ligand, no other atoms present
c***********************************************************************
c***********************************************************************
c     initialize working arrays, move input into working arrays
c***********************************************************************

         do i = 1, MAXAT
            itype(i) = 0
            n12(i) = 0
            do j = 1, MAXVAL
               i12(j,i) = 0
               bo12(j,i) = 0
            end do
            x(i) = 0.0e0
            y(i) = 0.0e0
            z(i) = 0.0e0
            pichar(i) = ' '
            atlab(i) = '  '
         end do

         do i = 1, natomin
            itype(i) = itypein(i)
            n12(i) = n12in(i)
            do j = 1, n12(i)
               i12(j,i) = i12in(j,i)
               bo12(j,i) = bo12in(j,i)
            end do
            x(i) = xin(i)
            y(i) = yin(i)
            z(i) = zin(i)
            pichar(i) = picharin(i)
            atlab(i) = atlabin(i)
         end do

         n = natomin

         spot = 0
         do i = 1, n
            if(itype(i).ge.300) spot = i
         end do
         if(spot.eq.0) then
            write(6,'(a)') 'Error makedata.f'
            write(6,'(a)') 'Unable to find a metal atom type > 300'
            goto 1000
         end if
         itype(spot) = mtype
         atlab(spot) = matlab(ii)

c***********************************************************************
c      get superimposable atom serial numbers for metal-ligand complex
c         if n12(metal)=2, ifit(1,1) = metal, ifit(1,2) = donor 1
c         if n12(metal)=3, adds ifit(1,3) = donor 2
c      and for aquo ion complex, assumes metal=#1, donor1=#2, donor2=#3 
c***********************************************************************

         if(n12(spot).eq.1) then
            ndonor = 1
            nfit = 2
            ifit(1,1) = spot                         ! complex metal ion
            ifit(1,2) = i12(1,spot)                    ! ligand donor #1
            ifit(2,1) = 1                          ! lost aquo metal ion
            ifit(2,2) = 2                                  ! lost OH2 #1
         elseif(n12(spot).eq.2) then
            ndonor = 2
            nfit = 3
            ifit(1,1) = spot                         ! complex metal ion
            ifit(1,2) = i12(1,spot)                    ! ligand donor #1
            ifit(1,3) = i12(2,spot)                    ! ligand donor #2
            ifit(2,1) = 1                          ! lost aquo metal ion
            ifit(2,2) = 2                                  ! lost OH2 #1
            ifit(2,3) = 3                                  ! lost OH2 #2
         else
            write(6,'(a)') 'Exceeded maximum attachments to metal'
            write(6,'(a,i2)') 'allowed values 2 or 3, read ',n12(spot)
         end if 
         
c***********************************************************************
c      superimpose atoms
c***********************************************************************

         call superimpose(n,x,y,z,naq,xaq,yaq,zaq,nfit,ifit,rmsd)

c***********************************************************************
c      strip out extra aquo ion atoms, fix connectivity
c***********************************************************************

         start = n12(spot) + 2                           ! either 3 or 4
         do i = start, naq
            n = n + 1
            x(n) = xaq(i)
            y(n) = yaq(i)
            z(n) = zaq(i)
            itype(n) = itypeaq(i)
            n12(n) = 1
            i12(1,n) = spot
            bo12(1,n) = bo12aq(1,i) 
            pichar(n) = ' '
            atlab(n) = atlabaq(i)
            n12(spot) = n12(spot) + 1
            i12(n12(spot),spot) = n
            bo12(n12(spot),spot) = bo12aq(1,i)
         end do

c***********************************************************************
c     write pcmodel input file, optimize metal complex
c***********************************************************************

         filename = 'infile.pcm'
         title =' input for metal complex'
         call writepcm(filename,title,n,x,y,z,itype,n12,i12,
     &    bo12,pichar)

         if(matlab(ii)(2:2).eq.' ') then
            filename = matlab(ii)(1:1)//'__'//char1//'+_'//
     &      ligname(1:length)//'.pcm'
         else
            filename = matlab(ii)//'_'//char1//'+_'//
     &      ligname(1:length)//'.pcm'
         end if
         title = ' metal complex' 

         open(unit=11,file='conpcm')
            write(11,'(a)') 'mode opt'
            write(11,'(a)') 'infile infile.pcm'
            write(11,'(2a)') 'outfile ',filename 
            write(11,'(a)') 'forcefield mm3'
            write(11,'(a)') 'dielec 4.0'
         close(11)
         call system('mengine')

         call readpcm(filename,n,itype,x,y,z,e1,e2,n12,i12,bo12,
     &    miss,fail,pichar)

c***********************************************************************
c     store optimized metal complex
c***********************************************************************

         do i = 1, MAXAT
            itypem(i) = 0
            n12m(i) = 0
            do j = 1, MAXVAL
               i12m(j,i) = 0
               bo12m(j,i) = 0
            end do
            xm(i) = 0.0e0
            ym(i) = 0.0e0
            zm(i) = 0.0e0
            picharm(i) = ' '
            atlabm(i) = '  '
         end do

         nm = n
         do i = 1, nm
            itypem(i) = itype(i)
            n12m(i) = n12(i)
            do j = 1, n12m(i)
               i12m(j,i) = i12(j,i)
               bo12m(j,i) = bo12(j,i)
            end do
            xm(i) = x(i)
            ym(i) = y(i)
            zm(i) = z(i)
            picharm(i) = pichar(i)
            atlabm(i) = atlab(i)
         end do

c***********************************************************************
c     compute ligand strain:
c     remove everything except ligand atoms
c     do a single point, ebound is energy of bound form
c     eligand is energy of ligand global minimum 
c     ligstrain = ebound - eligand
c***********************************************************************

c***********************************************************************
c     strip metal from structure, assumes that structure is limited to
c     ligand and one metal ion, and write pcmodel input file for ligand
c***********************************************************************

         spot = 0
         do i = 1, n
            if(itype(i).gt.300) then
               mtype = itype(i)
               spot = i
            end if
         end do
         if(spot.eq.0) then
            write(6,'(a)') 'Error makedata.f'
            write(6,'(a)') 'Unable to find a metal atom type > 300'
            goto 1000
         end if

         call delatom(spot,n,itype,n12,i12,bo12,x,y,z,atlab,pichar)

c***********************************************************************
c     remove atoms with no attachments
c***********************************************************************

   20    continue
         spot = 0
         do i = 1, n
            if(n12(i).eq.0) spot = i
         end do
         if(spot.ne.0) then
            call delatom(spot,n,itype,n12,i12,bo12,x,y,z,atlab,pichar)
            goto 20 
         end if

c***********************************************************************
c     optimize bound ligand, read ebound
c***********************************************************************

         filename = 'infile.pcm'
         title =' bound ligand geometry'
         call writepcm(filename,title,n,x,y,z,itype,n12,i12,bo12,pichar)

         open(unit=11,file='conpcm')
            write(11,'(a)') 'mode opt'
            write(11,'(a)') 'infile infile.pcm'
            write(11,'(a)') 'outfile junk.pcm'
            write(11,'(a)') 'forcefield mm3'
            write(11,'(a)') 'dielec 4.0'
         close(11)
         call system('mengine')
         filename = 'junk.pcm'
         call readpcm(filename,n,itype,x,y,z,ebound,e2,n12,i12,bo12,
     &    miss,fail,pichar)
         call system('rm junk.pcm')

         ligstrain = ebound - eligand

         if(ligstrain.gt.30.0e0) then
            print*,'LIGAND STRAIN > 30 kcal/mol'
            print*,'probably something wrong, check input'
            stop
         end if

c***********************************************************************
c     compute inner sphere strain:
c     remove everything except atoms attached to metal ion
c     do an optimizaiton and get initial and final energy
c     innerstrain = estart - efinal
c***********************************************************************

c***********************************************************************
c     keep metal and atoms attached to metal
c***********************************************************************
      
      n = 0
      do i = 1, nm
         if(itypem(i).ge.300) then                         ! found metal
            n = n + 1
            x(n) = xm(i)
            y(n) = ym(i)
            z(n) = zm(i)
            itype(n) = itypem(i)
            atlab(n) = atlabm(i)
            pichar(n) = picharm(i)
            n12(n) = 0
            do j = 1, MAXVAL
               i12(j,i) = 0
               bo12(j,i) = 0
            end do
            do j = 1, n12m(i)
               spot = i12m(j,i)
               n = n + 1
               x(n) = xm(spot)
               y(n) = ym(spot)
               z(n) = zm(spot)
               itype(n) = itypem(spot)
               atlab(n) = atlabm(spot)
               pichar(n) = ' '
               n12(n) = 1
               i12(1,n) = 1
               bo12(1,n) = 1
               n12(1) = n - 1
               i12(n-1,1) = n
               bo12(n-1,1) = 1
            end do 
         end if
      end do
 
c***********************************************************************
c     optimize innersphere, read estart and efinal
c***********************************************************************

         filename = 'infile.pcm'
         title =' inner sphere from metal complex'
         call writepcm(filename,title,n,x,y,z,itype,n12,i12,
     &    bo12,pichar)

         open(unit=11,file='conpcm')
            write(11,'(a)') 'mode opt'
            write(11,'(a)') 'infile infile.pcm'
            write(11,'(a)') 'outfile junk.pcm'
            write(11,'(a)') 'forcefield mm3'
            write(11,'(a)') 'dielec 4.0'
         close(11)
         call system('mengine')
         filename = 'junk.pcm'
         call readpcm(filename,n,itype,x,y,z,estart,efinal,n12,
     &    i12,bo12,miss,fail,pichar)
         call system('rm junk.pcm')

         innerstrain = estart - efinal

         dele = ligstrain + innerstrain

c***********************************************************************
c     write datafile for iith entry
c***********************************************************************

         if(matlab(ii)(2:2).eq.' ') then
            filename = matlab(ii)(1:1)//'__'//char1//'+_'//
     &      ligname(1:length)//'.dat'
         else
            filename = matlab(ii)//'_'//char1//'+_'//
     &      ligname(1:length)//'.dat'
         end if

         open(15,file=filename)

            write(15,'(a)') ' logK(I=0.0)   logK_in   I_in   Z_lig   '//
     &      'Z_met   nrot    met_r   met_CN   E_strain   G_solv    '//
     &      'rdhE     rdhC'
            write(15,300) logK(ii),logKin(ii),istrength(ii),ql,qm(ii),
     &      nrot,msir,mcoord,dele,msolv,mrdhe,mrdhc

  300 FORMAT(X,F8.2,4x,2F8.2,I6,5x,I3,5x,I3,2x,F8.3,5x,I2,4x,F8.2,
     & 2x,F8.1,x,F8.3,x,F8.3)

c***********************************************************************
c     strip out all single atoms attached to metal 
c***********************************************************************

            jj = 1
   30       continue
            if(n12m(jj).eq.1) then
               if(itypem(i12m(1,jj)).gt.300) then
                  call delatom(jj,nm,itypem,n12m,i12m,bo12m,xm,ym,zm,
     &             atlabm,picharm)
                  jj = jj - 1
               end if
            end if
            jj = jj + 1
            if(jj.le.nm) goto 30

c***********************************************************************
c     write SDF mol formatted molecule
c***********************************************************************

            nbond = 0
            do i = 1, nm 
               do j = 1, n12m(i)
                  nbond = nbond + 1
                  bo(nbond) = bo12m(j,i)
                  if(i.lt.i12m(j,i)) then
                     ibond(nbond) = i
                     jbond(nbond) = i12m(j,i)
                  else
                     ibond(nbond) = i12m(j,i)
                     jbond(nbond) = i
                  end if
               end do
            end do
            nsb = 1
            isb(1) = ibond(1)
            jsb(1) = jbond(1)
            ibo(1) = bo(1)
            do i = 1, nbond
               do j = 1, nsb
                  if((isb(j).eq.ibond(i)).and.
     &               (jsb(j).eq.jbond(i))) then
                     goto 40
                  end if
               end do
               nsb = nsb + 1
               isb(nsb) = ibond(i)
               jsb(nsb) = jbond(i)
               ibo(nsb) = bo(i)
   40          continue
            end do

            if((n.gt.999).or.(nsb.gt.999)) then
               if(n.gt.999) then
                  write(6,'(A)') 'skipped write, number '//
     &             'of atoms > sdf format limit of 999'
               else if(nsb.gt.999) then
                  write(6,'(A)') 'skipped write, number '//
     &             'of bonds > sdf format limit of 999'
               end if
               goto 1000
            end if

            if(matlab(ii)(2:2).eq.' ') then
               filename = matlab(ii)(1:1)//'__'//char1//'+_'//
     &         ligname(1:length)
            else
               filename = matlab(ii)//'_'//char1//'+_'//
     &         ligname(1:length)
            end if

            write(15,'(2a)') 'STRUCTURE: ',filename
            write(15,'(a)') 'Comment: MM3 optimized geometry '//
     &       'in MOL format, inner sphere aquo O atoms deleted'
            write(15,'(a)') 'Comment: optimization performed '//
     &       'with MENGINE (HostDesigner, Version 4.3)'
            write(15,200) nm, nsb
            do i = 1, nm
               write(15,201) xm(i),ym(i),zm(i),atlabm(i)
            end do
            do i = 1, nsb
               write(15,202) isb(i),jsb(i),ibo(i)
            end do

c***********************************************************************
c     charge line
c***********************************************************************

            nchg = 0
            do i = 1, nm

c carboxylates
               if( (itypem(i).eq.147).or.(itypem(i).eq.247).or.
     &          (itypem(i).eq.47) ) then
                  do j = 1, n12m(i)
                     spot = i12m(j,i)
                     if(itypem(spot).eq.3) then
                        if(bo12m(j,i).eq.1) then
                           nchg = nchg + 1
                           zatom(nchg) = i
                           zchg(nchg) = -1
                        end if
                     end if
                  end do

c diketonates
               elseif( (itypem(i).eq.168).or.(itypem(i).eq.272) ) then
                  do j = 1, n12m(i)
                     spot = i12m(j,i)
                     if( (itypem(spot).eq.2) .or.
     &                   (itypem(spot).eq.202)) then
                        if(bo12m(j,i).eq.1) then
                           nchg = nchg + 1
                           zatom(nchg) = i
                           zchg(nchg) = -1
                        end if
                     end if
                  end do

c alkoxides
               elseif( (itypem(i).eq.165).or.(itypem(i).eq.166).or.
     &          (itypem(i).eq.270).or.(itypem(i).eq.271) ) then
                  nchg = nchg + 1
                  zatom(nchg) = i
                  zchg(nchg) = -1

c sulfonates
               elseif(itypem(i).eq.169) then
                  do j = 1, n12m(i)
                     spot = i12m(j,i)
                     if(itypem(spot).eq.18) then  
                        if(bo12m(j,i).eq.1) then
                           nchg = nchg + 1
                           zatom(nchg) = i
                           zchg(nchg) = -1
                        end if
                     end if
                  end do

c ammonium
               elseif(itypem(i).eq.39)  then
                   nchg = nchg + 1
                   zatom(nchg) = i
                   zchg(nchg) = 1

c metal ion
               elseif(itypem(i).ge.300) then
                  nchg = nchg + 1
                  zatom(nchg) = i
                  zchg(nchg) = qm(ii)
               end if
            end do

            if(nchg.le.10) then
               write(15,203) nchg,(zatom(i),zchg(i),i=1,nchg)
            else
               write(6,'(a)') 'Error too many charged atoms'
               stop
            end if

c***********************************************************************
c     bottom of mol file
c***********************************************************************

            write(15,'(a)') 'M  END'
            write(15,'(a)') '$$$$'

  200 FORMAT(2I3,'  0  0  0  0  0  0  0  0999 V2000')
  201 FORMAT(3F10.4,x,a2,'  0  0  0  0  0  0')
  202 FORMAT(3I3)
  203 FORMAT('M  CHG',i3,20i4)

         close(15)

c***********************************************************************
c     bottom of loop over complexes
c***********************************************************************

  500 continue

      end do

c***********************************************************************
c     bottom of loop over ligand subdirectories
c***********************************************************************

         call system('rm conpcm')
         call system('rm pcmod.out')
         call system('rm pcmod.bak')
         call system('rm pcmod.par')
         call system('rm infile.pcm')
         
         cmd = 'mv input '//subdir(nn)(1:length)
         call system(cmd)
         cmd = 'mv *pcm '//subdir(nn)(1:length)
         call system(cmd)
         cmd = 'mv *dat '//subdir(nn)(1:length)
         call system(cmd)
         cmd = 'mv *xyz '//subdir(nn)(1:length)
         call system(cmd)

      end do 

c***********************************************************************
c     bottom of code
c***********************************************************************

      goto 1000

  999 continue
      write(6,'(a)') 'Error reading file named input'

 1000 continue

      end
