      subroutine atomtyper(natom,itype,atlab)

      implicit none
      include 'sizes.i'

      integer i,natom,itype(MAXAT)
      character*2 atlab(MAXAT)

      do i = 1, natom
         atlab(i) = '  '
         if( (itype(i).eq.1) .or. (itype(i).eq.2) .or.
     &   (itype(i).eq.3) .or. (itype(i).eq.4) .or.
     &   (itype(i).eq.22) .or. (itype(i).eq.29) .or.
     &   (itype(i).eq.30) .or. (itype(i).eq.38) .or.
     &   (itype(i).eq.50) .or. (itype(i).eq.56) .or.
     &   (itype(i).eq.57) .or. (itype(i).eq.58) .or.
     &   (itype(i).eq.67) .or. (itype(i).eq.68) .or.
     &   (itype(i).eq.71) .or. (itype(i).eq.106) .or.
     &   (itype(i).eq.114) .or. (itype(i).eq.156) .or.
     &   (itype(i).eq.157) .or. (itype(i).eq.158) .or.
     &   (itype(i).eq.160) .or. (itype(i).eq.160) .or.
     &   (itype(i).eq.162) .or. (itype(i).eq.202) ) then
            atlab(i) = 'C '
         elseif ( (itype(i).eq.5) .or. (itype(i).eq.21) .or.
     &   (itype(i).eq.23) .or. (itype(i).eq.24) .or.
     &   (itype(i).eq.28) .or. (itype(i).eq.44) .or.
     &   (itype(i).eq.48) .or. (itype(i).eq.73) .or.
     &   (itype(i).eq.124) ) then
            atlab(i) = 'H '
         elseif ( (itype(i).eq.6) .or. (itype(i).eq.7) .or.
     &   (itype(i).eq.41) .or. (itype(i).eq.47) .or.
     &   (itype(i).eq.49) .or. (itype(i).eq.69) .or.
     &   (itype(i).eq.70) .or. (itype(i).eq.145) .or.
     &   ((itype(i).ge.75).and.(itype(i).le.103)) .or.
     &   ((itype(i).ge.115).and.(itype(i).le.121)) .or.
     &   (itype(i).eq.147) .or. (itype(i).eq.148) .or.
     &   (itype(i).eq.149) .or. (itype(i).eq.159) .or.
     &   ((itype(i).ge.163).and.(itype(i).le.187)) .or.
     &   ((itype(i).ge.241).and.(itype(i).le.278)) .or.
     &   (itype(i).eq.205) .or. (itype(i).eq.206) .or.
     &   (itype(i).eq.207) ) then
            atlab(i) = 'O '
         elseif ( (itype(i).eq.8) .or. (itype(i).eq.9) .or.
     &   (itype(i).eq.10) .or. (itype(i).eq.37) .or.
     &   (itype(i).eq.39) .or. (itype(i).eq.40) .or.
     &   (itype(i).eq.43) .or. (itype(i).eq.45) .or.
     &   (itype(i).eq.46) .or. (itype(i).eq.72) .or.
     &   ((itype(i).ge.107).and.(itype(i).le.111)) .or.
     &   (itype(i).eq.143) .or. (itype(i).eq.144) .or.
     &   (itype(i).eq.146) .or. (itype(i).eq.150) .or.
     &   (itype(i).eq.151) .or. (itype(i).eq.155) .or.
     &   ((itype(i).ge.208).and.(itype(i).le.239)) ) then
            atlab(i) = 'N '
         elseif (itype(i).eq.11) then
            atlab(i) = 'F '
         elseif (itype(i).eq.12) then
            atlab(i) = 'Cl'
         elseif (itype(i).eq.13) then
            atlab(i) = 'Br'
         elseif (itype(i).eq.14) then
            atlab(i) = 'I '
         elseif (itype(i).eq.19) then
            atlab(i) = 'Si'
         elseif ( (itype(i).eq.25) .or. (itype(i).eq.60) .or.
     &   (itype(i).eq.153) ) then
            atlab(i) = 'P '
         elseif ( ((itype(i).ge.15).and.(itype(i).le.18)) .or.
     &   (itype(i).eq.42) .or. (itype(i).eq.74) .or.
     &   (itype(i).eq.104) .or. (itype(i).eq.105) .or.
     &   (itype(i).eq.154) ) then
            atlab(i) = 'S '
         elseif ( (itype(i).eq.26) .or. (itype(i).eq.27) ) then
            atlab(i) = 'B '
         elseif (itype(i).eq.200) then
            atlab(i) = 'Du'
         elseif (itype(i).eq.201) then
            atlab(i) = 'Lp'
         elseif (itype(i).gt.300) then
            atlab(i) = 'Ni'
         else
            write(6,'(a)') 'ERROR in atomtyper.f'
            write(6,'(a,i4)') ' Do not recognize atom type ',
     &       itype(i)
            goto 1000
         end if
         if(atlab(i).eq.'  ') then
            write(6,'(a)') 'ERROR in atomtyper.f'
            write(6,'(a,i4)') 'Failed to assign label to atom # ',i
            goto 1000
         end if
      end do

      return
 1000 continue
      stop
      end
