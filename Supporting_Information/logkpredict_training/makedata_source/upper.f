c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'upper' -- forces string to be upper case
c
c     passed variables:
c     character  s1  input string
c     character  s2  output string
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine upper(s1,s2)
      
      implicit none

      integer            :: i
      integer, PARAMETER :: DUC = ichar('A') - ichar('a')
      character(*)       :: s1
      character(len(s1)) :: s2
      character          :: ch

      do i = 1, len(s1)
         ch = s1(i:i)
         if (ch >= 'a'.and.ch <= 'z') ch = char(ichar(ch)+DUC)
         s2(i:i) = ch
      enddo

      return
      end
