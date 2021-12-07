c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'writepcm' -- writes structure in PCModel format. 
c
c     passed variables:
c     char     filename  name of file to be written
c     char     title     title line
c     integer  n         number of atoms
c     real     x,y,z     coordinates
c     integer  type      force field atom type
c     integer  n12       number of attached atoms
c     integer  i12       serial numbers of attached atoms
c     integer  bo12      bondorder array
c     integer  pichar    pi atom list
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine writepcm(filename,title,n,x,y,z,type,n12,i12,bo12,
     &   pichar)
      
      implicit none
      include 'sizes.i'
            
      integer i,j,length
      integer n,type(*),n12(*),i12(MAXVAL,*),bo12(MAXVAL,*)
      
      real x(*),y(*),z(*)
      
      character*1 pichar(*)
      character*40 filename
      character*50 title
      character*150 aline
      
      open(unit=77,file=filename,status='unknown')

      write(77,'(a4,a50)') '{PCM',title
      
      write(77,'(a2,i5)') 'NA', n

      write(77,'(a9,i2)') 'ATOMTYPES', 3
              
      do i = 1, n
      
         write(aline,10) 'AT', i, type(i), x(i), y(i), z(i),
     &   'B', (i12(j,i),',',bo12(j,i),j=1,n12(i))    
          
            length = 150
            do while (aline(length:length).eq.' ')
               length=length-1
            end do

            if(pichar(i).eq.'P') then
               length=length+2
               aline(length:length)='P'
            end if
            
            write(77,'(a)') aline(1:length)
            
         end do
         
      write(77,*) '}'
         
   10 FORMAT(a2,1x,i3,1x,i3,3(1X,F9.5),1x,a1,1x,20(i5,a1,i2,1x))
   
      close(77)
            
      return
      end
