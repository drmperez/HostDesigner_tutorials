c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'readpcm' -- reads a PCModel file containing a single 
c     structure within it, stores the data into the stomm common blocks
c
c     passed variables:
c     char       file        up to 40 char name 
c     integer    n           number of atoms 
c     integer    type        ff atom type
c     real       x,y,z       coordinates
c     real       e1,e2       steric energies
c     integer    n12         number of attached atoms
c     integer    i12         serial numbers of attached atoms
c     integer    bo12        bond order to attached atoms
c     integer    miss        missing parameters
c     integer    fail        calculation did not complete
c     char       pichar      pi atom designators
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine readpcm(filename,n,type,x,y,z,e1,e2,n12,i12,bo12,
     &   miss,fail,pichar)  
                                                        
      implicit none
      include 'sizes.i'

      integer i,j,k,miss,fail
      integer n,type(*),n12(*),i12(MAXVAL,*),bo12(MAXVAL,*)
      integer state,head

      real x(*),y(*),z(*),charge(MAXAT)
      real e1,e2,h1,h2
    
      character*1 pichar(*)
      character*30 token
      character*40 filename
      character*60 junk
      character*200 inputline

      logical nocaps
      
c-----------------------------------------------------------------------
c     open PCModel formatted file
c-----------------------------------------------------------------------

      open(unit=35,file=filename,status='old',err=500)
      
      do i = 1, MAXAT
         pichar(i) = ' '
         charge(i) = 0.0
         n12(MAXAT) = 0
         do j = 1, MAXVAL
            i12(j,i) = 0
            bo12(j,i) = 0
         end do
      end do

      miss = 0
      fail = 0
      e1 = 0.0e0
      e2 = 0.0e0
      h1 = 0.0e0
      h2 = 0.0e0
      
c-----------------------------------------------------------------------
c     read the header of the file
c-----------------------------------------------------------------------

      read(35,'(a)',end=200) inputline 
      if(inputline(1:4).eq.'{PCM') read(inputline,'(5x,a)') junk

   20 read(35,'(a)',end=200) inputline
      head=1
      call gettoken(state,head,nocaps,inputline,token)
      if(head.eq.1000) goto 20          

c-----------------------------------------------------------------------
c     read the energy line and assign relative energies
c-----------------------------------------------------------------------
   
      if(token.eq.'OPT') then
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
   26    if(token.eq.'PARA') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) miss
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'CONV') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) fail
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'IE') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) e1
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'FE') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) e2
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'IH') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) h1
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'FH') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) h2
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         end if
         goto 20

c-----------------------------------------------------------------------
c     read the number of atoms - n
c-----------------------------------------------------------------------

      else if(token.eq.'NA') then
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20 
         read(token,*) n                                                
         call gettoken(state,head,nocaps,inputline,token) 
         if(head.eq.1000) goto 20 
         goto 20 

c-----------------------------------------------------------------------
c     read the atom types
c-----------------------------------------------------------------------
 
      else if(token.eq.'AT') then
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
         read(token,*) i 
         call gettoken(state,head,nocaps,inputline,token) 
         if(head.eq.1000) goto 20
         read(token,'(i5)') type(i)
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20  
         read(token,*) x(i) 
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20 
         read(token,*) y(i)
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
         read(token,*) z(i)
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
   80    if(token.eq.'B') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20 
            j = 0
   90       j = j + 1 
            if(state.eq.1) then 
               read(token,*) i12(j,i) 
               n12(i) = j                                         
            else
              goto 80
            end if  
            call gettoken(state,head,nocaps,inputline,token) 
            if(head.eq.1000) goto 20
            if(state.eq.1) then
               read(token,*) bo12(j,i) 
            else 
               write(*,*) 'Subroutine readpcm - missing bondorder' 
               stop 
            end if 
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20 
            goto 90     
         else if(token.eq.'P') then 
            pichar(i) = 'P' 
            call gettoken(state,head,nocaps,inputline,token) 
            if(head.eq.1000) goto 20
            goto 80 
         else if(token.eq.'C') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) charge(i)
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 80        
         else if(token.eq.'S') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) k
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 80
         else if(token.eq.'N') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) k
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 80
         else if(token.eq.'D') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) k
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 80
         else if(token.eq.'R') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) k
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 80
         end if 

c-----------------------------------------------------------------------
c     end of structure
c-----------------------------------------------------------------------

      else if(token.eq.'}') then 
        goto 200          
      else   
        goto 20
      end if
      
  200 close(25)   
      return   
  500 continue
      write(*,'(a,a)') 'Subroutine readpcm - no file named ',
     & filename
      stop 
      end 
