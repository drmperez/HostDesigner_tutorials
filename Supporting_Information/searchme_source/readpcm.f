c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'readpcm' -- reads PCModel file containing a single structure 
c
c     passed variables:
c     char       filename    80 character string 
c     integer    ifield      mm3=3, mmff94=7
c     integer    n           number of atoms 
c     real       x,y,z       coordinates
c     integer    itype       atom type
c     integer    n12         number of attached atoms
c     integer    i12         serial numbers of attached atoms
c     integer    bo12        bond order to attached atoms
c     char       pichar      pi atom designators
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine readpcm(filename,ifield,n,x,y,z,itype,n12,i12,bo12,
     & pichar)
                                                        
      implicit none
      include 'params.i'

      integer i,j,k,ifield,ijunk
      integer n,itype(MAXATOM),n12(MAXATOM),i12(MAXVAL,MAXATOM)
      integer bo12(MAXVAL,MAXATOM)
      integer state,head

      real x(*),y(*),z(*),charge(MAXATOM),rjunk
    
      character*1 pichar(MAXATOM)

      character*30 token
      character*60 junk
      character*80 filename
      character*200 inputline

      logical nocaps
      
c-----------------------------------------------------------------------
c     Open PCModel formatted file
c-----------------------------------------------------------------------
        
      open(unit=25,file=filename,status='old',err=500)
      
      ifield = 0
      n = 0
      do i = 1, MAXATOM
         pichar(i) = ' '
         charge(i) = 0.0e0
         itype(i) = 0
         x(i) = 0.0e0
         y(i) = 0.0e0
         z(i) = 0.0e0
         n12(MAXATOM) = 0
         do j = 1, MAXVAL
            i12(j,i) = 0
            bo12(j,i) = 0
         end do
      end do
      
c-----------------------------------------------------------------------
c     Read header
c-----------------------------------------------------------------------

      read(25,'(a)',end=200) inputline 
      if(inputline(1:4).eq.'{PCM') read(inputline,'(5x,a)') junk

   20 read(25,'(a)',end=200) inputline
      head=1
      call gettoken(state,head,nocaps,inputline,token)
      if(head.eq.1000) goto 20          

c-----------------------------------------------------------------------
c     Read energy line
c-----------------------------------------------------------------------
   
      if(token.eq.'OPT') then
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
   26    if(token.eq.'PARA') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) ijunk
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'CONV') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) ijunk
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'IE') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) rjunk
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'FE') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) rjunk
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'IH') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) rjunk
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'FH') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) rjunk
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         end if
         goto 20

c-----------------------------------------------------------------------
c     Read number of atoms, n
c-----------------------------------------------------------------------

      else if(token.eq.'NA') then
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20 
         read(token,*) n                                                
         call gettoken(state,head,nocaps,inputline,token) 
         if(head.eq.1000) goto 20 
         goto 20 

c-----------------------------------------------------------------------
c     Read force field model
c-----------------------------------------------------------------------

      else if(token.eq.'ATOMTYPES') then
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20 
         read(token,*) ifield                                                
         call gettoken(state,head,nocaps,inputline,token) 
         if(head.eq.1000) goto 20 
         goto 20 

c-----------------------------------------------------------------------
c     Read atom types
c-----------------------------------------------------------------------
 
      else if(token.eq.'AT') then
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
         read(token,*) i 
         call gettoken(state,head,nocaps,inputline,token) 
         if(head.eq.1000) goto 20
         read(token,'(i5)') itype(i)
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
c     End of structure
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
