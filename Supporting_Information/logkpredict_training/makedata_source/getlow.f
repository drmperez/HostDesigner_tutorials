c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'getlow' -- reads search file created by PCModel, passes back
c     coordinates for lowest energy structure and its energy
c
c     passed variables:
c     char     outfile   name of the output file
c     real     x,y,z     coordinates
c     real     elow      energy of lowest energy conformer
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine getlow(outfile,x,y,z,elow)  
                                                        
      implicit none
      include 'sizes.i'

      integer i,j,k,state,head

      real a,x(MAXAT),y(MAXAT),z(MAXAT),elow
    
      character*1 char1
      character*30 token
      character*60 junk
      character*40 outfile
      character*200 inputline

      logical nocaps
      
c-----------------------------------------------------------------------
c     Open PCModel formatted file
c-----------------------------------------------------------------------
        
      elow = 0.00000e0

      open(unit=25,file=outfile,status='old',err=500)
      
c-----------------------------------------------------------------------
c     Start reading next hit
c-----------------------------------------------------------------------
            
      read(25,'(a)',end=200) inputline                                          
      if(inputline(1:4).eq.'{PCM') then
         read(inputline,'(5x,a)') junk
      else
         write(6,'(a)') 'Do not find {PCM marker, input problem'
         stop
      end if

   20 read(25,'(a)',end=200) inputline
      head=1
      call gettoken(state,head,nocaps,inputline,token)
      if(head.eq.1000) goto 20                                

c-----------------------------------------------------------------------
c     Read energy line and assign relative energies
c-----------------------------------------------------------------------

      if(token.eq.'OPT') then
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
   26    if(token.eq.'PARA') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) k
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'CONV') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) k
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'IE') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) a
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'FE') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) elow
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'IH') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) a
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         else if(token.eq.'FH') then         
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) a
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 26
         end if
         goto 20
         
c-----------------------------------------------------------------------
c     Read number of atoms
c-----------------------------------------------------------------------

      else if(token.eq.'NA') then
         call gettoken(state,head,nocaps,inputline,token) 
         if(head.eq.1000) goto 20 
         read(token,'(i5)') k
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
         goto 20                
                                                         
c-----------------------------------------------------------------------
c     Read atoms
c-----------------------------------------------------------------------
  
      else if(token.eq.'AT') then 
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
         read(token,*) i 
         call gettoken(state,head,nocaps,inputline,token)
         if(head.eq.1000) goto 20
         read(token,*) k
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
               read(token,*) k
            else
               goto 80
            end if
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            if(state.eq.1) then
               read(token,*) k
            else
               write(*,*) 'Subroutine readsearch - missing bond order'
               stop
            end if 
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 90
         else if(token.eq.'P') then
            char1 = 'P'
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            goto 80
         else if(token.eq.'C') then
            call gettoken(state,head,nocaps,inputline,token)
            if(head.eq.1000) goto 20
            read(token,*) a
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
      if(elow.eq.0.00000e0) then
         write(*,'(a)') 'ELOW = 0.0000e0, possible problem'
         stop
      endif
      return   
                                                           
  500 continue
      write(*,'(a,a)') 'Subroutine readsearch - no file named ',
     &   outfile
      close(25)
      stop                                                                   
      end                                                                       
