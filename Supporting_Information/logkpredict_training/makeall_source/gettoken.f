c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'gettoken' -- reads a variable from a string
c      adapted from Kevin Gilbert's routine used to read PCModel files
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine gettoken(state,head,nocaps,string,token)      
      
      implicit none                                               

      integer i, ic, foundtoken

      integer state, head                                             
      logical nocaps                                                            
      character*1 string(200),token(30)                                         
                           
      state=0                                                                   
      foundtoken=0                                                              
      do 10 i=1,30                                                              
        token(i)=' '                                                            
   10 continue                                                                  
      do 20 i=head,200                                                          
        if(string(i).ne.':'.and.string(i).ne.' '.and.string(i).ne.';'           
     *    .and.string(i).ne.','.and.string(i).ne.'    ') then                   
          ic=ichar(string(i))                                                   
          if(foundtoken.eq.0) then                                              
            foundtoken=i                                                        
            if(ic.ge.48.and.ic.le.57) then                                      
c             *** beginning to read a number                                    
              state=1                                                           
            elseif(string(i).eq.'.'.or.string(i).eq.'-') then                   
c             *** period or dash can't yet tell if word or number               
              state=3                                                           
            else                                                                
c             *** must be a word                                                
              state=2                                                           
            endif                                                               
          else                                                                  
            if(((state.eq.2).and.(ic.eq.45.or.ic.eq.46.or.ic.ge.48.and.         
     *        ic.le.57)).or.((state.eq.1).and.(ic.ne.46.and.(ic.lt.48)          
     *        .or.(ic.gt.57)))) then                                            
              foundtoken=i                                                      
              goto 30                                                           
            elseif(state.eq.3) then                                             
              if(ic.ge.48.and.ic.le.57.or.ic.eq.46) then                        
                state=1                                                         
              else                                                              
                state=2                                                         
              endif                                                             
            endif                                                               
          endif                                                                 
          token(i-foundtoken+1)=string(i)                                       
        elseif(foundtoken.ne.0) then                                            
          foundtoken=i                                                          
          goto 30                                                               
        endif                                                                   
   20 continue                                                                  
      if(foundtoken.eq.0) then                                                  
        head=1000                                                               
        return                                                                  
      endif                                                                     
   30 head=foundtoken                                                           
      if(state.eq.2.and..not.nocaps) then                                       
        do i=1,10                                                               
          if(lge(token(i),'a').and.lle(token(i),'z')) then                      
            token(i)=char(ichar(token(i))-32)                                   
          endif                                                                 
        enddo                                                                   
      endif                                                                     
      return                                                                    
      end                                                                       
