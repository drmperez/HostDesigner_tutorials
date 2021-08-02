c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                     HostDesigner, Version 4
c       General purpose, structure-based molecular design code
c             Copyright (C) 2015-2020 by Benjamin P Hay
c 
c     This program is free software: you can redistribute it and/or
c     modify it under the terms of the GNU General Public License as 
c     published by the Free Software Foundation, either version 3 of the
c     License or (at your option) any later version.
c
c     This program is distributed in the hope that it will be useful, 
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See GNU 
c     General Public License for more details.  You should have received
c     a copy of the GNU General Public License along with this program.
c     If not, see <http://www.gnu.org/licenses/>.
c
c     Benjamin P. Hay
c     Supramolecular Design Institute
c     127 Chestnut Hill Road
c     Oak Ridge, TN 37830
c
c     hayben@comcast.net
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     'gettoken' -- reads variable from string
c      adapted from Kevin Gilbert's routine used to read PCModel files
c
c     passed variables:
c     integer     state
c     integer     head
c     logical     nocaps
c     character   string
c     character   token
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine gettoken(state,head,nocaps,string,token)      
      
      implicit none 

      integer i,ic,foundtoken
      integer state,head                                             
      logical nocaps
      character*1 string(200),token(30)
                           
      state = 0 
      foundtoken = 0 
      do i = 1, 30 
        token(i) = ' '
      end do

      do i = head, 200
         if(string(i).ne.':'.and.string(i).ne.' '.and.string(i).ne.';'
     &      .and.string(i).ne.','.and.string(i).ne.'    ') then 
            ic = ichar(string(i))
            if(foundtoken.eq.0) then
               foundtoken = i
               if(ic.ge.48.and.ic.le.57) then 
                  state = 1                             ! must be number
               else if(string(i).eq.'.'.or.string(i).eq.'-') then 
                  state = 3                    ! might be number or word
               else 
                  state = 2                               ! must be word
               end if
            else
               if(((state.eq.2).and.(ic.eq.45.or.ic.eq.46.or.ic.ge.48
     &            .and.ic.le.57)).or.((state.eq.1).and.(ic.ne.46.and.
     &            (ic.lt.48).or.(ic.gt.57)))) then
                  foundtoken=i
                  goto 30
               else if(state.eq.3) then
                   if(ic.ge.48.and.ic.le.57.or.ic.eq.46) then
                      state=1 
                   else 
                      state=2
                   end if
               end if
            end if
            token(i-foundtoken+1) = string(i)
         else if(foundtoken.ne.0) then
            foundtoken=i
            goto 30
         end if
      end do

      if(foundtoken.eq.0) then
         head = 1000
         return
      end if

   30 head = foundtoken 
      if(state.eq.2.and..not.nocaps) then
         do i = 1, 10
            if(lge(token(i),'a').and.lle(token(i),'z')) then
               token(i) = char(ichar(token(i))-32)
            end if
         end do
      end if

      return
      end
