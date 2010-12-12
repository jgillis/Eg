!
! Copyright 2006-2007 Henk Krus, Cyclone Fluid Dynamics BV
! All Rights Reserved.
!
! Licensed under the Apache License, Version 2.0 (the "License"); 
! you may not use this file except in compliance with the License. 
! You may obtain a copy of the License at
!
! http://www.dolfyn.net/license.html
! 
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an 
! "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, 
! either express or implied. See the License for the specific 
! language governing permissions and limitations under the License.
!
! With contributions by Shibo (Harry) Kuang
!
subroutine CollectCells(level,ip,nc,CellList)
!
! collect the neighbors of each cell
!
! for now a quick and dirty subroutine; not very
! suitable for big models but useful for now.
!
!
   use constants
   use geometry
   use variables

   integer, dimension(100) :: CellList                 ! <= HARDWIRED
   integer, dimension(100) :: VertList                 ! <= HARDWIRED

   logical, save :: Done = .false.

   type CellNeighborData
     integer ic
     type(CellNeighborData), pointer :: next => NULL()
   end type
   
   type CellListData
     type(CellNeighborData), pointer :: head => NULL()
     type(CellNeighborData), pointer :: tail => NULL()
     integer                         :: n
   end type CellListData

   type(CellListData), allocatable, save :: CellListArray(:)

   type(CellNeighborData), pointer   :: CellPtr, CellEntry


   type VertNeighborData
     integer ic
     type(VertNeighborData), pointer :: next => NULL()
   end type
   
   type VertListData
     type(VertNeighborData), pointer :: head => NULL()
     type(VertNeighborData), pointer :: tail => NULL()
     integer                         :: n
   end type VertListData

   type(VertListData), allocatable   :: VertListArray(:)

   type(VertNeighborData), pointer   :: VertPtr, VertEntry

   logical :: Skip

   if( level == -1 )then
     if( allocated( CellListArray ) )then
       write(*,*)'CollectCells cleanup'

       do i=1,Ncel
         CellEntry => CellListArray(i)%head%next
         do while( associated(CellEntry) )
           CellPtr => CellEntry%next
           deallocate( CellEntry )
           CellEntry => CellPtr
         end do
         deallocate(CellListArray(i)%head)
         deallocate(CellListArray(i)%tail)
       end do

       deallocate(CellListArray)

       Done = .false.
     endif
     return
   endif

   write(IOdbg,*)'CollectCells'
   if( .not. Done)then
     !
     ! collect a list
     !
     write(*,*)'CollectCells allocating arrays'

     allocate(VertListArray(Nvrt))

     !
     ! now simple only one layer, later using 'level'
     ! also multiple layers around the cell. note
     ! for 3d cases it requires memory to grow with n**3.
     !
     write(*,*)'List of vertices with cells'
     VertList = 0
     do i=1,NCel
       if( mod(i,50000) == 0 ) write(*,*) ' ',i

       nv = 0       
       do j=1,Nfaces(i)
         k   = CFace(i,j)
         ivo = -1
         do l=1,4
           iv = Face(k)%vertices(l)
           if( iv == ivo ) cycle
           do m=1,nv
             if( iv == VertList(m) ) goto 1
           end do
           nv  = nv + 1
           VertList(nv) = iv
           ivo = iv
         1 continue   
         end do
       end do

       do j=1,nv
         iv = VertList(j)
         if( .not. associated(VertListArray(iv)%head) )then 
           allocate(VertListArray(iv)%head)
           allocate(VertListArray(iv)%tail)
           
           VertListArray(iv)%head%ic   = i
           VertListArray(iv)%n = 1
           VertListArray(iv)%tail%next => VertListArray(iv)%head
         else
         
           Skip = .false.
           VertPtr => VertListArray(iv)%head
           do m=1,VertListArray(iv)%n
             icold = VertPtr%ic
             if( icold == i )then
               Skip = .true.
               exit
             endif
             VertPtr => VertPtr%next  
           end do

           if( .not. Skip)then 
             VertPtr => VertListArray(iv)%tail%next

             allocate( VertEntry)
           
             VertListArray(iv)%n = VertListArray(iv)%n + 1
             VertEntry%ic = i
           
             VertPtr%next => VertEntry
           
             VertListArray(iv)%tail%next => VertEntry
           endif
         endif
       end do
     end do
     
     !do i=1,Nvrt
     !    write(88,'('' v '',i6,i6,'':'',$)') i,VertListArray(i)%n
     !    VertPtr => VertListArray(i)%head
     !    do j=1,VertListArray(i)%n
     !      write(88,'(i6,$)') VertPtr%ic
     !      VertPtr => VertPtr%next  
     !    end do
     !    write(88,*)
     !end do

     !
     ! now the relationship vertex -> cell is setup
     !
     ! loop over cells and faces again
     ! but relate cell to the cells of the vertices
     !
     allocate(CellListArray(Ncel))

     VertList = 0
     CellList = 0
     write(*,*)'List of cell neighbors'
     do i=1,NCel
       if( mod(i,50000) == 0 ) write(*,*) ' ',i

       nv = 0       
       do j=1,Nfaces(i)
         k   = CFace(i,j)
         ivo = -1
         do l=1,4
           iv = Face(k)%vertices(l)
           if( iv == ivo ) cycle
           if( iv ==  -1 ) cycle
           do m=1,nv
             if( iv == VertList(m) ) goto 2
           end do
           nv  = nv + 1
           if( nv > 100 )write(*,*)'*** F O U T!***'
           VertList(nv) = iv
           ivo = iv
         2 continue   
         end do
       end do

       do j=1,nv
         iv = VertList(j)

         VertPtr => VertListArray(iv)%head
         do m=1,VertListArray(iv)%n
           in = VertPtr%ic
           
           if( in /= i )then

             if( .not. associated(CellListArray(i)%head) )then 
               allocate(CellListArray(i)%head)
               allocate(CellListArray(i)%tail)

               CellListArray(i)%head%ic   = in
               CellListArray(i)%n = 1
               CellListArray(i)%tail%next => CellListArray(i)%head
             else

               Skip = .false.
               CellPtr => CellListArray(i)%head
               do n=1,CellListArray(i)%n
                 icold = CellPtr%ic
                 if( icold == in )then
                   Skip = .true.
                   exit
                 endif
                 CellPtr => CellPtr%next  
               end do

               if( .not. Skip)then 
                 CellPtr => CellListArray(i)%tail%next

                 allocate( CellEntry)

                 CellListArray(i)%n = CellListArray(i)%n + 1
                 CellEntry%ic = in

                 CellPtr%next => CellEntry

                 CellListArray(i)%tail%next => CellEntry
               endif
             endif

           endif
           VertPtr => VertPtr%next  
         end do
       end do

       !if( CellListArray(i)%n > 32 )then
       !  write(88,'('' c '',i8,i4,'':'',$)') i,CellListArray(i)%n
       !  CellPtr => CellListArray(i)%head
       !  do j=1,CellListArray(i)%n
       !    write(88,'(i8,$)') CellPtr%ic
       !    CellPtr => CellPtr%next  
       !  end do
       !  write(88,*)
       !endif
     end do

     !do i=1,Ncel
     !  write(*,'('' c '',i6,'':'',$)') i 
     !  CellPtr => CellListArray(i)%head
     !  do j=1,CellListArray(i)%n
     !    write(*,'(i6,$)') CellPtr%ic
     !    CellPtr => CellPtr%next  
     !  end do
     !  write(*,*)
     !end do

     !
     ! clean up; keep CellListArray
     !
     write(*,*)'Clean up'
     do i=1,Nvrt
       if( associated(VertListArray(i)%head) )then
         VertEntry => VertListArray(i)%head%next
         do while( associated(VertEntry) )
           VertPtr => VertEntry%next
           deallocate( VertEntry )
           VertEntry => VertPtr
         end do
         deallocate(VertListArray(i)%head)
         deallocate(VertListArray(i)%tail)
       endif
     end do

     deallocate(VertListArray)
     write(*,*)'Done clean up'

     Done = .true.
   endif
   
   !
   ! return result for cell ic
   !
   CellList = 0
   
   nc = CellListArray(ip)%n
   
   CellPtr => CellListArray(ip)%head
   do i=1,min(nc,100)                                      ! <= HARDWIRED
     CellList(i) = CellPtr%ic
     CellPtr => CellPtr%next  
   end do

   write(IOdbg,*)'CollectCells for cell',ip,'>',nc,':',CellList(1:nc)
   write(IOdbg,*)'CollectCells Done'

end subroutine CollectCells


