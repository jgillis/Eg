!
! Copyright 2003-2007 Henk Krus, Cyclone Fluid Dynamics BV
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
subroutine dolfyn2gmv

   use constants
   use geometry
   use variables
   
   integer :: indx(8)
   
   write(*,*)'Writing GMV data file...'
   
   call openfile(IOpst,casename,'.gmv','FORMATTED', &
                                       'SEQUENTIAL','UNKNOWN',debug)

   write(IOpst,'(A)')      'gmvinput ascii'
   write(IOpst,'(A)')      'comments'
   write(IOpst,'(A,f5.3)') ' Produced by Dolfyn vs ',float(version)*0.001
   write(IOpst,'(A)')      ' Copyright(C) 2003-2004 Cyclone Fluid Dynamics BV'
   write(IOpst,'(A)')      'endcomm'

   write(IOpst,'(A,i8)')  'nodev ',Nvrt
   do i=1,Nvrt
     write(IOpst,'(3(1x,f9.4))') Vert(i,:)
   end do
   write(*,*) 'GMV nodes written'
   !
   ! cells
   !
   !write(IOpst,'(A,i8,A)') 'cells ',Ncel+Nbnd
   write(IOpst,'(A,i8,A)') 'cells ',Ncel 

   call openfile(IOcel,casename,'.cel','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)

   NHex  = 0
   NPen  = 0
   NPyr  = 0
   NTet  = 0
   NQuad = 0
   NTri  = 0
   
   do i=1,Ncel
     read(IOcel,*) idummy,(indx(j),j=1,8)
     i1 = indx(1)
     i2 = indx(2)
     i3 = indx(3)
     i4 = indx(4)
     i5 = indx(5)
     i6 = indx(6)
     i7 = indx(7)
     i8 = indx(8)

     if( i7 /= 0 )then
       if( i7 /= i8 )then
         Nhex   = Nhex + 1
         write(IOpst,'(A)') 'hex 8'
         write(IOpst,'(8(1x,i8))') indx(1),indx(2),indx(3),indx(4), &
                                   indx(5),indx(6),indx(7),indx(8)
       else if( i5 /= i6 )then
         NPen = NPen + 1
         write(IOpst,'(A)') 'prism 6'
         write(IOpst,'(8(1x,i8))') indx(1),indx(2),indx(3), &
                                   indx(5),indx(6),indx(7)
       else if( i3 /= i4 )then
         NPyr  = NPyr + 1
         write(IOpst,'(A)') 'pyramid 5'
         write(IOpst,'(8(1x,i8))') indx(1),indx(2),indx(3),indx(4), &
                                   indx(5) 
       else if( i3 == i4 .and. i5 == i6 )then
         NTet   = NTet + 1
         write(IOpst,'(A)') 'tet 4'
         write(IOpst,'(8(1x,i8))') indx(1),indx(2),indx(3),indx(5)
       else
         write(*,*)'Unknown GMV shape, cell ',i,' : ',i1,i2,i3,i4,' ...'
       endif
     else
       !
       ! Reserved for baffles 
       !
       write(*,*)'Error: Unsuppported feature'
     
     endif     
   end do
   
   close(IOcel)

   do ib=1,Nbnd
     indx(1:4) = bnd(ib)%vertices
     if( indx(3) /= indx(4) )then
       NQuad   = NQuad + 1
   !    write(IOpst,'(A)') 'quad 4'
   !    write(IOpst,'(8(1x,i8))') indx(1),indx(2),indx(3),indx(4)
     elseif( i3 /= 0 )then
       NTri    = NTri + 1
   !    write(IOpst,'(A)') 'tri 3'
   !    write(IOpst,'(8(1x,i8))') indx(1),indx(2),indx(3) 
     else
       write(*,*)'*** Error: internal GMV boundary cell format (1)'
     endif
   end do

   
   if( NHex  > 0 )write(*,*)'Number of GMV hexa''s     :',NHex 
   if( NPen  > 0 )write(*,*)'Number of GMV prism''s    :',NPen 
   if( NPyr  > 0 )write(*,*)'Number of GMV pyramids''s :',NPyr 
   if( NTet  > 0 )write(*,*)'Number of GMV tetra''s    :',NTet 
   if( NQuad > 0 )write(*,*)'Number of GMV quad''s     :',NQuad
   if( NTri  > 0 )write(*,*)'Number of GMV triangles''s:',NTri 

   close(IOcel)
   write(*,*)'GMV cells written'

 !  write(*,'(1x,A)') 'cellids'
 !  write(IOpst,'(A)') 'cellids'
 !  do i=1,Ncel
 !    write(IOpst,'(8(1x,i4))') 1
 !  end do

 !  do i=1,Nbnd
 !    write(IOpst,'(8(1x,i4))') bnd(i)%rid+2
 !  end do

   write(*,'(1x,A)') 'velocities'
   write(IOpst,'(A)') 'velocity 0'
   do i=1,Ncel 
     write(IOpst,'(3(1x,1pe10.3))') U(i),V(i),W(i)
   end do

 !  do i=1,Nbnd 
 !    write(IOpst,'(3(1x,1pe10.3))') U(Ncel+i),V(Ncel+i),W(Ncel+i)
 !  end do

   write(IOpst,'(A)') 'variable'
   
   write(*,'(1x,A)') 'pressure'
   write(IOpst,'(A)') 'pressure 0'
   do i=1,Ncel 
     write(IOpst,'(3(1x,1pe10.3))') p(i)
   end do

 !  do i=1,Nbnd 
 !    write(IOpst,'(3(1x,1pe10.3))') p(Ncel+i)
 !  end do

   write(IOpst,'(A)') 'endvars'


   write(IOpst,'(A)') 'endgmv'

   close(IOpst)
   write(*,*)'Done'
   
end subroutine dolfyn2gmv


