!
! Copyright 2004-2007 Henk Krus, Cyclone Fluid Dynamics BV
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
! This is vtk.f90, the vtk interface to paraview, visit or mayavi
!
! Particles added by Johan Jacobs, Simex-Technology, summer 2007 
!
subroutine dolfyn2vtk
   use constants
   use geometry
   use variables
   use particles
   
   character(len=12) :: string1, string2
   integer :: indx(8)

   integer            :: datum(8)
   character (len=12) :: clock(3)
   character (len=10) :: datestring 

   character (len=96) :: CaseNameTransient = ' '

   real, dimension(3)   :: Xp, X, ds 
   real, dimension(3)   :: Xpn, Normal 
   real, dimension(10)  :: tmp              ! temp array 
   
   real, allocatable    :: NodalData1(:)    ! array of data at vertices
   real, allocatable    :: NodalData2(:)    
   real, allocatable    :: NodalData3(:)    
   
   integer, allocatable :: NodalCounter(:)  

   integer, allocatable :: VTKcells(:)  

   integer, save        :: ICounter = 0     ! counter for transient data


   call date_and_time(clock(1),clock(2),clock(3),datum)
   write(datestring,'(i2.2,''/'',i2.2,''/'',i4)') datum(3),datum(2),datum(1)

   allocate( VTKcells(Ncel+Nbnd),stat=istat)    
   allocate( NodalData1(Nvrt),stat=istat)    
   allocate( NodalCounter(Nvrt),stat=istat)    
   call TrackMemory(istat,Ncel+Nbnd+2*Nvrt,'Work arrays for NodalData allocated')
   
   write(*,*)        'Writing VTK data file...'
   write(IOdbg,*)    'Writing VTK data file...'
   write(IOdbg,'(A,A)')  ' Case:  ',casename(1:lens(casename)) 
   write(IOdbg,'(A,A)')  ' Title: ',title(1:lens(title)) 
   write(IOdbg,'(A,A)')  ' Date : ',datestring
   
   if( Transient )then
     n = lens(CaseName)
     ICounter = ICounter + 1
     if( ICounter > 9999 )then 
       write(*,*) '+++ Warning: Transient VTK output data file counter reset'
       ICounter = 1
     endif
     write(string1,'(i4.4)') ICounter
     CaseNameTransient( 1 : n ) = CaseName(1:n)
     CaseNameTransient(n+1:n+1) = '_'
     CaseNameTransient(n+2:n+5) = string1

     write(IOdbg,'(A,A)')  ' File:  ',CaseNameTransient(1:n+5)//'.vtk'

     call openfile(IOpst,CaseNameTransient,'.vtk','FORMATTED', &
                                         'SEQUENTIAL','UNKNOWN',debug)
   
   else
     call openfile(IOpst,casename,'.vtk','FORMATTED', &
                                         'SEQUENTIAL','UNKNOWN',debug)
   endif
   
   write(IOpst,'(A)')      '# vtk DataFile Version 3.0'
   write(IOpst,'(A,A,A)')   casename(1:lens(casename)),': ',title(1:lens(title)) 
   write(IOpst,'(A)')      'ASCII'

   !
   ! algemene info
   !
   write(IOpst,'(A)')      'DATASET UNSTRUCTURED_GRID'
     
   !
   ! vertices
   !
   write(IOpst,'(A,i8,A)') 'POINTS ',Nvrt,' float'
   do i=1,Nvrt
     write(IOpst,'(3(1x,1pe12.5))') Vert(i,:)
   end do
   !
   ! the VTK file format needs to know in advance what is coming
   !
   call openfile(IOcel,casename,'.cel','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)

   !  first walk through the cells establish the type
   !
   !    6 hexaeder:  1,2,3,4 5,6,7,8  
   !    5 pentaeder: 1,2,3,3 5,6,7,7  
   !    5 pyramid:   1,2,3,4 5,5,5,5  
   !    4 tetraeder: 1,2,3,3 4,4,4,4  
   !    1 quad:      1,2,3,4          
   !    1 triangle:  1,2,3,3          
   !
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
         VTKcells(i) = 12
       else if( i5 /= i6 )then
         NPen = NPen + 1
         VTKcells(i) = 13
       else if( i3 /= i4 )then
         NPyr  = NPyr + 1
         VTKcells(i) = 14
       else if( i3 == i4 .and. i5 == i6 )then
         NTet   = NTet + 1
         VTKcells(i) = 10
       else
         write(*,*)'Unknown VTK shape, cell ',i,' : ',i1,i2,i3,i4,' ...'
       endif
     else
       !
       ! Reserved for baffles 
       !
       write(*,*)'Error: Unsuppported feature'
              
       if( i4 /= -1 )then
         NQuad   = NQuad + 1
         VTKcells(i) = 9
       else 
         NTri    = NTri + 1
         VTKcells(i) = 5
       endif
     endif
   end do
   
   close(IOcel)

   do ib=1,Nbnd
     indx(1:4) = bnd(ib)%vertices
     !write(*,*)'>>',indx(1:4) 
     if( indx(4) /= -1 )then
       NQuad   = NQuad + 1
       VTKcells(Ncel+ib) = 9    
     elseif( indx(4) == -1 )then
       NTri    = NTri + 1
       VTKcells(Ncel+ib) = 5    
     else
       write(*,*)'*** Error: internal VTK boundary cell format (1)'
     endif
   end do

   
   if( NHex  > 0 )write(IOdbg,*)'Number of VTK hexa''s     :',NHex 
   if( NPen  > 0 )write(IOdbg,*)'Number of VTK prism''s    :',NPen 
   if( NPyr  > 0 )write(IOdbg,*)'Number of VTK pyramids''s :',NPyr 
   if( NTet  > 0 )write(IOdbg,*)'Number of VTK tetra''s    :',NTet 
   if( NQuad > 0 )write(IOdbg,*)'Number of VTK quad''s     :',NQuad
   if( NTri  > 0 )write(IOdbg,*)'Number of VTK triangles''s:',NTri 

   NVTKc = NHex * 9  + NPen * 7 + NPyr * 6 + NTet * 5
   NVTKb = NQuad * 5 + NTri * 4
   !
   ! cells
   !
   write(IOpst,'(A,i8,i8)') 'CELLS ',Ncel+Nbnd,NVTKc + NVTKb

   call openfile(IOcel,casename,'.cel','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)
   rewind(IOcel)
   
 1 format(i2,8(1x,i4))
 2 format(i2,8(1x,i5))   
 3 format(i2,8(1x,i6))   
 4 format(i2,8(1x,i7))  
 5 format(i2,8(1x,i8))  
 
   if( Nvrt <= 99999999 ) assign 5 to ifmt 
   if( Nvrt <=  9999999 ) assign 4 to ifmt 
   if( Nvrt <=   999999 ) assign 3 to ifmt 
   if( Nvrt <=    99999 ) assign 2 to ifmt 
   if( Nvrt <=     9999 ) assign 1 to ifmt 
      
   do i=1,Ncel
     read(IOcel,*) idummy,(indx(j),j=1,8)

     if( VTKcells(i) == 12 )then
       write(IOpst,ifmt) 8, &
                                    indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1, &
                                    indx(5)-1,indx(6)-1,indx(7)-1,indx(8)-1
     elseif( VTKcells(i) == 13 )then
       write(IOpst,ifmt) 6, &
                                    indx(1)-1,indx(2)-1,indx(3)-1, &
                                    indx(5)-1,indx(6)-1,indx(7)-1
     elseif( VTKcells(i) == 14 )then
       write(IOpst,ifmt) 5, &
                                    indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1, &
                                    indx(5)-1 
     elseif( VTKcells(i) == 10 )then
       write(IOpst,ifmt) 4, &
                                    indx(1)-1,indx(2)-1,indx(3)-1, &
                                    indx(5)-1 
     else
       write(*,*)'*** Error: internal VTK cell format'
     endif
   end do
   !
   ! next loop over boundaries
   !
   do ib=1,Nbnd
     indx(1:4) = bnd(ib)%vertices
     if( indx(4) /= -1 )then
       write(IOpst,ifmt) 4,indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1
     else 
       write(IOpst,ifmt) 3,indx(1)-1,indx(2)-1,indx(3)-1
     endif
   end do

   close(IOcel)
   if( Debug > 2 ) write(*,*)'Cells done'
   
   write(IOpst,'(A,i8,i8)') 'CELL_TYPES ',Ncel+Nbnd

   do i=1,Ncel,10
     k = min(i+9,Ncel)
     write(IOpst,'(i2,9(1x,i2))') (VTKcells(j),j=i,k)
   end do
   do ib=1,Nbnd,10
     k = min(ib+9,Nbnd)
     write(IOpst,'(i2,9(1x,i2))') (VTKcells(Ncel+j),j=ib,k)
   end do
   !
   ! cell id's
   !
   write(IOdbg,*)'Writing VTK cell ids'
   
   write(IOpst,'(A,i8)')  'CELL_DATA ',Ncel+Nbnd
   write(IOpst,'(A)')     'SCALARS type_id int'
   write(IOpst,'(A)')     'LOOKUP_TABLE default'

   IBoffset = maxval(cell(:)%ctid)
   !if( Debug > -1 ) write(*,*)'IBoffset:',IBoffset

   if( IBoffset <= 9 )then
     do i=1,Ncel,10
       k = min(i+9,Ncel)
       write(IOpst,'(10(i2))') (cell(j)%ctid,j=i,k)
     end do
   else if( IBoffset <= 99 )then
     do i=1,Ncel,10
       k = min(i+9,Ncel)
       write(IOpst,'(10(i3))') (cell(j)%ctid,j=i,k)
     end do
   else
     do i=1,Ncel,10
       k = min(i+9,Ncel)
       write(IOpst,'(10(i4))') (cell(j)%ctid,j=i,k)
     end do
   endif

   !if( Debug > -1 ) write(*,*)'IBoffset+Nreg:',IBoffset+Nreg
   if( IBoffset+Nreg+1 <= 9 )then
     do ib=1,Nbnd,10
       k = min(ib+9,Nbnd)
       write(IOpst,'(10(i2))') (bnd(j)%rid+(IBoffset+1),j=ib,k)
     end do
   else if( IBoffset+Nreg+1 <= 99 )then
     do ib=1,Nbnd,10
       k = min(ib+9,Nbnd)
       write(IOpst,'(10(i3))') (bnd(j)%rid+(IBoffset+1),j=ib,k)
     end do
   else
     do ib=1,Nbnd,10
       k = min(ib+9,Nbnd)
       write(IOpst,'(10(i4))') (bnd(j)%rid+(IBoffset+1),j=ib,k)
     end do
   endif

   !
   ! data
   !
   write(IOdbg,*)'Writing VTK vectors'
   write(IOpst,'(A)') 'VECTORS velocity float'

   do i=1,Ncel
     write(IOpst,'(3(1x,1pe10.3))') U(i),V(i),W(i)
   end do
   do i=1,Nbnd
     write(IOpst,'(3(1x,1pe10.3))') U(Ncel+i),V(Ncel+i),W(Ncel+i)
   end do
   
   !
   ! SCALARS CELL DATA
   !
   write(IOdbg,*)'Writing VTK magnitude'
   write(IOpst,'(A)') 'SCALARS velmag float 1'
   write(IOpst,'(A)') 'LOOKUP_TABLE default'
   
   do i=1,Ncel,10
     k = min(i+9,Ncel)
     write(IOpst,'(10(1x,1pe10.3))') &
       (sqrt( U(j)**2 + V(j)**2 + W(j)**2 ),j=i,k)
   end do

   do i=1,Nbnd,10
     k = min(i+9,Nbnd)
     write(IOpst,'(10(1x,1pe10.3))') &
       (sqrt( U(Ncel+j)**2 + V(Ncel+j)**2 + W(Ncel+j)**2 ),j=i,k)
   end do

   if( PostC(VarP) )then
     write(IOdbg,*)'Writing VTK pressure'
     write(IOpst,'(A)') 'SCALARS pressure float 1'
     write(IOpst,'(A)') 'LOOKUP_TABLE default'
     
     do i=1,Ncel,10
       k = min(i+9,Ncel)
       write(IOpst,'(10(1x,1pe10.3))') (p(j),j=i,k)
     end do

     do i=1,Nbnd,10
       k = min(i+9,Nbnd)
       write(IOpst,'(10(1x,1pe10.3))') (p(Ncel+j),j=i,k)
     end do
   endif

   if( SolveTurb )then
     if( PostC(VarTE) )then
       write(IOdbg,*)'Writing VTK turbulent kinetic energy'
       write(IOpst,'(A)') 'SCALARS k float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (max(0.0,TE(j)),j=i,k)
       end do

       !do ib=1,Nbnd
       !  ir = bnd(ib)%rid 
       !  if( Reg(ir)%typ == RWall )then
       !    TE(Ncel+ib) = 0.0
       !  endif
       ! end do

       do i=1,Nbnd,10
         k = min(i+9,Nbnd)
         write(IOpst,'(10(1x,1pe10.3))') (max(0.0,TE(Ncel+j)),j=i,k)
       end do
     endif

     if( PostC(VarED) )then
       write(IOdbg,*)'Writing VTK turbulent dissipation'
       write(IOpst,'(A)') 'SCALARS eps float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (max(0.0,ED(j)),j=i,k)
       end do
   
       !do ib=1,Nbnd
       !  ir = bnd(ib)%rid 
       !  if( Reg(ir)%typ == RWall )then
       !    ED(Ncel+ib) = 0.0
       !  endif
       !end do
       
       do i=1,Nbnd,10
         k = min(i+9,Nbnd)
         write(IOpst,'(10(1x,1pe10.3))') (max(0.0,ED(Ncel+j)),j=i,k)
       end do
     endif

   endif
   
   if( SolveEnthalpy )then
     if( PostC(VarT) )then
       write(IOdbg,*)'Writing VTK temperature'
       write(IOpst,'(A)') 'SCALARS temperature float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (T(j)-Tref,j=i,k)
       end do

       do i=1,Nbnd,10
         k = min(i+9,Nbnd)
         write(IOpst,'(10(1x,1pe10.3))') (T(Ncel+j)-Tref,j=i,k)
       end do
     endif 
   endif

   if( PostC(VarVIS) )then
     write(IOdbg,*)'Writing VTK effective viscosity'
     write(IOpst,'(A)') 'SCALARS viscosity float 1'
     write(IOpst,'(A)') 'LOOKUP_TABLE default'

     do i=1,Ncel,10
       k = min(i+9,Ncel)
       write(IOpst,'(10(1x,1pe10.3))') (VisEff(j),j=i,k)
     end do

     do i=1,Nbnd,10
       k = min(i+9,Nbnd)
       write(IOpst,'(10(1x,1pe10.3))') (VisEff(j),j=i,k)
     end do
   endif 

   if( SolveScalars )then
     if( PostC(VarSC) )then
         do is=1,NScal
           write(IOdbg,*)'Writing VTK scalar(s) ',Variable(NVar+is)
           write(IOpst,'(A,A,A)') 'SCALARS ',&
             Variable(NVar+is)(1:lens(Variable(NVar+is))),' float 1'
           write(IOpst,'(A)') 'LOOKUP_TABLE default'

           do i=1,Ncel,10
             k = min(i+9,Ncel)
             write(IOpst,'(10(1x,1pe10.3))') (max(0.0,SC(j,is)),j=i,k)
           end do

           do i=1,Nbnd,10
             k = min(i+9,Nbnd)
             write(IOpst,'(10(1x,1pe10.3))') (max(0.0,SC(j,is)),j=i,k)
           end do
         
         end do
     endif
   endif

   !
   ! SCALARS VERTEX DATA
   !
   if( PostV(VarU) .or. PostV(VarV) .or. PostV(VarW) .or. & 
       PostV(VarP) .or. PostV(VarT) .or.                  &
       PostV(VarSC) )                                     &
       write(IOpst,'(A,i8)') 'POINT_DATA ',Nvrt

   if( PostV(VarU) .or. PostV(VarV) .or. PostV(VarW) )then

     allocate( NodalData2(Nvrt),stat=istat)    
     allocate( NodalData3(Nvrt),stat=istat)    

     write(IOdbg,*)'Writing VTK data on nodes'
     call TrackMemory(istat,2*Nvrt,'Work arrays for extra NodalData allocated')

     call InterpolateData(0,VarU,U,dUdX,NodalData1,NodalCounter)
     call InterpolateData(0,VarV,V,dVdX,NodalData2,NodalCounter)
     call InterpolateData(0,VarW,W,dWdX,NodalData3,NodalCounter)
 
     !
     ! dump vectors on vertices
     !
     write(IOdbg,*)'Writing VTK vectors on nodes ',Nvrt
     write(IOpst,'(A)')    'VECTORS velocity float'

     do i=1,Nvrt
       if( abs(Nodaldata1(i)) < Small ) Nodaldata1(i) = 0.0
       if( abs(Nodaldata2(i)) < Small ) Nodaldata2(i) = 0.0
       if( abs(Nodaldata3(i)) < Small ) Nodaldata3(i) = 0.0

       write(IOpst,'(3(1x,1pe10.3))') Nodaldata1(i),NodalData2(i),NodalData3(i)
     end do

     !
     ! dump velmagnitude on vertices
     !
     write(IOdbg,*)'Writing VTK velocity magnitude on vertices'
     write(IOpst,'(A)')    'SCALARS velmag float 1'
     write(IOpst,'(A)')    'LOOKUP_TABLE default'

     do i=1,Nvrt
       if( abs(Nodaldata1(i)) < Small ) Nodaldata1(i) = 0.0
       if( abs(Nodaldata2(i)) < Small ) Nodaldata2(i) = 0.0
       if( abs(Nodaldata3(i)) < Small ) Nodaldata3(i) = 0.0

       tmpmag = sqrt( Nodaldata1(i)**2 + NodalData2(i)**2 + NodalData3(i)**2 )

       write(IOpst,'(1x,1pe10.3)') tmpmag
     end do

   endif

   if( PostV(VarP) )then   
     write(IOdbg,*)'Writing VTK pressure on vertices'
     write(IOpst,'(A)')    'SCALARS pressure float 1'
     write(IOpst,'(A)')    'LOOKUP_TABLE default'

     call InterpolateData(0,VarP,P,dPdX,NodalData1,NodalCounter)
 
     do i=1,Nvrt,10
       k = min(i+9,Nvrt)
       write(IOpst,'(10(1x,1pe10.3))') (Nodaldata1(j),j=i,k)
     end do

   endif

   if( PostV(VarT) )then

     write(IOdbg,*)'Writing VTK temperature on vertices'
     write(IOpst,'(A)')    'SCALARS temperature float 1'
     write(IOpst,'(A)')    'LOOKUP_TABLE default'

     call InterpolateData(0,VarT,T,dPdX,NodalData1,NodalCounter)

     do i=1,Nvrt,10
       k = min(i+9,Nvrt)
       write(IOpst,'(10(1x,1pe11.4))')  (Nodaldata1(j)-Tref,j=i,k)
     end do

   endif

   if( PostV(VarTE) )then

     write(IOdbg,*)'Writing VTK turbulent kinetic energy on vertices'
     write(IOpst,'(A)')    'SCALARS k float 1'
     write(IOpst,'(A)')    'LOOKUP_TABLE default'

     call InterpolateData(0,VarTE,TE,dPdX,NodalData1,NodalCounter)

     do i=1,Nvrt,10
       k = min(i+9,Nvrt)
       write(IOpst,'(10(1x,1pe11.4))')  (max(0.0,Nodaldata1(j)),j=i,k)
     end do

   endif

   if( PostV(VarED) )then

     write(IOdbg,*)'Writing VTK turbulent dissipation on vertices'
     write(IOpst,'(A)')    'SCALARS eps float 1'
     write(IOpst,'(A)')    'LOOKUP_TABLE default'

     call InterpolateData(0,VarED,ED,dPdX,NodalData1,NodalCounter)

     do i=1,Nvrt,10
       k = min(i+9,Nvrt)
       write(IOpst,'(10(1x,1pe11.4))')  (max(0.0,Nodaldata1(j)),j=i,k)
     end do

   endif

   if( PostV(VarVIS) )then

     write(IOdbg,*)'Writing VTK effective viscosity on vertices'
     write(IOpst,'(A)')    'SCALARS effvis float 1'
     write(IOpst,'(A)')    'LOOKUP_TABLE default'

     call InterpolateData(0,VarVIS,VISEFF,dPdX,NodalData1,NodalCounter)

     do i=1,Nvrt,10
       k = min(i+9,Nvrt)
       write(IOpst,'(10(1x,1pe11.4))')  (max(0.0,Nodaldata1(j)),j=i,k)
     end do

   endif

   if( SolveScalars )then
     if( PostV(VarSC) )then
         do is=1,NScal
         
           call InterpolateData(0,VarS(is),SC(:,is),dPdX,NodalData1,NodalCounter)

           write(IOdbg,*)'Writing VTK scalar(s) ',&
                                     Variable(NVar+is),' on vertices'
           write(IOpst,'(A,A,A)') 'SCALARS ',&
             Variable(NVar+is)(1:lens(Variable(NVar+is))),' float 1'
           write(IOpst,'(A)') 'LOOKUP_TABLE default'

           do i=1,NVrt,10
             k = min(i+9,Nvrt)
             write(IOpst,'(10(1x,1pe10.3))') (max(0.0,Nodaldata1(j)),j=i,k)
           end do
        
         end do
     endif
   endif

   if( allocated( NodalData1   ) ) deallocate ( NodalData1   )
   if( allocated( NodalData2   ) ) deallocate ( NodalData2   )
   if( allocated( NodalData3   ) ) deallocate ( NodalData3   )
   if( allocated( NodalCounter ) ) deallocate ( NodalCounter )
   if( allocated( VTKcells     ) ) deallocate ( VTKcells     )

   !
   ! do not forget to close, the unit number IOpst might be reused
   !
   close(IOpst)
   
   !
   ! VTK particle track plot, added by Johan Jacobs, Simex-Technology
   !
   if( UseParticles )then
     if( Npart >= 1000 )then 
       write(*,*) '+++ Warning: Number of particles must be lower then 1000'
       write(*,*) 'File tracks.vtk not written'
     else
       call openfile(IOpst,'tracks','.vtk','FORMATTED', &
                                      'SEQUENTIAL','UNKNOWN',debug)
       write(IOpst,'(A)')      '# vtk DataFile Version 3.0'
       write(IOpst,'(A,A,A)')   casename(1:lens(casename)),': ',title(1:lens(title)) 
       write(IOpst,'(A)')      'ASCII'
       write(IOpst,'(A)')      'DATASET POLYDATA'

       n_part_points = 0
       do i=1,Npart
         if( .not. associated( Tracks(i)%head ) )then
           write(*,*)'ParticlePrint: error head not associated',i
           cycle
         endif
         n_part_points = n_part_points + Tracks(i)%n      
       end do

       write(IOpst,'(A,i8,A)') 'POINTS ',n_part_points,' float'
       do i=1,Npart
         Track => Tracks(i)%head
         do while( associated(Track) )
           write(IOpst,'(3(1x,1pe12.5))') Track%x
           Track => Track%next
         end do
       end do
       
       write(IOpst,'(A,i8,i8)') 'LINES ',Npart,n_part_points+Npart
       i_part_point=0
       do i=1,Npart
         write(IOpst,'(i8)') Tracks(i)%n
         do j=1,Tracks(i)%n
           write(IOpst,'(i8)') i_part_point+j-1
         end do
         i_part_point=i_part_point+Tracks(i)%n   
       end do
       write(IOpst,'(A,i8)') 'POINT_DATA', n_part_points
       write(IOpst,'(A,i8)') 'SCALARS part_inr int'
       write(IOpst,'(A,i8)') 'LOOKUP_TABLE default'
       do i=1,Npart
         do j=1,Tracks(i)%n 
           write(IOpst,'(i8)') i
         end do
       end do
       close(IOpst)
     endif
   endif
   !
   ! end addition Johan Jacobs, Simex-Technology
   !
   
   if( Debug > 2 ) write(*,*)'Done'

end subroutine dolfyn2vtk


