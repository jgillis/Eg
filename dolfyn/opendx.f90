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
subroutine dolfyn2opendx

   use constants
   use geometry
   use variables
   use scalars
   use particles
   
   character(len=12) :: string1, string2, string3, string4
   integer :: indx(8)

   integer            :: datum(8)
   character (len=12) :: clock(3)
   character (len=10) :: datestring 
   character (len= 9) :: scstring1, scstring2
   character (len=96) :: CaseNameTransient = ' '

   real, dimension(3)   :: Xp, X, ds, Xpn, Normal 
   real, allocatable    :: NodalData1(:)     ! array of data at vertices
   real, allocatable    :: NodalData2(:)    
   real, allocatable    :: NodalData3(:)    
   
   integer, allocatable :: NodalCounter(:)  

   integer, save        :: ICounter = 0     ! counter for transient data

   call date_and_time(clock(1),clock(2),clock(3),datum)
   write(datestring,'(i2.2,''/'',i2.2,''/'',i4)') datum(3),datum(2),datum(1)

   allocate( NodalData1(Nvrt),stat=istat)    
   allocate( NodalCounter(Nvrt),stat=istat)    
   call TrackMemory(istat,2*Nvrt,'Work arrays for NodalData allocated')
   
   write(*,*)            'Writing OpenDX data file...'
   write(IOdbg,*)        'Writing OpenDX data file...'
   write(IOdbg,'(A,A)')  ' Case:  ',casename(1:lens(casename)) 
   write(IOdbg,'(A,A)')  ' Title: ',title(1:lens(title)) 
   
   if( Transient )then 
     n = lens(CaseName)
     ICounter = ICounter + 1
     if( ICounter > 9999 )then 
       write(*,*) '+++ Warning: Transient OpenDX output data file counter reset'
       ICounter = 1
     endif
     write(string1,'(i4.4)') ICounter
     CaseNameTransient( 1 : n ) = CaseName(1:n)
     CaseNameTransient(n+1:n+1) = '_'
     CaseNameTransient(n+2:n+5) = string1

     write(IOdbg,'(A,A)')  ' File:  ',CaseNameTransient(1:n+5)//'.odx'

     call openfile(IOpst,CaseNameTransient,'.odx','FORMATTED', &
                                       'SEQUENTIAL','UNKNOWN',debug)
     
   else 
     call openfile(IOpst,casename,'.odx','FORMATTED', &
                                       'SEQUENTIAL','UNKNOWN',debug)
   endif

   write(IOpst,'(A)')      '#'
   write(IOpst,'(A)')      '# Dolfyn => OpenDX'
   write(IOpst,'(A,f5.3)') '# Produced by Dolfyn vs ',float(version)*0.001
   write(IOpst,'(A)')      '# Copyright(C) 2007 Cyclone Fluid Dynamics BV'
   write(IOpst,'(A)')      '#'
   write(IOpst,'(A,A)')    '# Case:  ',casename(1:lens(casename)) 
   write(IOpst,'(A,A)')    '# Title: ',title(1:lens(title)) 
   write(IOpst,'(A,A)')    '# Date : ',datestring
   if( Transient )then
     write(IOpst,'(A,i8,A,1pe9.3)') '# Step/Time:  ',Iter,'  ',Time
   else
     write(IOpst,'(A,i8)')     '# Step:  ',Iter
   endif
   write(IOpst,'(A)')      '#'

   !
   ! algemene info
   !
   write(IOpst,'(A,f5.3,A,A,A,A,A)') &
                 'object "dolfyn" class string "dolfyn ', &
                              float(version)*0.001,' ',datestring,'"'                      
   write(IOpst,'(A,A,A,A,A)') 'object "casename" class string "', &
     casename(1:lens(casename)),': ',title(1:lens(title)),'"'

   TmpTime = Time
   if( .not. Transient ) TmpTime = -1.0
   write(IOpst,'(A)') &
     'object "time" class array type float rank 0 items 1 data follows'
   write(IOpst,'(1pe13.6)')   TmpTime
     
   !
   ! vertices
   !
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') '# vertices'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A,i8,A)') &
     'object "vertices" class array type float rank 1 shape 3 items ',&
     Nvrt,' data follows'
   rscalefactor = 1.0  ! / ScaleFactor ! terug naar mm?
   do i=1,Nvrt
     write(IOpst,'(3(1x,1pe11.4))') Vert(i,:) * rscalefactor
   end do
   write(IOpst,'(A)') 'attribute "dep" string "positions"'
   !
   ! cells
   !
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') '# cells'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A,i8,A)') &
     'object "cells" class array type int rank 1 shape 8 items ',&
     Ncel,' data follows'

   call openfile(IOcel,casename,'.cel','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)

   do i=1,Ncel
     read(IOcel,*) idummy,(indx(j),j=1,8)
     write(IOpst,'(8(1x,i8))') indx(1)-1,indx(2)-1,indx(5)-1,indx(6)-1, &
                               indx(4)-1,indx(3)-1,indx(8)-1,indx(7)-1
   end do

   close(IOcel)
   if( Debug > 2 ) write(*,*)'Cells done'
   
   write(IOpst,'(A)') 'attribute "element type" string "cubes"'
   write(IOpst,'(A)') 'attribute "dep" string "connections"'

   !
   ! wall/boundary cells
   !
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') '# section with wall/boundary cells'
   write(IOpst,'(A)') '#'

   do ir=0,Nreg

     nwalls = count( Bnd(:)%rid == ir )

     if( nwalls == 0 )then
       if( Debug > 2 ) write(*,*)'Empty region:',ir
       write(IOpst,'(A,i8)') '# note: empty region ',ir
       cycle
     else
       if( Debug > 2 ) write(*,*)'Region ',ir,'(',nwalls,')'
     endif
     
     if( ir >= 99 ) STOP 'pas opendx.f90 aan'
     write(string1,'(''walls'',i2.2)') ir
     write(string2,'(''w'',i2.2)') ir
     
     write(IOpst,'(A,A,A,i8,A)') 'object "',string1(1:7), &
       '" class array type int rank 1 shape 4 items ',&
        nwalls,' data follows'

     icnt = 0
     do ib=1,Nbnd
       if( bnd(ib)%rid == ir )then
         indx(1:4) = bnd(ib)%vertices
         if( indx(4) /= -1 )then
           write(IOpst,'(8(1x,i8))') indx(1)-1,indx(4)-1,indx(2)-1,indx(3)-1
           icnt = icnt + 1
         else
           write(IOpst,'(8(1x,i8))') indx(1)-1,indx(3)-1,indx(2)-1,indx(3)-1
           icnt = icnt + 1
         endif
       endif
     end do
     write(IOpst,'(A)') 'attribute "element type" string "quads"'
     write(IOpst,'(A)') 'attribute "dep" string "connections"'

     write(IOpst,'(A,A,A)') 'object "',string2(1:3),'" class field'
     write(IOpst,'(A)')     'component "positions" "vertices"'
     write(IOpst,'(A,A,A)')  'component "connections" "',string1(1:7),'"'
     write(IOpst,'(A)')     '###'

     if( icnt /= nwalls )then
       write(*,*)'Warning walls mismatch'
     endif
   end do
   if( Debug > 2 ) write(*,*)'Walls done'
   

   if( Debug > 2 ) write(*,*)'Writing vectors'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') '# velocity vectors'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A,i8,A)') &
     'object "velocity" class array type float rank 1 shape 3 items ',&
     Ncel,' data follows'

   do i=1,Ncel
     if( abs(U(i)) < Small ) U(i) = 0.0
     if( abs(V(i)) < Small ) V(i) = 0.0
     if( abs(W(i)) < Small ) W(i) = 0.0
     
     write(IOpst,'(3(1x,1pe10.3))') U(i),V(i),W(i)
   end do
   write(IOpst,'(A)') 'attribute "dep" string "connections"'
   if( Debug > 2 ) write(*,*)'Vectors done'
   
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') '# velocity vector object (cell data)'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') 'object "V" class field'
   write(IOpst,'(A)') 'component "positions" "vertices"'
   write(IOpst,'(A)') 'component "connections" "cells"'
   write(IOpst,'(A)') 'component "data" "velocity"'

   if( PostV(VarU) .or. PostV(VarV) .or. PostV(VarW) )then

     allocate( NodalData2(Nvrt),stat=istat)    
     allocate( NodalData3(Nvrt),stat=istat)    
     call TrackMemory(istat,2*Nvrt,'Work arrays for extra NodalData allocated')
     !
     ! dump velocity vectors on the nodes as well
     !
     if( Debug > 2 ) write(*,*)'Writing vectors on nodes'
      
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# velocity vectors at nodes'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A,i8,A)') &
       'object "velocity2" class array type float rank 1 shape 3 items ',&
       Nvrt,' data follows'

     call InterpolateData(0,VarU,U,dUdX,NodalData1,NodalCounter)
     call InterpolateData(0,VarV,V,dVdX,NodalData2,NodalCounter)
     call InterpolateData(0,VarW,W,dWdX,NodalData3,NodalCounter)
 
     do i=1,Nvrt
       if( abs(Nodaldata1(i)) < Small ) Nodaldata1(i) = 0.0
       if( abs(Nodaldata2(i)) < Small ) Nodaldata2(i) = 0.0
       if( abs(Nodaldata3(i)) < Small ) Nodaldata3(i) = 0.0

       write(IOpst,'(3(1x,1pe10.3))') Nodaldata1(i),NodalData2(i),NodalData3(i)
     end do
     write(IOpst,'(A)') 'attribute "dep" string "positions"'

     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# velocity vector object (nodal data)'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') 'object "V2" class field'
     write(IOpst,'(A)') 'component "positions" "vertices"'
     write(IOpst,'(A)') 'component "connections" "cells"'
     write(IOpst,'(A)') 'component "data" "velocity2"'
   
     if( allocated( NodalData2   ) ) deallocate ( NodalData2   )
     if( allocated( NodalData3   ) ) deallocate ( NodalData3   )
   endif
   !
   ! SCALARS
   !
   if( Debug > 2 ) write(*,*)'Writing cell type id'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') '# cell type id'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A,i8,A)') &
     'object "ctid" class array type integer rank 0 items ',&
     Ncel,' data follows'

   do i=1,Ncel
     !if( cell(i)%ctid < 99 )then
     if( i < 99 )then
       write(IOpst,'(i2)') cell(i)%ctid
       !write(IOpst,'(i2)') i
     else
       write(IOpst,'(i2)') cell(i)%ctid
       !write(IOpst,'(i4)') i
     endif
   end do
   write(IOpst,'(A)') 'attribute "dep" string "connections"'

   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') '# cell type id object (cell data)'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') 'object "id" class field'
   write(IOpst,'(A)') 'component "positions" "vertices"'
   write(IOpst,'(A)') 'component "connections" "cells"'
   write(IOpst,'(A)') 'component "data" "ctid"'

   if( PostC(VarP) )then
     if( Debug > 2 ) write(*,*)'Writing pressure'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# pressure data'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A,i8,A)') &
       'object "pressure" class array type float rank 0 items ',&
       Ncel,' data follows'

     do i=1,Ncel
       if( abs(P(i)) < Small ) P(i) = 0.0

       write(IOpst,'(3(1x,1pe10.3))') p(i)
     end do
     write(IOpst,'(A)') 'attribute "dep" string "connections"'

     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# pressure object (cell data)'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') 'object "P" class field'
     write(IOpst,'(A)') 'component "positions" "vertices"'
     write(IOpst,'(A)') 'component "connections" "cells"'
     write(IOpst,'(A)') 'component "data" "pressure"'
   endif
   
   if( PostV(VarP) )then
   
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# pressure data at nodes'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A,i8,A)') &
       'object "pressure2" class array type float rank 0 items ',&
       Nvrt,' data follows'

     call InterpolateData(0,VarP,P,dPdX,NodalData1,NodalCounter)
 
     do i=1,Nvrt
       write(IOpst,'(1x,1pe11.4)') Nodaldata1(i)
     end do

     write(IOpst,'(A)') 'attribute "dep" string "positions"'

     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# pressure object (nodal data)'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') 'object "P2" class field'
     write(IOpst,'(A)') 'component "positions" "vertices"'
     write(IOpst,'(A)') 'component "connections" "cells"'
     write(IOpst,'(A)') 'component "data" "pressure2"'
   endif

   if( SolveEnthalpy )then
     if( PostC(VarT) )then
       if( Debug > 2 ) write(*,*)'Writing temperature'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# temperature data'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,i8,A)') &
         'object "temperature" class array type float rank 0 items ',&
         Ncel,' data follows'

       do i=1,Ncel
         write(IOpst,'(1x,1pe11.4))')  T(i)-Tref
       end do

       write(IOpst,'(A)') 'attribute "dep" string "connections"'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# temperature object (cell data)'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') 'object "T" class field'
       write(IOpst,'(A)') 'component "positions" "vertices"'
       write(IOpst,'(A)') 'component "connections" "cells"'
       write(IOpst,'(A)') 'component "data" "temperature"'
     endif 
   
     if( PostV(VarT) )then
       if( Debug > 2 ) write(*,*)'Writing nodal temperature ',Tref
   
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# temperature data at nodes'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,i8,A)') &
         'object "temperature2" class array type float rank 0 items ',&
         Nvrt,' data follows'

       call InterpolateData(0,VarT,T,dPdX,NodalData1,NodalCounter)

       do i=1,Nvrt
         write(IOpst,'(1x,1pe11.4)')  (Nodaldata1(i)-Tref)
       end do

       write(IOpst,'(A)') 'attribute "dep" string "positions"'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# temperature object (nodal data)'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') 'object "T2" class field'
       write(IOpst,'(A)') 'component "positions" "vertices"'
       write(IOpst,'(A)') 'component "connections" "cells"'
       write(IOpst,'(A)') 'component "data" "temperature2"'
     endif
   else
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# NO temperature data'
     write(IOpst,'(A)') '#'
   endif

   if( SolveTurb )then
     if( PostC(VarTE) )then
       if( Debug > 2 ) write(*,*)'Writing turbulent kinetic energy'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent kinetic energy data'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,i8,A)') &
         'object "turbke" class array type float rank 0 items ',&
         Ncel,' data follows'

       do i=1,Ncel
         write(IOpst,'(3(1x,1pe9.3))') max(0.0,TE(i))
       end do
       write(IOpst,'(A)') 'attribute "dep" string "connections"'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent kinetic energy object (cell data)'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') 'object "K" class field'
       write(IOpst,'(A)') 'component "positions" "vertices"'
       write(IOpst,'(A)') 'component "connections" "cells"'
       write(IOpst,'(A)') 'component "data" "turbke"'
     endif
     if( PostV(VarTE) )then
       if( Debug > 2 ) write(*,*)'Writing turbulent kinetic energy'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent kinetic energy (nodal) data'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,i8,A)') &
         'object "turbke2" class array type float rank 0 items ',&
         Nvrt,' data follows'

       call InterpolateData(0,VarTE,TE,dPdX,NodalData1,NodalCounter)

       do i=1,Nvrt
         write(IOpst,'(1x,1pe10.3)') max(0.0,Nodaldata1(i))
       end do

       write(IOpst,'(A)') 'attribute "dep" string "positions"'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent kinetic energy object (nodal data)'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') 'object "K2" class field'
       write(IOpst,'(A)') 'component "positions" "vertices"'
       write(IOpst,'(A)') 'component "connections" "cells"'
       write(IOpst,'(A)') 'component "data" "turbke2"'
     endif

     if( PostC(VarED) )then
       if( Debug > 2 ) write(*,*)'Writing turbulent dissipation'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent dissipation data'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,i8,A)') &
         'object "turbeps" class array type float rank 0 items ',&
         Ncel,' data follows'

       do i=1,Ncel
         write(IOpst,'(3(1x,1pe9.3))') max(0.0,ED(i))
       end do
       write(IOpst,'(A)') 'attribute "dep" string "connections"'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent dissipation object (cell data)'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') 'object "E" class field'
       write(IOpst,'(A)') 'component "positions" "vertices"'
       write(IOpst,'(A)') 'component "connections" "cells"'
       write(IOpst,'(A)') 'component "data" "turbeps"'
     endif

     if( PostV(VarED) )then
       if( Debug > 2 ) write(*,*)'Writing turbulent dissipation'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent dissipation (nodal) data'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,i8,A)') &
         'object "turbeps2" class array type float rank 0 items ',&
         Nvrt,' data follows'

       call InterpolateData(0,VarED,ED,dPdX,NodalData1,NodalCounter)
 
       do i=1,Nvrt
         write(IOpst,'(3(1x,1pe10.3))') max(0.0,Nodaldata1(i))
       end do
       write(IOpst,'(A)') 'attribute "dep" string "positions"'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent dissipation object (nodal data)'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') 'object "E2" class field'
       write(IOpst,'(A)') 'component "positions" "vertices"'
       write(IOpst,'(A)') 'component "connections" "cells"'
       write(IOpst,'(A)') 'component "data" "turbeps2"'
     endif

     if( PostC(VarVis) )then
       if( Debug > 2 ) write(*,*)'Writing turbulent viscosity'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent viscosity data'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,i8,A)') &
         'object "turbvis" class array type float rank 0 items ',&
         Ncel,' data follows'

       do i=1,Ncel
         write(IOpst,'(3(1x,1pe9.3))') max(0.0,VisEff(i) - VisLam)
       end do
       write(IOpst,'(A)') 'attribute "dep" string "connections"'

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent viscosity object (cell data)'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') 'object "VIS" class field'
       write(IOpst,'(A)') 'component "positions" "vertices"'
       write(IOpst,'(A)') 'component "connections" "cells"'
       write(IOpst,'(A)') 'component "data" "turbvis"'
     endif

     if( PostV(VarVis) )then
       if( Debug > 2 ) write(*,*)'Writing turbulent viscosity'

       call GradientPhi(VarVis,VisEff,dPdX)      

       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent viscosity nodal data'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,i8,A)') &
         'object "turbvis2" class array type float rank 0 items ',&
         Nvrt,' data follows'

       call InterpolateData(0,VarVis,VisEff,dPdX,NodalData1,NodalCounter)
 
       do i=1,Nvrt
         write(IOpst,'(3(1x,1pe10.3))') max(0.0,Nodaldata1(i) - VisLam)
       end do

       write(IOpst,'(A)') 'attribute "dep" string "positions"'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# turbulent viscosity object (cell data)'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') 'object "VIS2" class field'
       write(IOpst,'(A)') 'component "positions" "vertices"'
       write(IOpst,'(A)') 'component "connections" "cells"'
       write(IOpst,'(A)') 'component "data" "turbvis2"'
     endif
   else
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# NO turbulence data'
     write(IOpst,'(A)') '#'
   endif
   
   if( SolveScalars )then

     do is=1,NScal
       write(scstring1,'(''scalar'',i3.3)') is
      !write(scstring2,'(''SC'',i3.3)') is
       if( Debug > 1 ) write(*,*)'Writing scalar: ',Variable(NVar+is),scstring1,scstring2
       if( PostC(VarSC) )then
         if( Debug > 2 ) write(*,*)'Writing scalar(s)'
         write(IOpst,'(A)')   '#'
         write(IOpst,'(A,A)') '# scalar data ',Variable(NVar+is)
         write(IOpst,'(A)')   '#'
         write(IOpst,'(A,A,A,i8,A)') &
           'object "',scstring1,'" class array type float rank 0 items ',&
           Ncel,' data follows'

         do i=1,Ncel
           write(IOpst,'(1x,1pe9.3)') max(SC(i,is),0.0)            ! <=== !! SCALARS POSITIVE?
         end do
         
         write(IOpst,'(A)') 'attribute "dep" string "connections"'

         write(IOpst,'(A)') '#'
         write(IOpst,'(A)') '# scalar object (cell data)'
         write(IOpst,'(A)') '#'
         write(IOpst,'(A,A,A)') 'object "', &
           Variable(NVar+is)(1:lens(Variable(NVar+is))),'" class field'
         write(IOpst,'(A)') 'component "positions" "vertices"'
         write(IOpst,'(A)') 'component "connections" "cells"'
         write(IOpst,'(A,A,A)') 'component "data" "',scstring1,'"'
       endif
       if( PostV(VarSC) )then
         if( Debug > 2 ) write(*,*)'Writing scalar(s)'

         write(IOpst,'(A)') '#'
         write(IOpst,'(A,A)') '# scalar (nodal) data ',Variable(NVar+is)
         write(IOpst,'(A)') '#'
         write(IOpst,'(A,A,A,i8,A)') &
           'object "V',scstring1,'" class array type float rank 0 items ',&
           Nvrt,' data follows'

         call InterpolateData(0,VarS(is),SC(1,IS),dPdX,NodalData1,NodalCounter)
 
         do i=1,Nvrt
           write(IOpst,'(3(1x,1pe10.3))') max(0.0,Nodaldata1(i)) ! <=== !! SCALARS POSITIVE?
         end do
         write(IOpst,'(A)') 'attribute "dep" string "positions"'

         write(IOpst,'(A)') '#'
         write(IOpst,'(A)') '# scalar object (nodal data)'
         write(IOpst,'(A)') '#'
         write(IOpst,'(A,A,A)') 'object "', &
           Variable(NVar+is)(1:lens(Variable(NVar+is))),'2" class field'
         write(IOpst,'(A)') 'component "positions" "vertices"'
         write(IOpst,'(A)') 'component "connections" "cells"'
         write(IOpst,'(A,A,A)') 'component "data" "V',scstring1,'"'
       endif

     end do
   endif
   !
   ! wall data stuff
   !
   if( SolveEnthalpy )then
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# section with thermal wall data'
     write(IOpst,'(A)') '#'

     do ir=0,Nreg

       nwalls = count( Bnd(:)%rid == ir )

       if( Reg(ir)%typ == RWall )then

         if( nwalls == 0 )then
           if( Debug > 2 ) write(*,*)'Empty region:',ir
           write(IOpst,'(A,i8)') '# note: empty region ',ir
           cycle
         else
           if( Debug > -1 ) write(*,*)'Wall region ',ir,'(',nwalls,')'
         endif

         !
         ! heat transfer coefficient
         !
         write(string1,'(''wallh'',i2.2)') ir
         write(string2,'(''wh'',i2.2)') ir
         write(string3,'(''walls'',i2.2)') ir

         write(IOpst,'(A,A,A,i8,A)') 'object "',string1(1:7), &
           '" class array type float rank 0 items ',&
            nwalls,' data follows'

         do ib=1,Nbnd
           if( bnd(ib)%rid == ir ) write(IOpst,'(1x,1pe10.3)') bnd(ib)%h           
         end do
         
         write(IOpst,'(A)') 'attribute "dep" string "connections"'
         write(IOpst,'(A,A,A)') 'object "',string2(1:4),'" class field'
         write(IOpst,'(A)')     'component "positions" "vertices"'
         write(IOpst,'(A,A,A)') 'component "connections" "',string3(1:7),'"'
         write(IOpst,'(A,A,A)') 'component "data" "',string1(1:7),'"'
         write(IOpst,'(A)')     '###'
         !
         ! heat flux q
         !
         write(string1,'(''wallq'',i2.2)') ir
         write(string2,'(''wq'',i2.2)') ir
         write(string3,'(''walls'',i2.2)') ir

         write(IOpst,'(A,A,A,i8,A)') 'object "',string1(1:7), &
           '" class array type float rank 0 items ',&
            nwalls,' data follows'

         do ib=1,Nbnd
         
           if( bnd(ib)%rid == ir ) write(IOpst,'(1x,1pe10.3)') bnd(ib)%q
           
         end do
         write(IOpst,'(A)') 'attribute "dep" string "connections"'

         write(IOpst,'(A,A,A)') 'object "',string2(1:4),'" class field'
         write(IOpst,'(A)')     'component "positions" "vertices"'
         write(IOpst,'(A,A,A)') 'component "connections" "',string3(1:7),'"'
         write(IOpst,'(A,A,A)') 'component "data" "',string1(1:7),'"'
         write(IOpst,'(A)')     '###'
         !
         ! yplus
         !
         write(string1,'(''wally'',i2.2)') ir
         write(string2,'(''wy'',i2.2)') ir
         write(string3,'(''walls'',i2.2)') ir

         write(IOpst,'(A,A,A,i8,A)') 'object "',string1(1:7), &
           '" class array type float rank 0 items ',&
            nwalls,' data follows'

         do ib=1,Nbnd
         
           if( bnd(ib)%rid == ir ) write(IOpst,'(1x,1pe10.3)') bnd(ib)%yplus
           
         end do
         write(IOpst,'(A)') 'attribute "dep" string "connections"'

         write(IOpst,'(A,A,A)') 'object "',string2(1:4),'" class field'
         write(IOpst,'(A)')     'component "positions" "vertices"'
         write(IOpst,'(A,A,A)') 'component "connections" "',string3(1:7),'"'
         write(IOpst,'(A,A,A)') 'component "data" "',string1(1:7),'"'
         write(IOpst,'(A)')     '###'
         !
         ! temperature
         !
         write(string1,'(''wallt'',i2.2)') ir
         write(string2,'(''wt'',i2.2)') ir
         write(string3,'(''walls'',i2.2)') ir

         write(IOpst,'(A,A,A,i8,A)') 'object "',string1(1:7), &
           '" class array type float rank 0 items ',&
            nwalls,' data follows'

         do ib=1,Nbnd
           if( bnd(ib)%rid == ir )then
             if(  T(Ncel+ib) > Small )then
               write(IOpst,'(1x,1pe10.3)') T(Ncel+ib)-Tref
             else
               write(IOpst,'(1x,1pe10.3)') 0.0             
             endif
           endif
         end do
         write(IOpst,'(A)') 'attribute "dep" string "connections"'

         write(IOpst,'(A,A,A)') 'object "',string2(1:4),'" class field'
         write(IOpst,'(A)')     'component "positions" "vertices"'
         write(IOpst,'(A,A,A)') 'component "connections" "',string3(1:7),'"'
         write(IOpst,'(A,A,A)') 'component "data" "',string1(1:7),'"'
         write(IOpst,'(A)')     '###'



       endif
     end do





     if( Debug > -1 ) write(*,*)'Thermal wall data done'
   endif

   !
   ! PARTICLES
   !
   if( UseParticles )then
     if( Npart >= 1000 ) STOP 'pas opendx.f90 aan'

     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# particles'
     write(IOpst,'(A)') '#'
   
     rscalefactor = 1.0  ! / ScaleFactor ! terug naar mm?
     do i=1,Npart
       if( .not. associated( Tracks(i)%head ) )then
         write(*,*)'ParticlePrint: error head not associated',i
         cycle
       endif

       write(string1,'(''partxyz'',i4.4)') i
       write(string2,'(''track'',i4.4)')   i
       write(string3,'(''partcon'',i4.4)') i
       write(string4,'(''partdat'',i4.4)') i

       write(IOpst,'(A,A11,A,i8,A)') &
       'object "',string1,'" class array type float rank 1 shape 3 items ',&
       Tracks(i)%n,' data follows'                                     

       Track => Tracks(i)%head
       do while( associated(Track) )
         write(IOpst,'(3(1x,1pe12.5))') Track%x * rscalefactor
         Track => Track%next
       end do

       write(IOpst,'(A)') 'attribute "dep" string "positions"'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# define as a line'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,A11,A,i8)') 'object "',string3, &
                                   '" gridconnections counts ', Tracks(i)%n
       write(IOpst,'(A)') 'attribute "element type" string "lines"'
       write(IOpst,'(A)') 'attribute "dep" string "connections"'
       write(IOpst,'(A)') 'attribute "ref" string "positions"'
   
       write(IOpst,'(A)') '#'
       write(IOpst,'(A)') '# particle object'
       write(IOpst,'(A)') '#'
       write(IOpst,'(A,A11,A,i8)') 'object "',string4,&
                                  '" class regulararray count ',Tracks(i)%n
       write(IOpst,'(A,i8,A)') 'origin ',i,' delta 1'
       write(IOpst,'(A)') 'attribute "dep" string "positions"'
       write(IOpst,'(A)') '#'

       write(IOpst,'(A,A9,A)') 'object "',string2,'" class field'
       write(IOpst,'(A,A11,A)') 'component "positions" "',string1,'"'
       write(IOpst,'(A,A11,A)') 'component "connections" "',string3,'"'
       write(IOpst,'(A,A11,A)') 'component "data" "',string4,'"'



       write(IOpst,'(A)')    '#'
       write(IOpst,'(A,i8)') '# end particle object, particle ',i
       write(IOpst,'(A)')    '#'

     end do

     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# tracks'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') 'object "tracks" class group'
     icnt = 0
     do i=1,Npart
       if( .not. associated( Tracks(i)%head ) )then
         write(*,*)'ParticlePrint: error head not associated',i
         cycle
       endif
       write(string2,'(''track'',i4.4)') i
       write(IOpst,'(A,i4,A,A9,A)') 'member ',icnt,' value "',string2,'"'
       icnt = icnt + 1       
     end do
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# end tracks'
     write(IOpst,'(A)') '#'
   endif

   ! 
   ! =======================================================================
   !
   ! debug stuff
   !
   if( DXnormals )then
     if( Debug > 2 ) write(*,*)'Writing DXnormals'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# face normals'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A,i8,A)') &
       'object "normalnodes" class array type float rank 1 shape 3 items ',&
       Nfac*2,' data follows'
     do i=1,Nfac
       write(IOpst,'(3(1x,f9.4))') Face(i)%x,Face(i)%x+1.0*Face(i)%n
       
       ! vector Xpn
       !ip = Face(i)%cell1
       !in = Face(i)%cell2
       !if( in > 0 )then
       !  Xpn    = Cell(in)%x - Cell(ip)%x
       !  write(IOpst,'(3(1x,f9.4))') Cell(ip)%x,Cell(ip)%x+0.9*Xpn
       !else
       !  write(IOpst,'(3(1x,f9.4))') Cell(ip)%x,Cell(ip)%x
       !endif
     end do
     write(IOpst,'(A)') 'attribute "dep" string "positions"'

     write(IOpst,'(A,i8,A)') &
       'object "normallines" class array type int rank 1 shape 2 items ',&
       Nfac,' data follows'
     do i=1,Nfac*2,2
       write(IOpst,'(2(1x,i8))') i-1,i
     end do
     write(IOpst,'(A)') 'attribute "element type" string "lines"'
     write(IOpst,'(A)') 'attribute "dep" string "connections"'
     write(IOpst,'(A)') 'object "normals" class field'
     write(IOpst,'(A)') 'component "positions" value "normalnodes"'
     write(IOpst,'(A)') 'component "connections" value "normallines"'
   endif

   if( DXcenters )then
     if( Debug > 2 ) write(*,*)'Writing DXcenters'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# cell centers'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A,i8,A)') &
       'object "centernodes" class array type float rank 1 shape 3 items ',&
       Ncel,' data follows'
     do i=1,Ncel
       write(IOpst,'(3(1x,f9.4))') Cell(i)%x
     end do
     write(IOpst,'(A)') 'attribute "dep" string "positions"'
     write(IOpst,'(A)') 'object "centers" class field'
     write(IOpst,'(A)') 'component "positions" value "centernodes"'
   endif
     
   if( DXmassflux )then
     if( Debug > 2 ) write(*,*)'Writing DXmassflux'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# face mass fluxes'
     write(IOpst,'(A)') '#' 
     write(IOpst,'(A,i8,A)') &
       'object "massfluxnodes" class array type float rank 1 shape 3 items ',&
       Nfac*2,' data follows'
     do i=1,Nfac
       normal = face(i)%n
       call normalise(normal)
       write(IOpst,'(6(1x,1pe12.3))') Face(i)%x, &
           Face(i)%x + 0.5*MassFlux(i)*normal
     end do
     write(IOpst,'(A)') 'attribute "dep" string "positions"'
   
     write(IOpst,'(A,i8,A)') &
       'object "massfluxlines" class array type int rank 1 shape 2 items ',&
       Nfac,' data follows'
     do i=1,Nfac*2,2
       write(IOpst,'(2(1x,i8))') i-1,i
     end do
    write(IOpst,'(A)') 'attribute "element type" string "lines"'
    write(IOpst,'(A)') 'attribute "dep" string "connections"'

     write(IOpst,'(A,i8,A)') &
       'object "massfluxdata" class array type float rank 1 shape 3 items ',&
       Nfac,' data follows'
     do i=1,Nfac
       normal = face(i)%n
       call normalise(normal)
       write(IOpst,'(3(1x,1pe12.3))') MassFlux(i)*normal
     end do
     write(IOpst,'(A)') 'attribute "dep" string "connections"'

     write(IOpst,'(A)') 'object "massflux" class field'
     write(IOpst,'(A)') 'component "positions" value "massfluxnodes"'
     write(IOpst,'(A)') 'component "connections" value "massfluxlines"'
     write(IOpst,'(A)') 'component "data" value "massfluxdata"'
   endif

   if( DXdebugdata )then
     if( Debug > 2 ) write(*,*)'Writing DX debug cell data'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# cell debug data'
     write(IOpst,'(A)') '#' 

     write(IOpst,'(A,i8,A)') &
       'object "dxdebugdata" class array type float rank 0 items ',&
       Ncel,' data follows'
     do i=1,Ncel
       write(IOpst,'(3(1x,e13.7))') DXdebug(i)
 !      write(IOpst,'(3(1x,e13.7))') Cell(i)%Vol
     end do
     write(IOpst,'(A)') 'attribute "dep" string "connections"'

     write(IOpst,'(A)') 'object "debug" class field'
     write(IOpst,'(A)') 'component "positions" "vertices"'
     write(IOpst,'(A)') 'component "connections" "cells"'
     write(IOpst,'(A)') 'component "data" "dxdebugdata"'

     if( Debug > 2 ) write(*,*)'Writing DX debug gradient data'
     write(IOpst,'(A)') '#'
     write(IOpst,'(A)') '# cell debug gradient data'
     write(IOpst,'(A)') '#' 

     write(IOpst,'(A,i8,A)') &
       'object "dxdebugdata2" class array type float rank 1 shape 3 items ',&
       Ncel,' data follows'
     do i=1,Ncel
       write(IOpst,'(3(1x,e14.7))') DXgrad(i,:)
     end do
     write(IOpst,'(A)') 'attribute "dep" string "connections"'

     write(IOpst,'(A)') 'object "debug2" class field'
     write(IOpst,'(A)') 'component "positions" "vertices"'
     write(IOpst,'(A)') 'component "connections" "cells"'
     write(IOpst,'(A)') 'component "data" "dxdebugdata2"'
   endif

   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') '# data object in order to address ids'
   write(IOpst,'(A)') '#'
   write(IOpst,'(A)') 'object "DATA" class group'
   write(IOpst,'(A)') 'member  "id" value "id"'
   write(IOpst,'(A)') 'member  "V"  value "V"'
   if( PostC(VarP) )then
     write(IOpst,'(A)') 'member "P" value "P"'
   endif
   if( PostC(VarT) )then
     write(IOpst,'(A)') 'member "T" value "T"'
   endif
   if( SolveTurb .and. PostC(VarTE) )then
     write(IOpst,'(A)') 'member "K" value "K"'
   endif
   if( SolveTurb .and. PostC(VarED) )then
     write(IOpst,'(A)') 'member "E" value "E"'
   endif
   if( SolveTurb .and. PostC(VarVis))then
     write(IOpst,'(A)') 'member "VIS" value "VIS"'
   endif


   if( allocated( NodalData1   ) ) deallocate ( NodalData1    )
   if( allocated( NodalCounter ) ) deallocate ( NodalCounter )

   !
   ! do not forget to close, the unit number IOpst might be reused
   !
   close(IOpst)
   if( Debug > 2 ) write(*,*)'Done'
   
end subroutine dolfyn2opendx
subroutine InterpolateData(Mode,IVar,Phi,dPhidX,NodalData,NodalCounter)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX   
   real, dimension(Nvrt)        :: NodalData
 
   integer, dimension(Nvrt)     :: NodalCounter  

   real, dimension(3)   :: Xp, X, ds 

   if( Debug > 2 ) write(*,*)'*** InterpolateData',Variable(IVar)

   call GradientPhi(IVar,Phi,dPhidX)     

   NodalData    = 0.0
   NodalCounter =  0
   
   do ip=1,Ncel
     Xp   = Cell(ip)%x
     tmpp = Phi(ip)
     do j=1,NFaces(ip)
       k  = CFace(ip,j)

       if( Face(k)%vertices(4) > 0 )then                
         do i1=1,4                                      
           i2  = Face(k)%vertices(i1)                   
           X   = Vert(i2,:)                             
           ds  = X - Xp                                 
           tmp = tmpp + dot_product( dPhidX(ip,:), ds )   
           NodalData(i2)  = NodalData(i2)  + tmp        
           NodalCounter(i2) = NodalCounter(i2) + 1      
         end do                                         
       else                                             
         do i1=1,3                                      
           i2  = Face(k)%vertices(i1)                   
           X   = Vert(i2,:)                             
           ds  = X - Xp                                 
           tmp = tmpp + dot_product( dPhidX(ip,:), ds )   
           NodalData(i2)  = NodalData(i2)  + tmp       
           NodalCounter(i2) = NodalCounter(i2) + 1      
         end do                                         
       endif                                            
     end do    
   end do

   do i=1,Nvrt
     if( NodalCounter(i) > 0 )then
       NodalData(i)  = NodalData(i) /float(NodalCounter(i))
     endif
   end do

   !
   ! mode = 1: do not use boundary values, cell centered data only
   ! useful for particle tracking
   !
   if( Mode == 1 ) return
   
   !
   ! boundary faces overrule! (no gradient, simple average)
   !
   ! first loop over the nodes to reset them, reset nodalcounter too
   !
   if( Debug > 2 ) write(*,*)'InterpolateData: boundaries of ',Variable(IVar)

   NodalCounter = 0
   
   do ib=1,Nbnd           
     i  = Bnd(ib)%face
     i4 = 3
     if( Face(i)%vertices(4) > 0 ) i4 = 4
     
     do i1=1,i4
       i2  = Face(i)%vertices(i1)
       NodalData(i2)  = 0.0
     end do 
   end do
   
   do ib=1,Nbnd           
     i  = Bnd(ib)%face
     i4 = 3
     if( Face(i)%vertices(4) > 0 ) i4 = 4
     
     do i1=1,i4
       i2  = Face(i)%vertices(i1)
       NodalData(i2)  = NodalData(i2)  + Phi(Ncel+ib)
       NodalCounter(i2) = NodalCounter(i2) + 1
     end do 
   end do

   do i=1,Nvrt
     if( NodalCounter(i) > 0 )then
       NodalData(i)  = NodalData(i) /float(NodalCounter(i))
     endif
   end do

   if( Debug > 2 ) write(*,*)'=== InterpolateData'

end subroutine InterpolateData

