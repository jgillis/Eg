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
program dolfyn

   use constants
   use geometry
   use variables
   use scalars

   real    :: TimeElapsed, TimeNew
   logical :: Ready, Exists, Divergence

   Exists      = .false.
   Divergence  = .false.
   
   call cpu_time( TimeElapsed )
   TimeNew = -1.0
      
   call ReadCase   
   !
   ! first read geometry (needed to set certain constants)
   !
   call ReadGeometry

   call ReadControlFile

   if( UsePatches ) call PatchesSetUp(1)

   !
   ! initialisation and calculate geometry data
   !
   call InitialiseVariables

   !call InitializeProperties

   if( UsePatches ) call PatchesSetUp(2)  ! stage 2, variables allocated

   if( UseParticles ) call ParticleSetUp 

   call SetBoundaryConditions
      
   call PrettyPrint(IOdef)
   call PrettyPrint(IOdbg)
      
   if( count( Solver(:) == SparseKit2 ) > 0 )then
     write(*,*)'Allocating SparsKit2 Work array'
     allocate( Work(NNZ,8),stat=istat)    
     call TrackMemory(istat,(NNZ)*8*2,'Work array allocated')
   endif

   IterStart = 0
   Time      = 0.0
   dppmax    = 0.0
   
   if( Restart > 0 )then

     write(*,*)'Reading restart data, ',Restart
     write(IOdbg,*)'Reading restart data'
     call ReadRestartField(IterStart)
     
     if( Restart == 2 .or. Restart == 3 )then
       write(*,*) 'Overriding boundary conditions'
       write(IOdbg,*) 'Overriding boundary conditions'
       
       IterStart = 0
       Time      = 0.0
       call SetBoundaryConditions
       call Set_Normalisation_Factors
       if( Debug > 0 )write(*,*) 'Opening residuals file (1)'
       call openfile(IOres,Casename,'.res','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)
       call Set_Normalisation_Factors

       call InitialField(dppmax)
  
       if( Debug > 0 )write(*,*) 'Opening monitor file (1)'
       call openfile(IOmon,Casename,'.mon','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)

     else
       if( Debug > 0 )write(*,*) 'Opening residuals file (2)'
       call openfile(IOres,Casename,'.res','FORMATTED','APPEND','UNKNOWN',Debug)
       call Set_Normalisation_Factors

       if( Debug > 0 )write(*,*) 'Opening monitor file (2)'
       call openfile(IOmon,Casename,'.mon','FORMATTED','APPEND','UNKNOWN',Debug)
     endif

   else
     if( Debug > 0 )write(*,*) 'Opening residuals file (3) ',Restart
     call openfile(IOres,Casename,'.res','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)
     call Set_Normalisation_Factors
   
     call InitialField(dppmax)
   
     if( Debug > 0 )write(*,*) 'Opening monitor file (3)'
     call openfile(IOmon,Casename,'.mon','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)
   endif
      
   !
   ! the big time loop
   !
   if( Debug > -1 ) write(*,'(//,'' Go! '',/)')

   TimeSteps: do iter=IterStart+1,IterStart+Niter
     
     !
     ! reset the array of konsole output flags
     !
     Flags(1:NFlags) = ' '
     
     Time = Time + dt
     
     if( Debug > 1 )then
       write(*,'(1x,A)')     '****************************'
       if( Transient )then
         write(*,'(1x,A,i5,1x,f12.5)')  '*** Time step:  ',iter,time
       else
         write(*,'(1x,A,i5)')  '*** Iteration:  ',iter
       endif
       write(*,'(1x,A)')     '****************************'
     endif
     write(IOdbg,'(1x,A)')   '**********************************************'
     write(IOdbg,'(1x,A,i5,1x,f12.5)')  '*** Iteration/Time step:  ',iter,time
     write(IOdbg,'(1x,A)')   '**********************************************'
          
     ! 
     ! store old solutions etc
     !
     if( Transient )then
       if( Debug > 2 )write(*,*) 'Storing old values'
       if( QuadTime )then
         Uold2 = Uold
         Vold2 = Vold
         Wold2 = Wold
       endif
       Uold = U
       Vold = V
       Wold = W
       if( SolveEnthalpy .and. QuadTime ) Told2 = Told
       if( SolveEnthalpy ) Told = T
       if( SolveTurb )then
         if( QuadTime )then
           TEold2 = TEold
           EDold2 = EDold
         endif
         TEold = TE
         EDold = ED
       endif
       if( SolveScalars )then
         if( QuadTime ) SCold2 = SCold
         SCold = SC
       endif
       if( MovingGrid )then
         write(*,*)'Not implemented yet'
       endif
       call Set_Normalisation_Factors  ! inlet might be changed
     endif

     MaxOUTER = 40                     ! <= input deck
     iouter   =  0
     ready    = .false.
     OuterSteps: do while( .not. ready )
       
       iouter = iouter + 1
       !
       ! velocity gradients 
       !                                      
       call GradientPhi(VarU,U,dUdX)      !
       call GradientPhi(VarV,V,dVdX)      ! supporting stuff
       call GradientPhi(VarW,W,dWdX)      ! needed at various places
       call GradientPhi(VarP,P,dPdX)      !

       if( SolveUVW) call CalculateUVW 

       if( SolveP )  call CalculatePressure(dppmax)

       call UpdateBC

       if( SolveTurbEnergy ) call CalculateScalar(VarTE,TE,dpdx,TEold,TEold2)
       if( SolveTurbDiss )   call CalculateScalar(VarED,ED,dpdx,EDold,EDold2)
       if( SolveVisc )       call CalculateViscosity

       if( SolveEnthalpy )   call CalculateScalar(VarT,T,dpdx,Told,Told2)

       if( SolveScalars )then 
         do is=1,Nscal        
           call CalculateScalar(VarS(is),SC(:,is),dpdx,SCold(:,is),SCold2(:,is))
         end do
       endif

       !
       ! fluid properties (den,cp,vis) will be set here
       !
       
       !if( UseParticles )    call ParticleInterface ! only when coupled

       call PrintSummary(IOdbg,iouter)

       !
       ! stopping criteria
       !          
       if( iouter >= MaxOUTER .or. (.not. Transient) ) ready = .true.

       tmpres = Residual(6)   ! do not use the epsilon-residual
       Residual(6) = Small    ! because its nature is arbitrary

       if( maxval(Residual(1:8)) < ResMax )then
         if( Debug > 1 .or. .not. Transient )then
           write(*,'(/)')
           write(*,*)'*** CONVERGENCE ***',maxval(Residual(1:8))
         endif
         Residual(6) = tmpres
         if( Transient )then
           exit OuterSteps
         else
           exit TimeSteps
         endif
       endif
     
       if( maxval(Residual(1:8)) > 1.e+18 )then
         write(*,'(/)')
         write(*,*)'*** DIVERGENCE ***',maxval(Residual(1:8))
         Divergence = .true.
         exit TimeSteps
       endif

       Residual(6) = tmpres

     end do OuterSteps

     call PrintSummary(IOdef,iouter)
     
     !write(IOext,*) TIME, Qtransfer
     
     if( Debug > 2 ) write(*,*)'Maximum change in PP:', &
                                dppmax,maxval(Residual(1:8)),ResMax
     
     !if( dppmax < Small .and. iter > 20 .and. .not. Transient) EXIT timesteps

     !
     ! intermediate opendx debug dump
     !
     !if( mod(Iter,DXdump) == 0 .and. Iter /= IterStart+Niter )then
     !
     !  call dolfyn2opendx
     !
     !endif
     !
     ! intermediate output data options
     !
     if( mod(Iter,NOutput) == 0 .and. Iter /= IterStart+Niter )then

       write(*,*)'Writing postprocessing file(s)' 
       if( UseOpenDX  ) call dolfyn2opendx
       if( UseVTK     ) call dolfyn2vtk
       if( UseTECPLOT ) call dolfyn2tecplt(0)
       if( UseGMV     ) call dolfyn2gmv

     endif

     if( Iter == IOutput)then
       write(*,*)'Writing postprocessing file(s) at iteration ',Iter 
       if( UseOpenDX  ) call dolfyn2opendx
       if( UseVTK     ) call dolfyn2vtk
       if( UseTECPLOT ) call dolfyn2tecplt(0)
       if( UseGMV     ) call dolfyn2gmv
     endif
     
     if( TOutput /= -1. .and. Time >= TOutput .and. .not. TOutputDone )then
       write(*,*)'Writing postprocessing file(s) at time ',Time
       if( UseOpenDX  ) call dolfyn2opendx
       if( UseVTK     ) call dolfyn2vtk
       if( UseTECPLOT ) call dolfyn2tecplt(0)
       if( UseGMV     ) call dolfyn2gmv
       TOutputDone = .true.
     endif 
     !
     ! saving restart data options
     !
     if( mod(Iter,NSave) == 0 .and. Iter /= IterStart+Niter )then

       call WriteRestartField
     
     endif
     
     if( Iter == ISave )then
       write(*,*)'Saving restart data at iteration ',Iter 
       call WriteRestartField
     endif
     
     if( TSave > 0.0 .and. Time >= Tsave .and. .not. TSaveDone )then
       write(*,*)'Saving restart data at time ',Time 
       call WriteRestartField
       TSaveDone = .true.
     endif 

     if( TCPU > 0.0 )then
       if( TimeNew < 0.0 ) TimeNew = TCPU              ! if not set, do it
       call cpu_time( TimeElapsed )                    ! get current time
     
       if( TimeElapsed >= TimeNew )then
         if( CPUstop )then
           write(*,*)'*** CPU time exceeded ***'
           write(IOdbg,*)'*** CPU time exceeded ***'
           exit timesteps
         else
           TimeNew = TimeNew + TCPU                      ! set new save time
           write(*,*)'Saving restart data CPU time ',TimeElapsed 
           call WriteRestartField
         endif
       endif
     endif 
     !
     ! trick to stop a running job 
     !
     inquire(file='STOP',exist=Exists)
     if( Exists )then
       write(*,*)'*** OK stop ***'
       write(IOdbg,*)'*** OK stop ***'
       exit timesteps
     endif
     
     !
     ! option to fix values for ABL
     !
     if( UseFixABL )then
       !write(*,*) 'Fixing values in ABL citd ',IdFixABL     
       write(IOdbg,*) 'Fixing values in ABL citd ',IdFixABL     
       do i=1,Ncel
         if( cell(i)%ctid == IdFixABL )then
           U(i)  = UFixABL 
           V(i)  = VFixABL 
           W(i)  = WFixABL 
           TE(i) = TeFixABL   
           ED(i) = EdFixABL  
         endif
       end do
     endif
     !
     ! warning: not a standard call
     !
     call flush(IOdbg)
     
   end do timesteps

   if( Iter == IterStart+Niter+1 )then
     iter = iter - 1
     write(*,'(/)')
     write(*,*)'Number of requested steps done. ',Niter
     write(*,'(/)')
     write(IOdbg,'(/)')
     write(IOdbg,*)'Number of requested steps done. ',Niter
     write(IOdbg,'(/)')
   endif
      
   if( .not. Divergence ) call WriteRestartField
      
   if( PrintCellVar .or. PrintWallVar) call ShowCells
   
   !
   ! clean up a bit before the external interfaces are called
   !   
   if( allocated( Work ) )then
     deallocate( Work )
     call TrackMemory(istat,-NNZ*8*2,'SparsKit2 work array deallocated')
   endif

   if( UseParticles ) call ParticleInterface      ! only passive particles!


   if( UseOpenDX  ) call dolfyn2opendx
   if( UseVTK     ) call dolfyn2vtk
   if( UseTECPLOT ) call dolfyn2tecplt(1)
   if( UseGMV     ) call dolfyn2gmv

   if( Debug > -1 ) write(*,*)'Maximum change in PP:', &
                           dppmax,maxval(Residual(1:8)),ResMax
   
   call CleanUp
   
   write(*,*)'Done ',casename(1:lens(casename))
   write(IOdbg,*)'Done.'

   
end program Dolfyn
subroutine ReadGeometry

   use constants
   use geometry

   real, dimension(3)     :: Xp, Xn, X 
   real, dimension(3)     :: X1, X2, X3, Xcg 

   integer, parameter     :: Ncg = 30
   real, dimension(Ncg)   :: TetVol
   real, dimension(Ncg,3) :: TetSurf, TetNorm, TetCG

   real :: pi = 3.1415927

   character(len=32) string
   integer      geoversion
   logical      bin  

   logical      warning
   !
   ! extra information will be dumped in this file:
   !
   call openfile(IOdbg,casename,'.dbg','FORMATTED','APPEND','UNKNOWN',debug)

   call Disclaimer
   
   !
   ! inquire for the type of dolfyn geometry file
   !
   call feelfile(IOgeo,IFORM,casename,'.geo','FORMATTED','SEQUENTIAL','OLD',debug)
   
   if( IFORM == 0 )then
     bin = .true.
   else
     bin = .false.
   endif
   !
   ! read dolfyn files 
   !
   write(*,*) 'Opening geometry file'
   
   if( bin)then
     call openfile(IOgeo,casename,'.geo','UNFORMATTED','SEQUENTIAL','OLD',debug)
   else
     call openfile(IOgeo,casename,'.geo','FORMATTED','SEQUENTIAL','OLD',debug)
   endif

   if( bin ) then
     !
     ! binary
     !
     read(IOgeo) string(1:12)
     if( string(1:12) /= 'dolfyn bin g' )then
       write(*,*)'*** Error: Not a dolfyn geometry file'
       stop
     endif
     read(IOgeo) geoversion
     if( geoversion /= version )then
       write(*,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOdbg,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
     endif
     read(IOgeo) ScaleFactor
     write(IOdbg,*) 'Using ScaleFactor = ',ScaleFactor

     read(IOgeo) string(1:4)
     if( string(1:3) /= 'ce:' )then
       write(*,*)'*** Error: Cell definitions expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     read(IOgeo) icel
     write(IOdbg,*)'Reading cells:',icel,' BINARY'

     Ncel = icel
     allocate(Cell(Ncel),stat=istat)          
     call TrackMemory(istat,Ncel,'Cell array allocated')

     do i=1,Ncel
       read(IOgeo) j,cell(j)%ctid,cell(j)%x,cell(j)%vol
       if( j /= i ) write(*,*)'*** Warning: cell array corrupted'
     end do   
   else
     !
     ! ascii
     !
     read(IOgeo,'(a12)') string(1:12)
     if( string(1:12) /= 'dolfyn asc g' )then
       write(*,*)'*** Error: Not a dolfyn geometry file'
       stop
     endif
     read(IOgeo,'(8x,i6)') geoversion
     if( geoversion /= version )then
       write(*,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOdbg,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
     endif
     read(IOgeo,'(6x,e12.5)') ScaleFactor
     write(IOdbg,*) 'Using ScaleFactor = ',ScaleFactor

     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'ce:' )then
       write(*,*)'*** Error: Cell definitions expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     read(IOgeo,*) icel
     write(IOdbg,*)'Reading cells:',icel

     Ncel = icel
     allocate(Cell(Ncel),stat=istat)          
     call TrackMemory(istat,Ncel,'Cell array allocated')

     do i=1,Ncel
       read(IOgeo,'(2(i8,1x),4(1pe16.9,1x))') &
                                    j,cell(j)%ctid,cell(j)%x,cell(j)%vol
       if( j /= i ) write(*,*)'*** Warning: cell array corrupted'
     end do   
   endif
   
   !
   ! apply scalefactor
   !
   Cell(1:Ncel)%x(1) = ScaleFactor * Cell(1:Ncel)%x(1)
   Cell(1:Ncel)%x(2) = ScaleFactor * Cell(1:Ncel)%x(2)
   Cell(1:Ncel)%x(3) = ScaleFactor * Cell(1:Ncel)%x(3)
   
   ScaleFactor3 = ScaleFactor**3

   Cell(1:Ncel)%vol = ScaleFactor3 * Cell(1:Ncel)%vol
   
   !
   ! some output
   !
   if( debug > 0 )then
     write(*,'(1x,A,i2,A,i2)') 'Cell types from ', &
                         minval(cell(:)%ctid),' to ',maxval(cell(:)%ctid)
     write(*,'(1x,A)') 'Cell centers are:'
     write(*,'(1x,A,1pe10.3,A,e10.3)') &
       'Dimension ',minval(cell(:)%x(1)),' <  x  < ',maxval(cell(:)%x(1))
     write(*,'(1x,A,1pe10.3,A,e10.3)') &
       'Dimension ',minval(cell(:)%x(2)),' <  y  < ',maxval(cell(:)%x(2))
     write(*,'(1x,A,1pe10.3,A,e10.3)') &
       'Dimension ',minval(cell(:)%x(3)),' <  z  < ',maxval(cell(:)%x(3))
     write(*,'(1x,A,1pe10.3,A,e10.3)') &
       'Volumes   ',minval(cell(:)%vol), ' < Vol < ',maxval(cell(:)%vol)  
     write(*,'(1x,A,1pe10.3,A,e10.3)') &
       'Total volume: ',sum(cell(:)%vol)  
   endif
   write(IOdbg,'(1x,A,i2,A,i2)') 'Cell types from ', &
                       minval(cell(:)%ctid),' to ',maxval(cell(:)%ctid)
   write(IOdbg,'(1x,A)') 'Cell centers are:'
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(cell(:)%x(1)),' <  x  < ',maxval(cell(:)%x(1))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(cell(:)%x(2)),' <  y  < ',maxval(cell(:)%x(2))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(cell(:)%x(3)),' <  z  < ',maxval(cell(:)%x(3))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Volumes   ',minval(cell(:)%vol), ' < Vol < ',maxval(cell(:)%vol)  
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Total volume: ',sum(cell(:)%vol)  


   !
   ! list of cell faces
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:4) /= 'cf: ' )then
       write(*,*)'*** Error: List of cell faces expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading list of cell faces'

     read(IOgeo) maxfaces
     if( maxfaces > 6 )then
       write(*,*)'*** Warning: Unsupported feature'
     endif
     allocate(NFaces(Ncel),stat=istat)          
     call TrackMemory(istat,Ncel,'NFaces array allocated')

     allocate(CFace(Ncel,12),stat=istat)          
     call TrackMemory(istat,Ncel*12,'CFace array allocated')

     do i=1,Ncel 
       read(IOgeo) j,NFaces(j)
       do k=1,Nfaces(j)
         read(IOgeo) CFace(j,k)
       end do
       if( j /= i ) write(*,*)'*** Warning: cell face array corrupted'
     end do   
     
   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'cf:' )then
       write(*,*)'*** Error: List of cell faces expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading list of cell faces'

     read(IOgeo,*) maxfaces
     if( maxfaces > 6 )then
       write(*,*)'*** Warning: Unsupported feature'
     endif
     allocate(NFaces(Ncel),stat=istat)          
     call TrackMemory(istat,Ncel,'NFaces array allocated')

     allocate(CFace(Ncel,12),stat=istat)          
     call TrackMemory(istat,Ncel*12,'CFace array allocated')

     do i=1,Ncel
       read(IOgeo,'(12(i8,1x))') j,NFaces(j),(CFace(j,k),k=1,NFaces(j))
       if( j /= i ) write(*,*)'*** Warning: cell face array corrupted'
     end do   
   endif

   !
   ! the faces
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:3) /= 'fc: ' )then
       write(*,*)'*** Error: Cell faces expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading cell faces'

     read(IOgeo) Nfac
     allocate(Face(Nfac),stat=istat)          
     call TrackMemory(istat,Nfac*16,'Faces array allocated')

     allocate(RFace(Nfac,2),stat=istat)          
     call TrackMemory(istat,Nfac*2,'RFace array allocated')

     do i=1,Nfac
       read(IOgeo) j, face(j)
     end do
   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'fc:' )then
       write(*,*)'*** Error: Cell faces expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading cell faces'
   
     read(IOgeo,*) Nfac
     allocate(Face(Nfac),stat=istat)          
     call TrackMemory(istat,Nfac*16,'Faces array allocated')

     allocate(RFace(Nfac,2),stat=istat)          
     call TrackMemory(istat,Nfac*2,'RFace array allocated')

     do i=1,Nfac
       read(IOgeo,'(8(i8,1x),8(1pe16.9,1x))') j, face(j)
     end do
   endif

   !
   ! apply scalefactor
   !
   ScaleFactor2 = ScaleFactor**2

   Face(1:Nfac)%area = ScaleFactor2 * Face(1:Nfac)%area
   
   Face(1:Nfac)%n(1) = ScaleFactor2 * Face(1:Nfac)%n(1)
   Face(1:Nfac)%n(2) = ScaleFactor2 * Face(1:Nfac)%n(2)
   Face(1:Nfac)%n(3) = ScaleFactor2 * Face(1:Nfac)%n(3)

   Face(1:Nfac)%x(1) = ScaleFactor * Face(1:Nfac)%x(1)
   Face(1:Nfac)%x(2) = ScaleFactor * Face(1:Nfac)%x(2)
   Face(1:Nfac)%x(3) = ScaleFactor * Face(1:Nfac)%x(3)
   
   !
   ! boundaries
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:4) /= 'bn: ' )then
       write(*,*)'*** Error: Boundary faces expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading boundaries'

     read(IOgeo) Nbnd                           ! , Nreg <=hoeft niet!
     allocate(Bnd(Nbnd),stat=istat)          
     call TrackMemory(istat,Nbnd,'Boundary array allocated')

     do i=1,Nbnd
       read(IOgeo) j,bnd(j)%face, &
                     bnd(j)%vertices, bnd(j)%rid, bnd(j)%distance
       bnd(j)%yplus =  0.0
       bnd(j)%uplus =  0.0
       bnd(j)%shear =  0.0
       bnd(j)%h     =  0.0
       bnd(j)%q     =  0.0
     end do
   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'bn:' )then
       write(*,*)'*** Error: Boundary faces expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading boundaries'

     read(IOgeo,*) Nbnd                           ! , Nreg <=hoeft niet!
     allocate(Bnd(Nbnd),stat=istat)          
     call TrackMemory(istat,Nbnd,'Boundary array allocated')

     do i=1,Nbnd
       read(IOgeo,'(7(i8,1x),1(e12.5,1x))') j,bnd(j)%face, &
                     bnd(j)%vertices, bnd(j)%rid, bnd(j)%distance
       bnd(j)%yplus =  0.0
       bnd(j)%uplus =  0.0
       bnd(j)%shear =  0.0
       bnd(j)%h     =  0.0
       bnd(j)%q     =  0.0
     end do   
   endif
   !
   ! apply scalefactor
   !
   Bnd(1:Nbnd)%distance = ScaleFactor * Bnd(1:Nbnd)%distance
   
   
   Nreg = maxval( Bnd(:)%rid )
   write(*,*) 'Maximum region number found: ',Nreg
   write(IOdbg,*) 'Maximum region number found: ',Nreg
   allocate(Reg(0:Nreg),stat=istat)          
   call TrackMemory(istat,(1+Nreg)*17, &
     'Boundary region definitions array allocated')

   !
   ! finally the vertices
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:4) /= 'vr: ' )then
       write(*,*)'*** Error: Vertices expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading Vertices'

     read(IOgeo) Nvrt
     allocate(Vert(Nvrt,3),stat=istat)          
     call TrackMemory(istat,Nvrt,'Vertices array allocated')

     do i=1,Nvrt
       read(IOgeo) j,(vert(j,k),k=1,3)
     end do
   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'vr:' )then
       write(*,*)'*** Error: Vertices expected'
       write(*,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading Vertices'

     read(IOgeo,*) Nvrt
     allocate(Vert(Nvrt,3),stat=istat)          
     call TrackMemory(istat,Nvrt,'Vertices array allocated')

     do i=1,Nvrt
       read(IOgeo,'(1(i8,1x),3(1pe16.9,1x))') j,(vert(j,k),k=1,3)
     end do
   endif
   !
   ! apply scalefactor
   !
   Vert = ScaleFactor * Vert
   
   !
   ! *** consistency checks ***
   !
   write(*,*)    'Geometry checks'
   write(IOdbg,*)'Geometry checks'
   !
   ! check volume using Gauss (see eq. 8.42)
   !
   vol1 = 0.0
   vol2 = 0.0
   vol3 = 0.0
   do i=1,Ncel
   
     Xp = Cell(i)%x
     v1 = 0.0
     v2 = 0.0
     v3 = 0.0

     do j=1,NFaces(i)
       k  = CFace(i,j)
       ip = Face(k)%cell1
       in = Face(k)%cell2

       if( ip == i )then 
         v1 = v1 + Face(k)%x(1) * Face(k)%n(1)
         v2 = v2 + Face(k)%x(2) * Face(k)%n(2)
         v3 = v3 + Face(k)%x(3) * Face(k)%n(3)
       else if( in == i )then
         v1 = v1 - Face(k)%x(1) * Face(k)%n(1)
         v2 = v2 - Face(k)%x(2) * Face(k)%n(2)
         v3 = v3 - Face(k)%x(3) * Face(k)%n(3)
       else
         write(*,*)'Error in check volume consistency'
       endif
     end do
     vol1 = vol1 + v1
     vol2 = vol2 + v2
     vol3 = vol3 + v3
   end do
   if( abs(vol1-vol2) > 1.e-3 * vol1 .or. &
       abs(vol1-vol3) > 1.e-3 * vol1 .or. & 
       abs(vol2-vol3) > 1.e-3 * vol1 )then
     write(*,*)    'WARNING: Check volumes'
     write(IOdbg,*)'WARNING: Check volumes'
   endif 

   write(IOdbg,*)'Volume based on xi:',vol1 
   write(IOdbg,*)'Volume based on yj:',vol2
   write(IOdbg,*)'Volume based on zk:',vol3

   !write(IOdbg,*)'Checking centroids'
   !
   !warning = .false.
   !do i=1,Ncel
   !
   !  Xp = Cell(i)%x
   !  X  = 0.0
   !
   !  do j=1,NFaces(i)
   !    k  = CFace(i,j)
   !    X = X + (Face(k)%x - Xp)*Face(i)%n
   !    !X = X + (Face(k)%x - Xp)*Face(k)%n ! bug found by Jojo
   !  end do
   !  if( abs(X(1)) > 1.e-2 .or. &
   !      abs(X(2)) > 1.e-2 .or. &
   !      abs(X(3)) > 1.e-2 )then
   !    write(IOdbg,*) ' '
   !    write(IOdbg,*) 'centroid',i,X
   !    write(IOdbg,*) '  nfaces',NFaces(i),Xp
   !    warning = .true.
   !  endif 
   !end do
   !
   ! alternative
   !
   !do i=1,Ncel
   !
   !  Xp = Cell(i)%x
   !  X  = 0.0
   ! 
   !  itet = 0
   !  do j=1,NFaces(i)
   !    k  = CFace(i,j)
   !    !write(*,*)'cell/face:',i,j,k
   !    !
   !    ! face still has a maximum of 4 vertices !!!
   !    ! build first tetrahedron
   !    !
   !    i1   = Face(k)%vertices(1)
   !    i2   = Face(k)%vertices(2)
   !    i3   = Face(k)%vertices(3)
   !    i4   = Face(k)%vertices(4)
   !
   !    x1   = Vert(i1,:)
   !    x2   = Vert(i2,:)
   !    x3   = Vert(i3,:) 
   !
   !    vol  = determinant( (x1(1)-xp(1)),(x2(1)-xp(1)),(x3(1)-xp(1)), &
   !                        (x1(2)-xp(2)),(x2(2)-xp(2)),(x3(2)-xp(2)), &
   !                        (x1(3)-xp(3)),(x2(3)-xp(3)),(x3(3)-xp(3))  )
   !
   !    vol  = 1./6. * vol
   !    
   !    call cross_product((x2-x1),(x3-x1),Xn)
   !
   !    itet = itet + 1
   !
   !    TetVol(itet)    = abs(vol)
   !    TetSurf(itet,:) = Xn
   !
   !    call normalise( Xn )       
   !    TetNorm(itet,:) = Xn
   !     
   !    TetCG(itet,:)   = 0.25*( Xp + x1 + x2 + x3 ) 
   !    !
   !    ! second tetrahedron
   !    !
   !    if( i4 > 0 )then
   !
   !      itmp = i1
   !      i1   = i3
   !      i2   = i4
   !      i3   = itmp
   !
   !      x1   = Vert(i1,:)
   !      x2   = Vert(i2,:)
   !      x3   = Vert(i3,:)
   !
   !      vol  = determinant( (x1(1)-xp(1)),(x2(1)-xp(1)),(x3(1)-xp(1)), &
   !                          (x1(2)-xp(2)),(x2(2)-xp(2)),(x3(2)-xp(2)), &
   !                          (x1(3)-xp(3)),(x2(3)-xp(3)),(x3(3)-xp(3))  )
   !
   !      vol  = 1./6. * vol
   !
   !      call cross_product((x2-x1),(x3-x1),Xn)
   !
   !      itet = itet + 1
   !
   !      TetVol(itet)    = abs(vol)
   !      TetSurf(itet,:) = Xn
   !
   !      call normalise( Xn )       
   !      TetNorm(itet,:) = Xn
   !
   !      TetCG(itet,:)   = 0.25*( Xp + x1 + x2 + x3 ) 
   !
   !    endif
   !    if( itet > Ncg ) stop 'Increase parameter Ncg'
   !  end do
   !
   !  vol  = sum( TetVol(1:itet) )          
   !  dvol = abs( Cell(i)%vol - vol )              ! absoluut, moet relatief
   !  if( dvol > 1.e-3 )then
   !    write(*,*)'volume:',i,Cell(i)%vol,':',vol
   !  endif
   !  
   !  Xcg = 0.0
   !  do j=1,itet
   !    Xcg = Xg + TetVol(j)*( TetCG(j,:) - Xp )
   !  end do
   !  Xcg = Xcg / vol
   !  
   !  dcg = vector_length( Xcg )                   ! absoluut, moet relatief
   !  if( dcg > 1.e-2 )then
   !    write(*,*)'centroid:',i,itet,dcg,vol
   !    write(*,*)'      Xp:',Xp
   !    write(*,*)'    dXcg:',Xcg
   !  endif
   !
   !end do
   !
   !if( warning ) write(*,*)'WARNING: Check centroids'
   
   !
   ! check for a closed boundary
   !

   !
   ! check non-orthogonality
   !
   anglemx  =  0.0
   anglemn  = 90.0
   anglelm  = 70.0
   iangle   =   0

   do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )then 
       ! internal face
       Xp = Cell(ip)%x   ! cell centre face owner
       Xn = Cell(in)%x   ! cell centre neigbor

       Xp = Xn - Xp      ! vector from P to N
       Xn = Face(i)%n    ! surface vector

       call normalise(Xp)
       call normalise(Xn)

       dotp  = dot_product( Xp , Xn )
       ! BT, Visual Fortran does not like acos(1.0)
       if(dotp > 0.99999 ) dotp = 0.99999
       angle = acos( dotp )*180./pi     
     else
       Xp = Cell(ip)%x   ! cell centre face owner
       Xn = Face(i)%x    ! face centre  

       Xp = Xn - Xp      ! vector from P to N
       Xn = Face(i)%n    ! surface vector

       call normalise(Xp)
       call normalise(Xn)

       dotp  = dot_product( Xp , Xn )
       ! BT, Visual Fortran does not like acos(1.0)
       if(dotp > 0.99999 ) dotp = 0.99999
       angle = acos( dotp )*180./pi     
     endif
     anglemx =  max(anglemx,angle)
     anglemn =  min(anglemn,angle)
     if( angle > anglelm ) iangle = iangle + 1
   end do
   
   write(*,*) 'Angles:',anglemn,anglemx,iangle
   
   
   write(IOdbg,'(1x,A)') 'Vertices are:'
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,1)),' <  x  < ',maxval(vert(:,1))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,2)),' <  y  < ',maxval(vert(:,2))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,3)),' <  z  < ',maxval(vert(:,3))
   !
   ! final stuff
   ! 
   close(IOgeo)
   write(IOdbg,*)'End reading geometry data'
   write(IOdbg,*)'Current allocated memory (words): ',memory
   write(IOdbg,*)'Status:'
   write(IOdbg,'(1x,A,i8)') 'Number of cells      =',Ncel
   write(IOdbg,'(1x,A,i8)') '          Special''s  =',count( Nfaces > 6 )
   write(IOdbg,'(1x,A,i8)') '             Hexa''s  =',count( Nfaces == 6 )
   write(IOdbg,'(1x,A,i8)') '    Wedges, Pyramids =',count( Nfaces == 5 )
   write(IOdbg,'(1x,A,i8)') '             Tetra''s =',count( Nfaces == 4 )
   write(IOdbg,'(1x,A,i8)') 'Number of vertices   =',Nvrt
   write(IOdbg,'(1x,A,i8)') 'Number of faces      =',Nfac
   Nint = count( Face(:)%cell2 /= 0)
   write(IOdbg,'(1x,A,i8)') '      Internal faces =',Nint
   if( (Nfac-Nint) /= Nbnd ) write(IOdbg,*)'Error: inconsistent geometry data'
   write(IOdbg,'(1x,A,i8)') 'Number of boundaries =',Nbnd
   do i=0,Nreg
     write(IOdbg,'(1x,A,i2,A,i8,A)') &
               '           Region ',i,' =',count( Bnd(:)%rid == i)
   end do
   write(*,'(1x,A)') 'End reading geometry'
   
end subroutine ReadGeometry
subroutine InitialiseVariables
!========================================================================

   use constants
   use geometry
   use variables
   use scalars

   Small = epsilon(U(1)) * 0.1
   Large = huge(large) * 1.e-3
   
   SMALL = 1.e-12
   
   if( Debug > -1 )then
     write(IOdbg,*)'Constants small/large:',Small,Large
     write(IOdbg,*)'Tiny     :',tiny(small) 
     write(IOdbg,*)'Epsilon  :',epsilon(U(1))
     write(IOdbg,*)'Precision:',precision(small),digits(small), &
                            minexponent(small),maxexponent(small)
   endif
   
   i = 0
   allocate( U(Ncel+Nbnd),stat=istat)
   i = i + istat 
   allocate( V(Ncel+Nbnd),stat=istat )
   i = i + istat 
   allocate( W(Ncel+Nbnd),stat=istat )
   istat = i + istat 
   call TrackMemory(istat,(Ncel+Nbnd)*3,&
                   'Velocity components array allocated')

   U  = Guess(VarU)
   V  = Guess(VarV)
   W  = Guess(VarW)

   i = 0
   allocate( P(Ncel+Nbnd),stat=istat)
   i = i + istat 
   allocate( PP(Ncel+Nbnd),stat=istat)
   istat = i + istat 
   call TrackMemory(istat,(Ncel+Nbnd)*2,'Pressure array allocated')

   P  = Guess(VarP)
   PP = 0.0
   
   i = 0
   allocate( dudx(Ncel+Nbnd,3),stat=istat )
   i = i + istat 
   allocate( dvdx(Ncel+Nbnd,3),stat=istat )
   i = i + istat 
   allocate( dwdx(Ncel+Nbnd,3),stat=istat )
   i = i + istat 
   allocate( dpdx(Ncel+Nbnd,3),stat=istat )
   istat = i + istat 
   call TrackMemory(istat,(Ncel+Nbnd)*12,'Gradient arrays allocated')

   dudx = 0.0
   dvdx = 0.0
   dwdx = 0.0

   dpdx = 0.0

   allocate( Den(Ncel+Nbnd),stat=istat)
   call TrackMemory(istat,Ncel+Nbnd,'Density array allocated')

   Den = DensRef                                                  !<=== let op material 1 

   allocate( Cp(Ncel+Nbnd),stat=istat)
   call TrackMemory(istat,Ncel+Nbnd,'Cp array allocated')

   Cp  = CpStd                                                    !<=== let op material 1 

   allocate( VisEff(Ncel+Nbnd),stat=istat)
   call TrackMemory(istat,Ncel+Nbnd,'Effective viscosity array allocated')

   VisEff  = VisLam                                               !<=== let op material 1 

   if( SolveEnthalpy )then
     allocate( T(Ncel+Nbnd),stat=istat)
     call TrackMemory(istat,Ncel+Nbnd,'Temperature array allocated')

     T = Tref          !Tref wordt 0! we werken straks relatief en per materiaal!
   
!   else
!     
!     allocate( T(1),stat=istat)     ! temp. workaround of g95-bug!!!
!     
   endif

   if( SolveTurb )then
     i = 0
     allocate( TE(Ncel+Nbnd),stat=istat )
     i = i + istat 
     allocate( TurbP(Ncel+Nbnd),stat=istat )
     i = i + istat 
     allocate( ED(Ncel+Nbnd),stat=istat )
     i = i + istat 
     call TrackMemory(istat,(Ncel+Nbnd)*3,'K-eps model arrays allocated')

     TE    = Guess(VarTE)
     TurbP = 0.0
     
     ED    = Guess(VarED)

!   else
!     
!     allocate( TE(1),stat=istat)     ! temp. workaround of g95-bug!!!
!     allocate( ED(1),stat=istat)     ! temp. workaround of g95-bug!!!
!          
   endif
 
   if( SolveScalars )then

     write(*,*)'Allocating ',Nscal,' scalar arrays'
     allocate( SC(Ncel+Nbnd,Nscal),stat=istat)
     call TrackMemory(istat,(Ncel+Nbnd)*Nscal,&
                     'Scalar array(s) allocated')

     SC = 0.0
   
   endif
 
   if( DXdebugdata )then
     allocate( DXdebug(Ncel),stat=istat)
     dxdebug = 0.0
     
     allocate( DXgrad(Ncel,3),stat=istat)
     call TrackMemory(istat,Ncel,'DX debug arrays allocated')
     dxgrad = 0.0
   endif

   allocate( dPhidXo(Ncel+Nbnd,3),stat=istat)  ! vraag? moet Nbnd???
   call TrackMemory(istat,Ncel*6,'Old gradient dPhi/dX array allocated')
   dPhidXo = 0.0
   
   allocate( MassFlux(Nfac),stat=istat)
   call TrackMemory(istat,Nfac,'Mass fluxes array allocated')
   MassFlux = 0.0

   !
   ! linear solver section
   !
   i = 0
   allocate( Ar(Ncel),stat=istat)
   i = i + istat   
   allocate( Au(Ncel),stat=istat)
   i = i + istat
   allocate( Av(Ncel),stat=istat)
   i = i + istat
   allocate( Aw(Ncel),stat=istat)
   i = i + istat
   allocate( Su(Ncel),stat=istat)
   i = i + istat
   allocate( Sv(Ncel),stat=istat)
   i = i + istat
   allocate( Sw(Ncel),stat=istat)
   i = i + istat
   allocate( Res(Ncel),stat=istat)
   i = i + istat
   call TrackMemory(istat,Ncel*8,'Cell matrices allocated')

   Ar  = 0.0
   Au  = 0.0
   Su  = 0.0
   Av  = 0.0
   Sv  = 0.0
   Aw  = 0.0
   Sw  = 0.0
   Res = 0.0
      
   NNZ = Ncel + 2*Nint   
   i = 0
   allocate( Acoo(NNZ),stat=istat)     ! double prec.
   i = i + istat
   allocate( Arow(NNZ),stat=istat)
   i = i + istat
   allocate( Acol(NNZ),stat=istat)
   i = i + istat
   allocate( Acsr(NNZ),stat=istat)     ! double prec.
   i = i + istat
   allocate( Arwc(Ncel+1),stat=istat)
   i = i + istat
   allocate( Aclc(NNZ),stat=istat)
   istat = i + istat

   allocate( RHS(Ncel),stat=istat)     ! righthand side in dble. prec.
   istat = i + istat
   allocate( SOL(Ncel),stat=istat)     ! solution in dble. prec.
   istat = i + istat

   call TrackMemory(istat,(NNZ)*10,'Matrix A in COO/CSR format allocated')
   if( debug > 2 ) write(*,*)'NNZ: ',NNZ

   !
   ! Preconditioner arrays
   !
   ! allocate( ALU(NNZ),stat=istat)
   ! istat = i + istat
   ! allocate( JAU(NNZ),stat=istat)
   ! istat = i + istat
   ! allocate( JU(NNZ),stat=istat)
   ! istat = i + istat

   !
   ! store arrays for transient simulations
   !
   if( Transient )then
     i = 0
     allocate( Uold(Ncel+Nbnd),stat=istat)
     i = i + istat 
     allocate( Vold(Ncel+Nbnd),stat=istat )
     i = i + istat 
     allocate( Wold(Ncel+Nbnd),stat=istat )
     istat = i + istat 
     call TrackMemory(istat,(Ncel+Nbnd)*3,&
                     'Old velocity components arrays allocated')

     Uold = Guess(VarU)
     Vold = Guess(VarV)
     Wold = Guess(VarW)

     if( QuadTime )then
       i = 0
       allocate( Uold2(Ncel+Nbnd),stat=istat)
       i = i + istat 
       allocate( Vold2(Ncel+Nbnd),stat=istat )
       i = i + istat 
       allocate( Wold2(Ncel+Nbnd),stat=istat )
       istat = i + istat 
       call TrackMemory(istat,(Ncel+Nbnd)*3,&
                       'Old2 velocity components arrays allocated')

       Uold2 = Guess(VarU)
       Vold2 = Guess(VarV)
       Wold2 = Guess(VarW)
     endif

     if( SolveEnthalpy )then
       allocate( Told(Ncel+Nbnd),stat=istat)
       call TrackMemory(istat,(Ncel+Nbnd)*3,&
                       'Old enthalpy array allocated')
       Told = Guess(VarT)
       
       if( QuadTime )then
         allocate( Told2(Ncel+Nbnd),stat=istat)
         call TrackMemory(istat,(Ncel+Nbnd)*3,&
                       'Old2 enthalpy array allocated')
         Told2 = Guess(VarT)
       endif
     endif
     
     if( SolveTurb )then
       i = 0
       allocate( TEold(Ncel+Nbnd),stat=istat)
       i = i + istat 
       allocate( EDold(Ncel+Nbnd),stat=istat )
       istat = i + istat 
       call TrackMemory(istat,(Ncel+Nbnd)*3,&
                     'Old turbulence model arrays allocated')
       TEold = Guess(VarTE)
       EDold = Guess(VarED)

       if( QuadTime )then
         i = 0
         allocate( TEold2(Ncel+Nbnd),stat=istat)
         i = i + istat 
         allocate( EDold2(Ncel+Nbnd),stat=istat )
         istat = i + istat 
         call TrackMemory(istat,(Ncel+Nbnd)*3,&
                     'Old2 turbulence model arrays allocated')
         TEold2 = Guess(VarTE)
         EDold2 = Guess(VarED)
       endif

     endif
     
     if( SolveScalars )then

       allocate( SCold(Ncel+Nbnd,Nscal),stat=istat)
       call TrackMemory(istat,(Ncel+Nbnd)*Nscal,&
                       'Old scalar array(s) allocated')       
       SCold = 0.0
       
       if( QuadTime )then
         allocate( SCold2(Ncel+Nbnd,Nscal),stat=istat)
         call TrackMemory(istat,(Ncel+Nbnd)*Nscal,&
                       'Old2 scalar array(s) allocated')
         SCold2 = 0.0
       endif
     endif
!   else
!     
!     allocate( TEold(1),stat=istat)     ! temp. workaround of g95-bug!!!
!     allocate( EDold(1),stat=istat)     ! temp. workaround of g95-bug!!!
!     allocate( TEold2(1),stat=istat)    ! temp. workaround of g95-bug!!!
!     allocate( EDold2(1),stat=istat)    ! temp. workaround of g95-bug!!!
!          
   endif

end subroutine InitialiseVariables
subroutine SetBoundaryConditions
!========================================================================
!
!   this routine calculates the convective and diffusive fluxes  
!   through the cell faces.
!

   use constants
   use geometry
   use variables
   use scalars

   real Xn(3)

   if( Debug > 3 ) write(*,*)'*** SetBoundaryConditions'
   !
   ! first collect areas and inlet fluxes
   !
   areain  = 0.0
   fluxin  = 0.0
   areaout = 0.0
   icntin  = 0
   icntout = 0
   do ib=1,Nbnd
     
     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == RInlet )then
       icntin = icntin + 1
       areain = areain + Face(i)%area
       if( Reg(ir)%User )then

         ip = Face(i)%cell1
         
         Utmp  = Reg(ir)%uvw(1)
         Vtmp  = Reg(ir)%uvw(2)
         Wtmp  = Reg(ir)%uvw(3)
         TEtmp = Reg(ir)%k
         EDtmp = Reg(ir)%e
         Ttmp  = Reg(ir)%T
         Dtmp  = Reg(ir)%den
         
         Xtmp  = Face(i)%x(1)
         Ytmp  = Face(i)%x(2)
         Ztmp  = Face(i)%x(3)
         Atmp  = Face(i)%area
 
         call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                        Pref,Tref,Iter,Time,                  &
                        Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                        Reg(ir)%name1,Reg(ir)%name2)
         
         fluxin = fluxin +   Dtmp * ( Utmp*Face(i)%n(1)   &
                + Vtmp*Face(i)%n(2) + Wtmp*Face(i)%n(3) )

       else
         fluxin = fluxin + Reg(ir)%den * dot_product( Reg(ir)%uvw , Face(i)%n )
       endif
     else if( it == ROutlet )then
       icntout = icntout + 1
       areaout = areaout + Face(i)%area
     endif
   end do
   !
   ! check if valid
   !
   if( icntin > 0 .and. icntout == 0 )then
     write(*,*) 'Invalid case: inflow without outflow. Stop'
     stop
   endif
   if( icntin == 0 .and. icntout > 0 )then
     write(*,*) 'Invalid case: outflow without inflow. Stop'
     stop
   endif
   if( icntin == 0 .and. icntout == 0 )then
     write(*,*) 'Case without in and outflow.'
   endif
   
   if( icntin /= 0 )then
     rate = -1.0 * fluxin/(areain+Small)        ! mf < 0 bij inlet
     if( Debug > 1 )then
       write(*,*)'Mass flow in:',fluxin,fluxin/areain
       write(*,*)'     Area in:',areain
       write(*,*)'    Area out:',areaout
     endif
     write(IOdbg,*)'Mass flow in:',fluxin,fluxin/areain
     write(IOdbg,*)'     Area in:',areain
     write(IOdbg,*)'    Area out:',areaout
     
     rate = areain / (areaout+Small) * rate     ! corr. voor opp. verschil
     if( Debug > 1 ) write(*,*)' Rate in/out:',rate
     write(IOdbg,*)' Rate in/out:',rate
     if( rate <= 0.0 )then
       write(*,*)''
       write(*,*)'+++ Warning: Check inlet conditions'
       write(*,*)''
     endif
   else
     write(*,*)'Mass flow in: absent'
     write(IOdbg,*)'Mass flow in: absent'
     rate = 0.0
   endif
   
   do ib=1,Nbnd
     
     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in /= 0 )write(*,*)'+ internal error: corrupt bc ',ib,i,ip,in
          
     if( it == RInlet )then
       U(Ncel+ib)       = Reg(ir)%uvw(1) 
       V(Ncel+ib)       = Reg(ir)%uvw(2) 
       W(Ncel+ib)       = Reg(ir)%uvw(3) 
       P(Ncel+ib)       = P(ip)
       Den(Ncel+ib)     = Reg(ir)%den
       if( SolveEnthalpy) T(Ncel+ib) = Reg(ir)%T

       if( SolveTurbEnergy ) TE(Ncel+ib) = Reg(ir)%k   
       if( SolveTurbDiss )   ED(Ncel+ib) = Reg(ir)%e    
       if( SolveTurb )then
         VisEff(Ncel+ib) = VisLam + &
           Den(Ncel+ib)*TMCmu*TE(Ncel+ib)/(ED(Ncel+ib)+Small)
       else
         VisEff(Ncel+ib) = VisLam
       endif

       if( Restart == 0 )then
         U(ip)  = U(Ncel+ib)
         V(ip)  = V(Ncel+ib)                     ! initial guess
         W(ip)  = W(Ncel+ib)
       endif
       !
       ! extra
       !
       MassFlux(i) = Reg(ir)%den * dot_product( Reg(ir)%uvw , Face(i)%n )
       !
       if( Reg(ir)%user )then
       
         ip = Face(i)%cell1
         
         Utmp  = Reg(ir)%uvw(1)
         Vtmp  = Reg(ir)%uvw(2)
         Wtmp  = Reg(ir)%uvw(3)
         TEtmp = Reg(ir)%k
         EDtmp = Reg(ir)%e
         Ttmp  = Reg(ir)%T
         Dtmp  = Reg(ir)%den
         
         Xtmp  = Face(i)%x(1)
         Ytmp  = Face(i)%x(2)
         Ztmp  = Face(i)%x(3)
         Atmp  = Face(i)%area
         
         call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                        Pref,Tref,Iter,Time,                  &
                        Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                        Reg(ir)%name1,Reg(ir)%name2)
         
         MassFlux(i) = Dtmp * ( Utmp*Face(i)%n(1)   &
                              + Vtmp*Face(i)%n(2)   &
                              + Wtmp*Face(i)%n(3) )

         U(Ncel+ib)       = Utmp
         V(Ncel+ib)       = Vtmp
         W(Ncel+ib)       = Wtmp
         P(Ncel+ib)       = P(ip)
         if( SolveTurbEnergy ) TE(Ncel+ib) = TEtmp
         if( SolveTurbDiss )   ED(Ncel+ib) = EDtmp
         if( SolveEnthalpy)    T(Ncel+ib) = Ttmp
         Den(Ncel+ib)     = Dtmp
         
         if( SolveTurb )then
           VisEff(Ncel+ib) = VisLam + Dtmp*TMCmu*TEtmp/(EDtmp+Small)
         else
           VisEff(Ncel+ib) = VisLam
         endif

         if( Restart == 0 )then
           U(ip)  = U(Ncel+ib)
           V(ip)  = V(Ncel+ib)                  ! initial guess
           W(ip)  = W(Ncel+ib)
         endif
            
       endif
       
     else if( it == ROutlet .and. (Restart == 0) )then
     
       Xn     = Face(i)%n                      ! de surface vector
       call normalise(Xn)                      ! normaal vector met lengte 1
       
       U(Ncel+ib) = Xn(1) * rate / DensRef !
       V(Ncel+ib) = Xn(2) * rate / DensRef !
       W(Ncel+ib) = Xn(3) * rate / DensRef !

       U(ip)  = U(Ncel+ib)
       V(ip)  = V(Ncel+ib)                     ! initial guess
       W(ip)  = W(Ncel+ib)

       P(Ncel+ib)       = P(ip)
       Den(Ncel+ib)     = Den(ip)
       if( SolveEnthalpy) T(Ncel+ib) = T(ip)
       VisEff(Ncel+ib)  = VisEff(ip)

     else if( it == RSymp .and. (Restart == 0 ) )then
       U(Ncel+ib)       = U(ip)
       V(Ncel+ib)       = V(ip)
       W(Ncel+ib)       = W(ip)
       P(Ncel+ib)       = P(ip)
       Den(Ncel+ib)     = Den(ip)
       if( SolveEnthalpy) T(Ncel+ib) = T(ip)
       VisEff(Ncel+ib)  = VisEff(ip)
     else if( it == RWall .and. (Restart == 0 ) )then
       U(Ncel+ib)       = Reg(ir)%uvw(1)
       V(Ncel+ib)       = Reg(ir)%uvw(2)
       W(Ncel+ib)       = Reg(ir)%uvw(3)
       P(Ncel+ib)       = P(ip)
       Den(Ncel+ib)     = Den(ip)
       if( SolveEnthalpy) T(Ncel+ib) = Reg(ir)%T
       VisEff(Ncel+ib)  = VisEff(ip)
     else if( Restart == 0 )then
       write(*,*)'SetBoundaryConditions, unknown Boundary Condition'
       write(*,*)'i,ir,it: ',i,ir,it 
     endif

   end do

   if( Debug > 3 ) write(*,*)'=== SetBoundaryConditions'

end subroutine SetBoundaryConditions
subroutine CalculateScalar(ivar,Phi,dPhidX,PhiOld,PhiOld2)
!========================================================================

   use constants
   use geometry
   use variables
   use scalars 

   !use hypre_mod, only: solve_type                        ! === Hypre ===

   integer, intent(in)          :: iVar      ! variable/scalar id 

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd)   :: PhiOld
   real, dimension(Ncel+Nbnd)   :: PhiOld2
   real, dimension(Ncel+Nbnd,3) :: dPhidX   
   
   logical :: Warning = .false.

itmpdebug = Debug
!Debug = 0

   if( Debug > 3 ) write(*,*)'*** CalculateScalar: ',Variable(ivar)
   !
   ! calculate gradients of phi
   !
!debug=8
   call GradientPhi(iVar,Phi,dPhidX)
!debug=0

   Au  = 0.0  ! 
   Su  = 0.0  ! use the arrays for the u-velocity component 
   Acoo= 0.0  !

   !
   ! calculate fluxes through all inner faces
   !
   call FluxScalar(ivar,Phi,dPhidX)

   !
   ! unsteady term
   !
   if( Transient )then 
     rdt = 1.0 / dt
     if( Debug > 2 ) write(*,*)'Instat. Scalar term, rdt =',rdt
     if( Euler )then 
       do i=1,Ncel
         f = Den(i) * Cell(i)%vol * rdt 
         Su(i) = Su(i) + f * PhiOld(i)        
         Au(i) = Au(i) + f
       end do 
     else
       do i=1,Ncel
         f = Den(i) * Cell(i)%vol * rdt 
         Su(i) = Su(i) + f *( (1.0 + GammaTime )*PhiOld(i) - 0.5*GammaTime*PhiOld2(i) )
         Au(i) = Au(i) + f *( 1.0 + 0.5*GammaTime )
       end do 
     endif
   endif 

   !
   ! de term P*div(Vol) ?
   !

   ! 
   ! moleculaire dissipatie
   !
   if( SolveEnthalpy .and. iVar == VarT )then
     !write(*,*)'Molecular dissipation term'
     !do i=1,Ncel
     !  Dissip = GenP * Vism * vol
     !  Su(i) = Su(i) + Dissip
     !end do
     !QMdis(mat) = QMdis(mat) + Dissip(1:Ncel)
   endif
   
   !
   ! source and wall terms for turbulence models
   !
   if( iVar == VarTE .or. ivar == VarED )then

     call TurbulenceModels(iVar,Phi,dPhidX)

   endif
   
   if( UsePatches ) call PatchesSource(iVar, Phi, Au, Su)    !<====!!!!!!

   !
   ! pick underrelaxation factor
   !
   if( iVar > NVar )then
     URFScalar  = URF(VarSC)                      
   else
     URFScalar  = URF(iVar)                      
   endif

   call SetUpMatrixA( iVar,URFScalar, Au,Phi,Su  )

   !solve_type = BiCGstab2                                 ! === Hypre ===

   call SolveMatrixA( iVar,DoubleS, ResS0,ResS1, Au,Phi,Su )
 
   Residual(iVar)    = ResS1

   if( Debug > 2 )then
       write(*,*)'S:',minval(phi(:)),'< ',Variable(ivar),' <',maxval(phi(:))
       write(*,*)'Residu:',ResS0,' -> ',ResS1
   endif

   !
   ! final feasibility check?
   !
   !ic = 0
   !if( iVar == VarT )then
   !  do ip=1,Ncel
   !    if( phi(ip) < 293. .or. phi(ip) > 294. )then
   !      ic = ic + 1 
   !      phi(ip) = min(max(293.0,phi(ip)),294.0)
   !    endif
   !  end do   
   !endif
   !if( ic > 0 ) write(*,*)'Ingrepen in T:',ic

   if( iVar == VarTE .or. ivar == VarED )then
   
     Warning = .false.
     icnt    = 0
   
     do ip=1,Ncel
       if( Phi(ip) <= 0.0 )then
         icnt = icnt + 1
         if( icnt <= 10 ) write(IOdbg,*)'+ neg. value ',Phi(ip), &
                                       ' for cell ',ip,Variable(iVar)
         Warning = .true.  

         PhiAv = 0.0
         icells = 0 
         do j=1,NFaces(ip)
           k   = CFace(ip,j)
           ipp = Face(k)%cell1
           inn = Face(k)%cell2
           if( inn > 0 )then
             ! internal
             icells = icells + 1
             if( ipp == ip )then
               PhiAv = PhiAv + Phi(inn)
             elseif( inn == ip )then
               PhiAv = PhiAv + Phi(ipp)
             endif
           endif
         end do
         if( icells > 0 )then
           PhiAv = PhiAv / float(icells )
           !Phi(ip) = max(PhiAv,Small)            ! must be positive
           if( iVar == VarTE )then
             Phi(ip) = max(PhiAv,1.e-12)           ! Small is too small!        
           elseif( iVar == VarED )then
             Phi(ip) = max(PhiAv,1.e-12)           
           else
             Phi(ip) = max(PhiAv,Small)             
           endif
           if( icnt <= 10 ) write(IOdbg,*) '+ replaced by:',Phi(ip)
         else
           Phi(ip) = Small
         endif
       endif
     end do
     
     if( Warning )then
       write(IOdbg,*) '+ neg. values found for ',icnt,' cells'
       Flags(IFlagTE) = 't'
     endif
   endif

   !
   ! update sym. boundaries
   !
   do ib=1,Nbnd    
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(*,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == Rsymp .or. it == Routlet)then

       Phi(Ncel+ib) = Phi(ip)

     endif
   end do    
   !call Update_Scalars_at_boundaries2(ivar,Phi)

   if( Debug > 3 ) write(*,*)'=== CalculateScalar'

 Debug = itmpdebug   
end subroutine CalculateScalar
subroutine GradientPhiGauss(ivar,Phi,dPhidX)
!========================================================================

   use constants 
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX   

   real :: facn, facp
   real :: Xac(3), dPhidXac(3) 

   real, dimension(3) :: Xp, Xn, Xf, Xnorm, Xpac, Xnac, delp, deln

   if( Debug > 3 ) write(*,*)'*** GradientPhiGauss: ',Variable(ivar)
   write(IOdbg,*)'*** GradientPhiGauss: ',Variable(ivar),Ngradient
   
   dPhidXo = 0.0
   !
   ! iterative calculation of gradients
   !
   gradientloop: do igrad=1,Ngradient
     dPhidX  = 0.0
     faceloop: do i=1,Nfac
       ip = Face(i)%cell1
       in = Face(i)%cell2

       if( in > 0 )then 
         !
         ! internal cell face which points to two cells
         !
         facn = Face(i)%lambda
         facp = 1.0 - facn
         !
         ! for the next section see the discussion in the text
         !
         Xac      =     Cell(in)%x * facn  +   Cell(ip)%x * facp
         dPhidXac =  dPhidXo(in,:) * facn + dPhidXo(ip,:) * facp
         !
         ! now the gradient at the shifted position is known
         ! correct the value at the cell face center
         !
         PhiFace = Phi(in) * facn + Phi(ip) * facp            ! the standard 

         delta   = dot_product( dPhidXac , Face(i)%x - Xac )  ! correction
         PhiFace = PhiFace + delta
         
         !
         ! alternative method
         !
         !PhiFace = 0.5*( Phi(ip) + Phi(in) + &
         !   dot_product( dPhidXo(ip,:) , Face(i)%x - Cell(ip)%x ) + &
         !   dot_product( dPhidXo(in,:) , Face(i)%x - Cell(in)%x ) )
         !
         !
         ! dit lijkt niet te helpen...
         !
         !Xp = Cell(ip)%x
         !Xn = Cell(in)%x
         !Xf = Face(i)%x 
         !Xnorm = Face(i)%n 
         !
         !call normalise(Xnorm)
         !Xpac = Xf - dot_product(Xf-Xp,0.25*Xnorm)*Xnorm
         !Xnac = Xf - dot_product(Xf-Xn,0.25*Xnorm)*Xnorm

         !delp = Xpac - Xp
         !PhiP = Phi(ip) + dot_product( dPhidXo(ip,:) , delp )
         !deln = Xnac - Xn
         !PhiN = Phi(in) + dot_product( dPhidXo(in,:) , deln )

         !delt1 = sqrt( dot_product(( Xf  - Xpac),( Xf  - Xpac)) )
         !delt2 = sqrt( dot_product((Xnac - Xpac),(Xnac - Xpac)) )
         !flam  = delt1/delt2
         
         !PhiFace2 = PhiN * flam + PhiP * (1.0-flam)
         
         !write(*,*)'CD:',ip,in,':',PhiFace,'->',Phiface2
         !phiface = phiface2
         !
         ! now only the value at the face center is known
         ! multiply it by the area-components
         ! this is basically Gauss' theorem 
         !
         dPhidX(ip,:) = dPhidX(ip,:) + PhiFace * Face(i)%n    ! normal p -> n
         dPhidX(in,:) = dPhidX(in,:) - PhiFace * Face(i)%n    ! reversed

       else
         !
         ! boundary face
         !
         ib = Face(i)%bnd
         
         PhiFace = Phi(Ncel+ib) 
         dPhidX(ip,:) = dPhidX(ip,:) + PhiFace * Face(i)%n

       endif
               
     end do faceloop

     dPhidX(1:Ncel,1) = dPhidX(1:Ncel,1) / Cell(1:Ncel)%Vol 
     dPhidX(1:Ncel,2) = dPhidX(1:Ncel,2) / Cell(1:Ncel)%Vol 
     dPhidX(1:Ncel,3) = dPhidX(1:Ncel,3) / Cell(1:Ncel)%Vol 

     dPhidXo = dPhidX 
     !dPhidXo = dPhidXo + 0.95*( dPhidX - dPhidXo) ! under relaxation

   end do gradientloop

   if( Debug > 4 )then
     write(*,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,1)),' <  dPhi/dX  < ',maxval(dPhidX(1:Ncel,1))
     write(*,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,2)),' <  dPhi/dY  < ',maxval(dPhidX(1:Ncel,2))
     write(*,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,3)),' <  dPhi/dZ  < ',maxval(dPhidX(1:Ncel,3))
   endif
     
   if( Debug > 3 ) write(*,*)'=== GradientPhi: ',Variable(ivar)

end subroutine GradientPhiGauss
subroutine FluxScalar(ivar,Phi,dPhidX)
!========================================================================
!
!   this routine calculates the convective and diffusive scalar fluxes  
!   through the cell faces.
!
   use constants
   use geometry
   use variables
   use scalars

   integer, intent(in)          :: iVar      ! variable/scalar id 

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX 
   
   real Xac(3), Xpn(3), dPhidXac(3), dPhiC(3), dPhiD(3), dPhiP(3), dPhiN(3)
   real d(3), d1(3), d2(3)
   real Xf(3), Xp(3), Xn(3), Xnorm(3), Xpac(3), Xnac(3), delp(3), deln(3)

   if( Debug > 3 ) write(*,*)'*** FluxScalar: ',Variable(ivar),iVar

   pe0 =  Large
   pe1 = -Large

   GammaBlend = Gamma(iVar)
   IScheme    = Scheme(ivar)

   is1c = 0
   is2c = 0
   is3c = 0
   is4c = 0
   
   QtransferIn  =   0.0
   QtransferOut =   0.0
   Atot         =   0.0
   qmin         =  Large
   qmax         = -Large
   qtot         =   0.0
   hmin         =  Large
   hmax         = -Large
   htot         =   0.0

! initialize flux summations for all regions to 0 for this variable
!   do ir=1, nReg
!           Reg(ir)%rInFlux(ivar) = 0.
!           Reg(ir)%rOutFlux(ivar) = 0.
!   end do

   faceloop: do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2

     if( in > 0 )then
       !
       ! internal cell face which points to two cells
       ! see eq. 4.13 and 4.14, note lambda is from Xf to Xp
       !
       facn = Face(i)%lambda
       facp = 1.0 - facn
       !
       ! Viscosity Uncorrected! to be done?
       !
       Xac      =   Cell(in)%x * facn +   Cell(ip)%x * facp
       Phiac    =      Phi(in) * facn +      Phi(ip) * facp != central diff. scheme!
       Visac    =   VisEff(in) * facn +   VisEff(ip) * facp !<= moet alg. diff. coef worden

       !
       ! dif. coef. is dependent on the kind of simulation
       ! (turb. on or off) and the variable
       !
       ! recall: VisEff = VisLam + VisTurb
       !
       if( SolveTurb )then
         !
         ! turbulent diffusion
         !
         Visac = Visac - VisLam 
         
         if( ivar == VarT )then
           Visac = ( VisLam + Visac / Sigma_T )/Prandtl
         else if( ivar == VarTE )then
           Visac = VisLam + Visac / Sigma_k 
         else if( ivar == VarED )then
           !
           ! of is het alleen: Visac / Sigma_e ????
           ! nee want anders kloppen de extreme diss. sim. niet
           !
           Visac = VisLam + Visac / Sigma_e 
         else
           Visac = ( VisLam + Visac / Sigma_s )/Schmidt
         endif
       else
         !
         ! only laminar scalars (incl. T)
         !
         ! note: thermal dif. is defined by: lambda/Cp
         !       using the definition of the Prandtl-number
         !       Pr = (vislam Cp)/lambda one can use
         !       vislam/Pr instead.
         !
         if( ivar == VarT )then
           Visac  = Visac / Prandtl
         else
           Visac  = Visac / Schmidt
         endif
       endif
 
       dPhidXac = dPhidX(in,:) * facn + dPhidX(ip,:) * facp
       delta    = dot_product( dPhidXac , Face(i)%x - Xac )  ! correction
       
       Xpn      = Cell(in)%x - Cell(ip)%x
       VisFace  = Visac*Face(i)%area/vector_length(Xpn)

       !
       ! test stukje voor expliciete convectie schema's
       !
       ! find direction of massflux and thus PhiC and PhiD
       ! (PhiD is the downwind value of Phi)
       !
       Xf    = Face(i)%x 
       PhiP  = Phi(ip)
       PhiN  = Phi(in)
       dPhiP = dPhidX(ip,:)
       dPhiN = dPhidX(in,:)
       Xp    = Cell(ip)%x
       Xn    = Cell(in)%x
       flam  = Face(i)%lambda

       if( MassFlux(i) >= 0.0 )then
         PhiC  = PhiP
         PhiD  = PhiN
         dPhiC = dPhidX(ip,:)
         dPhiD = dPhidX(in,:)

         dPhi  = PhiD - PhiC
         d     = Xn - Xp
         ddPhi = 2.0 * dot_product( dPhiC , d )
         PhiCt = 1.0 - dphi/(ddphi+Small)

         iCelC = ip
         iCelD = in                    
         inp   = 0
       else  
         PhiC  = PhiN
         PhiD  = PhiP
         dPhiC = dPhidX(in,:)
         dPhiD = dPhidX(ip,:)

         dPhi  = PhiD - PhiC
         d     = Xp - Xn
         ddPhi = 2.0 * dot_product( dPhiC , d )
         PhiCt = 1.0 - dphi/(ddphi+Small)

         iCelC = in
         iCelD = ip
         inp = 1
       endif

       !
       ! UDS 
       !
       PhiUD    = PhiC

       !
       ! CDS (1)
       !
       delta    = dot_product( dPhidXac , Xf - Xac ) 
       PhiCD1   = Phiac + delta

       if( PhiCD1 < min(PhiC,PhiD) .or. &
           PhiCD1 > max(PhiC,PhiD) ) PhiCD1 = Phiac 
       !
       ! CDS (2)
       !
       d1       = Xf - Cell(ip)%x
       d2       = Xf - Cell(in)%x 
       del1     = dot_product( dPhidX(ip,:) , d1 )
       del2     = dot_product( dPhidX(in,:) , d2 ) 
       PhiCD2   = 0.5*( Phi(in) + Phi(ip) + del1 + del2 ) 

       if( PhiCD2 < min(PhiC,PhiD) .or. &
           PhiCD2 > max(PhiC,PhiD) ) PhiCD2 = Phiac 
       !
       ! LUDS  
       !
       if( PhiCt <= 0.0 .or. PhiCt >= 1.0 )then
         PhiLUD = PhiC
       else 
         d      = Xf - Cell(iCelC)%x
         PhiLUD = PhiC + dot_product( dPhiC , d )   

         if( PhiLUD < min(PhiC,PhiD) .or. &
             PhiLUD > max(PhiC,PhiD) ) PhiLUD = PhiC 
       endif

       !
       ! MinMod
       !
       if( PhiCt <= 0.0 .or. PhiCt >= 1.0 )then
         PhiMMD = PhiC
       elseif( PhiCt < 0.5 )then 
         PhiMMD = PhiLUD
       else
         PhiMMD = PhiCD1
       endif

       !
       ! Gamma  (met 0.1 <= Beta_m <= 0.5)
       !
       Beta_m = 0.1                  
       if( PhiCt <= 0.0 .or. PhiCt >= 1.0 )then
         PhiGAM = PhiUD
       elseif( Beta_m <= PhiCt .and. PhiCt < 1.0 )then
         PhiGAM = PhiCD1
       elseif( 0.0 < PhiCt .and. PhiCt < Beta_m )then
         fg = PhiCt / Beta_m
         PhiGAM = fg*PhiCD1 + (1.0-fg)*PhiUD
       endif

       !
       ! select scheme
       !

       if( IScheme == DSud )then
         PhiFace = PhiUD
         is1c = is1c + 1
       elseif( IScheme == DScd1 )then
         PhiFace = PhiCD1
         is2c = is2c + 1
       elseif( IScheme == DScd2 )then         
         PhiFace = PhiCD2
         is2c = is2c + 1
       elseif( IScheme == DSlud )then         
         PhiFace = PhiLUD
         is4c = is4c + 1
       elseif( IScheme == DSmod )then         
         PhiFace = PhiMMD
       elseif( IScheme == DSgam )then         
         PhiFace = PhiGAM
         is3c = is3c + 1
       else
         write(*,*)'OEPS'
         PhiFace = PhiUD
       endif      
       
       if( PhiFace < min(PhiC,PhiD) .or. PhiFace > max(PhiC,PhiD) )then
  !      write(IOdbg,*)
  !      write(IOdbg,*)'v ',Variable(ivar),IScheme,GammaBlend
  !      write(IOdbg,*)'c:',ip,'(ip)  ',in,'(in)'
  !      write(IOdbg,*)'#:',PhiC,PhiFace,PhiD
  !      write(IOdbg,*)'#:',dPhi,PhiCt,PhiUD
  !      write(IOdbg,*)'#:',PhiCD1,PhiCD2,PhiLUD,PhiGAM
  !      write(IOdbg,*)'#:',PhiAc,delta
  !      write(IOdbg,*)'> ',Xf
  !      write(IOdbg,*)'> ',Xac 
  !      write(IOdbg,*)'> ',Xf - Xac 
  !      write(IOdbg,*)'g ',dPhidXac
         
         PhiFace = PhiUD
       endif     
       
       !
       ! test stukje voor expliciete gecorrigeerde diffusie
       ! Jasaks over-relaxed approach (thesis). nope nog niet ok
       !
       d    = Cell(in)%x - Cell(ip)%x
       dXpn = vector_length( d )
       S2   = Face(i)%area**2
              
       d1   = d *S2/dot_product( d , Face(i)%n )  ! orthogonaal deel       
       d2   = Face(i)%n - d1                      ! niet-orthogonaal deel
        
       fde2 = Visac*(vector_length(d1)*( PhiN - PhiP )/dXpn + &
              dot_product( d2 , dPhidXac ))
       
       !
       ! explicit higher order fluxes 
       !
       fce = MassFlux(i) * PhiFace
       
       fde1 = Visac * dot_product( dPhidXac , Face(i)%n )  
       fde  = fde1
       
!###############
!  VISFACE = 0.0  !######## OPLETTEN!!!!
!  FDE = 0.0
!###############       
       !
       ! implicit lower order (simple upwind) 
       ! convective and diffusive fluxes
       !
       fci = min( MassFlux(i) , 0.0 ) * Phi(in) &
           + max( MassFlux(i) , 0.0 ) * Phi(ip) 

       fdi = VisFace * dot_product( dPhidXac , Xpn )

       !
       ! convective coefficients with deferred correction with
       ! gamma as the blending factor (0.0 <= gamma <= 1.0)
       !
       !      low            high    low  OLD
       ! F = F    + gamma ( F     - F    )
       !     ----   -------------------------
       !      |                  |
       !  implicit           explicit (dump into source term)
       !
       Rface(i,1) = -VisFace - max( MassFlux(i) , 0.0 )  ! Ap, Ae
       Rface(i,2) = -VisFace + min( MassFlux(i) , 0.0 )  ! An, Aw 

       blend = GammaBlend*( fce - fci )
        
       !
       ! assemble the two source terms
       !
       Su(ip) = Su(ip) - blend + fde - fdi
       Su(in) = Su(in) + blend - fde + fdi     

       peclet = MassFlux(i)/Face(i)%area*vector_length(Xpn)/(Visac+Small)

       pe0 = min( pe0 , peclet )
       pe1 = max( pe1 , peclet )
       
     else
       !
       ! boundary face
       !
       ib = Face(i)%bnd
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       ip = Face(i)%cell1
       
       if( it == RInlet )then

         dPhidXac = dPhidX(ip,:)              
         Xac      = Face(i)%x
         
         if( Reg(ir)%user )then
           Utmp  = Reg(ir)%uvw(1)
           Vtmp  = Reg(ir)%uvw(2)
           Wtmp  = Reg(ir)%uvw(3)
           TEtmp = Reg(ir)%k
           EDtmp = Reg(ir)%e
           Ttmp  = Reg(ir)%T
           Dtmp  = Reg(ir)%den

           Xtmp  = Face(i)%x(1)
           Ytmp  = Face(i)%x(2)
           Ztmp  = Face(i)%x(3)
           Atmp  = Face(i)%area

           call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                          Pref,Tref,Iter,Time,                  &
                          Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                          Reg(ir)%name1,Reg(ir)%name2)

           if( ivar == VarT )then
             PhiFace  = Ttmp 
           else if( ivar == VarTE )then
             PhiFace  = TEtmp     
           else if( ivar == VarED )then
             PhiFace  = EDtmp      
           !else
           !  PhiFace  = ScReg(ir,(iVar-Nvar))%value 
           endif
         else
           if( ivar == VarT )then
             PhiFace  = Reg(ir)%T      
           else if( ivar == VarTE )then
             PhiFace  = Reg(ir)%k      
           else if( ivar == VarED )then
             PhiFace  = Reg(ir)%e      
           else if( ivar > NVar )then
             PhiFace  = ScReg(ir,(iVar-Nvar))%value
             !write(*,*)'FluxScalar: ScReg',ir,iVar-Nvar,PhiFace
           else
             write(*,*)'+++ internal error FluxScalar: inlet for ',&
                       Variable(ivar),ivar   
           endif
         endif

         Visac    = VisEff(Ncel+ib) 

         if( SolveTurb )then
           !
           ! turbulent diffusion
           !
           Visac = Visac - VisLam 

           if( ivar == VarT )then
             Visac = ( VisLam + Visac / Sigma_T )/Prandtl
           else if( ivar == VarTE )then
             Visac = VisLam + Visac / Sigma_k 
           else if( ivar == VarED )then
             Visac = VisLam + Visac / Sigma_e 
           else
             Visac = ( VisLam + Visac / Sigma_s )/Schmidt
           endif
         else
           !
           ! only laminar scalars (incl. T)
           !
           if( ivar == VarT )then
             Visac  = Visac / Prandtl
           else
             Visac  = Visac / Schmidt
           endif
         endif

         Xpn      = Xac - Cell(ip)%x
         VisFace  = Visac*Face(i)%area/vector_length(Xpn)

         fde1 = Visac * dot_product( dPhidXac , Face(i)%n )


         !
         ! expliciete diffusie nauwkeuriger???
         !
         !Xface = Face(i)%x 
         !Xnorm = Face(i)%n 
         !call normalise(Xnorm)
         !
         !Xpac = Xface - dot_product(Xface-Cell(ip)%x,Xnorm)*Xnorm
         !dXpn = vector_length(Xface-Xpac)
         !
         !delp = Xpac - Cell(ip)%x
         !PhiP = Phi(ip) + dot_product( dPhidX(ip,:) , delp )
         !
         !fde2  = Visac * ( PhiP - PhiFace )/dXpn * &
         !        dot_product( Xnorm , Face(i)%n )  

         
         fce = MassFlux(i) * PhiFace
         fde = fde1

         fci = min( MassFlux(i) , 0.0 ) * PhiFace &
             + max( MassFlux(i) , 0.0 ) * Phi(ip)   ! voor het geval van neg. inlet 

         fdi = VisFace * dot_product( dPhidXac , Xpn ) 
         
         f   = -VisFace + min( MassFlux(i) , 0.0 )
         
         Au(ip) = Au(ip) - f
         Su(ip) = Su(ip) + fde - fdi - f*PhiFace      ! TOCH EENS NALOPEN!
         Phi(Ncel+ib) = PhiFace 

         !if( ivar == VarT ) write(*,*)'inlet:',ir,f*PhiFace

       else if( it == ROutlet )then
         
         dPhidXac = dPhidX(ip,:)     
         
         Xac      = Face(i)%x
         Visac    = VisEff(ip) 
         if( SolveTurb )then
           !
           ! turbulent diffusion
           !
           Visac = Visac - VisLam 

           if( ivar == VarT )then
             Visac = ( VisLam + Visac / Sigma_T )/Prandtl
           else if( ivar == VarTE )then
             Visac = VisLam + Visac / Sigma_k 
           else if( ivar == VarED )then
             Visac = VisLam + Visac / Sigma_e 
           else
             Visac = ( VisLam + Visac / Sigma_s )/Schmidt
           endif
         else
           !
           ! only laminar scalars (incl. T)
           !
           if( ivar == VarT )then
             Visac  = Visac / Prandtl
           else
             Visac  = Visac / Schmidt
           endif
         endif

         Xpn      = Xac - Cell(ip)%x
         PhiFace  = Phi(ip) + dot_product( dPhidX(ip,:) , Xpn )
         VisFace  = Visac*Face(i)%area/vector_length(Xpn)

         fce = MassFlux(i) * PhiFace
         fde = Visac * dot_product( dPhidXac , Face(i)%n )

         fci = MassFlux(i) * Phi(ip)                    ! per. def naar buiten
         fdi = VisFace * dot_product( dPhidXac , Xpn ) 
         
         Su(ip) = Su(ip)  + fde - fdi 
         Phi(Ncel+ib) = PhiFace 

       else if( it == RSymp )then
         !
         ! do nothing just set PhiFace equal to Phi(ip)
         ! (te simpel maar ok nu)
         !
         Phi(Ncel+ib) = Phi(ip)    

       else if( it == RWall )then
         !write(*,*)'Scalar wand: ',SolveEnthalpy,Reg(ir)%adiab,Reg(ir)%flux
         if( SolveEnthalpy .and. ivar == VarT )then
           if( Reg(ir)%adiab )then
             !
             ! adiabatic wall, kan beter nu even simpel
             !
             Phi(Ncel+ib) = Phi(ip)    
             Bnd(ib)%q    = 0.0                     ! safety first
           else 
             ! 
             ! wall temp or heat flux
             !
             if( Reg(ir)%flux )then
               PhiFlux      = Reg(ir)%T
             else
               PhiFace      = Reg(ir)%T
               Phi(Ncel+ib) = PhiFace
             endif

             Visac  = VisEff(Ncel+ib)               ! = VisLam + VisTurb
             dn     = Bnd(ib)%distance              ! dist. wall center to P
             Resist = Reg(ir)%R                     ! extra thermal resistance
             
             if( .not. SolveTurb )then
               !
               ! laminar case; note lambda/Cp = VisLam/Pr
               !
               VisFace = VisLam / Prandtl / dn 
             
               Hcoef   = 1.0/( 1.0/VisFace + Resist*Cp(Ncel+ib) ) * Face(i)%area
             else
               !
               ! turbulent case 
               !
               ! the standard Jayatilleke extra sublayer resistance factor 
               !
               SLres = 9.24*( (Prandtl/Sigma_T)**(0.75) - 1.0 ) *  &
                      (1.0 + 0.28*exp(-0.007*Prandtl/Sigma_T) )

               Cmu25 = TMCmu**0.25

               Tplus = Sigma_T * ( Bnd(ib)%uplus + SLres )
               utau  = Cmu25 * sqrt( TE(ip) )
   
               if( Bnd(ib)%yplus < Reg(ir)%ylog  )then
                 !
                 ! in the viscous sublayer
                 !
                 VisFace = VisLam / Prandtl / dn 

                 Hcoef   = 1.0/( 1.0/VisFace + Resist*Cp(Ncel+ib) ) * Face(i)%area
               else
                 VisFace = Den(ip)*Utau/(Tplus + Small)  
                 
                 Hcoef   = 1.0/( 1.0/VisFace + Resist*Cp(Ncel+ib) ) * Face(i)%area
               endif               
             endif

             if( Reg(ir)%flux )then
               PhiFace = Phi(ip) + PhiFlux / (Hcoef*Cp(Ncel+ib)/Face(i)%area)
               Phi(Ncel+ib) = PhiFace
             endif

             Au(ip) = Au(ip) + Hcoef
             Su(ip) = Su(ip) + Hcoef * PhiFace

             Tdif      = ( PhiFace - Phi(ip) )     ! temperature difference
             Bnd(ib)%h = Hcoef * Cp(Ncel+ib) / &   ! take cp out of Hcoef
                                 Face(i)%area 
             Bnd(ib)%q = bnd(ib)%h  * Tdif         ! flux q in W/m2 (incl. R!)
             Bnd(ib)%h = Visface * Cp(Ncel+ib)     ! take out cp and resistance
             Bnd(ib)%T = PhiFace                   ! just record the wall temp.
               
             if( Bnd(ib)%q > 0.0 )then
               QtransferIn  = QtransferIn  + Bnd(ib)%q * Face(i)%area
             else
               QtransferOut = QtransferOut + Bnd(ib)%q * Face(i)%area
             endif

             qmin = min(qmin,Bnd(ib)%q)
             qmax = max(qmax,Bnd(ib)%q)
             hmin = min(hmin,Bnd(ib)%h)
             hmax = max(hmax,Bnd(ib)%h)

             Atot = Atot + Face(i)%area 
             qtot = qtot + Bnd(ib)%q * Face(i)%area
             htot = htot + Bnd(ib)%h * Face(i)%area

           endif
     
         elseif( SolveTurb )then
           !
           ! doe even niets
           !
         elseif( SolveScalars .and. iVar > Nvar )then
           ! 
           ! 'adiabatic' scalar only now
           !
           Phi(Ncel+ib) = Phi(ip)    

         else
         
           write(*,*)'+ internal error FluxScalar (KAN NOG NIET)'
           write(*,*)'  Variable:',ivar,Variable(ivar)

         endif
       endif

     endif
          
   end do faceloop 

   Qtransfer = QtransferIn + QtransferOut

   if( is1c+is2c+is3c+isc4 > 0 )then
     f = float(is1c+is2c+is3c+isc4)*100.
     write(IOdbg,*)'Schemes:  UD ',is1c,float(is1c)/f
     write(IOdbg,*)'Schemes:  CD ',is2c,float(is2c)/f
     write(IOdbg,*)'Schemes: gam ',is3c,float(is3c)/f
     write(IOdbg,*)'Schemes: LUD ',is4c,float(is4c)/f
   endif
   
   if( Debug > 0 .and. ivar == VarT ) &
     write(*,*)'Total heat transferred:',Qtransfer,' W'

   if( ivar == VarT )then 
     write(IOdbg,*)'Heat transferred IN    :',QtransferIN,' W'
     write(IOdbg,*)'Heat transferred OUT   :',QtransferOUT,' W'
     write(IOdbg,*)'Total heat transferred :',Qtransfer,' W'
     if( atot > 0.0 )then
       write(IOdbg,*) hmin,'<  h  <',hmax,' '
       write(IOdbg,*) qmin,'< q_w <',qmax,' W/m2'
       write(IOdbg,*)'Total heat area        :',atot,' m2'
       write(IOdbg,*)'Average heat flux      :',qtot/atot,' W/m2'
       write(IOdbg,*)'Average heat tran.coef.:',htot/atot,' W/m2K'
     endif
   endif

   !if( Debug > 0 )then
   !    write(*,*) 'Scl: ',pe0,' < Peclet < ',pe1,'(',Variable(iVAR),')'   
       write(IOdbg,*) 'Scl: ',pe0,' < Peclet < ',pe1,'(',Variable(iVAR),')'   
   !endif
     
   if( Debug > 3 ) write(*,*)'=== FluxScalar'   

end subroutine FluxScalar
subroutine CalculateUVW 
!========================================================================

   use constants
   use geometry
   use variables

   !use hypre_mod, only: solve_type                        ! === Hypre ===

   logical :: copied = .false.   ! is the Ar array copied?
   
   if( Debug > 3 ) write(*,*)'*** CalculateUVW'

   copied = .false.

   Ar  = 0.0      ! reciprocal of one of the A-matrices
   
   Au  = 0.0      ! A-matrix for U-component
   Su  = 0.0      ! source term

   Av  = 0.0      ! V-component
   Sv  = 0.0  

   Aw  = 0.0      ! W-component
   Sw  = 0.0  

   ResU0 = 0.0
   ResU1 = 0.0
   ResV0 = 0.0
   ResV1 = 0.0
   ResW0 = 0.0
   ResW1 = 0.0
   
   !
   ! calculate fluxes through all inner faces
   ! (including processing of bc's)
   !
   call FluxUVW

   !
   ! bouyancy 
   ! (later ook = (Den(i)-DenRef(imat))*Cell(i)%vol)
   ! algemeen: bodyforce(1-3,i)*Cell(i)%vol
   !
   if( SolveEnthalpy )then
     if( Debug > 2 )write(*,*)'BodyForce',Gravity
     do i=1,Ncel
       BodyForce  = -Beta*Den(i)*Cell(i)%vol*( T(i) - Tref )

       Su(i) = Su(i) + Gravity(1) * BodyForce
       Sv(i) = Sv(i) + Gravity(2) * BodyForce
       Sw(i) = Sw(i) + Gravity(3) * BodyForce
     end do
   endif

   if( UsePatches )then
    call PatchesSource(VarU, U, Au, Su)      
    call PatchesSource(VarV, V, Av, Sv)                      !<====!!!!!!
    call PatchesSource(VarW, W, Aw, Sw)     
   endif
   !
   ! pressure force
   !
   if( Debug > 2 )write(*,*)'PressureForce'
   do i=1,Ncel
     Su(i) = Su(i) - dPdX(i,1)*Cell(i)%vol
     Sv(i) = Sv(i) - dPdX(i,2)*Cell(i)%vol
     Sw(i) = Sw(i) - dPdX(i,3)*Cell(i)%vol
   end do
   !
   ! unsteady term
   !
   if( Transient )then 
     rdt = 1.0 / dt
     if( Debug > 2 ) write(*,*)'CalculateUVW: Instat. term, rdt =',rdt
     if( Euler )then
       ! 
       ! see section 6.3.2 Implicit Euler Method
       ! 
       do i=1,Ncel
         f = Den(i) * Cell(i)%vol * rdt 
       
         Su(i) = Su(i) + f * Uold(i) 
         Sv(i) = Sv(i) + f * Vold(i)   ! op een later tijdstip  
         Sw(i) = Sw(i) + f * Wold(i)   ! ook 2 oude tijden
       
         Au(i) = Au(i) + f
         Av(i) = Av(i) + f
         Aw(i) = Aw(i) + f
       end do 
     else
       !
       ! see chapter 6.3.2 the Three Time Level Method
       !
       do i=1,Ncel
         f = Den(i) * Cell(i)%vol * rdt 
       
         Su(i) = Su(i) + f *( (1.0 + GammaTime )*Uold(i) - 0.5*GammaTime*Uold2(i) )
         Sv(i) = Sv(i) + f *( (1.0 + GammaTime )*Vold(i) - 0.5*GammaTime*Vold2(i) )
         Sw(i) = Sw(i) + f *( (1.0 + GammaTime )*Wold(i) - 0.5*GammaTime*Wold2(i) )
       
         Au(i) = Au(i) + f *( 1.0 + 0.5*GammaTime )
         Av(i) = Av(i) + f *( 1.0 + 0.5*GammaTime )
         Aw(i) = Aw(i) + f *( 1.0 + 0.5*GammaTime )
       end do 
     endif
   endif 
   
   ! ====================
   ! Velocity component U
   ! ====================
   if( SolveU )then    

     !
     ! set under relaxation factor
     !   
     URFactor = URF(VarU)
     
     !
     ! According to Ferziger/Peric chap. 8.8
     ! "since A_P^{u_i} is the same for all vel. comp. in a
     ! given CV (exept near some boundaries), one can replace 
     ! A_P^{v_n} by this quantity."
     !
     ! see also 8.10.3 (wall /symm. treatment)

     call SetUpMatrixA( VarU,URFactor, Au,U,Su  )
     
     !
     ! later anders ivm solid/fluid en poreuze media  ?!?!
     !
     if( .not. copied )then
       where( Au /= 0.0 )
         Ar = 1.0 / Au 
       elsewhere
         Ar = 0.0
       end where
       copied = .true.
       if( Debug > 3 ) write(*,*)'Calculate UVW (U): Ar based on Au'
     endif

     !solve_type = GMRES                                   ! === Hypre ===

     call SolveMatrixA( VarU,DoubleU, ResU0,ResU1, Au,U,Su )

   endif
   ! ====================
   ! Velocity component V
   ! ====================
   if( SolveV )then

     URFactor = URF(VarV)

     call SetUpMatrixA( VarV,URFactor,  Av,V,Sv  )

     if( .not. copied )then
       where( Av /= 0.0 )
         Ar = 1.0 / Av
       elsewhere
         Ar = 0.0
       end where
       copied = .true.
       if( Debug > 3 ) write(*,*)'Calculate UVW (V): Ar based on Av'
     endif

     !solve_type = GMRES                                   ! === Hypre ===

     call SolveMatrixA( VarV,DoubleV, ResV0,ResV1, Av,V,Sv )

   endif
   ! ====================
   ! Velocity component W
   ! ====================
   if( SolveW )then

     URFactor = URF(VarW)

     call SetUpMatrixA( VarW,URFactor,  Aw,W,Sw  )

     if( .not. copied )then
       where( Aw /= 0.0 )
         Ar = 1.0 / Aw
       elsewhere
         Ar = 0.0
       end where
       copied = .true.
       if( Debug > 3 ) write(*,*)'Calculate UVW (W): Ar based on Aw'
     endif
 
     !solve_type = GMRES                                   ! === Hypre ===

     call SolveMatrixA( VarW,DoubleW, ResW0,ResW1, Aw,W,Sw )

  endif

  if( Debug > 2 )then
    write(*,*)'Overzicht residuen:'
    write(*,*)'U:',ResU0,ResU1,ResU1/(ResU0+Small)
    write(*,*)'V:',ResV0,ResV1,ResV1/(ResV0+Small)
    write(*,*)'W:',ResW0,ResW1,ResW1/(ResW0+Small)
  endif

  call PrintSummary(IOdbg,iter)

  if( Debug > 3 ) write(*,*)'=== CalculateUVW'

end subroutine CalculateUVW
subroutine FluxUVW
!========================================================================
!
!   this routine calculates the convective and diffusive fluxes  
!   through the cell faces.
!

   use constants
   use geometry
   use variables

   real Xac(3), Xpn(3), ds(3),   Up(3), Ub(3), Un(3), Ut(3), &
        Force(3), TotalForce(3), Uw(3), Us(3), Xn(3), Xp(3), Xt(3)

   real d(3), d1(3), d2(3)
   real Xf(3), Xpa(3), Xna(3), dsp(3), dsn(3), dX(3)
   real Tau_nn(3), Tau_nt(3)
   real dUdXac(3), dVdXac(3), dWdXac(3)

   real dUdXP(3),  dVdXP(3),  dWdXP(3), &
        dUdXN(3),  dVdXN(3),  dWdXN(3), &
        dUdXC(3),  dVdXC(3),  dWdXC(3), &
        dUdXD(3),  dVdXD(3),  dWdXD(3), & 
        dUdXe(3),  dVdXe(3),  dWdXe(3) 

   

   if( Debug > 3 ) write(*,*)'*** FluxUVW'

   GammaBlend = Gamma(VarU)
   IScheme    = Scheme(VarU)

   if( Debug > 3 ) write(*,*)'Variable ',VarU,GammaBlend,IScheme

   pe0 =  9999.
   pe1 = -9999. 
   TotalForce = 0.0

   faceloop: do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2
     it = -1
     if( in > 0 )then
       !
       ! internal cell face which points to two cells
       !
       facn = Face(i)%lambda
       facp = 1.0 - facn
       !
       ! Viscosity Uncorrected! to be done?
       !
       Xac      =   Cell(in)%x * facn +   Cell(ip)%x * facp

       Uac      =      U(in)   * facn +      U(ip)   * facp
       Vac      =      V(in)   * facn +      V(ip)   * facp
       Wac      =      W(in)   * facn +      W(ip)   * facp

       dUdXac   =   dUdX(in,:) * facn +   dUdX(ip,:) * facp 
       dVdXac   =   dVdX(in,:) * facn +   dVdX(ip,:) * facp 
       dWdXac   =   dWdX(in,:) * facn +   dWdX(ip,:) * facp 

       Visac    =   VisEff(in) * facn +   VisEff(ip) * facp 

       Xpn      = Cell(in)%x - Cell(ip)%x
       
       VisFace  = Visac*Face(i)%area/vector_length(Xpn)

       ds       = Face(i)%x - Xac 
       dU       = dot_product( dUdXac , ds )
       dV       = dot_product( dVdXac , ds )
       dW       = dot_product( dWdXac , ds )

       UFace    = Uac + dU
       VFace    = Vac + dV
       WFace    = Wac + dW
       
       !
       ! differentie schemas
       !
       Xf    = Face(i)%x 

       UPs   = U(ip)
       VPs   = V(ip)
       WPs   = W(ip)

       UNs   = U(in)
       VNs   = V(in)
       WNs   = W(in)

       dUdXP = dUdX(ip,:)
       dVdXP = dVdX(ip,:)
       dWdXP = dWdX(ip,:)

       dUdXN = dUdX(in,:)
       dVdXN = dVdX(in,:)
       dWdXN = dWdX(in,:)

       Xp    = Cell(ip)%x
       Xn    = Cell(in)%x
       flam  = Face(i)%lambda
       
       if( MassFlux(i) >= 0.0 )then
         UC  = UPs
         VC  = VPs
         WC  = WPs
         
         UD  = UNs
         VD  = VNs
         WD  = WNs
         
         dUdXC = dUdXP
         dVdXC = dVdXP
         dWdXC = dWdXP

         dUdXD = dUdXN
         dVdXD = dVdXN
         dWdXD = dWdXN

         dU    = UD - UC
         dV    = VD - VC
         dW    = WD - WC
         
         d     = Xn - Xp
         
         ddU   = 2.0 * dot_product( dUdXC , d )
         UCt   = 1.0 - dU/(ddU+Small)

         ddV   = 2.0 * dot_product( dVdXC , d )
         VCt   = 1.0 - dV/(ddV+Small)

         ddW   = 2.0 * dot_product( dWdXC , d )
         WCt   = 1.0 - dW/(ddW+Small)

         iCelC = ip
         iCelD = in                    
         inp   = 0
       else  
         UC  = UNs
         VC  = VNs
         WC  = WNs
         
         UD  = UPs
         VD  = VPs
         WD  = WPs
         
         dUdXC = dUdXN
         dVdXC = dVdXN
         dWdXC = dWdXN

         dUdXD = dUdXP
         dVdXD = dVdXP
         dWdXD = dWdXP

         dU    = UD - UC
         dV    = VD - VC
         dW    = WD - WC
         
         d     = Xn - Xp
         
         ddU   = 2.0 * dot_product( dUdXC , d )
         UCt   = 1.0 - dU/(ddU+Small)

         ddV   = 2.0 * dot_product( dVdXC , d )
         VCt   = 1.0 - dV/(ddV+Small)

         ddW   = 2.0 * dot_product( dWdXC , d )
         WCt   = 1.0 - dW/(ddW+Small)

         iCelC = in
         iCelD = ip                    
         inp   = 1
       endif
       
       !
       ! UDS 
       !
       UUD    = UC
       VUD    = VC
       WUD    = WC
 
       !
       ! CDS (1)
       !
       ds       = Face(i)%x - Xac 
       
       deltaU   = dot_product( dUdXac , ds ) 
       deltaV   = dot_product( dVdXac , ds ) 
       deltaW   = dot_product( dWdXac , ds ) 

       UCD1     = Uac + deltaU
       VCD1     = Vac + deltaV
       WCD1     = Wac + deltaW

       !
       ! CDS (2)
       !
       d1       = Face(i)%x - Cell(ip)%x
       d2       = Face(i)%x - Cell(in)%x 

       del1     = dot_product( dUdX(ip,:) , d1 )
       del2     = dot_product( dUdX(in,:) , d2 ) 
       UCD2     = 0.5*( U(in) + U(ip) + del1 + del2 ) 

       del1     = dot_product( dVdX(ip,:) , d1 )
       del2     = dot_product( dVdX(in,:) , d2 ) 
       VCD2     = 0.5*( V(in) + V(ip) + del1 + del2 ) 

       del1     = dot_product( dWdX(ip,:) , d1 )
       del2     = dot_product( dWdX(in,:) , d2 ) 
       WCD2     = 0.5*( W(in) + W(ip) + del1 + del2 ) 

       !
       ! LUDS  
       !
       d = Face(i)%x - Cell(iCelC)%x
       
       if( UCt <= 0.0 .or. UCt >= 1.0 )then
         ULUD = UC
       else 
         ULUD = UC + dot_product( dUdXC , d )   
       endif

       if( VCt <= 0.0 .or. VCt >= 1.0 )then
         VLUD = VC
       else 
         VLUD = VC + dot_product( dVdXC , d )   
       endif

       if( WCt <= 0.0 .or. WCt >= 1.0 )then
         WLUD = WC
       else 
         WLUD = WC + dot_product( dWdXC , d )   
       endif
       
       !
       ! MinMod
       !
       if( UCt <= 0.0 .or. UCt >= 1.0 )then
         UMMD = UC
       elseif( UCt < 0.5 )then 
         UMMD = ULUD
       else
         UMMD = UCD1
       endif

       if( VCt <= 0.0 .or. VCt >= 1.0 )then
         VMMD = VC
       elseif( VCt < 0.5 )then 
         VMMD = VLUD
       else
         VMMD = VCD1
       endif
       
       if( WCt <= 0.0 .or. WCt >= 1.0 )then
         WMMD = WC
       elseif( WCt < 0.5 )then 
         WMMD = WLUD
       else
         WMMD = WCD1
       endif
      
       !
       ! Gamma
       !
       Beta_m = 0.1                   ! 0.1 <= Beta_m <= 0.5

       if( UCt <= 0.0 .or. UCt >= 1.0 )then
         UGAM = UUD
       elseif( Beta_m <= UCt .and. UCt < 1.0 )then
         UGAM = UCD1
       elseif( 0.0 < UCt .and. UCt < Beta_m )then
         fg = UCt / Beta_m
         UGAM = fg*UCD1 + (1.0-fg)*UUD
       endif

       if( VCt <= 0.0 .or. VCt >= 1.0 )then
         VGAM = VUD
       elseif( Beta_m <= VCt .and. VCt < 1.0 )then
         VGAM = VCD1
       elseif( 0.0 < VCt .and. VCt < Beta_m )then
         fg = VCt / Beta_m
         VGAM = fg*VCD1 + (1.0-fg)*VUD
       endif

       if( WCt <= 0.0 .or. WCt >= 1.0 )then
         WGAM = WUD
       elseif( Beta_m <= WCt .and. WCt < 1.0 )then
         WGAM = WCD1
       elseif( 0.0 < WCt .and. WCt < Beta_m )then
         fg = WCt / Beta_m
         WGAM = fg*WCD1 + (1.0-fg)*WUD
       endif

       !
       ! select scheme
       !
       if( IScheme == DSud )then
         UFace = UUD
         VFace = VUD
         WFace = WUD
       elseif( IScheme == DScd1 )then
         UFace = UCD1
         VFace = VCD1
         WFace = WCD1
       elseif( IScheme == DScd2 )then         
         UFace = UCD2
         VFace = VCD2
         WFace = WCD2
       elseif( IScheme == DSlud )then         
         UFace = ULUD
         VFace = VLUD
         WFace = WLUD
       elseif( IScheme == DSmod )then         
         UFace = UMMD
         VFace = VMMD
         WFace = WMMD
       elseif( IScheme == DSgam )then         
         UFace = UGAM
         VFace = VGAM
         WFace = WGAM
       else
         UFace = UUD
         VFace = VUD
         WFace = WUD
       endif      
       
       !
       ! check should be on the vector not on the
       ! components as if they were scalars
       !
       if( UFace < min(UC,UD) .or. UFace > max(UC,UD) )then
        !write(IOdbg,*)
        !write(IOdbg,*)'--  U  --'
        !write(IOdbg,*)'v ',Variable(VarU),IScheme,GammaBlend
        !write(IOdbg,*)'c:',ip,'(ip)  ',in,'(in)'
        !write(IOdbg,*)'#:',UC,UFace,UD
        !write(IOdbg,*)'#:',ddU,UCt,UUD
        !write(IOdbg,*)'#:',UCD1,UCD2,ULUD,UGAM
        !write(IOdbg,*)'#:',UAc,deltaU
        !write(IOdbg,*)'> ',Xf
        !write(IOdbg,*)'> ',Xac 
        !write(IOdbg,*)'> ',Xf - Xac 
        !write(IOdbg,*)'g ',dUdXac
         
         UFace = UUD
       endif     
       if( VFace < min(VC,VD) .or. VFace > max(VC,VD) )then
        !write(IOdbg,*)
        !write(IOdbg,*)'--  V  --'
        !write(IOdbg,*)'v ',Variable(VarV),IScheme,GammaBlend
        !write(IOdbg,*)'c:',ip,'(ip)  ',in,'(in)'
        !write(IOdbg,*)'#:',VC,VFace,VD
        !write(IOdbg,*)'#:',dDV,VCt,VUD
        !write(IOdbg,*)'#:',VCD1,VCD2,VLUD,VGAM
        !write(IOdbg,*)'#:',VAc,deltaV
        !write(IOdbg,*)'> ',Xf
        !write(IOdbg,*)'> ',Xac 
        !write(IOdbg,*)'> ',Xf - Xac 
        !write(IOdbg,*)'g ',dVDXac
         
         VFace = VUD
       endif     
       if( WFace < min(WC,WD) .or. WFace > max(WC,WD) )then
        !write(IOdbg,*)
        !write(IOdbg,*)'--  W  --'
        !write(IOdbg,*)'v ',Variable(VarW),IScheme,GammaBlend
        !write(IOdbg,*)'c:',ip,'(ip)  ',in,'(in)'
        !write(IOdbg,*)'#:',WC,WFace,WD
        !write(IOdbg,*)'#:',dDW,WCt,WUD
        !write(IOdbg,*)'#:',WCD1,WCD2,WLUD,WGAM
        !write(IOdbg,*)'#:',WAc,deltaW
        !write(IOdbg,*)'> ',Xf
        !write(IOdbg,*)'> ',Xac 
        !write(IOdbg,*)'> ',Xf - Xac 
        !write(IOdbg,*)'g ',dWdXac
         
         WFace = WUD
       endif     
       !
       ! explicit higher order convective flux (see eg. eq. 8.16)
       !
       fuce = MassFlux(i) * UFace
       fvce = MassFlux(i) * VFace
       fwce = MassFlux(i) * WFace

       sx = Face(i)%n(1)
       sy = Face(i)%n(2)
       sz = Face(i)%n(3)
       
       !
       ! explicit higher order diffusive flux based on simple uncorrected
       ! interpolated cell centred gradients(see eg. eq. 8.19)
       !
       fude1 = (dUdXac(1)+dUdXac(1))*sx + (dUdXac(2)+dVdXac(1))*sy + (dUdXac(3)+dWdXac(1))*sz
       fvde1 = (dUdXac(2)+dVdXac(1))*sx + (dVdXac(2)+dVdXac(2))*sy + (dVdXac(3)+dWdXac(2))*sz
       fwde1 = (dUdXac(3)+dWdXac(1))*sx + (dWdXac(2)+dVdXac(3))*sy + (dWdXac(3)+dWdXac(3))*sz

       !
       ! expliciete diffusieve flux/gecorrigeerde diffusie
       ! Jasaks over-relaxed approach (thesis).
       !
       d     = Cell(in)%x - Cell(ip)%x
       dXpn  = vector_length( d )
       S2    = Face(i)%area**2

       d1    = d *S2/dot_product( d , Face(i)%n )  ! orthogonaal deel       
       d2    = Face(i)%n - d1                      ! niet-orthogonaal deel

       dUdXe = vector_length(d1)*( UNs - UPs )/dXpn + &
               dot_product( d2 , dUdXac )

       dVdXe = vector_length(d1)*( VNs - VPs )/dXpn + &
               dot_product( d2 , dVdXac )

       dWdXe = vector_length(d1)*( WNs - WPs )/dXpn + &
               dot_product( d2 , dWdXac )

       fude2 = (dUdXe(1)+dUdXe(1)) + (dUdXe(2)+dVdXe(1)) + (dUdXe(3)+dWdXe(1))
       fvde2 = (dUdXe(2)+dVdXe(1)) + (dVdXe(2)+dVdXe(2)) + (dVdXe(3)+dWdXe(2))
       fwde2 = (dUdXe(3)+dWdXe(1)) + (dWdXe(2)+dVdXe(3)) + (dWdXe(3)+dWdXe(3))

       fude = Visac * fude1
       fvde = Visac * fvde1
       fwde = Visac * fwde1
       
       if( ip == 9000000 )then
         write(*,*)' '
         write(*,*)'FDE:  ',fude1,fude2
         write(*,*)'FDE:  ',fvde1,fvde2
         write(*,*)'FDE:  ',fwde1,fwde2
         write(*,*)'d :   ',d
         write(*,*)'d1:   ',d1
         write(*,*)'d2:   ',d2
         write(*,*)'S :   ',Face(i)%n
         
         write(*,*)'dudx1:',dUdXe
         write(*,*)'dudx2:',dUdXac 
         write(*,*)'dudx3:',dUdXac * Face(i)%n
         write(*,*)'dvdx1:',dvdXe
         write(*,*)'dvdx2:',dvdXac 
         write(*,*)'dvdx3:',dvdXac * Face(i)%n
         write(*,*)'dwdx1:',dwdXe
         write(*,*)'dwdx2:',dwdXac 
         write(*,*)'dwdx3:',dwdXac * Face(i)%n
       endif
       !
       ! implicit lower order (simple upwind) 
       ! convective and diffusive fluxes
       !
       fmin = min(MassFlux(i),0.0)
       fmax = max(MassFlux(i),0.0)
       
       fuci = fmin * U(in) + fmax * U(ip) 
       fvci = fmin * V(in) + fmax * V(ip) 
       fwci = fmin * W(in) + fmax * W(ip) 

       fudi = VisFace * dot_product( dUdXac , Xpn )
       fvdi = VisFace * dot_product( dVdXac , Xpn )
       fwdi = VisFace * dot_product( dWdXac , Xpn )
       !
       ! convective coefficients with deferred correction with
       ! gamma as the blending factor (0.0 <= gamma <= 1.0)
       !
       !      low            high    low  OLD
       ! F = F    + gamma ( F     - F    )
       !     ----   -------------------------
       !      |                  |
       !  implicit           explicit (dump into source term)
       !
       !            diffusion       convection
       !                v               v
       Rface(i,1) = -VisFace - max( MassFlux(i) , 0.0 )  ! P (e)
       Rface(i,2) = -VisFace + min( MassFlux(i) , 0.0 )  ! N (w)

       blend_u = GammaBlend*( fuce - fuci )
       blend_v = GammaBlend*( fvce - fvci )
       blend_w = GammaBlend*( fwce - fwci )
       !
       ! assemble the two source terms
       !
       Su(ip) = Su(ip) - blend_u + fude - fudi
       Su(in) = Su(in) + blend_u - fude + fudi     

       Sv(ip) = Sv(ip) - blend_v + fvde - fvdi
       Sv(in) = Sv(in) + blend_v - fvde + fvdi     

       Sw(ip) = Sw(ip) - blend_w + fwde - fwdi
       Sw(in) = Sw(in) + blend_w - fwde + fwdi     


       peclet = MassFlux(i)/Face(i)%area*vector_length(Xpn)/(Visac+Small)
       pe0 = min( pe0 , peclet )
       pe1 = max( pe1 , peclet )

     else
       !
       ! boundary face               Ohne Gewaehr!  kritisch blijven !!!!
       !                             ====================================
       ib = Face(i)%bnd
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       ip = Face(i)%cell1

       if( it == RInlet )then
         dUdXac   = dUdX(ip,:)
         dVdXac   = dVdX(ip,:)
         dWdXac   = dWdX(ip,:)
         
         Xac      = Face(i)%x
         if( Reg(ir)%user )then
           Utmp  = Reg(ir)%uvw(1)
           Vtmp  = Reg(ir)%uvw(2)
           Wtmp  = Reg(ir)%uvw(3)
           TEtmp = Reg(ir)%k
           EDtmp = Reg(ir)%e
           Ttmp  = Reg(ir)%T
           Dtmp  = Reg(ir)%den

           Xtmp  = Face(i)%x(1)
           Ytmp  = Face(i)%x(2)
           Ztmp  = Face(i)%x(3)
           Atmp  = Face(i)%area

           call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                          Pref,Tref,Iter,Time,                  &
                          Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                          Reg(ir)%name1,Reg(ir)%name2)

           UFace    = Utmp
           VFace    = Vtmp
           WFace    = Wtmp
            
         else
           UFace    = Reg(ir)%uvw(1) 
           VFace    = Reg(ir)%uvw(2) 
           WFace    = Reg(ir)%uvw(3) 
         endif
         
         Visac    = VisEff(Ncel+ib) 

         ! piet's stukje
         !UFace    = Face(i)%x(2)*(1.0-Face(i)%x(2))
         !write(*,*)'Uin:',Uface 
         !
         ! stuwpunt
         !UFace    =  Face(i)%x(1)
         !VFace    = -Face(i)%x(2)
         !
         ! standaard test
         !Xac(3) = 0.0
         !Umag  = vector_length( Xac )
         !Uface = Umag * 0.866025
         !Vface = Umag * 0.5
         !Xac      = Face(i)%x

         Xpn      = Xac - Cell(ip)%x
         VisFace  = Visac*Face(i)%area/vector_length(Xpn)

         fuce = MassFlux(i) * UFace
         fvce = MassFlux(i) * VFace
         fwce = MassFlux(i) * WFace

         sx = Face(i)%n(1)
         sy = Face(i)%n(2)
         sz = Face(i)%n(3)

         fude = (dUdXac(1)+dUdXac(1))*sx + (dUdXac(2)+dVdXac(1))*sy + (dUdXac(3)+dWdXac(1))*sz
         fvde = (dUdXac(2)+dVdXac(1))*sx + (dVdXac(2)+dVdXac(2))*sy + (dVdXac(3)+dWdXac(2))*sz
         fwde = (dUdXac(3)+dWdXac(1))*sx + (dWdXac(2)+dVdXac(3))*sy + (dWdXac(3)+dWdXac(3))*sz

         fude = Visac * fude
         fvde = Visac * fvde
         fwde = Visac * fwde

         fmin = min(MassFlux(i),0.0)   ! in inlaat stroomt het naar binnen mf < 0
         fmax = max(MassFlux(i),0.0)   ! fmin /= 0.0 en fmax = 0.0 !

         fuci = fmin * UFace + fmax * U(ip)     !
         fvci = fmin * VFace + fmax * V(ip)     ! voor het geval van neg. inlet ?
         fwci = fmin * WFace + fmax * W(ip)     ! nog sjekken of ok!

         fudi = VisFace * dot_product( dUdXac , Xpn )
         fvdi = VisFace * dot_product( dVdXac , Xpn )
         fwdi = VisFace * dot_product( dWdXac , Xpn )
         !
         ! by definition points a boundary normal outwards 
         ! therefore an inlet results in a mass flux < 0.0
         !
         !write(*,*)'fluxuvw: inl:',ip,massflux(i)
         
         f   = -VisFace + min( MassFlux(i) , 0.0 )

         Au(ip) = Au(ip) - f
         Su(ip) = Su(ip) - f * UFace + fude - fudi  
         U(Ncel+ib) = UFace 

         Av(ip) = Av(ip) - f
         Sv(ip) = Sv(ip) - f * VFace + fvde - fvdi 
         V(Ncel+ib) = VFace 

         Aw(ip) = Aw(ip) - f
         Sw(ip) = Sw(ip) - f * WFace + fwde - fwdi
         W(Ncel+ib) = WFace 

        !write(*,*)'inlet:',uface,vface,wface
        !write(*,*)' dUdX:',dUdXac
        !write(*,*)' dVdX:',dVdXac
        !write(*,*)' dWdX:',dWdXac
        !write(*,*)' expl:',fude,fvde,fwde
        !write(*,*)' impl:',fudi,fvdi,fwdi

       else if( it == ROutlet )then
         
         dUdXac   = dUdX(ip,:)
         dVdXac   = dVdX(ip,:)
         dWdXac   = dWdX(ip,:)
         
         Xac      = Face(i)%x
         Visac    = VisEff(ip) 

         Xpn      = Xac - Cell(ip)%x

        !UFace    = U(ip) + dot_product( dUdX(ip,:) , Xpn )    
        !VFace    = V(ip) + dot_product( dVdX(ip,:) , Xpn )  
        !WFace    = W(ip) + dot_product( dWdX(ip,:) , Xpn )  
        !
        ! modification by Harry, discard the gradient
        !
         UFace    = U(ip)    
         VFace    = V(ip)  
         WFace    = W(ip)  

         VisFace  = Visac*Face(i)%area/vector_length(Xpn)

         fuce = MassFlux(i) * UFace
         fvce = MassFlux(i) * VFace
         fwce = MassFlux(i) * WFace

         sx = Face(i)%n(1)
         sy = Face(i)%n(2)
         sz = Face(i)%n(3)

         fude = (dUdXac(1)+dUdXac(1))*sx + (dUdXac(2)+dVdXac(1))*sy + (dUdXac(3)+dWdXac(1))*sz
         fvde = (dVdXac(1)+dUdXac(2))*sx + (dVdXac(2)+dVdXac(2))*sy + (dVdXac(3)+dWdXac(2))*sz
         fwde = (dWdXac(1)+dUdXac(3))*sx + (dWdXac(2)+dVdXac(3))*sy + (dWdXac(3)+dWdXac(3))*sz

         fude = Visac * fude
         fvde = Visac * fvde
         fwde = Visac * fwde

         fmin = min(MassFlux(i),0.0)   ! in uitlaat stroomt het naar buiten mf >=0         
         fmax = max(MassFlux(i),0.0)   ! fmin = 0.0 en fmax > 0 !
                                                
         fuci = fmin * UFace + fmax * U(ip)     
         fvci = fmin * VFace + fmax * V(ip)     
         fwci = fmin * WFace + fmax * W(ip)     

         fudi = VisFace * dot_product( dUdXac , Xpn )
         fvdi = VisFace * dot_product( dVdXac , Xpn )
         fwdi = VisFace * dot_product( dWdXac , Xpn )
         !
         ! by definition points a boundary normal outwards
         ! therefore an outlet results in a mass flux >= 0.0
         !
         if( MassFlux(i) < 0.0 )then
           write(*,*)'+ internal error: neg. massflux in outlet',massflux(i)
           massflux(i) = 1.e-6
         endif
         f    = -VisFace + min( MassFlux(i) , 0.0 )

         Au(ip) = Au(ip) - f
         Su(ip) = Su(ip) - f * UFace + fude - fudi  
         U(Ncel+ib) = UFace 

         Av(ip) = Av(ip) - f
         Sv(ip) = Sv(ip) - f * VFace + fvde - fvdi 
         V(Ncel+ib) = VFace 

         Aw(ip) = Aw(ip) - f
         Sw(ip) = Sw(ip) - f * WFace + fwde - fwdi
         W(Ncel+ib) = WFace 

        ! if( massflux(i) == 0.0 )then
        !   write(*,*)'fluxuvw:',i,MassFlux(i),f
        !   write(*,*)'    uvw:',Uface,Vface,Wface
        !   write(*,*)'  uvw_P:',U(ip),V(ip),W(ip)
        !   write(*,*)'     sx:',sx,sy,sz
        ! endif   

       else if( it == RSymp )then
         !
         ! er is alleen een schuifspanning _|_ op het vlak
         ! zie hoofdstuk 8.10.4
         !          
         Xac     = Face(i)%x                 ! face centre
         Xpn     = Xac - Cell(ip)%x          ! dist. Face to P

         dUdXac  = dUdX(ip,:)                !
         dVdXac  = dVdX(ip,:)                ! use values in symm.
         dWdXac  = dWdX(ip,:)                ! cell center

         Visac   = VisEff(ip)                !

         VisFace = Visac*Face(i)%area/vector_length(Xpn)

        !
        ! modification by Harry, discard the gradient
        !
         dU      = dot_product( dUdXac , Xpn )
         dV      = dot_product( dVdXac , Xpn )
         dW      = dot_product( dWdXac , Xpn )
        
         Us(1)   = U(ip) + dU
         Us(2)   = V(ip) + dV
         Us(3)   = W(ip) + dW

        !Us(1)   = U(ip) 
        !Us(2)   = V(ip) 
        !Us(3)   = W(ip) 

         Xn = -Face(i)%n                     ! flip the normal inwards
         call normalise(Xn)                  ! normaal vector met lengte 1
         
         dp = dot_product( Us , Xn )
         Un = dp * Xn                        ! normale snelh. component

         dn     = Bnd(ib)%distance 

         Tau_nn = 2.0 * Visac * Un / dn      ! eq. (8.75)
         Force  = Tau_nn * Face(i)%area !&   ! tau_nn zorgt voor de richting
                                        !* abs(Xn) ! dus |Xn| ontbindt juist

         !if( ip == 1 )then
         !  write(*,*)' '
         !  write(*,*)'S ip xac:',ip,Xac
         !  write(*,*)'S ir xn :',ir,Xn
         !  write(*,*)'S ib un :',ib,un 
         !  write(*,*)'S i  tau:',i,tau_nn
         !  write(*,'(A,i3,5(1x1pe10.3))')' S force:',ip,Force(1:3),Bnd(ib)%yplus,dn
         !endif
         !               
         !               
         !
         if( .not. Initialisation )then               
           Au(ip) = Au(ip) + VisFace
           Av(ip) = Av(ip) + VisFace 
           Aw(ip) = Aw(ip) + VisFace

           Su(ip) = Su(ip) + VisFace*U(ip) - Force(1) 
           Sv(ip) = Sv(ip) + VisFace*V(ip) - Force(2)
           Sw(ip) = Sw(ip) + VisFace*W(ip) - Force(3)
         endif
                  
         U(Ncel+ib) = Us(1) 
         V(Ncel+ib) = Us(2) 
         W(Ncel+ib) = Us(3) 

       else if( it == RWall )then
         !
         ! schuifspanning in het vlak; hoofdstuk 8.10.3
         !
         Xac      = Face(i)%x                 ! face centre
         Uw(1)    = Reg(ir)%uvw(1)            ! 
         Uw(2)    = Reg(ir)%uvw(2)            ! wall bc's
         Uw(3)    = Reg(ir)%uvw(3)            !

         dUdXac   = dUdX(ip,:)                !
         dVdXac   = dVdX(ip,:)                ! use values in wall   
         dWdXac   = dWdX(ip,:)                ! cell center

         Visac    = VisEff(Ncel+ib)           !
         
         Xpn      = Face(i)%x - Cell(ip)%x    ! dist. wall center to P

         coef = Visac * Face(i)%area / vector_length(Xpn)

         Xn = Face(i)%n                       ! de surface vector
         call normalise(Xn)                   ! normaal vector met lengte 1

         Up(1) = U(ip)
         Up(2) = V(ip)
         Up(3) = W(ip)
         
         Up = Up - Uw                         ! snelh. verschil

         dp = dot_product( Up , Xn )
         Un = dp * Xn                         ! normale snelh. component
         Ut = Up - Un                         ! tang. snelh. component

         Uvel = abs(Ut(1)) + abs(Ut(2)) + abs(Ut(3))
         if( Uvel > Small )then               ! simpele test if |Ut| > 0

           dn     = Bnd(ib)%distance          
           Tau_nt = Visac * Ut / dn
           Force  = Tau_nt * Face(i)%area     ! tau_nt zorgt voor de richting
           
           Bnd(ib)%shear = Force              ! absolute force (in N)
           
         else
         
           Force = 0.0                        ! voor het geval dat...
           Bnd(ib)%shear = Force
         
         endif

         if( .not. Initialisation )then               
           
           TotalForce = TotalForce + Force

           !               standard
           !               implicit
           !                  V
           Au(ip) = Au(ip) + coef       !
           Av(ip) = Av(ip) + coef       ! deze 3 MOETEN gelijk 
           Aw(ip) = Aw(ip) + coef       ! zijn ivm druk iteratie
           !
           !                    corr.     expliciet
           !                  impliciet  
           !                     V           V
           Su(ip) = Su(ip) + coef*U(ip) - Force(1)
           Sv(ip) = Sv(ip) + coef*V(ip) - Force(2)
           Sw(ip) = Sw(ip) + coef*W(ip) - Force(3)

         endif
                
         U(Ncel+ib) = Uw(1)
         V(Ncel+ib) = Uw(2) 
         W(Ncel+ib) = Uw(3) 


       endif

     endif
      
     !write(*,*)'f    :',i,ip,it,V(ip)
     !write(*,*)'Force:',ip,force
     !write(*,*)'Ai   :',Au(ip),Av(ip),Aw(ip)
     !write(*,*)'Si   :',Su(ip),Sv(ip),Sw(ip)
     !write(*,*)'uvw  :',Uw(1:3)
     
   end do faceloop 

   !if( Debug > 0 )then
   !  write(*,*) 'UVW: ',pe0,' < Peclet < ',pe1   
   !  write(*,*) 'Total shear force:',TotalForce
   !endif
   write(IOdbg,*) 'UVW: ',pe0,' < Peclet < ',pe1   
   write(IOdbg,*) 'Total shear force:',TotalForce

   !write(*,*)'Au:',Au
   !write(*,*)'Su:',Su

   !write(*,*)'Av:',Av
   !write(*,*)'Sv:',Sv

   !write(*,*)'Aw:',Aw
   !write(*,*)'Sw:',Su
   
   
   if( Debug > 3 ) write(*,*)'=== FluxUVW'   

end subroutine FluxUVW
subroutine CalculatePressure(dppmax)
!========================================================================
!     This routine assembles and solves the pressure-correction
!     equation using colocated grid. SIMPLE algorithm with one
!     or more corrector steps (non-orthogonality effects taken
!     into account as described in Sect. 8.8.
!
   use constants
   use geometry
   use variables
   
   !use hypre_mod, only: new_matrix, solve_type            ! === Hypre ===
   use artificial_compressibility                         ! === AC ===

   real, dimension(3) :: PressureForce, ds, Xn

   logical, save :: ready  = .false., break = .false.
   integer, save :: idumpA = 0
   integer, save :: iprint = 0
   
   integer ipar(16)
   real    fpar(16)
   
   if( Debug > 3 ) write(*,*)'*** CalculatePressure'

   if( UseArtificialComp )then                            ! === AC ===
     if( Debug > 3 ) write(*,*)'Allocating art. compr. arrays'
     allocate( art_snd_speed(Ncel) )
     allocate( delta_t(Ncel) )
     allocate( density_change(Ncel) )
   endif

   dppmax = Large

   iloop = 0  
 1 continue
   iloop = iloop + 1
   
   Au  = 0.0         ! re-use of the U-velocity component arrays
   Su  = 0.0         !  
   MassFlux = 0.0
   PressureForce = 0.0
   
   Rface = 0.0       ! safety first; op nul zetten kan straks weg!
   !
   ! calculate fluxes through all inner faces
   ! in this routine the coef for the p'-equation are set
   ! Au and Su are set here
   !
   call FluxMass
   
   !
   ! assemble p' equations
   !
   ready = .false.
   break = .false.
   PP(:) = 0.0
   ia    = 0                     ! index for A-matrix
   
   reference = 0.0
   if( UseArtificialComp )then                            ! === AC ===

      call artificial_comp()                                 

    endif
      
   ! opendx cityplot
   !open(31,file='matrix.dx')
   
   do i=1,Ncel

     app = 0.0                   ! <= *** LET OP *** AANGEPAST
     
     do j=1,NFaces(i)
       k  = CFace(i,j)
       ip = Face(k)%cell1
       in = Face(k)%cell2
       if( in > 0 )then
       
         ! internal faces only
         if( ip == i )then
           aface   = RFace(k,2)
           ia = ia + 1
           Acoo(ia) = DBLE( aface )
           Arow(ia) = i
           Acol(ia) = in

           Su(i) = Su(i) - MassFlux(k)
          
           !write(31,'(i6,1x,i6,1x,e12.5)') i,in,aface
           
         else if( in == i )then
           aface = RFace(k,1)
           ia = ia + 1
           Acoo(ia) = DBLE( aface )
           Arow(ia) = i
           Acol(ia) = ip

           Su(i) = Su(i) + MassFlux(k)
           
           !write(31,'(i6,1x,i6,1x,e12.5)') i,ip,aface
           
         else
           write(*,*)'+ internal error: assembly A-matrix. in=',in
         endif

         app = app - aface
         
       endif
     end do

     !
     ! andere solver nog voor echt perfect symm. A-matrix?
     ! lijkt niet per se nodig; scheelt zeker in cpu-tijd
     !
     if( UseArtificialComp )then                          ! === AC ===
       Au(i)  = app + (Cell(i)%vol)/(art_snd_speed(i)**2*delta_t(i))
     else
       Au(i)  = app     
     endif
     
     ia = ia + 1
     Acoo(ia) = DBLE( Au(i) )
     Arow(ia) = i
     Acol(ia) = i

     !write(31,'(i6,1x,i6,1x,e12.5)') i,i,app

   end do
     
   if( Debug > 2 )then  
     sumsu = sum( Su(:) )
     write(*,*)'Starting residual P:',sumsu
   endif

   ! opendx dump
   !close(31)
    
   if( ia /= NNZ ) write(*,*)'+ error: Calculate P: NNZ =',ia,' =/=',NNZ
   !     
   ! solving equation system Ax=B with x=Phi
   ! even makkelijk omgooien later direct in CSR-formaat
   !
   call coocsr(Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc)

   !new_matrix = .true.                                   ! === Hypre ===
   
   ipcor = 0
   
   do while( .not. ready )

     ipcor = ipcor + 1

     if( Debug > 3 )then
       write(*,*) 'Ervoor:'
       write(*,*) minval(PP(1:Ncel)),'< PP <',maxval(PP(1:Ncel))
     endif
     
     RHS(1:Ncel) = DBLE( Su(1:Ncel) )
     SOL(1:Ncel) = DBLE( PP(1:Ncel) )

     if( Solver(VarP) == SparseKit2 )then

       ipar  =  0
       fpar  = 0.0

       ipar(2) = 0           ! 0 == no preconditioning
       ipar(3) = 0
       ipar(4) = NNZ*8       ! workspace for BGStab
       ipar(6) = 800         ! maximum number of matrix-vector multiplies

       fpar(1) = RTOL(VarP)  ! relative tolerance, must be between (0, 1)
       fpar(2) = ATOL(VarP)  ! absolute tolerance, must be positive

       write(IOdbg,*)'SolveMatrixA: SparsKit2 pressure'

       call solve_spars(debug,IOdbg,Ncel,RHS,SOL,ipar,fpar,&
                              WORK,Acsr,Aclc,Arwc)

       if( ipar(1) ==  -3 ) break = .true.

       if( ipar(1) ==  -1 ) Flags(IFlagP) = 'p'

!    elseif( UseHypre )then                               ! === Hypre ===
!    
!      solve_type = BiCGstab1
!
!      call solve_hypre(RHS,SOL)
!      new_matrix = .false.
     
     else if( Solver(VarP) == SolverUser )then

       write(*,*)'+ User solver for pressure'
       !call solve_user(debug,IOdbg,Ncel,NNZ,Acoo,Arow,Acol,S,Phi)

     else

       write(*,*)'+ internal error: solver not set'

     endif

     PP(1:Ncel) = REAL( SOL(1:Ncel) )

     if( Debug > 2 )then
       write(*,*) minval(PP(1:Ncel)),'< PP <',maxval(PP(1:Ncel))
       write(*,*) 'Using tol.:',RTOL(VarP),ATOL(VarP)
     endif
     dppmax = maxval( abs(pp) )
     if( ipcor == 1 ) dppref = dppmax
     iprint = iprint + 1

     !
     ! check if residual get's bigger
     !
     if( .not. break )then
       if( ipcor == 2 )then
         reference = fpar(6)
         if( Debug > 2 ) write(*,*)'Using Residual:',reference
       else if( ipcor > 2 )then
         if( fpar(6) > reference )then
           break = .true.
         else
           reference = fpar(6)
         if( Debug > 2 ) write(*,*)'Using new Residual:',reference
         endif
       endif
     endif
     !
     ! update P' (PP) along boundaries, needed for gradient
     !
     call Update_P_at_boundaries(VarPP,PP)
     !call Update_P_at_boundaries2(VarPP,PP)

     call GradientPhi(VarPP,PP,dPdX) 

     do i=1,Nfac
       ip = Face(i)%cell1
       in = Face(i)%cell2
       if( in > 0 )then 
    
         MassFlux(i) = MassFlux(i) + RFace(i,1)*( PP(in) - PP(ip) )

       endif      
     end do
     
     PPref = pp(IPref)             ! <= IPref voor mat. 1
     URFp  = URF(VarP)
     
     if( UseArtificialComp )then                          ! === AC ===
       do i=1,NCel
         P(i) = P(i) + URFp*( PP(i) - PPref )  
         density_change(i) = PP(i)/art_snd_speed(i)

         fact = Cell(i)%vol * Ar(i)
         if( SolveU ) U(i) = U(i) - dPdX(i,1) * fact + &
                            (density_change(i) * U(i))/Densref 
         if( SolveV ) V(i) = V(i) - dPdX(i,2) * fact + &
                            (density_change(i) * V(i))/Densref
         if( SolveW ) W(i) = W(i) - dPdX(i,3) * fact + &
                            (density_change(i) * W(i))/Densref
       end do
     else
       do i=1,NCel
         P(i) = P(i) + URFp*( PP(i) - PPref )  

         fact = Cell(i)%vol * Ar(i)
         if( SolveU ) U(i) = U(i) - dPdX(i,1) * fact
         if( SolveV ) V(i) = V(i) - dPdX(i,2) * fact
         if( SolveW ) W(i) = W(i) - dPdX(i,3) * fact

         !write(*,*)'dp:',i,dPdX(i,:)
       end do
     endif
     !
     ! een residu bouwen
     !
     Res(:) = 0.0
     do i=1,Ncel
       do j=1,NFaces(i)
         k  = CFace(i,j)
         ip = Face(k)%cell1
         in = Face(k)%cell2
         if( in > 0 )then
           ! internal faces 
           if( ip == i )then
             Res(i) = Res(i) - MassFlux(k)
           else if( in == i )then
             Res(i) = Res(i) + MassFlux(k)
           endif
         else
           ! boundary face
           Res(i) = Res(i) - MassFlux(k)
         endif
       end do
     end do
     Residual(VarP) = sum(abs(Res)) * ResiNorm(VarP)

     Su = 0.0
     P(IPref) = 0.0               ! <= let op mat(1) en in input deck!

     ! 
     ! update the source term Su for the second step
     ! 
     call FluxMass2

     PP = 0.0                     ! reset p' (=p''), is this really needed?

     sumsu = sum( Su(:) )

     if( debug > 2 )write(*,*)'P: ',ipcor,reference,Residual(VarP),sumsu,dppmax

     !
     ! stopping criteria
     !          
     if( ipcor >= MaxPCOR .or. dppmax <= FactDPP*dppref .or. break ) ready = .true.
     
   end do

   if( Debug > 1 )then
     write(*,*)'Number of pressure corrections:',ipcor,dppmax
     write(IOdbg,*)'Number of pressure corrections:',ipcor
   endif

   if( dppmax > Small )then
     dppmax = dppref
   else
     dppmax = 0.0
   endif
   !
   ! update pressure P along boundaries
   !
   !call Update_P_at_boundaries(VarP,P)
   call Update_P_at_boundaries(VarP,P)
   call GradientPhi(VarP,P,dPdX)
   
   call UpdateP(P,dPdX)

   !DXdebug(1:Ncel) = su(1:Ncel)
   !DXdebug(1:Ncel) = dpdX(1:Ncel,2)
   !call CheckMassFlux
   !call ShowCells

   if( UseArtificialComp )then                            ! === AC ===

     deallocate( art_snd_speed )
     deallocate( delta_t )
     deallocate( density_change )

   endif
!  if( UseHypre )then                                     ! === Hypre ===
!
!    call hypre_destroy()
!
!  endif
   !
   ! extract pressure force on solid walls
   !
   PressureForce = 0.0
   
   do ib=1,Nbnd           
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(*,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == Rwall)then

       ds = Face(i)%x - Cell(ip)%x
          
      !
      ! modification by Harry, discard the gradient
      !
       Pwall = P(ip) + dot_product( dPdX(ip,:) , ds ) 
      !
      !Pwall = P(ip)  

       Xn    = Face(i)%n                       ! de surface vector

       PressureForce = PressureForce + Pwall * Xn
       
     endif
   end do

   if( (PressureForce(1)+PressureForce(2)+PressureForce(3)) > Small )then
   
     write(IOdbg,*) 'Total pressure force:',PressureForce
   
   endif


   if( Debug > 3 ) write(*,*)'=== CalculatePressure'
   
end subroutine CalculatePressure
subroutine FluxMass
!========================================================================
   use constants
   use geometry
   use variables
   use scalars

   real :: facn, facp
   real :: Xac(3), Xpn(3), Xn(3), Xface(3), Xnorm(3)
   real :: Xpac(3), Xnac(3), Xpac2(3), Xnac2(3), Xnp(3), Xnn(3)
   real :: deln(3), delp(3), Xpn2(3)
   real :: dUdXac(3), dVdXac(3), dWdXac(3) ,delta(3)
   real :: Uin(3), Uout(3), Uf(3)

   real :: FlowRegion(30), FlowFact(30)   ! MOET BETER

   logical :: WarningOutlet

   if( Debug > 3 )write(*,*)'*** FluxMass'

   WarningOutlet = .false.

   icinl = 0
   icout = 0
   icsym = 0
   icwal = 0   

   faceloop: do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2
     
     if( in > 0 )then

       facn = Face(i)%lambda
       facp = 1.0 - facn
       
       dUdXac = dUdX(in,:) * facn + dUdX(ip,:) * facp
       dVdXac = dVdX(in,:) * facn + dVdX(ip,:) * facp
       dWdXac = dWdX(in,:) * facn + dWdX(ip,:) * facp

       !
       ! density here with a CD-like scheme; 
       ! to be done introduce a blend with UD (small amount, 5-10%)
       !
       Denf   = Den(in)* facn + Den(ip)* facp    
        
       Xac    = Cell(in)%x * facn + Cell(ip)%x * facp

       Xface  = Face(i)%x                        
       delta  = Xface - Xac 
       
       Uface  = U(in)*facn + U(ip)*facp + dot_product( dUdXac , delta )
       Vface  = V(in)*facn + V(ip)*facp + dot_product( dVdXac , delta )
       Wface  = W(in)*facn + W(ip)*facp + dot_product( dWdXac , delta )

       !
       ! this alternative method seems nicer, simpler, even better?
       ! disadvantage: not the same method as density (inconsistent?)
       !
       ! delta  = Face(i)%x - Cell(ip)%x
       ! UfaceP = U(ip) + dot_product( dUdX(ip,:) , delta )
       ! VfaceP = V(ip) + dot_product( dVdX(ip,:) , delta )
       ! WfaceP = W(ip) + dot_product( dWdX(ip,:) , delta )
       ! 
       ! delta  = Face(i)%x - Cell(in)%x
       ! UfaceN = U(in) + dot_product( dUdX(in,:) , delta )
       ! VfaceN = V(in) + dot_product( dVdX(in,:) , delta )
       ! WfaceN = W(in) + dot_product( dWdX(in,:) , delta )
       !
       ! Uface  = 0.5*( UfaceP + UfaceN )     !
       ! Vface  = 0.5*( VfaceP + VfaceN )     ! simple average
       ! Wface  = 0.5*( WfaceP + WfaceN )     !
       !
       
       MassFlux(i) = Denf * ( Uface*Face(i)%n(1) + &   !
                              Vface*Face(i)%n(2) + &   ! = rho.v_n
                              Wface*Face(i)%n(3) )     !

       !
       ! this was the easy part
       !

       Xnorm = Face(i)%n 
       call normalise(Xnorm)

       !
       ! find auxiliary nodes  (eq. 8.53)
       ! _    _       _   _   _  _
       ! r ,= r  - [( r - r ).n] n
       !  P    e       e   P
       !
       Xpac = Xface - dot_product(Xface-Cell(ip)%x,Xnorm)*Xnorm

       Xnp  = Xpac - Xface 

       ! _    _       _   _   _  _
       ! r ,= r  - [( r - r ).n] n
       !  E    e       e   E
       !
       Xnac = Xface - dot_product(Xface-Cell(in)%x,Xnorm)*Xnorm

       Xnn  = Xnac - Xface 

       !
       ! use the shortest
       !
       dp   = vector_length(Xnp)
       dn   = vector_length(Xnn)
       
       if( dn >= dp )then
         Xpac2 = Xface + Xnp
         Xnac2 = Xface - Xnp
         dl = 2 * dp
       else
         Xpac2 = Xface - Xnn 
         Xnac2 = Xface + Xnn
         dl = 2 * dn
       endif

       !
       ! set p ,and p , (see eqn. (8.54))
       !      P      E
       !
       delp = Xpac2 - Cell(ip)%x
       Pip  = P(ip) + dot_product( dPdX(ip,:) , delp )
       
       deln = Xnac2 - Cell(in)%x
       Pin  = P(in) + dot_product( dPdX(in,:) , deln )
       
       Xpn  = Xnac2 - Xpac2
       Xpn2 = Cell(in)%x - Cell(ip)%x

       apv1 = Den(ip)*Ar(ip)          !   * dp   ! stom! dl en delen door dl!
       apv2 = Den(in)*Ar(in)          !   * dn
       !ApV  = ( apv1 + apv2 )                     ! average
       
       fc = vector_length(Xpn)/vector_length(Xpn2)
       
       ApV  = 0.5 *( apv1 + apv2 )    !* fc                     ! average
       
       ApV  = ApV*Face(i)%area**2     !/ dl                     ! times surf. S^2

       !dpx  = 0.5*( dPdX(in,1) + dPdX(ip,1) ) * Xpn2(1) !
       !dpy  = 0.5*( dPdX(in,2) + dPdX(ip,2) ) * Xpn2(2) ! (gradP).Xpn term
       !dpz  = 0.5*( dPdX(in,3) + dPdX(ip,3) ) * Xpn2(3) !
       
       dpx  = ( dPdX(in,1) * facn + dPdX(ip,1) * facp ) * Xpn(1) !
       dpy  = ( dPdX(in,2) * facn + dPdX(ip,2) * facp ) * Xpn(2) ! (gradP).Xpn term
       dpz  = ( dPdX(in,3) * facn + dPdX(ip,3) * facp ) * Xpn(3) !

       fact = ApV 
                     
       Rface(i,1) = -fact                              ! the coefs for the
       Rface(i,2) = -fact                              ! P'-equation see 8.56


       !
       ! See eq. 8.56 FP 3rd edition
       !
       !            first part           the term in square brackets
       !
       !               m*           /delta Omega     1
       !              v          - --------------- (--)[ .. ]
       !               n           ( r   - r  ).n   Ap
       !                              E     P
       !
       MassFlux(i) = MassFlux(i) - fact*( (Pin-Pip) - dpx - dpy - dpz )

     else
       !
       ! boundary face
       !
       ib = Face(i)%bnd
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       ip = Face(i)%cell1
       
       if( it == RInlet )then
         icinl = icinl + 1
         !
         ! inlet simply has prescribed velocities
         ! get them from the input (or use user data)
         !

         Xac    = Face(i)%x
         
         if( Reg(ir)%user )then
           Utmp  = Reg(ir)%uvw(1)
           Vtmp  = Reg(ir)%uvw(2)
           Wtmp  = Reg(ir)%uvw(3)
           TEtmp = Reg(ir)%k
           EDtmp = Reg(ir)%e
           Ttmp  = Reg(ir)%T
           Dtmp  = Reg(ir)%den

           Xtmp  = Face(i)%x(1)
           Ytmp  = Face(i)%x(2)
           Ztmp  = Face(i)%x(3)
           Atmp  = Face(i)%area

           call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                          Pref,Tref,Iter,Time,                  &
                          Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                          Reg(ir)%name1,Reg(ir)%name2)

           Uin(1) = Utmp
           Uin(2) = Vtmp
           Uin(3) = Wtmp 
           
           dens   = Dtmp
         else
           Uin(1) = Reg(ir)%uvw(1) 
           Uin(2) = Reg(ir)%uvw(2) 
           Uin(3) = Reg(ir)%uvw(3) 
           dens   = Reg(ir)%den 
         endif
         
         ! piet's stukje
         !Uin(1) = Face(i)%x(2)*(1.0-Face(i)%x(2))
         !
         ! stuwpunt
         !Uin(1)    =  Face(i)%x(1)
         !Uin(2)    = -Face(i)%x(2)
         !
         ! standaard test
         !Xac(3) = 0.0
         !Umag   = vector_length( Xac )
         !Uin(1) = Umag * 0.866025
         !Uin(2) = Umag * 0.5
         !write(*,*)'Uin:',umag,Uin

         MassFlux(i) = dens * dot_product( Uin , Face(i)%n )

         su(ip) = su(ip) - MassFlux(i)    ! mf < 0 => source

         !write(*,*)'inlet:',ip,Uin,MassFlux(i)

       else if( it == ROutlet )then
         icout = icout + 1
         !
         ! outlet simply looks at node P, correct using gradient
         !
         delta  = Face(i)%x - Cell(ip)%x 

        !
        ! modification by Harry, discard the gradient
        !
        !Uface  = U(ip) + dot_product( dUdX(ip,:) , delta )
        !Vface  = V(ip) + dot_product( dVdX(ip,:) , delta )
        !Wface  = W(ip) + dot_product( dWdX(ip,:) , delta )
        !
         Uface  = U(ip) 
         Vface  = V(ip) 
         Wface  = W(ip) 

         Denf   = Den(ip)                            ! NIET GECORR.

         MassFlux(i) = Denf * ( Uface*Face(i)%n(1) + &   !
                                Vface*Face(i)%n(2) + &   ! = rho.v_n
                                Wface*Face(i)%n(3) )     !  
         !
         ! For an outlet MassFlux must be 0.0 or positive
         !
         if( MassFlux(i) < 0.0 )then
           if( .not. WarningOutlet )then
             write(IOdbg,*)'+ Warning: Inflow detected on outflow boundary'
             Flags(IFlagInflow) = 'i'
             WarningOutlet      = .true.
           endif
           MassFlux(i) = 1.e-6
           
           !
           ! to be sure reset add. variables too
           !
           if( SolveTurbEnergy ) TE(Ncel+ib) = TE(ip)
           if( SolveTurbDiss   ) ED(Ncel+ib) = ED(ip)
           if( SolveVisc       ) VisEff(Ncel+ib) = VisEff(ip)
           if( SolveEnthalpy   ) T(Ncel+ib) = T(ip)
           if( SolveScalars    ) SC(Ncel+ib,1:Nscal) = SC(ip,1:Nscal)
         endif
         !
         ! set Su below after the mass flux into the domain is known
         !
       else if( it == RSymp )then
         icsym = icsym + 1
         MassFlux(i) = 0.0
       else if( it == RWall )then
         icwal = icwal + 1
         MassFlux(i) = 0.0
       else
         write(*,*)'MassFlux: unknown bc, type ',it
       endif
     endif
   end do faceloop
   if( Debug > 2 )then
     write(*,*)'MassFlux: bc inlet :',icinl
     write(*,*)'          bc outlet:',icout
     write(*,*)'          bc symp. :',icsym
     write(*,*)'          bc wall  :',icwal
   endif
 
   if( Debug > 2 )write(*,*)'Checking global mass flux'
   i1 = count( Reg%typ == RInlet  )
   i2 = count( Reg%typ == ROutlet )    

   if( i1 > 0 .and. i2 == 0 )write(*,*)'FOUT IN BC'
   
   if( i1 > 0 )then                        ! flow with in/out
     FlowRegion = 0.0                      ! reset to zero

     flowin = 0.0
     do ib=1,Nbnd
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       if( it == RInlet )then
         i  = Bnd(ib)%face
         flowin = flowin + MassFlux(i)
       endif
     end do
     write(IOdbg,*)'Mass flow rate in: ',FlowIn

     flowout = 0.0          
     do ib=1,Nbnd                         ! flowsplit moet nog!
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       if( it == ROutlet )then
         i  = Bnd(ib)%face
         FlowOut = FlowOut + MassFlux(i)
         FlowRegion(ir) = FlowRegion(ir) + MassFlux(i)
       endif
     end do
     write(IOdbg,*)'Mass flow rate out:',FlowOut
     do ir=1,Nreg
       it = Reg(ir)%typ
       if( it == ROutlet ) &
         write(IOdbg,*)'Region ',ir,FlowRegion(ir),FlowRegion(ir)/FlowOut 
     end do
     do ir=1,Nreg
       it = Reg(ir)%typ
       if( it == ROutlet )then
         split = Reg(ir)%splvl
         FlowFact(ir) = -(split * FlowIn) / FlowRegion(ir)
         write(IOdbg,*)'Fact Region ',ir,split,FlowFact(ir),split * FlowIn
       else
         FlowFact(ir) = 0.0
       endif
     end do

     if( flowout < Small )then
       ! nothing set yet, thus an initial guess,
       ! start with calculating area
       areaout = 0.0
       do ib=1,Nbnd                        
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         if( it == ROutlet )then
           i  = Bnd(ib)%face
           areaout = areaout + Face(i)%area
         endif
       end do
       ratearea = -flowin/areaout        
       write(IOdbg,*)'Flowrate per area:',areaout,' m2',ratearea,' kg/s/m2'

       flowout = 0.0         
       do ib=1,Nbnd           
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         if( it == ROutlet )then
           i     = Bnd(ib)%face
           split = Reg(ir)%splvl
           MassFlux(i) = ratearea*Face(i)%area * split   ! flowsplit
           if( Face(i)%n(1) > Small )then
             U(Ncel+ib) = MassFlux(i)/Face(i)%n(1)
           else
             U(Ncel+ib) = 0.0
           endif
           if( Face(i)%n(2) > Small )then
             V(Ncel+ib) = MassFlux(i)/Face(i)%n(2)
           else
             V(Ncel+ib) = 0.0
           endif
           if( Face(i)%n(3) > Small )then
             W(Ncel+ib) = MassFlux(i)/Face(i)%n(3)
           else
             W(Ncel+ib) = 0.0
           endif
           flowout = flowout + MassFlux(i)
         endif
       end do
       write(IOdbg,*)'Initial flow rate out: ',flowout
     endif

     fact = -flowin /( flowout + Small )         ! amount of unbalance
     write(IOdbg,*)'Correct flow rate:',fact

     if( abs(1.0-fact) > 0.001 ) Flags(IFlagMass) = 'm'
      
     flowout2 = 0.0          
     do ib=1,Nbnd                                ! flowsplit moet nog!
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       if( it == ROutlet )then                   ! *** MOET DIT ANDERS? ***
         i     = Bnd(ib)%face
         
         MassFlux(i) = MassFlux(i) * FlowFact(ir)
         flowout2 = flowout2 + MassFlux(i)
         
         if( SolveU ) U(Ncel+ib) = U(Ncel+ib) * fact    ! 
         if( SolveV ) V(Ncel+ib) = V(Ncel+ib) * fact    ! dus toch wel?
         if( SolveW ) W(Ncel+ib) = W(Ncel+ib) * fact    !
         
         ip = Face(i)%cell1
         su(ip) = su(ip) - MassFlux(i)

         !write(*,*)'outlet:',ip,U(Ncel+ib),V(Ncel+ib),W(Ncel+ib),MassFlux(i)

       endif
     end do
     if( Debug > 2 )write(*,*)'Adapted mass flow rate out:',flowout2
     
   endif
      
   if( Debug > 3 )write(*,*)'=== FluxMass'
   
end subroutine FluxMass
subroutine FluxMass2
!========================================================================
!
! this routine should be consistent with FluxMass
! internal faces only
!
   use constants
   use geometry
   use variables

   real :: facn, facp
   real :: Xface(3), Xnorm(3), Xp(3), Xn(3), Xpn(3) ,delta(3), Xpa(3), Xna(3)
   real :: Xpac(3), Xnac(3), Xpac2(3), Xnac2(3), Xnp(3), Xnn(3)

   if( Debug > 3 )write(*,*)'*** FluxMass2'

   faceloop: do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )then

       facn = Face(i)%lambda
       facp = 1.0 - facn

       Xface = Face(i)%x
       
       Xnorm = Face(i)%n 
       call normalise(Xnorm)

       !
       ! find auxiliary nodes  (eq. 8.53)
       ! _    _       _   _   _  _
       ! r ,= r  - [( r - r ).n] n
       !  P    e       e   P
       !
       Xpac = Xface - dot_product(Xface-Cell(ip)%x,Xnorm)*Xnorm

       Xnp  = Xpac - Xface 

       ! _    _       _   _   _  _
       ! r ,= r  - [( r - r ).n] n
       !  E    e       e   E
       !
       Xnac = Xface - dot_product(Xface-Cell(in)%x,Xnorm)*Xnorm

       Xnn  = Xnac - Xface 

       !
       ! use the shortest
       !
       dp   = vector_length(Xnp)
       dn   = vector_length(Xnn)
       
       if( dn >= dp )then
         Xpac2 = Xface + Xnp
         Xnac2 = Xface - Xnp
         dl = 2 * dp
       else
         Xpac2 = Xface - Xnn 
         Xnac2 = Xface + Xnn
         dl = 2 * dn
       endif

       Xn = Xnac2 - Cell(in)%x
       Xp = Xpac2 - Cell(ip)%x

       !if( ip == 4 )then
       !  write(*,*)'*** Cel:',ip,in
       !  write(*,*)'  dp,dn:',dp,dn
       !  write(*,*)'  Xnac :',Xnac
       !  write(*,*)'  Xpac :',Xpac
       !  write(*,*)'    Xn :',Xn
       !  write(*,*)'    Xp :',Xp
       !endif
       !
       ! the coeff. is the same as in FluxMass
       !
       fact = Rface(i,1)

       !          ,    
       !   (grad p ).(r ,-r )   
       !          E    E   E
       !
       dpx = dPdX(in,1)*Xn(1) - dpdx(ip,1)*Xp(1)
       dpy = dPdX(in,2)*Xn(2) - dpdx(ip,2)*Xp(2)
       dpz = dPdX(in,3)*Xn(3) - dpdx(ip,3)*Xp(3)

       fc  = fact*( dpx + dpy + dpz ) 

       MassFlux(i) = MassFlux(i) + fc
       
       Su(ip) = Su(ip) - fc
       Su(in) = Su(in) + fc
 
     endif
   end do faceloop
  
   if( Debug > 3 )write(*,*)'=== FluxMass2'
   
end subroutine FluxMass2
subroutine UpdateP(Phi,dPhiDx)
!========================================================================
   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX   

   real :: ds(3)

   if( Debug > 3 )write(*,*)'*** UpdateP'

   do ib=1,Nbnd           
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(*,*)'*** Error: Bnd array corrupt'
     
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == ROutlet )then
       !
       ! the pressure gradient in the outlet should be 0
       !
       Phi(Ncel+ib) = Phi(ip) 
        
     else
       !
       ! modification by Harry, discard the gradient
       !
       ds = Face(i)%x - Cell(ip)%x
       Phi(Ncel+ib) = Phi(ip) + dot_product( dPhiDx(ip,:) , ds ) 
       
       !Phi(Ncel+ib) = Phi(ip) 

     endif

   end do

   if( Debug > 3 )write(*,*)'=== UpdateP'
end subroutine UpdateP
subroutine Update_P_at_boundaries(ivar,Phi)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(3) :: GradPhi, Xp, dX, ds 

   if( Debug > 3 ) write(*,*)'*** Update_P_at_boundaries',variable(ivar)
    
   do ib=1,Nbnd           
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(*,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ
          
     Xp = Cell(ip)%x 
     GradPhi = 0.0
     icnt  = 0
     icnt1 = 0
     icnt2 = 0
     icnt3 = 0
     
     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2

       if( ipn > 0 )then
         ! internal face
         if( ipp == ip )then
           dPhi = Phi(ipn) - Phi(ipp)
           dX   = Cell(ipn)%x - Cell(ipp)%x
         else
           dPhi = Phi(ipp) - Phi(ipn)
           dX   = Cell(ipp)%x - Cell(ipn)%x
         endif

         icnt = icnt + 1
         if( dX(1) /= 0.0 )then
           icnt1      = icnt1 + 1
           GradPhi(1) = GradPhi(1) + dPhi/dX(1)
         endif
         if( dX(2) /= 0.0 )then
           icnt2      = icnt2 + 1
           GradPhi(2) = GradPhi(2) + dPhi/dX(2)
         endif
         if( dX(3) /= 0.0 )then
           icnt3      = icnt3 + 1
           GradPhi(3) = GradPhi(3) + dPhi/dX(3)
         endif
                
       endif 
     end do   

     if( icnt > 0 ) GradPhi = GradPhi / float(icnt)

     ds = Face(i)%x - Xp

     !Phi(Ncel+ib) = Phi(ip) + dot_product( ds , GradPhi )

     Phi(Ncel+ib) = Phi(ip) 

   end do

   if( Debug > 3 )write(*,*)'=== Update_P_at_boundaries'
end subroutine Update_P_at_boundaries
subroutine UpdateBC
!========================================================================
   use constants
   use geometry
   use variables

   real :: ds(3)

   if( Debug > 3 )write(*,*)'*** UpdateBC (update symm. bc''s)'
   
   do ib=1,Nbnd           
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ
     if( in > 0 )write(*,*)'*** Error: Bnd array corrupt'
     
     if( it == RSymp )then
       !
       ! modification by Harry, discard the gradient
       !
       ds = Face(i)%x - Cell(ip)%x
       !   
       if( SolveU ) U(Ncel+ib) = U(ip) + dot_product( dUdX(ip,:) , ds ) 
       if( SolveV ) V(Ncel+ib) = V(ip) + dot_product( dVdX(ip,:) , ds )     
       if( SolveW ) W(Ncel+ib) = W(ip) + dot_product( dWdX(ip,:) , ds )     
       ! (P already done)

       !if( SolveU ) U(Ncel+ib) = U(ip)  
       !if( SolveV ) V(Ncel+ib) = V(ip)      
       !if( SolveW ) W(Ncel+ib) = W(ip)      
       
     endif
   end do

   if( Debug > 3 ) write(*,*)'=== UpdateBC'
end subroutine UpdateBC
subroutine ShowMass(btype)
!========================================================================

   use constants
   use geometry
   use variables

   integer btype

   sumtot = 0.0
   do ib=1,Nbnd
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ
     if( it == btype )then
       i   = Bnd(ib)%face
       ic  = Face(i)%cell1
       sum = 0.0       
       do j=1,NFaces(ic)
         k = CFace(ic,j)
         ip = Face(k)%cell1
         in = Face(k)%cell2
         if( ip == ic )then
           sum = sum - MassFlux(k)
         else
           sum = sum + MassFlux(k)
         endif
       end do
       sumtot = sumtot + sum
       !write(*,*)'Cell:',ip,sum,sumtot
     endif
   end do

end subroutine ShowMass
subroutine ShowCells
!========================================================================

   use constants
   use geometry
   use variables

   IOloc = IOdef
   if( PrintCellVar .or. PrintWallVar )then
     if( PrintFile(1:1) /= '-' )then
       call openfile(IOprt,PrintFile,'.dat','FORMATTED','SEQUENTIAL', &
                     'UNKNOWN',Debug)
       IOloc = IOprt
       write(*,*) 'Printing to ',PrintFile(1:lens(PrintFile)),'.dat'
     endif
   endif
   
   if( PrintCellVar )then
     if( PrintCellVarUser )then
     
       call UserPrintCell(IOloc)
     
     else
       write(IOloc,'(/,'' Cell data'')')
       write(IOloc,'(A,A)')'    Cell      U          V          W          P    ',&
                       '      k       epsilon    T_rel      Vis_t'

       do i=IPrintCellVarStart,min(IPrintCellVarEnd,NCel),IPrintCellVarInc

         TEtmp = 0.0
         EDtmp = 0.0
         Ttmp  = Tref

         if( SolveTurbEnergy )TEtmp = TE(i)
         if( SolveTurbDiss )  EDtmp = ED(i)
         if( SolveEnthalpy )  Ttmp  = T(i)

         write(IOloc,'(i8,8(1x,1pe10.3))') i,u(i),v(i),w(i),p(i), &
                                     TEtmp,EDtmp,Ttmp-Tref,VisEff(i)-Vislam

       end do
     endif
   endif
   
   if( PrintWallVar )then
     if( PrintWallVarUser )then
     
       call UserPrintWall(IOloc)
     
     else
       write(IOloc,'(/,'' Wall data'')')
       write(IOloc,'(A,A)')'    Cell   Reg  X-shear    Y-shear    Z-shear',&
        '      Y+         U+         dn        T_w        h_w        q_w'

       do ib=IPrintWallVarStart,min(IPrintWallVarEnd,NBnd),IPrintWallVarInc

         i  = Bnd(ib)%face
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ

         ip = Face(i)%cell1
         in = Face(i)%cell2
         if( in /= 0 )write(*,*)'+ internal error: corrupt bc ',ib,i,ip,in

         if( it == RWall )then

           Ttmp = Tref
           if( SolveEnthalpy )  Ttmp  = T(Ncel+ib)

           write(IOloc,'(i8,1x,i4,9(1x,1pe10.3))') ip,ir,bnd(ib)%shear, bnd(ib)%yplus,&
                      bnd(ib)%uplus,bnd(ib)%distance,Ttmp-Tref,bnd(ib)%h,bnd(ib)%q

         endif
       end do
     endif
   endif
   !write(*,'(/,'' Face data'')')
   !write(*,'(A,A)')'    Face    ip       in   massflux' 

   !do i=1,Nfac
   !  write(*,'(i8,2(1x,i4),9(1x,1pe10.3))') i, &
   !           Face(i)%cell1,Face(i)%cell2, MassFlux(i)
   !end do

   if( IOprt /= IOdef ) close(IOprt)
   
end subroutine ShowCells
subroutine Update_P_at_boundaries2(ivar,Phi)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(3)   :: GradPhi, Xp, dX, ds 

   real, dimension(3,3) :: A 
   real, dimension(3)   :: RHS_A 

   real, dimension(2,2) :: B
   real, dimension(2)   :: RHS_B

   real :: p1, p2, p3, p4, p5, p6 
    
   integer, dimension(3):: IPIV  ! LAPACK

   logical :: UseX, UseY, UseZ

   if( Debug > -1 )write(*,*)'*** Update_P_at_boundaries2: ',variable(ivar)
    
   do ib=1,Nbnd           
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(*,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ
          
     Xp = Cell(ip)%x 
     GradPhi = 0.0
     A       = 0.0       ! reset A
     B       = 0.0       ! reset B
     RHS_A   = 0.0
     RHS_B   = 0.0
     
     icnt = 0
     UseX = .false.
     UseY = .false.
     UseZ = .false.

     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2

       if( ipn > 0 )then
         ! internal faces only
         icnt = icnt + 1

         if( ipp == ip )then
           dPhi = Phi(ipn) - Phi(ipp)
           dX   = Cell(ipn)%x - Cell(ipp)%x
         else
           dPhi = Phi(ipp) - Phi(ipn)
           dX   = Cell(ipp)%x - Cell(ipn)%x
         endif
      
         A(1,1) = A(1,1) + dX(1)*dX(1)
         A(2,1) = A(2,1) + dX(1)*dX(2)
         A(3,1) = A(3,1) + dX(1)*dX(3)

         A(1,2) = A(2,1)
         A(2,2) = A(2,2) + dX(2)*dX(2)
         A(3,2) = A(3,2) + dX(2)*dX(3)

         A(1,3) = A(3,1)
         A(2,3) = A(3,2) 
         A(3,3) = A(3,3) + dX(3)*dX(3)

         RHS_A(1) = RHS_A(1) + dX(1)*dPhi
         RHS_A(2) = RHS_A(2) + dX(2)*dPhi
         RHS_A(3) = RHS_A(3) + dX(3)*dPhi
         
         !if( ib == 16 )then
         !  write(*,*)'dphi,dx:',dphi,dx
         !  write(*,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
         !  write(*,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
         !  write(*,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
         !endif

       endif 
     end do   

     if( icnt <= 0 )then
       write(*,*)'Foutje!'
     !else
     !  write(*,*)'Aantal relaties:',ip,j,icnt 
     endif
  
     if( UseLapack )then
       
       call SGESV( 3, 1, A, 3, IPIV, RHS_A, 3, INFO )

       if( info /= 0 ) write(*,*)'Lapack (sgesv) info(1):',INFO

       !if( info /= 0 .and. ib ==16)then
       !  write(*,*)'Lapack (sgesv) info(2):',INFO
       !  write(*,*)'ib,ip,ir,it:',ib,ip,ir,it
       !  write(*,*)'Xp:',Xp
       !  write(*,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
       !  write(*,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
       !  write(*,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
       !endif
       
       GradPhi(1) = RHS_A(1)
       GradPhi(2) = RHS_A(2)
       GradPhi(3) = RHS_A(3)
     
     else

       irang = 0
       if( A(1,1) > Small )then  ! altijd >= 0
         UseX = .true.
         irang = irang + 1
       else
         GradPhi(1) = 0.0
       endif    
       if( A(2,2) > Small )then  ! altijd >= 0
         UseY = .true.
         irang = irang + 1
       else
         GradPhi(2) = 0.0
       endif  
       if( A(3,3) > Small )then  ! altijd >= 0
         UseZ = .true.
         irang = irang + 1
       else
         GradPhi(3) = 0.0
       endif  

       if( irang == 1 )then
         !
         ! too simple, do nothing
         !

         GradPhi = 0.0

       else if( irang == 2 )then
         !
         ! only two components can be calculated
         ! assume the third to be zero
         !
         if( UseX .and. UseY )then
           B(1,1) = A(1,1)
           B(1,2) = A(1,2)
           B(2,1) = A(2,1)
           B(2,2) = A(2,2)

           RHS_B(1) = RHS_A(1)
           RHS_B(2) = RHS_A(2)
         else if( UseX .and. UseZ )then
           B(1,1) = A(1,1)
           B(1,2) = A(1,3)
           B(2,1) = A(3,1)
           B(2,2) = A(3,3)

           RHS_B(1) = RHS_A(1)
           RHS_B(2) = RHS_A(3)
         else if( UseY .and. UseZ )then
           B(1,1) = A(2,2)
           B(1,2) = A(2,3)
           B(2,1) = A(3,2)
           B(2,2) = A(3,3)

           RHS_B(1) = RHS_A(2)
           RHS_B(2) = RHS_A(3)
         endif
         !
         ! B inverteren
         !
         if( B(1,1) /= 0.0 )then
           p1 = - B(2,1)/B(1,1)
           B(2,1) = B(2,1) + p1 * B(1,1)
           B(2,2) = B(2,2) + p1 * B(1,2)

           RHS_B(2) = RHS_B(2) + p1 * RHS_B(1)

         else
           write(*,1) 1
         endif
         if( B(2,2) /= 0.0 )then
           p2 = - B(1,2)/B(2,2)
           B(1,1) = B(1,1) + p2 * B(2,1)
           B(1,2) = B(1,2) + p2 * B(2,2)

           RHS_B(1) = RHS_B(1) + p2 * RHS_B(2)

         else
           write(*,1) 2
         endif
         !
         ! normeren
         !
         if( B(1,1) /= 0.0 )then
           p3 = 1.0/B(1,1)
           B(1,1) = p3 * B(1,1)
           B(1,2) = p3 * B(1,2)

           RHS_B(1) = p3 * RHS_B(1)

         else
           write(*,1) 3,B(1,1)
         endif
         if( B(2,2) /= 0.0 )then
           p4 = 1.0/B(2,2)
           B(2,1) = p4 * B(2,1)
           B(2,2) = p4 * B(2,2)

           RHS_B(2) = p4 * RHS_B(2)

         else
           write(*,1) 4,B(2,2)
         endif
         !
         ! de gradient
         !
         if( UseX .and. UseY )then
           GradPhi(1) = RHS_B(1)
           GradPhi(2) = RHS_B(2)
           GradPhi(3) = 0.0
         else if( UseX .and. UseZ )then
           GradPhi(1) = RHS_B(1)
           GradPhi(2) = 0.0
           GradPhi(3) = RHS_B(2)
         else if( UseY .and. UseZ )then
           GradPhi(1) = 0.0
           GradPhi(2) = RHS_B(1)
           GradPhi(3) = RHS_B(2)
         endif

       else if( irang == 3 )then
         !
         ! the case were all three components
         ! can be calculated
         !
         if( A(1,1) /= 0.0 )then
           p1 = - A(2,1)/A(1,1)
           A(2,1)   = A(2,1) + p1 * A(1,1)
           A(2,2)   = A(2,2) + p1 * A(1,2)
           A(2,3)   = A(2,3) + p1 * A(1,3)
           RHS_A(2) = RHS_A(2) + p1 * RHS_A(1)

           p2 = - A(3,1)/A(1,1)
           A(3,1)   = A(3,1) + p2 * A(1,1)
           A(3,2)   = A(3,2) + p2 * A(1,2)
           A(3,3)   = A(3,3) + p2 * A(1,3)
           RHS_A(3) = RHS_A(3) + p2 * RHS_A(1)
         else
           write(*,1) 5
         endif
         if( A(2,2) /= 0.0 )then
           p3 = - A(1,2)/A(2,2)
           A(1,1)   = A(1,1) + p3 * A(2,1)
           A(1,2)   = A(1,2) + p3 * A(2,2)
           A(1,3)   = A(1,3) + p3 * A(2,3)
           RHS_A(1) = RHS_A(1) + p3 * RHS_A(2)

           p4 = - A(3,2)/A(2,2)
           A(3,1)   = A(3,1) + p4 * A(2,1)
           A(3,2)   = A(3,2) + p4 * A(2,2)
           A(3,3)   = A(3,3) + p4 * A(2,3)
           RHS_A(3) = RHS_A(3) + p4 * RHS_A(2)
         else
           write(*,1) 6
         endif
         if( A(3,3) /= 0.0 )then
           p5 = - A(1,3)/A(3,3)
           A(1,1)   = A(1,1) + p5 * A(3,1)
           A(1,2)   = A(1,2) + p5 * A(3,2)
           A(1,3)   = A(1,3) + p5 * A(3,3)
           RHS_A(1) = RHS_A(1) + p5 * RHS_A(3)

           p6 = - A(2,3)/A(3,3)
           A(2,1)   = A(2,1) + p6 * A(3,1)
           A(2,2)   = A(2,2) + p6 * A(3,2)
           A(2,3)   = A(2,3) + p6 * A(3,3)
           RHS_A(2) = RHS_A(2) + p6 * RHS_A(3)
         else
           write(*,1) 7
           !write(*,*)'c : ',ip,Nfaces(ip)
           !write(*,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
           !write(*,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
           !write(*,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
           !write(*,*)'Use:',UseX,UseY,UseZ,irang
         endif
         !
         ! normeren
         !
         if( abs(A(1,1)) > Small )then
           p7 = 1.0/A(1,1)
           A(1,1)   = p7 * A(1,1)
           A(1,2)   = p7 * A(1,2)
           A(1,3)   = p7 * A(1,3)
           RHS_A(1) = p7 * RHS_A(1)
         else
           write(*,1) 8
         endif
         if( abs(A(2,2)) > Small )then
           p8 = 1.0/A(2,2)
           A(2,1)   = p8 * A(2,1)
           A(2,2)   = p8 * A(2,2)
           A(2,3)   = p8 * A(2,3)
           RHS_A(2) = p8 * RHS_A(2)
         else
           write(*,1) 9
         endif
         if( abs(A(3,3)) > Small )then
           p9 = 1.0/A(3,3)
           A(3,1)   = p9 * A(3,1)
           A(3,2)   = p9 * A(3,2)
           A(3,3)   = p9 * A(3,3)
           RHS_A(3) = p9 * RHS_A(3)
         else
           write(*,1) 10
           write(*,*)'c : ',ip,Nfaces(ip)
           write(*,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_a(1)
           write(*,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_a(2)
           write(*,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_a(3)
           write(*,*)'p13:',p1,p2,p3
           write(*,*)'p46:',p4,p5,p6
           write(*,*)'p79:',p7,p8,p9
           write(*,*)'Use:',UseX,UseY,UseZ,irang

           RHS_A(3) = 0.0
         endif
         !
         ! de gradient
         !

         GradPhi(1) = RHS_A(1)
         GradPhi(2) = RHS_A(2)
         GradPhi(3) = RHS_A(3)

       else
         write(*,*)'+ internal error catch 22'
       endif

     2 continue

       !write(*,*)'achter de rug'
       !write(*,*)'B1:',B(1,1),B(1,2),'|',RHS_B(1)
       !write(*,*)'B2:',B(2,1),B(2,2),'|',RHS_B(2)
       !write(*,*)'---'
       !write(*,*)'Cel:',ip,GradPhi

     endif
     
     ds = Face(i)%x - Xp

     Phi(Ncel+ib) = Phi(ip) + dot_product( ds , GradPhi )

   end do

 1 format('+ internal error: Update_P_at_boundaries2 ',i2,e14.4)
   if( Debug > 3 )write(*,*)'=== Update_P_at_boundaries2'

end subroutine Update_P_at_boundaries2
subroutine Update_Scalars_at_boundaries2(ivar,Phi)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(3)   :: GradPhi, Xp, dX, ds 
   real, dimension(3,3) :: A
   real, dimension(2,2) :: B
   real, dimension(3)   :: RHS_A
   real, dimension(2)   :: RHS_B
    
   logical :: UseX, UseY, UseZ

   if( Debug > 3 )write(*,*)'*** Update_Scalars_at_boundaries2: ',variable(ivar)
    
   do ib=1,Nbnd           
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(*,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     !
     ! LET OP: alleen voor k en epsilon ook de wanden
     !         bij T wordt de wandtemp. elders gezet
     !         en mag niet overschreven worden.
     !
     if( it == Rsymp .or. it == Routlet)then
    
       Xp = Cell(ip)%x 
       GradPhi = 0.0
       A     = 0.0           ! reset A
       B     = 0.0           ! reset B
       RHS_A = 0.0
       RHS_B = 0.0

       icnt = 0
       UseX = .false.
       UseY = .false.
       UseZ = .false.

       do j=1,Nfaces(ip)
         k  = CFace(ip,j)
         ipp = Face(k)%cell1
         ipn = Face(k)%cell2

         if( ipn > 0 )then
           ! internal faces only
           icnt = icnt + 1

           if( ipp == ip )then
             dPhi = Phi(ipn) - Phi(ipp)
             dX   = Cell(ipn)%x - Cell(ipp)%x
           else
             dPhi = Phi(ipp) - Phi(ipn)
             dX   = Cell(ipp)%x - Cell(ipn)%x
           endif

           A(1,1) = A(1,1) + dX(1)*dX(1)
           A(2,1) = A(2,1) + dX(1)*dX(2)
           A(3,1) = A(3,1) + dX(1)*dX(3)

           A(1,2) = A(2,1)
           A(2,2) = A(2,2) + dX(2)*dX(2)
           A(3,2) = A(3,2) + dX(2)*dX(3)

           A(1,3) = A(3,1)
           A(2,3) = A(3,2) 
           A(3,3) = A(3,3) + dX(3)*dX(3)

           RHS_A(1) = RHS_A(1) + dX(1)*dPhi
           RHS_A(2) = RHS_A(2) + dX(2)*dPhi
           RHS_A(3) = RHS_A(3) + dX(3)*dPhi

         endif 
       end do   

       if( icnt <= 0 )then
         write(*,*)'Foutje!'
       !else
       !  write(*,*)'Aantal relaties:',ip,j,icnt
       endif


       irang = 0
       if( A(1,1) > Small )then  ! altijd >= 0
         UseX = .true.
         irang = irang + 1
       else
         A(1,1) = 0.0
         GradPhi(1) = 0.0
       endif    
       if( A(2,2) > Small )then  ! altijd >= 0
         UseY = .true.
         irang = irang + 1
       else
         A(2,2) = 0.0
         GradPhi(2) = 0.0
       endif  
       if( A(3,3) > Small )then  ! altijd >= 0
         UseZ = .true.
         irang = irang + 1
       else
         A(3,3) = 0.0
         GradPhi(3) = 0.0
       endif  

       if( irang == 1 )then
         !
         ! too simple, do nothing
         !

         GradPhi = 0.0

       else if( irang == 2 )then
         !
         ! only two components can be calculated
         ! assume the third to be zero
         !
         if( UseX .and. UseY )then
           B(1,1) = A(1,1)
           B(1,2) = A(1,2)
           B(2,1) = A(2,1)
           B(2,2) = A(2,2)

           RHS_B(1) = RHS_A(1)
           RHS_B(2) = RHS_A(2)
         else if( UseX .and. UseZ )then
           B(1,1) = A(1,1)
           B(1,2) = A(1,3)
           B(2,1) = A(3,1)
           B(2,2) = A(3,3)

           RHS_B(1) = RHS_A(1)
           RHS_B(2) = RHS_A(3)
         else if( UseY .and. UseZ )then
           B(1,1) = A(2,2)
           B(1,2) = A(2,3)
           B(2,1) = A(3,2)
           B(2,2) = A(3,3)

           RHS_B(1) = RHS_A(2)
           RHS_B(2) = RHS_A(3)
         endif
         !
         ! B inverteren
         !
         if( B(1,1) /= 0.0 )then
           p1 = - B(2,1)/B(1,1)
           B(2,1) = B(2,1) + p1 * B(1,1)
           B(2,2) = B(2,2) + p1 * B(1,2)

           RHS_B(2) = RHS_B(2) + p1 * RHS_B(1)

         else
           write(*,1) 1
         endif
         if( B(2,2) /= 0.0 )then
           p2 = - B(1,2)/B(2,2)
           B(1,1) = B(1,1) + p2 * B(1,2)
           B(1,2) = B(1,2) + p2 * B(2,2)

           RHS_B(1) = RHS_B(1) + p2 * RHS_B(2)

         else
           write(*,1) 2
         endif
         !
         ! normeren
         !
         if( abs(B(1,1)) > Small )then
           p3 = 1.0/B(1,1)
           B(1,1) = p3 * B(1,1)
           B(1,2) = p3 * B(1,2)

           RHS_B(1) = p3 * RHS_B(1)

         else
           write(*,1) 3
         endif
         if( abs(B(2,2)) > Small )then
           p4 = 1.0/B(2,2)
           B(2,1) = p4 * B(2,1)
           B(2,2) = p4 * B(2,2)

           RHS_B(2) = p4 * RHS_B(2)

         else
           write(*,1) 4
         endif
         !
         ! de gradient
         !
         if( UseX .and. UseY )then
           GradPhi(1) = RHS_B(1)
           GradPhi(2) = RHS_B(2)
           GradPhi(3) = 0.0
         else if( UseX .and. UseZ )then
           GradPhi(1) = RHS_B(1)
           GradPhi(2) = 0.0
           GradPhi(3) = RHS_B(2)
         else if( UseY .and. UseZ )then
           GradPhi(1) = 0.0
           GradPhi(2) = RHS_B(1)
           GradPhi(3) = RHS_B(2)
         endif

       else if( irang == 3 )then
         !
         ! the case were all three components
         ! can be calculated
         !
         if( A(1,1) /= 0.0 )then
           p1 = - A(2,1)/A(1,1)
           A(2,1) = A(2,1) + p1 * A(1,1)
           A(2,2) = A(2,2) + p1 * A(1,2)
           A(2,3) = A(2,3) + p1 * A(1,3)
           RHS_A(2) = RHS_A(2) + p1 * RHS_A(1)

           p2 = - A(3,1)/A(1,1)
           A(3,1) = A(3,1) + p2 * A(1,1)
           A(3,2) = A(3,2) + p2 * A(1,2)
           A(3,3) = A(3,3) + p2 * A(1,3)
           RHS_A(3) = RHS_A(3) + p2 * RHS_A(1)
         else
           write(*,1) 5
         endif
         if( A(2,2) /= 0.0 )then
           p3 = - A(1,2)/A(2,2)
           A(1,1) = A(1,1) + p3 * A(2,1)
           A(1,2) = A(1,2) + p3 * A(2,2)
           A(1,3) = A(1,3) + p3 * A(2,3)
           RHS_A(1) = RHS_A(1) + p3 * RHS_A(2)

           p4 = - A(3,2)/A(2,2)
           A(3,1) = A(3,1) + p4 * A(2,1)
           A(3,2) = A(3,2) + p4 * A(2,2)
           A(3,3) = A(3,3) + p4 * A(2,3)
           RHS_A(3) = RHS_A(3) + p4 * RHS_A(2)
         else
           write(*,1) 6
         endif
         if( A(3,3) /= 0.0 )then
           p5 = - A(1,3)/A(3,3)
           A(1,1) = A(1,1) + p5 * A(3,1)
           A(1,2) = A(1,2) + p5 * A(3,2)
           A(1,3) = A(1,3) + p5 * A(3,3)
           RHS_A(1) = RHS_A(1) + p5 * RHS_A(3)

           p6 = - A(2,3)/A(3,3)
           A(2,1) = A(2,1) + p6 * A(3,1)
           A(2,2) = A(2,2) + p6 * A(3,2)
           A(2,3) = A(2,3) + p6 * A(3,3)
           RHS_A(2) = RHS_A(2) + p6 * RHS_A(3)
         else
           write(*,1) 7
           !write(*,*)'c : ',ip,Nfaces(ip)
           !write(*,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
           !write(*,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
           !write(*,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
           !write(*,*)'Use:',UseX,UseY,UseZ,irang
         endif
         !
         ! normeren
         !
         if( abs(A(1,1)) > Small )then
           p7 = 1.0/A(1,1)
           A(1,1) = p7 * A(1,1)
           A(1,2) = p7 * A(1,2)
           A(1,3) = p7 * A(1,3)
           RHS_A(1) = p7 * RHS_A(1)
         else
           write(*,1) 8
         endif
         if( abs(A(2,2)) > Small )then
           p8 = 1.0/A(2,2)
           A(2,1) = p8 * A(2,1)
           A(2,2) = p8 * A(2,2)
           A(2,3) = p8 * A(2,3)
           RHS_A(2) = p8 * RHS_A(2)
         else
           write(*,1) 9
         endif
         if( abs(A(3,3)) > Small )then
           p9 = 1.0/A(3,3)
           A(3,1) = p9 * A(3,1)
           A(3,2) = p9 * A(3,2)
           A(3,3) = p9 * A(3,3)
           RHS_A(3) = p9 * RHS_A(3)
         else
           write(*,1) 10
           write(*,*)'c : ',ip,Nfaces(ip)
           write(*,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
           write(*,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
           write(*,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
           write(*,*)'Use:',UseX,UseY,UseZ,irang

           rhs_A(3) = 0.0
         endif
         !
         ! de gradient
         !

         GradPhi(1) = RHS_A(1)
         GradPhi(2) = RHS_A(2)
         GradPhi(3) = RHS_A(3)

       else
         write(*,*)'+ internal error catch 22'
       endif

       !write(*,*)'achter de rug'
       !write(*,*)'B1:',B(1,1),B(1,2),'|',RHS_B(1)
       !write(*,*)'B2:',B(2,1),B(2,2),'|',RHS_B(2)
       !write(*,*)'---'
       !write(*,*)'Cel:',ip,GradPhi

       ds = Face(i)%x - Xp

       Phi(Ncel+ib) = Phi(ip) + dot_product( ds , GradPhi )
     
     endif
   end do

 1 format('+ internal error: Update_Scalars_at_boundaries2 ',i2,e14.4)
   if( debug > 3 )write(*,*)'=== Update_Scalars_at_boundaries2'

end subroutine Update_Scalars_at_boundaries2
subroutine GradientPhi(ivar,Phi,dPhidX)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd), intent(IN)    :: Phi
   real, dimension(Ncel+Nbnd,3), intent(OUT) :: dPhidX   

   integer, intent(IN) :: ivar

!   if( GradAlg == GradLS .and. ivar == VarT  )then
   if( GradAlg == GradLS )then
     call GradientPhiLeastSquares(ivar,Phi,dPhidX)
!   else if( GradAlg == GradGauss .and. ivar == VarT )then
   else if( GradAlg == GradGauss )then
     call GradientPhiGauss(ivar,Phi,dPhidX)
     !call GradientPhiLeastSquaresN(ivar,Phi,dPhidX)
     !call GradientPhiNodes(ivar,Phi,dPhidX)
     !call GradientPhiLeastSquaresFaces(ivar,Phi,dPhidX)
   else
     call GradientPhiLeastSquares(ivar,Phi,dPhidX)
    !write(*,*)'+ internal error invalid gradient algorithm'
   endif

end subroutine GradientPhi
subroutine GradientPhiLeastSquares(ivar,Phi,dPhidX)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX   

   real, dimension(3)   :: GradPhi, Xf, Xac, Xp, dX, ds 
   real, dimension(3,3) :: A
   real, dimension(3)   :: RHS_A
    
   integer, dimension(3):: IPIV    ! LAPACK

   if( Debug > 3 ) write(*,*)'*** GradientPhiLeastSquares ',Variable(ivar)
   write(IOdbg,*)'*** GradientPhiLeastSquares ',Variable(ivar)

   dPhidX = 0.0
      
   do ip=1,Ncel         
          
     Xp = Cell(ip)%x 
     A     = 0.0
     RHS_A = 0.0     
     
     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2

       if( ipn > 0 )then
         !
         ! original method Sum from P to N
         !
         if( ipp == ip )then
           dPhi = Phi(ipn) - Phi(ipp)
           dX   = Cell(ipn)%x - Cell(ipp)%x
         else
           dPhi = Phi(ipp) - Phi(ipn)
           dX   = Cell(ipp)%x - Cell(ipn)%x 
         endif

         !if( ip == 8)then
         !  write(*,*) 'ls :',ip,ipp,':',dphi
         !  write(*,*) 'ls :',dX         
         !endif
         !
         ! alternative use face: P to Xf
         ! in case of an orthogonal mesh nothing changes
         ! however
         ! THE disadvantage is the unknown correcion
         ! since dPhidX is unknown
         !
         !facn = Face(k)%lambda
         !facp = 1.0 - facn
         !Xac  = Cell(ipn)%x * facn + Cell(ipp)%x * facp
         !Xf   = Face(k)%x 
         
         !Phiac   = Phi(ipn) * facn + Phi(ipp) * facp 
         !GradPhi = dPhidX(ipn,:) * facn + dPhidX(ipp,:) * facp
         !delta   = dot_product( GradPhi , Xf - Xac )  ! correction

         !PhiFace = Phiac !+ delta

         !if( ipp == ip )then
         !  dPhi = PhiFace - Phi(ipp)
         !  dX   = Xf - Cell(ipp)%x
         !else
         !  dPhi = Phi(ipp) - PhiFace
         !  dX   = Cell(ipp)%x - Xf
         !endif
         !

         A(1,1) = A(1,1) + dX(1)*dX(1)    
         A(2,1) = A(2,1) + dX(1)*dX(2)
         A(3,1) = A(3,1) + dX(1)*dX(3)

         A(1,2) = A(2,1)
         A(2,2) = A(2,2) + dX(2)*dX(2)
         A(3,2) = A(3,2) + dX(2)*dX(3)

         A(1,3) = A(3,1)
         A(2,3) = A(3,2) 
         A(3,3) = A(3,3) + dX(3)*dX(3)

         RHS_A(1) = RHS_A(1) + dX(1)*dPhi
         RHS_A(2) = RHS_A(2) + dX(2)*dPhi
         RHS_A(3) = RHS_A(3) + dX(3)*dPhi
       
       else
         ! boundary face 

         ib = Face(k)%bnd
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         
         dPhi = Phi(Ncel+ib) - Phi(ipp)
         dX   = Face(k)%x - Cell(ipp)%x
   	 
         !if( ip == 8)then
         !  write(*,*) 'lsb:',ip,ipp,':',dphi
         !  write(*,*) 'lsb:',dX         
         !endif

         A(1,1) = A(1,1) + dX(1)*dX(1)       
         A(2,1) = A(2,1) + dX(1)*dX(2)
         A(3,1) = A(3,1) + dX(1)*dX(3)

         A(1,2) = A(2,1)
         A(2,2) = A(2,2) + dX(2)*dX(2)
         A(3,2) = A(3,2) + dX(2)*dX(3)

         A(1,3) = A(3,1)
         A(2,3) = A(3,2) 
         A(3,3) = A(3,3) + dX(3)*dX(3)

         RHS_A(1) = RHS_A(1) + dX(1)*dPhi
         RHS_A(2) = RHS_A(2) + dX(2)*dPhi
         RHS_A(3) = RHS_A(3) + dX(3)*dPhi
         
       endif 
     end do   
  
     if( UseLapack )then

       call SGESV( 3, 1, A, 3, IPIV, RHS_A, 3, INFO )

       if( info /= 0 ) write(*,*)'Lapack (sgesv) info(1):',INFO
     
       GradPhi(1) = RHS_A(1)
       GradPhi(2) = RHS_A(2)
       GradPhi(3) = RHS_A(3)

       dPhidX(ip,:) = GradPhi 
     
     else     

       call A33xB3(A,RHS_A)

       GradPhi(1) = RHS_A(1)
       GradPhi(2) = RHS_A(2)
       GradPhi(3) = RHS_A(3)

       dPhidX(ip,:) = GradPhi 

     endif

   end do

   if( allocated(DXdebug) .and. allocated(T) ) DXdebug =T(1:Ncel)
   if( allocated(DXgrad) )                     DXgrad = dPhidX

   if( Debug > 2 )then
     write(*,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,1)),' <  dPhi/dX  < ',maxval(dPhidX(1:Ncel,1))
     write(*,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,2)),' <  dPhi/dY  < ',maxval(dPhidX(1:Ncel,2))
     write(*,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,3)),' <  dPhi/dZ  < ',maxval(dPhidX(1:Ncel,3))
   endif
   if( Debug > 3 ) write(*,*)'=== GradientPhiLeastSquares'

end subroutine GradientPhiLeastSquares
subroutine GradientPhiLeastSquaresN(ivar,Phi,dPhidX)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX   

   real, dimension(3)   :: GradPhi, Xf, Xac, Xp, dX, ds 
   real, dimension(3,3) :: A
   real, dimension(3)   :: RHS_A
    
   integer, dimension(3):: IPIV    ! LAPACK

   if( Debug > 1 ) write(*,*)'*** GradientPhiLeastSquaresN ',Variable(ivar)
   
   dPhidXo = 0.0
   !
   ! iterative calculation of gradients
   !
   do igrad=1,Ngradient
     write(*,*)'loop ',igrad
     dPhidX  = 0.0
     do ip=1,Ncel         

       Xp = Cell(ip)%x 
       A  = 0.0
       RHS_A = 0.0     

       do j=1,Nfaces(ip)
         k  = CFace(ip,j)
         ipp = Face(k)%cell1
         ipn = Face(k)%cell2

         if( ipn > 0 )then
           
           facn = Face(k)%lambda
           facp = 1.0 - facn
           Xac  = Cell(ipn)%x * facn + Cell(ipp)%x * facp
           Xf   = Face(k)%x 
         
           Phiac   = Phi(ipn) * facn + Phi(ipp) * facp 
           GradPhi = dPhidXo(ipn,:) * facn + dPhidXo(ipp,:) * facp
           delta   = dot_product( GradPhi , Xf - Xac )  ! correction
           PhiFace = Phiac + delta

           !PhiFace = 0.5*( Phi(ipp) + Phi(ipn) + &
           !dot_product( dPhidXo(ipp,:) , Xf - Cell(ipp)%x ) + &
           !dot_product( dPhidXo(ipn,:) , Xf - Cell(ipn)%x ) )

           if( ipp == ip )then
             dPhi = PhiFace - Phi(ipp)
             dX   = Xf - Cell(ipp)%x
           else
             dPhi = Phi(ipp) - PhiFace
             dX   = Cell(ipp)%x - Xf
           endif

           A(1,1) = A(1,1) + dX(1)*dX(1)    
           A(2,1) = A(2,1) + dX(1)*dX(2)
           A(3,1) = A(3,1) + dX(1)*dX(3)

           A(1,2) = A(2,1)
           A(2,2) = A(2,2) + dX(2)*dX(2)
           A(3,2) = A(3,2) + dX(2)*dX(3)

           A(1,3) = A(3,1)
           A(2,3) = A(3,2) 
           A(3,3) = A(3,3) + dX(3)*dX(3)

           RHS_A(1) = RHS_A(1) + dX(1)*dPhi
           RHS_A(2) = RHS_A(2) + dX(2)*dPhi
           RHS_A(3) = RHS_A(3) + dX(3)*dPhi

         else
           ! boundary face 

           ib = Face(k)%bnd
           ir = Bnd(ib)%rid
           it = Reg(ir)%typ

           dPhi = Phi(Ncel+ib) - Phi(ipp)
           dX   = Face(k)%x - Cell(ipp)%x

           A(1,1) = A(1,1) + dX(1)*dX(1)       
           A(2,1) = A(2,1) + dX(1)*dX(2)
           A(3,1) = A(3,1) + dX(1)*dX(3)

           A(1,2) = A(2,1)
           A(2,2) = A(2,2) + dX(2)*dX(2)
           A(3,2) = A(3,2) + dX(2)*dX(3)

           A(1,3) = A(3,1)
           A(2,3) = A(3,2) 
           A(3,3) = A(3,3) + dX(3)*dX(3)

           RHS_A(1) = RHS_A(1) + dX(1)*dPhi
           RHS_A(2) = RHS_A(2) + dX(2)*dPhi
           RHS_A(3) = RHS_A(3) + dX(3)*dPhi

         endif 
       end do   

       if( UseLapack )then
         call SGESV( 3, 1, A, 3, IPIV, RHS_A, 3, INFO )

         if( info /= 0 ) write(*,*)'Lapack (sgesv) info(1):',INFO

         GradPhi(1) = RHS_A(1)
         GradPhi(2) = RHS_A(2)
         GradPhi(3) = RHS_A(3)

         dPhidX(ip,:) = GradPhi 

       else     
         call A33xB3(A,RHS_A)

         GradPhi(1) = RHS_A(1)
         GradPhi(2) = RHS_A(2)
         GradPhi(3) = RHS_A(3)

         dPhidX(ip,:) = GradPhi 

       endif

     end do ! cell loop

     dPhidXo = dPhidX 
     !dPhidXo = dPhidXo + 0.75*( dPhidX - dPhidXo) ! under relaxation

     if( Debug > 2 )then
       write(*,'(1x,1pe10.3,A,e10.3)') &
         minval(dPhidX(1:Ncel,1)),' <  dPhi/dX  < ',maxval(dPhidX(1:Ncel,1))
       write(*,'(1x,1pe10.3,A,e10.3)') &
         minval(dPhidX(1:Ncel,2)),' <  dPhi/dY  < ',maxval(dPhidX(1:Ncel,2))
       write(*,'(1x,1pe10.3,A,e10.3)') &
         minval(dPhidX(1:Ncel,3)),' <  dPhi/dZ  < ',maxval(dPhidX(1:Ncel,3))
     endif

   end do   ! gradient loop

   if( Debug > 3 ) write(*,*)'=== GradientPhiLeastSquaresN'

end subroutine GradientPhiLeastSquaresN
subroutine Set_Normalisation_Factors
!========================================================================

   use constants
   use geometry
   use variables
   use scalars

   if( Debug > 3 ) write(*,*)'*** Set_Normalisation_Factors'

   areain  = 0.0
   fluxin  = 0.0

   sum0 = 0.0
   sum1 = 0.0    
   sum2 = 0.0
   sum3 = 0.0
   
   areaout = 0.0

   icntin  = 0
   icntout = 0

   do ib=1,Nbnd
     
     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == RInlet )then
       icntin = icntin + 1
       areain = areain + Face(i)%area
       
       if( Reg(ir)%user )then
       
         ip = Face(i)%cell1
         
         Utmp  = Reg(ir)%uvw(1)
         Vtmp  = Reg(ir)%uvw(2)
         Wtmp  = Reg(ir)%uvw(3)
         TEtmp = Reg(ir)%k
         EDtmp = Reg(ir)%e
         Ttmp  = Reg(ir)%T
         Dtmp  = Reg(ir)%den
         
         Xtmp  = Face(i)%x(1)
         Ytmp  = Face(i)%x(2)
         Ztmp  = Face(i)%x(3)
         Atmp  = Face(i)%area
         
         call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                        Pref,Tref,Iter,Time,                  &
                        Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                        Reg(ir)%name1,Reg(ir)%name2)

         fluxin =  Dtmp * ( Utmp*Face(i)%n(1)   &
                          + Vtmp*Face(i)%n(2)   &
                          + Wtmp*Face(i)%n(3) )
         sum0   = sum0 + fluxin

         sum1   = sum1 + &
                  fluxin *sqrt( Utmp*Utmp + Vtmp*Vtmp + Wtmp*Wtmp )

         sum2   = sum2 + &
                  fluxin * ( Utmp*Utmp + Vtmp*Vtmp + Wtmp*Wtmp )

         if( SolveEnthalpy )then

           sum3 = sum3 + &
                     Cp(Ncel+ib) * Ttmp * Dtmp * &
                       ( Utmp*Face(i)%n(1)       &
                       + Vtmp*Face(i)%n(2)       &
                       + Wtmp*Face(i)%n(3) )
         
         endif
       else
       
         fluxin = Reg(ir)%den * dot_product( Reg(ir)%uvw , Face(i)%n )
       
         sum0   = sum0 + fluxin

         sum1   = sum1 + &
                  fluxin *sqrt( dot_product( Reg(ir)%uvw , Reg(ir)%uvw ) )

         sum2   = sum2 + &
                  fluxin * dot_product( Reg(ir)%uvw , Reg(ir)%uvw ) 

         if( SolveEnthalpy )then

           sum3 = sum3 + &
                     Cp(Ncel+ib) * Reg(ir)%T * Reg(ir)%den * &
                     dot_product( Reg(ir)%uvw , Face(i)%n )

         endif
       endif

     else if( it == ROutlet )then
       icntout = icntout + 1
       areaout = areaout + Face(i)%area
     endif
   end do

   write(IOdbg,*)'Residuals with:'
   write(IOdbg,*)'   fluxin: ',fluxin
   write(IOdbg,*)'   areain: ',areain
   write(IOdbg,*)'  areaout: ',areaout
   write(IOdbg,*)'     sum0: ',sum0
   write(IOdbg,*)'     sum1: ',sum1
   write(IOdbg,*)'     sum2: ',sum2
   if( SolveEnthalpy ) write(IOdbg,*)'  sum3: ',sum3   

   ResiNorm(:) = 0.0
   
   if( icntin /= 0 )then
     
     if( SolveU ) ResiNorm(VarU) = 1.0 / (sum1 + Small)
     if( SolveV ) ResiNorm(VarV) = 1.0 / (sum1 + Small)
     if( SolveW ) ResiNorm(VarW) = 1.0 / (sum1 + Small)
     
     if( SolveP ) ResiNorm(VarP) = 1.0 / (sum0 + Small)
  
     if( SolveTurbEnergy ) ResiNorm(VarTE) = 1.0 / (sum2 + Small)
     
     if( SolveTurbDiss )then
        Uin = fluxin / areain              ! mean velocity
        dL  = 1.0                          ! some length scale
        ResiNorm(VarED) = 1.0 / (sum2 + Small) * dL / Uin
     endif
     
     if( SolveEnthalpy ) ResiNorm(VarT) = 1.0 / (sum3 + Small)
     
     if( SolveScalars ) ResiNorm(VarSC) = 1.0   ! not used!              !<= MOET BETER!
     
     ResiNorm(:) = abs( ResiNorm )
      
   else
     
     ResiNorm(:) = 1.0

     if( SolveEnthalpy .and. sum3 > Small )then
       ResiNorm(VarT) = 1.0/sum3
     endif

   endif

   write(IOdbg,*)'Using ResiNorm 1-3:',ResiNorm(1:3)
   write(IOdbg,*)'      ResiNorm 4-6:',ResiNorm(4:6)
   write(IOdbg,*)'      ResiNorm 7-8:',ResiNorm(7:8)

   if( SolveScalars ) call Set_Scalar_Norm()


   if( Debug > 3 ) write(*,*)'=== Set_Normalisation_Factors'

end subroutine Set_Normalisation_Factors
subroutine CalculateViscosity
!========================================================================
!                                     2
!  Eddy viscosity: VisTurb = rho Cmu k / epsilon (eq. 9.43)
!
!  Effective viscosity: VisEff = VisLam + VisTurb
!
!========================================================================

   use constants
   use geometry
   use variables

   logical warning

   if( Debug > 3 ) write(*,*)'*** CalculateViscosity'

   warning = .false.
   Cmu    = TMCmu
   Cmu25  = Cmu**(0.25)
   VisURF = 1.0   ! 0.99

   do ip=1,Ncel
     VisOld = VisEff(ip)
     if( ED(ip) > 1.e-8 )then
       VisNew = VisLam + Cmu*Den(ip)*TE(ip)**2/(ED(ip)+Small)
     else
       VisNew = VisLam 
     endif
     VisEff(ip) = VisOld + VisURF*( VisNew - VisOld )
   end do

   !
   ! that was the easy part... now modify the viscosities
   ! at the boundaries
   !
   
   do ir=1,Nreg
   
     if( Reg(ir)%typ == RWall .and. Reg(ir)%std )then
       !
       ! standard smooth wall: u+ = 1/k ln(Ey+)
       !
       yplus = 11.
       elog  = Reg(ir)%elog
       
       i = 0
       do while( abs( yplus - yplustmp ) > 0.001 .and. i < 50 ) 
         yplustmp = yplus
         yplus    = log( elog * yplustmp )/kappa
         i = i + 1
       end do
       
       Reg(ir)%ylog = yplus
       
       write(IOdbg,'(1x,''Region'',i2,'' y+m ='',f5.1,i3)') ir,yplus,i

     else if( Reg(ir)%typ == RWall .and. .not. Reg(ir)%std )then
       !
       ! rough wall:  u+ = B_rough + 1/k ln( (y-do)/z0 )
       !
       
       Reg(ir)%ylog = 0.0              ! y+m = 0.0 by definition

     endif
   end do
   
   YplusMin = Large
   YplusMax = Small
   
   do ib=1,Nbnd
   
     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in /= 0 )write(*,*)'+ internal error: corrupt bc ',ib,i,ip,in
          
     if( it == RInlet )then
       
       if( Reg(ir)%user )then
       
         ip = Face(i)%cell1
         
         Utmp  = Reg(ir)%uvw(1)
         Vtmp  = Reg(ir)%uvw(2)
         Wtmp  = Reg(ir)%uvw(3)
         TEtmp = Reg(ir)%k
         EDtmp = Reg(ir)%e
         Ttmp  = Reg(ir)%T
         Dtmp  = Reg(ir)%den
         
         Xtmp  = Face(i)%x(1)
         Ytmp  = Face(i)%x(2)
         Ztmp  = Face(i)%x(3)
         Atmp  = Face(i)%area
         
         call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                        Pref,Tref,Iter,Time,                  &
                        Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                        Reg(ir)%name1,Reg(ir)%name2)
       
         VisEff(Ncel+ib)  = VisLam + &
           Dtmp * Cmu * TEtmp**2 / (EDtmp + Small)

       else
       
         VisEff(Ncel+ib)  = VisLam + &
           Reg(ir)%den * Cmu * Reg(ir)%k**2 / (Reg(ir)%e + Small)
       
       endif
     else if( it == ROutlet )then
       !
       ! simpel
       !
       VisEff(Ncel+ib)  = VisEff(ip)

     else if( it == RSymp )then

       VisEff(Ncel+ib)  = VisEff(ip)                ! uitvoeriger zoals bij Rwall?

     else if( it == RWall )then

       dist = Bnd(ib)%distance                      ! wall distance
       turb = TE(ip)                                ! k of wall cel
       
       if( turb < Zero )then
         turb = Zero                                ! don't let k become negative
         write(*,*)'+ CalculateViscosity: Found negative k, set to zero.'
       endif
       
       Utau  = Cmu25*sqrt( turb )                   ! schuifspanningssnelheid (9.48)
                
       Yplus = dist * Utau * Den(ip) /        &     ! y+ (9.47)
                               (VisLam + Small)     
       Bnd(ib)%yplus = Yplus

       if( Reg(ir)%std  )then      !<<===== LET OPPPPP!!!!
         !
         ! standard wall function based on:
         !
         !   u+ = 1/kappa ln( E y+ )
         !
         ! test if y+ < Y+m; then in viscous layer => u+ = y+
         !
         ! the product of E y+ may not become less or equal than 1
         ! and the natural logarithm becomes negative
         ! (which means that the central node submerges into the
         !  sandrougness). Use a slightly higher value of say 1.1 
         !
         if( Yplus < Reg(ir)%ylog  )then
           
           Uplus = Yplus
           
           if( Yplus < Small )then
             Yplus = 0.0
             Uplus = Small
           endif
         else
           tmp = Reg(ir)%elog * Bnd(ib)%yplus 
           if( tmp < 1.1 )then
             write(*,*)'+ CalculateViscosity: Warning: E y+ < 1.1 ',tmp
             tmp = 1.1
           endif
           !
           !  u+ = 1/kappa ln( E y+ )
           !              
           
           Uplus = log( tmp ) / kappa
         
         endif

         Bnd(ib)%uplus = Uplus
         
         YplusMin = min( YplusMin, Yplus )
         YplusMax = max( YplusMax, Yplus )
         !
         ! adhere to the workflow of sub. FluxUVW in which
         !
         !       Tau_(wall) = VisEff * U_(tan) / n
         !
         ! now Tau_(wall) is known and therefore:
         !
         !     VisEff = Tau_(wall) / ( U_(tan) / n )
         !
         !                        2
         !            = rho * Utau / ( U_(tan) / n )
         !
         !            = ( rho * Utau * n / lamvisc )*
         !                               ( lamvisc / ( U_(tan) / Utau )
         !
         !            = Y+ * lamvisc / U+
         !
         !            = Y+/U+ * lamvisc     (ie a factor for the viscosity) 
         !
         VisEff(Ncel+ib) = &                       ! a viscosity lower than
               max( 1.0 , Yplus/Uplus ) * VisLam   ! VisLam is quite unlikely
                  
       else
         !                           y - d0
         ! rough wall:  u+ = 1/k ln( ------ )
         !                             z0
         if( .not. warning )then
           write(*,*)'>> RUWE WAND VOOR ',ir
           warning = .true.
         endif
         
         d0 = 0.0                                  ! displacement thickness
         z0 = 0.1                                  ! "aerodynamic" roughness
         
         Uplus = log( (dist-d0)/z0 ) / kappa

         Bnd(ib)%uplus = Uplus
         
         VisEff(Ncel+ib) = &                       
               max( 1.0 , Yplus/Uplus ) * VisLam    
         
         YplusMin = min( YplusMin, Yplus )
         YplusMax = max( YplusMax, Yplus )
       endif
     else
       write(*,*)'+ CalculateViscosity: Unknown boundary condition'
       write(*,*)'i,ir,it: ',i,ir,it
     endif

   end do       
      
   write(IOdbg,*) YplusMin,'< y+ <',YplusMax

   !
   ! at the start of a turb. simulation k and epsilon may produce
   ! curious results and as a result viseff may increase to unrealistically
   ! high values. this is a kind of limiter (6 orders is sufficient?)
   !
   VisEffMax = maxval( VisEff )
!   if( VisEffMax > 20000000.* VisLam ) Flags(IFlagVis) = 'v'
   if( VisEffMax > 2000000.* VisLam ) Flags(IFlagVis) = 'v'
   
   do ip=1,Ncel
     VisEff(ip) = min( VisEff(ip) , 2000000.* VisLam )         
   end do


   if( Debug > 3 ) write(*,*)'=== CalculateViscosity'

end subroutine CalculateViscosity
subroutine TurbulenceModels(iVar,Phi,dPhidX)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX   
   real, dimension(3)           :: dUdXp, dVdXp, dWdXp
   real, dimension(3)           :: Uw, Up, Xpn, Xn, Un, Ut

   if( Debug > 3 ) write(*,*)'*** TurbulenceModels',TurbModel,Variable(iVar)


   if( TurbModel == TMkeps .or. TurbModel == TMrng )then
     !
     ! standard k-epsilon model
     !
     if( Debug > 3 ) write(*,*)'k-epsilon model'
     if( iVar == VarTE )then
       if( Debug > 3 ) write(*,*)'Turb. kinetic energy'
       
       !
       ! LHS = VisT * Production - rho * epsilon 
       !     - compr. effects + non lin. term
       !
       ! the diffusion term is absorbed in the transp. eq.
       !
       ! source term
       !
       do ip=1,Ncel
         
         dUdXp = dUdX(ip,:)
         dVdXp = dVdX(ip,:)
         dWdXp = dWdX(ip,:)
       
         !
         ! rate of production of turbulent energy (eq. 9.40)
         !
       
         s1   = (dUdXp(1)+dUdXp(1))*dUdXp(1) + (dUdXp(2)+dVdXp(1))*dUdXp(2) + (dUdXp(3)+dWdXp(1))*dUdXp(3)
         s2   = (dVdXp(1)+dUdXp(2))*dVdXp(1) + (dVdXp(2)+dVdXp(2))*dVdXp(2) + (dVdXp(3)+dWdXp(2))*dVdXp(3)
         s3   = (dWdXp(1)+dUdXp(3))*dWdXp(1) + (dWdXp(2)+dVdXp(3))*dWdXp(2) + (dWdXp(3)+dWdXp(3))*dWdXp(3)
       
         VisT = VisEff(ip) - VisLam
         
         Pk   = VisT * ( s1 + s2 + s3 )  
         
         !
         ! bouyancy production term: - Gi/(sigma_h,t rho) drho/dx
         !
         !Pbouy =
         
         TurbP(ip) = Pk
         
         !
         ! compres. amplification term
         !
         !Comp =

         !
         ! dissipation
         !
         
         Dis = Den(ip) * ED(ip)
         
         !
         ! finally the coefficients
         !
         ! pull -rho.eps on the RHS to the LHS by
         ! multiplying with k/k (=1.0) and using
         ! rho.eps/k as the coefficient.
         !
         
         Su(ip) = Su(ip) + TurbP(ip) * Cell(ip)%vol
         Au(ip) = Au(ip) + Dis /(TE(ip) + Small) * Cell(ip)%vol
          
       end do
       !
       ! solid walls
       !
       Cmu   = TMCmu
       Cmu25 = Cmu**0.25

       do ib=1,Nbnd
         i  = Bnd(ib)%face
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         ip = Face(i)%cell1
         in = Face(i)%cell2
         if( in /= 0 )write(*,*)'*** Error: Bnd array corrupt'

         if( it == RWall .and. .not. Initialisation )then
         
           !
           ! in the preceding loop the source term for ip has been
           ! calculated; it needs to be replaced next to a solid
           ! wall; thus first subtract the 'standard' production term
           !
           Su(ip) = Su(ip) - TurbP(ip) * Cell(ip)%vol
           
           !pk_oud = TurbP(ip)
           
           Uw(1)    = Reg(ir)%uvw(1)            ! 
           Uw(2)    = Reg(ir)%uvw(2)            ! wall bc's
           Uw(3)    = Reg(ir)%uvw(3)            !

           Up(1) = U(ip)
           Up(2) = V(ip)
           Up(3) = W(ip)

           Visac    =  VisEff(Ncel+ib)
           
           Xpn      = Face(i)%x - Cell(ip)%x    ! dist. wall center to P

           !coef = Visac * Face(i)%area / vector_length(Xpn)

           Xn = Face(i)%n                       ! de surface vector
           call normalise(Xn)                   ! normaal vector met lengte 1
        
           Up = Up - Uw                         ! snelh. verschil

           dp = dot_product( Up , Xn )
           Un = dp * Xn                         ! normale snelh. component
           Ut = Up - Un                         ! tang. snelh. component

           Uvel  = vector_length( Ut )          !
           dn    = Bnd(ib)%distance             ! alleen de grootte nodig
           Tau_w = Visac * Uvel / dn            ! niet de richting

           if(  dn   <= Zero ) write(*,*)'Foute boel 1!'
           if( Tau_w <= Zero ) write(*,*)'Foute boel 2!'
           
           !
           ! the production in the wall region is (eq. 9.50)
           !
           !    Pk = Tau_w d( Ut )/dn
           !
           ! for a node in the law of the wall region is (eq. 9.51)
           !
           !   d( Ut )/dn = u_tau/( kappa dn )
           ! 

           if( Bnd(ib)%yplus > Reg(ir)%ylog  )then
           
             turb = TE(ip)
             if( turb < Zero )then
                turb = 0.001                          ! don't let k become negative
                write(*,*)'+ TurbulenceModels: Found negative k, set to zero.'
             endif
         
             Utau  = Cmu25*sqrt(turb)                     ! (eq. 9.48)

             !Tau_w = Den(ip)*Utau * kappa * Uvel / &      ! (eq. 9.49) BUG!!!
             !        log( Bnd(ib)%yplus * Reg(ir)%elog )  ! stond ylog ipv elog
             !TurbP(iP) = Tau_w * Utau /( kappa * dn )     ! (eq. 9.50/9.51)

             Tau_w     = Visac * Uvel / dn            
             TurbP(iP) = Tau_w * Utau /( kappa * dn )     ! (eq. 9.50/9.51)


           else
           
             Tau_w     = Visac * Uvel / dn            
             Utau      = sqrt( Tau_w / Den(ip) )
             
             TurbP(iP) = Tau_w * Utau /( kappa * dn )  
           
           endif       

           TE(Ncel+ib) = TE(ip)
         
           !if( ip < 5 )then
           !  write(*,*)'Pk:',ip,Te(ip),TurbP(iP),pk_oud,Tau_w
           !  write(*,*)'  :',ip,Uvel,Bnd(ib)%yplus,Reg(ir)%elog,kappa
           !  write(*,*)'  :',ip,Utau,Tau_w*face(i)%area
           !  write(*,*)' c:',U(ip),V(ip),TE(ip),ED(ip)
           !  write(*,*)' y:',bnd(ib)%yplus,bnd(ib)%uplus,Reg(ir)%ylog
           !endif
           !
           ! finally set the source term for the wall cel
           !
           
           Su(ip) = Su(ip) + TurbP(ip) * Cell(ip)%vol
           
         else if( it == RWall .and. Initialisation )then

           TE(Ncel+ib) = TE(ip)
         
         endif 
       end do
        
     elseif( iVar == VarED )then
       if( Debug > 3 ) write(*,*)'Turb. dissipation'
    
       !
       ! see eq. 9.42 for:
       !
       ! RHS = Ce1 Prod eps/k - rho Ce2 eps^2/k
       !                        
       !     = Ce1 Prod eps/k - rho Ce2 eps/k * eps
       !                        -------------
       !                              v
       !                   in the diagonal A-term * eps
       ! source term
       !
       do ip=1,Ncel

         fact = ED(ip)/(TE(ip)+Small) * Cell(ip)%vol 
       
         Su(ip) = Su(ip) + TMCeps1 * fact * TurbP(ip)
         Au(ip) = Au(ip) + TMCeps2 * fact * Den(ip)
       
       end do
    
       !
       ! solid walls
       !
       Cmu   = TMCmu
       Cmu75 = Cmu**0.75

       do ib=1,Nbnd
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         if( it == RWall .and. .not. Initialisation )then
           i  = Bnd(ib)%face
           ip = Face(i)%cell1
           in = Face(i)%cell2
           if( in > 0 )write(*,*)'*** Error: Bnd array corrupt'
         
           !
           ! due to the requirement Prod = Diss just
           ! force/fix the value in the wall cell (eq. 9.52)
           !            
           
           turb = TE(ip)
           if( turb < Zero )then
             turb = 0.001                       ! don't let k become negative
             write(*,*)'+ TurbulenceModels: Found negative k, set to zero (2).'
           endif
           
           dn  = Bnd(ib)%distance
         
           Dis = Cmu75 * turb**1.5 /( dn * kappa ) 
           
           !
           ! to brutally fix a value enforce Ap * phi = Su
           ! (more efficient than bumping Su and Au?)
           !
           do j=1,NFaces(ip)
             k   = CFace(ip,j)
             ipf = Face(k)%cell1
             inf = Face(k)%cell2
             if( inf > 0 )then
               ! internal
               if( ipf == ip )then
                 RFace(k,2) = 0.0
               elseif( inf == ip )then
                 RFace(k,1) = 0.0
               else
                 write(*,*)'+ internal error: No fix'
               endif
             endif
           end do

           if( .not. Initialisation )then               
             ED(ip) = Dis
             Su(ip) = Dis
             Au(ip) = 1.0 
           endif

           ED(Ncel+ib) = ED(ip)
           
         else if( it == RWall .and. Initialisation )then

           ED(Ncel+ib) = ED(ip)

         endif
       end do
     else
       write(*,*)'+ internal error: invalid turb. model scalar:',Variable(iVar)
     endif
   !elseif( TurbModel == TMrng )then
   !  !
   !  ! RNG k-epsilon model
   !  !
   !  if( debug > 0 ) write(*,*)'RNG k-epsilon model'
     
   else
     !
     ! unknown model
     !
     write(*,*) 'Unknown turb. model / not implemented'
   endif

   if( Debug > 3 ) write(*,*)'=== TurbulenceModels ',Variable(iVar)

end subroutine TurbulenceModels

