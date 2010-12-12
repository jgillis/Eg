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
subroutine WriteRestartField
!
!  simple wrapper subroutine
!

   call WriteRestartField_400 

end subroutine WriteRestartField
subroutine ReadRestartField(Iteration)
!
!  simple wrapper subroutine
!
   use constants

   if( Debug > 0 ) write(*,*) 'Opening restart file'
   call openfile(IOrst,Casename,'.rst','UNFORMATTED', &
                                       'SEQUENTIAL','UNKNOWN',Debug)

   read(IOrst) IRVersion
   if( IRVersion <= 300 )then
     write(*,*)'+ Restart file version < 3.00'
     call ReadRestartField_300(Iteration)
   else
     call ReadRestartField_400(Iteration)
   endif
   
end subroutine ReadRestartField
integer function IntStore(i,dummy)

   use constants
   use variables

   logical :: dummy

!write(*,*) 'test'
   if( dummy )then
!     if( Debug > 0 ) write(*,*) 'Storing ',Variable(i)
     IntStore = 1
   else
     IntStore = 0
   endif
!write(*,*) 'klaar'

end function IntStore
logical function IsStored(i,idummy)

   use constants
   use variables

   integer :: idummy

   if( idummy == 1 )then
     IsStored = .true.
     if( Debug > 0 ) write(*,*) 'Reading ',Variable(i)
   else
     IsStored = .false.
     if( Debug > 0 ) write(*,*) 'Header ',Variable(i)
   endif


end function IsStored
subroutine WriteRestartField_400
!
!  First run: when a variable is solved then store it
!  Additional (restart/initial) runs: when stored read it
!  and store it (even if solve has been turned off).
!
!  So 'solve' is always 'store', when restarted then
!  see if a variable is present then read and set 'store'.
!  However when restarted as an initial guess and a 
!  variable is not present (eg. when a flow has been solved 
!  and a scalar is later included) this routine should not
!  abort.
!
!  In a process like steady state => transient => 
!  back to steady state; during the last step we do not need
!  to store the old transient arrays
!
   use constants
   use geometry
   use variables

   integer IntStore

   if( Debug > -1 ) write(*,*) 'Opening restart file'
   call openfile(IOrst,Casename,'.rst','UNFORMATTED', &
                                       'SEQUENTIAL','UNKNOWN',Debug)

   write(IOrst) Version
   write(IOrst) Iter,Time
   write(IOrst) IPref, Pref, Tref
   
   write(IOrst) Ncel, Nbnd, Nfac, Nscal, Nmat
   write(IOrst) NVar, MaxVar, MaxSC, MaxMat   

   write(IOdbg,*) 'Write restart file'
   write(IOdbg,*) 'Version  :',Version
   write(IOdbg,*) 'Iter     :',Iter,Time
   write(IOdbg,*) 'Reference:',IPref, Pref, Tref
   write(IOdbg,*) 'Set1: ',Ncel, Nbnd, Nfac, Nscal, Nmat
   write(IOdbg,*) 'Set2: ',NVar, MaxVar, MaxSC, MaxMat   

   where( Solve ) Store = .true.

   write(IOdbg,*) 'Store:',Store(1:4),':',Store(5:8),'x',Store(9:MaxVar)
   if( SolveScalars ) write(IOdbg,*) 'Solve scalar:',Solve(Nvar+1:MaxVar)
   if( SolveScalars ) write(IOdbg,*) 'Store scalar:',Store(Nvar+1:MaxVar)      

   if( Debug > -1 ) write(*,*) 'Data'

   write(IOrst) VarU,IntStore(VarU,Store(VarU))
   if( Debug > -1 ) write(*,*) 'U'
   if( Store(VarU)    ) write(IOrst) U(:)
   
   if( Debug > -1 ) write(*,*) 'V'
   write(IOrst) VarV,IntStore(VarV,Store(VarV))
   if( Store(VarV)    ) write(IOrst) V(:)
   
   if( Debug > -1 ) write(*,*) 'W'
   write(IOrst) VarW,IntStore(VarW,Store(VarW))
   if( Store(VarW)    ) write(IOrst) W(:)
   
   if( Debug > -1 ) write(*,*) 'P'
   write(IOrst) VarP,IntStore(VarP,Store(VarP))
   if( Store(VarP)    ) write(IOrst) P(:)

   if( Debug > -1 ) write(*,*) 'KE'
   write(IOrst) VarTE,IntStore(VarTE,Store(VarTE))
   if( Store(VarTE)   ) write(IOrst) TE(:)

   if( Debug > -1 ) write(*,*) 'ED'
   write(IOrst) VarED,IntStore(VarED,Store(VarED))
   if( Store(VarED)   ) write(IOrst) ED(:)

   if( Debug > -1 ) write(*,*) 'T'
   write(IOrst) VarT,IntStore(VarT,Store(VarT))
   if( Store(VarT)    ) write(IOrst) T(:)

   if( Debug > -1 ) write(*,*) 'Density'
   write(IOrst) VarDEN,IntStore(VarDEN,Store(VarDEN))
   if( Store(VarDEN)  ) write(IOrst) Den(1:NCEL)

   if( Debug > -1 ) write(*,*) 'VisTurb',VarVIS,IntStore(VarVIS,Store(VarVIS))
   write(IOrst) VarVIS,IntStore(VarVIS,Store(VarVIS))
   if( Store(VarVis)  ) write(IOrst) VisEff(:)

!   if( Store(VarLvis) ) write(IOrst) VarLvis 
!   if( Store(VarLvis) ) write(IOrst) LamVisc(:)   
   
   if( Debug > -1 ) write(*,*) 'CP'
   write(IOrst) VarCP,IntStore(VarCP,Store(VarCP))
   if( Store(VarCP)   ) write(IOrst) CP(:)   

   do is=1,Nscal
     if( Debug > -1 ) write(*,*) 'is',is
     write(IOrst) VarS(is),IntStore(VarS(is),Store(Nvar+is))
     if( Store(Nvar+is) )then
       write(IOrst) Scalar(is)             ! store the name as well
       write(IOrst) SC(:,is)
     endif
   end do

   if( Debug > -1 ) write(*,*) 'Mass flux'
   write(IOrst) -3
   write(IOrst) MassFlux(:)
   write(IOrst) -4
   if( Debug > -1 ) write(*,*) 'Boundaries'
   write(IOrst) Bnd(:)
   

   if( Transient )then
     if( Debug > 0 ) write(*,*)'Restartfile (transient data)'
     write(IOrst) VarU,IntStore(VarU,Store(VarU))
     if( Store(VarU)  ) write(IOrst) Uold(:)
     write(IOrst) VarV,IntStore(VarV,Store(VarV))
     if( Store(VarV)  ) write(IOrst) Vold(:)
     write(IOrst) VarW,IntStore(VarW,Store(VarW))
     if( Store(VarW)  ) write(IOrst) Wold(:)

     write(IOrst) VarTE,IntStore(VarTE,Store(VarTE))
     if( Store(VarTE) ) write(IOrst) TEold(:)
     write(IOrst) VarED,IntStore(VarED,Store(VarED))
     if( Store(VarED) ) write(IOrst) EDold(:)  
     write(IOrst) VarT,IntStore(VarT,Store(VarT))
     if( Store(VarT)  ) write(IOrst) Told(:)  

!     if( Store(VarDEN)) write(IOrst) VarDen
!     if( Store(VarDEN)) write(IOrst) DENold(:)

     do is=1,Nscal
       if( Store(Nvar+is) )then
         write(IOrst) VarS(is),IntStore(VarS(is),Store(Nvar+is))
         write(IOrst) SCold(:,is)
       endif
     end do
     
     if( .not. Euler )then
       write(IOrst) .true.
       write(IOrst) VarU,IntStore(VarU,Store(VarU))
       if( Store(VarU)  ) write(IOrst) Uold2(:)
       write(IOrst) VarV,IntStore(VarV,Store(VarV))
       if( Store(VarV)  ) write(IOrst) Vold2(:)
       write(IOrst) VarW,IntStore(VarW,Store(VarW))
       if( Store(VarW)  ) write(IOrst) Wold2(:)

       write(IOrst) VarTE,IntStore(VarTE,Store(VarTE))
       if( Store(VarTE) ) write(IOrst) TEold2(:)
       write(IOrst) VarED,IntStore(VarED,Store(VarED))
       if( Store(VarED) ) write(IOrst) EDold2(:)  
       write(IOrst) VarT,IntStore(VarT,Store(VarT))
       if( Store(VarT)  ) write(IOrst) Told2(:)  

       do is=1,Nscal
         if( Store(Nvar+is) )then
           write(*,*) 'Storing very old scalar ',Scalar(is),is
           write(IOrst) VarS(is),IntStore(VarS(is),Store(Nvar+is))
           write(IOrst) SCold2(:,is)
         endif
       end do
     endif

   endif

   close(IOrst)
  
   if( Debug > -1 ) write(*,*)'Restartfile 400 written'
  
end subroutine WriteRestartField_400
subroutine ReadRestartField_400(IterStart)

   use constants
   use geometry
   use variables

   logical :: Ldummy, OK, IsStored

   logical, allocatable :: Stored(:), TmpLogical(:)

   character(len=12)    :: string
   
   !if( Debug > 0 ) write(*,*) 'Opening restart file'
   !call openfile(IOrst,Casename,'.rst','UNFORMATTED', &
   !                                    'SEQUENTIAL','UNKNOWN',Debug)
   !
   !read(IOrst) IRVersion
   !if( IRVersion /= Version )write(*,*)'+ Version mismatch'

 itmpdebug = Debug
 Debug = 6
   
   write(IOdbg,*) 'ReadRestartField vs400'

   read(IOrst) IterStart,Time
   read(IOrst) IPref, Pref, Tref
   read(IOrst) INcel,INbnd,INfac,INscal,INmat
   if( INcel  /= Ncel  ) write(*,*)'+ restartdata, Ncel  = ',INcel
   if( INbnd  /= Nbnd  ) write(*,*)'+ restartdata, Nbnd  = ',INbnd
   if( INbnd  /= Nbnd  ) Restart = 4
   if( INfac  /= Nfac  ) write(*,*)'+ restartdata, Nfac  = ',INfac
   if( INscal /= Nscal ) write(*,*)'+ restartdata, Nscal = ',INscal
   if( INmat  /= Nmat  ) write(*,*)'+ restartdata, Nmat  = ',INmat

   read(IOrst) INVar, IMaxVar, IMaxSC, IMaxMat
   if( INvar   /= Nvar  ) write(*,*)'+ restartdata, Nvar  = ',INvar  
   if( IMaxVar /= MaxVar) write(*,*)'+ restartdata, MaxVar= ',IMaxVar
   if( IMaxSC  /= MaxSC ) write(*,*)'+ restartdata, MaxSC = ',IMaxSC 
   if( IMaxMat /= MaxMat) write(*,*)'+ restartdata, MaxMat= ',IMaxMat

   allocate(Stored(IMaxVar),stat=istat)
   call TrackMemory(istat,IMaxVar,'Stored array allocated')

   Stored(:) = .false.

   !
   ! the data
   !
   read(IOrst) idata,iflag
   if( idata /= VarU ) write(*,*)'+ Corrupt restartdata, U = ',idata
   Stored(VarU) = IsStored(VarU,iflag)
   
   if( Stored(VarU) )then
     if( allocated( U ) .and. Restart == 4 )then
       read(IOrst) U(1:Ncel)       
     else if( allocated( U ) )then
       read(IOrst) U(:)       
     else
       write(*,*)'+ Internal error restartdata U' ! U should be alloc.
     endif
   else if( Solve(VarU) )then
     write(*,*)'+ U not found on restartdata'
     U(:) = Guess(VarU)
   endif

   read(IOrst) idata,iflag
   if( idata /= VarV ) write(*,*)'+ Corrupt restartdata, V = ',idata
   Stored(VarV) = IsStored(VarV,iflag)

   if( Stored(VarV) )then
     if( allocated( V ) .and. Restart == 4 )then
       read(IOrst) V(1:Ncel)       
     else if( allocated( V ) )then
       read(IOrst) V(:)   
     else    
       write(*,*)'+ Internal error restartdata V' ! V should be alloc.
     endif
   else if( Solve(VarV) )then
     write(*,*)'+ V not found on restartdata'
     V(:) = Guess(VarV)
   endif

   read(IOrst) idata,iflag
   if( idata /= VarW ) write(*,*)'+ Corrupt restartdata, W = ',idata
   Stored(VarW) = IsStored(VarW,iflag)

   if( Stored(VarW) )then
     if( allocated( W ) .and. Restart == 4 )then
       read(IOrst) W(1:Ncel)       
     else if( allocated( W ) )then
       read(IOrst) W(:)       
     else
       write(*,*)'+ Internal error restartdata W' ! W should be alloc.
     endif
   else if( Solve(VarW) )then
     write(*,*)'+ W not found on restartdata'
     W(:) = Guess(VarW)
   endif

   read(IOrst) idata,iflag
   if( idata /= VarP ) write(*,*)'+ Corrupt restartdata, P = ',idata
   Stored(VarP) = IsStored(VarP,iflag)

   if( Stored(VarP) )then
     if( allocated( P ) .and. Restart == 4 )then
       read(IOrst) P(1:Ncel)       
     else if( allocated( P ) )then
       read(IOrst) P(:)       
     else
       write(*,*)'+ Internal error restartdata P' ! P should be alloc.
     endif
   else if( Solve(VarP) )then
     write(*,*)'+ P not found on restartdata'
     P(:) = Guess(VarP)
   endif

   read(IOrst) idata,iflag
   if( idata /= VarTE ) write(*,*)'+ Corrupt restartdata, TE = ',idata
   Stored(VarTE) = IsStored(VarTE,iflag)

   if( Stored(VarTE) )then
     if( allocated( TE ) .and. Restart == 4 )then
       read(IOrst) TE(1:Ncel)       
     else if( allocated( TE ) )then
       read(IOrst) TE(:)       
     else
       write(*,*)'+ TE found on restartdata, TE not allocated'
       Store(VarTE) = .true.
       allocate(TE(Ncel+Nbnd),stat=istat)
       call TrackMemory(istat,Ncel+Nbnd,'Stored array TE allocated')
       read(IOrst) TE(:)       
     endif
   else if( Solve(VarTE) )then
     write(*,*)'+ TE not found on restartdata'
     TE(:) = Guess(VarTE)
   endif

   read(IOrst) idata,iflag
   if( idata /= VarED ) write(*,*)'+ Corrupt restartdata, ED = ',idata
   Stored(VarED) = IsStored(VarED,iflag)

   if( Stored(VarED) )then
     if( allocated( ED ) .and. Restart == 4 )then
       read(IOrst) ED(1:Ncel)       
     else if( allocated( ED ) )then
       read(IOrst) ED(:)       
     else
       write(*,*)'+ ED found on restartdata, ED not allocated'
       Store(VarED) = .true.
       allocate(ED(Ncel+Nbnd),stat=istat)
       call TrackMemory(istat,Ncel+Nbnd,'Stored array ED allocated')
       read(IOrst) ED(:)       
     endif
   else if( Solve(VarED) )then
     write(*,*)'+ ED not found on restartdata'
     ED(:) = Guess(VarED)
   endif

   read(IOrst) idata,iflag
   if( idata /= VarT ) write(*,*)'+ Corrupt restartdata, T = ',idata
   Stored(VarT) = IsStored(VarT,iflag)

   if( Stored(VarT) )then
     if( allocated( T ) .and. Restart == 4 )then
       read(IOrst) T(1:Ncel)       
     else if( allocated( T ) )then
       read(IOrst) T(:)       
     else
       write(*,*)'+ T found on restartdata, T not allocated'
       Store(VarT) = .true.
       allocate(T(Ncel+Nbnd),stat=istat)
       call TrackMemory(istat,Ncel+Nbnd,'Stored array T allocated')
       read(IOrst) T(:)       
     endif
   else if( Solve(VarT) )then
     write(*,*)'+ T not found on restartdata'
     T(:) = Guess(VarT)
   endif

   read(IOrst) idata,iflag
   if( idata /= VarDEN ) write(*,*)'+ Corrupt restartdata, DEN = ',idata
   Stored(VarDEN) = IsStored(VarDEN,iflag)

   if( Stored(VarDEN) )then
     if( allocated( DEN ) .and. Restart == 4 )then
       read(IOrst) DEN(1:Ncel)       
     else if( allocated( DEN ) )then
       read(IOrst) DEN(:)       
     else
       write(*,*)'+ Internal error restartdata DEN' ! DEN should be alloc.
     endif
   else if( Solve(VarDEN) )then
     write(*,*)'+ DEN not found on restartdata'
     DEN(:) = DensRef
   endif

   read(IOrst) idata,iflag
   if( idata /= VarVIS ) write(*,*)'+ Corrupt restartdata, VisEff = ',idata
   Stored(VarVIS) = IsStored(VarVIS,iflag)

   if( Stored(VarVis) )then
     if( allocated( VisEff ) .and. Restart == 4 )then
       Store(VarVis) = .true.          ! when stored then always store
       read(IOrst) VisEff(1:Ncel)       
     else if( allocated( VisEff ) )then
        Store(VarVis) = .true.          ! when stored then always store
       read(IOrst) VisEff(:)       
    else
       Store(VarVis) = .true.          ! when stored then always store
       read(IOrst) VisEff(:)       

       write(*,*)'+ Internal error restartdata VisEff' ! VisEff should be alloc.
     endif
   else if( Solve(VarVis) )then
     write(*,*)'+ VisEff not found on restartdata'
     VisEff(:) = VisLam
   endif

   read(IOrst) idata,iflag
   if( idata /= VarCP ) write(*,*)'+ Corrupt restartdata, CP = ',idata
   Stored(VarCP) = IsStored(VarCP,iflag)

   if( Stored(VarCP) )then
     if( allocated( CP ) .and. Restart == 4 )then
       read(IOrst) CP(1:Ncel)       
     else if( allocated( CP ) )then
       read(IOrst) CP(:)       
     else
       write(*,*)'+ Internal error restartdata CP' ! CP should be alloc?
       Store(VarCP) = .true.
       allocate(CP(Ncel+Nbnd),stat=istat)
       call TrackMemory(istat,Ncel+Nbnd,'Stored array CP allocated')
       read(IOrst) CP(:)       
     endif
   else if( Solve(VarCP) )then
     write(*,*)'+ CP not found on restartdata'
     CP(:) = CpStd
   endif

   !
   ! scalar arrays can grow in time...
   !
   if( INscal > Nscal )then
     if( allocated( SC ) )then
       deallocate( SC )
       call TrackMemory(istat,-(Ncel+Nbnd)*Nscal,'SC deallocated')

       allocate(SC(Ncel+Nbnd,INscal),stat=istat)
       call TrackMemory(istat,(Ncel+Nbnd)*INscal,'Stored array SC allocated (1)')

       write(*,*)'Found more scalars then expected!'
       write(*,*)'Growing arrays not implemented yet. Stop'
       stop

     else
       Nscal  = INscal
       MaxSC  = NScal
       MaxVar = NVar + MaxSC
       write(*,*)'Reset NScal, MaxSC, MaxVar to ',NScal,MaxSC,MaxVar
       
       allocate(SC(Ncel+Nbnd,Nscal),stat=istat)
       call TrackMemory(istat,(Ncel+Nbnd)*INscal,'Stored array SC allocated (2)')

       allocate(Scalar(MaxSC),stat=istat)
       call TrackMemory(istat,MaxSC,'Scalar array allocated')

       allocate(SolveScalar(MaxSC),stat=istat)
       call TrackMemory(istat,MaxSC,'SolveScalar array allocated')

       allocate(StoreScalar(MaxSC),stat=istat)
       call TrackMemory(istat,MaxSC,'StoreScalar array allocated')

       SolveScalar = .false.
       StoreScalar = .false.

       allocate(VarS(MaxSC),stat=istat)
       call TrackMemory(istat,MaxSC,'VarS array allocated')

       allocate(TmpLogical(MaxVar))
       
       TmpLogical(1:NVar) = Store(1:NVar)
       deallocate(Store)
       allocate(Store(MaxVar)) 
       Store(:) = .false. 
       Store(1:NVar) = TmpLogical(1:NVar) 

       TmpLogical(1:NVar) = Stored(1:NVar)
       deallocate(Stored)
       allocate(Stored(MaxVar))
       Stored(:) = .false. 
       Stored(1:NVar) = TmpLogical(1:NVar) 

       write(*,*) 'Found scalars, reading them anyway'
       write(*,*) 'Store:',Store(1:4),':',Store(5:8),'x',Store(9:MaxVar)
         
       do i=1,MaxSC
         indx    = Nvar + i
         VarS(i) = indx
       end do

     endif
     
   endif

   if( Nscal <= INscal )then
     do is=1,Nscal

       read(IOrst) idata,iflag
       if( idata /= VarS(is) ) write(*,*)'+ Corrupt restartdata, VarS = ',idata
       Stored(VarS(is)) = IsStored(VarS(is),iflag)

       if( Stored(Nvar+is) )then
         read(IOrst) string
         write(*,*)'Scalar found: ',string
         Scalar(is) = string

         read(IOrst) SC(:,is)       

         StoreScalar(is) = .true.
         Store(NVar+is)  = .true.
       else

         write(*,*)'+ Error scalar data not found ',is

         exit

       endif
     end do
   endif
   
   write(*,*)'Variables passed'
   
   read(IOrst) IN3
   if( IN3 /= -3 )    write(*,*)'+ Corrupt restartdata, MassFlux = ',IN3
   if( Restart == 4 )then
     read(IOrst) MassFlux(1)
     write(*,*)'+ massfluxes set to zero'
     MassFlux = 0.0
     call GradientPhi(VarU,U,dUdX)     
     call GradientPhi(VarV,V,dVdX)     
     call GradientPhi(VarW,W,dWdX)     
     call GradientPhi(VarP,P,dPdX)     
     do ib=1,Nbnd                              
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       if( it == ROutlet )then                 
         MassFlux(Bnd(ib)%face) = 1.0
       endif
     end do
     call FluxMass
   else
     read(IOrst) MassFlux(:)
   endif   
   write(*,*)'Massfluxes passed'

   read(IOrst) IN4
   if( IN4 /= -4 )    write(*,*)'+ Corrupt restartdata, Bnd = ',IN4
   if( Restart == 2 .or. Restart == 4 )then
     !
     ! skip the bnd-data
     !
     read(IOrst) Bnd(1)
   else
     read(IOrst) Bnd(:)
   endif
    
   write(*,*)'Boundaries passed'

   if( Transient .and. Restart == 2 )then
     !
     ! note Restart = 4 is de facto useless
     !
     write(*,*)'Restartfile as initial guess'

     Uold = U
     Vold = V
     Wold = W
     if( SolveTurbEnergy ) TEold = TE
     if( SolveTurbDiss )   EDold = ED 
     if( SolveEnthalpy )   Told  = T
     if( SolveScalars )then
       do is=1,Nscal
         SCold(:,is) = SC(:,is)
       end do
     endif

     if( .not. Euler )then
       write(*,*)'Setting quad arrays'
       Uold2 = U
       Vold2 = V
       Wold2 = W
       if( SolveTurbEnergy ) TEold2 = TE
       if( SolveTurbDiss )   EDold2 = ED 
       if( SolveEnthalpy )   Told2  = T
       if( SolveScalars )then
         do is=1,Nscal
           SCold2(:,is) = SC(:,is)
         end do
       endif
     endif
   else if( Transient .and. ( Restart == 1 .or. Restart == 3) )then
     write(*,*)'Restartfile (transient data)'

     read(IOrst) idata,iflag
     if( idata /= VarU ) write(*,*)'+ Corrupt restartdata, Uold = ',idata
     Stored(VarU) = IsStored(VarU,iflag)

     if( Stored(VarU) )then
       if( allocated( Uold ) )then
         read(IOrst) Uold(:)       
       else
         write(*,*)'+ Internal error restartdata Uold'
       endif
     else if( Solve(VarU) )then
       write(*,*)'+ Uold not found on restartdata'
       Uold(:) = U(:)
     endif

     read(IOrst) idata,iflag
     if( idata /= VarV ) write(*,*)'+ Corrupt restartdata, Vold = ',idata
     Stored(VarV) = IsStored(VarV,iflag)

     if( Stored(VarV) )then
       if( allocated( Vold ) )then
         read(IOrst) Vold(:)       
       else
         write(*,*)'+ Internal error restartdata Vold'
       endif
     else if( Solve(VarV) )then
       write(*,*)'+ Vold not found on restartdata'
       Vold(:) = V(:)
     endif

     read(IOrst) idata,iflag
     if( idata /= VarW ) write(*,*)'+ Corrupt restartdata, Wold = ',idata
     Stored(VarW) = IsStored(VarW,iflag)

     if( Stored(VarW) )then
       if( allocated( Wold ) )then
         read(IOrst) Wold(:)       
       else
         write(*,*)'+ Internal error restartdata Wold'
       endif
     else if( Solve(VarW) )then
       write(*,*)'+ Wold not found on restartdata'
       Wold(:) = W(:)
     endif

     read(IOrst) idata,iflag
     if( idata /= VarTE ) write(*,*)'+ Corrupt restartdata, TEold = ',idata
     Stored(VarTE) = IsStored(VarTE,iflag)

     if( Stored(VarTE) )then
       read(IOrst) TEold(:)
     endif 
      
     read(IOrst) idata,iflag
     if( idata /= VarED ) write(*,*)'+ Corrupt restartdata, EDold = ',idata
     Stored(VarED) = IsStored(VarED,iflag)

     if( Stored(VarED) )then
       read(IOrst) EDold(:)  
     endif 
     
     read(IOrst) idata,iflag
     if( idata /= VarT ) write(*,*)'+ Corrupt restartdata, Told = ',idata
     Stored(VarT) = IsStored(VarT,iflag)

     if( Stored(VarT) )then
       read(IOrst) Told(:)  
     endif  

     write(*,*)'Restartfile (transient data, old arays)'
     if( .not. Euler )then

       OK = .false.
       read(IOrst,end=4) Ldummy
       OK = .true.

 4     continue
       if( .not. OK )then
         write(*,*)'Old2 arrays not on restart file'
       endif

       if( Ldummy )then
         write(*,*) 'Restartfile (transient data, old2 arrays)'
         
         read(IOrst) idata,iflag
         if( idata /= VarU ) write(*,*)'+ Corrupt restartdata, Uold2 = ',idata
         Stored(VarU) = IsStored(VarU,iflag)

         if( Stored(VarU) )then
           if( allocated( Uold2 ) )then
             read(IOrst) Uold2(:)       
           else
             write(*,*)'+ Internal error restartdata Uold2'
           endif
         else if( Solve(VarU) )then
           write(*,*)'+ Uold2 not found on restartdata'
           Uold2(:) = Uold(:)
         endif

         read(IOrst) idata,iflag
         if( idata /= VarV ) write(*,*)'+ Corrupt restartdata, Vold2 = ',idata
         Stored(VarV) = IsStored(VarV,iflag)

         if( Stored(VarV) )then
           if( allocated( Vold2 ) )then
             read(IOrst) Vold2(:)       
           else
             write(*,*)'+ Internal error restartdata Vold2'
           endif
         else if( Solve(VarV) )then
           write(*,*)'+ Vold2 not found on restartdata'
           Vold2(:) = Vold(:)
         endif

         read(IOrst) idata,iflag
         if( idata /= VarW ) write(*,*)'+ Corrupt restartdata, Wold2 = ',idata
         Stored(VarW) = IsStored(VarW,iflag)

         if( Stored(VarW) )then
           if( allocated( Wold2 ) )then
             read(IOrst) Wold2(:)       
           else
             write(*,*)'+ Internal error restartdata Wold2'
           endif
         else if( Solve(VarW) )then
           write(*,*)'+ Wold2 not found on restartdata'
           Wold2(:) = Wold(:)
         endif

         read(IOrst) idata,iflag
         if( idata /= VarTE ) write(*,*)'+ Corrupt restartdata, TEold2 = ',idata
         Stored(VarTE) = IsStored(VarTE,iflag)

         if( Stored(VarTE) )then
           if( allocated( TEold2 ) )then
             read(IOrst) TEold2(:)       
           else
             write(*,*)'+ Internal error restartdata TEold2'
           endif
         else if( Solve(VarTE) )then
           write(*,*)'+ TEold2 not found on restartdata'
           TEold2(:) = TEold(:)
         endif

         read(IOrst) idata,iflag
         if( idata /= VarED ) write(*,*)'+ Corrupt restartdata, EDold2 = ',idata
         Stored(VarED) = IsStored(VarED,iflag)

         if( Stored(VarED) )then
           if( allocated( EDold2 ) )then
             read(IOrst) EDold2(:)       
           else
             write(*,*)'+ Internal error restartdata EDold2'
           endif
         else if( Solve(VarED) )then
           write(*,*)'+ EDold2 not found on restartdata'
           EDold2(:) = EDold(:)
         endif

         read(IOrst) idata,iflag
         if( idata /= VarT ) write(*,*)'+ Corrupt restartdata, Told2 = ',idata
         Stored(VarT) = IsStored(VarT,iflag)

         if( Stored(VarT) )then
           if( allocated( Told2 ) )then
             read(IOrst) Told2(:)       
           else
             write(*,*)'+ Internal error restartdata Told2'
           endif
         else if( Solve(VarT) )then
           write(*,*)'+ Told2 not found on restartdata'
           Told2(:) = Told(:)
         endif

       else

         Uold2 = Uold
         Vold2 = Vold
         Wold2 = Wold
         if( SolveTurbEnergy ) TEold2 = TEold
         if( SolveTurbDiss )   EDold2 = EDold 
         if( SolveEnthalpy )   Told2  = Told

       endif
     endif
   endif

   close(IOrst)
   
   !
   ! overview in debug file
   !
   write(IOdbg,*)'Restartfile read from ',IterStart,Time

   write(IOdbg,*)'U:',minval(U(1:Ncel)),'< U <',maxval(U(1:Ncel)), &
                                  sum( U(1:Ncel) )/float(Ncel)     
   write(IOdbg,*)'V:',minval(V(1:Ncel)),'< V <',maxval(V(1:Ncel)), &
                                  sum( V(1:Ncel) )/float(Ncel)     
   write(IOdbg,*)'W:',minval(W(1:Ncel)),'< W <',maxval(W(1:Ncel)), &
                                  sum( W(1:Ncel) )/float(Ncel)     
   write(IOdbg,*)'P:',minval(P(1:Ncel)),'< P <',maxval(P(1:Ncel)), &
                                  sum( P(1:Ncel) )/float(Ncel)     
   write(IOdbg,*)'D:',minval(Den(1:Ncel)),'< D <',maxval(Den(1:Ncel)), &
                                  sum( Den(1:Ncel) )/float(Ncel)     

   if( SolveTurb )then
     write(IOdbg,*)'k:',minval(TE(1:Ncel)),'< k <', &
                   maxval(TE(1:Ncel)), sum( TE(1:Ncel) )/float(Ncel)     
     write(IOdbg,*)'e:',minval(ED(1:Ncel)),'< e <', &
                   maxval(ED(1:Ncel)), sum( ED(1:Ncel) )/float(Ncel)     
   endif
   write(IOdbg,*)'v:',minval(VisEff(1:Ncel)),'< v <', &
                 maxval(VisEff(1:Ncel)), sum( VisEff(1:Ncel) )/float(Ncel)     

   if( SolveEnthalpy ) write(IOdbg,*)'T:',minval(T(1:Ncel)),'< T <', &
                      maxval(T(1:Ncel)), sum( T(1:Ncel) )/float(Ncel)     

   do i=1,Nscal
     write(IOdbg,*)'S:',minval(Sc(1:Ncel,i)),'< S <', &
           maxval(Sc(1:Ncel,i)), sum(Sc(1:Ncel,i))/float(Ncel),':',i    
   end do
  
   if( Debug > 0 ) write(*,*)'Restartfile read from ',IterStart,Time

  Debug = itmpdebug
 
end subroutine ReadRestartField_400
