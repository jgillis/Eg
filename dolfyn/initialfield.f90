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
subroutine InitialField(dppmax)
!========================================================================
!     This routine assembles and solves the pressure-correction
!     equation using colocated grid. SIMPLE algorithm with one
!     or more corrector steps (non-orthogonality effects taken
!     into account as described in Sect. 8.8 of the corrected
!     2nd printing.
!
   use constants
   use geometry
   use variables
   use scalars

!  logical :: Bernoulli = .false.

   if( Restart == 3 ) return 
   
   if( Debug > 3 ) write(*,*)'*** InitialField'

   GammaSet    = Gamma(VarU)  
   Gamma(VarU) = 0.0

   iDebugSet = debug
   debug = 0

   Initialisation = .true.

   if( Restart == 0 .and. .not. UserInitialisation )then
   
     write(*,*)'InitialField (1)'

     VisEff = VisInitialFactor * VisLam

     dppmax = 0.0
     do iter=1,InitialSteps

         call GradientPhi(VarU,U,dUdX)     
         call GradientPhi(VarV,V,dVdX)     
         call GradientPhi(VarW,W,dWdX)     
         call GradientPhi(VarP,P,dPdX)     

         if( SolveUVW) call CalculateUVW 

         if( SolveP )  call CalculatePressure(dppmax)
         write(*,*)'>> ',iter,dppmax 

         call UpdateBC

         !if( SolveTurbEnergy ) &
         !  call CalculateScalar(VarTE,TE,dpdx,TEold,TEold2)
         !if( SolveTurbDiss )   &
         !  call CalculateScalar(VarED,ED,dpdx,EDold,EDold2)

         call PrintSummary(IOdbg,-iter)

     end do

     VisEff = VisEff / VisInitialFactor 

     !
     ! now there is an initial velocity and pressure field
     ! set other variables if neccesary
     !
     ! the pressure is based on a higher viscosity, but as
     ! all the walls were frictionless the resulting pressure
     ! drop might be moderate
     !
     !if( Bernoulli )then
     !  Pref  = P(IPref)
     !  Vref2 = U(IPref)**2 + V(IPref)**2 + W(IPref)**2
     !
     !  do i=1,Ncel
     !    Vmag2 = U(i)**2 + V(i)**2 + W(i)**2
     !    P(i)  = Pref + 0.5*Den(i)*(Vref2 - Vmag2)
     !  end do
     !endif

     if( SolveTurb )then

       write(*,*)'InitialField (1), SolveTurb'
       !
       ! assume a low turbulence intensity
       !
       !ti2 = 0.01 **2
       
       !
       ! test if something is set as a guess
       !
       !sumte = sum( TE(1:Ncel) )

       do ip=1,Ncel

         !Vmag2 = U(ip)**2 + V(ip)**2 + W(ip)**2
         !TE(ip) = 1.5 * ti2 * Vmag2
         !ED(ip) = TE(ip)**(1.5) / TMLenSc + Small
       
         TE(ip) = Guess(VarTE)
         ED(ip) = Guess(VarED)
       
         VisEff(ip) = VisLam + DensRef * TMCmu * &
                               TE(ip)**2/(ED(ip)+Small)
       end do

       write(*,*)'k:',minval(TE(1:Ncel)),'< k <', &
                     maxval(TE(1:Ncel)), sum( TE(1:Ncel) )/float(Ncel)     
       write(*,*)'e:',minval(ED(1:Ncel)),'< e <', &
                     maxval(ED(1:Ncel)), sum( ED(1:Ncel) )/float(Ncel)     
       write(*,*)'v:',minval(VisEff(1:Ncel)),'< v <', &
                     maxval(VisEff(1:Ncel)), sum( VisEff(1:Ncel) )/float(Ncel)     
                     
     else

       VisEff = VisLam

     endif                  
   
     if( SolveEnthalpy ) T = Guess(VarT)
                    
 !  elseif( .not. SolveTurb .and. .not. UserInitialisation .and. Restart == 0 )then
 !  
 !    write(*,*)'InitialField (2)'
 !
 !    VisEff = VisLam  
 !    if( SolveEnthalpy ) T = Guess(VarT)
   
   elseif( SolveTurb .and. UserInitialisation .and. Restart == 0 )then

     do ip=1,Ncel

       Utmp   = Guess(VarU)
       Vtmp   = Guess(VarV)
       Wtmp   = Guess(VarW)
       Ptmp   = Guess(VarP)
       TEtmp  = Guess(VarTE)
       EDtmp  = Guess(VarED)
       Ttmp   = Guess(VarT)
       Dtmp   = DensRef
        
       Xtmp   = Cell(ip)%x(1)
       Ytmp   = Cell(ip)%x(2)
       Ztmp   = Cell(ip)%x(3)
       VOLtmp = Cell(ip)%vol
       
       ! write(*,*)'1:',utmp,vtmp,wtmp,ptmp
       ! write(*,*)'2:',tetmp,edtmp,ttmp,dtmp

       call UserInitialField(ip,Xtmp,Ytmp,Ztmp,VOLtmp, &
                       Pref,Tref,Iter,Time,            &
                     Utmp,Vtmp,Wtmp,Ptmp,TEtmp,EDtmp,Ttmp,Dtmp)
        
       ! write(*,*)'3:',utmp,vtmp,wtmp,ptmp
       ! write(*,*)'4:',tetmp,edtmp,ttmp,dtmp

       U(ip)  = Utmp
       V(ip)  = Vtmp
       W(ip)  = Wtmp
       P(ip)  = Ptmp
       if( SolveTurb )     TE(ip) = TEtmp
       if( SolveTurb )     ED(ip) = EDtmp
       if( SolveEnthalpy ) T(ip)  = Ttmp
       Den(ip)= Dtmp
       
       if( SolveTurb )then
         VisEff(ip) = VisLam + Dtmp * TMCmu * &
                              TEtmp**2/(EDtmp+Small)
       else
         VisEff(ip) = VisLam
       endif                  
                       
       ! write(*,*)'5:',utmp,vtmp,wtmp,ptmp
       ! write(*,*)'6:',tetmp,edtmp,ttmp,dtmp
     end do       
     
999  call FluxMass
     
     write(*,*)'User initialisation done'

     VisEff = VisInitialFactor * VisLam
     dppmax = 0.0
     do iter=1,InitialSteps
         call GradientPhi(VarU,U,dUdX)     
         call GradientPhi(VarV,V,dVdX)     
         call GradientPhi(VarW,W,dWdX)     
         call GradientPhi(VarP,P,dPdX)     
         if( SolveUVW) call CalculateUVW 
         if( SolveP )  call CalculatePressure(dppmax)
         write(*,*)'>> ',iter,dppmax 
         call UpdateBC
         call PrintSummary(IOdbg,-iter)
     end do
     VisEff = VisEff / VisInitialFactor 
   
   elseif( SolveTurb .and. .not. UserInitialisation .and. Restart == 2 )then
   
     write(*,*)'InitialField (4)'

     if( ED(1) < 0.0 )then
       !
       ! invalid epsilon value found (set in ReadRestartField)
       !
       TKguess  = Guess(VarTE)
       EDguess  = Guess(VarED)
       
       !EDguess  = TKguess**(1.5) / TMLenSc + Small
     
       VISguess = VisLam + DensRef * TMCmu * &
                  TKguess **2/(EDguess +Small)
   
       write(*,*)'Set initial turbulence variables'
       write(*,*)'Turb. length scale  : ',TMLenSc 
       write(*,*)'Turb. kinetic energy: ',TKguess
       write(*,*)'Turb. dissipation   : ',EDguess
       write(*,*)'Effective viscosity : ',VISguess

       TE(1:Ncel) = TKguess
       ED(1:Ncel) = EDguess
       VisEff(1:Ncel) = VISguess
       
     endif
     
     write(*,*)'Simple override temperature array???'    ! TO DO BETTER SOLUTION!

     !if( SolveEnthalpy ) T = Guess(VarT)

     write(*,*)'Simple override density array???'        ! TO DO BETTER SOLUTION!

     DEN = DensRef
   
   endif
   
   Gamma(VarU) = GammaSet
   debug = iDebugSet
      
   Initialisation = .false.

99 continue
   
   if( Debug > 3 ) write(*,*)'=== InitialField'
   
end subroutine InitialField
subroutine WriteRestartField_300

   use constants
   use geometry
   use variables

   if( Debug > -1 ) write(*,*) 'Opening restart file'
   call openfile(IOrst,Casename,'.rst','UNFORMATTED', &
                                       'SEQUENTIAL','UNKNOWN',Debug)

   write(IOrst) Version
   write(IOrst) Iter,Time
   write(IOrst) IPref, Pref, Tref
   write(IOrst) Ncel,Nbnd,Nfac

   write(IOrst) VarU
   write(IOrst) U(:)
   write(IOrst) VarV
   write(IOrst) V(:)
   write(IOrst) VarW
   write(IOrst) W(:)
   write(IOrst) VarP
   write(IOrst) P(:)
   write(IOrst) -1
   write(IOrst) Den(:)
   write(IOrst) -2
   write(IOrst) VisEff(:)
   write(IOrst) -3
   write(IOrst) MassFlux(:)
   write(IOrst) -4
   write(IOrst) Bnd(:)
   
   if( SolveTurbEnergy )then
     write(IOrst) VarTE
     write(IOrst) TE(:)
   endif 
   if( SolveTurbDiss )then
     write(IOrst) VarED
     write(IOrst) ED(:)  
   endif 
   if( SolveEnthalpy )then
     write(IOrst) VarT
     write(IOrst) T(:)  
   endif  

   if( Transient )then
     if( Debug > 0 ) write(*,*)'Restartfile (transient data)'
     write(IOrst) VarU
     write(IOrst) Uold(:)
     write(IOrst) VarV
     write(IOrst) Vold(:)
     write(IOrst) VarW
     write(IOrst) Wold(:)
     if( SolveTurbEnergy )then
       write(IOrst) VarTE
       write(IOrst) TEold(:)
     endif 
     if( SolveTurbDiss )then
       write(IOrst) VarED
       write(IOrst) EDold(:)  
     endif 
     if( SolveEnthalpy )then
       write(IOrst) VarT
       write(IOrst) Told(:)  
     endif  
     
     if( .not. Euler )then
       write(IOrst) .true.
       write(IOrst) VarU
       write(IOrst) Uold2(:)
       write(IOrst) VarV
       write(IOrst) Vold2(:)
       write(IOrst) VarW
       write(IOrst) Wold2(:)
       if( SolveTurbEnergy )then
         write(IOrst) VarTE
         write(IOrst) TEold2(:)
       endif 
       if( SolveTurbDiss )then
         write(IOrst) VarED
         write(IOrst) EDold2(:)  
       endif 
       if( SolveEnthalpy )then
         write(IOrst) VarT
         write(IOrst) Told2(:)  
       endif  
     endif

   endif

   close(IOrst)
  
   if( Debug > 0 ) write(*,*)'Restartfile (300) written'
  
end subroutine WriteRestartField_300
subroutine ReadRestartField_300(IterStart)

   use constants
   use geometry
   use variables

   logical :: Ldummy, OK

   !if( Debug > 0 ) write(*,*) 'Opening restart file'
   !call openfile(IOrst,Casename,'.rst','UNFORMATTED', &
   !                                   'SEQUENTIAL','UNKNOWN',Debug)
   !
   !read(IOrst) IRVersion
   !if( IRVersion /= Version )write(*,*)'+ Version mismatch'
   
   read(IOrst) IterStart,Time
   read(IOrst) IPref, Pref, Tref
   read(IOrst) INcel,INbnd,INfac
   if( INcel /= Ncel )write(*,*)'+ Corrupt restartdata, Ncel = ',INcel
   if( INbnd /= Nbnd )write(*,*)'+ Corrupt restartdata, Nbnd = ',INbnd
   if( INfac /= Nfac )write(*,*)'+ Corrupt restartdata, Nfac = ',INfac

   read(IOrst) IVarU
   if( IVarU /= VarU )write(*,*)'+ Corrupt restartdata, ck01 = ',IVarU
   read(IOrst) U(:)
   
   read(IOrst) IVarV
   if( IVarV /= VarV )write(*,*)'+ Corrupt restartdata, ck02 = ',IVarV
   read(IOrst) V(:)
   
   read(IOrst) IVarW
   if( IVarW /= VarW )write(*,*)'+ Corrupt restartdata, ck03 = ',IVarW
   read(IOrst) W(:)
   
   read(IOrst) IVarP
   if( IVarP /= VarP )write(*,*)'+ Corrupt restartdata, ck04 = ',IVarP
   read(IOrst) P(:)
   
   read(IOrst) IN1
   if( IN1 /= -1 )    write(*,*)'+ Corrupt restartdata, ck05 = ',IN1
   read(IOrst) Den(:)
   
   read(IOrst) IN2
   if( IN2 /= -2 )    write(*,*)'+ Corrupt restartdata, ck06 = ',IN2
   read(IOrst) VisEff(:)
   
   read(IOrst) IN3
   if( IN3 /= -3 )    write(*,*)'+ Corrupt restartdata, ck07 = ',IN3
   read(IOrst) MassFlux(:)
   
   read(IOrst) IN4
   if( IN4 /= -4 )    write(*,*)'+ Corrupt restartdata, ck08 = ',IN4
   if( Restart == 2 )then
     !
     ! skip the bnd-data
     !
     read(IOrst) Bnd(1)
   else
     read(IOrst) Bnd(:)
   endif
    
   if( SolveTurbEnergy )then
     OK = .false.
     read(IOrst,end=1) IVarTE
     if( IVarTE /= VarTE )write(*,*)'+ Corrupt restartdata, ck09 = ',IVarTE
     read(IOrst) TE(:)
     OK = .true.

 1   continue
     if( .not. OK )then
       write(*,*)'TE not on restart file'
       TE(:) = Guess(VarTE)
      endif
   endif 
   
   if( SolveTurbDiss )then
     OK = .false.
     read(IOrst,end=2) IVarED
     if( IVarED /= VarED )write(*,*)'+ Corrupt restartdata, ck10 = ',IVarED
     read(IOrst) ED(:)  
     OK = .true.

 2   continue
     if( .not. OK )then
       write(*,*)'ED not on restart file'
       
       ED(:) = 0.0
       
       ! flag the first entry INVALID
       !
       ED(1) = -1.0
       
     endif
   endif 

   if( SolveEnthalpy )then
     OK = .false.
     read(IOrst,end=3) IVarT
     if( IVarT /= VarT )write(*,*)'+ Corrupt restartdata, ck11 = ',IVarT
     read(IOrst) T(:)  
     OK = .true.

 3   continue
     if( .not. OK )then
       write(*,*)'T  not on restart file'
     endif
   endif  

   if( Transient .and. Restart == 2 )then
     write(*,*)'Restartfile as initial guess'

     Uold = U
     Vold = V
     Wold = W
     if( SolveTurbEnergy ) TEold = TE
     if( SolveTurbDiss )   EDold = ED 
     if( SolveEnthalpy )   Told  = T

     if( .not. Euler )then
       write(*,*)'Setting quad arrays'
       Uold2 = U
       Vold2 = V
       Wold2 = W
       if( SolveTurbEnergy ) TEold2 = TE
       if( SolveTurbDiss )   EDold2 = ED 
       if( SolveEnthalpy )   Told2  = T
     endif
   else if( Transient .and. ( Restart == 1 .or. Restart == 3) )then
     write(*,*)'Restartfile (transient data)'

     read(IOrst) IVarU
     if( IVarU /= VarU )write(*,*)'+ Corrupt restartdata, ck12 = ',IVarU
     read(IOrst) Uold(:)
     
     read(IOrst) IVarV
     if( IVarV /= VarV )write(*,*)'+ Corrupt restartdata, ck13 = ',IVarV
     read(IOrst) Vold(:)
     
     read(IOrst) IVarW
     if( IVarW /= VarW )write(*,*)'+ Corrupt restartdata, ck14 = ',IVarW
     read(IOrst) Wold(:)
     
     if( SolveTurbEnergy )then
       read(IOrst) IVarTE
       if( IVarTE /= VarTE )write(*,*)'+ Corrupt restartdata, ck15 = ',IVarTE
       read(IOrst) TEold(:)
      endif 
      
     if( SolveTurbDiss )then
       read(IOrst) IVarED
       if( IVarED /= VarED )write(*,*)'+ Corrupt restartdata, ck16 = ',IVarED
       read(IOrst) EDold(:)  
     endif 
     
     if( SolveEnthalpy )then
       read(IOrst) IVarT
       if( IVarT /= VarT )write(*,*)'+ Corrupt restartdata, ck17 = ',IVarT
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
         
         read(IOrst) IVarU
         if( IVarU /= VarU )write(*,*)'+ Corrupt restartdata, ck18 = ',IVarU
         read(IOrst) Uold2(:)

         read(IOrst) IVarV
         if( IVarV /= VarV )write(*,*)'+ Corrupt restartdata, ck19 = ',IVarV
         read(IOrst) Vold2(:)

         read(IOrst) IVarW
         if( IVarW /= VarW )write(*,*)'+ Corrupt restartdata, ck20 = ',IVarW
         read(IOrst) Wold2(:)

         if( SolveTurbEnergy )then
           read(IOrst) IVarTE
           if( IVarTE /= VarTE )write(*,*)'+ Corrupt restartdata, ck21 = ',IVarTE
           read(IOrst) TEold2(:)
         endif 

         if( SolveTurbDiss )then
           read(IOrst) IVarED
           if( IVarED /= VarED )write(*,*)'+ Corrupt restartdata, ck22 = ',IVarED
           read(IOrst) EDold2(:)  
         endif 

         if( SolveEnthalpy )then
           read(IOrst) IVarT
           if( IVarT /= VarT )write(*,*)'+ Corrupt restartdata, ck23 = ',IVarT
           read(IOrst) Told2(:)  
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
  
   if( Debug > 0 ) write(*,*)'Restartfile read from ',IterStart,Time
  
end subroutine ReadRestartField_300
subroutine PrintSummary(IO,iouter)

   use constants
   use geometry
   use variables
   use scalars

   real, dimension(20) :: SCMON, ZEROS = 0.0         !<===== tijdelijk!

   iprint = Iter
   
   TEmon = 0.0
   EDmon = 0.0
   Tmon  = Tref
   
   Umon  = U(IMoni)
   Vmon  = V(IMoni) 
   Wmon  = W(Imoni)
   Pmon  = P(IMoni)
   if( SolveTurb )     TEmon = TE(IMoni)
   if( SolveTurb )     EDmon = ED(Imoni)
   if( SolveEnthalpy ) Tmon  = T(IMoni)
   if( SolveScalars)   SCMON(1:MaxSC) = SC(IMoni,1:MaxSC)
   
   if( Debug > -1 )then
     !write(IO,'(i6,7(1x,1pe9.2))') iprint,Residual(1:7)
     !write(IO,'(6x,7(1x,1pe9.2))') Umon,Vmon,Wmon,Pmon,TEmon,EDmon,Tmon
     if( Transient .and. Debug < 1 )then
       write(IO,'(1x,i6,i4,1x,1pe10.3,'':'',7(1x,1pe9.2),3x,7(1x,1pe9.2))') &
         iprint,iouter,Time,Residual(1:7), &
         Umon,Vmon,Wmon,Pmon,TEmon,EDmon,Tmon-Tref
     else if( Transient )then
       write(IO,'(1x,i6,i4,1x,1pe10.3,'':'',7(1x,1pe9.2))') &
         iprint,iouter,Time,Residual(1:7)
       write(IO,'(23x,7(1x,1pe9.2))') &
         Umon,Vmon,Wmon,Pmon,TEmon,EDmon,Tmon-Tref
     else if( Debug < 1 )then
!      write(IO,'(i6,'':'',7(1x,1pe9.2),3x,7(1x,1pe9.2))') &
!        iprint,Residual(1:7),Umon,Vmon,Wmon,Pmon,TEmon,EDmon,Tmon-Tref
       write(IO,'(5a1,'':'',i6,'':'',7(1x,1pe9.2),3x,7(1x,1pe9.2))') &
         (Flags(i),i=1,NFlags), &
         iprint,Residual(1:7),Umon,Vmon,Wmon,Pmon,TEmon,EDmon,Tmon-Tref
       if( SolveScalars )then
       write(IO,'(7x,7(1x,1pe9.2),3x,7(1x,1pe9.2))') &
         Residual(NVar+1:NVar+MaxSC),Zeros(MaxSC+1:7), &     ! <== niet algemeen genoeg!
         SCMON(1:MaxSC),Zeros(MaxSC+1:7)
       endif
     else
       write(IO,'(i6,'':'',7(1x,1pe9.2))') &
         iprint,Residual(1:7) 
       write(IO,'(7x,7(1x,1pe9.2))') &
         Umon,Vmon,Wmon,Pmon,TEmon,EDmon,Tmon-Tref
       if( SolveScalars )then
         write(IO,'(5x,''s:'',7(1x,1pe9.2))') &
           Residual(NVar+1:NVar+MaxSC)                      ! <== niet algemeen genoeg!
         write(IO,'(5x,''s:'',7(1x,1pe9.2))') &
           SCMON(1:MaxSC)                      
       endif
     endif

   endif
   
   !
   ! residuals in a separate file (once!):
   !
   if( IO == IOdef )then
     if( .not. Transient )then
       write(IOres,'(i6,7(1x,1pe9.2))') iprint,Residual(1:7)
       write(IOmon,'(i6,7(1x,1pe9.2))') iprint, &
                                  Umon,Vmon,Wmon,Pmon,TEmon,EDmon,Tmon-Tref
     else
       write(IOres,'(1pe12.5,7(1x,1pe12.5))') Time,Residual(1:7)
       write(IOmon,'(1pe12.5,7(1x,1pe12.5))') Time, &
                                  Umon,Vmon,Wmon,Pmon,TEmon,EDmon,Tmon-Tref
     endif		  			 
   endif

   if( ( Debug > 1 .and. IO == IOdef ) .or. &
       ( IO == IOdbg ) )then
     write(IO,'(//)')
     write(IO,*)'U:',minval(U(1:Ncel)),'< U <',maxval(U(1:Ncel)), &
                                    sum( U(1:Ncel) )/float(Ncel)     
     write(IO,*)'V:',minval(V(1:Ncel)),'< V <',maxval(V(1:Ncel)), &
                                    sum( V(1:Ncel) )/float(Ncel)     
     write(IO,*)'W:',minval(W(1:Ncel)),'< W <',maxval(W(1:Ncel)), &
                                    sum( W(1:Ncel) )/float(Ncel)     
     write(IO,*)'P:',minval(P(1:Ncel)),'< P <',maxval(P(1:Ncel)), &
                                    sum( P(1:Ncel) )/float(Ncel)     
     write(IO,*)'D:',minval(Den(1:Ncel)),'< D <',maxval(Den(1:Ncel)), &
                                    sum( Den(1:Ncel) )/float(Ncel)     

     if( SolveTurb )then
       write(IO,*)'k:',minval(TE(1:Ncel)),'< k <', &
                     maxval(TE(1:Ncel)), sum( TE(1:Ncel) )/float(Ncel)     
       write(IO,*)'e:',minval(ED(1:Ncel)),'< e <', &
                     maxval(ED(1:Ncel)), sum( ED(1:Ncel) )/float(Ncel)     
     endif
     write(IO,*)'v:',minval(VisEff(1:Ncel)),'< v <', &
                   maxval(VisEff(1:Ncel)), sum( VisEff(1:Ncel) )/float(Ncel)     

     if( SolveEnthalpy ) write(IO,*)'T:',minval(T(1:Ncel)),'< T <', &
                        maxval(T(1:Ncel)), sum( T(1:Ncel) )/float(Ncel)     

     if( SolveScalars )then
       do i=1,Nscal
        write(IO,*)'S:',minval(Sc(1:Ncel,i)),'< S <', &
                     maxval(Sc(1:Ncel,i)), sum(Sc(1:Ncel,i))/float(Ncel),':',i    
       end do
     endif

   endif
   
end subroutine PrintSummary
subroutine Disclaimer

   use constants

   write(IOdbg,'(1x,A,f5.3)') 'This is dolfyn version ',float(Version)*0.001
   write(IOdbg,'(1x,A)') 'Copyright(C) 2002-2007 Cyclone Fluid Dynamics BV'
   write(IOdbg,'(1x,A)') 'NL-5583 XM, Waalre, The Netherlands'
   write(IOdbg,'(1x,A)') 'see http://www.cyclone.nl and http://www.dolfyn.net'
   write(IOdbg,'(/,1x,A)') 'Using Sparsekit2 by Yousef Saad'
   write(IOdbg,'(1x,A)') '(c) 2005, the Regents of the University of Minnesota'
   write(IOdbg,'(/,1x,A)') 'Modules with patches (C) 2004-2007 by B. Tuinstra'
   write(IOdbg,'(1x,A)') 'see http://www.home.zonnet.nl/bouke_1/dolfyn'
   write(IOdbg,'(/,1x,A)') 'Modules with particles (C) 2006-2007 by'
   write(IOdbg,'(1x,A)')   'Henk Krus and Shibo (Harry) Kuang'
   write(IOdbg,'(/,1x,A)') 'Tecplot interface (tecplt.f90) (c) 2006-2007 by'
   write(IOdbg,'(1x,A)')   'Shibo (Harry) Kuang'
   write(IOdbg,'(1x,A)')   '(University of New South Wales, Sydney, Australia)'
   write(IOdbg,'(/,1x,A)') 'Particles in VTK output (C) 2007 by J. Jacobs'

   write(IOdbg,'(/)')
   write(IOdbg,'(1x,A)')'----------------------------------------------------'
   write(IOdbg,'(1x,A)')'Disclaimer of Warranty. Unless required by          '
   write(IOdbg,'(1x,A)')'applicable law or agreed to in writing, Licensor    '
   write(IOdbg,'(1x,A)')'provides the Work (and each Contributor provides its'
   write(IOdbg,'(1x,A)')'Contributions) on an "AS IS" BASIS, WITHOUT         '
   write(IOdbg,'(1x,A)')'WARRANTIES OR CONDITIONS OF ANY KIND, either express'
   write(IOdbg,'(1x,A)')'or implied, including, without limitation, any      '
   write(IOdbg,'(1x,A)')'warranties or conditions of TITLE, NON-INFRINGEMENT,'
   write(IOdbg,'(1x,A)')'MERCHANTABILITY, or FITNESS FOR A PARTICULAR        '
   write(IOdbg,'(1x,A)')'PURPOSE. You are solely responsible for determining '
   write(IOdbg,'(1x,A)')'the appropriateness of using or redistributing the  '
   write(IOdbg,'(1x,A)')'Work and assume any risks associated with Your      '
   write(IOdbg,'(1x,A)')'exercise of permissions under this License.         '
   write(IOdbg,'(1x,A)')'                                                    '
   write(IOdbg,'(1x,A)')'Limitation of Liability. In no event and under no   '
   write(IOdbg,'(1x,A)')'legal theory, whether in tort (including            '
   write(IOdbg,'(1x,A)')'negligence), contract, or otherwise, unless required'
   write(IOdbg,'(1x,A)')'by applicable law (such as deliberate and grossly   '
   write(IOdbg,'(1x,A)')'negligent acts) or agreed to in writing, shall any  '
   write(IOdbg,'(1x,A)')'Contributor be liable to You for damages, including '
   write(IOdbg,'(1x,A)')'any direct, indirect, special, incidental, or       '
   write(IOdbg,'(1x,A)')'consequential damages of any character arising as a '
   write(IOdbg,'(1x,A)')'result of this License or out of the use or         '
   write(IOdbg,'(1x,A)')'inability to use the Work (including but not limited'
   write(IOdbg,'(1x,A)')'to damages for loss of goodwill, work stoppage,     '
   write(IOdbg,'(1x,A)')'computer failure or malfunction, or any and all     '
   write(IOdbg,'(1x,A)')'other commercial damages or losses), even if such   '
   write(IOdbg,'(1x,A)')'Contributor has been advised of the possibility of  '
   write(IOdbg,'(1x,A)')'such damages.                                       '
   write(IOdbg,'(1x,A)')'----------------------------------------------------'
   write(IOdbg,'(/)')

end subroutine Disclaimer


