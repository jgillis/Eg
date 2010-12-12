!
! Copyright 2006-2007 Bouke Tuinstra
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
! The subroutines in this file originally written by B. Tuinstra, 
! see www.home.zonnet.nl/bouke_1/dolfyn
!
! Modified/pulled into the main code, febr. 2006
!
subroutine PatchesScalarsSetUp(iStage)

   use constants
   use geometry    ! need NReg...
   use scalars 

   !
   ! arguments
   !
   integer  iStage

   select case (iStage)
     case (1)
       !
       ! allocate memory for scalar region extension
       !
       write(*,*)'Allocating ',NReg*MaxSC,' for ScReg'
       allocate(ScReg(NReg,MaxSC), stat=istat)
       call TrackMemory(istat, NReg*MaxSC, 'Scalar Region array allocated')

       write(*,*)'Allocating ',MaxSC,' for ScProp'
       allocate(ScProp(MaxSC), stat=istat)
       call TrackMemory(istat, MaxSC, 'Scalar Property array allocated')
       !
       ! initialize default scalar boundary values (symmetry)
       !
       do i=1,NReg
         do j=1,MaxSC
           ScReg(i,j)%value = 0.0
         end do
       end do
       
      case (2)
         !
         ! after allocation of variables
         !
!         call ScalarsInit
   end select

end subroutine PatchesScalarsSetUp
subroutine DeAllocateScalars

   use constants, only: MaxSC
   use geometry, only: Nreg
   use scalars

   deallocate(ScReg, stat=istat)
   call TrackMemory(istat,-NReg*MaxSC, 'Scalar Region array deallocated')

end subroutine DeAllocateScalars
subroutine PatchesInitializeScalars

   use constants
   use geometry
   use variables
   use scalars

   do i=1,NScal
     SolveScalar(i) = ScProp(i)%solve
     Solve(NVar+i)  = SolveScalar(i)
   end do

   do i=1,NScal
     StoreScalar(i) = .true.
     Solve(NVar+i)  = .true.
   end do

   do i=1,NScal
     Scalar(i)        = ScProp(i)%name
     Variable(NVar+i) = ScProp(i)%name
     Gamma(NVar+i)    = Gamma(VarSC)          ! all scalars the same scheme
   end do
   
end subroutine PatchesInitializeScalars
logical function LFReadScalarData(iaIO, saKey, saString)
   !
   ! read data from input or patch file
   ! uses parser from readcontrolfile.f90
   !
   use constants
   use variables, only: Nscal, Scalar
   use geometry    ! need Reg()
   use scalars

   ! arguments
   integer, intent(in)       :: iaIO      ! file handle
   character*(*), intent(in) :: saKey     ! first key
   character*(*), intent(in) :: saString  ! first line

   ! locals
   character(len=120)  :: string
   character(len=32)   :: s
   integer, parameter  :: MaxKeys = 20
   character(len=16)   :: keys(MaxKeys) 
   integer             :: keyi(MaxKeys)
   real                :: keyr(MaxKeys)
   logical             :: keyu(MaxKeys)
   integer             :: i, n

   LFReadScalarData = .false.

   if( saKey == 'scalar' )then
     ! SCALAR block:
     ! scalar, [name], [unit]       ! unit is optional?
     ! solved, [on|off]             ! if not solved: constant field (property, ...)
     ! stored, [on|off]             ! written to output file?
     ! prlam, [value]               ! laminar diffusion
     ! prtur, [value]               ! turbulent diffusion
     
     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,saString)

     n = NScal + 1

     if( n > MaxSC ) then
       write(*,*) 'Error: maximum number of scalars read: (', &
                  keys(2),') skipped.'
       return
     end if
     
     ScProp(n)%name = keys(2)
     ScProp(n)%unit = keys(3)

     if( iReadOne(iaIO, string) < 0 ) return
     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

     if( keys(1) .ne. 'name' ) &
       write(*,*) 'Error reading scalar (',ScProp(n)%name, &
                  '): name expected instead of ',keys(1)

     if( iReadOne(iaIO, string) < 0 ) return
     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
     if( keys(1) .ne. 'solved' ) &
       write(*,*) 'Error reading scalar (',ScProp(n)%name, &
       '): solved expected instead of ', keys(1)
       
     ScProp(n)%solve = ( keys(2) == 'on' .or. keys(2) == 'yes')

     if( iReadOne(iaIO, string) < 0 ) return
     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

     if( keys(1)(1:6) .ne. 'stored' ) &                           ! stond SAVED!
       write(*,*) 'Error reading scalar (',ScProp(n)%name, &
                  '): ''stored'' expected instead of ', keys(1)
                  
     ScProp(n)%save = ( keys(2) == 'on' .or. keys(2) == 'yes' )

     if( iReadOne(iaIO, string) < 0 ) return
     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
     if( keys(1) .ne. 'prlam' ) &
       write(*,*) 'Error reading scalar (',ScProp(n)%name, &
                  '): prlam expected instead of ', keys(1)
     
     ScProp(n)%PrL= keyr(2)

     if( iReadOne(iaIO, string) < 0 ) return
     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
     if( keys(1) .ne. 'prtur' ) &
       write(*,*) 'Error reading scalar (',ScProp(n)%name, &
                  '): prtur expected instead of ', keys(1)
                  
     ScProp(n)%PrT= keyr(2)

     NScal = n
     Scalar(n) = ScProp(n)%name
     
     !write(*,*)'>>>>',Scalar(n)

     LFReadScalarData = .true.

   elseif( saKey == 'scalarbounds' )then
     ! scalarbounds, iID
     ! number, nS
     ! => nS times: scalar, sNames(i), ['flux'|'value'], rValue(i)
     !
     ! type ScalarBounds
     !   integer          iId                  ! number of corresponding boundary
     !   integer          nS                   ! number of scalars defined
     !   character(len=4) sNames(npMaxScalars) ! names of defined scalars
     !   real             rValue(npMaxScalars) ! value of defined scalar
     !   logical          lFlux(npMaxScalars)

     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,saString)
     if( keyi(2) < 0 .or. keyi(2) > NReg ) &
       write(*,*) 'Error reading scalarbounds: region number ', &
                  keyi(2),' is invalid'
     ir = keyi(2)

     if( iReadOne(iaIO, string) < 0 ) return
     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
     if( keys(1) .ne. 'number' ) &
       write(*,*) 'Error reading scalarbounds (',ir, &
                  '): ''number'' expected instead of ', keys(1)
     if( keyi(2) .lt. 0 .or. keyi(2) > MaxSC ) &
       write(*,*) 'Error reading scalarbounds (',ir, &
                  '): number of scalar conditions ',keys(2),' is invalid'
     
     ns = keyi(2)

     ! read in scalar names & values.
     do i=1,ns
       if( iReadOne(iaIO, string) < 0 ) return
       call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
       if( keys(1) .ne. 'scalar' ) &
         write(*,*) 'Error reading scalarbounds (',ir, &
                    '): ''scalar'' expected instead of ', &
                    keys(1),' for scalar line ', i
       
       j = IFindVar(keys(2)) 

       if( j <= 0 ) write(*,*) 'Error reading scalarbounds (', &
                               ir,'): unknown variable name ',keys(2),&
                               ' for scalar line ',i

       if( keys(3) == 'flux' )then
         ScReg(ir,j)%flux = .true.
       elseif( keys(3) == 'value' )then
         ScReg(ir,j)%flux = .false.
       else 
         write(*,*) 'Error reading scalarbounds (',ir, &
                    '): ''flux'' or ''value'' expected instead of ', &
                    keys(3),' for scalar line ', i
       endif
       
       if( ScReg(ir,j)%flux .and. &
         ( Reg(ir)%typ == RInlet .or. Reg(ir)%typ == ROutlet) )then
         write(*,*) 'Error reading scalarbounds (',ir, &
         '): flux condition is not allowed for inlets/outlets; value assumed'
       endif
       
       ScReg(ir,j)%value = keyr(4)
       
     end do

     LFReadScalarData = .true.

   endif

end function LFReadScalarData
subroutine WriteScalarResiduals(io)

   use constants
   use variables
   use geometry
   use scalars

   integer, intent(in) :: io        ! file handle
   character(len=15)   :: string
   integer             :: i

   write(string,'(a4,i2,a9)') '(a6,',NScal,'(1pe9.2))'
!   write(io,string)'  sc->',(SCResidual(i),i=1,NScal)

end subroutine WriteScalarResiduals
function IFindVar(target)
   !
   ! find variable index of variable named target
   !
   use constants
   use variables
   use geometry 

   character(len=12), intent(in) :: target
   character(len=4)  :: key, name
   integer           :: is

   key = target(1:4)
   call lowercase(key)

   is = 0
   if( key == 'none' )then
     is = 0
   elseif( key(1:3) == 'all' )then
     is = -1
   else
     do i=1,nVar
       name = Variable(i)(1:4)
       call lowercase(name)
       if( key == name )then
         is = i
         goto 10
       endif
     end do

     do i=1,nScal
       name = Scalar(i)(1:4)
       call lowercase(name)
       if( key == name )then
         is = i
         goto 10
       endif
     end do
10   continue
   endif
  
   IFindVar = is
   
end function IFindVar
subroutine Set_Scalar_Norm()
   !
   ! Calculate residual normalization factor for scalars
   ! this code could be incorporated into routine Set_Normalisation_Factors
   ! but I keep it here to minimize code changes in Dolfyn.f90 
   ! (it is called from there)
   !
   use constants
   use geometry
   use variables
   use scalars

   if( Debug > 3 ) write(*,*)'*** Set_Scalar_Norm'
   !
   ! loop over all scalars
   !
   do is=1,NScal
     ResiNorm(NVar+is) = 1.0
     sum = 0.0
     
     if( ScProp(is)%solve )then
       !
       ! loop over all boundaries
       !
       do ib=1,Nbnd
         i  = Bnd(ib)%face
         ir = Bnd(ib)%rid
         select case( Reg(ir)%typ )
           case( RInlet )
             !
             ! Inlet condition: source is mass flux * value
             !
             fluxin = Reg(ir)%den * dot_product( Reg(ir)%uvw , Face(i)%n )
             sum = sum + abs( fluxin * ScReg(ir,is)%value )
           case( ROutlet )
             !
             ! Outlet condition: no contribution?
             !
           case( RSymp )
             !
             ! Symmetry condition: no contribution
             !
           case( RWall )
             if( ScReg(ir,is)%flux )then
               !
               ! for now, only flux wall conditions are taken into account...
               !
               sum = sum + abs( ScReg(ir,is)%value )*Face(i)%area
             endif
         end select
       end do
     end if 
     
     !
     ! scalar solved
     !
     if( sum > 0.0 ) ResiNorm(NVar+is) = 1.0/sum

     write(IOdbg,*)'ResiNorm for scalar ',is,':',ResiNorm(NVar+is)
     
   end do

   if( Debug > 3 ) write(*,*)'=== Set_Scalar_Norm'

end subroutine Set_Scalar_Norm
