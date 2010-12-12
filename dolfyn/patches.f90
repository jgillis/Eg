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
module PatchesModule

   ! constants
   integer, parameter :: IOpat = 100         ! file handle for patch file
   integer, parameter :: ipMaxPatData = 20   ! max. nr of data values per patch
   integer, parameter :: npMaxPatches = 100  ! max. nr of patches in case

   integer, parameter :: ipPatNoSrc=-1, ipPatInit=0, ipPatConVSrc=1, ipPatConMSrc=11
   integer, parameter :: ipPatLinVSrc=2, ipPatLinMSrc=12
   integer, parameter :: ipPatQuadVSrc=3, ipPatQuadMSrc=13
   integer, parameter :: ipPatCatBed=100

   integer, parameter :: ipPatTest=99
   !
   ! type definitions
   !
   type PatchData
      integer           iId                ! ID of patch (used to report sources etc.)
      character(len=15) sName              ! Name of patch (used to report sources etc.)
      integer           iType              ! Type of patch (constant, linear, special...)
      character(len=4)  sVar               ! Variable short name or 'none' or 'all'
      integer           iVar               ! Variable number for this patch (for quick ref) or 0 for none or -1 for all
      integer           iCellID            ! Restrict working to a cell type (or 0 for all types)
      real              xyz1(3)            ! low point of box
      real              xyz2(3)            ! high point of box
      integer           iDatlen            ! length of additional data (real numbers) that follow
      real              dat(ipMaxPatData)  ! additional data
   end type

   !
   ! variables
   !
   type(PatchData) :: Patch(npMaxPatches)    
   integer         :: nPatches             ! count number of patches

end module PatchesModule
subroutine PatchesSetUp(iStage)
   !
   ! general initialization for user routines
   !
   use constants
   use geometry
   use variables
   
   use PatchesModule

   !
   ! arguments
   !
   integer  iStage

   select case (iStage)
     case (0)
       !
       ! before anything else
       !
       write(*,*)"NOTE: this is DOLFYN with patches by B. Tuinstra"
       write(*,*)"DATE of this version: Febr. xx, 2006"
       write(*,*)"see http://www.home.zonnet.nl/bouke_1/dolfyn"
       write(*,*)"Use at your own risk"
     case (1)

       if( UseScalars )then
         write(*,*)'Setting up scalars for patches...'
         call PatchesScalarsSetUp(1)
       endif
       !
       ! standard input files have been read
       !
       write(*,*)'Reading patch file...'
       call ReadPatchFile(casename)

       !
       ! initialize scalar arrays 
       !
       call PatchesInitializeScalars
       
       !
       ! Get variable numbers from names in patches
       !
       do ip=1,NPatches
          Patch(ip)%iVar = IFindVar(Patch(ip)%sVar)
       end do

     case (2)
        !
        ! after allocation of variables
        !
        call PatchesInit

        write(*,*)     'Patches initialised'
        write(IOdbg,*) 'Patches initialised'
   end select

end subroutine PatchesSetUp
subroutine PatchesExit(iStage)
   !
   ! chance to clean up user routines
   !
   use constants
   use geometry
   use variables
   
   use PatchesModule

   !
   ! arguments
   !
   integer iStage

   select case (iStage)
     case (0)
       !
       ! At the very end, after results are written but before main prog cleanup
       !
       ! Clean up scalar regions (scalar variables are done in Dolfyn.f90)
       !
!       call DeAllocateScalars()

     case (1)
       !
       ! Cleanup before results are written
       !
       ! Write global balance to debug file
       ! 
       call WriteBalance(ioDbg)

   end select
   
end subroutine PatchesExit
subroutine PatchesInit
   !
   ! Chance to initialize field, patch-based
   !
   
   !
   ! modules
   !
   use constants
   use geometry
   use variables

   use scalars
   
   use PatchesModule

   !
   ! external function types
   !
   logical LCelInPat

   !
   ! locals
   !
   integer ip, ic, it, is

   do ip = 1, nPatches
     it = Patch(ip)%iType
     if( it == 0 )then
       do ic=1,Ncel
         ! if cell within patch box
         if( LCelInPat(ip,ic) )then
         
           select case (Patch(ip)%iVar)
             case(-1)       ! initialize 'all' variables
                    ! patch data format?
             case(varU)
                U(ic)  = Patch(ip)%dat(1)
             case(varV)
                V(ic)  = Patch(ip)%dat(1)
             case(varW)
                W(ic)  = Patch(ip)%dat(1)
             case(varP)
                P(ic)  = Patch(ip)%dat(1)
             case(varT)
                T(ic)  = Patch(ip)%dat(1)
             case(varTE)
                TE(ic) = Patch(ip)%dat(1)
             case(varED)
                ED(ic) = Patch(ip)%dat(1)                
           end select
           
           if( Patch(ip)%iVar > nVar ) then
             
!            is = Patch(ip)%iVar-npStdVars         !<=============!!!
             is = Patch(ip)%iVar-NVar          
             
             SC(ic, is) = Patch(ip)%dat(1)
             
           endif
         endif
       end do
     endif
   end do
end subroutine PatchesInit
subroutine PatchesSource(iVar, Phi, A, S)
   !
   ! calculate sources for all cells
   !
   use constants
   use geometry
   use variables
   
   use PatchesModule

   !
   ! arguments
   !
   integer,  intent(in)                   :: iVar   ! variable number
   real, dimension(Ncel+Nbnd), intent(in) :: Phi    ! variable field
   real, dimension(Ncel), intent(inout)   :: A      ! coefficient
   real, dimension(Ncel),intent(inout)    :: S      ! source [[Phi] kg/s]

   ! The source term S (as all other terms) 
   ! should have the unit of
   ! V d{rho phi}/dt, or [meter^3] [mass/meter^3] [phi]/[time] 
   !                   = [phi][mass]/[time]

   !
   ! external function types
   !
   logical LCelInPat

   ! locals
   integer ip, ic


   if( Debug > 3 ) write(*,*)'*** PatchesSource: ',Variable(iVar),nPatches 
   
   !
   ! loop over all patches
   !
   do ip=1,nPatches
     !
     ! if iVar is patch variable or patch variable =-1 (=> 'all')
     !
     if( iVar == Patch(ip)%iVar .OR. Patch(ip)%iVar == -1 )then
       !write(*,*)     'Apply PatchesSource to ',Variable(iVar)
       write(IOdbg,*) 'Apply PatchesSource to ',Variable(iVar)
       select case (Patch(ip)%iType)
         case(ipPatNoSrc)       
           !
           ! no source
           !
         case(ipPatConVSrc)     
           !
           ! constant source per volume: 
           ! source = dat(1) [[Phi] kg/s/m3]
           !
           do ic=1,Ncel
             !
             ! if cell within patch box
             !
             if( LCelInPat(ip,ic) )then
               !
               ! Set source for this patch & cell
               ! Below is correct for heat source term patch%dat(1)=(Q/Cp) 
               ! where Q in [W/m3] so Q/Cp in [kg K/m3/s]
               !
               SS = Patch(ip)%dat(1) * Cell(ic)%vol
               S(ic) = S(ic) + SS
             endif
           end do
           
         case(ipPatLinVSrc)     
           !
           ! linear source per volume:
           ! source = dat(1) - dat(2)*Phi [[Phi] kg/s/m3]
           !
           do ic=1,Ncel
             if( LCelInPat(ip, ic) )then
               !
               ! source = Sref - A*phi; dat(1) = Sref, dat(2) = A
               ! for heat source: Sref = Q[W/m3]/CP in [kg K/s/m3] 
               ! and A = (dQ/dT)/Cp in [kg/s/m3]
               !
               SS = Patch(ip)%dat(1) * Cell(ic)%vol
               SA = Patch(ip)%dat(2) * Cell(ic)%vol
               S(ic) = S(ic) + SS
               A(ic) = A(ic) + SA
             endif
           end do
           
         case(ipPatQuadVSrc)    
           ! 
           ! quadratic source per volume 
           ! (note the absolute value to prevent positive coef):
           ! source = dat(1) - dat(2)*Phi - dat(3)*Phi*abs(Phi) [[Phi] kg/s/m3]
           !
           !write(*,*)'pat3>',Patch(ip)%dat(1),Patch(ip)%dat(2),Patch(ip)%dat(3)

           do ic=1,nCel
             if( LCelInPat(ip,ic) )then
               SS = Patch(ip)%dat(1) * Cell(ic)%vol
               SA = (Patch(ip)%dat(2)+Patch(ip)%dat(3)*abs(Phi(ic))) * Cell(ic)%vol
               S(ic) = S(ic) + SS
               A(ic) = A(ic) + SA
             endif
           end do
           
         case(ipPatTest)    
           ! 
           ! quadratic source per volume 
           ! (note the absolute value to prevent positive coef):
           ! source = dat(1) - dat(2)*Phi - dat(3)*Phi*abs(Phi) [[Phi] kg/s/m3]
           !
           !write(*,*)'pat4>',Patch(ip)%dat(1),Patch(ip)%dat(2),Patch(ip)%dat(3)

           icnt = 0
           valmin =  1.e30
           valmax = -1.e30
           
           do ic=1,nCel
             if( LCelInPat(ip,ic) )then
               
               icnt = icnt + 1
               
               quitelarge = 1000000.0
               
               vmag = sqrt( U(ic)**2 + V(ic)**2 + W(ic)**2 )
               func = Patch(ip)%dat(1) +  Patch(ip)%dat(2) * vmag & 
                    + Patch(ip)%dat(3) * vmag**2
                              
               valmin = min( valmin, func )
               valmax = max( valmax, func )
               !write(*,*) '>>',ic,vmag,func
                              
               SS =   quitelarge * func  
               SA =   quitelarge 
               
               S(ic) = S(ic) + SS
               A(ic) = A(ic) + SA
             endif
           end do

           write(IOdbg,*) 'pat typ 99 ',Variable(iVar),':',icnt,valmin,valmax
           
         case(ipPatConMSrc)     
           !
           ! constant source per kg fluid: 
           ! source = dat(1) [[Phi] m3/s/m3]
           !
           do ic=1,Ncel
             if( LCelInPat(ip, ic) )then
               SS = Patch(ip)%dat(1) * Cell(ic)%vol * DEN(ic)
               S(ic) = S(ic) + SS
             endif
           end do
           
         case( ipPatLinMSrc )     
           !
           ! linear source per kg fluid:
           ! source = dat(1) - dat(2)*Phi [[Phi] kg/s/m3]
           !
           do ic=1,nCel
             if( LCelInPat(ip,ic) )then
               !
               ! see type 2 for details; this differs factor DEN(ic)
               !
               SS = Patch(ip)%dat(1) * Cell(ic)%vol * DEN(ic)
               SA = Patch(ip)%dat(2) * Cell(ic)%vol * DEN(ic)
               S(ic) = S(ic) + SS
               A(ic) = A(ic) + SA
             endif
           end do
           
         case( ipPatQuadMSrc )    
           !
           ! quadratic source per kg fluid 
           ! (note the absolute value to prevent positive coef):
           ! source = dat(1) - dat(2)*Phi - dat(3)*Phi*abs(Phi) [[Phi] kg/s/m3]
           !
           do ic=1,Ncel
             if( LCelInPat(ip,ic) )then
               SS = Patch(ip)%dat(1)*Cell(ic)%vol*Den(ic)
               SA = (Patch(ip)%dat(2) + Patch(ip)%dat(3)*abs(Phi(ic)))*Cell(ic)%vol*Den(ic)
               S(ic) = S(ic) + SS
               A(ic) = A(ic) + SA
             endif
           end do
           
       end select
     end if
   end do

   if( Debug > 3 ) write(*,*)'=== PatchesSource  ' 
   
end subroutine PatchesSource
subroutine PatchesDispersion()

   use constants
   use geometry
   use variables
!   use scalars
   
   use PatchesModule
   !
   ! external function types
   !
   logical LCelInPat

   !
   ! locals
   !
   integer iPat, iFac, iReg, iType
   logical lC1, lC2

!   do iFac=1,Nfac
!      do iPat=1,nPatches
!      
!      end do
!   end do

end subroutine PatchesDispersion
subroutine PatchesBnd(iaF, iaVar, raVal, raCo)

   use constants
   use geometry
   use variables
   
   use PatchesModule
   !
   ! arguments
   !
   integer, intent(IN) :: iaF    ! Face # fow which boundary data is needed
   integer, intent(IN) :: iaVar  ! Variable # for which boundary data is needed
   real, intent(OUT)   :: raVal  ! Face value of the variable
   real, intent(OUT)   :: raCO   ! Coefficient
   !
   ! locals
   !
   integer  ib, ir
   integer  iscal                ! scalar number

   ib = Face(iaF)%bnd
   ir = Bnd(ib)%rid

   write(*,*) 'Error: UserBND called but not implemented yet'
   raVal = 0.0
   raCo  = 0.0
    
end subroutine PatchesBnd
logical function LCelInPat(ip, ic)

   use constants
   use geometry
   use variables
   
   use PatchesModule
   !
   ! arguments
   !
   integer, intent(IN)             ::      ip
   integer, intent(IN)             ::      ic
   !
   ! locals
   !
   logical linpat

   !
   ! Restrict action to given fluid/cell type id BUT 
   ! the cell has to be IN the defined box!
   ! (or no restriction if iCellID=0)
   !
   
   linpat = (Patch(ip)%iCellID == 0 .or. Patch(ip)%iCellID == Cell(ic)%ctid )

   if( .NOT.linpat ) goto 10

   do i=1,3

     if( Cell(ic)%x(i) < Patch(ip)%xyz1(i) )then
       linpat = .FALSE.
       exit
     end if
     if( Cell(ic)%x(i) > Patch(ip)%xyz2(i) )then
       linpat = .FALSE.
       exit
     end if
   end do

10 continue

   LCelInPat = linpat 

end function LCelInPat
subroutine ReadPatchFile(sacasename)
   !
   ! read PATch file
   !
   use constants
   use variables
   use geometry
!   use scalars
   
   use PatchesModule

   !
   ! arguments
   !
   character*(*)       :: sacasename

   !
   ! external function types
   !
   logical             :: LFReadScalarData
   logical             :: LFReadMediumData

   !
   ! locals
   !
   character(len=120)  :: string
   character(len=32)   :: s

   integer, parameter  :: MaxKeys = 20
   character(len=16)   :: keys(MaxKeys) 
   integer             :: keyi(MaxKeys)
   real                :: keyr(MaxKeys)
   logical             :: keyu(MaxKeys)
   integer             :: n
   logical             :: llExists

   nPatches = 0
   
   !
   ! look if patch file exists
   !
   n = lens(saCaseName)
   string = saCaseName(1:n)//'.pat'
   n = lens(string)
    
   inquire( file=string(1:n), exist=llExists )

   ! if no patch file, do not try to read it
   if( .not. llExists )return

   !
   ! patch file exitst, open & read it.
   !
   call openfile(IOpat,sacasename,'.pat','FORMATTED', &
                 'SEQUENTIAL','OLD',debug)

   do while (iReadOne(IOpat, string)>=0)
     
     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

     write(ioDbg,*) 'Patch file patch data read: ', string

     if( Nkeys == 0 ) cycle           ! skip blank line

     if( keys(1) == 'patch' )then
       !
       ! note: within a patch, the fields must be in a fixed order
       !
       if( nPatches == npMaxPatches )then
          write(*,*) 'Error: maximum number of patches read'
          goto 90
       end if
       n = nPatches + 1
       Patch(n)%iId=keyi(2)

       if(iReadOne(IOpat, string) < 0 ) goto 90
       call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
       if (keys(1).ne.'name') &
          write(*,*) 'Error reading patch (',Patch(n)%iId, &
                     '): name expected instead of ', keys(1)
       Patch(n)%sName= keys(2)

       if(iReadOne(IOpat, string) < 0 ) goto 90
       call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
       if (keys(1).ne.'type') &
          write(*,*) 'Error reading patch (',Patch(n)%iId, &
                     '): type expected instead of ', keys(1)
       Patch(n)%iType= keyi(2)

       if (iReadOne(IOpat, string)<0) goto 90
       call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
       if (keys(1).ne.'variable') &
          write(*,*) 'Error reading patch (',Patch(n)%iId, &
                     '): variable expected instead of ', keys(1)
       Patch(n)%sVar = keys(2)

       if (iReadOne(IOpat, string)<0) goto 90
       call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
       if( keys(1).ne.'cellid') &
          write(*,*) 'Error reading patch (',Patch(n)%iId, &
                     '): cellid expected instead of ', keys(1)
       Patch(n)%iCellID = keyi(2)

       if (iReadOne(IOpat, string)<0) goto 90
       call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
       if (keys(1).ne.'box') & 
          write(*,*) 'Error reading patch (',Patch(n)%iId, &
                     '): box expected instead of ', keys(1)
       Patch(n)%xyz1(1:3) = keyr(2:4)
       Patch(n)%xyz2(1:3) = keyr(5:7)

       if (iReadOne(IOpat, string)<0) goto 90
       call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
       if (keys(1).ne.'datlen') &
          write(*,*) 'Error reading patch (',Patch(n)%iId, &
                     '): datlen expected instead of ', keys(1)
       Patch(n)%iDatLen= keyi(2)
       
       if( Patch(n)%iDatLen > ipMaxPatData )then
          Patch(n)%iDatLen = ipMaxPatData
          write(*,*) 'Error reading patch (',Patch(n)%iId, &
                     '): number of data items(', keyi(1), &
                     ') exceeds maximum allowed (',ipMaxPatData,')'
       end if
 
       do i=1,Patch(n)%iDatLen
          if (iReadOne(IOpat, string)<0) goto 90
          call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
          Patch(n)%dat(i)= keyr(2)
       end do

       nPatches = n
     
     else if (LFReadScalarData(ioPat, keys(1), string)) then
     
       write(ioDbg,*) 'Patch file scalar data read: ', string
     
     !else if (LFReadMediumData(ioPat, keys(1), string)) then
     
     !   write(ioDbg,*) 'Medium data read: ', string
     
     end if
   end do
90 continue
   
   close(IOPat)
   
end subroutine ReadPatchFile
integer function iReadOne(IO, string)
   !
   ! arguments
   !
   integer         IO
   character(*)    string
   !
   ! locals
   !
   integer    iRet
   integer    indx

   iRet = -1

   !
   ! Read line, skip comment lines (starting with '#')
   !
   string(1:1)='#'
   do while (string(1:1)=='#')
     
     read(IO,'(A)',end=90) string                    
   
   end do

   !
   ! Strip trailing comments
   !
   indx = index(string,'#')
   if( indx >= 2 )then
     do i=indx,len(string)
       string(i:i) = ' '
     end do
   endif
   
   iRet = lens(string)
   if (string(1:1)=='#') iRet=-1

90 continue
   
   iReadOne = iRet
   
end function iReadOne

!
! hoort hier niet!
!
subroutine WriteBalance(io)

   use Constants
   use Geometry
   use Variables
   use Scalars

   !
   ! arguments
   !
   integer, intent(IN) :: io     ! file handle

   !
   ! locals
   !
   integer j, it, iVar, iReg
   real    rin, rout

   do i=1,MaxSC
      if( SolveScalar(i) )then
         rin  = 0
         rout = 0
         write(io,*)
         write(io,*) 'Scalar: ', Scalar(i)
         write(io, '(A4,A8, 2A15)')'Bnd#','Type', 'Flux in', 'Flux out'
         
         do iReg=1,Nreg
            it = Reg(iReg)%typ
  !          write(io, '(I4,A8, 2E15.7)') &
  !             iReg, Region(it), Reg(iReg)%ScFluxIn(i), Reg(iReg)%ScFluxOut(i)
  !          rin  = rin + Reg(iReg)%ScFluxIn(i)
  !          rout = rout + Reg(iReg)%ScFluxOut(i)
         end do
         write(io,'(A4, A8, 2A15)')'+ --','------','---------------', '---------------'
         write(io,'(A4, A8, 3E15.7)')'Tot:','      ', rIn, rOut, rIn-rOut
      endif
   end do
end subroutine WriteBalance
