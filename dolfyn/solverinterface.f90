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
subroutine SetUpMatrixA(ivar,URFactor,A,Phi,S)

   use constants
   use geometry
   use variables
   use scalars 

   real, dimension(Ncel)      :: A, S 
   real, dimension(Ncel+Nbnd) :: Phi 

   if( Debug > 3 ) write(*,*)'*** SetUpMatrixA  ',Variable(ivar),URFactor

   RURF = 1./URFactor
  
   Res(1:Ncel) = 0.0
   
   ia = 0    
   do i=1,Ncel

     app = A(i)

     do j=1,NFaces(i)
       k  = CFace(i,j)
       ip = Face(k)%cell1
       in = Face(k)%cell2
       if( in > 0 )then
           ! internal
           if( ip == i )then
             aface   = RFace(k,2)
             ia = ia + 1
             Acoo(ia) = DBLE( aface )
             Arow(ia) = i
             Acol(ia) = in
             
             Res(i) = Res(i) - aface*Phi(in)
           elseif( in == i )then
             aface = RFace(k,1)
             ia = ia + 1
             Acoo(ia) = DBLE( aface )
             Arow(ia) = i
             Acol(ia) = ip

             Res(i) = Res(i) - aface*Phi(ip)
           else
             write(*,*)'+ internal error: assembly A-matrix. in=',in
             write(*,*)'+ in SetUpMatrixA for ',Variable(ivar)
           endif
           app = app - aface
       endif
     end do

     A(i)  = app * RURF  
     S(i)  = S(i) + (1.0-URFactor)*A(i)*Phi(i)

     ia = ia + 1
     Acoo(ia) = DBLE( A(i) )
     Arow(ia) = i
     Acol(ia) = i

     Res(i) = Res(i) + S(i) - A(i)*Phi(i)
   end do

   if( ia /= NNZ ) write(*,*)'+ error: SetUpMatrixA: NNZ =',ia,' =/=',NNZ

   Res0 = sum(abs(Res)) * ResiNorm(iVar)
    
   if( Debug > 2 ) write(*,*)'Res0:',Res0,ResiNorm(iVar),ivar,Variable(ivar)

   !do i=1,Ncel     ! row
   !  write(*,'(i2,'': '',$)') i
   !  do j=1,Ncel   ! col
   !    ic = 0
   !    do k=1,NNZ
   !      if( Arow(k) == i )then
   !        if( Acol(k) == j )then
   !          ic = 1
   !          write(*,'(1x,f5.2,$)') SNGL( Acoo(k) )
   !          exit
   !        endif
   !      endif
   !    end do
   !    if( ic == 0 ) write(*,'(''   0  '',$)')
   !  end do
   !  write(*,*) '.',Phi(i),'=',S(i)
   !end do
   
   if( Debug > 3 ) write(*,*)'=== SetUpMatrixA  ' 

end subroutine SetUpMatrixA
subroutine SolveMatrixA(ivar,DoubleVar,res0,res1,A,Phi,S)

   use constants
   use geometry
   use variables
   use scalars 

!   use hypre_mod, only: new_matrix, solve_type

   real, dimension(Ncel)      :: A, S 
   real, dimension(Ncel+Nbnd) :: Phi 
   logical :: DoubleVar

   integer ipar(16)
   real    fpar(16)
   
   if( Debug > 3 ) write(*,*)'*** SolveMatrixA  ',Variable(ivar)
   
   Res0    = sum(abs(Res)) * ResiNorm(iVar)
   ichoose = Solver(iVar)

   !        
   ! solving equation system Ax=B with x=Phi
   ! even makkelijk omgooien later direct in CSR-formaat
   !  
   !           nrow nnz  a    ir   jc   ao  jao  iao
   !
   call coocsr(Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc)
   
   RHS(1:Ncel) = DBLE( S(1:Ncel) )
   SOL(1:Ncel) = DBLE( Phi(1:Ncel) )

   if( ichoose == SparseKit2 )then
     !
     ! SparseKit2 - section
     !

     ipar =  0            ! reset ipar-array
     fpar = 0.0           ! reset fpar-array

     ipar(2) = 0          ! 0 == no preconditioning
     ipar(3) = 0          ! stopping criteria
     ipar(4) = NNZ*8      ! workspace for BGStab = 8n
                          ! gmres = (n+3)*(m+2) + (m+1)*m/2 
     ipar(5) = 0          ! size of the Krylov subspace m, 
                          ! default for gmres => m =15
     ipar(6) = 800        ! maximum number of matrix-vector multiplies

     if( iVar > NVar )then
       fpar(1) = RTOL(VarSC)
       fpar(2) = ATOL(VarSC)
     else
       fpar(1) = RTOL(ivar) ! relative tolerance, must be between (0, 1)
       fpar(2) = ATOL(ivar) ! absolute tolerance, must be positive
     endif

     write(IOdbg,*)'SolveMatrixA: SparsKit2 ',Variable(ivar)

     call solve_spars(debug,IOdbg,Ncel,RHS,SOL,ipar,fpar,&
                            WORK,Acsr,Aclc,Arwc)
     !
     ! and the solution is...
     !
     Phi(1:Ncel) = REAL( SOL )
   
   elseif( ichoose == SolverBiCGstabL   )then
     !
     ! BiCGstabL - section
     !
     write(IOdbg,*)'SolveMatrixA: BiCGstabL' 

     !mxmv  = 100
     !iwork = NNZ*8
     !
     !call solve_bicgstab2(.true.,1,Ncel,SOL,RHS,matvec,.false., &
     !                     1.D-03,'rel',mxmv,WORK,iwork,info)

     !Phi(1:Ncel) = REAL( SOL )
          
   elseif( ichoose == Hypre  )then                                 ! === Hypre ===
     !
     ! Hypre - section
     !
     write(IOdbg,*)'SolveMatrixA: Hypre' 

     ! new_matrix = .true.
     !
     ! call solve_hypre(RHS,SOL)
     !
     ! and the solution is...
     !
     !Phi(1:Ncel) = REAL( SOL )

   elseif( ichoose == SolverDirect  )then
     !
     ! Direct - section
     !
     write(IOdbg,*)'SolveMatrixA: Direct' 

     !
     ! hsl ma27 is only for a symmetric martix A (var. P or PP)
     !
     !call solve_hslma27(debug,IOdbg,Ncel,NNZ,Acoo,Arow,Acol,S,Phi)
     
   elseif( ichoose == SolverUser  )then
     !
     ! User - section
     !
     write(IOdbg,*)'SolveMatrixA: User' 

     !call solve_user(debug,IOdbg,Ncel,NNZ,Acoo,Arow,Acol,S,Phi)
     
   else
     write(*,*)'+ internal error: solver not set'
   endif

   if( debug > 2 ) write(*,*)'SolveMatrixA: residuals' 
   
   Res(1:Ncel) = 0.0

   do i=1,Ncel
     do j=1,NFaces(i)
       k  = CFace(i,j)
       ip = Face(k)%cell1
       in = Face(k)%cell2
       if( in > 0 )then
           ! internal
           if( ip == i )then
             Res(i) = Res(i) - RFace(k,2)*Phi(in)
           elseif( in == i )then
             Res(i) = Res(i) - RFace(k,1)*Phi(ip)
           endif
       endif
     end do
     Res(i) = Res(i) + S(i) - A(i)*Phi(i)
   end do

   Res1 = sum(abs(Res)) * ResiNorm(iVar)

   Residual(iVar) = Res1 

   write(IOdbg,*)'SolveMatrix:',Variable(ivar),':',Res0,'->',Res1
   
   if( Debug > 3 ) write(*,*)'=== SolveMatrixA ' 

end subroutine SolveMatrixA
