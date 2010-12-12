!
! Copyright 2006-2007 Shibo (Harry) Kuang
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
subroutine dolfyn2tecplt(Model)

   use constants
   use geometry
   use variables
   use particles

   
   character(len=40) :: string1, string2
   integer :: indx(8)

   integer            :: datum(8)
   integer            :: Iuse               ! Number of ploted variables

   character (len=12) :: clock(3)
   character (len=10) :: datestring 

   character (len=96) :: PltVar=' '
   character (len=10) :: CharTemp=''

   integer, intent(IN)  :: Model            ! Model=1 write particle data

   real, dimension(3)   :: Xp, X, ds 
   real, dimension(3)   :: Xpn, Normal 
   real, dimension(10)  :: tmp              ! temp array 
   real, dimension(12)  :: Ptmp=0.0
   
   real, allocatable    :: NodalData1(:)    ! array of data at vertices
   real, allocatable    :: NodalData2(:)    
   real, allocatable    :: NodalData3(:)    
   
   integer, allocatable :: NodalCounter(:)  

   integer, allocatable :: PLTcells(:)  

   integer, save        :: ICounter = 0     ! counter for transient data
   real, save           :: TimePrevious=0

   logical :: Exists

   call date_and_time(clock(1),clock(2),clock(3),datum)
   write(datestring,'(i2.2,''/'',i2.2,''/'',i4)') datum(3),datum(2),datum(1)

   allocate( PLTcells(Ncel+Nbnd),stat=istat)    
   allocate( NodalData1(Nvrt),stat=istat)    
   allocate( NodalCounter(Nvrt),stat=istat)   
   call TrackMemory(istat,Ncel+Nbnd+2*Nvrt,'Work arrays for NodalData allocated')
   
   write(*,*)        'Writing TECPLT data file...'
   write(IOdbg,*)    'Writing TECPLT data file...'
   write(IOdbg,'(A,A)')  ' Case:  ',casename(1:lens(casename)) 
   write(IOdbg,'(A,A)')  ' Title: ',title(1:lens(title)) 
   write(IOdbg,'(A,A)')  ' Date : ',datestring

   ICounter=ICounter+1
   Exists      = .false.

   if(ICounter==1)then
     inquire(file=casename(1:lens(casename))//'.dat',exist=Exists)
     if( Exists )then
       call openfile(IOpst,casename,'.dat','FORMATTED', &
                           'SEQUENTIAL','OLD',debug)
       write(*,'(2x,A)')'Remove existing file'
       close(IOpst,status='DELETE') 
     endif
   endif
   
   if( Transient )then
     write(string1,'(f10.4)') time
     if(time.le.TimePrevious)return
       TimePrevious=time
       string1='"Time='//string1(1:lens(string1))//'Sec"'
       call openfile(IOpst,casename,'.dat','FORMATTED', &
                                         'APPEND','UNKNOWN',debug)
     else
       string1='"Steady"'
       call openfile(IOpst,casename,'.dat','FORMATTED', &
                                         'SEQUENTIAL','UNKNOWN',debug)
   endif
   
! Postv(VarU)=.True.
! PostV(VarP)=.True.
! PostV(VarTE)=.false.
! PostV(VarVIS)=.false.

   !
   !  Note, Order of Variables corresponds to their data
   !
   PltVar='VARIABLES = "X", "Y", "Z", "Vmag" '
   Iuse=1
   
   if( PostV(VarU) .or. PostV(VarV) .or. PostV(VarW) )then
     PltVar=PltVar(1:lens(PltVar))//',"U","V","W"'
     Iuse=Iuse+3
   endif
   
   if( PostV(VarP) )then
     PltVar=PltVar(1:lens(PltVar))//',"P"'
     Iuse=Iuse+1
   endif
   
   if( PostV(VarT) )then
     PltVar=PltVar(1:lens(PltVar))//',"T"'
     Iuse=Iuse+1
   endif

   if( PostV(VarTE) .or.PostV(VarED) )then
     PltVar=PltVar(1:lens(PltVar))//',"TE","ED"'
     Iuse=Iuse+2
   endif
   
   if( PostV(VarVIS) )then
     PltVar=PltVar(1:lens(PltVar))//',"VIS"'
     Iuse=Iuse+1
   endif
   
   if( SolveScalars )then
     if( PostV(VarSC) )then
       do is=1,NScal
         write(Chartemp,'(I2)')is
         PltVar=PltVar(1:lens(PltVar))//',SCA'//Chartemp
         Iuse=Iuse+1
       end do
     endif
   endif

   write(IOpst,'(A)')'TITLE='//casename(1:lens(casename)) 
   write(IOpst,'(A)')PltVar(1:lens(PltVar)) 
   write(IOpst,'(A,1x,I8,1x,A,1x,I8,1x,A)')  &
         'ZONE T='//string1(1:lens(string1))//' N=',NVrt, &
         ', E=',Ncel,', F=FEBLOCK, ET=BRICK'
   !
   ! vertices
   !
   do i=1,Nvrt,10
     k = min(i+9,Nvrt)
     write(IOpst,'(10(1x,1pe12.5))')(Vert(j,1),j=i,k) 
   end do
   
   do i=1,Nvrt,10
     k = min(i+9,Nvrt)
     write(IOpst,'(10(1x,1pe12.5))')(Vert(j,2),j=i,k) 
   end do
   
   do i=1,Nvrt,10
     k = min(i+9,Nvrt)
     write(IOpst,'(10(1x,1pe12.5))')(Vert(j,3),j=i,k) 
   end do

   allocate( NodalData2(Nvrt),stat=istat)    
   allocate( NodalData3(Nvrt),stat=istat)    

   call TrackMemory(istat,2*Nvrt,'Work arrays for extra NodalData allocated')


   call InterpolateData(0,VarU,U,dUdX,NodalData1,NodalCounter)
   call InterpolateData(0,VarV,V,dVdX,NodalData2,NodalCounter)
   call InterpolateData(0,VarW,W,dWdX,NodalData3,NodalCounter)

   write(IOdbg,*)'Writing TECPLT Vmag on nodes ',Nvrt
   do i=1,Nvrt,10
     k = min(i+9,Nvrt)
     write(IOpst,'(10(1x,1pe10.3))') &
       (sqrt( Nodaldata1(j)**2 + Nodaldata2(j)**2 + Nodaldata3(j)**2 ),j=i,k)
   end do


   if( PostV(VarU) .or. PostV(VarV) .or. PostV(VarW) ) then
      write(IOdbg,*)'Writing TECPLT vectors on nodes ',Nvrt

      do i=1,Nvrt,10
        k = min(i+9,Nvrt)
        if( abs(Nodaldata1(i)) < Small ) Nodaldata1(i) = 0.0
        write(IOpst,'(10(1x,1pe10.3))') (Nodaldata1(j),j=i,k)
      end do

      do i=1,Nvrt,10
        k = min(i+9,Nvrt)
        if( abs(Nodaldata2(i)) < Small ) Nodaldata2(i) = 0.0
        write(IOpst,'(10(1x,1pe10.3))') (Nodaldata2(j),j=i,k)
      end do

      do i=1,Nvrt,10
        k = min(i+9,Nvrt)
        if( abs(Nodaldata3(i)) < Small ) Nodaldata3(i) = 0.0
        write(IOpst,'(10(1x,1pe10.3))') (Nodaldata3(j),j=i,k)
      end do

    endif

    if( PostV(VarP) )then
      write(IOdbg,*)'Writing tecplt pressure on vertices'

      call InterpolateData(0,VarP,P,dPdX,NodalData1,NodalCounter)

      do i=1,Nvrt,10
        k = min(i+9,Nvrt)
        write(IOpst,'(10(1x,1pe10.3))') (Nodaldata1(j),j=i,k)
      end do
    endif

    if( PostV(VarT) )then
      write(IOdbg,*)'Writing tecplt temperature on vertices'

      call InterpolateData(0,VarT,T,dPdX,NodalData1,NodalCounter)

      do i=1,Nvrt,10
        k = min(i+9,Nvrt)
        write(IOpst,'(10(1x,1pe11.4))')  (Nodaldata1(j)-Tref,j=i,k)
      end do
    endif

    if(PostV(VarTE).or.PostV(VarED))then
      write(IOdbg,*)'Writing tecplt turbulent kinetic energy on vertices'
 
      call InterpolateData(0,VarTE,TE,dPdX,NodalData1,NodalCounter)

      do i=1,Nvrt,10
        k = min(i+9,Nvrt)
        write(IOpst,'(10(1x,1pe11.4))')  (max(0.0,Nodaldata1(j)),j=i,k)
      end do

      write(IOdbg,*)'Writing tecplt turbulent dissipation on vertices'
  
      call InterpolateData(0,VarED,ED,dPdX,NodalData1,NodalCounter)

      do i=1,Nvrt,10
        k = min(i+9,Nvrt)
        write(IOpst,'(10(1x,1pe11.4))')  (max(0.0,Nodaldata1(j)),j=i,k)
      end do
    endif  

    if( PostV(VarVIS) )then
      write(IOdbg,*)'Writing tecplt effective viscosity on vertices'

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
 
          write(IOdbg,*)'Writing tecplt scalar(s) ',&
                                    Variable(NVar+is),' on vertices'

          do i=1,NVrt,10
            k = min(i+9,Nvrt)
            write(IOpst,'(10(1x,1pe10.3))') (max(0.0,Nodaldata1(j)),j=i,k)
          end do
        end do
      endif
    endif

   !
   ! the TECPLT file format needs the topology of vertices
   !
   call openfile(IOcel,casename,'.cel','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)
 1 format(8(1x,i4))
 2 format(8(1x,i5))   
 3 format(8(1x,i6))   
 4 format(8(1x,i7))  
 5 format(8(1x,i8))  
 
   if( Nvrt <= 99999999 ) assign 5 to ifmt 
   if( Nvrt <=  9999999 ) assign 4 to ifmt 
   if( Nvrt <=   999999 ) assign 3 to ifmt 
   if( Nvrt <=    99999 ) assign 2 to ifmt 
   if( Nvrt <=     9999 ) assign 1 to ifmt 

   do i=1,Ncel
     read(IOcel,*) idummy,(indx(j),j=1,8)
     write(IOpst,ifmt)indx(1),indx(2),indx(3),indx(4), &
                      indx(5),indx(6),indx(7),indx(8)
   end do
   
   close(IOcel)   

   !
   !  for particle
   !
   if( UseParticles .and. Model==1 )then
     write(*,*)        'Writing Particle data'
     write(IOdbg,*)    'Writing Particle data'
     do i=1,Npart
       if( .not. associated( Tracks(i)%head ) )then
         write(*,*)'ParticlePrint: error head not associated',i
         cycle
       endif
       write(string2,'(I12)')i
       string2='"ParticleId='//string2(1:lens(string2))//'"'
       write(IOpst,'(A)')PltVar(1:lens(PltVar)) 
       write(IOpst,'(A,1x,I8,1x,A,1x,I8,1x,A)')  &
             'ZONE T='//string2(1:lens(string2))//' I=',Tracks(i)%n, &
             ', J=1, K=1, F=POINT'
       Track => Tracks(i)%head
   
       do while( associated(Track) )
         write(IOpst,'(15(1x,1pe10.3))') Track%x,Ptmp(1:Iuse)
         Track => Track%next
       end do
     end do
   endif

   if( allocated( NodalData1   ) ) deallocate ( NodalData1   )
   if( allocated( NodalData2   ) ) deallocate ( NodalData2   )
   if( allocated( NodalData3   ) ) deallocate ( NodalData3   )
   if( allocated( NodalCounter ) ) deallocate ( NodalCounter )
   if( allocated( PLTcells     ) ) deallocate ( PLTcells     )

   !
   ! do not forget to close, the unit number IOpst might be reused
   !
   close(IOpst)

end subroutine dolfyn2tecplt
