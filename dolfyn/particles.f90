!
! Copyright 2006-2007 Henk Krus, Cyclone Fluid Dynamics BV
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
! With contributions by Shibo (Harry) Kuang
!
module particles

   integer :: Npart = 0

   type ParticleProp
     real    :: x0(3) = 0.0                        ! starting coordinates
     real    :: v0(3) = 0.0                        ! starting velocity
     real    :: x1(3) = 0.0                        ! current coordinates
     real    :: v1(3) = 0.0                        ! current velocity
     real    :: a1(3) = 0.0                        ! current acceleration
     integer :: cell0 =  0                         ! starting cell the particle is in
     integer :: cell1 =  0                         ! current cell the particle is in
     integer :: face  =  0                         ! current face the particle is on
     real    :: dens0 = -1.0  ! <= undefined flag  ! current density
     real    :: diam0 =  1.0
     real    :: mass0 =  1.0
     logical :: wall  = .false.
   end type

   type ParticleTrackData
     real    x(3)                                       ! coordinates
     real    v(3)                                       ! velocity
     integer cell                                       ! cell the particle is in
     real    time                                       ! travelling time since start
     type(ParticleTrackData), pointer :: next => NULL() ! list members
   end type

   type ParticleList
     type(ParticleTrackData), pointer :: head => NULL() ! start of list i
     type(ParticleTrackData), pointer :: last => NULL() ! last entry
     integer                          :: n              ! number of entries
   end type ParticleList

   type(ParticleList), allocatable    :: Tracks(:)      ! all the tracks
   
   type(ParticleTrackData), pointer   :: Track, TrackEntry

   type(ParticleProp), allocatable    :: Particle(:)

   type(ParticleProp), allocatable    :: ParticleTmp(:) ! temporary array to store

end module particles
!
! ***********************************************************************
!
real function Cd(Re)
!========================================================================

   real, intent(IN) :: Re

   if( Re <= 0.001 )then                     !
     Cd = 24000.                             ! some upper limit
   else if( Re <= 0.2 )then                  !
     Cd = 24./Re                             ! Stoke's
   else if( Re <= 1000. )then                !
     Cd = 24.*( 1.0 + 0.15*Re**(0.687))/Re   ! Newton
   else if( Re <= 200000. )then              !
     Cd = 0.44                               ! before, and
   else                                      !
     Cd = 0.10                               ! after supercritical
   endif                                     !

end function Cd
real function ParticleDTcourant(Vel,Nodes,Nnod,MaxNod)
!========================================================================
!
!  Courant number: Co = V dt / dX 
!
!  or: dt = || Co dX / V ||
!
   use constants, only: Large, Small, IOdbg, ParticleCourant 
   use geometry,  only: Vert

   !use particles, only: ParticleCourant
   
   integer, intent(IN)            :: Nodes(1:MaxNod) ! array of nodes of a cell
   real, dimension(3), intent(IN) :: Vel
   integer, intent(IN)            :: Nnod, MaxNod        

   real, dimension(3) :: DTc, DXc
   integer            :: i
   
   xmin =  Large
   xmax = -Large
   ymin =  Large
   ymax = -Large
   zmin =  Large
   zmax = -Large

   DTc = 0.0
   do i=1,Nnod
     xmin = min(xmin,Vert(Nodes(i),1))
     xmax = max(xmax,Vert(Nodes(i),1))
     ymin = min(ymin,Vert(Nodes(i),2))
     ymax = max(ymax,Vert(Nodes(i),2))
     zmin = min(zmin,Vert(Nodes(i),3))
     zmax = max(zmax,Vert(Nodes(i),3))
   end do

   DXc(1) = xmax - xmin
   DXc(2) = ymax - ymin
   DXc(3) = zmax - zmin

   Courant = ParticleCourant                         ! <= in input deck!
   
   DTc = DXc * Courant/( abs(Vel) + Small )

   ParticleDTcourant = minval( DTc )
   !write(*,*)'dT courant:',ParticleDTcourant   
 
end function ParticleDTcourant
subroutine ParticleSetup
!========================================================================

   use constants
   use geometry 
   use particles

   real    :: Xp(3)
   integer :: node(4)
   
   real, allocatable    :: CellXYZ(:,:)        ! cell bounding box 
   real, allocatable    :: DomXYZ(:,:)         ! subdomain bounding box 

   logical :: flag, found, ParticleInCell, ParticleOnFace
   character*10 :: note

   real, parameter    :: Pi = 3.1415927
   !
   ! linked list for cells mapped to subdomains
   !
   type DomEntry
     integer                 :: cell
     type(DomEntry), pointer :: next => NULL()
   end type

   type DomPtr
     type(DomEntry), pointer   :: head
     type(DomEntry), pointer   :: last
   end type
         
   type(DomPtr), allocatable   :: DomMap(:)
   type(DomEntry), pointer     :: cptr, tptr, new
   
   write(*,*)     'Setting up particle data'
   write(IOdbg,*) 'Setting up particle data'    

   !
   ! copy the initial input data
   ! scan the array first backwards for last undefined particle 
   ! (density = < 0.0)
   !
   do i=Npart,1,-1
     if( Particle(i)%dens0 > 0.0 ) exit
   end do
   LastParticle = i

   if( LastParticle > 0 )then
     allocate(ParticleTmp(LastParticle),stat=istat)
     call TrackMemory(istat,LastParticle,'ParticleTmp array allocated')

     ParticleTmp(1:LastParticle) = Particle(1:LastParticle)
     
     deallocate(Particle,stat=istat)
     call TrackMemory(istat,-Npart,'Particle array deallocated')

     allocate(Particle(NPart),stat=istat)
     call TrackMemory(istat,Npart,'New particle array allocated')
     
     Particle(1:LastParticle) = ParticleTmp(1:LastParticle)

     deallocate(ParticleTmp,stat=istat)
     call TrackMemory(istat,-LastParticle,'ParticleTmp array deallocated')
     
     Npart = LastParticle
   endif
   
   !
   ! check for holes??
   !
   icnt = count( Particle%dens0 < 0.0 )  
   if( icnt > 0 )then
     write(*,*)    'Particle array has undefined particles'
     write(IOdbg,*)'Particle array has undefined particles'
   endif

   do i=1,Npart
     if( Particle(i)%dens0 > 0.0 )then
       Particle(i)%x1 = Particle(i)%x0
       Particle(i)%v1 = Particle(i)%v0
       Particle(i)%mass0 = Particle(i)%dens0 * Pi/6.0 * Particle(i)%diam0**3 
     endif
   end do

   write(IOdbg,'(1x,3A)') &
                 '==========================================',&
                 '========= Particle data ===========',&
                 '==================================='
   write(IOdbg,'(1x,3A)') &
                 'particle   dens0      diam0      mass0   |',&
                 '     x0         y0         z0     |',&
                 '     u0         v0         w0     |'
   do i=1,Npart
     if( Particle(i)%dens0 > 0.0 )then
   write(IOdbg,'(i8,3(1x,1pe10.3),'' |'',3(1x,1pe10.3),'' |'',3(1x,1pe10.3),'' |'')') &
       i,Particle(i)%dens0,Particle(i)%diam0,Particle(i)%mass0, &
         Particle(i)%x0(1:3),Particle(i)%v0(1:3)
     else
   write(IOdbg,'(i8,3(1x,1pe10.3),'' |'',3(1x,1pe10.3),'' |'',3(1x,1pe10.3),'' |'')') &
       i
     endif
   end do
   !
   ! find bounding box of complete domain
   !
   XminDom   = minval(vert(:,1))
   XmaxDom   = maxval(vert(:,1))
   
   YminDom   = minval(vert(:,2))
   YmaxDom   = maxval(vert(:,2))

   ZminDom   = minval(vert(:,3))
   ZmaxDom   = maxval(vert(:,3))

  !VolDom    = (XmaxDom - XminDom) * &
  !            (YmaxDom - YminDom) * &
  !            (ZmaxDom - ZminDom) 

   VolModel  = sum(cell(:)%vol)
   
   VolCell   = VolModel / float(Ncel)   ! average cell volume
   DivLength = VolCell**(1./3.)         ! average cell dimension
   
  !Nsplit    = 20
   
   NxDom     = max(1,floor((XmaxDom - XminDom)/DivLength))
   NyDom     = max(1,floor((YmaxDom - YminDom)/DivLength))
   NzDom     = max(1,floor((ZmaxDom - ZminDom)/DivLength))
   
   NDom = NxDom * NyDom * NzDom

   !write(*,*) 'div:',divlength
   !write(*,*) 'nxyz:',nxdom,nydom,nzdom,':',Ndom

   !
   ! something clever here or set in input deck?
   !
   NxDom = 4
   NyDom = 4
   NzDom = 1

   !
   ! in order to be sure that all cells fit;
   ! expand the domain by a small amount
   !
   !dx    = 0.01*(XmaxDom - XminDom) 
   !dy    = 0.01*(YmaxDom - YminDom) 
   !dz    = 0.01*(ZmaxDom - ZminDom) 

   !XminDom = XminDom - dx
   !XmaxDom = XmaxDom + dx
   !YminDom = YminDom - dy     NIET NODIG: gebruik cell center
   !YmaxDom = YmaxDom + dy
   !ZminDom = ZminDom - dz
   !ZmaxDom = ZmaxDom + dz

   DxDom = (XmaxDom - XminDom)/float(NxDom)
   DyDom = (YmaxDom - YminDom)/float(NyDom)
   DzDom = (ZmaxDom - ZminDom)/float(NzDom)

   NDom = NxDom * NyDom * NzDom
   !
   ! create subdomains
   !
   allocate(DomXYZ(NDom,6),stat=istat)          
   call TrackMemory(istat,NDom*6,'Subdomain bounding box array allocated')

   id = 0
   do i=1,NxDom
     do j=1,NyDom
       do k=1,NzDom
         id = id + 1
         DomXYZ(id,1) = XminDom + float(i-1)*DxDom
         DomXYZ(id,2) = XminDom + float( i )*DxDom
         DomXYZ(id,3) = YminDom + float(j-1)*DyDom
         DomXYZ(id,4) = YminDom + float( j )*DyDom
         DomXYZ(id,5) = ZminDom + float(k-1)*DzDom
         DomXYZ(id,6) = ZminDom + float( k )*DzDom
       end do
     end do
   end do 

   write(IOdbg,*)'Particle search domains'
   do i=1,NDom
     write(IOdbg,1) i,DomXYZ(i,:)
   end do   
 1 format(i3,'   ',1pe10.3,' < x < ',e10.3,  &
             ' , ',1pe10.3,' < y < ',e10.3,  &
             ' , ',1pe10.3,' < z < ',e10.3)
   !
   ! start search for starting cell
   !
   allocate(CellXYZ(Ncel,6),stat=istat)          
   call TrackMemory(istat,Ncel*6,'Cell bounding box array allocated')

   allocate(DomMap(NDom),stat=istat)
   call TrackMemory(istat,NDom*2,'Domain map array allocated')
   
   do i=1,NDom
     nullify(DomMap(i)%head)
     nullify(DomMap(i)%last)
   end do

   do i=1,Ncel

     xmin =  Large
     xmax = -Large
     ymin =  Large
     ymax = -Large
     zmin =  Large
     zmax = -Large

     do j=1,Nfaces(i)
       k  = CFace(i,j)
       
       node(1) = Face(k)%vertices(1)
       node(2) = Face(k)%vertices(2)
       node(3) = Face(k)%vertices(3)
       node(4) = Face(k)%vertices(4)
   
       do indx=1,4
         xmin = min(xmin,Vert(node(indx),1))
         xmax = max(xmax,Vert(node(indx),1))
         ymin = min(ymin,Vert(node(indx),2))
         ymax = max(ymax,Vert(node(indx),2))
         zmin = min(zmin,Vert(node(indx),3))
         zmax = max(zmax,Vert(node(indx),3))       
       end do
     end do
     
     CellXYZ(i,1) = xmin
     CellXYZ(i,2) = xmax
     CellXYZ(i,3) = ymin
     CellXYZ(i,4) = ymax
     CellXYZ(i,5) = zmin
     CellXYZ(i,6) = zmax

     !
     ! allocate this cell to a subdomain
     ! according to position of cell center
     !
     do id=1,NDom
     
       if( Cell(i)%x(1) <  DomXYZ(id,1) .or. &
           Cell(i)%x(1) >= DomXYZ(id,2) ) cycle
       if( Cell(i)%x(2) <  DomXYZ(id,3) .or. &
           Cell(i)%x(2) >= DomXYZ(id,4) ) cycle
       if( Cell(i)%x(3) <  DomXYZ(id,5) .or. &
           Cell(i)%x(3) >= DomXYZ(id,6) ) cycle
       
       if( .not. associated( DomMap(id)%head ) )then
         allocate(DomMap(id)%head)                  ! create new 'head' of list
         allocate(DomMap(id)%last)                  ! create entry to store a 'previous' entry

         DomMap(id)%head%cell = i                   ! first cell in 'head'

         DomMap(id)%last%next => DomMap(id)%head    ! store address of this 'head'
         
!         cell(i)%ctid = id                   
       else
         cptr => DomMap(id)%last%next               ! retrieve last stored adress
         
         allocate(new)                              ! allocate 'new' storage

         new%cell = i                               ! store cell number in 'new'

         cptr%next => new                           ! store 'new' adress in previous 'next'
         DomMap(id)%last%next => new                ! remember address of 'new'
         
!         cell(i)%ctid = id                  
       endif
       
     end do
   end do
   !
   ! find starting cell 
   ! todo: check if it is in a FLUID cell!
   !
   do ipart=1,Npart

     if( Particle(ipart)%dens0 < 0.0 ) cycle      ! undefined particle

     Xp = Particle(ipart)%x0
     !write(IOdbg,*) 'Particle ',ipart,':',Xp

     found   = .false.                        ! flag
     factor  = 0.0                            ! factor of subdomain size

     !
     ! first test on bounding box of complete domain
     !
     if( Xp(1) <  (XminDom ) .or. &
         Xp(1) >= (XmaxDom ) ) cycle
     if( Xp(2) <  (YminDom ) .or. &
         Xp(2) >= (YmaxDom ) ) cycle
     if( Xp(3) <  (ZminDom ) .or. &
         Xp(3) >= (ZmaxDom ) ) cycle
 
     do while( .not. found )
       !
       ! first find the subdomain, initially start
       ! with a factor of 0.0 (just hope the particle
       ! is well in a subdomain)
       !     
       do id=1,NDom
         if( Xp(1) <  (DomXYZ(id,1) - factor*DxDom ) .or. &
             Xp(1) >= (DomXYZ(id,2) + factor*DxDom ) ) cycle
         if( Xp(2) <  (DomXYZ(id,3) - factor*DyDom ) .or. &
             Xp(2) >= (DomXYZ(id,4) + factor*DyDom ) ) cycle
         if( Xp(3) <  (DomXYZ(id,5) - factor*DzDom ) .or. &
             Xp(3) >= (DomXYZ(id,6) + factor*DzDom ) ) cycle

         !write(*,*)'current domain:',id
         cptr => DomMap(id)%head

         do while( associated(cptr) )           ! test for end of list

           ic = cptr%cell

           if( Xp(1) < CellXYZ(ic,1) )then
             cptr => cptr%next
             cycle
           elseif( Xp(1) > CellXYZ(ic,2) )then 
             cptr => cptr%next
             cycle
           elseif( Xp(2) < CellXYZ(ic,3) )then 
             cptr => cptr%next
             cycle
           elseif( Xp(2) > CellXYZ(ic,4) )then 
             cptr => cptr%next
             cycle
           elseif( Xp(3) < CellXYZ(ic,5) )then 
             cptr => cptr%next
             cycle
           elseif( Xp(3) > CellXYZ(ic,6) )then 
             cptr => cptr%next
             cycle
           endif
           !
           ! candidate found, check if particle is
           ! really in this cell. if so exit
           !
           !write(*,*) 'checking cell =',ic
           !write(*,*) ip,CellXYZ(ip,1:4)

           flag = ParticleInCell(Xp,ic,idummy,-1.e-6)     ! test if it is in cell ic

           if( flag )then
             found = .true.                       ! found it
             Particle(ipart)%cell0 = ic           ! store the starting cell
             Particle(ipart)%cell1 = ic           ! store the starting cell
             !
             ! check if it starts on a face
             !
             Particle(ipart)%face = 0             ! asume not on a face
             do j=1,NFaces(ic)
               k  = CFace(ic,j)
               if( ParticleOnFace(Xp,ic,k) )then
                 write(IOdbg,*) '*** Particle starts on a face: ',ipart,ic,k
                 Particle(ipart)%face = k
                 exit               
               endif
             end do

             exit                                 ! particle in this cell! stop search
           endif
           !
           ! falls through, keep on searching through the list
           !       
           cptr => cptr%next

         end do ! do-while loop over cell of current subdomain    
         if( found ) exit
       end do ! domain loop
       !
       ! if not found than retry with slightly larger
       ! subdomains
       !
       factor = factor + 0.26
       
       if( factor >= 0.5 )then
         write(IOdbg,*) '*** Particle not found! ',ipart
         found = .true.  ! in order to exit the while
       endif
     end do
   end do
   write(*,*) 'Search completed'

   do i=1,NPart
     if( Particle(i)%dens0 < 0.0 )then
       note = ' UNDEFINED'
     else if( particle(i)%cell0 == 0 )then
       note = ' NOT FOUND'
     else
       note = '          '
     endif
     write(IOdbg,*)'Particle ',i,' starts in cell ',Particle(i)%cell0,note 
     Particle(i)%x1 = Particle(i)%x0
     Particle(i)%v1 = 0.0
   end do

   !
   ! allocate memory for tracks
   !
   allocate(Tracks(NPart),stat=istat)
   
   call ParticleSetUpTracks
   
   do i=1,Npart
     call ParticleAddToTrack(i)
   end do
   !
   ! clean up search lists
   !
   do id=1,NDom
     cptr => DomMap(id)%head
     if( associated(DomMap(id)%last) )then
       deallocate( DomMap(id)%last)
     endif
     do while( associated(cptr) )      ! test for end of list
      !ip = cptr%cell
       tptr => cptr%next               ! first get the next one
       deallocate( cptr )              ! delete where cptr points to
       cptr => tptr
     end do
   end do
   
   deallocate(DomMap,stat=istat)          
   call TrackMemory(istat,NDom*2,'Domain map array allocated')
   deallocate(DomXYZ,stat=istat)          
   call TrackMemory(istat,-NDom*6,'Subdomain bounding box array deallocated')
   deallocate(CellXYZ,stat=istat)
   call TrackMemory(istat,-Ncel*6,'Cell bounding box array deallocated')
      
         
end subroutine ParticleSetup
subroutine ParticleSetUpTracks
!========================================================================

   use particles
   
   write(*,*) 'Setting up tracks'
   
   do i=1,NPart
     if( .not. associated( Tracks(i)%head ) )then
       allocate( Tracks(i)%head )
       allocate( Tracks(i)%last )
     
       Tracks(i)%n = 1

       Tracks(i)%head%x    = Particle(i)%x1
       Tracks(i)%head%v    = Particle(i)%v1
       Tracks(i)%head%cell = Particle(i)%cell0
       !Tracks(i)%head%time = Particle(i)%time

       Tracks(i)%last%next => Tracks(i)%head
     else
       write(*,*)'ParticleSetUpTracks: error track exists:',i
     endif
   end do

end subroutine ParticleSetUpTracks
subroutine ParticleAddToTrack(i)
!========================================================================

   use particles
   
   if( Particle(i)%cell0 == 0 ) return                ! nothing to be done

   if( .not. associated( Tracks(i)%head ) )then
     write(*,*)'ParticleAddToTrack: error head not associated',i
     return
   endif

   Track => Tracks(i)%last%next
   
   allocate( TrackEntry )

   Tracks(i)%n = Tracks(i)%n + 1
   
   TrackEntry%x     = Particle(i)%x1
   TrackEntry%v     = Particle(i)%v1
   TrackEntry%cell  = Particle(i)%cell1
   !TrackEntry%time  = Particle(i)%time

   Track%next => TrackEntry
   
   Tracks(i)%last%next => TrackEntry 
   
end subroutine ParticleAddToTrack
subroutine ParticleClearTrack(i)
!========================================================================

   use particles

   type(ParticleTrackData), pointer :: ptr

   if( .not. associated( Tracks(i)%head ) )then
     write(*,*)'ParticleClearTrack: error head not associated',i
     return
   endif

   !write(*,*)'+++ clearing track',i
   
   Track => Tracks(i)%head%next

   do while( associated(Track) )
     ptr => Track%next
     deallocate( Track )
     Track => ptr
   end do

   Tracks(i)%last%next => Tracks(i)%head
   Tracks(i)%n = 1
   
end subroutine ParticleClearTrack
subroutine ParticlePrint
!========================================================================

   use particles
                                                       ! HARDWIRED
   rewind(8)                                           !<===== NOTE!!!!

   do i=1,Npart
     if( .not. associated( Tracks(i)%head ) )then
       write(*,*)'ParticlePrint: error head not associated',i
       cycle
     endif

     Track => Tracks(i)%head
   
     write(8,*) 'n =',Tracks(i)%n

     do while( associated(Track) )

       write(8,'(i4,6(1x,1pe10.3),i8)') i,Track%x,Track%v,Track%cell

       Track => Track%next
     end do
     
     write(8,'(i4,6(1x,1pe10.3),i8)') -i,0.0,0.0,0.0,0.0,0.0,0.0,0

   end do
   
end subroutine ParticlePrint
logical function ParticleInCell(Xp,ic,ifound,Sloppy)
!========================================================================

   use geometry, only: Nfaces, CFace, Face  

   integer, intent(IN)  :: ic
   real, intent(IN)     :: Xp(3), Sloppy
   integer, intent(OUT) :: ifound

   integer              :: j, k, ip, in
   real                 :: Xf(3), Xn(3)
   logical              :: flag 

   flag   = .true.                          ! assume it is in cell ic
   icnt   = 0
   ifound = -1

   do j=1,NFaces(ic)
     k  = CFace(ic,j)
     ip = Face(k)%cell1
     in = Face(k)%cell2
     if( ip == ic )then 
       Xn =  Face(k)%n 
     else if( in == ic )then
       Xn = -Face(k)%n 
     else
       write(*,*)'Error in particle in cell'
     endif

     Xf = Face(k)%x - Xp                  ! vector from particle to face centre

     call normalise(Xf)
     call normalise(Xn)                   ! normal is not normalised normalize

     dotp  = dot_product( Xf , Xn )       ! dot product

!   if( ic == 338 )write(*,*)'x>',j,k,dotp,sloppy
!   if( ic == 218 )write(*,*)'x2>',j,k,dotp,sloppy
     if( dotp < Sloppy )then
       flag   = .false.                   ! Oops! outside of face! (0.0 is on the face)
       icnt   = icnt + 1
       ifound = k

       !if( ic == 225 )hoek = acos( dotp )*180./3.1415927  
       !if( ic == 225 )write(*,*)'f:',ic,j,dotp,' hoek:', hoek,flag,sloppy  
!   write(*,*)'f:',ic,'=>',j,flag,dotp,acos( dotp )*180./3.1415927,sloppy
     endif

   end do
!   if( flag == .false. )write(*,*)'f:',ic,'=>',flag,k

   !if( icnt == 1 )then
   !  write(*,*)'one face failed'
   !else if( icnt > 1 )then
   !  write(*,*)'multiple faces failed ',icnt
   !  ifound = -1
   !endif

   ParticleInCell = flag

end function ParticleInCell
logical function ParticleOnFace(Xp,ic,k)
!========================================================================

   use constants
   use geometry 

   integer, intent(IN) :: ic, k
   real, intent(IN)    :: Xp(3)

   real                :: Xf(3), Xn(3)
   logical             :: flag 

   flag = .true.                            

   ip = Face(k)%cell1
   in = Face(k)%cell2
   if( ip == ic )then 
     Xn =  Face(k)%n 
   else if( in == ic )then
     Xn = -Face(k)%n 
   endif

   Xf = Face(k)%x - Xp                  ! vector from particle to face centre

   call normalise(Xf)
   call normalise(Xn)                   ! normal is not normalised normalize

   dotp  = dot_product( Xf , Xn )       ! dot product

   !hoek = acos( dotp )*180./3.1415927  
   !write(*,*)'Particle on face:',k,dotp,' hoek:', hoek  

   if( abs(dotp) < 0.01 )then
     flag = .true.
   else
     flag = .false.
   endif
   
   ParticleOnFace = flag

end function ParticleOnFace
subroutine ParticleProjectToFace(Xp,k)
!========================================================================

   use geometry 

   integer, intent(IN) :: k 
   real, intent(INOUT) :: Xp(3)

   double precision    :: Xs(3), Xn(3), dS(3), dotp, s

   Xn = dble( Face(k)%n )
   s  = 1.0/sqrt( dot_product(Xn,Xn) )  ! normal is not normalised normalize
   Xn = s*Xn                            !

   Xs = dble( Xp - Face(k)%x )          ! vector from particle to face centre

   dotp  = dot_product( Xs , Xn )       ! dot product
   
   dS = dotp * Xn 
   Xp = Xp - real( dS )
   
end subroutine ParticleProjectToFace
subroutine ParticleIntersectsFace(found,factor,Xi,Xp,Up,k)
!========================================================================

   use constants
   use geometry 

   integer, intent(IN)   :: k
   real, intent(IN)      :: Xp(3), Up(3)
   logical, intent(OUT)  :: found 
   real, intent(OUT)     :: factor
   real, dimension(3), intent(OUT) :: Xi(3)

   real                  :: Xf(3), Xn(3)

   found  = .false.
   factor = -1.0
   !
   !  Alternative provided by Shibo 'Harry' Kuang (Australia)
   !
   !  eq. of line =>  Xi = Xp + lambda Up
   !
   !  for a point on a surface => ( Xi - Xc ) dot Xn = 0.0       
   !
   !  thus: (Xp + lambda Up - Xc ).Xn = 0.0
   !
   !   =>   (Xp - Xc).Xn + lambda Up.Xn = 0.0
   !
   !   =>   lambda = - ((Xp - Xc).Xn)/(Up.Xn)
   !
   found  = .false.
   factor = -1.0
   Xi     =  0.0

   Xn = Face(k)%n                    ! direction of Xn does not matter
   call normalise( Xn)               ! but normalise it

   dotp = dot_product( Xn , Up )     ! dot product

   if( abs(dotp) > 0.0 )then         !
     found = .true.                  ! an intersection somewhere
   else                              !
     return                          ! more or less parallel to face
   endif                             !

   Xf = Xp - Face(k)%x               ! vector from face centre to intersection

   dotp2 = dot_product( Xf , Xn )    ! dot product

   factor = - dotp2 / dotp           ! lambda and
   Xi     = Xp + factor * Up         ! intersects at Xi

   
end subroutine ParticleIntersectsFace
subroutine ParticleIntersectsFace2(found,factor,Xi,Xp,Up,k)
!========================================================================

   use constants
   use geometry 

   integer, intent(IN)   :: k
   real, intent(IN)      :: Xp(3), Up(3)
   logical, intent(OUT)  :: found 
   real, intent(OUT)     :: factor
   real, dimension(3), intent(OUT) :: Xi(3)

   real                  :: Xf(3), V1(3), V2(3)
   integer               :: node(4)

   real, dimension(3,3)  :: A 
   real, dimension(3)    :: RHSA 
   integer, dimension(3) :: IPIV  ! LAPACK

   found  = .false.
   factor = -1.0
   !
   !  eq. of line   =      eq. of surface
   ! _           _    _          _         _
   ! Xp + lambda Up = Xf + alpha Vf + beta Vf
   !                               1         2
   !
   node(1) = Face(k)%vertices(1)
   node(2) = Face(k)%vertices(2)
   node(4) = Face(k)%vertices(4)

   Xf = Vert(node(1),1:3)
   V1 = Vert(node(2),1:3)
   V2 = Vert(node(4),1:3)

   V1 = V1 - Xf
   V2 = V2 - Xf

   A(1,1) = Up(1)
   A(2,1) = Up(2)
   A(3,1) = Up(3)

   A(1,2) = -V1(1)
   A(2,2) = -V1(2)
   A(3,2) = -V1(3)

   A(1,3) = -V2(1)
   A(2,3) = -V2(2)
   A(3,3) = -V2(3)

   RHSA(1) = Xf(1) - Xp(1)
   RHSA(2) = Xf(2) - Xp(2)
   RHSA(3) = Xf(3) - Xp(3)
   
   call SGESV( 3, 1, A, 3, IPIV, RHSA, 3, INFO )

   !write(*,*)'IS:',k,':',RHSA,INFO
   
   if( INFO == 0 )then
     found  = .true.
     factor = RHSA(1)
     
     !Xi1 = Xp + RHSA(1)*Up
     !Xi2 = Xf + RHSA(2)*V1 + RHSA(3)*V2

     Xi  = Xp + RHSA(1)*Up
   endif

end subroutine ParticleIntersectsFace2
subroutine ParticleUpdateCell(ipart,icell,iface,Vpt,inewcell)
!========================================================================

   use constants
   use geometry 
   use particles

   integer, intent(IN)  :: ipart,icell,iface
   integer, intent(OUT) :: inewcell
   integer              :: ip, in
   real, dimension(3), intent(IN)  :: Vpt
   real, dimension(3)   :: Xn, Vn

   inewcell = 0
   Vn = Vpt

   ip = Face(iface)%cell1
   in = Face(iface)%cell2

   if( ip == icell )then 
     Xn =  Face(iface)%n 
   else if( in == icell )then
     Xn = -Face(iface)%n 
   else
     write(*,*) &
     'ParticleUpdateCell: internal error, face does not belong to cell'
   endif

   call normalise( Xn )               
   call normalise( Vn )               
   dotp = dot_product( Vn , Xn )      

   !if( ipart == 3 )then
   !  write(*,*)'cell=>',icell
   !  write(*,*)'Xn:',Xn
   !  write(*,*)'Vn:',Vn
   !  write(*,*)' .:',dotp
   !  hoek = acos( dotp )*180./3.1415927  
   !  write(*,*)'hoek:', hoek
   !endif
   
   if( in > 0 )then
     if( ip == icell )then 
       if( dotp > 0 )then
         Particle(ipart)%cell1 = in 
       else
         Particle(ipart)%cell1 = ip 
       endif
 !if( ipart == 3 )then
 !  write(IOdbg,*)'Cell:',icell,' ->',in,' (n)'
 !  write(*,*)'Cell:',icell,' ->',ip,in,' (n)'
 !endif
     else if( in == icell )then
       if( dotp > 0 )then
         Particle(ipart)%cell1 = ip
       else
         Particle(ipart)%cell1 = in 
       endif
 !if( ipart == 3 )then
 !  write(IOdbg,*)'Cell:',icell,' ->',ip,' (p)'
 !  write(*,*)'Cell:',icell,' ->',ip,in,' (p)'
 !endif
     else
       write(IOdbg,*)'CELL:',icell,' =>',ip,in,'<='
       stop 'ParticleUpdateCell: internal error 2'
     endif 

     inewcell = Particle(ipart)%cell1
   else
     Particle(ipart)%wall = .true.         

     ib = Face(iface)%bnd
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == RInlet )then
       write(IOdbg,*)'*** Inlet'
     else if( it == ROutlet )then
       write(IOdbg,*)'*** OUTLET'
     else if( it == RSymp )then
       write(IOdbg,*)'*** Symmetry plane'
     else if( it == RWall )then
       write(IOdbg,*)'*** WALL'
     else
       write(IOdbg,*)'+ internal error ParticleUpdateCell'     
     endif
  
   endif

end subroutine ParticleUpdateCell
subroutine ParticleGetVelocity(Xp,Up,Nodes,Wnode,Nnod,MaxNod,Vvrt)
!========================================================================

   use constants
   use geometry 
   use variables

   real, intent(IN)    :: Vvrt(1:Nvrt,1:3)           ! array of data at vertices
   integer, intent(IN) :: Nodes(1:MaxNod)            ! array of nodes of a cell
   integer, intent(IN) :: Nnod, MaxNod               !  

   real                :: Wnode(1:MaxNod)            ! array of nodal weights
   real, intent(IN)    :: Xp(3) 
   real, intent(OUT)   :: Up(3)
   real                :: Xc(3)
   
   Sum = 0.0
   do i=1,Nnod
     Xc   = Xp - Vert(Nodes(i),1:3)          ! vector particle to node
     Wnode(i) = dot_product(Xc,Xc)           ! length squared
     Wnode(i) = 1.0/(Wnode(i)+Small)         ! weight function
     Sum = Sum + Wnode(i) 
   end do 

   Wnode(1:Nnod) = Wnode(1:Nnod)/Sum      
   
   Up = 0.0
   do i=1,Nnod
     Up(1) = Up(1) + Wnode(i)*Vvrt( Nodes(i) ,1)
     Up(2) = Up(2) + Wnode(i)*Vvrt( Nodes(i) ,2)
     Up(3) = Up(3) + Wnode(i)*Vvrt( Nodes(i) ,3)
   end do

end subroutine ParticleGetVelocity 
subroutine CollectCellNodes(icell,Nodes,Nnod,MaxNod)
!========================================================================

   use constants
   use geometry 

   integer, intent(IN)  :: icell                    ! cell
   integer, intent(IN)  :: MaxNod                   ! still hard wired!!!!

   integer, intent(OUT) :: Nodes(1:MaxNod)          ! array of nodes of a cell
   integer, intent(OUT) :: Nnod                     ! number of nodes

   Nnod = 0
   Faces: do j=1,Nfaces(icell)
     k  = CFace(icell,j)
     FaceNodes: do l=1,4
       if( Face(k)%vertices(l) < 0 ) cycle
       do m=1,Nnod
         if( Nodes(m) == Face(k)%vertices(l) ) exit FaceNodes ! already found
       end do
       Nnod = Nnod + 1
       Nodes(Nnod) = Face(k)%vertices(l)              
     end do FaceNodes
   end do Faces

end subroutine CollectCellNodes
integer function ParticleInFace(Xp,Up,iface)
!========================================================================
!  to determine whether particle is inside a face with bound
!              B    --------------- 
!                  /|    <---     /   
!                 / | b   . p(t+dt)          b=B-p(t)
!                /  |    /      /            a=A-p(t)
!             A  ---------------             t=p(t+dt)-p(t)
!               \   |  / t
!              a \  | /              
!                   .  p(t)
!
! dotp=(aXb).t  if all dotp<0, then particle will intersect with face                  


   use geometry 

   integer, intent(IN) :: iface
   real, intent(IN)    :: Xp(3),Up(3)

   integer             :: i1,i21,i22,Npoint, flag
   real                :: Line1(3), Line2(3), Tmp(3)

   flag = 1                           ! assume it is inside face iface

   if( Face(iface)%vertices(4) > 0 )then                
     Npoint = 4
   else
     Npoint = 3
   endif

   do i1=1,Npoint                                      
     if( i1+1 <= Npoint )then
       i21  = Face(iface)%vertices(i1)
       i22  = Face(iface)%vertices(i1+1)
     else if( i1 == Npoint )then
       i21  = Face(iface)%vertices(i1)
       i22  = Face(iface)%vertices(1)
     endif                   

     Line2 = Vert(i22,:)-xp
     Line1 = Vert(i21,:)-xp
           
     call cross_product(line2,line1,tmp)
           
     dotp= dot_product(tmp,up)

     if( dotp > 0.0 )then
       flag = 0
       exit
     endif

   end do    
       
   ParticleInFace = flag

end function ParticleInFace
subroutine ParticleSearchCell(level,icell,Xp,inewcell,inewface)
!====================================================================================
 
   use constants, only: IOdbg
   use geometry, only:  Nfaces, CFace

   integer, intent(IN)  :: level, icell
   real, intent(IN)     :: Xp(3)
   integer,intent(OUT)  :: inewcell, inewface

   integer, dimension(100) :: CellList                     ! <= HARDWIRED

   logical :: flag, ParticleInCell, ParticleOnFace

   inewface = -1         ! flag as not found (yet)
   inewcell = -1
   !
   ! run over the neighboring cells only
   !
   call CollectCells(1,icell,nc,CellList)   

! write(*,*)'list:',icell,'>',CellList(1:nc)
!if( icell == 338 ) write(*,*)'list:',icell,'>',CellList(1:nc)
!if( icell == 210 ) write(*,*)'list:',icell,'>',CellList(1:nc)
!if( icell == 218 ) write(*,*)'list:',icell,'>',CellList(1:nc)

   do i=1,nc

!write(*,*)'xxxx:',i,nc,icell,'>',CellList(1:nc)

     ip = CellList(i)

!if( ip == 218 ) write(*,*)'list2:',icell,'>',CellList(1:nc)
     
     flag  = ParticleInCell(Xp,ip,idummy,0.0)        

!     write(*,*)'>>>>',icell,nc,i,'>',ip, flag,idummy
     !write(*,*)'....',Xp
     
     if( flag )then
       inewcell = ip
       !
       ! loop over faces 
       !
       do j=1,NFaces(inewcell)
         k  = CFace(inewcell,j)
!write(*,*)'....',j,k
         if( ParticleOnFace(Xp,inewcell,k) )then
!      write(*,*)'..HEBBES'
           inewface = k
           exit   
         endif            
       end do
       exit
     endif
   end do
!write(*,*)'hier'
   if( inewcell == -1 )then
     !
     ! using tight tolerance of 0.0 unsuccesful
     ! try a slight tolerance
     !
     write(IOdbg,*)'ParticleSearchCell using tolerance'
     write(*,*)'ParticleSearchCell using tolerance'
     tol = -1.e-6
     do while( inewcell == -1 )
       do i=1,nc

         ip = CellList(i)

         flag  = ParticleInCell(Xp,ip,idummy,tol)        

!         write(*,*)'>>>>',icell,nc,i,'>',ip, flag,idummy
         !write(*,*)'....',Xp

         if( flag )then
           inewcell = ip
           !
           ! loop over faces 
           !
           do j=1,NFaces(inewcell)
             k  = CFace(inewcell,j)
             if( ParticleOnFace(Xp,inewcell,k) )then
               inewface = k
               exit   
             endif            
           end do
           exit
         endif
       end do
       
       if( tol < -1.e-2 )then
         exit
       else
         tol = 10.0*tol
         write(IOdbg,*)'ParticleSearchCell set tolerance to ',tol
       endif
     end do
   endif
 write(*,*)'einde zoek'   
end subroutine ParticleSearchCell
!
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! 
! hack hoek
!
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!
subroutine ParticleInterface
!========================================================================

   use constants
   use geometry 
   use variables
   use particles
   
   !real    :: Pi = 3.1415927
   real    :: vector_length
   
   real    :: Xp(3), Xn(3)
   real    :: Up(3), Vp(3), Xi(3), Xtmp(3), Utmp(3)
   real    :: Upt(3), Vpt(3), Xpt(3), dX(3), X1(3), Vi(3), dXi(3)
    
   real    :: V1(3), V2(3), Vslip(3), Vf(3)

   logical :: Done, Found, flag, NewCell, OnFace, TinyStep
   logical :: ParticleInCell, ParticleOnFace

   real, allocatable    :: Vvrt(:,:)           ! array of data at vertices
   integer, allocatable :: NodalCounter(:)  
   integer, allocatable :: Nodes(:)            ! array of nodes of a cell
   real, allocatable    :: Wnode(:)            ! array of nodal weights
   
   write(*,*) 'Interpolating velocities'

   allocate( Vvrt(Nvrt,3),stat=istat)    
   allocate( NodalCounter(Nvrt),stat=istat)    

   MaxNod = 20                                   !<====== HARDWIRED

   allocate( Nodes(MaxNod),stat=istat)    
   allocate( Wnode(MaxNod),stat=istat)    
   call TrackMemory(istat,4*Nvrt+2*MaxNod,'Work arrays for nodal velocity components allocated')

   call InterpolateData(1,VarU,U,dUdX,Vvrt(:,1),NodalCounter)
   call InterpolateData(1,VarV,V,dVdX,Vvrt(:,2),NodalCounter)
   call InterpolateData(1,VarW,W,dWdX,Vvrt(:,3),NodalCounter)

   write(*,*) 'Done. Start tracking ',Npart

   DtMax   = 5.e-4
   TimeMax =  50

   IPDBG = -1

   do ipart=1,Npart
 
     if( Particle(ipart)%dens0 < 0.0   ) cycle     ! undefined particle

     if( Particle(ipart)%diam0 < 1.e-8 ) cycle     ! unphysical particle

     if( Particle(ipart)%cell0 == 0    ) cycle     ! no starting point

     write(*,*) '*** particle:',ipart
     write(IOdbg,*) '*** particle:',ipart
     call ParticleClearTrack(ipart)

     Particle(ipart)%x1    = Particle(ipart)%x0
     Particle(ipart)%v1    = Particle(ipart)%v0
     Particle(ipart)%cell1 = Particle(ipart)%cell0
     
     Done    = .false.
     NewCell = .true.
     if( Particle(ipart)%face == 0 )then
       OnFace  = .false.
     else
       OnFace  = .true.
     endif

     icnt        =  0
     iface       = -1
     Vpt         = 0.0
     lastcell    = -1
     LastFace    = -1

     do while( .not. done )
       icnt = icnt + 1
       ic = Particle(ipart)%cell1                ! is currently in cell ic
       Xp = Particle(ipart)%x1                   ! get current position
       Vp = Particle(ipart)%v1                   ! get current velocity

       DenF  = Den(ic)                           !
       DenP  = Particle(ipart)%dens0             ! some data
       Diam  = Particle(ipart)%diam0             !
       VisF  = VisLam                            ! has to be VisEff or not?


       if( NewCell )then

         if( ipart == IPDBG ) write(*,*)'cell:',ic,icnt
         icntcell = 0

         call CollectCellNodes(ic,Nodes,Nnod,MaxNod)

       else
       
         icntcell = icntcell + 1
         
         if( icntcell > 2 ) LastFace    = -1

       endif

       call ParticleGetVelocity(Xp,Up,Nodes,Wnode,Nnod,MaxNod,Vvrt)

       vel12 = dot_product(Up,Up)              ! starting fluid vel. squared
       vel22 = dot_product(Vp,Vp)              ! starting particle vel. squared
       if( vel12 < vel22 )then
         Utmp = Vp
       else
         Utmp = Up
       endif
       !
       ! establish DT courant 
       ! based on the highest velocity  
       !
       DTcourant = ParticleDTcourant(Utmp,Nodes,Nnod,MaxNod)

       call ParticleGetVelocity(Xp,Upt,Nodes, &
                                Wnode,Nnod,MaxNod,Vvrt) ! Upt tpv Xpt

       Vslip = Upt - Vp 
       Vmag  = sqrt( dot_product(Vslip,Vslip) )
       Re    = DenF * Vmag * Diam / VisF
   
       c1  = 3./4. *( Cd(Re) * DenF * Vmag )/ (DenP + 0.5*DenF) /Diam 
       c2  = (DenP - DenF) / (DenP + 0.5*DenF)
       c3  = c2/(c1 + Small)

       V1  = Upt + c3 * Gravity
       V2  = V1 - Vp

       if( Re < 0.2 )then
         TauM = (DenP + 0.5*DenF)*Diam**2 / (18. * VisF)
       else
         TauM = 1.0/( c1 + Small )
       endif

       V1mag  = sqrt( dot_product(V1,V1) )
       Vpmag  = sqrt( dot_product(Vp,Vp) )
       vrel   = vpmag/(v1mag+Small)
       
       if( vrel > 0.98 .and. vrel < 1.02 ) TauM = DTcourant
              
       dtmp = min( DTcourant, TauM)

       Vpt = V1 - V2*exp( -c1 * Dtmp )

       !dX  = Vpt * DTmp
       dX  = V1*Dtmp - V2/c1 * ( 1.0 - exp(-c1*Dtmp) )

       X1  = Xp + dX

       dS  = vector_length( dX ) 
       if( dS < 1.e-12 )then
         write(*,*)'Particle hardly moves. Stop tracking it.'
         write(*,*)'dS = ',dS
         Done = .true.
         cycle
       endif
       
       !
       ! nog in de cell? zal vaak het geval zijn
       !
       flag = ParticleInCell(X1,ic,idummy,0.0)

       if( flag )then
         !write(*,*)'In cell:',ic
         !
         ! still in this cell
         !
         icntcell = icntcell + 1
         
         if( icntcell > 10 )write(*,*)'=>', &
           Vpt(1),up(1),DS,dtmp,DTcourant,taum,vrel

         if( icntcell == 400 )then
           write(IOdbg,*)'Particle stuck in',ic
           write(*,*)'Stuck in ',ic,ipart,Vpt
           Done = .true.
           cycle
         endif
         
         Particle(ipart)%face = 0
         Particle(ipart)%x1   = X1
         Particle(ipart)%v1   = Vpt
        
         call ParticleAddToTrack(ipart)
        
         LastFace = -1
         NewCell  = .false.
         OnFace   = .false.
         
       else
         !
         ! waar is ie? waar gaat ie heen?
         !
         if( ipart == IPDBG ) &
           write(*,*)'NIET MEER in cell:',ic,' part:',ipart
         OnFace   = .false.
         TinyStep = .false.
         factor0  = 1.e12
         dS       = 1.e12
         
         Vi = dX/Dtmp
         
         do j=1,NFaces(ic)
           k  = CFace(ic,j)

           call ParticleIntersectsFace(found,factor1,Xtmp,Xp,Vi,k)
           !write(*,*)'N>',j,k,found,factor1,Face(k)%cell1,Face(k)%cell2
 
           if( k == LastFace .and. abs(factor1) < 1.e-6  ) cycle

           Epsilon = 1.e-12

           if( found .and. factor1 > Epsilon )then
             if( abs(factor1-factor0) < Epsilon )then
               dXi = Xtmp - Face(k)%x
               if( vector_length( dXi ) < dS )then
                 factor0 = factor1
                 iface   = k
                 Xi      = Xtmp
                 dS      = vector_length( dXi ) 
               endif
             elseif( factor1 < factor0 )then
               dXi = Xtmp - Face(k)%x
               factor0 = factor1
               iface   = k
               Xi      = Xtmp
               dS      = vector_length( dXi ) 
             endif
           endif
           
         end do

         !
         ! test if near a face
         !
         if( found .and. factor0 < 0.01 * DTmp  )then

           TinyStep = .true.
           if( ipart == IPDBG ) write(*,*)'PART on face:',factor0

         endif
         
         if( factor0 < Dtmp .and. .not. TinyStep )then
           if( ipart == IPDBG ) write(*,*)'End step',Dtmp,ic
           Xpt = Xi
           Vpt = V1 - V2*exp( -c1 * factor0 ) 
           call ParticleProjectToFace(Xpt,iface) 

           LastFace = iface

           Particle(ipart)%face = iface
           Particle(ipart)%x1   = Xpt
           Particle(ipart)%v1   = Vpt

           call ParticleAddToTrack(ipart)
           call ParticleUpdateCell(ipart,ic,iface,Vpt,icnew)

           NewCell = .true.
           OnFace  = .true.

         else if( TinyStep )then
           if( ipart == IPDBG ) write(*,*)'Tiny step',Dtmp,ic

           Xpt = Xp
           Vpt = Vp
           
           call ParticleProjectToFace(Xpt,iface) 
 
           LastFace = iface

           Particle(ipart)%face = iface
           Particle(ipart)%x1   = Xpt
           Particle(ipart)%v1   = Vpt

           call ParticleUpdateCell(ipart,ic,iface,Vpt,icnew)

           NewCell = .true.
           if( ic == icnew ) NewCell = .false.
           OnFace  = .true.

         else if( factor0 > 1.e10 )then

           write(IOdbg,*)'ABORT TRACK',ipart
           write(*,*)'ABORT TRACK',ipart
           Done = .true.
           cycle

         else
           write(IOdbg,*)'Particle lost',ipart
           write(*,*)'kwijt',ipart

           call ParticleSearchCell(0,ic,Xp,inewcell,inewface)

           write(*,*)'gevonden',ic,inewcell,inewface
           if( inewcell > 0 )then

             if( inewface == -1 ) inewface = 0

             Particle(ipart)%cell1 = inewcell
             Particle(ipart)%face  = inewface
             Particle(ipart)%x1    = X1
             Particle(ipart)%v1    = Vpt
        
             call ParticleAddToTrack(ipart)
           
             NewCell = .false.
             OnFace  = .false.

           write(*,*)'verder...',done,ipart
         
           else
             write(IOdbg,*)'LOST',ipart
             write(*,*)'LOST',ipart
             Done = .true.
             cycle           
           endif
         
         endif       

       endif

       if( Particle(ipart)%wall ) Done = .true.
       if( icnt > 1800 )then
         write(*,*)'STUCK', ipart     
         done = .true.
       endif

     end do ! track of one particle
     
     write(IOdbg,*)'Done particle:',ipart,icnt,ic
     write(*,*)'Done particle:',ipart,icnt,ic,Npart

   end do   ! all particles

   if( allocated( Vvrt         ) ) deallocate ( Vvrt         )
   if( allocated( NodalCounter ) ) deallocate ( NodalCounter )
   if( allocated( Nodes        ) ) deallocate ( Nodes        )
   if( allocated( Wnode        ) ) deallocate ( Wnode        )
   
   call ParticlePrint
end subroutine ParticleInterface
