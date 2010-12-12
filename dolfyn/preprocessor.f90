!
! Copyright 2002-2007 Henk Krus, Cyclone Fluid Dynamics BV
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
module constants
   !
   ! global constants
   !
   integer, parameter :: Solid = 1, Fluid = 2, Ghost = 3, Baffle = 4, Shell = 5

   integer, parameter :: Hexa = 1, Prism = 2, Pyramid  = 3, Tetra = 5
   integer, parameter :: Poly = 6, Quad  = 7, Triangle = 8, Polygon = 9

   integer, parameter :: NFhex  = 6, NFprism = 5, NFpyrmd = 5, NFtet = 4
   integer, parameter :: NFpoly = 8, NFquad  = 2, NFtri = 2

   integer, parameter :: IOgeo =  8, IOcel = 14, IOvrt = 15, IObnd = 23
   integer, parameter :: IOdbg = 63, IOcfg = 12

   integer :: debug = 0

   integer, parameter :: version = 0420 
end module constants
!========================================================================
module geometry
   !
   ! the basic geometry related features
   ! the rest of the program should rely on this geometry related
   ! module only for its data
   !
   use constants

   public 
   ! 
   ! general stuff
   !
   integer :: Nvrt = 0, Ncel = 0, Nbnd = 0
   
   !
   ! vertex section
   !
   real, allocatable :: Vert(:,:)   ! array of vertices

   !
   ! face section
   !
   type FaceData
   ! integer shape               ! Quad or Triangle => niet nodig?
     integer bnd                 ! internal (0), boundary (#) ...
     integer cell1               ! first cell
     integer cell2               ! second cell
     integer vertices(4)         ! 
     real    area                ! area
     real    n(3)                ! normal / surface area components
     real    x(3)                ! center coordinates
     real    lambda              ! interpolation factor cell1 -> cell2
   end type

   type(FaceData), allocatable   :: Face(:)    

   real, allocatable :: FaceNormal(:,:) ! Face normals

   !
   !--- for each computational cell there will be a list of faces
   !
   type CellFaceData
     integer face
     type(CellFaceData), pointer :: next
   end type
   
   type CFadresses
     type(CellFaceData), pointer :: ptr
   end type

   type(CFadresses), allocatable :: CFhead(:)
   type(CFadresses), allocatable :: CFtail(:)

   integer, allocatable :: NFaces(:) ! number of faces for a comp. cell
   integer              :: Nfac      ! total number of faces

   !
   !--- for each vertex there will be a list of faces
   !
   type VertFaceData
     integer face
     type(VertFaceData), pointer :: next
   end type
   
   type VFadresses
     type(VertFaceData), pointer :: ptr
   end type

   type(VFadresses), allocatable :: VFhead(:)
   type(VFadresses), allocatable :: VFtail(:)

   integer, allocatable :: NVertFaces(:) ! number of faces for a vertex

   !
   ! cell section
   !
   type :: CellData
     integer, dimension(8) :: vertices   ! the 8 vertices 
     real    x(3)                        ! center coordinates
     real    vol                         ! volume
     integer ctid                        ! fluid type id as set in fluid table
     integer corg                        ! original cell number
   end type
   type(CellData), allocatable :: Cell(:)   ! array of computational cells

   !
   ! boundary section (both set by input and default boundaries)
   !
   type :: BoundaryData
     integer :: face                     ! belongs to face...
     integer, dimension(4) :: vertices   ! the 4 vertices 
     integer :: rid                      ! region id as set in rtable
     real    :: distance                 ! normal distance from cell face 
                                         ! center to cell center
    end type
   type(BoundaryData), allocatable :: Bnd(:) ! array of computational cells
   
   contains
   
   subroutine set_up_add_face_to_cell(n)
   !=====================================
   integer n
   
   if( debug > 0 ) write(*,*)'Setting up list of faces for cells:',n
   
   allocate(CFhead(n))
   allocate(CFtail(n))
   
   do i=1,n
     nullify(CFhead(i)%ptr)
     nullify(CFtail(i)%ptr)
   end do

   end subroutine set_up_add_face_to_cell
   !=====================================
   subroutine add_face_to_cell(cell,face)
   !
   ! add a face to a cell
   !
   integer cell, face
      
   if( .not. associated( CFhead(cell)%ptr ) )then
     ! empty list, just add the face right here
     !if( debug > 0 ) write(*,*)'Initialise CellFace list for cell ',cell
     allocate(CFhead(cell)%ptr, stat=istat)
     CFtail(cell)%ptr => CFhead(cell)%ptr
     nullify( CFtail(cell)%ptr%next )
     CFtail(cell)%ptr%face = face
   else 
     !  write(*,*)'add vertex ',vertex
     allocate(CFtail(cell)%ptr%next, stat=istat)
     CFtail(cell)%ptr => CFtail(cell)%ptr%next
     nullify( CFtail(cell)%ptr%next )
     CFtail(cell)%ptr%face = face
   endif

   end subroutine add_face_to_cell
   !=====================================================================
   subroutine debug_face_to_cell(n)
   integer n   
   type(CFadresses) ptr
 
   do i=1,min(n,10)
       ptr%ptr => CFhead(i)%ptr
       write(*,'(/,a,i8,a,i3,a)',advance='no') ' c',i,' (',NFaces(i),')->'
       do while( associated(ptr%ptr) )
         j = ptr%ptr%face
         write(*,'(i8)',advance='no') j
         ptr%ptr => ptr%ptr%next
       end do
   end do

   if( n > 20 )then
     write(*,'(3(/,''   .''))')   
     do i=n-10,n
         ptr%ptr => CFhead(i)%ptr
         write(*,'(/,a,i8,a,i3,a)',advance='no') ' c',i,' (',NFaces(i),')->'
         do while( associated(ptr%ptr) )
           j = ptr%ptr%face
           write(*,'(i8)',advance='no') j
           ptr%ptr => ptr%ptr%next
         end do
     end do
   endif

   !do i=1,n
   !    ptr%ptr => CFhead(i)%ptr
   !    write(IOdbg,'(/,a,i8,a,i3,a)',advance='no') ' c',i,' (',NFaces(i),')->'
   !    do while( associated(ptr%ptr) )
   !      j = ptr%ptr%face
   !      write(IOdbg,'(i8)',advance='no') j
   !      ptr%ptr => ptr%ptr%next
   !    end do
   !end do

   write(*,'(/,'' Dump done'')')   
   end subroutine debug_face_to_cell   
   
   subroutine set_up_add_face_to_vertex(n)
   !=====================================
   integer n
   
   if( debug > 0 ) write(*,*)'Setting up list of faces for vertices:',n
   
   allocate(VFhead(n))
   allocate(VFtail(n))
   
   do i=1,n
     nullify(VFhead(i)%ptr)
     nullify(VFtail(i)%ptr)
   end do

   end subroutine set_up_add_face_to_vertex
   !=====================================
   subroutine add_face_to_vertex(iv,face)
   !
   ! add a face to a cell
   !
   integer iv, face
      
   !write(*,*)'add vertex/face ',iv,face
   
   if( .not. associated( VFhead(iv)%ptr ) )then
     allocate(VFhead(iv)%ptr, stat=istat)
     VFtail(iv)%ptr => VFhead(iv)%ptr
     nullify( VFtail(iv)%ptr%next )
     VFtail(iv)%ptr%face = face
   else 
     !write(*,*)'add vertex/face ',iv,face
     allocate(VFtail(iv)%ptr%next, stat=istat)
     VFtail(iv)%ptr => VFtail(iv)%ptr%next
     nullify( VFtail(iv)%ptr%next )
     VFtail(iv)%ptr%face = face
   endif

   end subroutine add_face_to_vertex

end module geometry
!========================================================================
module store_vertices
   !
   ! get a vertex, store it in a temporary linked list
   ! when the linked list is filled dump the final
   ! list into an allocated vertex-array
   ! in this way the vertex list mat come in in any
   ! order and may even include 'holes' 
   ! 
   use constants

   private vertex

   integer, private :: min_vertex =  99999999
   integer, private :: max_vertex = -99999999

   type :: vertex
     integer :: index
     real, dimension(3) :: xyz
     type (vertex), pointer :: next
   end type

   type (vertex), private, pointer :: head ! pointer to the head of the list
   type (vertex), private, pointer :: tail ! pointer to the new vertex

   public 

   contains
   subroutine add_vertex(vertex,xyz)
   !
   ! add a new vertex to the list
   !
   integer vertex
   real, dimension(3) :: xyz

   ! write(*,*)'add_vertex: ',vertex

   if( .not. associated(head) )then
     ! empty list, just add the vertex right here
     if( debug > 0 ) write(*,*)'Initialise vertex list'
     allocate(head, stat=istat)
     tail => head
     nullify( tail%next )
     tail%index = vertex
     tail%xyz   = xyz
   else 
     !  write(*,*)'add vertex ',vertex
     allocate(tail%next, stat=istat)
     tail => tail%next
     nullify( tail%next )
     tail%index = vertex
     tail%xyz   = xyz
   endif

   min_vertex = min(min_vertex, vertex)
   max_vertex = max(max_vertex, vertex)

   end subroutine add_vertex
   !======================================================================
   subroutine get_vertices(n,v)
   !
   ! transfer the vertices from the list to the vertex array
   !
   integer n,i,icnt
   real, dimension(n,3) :: v  
   real, dimension(3)   :: temp
   logical, allocatable :: used(:)
   integer :: errors = 0

   type (vertex), pointer :: ptr, tmp

   if( debug > 0 ) write(*,*)'Transfering vertex list into vertex array'
   if( .not. associated(head) )then
     write(*,*)'get_vertices: head of temp. vertex list not associated. stop'
     stop
   else
     ptr  => head                       ! point to start of list
     icnt = 0                           ! set counter to zero
     allocate( used(n) )                ! allocate temp. logical used array
     used = .false.                     ! initialise to .false.
     v    = 0.0                         ! initialise vertex array
     do
       if( .not. associated(ptr) )exit  ! test for end of list
       icnt    = icnt + 1               ! bump counter
       i       = ptr%index              ! get vertex index
       if( used(i) )then                ! test if already used
         errors = errors + 1            ! count the errors
         write(*,*)'+ Error: double index found:',i
       endif
       temp    = ptr%xyz                ! copy xyz-coords into temp
       v(i,1)  = temp(1)                ! set v(i,1), the x-coordinate
       v(i,2)  = temp(2)                ! set v(i,2), the y-coordinate
       v(i,3)  = temp(3)                ! set v(i,3), the z-coordinate
       used(i) = .true.                 ! flag this index as used
       tmp     => ptr                   ! store the current pointer
       ptr     => ptr%next              ! get next item in list
       deallocate(tmp)                  ! remove the current item
     end do    
     deallocate( used )                 ! used is not needed anymore
     if( errors > 0 )then
       write(*,*)'+ Error: number of double indices found:',errors,'. stop.'
       stop      
     endif
     if( debug > 0 ) write(*,*)'Vertices transfered: ',icnt    
   endif

   end subroutine get_vertices
   !======================================================================
   integer function get_min_vertex()
   get_min_vertex = min_vertex
   end function get_min_vertex
   !======================================================================
   integer function get_max_vertex()
   get_max_vertex = max_vertex
   end function get_max_vertex
end module store_vertices
!========================================================================
module store_original_cells
   !
   ! get a cell, store it in a temporary linked list
   ! when the linked list is filled dump the final
   ! list into an allocated cell-array
   ! in this way the cell list may come in any
   ! order and may even include 'holes'
   ! in contrast to the vertices the computational cell-list
   ! may not contain holes; use a map to store the original list
   ! 
   use constants

   private cell_org

   integer, private :: min_cell =  99999999
   integer, private :: max_cell = -99999999

   type :: cell_org
     integer :: index                    ! index as read
     integer, dimension(8) :: vertices   ! the 8 vertices 
     integer :: ctid                     ! ctable id
     type (cell_org), pointer :: next    ! next original cell definition
   end type

   type (cell_org), private, pointer :: head ! pointer to the head of the list
   type (cell_org), private, pointer :: tail ! pointer to the new orginal cell

   public 

   contains

   subroutine add_original_cell(icell,vertices,ctid)
   !
   ! just add a new cell to the list
   !
   integer icell, ctid
   integer, dimension(8) :: vertices

   ! write(*,*)'add_cell: ',icell

   if( .not. associated(head) )then
     ! empty list, just add the cell right here
     if( debug > 0 ) write(*,*)'Initialise original cell list'
     allocate(head, stat=istat)
     tail => head
     nullify( tail%next )
     tail%index    = icell
     tail%vertices = vertices
     tail%ctid     = ctid
   else 
     !  write(*,*)'add vertex ',vertex
     allocate(tail%next, stat=istat)
     tail => tail%next
     nullify( tail%next )
     tail%index    = icell
     tail%vertices = vertices
     tail%ctid     = ctid
   endif

   min_cell = min(min_cell, icell)
   max_cell = max(max_cell, icell)

   end subroutine add_original_cell
   !======================================================================
   subroutine get_original_cells(n,corg)
   !
   ! transfer the cells from the list to the original cell array
   !
   integer n,i,icnt, ctid
   integer, dimension(n,9) :: corg  
   integer, dimension(8)   :: temp
   logical, allocatable    :: used(:)
   integer :: errors = 0

   type (cell_org), pointer :: ptr, tmp

   if( debug > 0 ) write(*,*)'Transfering cell list into original cell array'
   if( .not. associated(head) )then
     write(*,*)'get_original_cells: head of temp. cell list not assiociated. stop'
     stop
   else
     ptr  => head                       ! point to start of list
     icnt = 0                           ! set counter to zero
     allocate( used(n) )                ! allocate temp. logical used array
     used = .false.                     ! initialise to .false.
     corg = 0                           ! initialise corg array
     do
       if( .not. associated(ptr) )exit  ! test for end of list
       icnt    = icnt + 1               ! bump counter
       i       = ptr%index              ! get vertex index
       if( used(i) )then                ! test if already used
         errors = errors + 1            ! count the errors
         write(*,*)'+ Error: double index found:',i
       endif
       temp    = ptr%vertices           ! copy vertices into temp
       corg(i,1) = temp(1)              ! set corg(i,1), first vertex
       corg(i,2) = temp(2)              ! set corg(i,2), 2nd vertex
       corg(i,3) = temp(3)              ! set corg(i,3)    .
       corg(i,4) = temp(4)              ! set corg(i,4)
       corg(i,5) = temp(5)              ! set corg(i,5)    .
       corg(i,6) = temp(6)              ! set corg(i,6)
       corg(i,7) = temp(7)              ! set corg(i,7)    .
       corg(i,8) = temp(8)              ! set corg(i,8), last vertex
       corg(i,9) = ptr%ctid             ! cell id
       used(i) = .true.                 ! flag this index as used
       tmp     => ptr                   ! store the current pointer
       ptr     => ptr%next              ! get next item in list
       deallocate(tmp)                  ! remove the current item
     end do    
     deallocate( used )                 ! used is not needed anymore
     if( errors > 0 )then
       write(*,*)'+ Error: number of double cells found:',errors,' stop.'
       stop      
     endif
     if( debug > 0 ) write(*,*)'Cells transfered: ',icnt    
   endif

   end subroutine get_original_cells

   !======================================================================
   integer function get_min_cell()
   get_min_cell = min_cell
   end function get_min_cell
   !======================================================================
   integer function get_max_cell()
   get_max_cell = max_cell
   end function get_max_cell

end module store_original_cells
!========================================================================
module store_boundaries
   !
   ! get a boundary, store it in a temporary linked list
   ! when the linked list is filled dump the final
   ! list into an allocated boundary-array
   ! in this way the boundary list may come in any
   ! order and will probably even include 'holes' (default bnds)
   ! in contrast to the vertices the computational boundary-list
   ! may not contain holes; use a map to store the original list
   ! 
   use constants

   private bnd_org

   integer, private :: min_bnd =  99999999
   integer, private :: max_bnd = -99999999

   type :: bnd_org
     integer :: index                    ! index as read
     integer, dimension(4) :: vertices   ! the 4 vertices 
     integer :: rid                      ! region id
     type (bnd_org), pointer :: next     ! next original cell definition
   end type

   type (bnd_org), private, pointer :: head ! pointer to the head of the list
   type (bnd_org), private, pointer :: tail ! pointer to the new orginal cell

   public 

   contains

   subroutine add_boundary(bnd,vertices,rid)
   !
   ! just add a new cell to the list
   !
   integer bnd, rid
   integer, dimension(4) :: vertices

   ! write(*,*)'add_boundary: ',bnd

   if( .not. associated(head) )then
     ! empty list, just add the cell right here
     if( debug > 0 ) write(*,*)'Initialise original boundary list'
     allocate(head, stat=istat)
     tail => head
     nullify( tail%next )
     tail%index    = bnd
     tail%vertices = vertices
     tail%rid      = rid
   else 
     !  write(*,*)'add vertex ',vertex
     allocate(tail%next, stat=istat)
     tail => tail%next
     nullify( tail%next )
     tail%index    = bnd
     tail%vertices = vertices
     tail%rid      = rid
   endif

   min_bnd = min(min_bnd, bnd)
   max_bnd = max(max_bnd, bnd)

   end subroutine add_boundary
   !======================================================================
   subroutine get_original_bnds(n,borg)
   !
   ! transfer the boundaries from the list to the boundary array
   !
   integer n,i,icnt, rid
   integer, dimension(n,5) :: borg  
   integer, dimension(4)   :: temp
   logical, allocatable    :: used(:)
   integer :: errors = 0

   type (bnd_org), pointer :: ptr, tmp

   if( debug > 0 ) write(*,*)'Transfering boundary list into original boundary array'
   if( .not. associated(head) )then
     write(*,*)'get_original_bnds: head of temp. boundary list not assiociated. stop'
     stop
   else
     ptr  => head                       ! point to start of list
     icnt = 0                           ! set counter to zero
     allocate( used(n) )                ! allocate temp. logical used array
     used = .false.                     ! initialise to .false.
     borg = 0                           ! initialise borg array
     do
       if( .not. associated(ptr) )exit  ! test for end of list
       icnt    = icnt + 1               ! bump counter
       i       = ptr%index              ! get boundary index
       if( used(i) )then                ! test if already used
         errors = errors + 1            ! count the errors
         write(*,*)'+ Error: double index found:',i
       endif
       temp    = ptr%vertices           ! copy vertices into temp
       borg(i,1) = temp(1)              ! set corg(i,1), first vertex
       borg(i,2) = temp(2)              ! set corg(i,2), 2nd vertex
       borg(i,3) = temp(3)              ! set corg(i,3)    .
       borg(i,4) = temp(4)              ! set corg(i,4)
       borg(i,5) = ptr%rid              ! region id
       used(i) = .true.                 ! flag this index as used
       tmp     => ptr                   ! store the current pointer
       ptr     => ptr%next              ! get next item in list
       deallocate(tmp)                  ! remove the current item
     end do    
     deallocate( used )                 ! used is not needed anymore
     if( errors > 0 )then
       write(*,*)'+ Error: number of double boundaries found:',errors,'. stop.'
       stop      
     endif
     if( debug > 0 ) write(*,*)'Boundaries transfered: ',icnt    
   endif

   end subroutine get_original_bnds

   !======================================================================
   integer function get_min_bnd()
   get_min_bnd = min_bnd
   end function get_min_bnd
   !======================================================================
   integer function get_max_bnd()
   get_max_bnd = max_bnd
   end function get_max_bnd

end module store_boundaries
module find_faces
   !
   ! loop over cells, store for every vertex the cell
   
   use constants

   type :: vertex_to_cell
     integer :: cell                    
     type (vertex_to_cell), pointer :: next
   end type

   type adresses
     type (vertex_to_cell), pointer :: p
   end type

   type (adresses), allocatable, private:: head(:)
   type (adresses), allocatable, private:: tail(:)

   integer, allocatable :: vc(:)
   
   contains

   !======================================================================
   subroutine set_up_find_faces(n)
   !
   if( debug > 0 ) write(*,*)'Setting up find faces for cells:',n
   allocate(head(n))
   allocate(tail(n))
   allocate(vc(n))

   do i=1,n
     nullify(head(i)%p)
     nullify(tail(i)%p)
   end do
   
   vc = 0
   
   if( debug > 0 ) write(*,*)'ok'
   
   end subroutine set_up_find_faces
   !======================================================================
   subroutine add_cell_to_vertex(vertex,icell)
   !
   ! add a new cell to the list
   !
   integer vertex, icell, last
   
   if( .not. associated(head(vertex)%p) )then
     allocate(head(vertex)%p,stat=istat)
     tail(vertex)%p => head(vertex)%p
     nullify( tail(vertex)%p%next )
     tail(vertex)%p%cell = icell
     vc(vertex) = 1
   else 
     last = tail(vertex)%p%cell
     if( icell /= last )then
       allocate(tail(vertex)%p%next,stat=istat)
       tail(vertex)%p => tail(vertex)%p%next
       nullify( tail(vertex)%p%next )
       tail(vertex)%p%cell  = icell
       vc(vertex) = vc(vertex) + 1
     endif
   endif
   
   end subroutine add_cell_to_vertex
   !======================================================================
   subroutine debug_cell_to_vertex(n)
   integer n          ! number of vertices
   type(adresses) pt

   if( debug <= 3 )then
     write(*,*) 'Partial list:',min(n,10)
     do i=1,min(n,10)
       pt%p => head(i)%p
       write(*,'(/,a,i6,a,i2,a)',advance='no') ' v',i,' (',vc(i),')->'
       do while( associated(pt%p) )
         j = pt%p%cell
         write(*,'(i7)',advance='no') j
         pt%p => pt%p%next
       end do
     end do
     if( n > 30 )then
       write(*,'(3(/,''  .''))')
       do i=(n-9),n
         pt%p => head(i)%p
         write(*,'(/,a,i6,a,i2,a)',advance='no') ' v',i,' (',vc(i),')->'
         do while( associated(pt%p) )
           j = pt%p%cell
           write(*,'(i8)',advance='no') j
           pt%p => pt%p%next
         end do
       end do
     endif
   else
     write(*,*) 'Full list:',n
     do i=1,n
       pt%p => head(i)%p
       write(*,'(/,a,i4,a,i2,a)',advance='no') ' v',i,' (',vc(i),')->'
       do while( associated(pt%p) )
         j = pt%p%cell
         write(*,'(i4)',advance='no') j
         pt%p => pt%p%next
       end do
     end do
   endif

   !write(IOdbg,*) 'Cell to vertex:',n
   !do i=1,n
   !  pt%p => head(i)%p
   !  write(IOdbg,'(/,a,i6,a,i3,a)',advance='no') ' v',i,' (',vc(i),')->'
   !  do while( associated(pt%p) )
   !    j = pt%p%cell
   !    write(IOdbg,'(i8)',advance='no') j
   !    pt%p => pt%p%next
   !  end do
   !end do
   write(*,'(/,'' Dump done'')')
   end subroutine debug_cell_to_vertex
   !=====================================================================
   integer function get_number_of_cells(vertex)
     integer vertex
     
     get_number_of_cells = vc(vertex)
   
   end function get_number_of_cells
   !=====================================================================
   integer function get_cell(vertex,index)
     integer vertex, index
     type(adresses) pt

     if( index > get_number_of_cells(vertex) )then
       write(*,*)'+ Error: Internal error in get_cell'
       get_cell = -1
     endif
     
     pt%p => head(vertex)%p
     
     i = 0
     do while( associated(pt%p) )
       i = i + 1
       if( i == index )then
         j = pt%p%cell
         get_cell = j
         exit
       else
         pt%p => pt%p%next
       endif
     end do
   
   end function get_cell
end module find_faces
!========================================================================
module cell_faces
   !
   ! the basic cell array with pointers to faces
   !
   use constants
   use geometry
   use find_faces

   type TmpFaceData
     integer i1,i2,i3,i4         ! 
     integer cell1,side1         ! first cell
     integer cell2,side2         ! second cell
     real    area                ! area
     real    x(3)                ! center
     real    n(3)                ! normal/cell face area components
   end type

   type(TmpFaceData), allocatable :: TmpFace(:) ! too large: 6*Ncel !
   integer :: TmpFaceAll = 0
   integer :: TmpFaceCnt = 0

   integer, allocatable :: TmpFaceCell(:)          
   integer, allocatable :: TmpFaceSide(:)          
   logical, allocatable :: TmpFaceUsed(:)
   
   integer :: FaceCnt   = 0

   integer :: NAllFaces = 0          ! total number of faces


   contains

   !======================================================================
   subroutine set_up_cell_faces(Ncel,Nvrt)
   !
   if( debug > 0 ) write(*,*)'Setting up list faces for cells:',Ncel,Nvrt


   allocate(TmpFace(6*Ncel))    ! temporary face list
   allocate(TmpFaceCell(6*Ncel))          
  !allocate(TmpFaceSide(6*Ncel))          
   allocate(TmpFaceUsed(6*Ncel))          

   TmpFaceCell = 0
  !TmpFaceSide = 0
   TmpFaceUsed = .false.
   
   if( debug > 0 ) write(*,*)'ok'

   end subroutine set_up_cell_faces
   !=====================================================================
   subroutine store_tmp_face(icell,side,i1,i2,i3,i4)
   !
   ! a face get's in and just store it
   !   
   integer icell, side, i1, i2, i3, i4

   TmpFaceCnt = TmpFaceCnt + 1
   TmpFaceAll = TmpFaceCnt
   
   TmpFaceCell(TmpFaceCnt) = icell

   TmpFace(TmpFaceCnt)%i1    = i1
   TmpFace(TmpFaceCnt)%i2    = i2
   TmpFace(TmpFaceCnt)%i3    = i3
   TmpFace(TmpFaceCnt)%i4    = i4
   TmpFace(TmpFaceCnt)%cell1 = icell
   TmpFace(TmpFaceCnt)%side1 = side
   TmpFace(TmpFaceCnt)%cell2 = 0
   TmpFace(TmpFaceCnt)%side2 = 0

   end subroutine store_tmp_face
   !=====================================================================
   subroutine monitor_face
   
     write(*,*)'Debug face'
     do i=1,TmpFaceAll
       if( (TmpFace(i)%cell1 > 8 .and. TmpFace(i)%cell1 < 13) .or. &
           (TmpFace(i)%cell2 > 8 .and. TmpFace(i)%cell2 < 13)  ) & 
       write(*,'('' f'',i4,'' c:'',2(1x,i4),''    v:'',4(1x,i4),1x,l1)') &
           i,TmpFace(i)%cell1,TmpFace(i)%cell2, &
           TmpFace(i)%i1,TmpFace(i)%i2,TmpFace(i)%i3,TmpFace(i)%i4,&
           tmpfaceused(i)
     end do
   
   end subroutine monitor_face 
   subroutine process_tmp_faces(CellMap_CO)
   !
   ! 
   !   
   integer CellMap_CO(:)
   integer icell, face, side
   integer cell1, cell2, face1, face2, side1, side2, i1, i2, i3
   logical, allocatable :: mask(:)
   logical :: hit

   type(VFadresses) ptr
   integer nodes(4)

   !
   ! the vertex -> face pointers...
   !
   call set_up_add_face_to_vertex(Nvrt)

   do i=1,TmpFaceAll

     i1 = TmpFace(i)%i1
     i2 = TmpFace(i)%i2
     i3 = TmpFace(i)%i3
     i4 = TmpFace(i)%i4

     !write(*,'(A,5(1x,i6))') ' f:',i,i1,i2,i3,i4
     call add_face_to_vertex(i1,i)
     call add_face_to_vertex(i2,i)
     call add_face_to_vertex(i3,i)
     if( i4 > 0 ) call add_face_to_vertex(i4,i)

   end do
   
   do i=1,TmpFaceAll
     if( mod(i,250000) == 0 ) write(*,*) ' ',i
     if( TmpFaceUsed(i) ) cycle
   
     ! pak de (nu nog 4) knooppunten
     !
     i1 = TmpFace(i)%i1
     i2 = TmpFace(i)%i2
     i3 = TmpFace(i)%i3
     i4 = TmpFace(i)%i4
   
     ichk = i1 + i2 + i3 + i4
     
     nodes(1) = i1     
     nodes(2) = i2     
     nodes(3) = i3     
     nodes(4) = i4     

     nodeloop: do inode=1,4
       ix = nodes(inode)
       if( ix < 0 ) cycle 
       !
       ! pak de faces van deze vertex
       !
       ptr%ptr => VFhead(ix)%ptr
       whileloop: do while( associated(ptr%ptr) )
         j  = ptr%ptr%face
         if( j == i )then
           ptr%ptr => ptr%ptr%next
           cycle
         endif
         j1 = TmpFace(j)%i1
         j2 = TmpFace(j)%i2
         j3 = TmpFace(j)%i3
         j4 = TmpFace(j)%i4
         
         jchk = j1 + j2 + j3 + j4
         !
         ! als ichk == jchk dan pas verder kijken
         !
         if( ichk == jchk )then
           hit = .false. 
           
           if( i1 == j1 .or. i1 == j2 .or. i1 == j3 .or. i1 == j4 )then
             if( i2 == j1 .or. i2 == j2 .or. i2 == j3 .or. i2 == j4 )then
               if( i3 == j1 .or. i3 == j2 .or. i3 == j3 .or. i3 == j4 )then
                 if( i4 == j1 .or. i4 == j2 .or. i4 == j3 .or. i4 == j4 )then
                   hit = .true.
                 endif
               endif
             endif
           endif
           
           if( hit )then
           
             TmpFaceUsed(i) = .true.
             TmpFaceUsed(j) = .true.

             TmpFace(i)%cell2 = TmpFace(j)%cell1
             TmpFace(i)%side2 = TmpFace(j)%side1 

             ! 'delete' 2nd reference
             TmpFaceCell(j)   = 0
             TmpFace(j)%cell1 = 0
             TmpFace(j)%side1 = 0
           
             exit nodeloop         
           endif
         endif
         ptr%ptr => ptr%ptr%next
       end do whileloop
     end do nodeloop
   
     !
     ! dus dit is een single face
     !
     TmpFaceUsed(i) = .true.
     !
     !if( i < 20 ) write(*,*)'f:',i,j,ichk,jchk,hit
   end do
   
   n1 = 0
   n2 = 0
   do i=1,TmpFaceAll
     if( TmpFaceUsed(i) )then
       i1 = TmpFace(i)%cell1
       i2 = TmpFace(i)%cell2
        !if( i1 > 0 .and. i < 10 )write(*,'(7(1x,i4))') i,  &
        ! CellMap_CO(i1),CellMap_CO(i2) &
        !,tmpface(i)%i1,tmpface(i)%i2,tmpface(i)%i3,tmpface(i)%i4
       if( i1 > 0 .and. i2 ==0 ) n1 = n1 + 1
       if( i1 > 0 .and. i2 > 0 ) n2 = n2 + 1
     endif
   end do   
   write(*,*)'Single faces found:           ',n1 
   write(*,*)'Internal faces found:         ',n2 

   if( debug > 0 ) write(*,*)'Deallocating temp. arrays'

   deallocate(TmpFaceCell)          
  !deallocate(TmpFaceSide)          
   deallocate(TmpFaceUsed)          

   end subroutine process_tmp_faces
end module cell_faces
!========================================================================
!========================================================================
!========================================================================
!========================================================================
integer function lens(string)

   character(len=*) string
   
   do i=len(string),0,-1
     if( string(i:i) .ne. ' ')goto 10
   end do
   i = 0
10 continue

   lens = i
    
end function lens
!========================================================================
subroutine openfile(iunit,casename,extension,reqform,status,idebug)

  character(len=*)  casename
  character(len=*)  extension
  character(len=*)  reqform
  character(len=*)  status
  character(len=48) filename
  character(len=11) form 

  logical exists


  filename = casename(1:lens(casename))//extension(1:lens(extension))
  length   = lens(filename)

  if( idebug > 2 )write(*,*) 'Opening ',filename(1:length)
  
  if( status(1:3) == 'OLD' )then
    inquire(file=filename(1:length),exist=exists,form=form)
    if( .not. exists )then
      write(*,*) '*** Error: File ',filename(1:length),' does not exist'
      stop
    endif
  endif

  open(iunit,file=filename(1:length),form=reqform,status=status)

  if( idebug >= 2 ) write(*,*) 'File ',filename(1:length),' opened'
  
end subroutine openfile
subroutine readgeometry(casename)
   use constants
   use store_vertices
   use store_original_cells
   use store_boundaries
   use find_faces
   use cell_faces
   use geometry

   real, dimension(3) :: dummy
   real    xmin, xmax, ymin, ymax, zmin, zmax
   integer tmin, tmax, rmin, rmax
   
   integer, allocatable :: corg(:,:), borg(:,:)
   integer, allocatable :: CellMap_CO(:), CellMap_OC(:)
   integer, allocatable :: cellshape(:)
   logical, allocatable :: mask(:)
   logical  skip, hack
   
   integer idummy1, idummy2
   integer, dimension(8) :: idummy
   logical  exists
   character cpfa*4

   character(len=32) casename

   real, dimension(3) :: x1, x2, x3, x4, x31, x42, &
                         s, a, b, c, normal, center

   !
   ! read files for now assume they are called zzz.xxx
   !
   write(*,*) 'Opening vertex file'
   call openfile(IOvrt,casename,'.vrt','FORMATTED','OLD',debug)

   icnt = 0
   do 
     !read(IOvrt,*,end=1) idummy1, (dummy(i),i=1,3)
     read(IOvrt,'(i9,6x,3g16.9)',end=1) idummy1, (dummy(i),i=1,3)
     
     call add_vertex(idummy1,dummy)
     icnt = icnt + 1
     !write(*,*) 'v :',icnt,(dummy(i),i=1,3)
   end do
 1 continue  
   Nvrt = get_max_vertex()   
   
   write(*,*) 'Number of vertices read:   ',icnt
   if( icnt /= Nvrt )then
     write(*,*) '+ Warning vertex list has holes'
     write(*,*) '+ Counted, found:',icnt,Nvrt
   endif
   if( get_min_vertex() <= 0 )then
     write(*,*) '+ Error negative vertex index found. stop.'
     stop
   endif
   if( get_min_vertex() > 1 )then
     write(*,*) '+ Warning minimum vertex index /= 1 : ',get_min_vertex()
     stop
   endif
   if( debug > 0 )then
     write(*,*) '  Min. vertex found: ',get_min_vertex()
     write(*,*) '  Max. vertex found: ',get_max_vertex()   
   endif

   close(IOvrt)
   !
   ! process vertices
   !
   allocate(Vert(Nvrt,3),stat=istat)   

   call get_vertices(Nvrt,Vert)

   xmin = minval(Vert(:,1))
   xmax = maxval(Vert(:,1))
   ymin = minval(Vert(:,2))
   ymax = maxval(Vert(:,2))
   zmin = minval(Vert(:,3))
   zmax = maxval(Vert(:,3))

   write(*,*) xmin,' <= x <=',xmax
   write(*,*) ymin,' <= y <=',ymax
   write(*,*) zmin,' <= z <=',zmax
   !
   ! === end vertex section, start cells ===
   !
   !open(IOcel,file="zzz.cel")
   call openfile(IOcel,casename,'.cel','FORMATTED','OLD',debug)
   icnt = 0
   do 
     read(IOcel,*,end=2) idummy1,(idummy(i),i=1,8),idummy2
     call add_original_cell(idummy1,idummy,idummy2)
     icnt = icnt + 1
   end do
 2 continue  
   Norgcel = get_max_cell()
   write(*,*) 'Number of cells read:      ',Norgcel
   if( icnt /= Norgcel )then
     write(*,*) '+ Warning cell list has holes'
     write(*,*) '+ Counted, found:',icnt,Norgcel
   endif
   if( get_min_cell() <= 0 )then
     write(*,*) '+ Error negative cell index found. stop.'
     stop
   endif
   if( get_min_cell() > 1 )then
     write(*,*) '+ Warning minimum cell index /= 1 : ',get_min_cell()
   endif
   if( debug > 0 )then
     write(*,*) '  Min cell found: ',get_min_cell()
     write(*,*) '  Max cell found: ',get_max_cell()   
   endif
   close(IOcel)   

   !
   ! process original cells
   !
   allocate(corg(Norgcel,9),stat=istat)            ! 8 vertices + ctid

   call get_original_cells(Norgcel,corg) 

   tmin = minval(corg(:,9),corg(:,9) > 0 )
   tmax = maxval(corg(:,9))

   if( tmin == tmax )then
     write(*,*) 'One cell type found:',tmin
   else
     write(*,*) '  Min. cell type defined:',tmin
     write(*,*) '  Max. cell type defined:',tmax
   endif

   !do i=1,Norgcel
   !  write(IOdbg,'(A,i6,8(1x,i6))') ' co: ',i,corg(i,1:8)
   !end do
   !
   ! === end cell section, start boundaries ===
   !
   !open(IObnd,file="zzz.bnd")
   call openfile(IObnd,casename,'.bnd','FORMATTED','OLD',debug)
   icnt = 0
   do 
     read(IObnd,*,end=3) idummy1,(idummy(i),i=1,4),idummy2
     call add_boundary(idummy1,idummy,idummy2)
     icnt = icnt + 1
   end do
 3 continue  
   Norgbnd = get_max_bnd()
   write(*,*) 'Number of boundaries read: ',Norgbnd
   if( icnt /= Norgbnd )then
     write(*,*) '+ Warning boundary list has holes'
     write(*,*) '+ Counted, found:',icnt,Norgbnd
   endif
   close(IObnd)

   allocate(borg(Norgbnd,5),stat=istat)          ! 4 vertices + rid
   call get_original_bnds(Norgbnd,borg) 

   rmin = minval(borg(:,5), borg(:,5) > 0 )
   rmax = maxval(borg(:,5))

   write(*,*) '  Min. region defined:',rmin
   write(*,*) '  Max. region defined:',rmax
   !====================================================================
   !
   !  find the faces
   !
   !  first walk through the cells establish the type
   !    6 hexaeder:  1,2,3,4 1,2,3,4  
   !    5 pentaeder: 1,2,3,3 5,6,7,7  
   !    5 pyramid:   1,2,3,4 5,5,5,5  
   !    4 tetraeder: 1,2,3,3 4,4,4,4  
   !    1 quad:      1,2,3,4          
   !    1 triangle:  1,2,3,3          
   !
   !  dump the geometry in a file
   !
   !  x,y,z vertices
   !  cells types with faces numbers
   !  faces with pointers to 1 or 2 cells,  boundary flag
   !  
   !         4 E        3 F 
   !          +----------+ 
   !         /    (bot) /|     face 1: 1 2 3 4   bottom
   !        /   F4(f1) / |          2: 5 6 7 8   top
   !   8 G /   north 7/H |          3: 1 2 6 5   south
   !      +----------+ F6|east      4: 4 3 8 7   north
   !      |  1+A     |   +2 B       5: 1 5 8 4   west
   !  f5  |    F2    |  /           6: 2 3 7 6   east
   !  west|   top    | /  
   !      |          |/   
   !      +----------+
   !     5 C   f3   6 D  
   !         south
   !
   call set_up_find_faces(Nvrt)

   icnt = 0
   do i=1,Norgcel
     if( corg(i,1) /= 0 .and. corg(i,5) /= 0 )then
       !write(*,'(A,9(1x,i8))')'+ c,v:',i,(corg(i,j),j=1,8)
       do j=1,8
         if( corg(i,j) <= Nvrt )then
           call add_cell_to_vertex(corg(i,j),i)     
         else
           write(*,*)'+ Error: Vertex out of range in cell:',i,corg(i,j)
           icnt = icnt + 1
         endif
       end do
     endif
   end do
   if( icnt > 0 )then
     write(*,*)'+ Error: Vertices out of range. stop.'
     stop
   endif
   if( debug > 1 )then
     write(*,*)'Dumping list...'
!     call debug_cell_to_vertex(Norgcel)
     call debug_cell_to_vertex(Nvrt)
   endif
   !
   ! cell map needs to be compressed without holes
   ! first list the active fluid or solid cells 
   ! followed by the baffles/walls/boundaries
   !
   write(*,*)'Set up maps...',Norgcel
   allocate(CellMap_CO(Norgcel),stat=istat)   
   allocate(CellMap_OC(Norgcel),stat=istat)   
   
   icnt = 0
   do i=1,Norgcel
     if( corg(i,1) /= 0 )then           ! skip empty cells
       !if( corg(i,5) /= 0 )then         ! skip baffles/shells
         !if( corg(i,9) /= 6)then        ! skip ghosts
           icnt = icnt + 1
           CellMap_CO(icnt) = i           ! 
           CellMap_OC(i)    = icnt        ! 
         !endif
       !endif
     endif
   end do
   Ncel = icnt
   write(*,*)'Mapped active cells to original cell array '
   write(*,*)'Number of active fluid/solid cells:',Ncel
   !
   ! get shape of cell 
   !
   allocate(cellshape(Ncel),stat=istat)   

   write(*,*)'Set up cell faces...'
   call set_up_cell_faces(Ncel,Nvrt)
   
   ihex   = 0
   iprism = 0
   ipyrm  = 0
   itet   = 0
   iquad  = 0
   itri   = 0
   do i=1,Ncel
     j  = CellMap_CO(i)
     i1 = corg(j,1)
     i2 = corg(j,2)
     i3 = corg(j,3)
     i4 = corg(j,4)
     i5 = corg(j,5)
     i6 = corg(j,6)
     i7 = corg(j,7)
     i8 = corg(j,8)
     
     i9 = corg(j,9)
     
     !if( i9 == 6 ) CYCLE    ! HACK   SOLID   TEST   

     if( i7 /= 0 )then
       if( i7 /= i8 )then
         ihex   = ihex + 1
         cellshape(i) = Hexa
       else if( i5 /= i6 )then
         iprism = iprism + 1
         cellshape(i) = Prism
       else if( i3 /= i4 )then
         ipyrm  = ipyrm + 1
         cellshape(i) = Pyramid
       else if( i3 == i4 .and. i5 == i6 )then
         itet   = itet + 1
         cellshape(i) = Tetra
       else
         write(*,*)'Unknown shape ',i,' : ',i1,i2,i3,' ...'
       endif

       if( cellshape(i) == Hexa )then
         call store_tmp_face(i,1,i1,i2,i3,i4)   !face 1
         call store_tmp_face(i,2,i5,i6,i7,i8)   !face 2
         call store_tmp_face(i,3,i1,i2,i6,i5)   !face 3
         call store_tmp_face(i,4,i4,i3,i7,i8)
         call store_tmp_face(i,5,i1,i5,i8,i4)
         call store_tmp_face(i,6,i2,i6,i7,i3)   !face 6
       else if( cellshape(i) == Prism )then
         call store_tmp_face(i,1,i1,i2,i3,-1)
         call store_tmp_face(i,2,i5,i6,i7,-1)
         call store_tmp_face(i,3,i1,i2,i6,i5)
         ! collapsed
         call store_tmp_face(i,5,i1,i5,i8,i4)
         call store_tmp_face(i,6,i2,i6,i7,i3)
       else if( cellshape(i) == Pyramid )then
         
         !write(*,*)'*** PYRAMIDE OK? ***'

         call store_tmp_face(i,1,i1,i2,i3,i4)   !face 1
         ! collapsed
         call store_tmp_face(i,3,i1,i2,i5,-1)   !face 3
         call store_tmp_face(i,4,i4,i3,i5,-1)
         call store_tmp_face(i,5,i1,i5,i4,-1)
         call store_tmp_face(i,6,i2,i3,i5,-1)   !face 6

       else if( cellshape(i) == Tetra )then
         call store_tmp_face(i,1,i1,i2,i3,-1)   !face 1
         ! collapsed
         call store_tmp_face(i,3,i1,i2,i5,-1)   !face 3
         ! collapsed
         call store_tmp_face(i,5,i1,i5,i3,-1)
         call store_tmp_face(i,6,i2,i3,i5,-1)   !face 6
       else
         write(*,*)'+ Internal error: cellshape: ',i,j
       endif
     else
       !
       ! these 'cells' should be on any fluid/solid
       ! do not store them they should be found below
       !
       if( i3 /= i4 )then
         iquad   = iquad + 1
         cellshape(i) = Quad
       else
         itri    = itri + 1
         cellshape(i) = Triangle
       endif
     endif
   end do
   if( ihex   > 0 )write(*,*)'Number of hexa''s     :',ihex
   if( iprism > 0 )write(*,*)'Number of prism''s    :',iprism
   if( ipyrm  > 0 )write(*,*)'Number of pyramids''s :',ipyrm
   if( itet   > 0 )write(*,*)'Number of tetra''s    :',itet
   if( iquad  > 0 )write(*,*)'Number of quad''s     :',iquad
   if( itri   > 0 )write(*,*)'Number of triangles''s:',itri
   !call monitor_face                          

   !
   ! show cell table
   !
   write(*,*)'Cell table from ',tmin,' to ',tmax
   
   do it=tmin,tmax
     itmp = count( corg(:,9) == it )
     write(*,*)'Cell table id ',it,itmp
   end do     

   !
   ! here the cell table should be read in
   !
   ! cells:   fluid, solid,  ghost
   ! shells:  baffle, plate, shell
   !

   !
   ! now work on the faces of each active fluid/solid cel
   ! objective is to sort them and assign one or two cells
   ! to them
   !
   call process_tmp_faces(CellMap_CO)
   !
   ! all faces known set up faces
   !
   if( debug > 0 ) write(*,*)'Total number of temporary faces:',TmpFaceAll
   

   f250 = 1./4.
   f333 = 1./3.
   do i=1,TmpFaceAll
     
     i1 = TmpFace(i)%cell1
     if( i1 <= 0 ) cycle
     
     !write(*,*) 'i:',i,Tmpface(i)%i1,Tmpface(i)%i2,Tmpface(i)%i3,Tmpface(i)%i4
     
     x1(1:3) = Vert(Tmpface(i)%i1,1:3) 
     x2(1:3) = Vert(Tmpface(i)%i2,1:3) 
     x3(1:3) = Vert(Tmpface(i)%i3,1:3) 
     
     if( TmpFace(i)%i4 /= -1 )then
       x4(1:3) = Vert(Tmpface(i)%i4,1:3) 
     else
       x4(1:3) = 0.0
     endif

     !                             _   _       _    _
     ! surface is defined by: 1/2( x - x ) X ( x  - x  )
     !                              3   1       4    2
     if( TmpFace(i)%i4 /= -1 )then
       a = x3 - x1    
       b = x4 - x2
       !
       ! the surface normal
       !
       normal(1) = 0.5*( a(2)*b(3) - a(3)*b(2) )
       normal(2) = 0.5*( a(3)*b(1) - a(1)*b(3) )
       normal(3) = 0.5*( a(1)*b(2) - a(2)*b(1) )
     else
       a = x3 - x1    
       b = x2 - x1
       !
       ! the surface normal
       !
       normal(1) = 0.5*( a(2)*b(3) - a(3)*b(2) )
       normal(2) = 0.5*( a(3)*b(1) - a(1)*b(3) )  !<= OK ????
       normal(3) = 0.5*( a(1)*b(2) - a(2)*b(1) )
     endif

     area = sqrt( normal(1)*normal(1) &
                + normal(2)*normal(2) + normal(3)*normal(3) )

     TmpFace(i)%area = area

  !   if( i <= 4 )then
  !     write(*,*)'A  ',i
  !     write(*,*)'  i   ',Tmpface(i)%i1,Tmpface(i)%i2,Tmpface(i)%i3,Tmpface(i)%i4
  !     write(*,*)'  a   ',area,normal
  !     write(*,*)'  x1  ',x1
  !     write(*,*)'  x2  ',x2
  !     write(*,*)'  x3  ',x3
  !     write(*,*)'  x4  ',x4
  !     write(*,*)'x3-x1 ',x3-x1
  !     write(*,*)'x4-x2 ',x4-x2
  !     write(*,*)' a ,b ',a,b
  !   endif

     TmpFace(i)%n    = normal   ! later check if direction is ok (p -> face)
     !
     ! cell face centre
     !
     s = 0
     if( TmpFace(i)%i4 > 0 )then
       s(1) = f250 * ( Vert(Tmpface(i)%i1,1) + Vert(Tmpface(i)%i2,1) &
                     + Vert(Tmpface(i)%i3,1) + Vert(Tmpface(i)%i4,1) ) 
       s(2) = f250 * ( Vert(Tmpface(i)%i1,2) + Vert(Tmpface(i)%i2,2) &
                     + Vert(Tmpface(i)%i3,2) + Vert(Tmpface(i)%i4,2) )
       s(3) = f250 * ( Vert(Tmpface(i)%i1,3) + Vert(Tmpface(i)%i2,3) &
                     + Vert(Tmpface(i)%i3,3) + Vert(Tmpface(i)%i4,3) )         
     else if( TmpFace(i)%i4 == -1 )then
       s(1) = f333 * ( Vert(Tmpface(i)%i1,1) + Vert(Tmpface(i)%i2,1) &
                     + Vert(Tmpface(i)%i3,1) ) 
       s(2) = f333 * ( Vert(Tmpface(i)%i1,2) + Vert(Tmpface(i)%i2,2) &
                     + Vert(Tmpface(i)%i3,2) )
       s(3) = f333 * ( Vert(Tmpface(i)%i1,3) + Vert(Tmpface(i)%i2,3) &
                     + Vert(Tmpface(i)%i3,3) )         
     else
       write(*,*)'+ Error: tmp. cell face data structure corrupted'
     endif

     TmpFace(i)%x = s

   end do

   !
   ! geometry caluculations
   !
   write(*,*)'Starting geometry calculations...'
   allocate( Cell(1:Ncel) )
   
   totvol = 0.0

   f125 = 1./8.
   f166 = 1./6.
   f200 = 1./5.
   f250 = 1./4.
   
   do i=1,Ncel
     j  = CellMap_CO(i)
     i1 = corg(j,1)
     i2 = corg(j,2)
     i3 = corg(j,3)
     i4 = corg(j,4)
     i5 = corg(j,5)
     i6 = corg(j,6)
     i7 = corg(j,7)
     i8 = corg(j,8)
     
     x = 0.0
     y = 0.0
     z = 0.0
     volume = 0.0
     skip = .false.
     
     if( CellShape(i) == Hexa )then
     
       x = f125 * ( Vert(i1,1) + Vert(i2,1) + Vert(i3,1) + Vert(i4,1) &
                  + Vert(i5,1) + Vert(i6,1) + Vert(i7,1) + Vert(i8,1) )
       y = f125 * ( Vert(i1,2) + Vert(i2,2) + Vert(i3,2) + Vert(i4,2) &
                  + Vert(i5,2) + Vert(i6,2) + Vert(i7,2) + Vert(i8,2) )
       z = f125 * ( Vert(i1,3) + Vert(i2,3) + Vert(i3,3) + Vert(i4,3) &
                  + Vert(i5,3) + Vert(i6,3) + Vert(i7,3) + Vert(i8,3) )       
                  
       call  calc_vol(volume,i1,i2,i3,i4,i5,i6,i7,i8)
       
     else if( CellShape(i) == Prism )then
     
       x = f166 * ( Vert(i1,1) + Vert(i2,1) + Vert(i3,1) + Vert(i5,1) &
                  + Vert(i6,1) + Vert(i7,1) )
       y = f166 * ( Vert(i1,2) + Vert(i2,2) + Vert(i3,2) + Vert(i5,2) &
                  + Vert(i6,2) + Vert(i7,2) )
       z = f166 * ( Vert(i1,3) + Vert(i2,3) + Vert(i3,3) + Vert(i5,3) &
                  + Vert(i6,3) + Vert(i7,3) )        
     
       !write(*,*)'>>',x,y,z
       
       call  calc_vol(volume,i1,i2,i3,i3,i5,i6,i7,i7)

     else if( CellShape(i) == Pyramid )then
     
       x = f200 * ( Vert(i1,1) + Vert(i2,1) + Vert(i3,1) + Vert(i4,1) &
                  + Vert(i5,1) )
       y = f200 * ( Vert(i1,2) + Vert(i2,2) + Vert(i3,2) + Vert(i4,2) &
                  + Vert(i5,2) )
       z = f200 * ( Vert(i1,3) + Vert(i2,3) + Vert(i3,3) + Vert(i4,3) &
                  + Vert(i5,3) )        
     
       call  calc_vol(volume,i1,i2,i3,i4,i5,i5,i5,i5)

     else if( CellShape(i) == Tetra )then
     
       x = f250 * ( Vert(i1,1) + Vert(i2,1) + Vert(i3,1) + Vert(i5,1) )
       y = f250 * ( Vert(i1,2) + Vert(i2,2) + Vert(i3,2) + Vert(i5,2) )
       z = f250 * ( Vert(i1,3) + Vert(i2,3) + Vert(i3,3) + Vert(i5,3) )        
     
       call  calc_vol(volume,i1,i2,i3,i3,i5,i5,i5,i5)

     else if( CellShape(i) == Quad .or. CellShape(i) ==Triangle )then

       ! ignore
       skip = .true.

     else

       write(*,*)'+ Info: Unexpectedly skipped cell'     
       cycle

     endif

     !if( corg(j,9) == 6 )then  
     !                            ! ***************************
     !  skip = .true.             ! HACK test for solids/ghosts
     !                            ! ***************************
     !endif

     if( .not. skip )then
       Cell(i)%vertices(1) = i1
       Cell(i)%vertices(2) = i2
       Cell(i)%vertices(3) = i3
       Cell(i)%vertices(4) = i4
       Cell(i)%vertices(5) = i5
       Cell(i)%vertices(6) = i6
       Cell(i)%vertices(7) = i7
       Cell(i)%vertices(8) = i8
       Cell(i)%x(1) = x
       Cell(i)%x(2) = y
       Cell(i)%x(3) = z
       Cell(i)%ctid = corg(j,9)
       Cell(i)%vol  = volume
       Cell(i)%corg = j

       totvol = totvol + volume
     endif

   end do
   write(*,*)'Total volume: ',totvol     

   !
   ! the normal of the cell faces should lie in the direction 
   ! from point p (first cell of face) to the face center (or neighbor n)
   !
   pi = 4.0*atan(1.0)
   iflip = 0
   do i=1,TmpFaceAll

     ip = TmpFace(i)%cell1
     if( ip > 0 )then

       a = TmpFace(i)%x - Cell(ip)%x
       b = TmpFace(i)%n
       
       !write(*,*)'i :',i
       !write(*,*)' a:',a,b

       call normalise( a )
       call normalise( b )

       
       dotp   = dot_product( a , b )  !<= test op 0!!
       dotp   = min( dotp ,  1.0 )
       dotp   = max( dotp , -1.0 )
       angle  = acos( dotp )*180./pi
       if( angle > 90.0 .and. angle < 270.0 )then
         ! flip normal
         if( i <= 12 )write(*,*)'flip ',i,angle
         iflip = iflip + 1
         TmpFace(i)%n = -TmpFace(i)%n

       !else
       endif
     endif
   end do
   write(*,*)'Normals checked'
   if( iflip > 0 )write(*,*)'Normal flipped:',iflip
   !
   ! read couples (nog niet algemeen!)
   !
   inquire(file=casename(1:lens(casename))//'.cpl',exist=exists)
   !exists = .false.
   if( exists )then
     if( debug > -1 )then
       icnt = count( TmpFace(:)%cell1 > 0 )
       write(*,*)'Total number of active faces before couples:',icnt
     endif
     icnt = 0
     !open(11,file='zzz.cpl')
     call openfile(11,casename,'.cpl','FORMATTED','OLD',debug)
     
     nivs1 = 0
     nivs2 = 0
     nivs3 = 0
     do
       islave1 = -1
       islave2 = -1
       islave3 = -1

        read(11,'(A4,1x,i8,2(i8,i2))',end=4,eor=13,advance='no') &
                     cpfa,icpl,imaster,imasterside, islave1,iside1
    
        read(11,'(i8,i2)',end=4,eor=13,advance='no') islave2, iside2
        !let op eor=13 weggehaald voor gfortran
        read(11,'(i8,i2)',end=4,advance='yes') islave3, iside3
    
       !read(11,'(i8,i6,i9)',end=4,eor=13,advance='yes') &
       !             icpl,itmpcpl,itypecpl
    
       !read(11,'(i9,i2)',end=4,eor=13,advance='no')  imaster,imasterside
       !read(11,'(i9,i2)',end=4,eor=13,advance='no')  islave1,iside1
       !read(11,'(i9,i2)',end=4,eor=13,advance='no')  islave2, iside2
       !read(11,'(i9,i2)',end=4,eor=13,advance='yes') islave3, iside3

    13 continue   
       if( islave1 > 0 ) ivs = 1
       if( islave2 > 0 ) ivs = 2
       if( islave3 > 0 ) ivs = 3

       if( ivs == 1 ) nivs1 = nivs1 + 1
       if( ivs == 2 ) nivs2 = nivs2 + 1
       if( ivs == 3 ) nivs3 = nivs3 + 1

       if( ivs >= 2 )then
         !
         ! one master, two slaves
         !
         im = CellMap_OC(imaster)
         is = CellMap_OC(imasterside)
         do i=1,TmpFaceAll
           if( TmpFace(i)%cell1 == im )then
             if( TmpFace(i)%side1 == is )then
                 ifaceM = i
             endif
           endif
         end do

         im = CellMap_OC(islave1)
         is = CellMap_OC(iside1)
         do i=1,TmpFaceAll
           if( TmpFace(i)%cell1 == im )then
             if( TmpFace(i)%side1 == is )then
                 ifaceS1 = i
             endif
           endif
         end do

         im = CellMap_OC(islave2)
         is = CellMap_OC(iside2)
         do i=1,TmpFaceAll
           if( TmpFace(i)%cell1 == im )then
             if( TmpFace(i)%side1 == is )then
                 ifaceS2 = i
             endif
           endif
         end do

         if( ivs > 2 )then
           im = CellMap_OC(islave3)
           is = CellMap_OC(iside3)
           do i=1,TmpFaceAll
             if( TmpFace(i)%cell1 == im )then
               if( TmpFace(i)%side1 == is )then
                   ifaceS3 = i
               endif
             endif
           end do
         else
           ifaceS3 = -1
         endif

         icnt = icnt + 1
         write(*,*) 'Couple ',icnt,ifaceM,'=>',ifaceS1,ifaceS2,ifaceS3
         write(*,*) '  Master:',TmpFace(ifacem)%cell1,TmpFace(ifaceM)%side1
         write(*,*) '  Slave1:',TmpFace(ifaces1)%cell1,TmpFace(ifaces1)%side1
         write(*,*) '  Slave2:',TmpFace(ifaces2)%cell1,TmpFace(ifaces2)%side1
         if( ivs > 2 )then
           write(*,*) '  Slave3:',TmpFace(ifaces3)%cell1,TmpFace(ifaces3)%side1
         endif

         TmpFace(ifaceS1)%cell2 = TmpFace(ifaceM)%cell1
         TmpFace(ifaceS1)%side2 = TmpFace(ifaceM)%side1

         TmpFace(ifaceS2)%cell2 = TmpFace(ifaceM)%cell1
         TmpFace(ifaceS2)%side2 = TmpFace(ifaceM)%side1

         if( ivs > 2 )then
           TmpFace(ifaceS3)%cell2 = TmpFace(ifaceM)%cell1
           TmpFace(ifaceS3)%side2 = TmpFace(ifaceM)%side1
         endif

         TmpFace(ifaceM)%cell1  = 0 ! kill the lonely master...
         TmpFace(ifaceM)%side1  = 0
       else if( ivs == 1 )then
         !
         ! one master, one slave  NIET OK!
         !
         im = CellMap_OC(imaster)
         is = CellMap_OC(imasterside)
         do i=1,TmpFaceAll
           if( TmpFace(i)%cell1 == im )then
             if( TmpFace(i)%side1 == is )then
                 ifaceM = i
             endif
           endif
         end do

         im = CellMap_OC(islave1)
         is = CellMap_OC(iside1)
         do i=1,TmpFaceAll
           if( TmpFace(i)%cell1 == im )then
             if( TmpFace(i)%side1 == is )then
                 ifaceS1 = i
             endif
           endif
         end do

         icnt = icnt + 1
         write(*,*) 'Couple ',ifaceM,'=>',ifaceS1
         write(*,*) '  Master:',TmpFace(ifacem)%cell1,TmpFace(ifaceM)%side1
         write(*,*) '  Slave1:',TmpFace(ifaces1)%cell1,TmpFace(ifaces1)%side1

         TmpFace(ifaceS1)%cell2 = TmpFace(ifaceM)%cell1
         TmpFace(ifaceS1)%side2 = TmpFace(ifaceM)%side1

         TmpFace(ifaceM)%cell1  = 0 ! kill the lonely master...
         TmpFace(ifaceM)%side1  = 0
       else
         write(*,*)'+ Internal error, ivs=',ivs
       endif

     end do
 4   continue
     write(*,*) 'Number of couples found: ',icnt
     write(*,*) '  M / S     : ',nivs1
     write(*,*) '  M / S+S   : ',nivs2
     write(*,*) '  M / S+S+S : ',nivs3

   endif

   !
   ! baffle, plate, shell section
   !
   ! loop over the original cells
   ! where a quad/triangle matches as a solid wall
   ! just cut the relationship between the cells
   ! by removing the face pointers
   !
   hack = .false.
   
   if( hack )then
     write(*,*)'Checking for shells'

     icnt   = 0
     iloose = 0

     do ic=1,Norgcel
                                 ! ***************************
       i1 = corg(ic,1)           ! HACK test for solids/ghosts
       i2 = corg(ic,2)           ! ***************************
       i3 = corg(ic,3)
       i4 = corg(ic,4)
       i5 = corg(ic,5)

       if( i5 /= 0 ) cycle       ! looking for shells only

       i9 = corg(ic,9)            ! cell table id

       if( i9 >= 7 )then
         ichk  = i1 + i2 + i3 +i4
         ihit  = 0
         ifnd  = 0
         ifnd1 = 0
         ifnd2 = 0

         do i=1,TmpFaceAll 
           j1 = TmpFace(i)%i1
           j2 = TmpFace(i)%i2
           j3 = TmpFace(i)%i3
           j4 = TmpFace(i)%i4

           jchk =j1 + j2 + j3 +j4

           if( ichk == jchk )then
             if( i4 > 0 )then
               if( i1 == j1 .or. i1 == j2 .or. i1 == j3 .or. i1 == j4 )then
                 if( i2 == j1 .or. i2 == j2 .or. i2 == j3 .or. i2 == j4 )then
                   if( i3 == j1 .or. i3 == j2 .or. i3 == j3 .or. i3 == j4 )then
                     if( i4 == j1 .or. i4 == j2 .or. i4 == j3 .or. i4 == j4 )then
                       ihit = i
                       ifnd = ifnd + 1
                       if( ifnd == 1 )then
                        ifnd1 = ihit
                       elseif( ifnd == 2 )then
                        ifnd2 = ihit
                       else
                         write(*,*)'bah! 1'
                       endif
                     endif
                   endif
                 endif
               endif    
             else     
               if( i1 == j1 .or. i1 == j2 .or. i1 == j3 )then
                 if( i2 == j1 .or. i2 == j2 .or. i2 == j3 )then
                   if( i3 == j1 .or. i3 == j2 .or. i3 == j3 )then
                     ihit = i
                     ifnd = ifnd + 1
                     if( ifnd == 1 )then
                      ifnd1 = ihit
                     elseif( ifnd == 2 )then
                      ifnd2 = ihit
                     else
                       write(*,*)'bah! 2'
                     endif
                   endif
                 endif
               endif        
             endif
           endif
         end do

         if( ihit > 0 )then

           if( ifnd2 == 0 )then
             !
             ! do nothing, just a shell at 
             ! domain boundary, ignore it
             !
             iloose = iloose + 1

             !write(*,*)' - ',ic,j,i,i9,ifnd1,ifnd2
             !write(*,*)' 1 ',tmpface(ifnd1)%cell1,tmpface(ifnd1)%cell2

           else

             !write(*,*)'hit',ic,j,i,i9,ifnd1,ifnd2
             !write(*,*)' 1 ',tmpface(ifnd1)%cell1,tmpface(ifnd1)%cell2
             !write(*,*)' 2 ',tmpface(ifnd2)%cell1,tmpface(ifnd2)%cell2
             !write(*,*)' 3 ',tmpface(ifnd1)%side1,tmpface(ifnd2)%side1

             tmpface(ifnd2)%cell1 = tmpface(ifnd1)%cell2
             tmpface(ifnd1)%cell2 = 0

             tmpface(ifnd2)%area  = tmpface(ifnd1)%area
             tmpface(ifnd2)%x     = tmpface(ifnd1)%x
             tmpface(ifnd2)%n     = tmpface(ifnd1)%n * (-1)

           endif
         else
           write(*,*) '+ Error: Internal error duplicating faces: ',i,j       
         endif

       endif
     end do
   endif
   !
   ! now the active faces are known assemble them in their final list 
   !
   Nfac = count( TmpFace(:)%cell1 > 0 )
   write(*,*)'Total number of active faces:',Nfac

   allocate( Face(1:Nfac) )
   !allocate( FaceNormal(1:Nfac,3) )

   icnt = 0
   
   do i=1,TmpFaceAll
     if( TmpFace(i)%cell1 > 0 )then
       icnt = icnt + 1
       Face(icnt)%bnd    = 0                  ! bnd region by definition
       Face(icnt)%cell1  = TmpFace(i)%cell1
       Face(icnt)%cell2  = TmpFace(i)%cell2
       Face(icnt)%area   = TmpFace(i)%area
       Face(icnt)%n      = TmpFace(i)%n
       Face(icnt)%x      = TmpFace(i)%x
       Face(icnt)%lambda = -1.0

       Face(icnt)%vertices(1) = TmpFace(i)%i1
       Face(icnt)%vertices(2) = TmpFace(i)%i2
       Face(icnt)%vertices(3) = TmpFace(i)%i3
       Face(icnt)%vertices(4) = TmpFace(i)%i4

       !write(IOdbg,'('' f'',i8,'' c:'',2(1x,i6),''  v:'',4(1x,i6),&
       !  &'' a: '',1pe9.3,'' x: '',3(1x,1pe10.3),'' n: '',3(1x,1pe10.3))') &
       !      icnt,Face(icnt)%cell1,Face(icnt)%cell2,Face(icnt)%vertices, &
       !      Face(icnt)%area,Face(icnt)%x,Face(icnt)%n

     endif
   end do
   if( icnt /= Nfac ) write(*,*) '+ Error: Internal error Nfac /= icnt ',Nfac,icnt

   deallocate( TmpFace )
   write(*,*)'Deallocated TmpFace'
   !
   ! assemble the final data structure
   !
   allocate( NFaces(1:Ncel) )
   NFaces = 0
   
   do i=1,Nfac
     if( Face(i)%cell1 > 0 )then
       j = Face(i)%cell1
       NFaces(j) = NFaces(j) + 1
       if( Face(i)%cell2 > 0 )then
         j = Face(i)%cell2
         NFaces(j) = NFaces(j) + 1
       else
         if( Face(i)%cell2 < 0 )then
           write(*,*) '+ Error: Internal error, face data structure corrupted 1',i
         endif   
       endif   
     else
       write(*,*) '+ Error: Internal error, face data structure corrupted 2',i
     endif   
   end do

   if( debug > 0 )then
     write(*,*)'Checking NFaces'
  
     ipoly2  = count( NFaces(:)  > 6  )
     ihex2   = count( NFaces(:) == 6  )
     iprism2 = count( NFaces(:) == 5  )
     itet2   = count( NFaces(:) == 4  )
     ierr2   = count( NFaces(:) <= 3  )

     if( ipoly2  > 0 )write(*,*)'Number of poly''s     :',ipoly2
     if( ihex2   > 0 )write(*,*)'Number of hexa''s     :',ihex2
     if( iprism2 > 0 )write(*,*)'Number of penta''s    :',iprism2
     if( itet2   > 0 )write(*,*)'Number of tetra''s    :',itet2
     if( ierr2   > 0 )write(*,*)'Number of weirdo''s   :',ierr2

   endif
   !
   ! the cell -> face pointers...
   !
   call set_up_add_face_to_cell(Ncel)

   do i=1,Nfac
     icell1 = Face(i)%cell1

     call add_face_to_cell(icell1,i)

     if( Face(i)%cell2 > 0 )then
       icell2 = Face(i)%cell2
       call add_face_to_cell(icell2,i)
     endif
   end do

   call debug_face_to_cell(Ncel)
   
   !
   ! next loop over all boundaries and find/flag them
   !
   Nbnd = count( Face(:)%cell2 == 0 )
   allocate( Bnd(1:Nbnd) )
   write(*,*)'Setting up boundary list: ',Nbnd
    
   icnt  = 0
   icnt0 = 0
   do i=1,Nfac
     if( mod(i,250000) == 0 ) write(*,*) ' ',i
     if( Face(i)%cell2 == 0 )then
       ! boundary face found
       i1 = Face(i)%vertices(1)
       i2 = Face(i)%vertices(2)
       i3 = Face(i)%vertices(3)
       i4 = Face(i)%vertices(4)

       ichk = i1 + i2 + i3 + i4
       ihit = 0
              
       do j=1,Norgbnd
         j1 = borg(j,1)
         j2 = borg(j,2)
         j3 = borg(j,3)
         j4 = borg(j,4)
       
         if( j4 == j3 ) j4 = -1
         
         jchk =j1 + j2 + j3 +j4

         if( ichk == jchk )then
           if( i4 > 0 )then
             if( i1 == j1 .or. i1 == j2 .or. i1 == j3 .or. i1 == j4 )then
               if( i2 == j1 .or. i2 == j2 .or. i2 == j3 .or. i2 == j4 )then
                 if( i3 == j1 .or. i3 == j2 .or. i3 == j3 .or. i3 == j4 )then
                   if( i4 == j1 .or. i4 == j2 .or. i4 == j3 .or. i4 == j4 )then
                     ihit = j
                   endif
                 endif
               endif
             endif    
           else     
             if( i1 == j1 .or. i1 == j2 .or. i1 == j3 )then
               if( i2 == j1 .or. i2 == j2 .or. i2 == j3 )then
                 if( i3 == j1 .or. i3 == j2 .or. i3 == j3 )then
                   ihit = j
                 endif
               endif
             endif        
           endif
         endif
       end do
       if( ihit == 0 )then
         icnt  = icnt  + 1
         icnt0 = icnt0 + 1
         !write(*,*) 'Default Bound: ',i,ichk,jchk
         Face(i)%bnd = icnt
         Bnd(icnt)%face        = i
         Bnd(icnt)%vertices(1) = i1
         Bnd(icnt)%vertices(2) = i2
         Bnd(icnt)%vertices(3) = i3
         Bnd(icnt)%vertices(4) = i4
         Bnd(icnt)%rid         =  0
       else if( ihit > 0 )then
         icnt = icnt + 1
         !write(*,*) 'Bound: ',i,j
         Face(i)%bnd = icnt
         Bnd(icnt)%face        = i
         Bnd(icnt)%vertices(1) = i1
         Bnd(icnt)%vertices(2) = i2
         Bnd(icnt)%vertices(3) = i3
         Bnd(icnt)%vertices(4) = i4
         Bnd(icnt)%rid         =  borg(ihit,5)
       else
         write(*,*) '+ Error: Internal error finding boundaries: ',i,j       
       endif
     endif
   end do
   write(*,*)'Boundary list set up'
   deallocate(borg)

   if( debug > 3 )then
     write(IOdbg,'(//,'' Boundaries 1'')') 
     icnt = 0
     do i=1,Nbnd
       if( Bnd(i)%rid > 0 )then
         icnt = icnt + 1
         write(IOdbg,'('' b a:'',i4,''=>'',4(1x,i4),1x,i4)') &
                                          i,Bnd(i)%vertices,icnt
       endif 
     end do
     write(IOdbg,'(/,'' Boundaries 2'')') 
     icnt0 = 0
     do i=1,Nbnd
       if( Bnd(i)%rid == 0 )then
         icnt0 = icnt0 + 1
         write(IOdbg,'('' b x:'',i4,''=>'',4(1x,i4),1x,i4)') &
                                          i,Bnd(i)%vertices,icnt0
       endif 
     end do
   endif
   
   write(*,*)'Total number of boundaries: ',icnt
   write(*,*)'Default boundaries:         ',icnt0
   write(*,*)'Assigned boundaries:        ',icnt-icnt0
   
   write(*,*)'Done reading geometry'

   call calc_interpolation_factor

   call calc_normal_distance
   
end subroutine readgeometry
subroutine calc_vol(volume,i1,i2,i3,i4,i5,i6,i7,i8)

   use geometry
   
   integer i1, i2, i3, i4, i5, i6, i7, i8
   integer isbw, isbe, inbe, inbw, istw, iste, inte, intw

   real    volume, vol(6)
   
   isbw = i1
   isbe = i2
   inbe = i3
   inbw = i4
   istw = i5
   iste = i6
   inte = i7
   intw = i8

   xsbw = Vert(isbw,1)
   xsbe = Vert(isbe,1)
   xnbe = Vert(inbe,1)
   xnbw = Vert(inbw,1)
   xstw = Vert(istw,1)
   xste = Vert(iste,1)
   xnte = Vert(inte,1)
   xntw = Vert(intw,1)

   ysbw = Vert(isbw,2)
   ysbe = Vert(isbe,2)
   ynbe = Vert(inbe,2)
   ynbw = Vert(inbw,2)
   ystw = Vert(istw,2)
   yste = Vert(iste,2)
   ynte = Vert(inte,2)
   yntw = Vert(intw,2)

   zsbw = Vert(isbw,3)
   zsbe = Vert(isbe,3)
   znbe = Vert(inbe,3)
   znbw = Vert(inbw,3)
   zstw = Vert(istw,3)
   zste = Vert(iste,3)
   znte = Vert(inte,3)
   zntw = Vert(intw,3)
   !
   ! define delta's relative to some corner say eg 1 as origin
   !
   xo   = xsbw
   yo   = ysbw
   zo   = zsbw

   dx12 = xsbe - xo
   dy12 = ysbe - yo
   dz12 = zsbe - zo

   dx13 = xnbe - xo
   dy13 = ynbe - yo
   dz13 = znbe - zo

   dx14 = xnbw - xo
   dy14 = ynbw - yo
   dz14 = znbw - zo

   dx15 = xstw - xo
   dy15 = ystw - yo
   dz15 = zstw - zo

   dx16 = xste - xo
   dy16 = yste - yo
   dz16 = zste - zo

   dx17 = xnte - xo
   dy17 = ynte - yo
   dz17 = znte - zo

   dx18 = xntw - xo
   dy18 = yntw - yo
   dz18 = zntw - zo

   vol(1:6) = 0.0
   !   
   ! each cell can be decomposed into 6 'tetra's' 
   ! so the total volume is the average of these enclosed parallepipeds
   ! the procedure works also for the degenerated hexa's
   !
   ! a better form is needed see eg. Wesseling 11.4
   ! or Ferziger 8.6.4, eg. based on 8 tetra's
   !
   vol(1) = determinant(dx16,dy16,dz16,dx12,dy12,dz12,dx17,dy17,dz17)
   vol(2) = determinant(dx15,dy15,dz15,dx16,dy16,dz16,dx17,dy17,dz17)
   vol(3) = determinant(dx18,dy18,dz18,dx15,dy15,dz15,dx17,dy17,dz17)
   vol(4) = determinant(dx12,dy12,dz12,dx13,dy13,dz13,dx17,dy17,dz17)
   vol(5) = determinant(dx13,dy13,dz13,dx14,dy14,dz14,dx17,dy17,dz17)
   vol(6) = determinant(dx14,dy14,dz14,dx18,dy18,dz18,dx17,dy17,dz17)

   f166 = 1./6.
   volume = abs( f166*sum( vol(1:6) ) )

end subroutine calc_vol
real function determinant(a1,a2,a3,b1,b2,b3,c1,c2,c3) 

   determinant = ( a2*b3 - b2*a3 )*c1 + ( b1*a3 - a1*b3 )*c2 + &
                 ( a1*b2 - a2*b1 )*c3

end function determinant
subroutine calc_interpolation_factor
   !
   ! interpolation factor for node 'P' to neighbor
   !
   use geometry
   
   real, dimension(3) :: xp, xf, xn, xnorm
   
   eps = epsilon(1.0)
   do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )then
       xp = Cell(ip)%x
       xf = Face(i)%x
       xn = Cell(in)%x
       
      !delt1 = sqrt( dot_product((xf - xp),(xf - xp)) )
      !delt2 = sqrt( dot_product((xn - xf),(xn - xf)) )
      !
      ! new:
      ! interpolation along the projection on the face normal
      !
       xnorm = Face(i)%n
       delt1 = dot_product( xnorm , (xf - xp) )/Face(i)%area
       delt2 = dot_product( xnorm , (xn - xf) )/Face(i)%area
           
       if( delt1+delt2 > eps )then
         Face(i)%lambda = delt1/(delt1+delt2)

       else
         Face(i)%lambda = huge(1.0)
       endif
     else
       Face(i)%lambda = -1.0       
     endif
   
   end do
   
end subroutine calc_interpolation_factor
subroutine calc_normal_distance
   !
   ! interpolation factor for node 'P' to neighbor
   ! see Ferziger 8.10.3
   ! 
   use geometry
   
   real, dimension(3) :: xp, xf, xv1, xv2, xv3, xv4
   real, dimension(3) :: x1, x2, n1, n2, n3, n4, n

   real :: length
   real :: pi = 3.1415927

   distmin  =  999999.
   distmax  = -999999.
   anglemin =  999999.
   anglemax = -999999.

   ierr = 0
   do i=1,Nbnd
     iface = Bnd(i)%face
     xf = Face(iface)%x

     ic1 = Face(iface)%cell1
     ic2 = Face(iface)%cell2
     if( ic2 > 0 )then
       write(*,*)'+ Error: Internal error boundary faces corrupted'
       ierr = ierr + 1
     endif
     xp = Cell(ic1)%x

     i1 = Bnd(i)%vertices(1) 
     i2 = Bnd(i)%vertices(2) 
     i3 = Bnd(i)%vertices(3) 
     i4 = Bnd(i)%vertices(4) 

     xv1 = Vert(i1,1:3)
     xv2 = Vert(i2,1:3)
     xv3 = Vert(i3,1:3)
     if( i4 > 0 ) xv4 = Vert(i4,1:3)
          
     if( i4 > 0 )then
       ! quad average of 4 triangles from center to corners

       x1 = xv1 - xf
       x2 = xv2 - xf
       call cross_product(x1,x2,n1)
       
       x1 = xv2 - xf
       x2 = xv3 - xf
       call cross_product(x1,x2,n2)

       x1 = xv3 - xf
       x2 = xv4 - xf
       call cross_product(x1,x2,n3)

       x1 = xv3 - xf
       x2 = xv1 - xf
       call cross_product(x1,x2,n4)
       
       n = 0.25*( n1 + n2 + n3 + n4 )

       x1 = xp - xf
       x2 = n

       call normalise(x1)
       call normalise(x2)
       
       dot = (x1(1)*x2(1)+x1(2)*x2(2)+x1(3)*x2(3) )
       if( dot <= 0.0 )then
         dot = max(-1.0,dot)
       else
         dot = min(dot,1.0)
       endif
       angle = acos( dot )
       
       if( angle > pi/2.0 )then
         n = -n
       endif
     else
       ! triangle is just one triangle
       x1 = xv2 - xv1
       x2 = xv3 - xv1
       call cross_product(x1,x2,n)

       x1 = xp - xf
       x2 = n
            
       call normalise(x1)
       call normalise(x2)
       
       dot = (x1(1)*x2(1)+x1(2)*x2(2)+x1(3)*x2(3) )
       if( dot <= 0.0 )then
         dot = max(-1.0,dot)
       else
         dot = min(dot,1.0)
       endif
       angle = acos( dot )
       
       if( angle > pi/2.0 )then
         n = -n
       endif
              
     endif

     x1 = xp - xf
     length = vector_length( x1 )     
     x2 = n
          
     call normalise(x1)
     call normalise(x2)
     
     dot = (x1(1)*x2(1)+x1(2)*x2(2)+x1(3)*x2(3) )
     if( dot <= 0.0 )then
       dot = max(-1.0,dot)
     else
       dot = min(dot,1.0)
     endif
     angle = acos( dot )

     !    n|   + xp   
     !     |  /.
     !     | / .  distance: |xp-xf| cos(angle)
     !     |/  .
     !    -+--------------
     !    xf
     !
     
     distance = length * cos( angle )

     if( distance <= 0.0 )then
       if( ierr <= 10 )then
         write(*,*)'+ Error: Internal error distance to wall'
         write(*,*)'  Boundary ',i,' : ',Bnd(i)%vertices
         write(*,*)'  In cell  ',ic1
         write(*,*)'  Xp : ',xp
         write(*,*)'  Xf : ',xf
         write(*,*)'  x1 : ',x1
         write(*,*)'  x2 : ',x2
         write(*,*)'  L  : ',length
         ierr = ierr + 1
       else if( ierr == 11 )then
         ierr = ierr + 1
         write(*,*)'+ No more shown.'
       endif
     endif

     Bnd(i)%distance = distance

     distmin = min(distmin,distance)
     distmax = max(distmax,distance)
     anglemin = min(anglemin,angle)
     anglemax = max(anglemax,angle)
     
  !   write(*,*)'',i,distance,angle*180./pi,i1,i2,i3,i4


   end do

   write(*,*)'Minimum distance to a wall: ',distmin
   write(*,*)'Maximum distance to a wall: ',distmax
   write(*,*)'Minimum wall angle (degs) : ',anglemin*180./pi
   write(*,*)'Maximum wall angle (degs) : ',anglemax*180./pi
      
end subroutine calc_normal_distance
subroutine cross_product(a,b,c)

   real, dimension(3) :: a, b, c
   
   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine cross_product
subroutine normalise(a)

   real, dimension(3) :: a
   real :: rlength
   
   rlength = 1.0/sqrt( dot_product(a,a) )
   
   a(1) = a(1)*rlength
   a(2) = a(2)*rlength
   a(3) = a(3)*rlength

end subroutine normalise
real function vector_length(a)

   real, dimension(3) :: a

   vector_length = sqrt( dot_product(a,a) )
   
end function vector_length
program preprocessor

   use constants
   use geometry

   type(CFadresses) ptr

   character(len=32) :: casename = 'dolfyn'
   character(len= 3) :: string   = 'bin'
   logical           :: bin      = .true.
   logical           :: exists   = .false.
   logical           :: CFGused  = .false.

   write(*,*)'Dolfyn PreProcessor'

   !
   ! dolfyn.cfg used here as well
   ! (fixed format file)
   !
   inquire(file='dolfyn.cfg',exist=exists)
   if( exists )then
     call openfile(IOcfg,'dolfyn','.cfg','FORMATTED', &
                         'OLD',debug)
     write(*,*) 'Reading case from dolfyn.cfg'
     read(IOcfg,'(A)') Casename 
     read(IOcfg,*)     ScaleFactor
     read(IOcfg,'(A)') String
     close(IOcfg)                         
     CFGused = .true.
   else
     write(*,*)'Input casename:'
     read(*,'(A)') Casename 
   endif

   write(*,*)'Using ',casename(1:lens(casename)),' as input'

   if( debug > 1 ) &
     call openfile(IOdbg,casename,'.dbg','FORMATTED','UNKNOWN',debug)
      
   call readgeometry(casename)

   if( .not. CFGused )then
 1   continue
     write(*,*)'Enter scaling factor (1.0):'
     read(*,*) ScaleFactor
     if( ScaleFactor <= 0.0 ) goto 1
     write(*,*)'Using: ',ScaleFactor

 2   continue
     write(*,*)'Enter format of geometry file (bin|ascii):'
     read(*,*) String
     if( String == 'BIN' ) string = 'bin'
     if( String == 'ASC' ) string = 'asc'
     if( String == 'Bin' ) string = 'bin'
     if( String == 'Asc' ) string = 'asc'
     if( String /= 'bin' )then
       if( String /= 'asc' ) goto 2
     endif
   endif
   
   if( String == 'bin' ) bin = .true.
   if( String == 'asc' ) bin = .false.

   if( bin )then
     write(*,*)'Dump binary geometry file'
     call openfile(IOgeo,casename,'.geo','UNFORMATTED','UNKNOWN',debug)
   else
     write(*,*)'Dump ascii geometry file'
     call openfile(IOgeo,casename,'.geo','FORMATTED','UNKNOWN',debug)
   endif

   if( bin )then  !1234567812345678
     write(IOgeo) 'dolfyn bin g'
     write(IOgeo) version
     write(IOgeo) ScaleFactor
     !    
     ! write active cells as a list of faces, a pointer to the original cell number
     ! and the fluid type id
     !
     write(IOgeo) 'ce: '
     write(IOgeo) Ncel
     do i=1,Ncel
       write(IOgeo) i,cell(i)%ctid,cell(i)%x,cell(i)%vol
     end do   

     write(IOgeo) 'cf: '
     write(IOgeo) maxval(NFaces)
     do i=1,Ncel
       write(IOgeo) i,NFaces(i)
       ptr%ptr => CFhead(i)%ptr
       do while( associated(ptr%ptr) )
         j = ptr%ptr%face
         write(IOgeo) j
         ptr%ptr => ptr%ptr%next
       end do
     end do   
     !
     ! write the faces
     !
     write(IOgeo) 'fc: '
     write(IOgeo) Nfac
     do i=1,Nfac
       write(IOgeo) i,face(i)
     end do
     !
     ! write the boundaries
     !
     write(IOgeo) 'bn: '
     write(IOgeo) Nbnd,maxval(bnd%rid)
     do i=1,Nbnd
       write(IOgeo) i,bnd(i)
     end do
     !
     ! write the vertices
     !
     write(IOgeo) 'vr: '
     write(IOgeo) Nvrt
     do i=1,Nvrt
       write(IOgeo) i,(vert(i,j),j=1,3)
     end do
     !
     ! write active cells in standard format as they came in
     ! (should not be needed, but is useful for postprocessors!)
     !
     write(IOgeo) 'oc: '
     write(IOgeo) Ncel
     do i=1,Ncel
       write(IOgeo) i,Cell(i)%corg,cell(i)%vertices 
     end do  
     write(IOgeo) 'eof '
  
   else
     !
     ! ascii/text/debug format
     !
     write(IOgeo,'(a)')  'dolfyn asc g'
     write(IOgeo,'(a,1x,i6)') 'version',version
     write(IOgeo,'(a,e12.5)') 'scale:', ScaleFactor
     !    
     ! write active cells as a list of faces, a pointer to the original cell number
     ! and the fluid type id
     !
     write(IOgeo,'(a)') 'ce: cells'
     write(IOgeo,*) Ncel
     do i=1,Ncel
       write(IOgeo,'(2(i8,1x),4(1pe16.9,1x))') i,cell(i)%ctid,cell(i)%x,cell(i)%vol
     end do   

     write(IOgeo,'(a)') 'cf: cell faces'
     write(IOgeo,'(i8)',advance='no') maxval(NFaces)
     do i=1,Ncel
       write(IOgeo,'(/,2(i8,1x))',advance='no') i,NFaces(i)
       ptr%ptr => CFhead(i)%ptr
       do while( associated(ptr%ptr) )
         j = ptr%ptr%face
         write(IOgeo,'(i8,1x)',advance='no') j
         ptr%ptr => ptr%ptr%next
       end do
     end do   
     write(IOgeo,'(/)',advance='no') 
     !
     ! write the faces
     !
     write(IOgeo,'(a)') 'fc: faces'
     write(IOgeo,*) Nfac
     do i=1,Nfac
       write(IOgeo,'(8(i8,1x),8(1pe16.9,1x))') i,face(i)
     end do
     !
     ! write the boundaries
     !
     write(IOgeo,'(a)') 'bn: boundaries'
     write(IOgeo,*) Nbnd,maxval(bnd%rid)
     do i=1,Nbnd
       write(IOgeo,'(7(i8,1x),1(e12.5,1x))') i,bnd(i)
     end do
     !
     ! write the vertices
     !
     write(IOgeo,'(a)') 'vr: vertices'
     write(IOgeo,*) Nvrt
     do i=1,Nvrt
       write(IOgeo,'(1(i8,1x),3(1pe16.9,1x))') i,(vert(i,j),j=1,3)
     end do
     !
     ! write active cells in standard format as they came in
     ! (should not be needed, but is useful for postprocessors!)
     !
     write(IOgeo,'(a)') 'oc: original cells'
     write(IOgeo,*) Ncel
     do i=1,Ncel
       write(IOgeo,'(10(i8,1x))') i,Cell(i)%corg,cell(i)%vertices 
     end do  
     write(IOgeo,'(a)') 'eof '

   endif
   
end 
