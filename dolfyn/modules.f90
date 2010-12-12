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
! With contributions by B. Tuinstra  
! (www.home.zonnet.nl/bouke_1/dolfyn)
!
module constants
   !
   ! global constants
   !
   integer, parameter :: version = 0420 

   integer, parameter :: Solid = 1, Fluid = 2, Ghost = 3, Baffle = 4, Shell = 5

   integer, parameter :: Hexa = 1, Prism = 2, Pyramid  = 3, Tetra = 5,     &
                         Poly = 6, Quad  = 7, Triangle = 8, Polygon = 9

   integer, parameter :: NFhex  = 6, NFprism = 5, NFpyrmd = 5, NFtet = 4,  &
                         NFpoly = 8, NFquad  = 2, NFtri = 2

   integer, parameter :: IOgeo = 8,  IOcel = 14, IOvrt = 15, IObnd = 23,   &
                         IOdbg = 63, IOinp = 10, IOpst = 13, IOrst =  9,   &
                         IOmtx = 11, IOres =  7, IOdef =  6, IOcfg = 12,   &
                         IOext = 17, IOmon = 12, IOprt = 18

   integer, parameter :: RWall = 1, RInlet = 2, ROutlet = 3, RPressure = 4,&
                         RSymp = 5

   integer, parameter :: TMnone = 0, TMkeps = 1, TMrng = 2, TMcube = 3

   integer, parameter :: GradLS = 0, GradGauss = 1
   integer            :: GradAlg = GradLS

   character(len=8)   :: Region(5)

   data Region / 'Wall    ', 'Inlet   ', 'Outlet  ', 'Pressure',&
                 'Symplane' /

   real               :: Small, Large

   real, parameter    :: Zero = 0.0, One = 1.0

   integer, parameter :: DSud  = 1, DScd1 = 2, DScd2 = 3, &
                         DSgam = 4, DSlud = 5, DSmod = 6

   character(len=3)   :: DSch(6)
   
   data DSch / 'UD ','CD ','CD2','GAM','LUD','MMD' /
   !
   ! character array with konsole output flags
   !
   integer, parameter    :: NFlags = 5
   character*1, dimension(NFlags) :: Flags
 
   integer, parameter    :: IFlagP = 1, IFlagInflow = 2, IFlagMass = 3, &
                            IFlagTE = 4, IFlagVis = 5

   !
   ! keywords for variable density, viscosity, Cp
   !
   integer, parameter :: Constant = 1, Ideal = 2, Isobaric = 3, Inviscid = 4, &
                         Multicomponent = 4, Nonnewtonian = 5, User = 6,      &
                         Polynomial = 7
   !
   ! materials 
   !                         
   integer            :: Density = Constant     
   integer            :: LaminarViscosity = Constant

   !
   ! symbolic constants for number of variables
   !
   integer, parameter :: NVar  = 13               ! number of standard variables
   integer, parameter :: MinSC = NVar + 1         ! start of scalars
   integer            :: MaxVar, MaxSC, MaxMat    ! no parameter, will be set

   integer, parameter :: VarU   =  1, VarV   =  2, VarW   =  3, VarP   =  4,  &
                         VarTE  =  5, VarED  =  6, VarT   =  7, VarSC  =  8,  &
                         VarDen =  9, VarPP  = 10, VarVis = 11, VarLvis= 12,  &
                         VarCP  = 13      

   integer, allocatable :: VarS(:)

   character(len=64) :: casename = 'dolfyn'
   character(len=96) :: title    = 'dolfyn'
   
   integer   :: memory     =   0
   integer   :: Debug      =   0

   integer   :: NSave      = 500                 ! save every
   integer   :: ISave      =  -1                 ! save just iter i
   real      :: TSave      =  -1.                ! save only at time
   logical   :: TSaveDone  = .false.
   real      :: TCPU       =  -1.                ! save every elapsed  cpu time
   logical   :: CPUStop    = .false.             ! stop or continue

   integer   :: NOutput    = 500                 ! write every
   integer   :: IOutput    =  -1                 ! write just iter i
   real      :: TOutput    =  -1.                ! write only at time
   logical   :: TOutputDone= .false.


   logical   :: Transient  = .false.
   real      :: dt         = 0.001
   logical   :: Euler      = .true.
   logical   :: QuadTime   = .false.
   real      :: GammaTime  = 1.0
   
   integer   :: Restart    = 0
   logical   :: MovingGrid = .false.

   logical   :: SolveU          = .true.
   logical   :: SolveV          = .true.
   logical   :: SolveW          = .true.
   logical   :: SolveUVW        = .true.
   
   logical   :: SolveP          = .true.
   logical   :: SolveEnthalpy   = .false.
   logical   :: SolveTurb       = .false.
   logical   :: SolveVisc       = .false.
   logical   :: SolveScalars    = .false.
   logical   :: SolveTurbEnergy = .false.
   logical   :: SolveTurbDiss   = .false.

   logical   :: SolveDensity    = .false.
   logical   :: SolveCP         = .false.

   logical   :: Initialisation     = .false.
   logical   :: UserInitialisation = .false.
   integer   :: InitialSteps       =  10
   real      :: VisInitialFactor   =   1.0

   !
   ! data print options
   !
   character(len=64) :: PrintFile  = '-'

   logical   :: PrintCellVar       = .false.
   logical   :: PrintCellVarUser   = .false.

   integer   :: IPrintCellVarStart = 0
   integer   :: IPrintCellVarEnd   = 0
   integer   :: IPrintCellVarInc   = 1

   logical   :: PrintWallVar       = .false.
   logical   :: PrintWallVarUser   = .false.

   integer   :: IPrintWallVarStart = 0
   integer   :: IPrintWallVarEnd   = 0
   integer   :: IPrintWallVarInc   = 1
   !
   ! general options
   !
   logical   :: UsePatches         = .false.
   logical   :: UseScalars         = .false.
   logical   :: UseParticles       = .false.

   logical   :: UseArtificialComp  = .false.
   logical   :: UseGMV             = .false.
   logical   :: UseVTK             = .false.
   logical   :: UseTECPLOT         = .false.
   logical   :: UseLapack          = .false. 
   logical   :: UseOpenDX          = .true.
   !
   ! special ABL options
   !
   logical   :: UseFixABL          = .false.
   integer   :: IdFixABL           = 0
   real      :: UFixABL            = 0.0
   real      :: VFixABL            = 0.0
   real      :: WFixABL            = 0.0
   real      :: TeFixABL           = 0.0
   real      :: EdFixABL           = 0.0

   integer, parameter :: SparseKit2      =  1  ! Sparsekit2 BiCGstab by Youcef Saad
   integer, parameter :: Hypre           =  2  ! Hypre 
   integer, parameter :: SolverBiCGstabL =  3  ! van der Vorst BiCGstab Ell
   integer, parameter :: SolverGMRESR    =  4  ! van der Vorst GMRESr
   integer, parameter :: SolverICCG      =  5  !  
   integer, parameter :: SolverDirect    =  6  ! Harwell HSL MA27 (testing only)
   integer, parameter :: SolverUser      =  7  ! easy way to allow for your own solver

   logical   :: DoubleU   = .false.         ! obsolete
   logical   :: DoubleV   = .false.
   logical   :: DoubleW   = .false.
   logical   :: DoubleP   = .false. 
   logical   :: DoubleE   = .false. 
   logical   :: DoubleS   = .false.

   integer   :: TurbModel = TMnone
   integer   :: Nouter    =   1 
   integer   :: Niter     = 100 
   integer   :: Ngradient =   2

   !
   ! stand. k-epsilon turb. model constants
   !
   real      :: Kappa   = 0.419
   real      :: TMCmu   = 0.09
   real      :: TMCeps1 = 1.44
   real      :: TMCeps2 = 1.92
   real      :: TMLenSc = 0.1
   
   !
   ! rarely changed switches/constants
   !
   real      :: ParticleCourant = 0.35         ! default Courant number

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
   integer :: Nvrt = 0, Ncel = 0, Nbnd = 0, Nreg  = 0, &
              Nfac = 0, Nint = 0, NNZ = 0 
   !
   ! vertex section
   !
   real, allocatable :: Vert(:,:)   ! array of vertices

   !
   ! scalefactor
   !
   real :: ScaleFactor = 1.0
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
   integer, allocatable :: NFaces(:)     ! number of faces for a comp. cell
   integer, allocatable :: CFace(:,:)    ! list of faces
   real, allocatable    :: RFace(:,:)    ! Ae, Aw coef. list (face -> cell)

   real, allocatable    :: Ar(:)         ! reciprocal coef. of A (pres. cor.)

   real, allocatable    :: Au(:)         ! coef. for the cell (Au) 
   real, allocatable    :: Su(:)         ! cell source term (Su)       

   real, allocatable    :: Av(:)         ! coef. for the cell (Av) 
   real, allocatable    :: Sv(:)         ! cell source term (Sv)       

   real, allocatable    :: Aw(:)         ! coef. for the cell (Aw) 
   real, allocatable    :: Sw(:)         ! cell source term (Sw)       

   !
   ! SparsKit2 / Hypre solver work arrays
   !
   double precision, allocatable :: Work(:,:)      

   double precision, allocatable :: Acoo(:)   ! A-matrix in COO-format      
   integer, allocatable          :: Arow(:)   ! A-matrix row entries
   integer, allocatable          :: Acol(:)   ! A-matrix column entries

   double precision, allocatable :: Acsr(:)   ! A-matrix in CSR-format      
   integer, allocatable          :: Arwc(:)   ! A-matrix row entries
   integer, allocatable          :: Aclc(:)   ! A-matrix column entries

   double precision, allocatable :: ALU(:)    ! U-matrix in CSR-format      
   integer, allocatable          :: JLU(:)    ! U-matrix row entries
   integer, allocatable          :: JU(:)     ! U-matrix column entries
   integer, allocatable          :: JW(:)     ! ILU-matrix  
   integer, allocatable          :: IPERM(:)  ! ILU-matrix  

   double precision, allocatable :: RHS(:)    ! righthand side    
   double precision, allocatable :: SOL(:)    ! solution    

   !double precision, allocatable :: Acsr_DP(:)      
   !double precision, allocatable :: Acoo_DP(:)      
   !
   ! cell section
   !
   type :: CellData
     real    x(3)                        ! center coordinates
     real    vol                         ! volume
     integer ctid                        ! fluid type id as set in fluid table
   end type
   
   type(CellData), allocatable :: Cell(:) ! array of computational cells
   
   !
   ! boundary section (both set by input and default boundaries)
   !
   type :: BoundaryData
     integer :: face                     ! belongs to face...
     integer, dimension(4) :: vertices   ! the 4 vertices 
     integer :: rid                      ! region id as set in rtable
     real    :: distance                 ! normal distance from cell face 
                                         ! center to cell center
     real    :: yplus                    ! y+
     real    :: uplus                    ! u+
     real    :: shear(3)                 ! shearstress components
     real    :: h                        ! local heattransfer coef.
     real    :: q                        ! local heat flux (in W/m2)
     real    :: T                        ! local wall temperature
   end type
   
   type(BoundaryData), allocatable :: Bnd(:) ! array of boundaries

   !
   ! boundary region definitions
   !
   type :: RegionData
     integer :: typ                      ! type
     real    :: uvw(3)                   ! place for velocities
     real    :: den = 1.20500            ! place for density
     real    :: T = 293.0                ! place for temperature
     real    :: R = 0.0                  ! place for resistance   
     real    :: k = 0.0001               ! place for turb. kin. energy
     real    :: e = 1.0                  ! place for turb. diss.
     logical :: noslip = .true.          ! slip/noslip switch
     logical :: std  = .true.            ! standard/roughness switch
     real    :: elog = 9.0               ! E-parameter
     real    :: ylog = 11.0              ! y+ value where u+=y+ matches wall law
     real    :: z0   = 0.03              ! roughness parameter
     logical :: split = .true.           ! split/fixed flow rate switch
     real    :: splvl = 1.0 
     logical :: adiab = .true.           ! adiab./fixed temperature switch
     logical :: flux  = .false.          ! fixed flux temperature switch
     logical :: user  = .false.          ! call special user subroutine?     
     character(len=4) :: name1 = '    '  ! IDs of special user subroutine     
     character(len=4) :: name2 = '    '  ! IDs of special user subroutine     
   end type 
   
   type(RegionData), allocatable :: Reg(:) ! array of regions

end module geometry
module scalars 
   !
   ! scalar boundary region definitions
   !
   type :: ScalarRegionData
     integer :: typ      =  -1            ! type
     real    :: fraction =  0.0           ! boundary mass fraction
     logical :: adiab    = .true.         ! adiab./fixed temperature switch
     logical :: flux     = .false.        ! fixed flux temperature switch
     real    :: value    =  0.0           ! value
     logical :: user     = .false.        ! call special user subroutine?     
   end type 
   
   type(ScalarRegionData), allocatable :: ScReg(:,:)  ! array of scalar regions

   real, allocatable :: ScFluxIn(:,:)     ! store influx for each scalar
   real, allocatable :: ScFluxOut(:,:)    ! store outflux for each scalar
   
   !
   ! molecular scalar properties
   !
   type :: ScalarData
     character(len=12) name              ! name
     character(len=12) unit              ! unit

     logical :: active                   ! passive (simple tranport) or active 
     logical :: solve                    ! use it or not
     logical :: save                     ! if true, variable field will be written on restart file

     real    :: PrL = 1.0                ! laminar Prandtl number
     real    :: PrT = 1.0                ! turbulent Prandtl number

     real    :: den                      ! boundary mass fraction (kg/m3)
     real    :: mol                      ! molecular weight (kg/kmol)
     real    :: beta                     ! thermal expansion coefficient (1/K)
     real    :: vis                      ! molecular viscosity (kg/ms)
     real    :: cp                       ! specific heat (J/kgK)
     real    :: lambda                   ! conductivity (W/mK)
     real    :: heat                     ! heat of formation (J/kg)
     real    :: temp                     ! temp. of formation (K)
   end type 
   
   type(ScalarData), allocatable :: ScProp(:)  ! array of scalar properties

end module scalars 
!========================================================================
module artificial_compressibility

    !double precision, save :: delta_t !the pseudo timestep
    double precision, dimension(3) :: sum_VS          !factor for pseudo time step
    double precision, dimension(3) :: sum_area        !the sum of the areas in the direction of the flow
    double precision, dimension(3) :: U_P, V_P, W_P   !local velocity vectors
    double precision :: cfl = 1d50                     !courant criterium
    double precision, allocatable :: art_snd_speed(:) !artificial sound speed
    double precision :: ac_beta = 1d50                 !art snd speed prefactor
    double precision, allocatable :: delta_t(:)       !the pseudo timestep
    double precision, allocatable :: density_change(:)!artificial change in density
    integer :: iter

end module artificial_compressibility
!========================================================================
module variables
   !
   ! main variables defined at the cell centers 
   !   
   integer :: Nscal = 0, Nmat = 0


   character(len=12), allocatable  :: Variable(:)
   character(len=12), allocatable  :: Scalar(:)

   real, allocatable :: U(:),V(:),W(:), TE(:), ED(:), &
                        VisEff(:), TurbP(:), T(:), Res(:) 

   real, allocatable :: Cp(:),DEN(:), DENold(:), DENdp(:) 

   real, allocatable :: Uold(:), Vold(:), Wold(:),  TEold(:),  EDold(:),  Told(:) 
   real, allocatable :: Uold2(:),Vold2(:),Wold2(:), TEold2(:), EDold2(:), Told2(:) 

   real, allocatable :: SC(:,:), SCold(:,:), SCold2(:,:)

   real, allocatable :: P(:), PP(:)
   !double precision, allocatable  :: P(:) <=moet op een gegeven moment

   real, allocatable    :: Residual(:)
   real, allocatable    :: ResiNorm(:)
   
   integer, allocatable :: Solver(:)           ! blending factor UDS en HigherOrder DS.
   real, allocatable    :: Gamma(:)            ! blending factor UDS en HigherOrder DS.
   integer, allocatable :: Scheme(:)           ! choice of DS

   logical, allocatable :: Solve(:)            ! simple array filled in after readcontrolfile
   logical, allocatable :: Store(:)            ! simple array filled in after readcontrolfile

   logical, allocatable :: SolveScalar(:)      ! solve extra scalars
   logical, allocatable :: StoreScalar(:)      ! solve extra scalars

   logical, allocatable :: PostC(:)            ! postprocessing arrays, cells
   logical, allocatable :: PostV(:)            ! vertices

   !
   ! variables defined at the cell faces only
   !
   real, allocatable :: MassFlux(:)
   real, allocatable :: MassFluxDebug(:,:)

   real, allocatable :: DXdebug(:), DXgrad(:,:)
   !
   ! gradients 
   !
   real, allocatable :: dudx(:,:), dvdx(:,:), dwdx(:,:), &
                        dpdx(:,:)

   real, allocatable :: dPhidXo(:,:) 
   !
   ! materials section
   !
   !real, allocatable :: Tref(:), Pref(:), DensRef(:), VisMol(:)

   integer :: IPref  =   1        ! ref. cell for pressure (material 1)
   integer :: IMoni  =   1        ! monitor cell for pressure (material 1)
   integer :: Iter   =   0        ! iteration counter
   real    :: Time   =   0.0      ! current time (if applicable)

   real :: Tref      = 273.       ! ref. cell for temperature (material 1)
   real :: Pref      =   0.0      ! ref. cell for pressure (material 1)
   
   real :: DensRef   =   1.2      ! <= later lucht op 20C
   real :: VisLam    =   0.001    ! <= later lucht op 20C

   real :: Prandtl   =   0.6905   ! lucht op 20C  
   real :: Schmidt   =   0.9 

   real :: Sigma_T   =   0.9      !
   real :: Sigma_k   =   1.0      ! turbulence diff. coef. factors
   real :: Sigma_e   =   1.219    ! ie. turbulent Prandtl numbers
   real :: Sigma_s   =   0.9      ! 

   real :: Gravity(3) =     0.0   ! <= gravity vector
   real :: Beta       =     0.001 ! expansie coef.
   real :: CpStd      =  1006.0   !      
   real :: CvStd      =  1006.0   !
   real :: Lambda     =     0.02637 ! warmtegeleiding lucht

   real :: Qtransfer  = 0.0
   !                      u    v    w    p     k    e      T   Sc 
   real :: URF(8)   = (/ 0.5, 0.5, 0.5, 0.2,  0.5, 0.5,   0.5, 0.5 /) 
   real :: RTOL(8)  = (/ 0.1, 0.1, 0.1, 0.05, 0.1, 0.1,   0.1, 0.1 /) 
   real :: ATOL(8)  = (/ 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,   0.0, 0.0 /) 
   real :: Guess(8) = (/ 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 293.0, 0.0 /)  !<==!!!!

   real    :: ResMax   = 1.e-4
   integer :: MaxPCOR  =    4
   real    :: FactDPP  = 0.25 
   !
   ! special post processing section
   !      
   logical :: DXnormals     = .false.
   logical :: DXcenters     = .false.
   logical :: DXmassflux    = .false.
   logical :: DXtemperature = .false. 
   logical :: DXdebugdata   = .false.
   integer :: DXdump        =    100
   
end module variables
