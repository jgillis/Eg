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
subroutine ScanControlFile

   use constants
   use geometry
   use variables
   use scalars
   use particles
   
   character(len=120) string

   character(len=32)  :: key, key2, key3

   integer, parameter :: MaxKeys = 20
   character(len=16)  :: keys(MaxKeys) 
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   integer, save      :: MathMode = 0

   write(IOdbg,*)'Start scan of input deck'

   call openfile(IOinp,casename,'.din','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)

   MaxSC  = 0
   MaxMat = 1
   
   do 
     read(IOinp,'(A)',end=20) string   ! read string

     if( string(1:1) == '#' ) cycle    ! skip if comment line

     indx = index(string,'#')
     if( indx >= 2 )then               ! strip trailing comments
       do i=indx,len(string)
         string(i:i) = ' '
       end do
     endif

     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

     if( Nkeys == 0 ) cycle            ! skip blank line

     key = keys(1)
     call lowercase(key)

     select case (key) 

       case ('dimension')
         key = keys(2)
         call lowercase(key)
         select case (key(1:3))
           case ('scalars')
             idum = keyi(3)
             if( idum < 1 )then
               write(*,*)'+ error number of scalars should be > 0 :',idum
             endif

             MaxSC = max(MaxSC,abs(idum))
             write(IOdbg,'(1x,A,i3)')'+ scalars: MaxSC=',MaxSC 
           
           case ('materials')
             idum = keyi(3)
             if( idum < 1 )then
               write(*,*)'+ error number of materials should be > 0 :',idum
             endif

             MaxMat = max(MaxMat,abs(idum))
             write(IOdbg,'(1x,A,i3)')'+ materials: MaxMat=',MaxMat
           
           case default
             write(*,*)'+ error unknown dimension command: ',string
           
         end select
         
       case ('use')
         key = keys(2)
         call lowercase(key)
         select case (key(1:7))
           case ('patches')
             idum = keyi(3)
             if( idum < 0 )then
               write(*,*)'+ error number of scalars should be >= 0 :',idum
             endif 

             MaxSC = max(MaxSC,idum)
             write(IOdbg,'(1x,A,i3)')'+ scalars: MaxSC=',MaxSC 
           case ('particl')
             idum = keyi(3)
             if( idum < 0 )then
               write(*,*)'+ error number of particles should be >= 0 :',idum
             endif 

             Npart = idum
             write(IOdbg,'(1x,A,i3)')'+ particles: Npart=',Npart 
         end select

       case ('set')
         key = keys(2)
         call lowercase(key)
         !write(*,*)'set>',key
         r1 = keyr(3)
         r2 = keyr(4)

         call setvariable(key(1:8),r1,r2)

       case ('math')
         key2 = keys(2)
         call lowercase(key2)

         if( key2(1:3) == 'deg' )then
           write(IOdbg,*)'Mathmode set to degrees'
           MathMode = 0
         else if( key2(1:3) == 'rad' )then
           write(IOdbg,*)'Mathmode set to radians'
           MathMode = 1         
         else
         
           key3 = keys(3)
           call lowercase(key3)

           r1  = keyr(4)
           r2  = keyr(5)
           r3  = keyr(6)
           r4  = keyr(7)

           call mathvariable(key3,MathMode,r1,r2,r3,r4)

           r2  = 0.0        
           call setvariable(key2(1:8),r1,r2)

         endif

       case default
         !write(*,*) '+ skipped keyword: ',key(1:lens(key))

     end select

   end do
20 continue
   write(IOdbg,*)'End scan of input deck'

   close(IOinp)
   rewind(IOinp)

   !
   ! dimension array's !!!
   !
   ! scalars
   !
   MaxVar = NVar + MaxSC

   write(IOdbg,*)'MaxSc  = ',MaxSC
   write(IOdbg,*)'MaxVar = ',MaxVar
   
   allocate(Variable(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Variable array allocated')

   allocate(Residual(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Residual array allocated')

   allocate(ResiNorm(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'ResiNorm array allocated')

   Residual = 0.0
   ResiNorm = 1.0

   allocate(Solver(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Solver array allocated')

   Solver   = SparseKit2

   allocate(Gamma(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Gamma array allocated')

   Gamma    = 0.0

   allocate(Scheme(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Scheme array allocated')

   Scheme   = DScd1

   allocate(Solve(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Solve array allocated')

   Solve(:)        = .false.
   Solve(1:4)      = .true. 

   allocate(Store(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Store array allocated')

   Store(:)        = .false.
   Store(1:4)      = .true. 

   allocate(PostC(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'PostC array allocated')

   PostC(:)        = .false. 
   PostC(1:4)      = .true. 

   allocate(PostV(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'PostV array allocated')

   PostV(:)        = .false.

   if( MaxSC > 0 )then
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

     do i=1,MaxSC
       indx    = Nvar + i
       VarS(i) = indx
       write(key,'(''Scalar'',i3.3,''   '')') i 
       Variable(indx) = key(1:12)
       Scalar(i)      = key(1:12)
     end do
     
     NScal = 0
     
   endif  
   
   Variable(VarU )   = 'U       '
   Variable(VarV )   = 'V       '
   Variable(VarW )   = 'W       '
   Variable(VarP )   = 'Pressure'

   Variable(VarTE)   = 'KE      '
   Variable(VarED)   = 'ED      '
   Variable(VarT )   = 'T       '
   Variable(VarSC)   = 'Scalar  '

   Variable(VarDen)  = 'Density '
   Variable(VarPP)   = 'PressCor'
   Variable(VarVis)  = 'VisTurb '
   Variable(VarLVis) = 'VisLam  '
   Variable(VarCP)   = 'CP      '

   write(IOdbg,*)'Variable control arrays dimensioned'

   !
   ! particles
   !
   if( Npart > 0 )then
     allocate(Particle(NPart),stat=istat)
     call TrackMemory(istat,Npart,'Particle array allocated')
   endif

end subroutine ScanControlFile
subroutine ReadControlFile

   use constants
   use geometry
   use variables
   use scalars
   use particles
   
   character(len=120) string

   character(len=32)  :: key, key2, key3, key4, blank

   integer, parameter :: MaxKeys = 20
   character(len=16)  :: keys(MaxKeys) 
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   integer            :: ideltas(MaxKeys)
   real               :: rdeltas(MaxKeys)

   integer, save      :: MathMode = 0

   logical :: Messages = .false.
   logical :: Echo     = .false.
   
   logical :: LGenerator
   
   Noutlets = 0
   Split    = 0.0
   
   blank      = '                                '

   Small = tiny(small) * 1.e+3
   Large = huge(large) * 1.e-3
   
   SMALL = 1.e-12

   call SetUpVariables
   call ScanControlFile

   call SetUpVariables
   call openfile(IOinp,casename,'.din','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)
   
   do 
     read(IOinp,'(A)',end=20)    string   ! read string
     !write(IOdbg,'(1x,''>'',A)') string

     if( string(1:1) == '#' ) cycle    ! skip if comment line

     if( string(1:5) == 'title' .and.                     &
       ( string(6:6) == ',' .or. string(6:6) == ' ' .or.  &
         string(6:6) == '=') )then
       title = string(7:len(string))
       !write(*,'(1x,A,A)')'+ title ',title(1:lens(title))
       cycle
     endif

     indx = index(string,'#')
     if( indx >= 2 )then               ! strip trailing comments
       do i=indx,len(string)
         string(i:i) = ' '
       end do
     endif

     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

     if( Nkeys == 0 ) cycle           ! skip blank line

     key = keys(1)
     call lowercase(key)

     !write(*,*)'using key:',key

     select case (key) 

       case ('debug')
         i = keyi(2)
         debug = i
         if( debug > 1 ) Echo = .true.
         if( Echo ) write(*,'(1x,A,i2)')'+ debug level = ',i

       case ('echo')
         key = keys(2)
         call lowercase(key)

         select case (key(1:3))
           case ('no')
             if( Echo ) write(*,'(1x,A,i3)')'+ echo input off'
             Echo = .false.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ echo input off'
             Echo = .false.
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ echo input on'
             Echo = .true.
           case default
             if( Echo ) write(*,'(1x,A,i3)')'+ echo input on'
             Echo = .true.
         end select

       case ('steps')
         i = keyi(2)
         if( Echo ) write(*,*)'+ iterations/time steps = ',i
         Niter = i
  
         if( keyr(3) > 0.0 )then
           ResMax = keyr(3)
           if( Echo ) write(*,'(1x,A,1pe9.3)')'+ with target residual ',ResMax
         endif

       case ('pref')
         !
         ! nog chekken of in de juiste range/materiaal is enz.
         !
         i = keyi(2)
         if( i > Ncel )then
           write(*,*)'+++ Warning: Pressure reference cel out of range'
           i = Ncel 
         endif
         
         if( Echo ) write(*,*)'+ reference pressure in cell = ',i
         IPref = i

       case ('monitor')
         !
         ! nog chekken of in de juiste range/materiaal is enz.
         !
         i = keyi(2)
         if( i > Ncel )then
           write(*,*)'+++ Warning: Monitor cel out of range'
           i = Ncel 
         endif

         if( Echo ) write(*,*)'+ monitor cell = ',i
         IMoni = i

       case ('save')
         !
         ! save restart data
         !
         key = keys(2)
         call lowercase(key)

         select case (key(1:4))
           case ('ever')
             i = keyi(3)
             if( i < 0 )then
               write(*,*)'+++ Error: invalid entry in "save,every,":',i
               messages = .true.
             else
               if( Echo ) write(*,'(1x,A,i3)')'+ saving every ',i,' steps'
               NSave = i
             endif
           case ('time')
             r = keyr(3)
             if( r < 0.0 )then
               write(*,*)'+++ Error: invalid entry in "save,time,":',r
               messages = .true.
             else
               if( Echo ) write(*,'(1x,A,1pe10.3)')'+ saving at time =',r
               TSave = r
             endif
           case ('iter')
             i = keyi(3)
             if( i < 0 )then
               write(*,*)'+++ Error: invalid entry in "save,iteration,":',i
               messages = .true.
             else
               if( Echo ) write(*,'(1x,A,i3)')'+ saving iteration ',i
               ISave = i
             endif
           case ('cpu')
             r = keyr(3)
             key2 = keys(4)
             call lowercase(key2)
             key3 = keys(5)
             call lowercase(key3)

             select case (key2(1:1))
               case ('s')
                 fact =    1.0
               case ('m')
                 fact =   60.0
               case ('h')
                 fact = 3600.0
               case default
                 fact =   1.0
             end select

             select case (key3(1:1))
               case ('s')
                 CPUStop = .true.
               case default
                 CPUStop = .false.
             end select             
             
             if( r < 0.0 )then
               write(*,*)'+++ Error: invalid entry in "save,cpu,":',r
               messages = .true.
             else
               r = fact * r
               if( Echo )  write(*,'(1x,A,f6.0)') &
                         '+ saving after elapsed CPU time seconds:',r
               TCPU = r
             endif
           case default
             write(*,*)'+++ Error: unknown option in save command.'
             messages = .true.
         end select

       case ('output')
         !
         ! write output data...
         !
         key = keys(2)
         call lowercase(key)

         select case (key(1:4))
           case ('ever')
             i = keyi(3)
             if( i < 0 )then
               write(*,*)'+++ Error: invalid entry in "output,every,":',i
               messages = .true.
             else
               if( Echo ) write(*,'(1x,A,i3)')'+ output every ',i,' steps'
               NOutput = i
             endif
           case ('time')
             r = keyr(3)
             if( r < 0.0 )then
               write(*,*)'+++ Error: invalid entry in "output,time,":',r
               messages = .true.
             else
               if( Echo ) write(*,'(1x,A,1pe10.3)')'+ output at time =',r
               TOutput = r
             endif
           case ('iter')
             i = keyi(3)
             if( i < 0 )then
               write(*,*)'+++ Error: invalid entry in "output,iteration,":',i
               messages = .true.
             else
               if( Echo ) write(*,'(1x,A,i3)')'+ output iteration ',i
               IOutput = i
             endif
           case default
             write(*,*)'+++ Error: unknown option in output command.'
             messages = .true.
         end select

       case ('boundary')
         !
         ! boundaries section
         !
         i = keyi(2)
         if( Echo ) write(*,'(1x,A,i3)')'+ boundary category = ',i
         if( i < 0 .or. i > Nreg )then
           write(*,*)'Boundary out of range. Skipped.'
           messages = .true.
           goto 3
         endif
         read(IOinp,'(A)',err=10,end=20) string
         key = blank
         call getkeyword(string,key,jkey)
         call lowercase(key)
         if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)
         
         select case (key(1:3))
           case ('sym')
             write(*,'(1x,A,i3,A)') 'Boundary ',i,' sym. conditions'
             Reg(i)%typ = RSymp
           case ('inl')
             !
             ! inlet
             !
             write(*,'(1x,A,i3,A)') 'Boundary ',i,' inlet conditions'
             Reg(i)%typ = RInlet
             if( keys(3) == 'user' ) Reg(i)%user  = .true. 
             if( keys(4) /= '    ' ) Reg(i)%name1 = keys(4)
             if( keys(5) /= '    ' ) Reg(i)%name2 = keys(5)
             if( Reg(i)%user ) write(*,*) '  *** User inlet ',&
                                          Reg(i)%name1,' ',Reg(i)%name2

!            read(IOinp,*,err=11) Uin,Vin,Win
!            Reg(i)%uvw(1) = Uin
!            Reg(i)%uvw(2) = Vin
!            Reg(i)%uvw(3) = Win

             read(IOinp,'(A)',end=20) string    
             call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
             Uin = keyr(1)
             Vin = keyr(2)
             Win = keyr(3)
             Reg(i)%uvw(1) = Uin 
             Reg(i)%uvw(2) = Vin 
             Reg(i)%uvw(3) = Win 
             if( Echo ) write(*,*) '  UVW:',Reg(i)%uvw,' m/s'

!            read(IOinp,*,err=11) DENin
!            Reg(i)%den    = DENin
             read(IOinp,'(A)',end=20) string    
             call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
             Reg(i)%den = keyr(1)
             if( Echo ) write(*,*) '  DEN:',Reg(i)%den,' kg/m3'

!            read(IOinp,*,err=11) Tin
!            Reg(i)%T      = Tin
             read(IOinp,'(A)',end=20) string    
             call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
             Reg(i)%T= keyr(1)
             if( Echo ) write(*,*) '  T  :',Reg(i)%T,' K'

             Vmag2 = Uin**2 + Vin**2 + Win**2

             read(IOinp,'(A)',end=20) string
             key = blank
             call getkeyword(string,key,jkey)
             call lowercase(key)
             if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:4)
             if( key(1:4) == 'keps' )then
               if( Echo ) write(*,*) '  Using k and epsilon'
!              read(IOinp,*,err=11) TKin,EPSin
               read(IOinp,'(A)',end=20) string    
               call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
               TKin  = keyr(1)
               EPSin = keyr(2)

               Reg(i)%k = TKin
               Reg(i)%e = EPSin
               if( TKin <= 0.0 )then
                 write(*,*)'Error: k should be positive! Reseting.'
                 ti = 0.01
                 Reg(i)%k = 3./2.*( ti**2 * vmag2 ) 
                 messages = .true.
               endif
               if( EPSin <= 0.0 )then
                 write(*,*)'Error: epsilon should be positive! Reseting.'
                 tl = 1.0
                 Reg(i)%e = (TMCmu)**(0.75) * Reg(i)%k**(1.5) / tl
                 messages = .true.
               endif
               TKin  = sqrt( 2./3.*Reg(i)%k / (Vmag2 + Small) )
               EPSin = (TMCmu)**(0.75) * Reg(i)%k**(1.5) / (Reg(i)%e + Small)
               write(*,*) '  TE :',Reg(i)%k,'( i',TKin,')'
               write(*,*) '  ED :',Reg(i)%e,'( l',EPSin,')'
             elseif( key(1:4) == 'inle' )then
               if( Echo ) write(*,*) '  Using intensity and length scale'
!              read(IOinp,*,err=11) TKin,EPSin
               read(IOinp,'(A)',end=20) string    
               call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
               TKin  = keyr(1)
               EPSin = keyr(2)
               if( TKin <= 0.0 )then
                 write(*,*)'Error: intensity should be positive! Reseting.'
                 ti   = 0.01
                 TKin = ti
                 Reg(i)%k = 3./2.*( ti**2 * vmag2 ) 
                 messages = .true.
               else
                 Reg(i)%k = 3./2.*( TKin**2 * vmag2 ) 
               endif
               if( EPSin <= 0.0 )then
                 write(*,*)'Error: length scale should be positive! Reseting.'
                 tl    = 1.0
                 EPSin = tl
                 Reg(i)%e = (TMCmu)**(0.75) * Reg(i)%k**(1.5) / tl
                 messages = .true.
               else
                 Reg(i)%e = (TMCmu)**(0.75) * Reg(i)%k**(1.5) / EPSin
               endif
               if( Echo ) write(*,*) '  TE :',TKin,'( k',Reg(i)%k,')'
               if( Echo ) write(*,*) '  ED :',EPSin,'( e',Reg(i)%e,')'
             else
               write(*,*)'Error: [keps*|inlen] expected'
               messages = .true.
             endif
           case ('out')
             ! 
             ! outlet 
             !
             write(*,'(1x,A,i3,A)') 'Boundary ',i,' outlet conditions'
             Reg(i)%typ    = ROutlet
             read(IOinp,*,err=11) SplitIn
             Reg(i)%splvl  = SplitIn  
             
             Noutlets = Noutlets + 1
             Split    = Split + SplitIn
             
             !write(*,*)'>> ',Noutlets,SplitIn,Split
             if( Split > 1.0 )then
               write(*,*)'Error: Sum of splits greater than 1.0'
               messages = .true.
             endif
           case ('wal')
             !
             ! wall
             !
             write(*,'(1x,A,i3,A)') 'Boundary ',i,' wall conditions'
             Reg(i)%typ = RWall
             read(IOinp,'(A)',end=20) string
             key = blank
             call getkeyword(string,key,jkey)
             call lowercase(key)
             if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:4)
             if( key(1:4) == 'slip' )then
               if( Echo ) write(*,*) '  Free slip wall'
               Reg(i)%noslip = .false.
               Reg(i)%adiab  = .true.
               Reg(i)%flux   = .false.
               goto 3
             else if( key(1:4) == 'nosl' )then
               Reg(i)%noslip = .true.
             else
               write(*,*) '*** Error: [noslip*|slip] expected'
               write(*,*) '*** Assuming noslip wall'
               Reg(i)%noslip = .false.
               messages = .true.
             endif
             !
             ! noslip wall
             !
             if( Echo ) write(*,*) '  No slip wall'
             read(IOinp,'(A)',end=20) string    
             call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
             Reg(i)%uvw(1) = keyr(1) 
             Reg(i)%uvw(2) = keyr(2) 
             Reg(i)%uvw(3) = keyr(3) 
             if( Echo ) write(*,*) '  UVW:',Reg(i)%uvw

             read(IOinp,'(A)',end=20) string
             key = blank
             call getkeyword(string,key,jkey)
             call lowercase(key)
             if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:4)
             if( Echo ) write(*,*)'+ note: resistance R added after T or q!'
             if( key(1:4) == 'adia' )then
               if( Echo ) write(*,*) '  Adiabatic wall'
               Reg(i)%adiab  = .true.
               Reg(i)%flux   = .false.
               goto 3
             elseif( key(1:4) == 'fixe' )then
               Reg(i)%adiab  = .false.
               Reg(i)%flux   = .false.
               read(IOinp,'(A)',end=20) string    
               call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
               Reg(i)%T = keyr(1) 
               Reg(i)%R = keyr(2) 

               if( Echo ) write(*,*) '  T,R:',Reg(i)%T,Reg(i)%R
               goto 3
             elseif( key(1:4) == 'flux' )then
               Reg(i)%adiab  = .false.
               Reg(i)%flux   = .true.
               read(IOinp,'(A)',end=20) string    
               call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
               Reg(i)%T = keyr(1) 
               Reg(i)%R = keyr(2) 

               if( Echo ) write(*,*) '  T  :',Reg(i)%T,Reg(i)%R
               goto 3
             else
               write(*,*)'Error: [adiabatic*|fixed|flux] expected'
               messages = .true.
               goto 3
             endif
         end select
 3       continue
       case ('ngradient')
         i = keyi(2)
         if( Echo ) write(*,'(1x,A,i3)')'+ gradient iterations = ',i
         Ngradient = i

       case ('thermal')
         if( Echo ) write(*,'(1x,A,i3)')'+ thermal options'
         key = keys(2)
         call lowercase(key)
         if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)

         select case (key(1:3))
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ thermal mode switched off'
             SolveEnthalpy = .false.
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ thermal mode switched on'
             SolveEnthalpy = .true.
           case default
             SolveEnthalpy = .true.
         end select

       case ('turbulence')
         if( Echo ) write(*,'(1x,A,i3)')'+ turbulence options'
         key = keys(2)
         call lowercase(key)
         if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)
         select case (key(1:3))
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ turbulence modelling mode switched off'
             SolveTurb       = .false.
             SolveVisc       = .false.
             SolveTurbEnergy = .false.
             SolveTurbDiss   = .false.
             TurbModel       =  TMnone
           case ('ke')
             if( Echo ) write(*,'(1x,A,i3)')'+ standard k-e turbulence model'
             SolveTurb       = .true.
             SolveVisc       = .true.
             SolveTurbEnergy = .true.
             SolveTurbDiss   = .true.
             TurbModel       =  TMkeps
             if( keyr(3) > Small )then
               TMLenSc = keyr(3)
             endif
           case ('rng')
             if( Echo ) write(*,'(1x,A,i3)')'+ rng k-e turbulence model'
             SolveTurb       = .true.
             SolveVisc       = .true.
             SolveTurbEnergy = .true.
             SolveTurbDiss   = .true.
             TurbModel       =  TMrng
         end select

       case ('scalars')
         !
         write(*,*)'*** Scalars: experimental feature ***'
         !
         if( Echo ) write(*,'(1x,A,i3)')'+ scalar options'
         key = keys(2)
         call lowercase(key)
         if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)
         select case (key(1:3))
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ scalars switched off'
             SolveScalars = .false.
           case ('on ')
             if( Echo ) write(*,'(1x,A,i3)')'+ scalars switched on'
             SolveScalars = .true.
             UseScalars   = .true.
         end select

       case ('restart')
         key = keys(2)
         call lowercase(key)
         if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)
         select case (key(1:3))
           case ('no')
             if( Echo ) write(*,'(1x,A,i3)')'+ no restart file used'
             Restart       = 0
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ no restart file used'
             Restart       = 0
           case ('ini')
             !
             ! if ReadRestartField discovers that Nbnd and Nfaces
             ! have changed Restart is set there to 4
             !
             if( Echo ) write(*,'(1x,A,i3)')'+ use old results as initial geuss'
             Restart       = 2
           case ('res')
             if( Echo ) write(*,'(1x,A,i3)')'+ just reset counters'
             Restart       = 3
           case default
             if( Echo ) write(*,'(1x,A,i3)')'+ run restarted'
             Restart       = 1
         end select

       case ('gravity')
         Gravity(1) = keyr(2)
         Gravity(2) = keyr(3)
         Gravity(3) = keyr(4)
         if( Echo ) write(*,'(1x,A,3(1x,f6.3))')'+ gravity vector = ',Gravity(1:3)

       case ('beta')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ volumetric expansion coef = ',dummy
         beta = dummy
         
       case ('transient')
         dummy = keyr(2)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ transient, time step (dt) = ',dummy
         dt = dummy
         Transient = .true.
         Euler     = .true.
         QuadTime  = .false.
         if( keys(3) == 'quad' )then
           Euler    = .false.
           QuadTime = .true.
           GammaTime = keyr(4)
         endif
         !
         ! RESET RELAXATION FACTORS!
         !
         !URF(VarU ) = 1.0
         !URF(VarV ) = URF(VarU)
         !URF(VarW ) = URF(VarU)
         !URF(VarP ) = 0.5
         !URF(VarTE) = 1.0
         !URF(VarED) = URF(5)
         !URF(VarT ) = URF(5)
         !URF(VarSC) = URF(5)

       case ('vislam')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ laminar viscosity = ',dummy
         vislam = dummy

       case ('density')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ density = ',dummy
         DensRef = dummy

       case ('cp')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ specific heat Cp = ',dummy
         CpStd = dummy

       case ('cv')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ specific heat Cv = ',dummy
         CvStd = dummy

       case ('prandtl')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ Prandtl number = ',dummy
         Prandtl = dummy
         Lambda  = VisLam * CpStd / Prandtl

       case ('conductivity')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ conductivity = ',dummy
         Lambda = dummy
         Prandtl = VisLam * CpStd / Lambda

       case ('schmidt')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ Schmidt number = ',dummy
         Schmidt = dummy

       case ('relax*')
         URF(VarU ) = keyr(2)
         URF(VarV ) = keyr(3)
         URF(VarW ) = keyr(4)
         URF(VarP ) = keyr(5)
         URF(VarTE) = keyr(6)
         URF(VarED) = keyr(7)
         URF(VarT ) = keyr(8)
         URF(VarSC) = keyr(9)

         if( URF(VarU ) == 0.0 ) URF(VarU ) = 0.7   ! set defaults
         if( URF(VarV ) == 0.0 ) URF(VarV ) = 0.7
         if( URF(VarW ) == 0.0 ) URF(VarW ) = 0.7
         if( URF(VarP ) == 0.0 ) URF(VarP ) = 0.35
         if( URF(VarTE) == 0.0 ) URF(VarTE) = 0.7
         if( URF(VarED) == 0.0 ) URF(VarED) = 0.7
         if( URF(VarT ) == 0.0 ) URF(VarT ) = 0.9
         if( URF(VarSC) == 0.0 ) URF(VarSC) = 0.9

         call check_minmax(URF(VarU ),Small,1.0)
         call check_minmax(URF(VarV ),Small,1.0)
         call check_minmax(URF(VarW ),Small,1.0)
         call check_minmax(URF(VarP ),Small,1.0)
         call check_minmax(URF(VarTE),Small,1.0)
         call check_minmax(URF(VarED),Small,1.0)
         call check_minmax(URF(VarT ),Small,1.0)
         call check_minmax(URF(VarSC),Small,1.0)
         if( Echo )then
           write(*,'(1x,A,3(1x,f6.3))')'+ under relaxation factors (*) => '
           write(*,*) '+ uvwp: ',(URF(iurf),iurf=1,4)
           write(*,*) '+ kets: ',(URF(iurf),iurf=5,8)
         endif

       case ('relax')

         URF(VarU ) = keyr(2)
         URF(VarV ) = URF(VarU)
         URF(VarW ) = URF(VarU)
         URF(VarP ) = keyr(3)
         URF(VarTE) = keyr(4)
         URF(VarED) = URF(VarTE)
         URF(VarT ) = keyr(5)
         URF(VarSC) = keyr(6)

         if( URF(VarU ) == 0.0 ) URF(VarU ) = 0.7   ! set defaults
         if( URF(VarV ) == 0.0 ) URF(VarV ) = 0.7
         if( URF(VarW ) == 0.0 ) URF(VarW ) = 0.7
         if( URF(VarP ) == 0.0 ) URF(VarP ) = 0.35
         if( URF(VarTE) == 0.0 ) URF(VarTE) = 0.7
         if( URF(VarED) == 0.0 ) URF(VarED) = 0.7
         if( URF(VarT ) == 0.0 ) URF(VarT ) = 0.9
         if( URF(VarSC) == 0.0 ) URF(VarSC) = 0.9
         
         call check_minmax(URF(VarU) ,Small,1.0)
         call check_minmax(URF(VarP) ,Small,1.0)
         call check_minmax(URF(VarTE),Small,1.0)
         call check_minmax(URF(VarT) ,Small,1.0)
         call check_minmax(URF(VarSC),Small,1.0)
         if( Echo ) &
           write(*,'(1x,A,5(1x,f6.3))')'+ under relaxation factors  = ',&
           URF(VarU),URF(VarP),URF(VarTE),URF(VarT),URF(VarSC)
         
       case ('rtol') 
         RTOL(VarU ) = keyr(2)
         RTOL(VarV ) = RTOL(VarU)
         RTOL(VarW ) = RTOL(VarU)
         RTOL(VarP ) = keyr(3)
         RTOL(VarTE) = keyr(4)
         RTOL(VarED) = RTOL(VarTE)
         RTOL(VarT ) = keyr(5)
         RTOL(VarSC) = keyr(6)         

         if( RTOL(VarU ) == 0.0 ) RTOL(VarU ) = 0.1   ! set defaults
         if( RTOL(VarV ) == 0.0 ) RTOL(VarV ) = 0.1
         if( RTOL(VarW ) == 0.0 ) RTOL(VarW ) = 0.1
         if( RTOL(VarP ) == 0.0 ) RTOL(VarP ) = 0.05
         if( RTOL(VarTE) == 0.0 ) RTOL(VarTE) = 0.1
         if( RTOL(VarED) == 0.0 ) RTOL(VarED) = 0.1
         if( RTOL(VarT ) == 0.0 ) RTOL(VarT ) = 0.1
         if( RTOL(VarSC) == 0.0 ) RTOL(VarSC) = 0.1

         call check_minmax(RTOL(VarU) ,Small,1.0)
         call check_minmax(RTOL(VarP) ,Small,1.0)
         call check_minmax(RTOL(VarTE),Small,1.0)
         call check_minmax(RTOL(VarT) ,Small,1.0)
         call check_minmax(RTOL(VarSC),Small,1.0)
         if( Echo ) &
           write(*,'(1x,A,5(1x,f6.3))')'+ relative solver tolerance = ',&
           RTOL(VarU),RTOL(VarP),RTOL(VarTE),RTOL(VarT),RTOL(VarSC)
         
       case ('atol')
         ATOL(VarU ) = keyr(2)
         ATOL(VarV ) = ATOL(VarU)
         ATOL(VarW ) = ATOL(VarU)
         ATOL(VarP ) = keyr(3)
         ATOL(VarTE) = keyr(4)
         ATOL(VarED) = ATOL(5)
         ATOL(VarT ) = keyr(5)
         ATOL(VarSC) = keyr(6)       
         call check_minmax(ATOL(1),-Small,1.0)
         call check_minmax(ATOL(4),-Small,1.0)
         call check_minmax(ATOL(5),-Small,1.0)
         if( Echo ) &
           write(*,'(1x,A,5(1x,1pe7.1))')'+ absolute solver tolerance = ',&
           ATOL(VarU),ATOL(VarP),ATOL(VarTE),ATOL(VarT),ATOL(VarSC)

       case ('init')
         key2 = keys(2)
         call lowercase(key2)
         key3 = keys(3)
         call lowercase(key3)
         if( Echo ) &
           write(*,'(1x,A,A,A)')'+ init keyword with strings = ',key2,key3
         select case (key2(1:5))
           case ('field')
             Guess(VarU ) = keyr(3)
             Guess(VarV ) = keyr(4)
             Guess(VarW ) = keyr(5)
             Guess(VarP ) = keyr(6)
             Guess(VarTE) = keyr(7)
             Guess(VarED) = keyr(8)
             if( keyr(9) > Small ) Guess(VarT ) = keyr(9)
             Guess(VarSC) = keyr(10)
             if( Echo ) write(*,'(1x,A,4(1x,1pe8.1))')'+ initial fields = ',Guess(1:4)
             if( Echo ) write(*,'(1x,A,4(1x,1pe8.1))')'                   ',Guess(5:8)
           case ('user')
             if( Echo ) write(*,'(1x,A)') '+ use user subroutine'
             UserInitialisation = .true.
           case ('fact')
             VisInitialFactor = keyr(3)
             if( Echo ) write(*,'(1x,A,1pe8.1)')'+ initial viscosity factor  = ',&
                        VisInitialFactor
           case ('steps')
             InitialSteps = keyi(3)
             if( Echo ) write(*,'(1x,A,i3)')'+ initial steps = ',InitialSteps
           case default
             write(*,'(1x,A,i3)')'+ error in initial command'
             messages = .true.
         end select
             
       case ('pcor')
         key2 = keys(2)
         call lowercase(key2)
         key3 = keys(3)
         call lowercase(key3)
         if( Echo ) &
           write(*,'(1x,A,A,A)')'+ pcor keyword with strings = ',key2,key3
         select case (key2(1:3))
           case ('max')
             MaxPCOR = keyi(3)
             if( Echo ) write(*,'(1x,A,i3)')'+ maximum corr. steps = ',MaxPCOR
           case ('fac')
             FactDPP = keyr(3)
             if( Echo ) write(*,'(1x,A,f6.3)')'+ target factor = ',FactDPP
           case default
             write(*,'(1x,A,i3)')'+ error in pcor command'
             messages = .true.
         end select

       case ('gamma*')
         Gamma(VarU ) = keyr(2)                                      
         Gamma(VarV ) = keyr(3)                                      
         Gamma(VarW ) = keyr(4)                                      
         Gamma(VarP ) = keyr(5)                                      
         Gamma(VarTE) = keyr(6)                                      
         Gamma(VarED) = keyr(7)                                        
         Gamma(VarT ) = keyr(8)                                      
         Gamma(VarSC) = keyr(9)                                     
         if( Echo ) write(*,'(1x,A,4(1x,1pe7.1))')'+ blending (*) = ',Gamma(1:4)    
         if( Echo ) write(*,'(1x,A,4(1x,1pe7.1))')'                 ',Gamma(5:8)    
         
       case ('gamma')
         Gamma(VarU ) = keyr(2)                                      
         Gamma(VarV ) = Gamma(VarU )                                      
         Gamma(VarW ) = Gamma(VarU )                                      
         Gamma(VarP ) = keyr(3)         ! not used!                                  
         Gamma(VarTE) = keyr(4)                                      
         Gamma(VarED) = Gamma(VarTE)                                      
         Gamma(VarT ) = keyr(5)                                      
         Gamma(VarSC) = keyr(6)                                     
         if( Echo ) write(*,'(1x,A,4(1x,1pe7.1))')'+ blending = ',Gamma(1:4)    
         if( Echo ) write(*,'(1x,A,4(1x,1pe7.1))')'             ',Gamma(5:8)    

       case ('double')
         !
         ! obsolete!!!
         !
         key2 = keys(2)
         call lowercase(key2)
         key3 = keys(3)
         call lowercase(key3)
         
         if( Echo ) write(*,'(1x,A,A)')'+ double keyword with string = ',key2
         if( key2(1:3) == 'uvw' )then
           if( Echo ) write(*,'(1x,A)')'+ velocities with double precision'
           DoubleU   = .true.
           DoubleV   = .true.
           DoubleW   = .true. 
         endif 
         if( key2(1:3) == 'u  ' )then
           if( Echo ) write(*,'(1x,A)')'+ u velocity comp. with double precision'
           DoubleU   = .true.
         endif 
         if( key2(1:3) == 'v  ' )then
           if( Echo ) write(*,'(1x,A)')'+ v velocity comp. with double precision'
           DoubleV   = .true.
         endif 
         if( key2(1:3) == 'w  ' )then
           if( Echo ) write(*,'(1x,A)')'+ w velocity comp. with double precision'
           DoubleW   = .true.
         endif 
         if( key2(1:3) == 'p  ' )then
           if( Echo ) write(*,'(1x,A)')'+ pressure with double precision'
           DoubleP   = .true.
         endif 
         if( key2(1:3) == 'ke ' )then
           if( Echo ) write(*,'(1x,A)')'+ turbulence/scalars with double precision'
           DoubleS   = .true.
         endif 
         if( key2(1:3) == 't  ' )then
           if( Echo ) write(*,'(1x,A)')'+ temperature with double precision'
           DoubleE   = .true.
         endif 

       case ('solver')
         key2 = keys(2)
         call lowercase(key2)
         key3 = keys(3)
         call lowercase(key3)

         select case (key2)
           case ('u')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving U component'
                 SolveU      = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ U component not solved'
                 SolveU      = .false.
               case ('sparse')
                 if( Echo ) write(*,'(1x,A,i3)')'+ U component SparseKit2'
                 solver(VarU) = SparseKit2 
               case ('bicg')
                 if( Echo ) write(*,'(1x,A,i3)')'+ U component BiCGstabL'
                 solver(VarU) = SolverBiCGstabL 
               case ('direct')
                 if( Echo ) write(*,'(1x,A,i3)')'+ U component Direct'
                 solver(VarU) = SolverDirect 
               case ('user')
                 if( Echo ) write(*,'(1x,A,i3)')'+ U component User'
                 solver(VarU) = SolverUser
               case default
                 write(*,'(1x,A,i3)')'+ error in command solver'
                 write(*,'(1x,A,i3)')'+ solving U component'
                 SolveU      = .true.
                 messages = .true.
             end select
           case ('v')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving V component'
                 SolveV      = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ V component not solved'
                 SolveV      = .false.
               case ('sparse')
                 if( Echo ) write(*,'(1x,A,i3)')'+ V component SparseKit2'
                 solver(VarV) = SparseKit2 
               case ('bicg')
                 if( Echo ) write(*,'(1x,A,i3)')'+ V component BiCGstabL'
                 solver(VarV) = SolverBiCGstabL 
               case ('direct')
                 if( Echo ) write(*,'(1x,A,i3)')'+ V component Direct'
                 solver(VarV) = SolverDirect 
               case ('user')
                 if( Echo ) write(*,'(1x,A,i3)')'+ V component User'
                 solver(VarV) = SolverUser 
               case default
                 if( Echo ) write(*,'(1x,A,i3)')'+ error in command solver'
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving V component'
                 SolveV      = .true.
                 messages = .true.
             end select
           case ('w')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving W component'
                 SolveW      = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ W component not solved'
                 SolveW      = .false.
               case ('sparse')
                 if( Echo ) write(*,'(1x,A,i3)')'+ W component SparseKit2'
                 solver(VarW) = SparseKit2 
               case ('bicg')
                 if( Echo ) write(*,'(1x,A,i3)')'+ W component BiCGstabL'
                 solver(VarW) = SolverBiCGstabL 
               case ('direct')
                 if( Echo ) write(*,'(1x,A,i3)')'+ W component Direct'
                 solver(VarW) = SolverDirect 
               case ('user')
                 if( Echo ) write(*,'(1x,A,i3)')'+ W component User'
                 solver(VarW) = SolverUser
               case default
                 if( Echo ) write(*,'(1x,A,i3)')'+ error in command solver'
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving W component'
                 SolveW      = .true.
                 messages = .true.
             end select
           case ('p')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving pressure'
                 SolveP      = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ pressure not solved'
                 SolveP      = .false.
               case ('sparse')
                 if( Echo ) write(*,'(1x,A,i3)')'+ pressure SparseKit2'
                 solver(VarP) = SparseKit2 
               case ('bicg')
                 if( Echo ) write(*,'(1x,A,i3)')'+ pressure BiCGstabL'
                 solver(VarP) = SolverBiCGstabL 
               case ('direct')
                 if( Echo ) write(*,'(1x,A,i3)')'+ pressure Direct'
                 solver(VarP) = SolverDirect 
               case ('user')
                 if( Echo ) write(*,'(1x,A,i3)')'+ pressure User'
                 solver(VarP) = SolverUser 
               case default
                 if( Echo ) write(*,'(1x,A,i3)')'+ error in command solver'
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving pressure'
                 SolveP      = .true.
                 messages = .true.
             end select
           case ('t')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving temperature'
                 SolveEnthalpy = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature not solved'
                 SolveEnthalpy = .false.
               case ('sparse')
                 if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature SparseKit2'
                 solver(VarT) = SparseKit2 
               case ('bicg')
                 if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature BiCGstabL'
                 solver(VarT) = SolverBiCGstabL 
               case ('direct')
                 if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature Direct'
                 solver(VarT) = SolverDirect 
               case ('user')
                 if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature User'
                 solver(VarT) = SolverUser 
               case default
                 if( Echo ) write(*,'(1x,A,i3)')'+ error in command solver'
                 if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature not solved'
                 SolveEnthalpy = .false.
                 messages = .true.
             end select
           case ('k')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving turb. kin. energy k'
                 SolveTurbEnergy = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ turb. kin. energy k not solved'
                 SolveTurbEnergy = .false.
               case default
                 write(*,'(1x,A,i3)')'+ error in command solver'
                 write(*,'(1x,A,i3)')'+ turb. kin. energy k not solved'
                 SolveTurbEnergy = .false.
                 messages = .true.
             end select
           case ('eps')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ solving turb. dissipation epsilon'
                 SolveTurbDiss = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ turb. dissipation not solved'
                 SolveTurbDiss = .false.
               case default
                 write(*,'(1x,A,i3)')'+ error in command solver'
                 write(*,'(1x,A,i3)')'+ turb. dissipation not solved'
                 SolveTurbDiss = .false.
                 messages = .true.
             end select
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
             messages = .true.
         end select
!
! choice of differencing schemes
!
       case ('scheme*')
         Scheme(VarU ) = keyi(2)                                      
         Scheme(VarV ) = keyi(3)                                      
         Scheme(VarW ) = keyi(4)                                      
        !Scheme(VarP ) = keyi(5)                                      
         Scheme(VarTE) = keyi(6)                                      
         Scheme(VarED) = keyi(7)                                        
         Scheme(VarT ) = keyi(8)                                      
         Scheme(VarSC) = keyi(9)                                     
         if( Echo ) write(*,'(1x,A,4(1x,i2))')'+ scheme (*) = ',Scheme(1:4)    
         if( Echo ) write(*,'(1x,A,4(1x,i2))')'               ',Scheme(5:8)    

       case ('scheme')
         key2 = keys(2)
         call lowercase(key2)
         key3 = keys(3)
         call lowercase(key3)
         dum1 = keyr(4)
         dum2 = keyr(5)
         call check_minmax(dum1,0.0,1.0)
         call check_minmax(dum2,0.0,0.5)

         select case (key2)
           case ('uvw')
             select case (key3(1:3))
               case ('ud ')
                 Scheme(VarU) = DSud
                 Scheme(VarV) = DSud
                 Scheme(VarW) = DSud
                 if( Echo ) write(*,*)'+ velocity components: UD'
               case ('cd1')
                 Scheme(VarU) = DScd1
                 Scheme(VarV) = DScd1
                 Scheme(VarW) = DScd1
                 if( dum1 /= 0.0 )then
                   Gamma(VarU) = dum1
                   Gamma(VarV) = dum1
                   Gamma(VarW) = dum1
                 endif
                 if( Echo ) write(*,*)'+ velocity components: CD, blend ',dum1
               case ('cd2')
                 Scheme(VarU) = DScd2
                 Scheme(VarV) = DScd2
                 Scheme(VarW) = DScd2
                 if( dum1 /= 0.0 )then
                   Gamma(VarU) = dum1
                   Gamma(VarV) = dum1
                   Gamma(VarW) = dum1
                 endif
                 if( Echo ) write(*,*)'+ velocity components: CD2, blend ',dum1
               case ('cd ')
                 Scheme(VarU) = DScd1
                 Scheme(VarV) = DScd1
                 Scheme(VarW) = DScd1
                 if( dum1 /= 0.0 )then
                   Gamma(VarU) = dum1
                   Gamma(VarV) = dum1
                   Gamma(VarW) = dum1
                 endif
                 if( Echo ) write(*,*)'+ velocity components: CD, blend ',dum1
               case ('lud')
                 Scheme(VarU) = DSlud
                 Scheme(VarV) = DSlud
                 Scheme(VarW) = DSlud
                 if( dum1 /= 0.0 )then
                   Gamma(VarU) = dum1
                   Gamma(VarV) = dum1
                   Gamma(VarW) = dum1
                 endif
                 if( Echo ) write(*,*)'+ velocity components: LUD, blend ',dum1 
               case ('gam')
                 Scheme(VarU) = DSgam
                 Scheme(VarV) = DSgam
                 Scheme(VarW) = DSgam
                 if( dum1 /= 0.0 )then
                   Gamma(VarU) = dum1
                   Gamma(VarV) = dum1
                   Gamma(VarW) = dum1
                 endif
                !if( dum2 /= 0.0 ) Beta_m = dum2
                 if( Echo ) write(*,*)'+ velocity components: Gamma, blend ',dum1!,dum2
               case ('min')
                 Scheme(VarU) = DSmod
                 Scheme(VarV) = DSmod
                 Scheme(VarW) = DSmod
                 if( dum1 /= 0.0 )then
                   Gamma(VarU) = dum1
                   Gamma(VarV) = dum1
                   Gamma(VarW) = dum1
                 endif
                 if( Echo ) write(*,*)'+ velocity components: MinMod, blend ',dum1 
               case default
                 Scheme(VarU) = DSud
                 Scheme(VarV) = DSud
                 Scheme(VarW) = DSud
                 Gamma(VarU)  = 0.0
                 Gamma(VarV)  = 0.0
                 Gamma(VarW)  = 0.0
                 if( Echo ) write(*,*)'+ velocity components: default (UD)'
             end select
           case ('keps')
             select case (key3(1:3))
               case ('ud ')
                 Scheme(VarTE) = DSud
                 Scheme(VarED) = DSud
                 if( Echo ) write(*,*)'+ turb. model: UD'
               case ('cd1')
                 Scheme(VarTE) = DScd1
                 Scheme(VarED) = DScd1
                 if( dum1 /= 0.0 )then
                   Gamma(VarTE) = dum1
                   Gamma(VarED) = dum1
                 endif
                 if( Echo ) write(*,*)'+ turb. model: CD, blend ',dum1
               case ('cd2')
                 Scheme(VarTE) = DScd2
                 Scheme(VarED) = DScd2
                 if( dum1 /= 0.0 )then
                   Gamma(VarTE) = dum1
                   Gamma(VarED) = dum1
                 endif
                 if( Echo ) write(*,*)'+ turb. model: CD2, blend ',dum1
               case ('cd ')
                 Scheme(VarTE) = DScd1
                 Scheme(VarED) = DScd1
                 if( dum1 /= 0.0 )then
                   Gamma(VarTE) = dum1
                   Gamma(VarED) = dum1
                 endif
                 if( Echo ) write(*,*)'+ turb. model: CD, blend ',dum1
               case ('lud')
                 Scheme(VarTE) = DSlud
                 Scheme(VarED) = DSlud
                 if( dum1 /= 0.0 )then
                   Gamma(VarTE) = dum1
                   Gamma(VarED) = dum1
                 endif
                 if( Echo ) write(*,*)'+ turb. model: LUD, blend ',dum1 
               case ('gam')
                 Scheme(VarTE) = DSgam
                 Scheme(VarED) = DSgam
                 if( dum1 /= 0.0 )then
                   Gamma(VarTE) = dum1
                   Gamma(VarED) = dum1
                 endif
                !if( dum2 /= 0.0 ) Beta_m = dum2
                 if( Echo ) write(*,*)'+ turb. model: Gamma, blend ',dum1!,dum2
               case ('min')
                 Scheme(VarTE) = DSmod
                 Scheme(VarED) = DSmod
                 if( dum1 /= 0.0 )then
                   Gamma(VarTE) = dum1
                   Gamma(VarED) = dum1
                 endif
                 if( Echo ) write(*,*)'+ turb. model: MinMod, blend ',dum1 
               case default
                 Scheme(VarTE) = DSud
                 Scheme(VarED) = DSud
                 Gamma(VarTE)  = 0.0
                 Gamma(VarED)  = 0.0
                 if( Echo ) write(*,*)'+ turb. model: default (UD)'
             end select
           case ('t')
             select case (key3(1:3))
               case ('ud ')
                 Scheme(VarT)  = DSud
                 Scheme(VarSC) = DSud
                 if( Echo ) write(*,*)'+ temperature/scalars: UD'
               case ('cd1')
                 Scheme(VarT)  = DScd1
                 Scheme(VarSC) = DScd1
                 if( dum1 /= 0.0 )then
                   Gamma(VarT)  = dum1
                   Gamma(VarSC) = dum1
                 endif
                 if( Echo ) write(*,*)'+ temperature/scalars: CD, blend ',dum1
               case ('cd2')
                 Scheme(VarT)  = DScd2
                 Scheme(VarSC) = DScd2
                 if( dum1 /= 0.0 )then
                   Gamma(VarT)  = dum1
                   Gamma(VarSC) = dum1
                 endif
                 if( Echo ) write(*,*)'+ temperature/scalars: CD2, blend ',dum1
               case ('cd ')
                 Scheme(VarT)  = DScd1
                 Scheme(VarSC) = DScd1
                 if( dum1 /= 0.0 )then
                   Gamma(VarT)  = dum1
                   Gamma(VarSC) = dum1
                 endif
                 if( Echo ) write(*,*)'+ temperature/scalars: CD, blend ',dum1
               case ('lud')
                 Scheme(VarT)  = DSlud
                 Scheme(VarSC) = DSlud
                 if( dum1 /= 0.0 )then
                   Gamma(VarT)  = dum1
                   Gamma(VarSC) = dum1
                 endif
                 if( Echo ) write(*,*)'+ temperature/scalars: LUD, blend ',dum1 
               case ('gam')
                 Scheme(VarT)  = DSgam
                 Scheme(VarSC) = DSgam
                 if( dum1 /= 0.0 )then
                   Gamma(VarT)  = dum1
                   Gamma(VarSC) = dum1
                 endif
                !if( dum2 /= 0.0 ) Beta_m = dum2
                 if( Echo ) write(*,*)'+ temperature/scalars: Gamma, blend ',dum1!,dum2
               case ('min')
                 Scheme(VarT)  = DSmod
                 Scheme(VarSC) = DSmod
                 if( dum1 /= 0.0 )then
                   Gamma(VarT)  = dum1
                   Gamma(VarSC) = dum1
                 endif
                 if( Echo ) write(*,*)'+ temperature/scalars: MinMod, blend ',dum1 
               case default
                 Scheme(VarT)  = DSud
                 Scheme(VarSC) = DSud
                 Gamma(VarT)   = 0.0
                 Gamma(VarSC)  = 0.0
                 if( Echo ) write(*,*)'+ temperature/scalars: default (UD)'
             end select
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
             messages = .true.
         end select
!
! special models, commands, switches
!
       case ('use')
         key2 = keys(2)
         call lowercase(key2)

         select case (key2(1:5))
           case ('fixab')
             i = keyi(3)
             if( i == 0 ) write(*,*)'WARNING ctid incorrect!'
             
             tmpu = keyr(4)
             tmpv = keyr(5)
             tmpw = keyr(6)
             tmpk = keyr(7)
             tmpe = keyr(8)
             
             write(*,'(1x,A,i3)')'Fixing ABL top layer, ctid =',i
             write(IOdbg,'(1x,A,i3)')'Fixing ABL top layer, ctid =',i
             
             UseFixABL = .true.
             IdFixABL  = i
             UFixABL   = tmpu
             VFixABL   = tmpv
             WFixABL   = tmpw
             TeFixABL  = tmpk
             EdFixABL  = tmpe
             
             write(*,*)'Using:',tmpu,tmpv,tmpw,tmpk,tmpe
             write(IOdbg,*)'Using:',tmpu,tmpv,tmpw,tmpk,tmpe
           case ('artif')
             key3 = keys(3)
             call lowercase(key3)
             select case (key3(1:3))
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ artificial compressibility off'
                 UseArtificialComp = .false.
               case ('on')     
                 if( Echo ) write(*,'(1x,A,i3)')'+ artificial compressibility on'
                 UseArtificialComp = .true.
               case ('   ')     
                 if( Echo ) write(*,'(1x,A,i3)')'+ artificial compressibility on'
                 UseArtificialComp = .true.
               case default
                 write(*,'(1x,A,A)')'+ use,artif: unknown key: ',key3(1:3)
             end select
           case ('vtk')
             if( Echo ) write(*,'(1x,A,i3)')'+ using VTK'
             UseVTK = .true.
           case ('tecpl')
             if( Echo ) write(*,'(1x,A,i3)')'+ using Tecplot'
             UseTECPLOT = .true.
           case ('patch')
             if( Echo ) write(*,'(1x,A,i3)')'+ using patches by. B. Tuinstra'
             UsePatches   = .true.  ! switch on
             
             if( keyi(3) > 0 )then
               UseScalars   = .true.  ! patches relies on scalars, so switch on
               SolveScalars = .true.  ! 
             endif
             
             ! from ScanControlfile:
             if( Echo ) write(*,'(1x,A,i3)')'+ max. scalar dimension set to',MaxSC
             write(IOdbg,'(1x,A,i3)')'+ scalars: MaxSC=',MaxSC              
             
           case ('gmv')
             if( Echo ) write(*,'(1x,A,i3)')'+ using GMV'
             UseGMV = .true.
           case ('lapac')
             if( Echo ) write(*,'(1x,A,i3)') &
                    '+ using LAPACK in GradientPhiLeastSquares'
             UseLapack = .true.
           case ('parti')
             if( Echo ) write(*,'(1x,A,i3)')'+ using particles'
             UseParticles = .true.
           case ('gauss')
             i = keyi(3)
             if( i == 0 ) i = 2
             Ngradient = i

             if( Echo ) write(*,'(1x,A,i3)')'+ using Gauss for gradient',&
                                            Ngradient
             GradAlg = GradGauss 
           case default
             write(*,*)'+ use: unknown key: ',key2(1:5)
         end select
!
!===== postprocessors ==================================================
!
       case ('opendx')
         key2 = keys(2)
         call lowercase(key2)
         key3 = keys(3)
         call lowercase(key3)

         if( Echo ) write(*,'(1x,A,i3)')'+ opendx...'
         if( Echo ) write(*,*)'+ opendx keyword with strings: >',&
                      key2(1:lens(key2)),'< and >',&
                      key3(1:lens(key3)),'<'

         select case (key2)
           case ('normals')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ writing DXnormals'
                 DXnormals      = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ not writing DXnormals'
                 DXnormals      = .false.
               case default
                 write(*,'(1x,A,i3)')'+ error in command opendx'
                 messages = .true.
             end select
           case ('centers')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ writing DXcenters'
                 DXcenters      = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ not writing DXcenters'
                 DXcenters      = .false.
               case default
                 write(*,'(1x,A,i3)')'+ error in command opendx'
                 messages = .true.
             end select
           case ('massflux')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ writing DXmassflux'
                 DXmassflux      = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ not writing DXmassflux'
                 DXmassflux      = .false.
               case default
                 if( Echo ) write(*,'(1x,A,i3)')'+ error in command opendx'
                 messages = .true.
             end select
           case ('debug')
             select case (key3(1:lens(key3)))
               case ('on')
                 if( Echo ) write(*,'(1x,A,i3)')'+ writing DXdebugdata'
                 DXdebugdata      = .true.
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ not writing DXdebugdata'
                 DXdebugdata      = .false.
               case default
                 write(*,'(1x,A,i3)')'+ error in command opendx'
                 messages = .true.
             end select
           case ('dump')
             DXdump = keyi(3)
             if( Echo ) write(*,'(1x,A,i3,A)')'+ call opendx every ',&
                        DXdump,' steps'
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
             messages = .true.
         end select

       case ('post')
         key2 = keys(2)        ! variable
         call lowercase(key2)
         key3 = keys(3)        ! cell or vertex
         call lowercase(key3)
         key4 = keys(4)        ! yes or no
         call lowercase(key3)

         if( Echo ) &
           write(*,'(1x,A,A,A)')'+ post keyword with strings = ',key2,key3,key4
         select case (key2(1:3))
           case ('u  ')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarU) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarU) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarU) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarU) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('v  ')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarV) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarV) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarV) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarV) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('w  ')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarW) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarW) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarW) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarW) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('p  ')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarP) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarP) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarP) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarP) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('k  ')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarTE) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarTE) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarTE) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarTE) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('eps')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarED) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarED) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarED) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarED) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('t  ')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarT) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarT) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarT) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarT) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('sca')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarSC) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarSC) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarSC) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarSC) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('den')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarDen) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarDen) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarDen) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarDen) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('vis')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarVis) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarVis) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarVis) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarVis) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
           case ('lvi')
             if( key3(1:1) == 'c' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostC(VarLVis) = .true.
               else if( key4(1:1) == 'n' )then
                 PostC(VarLVis) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif
             else if( key3(1:1) == 'v' )then
               if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
                 PostV(VarLVis) = .true.
               else if( key4(1:1) == 'n' )then
                 PostV(VarLVis) = .false.
               else
                 write(*,'(1x,A,i3)')'+ error in post command:',key2           
               endif               
             else
               write(*,'(1x,A,i3)')'+ error in post command:',key3
               messages = .true.
             endif
               
           case default
             write(*,'(1x,A,i3)')'+ error in post command'
             messages = .true.
         end select
!
!===== particles =======================================================
!
       case ('gene')
         write(*,'(1x,A,i3)')'+ generation ignored here'
         messages = .true.
       
       case ('part')
         if( .not. UseParticles )then
           write(*,*)'+ error use particles not set'                         
           messages = .true.
           cycle
         endif
         !
         ! two forms:
         !
         ! PARTicle,#,PROP,dens,diam,COUPled etc
         ! 
         ! PARTicle,#,INIT,icoor,x,y,z,u,v,w ; icoor = 1 (global cartesian)
         ! 
         ! followed by a generation command. example:
         !
         ! PART,   1,PROP,1000.0,100e-6    define 1 particle
         ! GENE 19 1,,       0.0, 10e-6    and generate 19 others 
         !
         
         ipart = keyi(2) 
         key2  = keys(3)         
         call lowercase(key2)
         
         select case (key2(1:4))

           case ('prop')
             r1 = keyr(4) ! dens     
             r2 = keyr(5) ! diam

             Particle(ipart)%dens0 = r1
             Particle(ipart)%diam0 = r2
             
             if( LGenerator() )then
               call GetGenerators(ngen,ideltas,rdeltas)

               if( (ipart+ngen*ideltas(1)) > Npart )then
                 write(*,'(1x,A,4(1x,i6))')'+ too many particles generated (1):', &
                                        ipart,ngen,ideltas(1),ipart+ngen*ideltas(1)
                 messages = .true.                
               else
                 do i=1,ngen
                   j = i*ideltas(1)
                   Particle(ipart+j)%dens0 = Particle(ipart)%dens0 + float(i)*rdeltas(3) 
                   Particle(ipart+j)%diam0 = Particle(ipart)%diam0 + float(i)*rdeltas(4) 
                 end do
               endif
             endif

           case ('init')
             i1 = keyi(4) ! coordinate system (still to be done)
             if( i1 <= 0 .or. i1 >= 3 )then
               write(*,'(1x,A,i3)')'+ coordinate system still 1 or missing'
               messages = .true.                
             endif

             if( i1 == 1 )then
               r1 = keyr( 5) ! x     
               r2 = keyr( 6) ! y  
               r3 = keyr( 7) ! z  
               r4 = keyr( 8) ! u  
               r5 = keyr( 9) ! v  
               r6 = keyr(10) ! w  
             else
               x0 = 0.5
               y0 = 0.5       ! simpel test hack
               z0 = 0.125
               u0 = 0.0
               v0 = 0.0       
               w0 = 0.0

               !r1 = x0 + keyr( 5)*cosd(keyr( 6))      
               !r2 = y0 + keyr( 5)*sind(keyr( 6))     
               !r3 = keyr( 7)
               !r4 =      keyr( 8)*cosd(keyr( 9))   
               !r5 =      keyr( 8)*sind(keyr( 9))   
               !r6 = keyr(10)  
               write(*,*)'C>',r1,r2,r3,r4,r5,r6             
             endif
             
             Particle(ipart)%x0(1) = r1
             Particle(ipart)%x0(2) = r2
             Particle(ipart)%x0(3) = r3

             Particle(ipart)%v0(1) = r4
             Particle(ipart)%v0(2) = r5
             Particle(ipart)%v0(3) = r6


             if( LGenerator() )then
               call GetGenerators(ngen,ideltas,rdeltas)

               if( (ipart+ngen*ideltas(1)) > Npart )then
                 write(*,'(1x,A,4(1x,i6))')'+ too many particles generated (1):', &
                                        ipart,ngen,ideltas(1),ipart+ngen*ideltas(1)
                 messages = .true.                
               else
                 if( i1 == 1 )then
                   do i=1,ngen
                     j = i*ideltas(1)
                     Particle(ipart+j)%x0(1) = Particle(ipart)%x0(1) + float(i)*rdeltas(4) 
                     Particle(ipart+j)%x0(2) = Particle(ipart)%x0(2) + float(i)*rdeltas(5) 
                     Particle(ipart+j)%x0(3) = Particle(ipart)%x0(3) + float(i)*rdeltas(6) 
                     Particle(ipart+j)%v0(1) = Particle(ipart)%v0(1) + float(i)*rdeltas(7) 
                     Particle(ipart+j)%v0(2) = Particle(ipart)%v0(2) + float(i)*rdeltas(8) 
                     Particle(ipart+j)%v0(3) = Particle(ipart)%v0(3) + float(i)*rdeltas(9) 
                   end do
                 else
                   do i=1,ngen
                     j = i*ideltas(1)
                     !Particle(ipart+j)%x0(1) = x0 + keyr(5)*cosd(keyr(6)+float(i)*rdeltas(5))
                     !Particle(ipart+j)%x0(2) = y0 + keyr(5)*sind(keyr(6)+float(i)*rdeltas(5))
                     !Particle(ipart+j)%x0(3) = keyr( 7) + float(i)*rdeltas(6) 
                     !Particle(ipart+j)%v0(1) =      keyr(7)*cosd(keyr(8)+float(i)*rdeltas(8))
                     !Particle(ipart+j)%v0(2) =      keyr(7)*sind(keyr(8)+float(i)*rdeltas(8))
                     !Particle(ipart+j)%v0(3) = keyr(10) + float(i)*rdeltas(9) 
                   end do
                 endif
                 
               endif
             endif
                      
           case default
             write(*,'(1x,A,i3)')'+ error in particle command'
             messages = .true.
         end select
!
!===== special print options ===========================================
!
       case ('print')
         key2 = keys(2)
         call lowercase(key2)
         key3 = keys(3)
         call lowercase(key3)
         
         select case (key2(1:4))
           case ('cell')
             PrintCellVar     = .true.
             if( key3(1:4) == 'user' )then

               PrintCellVarUser = .true.             

             elseif( key3(1:3) == 'all' )then

               IPrintCellVarStart = 1
               IPrintCellVarEnd   = NCel
               IPrintCellVarInc   = 1

             else
             
               istart = keyi(3)
               iend   = keyi(4)
               iinc   = keyi(5)

               if( iend < istart )then
                 write(*,'(1x,A)')'+ error in print cell command: end < start'
                 messages = .true.
               elseif( iinc <= 0 )then
                 write(*,'(1x,A)')'+ error in print cell command: invalid increment'
                 messages = .true.
               elseif( istart <= 0 )then
                 write(*,'(1x,A)')'+ error in print cell command: invalid start'
                 messages = .true.
               elseif( iend <= 0 )then
                 write(*,'(1x,A)')'+ error in print cell command: invalid end'
                 messages = .true.
               endif

               IPrintCellVarStart = istart
               IPrintCellVarEnd   = iend
               IPrintCellVarInc   = iinc

               write(*,'(1x,A,3(i6))')'Printing cells ',istart,iend,iinc

             endif
           case ('wall')
             PrintWallVar     = .true.
             if( key3(1:4) == 'user' )then

               PrintWallVarUser = .true.             

             elseif( key3(1:3) == 'all' )then

               IPrintWallVarStart = 1
               IPrintWallVarEnd   = NBnd
               IPrintWallVarInc   = 1

             else
             
               istart = keyi(3)
               iend   = keyi(4)
               iinc   = keyi(5)

               if( iend < istart )then
                 write(*,'(1x,A)')'+ error in print wall command: end < start'
                 messages = .true.
               elseif( iinc <= 0 )then
                 write(*,'(1x,A)')'+ error in print wall command: invalid increment'
                 messages = .true.
               elseif( istart <= 0 )then
                 write(*,'(1x,A)')'+ error in print wall command: invalid start'
                 messages = .true.
               elseif( iend <= 0 )then
                 write(*,'(1x,A)')'+ error in print wall command: invalid end'
                 messages = .true.
               endif

               IPrintWallVarStart = istart
               IPrintWallVarEnd   = iend
               IPrintWallVarInc   = iinc

               write(*,'(1x,A,3(i6))')'Printing Walls ',istart,iend,iinc

             endif
           case ('file')
             PrintFile = keys(3)
             write(*,'(1x,A,A)')'Printing to file ',PrintFile
           case default
             write(*,'(1x,A)')'+ error in print command'
             messages = .true.
         end select
!
!===== rarely used changed constants/switches ==========================
!
       case ('switch')
         key2 = keys(2)
         call lowercase(key2)
         !write(*,*)'switch> ',key2 
         select case (key2(1:4))
           case ('cour')
             r = keyr(3)
             if( r < 1e-6 .or. r > 2.0 )then
               write(*,*)'Invalid particle Courant number:',r   
               messages = .true.
             else    
               ParticleCourant = r                          
               write(*,*)'Particle Courant number:',r
               write(IOdbg,*)'Switch: Particle Courant number:',r 
             endif
           case default
             write(*,'(1x,A)')'+ error in switch command'
             messages = .true.
         end select
!
!===== variables / library =============================================
!
       case ('set')
         key2 = keys(2)
         call lowercase(key2)

         r1  = keyr(3)
         r2  = keyr(4)
        
         call setvariable(key2(1:8),r1,r2)

       case ('math')
         key2 = keys(2)
         call lowercase(key2)

         if( key2(1:3) == 'deg' )then
           write(*,*)'Mathmode set to degrees'
           MathMode = 0
         else if( key2(1:3) == 'rad' )then
           write(*,*)'Mathmode set to radians'
           MathMode = 1         
         else
         
           key3 = keys(3)
           call lowercase(key3)

           r1  = keyr(4)
           r2  = keyr(5)
           r3  = keyr(6)
           r4  = keyr(7)

           call mathvariable(key3,MathMode,r1,r2,r3,r4)

           r2  = 0.0        
           call setvariable(key2(1:8),r1,r2)

         endif
!
!===== default =========================================================
!
       case default
         write(*,*) '+ unknown/unexpected keyword: ',key(1:lens(key))
         write(*,*) '+ string: ',string
         messages = .true.

     end select
     cycle 

10   continue
     write(*,*)'+ error in: ',string(1:lens(string))
     messages = .true.
     cycle

11   continue
     write(*,*)'+ read error: something else expected'
     write(*,*)'+ trying to continue'
     messages = .true.

   end do
20 continue
   write(*,*)'End of input deck'
   close(IOinp)

   !
   ! some final checks
   !
   if( Noutlets == 0 )then
     write(*,*) 'No outlets'
   elseif( Noutlets == 1 )then
     write(*,*) 'One outlet'
     if( Split /= 1.0 )then
       write(*,*)'Error: Split not equal to 1.0'
       messages = .true.       
     endif     
   elseif( Noutlets > 1 )then
     write(*,*) 'Multiple outlets'
     if( abs(Split-1.0) > 1.e-6 )then
       write(*,*)'Error: Sum of splits greater than 1.0', Split
       write(*,*)' > > : ',abs(Split-1.0) 
       messages = .true.       
     else
       write(*,*)'Sum split ',Split     
     endif
   endif

   !
   ! fill in the Solve-array 
   !
   Solve(VarU)    = SolveU
   Solve(VarV)    = SolveV
   Solve(VarW)    = SolveW
   Solve(VarP)    = SolveP
   Solve(VarTE)   = SolveTurbEnergy 
   Solve(VarED)   = SolveTurbDiss
   Solve(VarT)    = SolveEnthalpy
   Solve(VarSC)   = SolveScalars
   Solve(VarDEN)  = .false.
   Solve(VarVis)  = SolveVisc
   Solve(VarLVis) = .false.
   Solve(VarCP)   = .false.

   if( SolveScalars )then
     do is=1,Nscal
       Solve(NVar+is) = .true.
     end do
   endif

   !
   ! nicer/cleaner
   !
   if( Scheme(VarU)  == DScd1 .and. Gamma(VarU)  <= 0.001 ) Scheme(VarU)  = DSud
   if( Scheme(VarV)  == DScd1 .and. Gamma(VarV)  <= 0.001 ) Scheme(VarV)  = DSud
   if( Scheme(VarW)  == DScd1 .and. Gamma(VarW)  <= 0.001 ) Scheme(VarW)  = DSud
   if( Scheme(VarTE) == DScd1 .and. Gamma(VarTE) <= 0.001 ) Scheme(VarTE) = DSud
   if( Scheme(VarED) == DScd1 .and. Gamma(VarED) <= 0.001 ) Scheme(VarED) = DSud
   if( Scheme(VarT)  == DScd1 .and. Gamma(VarT)  <= 0.001 ) Scheme(VarT)  = DSud
   if( Scheme(VarSC) == DScd1 .and. Gamma(VarSC) <= 0.001 ) Scheme(VarSC) = DSud

   if( Scheme(VarU)  == DSud ) Gamma(VarU)  = 0.0
   if( Scheme(VarV)  == DSud ) Gamma(VarV)  = 0.0
   if( Scheme(VarW)  == DSud ) Gamma(VarW)  = 0.0
   if( Scheme(VarTE) == DSud ) Gamma(VarTE) = 0.0
   if( Scheme(VarED) == DSud ) Gamma(VarED) = 0.0
   if( Scheme(VarT)  == DSud ) Gamma(VarT)  = 0.0
   if( Scheme(VarSC) == DSud ) Gamma(VarSC) = 0.0
   
   if( messages )then
     write(*,*)'***'
     write(*,*)'*** Errors and/or warnings in input deck.'
     write(*,*)'***'
     stop
   endif

end subroutine ReadControlFile
subroutine parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

   character(len=*)   :: string
   character(len=120) :: tmpstring
   character(len=8)   :: name, value

   character(len=16)  :: keys(MaxKeys), blank
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   logical            :: math, found
   character(len=1)   :: s,s1 
 
   keyu(1:MaxKeys) = .false.
   math = .false.
   
   !
   ! test for blank line
   !
   if( lens(string) == 0 )then
     !write(*,*)'+ blank line'
     Nkeys = 0
     return
   endif

   !
   ! change first occurance of = into ,
   !
   indx = index(string,'=')
   if( indx > 0 ) string(indx:indx) = ' '
   !
   ! first remove extra spaces
   !  
   do i=1,len(string)
     s = string(i:i)
     if( s .eq. ' ' )then
       do j=i+1,len(string)
         s1 = string(j:j)
         if( s1 /= ' ' ) exit
       end do
       if( j /= (i+1) )then
         tmpstring = string(j:len(string)) 
         string(i+1:) = tmpstring
       endif
     endif
   end do

   !
   ! remove trailing komma's if any
   !
   do while( .true. )
     ilast = lens(string)
     if( string(ilast:ilast) == ',' )then
       string(ilast:ilast) = ' '
     else
       exit
     endif 
   end do

   !
   ! check for combination ', ' and replace by ',' 
   !
   do while( .true. )
     indx = index(string,', ')
     if( indx > 0 )then
       tmpstring = string(indx+2:len(string)) 
       string(indx:) = ','//tmpstring
     else
       exit
     endif 
   end do
   !
   ! test for leading space
   !
   if( string(1:1) == ' ' )then
     tmpstring = string(2:len(string))
     string(1:) = tmpstring
   endif

   !
   ! replace spaces by komma's
   !
   ilen = lens(string)
   do i=1,ilen
     if( string(i:i) == ' ' ) string(i:i) =','
   end do
   Nkeys = icnt

   !
   ! the line is filled with komma's as separators
   !
   ! count keywords
   !
   ilen = lens(string)
   icnt = 1
   do i=1,ilen
     if( string(i:i) == ',' ) icnt = icnt + 1
   end do
   Nkeys = icnt
   !write(*,*)'p5: number of keys:',Nkeys
      
   !
   ! fill in the strings
   !
   blank = '                '
   do i=1,MaxKeys
     keys(i) = blank
     keyi(i) =  0
     keyr(i) = 0.0
   end do
   
   i0 = 1
   i1 = 1
   do i=1,Nkeys
     indx    = index(string(i0:ilen),',')
     i1      = indx - 1
     !write(*,*) 'ip:',i,i0,i1
     if( i1 > 0 .and. i < Nkeys )then
       keyu(i) = .true.
       keys(i) = string(i0:i0+i1-1)
       i0      = i0 + i1 + 1
     elseif( i1 < 0 .and. i == Nkeys )then
       keyu(i) = .true.
       keys(i) = string(i0:ilen)
     elseif( i1 == 0 )then
       ! placeholder
       i0 = i0 + 1
       !write(*,*)'px:',i,'- - -',i0,i1
       keyu(i) = .false.
     else
       write(*,*)'+ internal error parser'
     endif
     !write(*,*)'px:',i,keys(i),i0,i1
   end do
  
   do i=1,Nkeys
     if( keyu(i) )read(keys(i),*,err=1) keyi(i) 
 1   continue
     !write(*,*)'pi:',i,keyi(i)
   end do 

   do i=1,Nkeys
     if( keyu(i) )read(keys(i),*,err=2) keyr(i)
 2   continue
     !write(*,*)'pr:',i,keyr(i)
   end do 

   !
   ! check for variables and insert them from library
   ! (set by 'set' command)
   !
   !write(*,*)'>>>',string(1:lens(string))
   do i=1,Nkeys
     if( keys(i)(1:1) == '$' )then
       name    = keys(i)(2:9)
       call getvariable(name,r1,r2)
       keyu(i) = .true.
       keyr(i) = r1
       keyi(i) = int(r1)     
       !write(*,*)'ps:',i,keys(i),'>',r1
     endif
   end do 
    
   !
   ! check for simple math (*/+-)
   !
   icnt = 0
   do i=2,Nkeys
     if( keys(i)(1:4) == '-   ' ) icnt = icnt + 1
     if( keys(i)(1:4) == '+   ' ) icnt = icnt + 1
     if( keys(i)(1:4) == '/   ' ) icnt = icnt + 1
     if( keys(i)(1:4) == '*   ' ) icnt = icnt + 1
     if( keys(i)(1:4) == '**  ' ) icnt = icnt + 1
   end do
   
   !if( icnt > 0 ) write(*,*) 'Math mode',icnt
   if( icnt > 0 ) math = .true.
   
   do while( icnt > 0 )

     do i=2,Nkeys-1
       found = .false. 
       if( keys(i)(1:4) == '-   ' )then
         r1 = keyr(i-1) - keyr(i+1)
         i1 = keyi(i-1) - keyi(i+1)       
         found = .true.
       else if( keys(i)(1:4) == '+   ' )then
         r1 = keyr(i-1) + keyr(i+1)
         i1 = keyi(i-1) + keyi(i+1)
         found = .true.
       else if( keys(i)(1:4) == '*   ' )then
         r1 = keyr(i-1) * keyr(i+1)
         i1 = keyi(i-1) * keyi(i+1)
         found = .true.
       else if( keys(i)(1:4) == '/   ' )then
         if( keyr(i+1) /= 0.0 ) r1 = keyr(i-1) / keyr(i+1)
         if( keyi(i+1) /=  0  )i1 = keyi(i-1) / keyi(i+1)
         found = .true.
       else if( keys(i)(1:4) == '**  ' )then
         r1 = keyr(i-1) ** keyr(i+1)
         i1 = keyi(i-1) ** keyi(i+1)
         found = .true.
       endif

       if( found )then
           
         keyu(i-1) = .true.
         keyr(i-1) = r1
         keyi(i-1) = i1 
         keys(i-1) = '(var)'

         do j=i,Nkeys-2
           keyu(j) = keyu(j+2)
           keyr(j) = keyr(j+2)
           keyi(j) = keyi(j+2)
           keys(j) = keys(j+2)
         end do

         keyu(Nkeys-1:MaxKeys) = .false.
         keyr(Nkeys-1:MaxKeys) =  0.0
         keyi(Nkeys-1:MaxKeys) =   0
         keys(Nkeys-1:MaxKeys) =  blank

         Nkeys = Nkeys - 2
         icnt  = icnt - 1 
         exit  
       endif 
     end do
   end do 

   !if( math )then
   !  write(*,*)'String:',string(1:lens(string))
   !  write(*,*)'New Nkeys:',Nkeys
   !  write(*,*)'keys:',keys(1:Nkeys)(1:4)
   !  write(*,*)'ints:',keyi(1:Nkeys)
   !  write(*,*)'rals:',keyr(1:Nkeys)
   !  write(*,*)'used:',keyu(1:Nkeys)
   !endif


   !write(*,*)'=== parser'

   
end subroutine parser
subroutine getkeyword(string,key,ikey)

   character(len=*) :: string, key
   character(len=1) :: s
   
   do i=1,len(string)
     s = string(i:i)
     if( s .eq. ' ' .or. s .eq. ',' .or. s .eq. '=') exit
     key(i:i) = string(i:i)
   end do
   ikey = i
      
end subroutine getkeyword
subroutine check_minmax(variable,varmin,varmax)

   if( variable < varmin )then
     write(*,*) '+ error: invalid entry'
     variable = varmin
   endif

   if( variable > varmax )then
     write(*,*) '+ error: invalid entry'
     variable = varmax
   endif

end subroutine check_minmax
logical function LGenerator()

   use constants
   character(len=120) string

   read(IOinp,'(A)',end=1) string   

   if( string(1:4) == 'gene' )then   
     !write(*,*)'Generator'
     backspace(IOinp)
     LGenerator = .true.
   else
     !write(*,*)'No generator'
     backspace(IOinp)
     LGenerator = .false.
   endif
   
   return
   
 1 continue
   write(*,*)'+ read error: premature eof (1)'
   stop
   
end function LGenerator
subroutine GetGenerators(ngen,ideltas,rdeltas)

   use constants
   
   character(len=120) string

   integer, parameter :: MaxKeys = 20
   character(len=16)  :: keys(MaxKeys) 
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   integer            :: ideltas(MaxKeys)
   real               :: rdeltas(MaxKeys)

   read(IOinp,'(A)',end=1) string    

   if( string(1:4) == 'gene' )then   
     !write(*,*)'Generator (2)'

     indx = index(string,'#')
     if( indx >= 2 )then               ! strip trailing comments
       do i=indx,len(string)
         string(i:i) = ' '
       end do
     endif

     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
 
     ngen    = keyi(2)
     if( Nkeys > 2 )then
       ideltas(1:Nkeys-2) = keyi(3:Nkeys) 
       rdeltas(1:Nkeys-2) = keyr(3:Nkeys) 
     else
       ideltas =  0 
       rdeltas = 0.0
     endif
   else
     write(*,*)'Internal error'
   endif

   return
   
 1 continue
   write(*,*)'+ read error: premature eof (2)'
   stop
   
end subroutine GetGenerators
subroutine dinlibrary

  use constants, only: IOdbg

  integer, parameter :: MaxDinVar = 80
  integer, save      :: NDinVar   =  0

  character(len=8)   :: variable
  logical            :: found

  character(len=8),dimension(MaxDinVar), save :: name
  real, dimension(MaxDinVar), save            :: rval1, rval2

  !=====================================================================
  entry setupvariables
  
  NDinVar = 0

  name(1)  = 'pi      '
  rval1(1) = 3.14159265
  rval2(1) = 0.0

  name(2)  = 'e       '
  rval1(2) = 2.71828183
  rval2(2) = 0.0

  do i=3,MaxDinVar
    name(i)  = '        '
    rval1(i) = 0.0
    rval2(i) = 0.0
  end do

  write(IOdbg,*)'DinLibrary variables set up:',MaxDinVar

  return

  !=====================================================================
  entry setvariable(variable,r1,r2)

  write(*,*)'set ',variable,r1,r2

  call lowercase(variable)

  do i=1,NDinVar
    if( name(i) == variable )then
      rval1(i) = r1
      rval2(i) = r2
      
      write(IOdbg,*)'DinLibrary modified:',variable,r1,r2      
      return
    end if
  end do
  
  if( NDinVar == MaxDinVar )then
    write(*,*)'+ error maximum variables reached'
    return
  endif
 
  
  NDinVar = NDinVar + 1
  name(NDinVar) = variable
  rval1(NDinVar) = r1
  rval2(NDinVar) = r2

  write(IOdbg,*)'DinLibrary stored:',variable,r1,r2
   
  return
  
  !=====================================================================
  entry getvariable(variable,r1,r2)

  found = .false.

  call lowercase(variable)

  do i=1,MaxDinVar
    if( variable == name(i) )then
      r1 = rval1(i)
      r2 = rval2(i)
      found = .true.
      exit
    endif
  end do

  if( .not. found )then
    write(*,*) '+ error: variable not found/set: ',variable
    r1 = 0.0
    r2 = 0.0
    stop
  endif 

  return
  
end subroutine dinlibrary
subroutine mathvariable(funct,mode,arg1,arg2,arg3,arg4)

   character(len=*), intent(IN)    :: funct
   integer, intent(IN)             :: mode

   real, intent(INOUT)             :: arg1
   real, intent(IN)                :: arg2
   real, intent(IN)                :: arg3
   real, intent(IN)                :: arg4
   
   real, parameter                 :: PI = 3.14159265
   real, parameter                 :: E  = 2.71828183 

   select case (funct) 

     case ('cos')
       if( mode == 0 )then
         arg1 = cos(arg1 * PI/180.)
       else if( mode == 1 )then
         arg1 = cos(arg1)
       endif
     case ('sin')
       if( mode == 0 )then
         arg1 = sin(arg1 * PI/180.)
       else if( mode == 1 )then
         arg1 = sin(arg1)
       endif
     case ('tan')
       if( mode == 0 )then
         arg1 = tan(arg1 * PI/180.)
       else if( mode == 1 )then
         arg1 = tan(arg1)
       endif
     case ('abs')
       arg1 = abs(arg1)
     case ('exp')
       arg1 = exp(arg1)
     case ('ln')
       arg1 = log(arg1)
     case ('log10')
       arg1 = log10(arg1)
     case ('sqrt')
       arg1 = sqrt(arg1)

     case default

       write(*,*)'Unknown math subcommand:',funct,mode
       arg1 = 0.0

   end select

end subroutine mathvariable
subroutine PrettyPrint(IO)

   use constants
   use geometry
   use variables
   use scalars

   use PatchesModule

   logical            :: warning
   integer            :: datum(8)
   character (len=12) :: clock(3)

   call date_and_time(clock(1),clock(2),clock(3),datum)

   write(IO,'(/)')
   write(IO,'(1x,A)')    '-------------------------------------------------------'
   write(IO,'(1x,A,A)')  'Case                 : ',Casename(1:lens(Casename))
   write(IO,'(1x,A,A)')  'Title                : ',title(1:lens(title))
   write(IO,'(1x,A,i2.2,A1,i2.2,A1,i4,A2,i2.2,A1,i2.2,A1)')   &
       'Date and time        : ',datum(3),'/',datum(2),'/',datum(1), &
       ' (',datum(5),':',datum(6),')'  
   write(IO,'(1x,A,i10)') 'Number of cells      : ',Ncel
   write(IO,'(1x,A,i10)') 'Number of vertices   : ',Nvrt
   write(IO,'(1x,A,i10)') 'Number of boundaries : ',Nbnd
   write(IO,'(1x,A,i10)') 'Number of faces      : ',Nfac
   write(IO,'(1x,A,i10)') 'Number of regions    : ',Nreg
   write(IO,'(1x,A,i10)') 'Number of scalars    : ',Nscal
   
   write(IO,'(1x,A)')    '-------------------------------------------------------'
   write(IO,'(1x,A,i10)') 'Number of iterations : ',Niter
   if( Transient )then
     if( Euler )then
       write(IO,'(1x,A,1pe10.3,A)') 'Transient run, dt    : ',dt,' (Euler)'
     else
       write(IO,'(1x,A,1pe10.3,A)') 'Transient run, dt    : ',dt,' (Quadratic)'
       write(IO,'(1x,A,1pe10.3)')   'Euler blending factor: ',GammaTime
     endif
   else
     write(IO,'(1x,A,1pe10.3)')   'Steady state run'
   endif
   if( UsePatches )then
   write(IO,'(1x,A,1pe10.3,A)') '*** Using patches'
   endif
   write(IO,'(1x,A,1pe10.3)')   'Target residual tol. : ',ResMax
   write(IO,'(1x,A,i10)')       'Monitor cell         : ',Imoni
   write(IO,'(1x,A,i10)')       'Pressure ref. cell   : ',IPref
   write(IO,'(1x,A,1pe10.3,A)') 'Reference pressure   : ',Pref,' Pa'
   if( SolveEnthalpy )then
   write(IO,'(1x,A,i10)')       'Temperature ref. cell: ',IPref
   write(IO,'(1x,A,1pe10.3,A)') 'Reference Temperature: ',Tref,' K'
   endif
   write(IO,'(1x,A,1pe10.3,A)') 'Reference density    : ',DensRef,' kg/m3'
   write(IO,'(1x,A,1pe10.3,A)') 'Molecular viscosity  : ',VisLam,' Pa s'
   if( SolveEnthalpy )then
   write(IO,'(1x,A,1pe10.3,A)') 'Prandtl number       : ',Prandtl,' -'
   write(IO,'(1x,A,1pe10.3,A)') 'Specific heat coef.  : ',CpStd,' J/kgK'
   write(IO,'(1x,A,1pe10.3,A)') 'Conductivity         : ',Lambda,' W/mK'
   endif
   write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Gravity vector       : ',Gravity,' '
   write(IO,'(1x,A)') '-------------------------------------------------------'
   write(IO,'(1x,A)') 'Var  SCV sch   blend      urf        rtol       guess'
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'U  : ',&
     SolveU ,PostC(VarU) ,PostV(VarU) ,DSch(Scheme(VarU)),Gamma(VarU) ,&
     URF(VarU) ,RTOL(VarU) ,Guess(VarU) 
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'V  : ',&
     SolveV ,PostC(VarV) ,PostV(VarV) ,DSch(Scheme(VarV)),Gamma(VarV) ,&
     URF(VarV) ,RTOL(VarV) ,Guess(VarV) 
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'W  : ',&
     SolveW ,PostC(VarW) ,PostV(VarW) ,DSch(Scheme(VarW)),Gamma(VarW) ,&
     URF(VarW) ,RTOL(VarW) ,Guess(VarW) 
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'P  : ',&
     SolveP ,PostC(VarP) ,PostV(VarP) ,'   ',Gamma(VarP) ,&
     URF(VarP) ,RTOL(VarP) ,Guess(VarP) 
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'TE : ',&
     SolveTurbEnergy,PostC(VarTE),PostV(VarTE),DSch(Scheme(VarTE)),Gamma(VarTE),&
     URF(VarTE),RTOL(VarTE),Guess(VarTE)
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'ED : ',&
     SolveTurbDiss,PostC(VarED),PostV(VarED),DSch(Scheme(VarED)),Gamma(VarED),&
     URF(VarED),RTOL(VarED),Guess(VarED)
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'T  : ',&
     SolveEnthalpy,PostC(VarT) ,PostV(VarT) ,DSch(Scheme(VarT)),Gamma(VarT) ,&
     URF(VarT) ,RTOL(VarT) ,Guess(VarT) 
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'SC : ',&
     SolveScalars,PostC(VarSC),PostV(VarSC),DSch(Scheme(VarSC)),Gamma(VarSC),&
     URF(VarSC),RTOL(VarSC),Guess(VarSC) 
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'DEN: ',&
     SolveDensity,PostC(VarDEN),PostV(VarDEN)
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'VIS: ',&
     SolveVISC,PostC(VarVIS),PostV(VarVIS)
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'CP : ',&
     SolveCP,PostC(VarCP),PostV(VarCP)

   if( Transient )then
     warning = .false.
     if( URF(VarU)  <= 0.99 ) warning = .true.
     if( URF(VarTE) <= 0.99 ) warning = .true.
     if( URF(VarT)  <= 0.99 ) warning = .true.
     if( URF(VarSC) <= 0.99 ) warning = .true.
     if( warning ) write(IO,'(/,1x,A)') &
        '  * Transient run with underrelaxtion factors < 1 '
   endif     

   write(IO,'(1x,A)')    '-------------------------------------------------------'
   write(IO,'(1x,A)')    '--- boundary conditions -------------------------------'
   write(IO,'(1x,A)')    '-------------------------------------------------------'
   write(IO,'(1x,A,i3)') 'Number of regions',NReg
   do i=0,NReg
     if( reg(i)%typ > 0 )then
       write(IO,'(1x,A,i3,A)')      '-  r e g i o n ',i,' -'
       write(IO,'(1x,A,A,i3)')      'Type         :  ',Region(reg(i)%typ)
       if( reg(i)%typ == RWall )then
         if( reg(i)%noslip )then
           write(IO,'(1x,A,3(1pe10.3,1x),A)') 'U, V, W      : ',reg(i)%uvw 
           write(IO,'(1x,A,3(1pe10.3,1x),A)') 'E log        : ',reg(i)%elog
           if( Reg(i)%adiab )then
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Adiabatic      ' 
           elseif( Reg(i)%flux )then
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'T flux, R    : ',reg(i)%T,reg(i)%R           
           else
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'T fixed, R   : ',reg(i)%T,reg(i)%R
           endif
         else
           write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Free slip wall'
         endif
       elseif( reg(i)%typ == RInlet )then
         write(IO,'(1x,A,3(1pe10.3,1x),L1)')'U, V, W      : ',reg(i)%uvw,reg(i)%user
         write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Density, T   : ',reg(i)%den,reg(i)%T
         write(IO,'(1x,A,3(1pe10.3,1x),A)') 'k, e         : ',reg(i)%k,reg(i)%e
       elseif( reg(i)%typ == ROutlet )then
         write(IO,'(1x,A,3(1pe10.3,1x),L1)')'Split        : ',reg(i)%splvl
       else
       endif
     endif
   end do

   if( UsePatches )then
     write(IO,'(1x,A)')    '-------------------------------------------------------'
     write(IO,'(1x,A)')    '--- patches -------------------------------------------'
     write(IO,'(1x,A)')    '-------------------------------------------------------'
     write(IO,'(1x,A,i3)') 'Number of patches',NPatches
   
     do i=1,NPatches
       write(IO,'(1x,A,i3,A)')      '- p a t c h ',i,' -'
       write(IO,'(1x,A,A,A)')       'Name         : ',Patch(i)%sName
       write(IO,'(1x,A,i3)')        'Type         : ',Patch(i)%itype
       write(IO,'(1x,A,A,i3)')      'Variable     : ',Patch(i)%sVar,Patch(i)%iVar
       write(IO,'(1x,A,i3)')              'Cell ID      : ',Patch(i)%iCellID
       write(IO,'(1x,A,3(1pe10.3,1x),A)') 'x,y,z low    : ',Patch(i)%xyz1
       write(IO,'(1x,A,3(1pe10.3,1x),A)') 'x,y,z high   : ',Patch(i)%xyz2
       write(IO,'(1x,A,i3)')        'Data length  : ',Patch(i)%iDatLen
       do j=1,Patch(i)%iDatLen
         write(IO,'(1x,A,i3,A,1pe10.3)')  '   data ',j,'  : ',Patch(i)%dat(j)
       end do
     end do
   endif

   if( SolveScalars )then
     write(IO,'(1x,A)')    '-------------------------------------------------------'
     write(IO,'(1x,A)')    '--- scalars -------------------------------------------'
     write(IO,'(1x,A)')    '-------------------------------------------------------'
     write(IO,'(1x,A,i3)') 'Number of scalars',NScal
     do i=1,NScal
       write(IO,'(1x,A,i3,A)')      '- s c a l a r ',i,' -'
       write(IO,'(1x,A,A,A)')       'Name         :  ',ScProp(i)%name,ScProp(i)%unit
       write(IO,'(1x,A,3L1)')       'Act/Sol/Sav  :  ',ScProp(i)%active,ScProp(i)%solve,ScProp(i)%save
       write(IO,'(1x,A,1pe10.3,A)') 'Lam. Prandtl : ',ScProp(i)%PrL,' -'
       write(IO,'(1x,A,1pe10.3,A)') 'Turb. Prandtl: ',ScProp(i)%PrT,' -'
     end do
     
     do i=1,NReg
       write(IO,'(1x,A,i3,A)')      '- r e g i o n ',i,' -'
       do j=1,NScal
         write(IO,'(i3,1x,A,3(1pe10.3,1x))')  j,ScProp(j)%name,ScReg(i,j)%fraction,ScReg(i,j)%value
       end do
     end do
   endif
   write(IO,'(1x,A)')    '-------------------------------------------------------'
   
end subroutine PrettyPrint
subroutine ReadCase

   use constants
   use geometry
   use variables

   logical :: Exists = .false.

   write(*,'(1x,A,f5.3)') 'This is dolfyn version ',float(Version)*0.001
   write(*,'(1x,A)') 'Copyright(C) 2002-2007 Cyclone Fluid Dynamics BV'
   write(*,'(1x,A)') 'NL-5583 XM, Waalre, The Netherlands'
   write(*,'(1x,A)') 'see http://www.cyclone.nl and http://www.dolfyn.net'
   write(*,'(1x,A)') 'Using Sparsekit2 by Yousef Saad'
   write(*,'(1x,A)') '(C) 2005, the Regents of the University of Minnesota'
   write(*,'(1x,A)') 'Modules with patches (C) 2004-2007 by B. Tuinstra'
   write(*,'(1x,A)') 'see http://www.home.zonnet.nl/bouke_1/dolfyn'
   write(*,'(1x,A)') 'Tecplot interface (C) 2006 by S.B. Kuang'
   write(*,'(1x,A)') 'VTK interface upated (C) 2007 by J. Jacobs'
   write(*,'(/)')  

   inquire(file='dolfyn.cfg',exist=exists)
   if( exists )then
     call openfile(IOcfg,'dolfyn','.cfg','FORMATTED', &
                         'SEQUENTIAL','OLD',debug)
     write(*,*) 'Reading case from dolfyn.cfg'
     read(IOcfg,'(A)') Casename 
     close(IOcfg)                         
   else
     write(*,*) 'Enter case:'
     read(*,'(A)') Casename 
   endif
   write(*,*) 'Using case: ',Casename(1:lens(Casename))

   !
   ! check if STOP-file is present; if so remove it
   ! (note: IOcfg used again)
   !
   inquire(file='STOP',exist=Exists)
   if( Exists )then
     call openfile(IOcfg,'STOP','','FORMATTED', &
                         'SEQUENTIAL','OLD',debug)
     write(*,*)'Removing STOP'
     close(IOcfg,status='DELETE') 
   endif

end subroutine ReadCase


