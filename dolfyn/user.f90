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
! This is user.f90; the central place for user subroutines
!
! Contains conrtibutions by:
! - B. Tuinstra, see www.home.zonnet.nl/bouke_1/dolfyn
!
subroutine UserInitialField(IP,X,Y,Z,Volume,&
                            Pref,Tref,Iteration,Time, &
                            U,V,W,P,TE,ED,T,DEN)
!=======================================================
!
!   This routine set the special user-defined initial  
!   conditions
!
!   All values are preset based on the already
!   provided input. Only deviations need to be
!   specified.
!   return values like: U,V,W,TE,ED,T,DEN,Scalars
!
!   in case other data is need, use a module
!
   use geometry   

   integer, intent(IN)  :: IP, Iteration
   
   real, intent(IN)     :: X,Y,Z, Pref,Tref, Volume, Time 
   
   real, intent (INOUT) :: U,V,W,P,TE,ED,T,DEN

   logical :: done = .false. 
   save done

   if( cell(ip)%ctid == 1 .or. cell(ip)%ctid == 6 )then
     U  = 3.907E-01
     V  = 9.205E-01
     W  = 0.0
     TE = 1.5E-04 
     ED = 3.0E-05
   else if( cell(ip)%ctid == 2 .or. cell(ip)%ctid == 3 )then
     U  = 0.0
     V  = 0.0
     W  = 0.1
     TE = 10.E-03  
     ED =  1.E-03
   else
     U  = 0.0
     V  = 0.0
     W  = 0.1
     TE = 10.E-03  
     ED =  1.E-03
   
   endif

   !
   ! wind
   !
   !call SetAbl(Z,Uabl,TEabl,EDabl,10.00,0.03,0.03)
   !
   !U  = Uabl
   !V  = 0.0
   !W  = 0.0
   !TE = TEabl
   !ED = EDabl
   !
   !if( cell(ip)%ctid == 4  .and. .not. done )then
   !  done = .true.
   !  
   !  write(*,*)'Bovenin: ',u,te,ed
   !
   !endif

   !
   ! stagnation flow
   !
   !U =  X
   !V = -Y
   !W = 0.0
   
end subroutine UserInitialField
subroutine UserInlet(IRegion,Ibnd,IP,X,Y,Z,Area,&
                     Pref,Tref,Iteration,Time,  &
                     U,V,W,TE,ED,T,DEN,         &
                     NAME1,NAME2)
!=======================================================
!
!   This routine set the special user-defined inlet  
!   boundary conditions
!
!   All values are preset based on the already
!   provided input. Only deviations need to be
!   specified.
!   return values like: U,V,W,TE,ED,T,DEN,Scalars
!
!   This example provides a special test cases
!
!
   use geometry, only  : Face, Bnd
   real, dimension(3) :: normaal, vradiaal

   integer, intent(IN)  :: IRegion, Ibnd, IP, Iteration
   
   real, intent(IN)     :: X,Y,Z, Pref,Tref, Area, Time 
   
   character(len=4), intent(IN)  :: NAME1, NAME2

   real, intent (INOUT) :: U,V,W,TE,ED,T,DEN

   if( IRegion == 1 .and. NAME1 == 'piet' )then
     !
     ! piet's stukje:
     !
     yloc = sqrt( X**2 + Y**2 )
     U  = yloc *( 1.0 - yloc )

   elseif( IRegion == 1 .and. NAME1 == 'abl ' )then
     !
     ! atmosferic boundary layer:
     !
     !call SetAbl(Z,Uabl,TEabl,EDabl,10.0,0.03,0.03)
     !
     !U  =  Uabl * ( 0.04728)
     !V  =  Uabl * (-0.99888)
     !W  =  0.0
     !TE = TEabl
     !ED = EDabl
     
   elseif( NAME1 == 'leon' .and. NAME2 == 'step' )then
     !
     ! Leonard's step profile test
     !
     if( Y < 0.16666666 )then
       T = 0.0
     else
       T = 1.0   
     endif
     
     T = T + 293.0
     
   elseif( NAME1 == 'leon' .and. NAME2 == 'sin2' )then
     !
     ! Leonard's sin^2 profile: sin^2(3pi(y-1/6)) 
     !
     T = 0.0
     if( 0.16666666 <= Y .and. Y < 0.5 ) &
       T = ( sin( 3.*3.14159265*(Y-0.16666666) ) )**2

     T = T + 293.0

   elseif( NAME1 == 'leon' .and. NAME2 == 'semi' )then
     !
     ! Leonard's semi-ellipse profile: sqrt(1-((x-1/3)/(1/6))^2) 
     !
     T = 0.0
     if( 0.16666666 <= Y .and. Y < 0.5 ) &
       T = sqrt(1.0-((Y-0.33333333)/(0.16666666))**2)

     T = T + 293.0
     
   else

     !call SetAbl(Z,Uabl,TEabl,EDabl,10.0,0.03,0.03)
     !
     !U  = Uabl
     !V  = 0.0
     !W  = 0.0
     !TE = TEabl
     !ED = EDabl

     !
     ! *** *** ***
     !
     ! read from table
     !
     !call ReadTable(x,y,z,ut,vt,wt)
     !
     !U = Ut
     !V = Vt
     !W = Wt
     !
     !TE = 1.5*(0.1*0.1) * (U*U+V*V+W*W)    ! intensiteit 10%
     !ED = 0.164 * TE**1.5 / 0.15           ! 10% van 1.5  
     !
     ! *** *** ***
     !
     ! radiaal instroom van een fan
     !
     !rpm   = -4000.0
     !omega = (2.*3.1415927)/60. * rpm 
     !r     = sqrt(X**2+Y**2)
     !if( x >= 0.0 )then
     !  alpha = atan(y/x)
     !else if( x < 0.0 .and. y >= 0.0 )then
     !  alpha = 3.1415927 - atan(y/(-x))
     !else 
     !  alpha = -3.1415927 + atan(y/x)
     !endif
     !
     !Vrad = 2.195445523171375 * 1.2058549321546825 ! m3
     !Vrad = 2.195445523171375 * 1.118998181413545  ! m2
     !Vtan = omega * r
     !
     !u1 =  Vrad*cos(alpha)
     !v1 =  Vrad*sin(alpha)
     !
     !u2 = -Vtan*sin(alpha)
     !v2 =  Vtan*cos(alpha)
     !
     !u  = u1 + u2
     !v  = v1 + v2
     !w  = 0.0

     !
     ! k = 3/2( i |U| )**2
     !
     ! e = Cmu**(3/4) k**(3/2) / l
     !
     ! *** *** ***
     !
     ! instroom normaal op het inlaat vlakje (een bol)
     ! 
     Vinl = 0.0001 * 18.902646 * 0.3333 * 2.579311581645777
     normaal = Face(Bnd(ibnd)%face)%n      ! positief naar buiten
     call normalise(normaal)
     Vradiaal = -Vinl * normaal 
     
     U  = Vradiaal(1)
     V  = Vradiaal(2)
     W  = Vradiaal(3)
     
     TE = 1.5*(0.001*0.001) * (U*U+V*V+W*W)   ! intensiteit 0.1%
     ED = 0.164 * TE**1.5 / 0.001             ! lengte schaal 1 m 
 
   endif

   !write(*,'(8(1x,f8.3))') alpha*180/3.1415927,u,v,sqrt(u**2+v**2)    

end subroutine UserInlet
subroutine SetAbl(Z,Uabl,TEabl,EDabl,Umet,Zomet,Zoloc)

   ! voorbeeld grenslaag
   !
   ! meteostation data
   !
   z0s  = Zomet
   Uoo  = Umet
   
   uss  = 0.419 * Uoo /log( (10.+z0s)/z0s )

   f1   = (1000.0 + z0s) / z0s
   Qref = z0s * uss / 0.419 * ( f1*log(f1) - f1)
   !
   ! doel
   !
   z0r  = Zoloc
   
   f2   = (1000.0 + z0r) / z0r
   usr  = Qref * 0.419 / z0r / ( f2*log(f2) - f2 )
   
   fs   = usr / 0.419
   fmu  = usr*usr / sqrt( 0.09 )
   f3   = usr**3 / 0.419
   
   Uabl  = fs * log( (Z + z0r)/z0r )
   TEabl = fmu
   EDabl = f3 / ( Z + z0r )
  
end subroutine SetAbl
subroutine ReadTable(x,y,z,u,v,w)

   real, intent(IN)  :: x,y,z 
   real, intent(OUT) :: u,v,w 

   parameter(N=12)

   character(len=32)          :: string
   real, save, dimension(N)   :: Xin, Yin, Zin
   real, save, dimension(N,N) :: Uin, Vin, Win

   logical, save :: Done = .false.

   if( .not. done )then
     !
     ! read file
     !
     write(*,*)'Open tabel'
     open(91,file='tabel.dat')
     
     read(91,*) string                ! skip kopregel

     icnt = 0
     do j=1,N
       do i=1,N
         read(91,*) Xtmp, Ztmp, Win(i,j), Vin(i,j), Uin(i,j)
         
         Xin(i)   = Xtmp - 0.751560   ! translatie
         Zin(j) = -Ztmp - 0.739853
         
         Uin(i,j) = -Uin(i,j)         !
         Vin(i,j) = -Vin(i,j)         ! omdraaien
         Win(i,j) = -Win(i,j)         !
         
         icnt = icnt + 1
       end do
     end do
     write(*,*)'Lines read:',icnt
     
     write(*,*)'Data tabel:'
     write(*,'('' X:'',12(1x,f6.3))') Xin(1:N)
     write(*,'('' Z:'',12(1x,f6.3))') Zin(1:N)
     write(*,*) minval(Uin),' < U < ',maxval(Uin)
     write(*,*) minval(Vin),' < V < ',maxval(Vin)
     write(*,*) minval(Win),' < W < ',maxval(Win)
     close(91)
     
     do j=1,N
       write(*,'(12(''('',f5.2,'','',f5.2,'','',f5.2,'')''))') &
         (Uin(i,j),Vin(i,j),Win(i,j),i=1,N)
     end do
     
     Done = .true.
   endif
   
   !
   ! find x(i) and x(i+1)
   !
   ix = -1
   do i=1,N-1
     if( X > Xin(i) .and. X <= Xin(i+1) )then
       ix = i
       goto 10
     endif
   end do

10 jz = -1
   do j=1,N-1
     if( Z < Zin(j) .and. Z >= Zin(j+1) )then
       jz = j
       goto 20
      endif
   end do

20 if( ix == -1 ) write(*,*)'FOUT: X niet gevonden'
   if( jz == -1 ) write(*,*)'FOUT: Z niet gevonden'

   xtmp = X - Xin(ix)
   ztmp = Z - Zin(jz)
   
   u1   = Uin(ix,jz)
   u2   = Uin(ix+1,jz)
   u3   = Uin(ix,jz+1)
   u4   = Uin(ix+1,jz+1)

   v1   = Vin(ix,jz)
   v2   = Vin(ix+1,jz)
   v3   = Vin(ix,jz+1)
   v4   = Vin(ix+1,jz+1)

   w1   = Win(ix,jz)
   w2   = Win(ix+1,jz)
   w3   = Win(ix,jz+1)
   w4   = Win(ix+1,jz+1)

   x1   = Xin(ix)
   x2   = Xin(ix+1)
   x3   = Xin(ix)
   x4   = Xin(ix+1)
   
   z1   = Zin(jz)
   z2   = Zin(jz)
   z3   = Zin(jz+1)
   z4   = Zin(jz+1)
   
   ugrad1 = (u2-u1)/(x2-x1)
   ugrad2 = (u4-u3)/(x4-x3)
   ugradx = 0.5*( ugrad1 + ugrad2 )
   
   ugrad1 = (u3-u1)/(z3-z1)
   ugrad2 = (u4-u2)/(z4-z2)
   ugradz = 0.5*( ugrad1 + ugrad2 )

   vgrad1 = (v2-v1)/(x2-x1)
   vgrad2 = (v4-v3)/(x4-x3)
   vgradx = 0.5*( vgrad1 + vgrad2 )
   
   vgrad1 = (v3-v1)/(z3-z1)
   vgrad2 = (v4-v2)/(z4-z2)
   vgradz = 0.5*( vgrad1 + vgrad2 )

   wgrad1 = (w2-w1)/(x2-x1)
   wgrad2 = (w4-w3)/(x4-x3)
   wgradx = 0.5*( wgrad1 + wgrad2 )
   
   wgrad1 = (w3-w1)/(z3-z1)
   wgrad2 = (w4-w2)/(z4-z2)
   wgradz = 0.5*( wgrad1 + wgrad2 )

   u = u1 + xtmp*ugradx + ztmp*ugradz 
   v = v1 + xtmp*vgradx + ztmp*vgradz 
   w = w1 + xtmp*wgradx + ztmp*wgradz 

end subroutine ReadTable
subroutine UserPrintCell(IOunit)

   use constants
   use geometry
   use variables

   !write(IOunit,'(/,'' User cell data'')')
   
   do i=IPrintCellVarStart,min(IPrintCellVarEnd,NCel),IPrintCellVarInc

     write(IOunit,*) Cell(i)%x(2),T(i)-Tref-20.0
     
   end do

end subroutine UserPrintCell
subroutine UserPrintWall(IOunit)

   use constants
   use geometry
   use variables

   real, dimension(20)    :: ab, pb, xw, yw, zw
   integer, dimension(20) :: ic

   write(IOunit,'(/,'' User wall data'')')

   ic =  0 
   ab = 0.0
   pb = 0.0
   xw = 0.0
   yw = 0.0
   zw = 0.0
   
   do i=1,Nbnd
     ir = Bnd(i)%rid
     if( ir >= 9 )then
       ib = Bnd(i)%face
       pf = P(Ncel+i)
       
       xw(ir) = xw(ir) + Face(ib)%x(1)
       yw(ir) = yw(ir) + Face(ib)%x(2)
       zw(ir) = zw(ir) + Face(ib)%x(3)
       
       ic(ir) = ic(ir) + 1
       ab(ir) = ab(ir) + Face(ib)%area
       pb(ir) = pb(ir) + Face(ib)%area * pf
     endif
   end do

   write(IOunit,*)'Wand'
   do ir=8,9
     c = ic(ir)
     write(IOunit,*) ir,ic(ir),xw(ir)/c,yw(ir)/c,zw(ir)/c, &
                ' a:',ab(ir),' p:',pb(ir)/ab(ir)
   end do

   write(IOunit,*)'Done wall'

end subroutine UserPrintWall

