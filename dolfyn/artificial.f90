!
! Copyright 2004-2005 Thomas Geenen 
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
subroutine artificial_comp()
    
    use constants
    use geometry
    use variables
    use artificial_compressibility
    
    !implicit none
    !integer ::   i, j, k, in, ip
    double precision :: U_ref, V_ref, W_ref 
!   ************************************** 
!   section for artificial compressibility

    if( Debug > 3 ) write(*,*)'*** ArtificialComp'

    write(IOdbg,*) 'Using artificial compressibility'

    U_ref = maxval(U)
    V_ref = maxval(V)
    W_ref = maxval(W)
   
    do i=1,Ncel
        art_snd_speed(i) = sqrt(U(i)**2 + V(i)**2 + W(i)**2)
        art_snd_speed(i) =ac_beta*max(art_snd_speed(i),(sqrt(0.5*(U_ref**2))))
        U_P(1:3) = (/U(i),0.,0./)
        V_P(1:3) = (/0.,V(i),0./)
        W_P(1:3) = (/0.,0.,W(i)/)
        sum_VS(:) = 0.0
        sum_area(:) = 0.0
 
     do j=1,NFaces(i)
        k  = CFace(i,j)
        ip = Face(k)%cell1
        in = Face(k)%cell2
 
        if( in > 0 )then
          ! interne wand (omdat in =/= 0 is)
          if( ip == i )then
              if (U_P(1) * Face(k)%n(1) .gt. 0.0) then
                  sum_VS(1) = sum_VS(1) + dot_product( U_P , Face(k)%n )!*Face(i)%area
                  sum_area(1) = sum_VS(1)/U_P(1)
              endif
              if (V_P(2) * Face(k)%n(2) .gt. 0.0) then
                  sum_VS(2) = sum_VS(2) + dot_product( V_P , Face(k)%n )!*Face(i)%area
                  sum_area(2) = sum_VS(2)/V_P(2)
              endif
              if (W_P(3) * Face(k)%n(3) .gt. 0.0) then
                  sum_VS(3) = sum_VS(3) + dot_product( W_P , Face(k)%n )!*Face(i)%area
                  sum_area(3) = sum_VS(3)/W_P(3)
              endif
          else if( in == i )then
              if (U_P(1) * Face(k)%n(1) .lt. 0.0) then
              sum_VS(1) = sum_VS(1) - dot_product( U_P , Face(k)%n )!*Face(i)%area
              sum_area(1) = sum_VS(1)/U_P(1)
              endif
              if (V_P(2) * Face(k)%n(2) .lt. 0.0) then
                  sum_VS(2) = sum_VS(2) - dot_product( V_P , Face(k)%n )!*Face(i)%area
                  sum_area(2) = sum_VS(2)/V_P(2)
              endif
              if (W_P(3) * Face(k)%n(3) .lt. 0.0) then
                  sum_VS(3) = sum_VS(3) - dot_product( W_P , Face(k)%n )!*Face(i)%area
                  sum_area(3) = sum_VS(3)/W_P(3)
              endif

          endif
        else
            ! externe wand (omdat in == 0 is)
            if (U_P(1) * Face(k)%n(1) .gt. 0.0) then
                sum_VS(1) = sum_VS(1) + dot_product( U_P , Face(k)%n )!*Face(i)%area
                sum_area(1) = sum_VS(1)/U_P(1)
            endif
            if (V_P(2) * Face(k)%n(2) .gt. 0.0) then
                sum_VS(2) = sum_VS(2) + dot_product( V_P , Face(k)%n )!*Face(i)%area
                sum_area(2) = sum_VS(2)/V_P(2)
            endif
            if (W_P(3) * Face(k)%n(3) .gt. 0.0) then
                sum_VS(3) = sum_VS(3) + dot_product( W_P , Face(k)%n )!*Face(i)%area
                sum_area(3) = sum_VS(3)/W_P(3)
            endif


        endif
     end do
     delta_t(i) = cfl*Cell(i)%vol/(abs(sum_VS(1)) &
                  + abs(sum_VS(2)) + abs(sum_VS(3))&
                  + art_snd_speed(i) * (sqrt(sum_area(1)**2 + sum_area(2)**2 &
                  + sum_area(3)**2)))
     end do

    if( Debug > 3 ) write(*,*)'=== ArtificialComp'

end subroutine artificial_comp
