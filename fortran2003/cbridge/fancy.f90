SUBROUTINE fancy(a, b,c,d,e) bind ( C, name="fancy" )
use iso_c_binding
IMPLICIT NONE
INTEGER, INTENT(IN) :: a
INTEGER, INTENT(OUT) :: b
character (kind=c_char, len=1), dimension (10), intent (in) :: c
character (len=10) :: regular_string
REAL*8, INTENT(IN) :: d
LOGICAL*1, INTENT(IN) :: e

! String conversion
integer :: i
regular_string = " "
loop_string: do i=1, 10
  if ( c (i) == c_null_char ) then
     exit loop_string
  else
     regular_string (i:i) = c (i)
  end if
end do loop_string
!end

b = a*a

print *, c
print *, d
print *, e

END SUBROUTINE fancy
