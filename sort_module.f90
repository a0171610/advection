module sort_module
  implicit none
  private

  public :: bubblesort

  contains
  subroutine bubblesort(N, array)
    !sikinote, 2016/08/08
    implicit none
    integer, intent(in) :: N
    real(8), intent(inout) :: array(1 : N)
    integer :: i,j
    real(8) :: t
     
    do i = 1,N-1
       do j = i+1,N
          if(array(i) .gt. array(j))then
             t = array(i)
             array(i) = array(j)
             array(j) = t
          end if
       end do
    end do
  
    return
  end subroutine bubblesort
  end module sort_module