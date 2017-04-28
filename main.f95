program main
    use spectrum_method
    implicit none
    real(8) a(10),a_next(10)
    real(8) u(11)
    integer :: n = 11, m = 3,i, k
    real(8) x(11), t1
    
    open(20, file='output.d')
    
    call init_am(a)
    do k = 1, (m+1)
        call make_next_a(a,a_next)
        write(*,*) a(:)
        do i = 1, n
            x(i) = (i-1)*0.1d0
            call make_u(x(i), a, u(i))
        enddo
        a(:) = a_next(:)
    enddo
    t1 = m*T
    do i = 1, 11
        write(20,*) x(i), u(i), exact_f(x(i),t1)
        write(*,*) x(i), u(i), exact_f(x(i),t1)
    enddo
   
    close(20)
contains
    real(8) function exact_f(x,t)
        real(8), intent(in) :: x,t
        exact_f = exp(-pi*pi*t)*sin(pi*x)
        return 
    end function exact_f
end program main