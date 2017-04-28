module spectrum_method
    use Simpui_module
    implicit none
    integer :: j = 10
    real(8) :: T = 0.02
contains    
    real(8) function make_am0(m)
        integer, intent(in) :: m
        real(8) :: x(21), y(21), s
        integer :: n = 21, err, i
        do i = 1, n
            x(i) = (i-1)*0.05
            y(i) = 2.0d0*sin(pi*x(i))*sin(x(i)*m*pi)
            !y(i) = 2.0d0*x(i)*(1.0d0-x(i))*sin(x(i)*m*pi)
        enddo
        call Simpui(x, y, n, s, err)
        make_am0 = s
        return 
    end function make_am0
    
    subroutine init_am(a)
        real(8), intent(out) :: a(:)
        integer m
        
        do m = 1, j
            a(m) = make_am0(m)
        enddo
    end subroutine init_am
    
    subroutine make_u(x, a, u)
        real(8), intent(in) :: x, a(:)
        real(8), intent(out) :: u
        integer i
        
        u = 0.0d0
        do i = 1, j
            u = u + a(i)*sin(i*pi*x)
        enddo
    end subroutine make_u
    
    subroutine make_next_a(a0, a_next)
        real(8), intent(in) :: a0(:)
        real(8), intent(out) :: a_next(:)
        real(8) cq(4), ckq(4), ck(4)
        real(8) q, ak, r, a
        integer m, i
        
        cq(:)=(/2.0d0,1.0d0,1.0d0,2.0d0/)
        ckq(:)=(/0.5d0,0.29289322d0,1.7071068d0,0.16666666d0/)
        ck(:)=(/0.5d0,0.29289322d0,1.7071068d0, 0.5d0/)
        do m = 1, j
!            a_next(m) = (1.0d0-T*(m*pi)*(m*pi))*a0(m)
            a = a0(m)
            do i = 1, 4
                ak = T*func(a, m)
                r = (ak-cq(i)*q)*ckq(i)
                a = a + r
                q = q + 3.0d0*r-ck(i)*ak
            enddo
            a_next(m) = a
        enddo
    end subroutine make_next_a

    real(8) function func(a,m)
        real(8), intent(in) :: a
        integer, intent(in) :: m
        func = -(m*pi)*(m*pi)*a
        return 
    end function func
end module spectrum_method