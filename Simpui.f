        module Simpui_module
            implicit none
            real(8) :: pi = 3.1415926
        contains
            subroutine integrate_simpui(x, Func, s)
                real(8), intent(in) :: x(:)
                real(8), allocatable :: y(:)
                real(8), intent(out) :: s
                real(8) Func
                integer :: m, err, i, n = 21
                
                m = size(x)
                allocate(y(m))
                do i = 1, n
                    y(i) = Func(x(i))
                enddo
                write(*,*) y(:)
                    
                call Simpui(x, y, n, s, err)
            end subroutine integrate_simpui
            
            subroutine Simpui(x, y, n, s, err)
                integer, intent(in) :: n
                real(8), intent(in) :: x(:), y(:)
                real(8), intent(out) :: s
                integer, intent(out) :: err
                real(8) tmp(3)
                real(8) h, hpd, hmd, d
                integer i, nm1
                if(n .lt. 3 .or. Mod(n,2) .eq. 0) then
                    err = 999
                    return
                else
                    nm1 = n -1
                    do i = 1, nm1
                        if(x(i+1) .le. x(i)) then
                            err = 999
                            return
                        endif
                    enddo
                    err = 0
                    s = 0.0d0
                    do i = 2, nm1, 2
                        tmp = (/x(i-1),x(i),x(i+1)/)
                        h = 0.5d0*(tmp(3)-tmp(1))
                        hpd = 1.0d0/(tmp(2)-tmp(1))
                        hmd = 1.0d0/(tmp(3)-tmp(2))
                        d = tmp(2) - tmp(1) - h
                        s = s + (1.0d0+2.0d0*d*hpd)*y(i-1)*h/3.0d0
                        s = s + 2.0d0*h*(hpd+hmd)*y(i)*h/3.0d0
                        s = s + (1.0d0-2.0d0*d*hmd)*y(i+1)*h/3.0d0
                    enddo
                endif
            end subroutine Simpui
        end module Simpui_module