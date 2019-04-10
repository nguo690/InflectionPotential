module wh_interpolate
    implicit none

    integer :: nx
    double precision, allocatable, dimension(:) :: x
    double precision, allocatable, dimension(:) :: y
    double precision, allocatable, dimension(:) :: ypp

    contains

    ! This function creates a new cubic_spline_task using the intel mkl function dfdnewtask1d, 
    ! and then initialises it with 
    subroutine initialise_cubic_spline(x_in,y_in,nx_in)
        implicit none

        ! Inputs
        ! ------
        ! number of interpolation points
        integer,intent(in) :: nx_in

        ! x positions of the interpolation points
        double precision, intent(in), dimension(nx) :: x_in

        ! y positions of the interpolation points
        double precision, intent(in), dimension(nx) :: y_in

        ! Local variables
        ! ---------------
        ! Matrix to solve
        double precision,dimension(nx_in,nx_in) :: M
        double precision,dimension(nx_in-1) :: sub_diag
        double precision,dimension(nx_in) :: diag
        double precision,dimension(nx_in-1) :: sup_diag
        double precision,dimension(nx_in-2) :: sup_diag_2
        ! Vector to solvve
        double precision,dimension(nx_in) :: b

        integer, dimension(nx_in) :: ipiv
        integer :: info

        if(allocated(x))   deallocate(x)
        if(allocated(y))   deallocate(y)
        if(allocated(ypp)) deallocate(ypp)

        nx = nx_in
        allocate(x(nx),y(nx),ypp(nx))
        ! Save the x and y values for later
        x = x_in
        y = y_in


        ! Initialise the sub diagonal, diagonal and super diagonal components of
        ! the matrix M, and the vector b
        ! The equation M * ypp = b is the set of equations for the second
        ! derivatives of the splines
        !
        ! We use 'Natural boundary conditions' i.e. y''=0 at the ends of the
        ! spline
        sub_diag = (/ (x(2:nx-1)-x(1:nx-2))/6, 0d0 /)
        diag = (/ 1d0, (x(3:nx)-x(1:nx-2))/3 , 1d0 /) 
        sup_diag = (/ 0d0, (x(3:nx)-x(2:nx-1))/6 /)

        ! Note the zeroes on either side for natural boundary conditions
        b = (/ 0d0 , (y(3:nx)-y(2:nx-1))/(x(3:nx)-x(2:nx-1)) - (y(2:nx-1)-y(1:nx-2))/(x(2:nx-1)-x(1:nx-2)) , 0d0 /)

        ! Solve the equation for the second derivatives
        ! (1) Compute LU factorisation
        call dgttrf(nx,sub_diag,diag,sup_diag,sup_diag_2,ipiv,info)
        ! (2) Solve the equations for ypp
        ypp=b
        call dgttrs('N',nx,1,sub_diag,diag,sup_diag,sup_diag_2,ipiv,ypp,nx,info)

        ! We now have our second derivatives and can compute splines


    end subroutine

    ! This function recieves the input value x, and the derivative order to
    ! compute, and outputs either y, y' or y'' depending on the value of dorder.
    function evaluate_cubic_spline(x_in, dorder) result(y_out)
        implicit none

        ! Inputs
        ! ------
        ! Inputs to spline 
        double precision, intent(in),target :: x_in
        ! Derivative order to compute
        integer,intent(in) :: dorder
        ! Output of spline
        double precision :: y_out

        ! Local variables
        ! ---------------
        double precision :: A,B,C,D

        integer :: j


        ! Find the bin that it is in
        if(x_in<=x(1)) then
            j = 1
            A = 1
            B = 0
            select case(dorder)
            case(0)
                y_out = y(1) + ( x_in - x(1) ) * spline_1()
            case(1)
                y_out = spline_1()
            case(2)
                y_out = 0
            case default
                y_out = 0d0
                write(*,'("Unrecognised dorder ", I5,  " in evaluate_cubic_spline")') dorder
            end select
            return
        else if(x_in>=x(nx)) then
            j = nx-1
            A = 0
            B = 1
            select case(dorder)
            case(0)
                y_out = y(nx) + ( x_in - x(nx) ) * spline_1()
            case(1)
                y_out = spline_1()
            case(2)
                y_out = 0
            case default
                y_out = 0d0
                write(*,'("Unrecognised dorder ", I5,  " in evaluate_cubic_spline")') dorder
            end select
            return
        else

            do j=1,nx-1

                if(x_in <= x(j+1) ) then

                    A = ( x(j+1)-x_in )  / ( x(j+1)-x(j) )
                    B = 1-A
                    C = (A**3 - A) * (x(j+1)-x(j))**2 / 6
                    D = (B**3 - B) * (x(j+1)-x(j))**2 / 6

                    select case(dorder)
                    case(0)
                        y_out = spline_0()
                    case(1)
                        y_out = spline_1()
                    case(2)
                        y_out = spline_2()
                    case default
                        y_out = 0d0
                        write(*,'("Unrecognised dorder ", I5,  " in evaluate_cubic_spline")') dorder
                    end select
                    return

                end if

            end do

        end if

        contains
        function spline_0
            implicit none
            double precision :: spline_0
            spline_0 = A*y(j) + B*y(j+1) + C*ypp(j) + D*ypp(j+1) 
        end function spline_0 

        function spline_1
            implicit none
            double precision :: spline_1
            spline_1 = ( y(j+1)-y(j) ) / ( x(j+1)-x(j) ) + ( (3*B**2-1)/6 * ypp(j+1) - (3*A**2-1)/6 * ypp(j) ) * (x(j+1)-x(j))  
        end function spline_1 

        function spline_2
            implicit none
            double precision :: spline_2
            spline_2 = A * ypp(j) + B * ypp(j+1)
        end function spline_2 


    end function evaluate_cubic_spline


end module wh_interpolate
