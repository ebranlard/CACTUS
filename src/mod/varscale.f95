module varscale

    ! Atmos props and scale parameters

    real :: rho             ! density
    real :: vis             ! viscosity
    real :: rpm             ! turbine rotation rate
    real :: Rmax            ! radius
    integer :: Incompr  ! Incompressibility flag, =1 (incompressible), =0 (compressible)

    real :: R_in_m   ! radius [m]
    real :: U_in_mps ! nominal wind velocity [m/s]
end module varscale
