subroutine CalcIndVel(NT,ntTerm,NBE,NB,NE,Px,Py,Pz,Vx,Vy,Vz)

    use wallsoln
    use blade
    use wake, only: InductionFlag

    ! Calculate wall and wake induced velocity (including bound vorticity component)

    real :: Px, Py, Pz, Vx, Vy, Vz, dUdX
    integer :: nt, ntTerm, nbe, nb, ne

    real :: Point(3), dVel(3), Vel(3)

    ! Calc wake induced velocity at wake locations
    if (InductionFlag==1) then
       CALL BladeIndVel(NT,ntTerm,NBE,NB,NE,Px,Py,Pz,Vx,Vy,Vz,dUdX,0,0)

       ! Calculate wall induced velocities at wake locations
       Point=[Px,Py,Pz]
       Call WallIndVel(Point,dVel)
       Vx=Vx+dVel(1)
       Vy=Vy+dVel(2)
       Vz=Vz+dVel(3)
    else
       Vx=0.0
       Vy=0.0
       Vz=0.0
    endif

    return
end subroutine CalcIndVel

