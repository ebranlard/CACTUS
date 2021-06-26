subroutine input(ErrFlag)

    use parameters

    use dystl
    use element
    use blade
    use wake
    use strut
    use varscale
    use shear
    use iecgust
    use airfoil
    use configr
    use pidef
    use vortex
    use wakedata
    use fielddata
    use time
    use wallgeom
    use wallsoln
    use wallsystem
    use regtest
    use output
    use tower
    use probesystem
    use fnames


    integer, parameter :: InBufferNumSectionTables = 1000
    integer, parameter :: InBufferNumWL = 1000
    integer, parameter :: MaxReadLine = 1000
    integer, parameter :: MaxTempAOA = 1000

    integer i, ii, jj, kk
    integer ErrFlag
    logical NotDone, NotBlank
    character(MaxReadLine) :: ReadLine
    character(MaxReadLine) :: GeomFilePath    ! path to geometry input file
    integer :: CI, EOF
    real :: temp, temp1(MaxTempAOA,4)


    ! Temp buffers

    character(MaxReadLine) :: AFDPath(InBufferNumSectionTables) ! Airfoil section data path
    character(MaxReadLine) :: WallMeshPath ! Wall mesh file path
    character(MaxReadLine) :: ProbeSpecPath ! Wall mesh file path


    integer :: WLI(InBufferNumWL)     ! wake line index buffer

    ! Namelist input file declaration
    NAMELIST/ConfigInputs/RegTFlag,GPFlag,WPFlag,FSFlag,nr,convrg,nti,iut,iWall,ivtxcor,VCRFB,VCRFT,VCRFS,vCutOffRad,ifc,convrgf,nric,ntif,iutf,ixterm,xstop, &
        Incompr,DSFlag,LBDynStallTp,PRFlag,AddMassFlag,InductionFlag,PrescribedGamma,&
        k1pos,k1neg,GPGridSF,GPGridExtent,FSGridSF,TSFilFlag,ntsf

    NAMELIST/CaseInputs/jbtitle,GeomFilePath,RPM,Ut,nSect,AFDPath, &
        hAG,dFS,rho,vis,tempr,hBLRef,slex,Cdpar,CTExcrM, &
        WLI,Igust,gustamp,gusttime,gustX0, &
        Itower,tower_Npts,tower_x,tower_ybot,tower_ytop,tower_D,tower_CD, &
        WallMeshPath

    NAMELIST/ConfigOutputs/OutputPath,BladeElemOutFlag,DynStallOutFlag,WallOutFlag,DiagOutFlag, &
        WakeElemOutFlag,WakeVTKOutFlag,WakeElemOutIntervalTimesteps,WakeElemOutStartTimestep,WakeElemOutEndTimestep, &
        FieldOutFlag,FieldOutIntervalTimesteps,FieldOutStartTimestep,FieldOutEndTimestep, &
        nxgrid,nygrid,nzgrid,xgridL,ygridL,zgridL,xgridU,ygridU,zgridU, &
        WallOutIntervalTimesteps,WallOutStartTimestep,WallOutEndTimestep, &
        ProbeFlag,ProbeOutIntervalTimesteps,ProbeOutStartTimestep,ProbeOutEndTimestep,ProbeSpecPath


    ! Default ConfigInputs
    RegTFlag   = 0
    GPFlag     = 0
    WPFlag   = 0 ! Flag for reading in arbitrary wall geometries
    FSFlag     = 0
    TSFilFlag  = 0
    ntsf       = 3
    nSect      = 1
    ifc        = 0
    nr         = 10
    convrg     = -1
    convrgf    = -1
    nti        = 20
    ntif       = -1
    ivtxcor    = 1
    vCutOffRad = 1e-7
    ixterm     = 0
    xstop      = 5.0
    iut        = 0
    iWall      = 0
    iutf       = 0
    nric       =-1
    VCRFB      = 1.0
    VCRFT      = 1.0
    VCRFS      = 1.0
    WLI(:)     = 0 ! all set to 0
    hAG        = 0.0
    dFS        = 0.0
    Incompr    = 0
    Cdpar      = 0.0
    CTExcrM    = 0.0
    DSFlag     = 1
    LBDynStallTp = 1.7
    PRFlag     = 1
    AddMassFlag= 1
    InductionFlag= 1
    PrescribedGamma= 0.0
    k1pos      = 1.0
    k1neg      = 0.5
    GPGridSF   = 1.0
    GPGridExtent = 10.0
    FSGridSF   = 1.0
    Igust      = 0
    gustamp    = 0.0
    gusttime   = 0.0
    gustX0     = 0.0
    Itower     = 0
    tower_Npts = 10
    tower_x    = 0.0
    tower_ybot = 0.0
    tower_ytop = 0.0
    tower_D    = 0.05
    tower_CD   = 1.0

    ! output options
    OutputPath         = 'output'
    BladeElemOutFlag      = 0
    DynStallOutFlag      = 0
    DiagOutFlag        = 0
    WakeElemOutFlag = 0
    WakeVTKOutFlag = 0
    FieldOutFlag    = 0
    WallOutFlag        = 0
    ProbeFlag          = 0

    ! field output default parameters
    nxgrid =    1
    nygrid =  100
    nzgrid =  100

    xgridL =  0.0
    xgridU =  0.0
    ygridL = -2.0     ! default grid is a ([-2.0,2.0],[-2.0,2.0]) 100x100 x-normal grid at x=0.0
    ygridU =  2.0
    zgridL = -2.0
    zgridU =  2.0

    ! Wake Output Frequency
    WakeElemOutIntervalTimesteps =  5       ! write wake element data every 5 timesteps
    WakeElemOutStartTimestep     =  1       ! write wake element data starting at first timestep
    WakeElemOutEndTimestep       = -1       ! stop writing wake element data at the last timestep

    ! Field Output Frequency
    FieldOutIntervalTimesteps    =  5       ! write field data every 5 timesteps
    FieldOutStartTimestep        =  1       ! write field data starting at first timestep
    FieldOutEndTimestep          = -1       ! stop writing field data at the last timestep

    ! Wall Output Frequency
    WallOutIntervalTimesteps        =  5       ! write wall data every 5 timesteps
    WallOutStartTimestep            =  1       ! write wall data starting at first timestep
    WallOutEndTimestep              = -1       ! stop writing wall data at the last timestep

    ! Probe Output Frequency
    ProbeOutIntervalTimesteps       =  1
    ProbeOutStartTimestep           =  1
    ProbeOutEndTimestep             = -1
    ! Namelist input
    read(4, nml=ConfigInputs)
    read(4, nml=CaseInputs)
    read(4, nml=ConfigOutputs)

    ! Write the geometry file to stdout (for bookkeeping)
    write(*,*) 'Geometry file'
    write(*,*) '--------------------------'
    Call file_to_stdout(GeomFilePath)
    write(*,*) ''

    ! Read geometry file
    Call InputGeom(GeomFilePath)

    ! Set array bounds based on inputs
    ! Geometry
    MaxBlades = nb
    MaxSegPerBlade = nbe
    MaxSegEndPerBlade = MaxSegPerBlade+1
    MaxSegEnds = MaxSegEndPerBlade*MaxBlades
    MaxSeg = MaxSegPerBlade*MaxBlades


    MaxStruts = NStrut
    ! Airfoil Data
    MaxAirfoilSect = nSect
    MaxReVals = 20
    MaxAOAVals = 1000
    ! Wake advancement
    MaxRevs = nr
    MaxTimeStepPerRev = nti
    MaxWakeNodes = MaxRevs * MaxTimeStepPerRev
    ! Non-linear convergence iteration
    MaxNLIters = 10
    ! Outputs
    MaxTimeSteps = MaxRevs * MaxTimeStepPerRev
    ! Wake outputs
    NWakeInd=0
    NotDone=.TRUE.
    i=1
    do while (NotDone .AND. i<=MaxSeg)
        if (WLI(i) == 0) then
            NWakeInd=NWakeInd+1
        else
            NotDone=.FALSE.
        end if
        i=i+1
    end do

    ! Array construction
    CALL blade_cns(MaxSegEnds)
    CALL wake_cns(MaxWakeNodes,MaxSegEnds)
    CALL element_cns(MaxSegEnds,MaxSegEndPerBlade)
    CALL airfoil_cns(MaxAOAVals,MaxReVals,MaxAirfoilSect)
    CALL wakedata_cns()
    CALL fielddata_cns()
    CALL dystl_cns(MaxAirfoilSect,MaxReVals,MaxSegEnds)
    CALL output_cns(MaxSeg,MaxBlades,MaxStruts,DSFlag)

    ! Write from buffer...

   !  WakeLineInd(1:NWakeInd)=WLI(1:NWakeInd)
    WakeLineInd = (/ (I, I = 1, MaxSeg) /)
    ! write(*,'(A)') WakeLineInd

    ! Set ground plane location for wall solution
    GPy=-hAG/Rmax

    ! Set depth and Froude number for free surface solution
    FSy=dFS/Rmax
    g=32.174  ! gravity, ft/s^2
    A=Rmax*(2.0*pi*rpm/60.0)**2/g   ! Accel ratio: w^2*R/g
    FnR=sqrt(A/ut**2)   ! Froude number based on turbine radius  FnR=Uinf/sqrt(g*R)

    ! Normalize ground shear inputs
    yref = hBLRef/Rmax  ! location of boundary layer edge (U/Uinf = 99% maybe) normalized to radius...
    ygc  = hAG/Rmax   ! Ground clearance normalized to radius

    ! Floor the tip speed ratio to the next lowest int
    ! and use this as the default update interval...
    if (iutf == 0) iutf = floor(ut)
    if (iut == 0) iut = floor(ut)
    if (iWall == 0) iWall = floor(ut)

    ! Default ntif to nti if nothing was input
    if (ntif .eq. -1) then
        ntif=nti
    end if

    ! Set number of RHS evaluations to average for the free surface calculation (should cover approx 1 revolution)
    NFSRHSAve=nti/iWall

    ne = (nbe+1)*nb ! Total number of blade segment ends (over all blades)

    ! Timestep filter setup
    KTF=1.0/real(ntsf)

    ! Read in wall geometry
    if (WPFlag == 1) then

        if (WPFlag == 1 .and. GPFlag == 1) then
            write(*,*) 'Error: WPFlag and GPFlag cannot both be set to 1!'
        end if

        write(*,'(A,A)') 'Reading walls from: ', WallMeshPath
        call read_p3d_walls(WallMeshPath)

        write(*,*) 'Summary of Wall Mesh'
        write(*,*) '--------------------------'
        write(*,*) 'Total walls:  ', Nwalls
        write(*,*) 'Total panels: ', NumWP_total

        do iw=1,NWalls
            write(*,*) 'Wall ', iw
            write(*,*) '    Dimensions:   ', Walls(iw)%NumWP1, ' x ', Walls(iw)%NumWP2
            write(*,*) '    Total Panels: ', Walls(iw)%NumWP
            write(*,*) ''
        end do

        write(*,*) ''

    end if

    ! Read in probe data
    if (ProbeFlag > 0) then
        ! Read probe data from file
        call read_probes(ProbeSpecPath)
    end if


    ! Airfoil Data Tables: Read CL, CD, CM vs AOA from data files
    ! Format Example:
    ! Title: AFTitle
    ! Thickness to chord ratio: 0.2
    ! Zero Lift AOA (deg): 0.0
    ! Reverse camber direction: 0
    !
    ! Reynolds Number: 1e6
    ! BV Dyn. Stall Model - Positive stall AOA (deg): 10
    ! BV Dyn. Stall Model - Negative stall AOA (deg): -10
    ! LB Dyn. Stall Model - Lift Coeff. Slope at Zero Lift AOA (per radian): 6.28
    ! LB Dyn. Stall Model - Positive Critical Lift Coeff.: 1.3
    ! LB Dyn. Stall Model - Negative Critical Lift Coeff.: -1.3
    ! AOA (deg) CL CD Cm
    ! ... ... ... ...
    !
    ! Reynolds Number: 5e6
    ! ...
    do kk = 1, nsect

        ! Open input file for this section
        open(15, file=AFDPath(kk))
        EOF=0

        ! Find title block
        NotDone=.TRUE.
        do while (NotDone)
            read(15,'(A)') ReadLine
            CI=index(ReadLine,':')
            if (CI>0) then
                NotDone=.FALSE.
            end if
        end do

        ! Read title and airfoil thickness
        if (len_trim(ReadLine)>CI) then
            aftitle(kk) = ReadLine(CI+1:len_trim(ReadLine))
        else
            aftitle(kk) = 'No Title'
        end if
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) tc(kk)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) alzer(kk)
        alzer(kk)=alzer(kk)*conrad
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) camb(kk)

        ! Reverse camber direction if desired
        if (camb(kk) == 1) then
            alzer(kk) = -alzer(kk)
        end if

        ! Find first Re block
        NotDone=.TRUE.
        do while (NotDone)
            read(15,'(A)',IOSTAT=EOF) ReadLine
            CI=index(ReadLine,':')
            if (CI>0 .OR. EOF<0) then
                NotDone=.FALSE.
            end if
        end do

        ! Read data for each Re value
        i=0
        do while (EOF >= 0  .AND. (i < MaxReVals))

            i=i+1
            ! Read Re and dyn. stall data
            read(ReadLine(index(ReadLine,':')+1:),*) TRE(i,kk)
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) alstlp(i,kk)
            alstlp(i,kk)=alstlp(i,kk)*conrad
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) alstln(i,kk)
            alstln(i,kk)=alstln(i,kk)*conrad
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) CLaData(i,kk)
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) CLCritPData(i,kk)
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) CLCritNData(i,kk)

            ! Reverse camber direction if desired
            if (camb(kk) == 1) then
                temp = alstlp(i,kk)
                alstlp(i,kk) = -alstln(i,kk)
                alstln(i,kk) = -temp
                temp = CLCritPData(i,kk)
                CLCritPData(i,kk) = -CLCritNData(i,kk)
                CLCritNData(i,kk) = -temp
            end if

            ! Read AOA data
            read(15,'(A)') ReadLine
            NotDone=.TRUE.
            ii=0
            do while (NotDone)
                read(15,'(A)',IOSTAT=EOF) ReadLine
                ! Check for carriage return (len_trim doesn't consider this a blank)
                NotBlank=.TRUE.
                if (len_trim(ReadLine)==0) then
                    NotBlank=.FALSE.
                else if (len_trim(ReadLine)==1) then
                    if (ichar(ReadLine(len_trim(ReadLine):len_trim(ReadLine))) == 13) then
                        NotBlank=.FALSE.
                    end if
                end if
                if (EOF>=0 .AND. NotBlank) then
                    if (ii == MaxAOAVals) then
                        write(6,'(A)') 'Max. allowed AOA values exceeded in airfoil data file: ', aftitle(kk)
                        ErrFlag=1
                        NotDone=.FALSE.
                    else
                        ii=ii+1
                        read(ReadLine,*) ta(ii,i,kk),tcl(ii,i,kk),tcd(ii,i,kk),tcm(ii,i,kk)
                    end if
                else
                    NotDone=.FALSE.
                end if
            end do
            ntbl(i,kk)=ii

            ! Check AOA limits
            if (ta(1,i,kk) > -180.0 .OR. ta(ntbl(i,kk),i,kk) < 180.0) then
                write(6,'(A)') 'AOA data needs to be +/-180 deg in airfoil data file: ', aftitle(kk)
                ErrFlag=1
            end if

            ! Reverse camber direction if desired
            if(camb(kk) == 1) then
                do ii = 1, ntbl(i,kk)
                    temp1(ii,1) = ta(ii,i,kk)
                    temp1(ii,2) = tcl(ii,i,kk)
                    temp1(ii,3) = tcd(ii,i,kk)
                    temp1(ii,4) = tcm(ii,i,kk)
                end do

                do ii = 1, ntbl(i,kk)
                    jj = ntbl(i,kk)-(ii-1)
                    ta(ii,i,kk) = -temp1(jj,1)
                    tcl(ii,i,kk) = -temp1(jj,2)
                    tcd(ii,i,kk) = temp1(jj,3)
                    tcm(ii,i,kk) = -temp1(jj,4)
                end do
            end if

            ! Find next Re block
            NotDone=.TRUE.
            if (EOF<0) then
                NotDone=.FALSE.
            end if
            do while (NotDone)
                read(15,'(A)',IOSTAT=EOF) ReadLine
                CI=index(ReadLine,':')
                if (CI>0 .OR. EOF<0) then
                    NotDone=.FALSE.
                end if
            end do

        end do
        ! Set number of Re vals for this section
        nRET(kk)=i

        ! Close input file for this section
        close(15)

        ! Check data
        if (i == 0) then
            write(6,'(A)') 'Error reading airfoil data file: ', aftitle(kk)
            ErrFlag=1
        end if
        if (EOF > 0) then
            write(6,'(A)') 'Warning: Max. allowed Re values exceeded in airfoil data file: ', aftitle(kk)
        end if

    end do

    ! ---
    R_in_m   = Rmax/3.280840
    U_in_mps = 2.0*pi*rpm/60.0*R_in_m/Ut ! U=Omega R/lambda
    PrescribedGamma =PrescribedGamma/(R_in_m*U_in_mps)

    return
601 format(' ','***airfoil section specified for blade segment ',i2,' is illegal. set to airfoil section 1***')
end subroutine input
