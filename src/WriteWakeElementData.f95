module FVW_VTK
    implicit none
    character(8), parameter :: RFMT='E17.8E3'
    character(8), parameter :: IFMT='I7'

   TYPE, PUBLIC :: FVW_VTK_Misc
      integer :: vtk_unit
      logical :: bFileOpen=.false.

      integer :: nData=0;
      integer :: nPoints=0;

      logical :: bBinary = .false.
      character(len=255) :: buffer

      ! Reference Frame
      real,dimension(3,3) :: T_g2b
      real,dimension(3)   :: PO_g
   END TYPE FVW_VTK_Misc

    character(1), parameter :: NL = char(10) ! New Line character

    interface vtk_dataset_structured_grid; module procedure &
            vtk_dataset_structured_grid_flat, &
            vtk_dataset_structured_grid_grid
    end interface

    interface vtk_point_data_vector; module procedure &
            vtk_point_data_vector_flat, &
            vtk_point_data_vector_grid2D,&
            vtk_point_data_vector_grid
    end interface
    interface vtk_point_data_scalar; module procedure &
            vtk_point_data_scalar_flat, &
            vtk_point_data_scalar_grid2D, &
            vtk_point_data_scalar_grid
    end interface
    interface vtk_cell_data_scalar; module procedure &
            vtk_cell_data_scalar_1d,&
            vtk_cell_data_scalar_2d 
    end interface

    public

contains
   subroutine WrVTK_Segments(filename, mvtk, SegPoints, SegConnct, SegGamma, SegEpsilon, bladeFrame) 
      character(len=*),intent(in)                 :: filename
      type(FVW_VTK_Misc),           intent(inout) :: mvtk       !< miscvars for VTK output
      real, dimension(:,:),      intent(in) :: SegPoints  !< 
      integer, dimension(:,:),  intent(in) :: SegConnct  !< 
      real,     dimension(:)  ,  intent(in) :: SegGamma   !< 
      real,     dimension(:)  ,  intent(in) :: SegEpsilon !< 
      logical,                      intent(in   ) :: bladeFrame !< Output in blade coordinate frame
      if ( vtk_new_ascii_file(filename,'Sgmt',mvtk) ) then
         call vtk_dataset_polydata(SegPoints(1:3,:),mvtk,bladeFrame)
         call vtk_lines(SegConnct(1:2,:)-1,mvtk) ! NOTE: VTK indexing at 0
         call vtk_cell_data_init(mvtk)
         call vtk_cell_data_scalar(SegGamma  ,'Gamma',mvtk)
         call vtk_cell_data_scalar(SegEpsilon,'Epsilon',mvtk)
   !       call vtk_cell_data_scalar(real(SegConnct(3,:), ReKi),'Age',mvtk)
         !call vtk_cell_data_scalar(real(SegConnct(4,:), ReKi),'Span',mvtk)
         call vtk_close_file(mvtk)
      endif
   end subroutine

   subroutine WrVTK_Lattice(filename, mvtk, LatticePoints, LatticeGamma, LatticeData3d, bladeFrame)
      character(len=*), intent(in)                         :: filename
      type(FVW_VTK_Misc),           intent(inout)          :: mvtk          !< miscvars for VTK output
      real, dimension(:,:,:), intent(in   )  :: LatticePoints  !< Array of points 3 x nSpan x nDepth
      real, dimension(:,:), intent(in  )             :: LatticeGamma  !< Array of            nSpan x nDepth
      real, dimension(:,:,:), intent(in  ), optional :: LatticeData3d !< Array of n x nSpan x nDepth KEEP ME
      logical,                      intent(in   )          :: bladeFrame    !< Output in blade coordinate frame
      !
      integer, dimension(:,:), allocatable :: Connectivity
      real, dimension(:,:), allocatable     :: Points

      CALL LatticeToPanlConnectivity(LatticePoints, Connectivity, Points)

      if ( vtk_new_ascii_file(filename,'',mvtk)) then
         call vtk_dataset_polydata(Points,mvtk,bladeFrame)
         call vtk_quad(Connectivity,mvtk)
         call vtk_cell_data_init(mvtk)
         call vtk_cell_data_scalar(LatticeGamma,'Gamma',mvtk)
         if (present(LatticeData3d)) then
            call vtk_point_data_init(mvtk)
            call vtk_point_data_vector(LatticeData3d,'Uconv',mvtk)
         endif
         call vtk_close_file(mvtk)
      endif

   end subroutine WrVTK_Lattice

   subroutine LatticeToPanlConnectivity(LatticePoints, Connectivity, Points)
      real, dimension(:,:,:), intent(in   )  :: LatticePoints  !< Array of points 3 x nSpan x nDepth
      integer, dimension(:,:), allocatable :: Connectivity
      real, dimension(:,:), allocatable     :: Points
      ! Local
      integer :: nSpan, nDepth
      integer :: iSpan, iDepth, k
      nSpan  = size(LatticePoints,2)
      nDepth = size(LatticePoints,3)

      if (allocated(Connectivity)) deallocate(Connectivity)
      allocate(Connectivity(1:4, 1:(nSpan-1)*(nDepth-1)))
      if (allocated(Points)) deallocate(Points)
      allocate(Points(1:3, 1:nSpan*nDepth))

      k=1
      do iDepth=1,nDepth-1; do iSpan=1,nSpan-1
         Connectivity(1,k)=(iDepth-1)*nSpan+(iSpan-1)
         Connectivity(2,k)=(iDepth-1)*nSpan+(iSpan )
         Connectivity(3,k)=(iDepth  )*nSpan+(iSpan)
         Connectivity(4,k)=(iDepth  )*nSpan+(iSpan-1)
         k=k+1
      enddo; enddo

      k=1
      do iDepth=1,nDepth; do iSpan=1,nSpan
         Points(1:3,k) = LatticePoints(1:3,iSpan,iDepth)
         k=k+1
      enddo; enddo
   end subroutine

   SUBROUTINE GetNewUnit(UnIn)
      INTEGER,        INTENT(OUT) :: UnIn           !< Logical unit for the file.
      INTEGER                     :: Un             !< Unit number
      LOGICAL                     :: Opened         !< Flag indicating whether or not a file is opened.
      INTEGER, PARAMETER          :: StartUnit = 10 !< Starting unit number to check (numbers less than 10 reserved)
      INTEGER, PARAMETER          :: MaxUnit   = 99 !< The maximum unit number available (or 10 less than the number of files you want to have open at a time)
      Un = StartUnit
      DO
         INQUIRE ( UNIT=Un , OPENED=Opened )
         IF ( .NOT. Opened )  EXIT
         Un = Un + 1
         IF ( Un > MaxUnit ) THEN
            Un=-1
            EXIT           ! stop searching now
         END IF
      END DO
      UnIn = Un
      RETURN
   END SUBROUTINE GetNewUnit


   subroutine vtk_misc_init(mvtk)
      type(FVW_VTK_Misc),intent(inout) :: mvtk
      mvtk%vtk_unit = -1           !< VTK output unit [-]
      mvtk%bFileOpen = .false.     !< binary file is open [-]
      mvtk%bBinary = .false.       !< write binary files [-]
      mvtk%nData = 0               !< number of data lines [-]
      mvtk%nPoints = 0             !< number of points [-]
   end subroutine

    !>
    subroutine set_vtk_binary_format(bBin,mvtk)
        logical, intent(in)::bBin
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        mvtk%bBinary=bBin
    end subroutine

    
    !> Save a coordinate transform
    ! ALL VTK Will be exported in this coordinate system!
    subroutine set_vtk_coordinate_transform(T_g2b_in,PO_g_in,mvtk)
        real,dimension(3,3), intent(in) :: T_g2b_in
        real,dimension(3)  , intent(in) :: PO_g_in
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        mvtk%T_g2b=T_g2b_in
        mvtk%PO_g=PO_g_in
    end subroutine

    logical function vtk_new_ascii_file(filename,label,mvtk)
        !use MainIO,     only: get_free_unit ,check_io
        !use MainIOData, only: bSTOP_ALLOWED
        !use FileSystem, only: file_exists
        !use Logging,    only: log_warning,log_error,log_info
        !
        character(len=*),intent(in)      :: filename
        character(len=*),intent(in)      :: label
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        !
        integer :: iostatvar
        logical :: b

        if (.not. mvtk%bFileOpen) then
            CALL GetNewUnit( mvtk%vtk_unit )   
            if (mvtk%bBinary) then
                ! Fortran 2003 stream, otherwise intel fortran !
                !form='UNFORMATTED',access='SEQUENTIAL',action='WRITE',convert='BIG_ENDIAN',recordtype='STREAM',buffered='YES',
               !print*,'Not available for this compiler' !COMPAQ-COMPILER
               !STOP !COMPAQ-COMPILER
!bjj: CONVERT is non-standard, so maybe this should be part of Sys*.f90? Like OpenUnfInpBEFile()?
                open(unit = mvtk%vtk_unit,file= trim(adjustl(filename)),form='UNFORMATTED',access = 'stream',& !OTHER-COMPILER
                    action = 'WRITE',convert= 'BIG_ENDIAN',iostat=iostatvar,status='replace') !OTHER-COMPILER
            else
                open(mvtk%vtk_unit,file=trim(adjustl(filename)),iostat=iostatvar,action="write",status='replace')
            endif
            if (iostatvar == 0) then
                if (mvtk%bBinary) then
                    write(mvtk%vtk_unit)'# vtk DataFile Version 3.0'//NL
                    write(mvtk%vtk_unit)trim(label)//NL
                    write(mvtk%vtk_unit)'BINARY'//NL
                else
                    write(mvtk%vtk_unit,'(a)') '# vtk DataFile Version 2.0'
                    write(mvtk%vtk_unit,'(a)') trim(label)
                    write(mvtk%vtk_unit,'(a)') 'ASCII'
                    write(mvtk%vtk_unit,'(a)') ' '
                endif

                mvtk%bFileOpen=.true.
                mvtk%nData=-1;
            endif
        else
            b=.false.
            !call log_error('VTK: Cannot open two vtk files at the same time, call vtk_close first')
        endif
        if (iostatvar ==0) then
           vtk_new_ascii_file=.true.
        else
           vtk_new_ascii_file=.false.
        endif
    end function

    subroutine vtk_close_file(mvtk)
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            close(mvtk%vtk_unit)
            mvtk%bFileOpen=.false.
        endif
    endsubroutine


    ! ------------------------------------------------------------------------- 
    ! --- POLYDATA STUFF
    ! ------------------------------------------------------------------------- 
    subroutine vtk_dataset_polydata(Points,mvtk,bladeFrame)
        real, dimension(:,:),intent(in) :: Points  !< 3 x n
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        logical, intent(in) :: bladeFrame
        integer :: i
        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=size(Points,2)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'DATASET POLYDATA'//NL
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints ,' double'
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
                if (bladeFrame)  then
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit)matmul(mvtk%T_g2b,Points(1:3,i)-mvtk%PO_g)
                    enddo
                else
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit)Points(1:3,i)
                    enddo
                endif
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET POLYDATA'
                write(mvtk%vtk_unit,'(A,I0,A)') 'POINTS ', mvtk%nPoints ,' double'
                if (bladeFrame)  then
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit,'(3'//RFMT//')') matmul(mvtk%T_g2b,Points(1:3,i)-mvtk%PO_g)
                    enddo
                else
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit,'(3'//RFMT//')') Points(1:3,i)
                    enddo
                endif
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine


    subroutine vtk_lines(L,mvtk)
        integer, dimension(:,:),intent(in) :: L    !< 2 x n
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        integer :: i

        if ( mvtk%bFileOpen ) then
            mvtk%nData=size(L,2)
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0,A,I0)')'LINES ',mvtk%nData,' ',3*mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit)2,L(1:2,i)
                enddo
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,I0,A,I0)')'LINES ',mvtk%nData,' ',3*mvtk%nData
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit,'(3'//IFMT//')') 2, L(1:2,i)
                enddo
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine

    subroutine vtk_quad(Q,mvtk)
        integer, dimension(:,:),intent(in) :: Q    !< 4 x n
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        integer :: i
        if ( mvtk%bFileOpen ) then
            mvtk%nData=size(Q,2)
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0,A,I0)')'POLYGONS ',mvtk%nData,' ',5*mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit)4,Q(1:4,i)
                enddo
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,I0,A,I0)') 'POLYGONS ', mvtk%nData,' ',5*mvtk%nData
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit,'(5'//IFMT//')') 4, Q(1:4,i)
                enddo
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine

    ! ------------------------------------------------------------------------- 
    ! --- RECTILINEAR
    ! ------------------------------------------------------------------------- 
    subroutine vtk_dataset_rectilinear(v1,v2,v3,mvtk)
        real, dimension(:),intent(in) :: v1,v2,v3  !<  n
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=size(v1)*size(v2)*size(v3)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET RECTILINEAR_GRID'//NL
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', size(v1),' ',size(v2),' ',size(v3)
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,I0,A)') 'X_COORDINATES ', size(v1), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)v1
                write(mvtk%vtk_unit)NL
                write(mvtk%buffer,'(A,I0,A)') 'Y_COORDINATES ', size(v2), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)v2
                write(mvtk%vtk_unit)NL
                write(mvtk%buffer,'(A,I0,A)') 'Z_COORDINATES ', size(v3), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)v3
                !write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET RECTILINEAR_GRID'
                write(mvtk%vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', size(v1),' ',size(v2),' ',size(v3)
                write(mvtk%vtk_unit,'(A,I0,A)') 'X_COORDINATES ', size(v1), ' double'
                write(mvtk%vtk_unit,'('//RFMT//')') v1
                write(mvtk%vtk_unit,'(A,I0,A)') 'Y_COORDINATES ', size(v2), ' double'
                write(mvtk%vtk_unit,'('//RFMT//')') v2
                write(mvtk%vtk_unit,'(A,I0,A)') 'Z_COORDINATES ', size(v3), ' double'
                write(mvtk%vtk_unit,'('//RFMT//')') v3
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine

    subroutine vtk_dataset_structured_points(x0,dx,n,mvtk)
        real,     dimension(3), intent(in) :: x0 !< origin
        real,     dimension(3), intent(in) :: dx !< spacing
        integer,        dimension(3), intent(in) :: n  !< length
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n(1)*n(2)*n(3)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_POINTS'//NL
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ',n(1),' ',n(2),' ',n(3)
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,3F16.8)') 'ORIGIN ', x0
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,3F16.8)') 'SPACING ', dx
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET STRUCTURED_POINTS'
                write(mvtk%vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n(1),' ',n(2),' ',n(3)
                write(mvtk%vtk_unit,'(A,3F16.8,A)') 'ORIGIN  ',x0
                write(mvtk%vtk_unit,'(A,3F16.8,A)') 'SPACING ',dx
            endif
        endif
    end subroutine


    ! ------------------------------------------------------------------------- 
    ! --- STRUCTURED GRID (Points dumped without for loop since memory is in proper order)
    ! ------------------------------------------------------------------------- 
    !> Subroutine using flat data as input (not in natural order)
    subroutine vtk_dataset_structured_grid_flat(D,n1,n2,n3,mvtk)
        integer , intent(in) :: n1,n2,n3
        real, dimension(:,:),intent(in)::D
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n1*n2*n3
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_GRID'//NL
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET STRUCTURED_GRID'
                write(mvtk%vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine

    !> Using Grid data 4d as input
    subroutine vtk_dataset_structured_grid_grid(D,n1,n2,n3,mvtk)
        integer , intent(in) :: n1,n2,n3
        real, dimension(:,:,:,:),intent(in)::D
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n1*n2*n3
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_GRID'//NL
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET STRUCTURED_GRID'
                write(mvtk%vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine
    
    
    
    ! ------------------------------------------------------------------------- 
    ! --- POINT DATA
    ! ------------------------------------------------------------------------- 
    subroutine vtk_point_data_init(mvtk)
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if(mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0)')'POINT_DATA ',mvtk%nPoints
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
            else
                write(mvtk%vtk_unit,'(A,I0)') 'POINT_DATA ', mvtk%nPoints
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_flat(D,sname,mvtk)
        real, dimension(:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_grid(D,sname,mvtk)
        real, dimension(:,:,:,:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_grid2D(D,sname,mvtk)
        real, dimension(:,:,:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    !>
    subroutine vtk_point_data_vector_flat(D,sname,mvtk)
        real, dimension(:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine
    !>
    subroutine vtk_point_data_vector_grid(D,sname,mvtk)
        real, dimension(:,:,:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine
    !>
    subroutine vtk_point_data_vector_grid2D(D,sname,mvtk)
        real, dimension(:,:,:),intent(in) :: D  !< 
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine


    ! ------------------------------------------------------------------------- 
    ! --- CELL DATA
    ! ------------------------------------------------------------------------- 
    subroutine vtk_cell_data_init(mvtk)
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0)')'CELL_DATA ',mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
            else
                write(mvtk%vtk_unit,'(A,I0)') 'CELL_DATA ', mvtk%nData
            endif
        endif
    end subroutine

    subroutine vtk_cell_data_scalar_1d(D,sname,mvtk)
        real, dimension(:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double 1'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,fmt='(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_cell_data_scalar_2d(D,sname,mvtk)
        real, dimension(:,:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double 1'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,fmt='(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine
    
    
    subroutine vtk_cell_data_vector(D,sname,mvtk)
        real, dimension(:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine

    ! --------------------------------------------------------------------------------}
    ! --- VTK Tools 
    ! --------------------------------------------------------------------------------{
    !> Exports a Plane From a mesh
    subroutine export_plane_grid3d(fname,v1,v2,v3,Values,mvtk)
        character(len=*),intent(in)             :: fname
        real,dimension(:), intent(in)       :: v1,v2,v3
        real,dimension(:,:,:,:), intent(in) :: Values
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        !  Variables
        integer :: nD

        ! Writting
        if ( vtk_new_ascii_file(trim(fname),'grid',mvtk)) then
            nD=size(Values,1)
            call vtk_dataset_rectilinear(v1,v2,v3,mvtk)
            ! Output as a structured grid, No need to reorder
            call vtk_point_data_init(mvtk)
            ! Could be a function of nDim, be careful
            if(nD==3) then
                call vtk_point_data_vector(Values(1:3,:,:,:),'Velocity',mvtk) ! Label...
            endif

            call vtk_close_file(mvtk)
        endif ! file opening
    end subroutine
    
    !> Exports a Plane From a mesh
    subroutine export_plane_grid2d(fname,v1,v2,v3,Values,mvtk)
        character(len=*),intent(in)             :: fname
        real,dimension(:), intent(in)       :: v1,v2,v3
        real,dimension(:,:,:), intent(in) :: Values
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        !  Variables
        integer :: nD

        ! Writting
        if ( vtk_new_ascii_file(trim(fname),'plane',mvtk) ) then
            nD=size(Values,1)
            call vtk_dataset_rectilinear(v1,v2,v3,mvtk)
            ! Output as a structured grid, No need to reorder
            call vtk_point_data_init(mvtk)
            ! Could be a function of nDim, be careful
            if(nD==3) then
                call vtk_point_data_vector(Values(1:3,:,:),'Velocity',mvtk) ! Label...
            endif

            call vtk_close_file(mvtk)
        endif ! file opening
    end subroutine
end module FVW_VTK

subroutine WriteWakeElemData()

    ! Write wake element positions and velocity

    use wakedata
    use blade
    use wake
    use wallsoln
    use configr
    use fnames
    use pathseparator

    implicit none

    integer :: tCount, tCountMax, wcount, node_id
    character(len=10) :: nt_str

    ! Optional wake element data output
    write(nt_str,'(I5.5)') nt
    WakeOutputFN=adjustl(trim(WakeElemOutputPath))//path_separator//trim(FNBase)//'_WakeElemData_'//trim(nt_str)//'.csv'
    OPEN(12, FILE=WakeOutputFN)
    write(12,'(A)') trim(WakeOutHead)

    tCountMax=nt
    do wcount=1,NWakeInd+nb
        do tCount=1,tCountMax
            ! compute the unique ID number of the current wake filament
            ! (higher numbers correspond to newer elements)
            node_id = (NWakeInd+nb)*(tcount-1) + wcount

            write(12,'(E13.7,",",$)') TimeN                ! Normalized simulation time (t*Uinf/Rmax)
            write(12,'(I0,",",$)') node_id                 ! Unique node ID
            write(12,'(I0,",",$)') wcount                  ! Node number that wake element originated from
            write(12,'(E13.7,",",$)') X(tCount,wcount)     ! Wake element X position
            write(12,'(E13.7,",",$)') Y(tCount,wcount)     ! Wake element Y position
            write(12,'(E13.7,",",$)') Z(tCount,wcount)     ! Wake element Z position
            write(12,'(E13.7,",",$)') U(tCount,wcount)     ! Wake element X velocity
            write(12,'(E13.7,",",$)') V(tCount,wcount)     ! Wake element Y velocity
            ! Dont suppress carriage return on last column
            write(12,'(E13.7)') W(tCount,wcount)           ! Wake element Z velocity
        end do
    end do

    ! close the output file
    CLOSE(12)

    return
end subroutine WriteWakeElemData


subroutine WriteWakeVTK()
   use FVW_VTK
    ! Write wake element positions and velocity
    use wakedata
    use blade
    use wake
    use wallsoln
    use configr
    use fnames
    use pathseparator
    implicit none

    integer :: tCount, tCountMax, wcount, node_id
    character(len=10) :: nt_str
    type(FVW_VTK_Misc) :: mvtk 
    real, dimension(:,:),     allocatable :: SegPoints  !< 
    integer, dimension(:,:),  allocatable:: SegConnct  !< 
    real,     dimension(:)  , allocatable :: SegGamma   !< 
    real,     dimension(:)  , allocatable :: SegEpsilon !< 
    real,     dimension(:,:) , allocatable :: LatticeGamma !< !< Array of            nSpan x nDepth
    real,     dimension(:,:,:) , allocatable :: LatticePoints !<  !< Array of points 3 x nSpan x nDepth
    character(1), dimension(26) :: I2ABC =(/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/)

    integer :: iSpan, iDepth
    integer :: nei, nej, i ,it,j

    allocate(LatticePoints(3,nbe+1,nt))
    allocate(LatticeGamma(nbe,nt-1))
    LatticeGamma=0.0

    ! Set new wake element positions
    print*,'>>> nb, nbe',nb,nbe
    !print*,'>>>x', x(1,:)
    !print*,'>>>y',-z(1,:)
    !print*,'>>>z', y(1,:)
    !print*,'>>>x', x(2,:)
    !print*,'>>>y',-z(2,:)
    !print*,'>>>z', y(2,:)
    !print*,'>>>x', x(3,:)
    !print*,'>>>y',-z(3,:)
    !print*,'>>>z', y(3,:)
    !print*,'GS',GS(1,:)
    !print*,'GS',GS(2,:)
    !print*,'GS',GS(3,:)
    !print*,'GT',GT(1,:)
    !print*,'GT',GT(2,:)
    !print*,'GT',GT(3,:)

    do i=1,nb ! Loop on blades
        nei=1+(i-1)*(nbe+1)
        do it =1,nt
           iDepth=nt-it+1 ! Elements are stored in reversed order
           ! --- Wake nodes
           do j=0,nbe
              nej=nei+j ! element index
              !iSpan = nbe+1-j
              iSpan = j+1
              LatticePoints(1, iSpan, iDepth) =   x(it,nej)
              LatticePoints(2, iSpan, iDepth) = - z(it,nej)
              LatticePoints(3, iSpan, iDepth) =   y(it,nej)
           enddo ! Loop on span
           ! --- Gamma ! TODO verify
           if (iDepth<nt) then
              j=0
              nej=nei+j ! element index
              iSpan = j+1
              LatticeGamma(iSpan,iDepth) = GT(it, nej)
              do j=1,nbe-1
                 nej=nei+j ! element index
                 !iSpan = nbe+1-j
                 iSpan = j+1
                 LatticeGamma(iSpan,iDepth)  = LatticeGamma(iSpan-1, iDepth) + GT(it,nej)
              enddo
              !j=nbe
              !nej=nei+j ! element index
              !iSpan = j+1
              !LatticeGamma(iSpan,iDepth) = - GT(it, nej)
           endif
           ! TODO TODO LatticeData3d Velocity
           ! write(12,'(E13.7,",",$)') U(tCount,wcount)     ! Wake element X velocity
           ! write(12,'(E13.7,",",$)') V(tCount,wcount)     ! Wake element Y velocity
           ! write(12,'(E13.7)') W(tCount,wcount)           ! Wake element Z velocity
        enddo ! Loop on depth

       write(nt_str,'(I9.9)') nt
       WakeOutputFN=adjustl(trim(WakeVTKOutputPath)//path_separator//trim(FNBase)//'_Bld'//I2ABC(i)//'.'//trim(nt_str)//'.vtk')
       print*,'>>>',trim(WakeOutputFN)
       ! NOTE GT GS trailed,shed TODO TODO
      call WrVTK_Lattice(WakeOutputFN, mvtk, LatticePoints, LatticeGamma, bladeFrame=.False.) !, LatticeData3d, bladeFrame)
    end do
    !print*,'>>>x',LatticePoints(1,:,1)
    !print*,'>>>y',LatticePoints(2,:,1)
    !print*,'>>>z',LatticePoints(3,:,1)

    !if (nt>1) then
    !   print*,'>>>x',LatticePoints(1,:,2)
    !   print*,'>>>y',LatticePoints(2,:,2)
    !   print*,'>>>z',LatticePoints(3,:,2)
    !endif

    !if (nt>2) then
    !   print*,'>>>x',LatticePoints(1,:,3)
    !   print*,'>>>y',LatticePoints(2,:,3)
    !   print*,'>>>z',LatticePoints(3,:,3)
    !endif
   if (allocated(LatticePoints)) deallocate(LatticePoints)
   if (allocated(LatticeGamma))  deallocate(LatticeGamma)
   if (allocated(SegConnct))     deallocate(SegConnct)
   if (allocated(SegPoints))     deallocate(SegPoints)
   if (allocated(SegGamma))      deallocate(SegGamma)
   if (allocated(SegEpsilon))    deallocate(SegEpsilon)

end subroutine

