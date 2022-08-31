program driver

    use read_parms
    use iso_fortran_env
   
    implicit none 
    integer, parameter :: dp = selected_real_kind(15,307)
    character(len=200) :: file_fchk_qm        ! G09/16 QM .fchk output file
    ! character(len=200) :: file_fchk_nb        ! G09/16 MM .fchk output file
    real(dp), dimension(:,:), allocatable :: WM, WH, WHW 
    real(dp), dimension(:,:), allocatable :: hess_cc, hess_ric
    real(dp), dimension(:,:), allocatable :: hess_weighted 
    real(dp), dimension(:,:), allocatable :: coords_xyz 
    real(dp), dimension(:), allocatable :: mass
    real(dp), dimension(:), allocatable :: lamb_cc, lamb_ric
    real(dp), dimension(:), allocatable :: freqs_cc, freqs_ric
    real(dp), dimension(:), allocatable :: vib_freqs
    real(dp), parameter :: bohr2m = 0.529177249d-10
    real(dp), parameter :: au2joule = 4.35974434d-18
    real(dp), parameter :: speed_of_light = 299792458
    real(dp), parameter :: avogadro = 6.0221413d+23
    real(dp), parameter :: pi = 4.d0*datan(1.d0) 
    real(dp) :: vib_constant  
    integer, dimension(:,:), allocatable :: ric_list 
    integer, dimension(4) :: ric_dim
    integer :: N_atoms
    integer :: N_vib
    integer :: N_ric
    integer :: i, j, ii, jj
    ! LAPACK
    integer, parameter :: LWMAX = 1000      
    integer :: LWORK, INFO                        
    real(dp), dimension(LWMAX) :: WORK 
    
    call get_command_argument(1, file_fchk_qm)
    call read_coords_hess(file_fchk_qm, N_atoms, coords_xyz, hess_cc)
    call read_coords_hess_ric(file_fchk_qm, N_atoms, ric_dim, & 
                              & hess_ric, ric_list)
    call  read_mass(file_fchk_qm, N_atoms, mass)

    allocate( hess_weighted(3*N_atoms, 3*N_atoms) )
    allocate( WM(3*N_atoms, 3*N_atoms) )
    allocate( WH(3*N_atoms, 3*N_atoms) )
    allocate( WHW(3*N_atoms, 3*N_atoms) )
    allocate( lamb_cc(3*N_atoms) )
    allocate( freqs_cc(3*N_atoms) )
    ! N_vib = 3*N_atoms + 5
    ! allocate( vib_freqs(N_vib) )

    N_ric = ric_dim(1)
    allocate( lamb_ric(N_ric))
    allocate( freqs_ric(N_ric) )

    hess_weighted = 0.d0
    WM = 0.d0
    WH = 0.d0
    WHW = 0.d0

    ! Build Diagonal Mass matrix WM
    do i = 1, N_atoms
        do j = 1, N_atoms
            ii = (i-1)*3+1
            jj = (j-1)*3+1
            WM(ii:3*i,ii:3*i) = 1/ sqrt( mass(i) )
        end do
    end do 
    do i = 1, ubound(WM, 2)
        do j = 1, ubound(WM,2)
            if ( i /= j ) then
                WM(i,j) = 0.d0
            end if
        end do
    end do

    ! Compute Mass-Weigthed Hessian
    WH = matmul(WM,hess_cc)
    WHW = matmul(WH,WM)

    ! Diagonalize  CC Hess
     LWORK = -1
     call dsyev('N', 'L', ubound(WHW,2), WHW,  ubound(WHW,2), lamb_cc, WORK, LWORK, INFO )
     LWORK = min(LWMAX, int( WORK(1) ) )
     call dsyev('N', 'L', ubound(WHW,2), WHW,  ubound(WHW,2), lamb_cc, WORK, LWORK, INFO )
     if ( INFO.gt.0 ) then
         write(*,*) 'The algorithm failed to compute eigenvalues.'
         stop
     end if

      ! Diagonalize  CC Hess
     LWORK = -1
     call dsyev('N', 'L', ubound(hess_ric,2), hess_ric, &
                &  ubound(hess_ric,2), lamb_ric, WORK, LWORK, INFO )
     LWORK = min(LWMAX, int( WORK(1) ) )
     call dsyev('N', 'L', ubound(hess_ric,2), hess_ric, &
                & ubound(hess_ric,2), lamb_ric, WORK, LWORK, INFO )
     if ( INFO.gt.0 ) then
         write(*,*) 'The algorithm failed to compute eigenvalues.'
         stop
     end if


     vib_constant = sqrt( ( avogadro * au2joule*1000) / (bohr2m*bohr2m) )
     vib_constant = vib_constant / (2*pi*speed_of_light*100)
     freqs_cc =  vib_constant * sqrt(lamb_cc)

     do i = 1, size(freqs_ric)
        freqs_ric(i) =  vib_constant * sqrt( abs(lamb_ric(i)) )
        if ( lamb_ric(i)  < 0 ) then
            freqs_ric(i) =  - freqs_ric(i)
        end if
    end do

    !  print*, freqs(3*N_atoms-6 + N_atoms + 1:)
     open(19, file = 'freq_cc.dat')
     write(19,'(A)') '# Frequencies from XYZ Hessian'
     write(19,'(A,2X,A)') '# nr. mode', 'Frequency (cm-1) '
     j = 0 
     do i = 1, size(freqs_cc)
        if ( freqs_cc(i) > 1.d0 ) then
            j = j + 1 
           write(19,'(I3,2X,F10.4)') j, freqs_cc(i) 
        end if
     end do
     close(19)

     open(29, file = 'freq_ric.dat')
     write(29,'(A)') '# Frequencies from RIC Hessian'
     write(29,'(A,2X,A)') '# nr. mode', 'Frequency (cm-1) '
    !  j = 0 
     do i = 1, size(freqs_ric)
        ! if ( freqs(i) > 1.d0 ) then
            ! j = j + 1 
           write(29,'(I3,2X,F10.4)') i, freqs_ric(i) 
        ! end if
     end do
     close(29)


end program driver 
