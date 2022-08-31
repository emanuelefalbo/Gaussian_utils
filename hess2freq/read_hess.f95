 module read_parms
  use iso_fortran_env
 
 implicit none
 contains

 subroutine read_coords_hess(file_input, N_atoms, coords_xyz, hess)
 implicit none
 integer, parameter :: dp = selected_real_kind(15, 307)
 real(dp), parameter :: pi = 4.d0*datan(1.d0) 
 real(dp), dimension(:,:), allocatable :: hess_tmp
 real(dp), dimension(:,:), allocatable :: hess
 real(dp), dimension(:,:), allocatable :: coords_xyz
 real(dp), dimension(:,:), allocatable :: coords_tmp
 real(dp), dimension(:), allocatable :: hess_1D
 real(dp), dimension(:), allocatable :: coords_1D
 character(len=*), parameter :: match_natoms = 'Number of atoms'
 character(len=*), parameter :: match_coords = 'Current cartesian coordinates'
 character(len=*), parameter :: match_hess_cc = 'Cartesian Force Constants'
 character(len=*) :: file_input
 character(len=200) :: text
 character(len=100) :: string
 character(len=100) :: filename
 integer :: i
 integer :: N_1D
 integer :: ios
 integer :: N_atoms, N_forces, N_coords, nlines
 integer :: fi = 18, fo = 38
 integer :: len_hess, len_tmp
 integer :: i_low, i_up
 
 open(fi, file = trim(file_input) )
 do 
    read(fi,'(A)', iostat=ios) text
    if (ios == iostat_end) exit
    if (text(1:15) == match_natoms) then       
       string = adjustl(trim(text(58:61)))
       read(string,'(I3)') N_atoms                  ! Convert N. of atoms into Integer
    else if (text(1:29) == match_coords) then
             N_coords = 3*N_atoms
             nlines = int ( N_coords/5 + 1 )
             allocate(coords_tmp(5,nlines))
    do i = 1, nlines
       read(fi, *, iostat=ios) coords_tmp(:,i)       ! Reading in Cartesian Coordinates
    end do 
    else if (text(1:25) == match_hess_cc) then
             string = adjustl(trim(text(50:61)))
             read(string,'(I5)') N_forces                 ! Convert String into Integer
             nlines =  int( N_forces/5 + 1 ) 
             allocate(hess_tmp(5,nlines))
    do i = 1, nlines
       read(fi, *, iostat=ios) hess_tmp(:,i)     ! Reading in Hessian matrix
    end do 
    end if
 end do
 close(fi)
 
 ! Change from Bohr to Ang
 coords_tmp = coords_tmp*0.529117
 
 ! Convert Cartesian Coords into xyz format
 len_tmp = 3*N_atoms                                  ! Temporary length to reshape coords_tmp
 allocate( coords_1D( len_tmp) )
 coords_1D = reshape( coords_tmp, (/ len_tmp /)   )
 allocate( coords_xyz(3, N_atoms) )
 do i = 1, N_atoms
 coords_xyz(:,i) =  coords_1D( (3*i - 2) : (3*i) )
 end do
 
 len_hess = 3*N_atoms                       ! Size RIC Hessian : N_ric x N_ric 
 allocate( hess(len_hess,len_hess) )
 hess = 0.d0
 
 ! Store unprocessed_hess into 1D array
 N_1D = size(hess_tmp)
 allocate( hess_1D(N_1D) )
 hess_1D = reshape( hess_tmp, (/ N_1D /)   )    
 
 ! Build full Symmetric Hessian:
 ! Adding n(n+1)/2 elements in up/low triangular matrix 
 do i = 1, len_hess
    i_low = int ( 0.5 * i * (i - 1) + 1 )  
    i_up = int ( 0.5 * i * (i + 1) + 1 )  
    hess(i,1:i) = hess_1D( i_low : i_up )
    hess(1:i,i) = hess_1D( i_low : i_up )
 end do
 
 ! Printout Hessian in au units
 filename = file_input( 1 : 2 )//"_hess_cc.dat"
 open(fo, file = trim(filename), status='replace') 
 do i = 1, len_hess
    write(fo,*) hess(i,:)
 end do 
 close(fo)
 
 end subroutine read_coords_hess

  
 subroutine read_coords_hess_ric(file_input, N_atoms, ric_dim, &
                               hess, ric_list)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: pi = 4.d0*datan(1.d0) 
    real(dp), dimension(:,:), allocatable :: hess_tmp
    real(dp), dimension(:), allocatable :: hess_1D
    real(dp), dimension(:,:), allocatable :: hess
    real(dp), dimension(:,:), allocatable :: coords_tmp
    real(dp), dimension(:), allocatable :: coords_1D
    real(dp), dimension(:,:), allocatable :: coords_xyz
    character(len=*), parameter :: match_natoms = 'Number of atoms'
    character(len=*), parameter :: match_coords = 'Current cartesian coordinates'
    character(len=*), parameter :: match_hess_ric = 'Internal Force Constants'
    character(len=*), parameter :: match_ric = 'Redundant internal dimensions'
    character(len=*), parameter :: match_ric_index = 'Redundant internal coordinate indices'
    character(len=200) :: file_input, text
    character(len=100) :: string
    integer, dimension(:,:), allocatable :: ric_list_tmp 
    integer, dimension(:,:), allocatable :: ric_list 
    integer, dimension(:), allocatable :: ric_list_1D 
    integer, dimension(4) :: ric_dim
    integer :: i, N_1D
    integer :: ios
    integer :: N_atoms, N_forces, N_coords, nlines
    integer :: fi = 18
    integer :: len_hess, len_tmp
    integer :: i_low, i_up
    integer :: N_ric, bond_ric, angle_ric, dihe_ric
    integer :: M_ric   ! # Atom indices in RIC list
   
    open(fi, file = trim(file_input) )
    do 
        read(fi,'(A)', iostat=ios) text
        if (ios == iostat_end) exit
        if (text(1:15) == match_natoms) then       
            string = adjustl(trim(text(58:61)))
            read(string,'(I3)') N_atoms                  ! Convert N. of atoms into Integer
          !  write(*,'(A,I3)') text(1:15) // ' = ', N_atoms
        else if (text(1:29) == match_coords) then
            N_coords = 3*N_atoms
            nlines = int ( N_coords/5 + 1 )
            allocate(coords_tmp(5,nlines))
   !         print*, text(1:29)
            do i = 1, nlines
               read(fi, *, iostat=ios) coords_tmp(:,i)       ! Reading in Cartesian Coordinates
   !            write(*,'(5ES15.6)') coords_tmp(:,i)
            end do 
        else if (text(1:29) == match_ric) then
             read(fi, *, iostat=ios) N_ric, bond_ric, angle_ric, dihe_ric 
             ric_dim = [ N_ric, bond_ric, angle_ric, dihe_ric]
          !    write(*, *, iostat=ios) N_ric, bond_ric, angle_ric, dihe_ric 
        else if (text(1:37) == match_ric_index) then
             string = adjustl(trim(text(50:61)))
            read(string,'(I5)') M_ric                 ! Convert String into Integer
            nlines =  int( M_ric/6 + 1 ) 
            allocate(ric_list_tmp(6,nlines))
            do i = 1, nlines
               read(fi, *, iostat=ios) ric_list_tmp(:,i)     ! Reading in RIC list
               ! write(*, *, iostat=ios) ric_list_tmp(:,i)     ! Reading in RIC list
            end do  
        else if (text(1:24) == match_hess_ric) then
            string = adjustl(trim(text(50:61)))
            read(string,'(I5)') N_forces                 ! Convert String into Integer
            nlines =  int( N_forces/5 + 1 ) 
            allocate(hess_tmp(5,nlines))
          !   write(*,'(A)') text(1:24)
            do i = 1, nlines
               read(fi, *, iostat=ios) hess_tmp(:,i)     ! Reading in Hessian matrix
            !   write(*,'(I4,5ES15.6)') i, hess_tmp(:,i)
            end do 
        end if
    end do
    close(fi)
   
    ! Change from Bohr to Ang
    coords_tmp = coords_tmp*0.529117
     
    ! Convert Cartesian Coords into xyz format
    len_tmp = 3*N_atoms                                  ! Temporary length to reshape coords_tmp
    allocate( coords_1D( len_tmp) )
    coords_1D = reshape( coords_tmp, (/ len_tmp /)   )
    allocate( coords_xyz(3, N_atoms) )
    do i = 1, N_atoms
       coords_xyz(:,i) =  coords_1D( (3*i - 2) : (3*i) )
    end do
   
    ! Buil RIC_list 
    allocate( ric_list_1D( M_ric) )
    ric_list_1D = reshape( ric_list_tmp, (/ M_ric /) )
    allocate( ric_list(N_ric, 4))
    ric_list = 0
    do i = 1, N_ric
      ric_list(i,:) = ric_list_1D( 2*(2*i) - 3 : 2*(2*i) )
    end do

    len_hess = N_ric                        ! Size RIC Hessian : N_ric x N_ric 
    allocate( hess(len_hess,len_hess) )
    hess = 0.d0
   
    ! Store unprocessed_hess into 1D array
    N_1D = size(hess_tmp)
    allocate( hess_1D(N_1D) )
    hess_1D = reshape( hess_tmp, (/ N_1D /)   )    
  
    ! Build full Symmetric Hessian:
    ! Adding n(n+1)/2 elements in up/low triangular matrix 
    do i = 1, len_hess
         i_low = int ( 0.5 * i * (i - 1) + 1 )  
         i_up = int ( 0.5 * i * (i + 1) + 1 )  
         hess(i,1:i) = hess_1D( i_low : i_up )
         hess(1:i,i) = hess_1D( i_low : i_up )
    end do

  end subroutine



 subroutine read_mass(file_input, N_atoms, mass_1D)
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 307)
   ! real(dp), dimension(:,:), allocatable :: mass
   real(dp), dimension(:,:), allocatable :: mass_tmp
   real(dp), dimension(:), allocatable :: mass_1D
   character(len=*), parameter :: match_natoms = 'Number of atoms'
   character(len=*), parameter :: match_mass = 'Real atomic weights'
   character(len=*) :: file_input
   character(len=200) :: text
   character(len=100) :: string
   integer :: i
   integer :: ios
   integer :: N_atoms,  nlines
   integer :: fi = 28
   integer :: len_tmp
   
   open(fi, file = trim(file_input) )
   do 
      read(fi,'(A)', iostat=ios) text
      if (ios == iostat_end) exit
      if (text(1:15) == match_natoms) then       
         string = adjustl(trim(text(58:61)))
         read(string,'(I3)') N_atoms                  ! Convert N. of atoms into Integer
      else if (text(1:25) == match_mass) then
               nlines =  int( N_atoms/5 + 1 ) 
               allocate(mass_tmp(5,nlines))
      do i = 1, nlines
         read(fi, *, iostat=ios) mass_tmp(:,i)     ! Reading in Masses
         ! write(*, *, iostat=ios) mass_tmp(:,i)     ! Reading in Hessian matrix
      end do 
      end if
   end do
   close(fi)
   
   ! Convert Masses into 1D-array
   len_tmp = N_atoms                                  ! Temporary length to reshape coords_tmp
   allocate( mass_1D( len_tmp) )
   mass_1D = reshape( mass_tmp, (/ len_tmp /)   )
      
   end subroutine read_mass

 end module
