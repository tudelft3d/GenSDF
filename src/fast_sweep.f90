module fast_sweep_mod
    !
    ! fast_sweep_mod — Godunov Fast Sweeping Method for Eikonal propagation.
    !
    use utils_io
    use mpi
    implicit none

    ! Large sentinel used for uninitialised cells (same sign convention as scalarvalue)
    real(dp), parameter :: INF_D = 1.0e30_dp

    ! Ghost plane buffers: allocated once on first call to fast_sweep_3d and reused.
    real(dp), allocatable :: send_buf_mod(:,:)
    real(dp), allocatable :: ghost_left_mod(:,:)
    real(dp), allocatable :: ghost_right_mod(:,:)

contains

    ! 
    !  Godunov upwind solution of |Grad phi| = 1 at a single cell.
    ! 
    pure function solve_eikonal(phix, phiy, phiz, hx, hy, hz) result(phi_new)
        !$acc routine seq
        implicit none
        real(dp), intent(in) :: phix, phiy, phiz, hx, hy, hz
        real(dp) :: phi_new

        ! Sort the three directional minima so p1 ≤ p2 ≤ p3
        real(dp) :: p(3), h(3), tmp
        real(dp) :: a, b, c, disc, sol

        p(1) = phix; h(1) = hx
        p(2) = phiy; h(2) = hy
        p(3) = phiz; h(3) = hz

        ! Simple bubble sort on 3 elements (by absolute value of p)
        if (abs(p(1)) > abs(p(2))) then
            tmp = p(1); p(1) = p(2); p(2) = tmp
            tmp = h(1); h(1) = h(2); h(2) = tmp
        end if
        if (abs(p(2)) > abs(p(3))) then
            tmp = p(2); p(2) = p(3); p(3) = tmp
            tmp = h(2); h(2) = h(3); h(3) = tmp
        end if
        if (abs(p(1)) > abs(p(2))) then
            tmp = p(1); p(1) = p(2); p(2) = tmp
            tmp = h(1); h(1) = h(2); h(2) = tmp
        end if

        ! Attempt 1-D solution along the smallest-value direction
        sol = p(1) + sign(1.0_dp, p(1)) * h(1)
        if (abs(sol) <= abs(p(2))) then
            phi_new = sol
            return
        end if

        ! Attempt 2-D solution
        a =  1.0_dp/(h(1)*h(1)) + 1.0_dp/(h(2)*h(2))
        b = -2.0_dp*(p(1)/(h(1)*h(1)) + p(2)/(h(2)*h(2)))
        c =  (p(1)/h(1))**2 + (p(2)/h(2))**2 - 1.0_dp
        disc = b*b - 4.0_dp*a*c
        if (disc >= -1.0e-14_dp) then
            disc = max(0.0_dp, disc)
            ! pick the root with same sign as p(1)
            if (p(1) >= 0.0_dp) then
                sol = (-b + sqrt(disc)) / (2.0_dp*a)
            else
                sol = (-b - sqrt(disc)) / (2.0_dp*a)
            end if
            if (abs(sol) <= abs(p(3))) then
                phi_new = sol
                return
            end if
        end if

        ! Attempt 3-D solution
        a =  1.0_dp/(h(1)*h(1)) + 1.0_dp/(h(2)*h(2)) + 1.0_dp/(h(3)*h(3))
        b = -2.0_dp*(p(1)/(h(1)*h(1)) + p(2)/(h(2)*h(2)) + p(3)/(h(3)*h(3)))
        c =  (p(1)/h(1))**2 + (p(2)/h(2))**2 + (p(3)/h(3))**2 - 1.0_dp
        disc = b*b - 4.0_dp*a*c
        if (disc >= -1.0e-14_dp) then
            disc = max(0.0_dp, disc)
            if (p(1) >= 0.0_dp) then
                phi_new = (-b + sqrt(disc)) / (2.0_dp*a)
            else
                phi_new = (-b - sqrt(disc)) / (2.0_dp*a)
            end if
        else
            ! Degenerate (nearly co-planar): fall back to 1-D
            phi_new = p(1) + sign(1.0_dp, p(1)) * h(1)
        end if
    end function solve_eikonal

    !  Single-cell Godunov update, shared by the CPU and GPU sweep drivers.
    subroutine godunov_update_cell(sdf, xstart, xend, ny_loc, nz_loc, ii, jj, kk, &
                                    ghost_left, ghost_right)
        !$acc routine seq
        implicit none
        integer,  intent(in)    :: xstart, xend, ny_loc, nz_loc, ii, jj, kk
        real(dp), intent(inout) :: sdf(xstart:xend, ny_loc, nz_loc)
        real(dp), intent(in)    :: ghost_left(ny_loc, nz_loc)   ! ghost col at xstart-1
        real(dp), intent(in)    :: ghost_right(ny_loc, nz_loc)  ! ghost col at xend+1
        real(dp) :: phix, phiy, phiz, phi_try, hx, hy, hz

        ! Skip cells already initialised by narrow-band phase (seeded |phi| << scalarvalue sentinel).
        if (abs(sdf(ii,jj,kk)) < scalarvalue * 0.5_dp) return

        hx = dx           ! uniform in x (xf array assumed uniform in this tool)
        hy = dyv(jj)      ! per-cell (spanwise stretching); equals dy on a uniform grid
        hz = dz(kk)       ! may be non-uniform

        ! Upwind neighbours in x: use ghost planes at slab boundaries
        if (ii == xstart) then
            phix = merge(ghost_left(jj,kk), sdf(ii+1,jj,kk), &
                         abs(ghost_left(jj,kk)) <= abs(sdf(ii+1,jj,kk)))
        else if (ii == xend) then
            phix = merge(sdf(ii-1,jj,kk), ghost_right(jj,kk), &
                         abs(sdf(ii-1,jj,kk)) <= abs(ghost_right(jj,kk)))
        else
            phix = merge(sdf(ii-1,jj,kk), sdf(ii+1,jj,kk), &
                         abs(sdf(ii-1,jj,kk)) <= abs(sdf(ii+1,jj,kk)))
        end if

        ! Upwind neighbours in y (local j index)
        if (jj == 1) then
            phiy = sdf(ii, 2, kk)
        else if (jj == ny_loc) then
            phiy = sdf(ii, ny_loc-1, kk)
        else
            phiy = merge(sdf(ii,jj-1,kk), sdf(ii,jj+1,kk), &
                         abs(sdf(ii,jj-1,kk)) <= abs(sdf(ii,jj+1,kk)))
        end if

        ! Upwind neighbours in z (local k index)
        if (kk == 1) then
            phiz = sdf(ii, jj, 2)
        else if (kk == nz_loc) then
            phiz = sdf(ii, jj, nz_loc-1)
        else
            phiz = merge(sdf(ii,jj,kk-1), sdf(ii,jj,kk+1), &
                         abs(sdf(ii,jj,kk-1)) <= abs(sdf(ii,jj,kk+1)))
        end if

        ! Sign of current cell (set by ray-crossing classification) is preserved during FSM propagation.
        phi_try = solve_eikonal(abs(phix), abs(phiy), abs(phiz), hx, hy, hz) &
                  * sign(1.0_dp, sdf(ii,jj,kk))
        if (abs(phi_try) < abs(sdf(ii,jj,kk))) sdf(ii,jj,kk) = phi_try
    end subroutine godunov_update_cell

    !
    !  One complete Godunov FSM sweep in the directions dictated by sx/sy/sz
    !  (each ±1, giving 8 combinations).
    !
    subroutine do_sweep(sdf, xstart, xend, ystart, yend, zstart, zend, &
                        sx_dir, sy_dir, sz_dir, xin, yin, zin, &
                        ghost_left, ghost_right)
        implicit none
        integer,  intent(in)    :: xstart, xend, ystart, yend, zstart, zend
        integer,  intent(in)    :: sx_dir, sy_dir, sz_dir
        real(dp), intent(in)    :: xin(:), yin(:), zin(:)
        real(dp), intent(inout) :: sdf(xstart:xend, 1:size(yin), 1:size(zin))
        real(dp), intent(in)    :: ghost_left(size(yin), size(zin))   ! ghost col at xstart-1
        real(dp), intent(in)    :: ghost_right(size(yin), size(zin))  ! ghost col at xend+1

        integer  :: ii, jj, kk, i_lo, i_hi, i_step
        integer  :: j_lo, j_hi, j_step, k_lo, k_hi, k_step
        integer  :: nx_loc, ny_loc, nz_loc

        nx_loc = xend - xstart + 1
        ny_loc = size(yin)
        nz_loc = size(zin)

        ! Build loop bounds from sweep direction
        if (sx_dir > 0) then; i_lo = xstart; i_hi = xend;  i_step =  1
        else;                  i_lo = xend;   i_hi = xstart; i_step = -1; end if
        if (sy_dir > 0) then; j_lo = 1;       j_hi = ny_loc; j_step =  1
        else;                  j_lo = ny_loc;  j_hi = 1;     j_step = -1; end if
        if (sz_dir > 0) then; k_lo = 1;       k_hi = nz_loc; k_step =  1
        else;                  k_lo = nz_loc;  k_hi = 1;     k_step = -1; end if

        do ii = i_lo, i_hi, i_step
            do jj = j_lo, j_hi, j_step
                do kk = k_lo, k_hi, k_step
                    call godunov_update_cell(sdf, xstart, xend, ny_loc, nz_loc, ii, jj, kk, &
                                              ghost_left, ghost_right)
                end do
            end do
        end do
    end subroutine do_sweep

#ifdef GPU_SDF
    !  GPU diagonal-wavefront reformulation of do_sweep (same dependency structure, reordered — see docs/external-lib-guidance/GenSDF_GPU_porting.md).
    subroutine do_sweep_wavefront_gpu(sdf, xstart, xend, ystart, yend, zstart, zend, &
                                       sx_dir, sy_dir, sz_dir, xin, yin, zin, &
                                       ghost_left, ghost_right)
        implicit none
        integer,  intent(in)    :: xstart, xend, ystart, yend, zstart, zend
        integer,  intent(in)    :: sx_dir, sy_dir, sz_dir
        real(dp), intent(in)    :: xin(:), yin(:), zin(:)
        real(dp), intent(inout) :: sdf(xstart:xend, 1:size(yin), 1:size(zin))
        real(dp), intent(in)    :: ghost_left(size(yin), size(zin))
        real(dp), intent(in)    :: ghost_right(size(yin), size(zin))

        integer :: nx_loc, ny_loc, nz_loc, level, level_max
        integer :: ip, jp, kp, ii, jj, kk

        nx_loc = xend - xstart + 1
        ny_loc = size(yin)
        nz_loc = size(zin)
        level_max = nx_loc + ny_loc + nz_loc

        do level = 3, level_max
            !$acc parallel loop collapse(2) default(none) &
            !$acc    private(ip, jp, kp, ii, jj, kk) &
            !$acc    present(sdf, ghost_left, ghost_right) &
            !$acc    firstprivate(level, nx_loc, ny_loc, nz_loc, xstart, xend, sx_dir, sy_dir, sz_dir)
            do jp = 1, ny_loc
                do ip = 1, nx_loc
                    kp = level - ip - jp
                    if (kp >= 1 .and. kp <= nz_loc) then
                        if (sx_dir > 0) then; ii = xstart + (ip-1); else; ii = xend - (ip-1);   end if
                        if (sy_dir > 0) then; jj = 1 + (jp-1);      else; jj = ny_loc - (jp-1); end if
                        if (sz_dir > 0) then; kk = 1 + (kp-1);      else; kk = nz_loc - (kp-1); end if
                        call godunov_update_cell(sdf, xstart, xend, ny_loc, nz_loc, ii, jj, kk, &
                                                  ghost_left, ghost_right)
                    end if
                end do
            end do
        end do
    end subroutine do_sweep_wavefront_gpu
#endif

    !
    !  Exchange a single yz ghost plane with a neighbour rank via MPI_SENDRECV.
    !  direction:  +1 → send right face (ii=xend),   fill ghost_out with right-neighbour data
    !              -1 → send left  face (ii=xstart),  fill ghost_out with left-neighbour data
    ! 
    subroutine exchange_ghost_plane(sdf, xstart, xend, ny_loc, nz_loc, &
                                    direction, myid, neighbour, &
                                    send_buf, ghost_out, ierror)
        use mpi
        implicit none
        integer,  intent(in)    :: xstart, xend, ny_loc, nz_loc
        integer,  intent(in)    :: direction, myid, neighbour
        real(dp), intent(in)    :: sdf(xstart:xend, 1:ny_loc, 1:nz_loc)
        real(dp), intent(inout) :: send_buf(ny_loc, nz_loc)
        real(dp), intent(out)   :: ghost_out(ny_loc, nz_loc)
        integer,  intent(out)   :: ierror
        integer  :: recv_status(MPI_STATUS_SIZE)
        integer  :: tag_send, tag_recv

        if (neighbour < 0 .or. neighbour >= nprocs) then
            ! No neighbour mainly for Neumann BC: ghost = boundary plane
            if (direction > 0) then
                ghost_out = sdf(xend, :, :)
            else
                ghost_out = sdf(xstart, :, :)
            end if
            ierror = 0
            return
        end if

        tag_send = 100 + direction
        tag_recv = 100 - direction

        if (direction > 0) then
            send_buf = sdf(xend, :, :)
        else
            send_buf = sdf(xstart, :, :)
        end if

        call MPI_SENDRECV(send_buf, ny_loc*nz_loc, MPI_DOUBLE_PRECISION, neighbour, tag_send, &
                          ghost_out, ny_loc*nz_loc, MPI_DOUBLE_PRECISION, neighbour, tag_recv, &
                          MPI_COMM_WORLD, recv_status, ierror)
    end subroutine exchange_ghost_plane

    ! 
    !  8-direction Godunov FSM sweeps over the domain slab owned by this rank.
    !  Ghost-plane MPI exchange with left_rank / right_rank after each
    !  x-direction reversal (16 exchanges total = 2 per sweep × 8 sweeps).
    ! 
    subroutine fast_sweep_3d(sdf, xstart, xend, ystart, yend, zstart, zend, &
                              xin, yin, zin, myid, left_rank, right_rank)
        implicit none
        integer,  intent(in)    :: xstart, xend, ystart, yend, zstart, zend
        integer,  intent(in)    :: myid, left_rank, right_rank
        real(dp), intent(in),  dimension(:) :: xin, yin, zin
        real(dp), intent(inout) :: sdf(xstart:xend, 1:size(yin), 1:size(zin))

        integer  :: sw, sx_dir, sy_dir, sz_dir, ierror
        integer  :: ny_loc, nz_loc

        ny_loc = size(yin)
        nz_loc = size(zin)

        ! Allocate ghost buffers on first call (or when grid size changes)
        if (.not. allocated(send_buf_mod) .or. &
            size(send_buf_mod,1) /= ny_loc .or. size(send_buf_mod,2) /= nz_loc) then
            if (allocated(send_buf_mod))  deallocate(send_buf_mod)
            if (allocated(ghost_left_mod))  deallocate(ghost_left_mod)
            if (allocated(ghost_right_mod)) deallocate(ghost_right_mod)
            allocate(send_buf_mod(ny_loc, nz_loc))
            allocate(ghost_left_mod(ny_loc, nz_loc), ghost_right_mod(ny_loc, nz_loc))
#ifdef GPU_SDF
            !$acc enter data create(ghost_left_mod, ghost_right_mod)
#endif
        end if
        ghost_left_mod  = INF_D
        ghost_right_mod = INF_D

#ifdef GPU_SDF
        !$acc data copy(sdf)
#endif
        ! 8 sweeps: all combinations of ±x, ±y, ±z.
        ! Each sweep processes cells in a fixed order along each axis, propagating
        ! the distance wave the full domain length in one pass per direction.
        ! With ray-crossing classification providing correct ±scalarvalue sentinels, a single
        ! round of 8 sweeps covers all octants.
        do sw = 0, 7
            sx_dir = merge( 1, -1, iand(sw, 1) == 0 )
            sy_dir = merge( 1, -1, iand(sw, 2) == 0 )
            sz_dir = merge( 1, -1, iand(sw, 4) == 0 )

#ifdef GPU_SDF
            ! Only the two boundary planes are needed on the host side for MPI.
            !$acc update host(sdf(xstart:xstart,:,:), sdf(xend:xend,:,:))
#endif
            ! Exchange ghost planes so boundary stencils use fresh neighbour data.
            call exchange_ghost_plane(sdf, xstart, xend, ny_loc, nz_loc, &
                                       -1, myid, left_rank,  send_buf_mod, ghost_left_mod,  ierror)
            call exchange_ghost_plane(sdf, xstart, xend, ny_loc, nz_loc, &
                                       +1, myid, right_rank, send_buf_mod, ghost_right_mod, ierror)

#ifdef GPU_SDF
            !$acc update device(ghost_left_mod, ghost_right_mod)
            call do_sweep_wavefront_gpu(sdf, xstart, xend, ystart, yend, zstart, zend, &
                          sx_dir, sy_dir, sz_dir, xin, yin, zin, ghost_left_mod, ghost_right_mod)
#else
            call do_sweep(sdf, xstart, xend, ystart, yend, zstart, zend, &
                          sx_dir, sy_dir, sz_dir, xin, yin, zin, ghost_left_mod, ghost_right_mod)
#endif
        end do
#ifdef GPU_SDF
        !$acc end data
#endif

    end subroutine fast_sweep_3d

end module fast_sweep_mod
