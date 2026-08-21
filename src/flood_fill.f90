module flood_fill_mod
    !
    ! flood_fill_mod — iterative flood-fill sign-determination.
    !
    use utils_io, only : dp, scalarvalue
    implicit none

contains

    subroutine fill_internal(grid, inx, iny, inz, sx, sy, sz, ex, ey, ez, large_negative_value)
        !
        ! Mark positive cells enclosed by the geometry surface with large_negative_value.
        !
        implicit none
        integer,  intent(in)    :: inx, iny, inz
        integer,  intent(in)    :: sx, sy, sz, ex, ey, ez
        real(dp), intent(in)    :: large_negative_value
        real(dp), intent(inout) :: grid(inx, iny, inz)

        integer,  allocatable :: stack(:,:), tmp_stack(:,:)
        ! flag: 0=unknown, 1=exterior, 2=interior
        integer,  allocatable :: flag(:,:,:)
        integer :: ii, jj, kk, ci, cj, ck, ni, nj, nk, top, i_dir
        integer, parameter :: n_dirs = 6
        integer :: dirs(3,n_dirs) = reshape([ &
            1, 0, 0,  -1, 0, 0, &
            0, 1, 0,   0,-1, 0, &
            0, 0, 1,   0, 0,-1], shape=[3,n_dirs])
        real(dp), parameter :: SENTINEL_FRAC = 0.5_dp
        integer :: stack_size

        ! Start with a modest stack size; double when needed
        stack_size = 1024
        allocate(flag(sx:ex, sy:ey, sz:ez))
        allocate(stack(3, stack_size))
        flag = 0

        !
        ! Pass 1: BFS from domain boundary through sentinel cells -> EXTERIOR
        !
        ! Only seed exterior from AABB faces that lie strictly inside the physical
        ! domain (the buffer zone). When an AABB face coincides with the physical
        ! domain boundary, the geometry may be cut there, so those cells never get
        ! a narrowband value and remain at +scalarvalue. Seeding them as exterior
        ! would let the BFS leak through the open cut face. Truly exterior cells on
        ! those clamped faces are still reached by BFS propagation from the
        ! buffer-zone seeds on the other faces.
        !
        top = 0
        do kk = sz, ez; do jj = sy, ey; do ii = sx, ex
            if ((ii==sx .and. sx>1)  .or. (ii==ex .and. ex<inx) .or. &
                (jj==sy .and. sy>1)  .or. (jj==ey .and. ey<iny) .or. &
                (kk==sz .and. sz>1)  .or. (kk==ez .and. ez<inz)) then
                if (grid(ii,jj,kk) > 0.0_dp .and. &
                    abs(grid(ii,jj,kk)) >= scalarvalue * SENTINEL_FRAC) then
                    flag(ii,jj,kk) = 1   ! exterior
                    top = top + 1
                    if (top > stack_size) then
                        stack_size = stack_size * 2
                        allocate(tmp_stack(3, stack_size))
                        tmp_stack(:, 1:top-1) = stack(:, 1:top-1)
                        call move_alloc(tmp_stack, stack)
                    end if
                    stack(:,top) = [ii,jj,kk]
                end if
            end if
        end do; end do; end do

        do while (top > 0)
            ci = stack(1,top); cj = stack(2,top); ck = stack(3,top)
            top = top - 1
            do i_dir = 1, n_dirs
                ni = ci + dirs(1,i_dir)
                nj = cj + dirs(2,i_dir)
                nk = ck + dirs(3,i_dir)
                if (ni>=sx .and. ni<=ex .and. nj>=sy .and. nj<=ey .and. nk>=sz .and. nk<=ez) then
                    if (flag(ni,nj,nk) == 0 .and. &
                        abs(grid(ni,nj,nk)) >= scalarvalue * SENTINEL_FRAC) then
                        flag(ni,nj,nk) = 1
                        top = top + 1
                        if (top > stack_size) then
                            stack_size = stack_size * 2
                            allocate(tmp_stack(3, stack_size))
                            tmp_stack(:, 1:top-1) = stack(:, 1:top-1)
                            call move_alloc(tmp_stack, stack)
                        end if
                        stack(:,top) = [ni,nj,nk]
                    end if
                end if
            end do
        end do

        !
        ! Pass 2: any sentinel cell not reached by the exterior BFS is enclosed
        ! by the narrowband surface, hence interior. This does not depend on the
        ! sign of narrowband values, so it is robust to inconsistent surface
        ! normals.
        !
        do kk = sz, ez; do jj = sy, ey; do ii = sx, ex
            if (flag(ii,jj,kk) == 0 .and. &
                abs(grid(ii,jj,kk)) >= scalarvalue * SENTINEL_FRAC) &
                flag(ii,jj,kk) = 2
        end do; end do; end do

        !
        ! Mark interior sentinel cells
        !
        do kk = sz, ez; do jj = sy, ey; do ii = sx, ex
            if (flag(ii,jj,kk) == 2 .and. grid(ii,jj,kk) > 0.0_dp) &
                grid(ii,jj,kk) = large_negative_value
        end do; end do; end do

        deallocate(flag, stack)
    end subroutine fill_internal

end module flood_fill_mod
