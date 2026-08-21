module spatial_hash_mod
    !
    ! spatial_hash_mod — uniform spatial hash over the triangle set for fast triangle queries.
    !
    use utils_io, only : dp
    implicit none

    ! Hash cell size h = hash_c * max(dx, dy, dz_min)
    real(dp), parameter :: hash_c = 6.0_dp

    ! Hash grid parameters (set by build_hash)
    integer  :: hcells(3)        ! number of hash cells in each direction
    real(dp) :: h_origin(3)      ! lower-left corner of hash domain (bbox_min - margin)
    real(dp) :: h_size(3)        ! hash cell size in each direction
    integer  :: total_hcells     ! hcells(1)*hcells(2)*hcells(3)

    ! CSR sparse storage: cell -> triangle list
    integer, allocatable :: csr_offset(:)   ! length total_hcells+1; offset into csr_triid
    integer, allocatable :: csr_triid(:)    ! triangle IDs, packed per cell

    ! Max triangles in any single hash cell (set by build_hash; used to bound query buffer)
    integer :: max_cell_count
    ! Number of faces stored in the hash (set by build_hash)
    integer :: nfaces_hash

    ! Per-thread GPU query_hash candidate buffer capacity
    integer, parameter :: max_cands_gpu = 32768

    !$acc declare create(hcells, h_origin, h_size, total_hcells)
    ! csr_offset/csr_triid are used inside query_hash_gpu (!$acc routine seq), so must be declared, not just enter-data'd.
    !$acc declare create(csr_offset, csr_triid)

contains

    !
    !  Inline helper: flat index of hash cell (hx,hy,hz)  (1-based)
    !
    pure integer function cell_flat(hx, hy, hz)
        !$acc routine seq
        integer, intent(in) :: hx, hy, hz
        cell_flat = (hx-1)*hcells(2)*hcells(3) + (hy-1)*hcells(3) + hz
    end function cell_flat

    !
    !  Hash cell index for a world coordinate (clipped to valid range)
    !
    pure integer function hash_ix(coord, dim)
        !$acc routine seq
        real(dp), intent(in) :: coord
        integer,  intent(in) :: dim
        hash_ix = max(1, min(hcells(dim), int((coord - h_origin(dim)) / h_size(dim)) + 1))
    end function hash_ix

    ! 
    subroutine build_hash(vertices, faces, nfaces, &
                          bbox_min, bbox_max, dx_in, dy_in, dz_min)
        !
        ! Build CSR spatial hash for nfaces triangles.
        !
        implicit none
        real(dp), intent(in) :: vertices(3,*)
        integer,  intent(in) :: faces(3,*)
        integer,  intent(in) :: nfaces
        real(dp), intent(in) :: bbox_min(3), bbox_max(3)
        real(dp), intent(in) :: dx_in, dy_in, dz_min

        integer  :: fi, ci, d
        real(dp) :: h, tri_min(3), tri_max(3)
        integer  :: lo(3), hi(3), hx, hy, hz, cflat
        integer, allocatable :: count(:)     ! per-cell triangle count (pass 1)
        integer, allocatable :: cursor(:)    ! write cursor (pass 2)
        real(dp), allocatable :: tri_aabbs(:,:,:)  ! cached AABBs: (fi, 1=min/2=max, dim)

        ! Choose hash cell size
        h = hash_c * max(dx_in, dy_in, dz_min)
        h_size = h

        ! Hash domain: geometry AABB with one cell margin
        h_origin = bbox_min - h
        do d = 1, 3
            hcells(d) = max(1, ceiling((bbox_max(d) - bbox_min(d) + 2.0_dp*h) / h))
        end do
        total_hcells = hcells(1) * hcells(2) * hcells(3)

        allocate(count(total_hcells))
        count = 0

        ! Cache all triangle AABBs once (reused in pass 2; avoids calling tri_aabb twice)
        allocate(tri_aabbs(nfaces, 2, 3))

        !  Pass 1: count (cell, triangle) pairs 
        do fi = 1, nfaces
            call tri_aabb(vertices, faces(:,fi), tri_min, tri_max)
            tri_aabbs(fi, 1, :) = tri_min
            tri_aabbs(fi, 2, :) = tri_max
            lo(1) = hash_ix(tri_min(1), 1);  hi(1) = hash_ix(tri_max(1), 1)
            lo(2) = hash_ix(tri_min(2), 2);  hi(2) = hash_ix(tri_max(2), 2)
            lo(3) = hash_ix(tri_min(3), 3);  hi(3) = hash_ix(tri_max(3), 3)
            do hx = lo(1), hi(1)
                do hy = lo(2), hi(2)
                    do hz = lo(3), hi(3)
                        cflat = cell_flat(hx, hy, hz)
                        count(cflat) = count(cflat) + 1
                    end do
                end do
            end do
        end do

        !  Prefix sum to build CSR offsets 
        allocate(csr_offset(total_hcells+1))
        csr_offset(1) = 1
        do ci = 1, total_hcells
            csr_offset(ci+1) = csr_offset(ci) + count(ci)
        end do
        allocate(csr_triid(csr_offset(total_hcells+1)-1))

        ! Record max triangles per cell (used to bound query_hash buffer)
        max_cell_count = 0
        do ci = 1, total_hcells
            max_cell_count = max(max_cell_count, csr_offset(ci+1) - csr_offset(ci))
        end do
        nfaces_hash = nfaces

        !  Pass 2: write triangle IDs (reuse cached AABBs)
        allocate(cursor(total_hcells))
        cursor = csr_offset(1:total_hcells)   ! start of each cell's write window

        do fi = 1, nfaces
            tri_min = tri_aabbs(fi, 1, :)
            tri_max = tri_aabbs(fi, 2, :)
            lo(1) = hash_ix(tri_min(1), 1);  hi(1) = hash_ix(tri_max(1), 1)
            lo(2) = hash_ix(tri_min(2), 2);  hi(2) = hash_ix(tri_max(2), 2)
            lo(3) = hash_ix(tri_min(3), 3);  hi(3) = hash_ix(tri_max(3), 3)
            do hx = lo(1), hi(1)
                do hy = lo(2), hi(2)
                    do hz = lo(3), hi(3)
                        cflat = cell_flat(hx, hy, hz)
                        csr_triid(cursor(cflat)) = fi
                        cursor(cflat) = cursor(cflat) + 1
                    end do
                end do
            end do
        end do

        deallocate(count, cursor, tri_aabbs)

#ifdef GPU_SDF
        if (27 * max_cell_count > max_cands_gpu) &
            error stop "GPU_SDF: hash cell too dense, raise max_cands_gpu (spatial_hash.f90)"
        !$acc update device(hcells, h_origin, h_size, total_hcells)
        !$acc update device(csr_offset, csr_triid)
#endif
    end subroutine build_hash

    subroutine query_hash(pt, candidate_tris, n_candidates)
        !
        ! Return all triangle IDs in the 3x3x3 neighbourhood of hash cells around pt.
        ! candidate_tris must be pre-allocated by the caller to at least
        ! 27 * max_cell_count elements.  n_candidates is set to the actual count.
        !
        implicit none
        real(dp), intent(in)    :: pt(3)
        integer,  intent(inout) :: candidate_tris(:)
        integer,  intent(out)   :: n_candidates
        integer :: hx, hy, hz, cx, cy, cz, cflat, lo, hi
        integer :: n_tmp

        hx = hash_ix(pt(1), 1)
        hy = hash_ix(pt(2), 2)
        hz = hash_ix(pt(3), 3)

        n_tmp = 0
        do cx = max(1, hx-1), min(hcells(1), hx+1)
            do cy = max(1, hy-1), min(hcells(2), hy+1)
                do cz = max(1, hz-1), min(hcells(3), hz+1)
                    cflat = cell_flat(cx, cy, cz)
                    lo = csr_offset(cflat)
                    hi = csr_offset(cflat+1) - 1
                    if (hi >= lo) then
                        candidate_tris(n_tmp+1 : n_tmp+hi-lo+1) = csr_triid(lo:hi)
                        n_tmp = n_tmp + hi - lo + 1
                    end if
                end do
            end do
        end do
        n_candidates = n_tmp
    end subroutine query_hash

#ifdef GPU_SDF
    !  Device-callable query_hash: fixed-capacity dummy array, private per GPU thread.
    subroutine query_hash_gpu(pt, candidate_tris, n_candidates)
        !$acc routine seq
        implicit none
        real(dp), intent(in)  :: pt(3)
        integer,  intent(out) :: candidate_tris(max_cands_gpu)
        integer,  intent(out) :: n_candidates
        integer :: hx, hy, hz, cx, cy, cz, cflat, lo, hi, kk
        integer :: n_tmp

        hx = hash_ix(pt(1), 1)
        hy = hash_ix(pt(2), 2)
        hz = hash_ix(pt(3), 3)

        n_tmp = 0
        do cx = max(1, hx-1), min(hcells(1), hx+1)
            do cy = max(1, hy-1), min(hcells(2), hy+1)
                do cz = max(1, hz-1), min(hcells(3), hz+1)
                    cflat = cell_flat(cx, cy, cz)
                    lo = csr_offset(cflat)
                    hi = csr_offset(cflat+1) - 1
                    do kk = lo, hi
                        n_tmp = n_tmp + 1
                        candidate_tris(n_tmp) = csr_triid(kk)
                    end do
                end do
            end do
        end do
        n_candidates = n_tmp
    end subroutine query_hash_gpu
#endif

    !
    subroutine free_hash()
        if (allocated(csr_offset)) deallocate(csr_offset)
        if (allocated(csr_triid))  deallocate(csr_triid)
    end subroutine free_hash

    ! 
    !  Internal helper: AABB of a single triangle
    ! 
    subroutine tri_aabb(vertices, tri_vids, tmin, tmax)
        implicit none
        real(dp), intent(in)  :: vertices(3,*)
        integer,  intent(in)  :: tri_vids(3)
        real(dp), intent(out) :: tmin(3), tmax(3)
        integer :: d
        do d = 1, 3
            tmin(d) = min(vertices(d,tri_vids(1)), vertices(d,tri_vids(2)), vertices(d,tri_vids(3)))
            tmax(d) = max(vertices(d,tri_vids(1)), vertices(d,tri_vids(2)), vertices(d,tri_vids(3)))
        end do
    end subroutine tri_aabb

end module spatial_hash_mod

