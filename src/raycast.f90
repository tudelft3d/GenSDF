module raycast_mod
    !
    ! raycast_mod — winding-independent, per-solid point-in-solid classification
    ! via triangle ray-crossing parity, used as the primary interior/exterior
    ! decision in place of (a) "nearest triangle wins" sign propagation and
    ! (b) flood_fill_mod's topological BFS.
    !
    ! Both of the replaced methods are fragile for real multi-object scenes:
    ! nearest-triangle sign silently hands a point to whichever solid happens
    ! to have the closest face, even an unrelated one (e.g. a small ground
    ! plane sitting under an open-bottomed building can "steal" the sign of
    ! points that are genuinely inside the building but happen to be nearer
    ! to the ground plane's top face); the topological BFS requires a fully
    ! closed shell of near-surface cells around every solid, which an open
    ! mesh (deliberately missing a bottom cap, say) does not guarantee.
    !
    ! Ray-crossing parity avoids both: it only asks "how many of THIS solid's
    ! own triangles does a ray from this point cross", which is correct
    ! regardless of what other geometry is nearby and regardless of triangle
    ! winding. A single axis is not enough on its own (a ray parallel to a
    ! missing face, e.g. straight up through an open-bottomed box, undercounts
    ! by exactly the missing crossing), so this casts along all three axes and
    ! takes a majority vote (>=2 of 3 agree), which is robust as long as a
    ! solid isn't degenerate along at least two of the three axes.
    !
    use utils_io,          only : dp
    use spatial_hash_mod,  only : hcells, csr_offset, csr_triid, cell_flat, hash_ix
    implicit none

contains

    ! -----------------------------------------------------------------------
    ! Gather every triangle ID whose hash bin lies in the 3x3 neighbourhood
    ! (in the two axes perpendicular to `axis`) of (coord_b, coord_c), across
    ! the full hash-cell range along `axis` — i.e. every candidate that could
    ! be crossed by the infinite line through (coord_b, coord_c) running
    ! along `axis`. axis/dim_b/dim_c follow the 1=x, 2=y, 3=z convention.
    subroutine query_hash_line(axis, coord_b, coord_c, candidate_tris, n_candidates)
        implicit none
        integer,  intent(in)    :: axis
        real(dp), intent(in)    :: coord_b, coord_c
        integer,  allocatable, intent(inout) :: candidate_tris(:)
        integer,  intent(out)   :: n_candidates
        integer :: dim_b, dim_c, hb, hc, cb, cc, i_along, cflat, lo, hi, n_tmp, need

        call perp_dims(axis, dim_b, dim_c)
        hb = hash_ix(coord_b, dim_b)
        hc = hash_ix(coord_c, dim_c)

        if (.not. allocated(candidate_tris)) allocate(candidate_tris(256))
        n_tmp = 0
        do cb = max(1, hb-1), min(hcells(dim_b), hb+1)
            do cc = max(1, hc-1), min(hcells(dim_c), hc+1)
                do i_along = 1, hcells(axis)
                    cflat = hcell_flat(axis, i_along, cb, cc)
                    lo = csr_offset(cflat); hi = csr_offset(cflat+1) - 1
                    if (hi >= lo) then
                        need = n_tmp + (hi - lo + 1)
                        if (need > size(candidate_tris)) call grow_int(candidate_tris, need)
                        candidate_tris(n_tmp+1:need) = csr_triid(lo:hi)
                        n_tmp = need
                    end if
                end do
            end do
        end do
        n_candidates = n_tmp
    end subroutine query_hash_line

    pure subroutine perp_dims(axis, dim_b, dim_c)
        integer, intent(in)  :: axis
        integer, intent(out) :: dim_b, dim_c
        select case (axis)
        case (1); dim_b = 2; dim_c = 3
        case (2); dim_b = 1; dim_c = 3
        case default; dim_b = 1; dim_c = 2
        end select
    end subroutine perp_dims

    pure integer function hcell_flat(axis, i_along, jb, jc)
        integer, intent(in) :: axis, i_along, jb, jc
        select case (axis)
        case (1); hcell_flat = cell_flat(i_along, jb, jc)
        case (2); hcell_flat = cell_flat(jb, i_along, jc)
        case default; hcell_flat = cell_flat(jb, jc, i_along)
        end select
    end function hcell_flat

    subroutine grow_int(arr, need)
        integer, allocatable, intent(inout) :: arr(:)
        integer, intent(in) :: need
        integer, allocatable :: tmp(:)
        allocate(tmp(max(need, 2*size(arr))))
        tmp(1:size(arr)) = arr
        call move_alloc(tmp, arr)
    end subroutine grow_int

    ! -----------------------------------------------------------------------
    ! Ray/triangle crossing test along `axis`, through (coord_b, coord_c) in
    ! the two perpendicular axes. Winding-independent: uses only vertex
    ! positions, never the stored/averaged triangle normal. Returns whether
    ! the infinite line crosses the triangle and, if so, the coordinate along
    ! `axis` at the crossing.
    pure subroutine ray_triangle_crossing(axis, coord_b, coord_c, v0, v1, v2, hit, coord_a)
        implicit none
        integer,  intent(in)  :: axis
        real(dp), intent(in)  :: coord_b, coord_c
        real(dp), intent(in)  :: v0(3), v1(3), v2(3)
        logical,  intent(out) :: hit
        real(dp), intent(out) :: coord_a
        integer  :: dim_b, dim_c
        real(dp) :: b0, c0, e1b, e1c, e2b, e2c, pb, pc
        real(dp) :: denom, s, t, u
        real(dp) :: edge1(3), edge2(3), nrm(3)

        hit = .false.
        coord_a = 0.0_dp
        call perp_dims(axis, dim_b, dim_c)

        b0 = v0(dim_b); c0 = v0(dim_c)
        e1b = v1(dim_b) - b0; e1c = v1(dim_c) - c0
        e2b = v2(dim_b) - b0; e2c = v2(dim_c) - c0
        pb  = coord_b - b0;   pc  = coord_c - c0

        ! 2-D barycentric test of (coord_b,coord_c) against the projected triangle
        denom = e1b*e2c - e2b*e1c
        if (abs(denom) < 1.0e-13_dp) return   ! projects to a degenerate line/point

        s = (pb*e2c - e2b*pc) / denom
        t = (e1b*pc - pb*e1c) / denom
        u = 1.0_dp - s - t
        if (s < 0.0_dp .or. t < 0.0_dp .or. u < 0.0_dp) return

        edge1 = v1 - v0
        edge2 = v2 - v0
        nrm(1) = edge1(2)*edge2(3) - edge1(3)*edge2(2)
        nrm(2) = edge1(3)*edge2(1) - edge1(1)*edge2(3)
        nrm(3) = edge1(1)*edge2(2) - edge1(2)*edge2(1)
        if (abs(nrm(axis)) < 1.0e-13_dp) return   ! triangle (near-)parallel to ray axis

        coord_a = v0(axis) - (nrm(dim_b)*(coord_b - v0(dim_b)) + nrm(dim_c)*(coord_c - v0(dim_c))) / nrm(axis)
        hit = .true.
    end subroutine ray_triangle_crossing

    ! -----------------------------------------------------------------------
    ! Classify every point (xin(sx:ex), yin(sy:ey), zin(sz:ez)) as inside or
    ! outside solid `sid`, via a 3-axis ray-crossing-parity majority vote.
    ! Processes one axis-perpendicular "line" of grid points at a time: for
    ! each line, gathers candidates once, buckets the crossing coordinates,
    ! sorts them, then sweeps the (already-sorted) grid coordinates along
    ! that line with a running pointer for an O(n log n + m) classification
    ! instead of re-testing every candidate at every grid point.
    subroutine classify_solid(vertices, faces, face_solid_id, nfaces, sid, &
                               xin, yin, zin, sx, ex, sy, ey, sz, ez, inside)
        implicit none
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: faces(:,:)
        integer,  intent(in) :: face_solid_id(:)
        integer,  intent(in) :: nfaces, sid
        real(dp), intent(in) :: xin(:), yin(:), zin(:)
        integer,  intent(in) :: sx, ex, sy, ey, sz, ez
        integer,  intent(out) :: inside(sx:ex, sy:ey, sz:ez)   ! vote count, 0-3

        integer, allocatable :: cands(:)
        real(dp), allocatable :: crossings(:)
        integer :: n_cand, n_cross, ii, jj, kk, fc, fid, cap
        real(dp) :: v0(3), v1(3), v2(3)
        logical  :: hit
        real(dp) :: coord_a
        real(dp) :: scale, eps_x, eps_y, eps_z
        integer, allocatable :: last_seen(:)
        integer :: query_gen

        inside = 0
        allocate(cands(256))

        ! build_hash bins a triangle into every hash cell its AABB overlaps
        ! (fine for nearest-distance search, since duplicates just get
        ! re-tested against the same running minimum). A crossing-parity
        ! count is not idempotent that way: a triangle whose AABB spans
        ! several hash cells within the queried neighbourhood — common for
        ! a big flat cap face relative to the hash cell size — would
        ! otherwise be counted as a crossing more than once per query,
        ! corrupting the parity. Tag each triangle with the query that last
        ! processed it so every candidate list is deduplicated in O(1) per
        ! candidate, without an O(nfaces) reset between columns.
        allocate(last_seen(nfaces))
        last_seen = 0
        query_gen = 0

        ! Deterministic per-axis jitter: a query point that lands exactly on
        ! the shared edge between two triangles (e.g. the diagonal of a
        ! quad face split into 2 triangles, common whenever a grid point
        ! coincides with a face's own centre — an axis-aligned box on an
        ! axis-aligned grid hits this constantly) would otherwise register
        ! as a crossing of BOTH triangles, doubling the count and flipping
        ! parity. Nudging the transverse coordinates by a tiny, distinct
        ! offset per axis makes an exact edge hit vanishingly unlikely
        ! without perturbing the classification of any real interior point.
        scale = max(abs(xin(ex)-xin(sx)), abs(yin(ey)-yin(sy)), abs(zin(ez)-zin(sz)), 1.0_dp)
        eps_x = scale * 7.5307258e-8_dp
        eps_y = scale * 3.1415927e-8_dp
        eps_z = scale * 4.1421356e-8_dp

        ! ---- axis = 3 (z): one line per (ii,jj), varying kk ----
        do ii = sx, ex
            do jj = sy, ey
                call query_hash_line(3, xin(ii)+eps_x, yin(jj)+eps_y, cands, n_cand)
                query_gen = query_gen + 1
                cap = max(n_cand, 1)
                allocate(crossings(cap))
                n_cross = 0
                do fc = 1, n_cand
                    fid = cands(fc)
                    if (last_seen(fid) == query_gen) cycle
                    last_seen(fid) = query_gen
                    if (face_solid_id(fid) /= sid) cycle
                    v0 = vertices(:, faces(1,fid)); v1 = vertices(:, faces(2,fid)); v2 = vertices(:, faces(3,fid))
                    call ray_triangle_crossing(3, xin(ii)+eps_x, yin(jj)+eps_y, v0, v1, v2, hit, coord_a)
                    if (hit) then
                        n_cross = n_cross + 1
                        crossings(n_cross) = coord_a
                    end if
                end do
                call sort_real(crossings, n_cross)
                call vote_along(crossings, n_cross, zin, sz, ez, inside(ii, jj, :))
                deallocate(crossings)
            end do
        end do

        ! ---- axis = 1 (x): one line per (jj,kk), varying ii ----
        do jj = sy, ey
            do kk = sz, ez
                call query_hash_line(1, yin(jj)+eps_y, zin(kk)+eps_z, cands, n_cand)
                query_gen = query_gen + 1
                cap = max(n_cand, 1)
                allocate(crossings(cap))
                n_cross = 0
                do fc = 1, n_cand
                    fid = cands(fc)
                    if (last_seen(fid) == query_gen) cycle
                    last_seen(fid) = query_gen
                    if (face_solid_id(fid) /= sid) cycle
                    v0 = vertices(:, faces(1,fid)); v1 = vertices(:, faces(2,fid)); v2 = vertices(:, faces(3,fid))
                    call ray_triangle_crossing(1, yin(jj)+eps_y, zin(kk)+eps_z, v0, v1, v2, hit, coord_a)
                    if (hit) then
                        n_cross = n_cross + 1
                        crossings(n_cross) = coord_a
                    end if
                end do
                call sort_real(crossings, n_cross)
                call vote_along(crossings, n_cross, xin, sx, ex, inside(:, jj, kk))
                deallocate(crossings)
            end do
        end do

        ! ---- axis = 2 (y): one line per (ii,kk), varying jj ----
        do ii = sx, ex
            do kk = sz, ez
                call query_hash_line(2, xin(ii)+eps_x, zin(kk)+eps_z, cands, n_cand)
                query_gen = query_gen + 1
                cap = max(n_cand, 1)
                allocate(crossings(cap))
                n_cross = 0
                do fc = 1, n_cand
                    fid = cands(fc)
                    if (last_seen(fid) == query_gen) cycle
                    last_seen(fid) = query_gen
                    if (face_solid_id(fid) /= sid) cycle
                    v0 = vertices(:, faces(1,fid)); v1 = vertices(:, faces(2,fid)); v2 = vertices(:, faces(3,fid))
                    call ray_triangle_crossing(2, xin(ii)+eps_x, zin(kk)+eps_z, v0, v1, v2, hit, coord_a)
                    if (hit) then
                        n_cross = n_cross + 1
                        crossings(n_cross) = coord_a
                    end if
                end do
                call sort_real(crossings, n_cross)
                call vote_along(crossings, n_cross, yin, sy, ey, inside(ii, :, kk))
                deallocate(crossings)
            end do
        end do

        deallocate(cands, last_seen)
    end subroutine classify_solid

    ! Given sorted crossing coordinates, sweep the (already-sorted) grid
    ! coordinates line(lo:hi) and accumulate a +1 vote wherever the parity of
    ! crossings strictly below that coordinate is odd (inside).
    subroutine vote_along(crossings, n_cross, line, lo, hi, votes)
        implicit none
        real(dp), intent(in)    :: crossings(:)
        integer,  intent(in)    :: n_cross, lo, hi
        real(dp), intent(in)    :: line(:)
        integer,  intent(inout) :: votes(lo:hi)
        integer :: idx, ptr

        ptr = 0
        do idx = lo, hi
            do while (ptr < n_cross)
                if (crossings(ptr+1) >= line(idx)) exit
                ptr = ptr + 1
            end do
            if (mod(ptr, 2) == 1) votes(idx) = votes(idx) + 1
        end do
    end subroutine vote_along

    ! -----------------------------------------------------------------------
    ! Top-level entry point: classify every point in [sx:ex,sy:ey,sz:ez] as
    ! solid/fluid via the 3-axis ray-crossing vote, one solid at a time (first
    ! solid to claim a point wins — solids are expected not to overlap), and
    ! write the result into `grid`/`objid_grid` (both full-domain (nx,ny,nz)
    ! arrays; only the [sx:ex,sy:ey,sz:ez] sub-box is touched, matching the
    ! AABB tagminmax already restricted the near-surface pass to). Cells that
    ! already carry a real near-surface distance (from compute_narrowband_sdf
    ! / compute_scalar_distance_face) keep their magnitude, with the sign
    ! forced to match the ray-cast vote; cells still at the sentinel value
    ! keep the sentinel magnitude (fast_sweep_3d fills those in afterwards).
    subroutine raycast_classify_domain(vertices, faces, face_solid_id, nfaces, nsolids, &
                                        xin, yin, zin, sx, ex, sy, ey, sz, ez, &
                                        scalarvalue_in, grid, objid_grid)
        implicit none
        real(dp), intent(in)    :: vertices(:,:)
        integer,  intent(in)    :: faces(:,:)
        integer,  intent(in)    :: face_solid_id(:)
        integer,  intent(in)    :: nfaces, nsolids
        real(dp), intent(in)    :: xin(:), yin(:), zin(:)
        integer,  intent(in)    :: sx, ex, sy, ey, sz, ez
        real(dp), intent(in)    :: scalarvalue_in
        real(dp), intent(inout) :: grid(:,:,:)
        real(dp), intent(inout) :: objid_grid(:,:,:)

        integer, allocatable :: votes(:,:,:)
        logical, allocatable :: any_inside(:,:,:)
        real(dp), allocatable :: winning_id(:,:,:)
        integer :: sid, ii, jj, kk

        allocate(votes(sx:ex, sy:ey, sz:ez))
        allocate(any_inside(sx:ex, sy:ey, sz:ez))
        allocate(winning_id(sx:ex, sy:ey, sz:ez))
        any_inside = .false.
        winning_id = 0.0_dp

        do sid = 1, nsolids
            call classify_solid(vertices, faces, face_solid_id, nfaces, sid, &
                                 xin, yin, zin, sx, ex, sy, ey, sz, ez, votes)
            do kk = sz, ez
                do jj = sy, ey
                    do ii = sx, ex
                        if (.not. any_inside(ii,jj,kk) .and. votes(ii,jj,kk) >= 2) then
                            any_inside(ii,jj,kk) = .true.
                            winning_id(ii,jj,kk) = real(sid, dp)
                        end if
                    end do
                end do
            end do
        end do

        do kk = sz, ez
            do jj = sy, ey
                do ii = sx, ex
                    if (any_inside(ii,jj,kk)) then
                        if (abs(grid(ii,jj,kk)) < scalarvalue_in * 0.5_dp) then
                            grid(ii,jj,kk) = -abs(grid(ii,jj,kk))
                        else
                            grid(ii,jj,kk) = -scalarvalue_in
                        end if
                        objid_grid(ii,jj,kk) = winning_id(ii,jj,kk)
                    else
                        if (abs(grid(ii,jj,kk)) < scalarvalue_in * 0.5_dp) then
                            grid(ii,jj,kk) = abs(grid(ii,jj,kk))
                        else
                            grid(ii,jj,kk) = scalarvalue_in
                        end if
                        objid_grid(ii,jj,kk) = 0.0_dp
                    end if
                end do
            end do
        end do

        deallocate(votes, any_inside, winning_id)
    end subroutine raycast_classify_domain

    ! Small-n insertion sort (crossing counts per line are typically tiny)
    pure subroutine sort_real(arr, n)
        real(dp), intent(inout) :: arr(:)
        integer,  intent(in)    :: n
        integer  :: i, j
        real(dp) :: key
        do i = 2, n
            key = arr(i); j = i - 1
            do while (j >= 1)
                if (arr(j) <= key) exit
                arr(j+1) = arr(j); j = j - 1
            end do
            arr(j+1) = key
        end do
    end subroutine sort_real

end module raycast_mod
