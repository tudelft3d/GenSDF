module narrowband_mod
    !
    ! narrowband_mod — triangle distance primitives and SDF computation.
    !
    use utils_io
#ifdef GPU_SDF
    use spatial_hash_mod, only : query_hash, max_cell_count, query_hash_gpu, max_cands_gpu
#else
    use spatial_hash_mod, only : query_hash, max_cell_count
#endif
    implicit none

contains

    ! -----------------------------------------------------------------------
    subroutine cross_product(res, u, v)
        !$acc routine seq
        implicit none
        real(dp), dimension(3), intent(in)  :: u, v
        real(dp), dimension(3), intent(out) :: res
        res(1) = u(2)*v(3) - u(3)*v(2)
        res(2) = u(3)*v(1) - u(1)*v(3)
        res(3) = u(1)*v(2) - u(2)*v(1)
    end subroutine cross_product

    ! -----------------------------------------------------------------------
    subroutine project_point_to_triangle(point, vert0, vert1, vert2, param_s, param_t)
        !$acc routine seq
        !
        ! Closest-point projection of point onto the triangle (vert0, vert1, vert2).
        ! Uses Ericson's 7-region algorithm (Real-Time Collision Detection §5.1.5),
        ! which correctly handles all Voronoi regions: 3 vertices, 3 edges, 1 face.
        !
        ! Returns (param_s, param_t) such that:
        !   closest_point = vert0 + param_s*(vert1-vert0) + param_t*(vert2-vert0)
        !   with param_s >= 0, param_t >= 0, param_s + param_t <= 1.
        !
        implicit none
        real(dp), intent(in)  :: point(3), vert0(3), vert1(3), vert2(3)
        real(dp), intent(out) :: param_s, param_t
        real(dp) :: ab(3), ac(3), ap(3)
        real(dp) :: d1, d2, d3, d4, d5, d6
        real(dp) :: va, vb, vc, denom

        ab = vert1 - vert0
        ac = vert2 - vert0
        ap = point  - vert0

        d1 = dot_product(ab, ap)
        d2 = dot_product(ac, ap)

        ! Vertex region of vert0
        if (d1 <= 0.0_dp .and. d2 <= 0.0_dp) then
            param_s = 0.0_dp; param_t = 0.0_dp; return
        end if

        d3 = dot_product(ab, ap - ab)
        d4 = dot_product(ac, ap - ab)

        ! Vertex region of vert1
        if (d3 >= 0.0_dp .and. d4 <= d3) then
            param_s = 1.0_dp; param_t = 0.0_dp; return
        end if

        ! Edge region of vert0-vert1 (t = 0)
        vc = d1*d4 - d3*d2
        if (vc <= 0.0_dp .and. d1 >= 0.0_dp .and. d3 <= 0.0_dp) then
            param_s = d1 / (d1 - d3)   ! d1-d3 = dot(ab,ab) > 0
            param_t = 0.0_dp; return
        end if

        d5 = dot_product(ab, ap - ac)
        d6 = dot_product(ac, ap - ac)

        ! Vertex region of vert2
        if (d6 >= 0.0_dp .and. d5 <= d6) then
            param_s = 0.0_dp; param_t = 1.0_dp; return
        end if

        ! Edge region of vert0-vert2 (s = 0)
        vb = d5*d2 - d1*d6
        if (vb <= 0.0_dp .and. d2 >= 0.0_dp .and. d6 <= 0.0_dp) then
            param_t = d2 / (d2 - d6)   ! d2-d6 = dot(ac,ac) > 0
            param_s = 0.0_dp; return
        end if

        ! Edge region of vert1-vert2 (s + t = 1)
        va = d3*d6 - d5*d4
        if (va <= 0.0_dp .and. (d4 - d3) >= 0.0_dp .and. (d5 - d6) >= 0.0_dp) then
            ! denominator = |V2-V1|^2, always > 0 when triangle is non-degenerate
            denom = (d4 - d3) + (d5 - d6)
            if (denom > 0.0_dp) then
                param_t = (d4 - d3) / denom
            else
                param_t = 0.5_dp
            end if
            param_s = 1.0_dp - param_t; return
        end if

        ! Interior (face) region
        denom = va + vb + vc
        if (abs(denom) < 1.0e-30_dp) then
            ! Degenerate triangle — project to vert0
            param_s = 0.0_dp; param_t = 0.0_dp; return
        end if
        denom = 1.0_dp / denom
        param_s = vb * denom
        param_t = vc * denom
    end subroutine project_point_to_triangle

    ! -----------------------------------------------------------------------
    subroutine distance_point_to_triangle(point, vert0, vert1, vert2, avg_normal, dist_to_face)
        !
        ! Signed distance from point to nearest point on triangle.
        ! Sign determined by dot product with the effective outward normal (positive = outside).
        !
        !$acc routine seq
        implicit none
        real(dp), intent(in)  :: point(3), vert0(3), vert1(3), vert2(3), avg_normal(3)
        real(dp), intent(out) :: dist_to_face
        real(dp) :: edge0(3), edge1(3), dist_vec(3), proj(3)
        real(dp) :: param_s, param_t, sign_indicator
        real(dp) :: geom_normal(3), geom_mag, sign_normal(3)

        edge0 = vert1 - vert0
        edge1 = vert2 - vert0

        ! Geometric (winding-based) face normal — no dependency on stored normals.
        call cross_product(geom_normal, edge0, edge1)
        geom_mag = sqrt(dot_product(geom_normal, geom_normal))

        if (geom_mag > 1.0e-14_dp) then
            geom_normal = geom_normal / geom_mag
            ! Use avg_normal only when it agrees with the winding-derived direction.
            if (dot_product(avg_normal, geom_normal) >= 0.0_dp) then
                sign_normal = avg_normal   ! consistent: avg_normal is reliable
            else
                sign_normal = geom_normal  ! stored normals are flipped → use geometric
            end if
        else
            ! Degenerate triangle: best-effort with whatever was passed in.
            sign_normal = avg_normal
        end if

        call project_point_to_triangle(point, vert0, vert1, vert2, param_s, param_t)
        proj = vert0 + param_s*edge0 + param_t*edge1
        dist_vec = point - proj
        sign_indicator = sign(1.0_dp, dot_product(dist_vec, sign_normal))
        dist_to_face   = sign_indicator * sqrt(dot_product(dist_vec, dist_vec))
    end subroutine distance_point_to_triangle

    ! -----------------------------------------------------------------------
    subroutine compute_scalar_distance_face(procid, xstart, xend, ystart, yend, zstart, zend, &
                                             xin, yin, zin, nfaces, input_faces, input_face_normals, &
                                             input_vertices, input_normals, buffer_point_size, distance)
        !
        ! Original O(Nf × buffer^3) brute-force distance loop.
        !
        use mpi, only : MPI_WTIME
        implicit none
        integer,  intent(in)  :: procid, xstart, xend, ystart, yend, zstart, zend
        real(dp), intent(in),  dimension(:) :: xin, yin, zin
        integer,  intent(in)  :: nfaces, buffer_point_size
        real(dp), intent(in),  dimension(:,:) :: input_vertices, input_normals
        integer,  intent(in),  dimension(:,:) :: input_faces, input_face_normals
        real(dp), intent(out) :: distance(xstart:xend, 1:size(yin), 1:size(zin))
        integer  :: face_id, ii, jj, kk
        real(dp) :: vertex_1(3), vertex_2(3), vertex_3(3)
        real(dp) :: norm_1(3), norm_2(3), norm_3(3)
        real(dp) :: min_query(3), max_query(3), query_point(3)
        integer  :: min_index(3), max_index(3)
        integer  :: ii_lo, ii_hi, jj_lo, jj_hi, kk_lo, kk_hi
        real(dp) :: avg_normal(3), avg_normal_mag_inv
        real(dp) :: temp_distance
        real(dp) :: deltax_inverse, deltay_inverse, stime
        integer  :: printrank

        printrank = nprocs / 2
        stime = MPI_WTIME()
        deltax_inverse = 1.0_dp / dx
        deltay_inverse = 1.0_dp / dy

#ifdef GPU_SDF
        !$acc data copy(distance) copyin(xin, yin, zin) &
        !$acc      present(input_vertices, input_normals, input_faces, input_face_normals)
#endif
        do face_id = 1, nfaces
            if (procid == printrank) call show_progress(face_id, nfaces, pbarwidth, stime)
            vertex_1 = input_vertices(:, input_faces(1,face_id))
            vertex_2 = input_vertices(:, input_faces(2,face_id))
            vertex_3 = input_vertices(:, input_faces(3,face_id))
            norm_1   = input_normals(:, input_face_normals(1,face_id))
            norm_2   = input_normals(:, input_face_normals(2,face_id))
            norm_3   = input_normals(:, input_face_normals(3,face_id))
            avg_normal = norm_1 + norm_2 + norm_3
            if (dot_product(avg_normal, avg_normal) > 1.0e-30_dp) then
                avg_normal_mag_inv = 1.0_dp / sqrt(dot_product(avg_normal, avg_normal))
                avg_normal = avg_normal * avg_normal_mag_inv
            else
                ! Normals cancel at sharp crease: fall back to geometric (winding) normal
                avg_normal(1) = (vertex_2(2)-vertex_1(2))*(vertex_3(3)-vertex_1(3)) &
                              - (vertex_2(3)-vertex_1(3))*(vertex_3(2)-vertex_1(2))
                avg_normal(2) = (vertex_2(3)-vertex_1(3))*(vertex_3(1)-vertex_1(1)) &
                              - (vertex_2(1)-vertex_1(1))*(vertex_3(3)-vertex_1(3))
                avg_normal(3) = (vertex_2(1)-vertex_1(1))*(vertex_3(2)-vertex_1(2)) &
                              - (vertex_2(2)-vertex_1(2))*(vertex_3(1)-vertex_1(1))
                avg_normal_mag_inv = 1.0_dp / sqrt(dot_product(avg_normal, avg_normal))
                avg_normal = avg_normal * avg_normal_mag_inv
            end if
            do ii = 1, 3
                min_query(ii) = min(vertex_1(ii), vertex_2(ii), vertex_3(ii))
                max_query(ii) = max(vertex_1(ii), vertex_2(ii), vertex_3(ii))
            end do
            min_index(1) = floor((min_query(1) - xin(1) + 0.5_dp*dx) * deltax_inverse) + 1 - buffer_point_size
            max_index(1) = floor((max_query(1) - xin(1) + 0.5_dp*dx) * deltax_inverse) + 1 + buffer_point_size
            min_index(2) = floor((min_query(2) - yin(1) + 0.5_dp*dy) * deltay_inverse) + 1 - buffer_point_size
            max_index(2) = floor((max_query(2) - yin(1) + 0.5_dp*dy) * deltay_inverse) + 1 + buffer_point_size
            min_index(3) = max(bsearch_le(zin, min_query(3), zstart, zend) - buffer_point_size, zstart)
            max_index(3) = min(bsearch_le(zin, max_query(3), zstart, zend) + buffer_point_size, zend)
            ii_lo = min_index(1); ii_hi = max_index(1)
            jj_lo = min_index(2); jj_hi = max_index(2)
            kk_lo = min_index(3); kk_hi = max_index(3)

#ifdef GPU_SDF
            ! Scalar loop bounds, not array elements — nvfortran mishandles firstprivate arrays used as collapsed loop bounds.
            !$acc parallel loop collapse(3) &
            !$acc    private(query_point, temp_distance) &
            !$acc    firstprivate(vertex_1, vertex_2, vertex_3, avg_normal, ii_lo, ii_hi, jj_lo, jj_hi, kk_lo, kk_hi) &
            !$acc    present(xin, yin, zin, distance)
#endif
            do kk = kk_lo, kk_hi
                do jj = jj_lo, jj_hi
                    do ii = ii_lo, ii_hi
                        if (ii>=xstart .and. ii<=xend .and. &
                            jj>=ystart .and. jj<=yend .and. &
                            kk>=zstart .and. kk<=zend) then
                            query_point = [xin(ii), yin(jj), zin(kk)]
                            call distance_point_to_triangle(query_point, vertex_1, vertex_2, vertex_3, &
                                                             avg_normal, temp_distance)
                            ! Equivalent to minloc(abs([temp_distance,distance])) — nvfortran mislowers minloc in this kernel.
                            if (abs(temp_distance) < abs(distance(ii,jj,kk))) distance(ii,jj,kk) = temp_distance
                        end if
                    end do
                end do
            end do
        end do
#ifdef GPU_SDF
        !$acc end data
#endif
    end subroutine compute_scalar_distance_face

    ! -----------------------------------------------------------------------
    subroutine compute_narrowband_sdf(procid, xstart, xend, ystart, yend, zstart, zend, &
                                       xin, yin, zin, nfaces, input_faces, input_face_normals, &
                                       input_vertices, input_normals, nbw, distance)
        !
        ! Narrow-band SDF computation using the spatial hash for candidate triangle queries.
        ! For each cell, we query the hash to get nearby triangles and compute the exact
        ! point-to-triangle distance, accepting the result only if it's within the narrow band.
        !
        use mpi, only : MPI_WTIME
        implicit none
        integer,  intent(in)  :: procid, xstart, xend, ystart, yend, zstart, zend
        real(dp), intent(in),  dimension(:) :: xin, yin, zin
        integer,  intent(in)  :: nfaces, nbw
        real(dp), intent(in),  dimension(:,:) :: input_vertices, input_normals
        integer,  intent(in),  dimension(:,:) :: input_faces, input_face_normals
        real(dp), intent(out) :: distance(xstart:xend, 1:size(yin), 1:size(zin))

        integer  :: ii, jj, kk, fc, fid, n_cand
#ifdef GPU_SDF
        integer  :: cands(max_cands_gpu)
#else
        integer, allocatable :: cands(:)
#endif
        real(dp) :: query_point(3)
        real(dp) :: v1(3), v2(3), v3(3)
        real(dp) :: n1(3), n2(3), n3(3), avg_n(3), avg_n_inv
        real(dp) :: d_try, best_d, nbw_dist
        real(dp) :: stime
        integer  :: printrank

        printrank = nprocs / 2
        stime     = MPI_WTIME()

        ! Initialise all cells to scalarvalue (large positive = far fluid)
        distance = scalarvalue

        ! Narrow-band distance estimate: nbw cells * max grid spacing
        nbw_dist = real(nbw, dp) * max(dx, dy, maxval(dz))

#ifdef GPU_SDF
        !$acc data copyin(xin, yin, zin) copy(distance) &
        !$acc      present(input_vertices, input_normals, input_faces, input_face_normals)
        !$acc parallel loop collapse(3) &
        !$acc    private(query_point, v1, v2, v3, n1, n2, n3, avg_n, avg_n_inv, d_try, best_d, &
        !$acc            fc, fid, n_cand, cands) default(none) &
        !$acc    present(xin, yin, zin, distance, input_vertices, input_normals, input_faces, input_face_normals) &
        !$acc    firstprivate(nbw_dist)
        do kk = zstart, zend
            do jj = ystart, yend
                do ii = xstart, xend
                    query_point = [xin(ii), yin(jj), zin(kk)]
                    call query_hash_gpu(query_point, cands, n_cand)

                    best_d = scalarvalue
                    do fc = 1, n_cand
                        fid = cands(fc)
                        v1 = input_vertices(:, input_faces(1,fid))
                        v2 = input_vertices(:, input_faces(2,fid))
                        v3 = input_vertices(:, input_faces(3,fid))
                        n1 = input_normals(:, input_face_normals(1,fid))
                        n2 = input_normals(:, input_face_normals(2,fid))
                        n3 = input_normals(:, input_face_normals(3,fid))
                        avg_n     = n1 + n2 + n3
                        if (dot_product(avg_n, avg_n) > 1.0e-30_dp) then
                            avg_n_inv = 1.0_dp / sqrt(dot_product(avg_n, avg_n))
                            avg_n     = avg_n * avg_n_inv
                        else
                            ! Normals cancel at sharp crease: fall back to geometric (winding) normal
                            avg_n(1) = (v2(2)-v1(2))*(v3(3)-v1(3)) - (v2(3)-v1(3))*(v3(2)-v1(2))
                            avg_n(2) = (v2(3)-v1(3))*(v3(1)-v1(1)) - (v2(1)-v1(1))*(v3(3)-v1(3))
                            avg_n(3) = (v2(1)-v1(1))*(v3(2)-v1(2)) - (v2(2)-v1(2))*(v3(1)-v1(1))
                            avg_n_inv = 1.0_dp / sqrt(dot_product(avg_n, avg_n))
                            avg_n     = avg_n * avg_n_inv
                        end if
                        call distance_point_to_triangle(query_point, v1, v2, v3, avg_n, d_try)
                        if (abs(d_try) < abs(best_d)) best_d = d_try
                    end do

                    ! Accept only cells inside the narrow band.
                    if (abs(best_d) < nbw_dist) distance(ii, jj, kk) = best_d
                end do
            end do
        end do
        !$acc end data
#else
        ! Pre-allocate candidate buffer once (27 hash-cell neighbourhood max)
        allocate(cands(27 * max(1, max_cell_count)))

        do kk = zstart, zend
            do jj = ystart, yend
                if (procid == printrank .and. jj == ystart) &
                    call show_progress(kk - zstart + 1, zend - zstart + 1, pbarwidth, stime)
                do ii = xstart, xend
                    query_point = [xin(ii), yin(jj), zin(kk)]
                    call query_hash(query_point, cands, n_cand)

                    best_d = scalarvalue
                    do fc = 1, n_cand
                        fid = cands(fc)
                        v1 = input_vertices(:, input_faces(1,fid))
                        v2 = input_vertices(:, input_faces(2,fid))
                        v3 = input_vertices(:, input_faces(3,fid))
                        n1 = input_normals(:, input_face_normals(1,fid))
                        n2 = input_normals(:, input_face_normals(2,fid))
                        n3 = input_normals(:, input_face_normals(3,fid))
                        avg_n     = n1 + n2 + n3
                        if (dot_product(avg_n, avg_n) > 1.0e-30_dp) then
                            avg_n_inv = 1.0_dp / sqrt(dot_product(avg_n, avg_n))
                            avg_n     = avg_n * avg_n_inv
                        else
                            ! Normals cancel at sharp crease: fall back to geometric (winding) normal
                            avg_n(1) = (v2(2)-v1(2))*(v3(3)-v1(3)) - (v2(3)-v1(3))*(v3(2)-v1(2))
                            avg_n(2) = (v2(3)-v1(3))*(v3(1)-v1(1)) - (v2(1)-v1(1))*(v3(3)-v1(3))
                            avg_n(3) = (v2(1)-v1(1))*(v3(2)-v1(2)) - (v2(2)-v1(2))*(v3(1)-v1(1))
                            avg_n_inv = 1.0_dp / sqrt(dot_product(avg_n, avg_n))
                            avg_n     = avg_n * avg_n_inv
                        end if
                        call distance_point_to_triangle(query_point, v1, v2, v3, avg_n, d_try)
                        if (abs(d_try) < abs(best_d)) best_d = d_try
                    end do
                    ! no deallocate — cands is pre-allocated across iterations

                    ! Accept only cells inside the narrow band.
                    if (abs(best_d) < nbw_dist) distance(ii, jj, kk) = best_d
                end do
            end do
        end do
        deallocate(cands)
#endif
    end subroutine compute_narrowband_sdf

    ! -----------------------------------------------------------------------
    ! Binary search: return largest index idx in [lo,hi] s.t. arr(idx) <= val.
    ! Returns lo-1 if all arr(lo:hi) > val.
    pure integer function bsearch_le(arr, val, lo, hi)
        real(dp), intent(in) :: arr(:), val
        integer,  intent(in) :: lo, hi
        integer :: l, r, mid
        l = lo; r = hi; bsearch_le = lo - 1
        do while (l <= r)
            mid = (l + r) / 2
            if (arr(mid) <= val) then
                bsearch_le = mid
                l = mid + 1
            else
                r = mid - 1
            end if
        end do
    end function bsearch_le

end module narrowband_mod
