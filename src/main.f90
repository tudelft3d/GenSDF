program generatesdf
    !
    ! GenSDF — Signed-Distance-Field Generator
    !
    use utils_io
    use narrowband_mod,   only : compute_scalar_distance_face, compute_narrowband_sdf
    use raycast_mod,      only : raycast_classify_domain
    use spatial_hash_mod, only : build_hash, free_hash
    use fast_sweep_mod,   only : fast_sweep_3d
    use mpi
#ifdef GPU_SDF
    use openacc
#endif

    implicit none

    !  geometry
    integer :: nfaces, nvertices, nnormals, nsolids
    real(dp), allocatable :: vertices(:,:), normals(:,:)
    integer,  allocatable :: faces(:,:), face_normals(:,:)
    integer,  allocatable :: face_solid_id(:)
    real(dp) :: bbox_min(3), bbox_max(3)
    integer  :: sx, ex, sy, ey, sz, ez
    integer  :: cx_lo, cx_hi   ! per-rank x-slab intersected with the AABB+buffer [sx,ex]

    !  SDF arrays
    real(dp), allocatable :: sdf(:,:,:)
    real(dp), allocatable :: objid(:,:,:)
    real(dp), allocatable :: flood_fill_arr(:,:,:)
    real(dp), allocatable :: objid_fill_arr(:,:,:)
    real(dp), allocatable :: fsm_local(:,:,:)    ! per-rank slab for parallel FSM

    !  MPI neighbours for FSM 
    integer :: left_rank, right_rank

    !  MPI 
    integer :: myid, ierror, mpi_dx
    integer, allocatable :: decomp_x_start(:), decomp_x_end(:), decomp_size(:)

    !  timing 
    real(dp) :: startTime, endTime, totalTime, time1, time2

    !  misc
    integer :: ii
    logical :: debug = .false.

#ifdef GPU_SDF
    !  multi-GPU device binding
    integer :: node_comm, node_rank, num_gpus
#endif

    ! MPI init and domain decomposition
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid,   ierror)
    allocate(decomp_x_start(0:nprocs-1), decomp_x_end(0:nprocs-1), decomp_size(0:nprocs-1))

#ifdef GPU_SDF
    ! One GPU per rank (node-local rank mod GPU count) — see docs/external-lib-guidance/GenSDF_GPU_porting.md.
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_comm, ierror)
    call MPI_COMM_RANK(node_comm, node_rank, ierror)
    num_gpus = acc_get_num_devices(acc_device_nvidia)
    if (num_gpus < 1) error stop "GPU_SDF: no NVIDIA GPU visible to this rank"
    call acc_set_device_num(mod(node_rank, num_gpus), acc_device_nvidia)
    call MPI_COMM_FREE(node_comm, ierror)
#endif

    if (myid == 0) then
        startTime = 0.0_dp
        call cpu_time(startTime)
        call printlogo()
        print *, "*** Starting with", nprocs, "MPI ranks ***"
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    call read_inputfile(myid)
#ifdef GPU_SDF
    !$acc update device(scalarvalue)
#endif
    call read_cans_grid(myid, 'data/', dp, ng, r0, non_uniform_grid, xp, yp, zp, xf, yf, zf, dz, lx, ly, lz)
    allocate(dz_inverse(ng(3)))
    call setup_grid_spacing(myid, xf, yf, zf, ng(3), dx, dy, dz, dx_inverse, dy_inverse, dz_inverse)
    ! dy above is the uniform-grid spacing; with a stretched spanwise grid keep the smallest cell width as the scalar dy
    if (y_stretched) then
        dy = minval(dyv)
        dy_inverse = 1.0_dp / dy
    end if
#ifdef GPU_SDF
    ! dz is allocated via dummy-argument aliasing in read_cans_grid, so declare create's automatic on-allocate hook never fires for it.
    !$acc enter data create(dz)
    !$acc enter data create(dyv)
    !$acc update device(dx, dy, dz, dyv)
#endif

    ! Load geometry: a single .stl/.obj file (one solid, ID 1) or a .list manifest
    ! of multiple .stl/.obj files, each assigned a solid ID by its manifest order.
    call load_multi_object_geometry(myid, trim(inputfilename), vertices, normals, faces, face_normals, &
                                    face_solid_id, nvertices, nnormals, nfaces, nsolids)

    ! Solver clamps solid IDs into src/global.f90's max_ibm_objects-sized BC arrays; keep both in sync.
    if (nsolids > 15) then
        print *, 'FATAL: manifest has ', nsolids, ' solids, exceeds solver max_ibm_objects=15'
        error stop
    end if

    call check_object_resolution(myid, vertices, faces, face_solid_id, nfaces, nsolids)

#ifdef GPU_SDF
    ! Geometry is read-only and reused unchanged across the p/u/v/w passes below.
    !$acc enter data copyin(vertices, normals, faces, face_normals, face_solid_id)
#endif

    call getbbox(myid, vertices, nvertices, bbox_min, bbox_max)
    call tagminmax(myid, xf, yp, zp, bbox_min, bbox_max, nx, ny, nz, dx, dy, dz(2), buffer_points, &
                   sx, ex, sy, ey, sz, ez)

    !  domain decomposition (x-slab over the FULL domain, not just the AABB).
    !  Restricting the slabs to [sx,ex] used to make fast_sweep_3d structurally
    !  unable to reach x-cells outside the AABB+buffer box, leaving them stuck
    !  at the raw sentinel scalarvalue forever (even with use_fast_sweep=.true.)
    !  while the y/z directions — already using the full 1:ny/1:nz range below —
    !  got a smooth Eikonal-propagated field all the way to the domain edges.
    !  That asymmetry produced a hard scalarvalue-magnitude cliff in x right at
    !  the AABB boundary. Decomposing over [1,nx] and only restricting the
    !  actual near-surface computation (via cx_lo/cx_hi below) keeps that work
    !  scoped to the AABB while letting the FSM pass fill the rest smoothly.
    mpi_dx = ceiling(real(nx, dp) / real(nprocs, dp))
    decomp_x_start(0) = 1
    decomp_x_end(0)   = min(mpi_dx, nx)
    decomp_size(0)     = decomp_x_end(0) - decomp_x_start(0) + 1
    do ii = 1, nprocs-2
        decomp_x_start(ii) = decomp_x_end(ii-1) + 1
        decomp_x_end(ii)   = min(decomp_x_start(ii) + mpi_dx - 1, nx)
        decomp_size(ii)     = decomp_x_end(ii) - decomp_x_start(ii) + 1
    end do
    if (nprocs > 1) then
        decomp_x_start(nprocs-1) = decomp_x_end(nprocs-2) + 1
        decomp_x_end(nprocs-1)   = nx
        decomp_size(nprocs-1)     = decomp_x_end(nprocs-1) - decomp_x_start(nprocs-1) + 1
    end if

    ! This rank's x-slab intersected with the AABB+buffer box: the only range
    ! where near-surface distance computation actually needs to run. Empty
    ! (cx_lo > cx_hi) when this rank's slab lies entirely outside the AABB.
    cx_lo = max(decomp_x_start(myid), sx)
    cx_hi = min(decomp_x_end(myid),   ex)

    if (debug .and. myid == 0) then
        do ii = 0, nprocs-1
            print *, "PID:", ii, "xstart:", decomp_x_start(ii), "xend:", decomp_x_end(ii)
        end do
    end if

    if (myid == 0) then
        call cpu_time(time1)
        print *, "-- Pre-processing done in", time1-startTime, "s"
        call estimated_memoryusage(nprocs, nfaces, nvertices, nx, ny, nz)
        print *, "*** Calculating SDF | cell-centres ***"
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    left_rank  = myid - 1
    right_rank = myid + 1
    if (myid == 0)        left_rank  = -1      ! Neumann BC at left boundary
    if (myid == nprocs-1) right_rank = nprocs  ! Neumann BC at right boundary

    ! Build spatial hash once (shared geometry, reused for all stagger locations).
    ! Needed unconditionally now: raycast_classify_domain (the primary sign/objid
    ! source, replacing flood_fill_mod's topological BFS) relies on it too, not
    ! just the fast-sweep narrow-band path.
    call build_hash(vertices, faces, nfaces, bbox_min, bbox_max, dx, dy, minval(dz))

    ! Allocate the full-domain array on all ranks (needed for broadcast before FSM)
    ! Only sdfp (cell-center) is calculated by default, compute_face_sdf = .true. handles face SDF 
    allocate(flood_fill_arr(nx, ny, nz))
    flood_fill_arr = scalarvalue
    allocate(objid_fill_arr(nx, ny, nz))
    objid_fill_arr = 0.0_dp

    allocate(sdf(decomp_x_start(myid):decomp_x_end(myid), ny, nz))
    allocate(objid(decomp_x_start(myid):decomp_x_end(myid), ny, nz))
    sdf = scalarvalue
    objid = 0.0_dp
    if (cx_lo <= cx_hi) then
        if (use_fast_sweep) then
            call compute_narrowband_sdf(myid, cx_lo, cx_hi, &
                                        sy, ey, sz, ez, xp, yp, zp, nfaces, faces, face_normals, &
                                        vertices, normals, face_solid_id, narrow_band_width, &
                                        sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
        else
            call compute_scalar_distance_face(myid, cx_lo, cx_hi, &
                                              sy, ey, sz, ez, xp, yp, zp, nfaces, faces, face_normals, &
                                              vertices, normals, face_solid_id, buffer_points, &
                                              sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
        end if
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                      1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
    call gather_array(myid, nprocs, objid, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                      1, ny, 1, nz, 0.0_dp, objid_fill_arr, ierror)
    ! gather_array fills untouched (outside-AABB) cells with the module sentinel `scalarvalue`,
    ! not the requested background — those cells are always far from any solid, so zero them here.
    if (myid == 0) where (objid_fill_arr == scalarvalue) objid_fill_arr = 0.0_dp
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    if (myid == 0) &
        call raycast_classify_domain(vertices, faces, face_solid_id, nfaces, nsolids, &
                                     xp, yp, zp, sx, ex, sy, ey, sz, ez, scalarvalue, &
                                     flood_fill_arr, objid_fill_arr)
    if (use_fast_sweep) then
        ! Broadcast sign-corrected grid to all ranks, run FSM in parallel on x-slabs
        call MPI_BCAST(flood_fill_arr, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
        allocate(fsm_local(decomp_x_start(myid):decomp_x_end(myid), ny, nz))
        fsm_local = flood_fill_arr(decomp_x_start(myid):decomp_x_end(myid), :, :)
        call fast_sweep_3d(fsm_local, decomp_x_start(myid), decomp_x_end(myid), 1, ny, 1, nz, &
                           xp, yp, zp, myid, left_rank, right_rank)
        call gather_array(myid, nprocs, fsm_local, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                          1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
        deallocate(fsm_local)
        ! objid is not touched by the fast-sweep propagation — it only needs to be
        ! correct in the narrow-band interface cells where ghost points actually live.
    end if
    if (myid == 0) then
        call cpu_time(time1)
        print *, "*** Writing | cell-centres ***"
        call write_sdf_padded('data/sdfp.bin', flood_fill_arr)
        call write_sdf_padded('data/sdfp_objid.bin', objid_fill_arr)
        call cpu_time(time2)
        print *, "-- Write done in", time2-time1, "s | cell-centres"
        if (compute_face_sdf) print *, "*** Calculating SDF | u-faces ***"
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    ! Face-staggered SDFs — only computed when compute_face_sdf = .true.
    ! (objid is recomputed but not written for these passes: the solver's thermal
    ! ghost list is cell-centre-only, so face-staggered object IDs are unused.)
    if (compute_face_sdf) then

        ! U-faces
        sdf = scalarvalue
        if (cx_lo <= cx_hi) then
            if (use_fast_sweep) then
                call compute_narrowband_sdf(myid, cx_lo, cx_hi, &
                                            sy, ey, sz, ez, xf, yp, zp, nfaces, faces, face_normals, &
                                            vertices, normals, face_solid_id, narrow_band_width, &
                                            sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
            else
                call compute_scalar_distance_face(myid, cx_lo, cx_hi, &
                                                  sy, ey, sz, ez, xf, yp, zp, nfaces, faces, face_normals, &
                                                  vertices, normals, face_solid_id, buffer_points, &
                                                  sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
            end if
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                          1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        if (myid == 0) &
            call raycast_classify_domain(vertices, faces, face_solid_id, nfaces, nsolids, &
                                         xf, yp, zp, sx, ex, sy, ey, sz, ez, scalarvalue, &
                                         flood_fill_arr, objid_fill_arr)
        if (use_fast_sweep) then
            call MPI_BCAST(flood_fill_arr, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
            allocate(fsm_local(decomp_x_start(myid):decomp_x_end(myid), ny, nz))
            fsm_local = flood_fill_arr(decomp_x_start(myid):decomp_x_end(myid), :, :)
            call fast_sweep_3d(fsm_local, decomp_x_start(myid), decomp_x_end(myid), 1, ny, 1, nz, &
                               xf, yp, zp, myid, left_rank, right_rank)
            call gather_array(myid, nprocs, fsm_local, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                              1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
            deallocate(fsm_local)
        end if
        if (myid == 0) then
            call cpu_time(time1)
            print *, "*** Writing | u-faces ***"
            call write_sdf_padded('data/sdfu.bin', flood_fill_arr)
            call cpu_time(time2)
            print *, "-- Write done in", time2-time1, "s | u-faces"
            print *, "*** Calculating SDF | v-faces ***"
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)

        ! V-faces
        sdf = scalarvalue
        if (cx_lo <= cx_hi) then
            if (use_fast_sweep) then
                call compute_narrowband_sdf(myid, cx_lo, cx_hi, &
                                            sy, ey, sz, ez, xp, yf, zp, nfaces, faces, face_normals, &
                                            vertices, normals, face_solid_id, narrow_band_width, &
                                            sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
            else
                call compute_scalar_distance_face(myid, cx_lo, cx_hi, &
                                                  sy, ey, sz, ez, xp, yf, zp, nfaces, faces, face_normals, &
                                                  vertices, normals, face_solid_id, buffer_points, &
                                                  sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
            end if
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                          1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        if (myid == 0) &
            call raycast_classify_domain(vertices, faces, face_solid_id, nfaces, nsolids, &
                                         xp, yf, zp, sx, ex, sy, ey, sz, ez, scalarvalue, &
                                         flood_fill_arr, objid_fill_arr)
        if (use_fast_sweep) then
            call MPI_BCAST(flood_fill_arr, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
            allocate(fsm_local(decomp_x_start(myid):decomp_x_end(myid), ny, nz))
            fsm_local = flood_fill_arr(decomp_x_start(myid):decomp_x_end(myid), :, :)
            call fast_sweep_3d(fsm_local, decomp_x_start(myid), decomp_x_end(myid), 1, ny, 1, nz, &
                               xp, yf, zp, myid, left_rank, right_rank)
            call gather_array(myid, nprocs, fsm_local, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                              1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
            deallocate(fsm_local)
        end if
        if (myid == 0) then
            call cpu_time(time1)
            print *, "*** Writing | v-faces ***"
            call write_sdf_padded('data/sdfv.bin', flood_fill_arr)
            call cpu_time(time2)
            print *, "-- Write done in", time2-time1, "s | v-faces"
            print *, "*** Calculating SDF | w-faces ***"
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)

        ! W-faces
        sdf = scalarvalue
        if (cx_lo <= cx_hi) then
            if (use_fast_sweep) then
                call compute_narrowband_sdf(myid, cx_lo, cx_hi, &
                                            sy, ey, sz, ez, xp, yp, zf, nfaces, faces, face_normals, &
                                            vertices, normals, face_solid_id, narrow_band_width, &
                                            sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
            else
                call compute_scalar_distance_face(myid, cx_lo, cx_hi, &
                                                  sy, ey, sz, ez, xp, yp, zf, nfaces, faces, face_normals, &
                                                  vertices, normals, face_solid_id, buffer_points, &
                                                  sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
            end if
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                          1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        if (myid == 0) &
            call raycast_classify_domain(vertices, faces, face_solid_id, nfaces, nsolids, &
                                         xp, yp, zf, sx, ex, sy, ey, sz, ez, scalarvalue, &
                                         flood_fill_arr, objid_fill_arr)
        if (use_fast_sweep) then
            call MPI_BCAST(flood_fill_arr, nx*ny*nz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
            allocate(fsm_local(decomp_x_start(myid):decomp_x_end(myid), ny, nz))
            fsm_local = flood_fill_arr(decomp_x_start(myid):decomp_x_end(myid), :, :)
            call fast_sweep_3d(fsm_local, decomp_x_start(myid), decomp_x_end(myid), 1, ny, 1, nz, &
                               xp, yp, zf, myid, left_rank, right_rank)
            call gather_array(myid, nprocs, fsm_local, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                              1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
            deallocate(fsm_local)
        end if
        if (myid == 0) then
            call cpu_time(time1)
            print *, "*** Writing | w-faces ***"
            call write_sdf_padded('data/sdfw.bin', flood_fill_arr)
            call cpu_time(time2)
            print *, "-- Write done in", time2-time1, "s | w-faces"
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    end if  ! compute_face_sdf

    if (myid == 0) then
        call cpu_time(endTime)
        totalTime = endTime - startTime
        print *, "*** SDF complete in", totalTime, "s ***"
    end if

    call free_hash()
    call MPI_FINALIZE(ierror)
end program generatesdf
