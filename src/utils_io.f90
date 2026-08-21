module utils_io
    !
    ! utils_io — data, grid I/O, geometry I/O, and output routines for GenSDF v2.
    !
    implicit none

    !  precision 
    integer, parameter, public :: dp = selected_real_kind(15,307), &
                                  sp = selected_real_kind(6 , 37)

    !  geometry 
    character(len=512) :: inputfilename

    !  Cartesian grid 
    integer, protected :: nx, ny, nz
    real(dp), dimension(3), protected :: r0
    integer, dimension(3) :: ng
    logical, protected :: non_uniform_grid
    real(dp), allocatable, dimension(:) :: xp, yp, zp, xf, yf, zf
    real(dp) :: lx, ly, lz, dx, dy, dx_inverse, dy_inverse
    real(dp), allocatable, dimension(:) :: dz, dz_inverse

    !  SDF / algorithm control 
    real(dp), protected :: scalarvalue
    integer,  protected :: buffer_points
    logical,  protected :: use_fast_sweep
    integer,  protected :: narrow_band_width
    !  2 = solver wall-normal is y;  3 = solver wall-normal is z (default, z is wall-normal)
    integer,  protected :: vertical_axis = 3
    ! Compute face SDFs (sdfu, sdfv, sdfw) in addition to cell-centre (sdfp).
    ! Default .false. — only sdfp is written.
    logical,  protected :: compute_face_sdf = .false.

    !  MPI 
    integer :: nprocs

    !  Auxiliary
    integer, protected :: pbarwidth

    !$acc declare create(dx, dy, scalarvalue)
    ! dz is referenced inside godunov_update_cell (!$acc routine seq), which requires
    ! module allocatables it touches to be declared, not just enter-data'd.
    !$acc declare create(dz)

contains

    ! Read input file 
    subroutine read_inputfile(procid)
        implicit none
        integer, intent(in) :: procid
        integer :: iunit, ierr
        character(len=512) :: dummyline

        iunit = 10
        open(newunit=iunit, file='parameters.in', status='old', action='read', iostat=ierr)
        if (ierr == 0) then
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) inputfilename
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) scalarvalue, buffer_points, pbarwidth
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) nx, ny, nz
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) r0(1), r0(2), r0(3)
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) non_uniform_grid
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) use_fast_sweep
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) narrow_band_width
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) vertical_axis
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) compute_face_sdf
            if (vertical_axis /= 2 .and. vertical_axis /= 3) then
                if (procid == 0) &
                    print *, "ERROR: vertical_axis must be 2 or 3, got:", vertical_axis
                error stop
            end if
            ! Permute so GenSDF's internal dir 3 (z) always = wall-normal direction
            if (vertical_axis == 2) then
                ng = [nx, nz, ny]
                nx = ng(1); ny = ng(2); nz = ng(3)
            end if
            ng(1) = nx
            ng(2) = ny
            ng(3) = nz
            if (procid == 0) print *, "*** Input file successfully read ***"
        else
            if (procid == 0) error stop "ERROR: parameters.in encountered a problem!"
        end if
        close(iunit)

    contains
        ! Skip blank lines and lines whose first non-blank character is '#', then
        ! return the first data line in `line`.
        subroutine skip_and_read(u, line, ios)
            integer,          intent(in)  :: u
            character(len=*), intent(out) :: line
            integer,          intent(out) :: ios
            character(len=1) :: fc
            integer :: pos
            do
                read(u,'(A)',iostat=ios) line
                if (ios /= 0) return
                pos = verify(line, ' ')           ! index of first non-space
                if (pos == 0) cycle               ! blank line → skip
                fc = line(pos:pos)
                if (fc /= '#') return             ! found a data line
                ! else: comment line → skip
            end do
        end subroutine skip_and_read
    end subroutine read_inputfile

    ! CaNS grid still default as it is easy and common 
    ! Read CaNS-style grid files: geometry.out (grid dimensions and lengths) and optionally grid.out (non-uniform z).
    subroutine read_cans_grid(procid, loc, iprecision, npoints, origin, non_uni_grid, &
                               xin_p, yin_p, zin_p, xin_f, yin_f, zin_f, dzin, &
                               inputlx, inputly, inputlz)
        implicit none
        integer,  intent(in)  :: procid, iprecision
        character(len=*), intent(in) :: loc
        real(dp), dimension(3), intent(in) :: origin
        logical,  intent(in)  :: non_uni_grid
        integer,  dimension(3), intent(out) :: npoints
        real(dp), intent(out), allocatable, dimension(:) :: xin_p, yin_p, zin_p
        real(dp), intent(out), allocatable, dimension(:) :: xin_f, yin_f, zin_f
        real(dp), intent(out), allocatable, dimension(:) :: dzin
        real(dp), intent(out) :: inputlx, inputly, inputlz
        logical :: fexists
        integer :: iter, ii, jj, unit
        real(dp) :: dl(3), l(3)
        real(sp), allocatable :: grid_z4(:,:)
        real(dp), allocatable :: grid_z8(:,:)
        character(len=512) :: geofile, grdfile, linebuf
        integer :: filerr

        inquire(file=trim(loc), exist=fexists)
        if (.not. fexists) then
            print *, "Input directory does not exist: ", trim(loc)
            stop
        end if

        geofile = trim(loc) // "geometry.out"
        open(newunit=unit, file=geofile, status='old')
        ! Skip comment lines starting with '#'
        do
            read(unit,'(A)',iostat=filerr) linebuf
            if (filerr /= 0) exit
            if (len_trim(adjustl(linebuf)) == 0) cycle
            if (linebuf(verify(linebuf,' '):verify(linebuf,' ')) /= '#') then
                read(linebuf, *) npoints; exit
            end if
        end do
        do
            read(unit,'(A)',iostat=filerr) linebuf
            if (filerr /= 0) exit
            if (len_trim(adjustl(linebuf)) == 0) cycle
            if (linebuf(verify(linebuf,' '):verify(linebuf,' ')) /= '#') then
                read(linebuf, *) l; exit
            end if
        end do
        close(unit)
        ! Permute so GenSDF's internal dir 3 (z) = wall-normal direction
        if (vertical_axis == 2) then
            npoints = [npoints(1), npoints(3), npoints(2)]
            l       = [l(1),       l(3),       l(2)      ]
        end if
        inputlx = l(1); inputly = l(2); inputlz = l(3)
        dl = l / real(npoints, dp)

        allocate(xin_p(npoints(1)), yin_p(npoints(2)), zin_p(npoints(3)))
        allocate(xin_f(npoints(1)), yin_f(npoints(2)), zin_f(npoints(3)))
        allocate(dzin(npoints(3)))

        do iter = 1, npoints(1)
            xin_p(iter) = origin(1) + (iter - 0.5_dp) * dl(1)
            xin_f(iter) = xin_p(iter) + dl(1) / 2.0_dp
        end do
        do iter = 1, npoints(2)
            yin_p(iter) = origin(2) + (iter - 0.5_dp) * dl(2)
            yin_f(iter) = yin_p(iter) + dl(2) / 2.0_dp
        end do
        do iter = 1, npoints(3)
            zin_p(iter) = origin(3) + (iter - 0.5_dp) * dl(3)
            zin_f(iter) = zin_p(iter) + dl(3) / 2.0_dp
        end do

        if (non_uni_grid) then
            grdfile = trim(loc) // "grid.out"
            if (iprecision == 4) then
                open(newunit=unit, file=grdfile, status='old', action='read', iostat=filerr)
                if (filerr /= 0) error stop "Error opening grid.out"
                allocate(grid_z4(npoints(3),5))
                do ii = 1, npoints(3)
                    read(unit,*,iostat=filerr) (grid_z4(ii,jj), jj=1,5)
                    if (filerr /= 0) error stop "Error reading grid.out"
                end do
                close(unit)
                zin_p = origin(3) + grid_z4(:,2)
                zin_f = origin(3) + grid_z4(:,3)
            else
                open(newunit=unit, file=grdfile, status='old', action='read', iostat=filerr)
                if (filerr /= 0) error stop "Error opening grid.out"
                allocate(grid_z8(npoints(3),5))
                do ii = 1, npoints(3)
                    read(unit,*,iostat=filerr) (grid_z8(ii,jj), jj=1,5)
                    if (filerr /= 0) error stop "Error reading grid.out"
                end do
                close(unit)
                zin_p = origin(3) + grid_z8(:,2)
                zin_f = origin(3) + grid_z8(:,3)
            end if
        end if

        if (procid == 0) print *, "*** Successfully read the grid ***"
    end subroutine read_cans_grid

    ! --
    subroutine setup_grid_spacing(procid, xin, yin, zin, nz_in, &
                                   dx_out, dy_out, dz_out, dx_inv_out, dy_inv_out, dz_inv_out)
        implicit none
        integer,  intent(in)  :: procid, nz_in
        real(dp), intent(in),  dimension(:) :: xin, yin, zin
        real(dp), intent(out) :: dx_out, dy_out, dx_inv_out, dy_inv_out
        real(dp), intent(out), dimension(:) :: dz_out, dz_inv_out
        integer :: ii

        dx_out = xin(2) - xin(1)
        dy_out = yin(2) - yin(1)
        dx_inv_out = 1.0_dp / dx_out
        dy_inv_out = 1.0_dp / dy_out
        ! z spacing: face(ii+1) - face(ii) gives exact width of cell ii
        ! Loop only to nz_in-1 to avoid reading zin(nz_in+1) past the array end;
        ! the last cell reuses the second-to-last spacing.
        do ii = 1, nz_in - 1
            dz_out(ii)     = zin(ii+1) - zin(ii)
            dz_inv_out(ii) = 1.0_dp    / dz_out(ii)
        end do
        dz_out(nz_in)     = dz_out(nz_in - 1)
        dz_inv_out(nz_in) = dz_inv_out(nz_in - 1)

        if (procid == 0) print *, "*** Successfully set up grid spacing ***"
    end subroutine setup_grid_spacing

    ! Read OBJ file (ASCII, triangulated, with vertex normals)
    subroutine read_obj(procid, filename, vertices, normals, faces, face_normals, &
                        num_vertices, num_normals, num_faces)
        implicit none
        integer,  intent(in)  :: procid
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: vertices(:,:), normals(:,:)
        integer,  allocatable, intent(out) :: faces(:,:), face_normals(:,:)
        integer,  intent(out) :: num_vertices, num_normals, num_faces
        integer :: unit, ios, v_count, vn_count, f_count
        character(len=512) :: line, dline
        real(dp) :: verx, very, verz, normx, normy, normz
        integer  :: verind1, verind2, verind3, vn1, vn2, vn3
        real(dp), allocatable :: tmp_row(:)
        integer,  allocatable :: tmp_irow(:)

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', trim(filename); stop
        end if

        v_count = 0; vn_count = 0; f_count = 0
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            select case (trim(adjustl(line(1:2))))
                case ('v '); v_count  = v_count  + 1
                case ('vn'); vn_count = vn_count + 1
                case ('f '); f_count  = f_count  + 1
            end select
        end do

        if (vn_count == 0 .and. procid == 0) then
            print *, "FATAL: Geometry has no normals. Use mesh.fix_normals() in trimesh."
            error stop
        end if

        allocate(vertices(3,v_count), normals(3,vn_count))
        allocate(faces(3,f_count), face_normals(3,f_count))
        num_vertices = v_count; num_normals = vn_count; num_faces = f_count

        rewind(unit)
        v_count = 0; vn_count = 0; f_count = 0
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (trim(adjustl(line(1:2))) == 'v ') then
                read(line(3:),*) verx, very, verz
                v_count = v_count + 1
                vertices(:, v_count) = [verx, very, verz]
            else if (trim(adjustl(line(1:3))) == 'vn ') then
                read(line(4:),*) normx, normy, normz
                vn_count = vn_count + 1
                normals(:, vn_count) = [normx, normy, normz]
            else if (trim(adjustl(line(1:2))) == 'f ') then
                dline = trim(line(3:))
                call parse_face_line(dline, verind1, vn1, verind2, vn2, verind3, vn3)
                f_count = f_count + 1
                faces(:, f_count)        = [verind1, verind2, verind3]
                face_normals(:, f_count) = [vn1, vn2, vn3]
            end if
        end do
        close(unit)

        if (vertical_axis == 2) then
            allocate(tmp_row(num_vertices))
            tmp_row = vertices(2,:); vertices(2,:) = vertices(3,:); vertices(3,:) = tmp_row
            deallocate(tmp_row)
            allocate(tmp_row(num_normals))
            tmp_row = normals(2,:); normals(2,:) = normals(3,:); normals(3,:) = tmp_row
            deallocate(tmp_row)
            ! y<->z is a reflection (det=-1): reverse winding so cross-product normals stay outward
            allocate(tmp_irow(num_faces))
            tmp_irow = faces(2,:); faces(2,:) = faces(3,:); faces(3,:) = tmp_irow
            tmp_irow = face_normals(2,:); face_normals(2,:) = face_normals(3,:); face_normals(3,:) = tmp_irow
            deallocate(tmp_irow)
        end if

        if (procid == 0) then
            print *, 'Successfully read OBJ: ', trim(filename)
            print *, 'Vertices: ', num_vertices, '  Normals: ', num_normals, '  Faces: ', num_faces
        end if
    end subroutine read_obj

    ! Read STL file (binary)
    subroutine read_stl_binary(procid, filename, vertices, normals, faces, face_normals, &
                               num_vertices, num_normals, num_faces)
        !
        ! Binary STL reader (§3).
        ! File format (little-endian, all ranks read independently):
        !   80 bytes  — header (ignored)
        !    4 bytes  — uint32 triangle count  (nf)
        !   Per triangle (50 bytes):
        !     3 × float32  — face normal
        !     3 × float32  — vertex 0
        !     3 × float32  — vertex 1
        !     3 × float32  — vertex 2
        !     2 bytes      — attribute (ignored)
        !
        ! Output arrays follow the same layout as read_obj:
        !   vertices(3, 3*nf)  — all vertex coordinates (3 per triangle)
        !   normals (3, nf)    — per-face normals (no averaging needed for STL)
        !   faces   (3, nf)    — vertex index triples (1-based, packed sequentially)
        !   face_normals(3,nf) — normal index triples (all equal to the face index)
        !
        use iso_fortran_env, only : int32
        implicit none
        integer,  intent(in)  :: procid
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: vertices(:,:), normals(:,:)
        integer,  allocatable, intent(out) :: faces(:,:), face_normals(:,:)
        integer,  intent(out) :: num_vertices, num_normals, num_faces
        ! Local
        integer            :: unit, ios, fi
        integer(int32)     :: nf_i32
        character(len=80)  :: header
        integer(2)         :: attr
        real(4)            :: fn(3), v0(3), v1(3), v2(3)  ! float32 reads
        real(dp), allocatable :: tmp_row(:)
        integer,  allocatable :: tmp_irow(:)

        open(newunit=unit, file=trim(filename), access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening STL file: ", trim(filename); stop
        end if

        read(unit) header   ! 80-byte header (ignored)
        read(unit) nf_i32   ! triangle count
        num_faces    = int(nf_i32)
        num_vertices = 3 * num_faces
        num_normals  = num_faces

        allocate(vertices(3, num_vertices))
        allocate(normals (3, num_normals))
        allocate(faces      (3, num_faces))
        allocate(face_normals(3, num_faces))

        do fi = 1, num_faces
            read(unit) fn(1), fn(2), fn(3)    ! face normal (float32)
            read(unit) v0(1), v0(2), v0(3)    ! vertex 0
            read(unit) v1(1), v1(2), v1(3)    ! vertex 1
            read(unit) v2(1), v2(2), v2(3)    ! vertex 2
            read(unit) attr                   ! attribute bytes (ignored)

            ! Store face normal (convert float32 -> dp)
            normals(1,fi) = real(fn(1),dp); normals(2,fi) = real(fn(2),dp); normals(3,fi) = real(fn(3),dp)

            ! Store vertices packed 3-per-face (1-based)
            vertices(:, 3*fi-2) = [real(v0(1),dp), real(v0(2),dp), real(v0(3),dp)]
            vertices(:, 3*fi-1) = [real(v1(1),dp), real(v1(2),dp), real(v1(3),dp)]
            vertices(:, 3*fi  ) = [real(v2(1),dp), real(v2(2),dp), real(v2(3),dp)]

            ! Faces point into the packed vertex array
            faces(:,fi)        = [3*fi-2, 3*fi-1, 3*fi]
            ! All three face-normal indices point to the same per-face normal
            face_normals(:,fi) = [fi, fi, fi]
        end do
        close(unit)

        if (vertical_axis == 2) then
            allocate(tmp_row(num_vertices))
            tmp_row = vertices(2,:); vertices(2,:) = vertices(3,:); vertices(3,:) = tmp_row
            deallocate(tmp_row)
            allocate(tmp_row(num_normals))
            tmp_row = normals(2,:); normals(2,:) = normals(3,:); normals(3,:) = tmp_row
            deallocate(tmp_row)
            ! y<->z is a reflection (det=-1): reverse winding so cross-product normals stay outward
            allocate(tmp_irow(num_faces))
            tmp_irow = faces(2,:); faces(2,:) = faces(3,:); faces(3,:) = tmp_irow
            deallocate(tmp_irow)
        end if

        if (procid == 0) then
            print *, "Successfully read binary STL: ", trim(filename)
            print *, "Faces: ", num_faces
        end if
    end subroutine read_stl_binary

    ! Parse face line (OBJ helper)
    subroutine parse_face_line(line, v1, vn1, v2, vn2, v3, vn3)
        implicit none
        character(len=512), intent(inout) :: line
        integer, intent(out) :: v1, vn1, v2, vn2, v3, vn3
        character(len=32) :: part
        integer :: pos

        pos = index(line,' '); part = trim(line(1:pos-1))
        call parse_vertex_normal_pair(part, v1, vn1); line = trim(line(pos+1:))
        pos = index(line,' '); part = trim(line(1:pos-1))
        call parse_vertex_normal_pair(part, v2, vn2); line = trim(line(pos+1:))
        part = trim(line)
        call parse_vertex_normal_pair(part, v3, vn3)
    end subroutine parse_face_line

    subroutine parse_vertex_normal_pair(pair, vertex, normal)
        !
        ! Parse a single OBJ face token.  Supported formats:
        !   v//vn      — vertex and normal, no texture (original format)
        !   v/vt/vn    — vertex, texture, and normal
        ! Vertex-only tokens (no slash) are rejected because normals are required.
        !
        implicit none
        character(len=*), intent(inout) :: pair
        integer, intent(out) :: vertex, normal
        integer :: dpos, dpos2
        dpos = index(pair,'//')
        if (dpos > 0) then
            ! Format v//vn
            read(pair(1:dpos-1),*) vertex
            read(pair(dpos+2:),  *) normal
        else
            ! Try format v/vt/vn: locate both slashes
            dpos  = index(pair,'/')
            dpos2 = index(pair(dpos+1:),'/') + dpos   ! position of second '/'
            if (dpos > 0 .and. dpos2 > dpos) then
                read(pair(1:dpos-1), *) vertex
                read(pair(dpos2+1:), *) normal
            else
                print *, "FATAL: unrecognised OBJ face token (expected v//vn or v/vt/vn): ", trim(pair)
                error stop
            end if
        end if
    end subroutine parse_vertex_normal_pair

    ! Get axis-aligned bounding box of the geometry
    subroutine getbbox(procid, vertices_in, num_vertices_in, bbox_min, bbox_max)
        implicit none
        integer,  intent(in)  :: procid, num_vertices_in
        real(dp), intent(in)  :: vertices_in(3, num_vertices_in)
        real(dp), intent(out) :: bbox_min(3), bbox_max(3)
        integer :: iter
        bbox_min = vertices_in(:,1); bbox_max = vertices_in(:,1)
        do iter = 2, num_vertices_in
            where (vertices_in(:,iter) < bbox_min) bbox_min = vertices_in(:,iter)
            where (vertices_in(:,iter) > bbox_max) bbox_max = vertices_in(:,iter)
        end do
        if (procid == 0) then
            print *, "Geometry bbox min:", bbox_min
            print *, "Geometry bbox max:", bbox_max
        end if
    end subroutine getbbox

    ! Tag grid points inside bounding box (with buffer) for focused SDF computation
    subroutine tagminmax(procid, xin, yin, zin, bbox_min, bbox_max, &
                          nx_in, ny_in, nz_in, dx_in, dy_in, dz_in, sdfresolution, &
                          sx, ex, sy, ey, sz, ez)
        implicit none
        integer,  intent(in)  :: procid, nx_in, ny_in, nz_in, sdfresolution
        real(dp), dimension(:), intent(in) :: xin, yin, zin
        real(dp), intent(in)  :: dx_in, dy_in, dz_in
        real(dp), dimension(:), intent(in) :: bbox_min, bbox_max
        integer,  intent(out) :: sx, ex, sy, ey, sz, ez
        integer  :: iteri, iterj, iterk
        real(dp) :: buf

        sx = 1; ex = nx_in; sy = 1; ey = ny_in; sz = 1; ez = nz_in
        buf = real(sdfresolution, dp)

        do iteri = 1, nx_in
            if (xin(iteri) <= bbox_min(1) - buf*dx_in) then; sx = iteri
            else if (xin(iteri) <= bbox_max(1) + buf*dx_in) then; ex = iteri
            else; exit; end if
        end do
        do iterj = 1, ny_in
            if (yin(iterj) <= bbox_min(2) - buf*dy_in) then; sy = iterj
            else if (yin(iterj) <= bbox_max(2) + buf*dy_in) then; ey = iterj
            else; exit; end if
        end do
        do iterk = 1, nz_in
            if (zin(iterk) <= bbox_min(3) - buf*dz_in) then; sz = iterk
            else if (zin(iterk) <= bbox_max(3) + buf*dz_in) then; ez = iterk
            else; exit; end if
        end do

        if (procid == 0) then
            print *, "AABB x:", sx, xin(sx), "|", ex, xin(ex)
            print *, "AABB y:", sy, yin(sy), "|", ey, yin(ey)
            print *, "AABB z:", sz, zin(sz), "|", ez, zin(ez)
        end if
    end subroutine tagminmax

    ! Gather local SDF arrays from all ranks into a global array on the root rank
    subroutine gather_array(proc_id, num_procs, local_sdf, global_nx, global_ny, global_nz, &
                             start_indices, end_indices, chunk_sizes, s_y, e_y, s_z, e_z, &
                             input_scalarvalue, assembled_array, ierr_mpi)
        use mpi
        implicit none
        integer,  intent(in)  :: proc_id, num_procs, global_nx, global_ny, global_nz
        integer,  intent(in)  :: s_y, e_y, s_z, e_z
        integer,  intent(in),  dimension(0:) :: start_indices, end_indices, chunk_sizes
        real(dp), intent(in),  dimension(:,:,:) :: local_sdf
        real(dp), intent(in)  :: input_scalarvalue
        real(dp), allocatable, intent(inout), dimension(:,:,:) :: assembled_array
        integer,  intent(out) :: ierr_mpi
        real(dp), allocatable, dimension(:) :: sendbuf, recvbuf
        integer,  allocatable, dimension(:)  :: recvcounts, displs
        integer :: rank, local_count

        ! Build per-rank send counts and receive displacements
        allocate(recvcounts(0:num_procs-1), displs(0:num_procs-1))
        do rank = 0, num_procs-1
            recvcounts(rank) = chunk_sizes(rank) * global_ny * global_nz
        end do
        displs(0) = 0
        do rank = 1, num_procs-1
            displs(rank) = displs(rank-1) + recvcounts(rank-1)
        end do

        ! Pack local slab into a contiguous send buffer
        local_count = chunk_sizes(proc_id) * global_ny * global_nz
        allocate(sendbuf(local_count))
        sendbuf = reshape(local_sdf, [local_count])

        ! Root receives all slabs; allocate assembled_array if not already done
        if (proc_id == 0) then
            if (.not. allocated(assembled_array)) then
                allocate(assembled_array(global_nx, global_ny, global_nz))
            end if
            assembled_array = scalarvalue
            allocate(recvbuf(sum(recvcounts)))
        else
            allocate(recvbuf(1))   ! unused on non-root ranks
        end if

        call MPI_GATHERV(sendbuf, local_count, MPI_DOUBLE_PRECISION, &
                         recvbuf, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                         0, MPI_COMM_WORLD, ierr_mpi)

        if (proc_id == 0) then
            do rank = 0, num_procs-1
                assembled_array(start_indices(rank):end_indices(rank), :, :) = &
                    reshape(recvbuf(displs(rank)+1 : displs(rank)+recvcounts(rank)), &
                            [chunk_sizes(rank), global_ny, global_nz])
            end do
            deallocate(recvbuf)
        end if

        deallocate(sendbuf, recvcounts, displs)
    end subroutine gather_array

    ! Write SDF to file, with ghost padding in x and wall-normal directions as expected by the solver's IBM reader.
    subroutine write_sdf_padded(outputfilename, interior)
        !
        ! Write SDF ghost-padded in x and the wall-normal direction, matching
        ! the solver IBM reader expectation of (nxm+2) x (n_wall+2) x n_span
        ! pages in solver coordinates.
        !
        implicit none
        character(len=*), intent(in) :: outputfilename
        real(dp), intent(in), dimension(:,:,:) :: interior
        integer :: d1, d2, d3, unit, j, k
        real(dp), allocatable :: padded(:,:,:)

        d1 = size(interior,1); d2 = size(interior,2); d3 = size(interior,3)

        ! Solver based padding on interior SDF
        allocate(padded(d1+2, d3+2, d2))
        padded = 0.0_dp
        do k = 1, d2           ! k loops spanwise = GenSDF dim 2
            do j = 1, d3       ! j loops wall-normal = GenSDF dim 3
                padded(2:d1+1, j+1, k) = interior(:, k, j)
            end do
        end do
        ! x ghost (dim 1)
        padded(1,    2:d3+1, :) = padded(2,    2:d3+1, :)
        padded(d1+2, 2:d3+1, :) = padded(d1+1, 2:d3+1, :)
        ! wall-normal ghost (dim 2)
        padded(:, 1,      :) = padded(:, 2,     :)
        padded(:, d3+2,   :) = padded(:, d3+1,  :)

        open(newunit=unit, file=trim(outputfilename), access='stream', status='replace', form='unformatted', &
             convert='big_endian')
        write(unit) padded
        close(unit)
        deallocate(padded)
    end subroutine write_sdf_padded

    ! Estimate memory usage based on input parameters
    subroutine estimated_memoryusage(n_procs, input_nfaces, input_nvertices, in_nx, in_ny, in_nz)
        implicit none
        integer, intent(in) :: n_procs, input_nfaces, input_nvertices, in_nx, in_ny, in_nz
        real(dp) :: total
        character(len=256) :: str
        total  = 6.0_dp*real(input_nfaces,    dp)*real(sizeof(in_nx),   dp)
        total  = total + 6.0_dp*real(input_nvertices, dp)*real(sizeof(total), dp)
        total  = total + 2.0_dp*real(in_nx*in_ny*in_nz, dp)*real(sizeof(total), dp)
        write(str,'(F5.2)') (n_procs*total)/1e9_dp
        print *, "-- Estimated minimum memory:", trim(str), " GiB(s)"
    end subroutine estimated_memoryusage

    ! Progress bar
    subroutine show_progress(current, total, width, start_time)
        use mpi, only : MPI_WTIME
        implicit none
        integer,  intent(in) :: current, total, width
        real(dp), intent(in) :: start_time
        integer  :: progress
        real(dp) :: percent, elapsed, remaining, now
        character(len=:), allocatable :: bar

        now = MPI_WTIME()
        elapsed = now - start_time
        if (current > 0) then
            remaining = elapsed * (real(total,dp) - real(current,dp)) / real(current,dp)
        else
            remaining = 0.0_dp
        end if
        if (total > 0) then
            percent  = real(current,dp) / real(total,dp) * 100.0_dp
            progress = int(real(current,dp) / real(total,dp) * width)
        else
            percent  = 0.0_dp
            progress = 0
        end if
        allocate(character(len=width) :: bar)
        bar = repeat('|', progress) // repeat(' ', width-progress)
        write(*,'(a,"|",a,"|",f6.2,"%",1x,i8,"/",i8," Elapsed:",f8.2,"s Remaining:",f8.2,"s")', &
              advance='no') char(13), bar, percent, current, total, elapsed, remaining
        flush(6)
        if (current == total) print *
        deallocate(bar)
    end subroutine show_progress

    ! LOGOOOO 
    subroutine printlogo()
        implicit none
        print *, " ░▒▓██████▓▒░░▒▓████████▓▒░▒▓███████▓▒░ ░▒▓███████▓▒░▒▓███████▓▒░░▒▓████████▓▒░ "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒▒▓███▓▒░▒▓██████▓▒░ ░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓██████▓▒░   "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, " ░▒▓██████▓▒░░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓███████▓▒░░▒▓█▓▒░        "
    end subroutine printlogo

end module utils_io
