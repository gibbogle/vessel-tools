! To estimate the probability distribution of distance from tissue to the nearest blood vessel.
! Uses points within the convex hull of a random set of points (?) 
!
! New method: use the close file from peel2 to determine if a point is inside or outside

module data_mod
use par_zig_mod
!use delaunay
use omp_lib

implicit none

type segment_type
	real :: end1(3), mid(3), end2(3)
	real :: r1, r2
end type

type segdist_type
	integer :: iseg
	real :: d2
end type

type edge_type
    integer :: npts
    integer, allocatable :: pt(:)
end type

real, allocatable :: point(:,:)

logical(2), allocatable :: in(:,:,:)
real, allocatable :: random_pt(:,:)

type(segment_type), allocatable :: segment(:)	! this is the list of all segments
type(edge_type), allocatable :: edge(:)
integer, allocatable :: seglist(:,:,:,:)		! for each (ib,jb,kb), a list of segment numbers
integer, allocatable :: nseglist(:,:,:)			! for each (ib,jb,kb), number of segments in the seglist
byte, allocatable :: imagedata(:,:,:)
byte, allocatable :: closedata(:)

real :: grid_dx, delp, pt_factor
real :: rmin(3), rmax(3), rmid(3), del(3), delb(3), voxelsize(3)
real :: sphere_centre(3), sphere_radius, constant_r, threshold_d
integer :: Nsegpts, N(3), NB(3), Nsegments, Mnodes, np_random, np_grid, offset(3)
integer :: nx, ny, nz, nx8, ny8, nz8, nmbytes, nxm, nym, nzm
character*(128) :: amfile, distfile, datafile, closefile
logical :: use_sphere, use_random, use_constant_radius, save_imagedata, use_close

logical :: dbug

integer, parameter :: nfam = 10, nfdist=11, nfout=12, nfdata=13, nfclose=14
real, parameter :: blocksize = 120
integer, parameter :: seed(2) = (/12345, 67891/)
integer, parameter :: MAX_NX = 900  ! a guess
integer, parameter :: test = 0
logical, parameter :: jiggle = .true.

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine input(res)
integer :: res
integer :: nlen, cnt, i, k, status
character*(2048) :: c, progname

call get_command (c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed to get command line with status = ', status
    write (nfout,*) 'get_command failed to get command line with status = ', status
    res = 1
    return
end if
write (*,*) 'command line = ', c(1:nlen)
write (nfout,*) 'command line = ', c(1:nlen)
call get_command_argument (0, c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'Getting command progname failed with status = ', status
    write (nfout,*) 'Getting command progname failed with status = ', status
    res = 1
    return
end if
progname = c(1:nlen)
cnt = command_argument_count ()
write (*,*) 'number of command arguments = ', cnt
write (nfout,*) 'number of command arguments = ', cnt
res = 0
if (cnt == 10) then
	use_sphere = .false.
	use_close = .false.
elseif (cnt == 11) then
	use_sphere = .false.
	use_close = .true.
elseif (cnt == 14) then
	use_sphere = .true.
	use_close = .false.
elseif (cnt == 15) then
	use_sphere = .true.
	use_close = .true.
else
    write(*,*) 'Bad command line argument count'
    write(nfout,*) 'Bad command line argument count'
	res = 2
    write(*,*) 'Use either: ',trim(progname), &
    ' amfile distfile grid_dx ncpu pt_factor threshold datafile vx vy vz' ! 10
    write(*,*) ' to analyze the whole network'
    write(*,*) 'or: ',trim(progname), &
    ' amfile distfile grid_dx ncpu pt_factor threshold datafile vx vy vz x0 y0 z0 R'  ! 14
    write(*,*) ' to analyze a spherical subregion with centre (x0,y0,z0), radius R'
    write(*,*) 'or: ',trim(progname), &
    ' amfile distfile grid_dx ncpu pt_factor threshold datafile vx vy vz close_file' ! 11
    write(*,*) ' to analyze the whole network using the close file to determine insideness'
    write(*,*) 'or: ',trim(progname), &
    ' amfile distfile grid_dx ncpu pt_factor threshold datafile vx vy vz x0 y0 z0 R close_file'  ! 15
    write(*,*) ' datafile is the temporary file used to pass the distribution data to maketiff which does the conversion'
    write(*,*) ' to analyze a spherical subregion with centre (x0,y0,z0), radius R, using the close file'
    write(*,*) ' If grid_dx = 0 the sampling points and randomly placed'
!   write(*,*) ' If constant_radius != 0 all vessels are given this radius'
    write(*,*) ' If threshold != 0 the sampling points are used to generate the image_file data, voxel=255 if distance > threshold'
    write(*,*) ' If the close files is used, voxel dimensions (um) must be specified: vx, vy, vz'
    return
endif

do i = 1, 10
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        write (nfout,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        res = 3
        return
    end if
    if (i == 1) then
        amfile = c(1:nlen)
        write(*,*) 'Network file: ',amfile
        write(nfout,*) 'Network file: ',amfile
    elseif (i == 2) then
        distfile = c(1:nlen)
        write(*,*) 'Distribution file: ',distfile
        write(nfout,*) 'Distribution file: ',distfile
    elseif (i == 3) then
        read(c(1:nlen),*) grid_dx																
        write(*,*) 'grid_dx: ',grid_dx
        write(nfout,*) 'grid_dx: ',grid_dx
    elseif (i == 4) then
        read(c(1:nlen),*) Mnodes
        write(*,*) 'Mnodes: ',Mnodes															
        write(nfout,*) 'Mnodes: ',Mnodes															
    elseif (i == 5) then
        read(c(1:nlen),*) pt_factor
        write(*,*) 'pt_factor: ',pt_factor
        write(nfout,*) 'pt_factor: ',pt_factor
    elseif (i == 6) then
        read(c(1:nlen),*) threshold_d
        write(*,*) 'threshold_d: ',threshold_d
        write(nfout,*) 'threshold_d: ',threshold_d
    elseif (i == 7) then
        read(c(1:nlen),*) datafile		
        write(*,*) 'datafile: ',datafile													
        write(nfout,*) 'datafile: ',datafile													
    elseif (i == 8) then
        read(c(1:nlen),*) voxelsize(1)																
        write(*,*) 'voxelsize(1): ',voxelsize(1)
        write(nfout,*) 'voxelsize(1): ',voxelsize(1)
    elseif (i == 9) then
        read(c(1:nlen),*) voxelsize(2)																
        write(*,*) 'voxelsize(2): ',voxelsize(2)
        write(nfout,*) 'voxelsize(2): ',voxelsize(2)
    elseif (i == 10) then
        read(c(1:nlen),*) voxelsize(3)	
        write(*,*) 'voxelsize(3): ',voxelsize(3)
        write(nfout,*) 'voxelsize(3): ',voxelsize(3)
    endif
enddo
k = 10
if (use_sphere) then
    do i = k+1, k+4
        call get_command_argument (i, c, nlen, status)
        if (status .ne. 0) then
            write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
            write (nfout,*) 'get_command_argument failed: status = ', status, ' arg = ', i
            res = 3
            return
        end if														
        if (i == k+1) then
            read(c(1:nlen),*) sphere_centre(1)
            write(*,*) 'sphere_centre(1): ',sphere_centre(1)																
            write(nfout,*) 'sphere_centre(1): ',sphere_centre(1)																
        elseif (i == k+2) then
            read(c(1:nlen),*) sphere_centre(2)																
            write(*,*) 'sphere_centre(2): ',sphere_centre(2)																
            write(nfout,*) 'sphere_centre(2): ',sphere_centre(2)																
        elseif (i == k+3) then
            read(c(1:nlen),*) sphere_centre(3)																
            write(*,*) 'sphere_centre(3): ',sphere_centre(3)																
            write(nfout,*) 'sphere_centre(3): ',sphere_centre(3)																
        elseif (i == k+4) then
            read(c(1:nlen),*) sphere_radius
            write(*,*) 'sphere_radius: ',sphere_radius
            write(nfout,*) 'sphere_radius: ',sphere_radius
        endif
    enddo
    k = k+4
endif
if (use_close) then
    i = k+1
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        write (nfout,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        res = 3
        return
    end if
    closefile = c(1:nlen)
    write(*,*) 'closefile: ',closefile
    write(nfout,*) 'closefile: ',closefile
endif										
if (grid_dx == 0) then
    use_random = .true.
else
    use_random = .false.
endif
constant_r = 0
if (constant_r > 0) then
    use_constant_radius = .true.
else
    use_constant_radius = .false.
endif
if (threshold_d > 0) then
    save_imagedata = .true.
else
    save_imagedata = .false.
endif
if (use_sphere) then
	write(*,'(a,4f8.2)') 'Sphere centre, radius: ',sphere_centre,sphere_radius
	write(nfout,'(a,4f8.2)') 'Sphere centre, radius: ',sphere_centre,sphere_radius
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
end subroutine

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(*,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    stop
endif
R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))

end function


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine setup(res)
integer res
integer :: k, i
real :: rng

call rng_initialisation
if (use_close) then
    call read_closedata(res)
    if (res /= 0) return
endif

Nsegpts = 2*nsegments
allocate(point(Nsegpts,3))
k = 0
do i = 1,nsegments
	k = k+1
	point(k,:) = segment(i)%end1
	k = k+1
	point(k,:) = segment(i)%end2
enddo
if (use_sphere) then
	write(*,'(a,4f8.2)') 'Using sphere: centre, radius: ',sphere_centre,sphere_radius
	rmin = sphere_centre - sphere_radius
	rmax = sphere_centre + sphere_radius
elseif (use_close) then     ! get ranges from the close file
    call getcloseranges
else
    rmin = 1.0e10
    rmax = -1.0e10
    do k = 1,Nsegpts
	    rmin = min(rmin,point(k,:))
	    rmax = max(rmax,point(k,:))
    enddo
endif
rmid = (rmin + rmax)/2
write(*,'(a,3f8.2)') 'rmin: ',rmin
write(*,'(a,3f8.2)') 'rmax: ',rmax
if (save_imagedata) then
    allocate(imagedata(MAX_NX,MAX_NX,MAX_NX))
    imagedata = 0
    nx = 0
    ny = 0
    nz = 0
endif
res = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine create_in(res)
integer :: res
integer :: i, j, k, ip(3), indx(3), ix, iy, iz, ierr, Ngridpts, kpar=0
real :: R(3), rsum, xyz(3), v(3), tri(3,3)
logical :: localize = .true.

write(*,*) 'create_in'
res = 0
!if (use_sphere) then
!	write(*,'(a,4f8.2)') 'Using sphere: centre, radius: ',sphere_centre,sphere_radius
!	rmin = sphere_centre - sphere_radius
!	rmax = sphere_centre + sphere_radius
!endif

!rmin = rmin - 1
!rmax = rmax + 1
del = grid_dx
N = (rmax - rmin)/del + 1
rmax = rmin + (N-1)*del
offset = rmin/del + 0.001
write(*,*) 'del: ',del
write(*,*) 'offset: ',offset
write(*,*) 'Dimensions of array in(:,:,:): ',N
write(nfdist,*) 'Dimensions of array in(:,:,:): ',N
allocate(in(N(1),N(2),N(3)),stat=ierr)
if (ierr /= 0) then
	write(*,*) 'Allocation failed on array in(:,:,:).  Increase grid_dx'
	write(nfout,*) 'Allocation failed on array in(:,:,:).  Increase grid_dx'
	res = 4
	return
endif
in = .false.
np_grid = 0
!del = (rmax-rmin)/(N-1)
if (use_sphere) then
	do ix = 1,N(1)
		xyz(1) = rmin(1) + (ix-1)*del(1)
		do iy = 1,N(2)
			xyz(2) = rmin(2) + (iy-1)*del(2)
			do iz = 1,N(3)
				xyz(3) = rmin(3) + (iz-1)*del(3)
				v = xyz - sphere_centre
				if (dot_product(v,v) < sphere_radius*sphere_radius) then
				    if (use_close) then
				        in(ix,iy,iz) = in_close(xyz)
				    else
    					in(ix,iy,iz) = .true.
                    endif
                    if (in(ix,iy,iz)) np_grid = np_grid + 1
				endif
			enddo
		enddo
	enddo
else
    if (use_close) then
	    do ix = 1,N(1)
		    xyz(1) = rmin(1) + (ix-1)*del(1)
		    do iy = 1,N(2)
			    xyz(2) = rmin(2) + (iy-1)*del(2)
			    do iz = 1,N(3)
				    xyz(3) = rmin(3) + (iz-1)*del(3)
                    in(ix,iy,iz) = in_close(xyz)
                    if (in(ix,iy,iz)) np_grid = np_grid + 1
                enddo
            enddo
        enddo
    else
        Ngridpts = pt_factor*N(1)*N(2)*N(3)
        write(*,*) 'Generating random interior grid points: ',Ngridpts
        write(nfdist,*) 'Generating random interior grid points: ',Ngridpts
	    do i = 1, Ngridpts/10      ! 10000000
    !	do i = 1, 10000000
		    j = 1
		    ! generate three random points from the list of vertices
		    ! this creates a triangle that is (almost certainly) completely within the tissue region
		    ! Improve this by requiring that the second two points are in blocks near the first point
    		
		    if (localize) then
		        call get_triangle(tri)
		    else
		        do
			        ip(j) = random_int(1,Nsegpts,kpar)
			        if (j > 1) then 
				        if (ip(j) == ip(1)) cycle
				        if (j == 3 .and. ip(j) == ip(2)) cycle
			        endif
			        j = j+1
			        if (j == 4) exit
		        enddo
		        do j = 1,3
		            tri(j,:) = point(ip(j),:)
		        enddo
		    endif
    		
		    do k = 1,10
			    xyz = 0
			    rsum = 0
			    do j = 1,3
				    R(j) = par_uni(kpar)
				    rsum = rsum + R(j)
    !				xyz = xyz + R(j)*point(ip(j),:)
				    xyz = xyz + R(j)*tri(j,:)
			    enddo
			    xyz = xyz/rsum
			    ! xyz is a random point inside the triangle
			    do j = 1,3
				    indx(j) = (xyz(j)-rmin(j))/del(j) + 1
				    if (indx(j) > N(j)) then
					    write(nfout,*) 'indx out of range:'
					    write(*,*) 'indx out of range:'
					    write(*,*) 'point: ',tri(j,:)
					    write(*,*) 'R: ',R/rsum
					    write(*,*) 'xyz: ',xyz
					    res = 5
					    return
				    endif
			    enddo
			    ! indx(:) gives the location of the nearest point on the grid to xyz
			    if (.not.in(indx(1),indx(2),indx(3))) np_grid = np_grid + 1
			    in(indx(1),indx(2),indx(3)) = .true.
		    enddo
	    enddo
	endif
endif
write(*,*) 'np_grid: ',np_grid
write(nfdist,*) 'np_grid: ',np_grid
end subroutine

!-----------------------------------------------------------------------------------------
! Test if the point xyz(:) corresponds to a lit voxel in closedata(:,:,:)
! Now using compressed close data, with each bit corresponding to a voxel
!-----------------------------------------------------------------------------------------
logical function in_close(xyz)
real :: xyz(3)
integer :: p(3), kbyte, kbit, nb
byte :: mbyte

p = xyz/voxelsize + 1
if (p(1) > nxm .or. p(2) > nym .or. p(3) > nzm) then
    in_close = .false.
!    write(*,*) 'in_close: p out of range: ',p,xyz
    return
endif
! Need to convert to the byte and bit
nb = p(1) + (p(2)-1)*nx8 + (p(3)-1)*nx8*ny8
!if (mod(nb,8) == 0) then
!    kbyte = nb/8
!    kbit = 0
!else
!    kbyte = nb/8 + 1
!    kbit = nb  - 8*(kbyte-1)
!endif
kbyte = (nb-1)/8 + 1
kbit = nb - 8*(kbyte-1) - 1
if (kbyte > nmbytes) then
    write(*,'(a,6i12)') 'Error: kbyte > nmbytes: ',p,nb,kbyte,nmbytes
    write(*,'(a,3f8.1)') 'xyz: ',xyz
    call exit(8)
endif
mbyte = closedata(kbyte)
if (btest(mbyte,kbit)) then
    in_close = .true.
else
    in_close = .false.
endif
    
!if (closedata(p(1),p(2),p(3)) == 0) then
!    in_close = .false.
!else
!    in_close = .true.
!endif
end function

!-----------------------------------------------------------------------------------------
! Deduce ranges from the closedata
!-----------------------------------------------------------------------------------------
subroutine getcloseranges
integer :: ix, iy, iz, k, kbyte, kbit
real :: xyz(3)
byte :: mbyte

rmin = 1.0e10
rmax = -1.0e10
k = 0
do iz = 1,nz8
    do iy = 1,ny8
        do ix = 1,nx8
            k = k + 1
            kbyte = (k-1)/8 + 1
            kbit = k-1 - 8*(kbyte-1)
            mbyte = closedata(kbyte)
            if (btest(mbyte,kbit)) then
                xyz(1) = (ix-1)*voxelsize(1)
                xyz(2) = (iy-1)*voxelsize(2)
                xyz(3) = (iz-1)*voxelsize(3)
                rmin = min(rmin,xyz)
                rmax = max(rmax,xyz)
            endif
        enddo
    enddo
enddo
write(*,*) 'getcloseranges: '
write(*,*) 'rmin: ',rmin
write(*,*) 'rmax: ',rmax
end subroutine

!-----------------------------------------------------------------------------------------
! After selecting the first point randomly from the list, the next two are chosen from
! the set of neighbour blocks (max 27)
!-----------------------------------------------------------------------------------------
subroutine get_triangle(tri)
real :: tri(3,3)
integer :: i, j, k, i1, ib(3), ibfrom(3), ibto(3), jb(3), ns, iseg, kpar = 0

i1 = random_int(1,Nsegpts,kpar)
tri(1,:) = point(i1,:)
ib(1) = (tri(1,1)-rmin(1))/delb(1) + 1
ib(2) = (tri(1,2)-rmin(2))/delb(2) + 1
ib(3) = (tri(1,3)-rmin(3))/delb(3) + 1
!ib = (tri(1,:)-rmin)/delb + 1
ibfrom = max(ib-1,1)
ibto = min(ib+1,NB)
!write(*,*) tri(1,:)
!write(*,*) rmin,delb
!write(*,*) ib
!write(*,*) ibfrom,ibto
do
    do j = 1,3
        jb(j) = random_int(ibfrom(j),ibto(j),kpar)
    enddo
    ns = nseglist(jb(1),jb(2),jb(3))
    if (ns == 0) cycle
    k = random_int(1,ns,kpar)
    exit
enddo
iseg = seglist(jb(1),jb(2),jb(3),k)
tri(2,:) = segment(iseg)%end1
tri(3,:) = segment(iseg)%end2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine create_random_pts(res)
integer :: res
integer :: kpar = 0
integer :: ntr, i, j, ip(3), k, kp, Nspherepts, nodpts, ix, iy, iz
real :: xyz(3), rsum, R(3), p(3), r2, sphr2, tri(3,3)
real, allocatable :: sphere_point(:,:)
logical :: hit, localize = .true.

if (use_sphere) then        ! Need to create a point list for the sphere
    sphr2 = sphere_radius*sphere_radius
    Nspherepts = 0
    do i = 1,Nsegpts
        p = point(i,:)
        R = p - sphere_centre
        r2 = R(1)*R(1) + R(2)*R(2) + R(3)*R(3)
        if (r2 < sphr2) then
            Nspherepts = Nspherepts + 1
        endif
    enddo
    allocate(sphere_point(Nspherepts,3))
    k = 0
    do i = 1,Nsegpts
        p = point(i,:)
        R = p - sphere_centre
        r2 = R(1)*R(1) + R(2)*R(2) + R(3)*R(3)
        if (r2 < sphr2) then
            k = k+1
            sphere_point(k,:) = p
        endif
    enddo
    nodpts = Nspherepts
    np_random = 100000
else
    nodpts = Nsegpts
    np_random = 1000000
endif
ntr = 4
allocate(random_pt(np_random,3))

hit = .false.
if (use_close) then
    hit = in_close(p)
else
    hit = .true.
endif

if (use_close) then
    kp = 0
    do
        ix = random_int(1,N(1),kpar)
        iy = random_int(1,N(2),kpar)
        iz = random_int(1,N(3),kpar)
	    xyz(1) = rmin(1) + (ix-1)*del(1)
        xyz(2) = rmin(2) + (iy-1)*del(2)
	    xyz(3) = rmin(3) + (iz-1)*del(3)
	    if (in_close(xyz)) then
	        kp = kp+1
	        random_pt(kp,:) = xyz
	    endif
	    if (kp == np_random) exit
	enddo
else
    kp = 0
    do i = 1,np_random/ntr
	    j = 1
	    ! generate three random points from the list of vertices
	    ! this creates a triangle that is (almost certainly) completely within the tissue region
    	
	    if (localize .and. .not.use_sphere) then
	        call get_triangle(tri)
	    else
	        do
		        ip(j) = random_int(1,nodpts,kpar)
		        if (j > 1) then 
			        if (ip(j) == ip(1)) cycle
			        if (j == 3 .and. ip(j) == ip(2)) cycle
		        endif
		        j = j+1
		        if (j == 4) exit
	        enddo
	    endif
	    do k = 1,ntr
		    xyz = 0
		    rsum = 0
		    do j = 1,3
			    R(j) = par_uni(kpar)
			    rsum = rsum + R(j)
			    if (use_sphere) then
    			    xyz = xyz + R(j)*sphere_point(ip(j),:)
			    elseif (localize) then
    			    xyz = xyz + R(j)*tri(j,:)
			    else
    			    xyz = xyz + R(j)*point(ip(j),:)
                endif
		    enddo
		    xyz = xyz/rsum
		    ! xyz is a random point inside the triangle
		    kp = kp+1
		    random_pt(kp,:) = xyz
	    enddo
    enddo
endif
res = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine show
integer :: x,y,z

x = N(1)/2
do y = 1,N(2)
	write(*,'(40L2)') (in(x,y,z),z=1,35)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Create a list of all segments, with end points and mid-point, and end radii.
!-----------------------------------------------------------------------------------------
subroutine read_network
integer :: i, j, k, nedges, npoints, iedge, iseg, ip, npts, imid
integer, allocatable :: edge_npts(:)
real, allocatable :: pt(:,:), diam(:)
character*(64) :: line
integer :: ns = 100
real :: test_centre(3) = (/500., 500., 500/)
real :: length, ds, p(3), dp(3)

!amfile = 'G00.am'
!grid_dx = 3
delp = 0.5
nsegments = 0
open(nfam,file=amfile,status='old')
do
	read(nfam,'(a)') line
	if (line(1:11) == 'define EDGE') then
		read(line(12:),'(i)') nedges
		allocate(edge(nedges))
		allocate(edge_npts(nedges))
	elseif (line(1:12) == 'define POINT') then
		read(line(13:),'(i)') npoints
		allocate(pt(npoints,3))
		allocate(diam(npoints))
	endif
	if (line(1:2) == '@3') then
		do iedge = 1,nedges
			read(nfam,*) edge_npts(iedge)
			edge(iedge)%npts = edge_npts(iedge)
			allocate(edge(iedge).pt(edge_npts(iedge)))
			nsegments = nsegments + edge_npts(iedge) - 1
		enddo
	endif
	if (line(1:2) == '@4') then
		do i = 1,npoints
			read(nfam,*) pt(i,:)
!			write(*,*) i,pt(i,:)
		enddo
	endif
	if (line(1:2) == '@5') then
		do i = 1,npoints
			read(nfam,*) diam(i)
!			write(*,*) i,diam(i)
		enddo
		exit
	endif
enddo
close(nfam)
allocate(segment(nsegments))
ip = 0
do iedge = 1,nedges
	do i = 1,edge_npts(iedge)
	    ip = ip+1
	    edge(iedge).pt(i) = ip
	enddo
enddo
if (ip /= npoints) then
    write(*,*) 'ReadAmiraFile: error: inconsistent # of points: ',ip,npoints
    stop
endif
open(nfdist,file=distfile,status='replace')
write(nfdist,'(a,a)') 'Spatialgraph file: ',amfile
if (use_sphere) then
    write(*,'(a,3f8.1,a,f6.1)') 'Sphere with centre: ',sphere_centre,'  radius: ',sphere_radius
    write(nfdist,'(a,3f8.1,a,f6.1)') 'Sphere with centre: ',sphere_centre,'  radius: ',sphere_radius
endif
write(nfdist,'(a,f6.2)') 'grid_dx: ',grid_dx
if (use_random) write(nfdist,'(a)') 'Random sampling points'
if (use_constant_radius) then
    write(*,*) '###########################################################'
    write(*,*) ' TESTING with constant radius !!!!!!!!!!!!!!!!!: ',constant_r
    write(*,*) '###########################################################'
    write(nfdist,*) '###########################################################'
    write(nfdist,*) ' TESTING with constant radius !!!!!!!!!!!!!!!!!: ',constant_r
    write(nfdist,*) '###########################################################'
endif
iseg = 0
do iedge = 1,nedges
    npts = edge_npts(iedge)
    ! testing!!!!!!!!!!!!
!    imid = (npts+1)/2
!    iseg = iseg + 1
!	ip = edge(iedge).pt(1)
!    segment(iseg)%end1 = pt(ip,:)
!    segment(iseg)%r1 = diam(ip)/2
!	ip = edge(iedge).pt(imid)
!    segment(iseg)%end2 = pt(ip,:)
!    segment(iseg)%r2 = diam(ip)/2
!    segment(iseg)%mid = (segment(iseg)%end1 + segment(iseg)%end2)/2
!!    write(*,'(2i6,3f8.1,4x,3f8.1)') iseg,iedge,segment(iseg)%end1,segment(iseg)%end2
!!    write(nfout,'(2i6,3f8.1,4x,3f8.1)') iseg,iedge,segment(iseg)%end1,segment(iseg)%end2
!    iseg = iseg + 1
!	ip = edge(iedge).pt(imid)
!    segment(iseg)%end1 = pt(ip,:)
!    segment(iseg)%r1 = diam(ip)/2
!	ip = edge(iedge).pt(npts)
!    segment(iseg)%end2 = pt(ip,:)
!    segment(iseg)%r2 = diam(ip)/2
!    segment(iseg)%mid = (segment(iseg)%end1 + segment(iseg)%end2)/2
!!    write(*,'(2i6,3f8.1,4x,3f8.1)') iseg,iedge,segment(iseg)%end1,segment(iseg)%end2
!!    write(nfout,'(2i6,3f8.1,4x,3f8.1)') iseg,iedge,segment(iseg)%end1,segment(iseg)%end2
!    cycle
	do i = 1,npts-1
	    ip = edge(iedge).pt(i)
	    iseg = iseg + 1
	    segment(iseg)%end1 = pt(ip,:)
	    segment(iseg)%r1 = diam(ip)/2
	    segment(iseg)%end2 = pt(ip+1,:)
	    segment(iseg)%r2 = diam(ip+1)/2
	    segment(iseg)%mid = (segment(iseg)%end1 + segment(iseg)%end2)/2
!		write(*,'(i6,3f8.2)') iseg,segment(iseg)%mid
	enddo
enddo
nsegments = iseg
write(*,*) 'nsegments: ',nsegments
		
end subroutine

!-----------------------------------------------------------------------------------------
! For each block (as the centre), make list of segments in the neighbourhood (max. 27 blocks),
! add midpoints, flag duplicates.
! The number of blocks in the three directions is NB(:)
! The size of a block in each direction is delb = (rmax - rmin)/NB
! The block index of point p(:) in direction i is (p(i)-rmin(i))/delb(i) + 1 
! The x-range of block(ib1,ib2,ib3) is rmin(1) + (ib1-1)*delb(1) -- rmin(1) + ib1*delb(1)
!-----------------------------------------------------------------------------------------
subroutine make_lists
integer :: iseg, ib(3), nmax

NB(1) = (rmax(1) - rmin(1))/blocksize + 1
if (NB(1) == 1) then
    NB(2) = 1
    NB(3) = 1
else
    NB(2) = NB(1)*(rmax(2) - rmin(2))/(rmax(1) - rmin(1))
    NB(3) = NB(1)*(rmax(3) - rmin(3))/(rmax(1) - rmin(1))
endif
write(*,*) 'NB: ',NB
allocate(nseglist(NB(1),NB(2),NB(3)))
if (NB(1) == 1) then
    delb = blocksize
else
    delb = (rmax - rmin)/NB
endif
write(*,*) 'Block size delb: ',delb

nseglist = 0
do iseg = 1,nsegments
	ib = (segment(iseg)%mid(:) - rmin)/delb + 1
	if (ib(1) < 1 .or. ib(2) < 1 .or. ib(3) < 1) cycle
	if (ib(1) > NB(1) .or. ib(2) > NB(2) .or. ib(3) > NB(3)) cycle
	nseglist(ib(1),ib(2),ib(3))  = nseglist(ib(1),ib(2),ib(3)) + 1
enddo
nmax = maxval(nseglist)
allocate(seglist(NB(1),NB(2),NB(3),nmax))
nseglist = 0
do iseg = 1,nsegments
	ib = (segment(iseg)%mid(:) - rmin)/delb + 1
	if (ib(1) < 1 .or. ib(2) < 1 .or. ib(3) < 1) cycle
	if (ib(1) > NB(1) .or. ib(2) > NB(2) .or. ib(3) > NB(3)) cycle
	nseglist(ib(1),ib(2),ib(3)) = nseglist(ib(1),ib(2),ib(3)) + 1
	seglist(ib(1),ib(2),ib(3),nseglist(ib(1),ib(2),ib(3))) = iseg
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! Get the distance from p(:) to the segment iseg, returning 0 if the point is inside the vessel.
!-----------------------------------------------------------------------------------------
subroutine get_segdist(p,iseg,d)
integer :: iseg
real :: p(3), d
real :: e1(3), e2(3), r1, r2, v(3), w(3), w_mod, v_mod, s, t, r, d2

e1 = segment(iseg)%end1
e2 = segment(iseg)%end2
r1 = segment(iseg)%r1
r2 = segment(iseg)%r2
v = e2 - e1
w = p - e1
w_mod = dot_product(w,w)
v_mod = dot_product(v,v)
if (v_mod == 0 .or. w_mod == 0) then
	d = 0
	return
endif
v_mod = sqrt(v_mod)
w_mod = sqrt(w_mod)
s = dot_product(w,v)/v_mod
if (s <= 0) then
!	d = 0
!	return
	d = w_mod
	r = r1
elseif (s >= v_mod) then
!	d = 0
!	return
	w = p - e2
	w_mod = dot_product(w,w)
	if (w_mod == 0) then
		d = 0
		return
	endif
	w_mod = sqrt(w_mod)
	d = w_mod
	r = r2
else
	t = s/v_mod
	d2 = w_mod*w_mod - s*s
	if (d2 <= 0) then
		d = 0
		return
	endif
	d = sqrt(d2)
	r = (1-t)*r1 + t*r2
endif
if (use_constant_radius) then
    r = constant_r
endif
if (d < r) then
	d = 0
else
	d = d - r
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine par_distribution
integer, parameter :: npdist = 500
integer :: ib(3), ix, iy, iz, np, id, idmax, i, nth, nth0, kpar
real :: p(3), dmin, R
real, allocatable :: pdist(:)
real, allocatable :: par_dmin(:)
logical :: drop, hit
real :: rad, theta, PI, d

allocate(pdist(npdist))
allocate(par_dmin(N(3)))
dbug = .false.
pdist = 0
np = 0
idmax = 0
      
do ix = 1,N(1)
	write(*,'(a,$)') '.'
	do iy = 1,N(2)
		par_dmin = 0
!$omp parallel do private (p,ib,dmin,hit,i,kpar)
        do iz = 1,N(3)
			if (.not.in(ix,iy,iz)) cycle 
            kpar = omp_get_thread_num()
        	p(1) = rmin(1) + (ix-0.5)*del(1)
	        ib(1) = (ix-0.5)*del(1)/delb(1) + 1
		    p(2) = rmin(2) + (iy-0.5)*del(2)
		    ib(2) = (iy-0.5)*del(2)/delb(2) + 1
			p(3) = rmin(3) + (iz-0.5)*del(3)
			ib(3) = (iz-0.5)*del(3)/delb(3) + 1
			! To get off the grid, in case voxel dimensions are integer
			if (jiggle) then
			    do i = 1,3
			        p(i) = p(i) + par_uni(kpar) - 0.5
			    enddo
			endif
			
            call get_mindist(p, ib, dmin, hit)
            
	        if (save_imagedata) then
	            if (dmin > threshold_d) then
!	                call add_voxel(p)
!                    write(*,*) '-----------------------voxel: ',ix,iy,iz,dmin
	                call add_ivoxel(ix,iy,iz)
	            endif
	        endif
            if (.not.hit .or. dmin == 0) cycle
			par_dmin(iz) = dmin
!			if (dbug) write(*,'(a,f6.1)') 'dmin: ',dmin
		enddo
!$omp end parallel do
				
		do iz = 1,N(3)
		    dmin = par_dmin(iz)
		    if (dmin == 0) cycle
			np = np+1
			id = min(npdist,int(dmin/delp + 1))
!			write(*,*) dmin,delp,id
			pdist(id) = pdist(id) + 1
			idmax = max(idmax,id)
		enddo

	enddo
enddo
pdist = pdist/np
write(nfdist,'(a,i10)') '# of segments: ',Nsegments
write(nfdist,'(a,i10)') '# of grid points: ',np_grid
write(nfdist,'(a,i10)') '# of tissue points: ',np
do i = 1,npdist
	write(nfdist,'(f6.2,e12.4)') (i-0.5)*delp,pdist(i)
enddo
close(nfdist)
write(*,*) 'Number of points: ',np
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_mindist(p, ib, dmin, hit)
integer :: ib(3)
real :: p(3), dmin
logical :: hit
integer :: iseg, ibfrom(3), ibto(3), nd, ib1, ib2, ib3, i, j, k
real :: d2, r(3), rad, r2, d, d2min

ibfrom = max(ib-1,1)
ibto = min(ib+1,NB)
!if (dbug) then
!write(*,'(3f6.0,3i6,2x,3i6,2x,3i6)') p,ib,ibfrom,ibto
!write(*,'(a,2f8.0)') 'block x range: ',rmin(1)+(ibfrom(1)-1)*delb(1),rmin(1)+ibto(1)*delb(1)
!write(*,'(a,2f8.0)') 'block y range: ',rmin(2)+(ibfrom(2)-1)*delb(2),rmin(2)+ibto(2)*delb(2)
!write(*,'(a,2f8.0)') 'block z range: ',rmin(3)+(ibfrom(3)-1)*delb(3),rmin(3)+ibto(3)*delb(3)
!endif
hit = .false.
dmin = 1.0e10
do ib1 = ibfrom(1),ibto(1)
	do ib2 = ibfrom(2),ibto(2)
		do ib3 = ibfrom(3),ibto(3)
!		    if (dbug) write(*,*) ib1,ib2,ib3,'   ',nseglist(ib1,ib2,ib3)
            if (nseglist(ib1,ib2,ib3) == 0) cycle
			do k = 1,nseglist(ib1,ib2,ib3)
				iseg = seglist(ib1,ib2,ib3,k)
	            hit = .true.
	            call get_segdist(p,iseg,d)
	            if (d == 0) then
	                dmin = 0
	                return
	            endif
	            dmin = min(dmin,d)
	        enddo
	    enddo
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_mindist1(p, ib, dmin, hit)
integer :: ib(3)
real :: p(3), dmin
logical :: hit
integer :: iseg, ibfrom(3), ibto(3), nd, ib1, ib2, ib3, i, j, k
real :: d2, r(3), rad, r2, d, d2min
integer, parameter :: Nbest = 10
type(segdist_type) :: segdist(Nbest)

ibfrom = max(ib-1,1)
ibto = min(ib+1,NB)
nd = 0
segdist%d2 = 1.0e10
do ib1 = ibfrom(1),ibto(1)
	do ib2 = ibfrom(2),ibto(2)
		do ib3 = ibfrom(3),ibto(3)
			do k = 1,nseglist(ib1,ib2,ib3)
				iseg = seglist(ib1,ib2,ib3,k)
				d2min = 1.0e10
				do j = 1,3
					if (j == 1) then
						r = p - segment(iseg)%end1
						rad = segment(iseg)%r1
					elseif (j == 3) then
						r = p - segment(iseg)%end2
						rad = segment(iseg)%r2
					else
						r = p - segment(iseg)%mid
						rad = (segment(iseg)%r1 + segment(iseg)%r2)/2
					endif
					d2 = dot_product(r,r)
					r2 = rad*rad
					if (d2 < r2) cycle
					d2min = min(d2min,d2)
				enddo
				d2 = d2min
				if (nd == 0) then
					nd = 1
					segdist(nd)%iseg = iseg
					segdist(nd)%d2 = d2
				else
					do i = 1,nd
						if (d2 < segdist(i)%d2) then
							if (i < Nbest) then
								do j = Nbest,i+1,-1
									segdist(j) = segdist(j-1)
								enddo
							endif
							segdist(i)%iseg = iseg
							segdist(i)%d2 = d2
							nd = min(nd+1,Nbest)
							exit
						endif
					enddo
				endif
			enddo
		enddo
	enddo
enddo

hit = .false.
if (nd == 0) then
	return
endif
! We have a list of the nd closest segments (mid-points)
dmin = 1.0e10
do i = 1,nd
	iseg = segdist(i)%iseg
	call get_segdist(p,iseg,d)
	if (d == 0) cycle
	dmin = min(dmin,d)
	hit = .true.
enddo
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine add_voxel(p)
real :: p(3)
integer :: ix, iy, iz

ix = p(1)/grid_dx + offset(1)
if (ix > MAX_NX) then
    write(*,*) 'ix > MAX_NX'
    call exit(9)
endif
nx = max(nx,ix)
if (use_random) nx = max(nx,ix+1)
iy = p(2)/grid_dx + offset(2)
if (iy > MAX_NX) then
    write(*,*) 'iy > MAX_NX'
    call exit(9)
endif
ny = max(ny,iy)
if (use_random) ny = max(ny,iy+1)
iz = p(3)/grid_dx + offset(3)
if (iz > MAX_NX) then
    write(*,*) 'iz > MAX_NX'
    call exit(9)
endif
nz = max(nz,iz)
if (use_random) nz = max(nz,iz+1)
if (use_random) then
    imagedata(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1) = 255
else
    imagedata(ix,iy,iz) = 255
endif
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine add_ivoxel(ix,iy,iz)
!real :: p(3)
integer :: ix, iy, iz
integer :: ixx, iyy, izz

ixx = ix + offset(1)
iyy = iy + offset(2)
izz = iz + offset(3)
!ix = p(1)/grid_dx
!if (ix > MAX_NX) then
!    write(*,*) 'ix > MAX_NX'
!    stop
!endif
nx = max(nx,ixx)
if (use_random) nx = max(nx,ixx+1)
!iy = p(2)/grid_dx
!if (iy > MAX_NX) then
!    write(*,*) 'iy > MAX_NX'
!    stop
!endif
ny = max(ny,iyy)
if (use_random) ny = max(ny,iyy+1)
!iz = p(3)/grid_dx
!if (iz > MAX_NX) then
!    write(*,*) 'iz > MAX_NX'
!    stop
!endif
nz = max(nz,izz)
if (use_random) nz = max(nz,izz+1)
if (use_random) then
    imagedata(ixx-1:ixx+1,iyy-1:iyy+1,izz-1:izz+1) = 255
else
    imagedata(ixx,iyy,izz) = 255
endif
end subroutine

!----------------------------------------------------------------------------------------- 
!#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)] 
!-----------------------------------------------------------------------------------------
subroutine write_imagedata
integer :: ix, iy, iz, offset(3)
integer :: x0, z0

x0 = 100
z0 = 100
do ix = 1,nx
    do iy = 1,ny
        do iz = 1,nz
            if (imagedata(ix,iy,iz) /= 0) then
                x0 = min(x0,ix)
                z0 = min(z0,iz)
            endif
        enddo
    enddo
enddo
write(*,*) 'write_imagedata: first ix, iy: ',x0,z0
nx = nx + 20
ny = ny + 20
nz = nz + 20
open(nfdata, file=datafile, status='replace', form = 'unformatted', access = 'stream')
write(nfdata) nx,ny,nz,(((imagedata(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
close(nfdata)
end subroutine

!----------------------------------------------------------------------------------------- 
!----------------------------------------------------------------------------------------- 
subroutine read_closedata(res)
integer res
integer :: ix, iy, iz, kbit, kbyte, k
byte :: readbyte, arraybyte

res = 0
open(nfclose, file=closefile, status='old', form = 'unformatted', access = 'stream', err=99)
read(nfclose) nxm,nym,nzm,nx8,ny8,nz8,nmbytes
write(*,*) 'closedata dimensions: ',nxm,nym,nzm,nx8,ny8,nz8,nmbytes
if (nmbytes /= nx8*ny8*nz8/8) then
    write(*,*) 'Error: this is not a compressed close data file'
    res = 7
    return
endif
!allocate(closedata(nxm,nym,nzm))
allocate(closedata(nmbytes))
!read(nfclose) (((closedata(ix,iy,iz),ix=1,nxm),iy=1,nym),iz=1,nzm)
read(nfclose) closedata
!write(*,'(20i4)') closedata(1:1000)
!kbit = 0
!kbyte = 0
!arraybyte = 0
!do iz = 1,nzm
!    write(*,*) iz
!    do iy = 1,nym
!        do ix = 1,nxm
!            read(nfclose) readbyte
!            kbit = kbit + 1
!            ! set bit kbit in arraybyte if readbyte != 0
!            if (kbit == 8) then
!                kbyte = kbyte + 1
!                closedata(kbyte) = arraybyte
!                kbit = 0
!                arraybyte = 0
!            endif
!        enddo
!    enddo
!enddo
close(nfclose)
return
99 continue
res = 6
!k = 0
!do kbyte = 1,nmbytes
!    if (closedata(kbyte) /= 0) k = k + 1
!enddo
!write(*,*) 'closedata lit voxels: ',k
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
!if (Mnodes == 1) return
#if defined(OPENMP) || defined(_OPENMP)
write(*,'(a,i2)') 'Requested Mnodes: ',Mnodes
npr = omp_get_num_procs()
write(*,'(a,i2)') 'Machine processors: ',npr

nth = omp_get_max_threads()
write(*,'(a,i2)') 'Max threads available: ',nth
if (nth < Mnodes) then
    Mnodes = nth
    write(*,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
endif

write(*,*) 'Mnodes: ',Mnodes
call omp_set_num_threads(Mnodes)
!!$omp parallel
!nth = omp_get_num_threads()
!write(*,*) 'Threads, max: ',nth,omp_get_max_threads()
!!$omp end parallel
#endif

write(*,*) 'Did omp_initialisation'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_omp
integer :: iz,kpar
real :: R

!$omp parallel do private (R,kpar)
do iz = 1,1000
    kpar = omp_get_thread_num()
    R = par_uni(kpar)
enddo
!$end parallel do
end subroutine

!-----------------------------------------------------------------------------------------
! Extract statistics of the number of vessels per mm^2 from 2D slices.
! For a given plane (parallel to XY, YZ or XZ planes) a vessel segment crosses the plane
! if the vessel ends are on opposite sides of the plane.
! Choose planes in a range of positions about the mean.
! E.g. for planes parallel to the XY plane (normal to the Z axis), the total range of z
! is rmin(3) < z < rmax(3), rmid(3) = (rmin(3) + rmax(3))/2 and rng(3) = rmax(3) - rmin(3)  
! We could set beta = 0.5 and use:
! rmid(3) - beta*rng(3)/2 < z < rmid(3) + beta*rng(3)/2
! Then choose a number of planes in this range, e.g. 100.
! The harder problem is to estimate the area of intersection of the plane with the network.
! This is solved by Delaunay triangulation.  After a set of triangles spanning all the vessel
! intersection points is created, the sum of their areas will be close to the intersection area.
!-----------------------------------------------------------------------------------------
subroutine histology
real :: x, y, z, pos1(3), pos2(3), area, ppos1(3), ppos2(3)
integer :: iseg, np, ntri, ip
real(kind=8) :: p(2)
real(kind=8), allocatable :: table(:,:)
integer, allocatable, dimension ( :, : ) :: triangle_node
logical :: repeat

z = rmid(3)
ppos1 = 0
ppos2 = 0
np = 0
do iseg = 1,Nsegments
    pos1 = segment(iseg)%end1
    pos2 = segment(iseg)%end2
    if ((pos1(3) < z  .and. z < pos2(3)) .or. (pos2(3) < z .and. z < pos1(3))) then
        np = np + 1
    endif
    if (pos1(1) == ppos1(1) .and. pos1(2) == ppos1(2) .and. pos1(3) == ppos1(3) .and. &
        pos2(1) == ppos2(1) .and. pos2(2) == ppos2(2) .and. pos2(3) == ppos2(3)) then
        write(*,'(a,i7,6f8.2)') 'Successive points are the same: iseg: ',iseg,pos1,pos2
        write(nfout,'(a,i7,6f8.2)') 'Successive points are the same: iseg: ',iseg,pos1,pos2
        stop
    endif
    ppos1 = pos1
    ppos2 = pos2
enddo
write(*,*) 'Nsegments, np: ',Nsegments,np
write(nfout,*) 'Nsegments, np: ',Nsegments,np
if (np < 3) return
allocate(table(2,np))
np = 0
do iseg = 1,Nsegments
    pos1 = segment(iseg)%end1
    pos2 = segment(iseg)%end2
    if ((pos1(3) < z  .and. z < pos2(3)) .or. (pos2(3) < z .and. z < pos1(3))) then
        p(1) = pos1(1) + (pos2(1)-pos1(1))*(z - pos1(3))/(pos2(3)-pos1(3))     
        p(2) = pos1(2) + (pos2(2)-pos1(2))*(z - pos1(3))/(pos2(3)-pos1(3))
        repeat = .false.
        do ip = 1,np
            if (p(1)==table(1,ip) .and. p(2)==table(2,ip)) then  
                repeat = .true.
                exit
            endif
        enddo
        if (repeat) cycle 
        np = np + 1
        table(:,np) = p
        write(nfout,'(i6,2f8.2)') np,table(1,np),table(2,np)
!        write(nfout,'(i6,6f8.2)') iseg,pos1,pos2
!        if (np == 1) cycle
!        if (table(1,np) == table(1,np-1) .and. table(2,np) == table(2,np-1)) then
!            write(*,*) 'repeated point: ',np,table(1,np),table(2,np)
!        endif
    endif
enddo
allocate ( triangle_node(3,3*np) )     ! node numbers for the triangles
!call triangulate(np, table, ntri, triangle_node)
!write(*,*) 'ntri: ',ntri
!call total_area(np, table, ntri, triangle_node, area)
deallocate ( table )
deallocate ( triangle_node )

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
logical function checkin(iseg,p)
integer :: iseg
real :: p(3), d

call get_segdist(p,iseg,d)  ! distance from the vessel wall, 0 if inside
if (d == 0) then
    checkin = .true.
else
    checkin = .false.
endif
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine random_distribution1
integer, parameter :: npdist = 500
integer, parameter :: Nbest = 10
integer :: iseg, ib(3), ip, ix, iy, iz, ibfrom(3), ibto(3), nd, ib1, ib2, ib3, i, j, k, np, id,idmax
real :: x, y, z, p(3), d2, v(3), rad, r2, d, dmin, d2min, theta, PI
real, allocatable :: pdist(:)
type(segdist_type) :: segdist(Nbest)
logical :: drop, inside
logical :: check_inside = .true.

allocate(pdist(npdist))
pdist = 0
np = 0
idmax = 0
do ip = 1,np_random
    if (mod(ip,1000) == 0) write(*,'(a,$)') '.'
    x = random_pt(ip,1)
    y = random_pt(ip,2)
    z = random_pt(ip,3)
	p(1) = x 
    p(2) = y
	p(3) = z
    ib(1) = (x-rmin(1))/delb(1) + 1
    ib(2) = (y-rmin(2))/delb(2) + 1
	ib(3) = (z-rmin(3))/delb(3) + 1
	ibfrom = max(ib-1,1)
	ibto = min(ib+1,NB)
	drop = .false.
	nd = 0
	segdist%d2 = 1.0e10
	do ib1 = ibfrom(1),ibto(1)
		if (drop) exit
		do ib2 = ibfrom(2),ibto(2)
			if (drop) exit
			do ib3 = ibfrom(3),ibto(3)
				if (drop) exit
				do k = 1,nseglist(ib1,ib2,ib3)
					iseg = seglist(ib1,ib2,ib3,k)
					if (check_inside) then
					    inside = checkin(iseg,p)
					    if (inside) then
					        drop = .true.
					        exit
					    endif
					endif
					! get min of squared distances to ends, middle
					d2min = 1.0e10
					do j = 1,3
						if (j == 1) then
							v = p - segment(iseg)%end1
							rad = segment(iseg)%r1
						elseif (j == 3) then
							v = p - segment(iseg)%end2
							rad = segment(iseg)%r2
						else
							v = p - segment(iseg)%mid
							rad = (segment(iseg)%r1 + segment(iseg)%r2)/2
						endif
						d2 = dot_product(v,v)
						r2 = rad*rad
						if (d2 < r2) then
						    drop = .true.   ! the point is inside this vessel
						    exit
						endif
						d2min = min(d2min,d2)
					enddo
					if (drop) exit
					d2 = d2min
					if (nd == 0) then
						nd = 1
						segdist(nd)%iseg = iseg
						segdist(nd)%d2 = d2
					else
						do i = 1,nd
							if (d2 < segdist(i)%d2) then
								if (i < Nbest) then
									do j = Nbest,i+1,-1
										segdist(j) = segdist(j-1)
									enddo
								endif
								segdist(i)%iseg = iseg
								segdist(i)%d2 = d2
								nd = min(nd+1,Nbest)
								exit
							endif
						enddo
					endif
				enddo
				if (drop) exit
			enddo
		enddo
	enddo
	if (nd == 0) then
		cycle
	endif
	! We have a list of the nd closest segments (based on ends and mid-points)
	! Note that the point may still lie inside a vessel
	dmin = 1.0e10
!	hit = .false.
	inside = .false.
	do i = 1,nd
		iseg = segdist(i)%iseg
		call get_segdist(p,iseg,d)  ! distance from the vessel wall, 0 if inside
		if (d == 0) then
		    inside = .true.
		    cycle
		endif
		dmin = min(dmin,d)
!		hit = .true.
	enddo
	if (inside) cycle
	np = np+1
	id = min(npdist,int(dmin/delp + 1))
	pdist(id) = pdist(id) + 1
	idmax = max(idmax,id)
enddo
write(*,*)
pdist = pdist/np
!open(nfdist,file=distfile,status='replace')
do i = 1,npdist
	write(nfdist,'(f6.2,e12.4)') (i-0.5)*delp,pdist(i)
enddo
close(nfdist)
deallocate(pdist)
write(*,*) 'Number of points: ',np
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine random_distribution
integer, parameter :: npdist = 500
integer, parameter :: Nbest = 10
integer :: iseg, ib(3), ip, ix, iy, iz, ibfrom(3), ibto(3), nd, ib1, ib2, ib3, i, j, k, np, id, kthread
real :: x, y, z, p(3), d, dmin
integer, allocatable :: np_th(:)
real, allocatable :: pdist(:), pdist_th(:,:)
logical :: drop, hit

grid_dx = 4     ! for save_imagedata
allocate(pdist(npdist))
allocate(pdist_th(npdist,Mnodes))
allocate(np_th(Mnodes))
pdist_th = 0
np_th = 0
!$omp parallel do private (x, y, z, p, ib, dmin, hit, id, kthread)
do ip = 1,np_random
    kthread = omp_get_thread_num()
    write(*,*) kthread
    kthread = kthread + 1
    if (kthread == 1 .and. mod(ip,1000) == 0) write(*,'(a,$)') '.'
    x = random_pt(ip,1)
    y = random_pt(ip,2)
    z = random_pt(ip,3)
	p(1) = x 
    p(2) = y
	p(3) = z
    ib(1) = (x-rmin(1))/delb(1) + 1
    ib(2) = (y-rmin(2))/delb(2) + 1
	ib(3) = (z-rmin(3))/delb(3) + 1

    call get_mindist(p, ib, dmin, hit)

!	ibfrom = max(ib-1,1)
!	ibto = min(ib+1,NB)
!	drop = .false.
!	hit = .false.
!	dmin = 1.0e10
!	do ib1 = ibfrom(1),ibto(1)
!		do ib2 = ibfrom(2),ibto(2)
!			do ib3 = ibfrom(3),ibto(3)
!				do k = 1,nseglist(ib1,ib2,ib3)
!					iseg = seglist(ib1,ib2,ib3,k)
!					call get_segdist(p,iseg,d)
!				    if (d == 0) then
!				        drop = .true.
!				        exit
!				    endif
!				    hit = .true.
!		            dmin = min(dmin,d)
!	            enddo
!	            if (drop) exit
!	        enddo
!            if (drop) exit
!	    enddo
!        if (drop) exit
!	enddo
!	if (drop .or. .not.hit) cycle
	if (save_imagedata) then
	    if (dmin > threshold_d) then
	        call add_voxel(p)
	    endif
	endif
    if (.not.hit .or. dmin == 0) cycle
	np_th(kthread) = np_th(kthread) + 1
	id = min(npdist,int(dmin/delp + 1))
	pdist_th(id,kthread) = pdist_th(id,kthread) + 1
enddo
np = 0
pdist = 0
do i = 1,Mnodes
    np = np + np_th(i)
    pdist(:) = pdist(:) + pdist_th(:,i)
enddo
write(*,*)
pdist = pdist/np
write(nfdist,'(a,i10)') '# of segments: ',Nsegments
write(nfdist,'(a,i10)') '# of random points: ',np_random
write(nfdist,'(a,i10)') '# of tissue points: ',np
do i = 1,npdist
	write(nfdist,'(f6.2,e12.4)') (i-0.5)*delp,pdist(i)
enddo
close(nfdist)
deallocate(pdist)
write(*,*) 'Number of points: ',np
end subroutine

end module

!-----------------------------------------------------------------------------------------
! Return error codes:
! 1  Command line error (1)
! 2  Command line error (2)
! 3  Command line error (3)
! 4  Allocation failure on array in(), increase grid_dx
! 5  Indx out of range
! 6  Open failure on close file
! 7  Supplied close file has incorrect format - did you specify .bin file?
!-----------------------------------------------------------------------------------------
program main
use data_mod
integer :: res
logical :: ok

res = 0

open(nfout, file='histo.out', status='replace')
call input(res)
if (res /= 0) then
	call exit(res)
endif
call omp_initialisation(ok)
call read_network
call setup(res)
if (res /= 0) then
    call exit(res)
endif
!call test_delaunay
!call histology
!call exit(res)

call make_lists
if (use_random) then
    call create_random_pts(res)
    write(*,*) 'Generated random points: ',np_random
else
    call create_in(res)
endif
if (res /= 0) then
	call exit(res)
endif
!call make_lists
!call show
if (use_random) then
    call random_distribution
else
    call par_distribution
    !call distribution
endif
if (save_imagedata) then
    call write_imagedata
endif
call exit(res)
end
