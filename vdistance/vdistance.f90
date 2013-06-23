! To estimate the probability distribution of distance from tissue to the nearest blood vessel.
! Uses points within the convex hull of a random set of points (?)

module data_mod
use par_zig_mod

implicit none

type segment_type
	real :: end1(3), mid(3), end2(3)
	real :: r1, r2
end type

type segdist_type
	integer :: iseg
	real :: d2
end type

real, allocatable :: point(:,:)

logical(2), allocatable :: in(:,:,:)

type(segment_type), allocatable :: segment(:)	! this is the list of all segments
integer, allocatable :: seglist(:,:,:,:)		! for each (ib,jb,kb), a list of segment numbers
integer, allocatable :: nseglist(:,:,:)			! for each (ib,jb,kb), number of segments in the seglist

real :: grid_dx, delp
real :: rmin(3), rmax(3), del(3), delb(3), centre(3), radius
integer :: Npts, N(3), NB(3), Nsegments
character*(128) :: amfile, distfile
logical :: use_sphere

integer, parameter :: nfam = 10, nfdist=11
real, parameter :: blocksize = 60
integer, parameter :: seed(2) = (/12345, 67891/)
integer, parameter :: test = 0

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine input(res)
integer :: res
integer :: nlen, cnt, i, status
character*(1024) :: c, progname

call get_command (c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    res = 1
    return
end if
!write (*,*) 'command line = ', c(1:len)
call get_command_argument (0, c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    res = 2
    return
end if
write (*,*) 'command name = ', c(1:nlen)
progname = c(1:nlen)
cnt = command_argument_count ()
write (*,*) 'number of command arguments = ', cnt
res = 0
if (cnt == 3) then
	use_sphere = .false.
elseif (cnt == 7) then
	use_sphere = .true.
else
	res = 3
    write(*,*) 'Use either: ',trim(progname),' amfile distfile grid_dx '
    write(*,*) ' to analyze the whole network'
    write(*,*) 'or: ',trim(progname),' amfile distfile grid_dx x0 y0 z0 R'
    write(*,*) ' to analyze a spherical subregion with centre (x0,y0,z0), radius R'
    return
endif

do i = 1, cnt
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        res = 1
        return
    end if
    if (i == 1) then
        amfile = c(1:nlen)
        write(*,*) 'Network file: ',amfile
    elseif (i == 2) then
        distfile = c(1:nlen)
        write(*,*) 'Distribution file: ',distfile
    elseif (i == 3) then
        read(c(1:nlen),*) grid_dx																
        write(*,*) 'grid_dx: ',grid_dx
    elseif (i == 4) then
        read(c(1:nlen),*) centre(1)																
    elseif (i == 5) then
        read(c(1:nlen),*) centre(2)																
    elseif (i == 6) then
        read(c(1:nlen),*) centre(3)																
    elseif (i == 7) then
        read(c(1:nlen),*) radius																
    endif
end do
if (use_sphere) then
	write(*,'(a,4f8.2)') 'centre, radius: ',centre,radius
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = 1
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
subroutine setup
integer :: k, i, ierr

write(*,*) "call rng_initialisation"
call rng_initialisation
write(*,*) "did rng_initialisation"

Npts = 2*nsegments
allocate(point(Npts,3))
rmin = 1.0e10
rmax = -1.0e10
k = 0
do i = 1,nsegments
	k = k+1
	point(k,:) = segment(i)%end1
	k = k+1
	point(k,:) = segment(i)%end2
enddo

do k = 1,Npts
	rmin = min(rmin,point(k,:))
	rmax = max(rmax,point(k,:))
enddo
if (test == 1) then
	rmin(2:3) = 450
	rmax(2:3) = 550
endif
write(*,'(a,3f8.2)') 'rmin: ',rmin
write(*,'(a,3f8.2)') 'rmax: ',rmax
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine create_in(res)
integer :: res
integer :: i, j, k, ip(3), indx(3), ix, iy, iz, ierr, kpar=0
real :: R(3), rsum, xyz(3), v(3)

res = 0
if (use_sphere) then
!	centre = (rmin + rmax)/2
!	radius = 1.0e10
!	do i = 1,3
!		radius = min(radius,centre(i)-rmin(i))
!	enddo
!	radius = radius/2
	write(*,'(a,4f8.2)') 'Using sphere: centre, radius: ',centre,radius
	rmin = centre - radius
	rmax = centre + radius
endif
rmin = rmin - 1
rmax = rmax + 1
del = grid_dx
N = (rmax - rmin)/del + 1
write(*,*) 'Dimensions of array in(:,:,:): ',N
allocate(in(N(1),N(2),N(3)),stat=ierr)
if (ierr /= 0) then
	write(*,*) 'Allocation failed on array in(:,:,:).  Increase grid_dx'
	res = 4
	return
endif
in = .false.
del = (rmax-rmin)/(N-1)
write(*,*) 'del: ',del
write(*,*) 'N: ',N
if (test == 1) then
	do ix = 1,N(1)
		do iy = 1,N(2)
			do iz = 1,N(3)
				in(ix,iy,iz) = .true.
			enddo
		enddo
	enddo
elseif (use_sphere) then
	do ix = 1,N(1)
		xyz(1) = rmin(1) + (ix-1)*del(1)
		do iy = 1,N(2)
			xyz(2) = rmin(2) + (iy-1)*del(2)
			do iz = 1,N(3)
				xyz(3) = rmin(3) + (iz-1)*del(3)
				v = xyz - centre
				if (dot_product(v,v) < radius*radius) then
					in(ix,iy,iz) = .true.
				else
					in(ix,iy,iz) = .false.
				endif
			enddo
		enddo
	enddo
else
	do i = 1,10000000
		j = 1
		do
			ip(j) = random_int(1,Npts,kpar)
			if (j > 1) then 
				if (ip(j) == ip(1)) cycle
				if (j == 3 .and. ip(j) == ip(2)) cycle
			endif
			j = j+1
			if (j == 4) exit
		enddo
		do k = 1,10
			xyz = 0
			rsum = 0
			do j = 1,3
				R(j) = par_uni(kpar)
				rsum = rsum + R(j)
				xyz = xyz + R(j)*point(ip(j),:)
			enddo
			xyz = xyz/rsum
			do j = 1,3
				indx(j) = (xyz(j)-rmin(j))/del(j) + 1
				if (indx(j) > N(j)) then
					write(*,*) 'indx out of range:'
					write(*,*) 'point: ',point(ip(j),:)
					write(*,*) 'R: ',R/rsum
					write(*,*) 'xyz: ',xyz
					res = 5
					return
				endif
			enddo
			in(indx(1),indx(2),indx(3)) = .true.
		enddo
	enddo
endif
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
integer :: i, j, k, nedges, npoints, iedge, iseg
integer, allocatable :: edgepts(:)
real, allocatable :: pt(:,:), diam(:)
character*(64) :: line
integer :: ns = 100
real :: test_radius = 5
real :: test_centre(3) = (/500., 500., 500/)
real :: length, ds, p(3), dp(3)

!amfile = 'G00.am'
!grid_dx = 3
delp = 0.5
if (test == 1) then		! create a simple "jack" network
	nsegments = 3*ns
	allocate(segment(nsegments))
	length = 2*centre(1)
	ds = length/ns
	k = 0
	do i = 1,1
		p = test_centre
		dp = 0
		dp(i) = ds
		do j = 1,ns
			p(i) = (j-1)*ds
			k = k+1
			segment(k)%end1 = p
			segment(k)%mid = p + dp/2
			segment(k)%end2 = p + dp
!			write(*,'(i4,3f6.1,2x,3f6.1)') k,segment(k)%end1,segment(k)%end2
		enddo
	enddo
	segment(:)%r1 = test_radius
	segment(:)%r2 = test_radius
	grid_dx = 0.42335
	return
endif

!write(*,*) 'Enter AM file name'
!read(*,*) amfile
!write(*,*) 'Enter grid dx in um (e.g. 3)'
!read(*,*) grid_dx
!write(*,*) 'Enter distance resolution required in um (e.g. 0.5)'
!read(*,*) delp
nsegments = 0
open(nfam,file=amfile,status='old')
do
	read(nfam,'(a)') line
	if (line(1:11) == 'define EDGE') then
		read(line(12:),'(i)') nedges
		allocate(edgepts(nedges))
	elseif (line(1:12) == 'define POINT') then
		read(line(13:),'(i)') npoints
		allocate(pt(npoints,3))
		allocate(diam(npoints))
	endif
	if (line(1:2) == '@3') then
		do i = 1,nedges
			read(nfam,*) edgepts(i)
			nsegments = nsegments + edgepts(i) - 1
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
iseg = 0
do iedge = 1,nedges
	do i = 1,edgepts(iedge)-1
		iseg = iseg + 1
		segment(iseg)%end1 = pt(iseg,:)
!		segment(iseg)%r1 = diam(iseg)/2
		segment(iseg)%r1 = 5
		segment(iseg)%end2 = pt(iseg+1,:)
!		segment(iseg)%r2 = diam(iseg+1)/2
		segment(iseg)%r2 = 5
		segment(iseg)%mid = (segment(iseg)%end1 + segment(iseg)%end2)/2
!		write(*,'(i6,3f8.2)') iseg,segment(iseg)%mid
	enddo
enddo
write(*,*) 'nsegments: ',nsegments
		
end subroutine

!-----------------------------------------------------------------------------------------
! For each block (as the centre), make list of segments in the neighbourhood (max. 27 blocks),
! add midpoints, flag duplicates.
! The number of blocks in the three directions is NB(:)
! The size of a block in each direction is delB = (rmax - rmin)/NB
! The block index of point p(:) in direction i is p(i)/delb(i) + 1 
!-----------------------------------------------------------------------------------------
subroutine make_lists
integer :: iseg, ib(3), nmax

!NB(1) = 40
NB(1) = (rmax(1) - rmin(1))/blocksize + 1
NB(2) = NB(1)*(rmax(2) - rmin(2))/(rmax(1) - rmin(1))
NB(3) = NB(1)*(rmax(3) - rmin(3))/(rmax(1) - rmin(1))
write(*,*) 'NB: ',NB
allocate(nseglist(NB(1),NB(2),NB(3)))
delb = (rmax - rmin)/NB
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
	nseglist(ib(1),ib(2),ib(3))  = nseglist(ib(1),ib(2),ib(3)) + 1
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
if (d < r) then
	d = 0
else
	d = d - r
endif
end subroutine

!-----------------------------------------------------------------------------------------
! For each inside point:
! 1. Determine the Nbest (10, say) closest segments, using the appropriate segment list.
! 2. Find the distance from each segment vessel wall.
! 3. Choose the minumum distance = centreline distance - vessel radius
!-----------------------------------------------------------------------------------------
subroutine distribution
integer, parameter :: npdist = 500
integer, parameter :: Nbest = 10
integer :: iseg, ib(3), ix, iy, iz, ibfrom(3), ibto(3), nd, ib1, ib2, ib3, i, j, k, np, id,idmax
real :: p(3), d2, r(3), rad, r2, d, dmin, d2min
real, allocatable :: pdist(:)
type(segdist_type) :: segdist(Nbest)
logical :: hit
!real, parameter :: delp = 0.5

allocate(pdist(npdist))
pdist = 0
np = 0
idmax = 0
do ix = 1,N(1)
	write(*,'(a,$)') '.'
	p(1) = rmin(1) + (ix-0.5)*del(1)
	ib(1) = (ix-0.5)*del(1)/delb(1) + 1
	do iy = 1,N(2)			
		p(2) = rmin(2) + (iy-0.5)*del(2)
		ib(2) = (iy-0.5)*del(2)/delb(2) + 1
		do iz = 1,N(3)
			if (.not.in(ix,iy,iz)) cycle
			p(3) = rmin(3) + (iz-0.5)*del(3)
			ib(3) = (iz-0.5)*del(3)/delb(3) + 1
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
			if (nd == 0) then
				cycle
			endif
			! We have a list of the nd closest segments (mid-points)
			dmin = 1.0e10
			hit = .false.
			do i = 1,nd
				iseg = segdist(i)%iseg
				call get_segdist(p,iseg,d)
!				do j = 1,3
!					if (j == 1) then
!						r = p - segment(iseg)%end1
!						radius = segment(iseg)%r1
!					elseif (j == 2) then
!						r = p - segment(iseg)%mid
!						radius = (segment(iseg)%r1 + segment(iseg)%r2)/2
!					else
!						r = p - segment(iseg)%end2
!						radius = segment(iseg)%r2
!					endif
!					d2 = dot_product(r,r)
!					r2 = radius*radius
!					if (d2 < r2) cycle
!					d = sqrt(d2) - radius
!					enddo
				if (d == 0) cycle
				dmin = min(dmin,d)
!				write(*,*) i,nd,d,dmin
				hit = .true.
			enddo
			if (.not.hit) cycle
			np = np+1
			id = min(npdist,int(dmin/delp + 1))
!			write(*,*) dmin,delp,id
			pdist(id) = pdist(id) + 1
			idmax = max(idmax,id)
!			if (dmin < 0.5) then
!				write(*,'(3i5,e12.3,2i4)') ix,iy,iz,dmin,id,idmax	
!			endif
		enddo
	enddo
enddo
write(*,*)
pdist = pdist/np
open(nfdist,file=distfile,status='replace')
do i = 1,npdist
	write(nfdist,'(f6.2,e12.4)') (i-0.5)*delp,pdist(i)
enddo
close(nfdist)
write(*,*) 'Number of points: ',np
end subroutine

end module

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
program main
use data_mod
integer :: res

res = 0
call input(res)
if (res /= 0) then
	call exit(res)
endif
call read_network
call setup
call create_in(res)
if (res /= 0) then
	call exit(res)
endif
!call show
call make_lists
call distribution
call exit(res)
end
