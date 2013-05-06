/*
 * To carry out fast averaging
 */

#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]

int width, height, depth;
int radius = 3;

// To keep a cube of side N (N = 2*R+1) completely within the image boundaries,
// R-1 <= x <= width-R, etc.
// One way to traverse the image is to start at the min x, y and z, then traverse the x-y 
// plane by incrementing x, then increment y and sweep again.  The first sum at the start of
// the sweep can be stored for reuse.
// When the plane is covered, z is incremented.  Again, the starting sum on the x-y plane
// can be saved for reuse.
// In this was at each (x,y,z) to update the sum requires only ne*ne voxel values to be accessed.
// Note that this assumes that the sum of the slice to be subtracted is also saved.
// Nomenclature:
//   cube = NxNxN block that we are averaging over
//   XYslice = slice of the cube in the xy plane
//   YZslice = slice of the cube in the yz plane
//   ZXslice = slice of the cube in the zx plane

N = 2*R + 1;

for (z=R; z<=depth-R-1; z++)
{
	for (y=R; y<=height-R-1; y++)
	{
		sumXYZ = 0;
		for (x=0; x<N; x++)
		{
			sumYZ[x] = 0;
			for (y=0; y<N; y++)
			{
				for (dz=-R; dz<=R; dz++)
				{
					sumYZ[x] += V(x,y,z+dz)
				}
			}
			sumXYZ += sumYZ[x];
		}
		for (x=R; x<=width-R-1; x++)
		{
			if (x > R)
			{
				sum = 0;
				for (dy=-R; dy<=R; dy++)
				{
					for (dz=-R; dz<=R; dz++)
					{
						sum += V(x,y+dy,z+dz)
					}
				}
				sumXYZ += sum - sumYZ[0];
				for (k=1; k<N; k++)
					sumYZ[k-1] = sumYZ[k];
				sumYZ[N-1] = sum;
			}
			A(x,y,z) = sumXYZ/(N*N*N);
		}
	}
}


