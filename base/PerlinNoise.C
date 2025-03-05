#include <cmath>
#include <iostream>

#include "PerlinNoise.h"
#include "Vector.h"

namespace lux
{

void PerlinNoise::init()
{
	int permutation[] = { 151,160,137,91,90,15,
		131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
		190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
		88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
		77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
		102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
		135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
		5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
		223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
		129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
		251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
		49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
		138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
	     };

	for (int i=0; i<256; i++)
	{
		p[256+i] = permutation[i];
		p[i] = permutation[i];
	}

}

const double PerlinNoise::fade(double t) const { return std::pow(t, 3) * (t * (t * 6 -15) + 10); }

const double PerlinNoise::lerp(double t, double a, double b) const { return a + t * (b - a); }

const double PerlinNoise::grad(int hash, double x, double y, double z) const
{
	int h = hash & 15;

	// Due to my stubborness in an ongoing debate over ternary operators, I will remain using the original implementation's ternary operators
	// I will provide comments explaining the "if statement" version though.
	
	double u = h<8 ? x : y;
	/* equivalent to:
 	 * 
 	 * if(h < 8) u = x;
 	 * else u = y;
 	 */

	double v = h<4 ? y : h==12||h==14 ? x : z;
	/* equivalent to:
	 * 
	 * if( h<4 ) 
	 * 	v = y;
	 * else
	 * {
	 * 	if(h==12 || h==14)
	 * 		v = x;
	 * 	else
	 * 		v = z;
	 * }
	 */

	return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
	/* equivalent to:
	 *
	 * float x1, x2;
	 *
	 * if((h&1) == 0) x1 = u;
	 * else		  x1 = -u;
	 *
	 * if((h&2) == 0) x2 = v;
	 * else		  x2 = -v;
	 *
	 * return x1 + x2;
	 */

}

const float PerlinNoise::eval( const Vector& vec ) const
{
	int X = (int)std::floor(vec.X()) & 255,
	    Y = (int)std::floor(vec.Y()) & 255,
	    Z = (int)std::floor(vec.Z()) & 255;
	
	float x = vec.X() - std::floor(vec.X());
	float y = vec.Y() - std::floor(vec.Y());
	float z = vec.Z() - std::floor(vec.Z());

	double	u = fade(x),
		v = fade(y),
		w = fade(z);

	int A = p[X  ]+Y, AA = p[A]+Z, AB = p[A+1]+Z,
	    B = p[X+1]+Y, BA = p[B]+Z, BB = p[B+1]+Z;

	return lerp(w, lerp(v, 	lerp(u, grad(p[AA ], x , y , z ),
					grad(p[BA ], x-1, y , z )), 
				lerp(u, grad(p[AB ], x , y-1, z ), 
					grad(p[BB ], x-1, y-1, z ))),
			lerp(v, lerp(u, grad(p[AA+1], x , y , z-1 ), 
					grad(p[BA+1], x-1, y , z-1 )), 
				lerp(u, grad(p[AB+1], x , y-1, z-1 ),
					grad(p[BB+1], x-1, y-1, z-1 ))));

}

const float PerlinNoise::eval( float arg ) const
{
	int X = (int)std::floor(arg) & 255;
	float x = std::floor(x);

	double	u = fade(x);

	int A = p[X  ], AA = p[A],
	    B = p[X+1], BA = p[B];

	return lerp(u,	grad(p[AA ], x , 0 , 0 ),
			grad(p[BA ], x-1, 0 , 0 )); 

}

const float PerlinNoiseGustavson::eval(float arg) const
{
	float x = arg;
	int ix0, ix1;
	float fx0, fx1;
	float s, n0, n1;

	ix0 = floor( x ); // Integer part of x
	fx0 = x - ix0;       // Fractional part of x
	fx1 = fx0 - 1.0f;
	ix1 = ( ix0+1 ) & 0xff;
	ix0 = ix0 & 0xff;    // Wrap to 0..255

	s = fade( fx0 );

	n0 = grad1( perm[ ix0 ], fx0 );
	n1 = grad1( perm[ ix1 ], fx1 );
	return 0.188f * ( lerp( s, n0, n1 ) );
}

const float PerlinNoiseGustavson::eval( const Vector& vec ) const
{

	float x = vec.X();
	float y = vec.Y();
	float z = vec.Z();
	float w = time;


	int ix0, iy0, ix1, iy1, iz0, iz1, iw0, iw1;
	float fx0, fy0, fz0, fx1, fy1, fz1, fw0, fw1;
	float s, t, r, q;
	float nxyz0, nxyz1, nxy0, nxy1, nx0, nx1, n0, n1;

	ix0 = floor( x ); // Integer part of x
	iy0 = floor( y ); // Integer part of y
	iz0 = floor( z ); // Integer part of z
	iw0 = floor( w ); // Integer part of w
	fx0 = x - ix0;        // Fractional part of x
	fy0 = y - iy0;        // Fractional part of y
	fz0 = z - iz0;        // Fractional part of z
	fw0 = w - iw0;        // Fractional part of w
	fx1 = fx0 - 1.0f;
	fy1 = fy0 - 1.0f;
	fz1 = fz0 - 1.0f;
	fw1 = fw0 - 1.0f;
	ix1 = ( ix0 + 1 ) & 0xff; // Wrap to 0..255
	iy1 = ( iy0 + 1 ) & 0xff;
	iz1 = ( iz0 + 1 ) & 0xff;
	iw1 = ( iw0 + 1 ) & 0xff;
	ix0 = ix0 & 0xff;
	iy0 = iy0 & 0xff;
	iz0 = iz0 & 0xff;
	iw0 = iw0 & 0xff;
	
	q = fade( fw0 );
	r = fade( fz0 );
	t = fade( fy0 );
	s = fade( fx0 );


	nxyz0 = grad4(perm[ix0 + perm[iy0 + perm[iz0 + perm[iw0]]]], fx0, fy0, fz0, fw0);
	nxyz1 = grad4(perm[ix0 + perm[iy0 + perm[iz0 + perm[iw1]]]], fx0, fy0, fz0, fw1);
	nxy0 = lerp( q, nxyz0, nxyz1 );
	    
	nxyz0 = grad4(perm[ix0 + perm[iy0 + perm[iz1 + perm[iw0]]]], fx0, fy0, fz1, fw0);
	nxyz1 = grad4(perm[ix0 + perm[iy0 + perm[iz1 + perm[iw1]]]], fx0, fy0, fz1, fw1);
	nxy1 = lerp( q, nxyz0, nxyz1 );
	    
	nx0 = lerp ( r, nxy0, nxy1 );

	nxyz0 = grad4(perm[ix0 + perm[iy1 + perm[iz0 + perm[iw0]]]], fx0, fy1, fz0, fw0);
	nxyz1 = grad4(perm[ix0 + perm[iy1 + perm[iz0 + perm[iw1]]]], fx0, fy1, fz0, fw1);
	nxy0 = lerp( q, nxyz0, nxyz1 );
	    
	nxyz0 = grad4(perm[ix0 + perm[iy1 + perm[iz1 + perm[iw0]]]], fx0, fy1, fz1, fw0);
	nxyz1 = grad4(perm[ix0 + perm[iy1 + perm[iz1 + perm[iw1]]]], fx0, fy1, fz1, fw1);
	nxy1 = lerp( q, nxyz0, nxyz1 );

	nx1 = lerp ( r, nxy0, nxy1 );

	n0 = lerp( t, nx0, nx1 );

	nxyz0 = grad4(perm[ix1 + perm[iy0 + perm[iz0 + perm[iw0]]]], fx1, fy0, fz0, fw0);
	nxyz1 = grad4(perm[ix1 + perm[iy0 + perm[iz0 + perm[iw1]]]], fx1, fy0, fz0, fw1);
	nxy0 = lerp( q, nxyz0, nxyz1 );
	    
	nxyz0 = grad4(perm[ix1 + perm[iy0 + perm[iz1 + perm[iw0]]]], fx1, fy0, fz1, fw0);
	nxyz1 = grad4(perm[ix1 + perm[iy0 + perm[iz1 + perm[iw1]]]], fx1, fy0, fz1, fw1);
	nxy1 = lerp( q, nxyz0, nxyz1 );

	nx0 = lerp ( r, nxy0, nxy1 );

	nxyz0 = grad4(perm[ix1 + perm[iy1 + perm[iz0 + perm[iw0]]]], fx1, fy1, fz0, fw0);
	nxyz1 = grad4(perm[ix1 + perm[iy1 + perm[iz0 + perm[iw1]]]], fx1, fy1, fz0, fw1);
	nxy0 = lerp( q, nxyz0, nxyz1 );
	    
	nxyz0 = grad4(perm[ix1 + perm[iy1 + perm[iz1 + perm[iw0]]]], fx1, fy1, fz1, fw0);
	nxyz1 = grad4(perm[ix1 + perm[iy1 + perm[iz1 + perm[iw1]]]], fx1, fy1, fz1, fw1);
	nxy1 = lerp( q, nxyz0, nxyz1 );

	nx1 = lerp ( r, nxy0, nxy1 );

	n1 = lerp( t, nx0, nx1 );

	return 0.87f * ( lerp( s, n0, n1 ) );	
}

void PerlinNoiseGustavson::init()
{
	int permutation[] = { 151,160,137,91,90,15,
		131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
		190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
		88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
		77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
		102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
		135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
		5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
		223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
		129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
		251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
		49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
		138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
	     };

	for (int i=0; i<256; i++)
	{
		perm[256+i] = permutation[i];
		perm[i] = permutation[i];
	}

}

}

