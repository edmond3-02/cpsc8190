#include <cstdio>
#include <cstdlib>

#include "Stamp.h"
#include "Volume.h"
#include "ProgressMeter.h"
#include "Noise.h"
#include "Fields.h"

namespace lux {

void StampField( VolumeGrid<float>& grid, Volume<float>* field )
{
	return;
}

void StampField( ScalarGrid& grid, ScalarField& field )
{
	//ProgressMeter prog((grid->nx()*grid->ny())*2, "Stamp Field" );
	#pragma omp parallel for
	for(int i=0; i<grid->nx();i++)
	{
		for(int j=0; j<grid->ny();j++)
		{
	//		prog.update();
			for(int k=0; k<grid->nz();k++)
			{
				Vector pos = grid->evalP(i,j,k);
				float val  = field->eval(pos);
				grid->set(i,j,k,val);
			}
		}
	}
	
}

void StampField( VectorGrid& grid, VectorField& field )
{
	ProgressMeter prog(grid->nx(), "Stamp Field" );
	for(int i=0; i<grid->nx();i++)
	{
		prog.update();
		#pragma omp parallel for
		for(int j=0; j<grid->ny();j++)
		{
			for(int k=0; k<grid->nz();k++)
			{
				Vector pos = grid->evalP(i,j,k);
				Vector val  = field->eval(pos);
				grid->set(i,j,k,val);
			}
		}
	}
}


void StampField( ScalarGrid& grid, ScalarField& field, int bandwidth)
{
	ProgressMeter prog((grid->nx()*grid->ny())*2, "Stamp Field" );

	float d = (grid->dx() + grid->dy() + grid->dz())/3.0;
	for(int i=0; i<grid->nx();i++)
	{
		for(int j=0; j<grid->ny();j++)
		{
			prog.update();
			for(int k=0; k<grid->nz();k++)
			{
				Vector pos = grid->evalP(i,j,k);
				float val  = field->eval(pos);
				if(val > (bandwidth * -d))
					grid->set(i,j,k,val);
			}
		}
	}

	for(int i=0; i<grid->nx();i++)
	{
		for(int j=0; j<grid->ny();j++)
		{
			prog.update();
			for(int k=0; k<grid->nz();k++)
			{
				Vector pos = grid->evalP(i,j,k);
				float val  = field->eval(pos);
				if( grid->goodBlock(i,j,k) )
					grid->set(i,j,k,val);
			}
		}
	}
	
}


void StampPyro( ScalarGrid& grid, const AnchorChain& particles )
{
	grid->setDefVal(0.0);
	for(int i=0; i<particles.size(); i++)
	{
		ScalarField pyro = RadialPyroclast(	particles[i].P, particles[i].radius, particles[i].amplitude,
							particles[i].octaves, particles[i].frequency, particles[i].roughness, 
							particles[i].offset, particles[i].time, particles[i].gamma);
		for(int x=0; x<grid->nx();x++)
		{
			for(int y=0; y<grid->ny();y++)
			{
				for(int z=0; z<grid->nz();z++)
				{
					Vector pos = grid->evalP(x,y,z);

					float val  = pyro->eval(pos);
					grid->set(x,y,z,val);
				}
			}
		}
	}
}

void StampNoise( ScalarGrid& grid, const AnchorChain& particles )
{
	for(int i=0; i<particles.size(); i++)
	{
		std::cout << "Particle " << i << " at " << particles[i].P.__str__() << std::endl;
		NoiseMachine noise = perlin(particles[i]);
		stamp_noise(grid, noise, particles[i].P, particles[i].radius, particles[i].falloff);
	}
}


void StampPointWisps( ScalarGrid& grid, const AnchorChain& particles )
{
	for(int i=0; i<particles.size(); i++)
	{
		// seed random
		srand(particles[i].seed);

		// FPSN 1
		NoiseMachine fspn1 = perlin(particles[i]);

		// FSPN 2
		Noise_t param;
		param.octaves = particles[i].wispOctaves;
		param.fjump = particles[i].wispFjump;
		param.roughness = particles[i].wispRoughness;
		param.frequency = particles[i].wispFreq;
		param.translate = particles[i].wispTranslate;
		param.offset = particles[i].wispOffset;
		NoiseMachine fspn2 = perlin(param);

		for(int j=0; j<particles[i].nbWisps; j++)
		{
			float x = (float(rand())/RAND_MAX)*2 - 1;
			float y = (float(rand())/RAND_MAX)*2 - 1;
			float z = (float(rand())/RAND_MAX)*2 - 1;

			Vector P0(x, y, z);
			Vector P1 = P0.unitvector();

			float R = std::pow(std::abs(fspn1->eval(P0)),particles[i].wispRadialGroup);

			Vector P3 = P1 * R;

			Vector tbn_adj = P3.X() * particles[i].tangent
				 	+P3.Y() * particles[i].normal
					+P3.Z() * particles[i].binormal;
			Vector P4 = particles[i].P + particles[i].pscale * tbn_adj;

			Vector offset	(1, 0, 0);
			Vector D	(fspn2->eval(P3), fspn2->eval(P3 + offset), fspn2->eval(P3 - offset));

			Vector P5 = P4 + D;

			int   ix, iix, iy, iiy, iz, iiz;
			float wx, wwx, wy, wwy, wz, wwz;
			
			grid->getLinearInterpolation(P5, ix, iix, wx, wwx, iy, iiy, wy, wwy, iz, iiz, wz, wwz);

			grid->getGridIndex(P5, ix, iy, iz);

			//std::cout << "weight sum = " << 	   wx*wy*wz +  wx*wwy*wz +  wx*wy*wwz +  wx*wwy*wwz
			//					+ wwx*wy*wz + wwx*wwy*wz + wwx*wy*wwz + wwx*wwy*wwz << std::endl;
//			std::cout << "Position: " << j << " " << P5.__str__() << std::endl;
			grid->set( ix,  iy,  iz, grid->get( ix,  iy,  iz) + (1.0 * wwx*wwy*wwz) );
			grid->set( ix, iiy,  iz, grid->get( ix, iiy,  iz) + (1.0 * wwx* wy*wwz) );
			grid->set( ix,  iy, iiz, grid->get( ix,  iy, iiz) + (1.0 * wwx*wwy* wz) );
			grid->set( ix, iiy, iiz, grid->get( ix, iiy, iiz) + (1.0 * wwx* wy* wz) );
			
			grid->set(iix,  iy,  iz, grid->get(iix,  iy,  iz) + (1.0 *  wx*wwy*wwz) );
			grid->set(iix, iiy,  iz, grid->get(iix, iiy,  iz) + (1.0 *  wx* wy*wwz) );
			grid->set(iix,  iy, iiz, grid->get(iix,  iy, iiz) + (1.0 *  wx*wwy* wz) );
			grid->set(iix, iiy, iiz, grid->get(iix, iiy, iiz) + (1.0 *  wx* wy* wz) );
		}
	}
}

void stamp_noise( ScalarGrid& grid, NoiseMachine& noise, const Vector& center, const float radius, const float fade )
{
	// iterate over grid points
	grid->setDefVal(0.0);
	#pragma omp parallel for
	for(int i=0; i<grid->nx();i++)
	{
		for(int j=0; j<grid->ny();j++)
		{
			for(int k=0; k<grid->nz();k++)
			{
				Vector pos = grid->evalP(i,j,k);
				// calculate q
				float q = (1.0 - (pos - center).magnitude()/radius)/fade;

				// calculate F
				float F;
				if( fade < 0)	F = 1;
				else if(q > 1)	F = 1;
				else if(q < 0)	F = 0;
				else 		F = q;


				if( q != 0) 
				{
					float cval = grid->get(i,j,k);
					float val  = noise->eval(pos);
					grid->set(i,j,k,std::max(cval, val * F));
				}
			}
		}
	}
}

}
