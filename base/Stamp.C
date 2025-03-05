#include <cstdio>

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
	std::cout << "dim: " << grid->nx() << " " << grid->ny() << " " << grid->nz() << std::endl;
	for(int i=0; i<grid->nx();i++)
	{
		for(int j=0; j<grid->ny();j++)
		{
			for(int k=0; k<grid->nz();k++)
			{;
				Vector pos = grid->evalP(i,j,k);
				float val  = field->eval(pos);
				//std::cout<<val << " ";
				grid->set(i,j,k,val);
			}
			//std::cout<<std::endl;
		}
		//std::cout << std::endl;
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
		NoiseMachine noise = perlin(particles[i]);
		stamp_noise(grid, noise, particles[i].P, particles[i].radius, particles[i].falloff);
	}
}

void stamp_noise( ScalarGrid& grid, NoiseMachine& noise, const Vector& center, const float radius, const float fade )
{
	// iterate over grid points
	grid->setDefVal(0.0);
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
					float val  = noise->eval(pos);
					grid->set(i,j,k,val * F);
				}
			}
		}
	}
}

}
