#include <cstdio>

#include "Stamp.h"
#include "Volume.h"
#include "ProgressMeter.h"


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


}
