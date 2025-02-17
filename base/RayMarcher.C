#include "RayMarcher.h"
#include <cmath>
#include "ProgressMeter.h"
#include "Vector.h"
#include "Fields.h"

namespace lux
{

void RenderFrame(const RenderData* d, ProgressMeter& pm, float* image, std::function<void(const RenderData*, const Vector&, Color&)> Render )
{

	//camera.setEyeViewUp(eye, -eye, up);

	for(int j=0;j<d->Ny();j++)
	{
		//progress.update();
		double y = (double)j/double(d->Ny());
	//#pragma omp parallel for
		for(int i=0;i<d->Nx();i++)
		{
			pm.update();

			double x = (double)i/double(d->Nx());
			Vector direction = d->camera.view(x,y);

			Color L(0,0,0,0);
			Render(d, direction, L);

			for(int c=0; c<d->Nc(); c++)
			{
				image[c + d->Nc()*(i + d->Nx()*j)] = L[c];
			}
		}
	}
	

}

void RayMarchEmission(const RenderData* d, const Vector& direction, Color& L)
{
	
	double T = 1;
	// Calculate snear and sfar ...
	double snear = d->snear;
	double sfar = d->sfar;
	double s = snear;
	Vector X = d->camera.eye() + direction*s;
	while( s < sfar && T > d->Tmin )
	{
		float den = d->densityField->eval(X);
		if( den > 0.0 )
		{
			float dT = std::exp( -d->ds * d->kappa * den );
			L += d->colorField->eval(X) * (1-dT) * T / d->kappa;
			T *= dT;
		}
		X += direction * d->ds;
		s += d->ds;
	}
	L[3] = 1;//-T; // set the alpha channel to the opacity
}

void RayMarchDSM(const RenderData* d, const Vector& direction, Color& L)
{
	double T = 1;
	// Calculate snear and sfar ...
	double snear = d->snear;
	double sfar = d->sfar;
	double s = snear;
	Vector X = d->camera.eye() + direction*s;
	while( s < sfar && T > d->Tmin )
	{
		float den = d->densityField->eval(X);
		if( den > 0.0 )
		{
			float dT = std::exp( -d->ds * d->kappa * den );
			// accumulate color
			Color Clights = Color();
			for(int a=0; a<d->lightColor.size(); a++)
			{
				Clights += d->lightColor[a] * d->dsmField[a]->eval(X);
			}

			L += d->colorField->eval(X) * Clights * (1-dT) * T / d->kappa;
			T *= dT;
		}
		X += direction * d->ds;
		s += d->ds;
	}
	L[3] = 1;//-T; // set the alpha channel to the opacity
}
	



ScalarField RayMarchDSMAccumulation( 	const RenderData* d, 
					const ScalarField& densityField, 
					const Vector& lightPosition, 
					float ds,
					ScalarGrid& dsmField)
{

	for(int i=0; i<dsmField->nx(); i++)
	{
		for(int j=0; j<dsmField->ny(); j++)
		{
			for(int k=0; k<dsmField->nz(); k++)
			{
				double arg = 0.0;
				Vector X = dsmField->evalP(i,j,k);

				if( densityField->eval(X) > 0.0)
				{
					float smax = (lightPosition-X).magnitude();
					Vector direction = (lightPosition-X).unitvector();
					double s = 0;
					while(s < smax)
					{
						arg += densityField->eval(X)*d->ds;
						X += direction * d->ds;
						s += d->ds;
					}
					dsmField->set(i,j,k,arg);
				}
			}
		}
	}

	return exp(volumize(dsmField) * -d->kappa );
}
	 
}
