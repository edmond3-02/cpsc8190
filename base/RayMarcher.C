#include "RayMarcher.h"
#include <cmath>

namespace lux
{

void RenderFrame(const RenderData* d, float* image)
{

	//camera.setEyeViewUp(eye, -eye, up);

	for(int j=0;j<d->Ny();j++)
	{
		//progress.update();
		double y = (double)j/double(d->Ny());
	//#pragma omp parallel for
		for(int i=0;i<d->Nx();i++)
		{

			double x = (double)i/double(d->Nx());
			Vector direction = d->camera.view(x,y);

			Color L = RayMarchEmission(d, direction);

			for(int c=0; c<d->Nc(); c++)
			{
				image[c + d->Nc()*(i + d->Nx()*j)] = L[c];
			}
		}
	}
	

}

Color RayMarchEmission(const RenderData* d, const Vector& direction)
{
	
	double T = 1;
	Color L(0,0,0,0);
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
	return L;
}


}
