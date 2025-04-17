#include "RayMarcher.h"
#include "ProgressMeter.h"
#include "Vector.h"
#include "Fields.h"

#include <cmath>
#include <limits>
#include <random>

namespace lux
{

void RenderFrame(const RenderData* d, ProgressMeter& pm, float* image, std::function<void(const RenderData*, const Vector&, Color&)> Render )
{

	//camera.setEyeViewUp(eye, -eye, up);

	for(int j=0;j<d->Ny();j++)
	{
		double y = (double)j/double(d->Ny());
	#pragma omp parallel for
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

void RenderFrame(const RenderData* d, ProgressMeter& pm, float* image)
{

	
	for(int j=0;j<d->Ny();j++)
	{

		pm.update();
		double y = (double)j/double(d->Ny());
//		#pragma omp parallel for
		for(int i=0;i<d->Nx();i++)
		{
			double x = (double)i/double(d->Nx());
			Vector direction = d->camera.view(x,y);
			
			// check ray against bounding boxes, if there are any
			bool hit = d->boundingBoxes.empty();
			Ray r(d->camera.eye(), direction);
			float near = d->snear;
			float far = d->sfar;
			for(int box=0; box<d->boundingBoxes.size(); box++)
			{
				if(d->boundingBoxes[box]->intersection(r, near, far))
				{
					hit = true;
					break;
				}
			}
			// if the ray misses the bounding boxes, skip to next pixel
			if(!hit) continue;

			// if ray hits, march it
			Color L(0,0,0,0);
			if(d->dsmField.size() == 0)
			{
				RayMarchEmission(d, direction, L);
			}
			else
				RayMarchDSM(d, direction, L, near, far);

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

void RayMarchDSM(const RenderData* d, const Vector& direction, Color& L, 
		 const float near, const float far)
{
	// Optimizations from "A Note on Computation of a Ray Bending Path", Section 4,
	// by Jerry Tessendorf
	// 1.
	int steps = 0;
	float ds = d->ds;
	std::uniform_real_distribution<> ran(0.0, 1.0);
	std::mt19937 gen(0);

	double T = 1;
	// Calculate snear and sfar ...
	double s = near;
	Vector X = d->camera.eye() + direction*s;
	Vector new_pos;
	while( s < far && T > d->Tmin )
	{

/*

		Vector grad = d->densityField->grad(X) / d->densityField->eval(X);// 3. (b)
		float epsilon = d->Fmax;
		do 
		{
			new_pos = X + direction * ds;// 3. (c)
			steps++;// 3. (d)
			// 3. (e)
			if(steps > d->maxSteps)
			{
				return;
			}

			// 3. (f)
			Vector new_grad = d->densityField->grad(new_pos) / d->densityField->eval(new_pos);
			// 3. (g)
			epsilon = (new_grad - grad).magnitude() / grad.magnitude();

			if(epsilon >= d->Fmax) ds *= ran(gen);
//		std::cout << "s = " << s << " || ds = " << ds << std::endl;

		} while(epsilon >= d->Fmax);
		
		X += direction * ds;
		s += ds;
		
		if(epsilon < d->Fmin) ds /= ran(gen);

// block off here
*/ 
	
		float den = d->densityField->eval(X);
		if( den > 0.0 )
		{
			float dT = std::exp( -ds * d->kappa * den );
			// accumulate color
			Color Clights = Color();
			for(int a=0; a<d->lightColor.size(); a++)
			{
				Clights += d->lightColor[a] * d->dsmField[a]->eval(X);
			}

			L += d->colorField->eval(X) * Clights * (1-dT) * T / d->kappa;
			T *= dT;
		}
//		X += direction * ds;
//		s += ds;
	}
	L[3] = 1-T; // set the alpha channel to the opacity
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
	L[3] = 1-T; // set the alpha channel to the opacity
}



ScalarField RayMarchDSMAccumulation( 	const RenderData* d, 
					const ScalarField& densityField, 
					const Vector& lightPosition, 
					float ds,
					ScalarGrid& dsmField)
{

	int  nb = dsmField->nx() * dsmField->ny();
	ProgressMeter progress(nb, "Deep Shadow Map for light: " + string(lightPosition.__str__()));

	for(int i=0; i<dsmField->nx(); i++)
	{
		for(int j=0; j<dsmField->ny(); j++)
		{
			progress.update();
			#pragma omp parallel for
			for(int k=0; k<dsmField->nz(); k++)
			{
				double arg = 0.0;
				Vector X = dsmField->evalP(i,j,k);
				if( densityField->eval(X) > 0.0)
				{
					double smax = (lightPosition-X).magnitude();
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
