#include "Matrix.h"
#include "ImplicitVectorShapes.h"

#include <iostream>
using namespace std;
using namespace lux;

// CONSTANT
ConstantVectorVolume::ConstantVectorVolume( const Vector& v ) { value = v; gradvalue = Matrix();}

const Vector ConstantVectorVolume::eval( const Vector& P ) const
{
	return value;
}

const Matrix ConstantVectorVolume::grad( const Vector& P ) const
{
	return gradvalue;
}

// ADD
AddVectorVolume::AddVectorVolume( const VectorField& v1, const VectorField& v2 ) :
      elem1(v1),
      elem2(v2)
    {}

const Vector AddVectorVolume::eval( const Vector& P ) const
{
	return  elem1->eval(P) + elem2->eval(P);
}

const Matrix AddVectorVolume::grad( const Vector& P ) const
{
	return  elem1->grad(P) + elem2->grad(P);
}

SubtractVectorVolume::SubtractVectorVolume( const VectorField& v1, const VectorField& v2 ) :
      elem1(v1),
      elem2(v2)
    {}

const Vector SubtractVectorVolume::eval( const Vector& P ) const
{
	return  elem1->eval(P) - elem2->eval(P);
}

const Matrix SubtractVectorVolume::grad( const Vector& P ) const
{
	return  elem1->grad(P) - elem2->grad(P);
}

IdentityVectorVolume::IdentityVectorVolume() :
	gradvalue(Matrix(1,0,0, 0,1,0, 0,0,1))
      {}

const Vector IdentityVectorVolume::eval( const Vector& P ) const { return P; }

const Matrix IdentityVectorVolume::grad( const Vector& P ) const { return gradvalue; }

// GRADIENT
GradientVectorVolume::GradientVectorVolume( Volume<float>* v, const float dx ) :
      elem(v),
      step(dx)
    {}

GradientVectorVolume::GradientVectorVolume( const ScalarField& v, const float dx ) :
      elem(v),
      step(dx)
    {}

const Vector GradientVectorVolume::eval( const Vector& P ) const { return elem->grad(P); }

// evaluate gradient by finite difference
// size of difference is controlled via the "dx" attribute
const Matrix GradientVectorVolume::grad( const Vector& P ) const 
{
   float dx = step;
   Vector exp = eval(P+Vector(dx,0,0));
   Vector eyp = eval(P+Vector(0,dx,0));
   Vector ezp = eval(P+Vector(0,0,dx));
   Vector exm = eval(P-Vector(dx,0,0));
   Vector eym = eval(P-Vector(0,dx,0));
   Vector ezm = eval(P-Vector(0,0,dx));
   Matrix result( exp-exm, eyp-eym, ezp-ezm );
   result /= 2.0*dx;
   return result.transpose();
}

// ITERATED NEAREST POINT
ImplicitPointVectorVolume::ImplicitPointVectorVolume( const ScalarField& v, const float dx, const int nb, const float thresh) :
	elem(v),
	step(dx),
	nbIterations(nb),
	threshold(thresh)
{}

const Vector ImplicitPointVectorVolume::eval( const Vector& P ) const
{
	Vector CPT = P;

	if(nbIterations == 0)
	{
		CPT -= elem->eval(CPT) * elem->grad(CPT);
	} else
	{
		for(int i=0; i<nbIterations; i++)
		{
			Vector g = elem->grad(CPT);
			CPT -= (elem->eval(CPT) * g) / (g * g);
		}
	}

	return CPT;
}

// ADVECT
AdvectVectorVolume::AdvectVectorVolume( const VectorField& v, const VectorField& u, const float delt ):
	elem(v),
	velocity(u),
	dt(delt)
      {}

const Vector AdvectVectorVolume::eval( const Vector& P ) const
{
	Vector X = P - (velocity->eval(P))*dt;
	return elem->eval(X);
}

const Matrix AdvectVectorVolume::grad( const Vector& P ) const 
{
	Vector X = P - (velocity->eval(P) ) * dt;
	Matrix M = unitMatrix() - ( velocity->grad(P) )*dt;
	return M * (elem->grad(X));
}

// WARP
WarpVectorVolume::WarpVectorVolume( Volume<Vector>* v, Volume<Vector>* u ) :
	elem(v),
	warp(u)
      {}

WarpVectorVolume::WarpVectorVolume( const VectorField& v, const VectorField& u ):
	elem(v),
	warp(u)
      {} 
const Vector WarpVectorVolume::eval( const Vector& P ) const
{
	Vector X = warp->eval(P);
	return elem->eval(X);
}

// NOISE
NoiseVectorVolume::NoiseVectorVolume( Noise* n, const float d ) :
	noise(n),
	dx(d)
{}

NoiseVectorVolume::NoiseVectorVolume( NoiseMachine n, const float d ) :
	noise(n),
	dx(d)
{}

const Vector NoiseVectorVolume::eval( const Vector& P ) const
{
	return Vector(noise->eval(P), noise->eval(P + Vector(dx,0,0)), noise->eval(P - Vector(dx,0,0)) );
}

// GRID
GriddedSGridVectorVolume::GriddedSGridVectorVolume( const VectorGrid& g ):
	vegrid(g),
	dx(vegrid->dx()),
	dy(vegrid->dy()),
	dz(vegrid->dz())
      {}

const Vector GriddedSGridVectorVolume::eval( const Vector& P ) const
{
	return vegrid->eval(P);
}

