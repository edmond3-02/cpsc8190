#include <cmath>
#include "ImplicitVolumeShapes.h"
#include <iostream>

using namespace lux;


ConstantVolume::ConstantVolume( const float v ) :
       value (v),
       gradvalue (Vector(0,0,0))
    {}


const float ConstantVolume::eval( const Vector& P ) const 
{
   return value; 
}

const Vector ConstantVolume::grad( const Vector& P ) const { return gradvalue; }


ScaleVolume::ScaleVolume(const ScalarField& field, const Vector val ) :
       e(field),
       v(val)
    {}


const float ScaleVolume::eval( const Vector& P ) const 
{
   Vector x = P;
   for(int i=0; i<3; i++) x[i] /= v[i];
   return e->eval(x);
}

TranslateVolume::TranslateVolume(const ScalarField& field, const Vector value ) :
       e(field),
       v(value)
    {}


const float TranslateVolume::eval( const Vector& P ) const 
{
   return e->eval(P - v);
}

RotateVolume::RotateVolume(const ScalarField& field, const float angle, const Vector axis ) :
       e(field),
       a(angle),
       a0(axis)
    {
	a0.normalize();
    }

const float RotateVolume::eval( const Vector& P ) const 
{
   Vector x = P;
   x = x * cos(a) + a0*(a0*x)*(1-cos(a)) + (x^a0)*sin(a);
   return e->eval(x);
}

NegateVolume::NegateVolume(const ScalarField& v ) :
       elem(v)
    {}

const float NegateVolume::eval( const Vector& P ) const
{
   return -elem->eval(P);
}

ExpVolume::ExpVolume( Volume<float>* v ) :
      elem (v)
    {}

ExpVolume::ExpVolume( const ScalarField& v ) :
      elem (v)
    {}

const float ExpVolume::eval( const Vector& P ) const 
{
   return std::exp( elem->eval(P) ); 
}

const Vector ExpVolume::grad( const Vector& P ) const { return eval(P) * elem->grad(P); }




SphereVolume::SphereVolume( const Vector& cen, const float rad ) :
       center (cen),
       radius (rad)
    {}


const float SphereVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   float rad = x.magnitude();
   return radius-rad;
}

const Vector SphereVolume::grad( const Vector& P ) const 
{
   Vector x = P-center;
   if( x.magnitude() != 0 ){ return -x.unitvector(); }
   return Vector(0,1,0);
}

EllipseVolume::EllipseVolume( const Vector& cen, const Vector& axs, const float majorrad, const float minorrad ) :
       center       (cen),
       axis       (axs),
       major_radius (majorrad),
       minor_radius (minorrad)
    {
	axis.normalize();
    }


const float EllipseVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;

   float f = 0.0;
   float Z = x*axis;
   Vector xt = x - Z*axis;
   f = 1 - std::pow(Z,2)/std::pow(major_radius,2) - std::pow(xt.magnitude(),2)/std::pow(minor_radius,2);
   return f;
}

AdvectVolume::AdvectVolume( const ScalarField& v, const VectorField& u, const float delt ):
	elem(v),
	velocity(u),
	dt(delt)
      {}

const float AdvectVolume::eval( const Vector& P ) const
{
	Vector X = P - ( velocity->eval(P) ) * dt;
	return elem->eval(X);
}


CsgBoxVolume::CsgBoxVolume( const Vector& cen, const float rad, float pwr) :
       center       (cen),
       radius       (rad),
       power	    (pwr)
    {}


const float CsgBoxVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;

   float f = radius - std::pow(x.X(),power) - std::pow(x.Y(),power) - std::pow(x.Z(),power);
   

   return f;
}


ConeVolume::ConeVolume( const Vector& cen, const Vector& axs, const float h, const float theta ) :
       center       (cen),
       axis       (axs),
       height     (h),
       angle      (theta)
    {}


const float ConeVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;

   float f = 0.0;
   float xdotn = x*axis;
   
   if(xdotn<0) f = xdotn;
   else if(xdotn>height) f = height - xdotn;
   else f = xdotn - x.magnitude()*cos(angle);

   return f;
}

PlaneVolume::PlaneVolume( const Vector& cen, const Vector& norm ) :
       center     (cen),
       normal	  (norm)
    {
	normal.normalize(); 
    }


const float PlaneVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   float f = -x * normal;

   return f;
}

TorusVolume::TorusVolume( const Vector& cen, const Vector& axs, const float majRad, const float minRad ) :
       center     (cen),
       axis	  (axs),
       major_radius(majRad),
       minor_radius(minRad)
    {
	axis.normalize(); 
    }


const float TorusVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   float f = 0.0;

   Vector xt = x - (x*axis)*axis;

   f = 4*std::pow(major_radius,2)*std::pow(xt.magnitude(),2) - std::pow(std::pow(x.magnitude(),2)+std::pow(major_radius,2)-std::pow(minor_radius,2),2);

   return f;
}

SteinerPatchVolume::SteinerPatchVolume( )   {}


const float SteinerPatchVolume::eval( const Vector& P ) const 
{
   Vector x = P;

   float x2 = std::pow(x.X(),2),
         y2 = std::pow(x.Y(),2),
         z2 = std::pow(x.Z(),2);

   float f = -(x2*y2 + x2*z2 + y2*z2 - x.X()*x.Y()*x.Z());

   return f;
}

IcosahedronVolume::IcosahedronVolume( )   {}


const float IcosahedronVolume::eval( const Vector& P ) const 
{
   Vector x = P;

   float f  = 0.0;
   float pi = 3.14159265359;
   float T  = 1.61803399;

   if(x.magnitude() > 1.8*pi)
	f = -1.8*pi;
   else if(x.magnitude() <= 1.8*pi)
	f = cos(x.X() + T*x.Y()) + cos(x.X() - T*x.Y())
	  + cos(x.Y() + T*x.Z()) + cos(x.Y() - T*x.Z())
	  + cos(x.Z() + T*x.X()) + cos(x.Z() - T*x.X())
	  - 2.0;

   return f;
}

AbsoluteVolume::AbsoluteVolume( const ScalarField& v ) :
     elem(v)
   {}

AbsoluteVolume::AbsoluteVolume( Volume<float> * v ) :
     elem(v)
   {}

const float AbsoluteVolume::eval( const Vector& P ) const
{
	return std::abs(elem->eval(P));
}

MultiplyVolume::MultiplyVolume( const ScalarField&  e, const float v ) :
      elem (e),
      factor(nullptr),
      constant (v)
    {}

MultiplyVolume::MultiplyVolume( const ScalarField&  e, const ScalarField& u) :
      elem (e),
      factor(u),
      constant (0)
    {}

const float MultiplyVolume::eval( const Vector& P ) const
{
   if(factor == nullptr)
	return  elem->eval(P) * constant;
   else
	return elem->eval(P) * factor->eval(P);
}

AddVolume::AddVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float AddVolume::eval( const Vector& P ) const
{
   return  e1->eval(P) + e2->eval(P);
}

const Vector AddVolume::grad( const Vector& P ) const
{
   return  e1->grad(P) + e2->grad(P);
}

MaskVolume::MaskVolume( const ScalarField&  v1 ) :
      e1 (v1)
    {}

const float MaskVolume::eval( const Vector& P ) const
{
   return  (e1->eval(P) > 0) ? 1 : 0;
}

UnionVolume::UnionVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

ClampVolume::ClampVolume( const ScalarField& v, float minv, float maxv ) :
      elem(v),
      vmin(minv),
      vmax(maxv)
    {}

const float ClampVolume::eval( const Vector& P ) const
{
	return std::min(std::max(elem->eval(P), vmin),vmax);
}
 
GammaVolume::GammaVolume( const ScalarField& v, float gam ) :
       elem(v),
       factor(nullptr),
       gamma(gam)
     {}
  
GammaVolume::GammaVolume( const ScalarField& v, const ScalarField& g ) :
       elem(v),
       factor(g),
       gamma(0)
     {}

const float GammaVolume::eval( const Vector& P ) const
{
	return (factor == nullptr) ? std::pow(elem->eval(P), gamma) : std::pow(elem->eval(P), factor->eval(P));
}

BlinnBlendVolume::BlinnBlendVolume( Volume<float> * v1, Volume<float> * v2, const float _alpha) :
        elem1(v1),
        elem2(v2),
        alpha(_alpha)
      {}
  
 BlinnBlendVolume::BlinnBlendVolume( const ScalarField& v1, const ScalarField& v2, const float _alpha) :
        elem1(v1),
        elem2(v2),
        alpha(_alpha)
      {}
 
const float BlinnBlendVolume::eval( const Vector& P ) const
{
	return std::exp(alpha * elem1->eval(P)) + std::exp(alpha * elem2->eval(P));
}
 
const float UnionVolume::eval( const Vector& P ) const
{
   return  std::max(e1->eval(P), e2->eval(P));
}

IntersectionVolume::IntersectionVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float IntersectionVolume::eval( const Vector& P ) const
{
   return  std::min(e1->eval(P), e2->eval(P));
}

CutoutVolume::CutoutVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float CutoutVolume::eval( const Vector& P ) const
{
   return  std::min(e1->eval(P), -e2->eval(P));
}

PyroclasticVolume::PyroclasticVolume( const NoiseMachine n, const ScalarField e, const float Amp, const int i, const float Gamma ) :
						noise(n),
						elem(e),
						amplitude(Amp),
						iter(i),
						gamma(Gamma)
	{}

const float PyroclasticVolume::eval(const Vector& P) const
{
	Vector Xc = P;
	if(iter == 0)
	{
		Xc -= elem->eval(Xc) * elem->grad(Xc);
	} else
	{
	// works?
	/////////////////////////////////////////////
		Vector g = elem->grad(Xc);
		for(int i=0; i<iter; i++)
		{
			if(g*g == 0.0) break;
			g = elem->grad(Xc);
			Xc -= (elem->eval(Xc) * g) / (g * g);
		}
	// Doesn't. why?
	//////////////////////////////////////////////	
	/*	for(int i=0; i<iter; i++)
		{
			Vector g = elem->grad(Xc);
			if(g*g == 0.0) break;
			Xc -= (elem->eval(Xc) * g) / (g * g);
		}*/
	}

	return elem->eval(P) + std::pow(std::abs(noise->eval(Xc)), gamma) * amplitude;

}
				      		

RadialPyroclasticVolume::RadialPyroclasticVolume( const Vector& Center, const float Radius, const float Amp, 
						  const float octaves,  const float freq,   const float rough, 
						  const float trans,    const float time,   const float Gamma ) :
						center   (Center),
						amplitude(Amp),
						gamma    (Gamma),
						translate(trans),
						radius   (Radius)
	{
		parms.octaves   = octaves;
		parms.frequency = freq;
		parms.roughness = rough;
		parms.time      = time;
		noise.setParameters(parms);
	}

const float RadialPyroclasticVolume::eval( const Vector& P ) const
{
	Vector X = P - center;
	Vector Xc = radius * (X/X.magnitude());

	return (radius - X.magnitude()) + std::pow(std::abs(noise.eval(Xc)), gamma) * amplitude;

}

ShellVolume::ShellVolume( const ScalarField&  v1, const float v2 ) :
      e (v1),
      d (v2)
    {}

const float ShellVolume::eval( const Vector& P ) const
{
   return  std::min( (e->eval(P) + d/2.0), -(e->eval(P) - d/2.0) );
}

// NOISE
NoiseVolume::NoiseVolume( Noise* n, const float d) :
	noise(n),
	dx   (d)
     {}


NoiseVolume::NoiseVolume( NoiseMachine n, const float d) :
	noise(n),
	dx   (d)
     {}

const float NoiseVolume::eval( const Vector& P ) const
{
	return noise->eval(P);	
} 

GriddedSGridVolume::GriddedSGridVolume( const ScalarGrid& g ):
	scgrid(g),
	dx(scgrid->dx()),
	dy(scgrid->dy()),
	dz(scgrid->dz())
     {}

const float GriddedSGridVolume::eval( const Vector& P ) const
{
	return scgrid->eval(P);
}


WarpVolume::WarpVolume( const ScalarField& v, VectorField& map ):
	elem(v),
	mapX(map)
      {}


const float WarpVolume::eval( const Vector& P ) const
{
	Vector X = mapX->eval(P);
	return elem->eval(X);
}

const Vector WarpVolume::grad(  const Vector& P ) const
{
	Vector X = mapX->eval(P);
	Matrix gX = mapX->grad(P);
	return gX * elem->grad(X);
}
