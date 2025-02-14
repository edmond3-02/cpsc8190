#include <cmath>
#include "ImplicitVolumeShapes.h"
#include <iostream>
using namespace std;

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

CsgBoxVolume::CsgBoxVolume( const Vector& cen, const float rad, float pwr) :
       center       (cen),
       radius       (rad),
       power	    (pwr)
    {}


const float CsgBoxVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;

   float f = std::pow(radius,power) - std::pow(x.X(),power) - std::pow(x.Y(),power) - std::pow(x.Z(),power);
   

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

MultiplyScalarVolume::MultiplyScalarVolume( const ScalarField&  e, const float v ) :
      e1 (e),
      val (v)
    {}

const float MultiplyScalarVolume::eval( const Vector& P ) const
{
   return  e1->eval(P) * val;
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

const float UnionVolume::eval( const Vector& P ) const
{
   return  max(e1->eval(P), e2->eval(P));
}

IntersectionVolume::IntersectionVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float IntersectionVolume::eval( const Vector& P ) const
{
   return  min(e1->eval(P), e2->eval(P));
}

CutoutVolume::CutoutVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float CutoutVolume::eval( const Vector& P ) const
{
   return  min(e1->eval(P), -e2->eval(P));
}

ShellVolume::ShellVolume( const ScalarField&  v1, const float v2 ) :
      e (v1),
      d (v2)
    {}

const float ShellVolume::eval( const Vector& P ) const
{
   return  min( (e->eval(P) + d/2.0), -(e->eval(P) - d/2.0) );
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
