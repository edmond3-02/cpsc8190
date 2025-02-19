#ifndef __GEOMETRYVOLUMESHAPES_H__
#define __GEOMETRYVOLUMESHAPES_H__
 
#include <cmath>
#include "Volume.h"
#include "LinearAlgebra.h"
//#include "Noise.h"
//#include "PerlinNoise.h"
 
namespace lux
{
 
 
class PlaneLevelSet : public Volume<float> 
{
  public:
  
    PlaneLevelSet( const Vector& p0, const Vector& N ) :
       P0     (p0),
       normal (N.unitvector())
    {}
 
   ~PlaneLevelSet(){}
 
    const float eval( const Vector& P ) const
    {
       return (P-P0)*normal;
    }
 
    const Vector grad( const Vector& P ) const
    {
       return normal;
    }
 
 
  protected:
 
    Vector P0;
    Vector normal;
};
 
 
 
 
class Triangle
{
  public:
 
    Triangle( const Vector& p0, const Vector& p1, const Vector& p2 ) :
       P0  (p0)
    {
       du = p1-p0;
       dv = p2-p0;
 
 
       double denom = (du*du) * (dv*dv) - (du*dv) * (du*dv);
       if( denom != 0 )
       {
          good = true;
          U = (du * (dv*dv) - dv * (du*dv))/denom;
          V = (dv * (du*du) - du * (du*dv))/denom;
          normal = du^dv/sqrt(denom);
       }
       else
       {
          good = false;
          normal = du^dv;
       }
    }
 
   ~Triangle(){}
 
    const Vector barycentric( const Vector& P ) const
    {
       Vector dP = P-P0;
       return Vector( dP*U, dP*V, dP*normal );
    }
 
    const bool Intersected( const Vector& P, const Vector& D ) const;
 
    const Vector& dU() const { return U; }
    const Vector& dV() const { return V; }
    const Vector& N()  const { return normal; }
    const Vector& P()  const { return P0; }
    const Vector& Eu()  const { return du; }
    const Vector& Ev()  const { return dv; }
 
    const bool isGood() const { return good; }
 
 
  private:
 
    Vector P0, U, V, normal;
    Vector du, dv;
    bool good;
 
 
};
 
 
class TexturedTriangle : public Triangle
{
  public:
 
    TexturedTriangle( const Vector& p0, const Vector& p1, const Vector& p2, 
                      const Vector& t0, const Vector& t1, const Vector& t2  ) :
    Triangle(p0,p1,p2),
    T0      (t0),
    dT1      (t1-t0),
    dT2      (t2-t0)
     {}
 
    ~TexturedTriangle(){}
 
    const Vector barycentricTexture( const Vector& P ) const
    {
       Vector U = barycentric(P);
       return T0 + U[0]*dT1 + U[1]*dT2;
    }
 
    const Vector& T() const { return T0; }
    const Vector& Tu() const { return dT1; }
    const Vector& Tv() const { return dT2; }
 
  protected:
 
    Vector T0,dT1,dT2;
};
 
 
 
//const Vector Intersection( const Triangle& t, const Vector& P, const Vector& D );
//const double ClosestDistance( const Triangle& t, const Vector& P );
//const Vector ClosestDistance( const TexturedTriangle& t, const Vector& P );
 
 
class TriangleLevelSet : public Volume<float> 
{
  public:
  
    TriangleLevelSet( const Vector& p0, const Vector& p1, const Vector& p2, const float bb );
    TriangleLevelSet( const Vector& p0, const Vector& p1, const Vector& p2 );
 
   ~TriangleLevelSet(){}
 
    const float eval( const Vector& P ) const;
    const Vector grad( const Vector& P ) const;
 
    const Vector& LLC() const { return llc; }
    const Vector& URC() const { return urc; }
 
    const bool isRejected() const { return reject; }
 
  protected:
 
    Triangle T;
    //Vector P0;
    //Vector edge10, edge20;
    //Vector normal;
 
    //Vector U,V;
 
    Vector llc, urc;
    bool useBB;
    bool reject;
 
};
 
 
 
 
}
#endif
