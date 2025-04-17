#ifndef __AARECTANGLE_H__
#define __AARECTANGLE_H__
 
#include "Vector.h"
#include "Ray.h"
//#include "RectangularGrid.h"
#include <memory>
 
namespace lux
{
 
class AARectangle
{
  public:
 
    AARectangle(); 
    AARectangle( const Vector& llc, const Vector& urc );
   ~AARectangle(){}
 
    const bool isInside( const Vector& P ) const;
 
    const double signedDistance( const Vector& P ) const;
 
    // returns true and the intersection point if there is an intersection
    // From: An Efficient and Robust Ray–Box Intersection Algorithm
    //       by Amy Williams Steve Barrus R. Keith Morley Peter Shirley
    const bool intersection( const Ray r, float& close, float& far ) const;
 
//    const double farIntersection( const Vector& P, const Vector& D ) const;
 
//    const double nearIntersection( const Vector& P, const Vector& D ) const;
 
//    const Vector normal( const Vector& P ) const;
 
    const Vector& llc() const { return _llc; }
    const Vector& urc() const { return _urc; }
    const Vector& length() const { return _length; }
    const Vector& center() const { return _center; }
 
//    void split( const int component, AARectangle& aabb1, AARectangle& aabb2 ) const;
 
//    const bool intersects( const AARectangle& aabb ) const;
 
//    AARectangle Union( const AARectangle& aabb ) const;
 
  private:
 
    Vector _llc, _urc;
    Vector _length, _center;
    Vector bounds[2];
 
};
 
 
typedef std::shared_ptr<AARectangle> AABBBase;
 
class AABB : public AABBBase
{
  public:
 
    AABB( const Vector& llc, const Vector& urc );
    //AABB( const AABB& aabb );
   //~AABB();
}; 

//AABB makeAABB( const Vector& llc, const Vector& urc );
//AABB expandAABB( const AABB& aabb, const double l );
//AABB shrinkAABB( const AABB& aabb, const double l );
//AABB expandAABB( const AABB& aabb, const Vector& l );
//AABB shrinkAABB( const AABB& aabb, const Vector& l );
//AABB translateAABB( const AABB& aabb, const Vector& t );
// 
//AABB Union( const AABB& aabb1, const AABB& aabb2 );
// 
//const bool isInside( const AABB& aabb, const Vector& P );
//const Vector& getAABBLLC( const AABB& aabb );
//const Vector& getAABBURC( const AABB& aabb );
//const Vector& getAABBCenter( const AABB& aabb );
//const Vector& getAABBLength( const AABB& aabb );
//const bool intersectingAABB(const AABB& aabb1, const AABB& aabb2 );
 

}
 
#endif
