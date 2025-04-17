#include "AARectangle.h"
#include "LinearAlgebra.h"

using namespace lux;

AARectangle::AARectangle():
	_llc(Vector(-.5,-.5,-.5)),
	_urc(Vector(.5,.5,.5)),
	_length(Vector(1,1,1)),
	_center(Vector(0,0,0))
      {
		bounds[0] = _llc;
		bounds[1] = _urc;
      }

AARectangle::AARectangle(const Vector& llc, const Vector& urc):
	_llc(llc),
	_urc(urc),
	_length(urc - llc),
	_center((urc+llc)/2)
      {
		bounds[0] = _llc;
		bounds[1] = _urc;
      }

const bool AARectangle::isInside( const Vector& P ) const
{
	if(P.X() < _llc.X()  ||  P.X() > _urc.X()) return false;
	if(P.Y() < _llc.Y()  ||  P.Y() > _urc.Y()) return false;
	if(P.Z() < _llc.Z()  ||  P.Z() > _urc.Z()) return false;

	return true;
}


const double AARectangle::signedDistance( const Vector& P ) const
{
	Vector d = max(_llc - P, P - _urc);
	float length = std::max(0.0, d.magnitude());
	length += std::min(0.0, std::max(std::max(d.X(), d.Y()), d.Z()) ); 
	return -length;
}
/*
const bool AARectangle::intersection( const Ray r, double& close, double& far ) const
{
	float tmin, tmax, tymin, tymax, tzmin, tzmax;
	tmin = (bounds[r.sign[0]].x() - r.origin.x()) * r.inv_direction.x();
	tmax = (bounds[1-r.sign[0]].x() - r.origin.x()) * r.inv_direction.x();
	tymin = (bounds[r.sign[1]].y() - r.origin.y()) * r.inv_direction.y();
	tymax = (bounds[1-r.sign[1]].y() - r.origin.y()) * r.inv_direction.y();
	if ( (tmin > tymax) || (tymin > tmax) )
		return false;
	if (tymin > tmin)
		 tmin = tymin;
	if (tymax < tmax)
		 tmax = tymax;
	tzmin = (bounds[r.sign[2]].z() - r.origin.z()) * r.inv_direction.z();
	tzmax = (bounds[1-r.sign[2]].z() - r.origin.z()) * r.inv_direction.z();
	if ( (tmin > tzmax) || (tzmin > tmax) )
		return false;
	if (tzmin > tmin)
		 tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;
	return ( (tmin < t1) && (tmax > t0) );

}
*/

const bool AARectangle::intersection( const Ray r, float& close, float& far ) const
{
	float t0=close, t1=far;
	float tmin, tmax, tymin, tymax, tzmin, tzmax;
	tmin = (bounds[r.sign(0)].X() - r.origin().X()) * r.inv_direction().X();
	tmax = (bounds[1-r.sign(0)].X() - r.origin().X()) * r.inv_direction().X();
	tymin = (bounds[r.sign(1)].Y() - r.origin().Y()) * r.inv_direction().Y();
	tymax = (bounds[1-r.sign(1)].Y() - r.origin().Y()) * r.inv_direction().Y();
	if ( (tmin > tymax) || (tymin > tmax) )
		return false;
	if (tymin > tmin)
		 tmin = tymin;
	if (tymax < tmax)
		 tmax = tymax;
	tzmin = (bounds[r.sign(2)].Z() - r.origin().Z()) * r.inv_direction().Z();
	tzmax = (bounds[1-r.sign(2)].Z() - r.origin().Z()) * r.inv_direction().Z();
	if ( (tmin > tzmax) || (tzmin > tmax) )
		return false;
	if (tzmin > tmin)
		 tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;

	close=tmin;
	far=tmax;
	return ( (tmin < tmax) && (tmax > t0) );

}


AABB::AABB( const Vector& llc, const Vector& urc ): 
	std::shared_ptr<AARectangle>( new AARectangle(llc, urc )) {}

