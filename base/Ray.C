#include "Ray.h"

#include "Vector.h"

using namespace lux;


Ray::Ray(const Vector &o, const Vector &d) 
{
	_origin = o;
	_direction = d;
	_inv_direction = Vector(1/d.X(), 1/d.Y(), 1/d.Z());
	_sign[0] = (_inv_direction.X() < 0);
	_sign[1] = (_inv_direction.Y() < 0);
	_sign[2] = (_inv_direction.Z() < 0);
}

const Vector Ray::origin() const { return _origin; }
const Vector Ray::direction() const { return _direction; }
const Vector Ray::inv_direction() const { return _inv_direction; }
const int    Ray::sign(int i) const { return _sign[i]; }
