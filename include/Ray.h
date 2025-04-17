#ifndef RAY_H
#define RAY_H

#include "Vector.h"

namespace lux{

class Ray 
{
	public:
		Ray(const Vector &o, const Vector &d);
		const Vector origin() const; 
		const Vector direction() const;
		const Vector inv_direction() const;
		const int sign(int i) const;


	private:	
		Vector _origin;
		Vector _direction;
		Vector _inv_direction;
		int _sign[3];
};

}

#endif
