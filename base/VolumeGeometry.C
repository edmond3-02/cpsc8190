
#include <cmath>
#include <climits>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include "VolumeGeometry.h"
#include "ProgressMeter.h"

namespace lux
{

Mesh createEmptyMesh() { return Mesh(new TriangleGeometry() ); }

void RayMarchLevelSet( const Mesh& geom, ScalarGrid& lsgrid, const int bandwidth )
{
	ProgressMeter  calcDistProg(geom->nbFaces(), "Calculate distances");

	lsgrid->setDefVal(100000);
	// ITERATE TRIS
	for(int t=0; t<geom->nbFaces(); t++)
	{
		calcDistProg.update();

		int i,j,k;
		geom->getFace(t, i,j,k);

		Vector v1 = geom->getVertex(i);
		Vector v2 = geom->getVertex(j);
		Vector v3 = geom->getVertex(k);


		Vector llc = min(v1, min(v2, v3) );
		Vector urc = max(v1, max(v2, v3) );
	

		int li, lj, lk;
		int ui, uj, uk;
		lsgrid->getGridIndex(llc, li, lj, lk);
		
		lsgrid->getGridIndex(urc, ui, uj, uk);
		// plus one to go to 'outer' index of urc
		ui++;
		uj++;
		uk++;

		for(int ci=li-bandwidth; ci<ui+bandwidth; ci++)
		{
			for(int cj=lj-bandwidth; cj<uj+bandwidth; cj++)
			{
				for(int ck=lk-bandwidth; ck<uk+bandwidth; ck++)
				{
					if(!lsgrid->isInside(ci,cj,ck)) continue;

					float val = lsgrid->get(ci, cj, ck);

					// CALCULATE MIN DISTANCE
					Vector point = lsgrid->evalP(ci,cj,ck);
					Vector closestPoint = closestPointTriangle(point, v1, v2, v3);
					
					float dist = std::sqrt(
							 std::pow( std::abs(point.X() - closestPoint.X()) , 2)
							+std::pow( std::abs(point.Y() - closestPoint.Y()) , 2)
							+std::pow( std::abs(point.Z() - closestPoint.Z()) , 2) );

					dist = std::abs(dist);

					if(dist < val)
						lsgrid->set(ci, cj, ck, dist);
				}
			}
		}

	}

	
	//ProgressMeter  SetSignProg(lsgrid->ny(), "Set signs");

	// ITERATE GRID POINTS
	for(int y=0; y<lsgrid->ny(); y++)
	{
		//SetSignProg.update();
		//#pragma omp parallel for
		for(int x=0; x<lsgrid->nx(); x++)
		{

			vector<float> intersections = FindAllIntersections(geom, lsgrid->evalP(x,y,0), Vector(0,0,1).unitvector());
			int intercount = 0;
// DEBUG
//				std::cout << intersections.size()<<" ";

			for(int z=0; z<lsgrid->nz(); z++)
			{
				float d = z * std::abs(lsgrid->dz());
				
				while( intercount < intersections.size() && d > intersections[intercount]) 
					intercount++;

				float val = lsgrid->get(x, y, z);
   				// removing the iff statement fixes it. WHY? and how to avoid?
				// WHY : because grid points within allocated blocks are initialized with positive values, and arent picked up by the if
				// HOW TO AVOID : use lsgrid->goodBlock?
				//if(val != lsgrid->getDefVal() || (intercount%2) == 1)   
				if(lsgrid->goodBlock(x,y,z) || (intercount%2) == 1)   
				lsgrid->set(x, y, z, val * std::pow(-1, 1+intercount) );


			}
		}

//		std::cout << std::endl;
	}

	lsgrid->setDefVal(-100000);

/*
	int x = lsgrid->nx()/2;
	for(int y=0; y<lsgrid->ny(); y++)
	{
		for(int z=0; z<lsgrid->nz(); z++)
		{
			if(z<10) continue;
			float val = lsgrid->get(x, y, z);
			if(val == std::abs(lsgrid->getDefVal())) std::cout << "++++ ";
			else if(val == lsgrid->getDefVal()) std::cout << "---- ";
			else std::cout << std::fixed << std::setprecision(3) << val << " ";


		}
		std::cout << std::endl;
	}
*/
}


vector<float> FindAllIntersections( const Mesh& g, const Vector P, const Vector D)
{
	vector<float> vec;

	for(int c=0; c<g->nbFaces(); c++)
	{
		int i,j,k;
		g->getFace(c, i,j,k);

		Vector x0 = g->getVertex(i);
		Vector x1 = g->getVertex(j);
		Vector x2 = g->getVertex(k);

		Vector e1 = (x1-x0);
		Vector e2 = (x2-x0);

		Vector nt = e2 ^ e1;
		nt.normalize();

		float dist = (nt * (x0 - P)) / (nt * D);
		Vector xxn = (P-x0) + (dist *D);
		Vector cross = e2 ^ e1;

		float u =  ((e2 ^ xxn) * (e2^e1))/std::pow((e2^e1).magnitude(), 2);
		float v = -((e1 ^ xxn) * (e2^e1))/std::pow((e2^e1).magnitude(), 2);

		if(u >= 0 && u <= 1)
		{
			if(v >= 0 && v <= 1)
			{
				if(u+v <= 1)
				{
					vec.push_back(std::abs(dist));
				}
			}
		}
	}
	if(!vec.empty())
	{
		sort(vec.begin(), vec.end());
	}


	return vec;


}


Vector closestPointTriangle(const Vector& p, const Vector& v1,const Vector& v2, const Vector& v3)
{
  const Vector e1 = v2 - v1;
  const Vector e2 = v3 - v1;

  // Check if first vertex is closest
  const Vector p1 = p - v1;
  const float d1 = e1 * p1;
  const float d2 = e2 * p1;
  if (d1 <= 0 && d2 <= 0) return v1; 

  // Check if second vertex is closest
  const Vector p2 = p - v2;
  const float d3 = e1 * p2;
  const float d4 = e2 * p2;
  if (d3 >= 0.0 && d4 <= d3) return v2; 

  // Check if third vertex is closest
  const Vector p3 = p - v3;
  const float d5 = e1 * p3;
  const float d6 = e2 * p3;
  if (d6 >= 0.0 && d5 <= d6) return v3; 

  // Check if v1-v2 edge is closest
  const float vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
  {
      const float v = d1 / (d1 - d3);
      return v1 + v * e1; 
  }
    
  // Check if v3-v1 edge is closest
  const float vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
  {
      const float v = d2 / (d2 - d6);
      return v1 + v * e2;
  }
    
  // Check if v2-v3 edge is closest
  const float va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
  {
      const float v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      return v2 + v * (v3 - v2);
  }

  // projected point is inside triangle
  const float denom = 1.0 / (va + vb + vc);
  const float v = vb * denom;
  const float w = vc * denom;
  return v1 + v * e1 + w * e2;
}


}
