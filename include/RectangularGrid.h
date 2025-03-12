#ifndef __RECTANGULARGRID_H__
#define __RECTANGULARGRID_H__
 
#include "Vector.h"
#include <vector>
#include <memory>
//#include "HighOrderInterpolator.h"
//#include "AARectangle.h"
#include <cstdio>
#include <iostream>
 
namespace lux
{
 
 
class RectangularGrid
{
  public:
 
    RectangularGrid() :
       nX   (0),
       nY   (0),
       nZ   (0),
       dX   (0),
       dY   (0),
       dZ   (0),
       lX   (0),
       lY   (0),
       lZ   (0),
       origin   (Vector(0,0,0)),
       isPeriodicX (false),
       isPeriodicY (false),
       isPeriodicZ (false)//,
       //ho_order (1)
    {}
 
   ~RectangularGrid(){}
 
 
    void init( int nx, int ny, int nz,
               double Lx, double Ly, double Lz,
               const Vector& Origin        )
    {
	origin = Origin;

	nX=nx;
	nY=ny;
	nZ=nz;

	lX=Lx;
	lY=Ly;
	lZ=Lz;
 
	dX = lX/nX;
	dY = lY/nY;
	dZ = lZ/nZ;

	topleft = Vector(origin.X()+lX, origin.Y()+lY, origin.Z()+lZ);
    }

    void init( const Vector&llc, const Vector& urc, const Vector& res )
    {
	origin = llc;
	
	nX=res.X(); 
	nY=res.Y(); 
	nZ=res.Z();
	
	lX = abs(urc.X()-llc.X());
	lY = abs(urc.Y()-llc.Y());
	lZ = abs(urc.Z()-llc.Z());

	dX = lX/nX;
	dY = lY/nY;
	dZ = lZ/nZ;
     }
 
    void reset( double Lx, double Ly, double Lz, const Vector& origin );
 
    int index( int i, int j, int k ) const;
    int index( const Vector& P ) const;
    void triple( const int indx, int& i, int& j, int& k ) const;
 
    const int& nx() const { return nX; }
    const int& ny() const { return nY; }
    const int& nz() const { return nZ; }
 
    const Vector& llc() const { return origin; }
    const Vector& urc() const { return topleft; }
 
    const double& dx() const { return dX; }
    const double& dy() const { return dY; }
    const double& dz() const { return dZ; }
    const double& Lx() const { return lX; }
    const double& Ly() const { return lY; }
    const double& Lz() const { return lZ; }
 
    const Vector evalP( int i, int j, int k ) const
    {
	float x=i*dX, y=j*dY, z=k*dZ;
	return Vector(x,y,z)+origin;
	
    }
 
    const bool isInside( const Vector& P ) const
    {

	if(P.X() > std::max(origin.X(), topleft.X()) || P.X() < std::min(origin.X(), topleft.X())) return false;
	if(P.Y() > std::max(origin.Y(), topleft.Y()) || P.Y() < std::min(origin.Y(), topleft.Y())) return false;
	if(P.Z() > std::max(origin.Z(), topleft.Z()) || P.Z() < std::min(origin.Z(), topleft.Z())) return false;

	return true;
    }

    const bool isInside( int i, int j, int k ) const
    {
	if(  i<0  || j<0  || k<0 ) return false;
	if( i>=nX || j>=nY || k>=nZ) return false;

	return true;
    }
 
    const bool getBox( const Vector& lower, const Vector& upper, int& ixl, int& ixu, int& iyl, int& iyu, int& izl, int& izu ) const;
    const bool getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const
    {
	Vector pp = P - origin;

	// X coords
	ix = int(floor(pp.X()/dX));
	
	// Y coords
	iy = int(floor(pp.Y()/dY));
	
	// Z coords
	iz = int(floor(pp.Z()/dZ));

	if(isInside(ix, iy, iz)) return true;
	else return false;

    }

    void getLinearInterpolation( const Vector& P,  
                                 int& ix, int& iix, float& wx, float& wwx, 
                 int& iy, int& iiy,  float& wy, float& wwy,
                 int& iz, int& iiz,  float& wz, float& wwz ) const
    {
	Vector pp = P - origin;
        
	// X coords
	ix = int(pp.X()/dX);
	iix = ix+1; 
	wx = (pp.X()/dX)-ix;
	wwx = 1.0-wx;
	
	// Y coords
	iy = int(pp.Y()/dY);
	iiy = iy+1; 
	wy = (pp.Y()/dY)-iy;
	wwy = 1.0-wy;
	
	// Z coords
	iz = int(pp.Z()/dZ);
	iiz = iz+1; 
	wz = (pp.Z()/dZ)-iz;
	wwz = 1.0-wz;

    }
 
    void getHighOrderInterpolation( const Vector& P, 
                                    std::vector<int>& indices_x, std::vector<int>& indices_y, std::vector<int>& indices_z, 
                                    std::vector<double>& weights_x,
                                    std::vector<double>& weights_y,
                                    std::vector<double>& weights_z ) const;
 
    virtual const Vector transform( const Vector& P ) const { return P-origin; }
 
    void setPeriodic() { isPeriodicX = isPeriodicY = isPeriodicZ = true; }
    void setPeriodicX() { isPeriodicX = true; }
    void setPeriodicY() { isPeriodicY = true; }
    void setPeriodicZ() { isPeriodicZ = true; }
 
    bool periodicX() const { return isPeriodicX; }
    bool periodicY() const { return isPeriodicY; }
    bool periodicZ() const { return isPeriodicZ; }
 
    //void setInterpolationOrder( const int o ) { ho_order = o; }
    //const int getInterpolationOrder() const { return ho_order; }
 
 
  protected:
 
    int nX, nY, nZ;
    double dX, dY, dZ;
    double lX, lY, lZ;
    Vector origin, topleft;
    bool isPeriodicX, isPeriodicY, isPeriodicZ;
 
  private:
 
    //HighOrderInterpolator ho_interp;
    //int ho_order;
 
};
 
 
RectangularGrid Union( const RectangularGrid& g1, const RectangularGrid& g2 );
 
 
typedef std::shared_ptr<RectangularGrid>  GridBoxBase;
 
 
class GridBox : public GridBoxBase
{
  public:
 
    GridBox( RectangularGrid* f );
   ~GridBox();
 
 
    GridBox operator+=( const GridBox& e2 ); // Union
 
    const Vector evalP( int i, int j, int k ) const;
    int index( const Vector& P ) const;
    char * __str__();
};
 
GridBox makeGridBox( const Vector& llc, const Vector& urc, const Vector& dx );
GridBox makeGridBox( const RectangularGrid& rg );

}
   
#endif
