#ifndef __VOLUMIZED_GRID_H__
#define __VOLUMIZED_GRID_H__
 
 
 
#include "Volume.h"
 
#include "VolumeGrid.h"
#include "SparseGrid.h"
 
namespace lux
{
 
template<typename U>
class VolumizedGrid : public Volume<U> 
{
  public:
 
 
   typedef U volumeDataType;
   typedef typename GradType<U>::GType volumeGradType;

    VolumizedGrid( const VolumeGrid<U>* g ) :
       grid (g),
       sparsegrid(0)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
    }
 
    VolumizedGrid( const SGrid<U>* g ) :
       grid (0),
       sparsegrid(g)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
    }
 
 
   ~VolumizedGrid(){}
 
    const volumeDataType eval( const Vector& P ) const { return ( grid != 0 ) ? grid->eval(P) : sparsegrid->eval(P); }
    const volumeGradType grad( const Vector& P ) const
    {
       float e = eval(P);
       Vector value( (eval(P+Vector(dx,0,0))-e)/dx,
                     (eval(P+Vector(0,dy,0))-e)/dy,
                     (eval(P+Vector(0,0,dz))-e)/dz    );
       return value;
    }
 
  private:
 
    const VolumeGrid<U>* grid;
    const SGrid<U>* sparsegrid;
    float dx, dy, dz;
};
 
}
  
#endif
