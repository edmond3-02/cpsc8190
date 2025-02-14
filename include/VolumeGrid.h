#ifndef __VOLUMEGRID_H__
#define __VOLUMEGRID_H__
 
#include "RectangularGrid.h"
#include "ProgressMeter.h"
#include "Vector.h"
#include "Color.h"
#include <iostream>
#include <fstream>
 
 
using namespace std;
 
 
namespace lux
{
 
 
template< typename DATA >
class VolumeGrid : public RectangularGrid
{
  public:
 
    VolumeGrid() :
       data (0)
    {}
 
   ~VolumeGrid()
    {
       if( data ){ delete[] data; }
    }
 
 
    void init( int nx, int ny, int nz, 
               float Lx, float Ly, float Lz,
           const Vector& Origin        )
    {
       if( data ){ delete[] data; }
       data = new DATA[ nx*ny*nz ];
       RectangularGrid::init( nx, ny, nz, Lx, Ly, Lz, Origin );
    }
 
    void setClearValue( const DATA value )
    {
       int N = nx()*ny()*nz();
       for( int i=0;i<N;i++ ){ data[i] = value; }
       clearValue = value;
    }
 
    const DATA getClearValue() const { return clearValue; }
 
    void setOutsideValue( const DATA value )
    {
       outsideValue = value;
    }
 
    const DATA getOutsideValue() const { return outsideValue; }
 
    DATA* rawPtr(){ return data; }
 
    const DATA& value( int i, int j, int k ) const { return data[ index(i,j,k) ]; }
 
          DATA& value( int i, int j, int k ) { return data[ index(i,j,k) ]; }
 
    void set( int i, int j, int k, const DATA& value ){ data[index(i,j,k)] = value; }
 
    void set( const Vector& P, const DATA& value )
    {
       if( !isInside(P) ){ return; }
       int ix, iix;
       float wx, wwx;
       int iy, iiy;
       float wy, wwy;
       int iz, iiz;
       float wz, wwz;
 
       getLinearInterpolation( P,  
                               ix, iix, wx,  wwx, 
                   iy, iiy, wy,  wwy,
                   iz, iiz, wz,  wwz );
 
       int i000 = ix  + nX*( iy + nY*iz );
       int i100 = iix  + nX*( iy + nY*iz );
       int i010 = ix  + nX*( iiy + nY*iz );
       int i110 = iix  + nX*( iiy + nY*iz );
       int i001 = ix  + nX*( iy + nY*iiz );
       int i101 = iix  + nX*( iy + nY*iiz );
       int i011 = ix  + nX*( iiy + nY*iiz );
       int i111 = iix  + nX*( iiy + nY*iiz );
       data[i000] += value * wwx * wwy * wwz;
       data[i100] += value * wx * wwy * wwz;
       data[i010] += value * wwx * wy * wwz;
       data[i110] += value * wx * wy * wwz;
       data[i001] += value * wwx * wwy * wz;
       data[i011] += value * wwx * wy * wz;
       data[i111] += value * wx * wy * wz;
       data[i101] += value * wx * wwy * wz;
    }
 
    const DATA eval( const Vector& P ) const
    {
       if( !isInside(P) ){ return outsideValue; }
       int ix, iix;
       float wx, wwx;
       int iy, iiy;
       float wy, wwy;
       int iz, iiz;
       float wz, wwz;
 
       getLinearInterpolation( P,  
                               ix, iix, wx,  wwx, 
                   iy, iiy, wy,  wwy,
                   iz, iiz, wz,  wwz );
 
       int i000 = ix  + nX*( iy + nY*iz );
       int i100 = iix  + nX*( iy + nY*iz );
       int i010 = ix  + nX*( iiy + nY*iz );
       int i110 = iix  + nX*( iiy + nY*iz );
       int i001 = ix  + nX*( iy + nY*iiz );
       int i101 = iix  + nX*( iy + nY*iiz );
       int i011 = ix  + nX*( iiy + nY*iiz );
       int i111 = iix  + nX*( iiy + nY*iiz );
       DATA density = data[i000] * wwx * wwy * wwz
            + data[i100] * wx * wwy * wwz
            + data[i010] * wwx * wy * wwz
            + data[i110] * wx * wy * wwz
            + data[i001] * wwx * wwy * wz
            + data[i011] * wwx * wy * wz
            + data[i111] * wx * wy * wz
            + data[i101] * wx * wwy * wz;
 
       return density;
   }
 
   void normalize( const VolumeGrid<float>& g )
   {
      for( int k=0;k<nz();k++ )
      {
         for( int j=0;j<ny();j++ )
     {
        for( int i=0;i<nx();i++ )
        {
           if( g.eval( evalP(i,j,k) ) > 0 )
           {
              DATA val = value(i,j,k);
          val /= g.eval( evalP(i,j,k) );
          value(i,j,k) = val;
           }
        }
     }
      }
   }
 
  protected:
 
    DATA* data;
    DATA clearValue;
    DATA outsideValue;
 
};
 
 
 
template<typename DATA>
void CopyVolumeGrid( VolumeGrid<DATA>& in, VolumeGrid<DATA>& copy )
{
   for( int k=0;k<copy.nz();k++ )
   {
      for( int j=0;j<copy.ny();j++ )
      {
         for( int i=0;i<copy.nx();i++ )
         {
        copy.value(i,j,k) = in.eval( copy.evalP(i,j,k) );
         }
      }
   }
}
 
template<typename DATA>
void WriteVolumeGrid( VolumeGrid<DATA>& grid, ofstream& out )
{
   int nx = grid.nx();
   int ny = grid.ny();
   int nz = grid.nz();
 
   Vector llc = grid.llc();
   float x = llc[0];
   float y = llc[1];
   float z = llc[2];
 
   float Lx = grid.Lx();
   float Ly = grid.Ly();
   float Lz = grid.Lz();
 
 
 
   out.write( (char*)&nx, sizeof(nx) );
   out.write( (char*)&ny, sizeof(ny) );
   out.write( (char*)&nz, sizeof(nz) );
 
   out.write( (char*)&Lx, sizeof(Lx) );
   out.write( (char*)&Ly, sizeof(Ly) );
   out.write( (char*)&Lz, sizeof(Lz) );
 
   out.write( (char*)&x, sizeof(x) );
   out.write( (char*)&y, sizeof(y) );
   out.write( (char*)&z, sizeof(z) );
 
   DATA cv = grid.getClearValue();
   out.write( (char*)&cv, sizeof(cv) );
   cv = grid.getOutsideValue();
   out.write( (char*)&cv, sizeof(cv) );
 
   DATA* value = new DATA[nx];
   int datasize = sizeof(DATA) * nx;
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
     {
        value[i] = grid.value(i,j,k);
     }
     out.write( (char*)&(*value), nx*(sizeof(DATA)) );
      }
   }
}
 
template<typename DATA>
void ReadVolumeGrid( VolumeGrid<DATA>& grid, ifstream& out )
{
   int nx;
   int ny;
   int nz;
 
   float x;
   float y;
   float z;
 
   float Lx;
   float Ly;
   float Lz;
 
   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
 
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );
 
   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );
 
   DATA cv, ov; 
   out.read( (char*)&cv, sizeof(cv) );
   out.read( (char*)&ov, sizeof(cv) );
   grid.init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );
   grid.setClearValue( cv );
   grid.setOutsideValue( ov );
 
   DATA* value = new DATA[nx];
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
     out.read( (char*)&(*value), nx*(sizeof(DATA)) );
         for( int i=0;i<nx;i++ )
     {
        grid.value(i,j,k) = value[i];
     }
 
      }
      
   }
   delete[] value;
 
}
 
 
 
void WriteFloatVolumeGrid( VolumeGrid<float>& grid, const string& fname );
void ReadFloatVolumeGrid( VolumeGrid<float>& grid, const string& fname );
 
void WriteVectorVolumeGrid( VolumeGrid<Vector>& grid, const string& fname );
void ReadVectorVolumeGrid( VolumeGrid<Vector>& grid, const string& fname );
 
void WriteColorVolumeGrid( VolumeGrid<Color>& grid, const string& fname );
void ReadColorVolumeGrid( VolumeGrid<Color>& grid, const string& fname );
 
 
 
 
 
 
/*
 
template<typename DATA>
class SparseMapGrid: public SparseRectangularGrid
{
  public:
    SparseMapGrid(){}
   ~SparseMapGrid(){}
 
    const DATA get(int i,int j,int k) const
    {
    map<int,map<int,map<int,float> > >::const_iterator iter1 = grid.find(i);
    if (iter1 == grid.end()) {
        return clearValue;
    }
 
    const map<int,map<int,float> > *holdMap1 = &(iter1->second);
    map<int,map<int,float> >::const_iterator iter2 = holdMap1->find(j);
    
    if (iter2 == holdMap1->end()) {
        return clearValue;
    }
    
    const map<int,float> *holdMap2 = &(iter2->second);
    map<int,float>::const_iterator iter3 = holdMap2->find(k);
    
    if (iter3 == holdMap2->end()) {
        return clearValue;
    }
    return iter3->second;
    }
 
    void set( DATA val, int i, int j, int k)
    {
       if( val != get(i,j,k) )
       {
          grid[i][j][k] = val;
       }
    }
 
    void setClearValue(DATA val) { clearValue = val; }
 
    
    const DATA eval( const Vector& P ) const
    {
       int ix, iix;
       float wx, wwx;
       int iy, iiy;
       float wy, wwy;
       int iz, iiz;
       float wz, wwz;
 
       getLinearInterpolation( P,  
                               ix, iix, wx,  wwx, 
                   iy, iiy, wy,  wwy,
                   iz, iiz, wz,  wwz );
 
       DATA density = get( ix, iy, iz) * wwx * wwy * wwz
            + get(iix, iy, iz) * wx * wwy * wwz
            + get( ix,iiy, iz) * wwx * wy * wwz
            + get(iix,iiy, iz) * wx * wy * wwz
            + get( ix, iy,iiz) * wwx * wwy * wz
            + get( ix,iiy,iiz) * wwx * wy * wz
            + get(iix,iiy,iiz) * wx * wy * wz
            + get(iix, iy,iiz) * wx * wwy * wz;
 
       return density;
   }
 
 
  private:
 
    map<int,map<int,map<int,DATA > > > grid;
    DATA clearValue;
};
 
*/
 
 
 
 
 
 
 
 
}
#endif
