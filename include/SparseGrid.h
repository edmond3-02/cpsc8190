#ifndef __SPARSEGRID_H__
#define __SPARSEGRID_H__
       
#include "RectangularGrid.h"
//#include "FrustumGrid.h"
#include "ProgressMeter.h"
#include "Color.h"
#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <map>
#include <boost/shared_ptr.hpp>
   
using namespace std;
    
          
namespace lux
{  

template <typename T>
class SGrid: public RectangularGrid
{
  public: 
   SGrid( int psize = 4 ) :
      partitionSize (psize),
      nbOccupied (0),
      data(NULL)
   {
       T empty;
       defVal = empty * 0.0;
   }
 
 
   void init(int Nx,int Ny,int Nz, float lx,float ly,float lz,const Vector & Origin)
   {
      deallocate();
      nnx = Nx/partitionSize;
      nny = Ny/partitionSize;
      nnz = Nz/partitionSize;
      if( nnx*partitionSize < Nx ){ nnx++; }
      if( nny*partitionSize < Ny ){ nny++; }
      if( nnz*partitionSize < Nz ){ nnz++; }
      int NNx = nnx * partitionSize;
      int NNy = nny * partitionSize;
      int NNz = nnz * partitionSize;
      RectangularGrid::init(NNx, NNy, NNz, lx, ly, lz, Origin);
      totalSize = nnx*nny*nnz;
      data = new T* [totalSize];
      for(int i = 0; i < totalSize; i++)
      {
         data[i] = NULL;
      }
      nbOccupied = 0;
   }
 
   void clear()
   {
      deallocate();
      data = new T* [totalSize];
      for(int i = 0; i < totalSize; i++)
      {
         data[i] = NULL;
      }
      nbOccupied = 0;
   }
 
   const T& get(int i,int j,int k) const
   {
      if(!isInside(i,j,k)){ return(defVal); }
      int myIndex = sindex(i, j, k);
      if(data[myIndex] == NULL)
      {  
         return(defVal); 
      }
      else
      {
         int ii = i % partitionSize;
         int jj = j % partitionSize;
         int kk = k % partitionSize;
         int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
         return data[myIndex][partitionIndex];
      } 
   }
 
   void setDefVal( const T& def ) { defVal = def; }
 
   const T& getDefVal() const { return defVal; }
 
   void set( int i,int j,int k, const T& val )
   {
         long myIndex = sindex(i, j, k);  
         if(data[myIndex] == NULL)
     {
            data[myIndex] = new T[partitionSize * partitionSize * partitionSize];
            for(int i = 0; i < partitionSize * partitionSize * partitionSize; i++)
        {
               data[myIndex][i] = defVal;
            }
        nbOccupied++;
         }
         int ii = i % partitionSize;
         int jj = j % partitionSize;
         int kk = k % partitionSize;
         int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
         data[myIndex][partitionIndex] = val;
   }
 
   const T eval( const Vector& P ) const
   {
       if( !isInside(P) ){ return defVal; }
 
/* 
       T accum = defVal * 0.0;
       vector<int> indices_x, indices_y, indices_z;
       vector<double> weights_x, weights_y, weights_z;

       getHighOrderInterpolation( P, indices_x, indices_y, indices_z, weights_x, weights_y, weights_z );
 
       if( indices_x.empty() || indices_y.empty() || indices_z.empty() ){ return defVal; }
 
       for( size_t k=0;k<indices_z.size();k++ )
       {
          for( size_t j=0;j<indices_y.size();j++ )
          {
             for( size_t i=0;i<indices_x.size();i++ )
             {
                accum = accum + get( indices_x[i], indices_y[j], indices_z[k] ) * weights_x[i] * weights_y[j] * weights_z[k];
             }
          }
       }
 	*/


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

	//std::cout << P.__str__();
	//std::cout<<ix<<" "<<iy<<" "<<iz<<"/"<<iix<<" "<<iiy<<" "<<iiz<<std::endl;
 
       T accum = get( ix, iy, iz ) * wwx * wwy * wwz;
       accum  += get( iix, iy, iz ) *  wx * wwy * wwz;
       accum  += get( ix, iiy, iz ) * wwx *  wy * wwz;
       accum  += get( iix, iiy, iz ) *  wx *  wy * wwz;
       accum  += get( ix, iy, iiz ) * wwx * wwy *  wz;
       accum  += get( iix, iy, iiz ) *  wx * wwy *  wz;
       accum  += get( ix, iiy, iiz ) * wwx *  wy *  wz;
       accum  += get( iix, iiy, iiz )*  wx *  wy *  wz;


       return accum;
   }
 
   ~SGrid()
   {
     deallocate();
   }
 
   void deallocate()
   {
      // Already deallocated.
      if(data == NULL)
      {
        return;
      }
 
      //totalSize = nnx*nny*nnz;
      for(long i = 0; i < totalSize; i++)
      {
         if(data[i] != NULL){ delete [] data[i]; }
      }
      delete [] data;
      data = NULL;
      nbOccupied = 0;
   }
 
   long Size() const { return nbOccupied; }
   long NbPartitions() const { return totalSize; }
 
   int blockSize() const { return partitionSize; }
   const bool goodBlock( const int i, const int j, const int k ) const { return ( data[sindex(i,j,k)] != NULL ); }
   const bool goodBlock( const int i ) const { return ( data[i] != NULL ); }
 
   void blockBounds( long block, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1 ) const
   {
      if ( block <0 || block >= totalSize )
      {
         i0 = j0 = k0 = i1 = j1 = k1 = -1;
         return;
      }
 
      k0 = block / (nnx*nny);
      long bblock = block - k0*nnx*nny;
      j0 = bblock / nnx;
      i0 = bblock - j0*nnx;
 
      i0 *= partitionSize;
      j0 *= partitionSize;
      k0 *= partitionSize;
 
      i1 = i0 + partitionSize;
      j1 = j0 + partitionSize;
      k1 = k0 + partitionSize;
 
      return;
   }
 
   const long sindex(int i,int j,int k) const
   {
      int ii = i/partitionSize;
      int jj = j/partitionSize;
      int kk = k/partitionSize;
      return (long)ii + (long)nnx*( (long)jj + (long)nny*(long)kk );
   }

 
  private:
      T **data;
      T defVal;
      int partitionSize;
      int nnx, nny, nnz;
      long nbOccupied;
      long totalSize;
 
 
/*
      int setUpSpace(int,int,int);
      int *clMap;
      int *clMapMap;
      float *clData;
      Vector4d* clColor;
      int numUsedMap;
      int numUsedMapMap;
*/
 
}; 
 
 
template<typename T>
void WriteVolumeGrid(SGrid<T>& grid, ofstream& out )
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
 
   T temp = grid.getDefVal();
   out.write( (char*)&temp, sizeof(T));
   
   for( int k=0;k<nz;k++ )   
   {    
      for( int j=0;j<ny;j++ )  
      {
         for( int i=0;i<nx;i++ )
         {
            temp = grid.get(i,j,k);
            if(temp != grid.getDefVal()){
               out.write((char*)&i, sizeof(int));
               out.write((char*)&j, sizeof(int));
               out.write((char*)&k, sizeof(int));
               out.write((char*)&temp, sizeof(T));
            }
         }
 
      }
   
   }
 
}
 
template<typename T>
void ReadVolumeGrid(SGrid<T>& grid, ifstream& out)
{
   int nx;
   int ny;
   int nz;
   
   int i, j, k;
   
   float x;
   float y;
   float z;
   
   float Lx;
   float Ly;
   float Lz;
   
   T value;
   
   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );
       
   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );
 
   out.read((char *)&value, sizeof(T));
   grid.setDefValue( value ); 
   grid.init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );
               
   while(out.read((char *)&i, sizeof(int)))
   {
       out.read((char *)&j, sizeof(int));
       out.read((char *)&k, sizeof(int));
       out.read((char *)&value, sizeof(T));
       grid.set(value, i, j, k);
   }
}
 
template<typename T>
void Blur( SGrid<T>& g )
{
   SGrid<T> temp;
   temp.setDefVal( g.getDefVal() );
   temp.init( g.nx(), g.ny(), g.nz(), g.Lx(), g.Ly(), g.Lz(), g.llc() ); 
   ProgressMeter meter( (g.nx())*(g.ny())*(g.nz()) * 2, "Blur" );
   for( int k=0;k<g.nz();k++ )
   {
      for( int j=0;j<g.ny();j++ )
      {
         for( int i=0;i<g.nx();i++ )
         {
        temp.set( g.get(i,j,k), i,j,k);
        meter.update();
         }
      }
   }
 
 
   for( int k=0;k<g.nz();k++ )
   {
      int kmin = k-1;
      if( kmin < 0 ){ kmin = k; }
      int kmax = k+1;
      if( kmax >= g.nz() ){ kmax = k;}
      for( int j=0;j<g.ny();j++ )
      {
         int jmin = j-1;
         if( jmin < 0 ){ jmin = j; }
         int jmax = j+1;
         if( jmax >= g.ny() ){ jmax = j;}
         for( int i=0;i<g.nx();i++ )
         {
            int imin = i-1;
            if( imin < 0 ){ imin = i; }
            int imax = i+1;
            if( imax >= g.nx() ){ imax = i;}
        T sum;
        sum *= 0.0;
        int nb = 0;
        for( int kk=kmin;kk<=kmax;kk++ )
        {
           for( int jj=jmin;jj<=jmax;jj++ )
           {
              for( int ii=imin;ii<=imax;ii++ )
              {
             sum += temp.get( ii, jj, kk );
             ++nb;
              }
           }
        }
        sum = sum/nb;
        g.set(sum, i,j,k);
        meter.update();
         }
      }
   }
}
 
       
 
 
class ScalarField;
class VectorField;
class ColorField;
class MatrixField;
 
 
 
 
typedef boost::shared_ptr<SGrid<float> > ScalarGridBase;
typedef boost::shared_ptr<SGrid<Vector> > VectorGridBase;
typedef boost::shared_ptr<SGrid<Color> > ColorGridBase;
typedef boost::shared_ptr<SGrid<Matrix> > MatrixGridBase;
 
class VectorGrid;
class ColorGrid;
class MatrixGrid;
 
class ScalarGrid : public ScalarGridBase
{
  public:
 
    ScalarGrid();
    ScalarGrid( SGrid<float>* f );
   ~ScalarGrid();
 
 
    ScalarGrid operator+=( const ScalarField& e2 );
    ScalarGrid operator-=( const ScalarField& e2 );
    ScalarGrid operator-();
    ScalarGrid operator*=( const ScalarField& e2 );
    ScalarGrid operator/=( const ScalarField& e2 );

};
 
 
class VectorGrid : public VectorGridBase
{
  public:
 
    VectorGrid( SGrid<Vector>* f );
   ~VectorGrid();
 
 
    VectorGrid operator+=( const VectorField& e2 );
    VectorGrid operator-=( const VectorField& e2 );
    VectorGrid operator-();
    VectorGrid operator*=( const ScalarField& e2 );
    VectorGrid operator/=( const ScalarField& e2 );
 
 
};
 
 
 
 
 
class ColorGrid : public ColorGridBase
{
  public:
 
    ColorGrid( SGrid<Color>* f );
   ~ColorGrid();
 
 
    ColorGrid operator+=( const ColorField& e2 );
    ColorGrid operator-=( const ColorField& e2 );
    ColorGrid operator-();
    ColorGrid operator*=( const ScalarField& e2 );
    ColorGrid operator/=( const ScalarField& e2 );
 
 
};
 
 
 
class MatrixGrid : public MatrixGridBase
{
  public:
 
    MatrixGrid( SGrid<Matrix>* f );
   ~MatrixGrid();
 
 
    MatrixGrid operator+=( const MatrixField& e2 );
    MatrixGrid operator-=( const MatrixField& e2 );
    MatrixGrid operator-();
    MatrixGrid operator*=( const ScalarField& e2 );
    MatrixGrid operator/=( const ScalarField& e2 );
 
 
}; 
 
 
 
 
 
 
 
 
 
 
 
/* 
template <typename T>
class FSGrid: public FrustumGrid
{
  public: 
   FSGrid( int psize = 4 ) :
      partitionSize (psize),
      nbOccupied (0)
   {
       T empty;
       defVal = empty * 0.0;
   }
 
 
   void init(int Nx,int Ny,int Nz, const Camera& cam)
   {
      nnx = Nx/partitionSize;
      nny = Ny/partitionSize;
      nnz = Nz/partitionSize;
      if( nnx*partitionSize > Nx ){ nnx++; }
      if( nny*partitionSize > Ny ){ nny++; }
      if( nnz*partitionSize > Nz ){ nnz++; }
      int NNx = nnx * partitionSize;
      int NNy = nny * partitionSize;
      int NNz = nnz * partitionSize;
      setEyeViewUp( cam.eye(), cam.view(), cam.up() );
      setFov( cam.fov() );
      setAspectRatio( cam.aspectRatio() );
      setNearPlane( cam.nearPlane() );
      setFarPlane( cam.farPlane() );
      FrustumGrid::init(NNx, NNy, NNz);
      totalSize = nnx*nny*nnz;
      data = new T* [totalSize];
      for(int i = 0; i < totalSize; i++)
      {
         data[i] = NULL;
      }
      nbOccupied = 0;
   }
 
   const T& get(int i,int j,int k) const
   {
      size_t myIndex = index(i, j, k);
      if(data[myIndex] == NULL)
      {  
         return(defVal); 
      }
      else
      {
         int ii = i % partitionSize;
         int jj = j % partitionSize;
         int kk = k % partitionSize;
         int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
         return data[myIndex][partitionIndex];
      } 
   }
 
   void setDefVal( const T& def ) { defVal = def; }
 
   const T& getDefVal() const { return defVal; }
 
   void set( int i,int j,int k, const T& val )
   {
      if(val != defVal)
      {
         size_t myIndex = index(i, j, k);  
         if(data[myIndex] == NULL)
     {
            data[myIndex] = new T[partitionSize * partitionSize * partitionSize];
            for(int i = 0; i < partitionSize * partitionSize * partitionSize; i++)
        {
               data[myIndex][i] = defVal;
            }
        nbOccupied++;
         }
         int ii = i % partitionSize;
         int jj = j % partitionSize;
         int kk = k % partitionSize;
         int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
         data[myIndex][partitionIndex] = val;
      }
   }
 
   const T eval( const Vector& P ) const
   {
       if( !isInside(P) ){ return defVal; }
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
 
       T accum = get( ix, iy, iz ) * wwx * wwy * wwz;
       accum  += get( iix, iy, iz ) *  wx * wwy * wwz;
       accum  += get( ix, iiy, iz ) * wwx *  wy * wwz;
       accum  += get( iix, iiy, iz ) *  wx *  wy * wwz;
       accum  += get( ix, iy, iiz ) * wwx * wwy *  wz;
       accum  += get( iix, iy, iiz ) *  wx * wwy *  wz;
       accum  += get( ix, iiy, iiz ) * wwx *  wy *  wz;
       accum  += get( iix, iiy, iiz )*  wx *  wy *  wz;
 
       return accum;
   }
 
   ~FSGrid()
   {
      long totalSize = nnx*nny*nnz;
      for(long i = 0; i < totalSize; i++)
      {
         if(data[i] != NULL){ delete [] data[i]; }
      }
      delete [] data;
   }
 
   long Size() const { return nbOccupied; }
   long NbPartitions() const { return nnx*nny*nnz; }
 
   int blockSize() const { return partitionSize; }
   const bool goodBlock( const int i, const int j, const int k ) const { return ( data[index(i,j,k)] != NULL ); }
   const bool goodBlock( const int i ) const { return ( data[i] != NULL ); }
 
   void blockBounds( long block, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1 ) const
   {
      if ( block <0 || block >= totalSize || nbOccupied == 0 )
      {
         i0 = j0 = k0 = i1 = j1 = k1 = -1;
         return;
      }
 
      k0 = block / (nnx*nny);
      long bblock = block - k0*nnx*nny;
      j0 = bblock / nnx;
      i0 = bblock - j0*nnx;
 
      i0 *= partitionSize;
      j0 *= partitionSize;
      k0 *= partitionSize;
 
      i1 = i0 + partitionSize;
      j1 = j0 + partitionSize;
      k1 = k0 + partitionSize;
 
      return;
   }
 
   const long index(int i,int j,int k) const
   {
      int ii = i/partitionSize;
      int jj = j/partitionSize;
      int kk = k/partitionSize;
      return ii + nnx*( jj + nny*kk );
   }
 
  private:
      T **data;
      T defVal;
      int partitionSize;
      int nnx, nny, nnz;
      long nbOccupied;
      long totalSize;
 
 
 
}; 
 
 
 
 
typedef boost::shared_ptr<FSGrid<float> > ScalarFrustumGridBase;
typedef boost::shared_ptr<FSGrid<Vector> > VectorFrustumGridBase;
typedef boost::shared_ptr<FSGrid<Color> > ColorFrustumGridBase;
typedef boost::shared_ptr<FSGrid<Matrix> > MatrixFrustumGridBase;
 
 
 
class ScalarFrustumGrid : public ScalarFrustumGridBase
{
  public:
 
    ScalarFrustumGrid( FSGrid<float>* f );
   ~ScalarFrustumGrid();
 
 
    ScalarFrustumGrid operator+=( const ScalarField& e2 );
    ScalarFrustumGrid operator-=( const ScalarField& e2 );
    ScalarFrustumGrid operator-();
    ScalarFrustumGrid operator*=( const ScalarField& e2 );
    ScalarFrustumGrid operator/=( const ScalarField& e2 );
 
 
};
 
 
class VectorFrustumGrid : public VectorFrustumGridBase
{
  public:
 
    VectorFrustumGrid( FSGrid<Vector>* f );
   ~VectorFrustumGrid();
 
 
    VectorFrustumGrid operator+=( const VectorField& e2 );
    VectorFrustumGrid operator-=( const VectorField& e2 );
    VectorFrustumGrid operator-();
    VectorFrustumGrid operator*=( const ScalarField& e2 );
    VectorFrustumGrid operator/=( const ScalarField& e2 );
};
 
 
 
 
 
class ColorFrustumGrid : public ColorFrustumGridBase
{
  public:
 
    ColorFrustumGrid( FSGrid<Color>* f );
   ~ColorFrustumGrid();
 
 
    ColorFrustumGrid operator+=( const ColorField& e2 );
    ColorFrustumGrid operator-=( const ColorField& e2 );
    ColorFrustumGrid operator-();
    ColorFrustumGrid operator*=( const ScalarField& e2 );
    ColorFrustumGrid operator/=( const ScalarField& e2 );
 
 
};
 
 
 
class MatrixFrustumGrid : public MatrixFrustumGridBase
{
  public:
 
    MatrixFrustumGrid( FSGrid<Matrix>* f );
   ~MatrixFrustumGrid();
 
 
    MatrixFrustumGrid operator+=( const MatrixField& e2 );
    MatrixFrustumGrid operator-=( const MatrixField& e2 );
    MatrixFrustumGrid operator-();
    MatrixFrustumGrid operator*=( const ScalarField& e2 );
    MatrixFrustumGrid operator/=( const ScalarField& e2 );
 
 
};
 
 
 
*/ 
 
 
}

#endif
