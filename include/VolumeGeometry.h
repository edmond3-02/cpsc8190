#ifndef __VOLUMEGEOMETRY_H__
#define __VOLUMEGEOMETRY_H__
 
#include <cmath>
#include <vector>
#include "Volume.h"
#include "VolumizedGrid.h"
#include "VolumeGrid.h"
//#include "Particle.h"
#include "ProgressMeter.h"
#include "GeometryVolumeShapes.h"
 
 
 
namespace lux
{
 
 
class TriangleGeometry
{
  public:
 
    TriangleGeometry():
      llc   (Vector(0,0,0)),
      urc   (Vector(0,0,0)),
      scal (1.0)
   {}
   ~TriangleGeometry(){}
 
    void merge( const TriangleGeometry& g );
 
    typedef std::vector<int> Face;
 
    void addVertex( const Vector& v ) 
    { 
	llc[0] = std::min(llc.X(), v.X());
	llc[1] = std::min(llc.Y(), v.Y());
	llc[2] = std::min(llc.Z(), v.Z());
    	
	urc[0] = std::max(urc.X(), v.X());
	urc[1] = std::max(urc.Y(), v.Y());
	urc[2] = std::max(urc.Z(), v.Z());
    	
    	vertices.push_back(v); 
    }

    void addNormal( const Vector& v ) { normals.push_back(v); }
    void addTextureCoordinate( const Vector& v );
    const Vector& getVertex( const size_t i ) const { return vertices[i]; }
    const Vector& getNormal( const size_t i ) const { return normals[i]; }
    const Vector faceNormal( const size_t i ) const;
    const double faceArea( const size_t i ) const;
    const Vector& getTextureCoordinate( const size_t i ) const { return textureCoordinates[i]; }
    const size_t nbVertices() const { return vertices.size(); }
    const size_t nbTextureCoordinates() const { return textureCoordinates.size(); }
    const size_t nbNormals() const { return normals.size(); }
 
    void setVertex( const size_t i, const Vector& v ) { vertices[i] = v; return; }
    void setNormal( const size_t i, const Vector& v ) { normals[i] = v; return; }
 
    void addFace( const int i, const int j, const int k )
    {
       Face f;
       f.push_back(i);
       f.push_back(j);
       f.push_back(k);
       faces.push_back(f);
    }
 
    void addFace( const Face& f )
    {
       faces.push_back(f);
    }
 
    void addTexturedFace( const int i, const int j, const int k, const int it, const int jt, const int kt )
    {
       Face f;
       f.push_back(i);
       f.push_back(j);
       f.push_back(k);
       faces.push_back(f);
       Face tf;
       tf.push_back(it);
       tf.push_back(jt);
       tf.push_back(kt);
       texturedfaces.push_back(tf);
 
    }
 
    void addTexturedNormaledFace( const int i, const int j, const int k, const int it, const int jt, const int kt, const int in, const int jn, const int kn )
    {
       Face f;
       f.push_back(i);
       f.push_back(j);
       f.push_back(k);
       faces.push_back(f);
       Face tf;
       tf.push_back(it);
       tf.push_back(jt);
       tf.push_back(kt);
       texturedfaces.push_back(tf);
       Face nf;
       nf.push_back(in);
       nf.push_back(jn);
       nf.push_back(kn);
       normaledfaces.push_back(nf);
    }
 
 
    void addTexturedFace( const Face& f, const Face& tf )
    {
       faces.push_back(f);
       texturedfaces.push_back(tf);
    }
 
 
 
    void getFace( const size_t f, int& i, int& j, int& k ) const
    {
       i = faces[f][0];
       j = faces[f][1];
       k = faces[f][2];
    }
 
    int getFaceIndex( const size_t f, int index ) const
    {
       return faces[f][index];
    }
 
 
    void getTexturedFace( const size_t f, int& i, int& j, int& k, int& it, int& jt, int& kt ) const
    {
       i = faces[f][0];
       j = faces[f][1];
       k = faces[f][2];
       it = texturedfaces[f][0];
       jt = texturedfaces[f][1];
       kt = texturedfaces[f][2];
    }
 
    void getTexturedNormaledFace( const size_t f, int& i, int& j, int& k, int& it, int& jt, int& kt, int& in, int& jn, int& kn ) const
    {
       i = faces[f][0];
       j = faces[f][1];
       k = faces[f][2];
       it = texturedfaces[f][0];
       jt = texturedfaces[f][1];
       kt = texturedfaces[f][2];
       in = normaledfaces[f][0];
       jn = normaledfaces[f][1];
       kn = normaledfaces[f][2];
    }
 
 
 
    const size_t nbFaces() const { return faces.size(); }
 
    void clear(){ vertices.clear(); faces.clear(); textureCoordinates.clear(); texturedfaces.clear(); connectivity.clear(); }
 
    const Vector& LLC() const { return llc; }
    const Vector& URC() const { return urc; }
 
    void setScaling( const float sc ) { scal = sc; }
    const float scaling() const { return scal; }
 
    void computeConnectivity();
 
    const std::vector<int>& connections( int i ) const { return connectivity[i]; }
 
    const float averageNeighborDistance( int i ) const;
 
    bool hasTextureCoordinates() const { return !textureCoordinates.empty(); }
    bool hasNormals() const { return !normals.empty(); }
 
 
    void translate( const Vector& t );
    void rotate( const Vector& t );
 
 
  protected:
 
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> textureCoordinates;
    std::vector<Face> faces;
    std::vector<Face> texturedfaces;
    std::vector<Face> normaledfaces;
    std::vector< std::vector<int> > connectivity;
 
    Vector llc, urc;
 
    float scal;
 
 
 
};
 
 
void writeObj( const string& fname, const TriangleGeometry& g, const bool use_normals = true, const bool compress = false );
 
 
 
 
typedef std::shared_ptr<TriangleGeometry> MeshBase;
 
 
class Mesh : public MeshBase
{
  public:
 
    Mesh()  :  std::shared_ptr<TriangleGeometry>() {}
    Mesh( TriangleGeometry* f ) :  std::shared_ptr<TriangleGeometry>(f) {}
   ~Mesh(){}
 
     char* __str__() 
     {
       static char typeLabel[2048];
       std::string lbl = "Mesh";
       size_t lbllength = lbl.size();
       if( lbllength > 2047 ){ lbllength = 2047; }
       lbllength = lbl.copy( typeLabel, lbllength);
       typeLabel[lbllength] = '\0';
       return typeLabel;
    }
 
 
 
};
 
Mesh createEmptyMesh();
const size_t nbMeshVertices(const Mesh& m);
const size_t nbMeshNormals(const Mesh& m);
const size_t nbMeshTextureCoordinates(const Mesh& m);
const size_t nbMeshFaces(const Mesh& m);
 
 
void writeObj( const string& fname, const Mesh& g, const bool use_normals = true, const bool compress = false );
 
 
class SignedDistance : public Volume<float>
{
  public:
 
    SignedDistance( const TriangleGeometry& geom );
    SignedDistance( const Triangle& geom );
   ~SignedDistance(){}
   
    void addTriangle( const Triangle& t ){ geometry.push_back(t); }
    const float eval( const Vector& P ) const;
 
 
  private:
 
    vector<Triangle> geometry;
};
 
 
class TexturedSignedDistance : public Volume<Vector>
{
  public:
 
    TexturedSignedDistance( const TriangleGeometry& geom );
    TexturedSignedDistance( const TexturedTriangle& geom );
   ~TexturedSignedDistance(){}
   
    void addTriangle( const TexturedTriangle& t ){ geometry.push_back(t); }
    const Vector eval( const Vector& P ) const;
 
 
  private:
 
    vector<TexturedTriangle> geometry;
};
 
 
 
 
 
Volume<float> * ProcessLevelSet( const TriangleGeometry& geom );
Volume<float> * ProcessLevelSet( const TriangleGeometry& geom, VolumeGrid<float>& lsgrid, const bool flip = false );
 
 
 
void RayMarchLevelSet( const TriangleGeometry& geom, VolumeGrid<float>& lsgrid, const float threshold = 0.1, int samps = 100 );
 
void RayMarchLevelSet( const TriangleGeometry& geom, ScalarGrid& lsgrid );
 
void RayMarchLevelSet( const TriangleGeometry& geom, VectorGrid& lsgrid );
 
 
 
 
 
 
 
Volume<float> * ProcessAABox( const Vector& llc, const Vector& urc );
Volume<float> * ProcessAABox( const Vector& llc, const Vector& urc, VolumeGrid<float>& lsgrid  );
 
Volume<float> * ProcessPyramid( const float length, const Vector center = Vector(0,0,0) );
Volume<float> * ProcessPyramid( const float length, const Vector center,  VolumeGrid<float>& lsgrid );
 
 
 
vector<float> FindAllIntersections( const vector<Triangle>& g, const Vector P, const Vector D, const float dmax, const float threshold );
vector<Vector> FindAllIntersections( const vector<TexturedTriangle>& g, const Vector P, const Vector D, const float dmax, const float threshold );
 
 
/* 
 
void Geometry2Particles( const TriangleGeometry& geom, ParticleGroupA& particles );
 
template<typename U>
void Blur( VolumeGrid<U>& g )
{
   VolumeGrid<U> temp;
   temp.init( g.nx(), g.ny(), g.nz(), g.Lx(), g.Ly(), g.Lz(), g.llc() ); 
   ProgressMeter meter( (g.nx())*(g.ny())*(g.nz()) * 2, "Blur" );
   for( int k=0;k<g.nz();k++ )
   {
      for( int j=0;j<g.ny();j++ )
      {
         for( int i=0;i<g.nx();i++ )
         {
        temp.value(i,j,k) = g.value(i,j,k);
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
        U sum;
        int nb = 0;
        for( int kk=kmin;kk<=kmax;kk++ )
        {
           for( int jj=jmin;jj<=jmax;jj++ )
           {
              for( int ii=imin;ii<=imax;ii++ )
              {
             sum += temp.value( ii, jj, kk );
             ++nb;
              }
           }
        }
        sum = sum/nb;
        g.value(i,j,k) = sum;
        meter.update();
         }
      }
   }
}
 
}
 
 */
 
#endif
