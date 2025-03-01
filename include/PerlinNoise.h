#ifndef __PERLINNOISE_H__
#define __PERLINNOISE_H__
 
#include "Vector.h"
 
 
namespace lux
{
 
#define PERLINB 0x100
class PerlinNoise
{
  public:
 
    PerlinNoise() { init(); }
   ~PerlinNoise(){}
 
    void setTime( const float time ) {}
 
    const float eval( float arg ) const;
    const float eval( const Vector& vec ) const;
 
  private:
 
   int p[PERLINB + PERLINB + 2];
   float g3[PERLINB + PERLINB + 2][3];
   float g2[PERLINB + PERLINB + 2][2];
   float g1[PERLINB + PERLINB + 2];
 
   void init();
   const double fade(double t) const;
   const double lerp(double t, double a, double b) const;
   const double grad(int hash, double x, double y, double z) const;
 
};
 
 
// This is Stefan Gustavson's DSOnoise implementation, ported to a C++ class.
class PerlinNoiseGustavson
{
  public:
 
    PerlinNoiseGustavson() : time (0) { init(); }
   ~PerlinNoiseGustavson(){}
 
 
    void setTime( const float t ) { time = t; }
 
    const float eval( float arg ) const;
 
 
    const float eval( const Vector& vec ) const;
 
  private:
 
   float time;
   int perm[PERLINB + PERLINB];

   float floor(float x) const { return ((x)>0) ? ((int)x) : ((int)x - 1); }
   float fade(float t) const { return t * t * t * ( t * ( t * 6 - 15 ) + 10); }
   float lerp(float t, float a, float b) const { return a + t * (b - a); }
   void init();

   float grad1( int hash, float x ) const {
    int h = hash & 15;
    float grad = 1.0 + (h & 7);  // Gradient value 1.0, 2.0, ..., 8.0
    if (h&8) grad = -grad;         // and a random sign for the gradient
    return ( grad * x );           // Multiply the gradient with the distance
}
 
float grad2( int hash, float x, float y ) const {
    int h = hash & 7;      // Convert low 3 bits of hash code
    float u = h<4 ? x : y;  // into 8 simple gradient directions,
    float v = h<4 ? y : x;  // and compute the dot product with (x,y).
    return ((h&1)? -u : u) + ((h&2)? -2.0*v : 2.0*v);
}
 
float grad3( int hash, float x, float y , float z ) const {
    int h = hash & 15;     // Convert low 4 bits of hash code into 12 simple
    float u = h<8 ? x : y; // gradient directions, and compute dot product.
    float v = h<4 ? y : h==12||h==14 ? x : z; // Fix repeats at h = 12 to 15
    return ((h&1)? -u : u) + ((h&2)? -v : v);
}
 
float grad4( int hash, float x, float y, float z, float t ) const {
    int h = hash & 31;      // Convert low 5 bits of hash code into 32 simple
    float u = h<24 ? x : y; // gradient directions, and compute dot product.
    float v = h<16 ? y : z;
    float w = h<8 ? z : t;
    return ((h&1)? -u : u) + ((h&2)? -v : v) + ((h&4)? -w : w);
}
 
 
};
 
 
 
 
 
 
}
#endif


