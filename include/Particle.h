#ifndef __PARTICLE_H__
#define __PARTICLE_H__
 
 
#include "Vector.h"
#include "Color.h"
 
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <memory>
 
using namespace std;
 
namespace lux
{
 
 
class Particle
{
  public:
 
    Particle() :
       _P           (Vector(0,0,0)),
       _Cd          (Color(1,1,1,1)),
       _id          (0),
       _pscale      (1.0),
       _octaves     (1.0),
       _roughness   (0.5),
       _freq        (1.0),
       _fjump       (2.0),
       _offset      (0.0),
       _translate   (Vector(0,0,0)),
       _fade        (0.0),
       _opacity     (1.0),
       _normal      (Vector(0,0,1)),
       _right       (Vector(1,0,0)),
       _up          (Vector(0,1,0)),
       _v              (Vector(0,0,0)),
       _accel          (Vector(0,0,0)),
       _shutter        (1.0),
       _framerate      (24.0),
       _lifetime       (0),
       _age            (0),
       _nbWisps        (1),
       _wisp_octaves   (1.0),
       _wisp_roughness (0.5),
       _wisp_freq      (1.0),
       _wisp_fjump     (2.0),
       _wisp_offset    (0.0),
       _wisp_translate (Vector(0,0,0)),
       _wisp_correlation(0.0),
       _wisp_radialgroup (1.0/3.0),
       _wisp_displacementScale (0.0),
       _pyroAmplitude  (1.0),
       _pyroGamma      (1.0),
       _pyroDensity    (1.0)
    {}
 
   ~Particle(){}
 
    const Vector& P() const { return _P; }
          Vector& P()       { return _P; }
 
    const Color& Cd() const { return _Cd; }
          Color& Cd()       { return _Cd; }
 
    const int& id() const { return _id; }
          int& id()       { return _id; }
 
    const float& pscale() const { return _pscale; }
          float& pscale()       { return _pscale; }
 
    const float& octaves() const { return _octaves; }
          float& octaves()       { return _octaves; }
 
    const float& roughness() const { return _roughness; }
          float& roughness()       { return _roughness; }
 
    const float& freq() const { return _freq; }
          float& freq()       { return _freq; }
 
    const float& fjump() const { return _fjump; }
          float& fjump()       { return _fjump; }
 
    const float& offset() const { return _offset; }
          float& offset()       { return _offset; }
 
    const Vector& translate() const { return _translate; }
          Vector& translate()       { return _translate; }
 
    const float& fade() const { return _fade; }
          float& fade()       { return _fade; }
 
    const float& opacity() const { return _opacity; }
          float& opacity()       { return _opacity; }
 
    const float& pyroAmplitude() const { return _pyroAmplitude; }
          float& pyroAmplitude()       { return _pyroAmplitude; }
 
    const float& pyroGamma() const { return _pyroGamma; }
          float& pyroGamma()       { return _pyroGamma; }
 
    const float& pyroDensity() const { return _pyroDensity; }
          float& pyroDensity()       { return _pyroDensity; }
 
 
    void addAttribute( const std::string& name, const float value )
    {
       std::map<std::string, float>::const_iterator check = attributes.find( name );
       if( check == attributes.end() )
       {
          attributes[ name ] = value;
       }
    }
 
 
    const float attribute( const std::string& name ) const
    {
       std::map<std::string,float>::const_iterator check = attributes.find( name );
       if( check != attributes.end() )
       {
          return (*check).second;
       }
       if( name == "PX" ){ return _P[0]; }
       if( name == "PY" ){ return _P[1]; }
       if( name == "PZ" ){ return _P[2]; }
       if( name == "VX" ){ return _v[0]; }
       if( name == "VY" ){ return _v[1]; }
       if( name == "VZ" ){ return _v[2]; }
       return 0;
    }
 
    // Wisp attributes
    const Vector& v() const { return _v; }
          Vector& v()       { return _v; }
 
    const Vector& accel() const { return _accel; }
          Vector& accel()       { return _accel; }
 
    const Vector& normal() const { return _normal; }
          Vector& normal()       { return _normal; }
 
    const Vector& right() const { return _right; }
          Vector& right()       { return _right; }
 
    const Vector& up() const { return _up; }
          Vector& up()       { return _up; }
 
    const int& nbWisps() const { return _nbWisps; }
          int& nbWisps()       { return _nbWisps; }
 
    const float& wispOctaves() const { return _wisp_octaves; }
          float& wispOctaves()       { return _wisp_octaves; }
 
    const float& wispRoughness() const { return _wisp_roughness; }
          float& wispRoughness()       { return _wisp_roughness; }
 
    const float& wispFreq() const { return _wisp_freq; }
          float& wispFreq()       { return _wisp_freq; }
 
    const float& wispFjump() const { return _wisp_fjump; }
          float& wispFjump()       { return _wisp_fjump; }
 
    const float& wispOffset() const { return _wisp_offset; }
          float& wispOffset()       { return _wisp_offset; }
 
    const Vector& wispTranslate() const { return _wisp_translate; }
          Vector& wispTranslate()       { return _wisp_translate; }
 
    const float& shutter() const { return _shutter; }
          float& shutter()       { return _shutter; }
 
    const float& framerate() const { return _framerate; }
          float& framerate()       { return _framerate; }
 
    const float& lifetime() const { return _lifetime; }
          float& lifetime()       { return _lifetime; }
 
    const float& age() const { return _age; }
          float& age()       { return _age; }
 
    const float& wispCorrelation() const { return _wisp_correlation; }
          float& wispCorrelation()       { return _wisp_correlation; }
 
    const float& wispRadialGroup() const { return _wisp_radialgroup; }
          float& wispRadialGroup()       { return _wisp_radialgroup; }
 
    const float& wispDisplacementScale() const { return _wisp_displacementScale; }
          float& wispDisplacementScale()       { return _wisp_displacementScale; }
 
    
 
  //protected:
 
    Vector _P;
    Color _Cd;
    int _id;
    float _pscale;
    float _octaves;
    float _roughness;
    float _freq;
    float _fjump;
    float _offset;
    Vector _translate;
    float _fade;
    float _opacity;
 
    Vector _normal, _right, _up;
    Vector _v;
    Vector _accel;
    float _shutter;
    float _framerate;
    float _lifetime;
    float _age;
 
    // additional parameters for point wisps
    int _nbWisps;
    float _wisp_octaves;
    float _wisp_roughness;
    float _wisp_freq;
    float _wisp_fjump;
    float _wisp_offset;
    Vector _wisp_translate;
    float _wisp_correlation;
    float _wisp_radialgroup;
    float _wisp_displacementScale;
 
    //additional parameters for pyroclastics
    float _pyroAmplitude;
    float _pyroGamma;
    float _pyroDensity;
 
    std::map<std::string, float> attributes;
 
 
};
 
void Interpolate( const Particle& p0, const Particle& p1, const float z, Particle& interp );
 
 
 
typedef vector<Particle> ParticleGroupA;
 
 
/*
typedef std::shared_ptr<Particle> Dot;
typedef std::shared_ptr<vector<Dot> > DotArray; 
 
class DotSet : public boost::shared_ptr<vector<Particle> >
{
  public:
 
    DotSet();
    DotSet( size_t nb );
    DotSet( const DotSet& ds );
 
   ~DotSet();
 
    const DotSet operator+( const DotSet& ds ) const;
    const DotSet subset( const size_t begin, const size_t end ) const;
 
    size_t nbDots() const;
 
};
*/
 
}
 
#endif
