#ifndef __NOISEMACHINE_H__
#define __NOISEMACHINE_H__
 
#include "Noise.h"
#include "PerlinNoise.h"
#include <boost/shared_ptr.hpp>
 
namespace lux
{
 
typedef boost::shared_ptr<Noise>  NoiseMachineBase;
 
 
class NoiseMachine : public NoiseMachineBase
{
  public:
 
    NoiseMachine() :  boost::shared_ptr<Noise>() {}
    NoiseMachine( Noise* f ) :  boost::shared_ptr<Noise>( f ) {}
   ~NoiseMachine(){}
 
};
 
 
NoiseMachine perlin( Noise_t n );
NoiseMachine terrain( Noise_t n );
 
 
 
 
 
 
}
 
#endif
