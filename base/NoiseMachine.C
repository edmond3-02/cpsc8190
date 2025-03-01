#include "NoiseMachine.h"
#include "PerlinNoise.h"


namespace lux
{

NoiseMachine perlin(Noise_t n)
{
	NoiseMachine noise = NoiseMachine(new FractalSum<PerlinNoiseGustavson>());
	noise->setParameters(n);
	return noise;
}

NoiseMachine terrain(Noise_t n)
{
	NoiseMachine noise = NoiseMachine(new FractalTerrain<PerlinNoiseGustavson>());
	noise->setParameters(n);
	return noise;
}

}
