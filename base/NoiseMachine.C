#include "NoiseMachine.h"
#include "PerlinNoise.h"


using namespace lux;


NoiseMachine lux::perlin(Noise_t n)
{
	NoiseMachine noise = NoiseMachine(new FractalSum<PerlinNoiseGustavson>());
	noise->setParameters(n);
	return noise;
}
NoiseMachine lux::terrain(Noise_t n)
{
	NoiseMachine noise = NoiseMachine(new FractalTerrain<PerlinNoiseGustavson>());
	noise->setParameters(n);
	return noise;
}

float lux::evaluate(const NoiseMachine& n, const Vector& X)
{
	return n->eval(X);
}
