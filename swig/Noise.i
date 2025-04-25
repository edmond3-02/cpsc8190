
%module bishop 
%{
#include "Noise.h"
%}

%include "Noise.h"
%template(FractalSumGustavson) lux::FractalSum<lux::PerlinNoiseGustavson>;
%template(FractalSumPerlin) lux::FractalSum<lux::PerlinNoise>;
%template(FractalTerrainGustavson) lux::FractalTerrain<lux::PerlinNoiseGustavson>;
%template(FractalTerrainPerlin) lux::FractalTerrain<lux::PerlinNoise>;
