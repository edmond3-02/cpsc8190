#ifndef __STAMP_H__
#define __STAMP_H__
 
#include "Volume.h"
#include "VolumeGrid.h"
#include "SparseGrid.h"
#include "Noise.h"
#include "NoiseMachine.h"
//#include "UniformPRN.h"
#include "PerlinNoise.h"
#include "Color.h"
#include "Particle.h"
 
namespace lux
{
 
//void StampNoise( VolumeGrid<float>& grid, Noise* noise, const Vector& center, const float radius );
//void StampNoise( VolumeGrid<float>& grid, Noise* noise, const Vector& center, const float radius, const Vector& vel, const Vector& accel, const float timestep, const int seed = 757584 );
//void StampNoiseAndColor( VolumeGrid<float>& grid, VolumeGrid<Color >& cgrid,  Noise* noise, const Vector& center, const float radius, const float fade, const Color& Cd );
//void StampNoiseAndColor( VolumeGrid<float>& grid, VolumeGrid<Color >& cgrid,  Noise* noise, const Vector& center, const float radius, const float fade, const Color& Cd, const Vector& vel, const Vector& accel, const float timestep, const int seed = 757584  );
// 
//void StampNoise( VolumeGrid<float>& grid, Noise* noise, const ParticleGroupA& particles );
//void StampNoiseAndColor( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid,  Noise* noise, const ParticleGroupA& particles );
// 
// 
//Volume<Vector>* StampSELMA( VolumeGrid<Vector>& Xgrid, Volume<Vector>* X, Volume<Vector>* velocity, float dt );
//Volume<Vector>* StampBFECCSELMA( VolumeGrid<Vector>& Xgrid, Volume<Vector>* X, Volume<Vector>* velocity, float dt );


/* 
class PointWispWanderer
{
  public:
 
    PointWispWanderer(){}
   ~PointWispWanderer(){}
 
    void reset( const Particle& p );
    void reset( const Anchor& p );
    const Vector step();
    const Vector& pos() const { return walkPosition; }
 
 
  private:
 
    Particle part;
    FractalSum<PerlinNoiseGustavson> displacementNoise;
    FractalSum<PerlinNoiseGustavson> transformNoise;
    UniformPRN walkNoise;
 
    Vector walkPosition;
    Vector d1, d2;
};
 */
 /*
class SplineWispWanderer
{
  public:
 
    SplineWispWanderer(){}
   ~SplineWispWanderer(){}
 
    void reset( const Particle& p0, const Particle& p1 );
    const Vector step();
    const Vector& pos() const { return walkPosition; }
    const Particle& interpolatedParticle(){ return interpolated; }
 
 
  private:
 
    Particle part0, part1;
    FractalSum<PerlinNoiseGustavson> displacementNoise;
    FractalSum<PerlinNoiseGustavson> transformNoise;
    UniformPRN walkNoise;
 
    Noise_t parms;
    Vector walkPosition;
    Vector d1, d2;
 
    Particle interpolated;
};
 */
 
//void StampPointWisps( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles );
//void StampPointWisps( VolumeGrid<float>& grid, const ParticleGroupA& particles );
// 
// 
//void StampPointWisps( SparseGrid& grid, SparseColorGrid& cgrid, const vector<string> files, const vector<Color> typeColors, vector<float>& opacity, float blurScale  );
// 
//void StampPointWisps( SparseGrid& grid, const ParticleGroupA& particles );
//void StampPointWisps( SparseGrid& grid, SparseColorGrid& cgrid, const ParticleGroupA& particles );
//void StampFlameWisps( SparseGrid& grid, SparseColorGrid& cgrid, const ParticleGroupA& particles, const Volume<Color>& baseColor );
// 
//void StampBlurredWisps( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& cd, int seed );
//void StampBlurredWisps( VolumeGrid<float>& grid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, int seed );
//void StampBlurredWisps( SparseGrid& grid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, int seed );
//void StampBlurredWisps( SparseGrid& grid, SparseColorGrid& cgrid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& cd, int seed );
 


//void StampField( VolumeGrid<float>& grid, Volume<float>* field );
void StampField( ScalarGrid& grid, ScalarField& field );
void StampField( ScalarGrid& grid, ScalarField& field, int bandwidth);
 
 
//void StampPyro(VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const Particle& part );
//void StampPyro( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles );
void StampPyro( ScalarGrid& grid, const AnchorChain& particles );
 
//void StampParticles( VolumeGrid<float>& grid, const ParticleGroupA& particles, const float timestep );
//void StampParticles( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles, const float timestep );
 
//void StampGridPattern( SparseGrid& grid, SparseColorGrid& cgrid, int spacing, const Color& col );
//void StampGridPattern( SparseGrid& grid, int spacing );
 
 
 
//void StampSplineWisps( ScalarGrid& grid, ColorGrid& cgrid, const ParticleGroupA& particles );
 
 
//void StampBlurredWisps( ScalarGrid& grid, ColorGrid& cgrid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& cd, int seed );
 
//void StampNoise( ScalarGrid& grid, Noise* noise, const Vector& center, const float radius, const float fade );
//void StampNoise( ScalarGrid& grid, Noise* noise, const Vector& center, const float radius, const Vector& vel, const Vector& accel, const float timestep, const float fade, const int seed = 757584 );
void StampNoise( ScalarGrid& grid, const AnchorChain& particles );
//void StampNoiseAndColor( ScalarGrid& grid, ColorGrid& cgrid,  Noise* noise, const AnchorChain& particles );
 
 
//void StampPointWisps( ScalarGrid& grid, const AnchorChain& particles );
//void StampPointWisps( ScalarGrid& grid, ColorGrid& cgrid, const AnchorChain& particles );
//void StampWisps( ScalarGrid& grid, const Vector& P, const float opacity );
//void StampBlurredWisps( ScalarGrid& grid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, int seed );
//void StampWisps( ScalarGrid& grid, ColorGrid& cgrid, const Vector& P, const float opacity, const Color& Cd );
//void StampBlurredWisps( ScalarGrid& grid, ColorGrid& cgrid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& Cd, int seed );
 
 
 
void stamp_noise( ScalarGrid& grid, NoiseMachine& noise, const Vector& center, const float radius, const float fade );

}
#endif
