#ifndef __RAYMARCHER_H__
#define __RAYMARCHER_H__
 
#include <vector>
#include <functional>
 
#include "Vector.h"
#include "Color.h"
#include "Volume.h"
#include "VolumeGrid.h"
//#include "AARectangle.h"
#include "SparseGrid.h"
//#include "PhaseFunction.h"
//#include "IntervalTree.h"
//#include "VolumeLight.h"
//#include "Image.h"
//#include "UniformPRN.h"
#include "Camera.h"
#include "ProgressMeter.h"
 
namespace lux
{
 
using namespace std;
 
 
class VShader;
 
class RenderData
{
  public:
   RenderData() :
      use_dsm_range      (false),
      dsm_range          (1.0),
      scatterCoefficient (Color(1,1,1,1)),
      snear		 (1),
      sfar		 (5),
      ds                 (0.01),
      maxPathlength      (100000000.0),
      Tmin		 (0.01),
      kappa		 (2.0),
      resolution	 {1920, 1080},
      channels		 (4)
   {}
 
   ~RenderData(){}
   int Nx() const {return resolution[0];}
   int Ny() const {return resolution[1];}
   int Nc() const {return channels;}
 
   ScalarField densityField;
   vector<ScalarField> dsmField;
   bool use_dsm_range;
   double dsm_range;
   vector<Color> lightColor;
   vector<Vector> lightPosition;
   ColorField colorField;  
   ColorField ambientColorField;  
   vector<Vector> startPosition;
   vector<Vector> startDirection;
   Color scatterCoefficient;

   float snear;
   float sfar;
   float ds;
   float maxPathlength;
   float Tmin;
   float kappa;
   
   int resolution[2];
   int channels;
   Camera camera;
 
};
 
void AddDSM( RenderData* d, const ScalarField& dsm );
void addLight( RenderData* d, Vector lightPosition, Color lightColor);
const ScalarField& GetDSM( RenderData* d, int i );
void SetDensityField(  RenderData* d, const ScalarField& field );
void SetColorField(  RenderData* d, const ColorField& field );
//void AddBoundingBox( RenderData *d, const Vector llc, const Vector urc );
//void AddBoundingBox( RenderData *d, const AABB& aabb );
//void AddBoundingBoxes( RenderData *d, const std::vector<AABB>& boxes );


void RenderFrame(const RenderData* d, ProgressMeter& pm, float* image, std::function<void(const RenderData*, const Vector&, Color&)> Render );

void RayMarchEmission(const RenderData* d, const Vector& direction, Color& L);

void RayMarchDSM(const RenderData* d, const Vector& direction, Color& L);
 
ScalarField RayMarchDSMAccumulation( 	const RenderData* d, 
					const ScalarField& densityField, 
					const Vector& lightPosition, 
					float ds,
					ScalarGrid& dsmField);

void RayMarchDSMAccumulation( const RenderData& input, 
                              const Vector& lightPosition, 
                              ScalarGrid& dsmField
                            );
 
 
void RayMarchVisibleDSMAccumulation( const RenderData& input, 
                                     const Vector& lightPosition, 
                                     const Camera& cam,
                                     ScalarGrid& dsmField
                                   );
 
 
void RayMarchDSMAccumulation( const ScalarField& densityField, 
                              const Vector& lightPosition, 
                              float ds,
                              ScalarGrid& dsmField,
                              ScalarField& holdout
                            );
 
 
 
 
void ssRayMarchAccumulation( const RenderData& input, vector<Color>& output_color, vector<Color>& output_opacity );

const Color gatherLight( const RenderData& data, const Vector P, const Vector Dir );
 
 
}
#endif
