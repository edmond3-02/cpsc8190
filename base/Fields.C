
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "ImplicitMatrixShapes.h"
#include "MoreImplicitVolumes.h"
#include "ImplicitColors.h"
#include "Fields.h"
#include "VolumizedGrid.h"
#include "SparseGrid.h"

using namespace std;
namespace lux
{


const float  evaluate( const ScalarField& v, const Vector& P ) { return v->eval(P); }
const Vector evaluate( const VectorField& v, const Vector& P ) { return v->eval(P); }
const Color  evaluate( const ColorField&  v, const Vector& P ) { return v->eval(P); }
const Matrix evaluate( const MatrixField& v, const Vector& P ) { return v->eval(P); }

ScalarField ScalarField::operator+( const ScalarField& e2 ) { return add(*this,e2); }
ScalarField ScalarField::operator-() { return negate(*this); }
ScalarField ScalarField::operator-( const ScalarField& e2 ) { return add(*this,negate(e2)); }
ScalarField ScalarField::operator*( const ScalarField& e2 ) { return multiply(*this, e2); }
ScalarField ScalarField::operator*( const float e2 ) { return multiply(*this, e2); }


ColorField ColorField::operator+( const ColorField& e2 ) { return add(*this,e2); }
ColorField ColorField::operator*( const ScalarField& e2 ) { return multiply(*this, e2); }

VectorField VectorField::operator+( const VectorField& e2 ) { return add(*this,e2); }
VectorField VectorField::operator-( const VectorField& e2 ) { return subtract(*this,e2); }

MatrixField MatrixField::operator+( const MatrixField& e2 ) { return MatrixField( new AddMatrixVolume(*this,e2) ); }



ScalarField constant( const float v ) { return ScalarField( new ConstantVolume(v) ); }
VectorField constant( const Vector& v ) { return VectorField(new ConstantVectorVolume(v)); }
ColorField  constant( const Color& v ) { return ColorField(new ConstantColor(v)); }
MatrixField constant( const Matrix& v ) { return MatrixField( new ConstantMatrixVolume(v) ); }

ScalarField scale( const ScalarField& f, const Vector& s ) { return ScalarField( new ScaleVolume(f, s) ) ; }
//ColorField  scale( const ColorField& c, const Vector& s ){ return ColorField( new ScaleColorVolume(c, s)); }

ScalarField translate( const ScalarField& f, const Vector& s ) { return ScalarField( new TranslateVolume(f, s) ) ; }

ScalarField rotate( const ScalarField& f, const float a , const Vector& s ) { return ScalarField( new RotateVolume(f, a, s) ) ; }

ScalarField exp( const ScalarField& v ) { return ScalarField(new ExpVolume(v)); }
MatrixField exp( const MatrixField& m ) { return MatrixField( new ExpMatrixVolume( m ) ); }

ScalarField negate( const ScalarField& v ) { return ScalarField(new NegateVolume(v)); }

ScalarField abs( const ScalarField& v ) { return ScalarField(new AbsoluteVolume(v)); }

ScalarField which( const ScalarField& v1, const ScalarField& v2, const ScalarField& swtch ) { return ScalarField(new SwitchVolume(v1, v2, swtch)); }
ColorField  which( const ColorField& v1, const ColorField& v2, const ScalarField& swtch ) { return ColorField(new SwitchColorVolume(v1, v2, swtch)); }

ScalarField multiply( const ScalarField& v, const float a ) { return ScalarField( new MultiplyVolume(v, a) ); }
ScalarField multiply( const ScalarField& v, const ScalarField& u ) { return ScalarField(new MultiplyVolume(v, u)); }
ColorField  multiply( const ColorField& v, const ScalarField& u ){ return ColorField(new MultiplyColorScalar(v, u) ); }

ScalarField add( const ScalarField&  v1, const ScalarField& v2 ) { return ScalarField(new AddVolume(v1,v2));}
VectorField add( const VectorField&  v1, const VectorField& v2 ) { return VectorField(new AddVectorVolume(v1,v2)); }
ColorField  add( const ColorField&  v1, const ColorField& v2 )   { return ColorField(new AddColor(v1,v2)); }
MatrixField add( const MatrixField&  v1, const MatrixField& v2 ) { return MatrixField( new AddMatrixVolume( v1, v2 )); }

VectorField subtract( const VectorField&  v1, const VectorField& v2 ) { return VectorField(new SubtractVectorVolume(v1, v2)); }


ScalarField Sphere( const Vector& cen, const float rad ) { return ScalarField( new SphereVolume(cen, rad) ); }
ScalarField Ellipse( const Vector& cen, const Vector& axs, const float majorrad, const float minorrad ) 
		{ return ScalarField( new EllipseVolume(cen, axs, majorrad, minorrad) ); }
ScalarField CsgBox( const Vector& cen, const float rad, const float pwr ){ return ScalarField( new CsgBoxVolume(cen, rad, pwr) ); }
ScalarField Cone( const Vector& cen, const Vector& ax, const float h, const float theta ){ return ScalarField( new ConeVolume(cen, ax, h, theta) ); }
ScalarField Plane( const Vector cen, const Vector norm ) { return ScalarField( new PlaneVolume(cen, norm) ); }
ScalarField Torus( const Vector& cen, const Vector& axis, const float majorRad, const float minorRad ) 
		{ return ScalarField(new TorusVolume(cen, axis, majorRad, minorRad) ); }
ScalarField SteinerPatch() { return ScalarField(new SteinerPatchVolume() ) ; }
ScalarField Icosahedron() { return ScalarField( new IcosahedronVolume() ) ; }
ScalarField Cylinder( const Vector axis, const float rad ) { return ScalarField( new InfiniteCylinder(axis, rad)); }
ScalarField CappedCylinder( const Vector cen, const Vector axis, const float length, const float radius ) 
		{ return ScalarField( new ImplicitCylinder(cen, axis, length, radius)); }


ScalarField mask( const ScalarField& v ) { return ScalarField( new MaskVolume(v));}
ScalarField clamp( const ScalarField& v, float minv, float maxv ) { return ScalarField( new ClampVolume(v, minv, maxv)); }
ScalarField pow( const ScalarField& v, float gam ) { return ScalarField( new GammaVolume(v, gam)); }
ScalarField pow( const ScalarField& v, const ScalarField& gam ) { return ScalarField( new GammaVolume(v, gam)); }
ScalarField BlinnBlend( const ScalarField& v1, const ScalarField& v2, const float a1, const float a2 ) 
		{ return ScalarField( new BlinnBlendVolume(v1, v2, a1, a2)); }
ScalarField Union( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new UnionVolume( v1, v2) ) ; }
ScalarField intersection( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new IntersectionVolume( v1, v2) ) ; }
ScalarField cutout( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new CutoutVolume(v1, v2) ) ;}
ScalarField Shell( const ScalarField& v, const float thickness ) { return ScalarField( new ShellVolume(v, thickness) ) ; }

ScalarField Pyroclast(  const ScalarField e, const NoiseMachine n, const float Amp, const int i, const float Gamma )
		{ return ScalarField(new PyroclasticVolume(n, e, Amp, i, Gamma)); }
ScalarField RadialPyroclast( const Vector& Center, const float Radius, const float Amp, 
                              const float octaves, const float freq, const float rough, 
                              const float trans, const float time, const float Gamma ) 
		{ return ScalarField(new RadialPyroclasticVolume( Center, Radius, Amp, octaves, freq, rough, trans, time, Gamma)); }

ScalarField gridded(const ScalarGrid& g) { return ScalarField( new GriddedSGridVolume( g ) ) ; }
VectorField gridded( const VectorGrid& g ){ return VectorField( new GriddedSGridVectorVolume( g ) ) ; }

ScalarField advect( const ScalarField& v, const VectorField& u, const float delt ) { return ScalarField( new AdvectVolume(v, u, delt)); }
VectorField advect( const VectorField& v, const VectorField& u, const float delt ) { return VectorField( new AdvectVectorVolume(v, u, delt)); }

VectorField advectBFECC( const VectorField& v, const VectorField& u, const float delt, int nb ) 
		{ return VectorField( new BFECCAdvectVectorVolume( v, u, delt, nb)); }

ScalarField warp( const ScalarField& v, VectorField& map ) { return ScalarField( new WarpVolume(v, map)); }
VectorField warp( const VectorField& v, VectorField& map ) { return VectorField( new WarpVectorVolume(v, map)); }

VectorField identity() { return VectorField( new IdentityVectorVolume()); }

ScalarField SFNoise( NoiseMachine n, const float d ) { return ScalarField( new NoiseVolume(n, d)); } 
VectorField VFNoise( NoiseMachine n, const float d ) { return VectorField( new NoiseVectorVolume(n, d)); } 


ScalarField volumize(const ScalarGrid& g) { return ScalarField( new GriddedSGridVolume( g ) ) ; }

}
