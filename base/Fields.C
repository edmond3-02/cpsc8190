
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "ImplicitMatrixShapes.h"
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

VectorField VectorField::operator+( const VectorField& e2 ) { return add(*this,e2); }

MatrixField MatrixField::operator+( const MatrixField& e2 ) { return MatrixField( new AddMatrixVolume(*this,e2) ); }



ScalarField constant( const float v ) { return ScalarField( new ConstantVolume(v) ); }
VectorField constant( const Vector& v ) { return VectorField(new ConstantVectorVolume(v)); }
ColorField  constant( const Color& v ) { return ColorField(new ConstantColor(v)); }
MatrixField constant( const Matrix& v ) { return MatrixField( new ConstantMatrixVolume(v) ); }

ScalarField scale( const ScalarField& f, const Vector& s ) { return ScalarField( new ScaleVolume(f, s) ) ; }

ScalarField translate( const ScalarField& f, const Vector& s ) { return ScalarField( new TranslateVolume(f, s) ) ; }

ScalarField rotate( const ScalarField& f, const float a , const Vector& s ) { return ScalarField( new RotateVolume(f, a, s) ) ; }

ScalarField exp( const ScalarField& v ) { return ScalarField(new ExpVolume(v)); }
MatrixField exp( const MatrixField& m ) { return MatrixField( new ExpMatrixVolume( m ) ); }

ScalarField negate( const ScalarField& v ) { return ScalarField(new NegateVolume(v)); }

ScalarField abs( const ScalarField& v ) { return ScalarField(new AbsoluteVolume(v)); }

ScalarField multiply( const ScalarField& v, const float a ) { return ScalarField( new MultiplyVolume(v, a) ); }
ScalarField multiply( const ScalarField& v, const ScalarField& u ) { return ScalarField(new MultiplyVolume(v, u)); }
ColorField  multiply( const ColorField& v, const ScalarField& u ){ return ColorField(new MultiplyColorScalar(v, u) ); }

ScalarField add( const ScalarField&  v1, const ScalarField& v2 ) { return ScalarField(new AddVolume(v1,v2));}
VectorField add( const VectorField&  v1, const VectorField& v2 ) { return VectorField(new AddVectorVolume(v1,v2)); }
ColorField  add( const ColorField&  v1, const ColorField& v2 )   { return ColorField(new AddColor(v1,v2)); }
MatrixField add( const MatrixField&  v1, const MatrixField& v2 ) { return MatrixField( new AddMatrixVolume( v1, v2 )); }



ScalarField Sphere( const Vector& cen, const float rad ) { return ScalarField( new SphereVolume(cen, rad) ); }
ScalarField Ellipse( const Vector& cen, const Vector& axs, const float majorrad, const float minorrad ) { return ScalarField( new EllipseVolume(cen, axs, majorrad, minorrad) ); }
ScalarField CsgBox( const Vector& cen, const float rad, const float pwr ){ return ScalarField( new CsgBoxVolume(cen, rad, pwr) ); }
ScalarField Cone( const Vector& cen, const Vector& ax, const float h, const float theta ){ return ScalarField( new ConeVolume(cen, ax, h, theta) ); }
ScalarField Plane( const Vector cen, const Vector norm ) { return ScalarField( new PlaneVolume(cen, norm) ); }
ScalarField Torus( const Vector& cen, const Vector& axis, const float majorRad, const float minorRad ) { return ScalarField(new TorusVolume(cen, axis, majorRad, minorRad) ); }
ScalarField SteinerPatch() { return ScalarField(new SteinerPatchVolume() ) ; }
ScalarField Icosahedron() { return ScalarField( new IcosahedronVolume() ) ; }


ScalarField mask( const ScalarField& v ) { return ScalarField( new MaskVolume(v));}
ScalarField clamp( const ScalarField& v, float minv, float maxv ) { return ScalarField( new ClampVolume(v, minv, maxv)); }
ScalarField pow( const ScalarField& v, float gam ) { return ScalarField( new GammaVolume(v, gam)); }
ScalarField pow( const ScalarField& v, const ScalarField& gam ) { return ScalarField( new GammaVolume(v, gam)); }
ScalarField Union( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new UnionVolume( v1, v2) ) ; }
ScalarField intersection( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new IntersectionVolume( v1, v2) ) ; }
ScalarField cutout( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField( new CutoutVolume(v1, v2) ) ;}
ScalarField Shell( const ScalarField& v, const float thickness ) { return ScalarField( new ShellVolume(v, thickness) ) ; }

ScalarField SFNoise( NoiseMachine n, const float d ) { return ScalarField( new NoiseVolume(n, d)); } 

ScalarField volumize(const ScalarGrid& g) { return ScalarField( new GriddedSGridVolume( g ) ) ; }

}
