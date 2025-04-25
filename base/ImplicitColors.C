
#include "ImplicitColors.h"

using namespace lux;


ConstantColor::ConstantColor( const Color& v ) :
       value (v)
    {}
    
const Color ConstantColor::eval( const Vector& P ) const 
{
	return value; 
}


MultiplyColorScalar::MultiplyColorScalar( Volume<Color> * e, Volume<float>* v ) :
  elem1(e),
  elem2(v)
{}

MultiplyColorScalar::MultiplyColorScalar( const ColorField& e, const ScalarField& v ) :
  elem1(e),
  elem2(v)
{}


const Color MultiplyColorScalar::eval( const Vector& P ) const 
{
   return ( elem1->eval(P) * elem2->eval(P) ); 
}


AddColor::AddColor( Volume<Color> * e, Volume<Color>* v ) :
  elem1(e),
  elem2(v)
{}

AddColor::AddColor( const ColorField& e, const ColorField& v ) :
  elem1(e),
  elem2(v)
{}


const Color AddColor::eval( const Vector& P ) const 
{
   return ( elem1->eval(P) + elem2->eval(P) ); 
}


SwitchColorVolume::SwitchColorVolume( Volume<Color>* v1, Volume<Color>* v2, Volume<float>* swtch ):
	elem1(v1),
	elem2(v2),
	swtchelem(swtch)
      {}

SwitchColorVolume::SwitchColorVolume( const ColorField& v1, const ColorField& v2, const ScalarField& swtch ):
	elem1(v1),
	elem2(v2),
	swtchelem(swtch)
      {}

const Color SwitchColorVolume::eval( const Vector& P ) const
{
	Color e1 = elem1->eval(P);
	Color e2 = elem2->eval(P);
	return e1 + (e2 - e1) * swtchelem->eval(P);
}
