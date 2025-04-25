

#ifndef __IMPLICITCOLORS_H__
#define __IMPLICITCOLORS_H__


#include "Volume.h"
#include "LinearAlgebra.h"

namespace lux
{

class ConstantColor : public Volume<Color> 
{
  public:

    ConstantColor( const Color& v );

   ~ConstantColor(){}


    const Color eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Constant";
       return lbl;
    }


  private:

    Color value;
};


class MultiplyColorScalar : public Volume<Color> 
{
  public:

    MultiplyColorScalar( Volume<Color> * e, Volume<float>* v );

    MultiplyColorScalar( const ColorField& e, const ScalarField& v );

   ~MultiplyColorScalar(){}


    const Color eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "MultiplyColorScalar";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem1;
    ScalarField elem2;
};






class AddColor : public Volume<Color> 
{
  public:

    AddColor( Volume<Color> * e, Volume<Color>* v );

    AddColor( const ColorField& e, const ColorField& v );

   ~AddColor(){}


    const Color eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Add";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem1, elem2;
};



class SwitchColorVolume : public Volume<Color>
 {
   public:
  
    SwitchColorVolume( Volume<Color>* v1, Volume<Color>* v2, Volume<float>* swtch );
  
    SwitchColorVolume( const ColorField& v1, const ColorField& v2, const ScalarField& swtch );
  
    ~SwitchColorVolume(){}
  
    const Color eval( const Vector& P ) const;
    
    //const Vector grad( const Vector& P ) const;
    
    virtual std::string typelabel() 
    { 
       std::string lbl = "Switch";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ",";
       lbl = lbl + swtchelem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
  
   private:
  
    const ColorField elem1;
    const ColorField elem2;
    const ScalarField swtchelem;
  
  
 };





}
#endif
