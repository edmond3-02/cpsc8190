
#ifndef __IMPLICITVOLUMESHAPES_H__
#define __IMPLICITVOLUMESHAPES_H__


#include <cmath>
#include <vector>
#include <iostream>
#include <string>


#include "Volume.h"
#include "LinearAlgebra.h"
#include "SparseGrid.h"

namespace lux
{


class ConstantVolume : public Volume<float> 
{
  public:

    ConstantVolume( const float v ) ;

   ~ConstantVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad( const Vector& P ) const; 

   virtual std::string typelabel() { return "Constant"; }


  private:

    float value;
    Vector gradvalue;
};

class ScaleVolume : public Volume<float> 
{
  public:

    ScaleVolume( const ScalarField& field, const Vector val ) ;

   ~ScaleVolume(){}


    const float eval( const Vector& P ) const; 

   virtual std::string typelabel() { return "Scale"; }


  private:
    ScalarField e;
    Vector v;
};

class TranslateVolume : public Volume<float> 
{
  public:

    TranslateVolume( const ScalarField& field, const Vector value) ;

   ~TranslateVolume(){}


    const float eval( const Vector& P ) const; 

   virtual std::string typelabel() { return "Translate"; }


  private:
    ScalarField e;
    Vector v;
};


class RotateVolume : public Volume<float> 
{
  public:

    RotateVolume( const ScalarField& field, const float angle , const Vector axis ) ;

   ~RotateVolume(){}


    const float eval( const Vector& P ) const; 

   virtual std::string typelabel() { return "Rotate"; }


  private:
    ScalarField e;
    float a;
    Vector a0;
};

class ExpVolume: public Volume<float> 
{
  public:

    ExpVolume( Volume<float>* v );

    ExpVolume( const ScalarField& v );

   ~ExpVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Exp";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;

};





class SphereVolume : public Volume<float> 
{
  public:

    SphereVolume( const Vector& cen, const float rad );

   ~SphereVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Sphere"); }


  private:

    Vector center;
    float radius;
};


class EllipseVolume : public Volume<float> 
{
  public:

    EllipseVolume( const Vector& cen, const Vector& axs, const float majorrad, const float minorrad );

   ~EllipseVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Ellipse"); }


  private:

    Vector center;
    Vector axis;
    float major_radius;
    float minor_radius;
};

class CsgBoxVolume : public Volume<float> 
{
  public:

    CsgBoxVolume( const Vector& cen, const float rad, const float pwr );

   ~CsgBoxVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("CsgBox"); }


  private:

    Vector center;
    float radius;
    float power;
};

class ConeVolume : public Volume<float> 
{
  public:

    ConeVolume( const Vector& cen, const Vector& axs, const float h, const float theta );

   ~ConeVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Cone"); }


  private:

    Vector center;
    Vector axis;
    float height;
    float angle;
};

class PlaneVolume : public Volume<float> 
{
  public:

    PlaneVolume( const Vector& cen, const Vector& norm );

   ~PlaneVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Plane"); }


  private:

    Vector center;
    Vector normal;
};

class TorusVolume : public Volume<float> 
{
  public:

    TorusVolume( const Vector& cen, const Vector& axs, const float majRad, const float minRad );

   ~TorusVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Torus"); }


  private:

    Vector center;
    Vector axis;
    float major_radius;
    float minor_radius;
};

class SteinerPatchVolume : public Volume<float> 
{
  public:

    SteinerPatchVolume();

   ~SteinerPatchVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("SteinerPatch"); }
 
};

class IcosahedronVolume : public Volume<float> 
{
  public:

    IcosahedronVolume();

   ~IcosahedronVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Icosahedron"); }
 
};

class MultiplyScalarVolume : public Volume<float> 
{
  public:

    MultiplyScalarVolume( const ScalarField&  e, const float v );

    ~MultiplyScalarVolume(){}


    const float eval( const Vector& P ) const;


   virtual std::string typelabel() 
   { 
      std::string lbl = "MultiplyScalar";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + std::to_string(val);
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1;
    const float val;
};




class AddVolume : public Volume<float> 
{
  public:

    AddVolume( Volume<float> * v1, Volume<float> * v2 );

    AddVolume( const ScalarField&  v1, const ScalarField& v2 );

    ~AddVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Add";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + e2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1, e2;
};


class MaskVolume : public Volume<float> 
{
  public:

    MaskVolume( Volume<float> * v1);

    MaskVolume( const ScalarField& v1);

    ~MaskVolume(){}


    const float eval( const Vector& P ) const;


   virtual std::string typelabel() 
   { 
      std::string lbl = "Mask";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1;
};


class UnionVolume : public Volume<float> 
{
  public:

    UnionVolume( const ScalarField&  v1, const ScalarField& v2 );

    ~UnionVolume(){}


    const float eval( const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Union";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + e2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1, e2;
};

class IntersectionVolume : public Volume<float> 
{
  public:
    
    IntersectionVolume( const ScalarField&  v1, const ScalarField& v2 );

    ~IntersectionVolume(){}


    const float eval( const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Intersection";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + e2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1, e2;
};

class CutoutVolume : public Volume<float> 
{
  public:
    
    CutoutVolume( const ScalarField&  v1, const ScalarField& v2 );

    ~CutoutVolume(){}


    const float eval( const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Cutout";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + e2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1, e2;
};

class ShellVolume : public Volume<float> 
{
  public:
    
    ShellVolume( const ScalarField&  v1, const float v2 );

    ~ShellVolume(){}


    const float eval( const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Shell";
      lbl = lbl + "(";
      lbl = lbl + e->typelabel();
      lbl = lbl + ",";
      lbl = lbl + std::to_string(d);
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e;
    const float  d;
};

/*
 *  PREFER VOLUMIZED GRID
 *  serious refactoring required, todo later
 */
class GriddedSGridVolume : public Volume<float> 
 {
   public:
  
     GriddedSGridVolume( const ScalarGrid& g );
  
    ~GriddedSGridVolume(){}
  
     const float eval( const Vector& P ) const;
     //const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }
  
   private:
  
     ScalarGrid scgrid;
     float dx, dy, dz;
 };


}




#endif



