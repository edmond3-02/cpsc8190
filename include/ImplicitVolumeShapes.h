
#ifndef __IMPLICITVOLUMESHAPES_H__
#define __IMPLICITVOLUMESHAPES_H__


#include <cmath>
#include <vector>
#include <iostream>
#include <string>


#include "Volume.h"
#include "LinearAlgebra.h"
#include "SparseGrid.h"
#include "NoiseMachine.h"

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

class NegateVolume : public Volume<float> 
 {
   public:
  
     NegateVolume( Volume<float> * v );
  
     NegateVolume( const ScalarField& v );
  
     ~NegateVolume(){}
  
  
     const float eval( const Vector& P ) const;
  
     //const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Negate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
   private:
   
     const ScalarField elem;
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

class AdvectVolume : public Volume<float> 
{
  public:
 
    AdvectVolume( Volume<float>* v, Volume<Vector>* u, const float delt ); 
 
    AdvectVolume( const ScalarField& v, const VectorField& u, const float delt ); 
 
   ~AdvectVolume(){}
 
    const float eval( const Vector& P ) const; 
    const Vector grad(  const Vector& P ) const;
 
    virtual std::string typelabel() 
    { 
       std::string lbl = "Advect";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + velocity->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
 
  private:
 
    const ScalarField elem;
    VectorField velocity;
    float dt;
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

class AbsoluteVolume : public Volume<float> 
{
  public:
 
    AbsoluteVolume( Volume<float> * v )  ; 
 
    AbsoluteVolume( const ScalarField& v )  ; 
 
    ~AbsoluteVolume(){}
 
 
    const float eval( const Vector& P ) const;
 
   //const Vector grad(  const Vector& P ) const;
 
   virtual std::string typelabel() 
   { 
      std::string lbl = "Absolute";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }
 
  private:
 
    const ScalarField elem;
 
};

class MultiplyVolume : public Volume<float> 
{
  public:

    MultiplyVolume( const ScalarField&  e, const float v );
    
    MultiplyVolume( const ScalarField& v, const ScalarField& u );

    ~MultiplyVolume(){}


    const float eval( const Vector& P ) const;


   virtual std::string typelabel() 
   { 
      std::string lbl = "Multiply";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ",";
      lbl = lbl + std::to_string(constant);
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;
    const ScalarField factor;
    float constant;
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

class ClampVolume : public Volume<float> 
 {
   public:
  
     //ClampVolume( Volume<float> * v, float minv, float maxv );
  
     ClampVolume( const ScalarField& v, float minv, float maxv );
  
  
     ~ClampVolume(){}
  
  
     const float eval( const Vector& P ) const;
  
     //const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Clamp";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
   private:
  
    const ScalarField elem;
     float vmin, vmax;
 };

class GammaVolume : public Volume<float> 
 {
   public:
  
     GammaVolume( Volume<float> * v, float gam );
  
     GammaVolume( const ScalarField& v, float gam );
  
     GammaVolume( const ScalarField& v, const ScalarField& g );
  
  
     ~GammaVolume(){}
  
  
     const float eval( const Vector& P ) const;
  
     //const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Gamma";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
   private:
  
     const ScalarField elem;
     const ScalarField factor;
     float gamma;
 };

class BlinnBlendVolume : public Volume<float> 
 {
   public:
  
     BlinnBlendVolume( Volume<float> * v1, Volume<float> * v2, const float _alpha = 1.0 );
  
     BlinnBlendVolume( const ScalarField& v1, const ScalarField& v2, const float _alpha = 1.0 );
  
  
     ~BlinnBlendVolume(){}
  
  
     const float eval( const Vector& P ) const;
 
     void setAlpha(const float a) { alpha = a; }
 
     //const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "BlinnBlend";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
  
   protected:
  
     const ScalarField elem1, elem2;
     float alpha;
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

class PyroclasticVolume : public Volume<float> 
 {
   public:
  
	PyroclasticVolume( const NoiseMachine n, const ScalarField e, const float Amp, const int i, const float Gamma = 1.0/3.0  );
  
    ~PyroclasticVolume(){}
  
	const float eval( const Vector& P ) const;
	//const Vector grad(  const Vector& P ) const;
  
    virtual std::string typelabel() 
    { 
	  std::string lbl = "Pyroclastic";
	  return lbl;
    }
  
   private:
  
	NoiseMachine noise;
	float amplitude;
	float gamma;
	int  iter;
	const ScalarField elem;
  
	float dx;
};

class RadialPyroclasticVolume : public Volume<float> 
{
  public:
 
    RadialPyroclasticVolume( const Vector& Center, const float Radius, const float Amp, 
                       const float octaves, const float freq, const float rough, const float trans, const float time, const float Gamma = 1.0/3.0 );
 
   ~RadialPyroclasticVolume(){}
 
    const float eval( const Vector& P ) const;
    //const Vector grad(  const Vector& P ) const;
 
   virtual std::string typelabel() 
   { 
      std::string lbl = "RadialPyroclastic";
      return lbl;
   }
 
  private:
 
    mutable FractalSum<PerlinNoiseGustavson> noise;
    float amplitude;
    Vector center;
    float gamma;
    float translate;
    float radius;
    mutable Noise_t parms;
 
    float dx;
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


class NoiseVolume : public Volume<float> 
 {
   public:
  
     NoiseVolume( Noise* n, const float d = 0.01 ); 
     NoiseVolume( NoiseMachine n, const float d = 0.01 ); 
  
    ~NoiseVolume(){}
  
     const float eval( const Vector& P ) const; 
     //const Vector grad(  const Vector& P ) const;  
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Noise";
       return lbl;
    }
  
   private:
  
     NoiseMachine noise;
     float dx;
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


class WarpVolume : public Volume<float> 
{
  public:
 
    WarpVolume( Volume<float>* v, Volume<Vector>* map );
 
    WarpVolume( const ScalarField& v, VectorField& map );
 
   ~WarpVolume(){}
 
    const float eval( const Vector& P ) const; 
    const Vector grad(  const Vector& P ) const;
 
    virtual std::string typelabel() 
    { 
       std::string lbl = "Warp";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + mapX->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
 
  private:
 
    const ScalarField elem;
    VectorField mapX;
};

}




#endif



