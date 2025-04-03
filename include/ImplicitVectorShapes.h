
#ifndef __IMPLICITVECTORSHAPES_H__
#define __IMPLICITVECTORSHAPES_H__

#include "Volume.h"
#include "LinearAlgebra.h"
#include <cmath>
#include <vector>
#include <iostream>

#include "Noise.h"
#include "NoiseMachine.h"
#include "VolumeGrid.h"
#include "SparseGrid.h"

namespace lux
{



class ConstantVectorVolume : public Volume<Vector> 
{
  public:

    ConstantVectorVolume( const Vector& v ) ;

   ~ConstantVectorVolume(){}


    const Vector eval( const Vector& P ) const; 

    const Matrix grad( const Vector& P ) const;

    virtual std::string typelabel() 
    { 
       std::string lbl = "Constant";
       return lbl;
    }


  private:

    Vector value;
    Matrix gradvalue;
};

class AddVectorVolume : public Volume<Vector> 
{
  public:

    AddVectorVolume( Volume<Vector> * v1, Volume<Vector> * v2 ) ;

    AddVectorVolume( const VectorField& v1, const VectorField& v2 ) ;

    ~AddVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


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
  
    VectorField elem1, elem2;
};


class GradientVectorVolume : public Volume<Vector> 
{
  public:

    GradientVectorVolume( Volume<float>* v, const float dx = 0.001 ) ;

    GradientVectorVolume( const ScalarField& v, const float dx = 0.001 ) ;

    ~GradientVectorVolume(){}

    const Vector eval( const Vector& P ) const;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Grad";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ScalarField elem;
    float step;
};


class ImplicitPointVectorVolume : public Volume<Vector> 
{
  public:
 
    ImplicitPointVectorVolume( Volume<float>* v, const float dx, const int nb = 3, const float thresh = 1.0e6 ) ;
 
    ImplicitPointVectorVolume( const ScalarField& v, const float dx, const int nb = 3, const float thresh = 1.0e6 ) ;
 
    ~ImplicitPointVectorVolume(){}
 
 
    const Vector eval( const Vector& P ) const;
 
    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
  //  const Matrix grad( const Vector& P ) const ;
 
 
    virtual std::string typelabel() 
    { 
       std::string lbl = "ImplicitPoint";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
 
  private:
 
    ScalarField elem;
    float step;
    int nbIterations;
    float threshold;
};

class AdvectVectorVolume : public Volume<Vector> 
{
  public:
 
    AdvectVectorVolume( Volume<Vector>* v, Volume<Vector>* u, const float delt ) ;
 
    AdvectVectorVolume( const VectorField& v, const VectorField& u, const float delt ) ;
 
   ~AdvectVectorVolume(){}
 
    const Vector eval( const Vector& P ) const ;
//    const Matrix grad( const Vector& P ) const ;
 
 
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
 
    VectorField elem, velocity;
    float dt;
};

class WarpVectorVolume : public Volume<Vector> 
{
  public:
 
    WarpVectorVolume( Volume<Vector>* v, Volume<Vector>* u ) ;
 
    WarpVectorVolume( const VectorField& v, const VectorField& u ) ;
 
   ~WarpVectorVolume(){}
 
    const Vector eval( const Vector& P ) const ;
//    const Matrix grad( const Vector& P ) const ;
 
 
    virtual std::string typelabel() 
    { 
       std::string lbl = "Warp";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + warp->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
 
 
  private:
 
    VectorField elem, warp;
};

class NoiseVectorVolume : public Volume<Vector> 
{
  public:
 
    NoiseVectorVolume( Noise* n, const float d = 0.01 ) ; 
    NoiseVectorVolume( NoiseMachine n, const float d = 0.01 ) ; 
 
   ~NoiseVectorVolume(){}
 
    const Vector eval( const Vector& P ) const;
//    const Matrix grad( const Vector& P ) const;  
 
 
    virtual std::string typelabel() 
    { 
       std::string lbl = "Noise";
       return lbl;
    }
 
  private:
 
    NoiseMachine noise;
    float dx;
};

class GriddedSGridVectorVolume : public Volume<Vector> 
{
  public:
 
    GriddedSGridVectorVolume( const VectorGrid& g ) ;
 
   ~GriddedSGridVectorVolume(){}
 
    const Vector eval( const Vector& P ) const ;
//    const Matrix grad( const Vector& P ) const ;
 
 
    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }
 
  private:
 
    VectorGrid elem;
    float dx, dy, dz;
};

}
#endif
