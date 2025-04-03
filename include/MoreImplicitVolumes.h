#ifndef __MOREIMPLICITVOLUMES_H__
 #define __MOREIMPLICITVOLUMES_H__
  
// #include "SpaceCurve.h"
 #include "Fields.h"
 #include "ImplicitVolumeShapes.h"
 #include "ImplicitVectorShapes.h"
// #include "GridVolumes.h"
// #include "AARectangle.h"
// #include "HermiteCurve.h"
 #include <cmath>
  
 namespace lux
 {
  
  
 class ImplicitFunctionDisplacement : public Volume<float>
 {
   public:
  
     ImplicitFunctionDisplacement( Volume<float>* f, Volume<float>* displacer, const float dx, const int nbsteps = 3 ) :
       elem(f),
       disp(displacer),
       ipvv( new ImplicitPointVectorVolume( elem, dx, nbsteps ) )
     {}
  
     ImplicitFunctionDisplacement( const ScalarField& f, const ScalarField& displacer, const float dx, const int nbsteps = 3 ) :
       elem(f),
       disp(displacer),
       ipvv( new ImplicitPointVectorVolume( elem, dx, nbsteps ) )
     {}
  
    ~ImplicitFunctionDisplacement(){}
  
     const float eval( const Vector& P ) const
     {
        float val0 =  elem->eval(P);
        Vector X =  ipvv->eval(P);
        float val1 =  disp->eval(X);
        return (val0 + val1);
     }
  
     const Vector grad( const Vector& P ) const
     {
        Vector val0 =  elem->grad(P);
        Vector X =  ipvv->eval(P);
        Vector val1 =  disp->grad(X);
        Matrix G = ipvv->grad(P);
        return val0 + G*val1;
     }
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "ImplicitFunctionDisplacement";
        lbl = lbl + "(";
        lbl = lbl + elem->typelabel();
        lbl = lbl + ",";
        lbl = lbl + disp->typelabel();
        lbl = lbl + ")";
        return lbl;
     }
  
   private:
  
     ScalarField elem, disp;
     VectorField ipvv;
  
 };
  
  
 class InfiniteCylinder : public Volume<float>
 {
   public:
  
     InfiniteCylinder( const Vector axis, const float rad ) :
        uaxis   (axis.unitvector()),
        radius  (rad)
     {}
  
   ~InfiniteCylinder(){}
  
  
    const float eval( const Vector& P ) const
    {
       Vector X = P - (P*uaxis)*uaxis;
       float rho = X.magnitude();
       return (radius - rho);
    }
  
    const Vector grad( const Vector& P ) const
    {
       Vector X = P - (P*uaxis)*uaxis;
       return X.unitvector();
    }
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "InfiniteCylinder";
        return lbl;
     }
  
  
   private:
  
     Vector uaxis;
     float radius;
  
  
 };
  
  
 class ImplicitCylinder : public Volume<float>
 {
   public:
  
     ImplicitCylinder( const Vector cen, const Vector axis, const float length, const float radius ) :
        elem ( new IntersectionVolume( 
                         ScalarField( new IntersectionVolume( 
                                ScalarField( ScalarField( new TranslateVolume( ScalarField(new InfiniteCylinder( axis.unitvector(), radius)), cen)) ), 
                                ScalarField( new PlaneVolume( cen+axis.unitvector()*length/2.0, -axis.unitvector() ) )  
                       ) 
                      ) ,   
                         ScalarField( new PlaneVolume( cen-axis.unitvector()*length/2.0, axis.unitvector()  ) )    
         )       
     )
     {}
  
   ~ImplicitCylinder(){}
  
  
    const float eval( const Vector& P ) const
    {
       return elem->eval(P);
    }
  
    const Vector grad( const Vector& P ) const
    {
       return elem->grad(P);
    }
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "ImplicitCylinder";
        return lbl;
     }
  
   private:
  
     ScalarField elem;
 };
  
  
  
  
 class ImplicitShell : public Volume<float>
 {
   public:
  
     ImplicitShell( Volume<float>* v, const float thickness ) :
       elem ( new IntersectionVolume(  
                            ScalarField( new AddVolume(
                                 ScalarField(v), 
                         ScalarField( new ConstantVolume(thickness*0.5)) 
                      ) 
                       ), 
                    ScalarField( new NegateVolume( ScalarField( new AddVolume( 
                                                                    ScalarField(v), 
                                            ScalarField(new ConstantVolume(-thickness*0.5) ) 
                                          ) 
                                          ) 
                                  ) 
                     ) 
              ) 
         ) 
     {}
  
     ImplicitShell( const ScalarField& v, const float thickness ) :
       elem ( new IntersectionVolume(  
                            ScalarField( new AddVolume(
                                 v, 
                         ScalarField( new ConstantVolume(thickness*0.5)) 
                      ) 
                       ), 
                    ScalarField( new NegateVolume( ScalarField( new AddVolume( 
                                                                    v, 
                                            ScalarField(new ConstantVolume(-thickness*0.5) ) 
                                          ) 
                                          ) 
                                  ) 
                     ) 
              ) 
         ) 
     {}
  
   ~ImplicitShell(){}
  
  
    const float eval( const Vector& P ) const
    {
       return elem->eval(P);
    }
  
    const Vector grad( const Vector& P ) const
    {
       return elem->grad(P);
    }
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "ImplicitShell";
        lbl = lbl + "(";
        lbl = lbl + elem->typelabel();
        lbl = lbl + ")";
        return lbl;
     }
  
   private:
     ScalarField elem;
 };
  
  
  
  
  
  
  
 class ProceduralDSM : public Volume<float>
 {
   public:
  
     ProceduralDSM( Volume<float>* density, const Vector lghtP, const float atten, const float step ) :
       elem(density),
       lightP (lghtP),
       attenuation (atten),
       ds (step)
     {}
  
    ~ProceduralDSM(){}
  
  
     const float eval( const Vector& P ) const
     {
        Vector D = (P-lightP);
        float distance = D.magnitude();
        D /= distance;
        float accum = 0;
        float s = 0;
        Vector X = P;
        D *= ds;
        while( s < distance )
        {
           accum += elem->eval(X);
       X += D;
       s += ds;
        }
        return accum * ds * attenuation;
     }
  
     const Vector grad( const Vector& P ) const { return Vector(0,0,0); }  // need to fix this
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "ProceduralDSM";
        lbl = lbl + "(";
        lbl = lbl + elem->typelabel();
        lbl = lbl + ")";
        return lbl;
     }
  
  
   private:
  
     ScalarField elem;
     Vector lightP;
     float attenuation;
     float ds;
  
 };
  
  
 class Bretzel2Volume : public Volume<float>
 {
  
   public:
  
     Bretzel2Volume( const float c = 1.0 ) : cutoff(c) {}
    ~Bretzel2Volume(){}
  
  
     const float eval( const Vector& P ) const
     {
        return cutoff - std::pow(  (double)((1.0-P[0]*P[0])*P[0]*P[0] - P[1]*P[1] ) ,(int)2 ) - P[2]*P[2];
     }
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "Bretzel2";
        return lbl;
     }
  
  
   private:
  
     float cutoff;
 };
  
  
 class CumuloVolume : public Volume<float>
 {
   public:
  
     CumuloVolume( Volume<float>* basefunction, vector<Volume<float>* >& dispArray, float step, int nbIterations ) 
     {
        cumuloed = ScalarField(basefunction);
        for( size_t i=0;i<dispArray.size();i++ )
        {
           VectorField lsX ( new ImplicitPointVectorVolume( cumuloed, step, nbIterations ) );
       ScalarField mappeddisp ( new WarpVolume( ScalarField(dispArray[i]), lsX ) );
           cumuloed =  cumuloed + mappeddisp;
        }
     }
  
     ~CumuloVolume(){}
  
     const float eval( const Vector& P ) const
     {
        return cumuloed->eval(P); 
     }
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "Cumulo";
        lbl = lbl + "(";
        lbl = lbl + cumuloed->typelabel();
        lbl = lbl + ")";
        return lbl;
     }
  
  
   protected:
  
     ScalarField cumuloed;
 };
  
  
 /* 
 class SpaceCurve;
  
 class SpaceCurveVolume : public Volume<float>
 {
   public:
  
     SpaceCurveVolume( const SpaceCurve& c, const float rad );
    ~SpaceCurveVolume();
  
     const float eval( const Vector& P ) const;
     const Vector grad( const Vector& P ) const;
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "SpaceCurve";
        return lbl;
     }
  
  
   protected:
  
     SpaceCurve curve;
     float radius;
 };
 */ 
  
 class PiecewiseCurveVolume : public Volume<float>
 {
   public:
  
     PiecewiseCurveVolume( const AnchorChain& c );
    ~PiecewiseCurveVolume();
  
     virtual const float eval( const Vector& P ) const;
     virtual const Vector grad( const Vector& P ) const;
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "PiecewiseCurve";
        return lbl;
     }
  
  
   protected:
  
     AnchorChain anchorList;
  
     void findClosestPoint( const Vector& P, int& anchor, double& position ) const;
     void interpolateAnchor( const int anchor, const double position, Anchor& value ) const;
     double distanceToSegment( const Vector& P, const size_t seg ) const;
     double positionOnSegment( const Vector& P, const size_t seg ) const;
 };
  
  
 class PiecewisePyroCurveVolume : public PiecewiseCurveVolume
 {
   public:
  
     PiecewisePyroCurveVolume( const AnchorChain& c );
    ~PiecewisePyroCurveVolume();
  
     const float eval( const Vector& P ) const;
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "PiecewisePyroCurve";
        return lbl;
     }
  
 };
  
  
 class PiecewiseNoiseCurveVolume : public PiecewiseCurveVolume
 {
   public:
  
     PiecewiseNoiseCurveVolume( const AnchorChain& c );
    ~PiecewiseNoiseCurveVolume();
  
     const float eval( const Vector& P ) const;
     const float fade( const float t ) const;
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "PiecewiseNoiseCurve";
        return lbl;
     }
 };
  
  
  
  
 /*
 class GriddedCumuloVolume : public Volume<float>
 {
   public:
  
     GriddedCumuloVolume( const Volume<float>* basefunction, const vector<Volume<float>* >& dispArray, float step, int nbIterations, int nx, int ny, int nz, const Vector& LLL, const Vector& origin ) 
     {
        addElement( basefunction );
        for( size_t i=0;i<dispArray.size();i++ ){ addElement( dispArray[i] ); }
        dx = step;
        iterations = nbIterations;
  
  
        const Volume<float>* cumuloed = basefunction;
        for( size_t i=0;i<dispArray.size();i++ )
        {
           Volume<Vector>* lsX =  new ImplicitPointVectorVolume( cumuloed, dx, iterations );
       // sample lsX onto a sparse grid
       VolumeGrid<Vector>* lsXgrid = new VolumeGrid<Vector>();
       lsXgrid->init( nx, ny, nz, LLL[0], LLL[1], LLL[2], origin );
       Sample( lsXgrid, lsX );
       delete lsX;
           lsX = new GriddedVectorVolume( lsXgrid );
       Volume<float>* mappeddisp = new WarpVolume( dispArray[i], lsX );
           cumuloed = new AddVolume( cumuloed, mappeddisp );
        }
        addElement( cumuloed );
     }
  
     ~GriddedCumuloVolume(){}
  
     const float eval( const Vector& P ) const
     {
        return getFloatElement( element.size()-1 )->eval(P); 
     }
  
  
  
   protected:
  
     float dx;
     int iterations;
  
 };
 */
  
 /* 
 class CompoundHermiteCurveVolume : public Volume<float>
 {
   public:
  
     CompoundHermiteCurveVolume( const lux::CompoundHermiteCurve& c, const float R );
    ~CompoundHermiteCurveVolume();
  
     const float eval( const Vector& P ) const;
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "CompoundHermiteCurve";
        return lbl;
     }
  
   private:
  
     lux::CompoundHermiteCurve elem;
     float radius;
 };
  
  
  */
  
  /*
  
 class BoxedVolume : public Volume<float>, AARectangle
 {
   public:
  
     BoxedVolume( const ScalarField& f, const Vector& llc, const Vector& urc, const float def );
    ~BoxedVolume(){}
  
     const float eval( const Vector& P ) const;
  
  
     virtual std::string typelabel() 
     { 
        std::string lbl = "BoxedVolume";
        lbl = lbl + "(";
        lbl = lbl + elem->typelabel();
        lbl = lbl + ")";
        return lbl;
     }
  
  
   private:
  
     const ScalarField& elem;
     const float defValue;
 };
  
  
  */ 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 }
  
 #endif
