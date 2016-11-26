#ifndef GEO_SHAPE_ISHAPE_H
#define GEO_SHAPE_ISHAPE_H

#include <Geo/Config.h>
#include <Geo/ShapeTypes.h>
#include <Geo/bv/BoundingVolume.h>
#include "IDomainSampler.h"

namespace geo {

/*! Shape Interface */
class IShape
{
public:
    //! \name Generic interface stuff
    //@{
    /*
    static const unsigned int cDimension = 0;
    typedef void transform_type;
    typedef void sdof_type;
    */
    //@}

public:
    IShape() {}
    virtual ~IShape() {}

//    virtual bool Create( const ShapeDef &shape_def ) { return false; }

    //! \name Dynamic interface
    //@{
    virtual EShapeType GetType() const = 0;
    virtual unsigned int GetNumSDOF() const = 0;
    virtual const Real* GetVecDefaultDOF() const = 0;

    // dimension-specific
    //void ComputeBV( bv::BoundingVolume &bv, const transform_type &transform, const sdof_type *vec_sdof ) const;
    //virtual const sdof_type *GetDefaultSDOF() const = 0;
    //@}

    /* \todo Offer Implicit/Parametric interfaces for each shape if suitable...
    virtual bool IsCurve() const;
    virtual bool IsSurface() const;
    virtual bool IsVolume() const;
    virtual const IImplicitCurve2 *GetImplicitCurve2() const {}
    virtual const IImplicitSurface3 *GetImplicitSurface3() const {}
    virtual const IParametricCurve2 *GetParametricCurve2() const {}
    virtual const IParametricSurface3 *GetParametricSurface3() const {}
    //same With surface, volume...
    */

    //\name Optional apis
    /*
    virtual IDomainSampler2 *CreateDomainSampler() const = 0;
    virtual const ISupportMap2 *GetSupportMap() const = 0;
    virtual const DistanceMap2 *GetDistanceMap() const = 0;
    //\todo virtual const IImplicit *GetImplicit() const = 0;
    //\todo MAYBE... probably not, virtual const ISAT *GetSAT() const = 0;
    */
};

template<unsigned D>
class GIShapeD: public IShape
{
public:
    static const unsigned int cDimension = 0;
    typedef mal::GTransform<Real,D> transform_type;

public:
    GIShapeD() {}
    ~GIShapeD() {}

    //\todo When DOF are generic and not SDOF, ComputeBV can be defined here
    //void ComputeBV( bv::GBoundingVolumeD<D> &bv, const mal::GTransform<Real,D> &transform, const DOF *vec_sdof ) const {};

    //\name Optional apis
    //@{
    virtual GIDomainSamplerD<D>* CreateDomainSampler() const = 0;
    //virtual const GISupportMap<D> *GetSupportMap() const = 0;
    //virtual const GDistanceMap<D> *GetDistanceMap() const = 0;
    //\todo virtual const GIImplicit<D> *GetImplicit() const = 0;
    //\todo MAYBE... probably not, virtual const ISAT *GetSAT() const = 0;
    //@}
};

typedef GIShapeD<2> IShape2;
typedef GIShapeD<3> IShape3;

#ifdef __DISABLED_TODO

//\todo These implicit/parametric curves/surfaces could be specified
//as GIFuncs?? or better not... there's no need to compose them, by
//now
// Implicit interface, Shapes that support it must return it when asked for
class IImplicitCurve2
{
    Real F( const &Vec2 point ) const;
    Vec2 dF_dp( const &Vec2 point ) const; //Gradient/Jacobian
    Mat2x2 d2F_dp2( const &Vec2 point ) const; //Hessian?
};

// Parametric interface
class IParametricCurve2
{
    Vec2 F( Real u ) const;
    Vec2 dF_du( Real u ) const;
    Vec2 d2F_du2( Real u ) const;
};
#endif

} //namespace geo

#endif // GEO_SHAPE_ISHAPE_H
