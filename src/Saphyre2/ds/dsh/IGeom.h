#ifndef S2_DS_DSH_IGEOM_H
#define S2_DS_DSH_IGEOM_H

#include "IEntity.h"
#include <Geo/IObject.h>
//#include <Mal/GIFunction.h> //Not by now...

namespace geo { namespace bp { struct Proxy; }}

namespace S2 {
namespace ds {

class IDynamicSystemHierarchy;

/*! Geometric Entity Interface

  IGeom models the dependency of a local geo::IObject on the DOF of a
  DSH or a static/kinematic Boundary geometry.

  Dependency on a parent DSH Transform/DOF will be implicit in the
  IGeom's ICoordinateTransform functions.
  - Boundary: Rel Transform to parent DSH and Rel Shape DOF can be animated externally.
  - DSH: Rel transform to own DSH is fixed, and Rel Shape DOF are bound to own DSH

  An IGeom object can have a relative DOF transformation wrt the DSH
  RefSys or DOF:
  - If bound to a DSH: IGeom = RelDOF * DSH.DOF * RefSys
    - Rel/AbsDOF change only due to DSH.DOF changes
  - If Boundary: IGeom = RelDOF * RefSys
    - RelDOF/AbsDOF can be changed to animate the Boundary geometry

  Subclasses will know the exact number of DOF of the associated
  Shape and DSH, and will allocate it statically if possible.
  - ParticleGeom + RigidShape
  - RigidGeom + RigidShape
  - ParticleSysGeom +
*/
class IGeom: public IEntity
{
public:
    IGeom( machine_uint_type uid, IDynamicSystemHierarchy *parent ) : IEntity(uid,parent) {}
    ~IGeom() {}

    finline unsigned int GetDimension() const { return GetGO()->GetDimension(); }

    //virtual bool Bind( transform_bf_type *p_transform_bf, dof_bf_type *p_dof_bf ) = 0;

    virtual geo::IObject *GetGO() = 0;
    virtual const geo::IObject *GetGO() const = 0;

    //!\name Broad-Phase related interface
    /*\todo Dimensionless API is ugly, forces IXXX root classes,
      consider moving them ALL to GIGeom<D> and forcing mandatory
      dimension-cast on IGeom to access dimension-specific stuff (eg:
      BV)... try to use the SAME PATTERN EVERYWHERE
      (geo/saphyre/etc...), see geo::IObject for a starting point
    */
    //@{
    virtual const geo::bv::IBoundingVolume *GetBV() const = 0;
    virtual const geo::bp::Proxy *GetProxyBP() const = 0;
    virtual void SetBV( geo::bv::IBoundingVolume *p_bv ) = 0;      //!< Only called from SS
    virtual void SetProxyBP( const geo::bp::Proxy *p_proxy ) = 0; //!< Only called from SS
    //@}
};

/*! IGeom basic partial implementation.
  Implements functionality and attributes common to all IGeom
  subclasses.
*/
class IBaseGeom: public IGeom
{

public:
    IBaseGeom( machine_uint_type uid, IDynamicSystemHierarchy *parent ) : IGeom(uid,parent), m_pBV(0), m_pProxyBP(0) {}
    ~IBaseGeom() {}

    /*
    bool Bind( IGeom::transform_bf_type *p_transform_bf, IGeom::dof_bf_type *p_dof_bf )
    {
        m_pTransformBF = p_transform_bf;
        m_pDOFBF = p_dof_bf;
    }
    */

    //!\name BP interface implementation
    //@{
    const geo::bv::IBoundingVolume *GetBV() const { return m_pBV; }
    const geo::bp::Proxy *GetProxyBP() const { return m_pProxyBP; }
    void SetBV( geo::bv::IBoundingVolume *p_bv ) { m_pBV = p_bv; RecomputeBV(); }
    void SetProxyBP( const geo::bp::Proxy *p_proxy ) { m_pProxyBP = p_proxy; }
    //@}

    //!\name IEntity implementation
    //@{
    EEntityType GetEntityType() const { return eEntity_Geom; }
    bool Destroy_Internal( ReturnStream &rets, bool b_notify_parent );
    void DoViz( util::VizStream &vs ) const;
    //@}

protected:
    //!\name Internal methods
    //@{
    inline void RecomputeBV() { if( 0 != m_pBV ) GetGO()->ComputeBV(*m_pBV); }
    //@}

protected:
    geo::bv::IBoundingVolume *m_pBV;
    const geo::bp::Proxy *m_pProxyBP;

    /*
    transform_bf_type *m_pTransformBF; //Ex PSys3D, mal::GCopy<Transform3,Transform3>, Rigid: mal::TransformGBF3f
    dof_bf_type *m_pDOFBF; //Ex PSys3D: mal::GCopy<Vec3,Vec3> , Rigid = NULL
    */
};

/*! IGeom dimension-specific interface */
template <unsigned D>
class GGeomD: public IBaseGeom
{
public:
    static const unsigned int cDimension = D;
    typedef geo::GIObjectD<D> go_type;
    typedef mal::GTransform<Real,D> transform_type;

public:
    GGeomD( machine_uint_type uid, IDynamicSystemHierarchy *parent, go_type *p_go )
    : IBaseGeom(uid,parent)
    , m_pGO(p_go) {}

    ~GGeomD() { delete m_pGO; }

    //!\name Bound geo::IObject interface
    //@{
    geo::IObject *GetGO() { return m_pGO; }
    const geo::IObject *GetGO() const { return m_pGO; }
    //@}

    //!\name IEntity Implementation
    //@{
    bool Create( const ParamIt &pit )
    {
        m_pGO->SetTransform( pit.Find("transform").Get<transform_type>() );
        //\todo Set Shape DOF, if specified
        RecomputeBV();
        return true;
    }
    bool Edit( const ParamIt &pit )
    {
        //\todo PROPAGATE EDIT TO SHAPE if required??
        ParamIt pit2 = pit.Find("transform");
        if( pit2.IsValid() )
        {
            m_pGO->SetTransform( pit2.Get<transform_type>() );
            RecomputeBV();
        }
        //\todo Set Shape DOF, if specified
        return true;
    }
    void QueryState( ReturnStream &rets ) const
    {
        // \todo Should return any changes in Transform and DOF, but,
        // as only BS KineXXX can change it, it's already guaranteed
        // to be in sync. This may change if ds-side kinematic
        // animation is supported (controllers, etc...)
    }
    //@}

private:
    go_type *m_pGO;
    //transform_type m_TransformRel; //geo::IObject to Parent DSH
};

//! \name Usual dimensions
//@{
typedef GGeomD<2> Geom2;
typedef GGeomD<3> Geom3;
//@}

}} // namespace S2::ds

#endif // S2_DS_DSH_IGEOM_H
