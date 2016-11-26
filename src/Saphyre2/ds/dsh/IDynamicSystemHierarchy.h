#ifndef S2_DS_DSH_IDSH_H
#define S2_DS_DSH_IDSH_H

#include "IEntity.h"
#include <Saphyre2/ds/ss/ISimulationScheme.h> //Req by IBaseDSH implementation
#include <util/GPointerContainer.h>
#include <Geo/IObject.h>

// Fwd declarations from Geo
namespace geo { namespace bp { struct Proxy; struct Pair; }}
namespace geo { namespace bv { class IBoundingVolume; }}

namespace S2 { namespace ds {

// Fwd declarations
class IConnector;
class IConstraint;
class IInteraction;
class IGeom;

typedef util::GPointerContainer<IDynamicSystemHierarchy> ContainerDSH;
typedef util::GPointerContainer<IConnector> ContainerConnectors;
typedef util::GPointerContainer<IConstraint> ContainerConstraints;
typedef util::GPointerContainer<IGeom> ContainerGeoms;

/*! Generic DynamicSystem Hierarchy node interface:
  - Manages children DSH and Entities
  - Manages update through an optional ISimulationScheme
*/
class IDynamicSystemHierarchy: public IEntity
{
public:
    inline IDynamicSystemHierarchy( machine_uint_type uid, IDynamicSystemHierarchy* p_parent ) : IEntity(uid,p_parent) {}
    virtual ~IDynamicSystemHierarchy() {}

    virtual EDSHType GetDSHType() const = 0;
    virtual unsigned int GetDimension() const { return 0; }

    virtual bool IsLeaf() const = 0;
    virtual bool IsAggregate() const = 0;

    virtual void Update( Real dt ) = 0; //!< Any implementation must at least: Step(dt), SyncGO() and RecomputeBV()
    void Update_Local( Real dt ) { Step(dt); SyncGO(); RecomputeBV(); }

    //! \name DOF interface
    //@{
    virtual unsigned int GetNumDOF() const { return 0; }
    virtual void GetDOF( Real* p_dof ) const {}
    //@}

    //! \name Simulation Scheme interface
    //@{
    virtual bool BindSS( ISimulationScheme* p_ss ) = 0;
    virtual ISimulationScheme* GetSS() const = 0;
    //@}

    //! \name Sub-entity interface
    //@{
    virtual unsigned int GetNumChildren() const = 0;

    virtual bool AddChild( IDynamicSystemHierarchy* p_child ) = 0;
    virtual bool AddConnector( IConnector* p_connector ) = 0;
    virtual bool AddConstraint( IConstraint* p_constraint ) = 0;
    virtual bool AddGeom( IGeom* p_geom ) = 0;

    virtual bool RemoveChild( IDynamicSystemHierarchy* p_child ) = 0;
    virtual bool RemoveConnector( IConnector* p_connector ) = 0;
    virtual bool RemoveConstraint( IConstraint* p_constraint ) = 0;
    virtual bool RemoveGeom( IGeom* p_geom ) = 0;

    virtual ContainerDSH::iterator GetChildrenIterator() const = 0;
    virtual ContainerConnectors::iterator GetConnectorIterator() const = 0;
    virtual ContainerConstraints::iterator GetConstraintIterator() const = 0;
    virtual ContainerGeoms::iterator GetGeomIterator() const = 0;
    //@}

    //! \name Bound-geometry interface
    //@{
    virtual bool BindGO( geo::IObject* p_go ) { return false; }
    virtual geo::IObject* GetGO() = 0;
    virtual const geo::IObject* GetGO() const = 0;
    //@}

    //!\name Broad-Phase related interface, only used by SS
    //@{
    virtual const geo::bv::IBoundingVolume* GetBV() const = 0;
    virtual const geo::bp::Proxy* GetProxyBP() const = 0;
    virtual void SetBV( geo::bv::IBoundingVolume* p_bv ) = 0;
    virtual void SetProxyBP( const geo::bp::Proxy* p_proxy ) = 0;
    virtual void NotifyPairBP( const geo::bp::Proxy* p_other, const geo::bp::Pair* p_pair ) = 0;
    //@}

    //!\name IEntity implementation
    //@{
    EEntityType GetEntityType() const { return eEntity_DSH; }
    virtual bool Destroy_Internal( ReturnStream& rets, bool b_notify_parent ) = 0; //!< Forwarded to IBaseDSH as public...
    //@}

protected:
    //!\name Local update sub-phases
    //@{
    virtual void Step( Real dt ) = 0; //!< Integrate DOF forward during dt
    virtual void SyncGO() = 0;        //!< Sync GO to current DOF
    virtual void RecomputeBV() = 0;   //!< Compute BV from current GO or DOF

    virtual void Predict( Real dt ) {}
    virtual void Unpredict( Real dt ) {}
    //@}

    friend class SS_AggregateDSH_Basic; //TEMPORAL
};

/*! IDynamicSystemHierarchy basic partial implementation.
  Implements functionality and attributes common to all
  IDynamicSystemHierarchy subclasses:
  - Manages children DSH and Entities
  - Manages update through an optional ISimulationScheme
*/
class IBaseDSH: public IDynamicSystemHierarchy
{
public:
    inline IBaseDSH( machine_uint_type uid, IDynamicSystemHierarchy* p_parent )
    : IDynamicSystemHierarchy(uid,p_parent)
    , m_pSS(0)
    , m_pBV(0), m_pProxyBP(0) {}
    ~IBaseDSH() {}

    bool IsLeaf() const { return false; }
    bool IsAggregate() const { return false; }

    void Update( Real dt ) { if(m_pSS) m_pSS->Update(dt); else Update_Local(dt); }

    //! \name Simulation Scheme management
    //@{
    bool BindSS( ISimulationScheme* p_ss ) { m_pSS = p_ss; return m_pSS->BindDSH(this); }
    ISimulationScheme* GetSS() const { return m_pSS; }
    //@}

    //! \name Sub-entity interface implementation
    //@{
    unsigned int GetNumChildren() const { return m_ChildrenDSH.Size(); }

    bool AddChild( IDynamicSystemHierarchy* p_child )
    { m_ChildrenDSH.Add(p_child); if(m_pSS) m_pSS->NotifyAddChild(p_child); return true; }
    bool AddConnector( IConnector* p_connector )
    { m_Connectors.Add(p_connector); if(m_pSS) m_pSS->NotifyAddConnector(p_connector); return true; }
    bool AddConstraint( IConstraint* p_constraint )
    { m_Constraints.Add(p_constraint); if(m_pSS) m_pSS->NotifyAddConstraint(p_constraint); return true; }
    bool AddGeom( IGeom* p_geom )
    { m_Geoms.Add(p_geom); if(m_pSS) m_pSS->NotifyAddGeom(p_geom); return true; }

    bool RemoveChild( IDynamicSystemHierarchy* p_child )
    { m_ChildrenDSH.Remove(p_child); if(m_pSS) m_pSS->NotifyRemoveChild(p_child); return true; }
    bool RemoveConnector( IConnector* p_connector )
    { m_Connectors.Remove(p_connector); if(m_pSS) m_pSS->NotifyRemoveConnector(p_connector); return true; }
    bool RemoveConstraint( IConstraint* p_constraint )
    { m_Constraints.Remove(p_constraint); if(m_pSS) m_pSS->NotifyRemoveConstraint(p_constraint); return true; }
    bool RemoveGeom( IGeom* p_geom )
    { m_Geoms.Remove(p_geom); if(m_pSS) m_pSS->NotifyRemoveGeom(p_geom); return true; }

    ContainerDSH::iterator GetChildrenIterator() const { return m_ChildrenDSH.Begin(); }
    ContainerConnectors::iterator GetConnectorIterator() const { return m_Connectors.Begin(); }
    ContainerConstraints::iterator GetConstraintIterator() const { return m_Constraints.Begin(); }
    ContainerGeoms::iterator GetGeomIterator() const { return m_Geoms.Begin(); }
    //@}

    //!\name BP interface implementation
    //@{
    const geo::bv::IBoundingVolume* GetBV() const { return m_pBV; }
    const geo::bp::Proxy* GetProxyBP() const { return m_pProxyBP; }
    void SetBV( geo::bv::IBoundingVolume* p_bv ) { DS_ASSERT(0!=p_bv); m_pBV = p_bv; RecomputeBV(); }
    void SetProxyBP( const geo::bp::Proxy* p_proxy ) { m_pProxyBP = p_proxy; }
    void NotifyPairBP( const geo::bp::Proxy* p_other, const geo::bp::Pair* p_pair ) { if(m_pSS) m_pSS->NotifyPairBP(p_other,p_pair); }
    //@}

    //!\name IEntity implementation
    //@{
    //! Recursive by default (overriden by Leafs)
    virtual void QueryState( ReturnStream& rets ) const;
    virtual void DoViz( ReturnStream& rets ) const;
    virtual bool Destroy_Internal( ReturnStream& rets, bool b_notify_parent );
    //@}

protected:
    ISimulationScheme* m_pSS;
    ContainerDSH m_ChildrenDSH;
    ContainerConnectors m_Connectors;
    ContainerConstraints m_Constraints;
    ContainerGeoms m_Geoms;

    geo::bv::IBoundingVolume* m_pBV;
    const geo::bp::Proxy* m_pProxyBP;

    //ICoordinateTransform* m_pLocalToParentCT;
};

class ILeafDSH: public IBaseDSH
{
public:
    inline ILeafDSH( machine_uint_type uid, IDynamicSystemHierarchy* p_parent ) : IBaseDSH(uid,p_parent) {}
    ~ILeafDSH() {}

    bool IsLeaf() const { return true; }

    //! \name Sub-entity interface implementation
    //@{
    bool AddChild( IDynamicSystemHierarchy* p_child ) { return false; }
    bool RemoveChild( IDynamicSystemHierarchy* p_child ) { return false; }
    //@}
};


class AggregateDSH: public IBaseDSH
{
public:
    inline AggregateDSH( machine_uint_type uid, IDynamicSystemHierarchy* p_parent ) : IBaseDSH(uid,p_parent) {}
    ~AggregateDSH() {}

    bool IsAggregate() const { return true; }

    void Step( Real dt )
    {
        for( ContainerDSH::iterator it_dsh=m_ChildrenDSH.Begin(); it_dsh.IsValid(); ++it_dsh )
            it_dsh->Update(dt);
    }

    void SyncGO()
    {
        // \todo Children GO is already Sync at Update(), here we ONLY
        // need to Sync AGGREGATE GO (if any) to AGGREGATE DOF...
    }

    void RecomputeBV()
    {
        // \todo Compute Aggregate BV from sub-dsh and sub-geom GO or BV
        /*
        for( ContainerDSH::iterator it_dsh=m_ChildrenDSH.Begin(); it_dsh.IsValid(); ++it_dsh )
            m_pBV->Merge( it_dsh->GetBV() );
        for( ContainerGeoms::iterator it_geom=m_Geoms.Begin(); it_geom.IsValid(); ++it_geom )
            m_pBV->Merge( it_geom->GetBV() );
        */
    }
};


template <unsigned D, class BaseT>
class GIDynamicSystemHierarchyD: public BaseT
{
public:
    static const unsigned int cDimension = D;
    typedef BaseT base_type;
    typedef GIDynamicSystemHierarchyD<D,base_type> self_type;
    typedef geo::GIObjectD<D> go_type;
    typedef mal::GTransform<Real,D> transform_type;

public:
    inline GIDynamicSystemHierarchyD( machine_uint_type uid, IDynamicSystemHierarchy* p_parent )
    : base_type(uid,p_parent)
    , m_pGO(0) {}
    ~GIDynamicSystemHierarchyD() {}

    unsigned int GetDimension() const { return D; }

    //! \name Bound-geometry interface
    //@{
    bool BindGO( geo::IObject* p_go )
        {
            DS_ASSERT( cDimension == p_go->GetDimension() );
            m_pGO = static_cast<go_type*>(p_go);
            return true;
        }
    geo::IObject* GetGO() { return m_pGO; }
    const geo::IObject* GetGO() const { return m_pGO; }
    //@}

    //!\name Refinements may set/get transform
    //@{
    finline const transform_type& GetTransformRel() const { return m_TransformRel; }
    finline void SetTransformRel( const transform_type& transform ) { m_TransformRel = transform; }
    void ComputeTransformAbs( transform_type& transform ) const { transform = m_TransformRel; } //\todo Should compose with GetParent()->GetTransformAbs()
    //@}

protected:
    go_type* m_pGO;
    transform_type m_TransformRel;
};

//! \name Some common typedefs
//@{
typedef GIDynamicSystemHierarchyD<2,ILeafDSH> ILeafDSH2;
typedef GIDynamicSystemHierarchyD<3,ILeafDSH> ILeafDSH3;
typedef GIDynamicSystemHierarchyD<2,AggregateDSH> AggregateDSH2;
typedef GIDynamicSystemHierarchyD<3,AggregateDSH> AggregateDSH3;
//@}

}} // namespace S2::ds

#endif // S2_DS_DSH_IDSH_H
