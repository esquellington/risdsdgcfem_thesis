#ifndef S2_DS_SS_SS_LEAFDSH_PARTICLE_SYSTEM_H
#define S2_DS_SS_SS_LEAFDSH_PARTICLE_SYSTEM_H

#include "ISimulationScheme.h"
#include <Saphyre2/ms/ParticleSystem.h>
#include <Geo/bv/GSphere.h>
#include <Geo/bp/IBroadPhase.h>

namespace S2 {
namespace ds {

// Fwd declarations
class IDynamicSystemHierarchy;
class IConnector;
class IConstraint;
class IInteraction;
class IGeom;

/*! Dimension-agnostic Leaf ParticleSystem Simulation Scheme.
  
IMPORTANT: The LeafDSH_ParticleSystem MUST be initialized (with
Create()) before binding to this SS.
  
*/
class SS_LeafDSH_ParticleSystem: public ISimulationScheme
{
public:
    SS_LeafDSH_ParticleSystem();
    ~SS_LeafDSH_ParticleSystem();
    
    bool BindDSH( IDynamicSystemHierarchy *p_dsh ); //!< Returns false if incompatible
    IDynamicSystemHierarchy *GetDSH() const { return m_pDSH; }

    void Update( Real dt );
    
    //! \name Bound DSH Notifications
    //@{
    void NotifyAddChild( IDynamicSystemHierarchy *p_child ) { DS_ASSERT(false); }
    void NotifyAddConnector( IConnector *p_connector ) {}
    void NotifyAddConstraint( IConstraint *p_constraint ) {}
    void NotifyAddGeom( IGeom *p_geom ) {}

    void NotifyRemoveChild( IDynamicSystemHierarchy *p_child ) { DS_ASSERT(false); }
    void NotifyRemoveConnector( IConnector *p_connector ) {}
    void NotifyRemoveConstraint( IConstraint *p_constraint ) {}
    void NotifyRemoveGeom( IGeom *p_geom ) {}
    void NotifyPairBP( const geo::bp::Proxy *p_other, const geo::bp::Pair *p_pair ) {}
    //@}

    //! \name BP access for external queries
    //@{
    const geo::bp::IBroadPhase *GetBP() const { return m_pBP; }
    //@}
    
    //! \name Debug
    //@{
    //virtual void DoLog( util::ItemStream &ls ) const {}//= 0;
    void DoViz( util::VizStream &vs ) const;
    //@}
    
private:
    void RecomputePerParticleBV();
    
private:
    IDynamicSystemHierarchy *m_pDSH;
    ms::IParticleSystem *m_pPS;
    union { ms::ParticleSystem2D *m_pPS2; ms::ParticleSystem3D *m_pPS3; };
    unsigned int m_Dimension;
    
    //!\name Per-particle annotations
    //@{
    const geo::bp::Proxy **m_vecProxy;
    union { geo::bv::Sphere2 *m_vecBV2; geo::bv::Sphere3 *m_vecBV3; };
    //@}

    //\name BroadPhase stuff
    //@{
    geo::bp::IBroadPhase *m_pBP;
    geo::bp::PairContainer m_PairsBP;
    //@}
    
    // specific ConstraintSets...
};

} } // namespace S2::ds

#endif // S2_DS_SS_SS_LEAFDSH_PARTICLE_SYSTEM_H
