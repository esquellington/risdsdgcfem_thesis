#ifndef S2_DS_SS_AGGREGATE_DSH_BASIC_H
#define S2_DS_SS_AGGREGATE_DSH_BASIC_H

#include "ISimulationScheme.h"
#include <Geo/bp/IBroadPhase.h>
#include <Saphyre2/ms/IConstraint.h>
#include <util/GPointerContainer.h>

#define __ENABLE_MULTI_SOLID3D
#ifdef __ENABLE_MULTI_SOLID3D
#  include <Geo/mp/TestContact.h>
#  include <util/GPoolDA.h>
#endif

namespace S2 { namespace ds {

/*! Simplest SS:
    - Children are updated locally
    - Collisions are detected only between direct children, and
      dispatched to child solvers.
    - Constraints:
      - User constraints are stored in the DSH
      - SS-created constraints are stored locally:
        - Contacts
        - "Translated" user-constraints notified from the DSH
*/
class SS_AggregateDSH_Basic: public ISimulationScheme
{
public:
    SS_AggregateDSH_Basic();
    ~SS_AggregateDSH_Basic();

    bool BindDSH( IDynamicSystemHierarchy* p_dsh );
    IDynamicSystemHierarchy* GetDSH() const { return m_pDSH; }

    void Update( Real dt );

    //! \name Bound DSH Notifications
    //@{
    void NotifyAddChild( IDynamicSystemHierarchy* p_child );
    void NotifyAddConnector( IConnector* p_connector ) {}
    void NotifyAddConstraint( IConstraint* p_constraint ) {}
    void NotifyAddGeom( IGeom* p_geom );
    void NotifyRemoveChild( IDynamicSystemHierarchy* p_child );
    void NotifyRemoveConnector( IConnector* p_connector ) {}
    void NotifyRemoveConstraint( IConstraint* p_constraint ) {}
    void NotifyRemoveGeom( IGeom* p_geom );
    void NotifyPairBP( const geo::bp::Proxy* p_other, const geo::bp::Pair* p_pair ) {}
    //@}

    //! \name BP access for external queries
    //@{
    const geo::bp::IBroadPhase* GetBP() const { return m_pBP; }
    //@}

    //! \name Debug
    //@{
    //virtual void DoLog( util::ItemStream& ls ) const {}//= 0;
    void DoViz( util::VizStream& vs ) const;
    //@}

private:
    IDynamicSystemHierarchy* m_pDSH;

    //\name BroadPhase stuff
    //@{
    geo::bp::IBroadPhase* m_pBP;
    geo::bp::PairContainer m_PairsBP;
    //@}

    //\name SS-created Constraints
    //@{
    typedef util::GPointerContainer<ms::IConstraint> ContainerConstraintImpl;
    ContainerConstraintImpl m_ConstraintsImpl;
    //@}

#ifdef __ENABLE_MULTI_SOLID3D
private:
    void MultiDSH_Solid3D_Update( Real dt );
    void MultiDSH_Solid3D_NotifyPairBP( const geo::bp::Pair* p_pair );
    void MultiDSH_Solid3D_ApplyContacts_PositionAlteration_Smoothed();
    void MultiDSH_Solid3D_ApplyContacts_Impulse();
    void MultiDSH_Solid3D_ApplyContacts_UpdateActiveSet();
    void MultiDSH_Solid3D_ApplyContacts_NormalImpulse();
    void MultiDSH_Solid3D_ApplyContacts_FrictionImpulse();

private:
    Real m_AccTime;
    //TEMP Modelled after ContactConstraintD2K in LeafDSH_Solid3D_FEM, see comments therein for further discussion
    struct ContactPointConstraintD2D
    {
        ContactPointConstraintD2D( const geo::PointOnFeature& pof1, const geo::PointOnFeature& pof2, const Vec3& normal, Real depth, Real radius )
        : m_POF1(pof1)
        , m_POF2(pof2)
        , m_Normal(normal)
        , m_Depth(depth)
        , m_Radius(radius)
        , m_LambdaN(0)
        , m_IsActive(true)
            {}
        geo::PointOnFeature m_POF1;
        geo::PointOnFeature m_POF2;
        Vec3 m_Normal;
        Real m_Depth;
        Real m_Radius;
        Real m_LambdaN; //Normal impulse
        bool m_IsActive;
    };
    struct ContactConstraintD2D
    {
        inline void Create() { Destroy(); m_TmpCD.Clear(); }
        inline void Destroy() { m_pDSH1 = 0; m_pDSH2 = 0; m_Age = 0; m_CC.Invalidate(); Clear(); }
        inline void Clear() { m_vecCPC.clear(); }
        inline void Reset( IDynamicSystemHierarchy* p_dsh1, IDynamicSystemHierarchy* p_dsh2,
                           const geo::np::ContactData3& cd, Real epsilon_length )
            {
                //\todo Use epsilon_length to merge close points, ONLY IF NOT already done by default with cd.Optimize(epsilon_length)
                DS_ASSERT( p_dsh1 != 0 && p_dsh2 != 0 );
                DS_ASSERT( cd.HasPOF() );
                m_pDSH1 = p_dsh1;
                m_pDSH2 = p_dsh2;
                m_Age = 0;
                m_vecCPC.clear();
                for( unsigned int it_cp=0; it_cp<cd.Size(); it_cp++ )
                    m_vecCPC.push_back( ContactPointConstraintD2D( cd.GetPOF1(it_cp),
                                                                   cd.GetPOF2(it_cp),
                                                                   cd.GetCP(it_cp).m_Normal,
                                                                   cd.GetCP(it_cp).m_Depth,
                                                                   cd.GetCP(it_cp).m_Radius ) );
            }
        inline void Persist( const geo::np::ContactData3& cd, Real epsilon_length, Real epsilon_direction )
            {
                // \todo Match and persist...
                m_Age++;
                Reset(m_pDSH1,m_pDSH2,cd,epsilon_length);
            }

        finline bool IsEmpty() const { return m_vecCPC.empty(); }
        finline unsigned int GetNumCPC() const { return m_vecCPC.size(); }
        finline ContactPointConstraintD2D& GetCPC( int i ) { return m_vecCPC[i]; }
        finline const ContactPointConstraintD2D& GetCPC( int i ) const { return m_vecCPC[i]; }

        IDynamicSystemHierarchy* m_pDSH1;
        IDynamicSystemHierarchy* m_pDSH2;
        uint32 m_Age;
        geo::np::ContactCache3 m_CC;
        std::vector<ContactPointConstraintD2D> m_vecCPC;

        geo::np::ContactData3 m_TmpCD; //TEMPORAL: Single CD stored to Viz later, otherwise local var
    };
    typedef util::GPoolDA<ContactConstraintD2D> PoolCCD2D;
    PoolCCD2D m_poolCCD2D;
#endif
};

} } // namespace S2::ds

#endif // S2_DS_SS_AGGREGATE_DSH_BASIC_H
