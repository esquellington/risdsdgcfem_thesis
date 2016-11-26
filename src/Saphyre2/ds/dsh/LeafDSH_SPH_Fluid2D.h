#ifndef S2_DS_DSH_LEAF_SPH_FLUID2D_H
#define S2_DS_DSH_LEAF_SPH_FLUID2D_H

#include <Saphyre2/ds/dsh/GLeafDSH_Model.h>
#include <Saphyre2/ms/SPH_Fluid2D.h>
#include <util/ItemStream.h>
#include <vector>

#include <Saphyre2/ds/DSG.h> //!\todo THIS BREAKS ENCAPSULATION, but better do it here than in <Saphyre2/ms/SPH_Fluid2D.h>

namespace S2 { namespace ds {

class LeafDSH_SPH_Fluid2D: public GLeafDSH_Model< ms::SPH_Fluid2D >
{
public:
    LeafDSH_SPH_Fluid2D( uint32 uid, IDynamicSystemHierarchy *p_parent )
    : GLeafDSH_Model< ms::SPH_Fluid2D >(uid,p_parent)
    {}

    EDSHType GetDSHType() const { return eDSH_Leaf_Fluid2D; }

    void SyncGO()
    {
        if( 0 != base_type::m_pGO )
        {
            //\note Assumes that transform == identity, as all
            //information is in the SDOF
            m_pGO->SetTransform( Transform2::Identity() );
            m_Model.GetDOF( static_cast< geo::GObjectSDOF<2,Vec2>* >(m_pGO)->GetVecDOF_WriteOnly() );
        }
    }

    void RecomputeBV()
    {
        /* MAYBE the Model already computes BV and there's no need to redo it
        if( 0 != base_type::m_pGO && 0 != base_type::m_pBV )
            base_type::m_pGO->ComputeBV( *base_type::m_pBV );
        */
    }

    bool Create( const ParamIt &pit )
    {
        ms::SPH_Fluid2D::EApproachType approach;
        switch( (int32)pit.Find("flags").Get<Flags32>() ) //TEMPORAL hack
        {
        case 0: approach = ms::SPH_Fluid2D::eApproachDG96; break;
        case 1: approach = ms::SPH_Fluid2D::eApproachKeiser; break;
        case 2: approach = ms::SPH_Fluid2D::eApproachClavet; break;
        default: approach = ms::SPH_Fluid2D::eApproachDG96; break;
        };

        GetModel().SetVizStream( DSG::GetVizIS() ); //TEMPORAL ugly, but necessary by now

        GetModel().SetDscr( pit.Find("num_particles").Get<uint32>(),
                            pit.Find("density").Get<Real>(),
                            pit.Find("thickness").Get<Real>(),
                            pit.Find("init_aabb_posmin").Get<Point2>(),
                            pit.Find("init_aabb_posmax").Get<Point2>(),
                            pit.Find("bounds_aabb_posmin").Get<Vec2>(),
                            pit.Find("bounds_aabb_posmax").Get<Vec2>(),
                            approach );
        return true;
    }

    bool Edit( const ParamIt &pit )
    {
        ParamIt pit2 = pit.Find("pressure");
        if( pit2.IsValid() )
        {
            const RadialPressureAtPoint2D &rpap = pit2.Get< RadialPressureAtPoint2D >();
            GetModel().ApplyPressure( rpap.m_Pos, rpap.m_Radius, rpap.m_Pressure );
        }
        return true;
    }
};

} } // namespace S2::ds

#endif // S2_DS_DSH_LEAF_SPH_FLUID2D_H
