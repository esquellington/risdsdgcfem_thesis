#ifndef S2_MS_CONSTRAINT_KINEMATIC_CONSTRAINT_SET_2_H
#define S2_MS_CONSTRAINT_KINEMATIC_CONSTRAINT_SET_2_H

#include "../Config.h"
#include <util/GPoolDA.h>

namespace S2 { namespace ms {

class KinematicNodeConstraintSet2
{
public:
    struct KNC;
    typedef KNC* knc_id_type;
    typedef uint16 node_index_type;
    enum EConstants { cInvalidKNCID = 0 };

    struct KNC
    {
        Vec2 m_Pos;
        Vec2 m_Vel;
        node_index_type m_NID;
        enum { eKNCT_Position, eKNCT_Velocity } m_Type;
    };
    typedef util::GPoolDA<KNC> PoolKNC;

public:
    KinematicNodeConstraintSet2()
    : m_poolKNC(16)
        {}
    ~KinematicNodeConstraintSet2() {}

    knc_id_type AddKNC( unsigned int nid, const Vec2 &pos, const Vec2 &vel )
        {
            KNC *pKNC = m_poolKNC.New();
            pKNC->m_Pos = pos;
            pKNC->m_Vel = vel;
            pKNC->m_NID = node_index_type(nid);
            pKNC->m_Type = KNC::eKNCT_Position;
            return knc_id_type(pKNC);
        }

    void RemoveKNC( knc_id_type kncid )
        {
            m_poolKNC.Delete(kncid);
        }

    void ResetKNC( knc_id_type kncid, const Vec2 &pos, const Vec2 &vel )
        {
            MS_ASSERT( m_poolKNC.IsValid(kncid) );
            kncid->m_Pos = pos;
            kncid->m_Vel = vel;
            kncid->m_Type = KNC::eKNCT_Position;
        }

    void SetPosKNC( knc_id_type kncid, const Vec2 &pos )
        {
            MS_ASSERT( m_poolKNC.IsValid(kncid) );
            kncid->m_Pos = pos;
            kncid->m_Type = KNC::eKNCT_Position;
        }

    void SetVelKNC( knc_id_type kncid, const Vec2 &vel )
        {
            MS_ASSERT( m_poolKNC.IsValid(kncid) );
            kncid->m_Vel = vel;
            kncid->m_Type = KNC::eKNCT_Velocity;
        }

    knc_id_type FindKNC( unsigned int nid ) const
        {
            for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
                if( nid == it_knc->m_NID )
                    return it_knc.GetPtr();
            return knc_id_type(0);
        }

    bool IsKNC( unsigned int nid ) const
        {
            for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
                if( nid == it_knc->m_NID )
                    return true;
            return false;
        }

    PoolKNC::iterator Begin() const
        {
            return m_poolKNC.Begin();
        }

private:
    PoolKNC m_poolKNC;
};

}} // namespace S2::ms

#endif //S2_MS_CONSTRAINT_KINEMATIC_CONSTRAINT_SET_2_H
