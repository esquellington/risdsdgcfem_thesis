#include <Saphyre2/bs/RayCastQuery.h>
#include <Saphyre2/bs/ISyncEntity.h>
#include <Saphyre2/ds/Commands.h>

#include <iostream>
#include <util/ItemStreamSerialization.h>

namespace S2
{

//----------------------------------------------------------------
// RC2D
//----------------------------------------------------------------

void RayCastQuery2D::Run( const Vec2& pos, const Vec2& dir,
                        const Interval& lambda_interval,
                        Real thickness,
                        Flags32 flags /*\todo filters....*/ )
{
    BS_BCMD("RayCastQuery2D::Run()");
    m_pChannel->BeginCommand( ds::eQuery_RayCast );
    {
        m_pChannel->Write( "qid", GetQID() );
        m_pChannel->Write( "eid", m_pEntity->GetEID() );
        m_pChannel->Write( "pos", pos );
        m_pChannel->Write( "dir", dir );
        m_pChannel->Write( "interval", lambda_interval );
        m_pChannel->Write( "thickness", thickness );
        m_pChannel->Write( "flags", flags );
    }
    m_pChannel->EndCommand();
    // Sync to obtain results
    m_NumHits = 0;
    bool bResult = m_pEntity->SyncMe( this ); //This is a "force-sync-me" request that always results in a universe force-sync by now
    //BS_LOG_WARNING( "RayCastQuery2D::Run #hits = %d", m_NumHits );
    BS_ECMD(bResult);
}

bool RayCastQuery2D::Closest( RayCastQuery2D::Hit2D& h2d ) const
{
    if( m_NumHits > 0 )
    {
        h2d = m_ClosestH2D;
        return true;
    }
    else
        return false;
}

bool RayCastQuery2D::ProcessReturn( const ds::ReturnIt& rit )
{
    m_NumHits = rit.GetSubItem().Find("num_hits").Get<uint32>();
    //BS_LOG_WARNING( "RayCastQuery2D::ProcessReturn #hits = %d", m_NumHits );
    m_ClosestH2D.m_RayHit.m_Interval = Intervalf::Empty(); //Clear closest
    uint32 num_hits(0);
    for( util::ItemStream::ItemIt it_hit = rit.GetSubItem().Find("array_hit").GetSubItem(); it_hit.IsValid(); ++it_hit )
    {
        //TEMPORAL: Hit reports are complex eType_Unknown by now
        BS_ASSERT( it_hit.IsComplex() && it_hit.GetType() == eType_Unknown );
        util::ItemStream::ItemIt it_hit_sub = it_hit.GetSubItem();
        Intervalf interval( it_hit_sub.Find("interval").Get<Intervalf>() );
        //if( m_ClosestH2D.m_RayHit2.m_Interval.IsEmpty() || m_ClosestH2D.m_RayHit2.m_Interval.Min() > interval.Min() )
        if( num_hits == 0 || m_ClosestH2D.m_RayHit.m_Interval.Min() > interval.Min() )
        {
            m_ClosestH2D.m_pEntity = reinterpret_cast<ISyncEntity*>( it_hit_sub.Find("uid").Get<machine_uint_type>() );
            m_ClosestH2D.m_RayHit.m_Interval = interval;
            m_ClosestH2D.m_RayHit.m_Point = it_hit_sub.Find("point").Get<Vec2f>();
            m_ClosestH2D.m_RayHit.m_Normal = it_hit_sub.Find("normal").Get<Vec2f>();
            m_ClosestH2D.m_RayHit.m_FeatureId.Set( (geo::EFeatureType)it_hit_sub.Find("feature_type").Get<uint32>(),
                                                   it_hit_sub.Find("feature_index").Get<geo::feature_index_type>() );
            m_ClosestH2D.m_RayHit.m_Extra_BarycentricCoords = it_hit_sub.Find("extra_bc").Get<Vec4f>();
        }
        num_hits++;
        /*
        BS_LOG_WARNING( "RayHit (%f,%f) entity %llx",
                        m_ClosestH2D.m_RayHit2.m_Interval.Min(),
                        m_ClosestH2D.m_RayHit2.m_Interval.Max(),
                        (machine_uint_type)m_ClosestH2D.m_pEntity );
        */
    }
    BS_ASSERT( num_hits == m_NumHits );
    return true;
}

//----------------------------------------------------------------
// RC3D
//----------------------------------------------------------------

void RayCastQuery3D::Run( const Vec3& pos, const Vec3& dir,
                          const Interval& lambda_interval,
                          Real thickness,
                          Flags32 flags /*\todo filters....*/ )
{
    BS_BCMD("RayCastQuery3D::Run()");
    m_pChannel->BeginCommand( ds::eQuery_RayCast );
    {
        m_pChannel->Write( "qid", GetQID() );
        m_pChannel->Write( "eid", m_pEntity->GetEID() );
        m_pChannel->Write( "pos", pos );
        m_pChannel->Write( "dir", dir );
        m_pChannel->Write( "interval", lambda_interval );
        m_pChannel->Write( "thickness", thickness );
        m_pChannel->Write( "flags", flags );
    }
    m_pChannel->EndCommand();
    // Sync to obtain results
    m_NumHits = 0;
    bool bResult = m_pEntity->SyncMe( this ); //This is a "force-sync-me" request that always results in a universe force-sync by now
    // BS_LOG_WARNING( "RayCastQuery3D::Run #hits = %d", m_NumHits );
    BS_ECMD(bResult);
}

bool RayCastQuery3D::Closest( RayCastQuery3D::Hit3D& h3d ) const
{
    if( m_NumHits > 0 )
    {
        h3d = m_ClosestH3D;
        return true;
    }
    else
        return false;
}

bool RayCastQuery3D::ProcessReturn( const ds::ReturnIt& rit )
{
    m_NumHits = rit.GetSubItem().Find("num_hits").Get<uint32>();
    //BS_LOG_WARNING( "RayCastQuery3D::ProcessReturn #hits = %d", m_NumHits );
    m_ClosestH3D.m_RayHit.m_Interval = Intervalf::Empty(); //Clear closest
    uint32 num_hits(0);
    for( util::ItemStream::ItemIt it_hit = rit.GetSubItem().Find("array_hit").GetSubItem(); it_hit.IsValid(); ++it_hit )
    {
        //TEMPORAL: Hit reports are complex eType_Unknown by now
        BS_ASSERT( it_hit.IsComplex() && it_hit.GetType() == eType_Unknown );
        util::ItemStream::ItemIt it_hit_sub = it_hit.GetSubItem();
        Intervalf interval( it_hit_sub.Find("interval").Get<Intervalf>() );
        //if( m_ClosestH3D.m_RayHit2.m_Interval.IsEmpty() || m_ClosestH3D.m_RayHit2.m_Interval.Min() > interval.Min() )
        if( num_hits == 0 || m_ClosestH3D.m_RayHit.m_Interval.Min() > interval.Min() )
        {
            m_ClosestH3D.m_pEntity = reinterpret_cast<ISyncEntity*>( it_hit_sub.Find("uid").Get<machine_uint_type>() );
            m_ClosestH3D.m_RayHit.m_Interval = interval;
            m_ClosestH3D.m_RayHit.m_Point = it_hit_sub.Find("point").Get<Vec3f>();
            m_ClosestH3D.m_RayHit.m_Normal = it_hit_sub.Find("normal").Get<Vec3f>();
            m_ClosestH3D.m_RayHit.m_FeatureId.Set( (geo::EFeatureType)it_hit_sub.Find("feature_type").Get<uint32>(),
                                                   it_hit_sub.Find("feature_index").Get<geo::feature_index_type>() );
            m_ClosestH3D.m_RayHit.m_Extra_BarycentricCoords = it_hit_sub.Find("extra_bc").Get<Vec4f>();
        }
        num_hits++;
        /*
        BS_LOG_WARNING( "RayHit (%f,%f) entity %llx",
                        m_ClosestH3D.m_RayHit.m_Interval.Min(),
                        m_ClosestH3D.m_RayHit.m_Interval.Max(),
                        (machine_uint_type)m_ClosestH3D.m_pEntity );
        */
    }
    BS_ASSERT( num_hits == m_NumHits );
    return true;
}


} //namespace S2
