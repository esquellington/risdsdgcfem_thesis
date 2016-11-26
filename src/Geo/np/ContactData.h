#ifndef GEO_NP_CONTACT_DATA_H
#define GEO_NP_CONTACT_DATA_H

#include <Geo/Config.h>
#include <Geo/np/PairwiseCache.h>
#include <vector>
#include <util/triad.h>

#define __ENABLE_GEO_NP_CONTACTDATA_VIZDATA

namespace geo {
namespace np {

template <unsigned D>
struct GContactPoint
{
    typedef mal::GVec<Real,D> vec_type;
    vec_type m_Pos1;
    vec_type m_Pos2;
    vec_type m_Normal;
    Real m_Depth;
    Real m_Radius;

    finline GContactPoint() {}
    finline GContactPoint( const vec_type& p1, const vec_type& p2, const vec_type& normal, Real depth, Real radius )
    : m_Pos1(p1), m_Pos2(p2), m_Normal(normal), m_Depth(depth), m_Radius(radius) {}

    finline void Flip() { vec_type tmp(m_Pos2); m_Pos2 = m_Pos1; m_Pos1 = tmp; m_Normal = -m_Normal; }
};

typedef GContactPoint<2> ContactPoint2;
typedef GContactPoint<3> ContactPoint3;

/*! Contact Data
 */
template <unsigned D>
struct GContactData
{
    typedef mal::GVec<Real,D> vec_type;
    vec_type m_AvgNormal;
    Real m_AvgDepth;
    unsigned int m_NumDisjointManifolds;
    std::vector< GContactPoint<D> > m_vecCP;
    std::vector< std::pair<PointOnFeature,PointOnFeature> > m_vecPOF;

    finline unsigned int Size() const { return m_vecCP.size(); }
    finline const GContactPoint<D>& GetCP( int cp_idx ) const { return m_vecCP[cp_idx]; }

    /*\todo May be more efficient if each POF array is stored
      separatedly, so that one/both POF array may be empty and
      therefore shapes with no valid POF do not need to store them
      \todo Also, consider storing POF / FP as type-specific contact
      data / annotations, similar to ISpecificPairwiseCache2, that
      only exists when used.
    */
    finline bool HasPOF() const { return m_vecPOF.size() > 0; }
    finline const PointOnFeature& GetPOF1( int cp_idx ) const { return m_vecPOF[cp_idx].first; }
    finline const PointOnFeature& GetPOF2( int cp_idx ) const { return m_vecPOF[cp_idx].second; }

    inline void Flip()
        {
            m_AvgNormal = -m_AvgNormal;
            for(unsigned int it_cp=0; it_cp<m_vecCP.size(); it_cp++ ) m_vecCP[it_cp].Flip();
            for(unsigned int it_cp=0; it_cp<m_vecPOF.size(); it_cp++ ) std::swap( m_vecPOF[it_cp].first, m_vecPOF[it_cp].second );
        }

    //\name Construction
    //@{
    finline void Clear()
        {
            m_AvgNormal = vec_type::Zero();
            m_AvgDepth = Real(0);
            m_NumDisjointManifolds = 0;
            m_vecCP.clear();
            m_vecPOF.clear();
#  ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
            m_VD.Clear();
#  endif
        }
    finline void Begin() { Clear(); }
    finline void AddCP( const vec_type& p1, const vec_type& p2, const vec_type& normal, Real depth, Real radius = Real(1) )
        {
            GEO_PARANOID_ASSERT( !mal::IsNaN(p1) && !mal::IsNaN(p2) && !mal::IsNaN(normal) && !mal::IsNaN(depth) && !mal::IsNaN(radius) );
            m_vecCP.push_back( GContactPoint<D>(p1,p2,normal,depth,radius) );
        }

    finline void AddPOF( const PointOnFeature& pof1, const PointOnFeature& pof2 )
        {
            m_vecPOF.push_back( std::make_pair(pof1,pof2) );
        }
    finline void SetNumDisjointManifolds( unsigned int num_dm ) { m_NumDisjointManifolds = num_dm; }

    // Closes CD, computes averages, returns true if contact is non-null
    inline bool End()
    {
        GEO_ASSERT( !HasPOF() || m_vecCP.size() == m_vecPOF.size() );
        m_AvgNormal = vec_type::Zero();
        m_AvgDepth = Real(0);
        for( unsigned int it_cp=0; it_cp<m_vecCP.size(); it_cp++ )
        {
            m_AvgNormal += m_vecCP[it_cp].m_Normal;
            m_AvgDepth += m_vecCP[it_cp].m_Depth;
        }
        if( m_vecCP.size() > 0 )
        {
            m_AvgNormal /= m_vecCP.size();
            m_AvgNormal.Normalize();
            m_AvgDepth /= m_vecCP.size();
            return true;
        }
        else return false;
    }
    //\todo inline void Optimize( Real epsilon_length ) Contact point reduction, based in position proximity and normals, clustering, convexity, whatever...
    //@}

#  ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
public:
    struct VizData
    {
        void Clear()
            {
                m_vecIP.clear(); m_vecCP.clear(); m_vecNP.clear(); m_vecFP.clear(); m_vecRNP.clear(); m_vecPCA.clear();
                m_vecIM_XSegment.clear(); m_vecIM_XPoint.clear(); m_vecIM_XLength.clear();
                m_vecIM_IB1_Points.clear(); m_vecIM_IB2_Points.clear();
                m_vecIM_IB1_Normals.clear(); m_vecIM_IB2_Normals.clear();
                m_vecIM_IB1_Radii.clear(); m_vecIM_IB2_Radii.clear();
                m_vecIM_IB_Discarded.clear();
                m_vecIM_Map.clear(); m_vecIM_RMap.clear();
                m_vecDCR3_Vs_Primitive_vecClippedT_Segments.clear();
                m_vecDCR2DCR_vecIC.clear(); m_vecDCR2DCR_vecUnconnectedCTP.clear();
                m_vecDCR2DCR_IB_vecSeedHE_Side[0].clear(); m_vecDCR2DCR_IB_vecSeedHE_Side[1].clear();
                m_vecDCR2DCR_IB_vecInternalT_Side[0].clear(); m_vecDCR2DCR_IB_vecInternalT_Side[1].clear();
                m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[0].clear(); m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[1].clear();
                m_vecDCR2DCR_IB_PosAndNormal_Side_x_IB[0].clear(); m_vecDCR2DCR_IB_PosAndNormal_Side_x_IB[1].clear();
                m_vecDCR2DCR_IB_vecP_Side_x_IB[0].clear(); m_vecDCR2DCR_IB_vecP_Side_x_IB[1].clear();
                m_vecDCR2DCR_IB_vecClippedT_Side[0].clear(); m_vecDCR2DCR_IB_vecClippedT_Side[1].clear();
                m_vecDCR2DCR_IB_vecCE_Side[0].clear(); m_vecDCR2DCR_IB_vecCE_Side[1].clear();
                m_vecDCR2DCR_IB_vecCT_Side[0].clear(); m_vecDCR2DCR_IB_vecCT_Side[1].clear();
            }
        std::vector< std::pair<vec_type,vec_type> > m_vecIP; //Intersecting
        std::vector< std::pair<vec_type,vec_type> > m_vecCP; //Closest
        std::vector< std::pair<vec_type,vec_type> > m_vecNP; //Near
        std::vector< std::pair<vec_type,vec_type> > m_vecFP; //Far
        std::vector< std::pair<vec_type,vec_type> > m_vecRNP; //Random Neighbours
        std::vector< std::pair<vec_type,vec_type> > m_vecPCA; //PCA
        //\todo OverlapBV...
        //\todo Clusters...

        // IM
        std::vector< std::pair<vec_type,vec_type> > m_vecIM_XSegment;
        std::vector< vec_type > m_vecIM_XPoint;
        std::vector< std::pair<float,float> > m_vecIM_XLength;
        std::vector< vec_type > m_vecIM_IB1_Points;
        std::vector< vec_type > m_vecIM_IB2_Points;
        std::vector< vec_type > m_vecIM_IB_Discarded;
        std::vector< std::pair<vec_type,vec_type> > m_vecIM_Map;
        std::vector< std::pair<vec_type,vec_type> > m_vecIM_RMap;
        // Boundariel xtra
        std::vector< vec_type > m_vecIM_IB1_Normals;
        std::vector< vec_type > m_vecIM_IB2_Normals;
        std::vector< Real > m_vecIM_IB1_Radii;
        std::vector< Real > m_vecIM_IB2_Radii;

        //-- DCR3 vs Primitive
        std::vector< std::pair<vec_type,vec_type> > m_vecDCR3_Vs_Primitive_vecClippedT_Segments;

        //-- DCR3 vs DCR3
        // Intersection Curve
        std::vector< std::vector< std::pair<vec_type,vec_type> > > m_vecDCR2DCR_vecIC;
        std::vector< vec_type > m_vecDCR2DCR_vecUnconnectedCTP;
        std::vector< std::pair<vec_type,vec_type> > m_vecDCR2DCR_IB_vecSeedHE_Side[2];
        std::vector< util::triad<vec_type,vec_type,vec_type> > m_vecDCR2DCR_IB_vecInternalT_Side[2];
        std::vector< std::vector< util::triad<vec_type,vec_type,vec_type> > > m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[2];
        std::vector< std::pair<vec_type,vec_type> > m_vecDCR2DCR_IB_PosAndNormal_Side_x_IB[2];
        std::vector< std::vector< std::pair<vec_type,vec_type> > > m_vecDCR2DCR_IB_vecP_Side_x_IB[2];
        std::vector< util::triad<vec_type,vec_type,vec_type> > m_vecDCR2DCR_IB_vecClippedT_Side[2];
        std::vector< std::pair<vec_type,vec_type> > m_vecDCR2DCR_IB_vecCE_Side[2];
        std::vector< std::pair<vec_type,vec_type> > m_vecDCR2DCR_IB_vecCT_Side[2];
    };
    VizData m_VD;
#  endif
};

typedef GContactData<2> ContactData2;
typedef GContactData<3> ContactData3;

}} //namespace geo::np

#endif // GEO_NP_CONTACT_DATA_H
