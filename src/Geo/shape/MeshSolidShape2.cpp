#include "MeshSolidShape2.h"
#include <memory.h> //TEMPORAL req by memcpy()

#include <Mal/GRandom.h> //for DomainSampler
#include <util/GIndexedRefCountedPoolDA.h> //for DomainSampler

#include "../util/GSimpleSpatialHash.h" //\todo USE SH in RebuildTopology()

// DCR \todo move elsewhere
#  include <Geo/shape/GPolygonalShape.h>
#  include <Geo/np/Intersection.h>

#ifdef __GEO_MSS_ENABLE_BVH
#  include <boost/bind.hpp> //for ILSS
#endif

#include <algorithm>

namespace geo
{

//-----------------------------------------------------------------------------
//---- MeshSolidShape2::DomainSampler Implementation
//-----------------------------------------------------------------------------

class MeshSolidShape2::DomainSampler: public GIDomainSamplerD<2>
{
public:
    //! Stores minimal local data to identify a POB and compute its geometry
    struct PointOnBoundary
    {
        PointOnBoundary( feature_index_type boundary_heid, Real lambda )
        : m_BoundaryHEID(boundary_heid), m_Lambda(lambda) {}
        ~PointOnBoundary() {}
        feature_index_type m_BoundaryHEID;
        Real m_Lambda;
    };

public:
    DomainSampler( const MeshSolidShape2 *p_shape ) : m_pShape(p_shape), m_vecSDOF(0), m_TotalArcLength(0) {}
    ~DomainSampler() {}

    unsigned int GetNumPOB() const { return m_poolPOB.Size(); }
    unsigned int GetNumAllocatedPOB() const { return m_poolPOB.NumAllocated(); }

    void UpdateDOF( const Real *vec_dof )
        {
            //TEMP: Compute actual_sdof (this is ugly... cast should be avoided, IShape should work with UNIFIED dof type (ether Real or GSRV)
            m_vecSDOF = reinterpret_cast<const MeshSolidShape2::sdof_type *>(vec_dof);
            //\todo recompute m_TotalArcLength
        }
    void ClearPOB() { m_poolPOB.Clear(); }

    //\name PointOnBoundary creation
    //@{
    point_on_boundary_id_type CreatePOB( const PointOnFeature& pof )
        {
            GEO_ASSERT( pof.m_FeatureId.IsSegment() );
            return m_poolPOB.New( PointOnBoundary( pof.m_FeatureId.AsSegment(), pof.m_BarycentricCoords[0] ) );
        }

    point_on_boundary_id_type StepPOB( point_on_boundary_id_type id,
                                       const vec_type& dir_local, Real step_length, Real eps_length )
        {
            GEO_ASSERT( m_poolPOB.IsValid(id) );
            GEO_ASSERT( mal::NormSq(dir_local) > Real(0) );
            GEO_ASSERT( step_length > Real(0) );
            GEO_ASSERT( eps_length >= Real(0) );

            point_on_boundary_id_type new_id(cInvalidPOB);

            /* If step_length is contained in the same HE, compute
               proper lambda, otherwise, advance/retreat to next/prev
               HE but set lambda=0 unconditionally.

               \note This stops the StepPOB() at vertices even if
               length<step_length which is probably the best decision
               for explicit meshes. However, notice that we ALWAYS
               progress, and the next step from the new POB will begin
               at either a larger/smaller lambda or the next/previous
               HE.

               \todo Consider input parameter b_stop_at_feature_change
               to give the user control over the stepping behaviour at
               feature changes.
            */
            const PointOnBoundary& seed_pob( m_poolPOB[id] );
            // Get POB BoundaryHE
            Vec2 he_p0( m_pShape->V_Pos( m_pShape->HE_OriginVID( seed_pob.m_BoundaryHEID ), m_vecSDOF ) );
            Vec2 he_p1( m_pShape->V_Pos( m_pShape->HE_FinalVID( seed_pob.m_BoundaryHEID ), m_vecSDOF ) );
            Real he_length( mal::Norm(he_p1-he_p0) );
            GEO_ASSERT( he_length > eps_length );
            Real inv_he_length( mal::Rcp(he_length) );
            Real pob_length( he_length * seed_pob.m_Lambda );

#define __USE_NORMALIZED_TANGENT
#ifdef __USE_NORMALIZED_TANGENT
            //\todo This is expensive, but fixes no-stepping bug
            //that arises too often when dir_local is very
            //small... CONSIDER requiring normalized dir_local input,
            //instead of normalizing here
            Vec2 tangent( mal::Normalized(dir_local) );
#else
            Vec2 tangent( dir_local );
#endif

            const Real dot_epsilon( 0.001f );
#define __USE_ROBUST_SELECT_DIRECTION
#ifdef __USE_ROBUST_SELECT_DIRECTION
            /*! Handle lambda=0 and lambda=1 cases specifically, as
              the dot_he_dir criterion FAILS for acute angles when
              the seed_pob is on a vertex.

              \todo Stepping direction computation is a bit risky. In
              2D it would be possible to just choose CW or CCW and
              pass it as a parameter to StePOB(), so that successive
              StepPOB() do NOT need to consider tangent nor compute
              stepping direction. However, in 3D this is not possible,
              and tangent is the only way to specify a preferred
              stepping direction... so we'll need to be VERY careful
              with acute angles in 3D...

              Also, simplified 2d stepping with CW/CCW parameter is
              discarded because it would require 2d-specific API in
              GIDomainSamplerD and we'd like to avoid losing
              genericity even at the expense of "slight" inefficiency.
            */
            // Select step direction
            enum { eForward, eBackwards, eNone } sd(eNone);
            if( seed_pob.m_Lambda == Real(0) )
            {
                // Get previous edge origin vertex
                Vec2 he_p_prev( m_pShape->V_Pos( m_pShape->HE_OriginVID( m_pShape->HE_Prev(seed_pob.m_BoundaryHEID) ), m_vecSDOF ) );
                // If tangent aligned with previous edge, we move forward, otherwise, backwards
                Real dot_he_prev_dir( mal::Dot( he_p0 - he_p_prev, tangent ) );
                if( dot_he_prev_dir > dot_epsilon ) sd = eForward;
                else if( dot_he_prev_dir < dot_epsilon ) sd = eBackwards;
                else sd = eNone;
            }
            else if( seed_pob.m_Lambda == Real(1) )
            {
                // Get next edge final vertex
                Vec2 he_p_next( m_pShape->V_Pos( m_pShape->HE_FinalVID( m_pShape->HE_Next(seed_pob.m_BoundaryHEID) ), m_vecSDOF ) );
                // If tangent aligned with previous edge, we move forward, otherwise, backwards
                Real dot_he_next_dir( mal::Dot( he_p_next - he_p1, tangent ) );
                if( dot_he_next_dir > dot_epsilon ) sd = eForward;
                else if( dot_he_next_dir < dot_epsilon ) sd = eBackwards;
                else sd = eNone;
            }
            else
            {
                //\todo THIS HAS VARIABLE SCALE BECAUSE tangent MAY NOT BE NORMALIZED
                Real dot_he_dir( mal::Dot( he_p1-he_p0, tangent ) * inv_he_length );
                if( dot_he_dir > dot_epsilon ) sd = eForward;
                else if( dot_he_dir < dot_epsilon ) sd = eBackwards;
                else sd = eNone;
            }
            // Step in the proper direction
            switch( sd )
            {
            case eForward:
                {
                    //GEO_LOG(" POB[%d] = (%d,%f) => Forward", id, m_poolPOB[id].m_BoundaryHEID, m_poolPOB[id].m_Lambda );
                    //\todo if( pob_length < he_length - eps_length ) THIS pob should be on a vtx... advance to next edge unconditionally!
                    Real new_lambda( (pob_length + step_length) * inv_he_length );
                    if( new_lambda < Real(1) - eps_length*inv_he_length  ) new_id = m_poolPOB.New( PointOnBoundary( seed_pob.m_BoundaryHEID, new_lambda ) );
                    else new_id = m_poolPOB.New( PointOnBoundary( m_pShape->HE_Next(seed_pob.m_BoundaryHEID), Real(0) ) );
                }
                break;
            case eBackwards:
                {
                    //GEO_LOG(" POB[%d] = (%d,%f) => Backwards", id, m_poolPOB[id].m_BoundaryHEID, m_poolPOB[id].m_Lambda );
                    //\todo if( pob_length > eps_length ) THIS pob should be on a vtx... retreat to prev edge unconditionally!
                    Real new_lambda( (pob_length - step_length) * inv_he_length );
                    if( new_lambda > eps_length*inv_he_length  ) new_id = m_poolPOB.New( PointOnBoundary( seed_pob.m_BoundaryHEID, new_lambda ) );
                    else new_id = m_poolPOB.New( PointOnBoundary( m_pShape->HE_Prev( seed_pob.m_BoundaryHEID ), Real(1) ) );
                }
                break;
            case eNone:
                {
                    GEO_LOG_ERROR("MeshSolidShape2::StepPOB() cannot advance, dir_local short OR ortogonal to edge_dir. Returning copy of POB[%d]", id );
                    //\todo IF THIS HAPPENS... consider passing 2 POB to
                    //STEP, so that the direction is implicit in their
                    //ordering... this wouldn't be so easy in 3d though.
                    new_id = m_poolPOB.New( seed_pob );
                    //m_poolPOB.AddRef(id); new_id = id; //??\todo No need to create a different POB...
                }
                break;
            default: GEO_LOG_ASSERT( false, "IMPOSSIBLE!" ); break;
            }
#else
            // Select direction
            Real dot_he_dir( mal::Dot( he_p1-he_p0, tangent ) * inv_he_length ); //\todo THIS HAS VARIABLE SCALE BECAUSE dir_local MAY NOT BE NORMALIZED
            const Real dot_epsilon( 0.001f );
            if( dot_he_dir > dot_epsilon ) //fwd
            {
                GEO_LOG(" POB[%d] = (%d,%f) => Forward", id, m_poolPOB[id].m_BoundaryHEID, m_poolPOB[id].m_Lambda );
                //\todo if( pob_length < he_length - eps_length ) THIS pob should be on a vtx... advance to next edge unconditionally!
                Real new_lambda( (pob_length + step_length) * inv_he_length );
                if( new_lambda < Real(1) - eps_length*inv_he_length  ) new_id = m_poolPOB.New( PointOnBoundary( seed_pob.m_BoundaryHEID, new_lambda ) );
                else new_id = m_poolPOB.New( PointOnBoundary( m_pShape->HE_Next(seed_pob.m_BoundaryHEID), Real(0) ) );
            }
            else if( dot_he_dir < dot_epsilon ) //bckwd
            {
                GEO_LOG(" POB[%d] = (%d,%f) => Backwards", id, m_poolPOB[id].m_BoundaryHEID, m_poolPOB[id].m_Lambda );
                //\todo if( pob_length > eps_length ) THIS pob should be on a vtx... retreat to prev edge unconditionally!
                Real new_lambda( (pob_length - step_length) * inv_he_length );
                if( new_lambda > eps_length*inv_he_length  ) new_id = m_poolPOB.New( PointOnBoundary( seed_pob.m_BoundaryHEID, new_lambda ) );
                else new_id = m_poolPOB.New( PointOnBoundary( m_pShape->HE_Prev( seed_pob.m_BoundaryHEID ), Real(1) ) );
            }
            else
            {
                GEO_LOG_ERROR("MeshSolidShape2::StepPOB() cannot advance, dir_local very short OR perpendicular to edge_dir. Returning copy of POB[%d]", id );
                //\todo IF THIS HAPPENS... consider passing 2 POB to
                //STEP, so that the direction is implicit in their
                //ordering... this wouldn't be so easy in 3d though.
                new_id = m_poolPOB.New( seed_pob );
                //m_poolPOB.AddRef(id); new_id = id; //??\todo No need to create a different POB...
            }
#endif
            //\todo IF boundary can be open, then some edges MAY NOT HAVE
            //A NEIGHBOUR!! ....This will happen with PathShape2, and
            //MeshSurfaceShape3, NOT with MeshSolidShape2,
            //Check actual stepping, ALWAYS possible in closed boundaries
            GEO_LOG_ASSERT( mal::NormSq( POB_Position(new_id) - POB_Position(id) ) > 0,
                            "[%d]=(%d,%f), [%d]=(%d,%f), sd = %s \n ep0=(%f,%f), ep1=(%f,%f), el=%f, t=(%f,%f)",
                            id, m_poolPOB[id].m_BoundaryHEID, m_poolPOB[id].m_Lambda,
                            new_id, m_poolPOB[new_id].m_BoundaryHEID, m_poolPOB[new_id].m_Lambda,
                            (sd==eForward) ? "Forward" : ( (sd==eBackwards) ? "Backwards" : "None" ),
                            he_p0[0], he_p0[1], he_p1[0], he_p1[1], he_length, tangent[0], tangent[1] );
            return new_id;
        }
    void CreateNeighbourPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                             point_on_boundary_id_type id,
                             Real neighbourhood_radius, Real neighbourhood_epsilon )
        {
            GEO_ASSERT( neighbourhood_radius > 0 );
            GEO_ASSERT( count == 2 );
            /* \note This code is adapted from StepPOB() to step both
               ways unconditionally, see comments there
            */
            const PointOnBoundary& seed_pob( m_poolPOB[id] );
            // Get POB BoundaryHE
            Vec2 he_p0( m_pShape->V_Pos( m_pShape->HE_OriginVID( seed_pob.m_BoundaryHEID ), m_vecSDOF ) );
            Vec2 he_p1( m_pShape->V_Pos( m_pShape->HE_FinalVID( seed_pob.m_BoundaryHEID ), m_vecSDOF ) );
            Real he_length( mal::Norm(he_p1-he_p0) );
            Real pob_length( he_length * seed_pob.m_Lambda );
            //bck
            Real bck_lambda( (pob_length - neighbourhood_radius) / he_length );
            if( bck_lambda*he_length > neighbourhood_epsilon ) vec_id[0] = m_poolPOB.New( PointOnBoundary( seed_pob.m_BoundaryHEID, bck_lambda ) );
            else vec_id[0] = m_poolPOB.New( PointOnBoundary( m_pShape->HE_Prev( seed_pob.m_BoundaryHEID ), Real(1) ) );
            //fwd
            Real fwd_lambda( (pob_length + neighbourhood_radius) / he_length );
            if( fwd_lambda*he_length < he_length-neighbourhood_epsilon ) vec_id[1] = m_poolPOB.New( PointOnBoundary( seed_pob.m_BoundaryHEID, fwd_lambda ) );
            else vec_id[1] = m_poolPOB.New( PointOnBoundary( m_pShape->HE_Next(seed_pob.m_BoundaryHEID), Real(0) ) );
            //\todo Check if the SAME POB is returned...

            //\todo IF MESH can be open, then some edges MAY NOT HAVE
            //A NEIGHBOUR!! ....This will happen with PathShape2, and
            //MeshSurfaceShape3, NOT with MeshSolidShape2.
        }

#define __ENABLE_FIND_CLOSEST_POB
#ifdef __ENABLE_FIND_CLOSEST_POB
    /*! Helper method to find the closest POB incrementally, stopping when exact point or boundary vertex found. */
    PointOnBoundary Find_ClosestPOB( const PointOnBoundary& pob, const vec_type& pos_local,
                                     Real neighbourhood_length, Real neighbourhood_epsilon,
                                     Real& sign, Real& remaining_length )
        {
            // Get POB BoundaryHE data
            Vec2 he_p0( m_pShape->V_Pos( m_pShape->HE_OriginVID( pob.m_BoundaryHEID ), m_vecSDOF ) );
            Vec2 he_p1( m_pShape->V_Pos( m_pShape->HE_FinalVID( pob.m_BoundaryHEID ), m_vecSDOF ) );
            Real he_length( mal::Norm(he_p1-he_p0) );
            Real inv_he_length( mal::Rcp(he_length) );
            Real pob_length( he_length * pob.m_Lambda );
            // Compute tangent and its step-length projection, clamped to neighbourhood radius
            Vec2 tangent( (he_p1 - he_p0) * inv_he_length );
            Vec2 p( he_p0 + pob.m_Lambda * (he_p1-he_p0) );
            Vec2 v( pos_local - p );
            Real step_length( mal::Clamp( mal::Dot( v, tangent ), -neighbourhood_length, neighbourhood_length ) );
            if( step_length > neighbourhood_epsilon ) //fwd
            {
                sign = +1;
                Real new_lambda( (pob_length + step_length) * inv_he_length );
                if( new_lambda < Real(1) - neighbourhood_epsilon*inv_he_length  )
                {
                    remaining_length = 0;
                    return PointOnBoundary( pob.m_BoundaryHEID, new_lambda );
                }
                else
                {
                    remaining_length = neighbourhood_length - (he_length - pob_length);
                    return PointOnBoundary( m_pShape->HE_Next(pob.m_BoundaryHEID), Real(0) );
                }
            }
            else if( step_length < -neighbourhood_epsilon ) //bckwd
            {
                sign = -1;
                Real new_lambda( (pob_length + step_length) * inv_he_length ); //step_length < 0!
                if( new_lambda > neighbourhood_epsilon*inv_he_length )
                {
                    remaining_length = 0;
                    return PointOnBoundary( pob.m_BoundaryHEID, new_lambda );
                }
                else
                {
                    remaining_length = neighbourhood_length - pob_length;
                    return PointOnBoundary( m_pShape->HE_Prev( pob.m_BoundaryHEID ), Real(1) );
                }
            }
            else
            {
                sign = 0;
                remaining_length = 0;
                return pob;
            }
        }
#endif //__ENABLE_FIND_CLOSEST_POB

    point_on_boundary_id_type ClosestPOB( point_on_boundary_id_type id,
                                          const vec_type& pos_local,
                                          Real neighbourhood_radius, Real neighbourhood_epsilon )
        {
            GEO_ASSERT( neighbourhood_radius > 0 );
            GEO_ASSERT( neighbourhood_epsilon >= 0 );

#ifdef __ENABLE_FIND_CLOSEST_POB
            const PointOnBoundary& seed_pob( m_poolPOB[id] );
            // Perform initial step
            Real seed_sign(0);
            Real remaining_length(0);
            PointOnBoundary pob = Find_ClosestPOB( seed_pob, pos_local,
                                                   neighbourhood_radius, neighbourhood_epsilon,
                                                   seed_sign, remaining_length );
            // Step while sign is not reversed and remaining_length > 0
            Real sign( seed_sign );
            if( seed_sign != 0 )
            {
                while( sign == seed_sign && remaining_length > 0 )
                    pob = Find_ClosestPOB( pob, pos_local,
                                           remaining_length, neighbourhood_epsilon,
                                           sign, remaining_length );
                /* If sign = 0, we found exact projection on a
                   segment, if sign reversed, the closest point is on
                   a vertex that we crossed in seed_sign direction and
                   that we hit again in the opposite direction due to
                   projection sign inversion. In both cases, the POB
                   is correct.

                   \note remaining_length reduction only effects
                   seed_sign advancing, not the possibility of sign
                   inversion that results in vertex-POB, as we never
                   need to actually ADVANCE backwards.
                */
                return m_poolPOB.New( pob );
            }
            else
                return id;
#else
            const PointOnBoundary& seed_pob( m_poolPOB[id] );
            // Get POB BoundaryHE data
            Vec2 he_p0( m_pShape->V_Pos( m_pShape->HE_OriginVID( seed_pob.m_BoundaryHEID ), m_vecSDOF ) );
            Vec2 he_p1( m_pShape->V_Pos( m_pShape->HE_FinalVID( seed_pob.m_BoundaryHEID ), m_vecSDOF ) );
            Real he_length( mal::Norm(he_p1-he_p0) );
            Real inv_he_length( mal::Rcp(he_length) );
            Real pob_length( he_length * seed_pob.m_Lambda );
            // Compute tangent and its step-length projection, clamped to neighbourhood radius
            Vec2 tangent( (he_p1 - he_p0) * inv_he_length );
            Vec2 p( he_p0 + seed_pob.m_Lambda * (he_p1-he_p0) );
            Vec2 v( pos_local - p );
            Real step_length( mal::Clamp( mal::Dot( v, tangent ), -neighbourhood_radius, neighbourhood_radius ) );
            // Advance step-length fwd/bckwd
            //\todo BY NOW WE NEVER continue on the next/previous edge, which makes CPOB usually stick to vertices incorrectly
            point_on_boundary_id_type new_id(cInvalidPOB);
            if( step_length > neighbourhood_epsilon ) //fwd
            {
                Real new_lambda( (pob_length + step_length) * inv_he_length );
                if( new_lambda < Real(1) - neighbourhood_epsilon*inv_he_length  ) new_id = m_poolPOB.New( PointOnBoundary( seed_pob.m_BoundaryHEID, new_lambda ) );
                else new_id = m_poolPOB.New( PointOnBoundary( m_pShape->HE_Next(seed_pob.m_BoundaryHEID), Real(0) ) );
            }
            else if( step_length < -neighbourhood_epsilon ) //bckwd
            {
                Real new_lambda( (pob_length + step_length) * inv_he_length ); //step_length < 0!
                if( new_lambda > neighbourhood_epsilon*inv_he_length ) new_id = m_poolPOB.New( PointOnBoundary( seed_pob.m_BoundaryHEID, new_lambda ) );
                else new_id = m_poolPOB.New( PointOnBoundary( m_pShape->HE_Prev( seed_pob.m_BoundaryHEID ), Real(1) ) );
            }
            else // id is already the closest POB
                new_id = id;
            return new_id;
#endif //__ENABLE_FIND_CLOSEST_POB
        }

    point_on_boundary_id_type CreateRandomPOB()
        {
            //\note Just use vec_id version... no point in repeating complex code
            point_on_boundary_id_type pob(0);
            CreateRandomPOB( &pob, 1 );
            return pob;
        }

    /* FAST_APPROX version can be considered mostly uniform if all
       boundary edges have the same length within a moderately small
       relative epsilon (ex 0.1...)
    */
    void CreateRandomPOB( point_on_boundary_id_type *vec_id, unsigned int count )
        {
#define __ENABLE_FAST_APPROX_CREATE_RANDOM_POB
#ifdef __ENABLE_FAST_APPROX_CREATE_RANDOM_POB
            /* FAST_APPROX version can be considered uniform if all boundary
               edges have the same length within a small epsilon.
            */
            // Generate count pob with random BHEID and lambda
            GEO_ASSERT( m_pShape->GetNumBoundaryP() == 1 );
            //\todo ASSUME SINGLE BP by now, if > 1, we should select randomly across ALL BP HE, not BP per BP...
            for( unsigned int it_bp=0; it_bp < m_pShape->GetNumBoundaryP(); it_bp++ )
            {
                unsigned int first_bheid( m_pShape->BP_FirstHEID(it_bp) );
                unsigned int num_bhe( m_pShape->BP_NumEdges(it_bp) );
                //GEO_LOG("BP[%d], FirstBHEID = %d, #BHE = %d", it_bp, first_bheid, num_bhe );
                for( unsigned int it_pob=0; it_pob < count; it_pob++ )
                    vec_id[it_pob] = m_poolPOB.New( PointOnBoundary( mal::RandomI<uint32>( first_bheid, first_bheid + num_bhe - 1 ),
                                                                     mal::RandomF<Real>( Real(0), Real(1) ) ) ) ;
            }
#else
            /* SLOW_UNIFORM version guarantees strictly uniform
             sampling, but at a huge O(#BE * #POB) cost
            */
            // Compute global arc length of ALL Boundary Polygons
            //\todo THIS COULD BE COMPUTED JUST ONCE AT UpdateDOF()!!
            Real total_arc_length(0);
            for( unsigned int it_bp=0; it_bp < m_pShape->GetNumBoundaryP(); it_bp++ )
            {
                unsigned int it_he( m_pShape->BP_FirstHEID(it_bp) );
                Vec2 p0( m_pShape->V_Pos( m_pShape->HE_OriginVID(it_he), m_vecSDOF ) );
                do
                {
                    it_he = m_pShape->HE_Next(it_he);
                    Vec2 p1( m_pShape->V_Pos( m_pShape->HE_OriginVID(it_he), m_vecSDOF ) );
                    total_arc_length += mal::Norm(p1-p0);
                    p0 = p1;
                } while ( it_he != m_pShape->BP_FirstHEID(it_bp) );
            }
            //GEO_LOG_WARNING( "total_arc_length = %f", total_arc_length );

            // Generate count pob with random arc length (\note we store random arc lengh as POB lambda temporarily)
            for( unsigned int it_pob=0; it_pob < count; it_pob++ )
                vec_id[it_pob] = m_poolPOB.New( PointOnBoundary( cInvalidFeatureIndex, mal::RandomF<Real>( 0, total_arc_length ) ) );

            //\todo CONSIDER: Sort pob with increasing global arc length

            // Iterate with increasing arc length and when matching pob is found, compute pob local data
            Real acc_arc_length = 0;
            for( unsigned int it_bp=0; it_bp < m_pShape->GetNumBoundaryP(); it_bp++ )
            {
                unsigned int it_he( m_pShape->BP_FirstHEID(it_bp) );
                Vec2 p0( m_pShape->V_Pos( m_pShape->HE_OriginVID(it_he), m_vecSDOF ) );
                do
                {
                    Vec2 p1( m_pShape->V_Pos( m_pShape->HE_FinalVID(it_he), m_vecSDOF ) );
                    Real length( mal::Norm(p1-p0) );
                    // Finish POB generation with local params
                    //\todo CONSIDER sorting generated POB by m_ArcLength to avoid having to check ALL generated POB for arc_length interval!!
                    for( unsigned int it_pob=0; it_pob<count; it_pob++ )
                    {
                        PointOnBoundary& pob( m_poolPOB[ vec_id[it_pob] ] );
                        if( pob.m_BoundaryHEID == cInvalidFeatureIndex
                            && pob.m_Lambda >= acc_arc_length
                            && pob.m_Lambda < acc_arc_length + length  )
                        {
                            pob.m_BoundaryHEID = it_he;
                            pob.m_Lambda = (m_poolPOB[ vec_id[it_pob] ].m_Lambda - acc_arc_length) / length;
                        }
                    }
                    // Advance
                    acc_arc_length += length;
                    it_he = m_pShape->HE_Next(it_he);
                    p0 = p1;
                } while ( it_he != m_pShape->BP_FirstHEID(it_bp) );
            }
#endif
        }

    void CreateRandomNeighbourPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                                   Real neighbourhood_radius, point_on_boundary_id_type id )
        {
            const PointOnBoundary& seed_pob( m_poolPOB[id] );
            unsigned int num_pob(0);
            // Generate count pob with random arc length in radius around POB[id]
            for( unsigned int it_pob=0; it_pob < count; it_pob++ )
                vec_id[it_pob] = m_poolPOB.New( PointOnBoundary( cInvalidFeatureIndex,
                                                                 mal::RandomF<Real>( -neighbourhood_radius,
                                                                                     neighbourhood_radius ) ) );
            // Generate neighbours on the same BoundaryHE
            Vec2 he_p0( m_pShape->V_Pos( m_pShape->HE_OriginVID( seed_pob.m_BoundaryHEID ), m_vecSDOF ) );
            Vec2 he_p1( m_pShape->V_Pos( m_pShape->HE_FinalVID( seed_pob.m_BoundaryHEID ), m_vecSDOF ) );
            Real he_length( mal::Norm(he_p1-he_p0) );
            Real pob_length( he_length * seed_pob.m_Lambda );
            for( unsigned int it_pob=0; it_pob < count; it_pob++ )
            {
                PointOnBoundary& pob( m_poolPOB[ vec_id[it_pob] ] );
                if( m_poolPOB[ vec_id[it_pob] ].m_BoundaryHEID == cInvalidFeatureIndex
                    && pob_length + m_poolPOB[ vec_id[it_pob] ].m_Lambda >= 0
                    && pob_length + m_poolPOB[ vec_id[it_pob] ].m_Lambda < he_length )
                {
                    m_poolPOB[ vec_id[it_pob] ].m_BoundaryHEID = seed_pob.m_BoundaryHEID;
                    m_poolPOB[ vec_id[it_pob] ].m_Lambda = (pob_length + m_poolPOB[ vec_id[it_pob] ].m_Lambda) / he_length;
                    if( m_poolPOB[ vec_id[it_pob] ].m_Lambda < Real(0) || m_poolPOB[ vec_id[it_pob] ].m_Lambda > Real(1) )
                        GEO_LOG_WARNING("Incorrect lambda mid %f", m_poolPOB[ vec_id[it_pob] ].m_Lambda );
                    num_pob++;
                }
            }
            // Generate forward: From NextHEID, iterate with increasing arc length and when matching pob is found, compute pob local data
            unsigned int num_pob_fwd(0);
            if( num_pob < count )
            {
                Real diff_arc_length( he_length - pob_length );
                unsigned int it_he( m_pShape->HE_Next( seed_pob.m_BoundaryHEID ) );
                Vec2 p0( m_pShape->V_Pos( m_pShape->HE_OriginVID(it_he), m_vecSDOF ) );
                do
                {
                    Vec2 p1( m_pShape->V_Pos( m_pShape->HE_FinalVID(it_he), m_vecSDOF ) );
                    Real length( mal::Norm(p1-p0) );
                    // Finish POB generation with local params
                    //\todo CONSIDER sorting generated POB by m_ArcLength to avoid having to check ALL generated POB for arc_length interval!!
                    for( unsigned int it_pob=0; it_pob<count; it_pob++ )
                    {
                        PointOnBoundary& pob( m_poolPOB[ vec_id[it_pob] ] );
                        if( pob.m_BoundaryHEID == cInvalidFeatureIndex
                            && pob.m_Lambda >= diff_arc_length
                            && pob.m_Lambda < diff_arc_length + length  )
                        {
                            pob.m_BoundaryHEID = it_he;
                            pob.m_Lambda = (pob.m_Lambda - diff_arc_length) / length;
                            if( pob.m_Lambda < Real(0) || pob.m_Lambda > Real(1) )
                                GEO_LOG_WARNING("Incorrect lambda fwd %f", pob.m_Lambda );
                            num_pob++;
                            num_pob_fwd++;
                        }
                    }
                    // Advance
                    diff_arc_length += length;
                    it_he = m_pShape->HE_Next(it_he);
                    p0 = p1;
                } while ( diff_arc_length <= neighbourhood_radius && it_he != seed_pob.m_BoundaryHEID );
            }
            // Generate backwards: From PredHEID, iterate with increasing arc length and when matching pob is found, compute pob local data
            //\note Backwards global arc_length lambdas MUST be negative, and local lambdas are computed as 1-x because p0->p1 points backwards!
            unsigned int num_pob_bwd(0);
            if( num_pob < count )
            {
                Real diff_arc_length( pob_length );
                unsigned int it_he( m_pShape->HE_Prev( seed_pob.m_BoundaryHEID ) );
                Vec2 p0( m_pShape->V_Pos( m_pShape->HE_FinalVID(it_he), m_vecSDOF ) );
                do
                {
                    Vec2 p1( m_pShape->V_Pos( m_pShape->HE_OriginVID(it_he), m_vecSDOF ) );
                    Real length( mal::Norm(p1-p0) );
                    // Finish POB generation with local params
                    //\todo CONSIDER sorting generated POB by m_ArcLength to avoid having to check ALL generated POB for arc_length interval!!
                    for( unsigned int it_pob=0; it_pob<count; it_pob++ )
                    {
                        PointOnBoundary& pob( m_poolPOB[ vec_id[it_pob] ] );
                        if( pob.m_BoundaryHEID == cInvalidFeatureIndex
                            && -pob.m_Lambda >= diff_arc_length
                            && -pob.m_Lambda < diff_arc_length + length  )
                        {
                            pob.m_BoundaryHEID = it_he;
                            pob.m_Lambda = Real(1) - (-pob.m_Lambda - diff_arc_length) / length;
                            if( pob.m_Lambda < Real(0) || pob.m_Lambda > Real(1) )
                                GEO_LOG_WARNING("Incorrect lambda bckwd %f", pob.m_Lambda );
                            num_pob++;
                            num_pob_bwd++;
                        }
                    }
                    // Advance
                    diff_arc_length += length;
                    it_he = m_pShape->HE_Prev(it_he);
                    p0 = p1;
                } while ( diff_arc_length <= neighbourhood_radius && it_he != seed_pob.m_BoundaryHEID );
            }
            //TEMP: Hack to FILL missing neighbours
            if( num_pob < count )
            {
                GEO_LOG_WARNING("MeshSolidShape2::DomainSampler::CreateRandomBoundaryNeighbours() Missing %d Neighbours!! NumPOB: Fwd %d, Bwd %d",
                                count - num_pob, (int32)num_pob_fwd, (int32)num_pob_bwd );
                for( unsigned int it_pob=0; it_pob<count; it_pob++ )
                {
                    PointOnBoundary& pob( m_poolPOB[ vec_id[it_pob] ] );
                    if( pob.m_BoundaryHEID == cInvalidFeatureIndex )
                    {
                        GEO_LOG_WARNING("\tUngenerated lambda: %f < %f, pobl = %f, hel = %f",
                                        m_poolPOB[vec_id[it_pob]].m_Lambda, neighbourhood_radius, pob_length, he_length );
                        pob = seed_pob;
                    }
                }
            }
        }

    //\name array-POB lifetime methods
    //@{
    void APOB_Destroy( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.Delete( vec_id[i] ); }
    void APOB_IncRef( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.IncRef( vec_id[i] ); }
    void APOB_DecRef( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.DecRef( vec_id[i] ); }
    //@}

    void AllPositionsAndNormals( vec_type *vec_pos, vec_type *vec_normal ) const
        {
            for( unsigned int it_pob=0; it_pob<m_poolPOB.NumAllocated(); it_pob++ ) //\note We iterate over ALL possible POB, but update only valid ones
                if( m_poolPOB.IsValid(it_pob) )
                {
                    Vec2 p0( m_pShape->V_Pos( m_pShape->HE_OriginVID( m_poolPOB[it_pob].m_BoundaryHEID ), m_vecSDOF ) );
                    Vec2 p1( m_pShape->V_Pos( m_pShape->HE_FinalVID( m_poolPOB[it_pob].m_BoundaryHEID ), m_vecSDOF ) );
                    vec_pos[it_pob] = p0 + m_poolPOB[it_pob].m_Lambda * (p1-p0);
                    vec_normal[it_pob] = mal::Normalized( mal::PerpendicularCW( p1-p0 ) );
                }
        }

    Vec2 POB_Position( point_on_boundary_id_type id ) const
        {
            Vec2 p0( m_pShape->V_Pos( m_pShape->HE_OriginVID( m_poolPOB[id].m_BoundaryHEID ), m_vecSDOF ) );
            Vec2 p1( m_pShape->V_Pos( m_pShape->HE_FinalVID( m_poolPOB[id].m_BoundaryHEID ), m_vecSDOF ) );
            return p0 + m_poolPOB[id].m_Lambda * (p1-p0);
        }
    Vec2 POB_Normal( point_on_boundary_id_type id ) const
        {
            Vec2 p0( m_pShape->V_Pos( m_pShape->HE_OriginVID( m_poolPOB[id].m_BoundaryHEID ), m_vecSDOF ) );
            Vec2 p1( m_pShape->V_Pos( m_pShape->HE_FinalVID( m_poolPOB[id].m_BoundaryHEID ), m_vecSDOF ) );
            return mal::Normalized( mal::PerpendicularCW( p1-p0 ) );
        }
    feature_id POB_FeatureId( point_on_boundary_id_type id ) const
        {
            //\todo MAY NEED TO ADD LAMBDA param to edge feature, and return eFT_Vertex lambda = 0,1
            return feature_id( eFT_Segment, m_poolPOB[id].m_BoundaryHEID );
        }

    PointOnFeature POB_PointOnFeature( point_on_boundary_id_type id ) const
        {
            return PointOnFeature( feature_id(eFT_Segment, m_poolPOB[id].m_BoundaryHEID),
                                   Vec4( m_poolPOB[id].m_Lambda, Real(1)-m_poolPOB[id].m_Lambda, 0, 0 ) );
        }

    const GIShapeD<2> *GetShape() const { return m_pShape; }

    //TEMP
    void Trace() const { m_poolPOB.Trace(); }

private:
    const MeshSolidShape2 *m_pShape;
    const MeshSolidShape2::sdof_type *m_vecSDOF;
    Real m_TotalArcLength;
    //util::GIndexedPoolDA<PointOnBoundary,point_on_boundary_id_type> m_poolPOB;
    util::GIndexedRefCountedPoolDA<PointOnBoundary,point_on_boundary_id_type> m_poolPOB;
};

//-----------------------------------------------------------------------------
//---- MeshSolidShape2 Implementation
//-----------------------------------------------------------------------------

MeshSolidShape2::MeshSolidShape2()
: m_NumV(0), m_NumP(0), m_NumHE(0)
, m_NumBoundaryP(0), m_NumBoundaryHE(0)
, m_NumL(0)
, m_NumAllocV(0), m_NumAllocP(0), m_NumAllocHE(0)
, m_vecPoints(0), m_vecV(0), m_vecP(0), m_vecHE(0)
, m_vecL(0)
, m_pBuffer(0)
, m_pDCR(0)
#ifdef __GEO_MSS_ENABLE_BVH
, m_pBVH(0)
#endif
{}

MeshSolidShape2::~MeshSolidShape2()
{
    ClearBakedData();
}

void MeshSolidShape2::ClearBakedData()
{
    if( m_pBuffer ) delete [] m_pBuffer;
    m_pBuffer = 0;
    m_NumV = 0; m_vecPoints = 0; m_vecV = 0;
    m_NumP = 0; m_vecP = 0;
    m_NumHE = 0; m_vecHE = 0;
    m_NumBoundaryP = 0; m_NumBoundaryHE = 0;
    m_NumAllocV = 0; m_NumAllocP = 0; m_NumAllocHE = 0;
    m_NumL = 0; m_vecL = 0;
    // DCR is strictly nonshared, by now
    if( m_pDCR ) delete m_pDCR;
    m_pDCR = 0;
#ifdef __GEO_MSS_ENABLE_BVH
    // BVH is strictly nonshared, by now
    if( m_pBVH ) delete m_pBVH;
    m_pBVH = 0;
#endif
}

/*\todo If the MSS has a BVH and it's up-to-date, here we should use
  it to Compute the required BV.
  OPTIMIZATION: MOREOVER, we could update the BVH here if it's not up
  to date and it is REQUIRED to compute the simple BV (it can easily
  happen for deformable objects if no conservative global BV is
  computed that includes the deformation... and in ALL CASES we'd like
  to avoid iterating all MSS geometry more than once to build the
  separate BV and the BVH...
  => NO! the BVH in an IShape is CONST, should contain only topology,
     the geometric/changing part should be in IObject, therefore,
     computebvd does NOT KNOW the actual BVH geometry unless IObject
     passes it as a param. This geometry part should be up-to-date
     when ComputeBVD() is called
*/
void MeshSolidShape2::ComputeBVD( bv::BoundingVolume2& bv, const transform_type& transform, const sdof_type *vec_sdof ) const
{
#ifdef __GEO_MSS_ENABLE_BVH_NOT_YEEEEEEEEEEEEEEET
    if( m_pBVH_Topology && p_bvh_geometry )
    {
        //\todo Use BVH[0] to compute BV
    }
#endif

    const sdof_type *actual_sdof( ( 0 != vec_sdof ) ? vec_sdof : m_vecPoints );
    switch( bv.GetType() )
    {
        /*\todo
    case bv::eBV_Sphere2:
        bv.As<bv::Sphere2>().SetPosRadius( transform.m_Pos, m_Radius );
        break;
        */
    case bv::eBV_AABB2:
        {
            bv::AABB2 aabb( transform*actual_sdof[0] );
            for( unsigned int i=1; i < m_NumV; i++ )
                aabb.Merge( transform*actual_sdof[i] );
            bv.As<bv::AABB2>() = aabb;
        }
        break;
    case bv::eBV_DOP2_K4:
        {
            bv::DOP2_K4 kdop( transform*actual_sdof[0] );
            for( unsigned int i=1; i < m_NumV; i++ )
                kdop.Merge( transform*actual_sdof[i] );
            bv.As<bv::DOP2_K4>() = kdop;
        }
        break;
    case bv::eBV_DOP2_K8:
        {
            bv::DOP2_K8 kdop( transform*actual_sdof[0] );
            for( unsigned int i=1; i < m_NumV; i++ )
                kdop.Merge( transform*actual_sdof[i] );
            bv.As<bv::DOP2_K8>() = kdop;
        }
        break;
        /*\todo
    case bv::eBV_LSS2:
        bv.As<bv::LSS2>().SetPosRadius( transform.m_Pos, transform.m_Pos, m_Radius );
        break;
        */
    case bv::eBV_Void: break;
    case bv::eBV_Infinite: break;
    default:
        GEO_LOG_ASSERT( false, "MeshSolidShape2::ComputeBVD() wrong BV type %d or dimension %d", (int32)bv.GetType(), (int32)bv.GetDimension() );
        break;
    }
}

IDomainSampler2 *MeshSolidShape2::CreateDomainSampler() const
{
    return new MeshSolidShape2::DomainSampler(this);
}

/* Init from external arrays

   If shared, arrays are NOT allocated and copied, only
   referenced. Otherwise all data allocated in m_pBuffer and copied,
   and the shape is self-contained.
*/
void MeshSolidShape2::SetBakedData( bool b_shared,
                                    uint32 num_v, uint32 num_p, uint32 num_he,
                                    uint32 num_boundary_p, uint32 num_boundary_he,
                                    uint32 num_l,
                                    const Vec2* vec_points,
                                    const vertex_type* vec_v,
                                    const polygon_type* vec_p,
                                    const half_edge_type* vec_he,
                                    const polygon_layer_type* vec_l )
{
    GEO_ASSERT( m_NumV == 0 && m_vecPoints == 0 && m_vecV == 0 && num_v > 0
                && m_NumP == 0 && m_vecP == 0 && num_p > 0
                && m_NumHE == 0 && m_vecHE == 0 && num_he > 0
                && m_NumBoundaryP == 0 && num_boundary_p > 0
                && m_NumBoundaryHE == 0 && num_boundary_he > 0
                && m_NumL == 0 && num_l > 0
                && m_pBuffer == 0 );
    m_NumV = num_v;
    m_NumP = num_p;
    m_NumHE = num_he;
    m_NumBoundaryP = num_boundary_p;
    m_NumBoundaryHE = num_boundary_he;
    m_NumL = num_l;
    m_NumAllocV = m_NumV;
    m_NumAllocP = m_NumP + m_NumBoundaryP;
    m_NumAllocHE = m_NumHE + m_NumBoundaryHE;
    if( b_shared )
    {
        m_vecPoints = vec_points;
        m_vecV = vec_v;
        m_vecP = vec_p;
        m_vecHE = vec_he;
        m_vecL = vec_l;
    }
    else
    {
        // Alloc and fill single buffer
        size_t size_points_4aligned = 4*((sizeof(Vec2)*m_NumAllocV+3)/4);
        size_t size_v_4aligned = 4*((sizeof(vertex_type)*m_NumAllocV+3)/4);
        size_t size_p_4aligned = 4*((sizeof(polygon_type)*m_NumAllocP+3)/4);
        size_t size_he_4aligned = 4*((sizeof(half_edge_type)*m_NumAllocHE+3)/4);
        size_t size_l_4aligned = 4*((sizeof(polygon_layer_type)*m_NumL+3)/4);
        size_t total_size_4aligned = size_points_4aligned + size_v_4aligned + size_p_4aligned + size_he_4aligned + size_l_4aligned;
        m_pBuffer = new uint32[ total_size_4aligned ];
        Vec2 *p_buffer_points = reinterpret_cast<Vec2*>( &m_pBuffer[0] );
        vertex_type *p_buffer_v = reinterpret_cast<vertex_type*>( &m_pBuffer[ size_points_4aligned ] );
        polygon_type *p_buffer_p = reinterpret_cast<polygon_type*>( &m_pBuffer[ size_points_4aligned + size_v_4aligned ] );
        half_edge_type *p_buffer_he = reinterpret_cast<half_edge_type*>( &m_pBuffer[ size_points_4aligned + size_v_4aligned + size_p_4aligned ] );
        polygon_layer_type *p_buffer_l = reinterpret_cast<polygon_layer_type*>( &m_pBuffer[ size_points_4aligned + size_v_4aligned + size_p_4aligned + size_he_4aligned ] );
        memcpy( p_buffer_points, vec_points, sizeof(Vec2)*m_NumAllocV );
        memcpy( p_buffer_v, vec_v, sizeof(vertex_type)*m_NumAllocV );
        memcpy( p_buffer_p, vec_p, sizeof(polygon_type)*m_NumAllocP );
        memcpy( p_buffer_he, vec_he, sizeof(half_edge_type)*m_NumAllocHE );
        memcpy( p_buffer_l, vec_l, sizeof(polygon_layer_type)*m_NumL );
        // Save const pointers
        m_vecPoints = p_buffer_points;
        m_vecV = p_buffer_v;
        m_vecP = p_buffer_p;
        m_vecHE = p_buffer_he;
        m_vecL = p_buffer_l;
    }

    /*TEMPORAL
    GEO_LOG_WARNING( "SetBakedData( %d, %d, %d; %d, %d )", num_v, num_p, num_he, num_boundary_p, num_boundary_he );
    for( unsigned int it_p=GetNumP(); it_p < GetNumP() + GetNumBoundaryP(); it_p++ )
    {
        GEO_LOG_WARNING( "BP[%d] = (%d,%d)", it_p, P_FirstHEID(it_p), P_NumEdges(it_p) );
    }
    */
}

Vec2 MeshSolidShape2::P_Barycenter_0( uint32 pid ) const
{
    Vec2 barycenter(0,0);
    unsigned int it_he( P_FirstHEID(pid) );
    do
    {
        barycenter += V_Pos_0( HE_OriginVID(it_he) );
        it_he = HE_Next(it_he);
    } while ( it_he != P_FirstHEID(pid) );
    return barycenter / Real( P_NumEdges(pid) );
}

Vec2 MeshSolidShape2::P_Barycenter( uint32 pid, const sdof_type *vec_sdof ) const
{
    Vec2 barycenter(0,0);
    unsigned int it_he( P_FirstHEID(pid) );
    do
    {
        barycenter += vec_sdof[ HE_OriginVID(it_he) ];
        it_he = HE_Next(it_he);
    } while ( it_he != P_FirstHEID(pid) );
    return barycenter / Real( P_NumEdges(pid) );
}

uint32 MeshSolidShape2::P_VecVID( uint32 pid, uint32* vec_vid, uint32 max_vid ) const
{
    uint32 num_vid( P_NumEdges(pid) );
    if( num_vid > max_vid )
    {
        GEO_LOG_ERROR( "P_VecVID(%u,XXX,%u) Polygon has too many vertices (%u)", pid, max_vid, num_vid );
        num_vid = max_vid;
    }
    unsigned int it_he( P_FirstHEID(pid) );
    for( uint32 i=0; i<num_vid; i++, it_he = HE_Next(it_he) )
        vec_vid[i] = HE_OriginVID(it_he);
    return num_vid;
}

//-----------------------------------------------------------------------------
//---- EditableMeshSolidShape2 Implementation
//-----------------------------------------------------------------------------
EditableMeshSolidShape2::EditableMeshSolidShape2()
: m_IsBeingEdited(false)
{
}

EditableMeshSolidShape2::~EditableMeshSolidShape2()
{
    Clear();
}

void EditableMeshSolidShape2::Clear()
{
    ClearEditData();
    MeshSolidShape2::ClearBakedData();
}

void EditableMeshSolidShape2::Set( const MeshSolidShape2& mss2 )
{
    Clear();
    // Copy other mesh baked data
    for( unsigned int i=0; i<mss2.GetNumV(); i++ )
        m_addV.push_back( editable_vertex_type( mss2.GetVecPoints()[i], mss2.GetVecV()[i].m_OutHEID, mss2.GetVecV()[i].m_NumEdges ) );
    for( unsigned int i=0; i<mss2.GetNumP(); i++ )
        m_addP.push_back( editable_polygon_type( mss2.GetVecP()[i].m_FirstHEID, mss2.GetVecP()[i].m_NumEdges ) );
    for( unsigned int i=0; i<mss2.GetNumHE(); i++ )
        m_addHE.push_back( editable_half_edge_type( mss2.GetVecHE()[i].m_OriginVID,
                                                    mss2.GetVecHE()[i].m_NextHEID,
                                                    mss2.GetVecHE()[i].m_SymHEID,
                                                    mss2.GetVecHE()[i].m_LeftPID ) );
    //\note Layers are always recomputed, not imported
    EndEdition();
}

void EditableMeshSolidShape2::BeginEdition()
{
    // Clear anything added outside Begin/End
    ClearEditData();
    // Copy baked data
    for( unsigned int i=0; i<m_NumV; i++ )
        m_addV.push_back( editable_vertex_type( m_vecPoints[i], m_vecV[i].m_OutHEID, m_vecV[i].m_NumEdges ) );
    for( unsigned int i=0; i<m_NumP; i++ )
        m_addP.push_back( editable_polygon_type( m_vecP[i].m_FirstHEID, m_vecP[i].m_NumEdges ) );
    for( unsigned int i=0; i<m_NumHE; i++ )
        m_addHE.push_back( editable_half_edge_type( m_vecHE[i].m_OriginVID,
                                                    m_vecHE[i].m_NextHEID,
                                                    m_vecHE[i].m_SymHEID,
                                                    m_vecHE[i].m_LeftPID ) );
    // Delete baked data
    MeshSolidShape2::ClearBakedData();
    m_IsBeingEdited = true;
}

bool EditableMeshSolidShape2::EndEdition()
{
    FixDegeneracies();
    RebuildTopology();
    RebuildLayers();
    RebuildBakedData();
    ClearEditData();
    m_IsBeingEdited = false;
    return true;
}

feature_index_type EditableMeshSolidShape2::AddVertex( const Vec2& point )
{
    GEO_ASSERT( m_addV.size() < cInvalidFeatureIndex-1 );
    /*IMPORTANT: NOOOOOO!!! This is wrong, it breaks Subdivide() if
      coincident vertices are generated, delay coincident vertex
      removal to FixDegenerateVertices()

    const Real epsilon_sq( mal::Sq( 0.001f ) ); //\todo should be parameter of constant somewhere
    feature_index_type vid = FindV( point, epsilon_sq );
    if( vid != cInvalidFeatureIndex )
        return vid;
    else
    {
    */
        m_addV.push_back( editable_vertex_type( point, cInvalidFeatureIndex, 0 ) );
        return feature_index_type(m_addV.size()-1);
        //}
}

feature_index_type EditableMeshSolidShape2::AddPolygon3( feature_index_type vid0, feature_index_type vid1, feature_index_type vid2 )
{
    feature_index_type vec_vid[3];
    vec_vid[0] = vid0;
    vec_vid[1] = vid1;
    vec_vid[2] = vid2;
    return AddPolygonN( 3, vec_vid );

    /* TEMP: Old triangle-specific code, kept as a reference, redundant but slightly faster.
    //\todo Consider looking for existing P
    GEO_ASSERT( m_addP.size() < cInvalidFeatureIndex-1 );
    GEO_ASSERT( cInvalidFeatureIndex == FindHE(vid0,vid1)
                && cInvalidFeatureIndex == FindHE(vid1,vid2)
                && cInvalidFeatureIndex == FindHE(vid2,vid0) );
    // Alloc new pid
    feature_index_type pid( m_addP.size() );
    // Find symmetrics
    feature_index_type sheid0( FindHE( vid1, vid0 ) );
    feature_index_type sheid1( FindHE( vid2, vid1 ) );
    feature_index_type sheid2( FindHE( vid0, vid2 ) );
    // Create 3 HE
    feature_index_type heid0( m_addHE.size() );
    feature_index_type heid1( heid0 + 1 );
    feature_index_type heid2( heid0 + 2 );
    m_addHE.push_back( editable_half_edge_type( vid0, heid1, sheid0, pid ) );
    m_addHE.push_back( editable_half_edge_type( vid1, heid2, sheid1, pid ) );
    m_addHE.push_back( editable_half_edge_type( vid2, heid0, sheid2, pid ) );
    // Create polygon
    m_addP.push_back( editable_polygon_type( heid0, 3 ) );
    // Link symmetric HE to this new polygon and its HE
    if( sheid0 != cInvalidFeatureIndex ) m_addHE[sheid0].m_SymHEID = heid0;
    if( sheid1 != cInvalidFeatureIndex ) m_addHE[sheid1].m_SymHEID = heid1;
    if( sheid2 != cInvalidFeatureIndex ) m_addHE[sheid2].m_SymHEID = heid2;
    return pid;
    */
}

feature_index_type EditableMeshSolidShape2::AddPolygonN( uint32 num_vid, feature_index_type *vec_vid )
{
    //\todo Consider looking for existing P
    GEO_ASSERT( m_addP.size() < cInvalidFeatureIndex-1 );
    for( unsigned int i=0; i<num_vid; i++ )
        GEO_ASSERT( cInvalidFeatureIndex == FindHE( vec_vid[i], vec_vid[(i+1)%num_vid] ) );
    // Alloc new pid
    feature_index_type pid( m_addP.size() );
    feature_index_type heid0( m_addHE.size() );
    // Add HE
    for( unsigned int i=0; i<num_vid; i++ )
    {
        unsigned int j( (i+1) % num_vid );
        feature_index_type vid0( vec_vid[i] );
        feature_index_type vid1( vec_vid[j] );
        // Find symmetric
        feature_index_type sheid0( FindHE( vid1, vid0 ) );
        // Link symmetric HE to this new polygon and its HE
        if( sheid0 != cInvalidFeatureIndex ) m_addHE[sheid0].m_SymHEID = m_addHE.size();
        // Add fully linked HE
        m_addHE.push_back( editable_half_edge_type( vid0, heid0 + j, sheid0, pid ) );
    }
    // Create polygon
    m_addP.push_back( editable_polygon_type( heid0, num_vid ) );
    return pid;
}

/* Move an existing vertex.
   \note This will NOT change topology... unless vertices are made
   coincident! therefore we MUST fully-rebuild in EndEdition()
   afterwards
*/
void EditableMeshSolidShape2::SetVertex( uint32 vid, const Vec2& point )
{
    GEO_ASSERT( IsBeingEdited() );
    GEO_ASSERT( vid < m_addV.size() );
    GEO_ASSERT( !IsNaN(point) );
    m_addV[vid].m_Pos = point;
}

/* Subdivision (Linear or Loop)
   - Linear: edge-split 1 Tri => 4 SubTri
   - Smoothing: Loop formulas, see "2000_SIGGRAPH_Course_SubdivisionForModelingAndAnimation.pdf" pg 68
*/
#define __ENABLE_LOOP_SUBD_PULL
bool EditableMeshSolidShape2::Subdivide( uint32 first_smooth_vid )
{
    GEO_ASSERT( !IsBeingEdited() );
    //GEO_LOG_WARNING("EditableMeshSolidShape2::Subdivide() BEGIN #P = %d, #V = %d, #HE = %d", GetNumP(), GetNumV(), GetNumHE() );

    // Ensure triangular faces
    for( unsigned int it_p=0; it_p<m_addP.size(); it_p++ )
        GEO_ASSERT( P_NumEdges(it_p) == 3 );

    BeginEdition();
    {
#ifdef __ENABLE_LOOP_SUBD_PULL //VERTICES
        // Smooth existing vertices (>= first_smooth_vid)
        std::vector<editable_vertex_type> vec_old_v(m_addV);
        std::swap( vec_old_v, m_addV );
        for( unsigned int it_v=first_smooth_vid; it_v<vec_old_v.size(); it_v++ )
        {
            uint32 k( vec_old_v[it_v].m_NumEdges );
            //GEO_LOG_WARNING("k[%d] = %d", it_v, k );
            GEO_ASSERT(k>1);
            Real w( mal::Rcp<Real>(k) * (Real(5)/8 - mal::Sq(Real(3)/8 + Real(0.25)*mal::Cos(mal::TwoPi<Real>()/k)) ) ); //\todo precompute w[k]!
            m_addV[it_v].m_Pos = (Real(1)-k*w)*vec_old_v[it_v].m_Pos;
            // We CANNOT use a baked-data iterator_polygons_around_vertex_ccw here
            uint32 it_heid_av( vec_old_v[it_v].m_OutHEID );
            do
            {
                // Get other vid ccw
                uint32 next_heid( m_addHE[it_heid_av].m_NextHEID );
                uint32 vid1( m_addHE[next_heid].m_OriginVID );
                m_addV[it_v].m_Pos += w*vec_old_v[vid1].m_Pos;
                // Get next edge ccw (see iterator_polygons_around_vertex_ccw)
                uint32 left_pid( m_addHE[it_heid_av].m_LeftPID );
                uint32 num_edges( m_addP[left_pid].m_NumEdges );
                uint32 heid0( m_addP[left_pid].m_FirstHEID );
                uint32 prev_heid( heid0 + ((it_heid_av - heid0 - 1 + num_edges) % num_edges ) );
                it_heid_av = m_addHE[prev_heid].m_SymHEID;
            } while( it_heid_av != vec_old_v[it_v].m_OutHEID     //internal V
                     && it_heid_av != cInvalidFeatureIndex ); //boundary V
            /*\todo Boundary vertices MUST be handled specifically, if
              I understand it correctly, they should behave as "crease
              vertices" as described in "Piecewise Smooth Surface
              Reconstruction". This requires

              By now we leave boundary vertices unchanged... which is
              what we want in subd-driven deformation
            */
            if( it_heid_av == cInvalidFeatureIndex )
                m_addV[it_v].m_Pos = vec_old_v[it_v].m_Pos;
        }
#else
        //TEMP; hack to avoid complicating edge subd code
        const std::vector<editable_vertex_type>& vec_old_v(m_addV);
#endif
        std::vector<feature_index_type> vec_hesd( m_addHE.size() );
        // Split all unique HE, storing mid-point VID
        for( unsigned int it_p=0; it_p<m_addP.size(); it_p++ )
        {
            editable_polygon_type& polygon = m_addP[it_p];
            unsigned int it_he0( polygon.m_FirstHEID );
            feature_index_type vid0( m_addHE[it_he0].m_OriginVID );
            do
            {
                unsigned int it_he1( m_addHE[it_he0].m_NextHEID );
                feature_index_type vid1( m_addHE[it_he1].m_OriginVID );
                feature_index_type sym_he0( m_addHE[it_he0].m_SymHEID );
#ifdef __ENABLE_LOOP_SUBD_PULL //EDGES
                if( cInvalidFeatureIndex == sym_he0 ) //Border edge
                    vec_hesd[it_he0] = AddVertex( Real(0.5f) * (vec_old_v[vid0].m_Pos + vec_old_v[vid1].m_Pos) );
                else if( it_he0 < sym_he0 ) //Regular edge
                {
                    feature_index_type vidL( m_addHE[ m_addHE[ it_he1 ].m_NextHEID ].m_OriginVID );
                    feature_index_type vidR( m_addHE[ m_addHE[ m_addHE[ sym_he0 ].m_NextHEID ].m_NextHEID ].m_OriginVID );
                    vec_hesd[it_he0] = AddVertex( Real(3.0f/8.0f) * (vec_old_v[vid0].m_Pos + vec_old_v[vid1].m_Pos)
                                                  + Real(1.0f/8.0f) * (vec_old_v[vidL].m_Pos + vec_old_v[vidR].m_Pos)  );
                    vec_hesd[sym_he0] = vec_hesd[it_he0];
                }
                //\todo We could ALSO displace vertices, but we do NOT
                //want to, as they will be given by the
                //simulation... we MAY BE ABLE TO skip vertex
                //refinement COMPLETELY for embedding purposes if we
                //ensure that the embedded Curve will be "far enough"
                //from any original vertex neigbourhood"
#else
                if( cInvalidFeatureIndex == sym_he0 || it_he0 < sym_he0 )
                {
                    vec_hesd[it_he0] = AddVertex( Real(0.5f) * (vec_old_v[vid0].m_Pos + vec_old_v[vid1].m_Pos) );
                    if( cInvalidFeatureIndex != sym_he0 ) vec_hesd[sym_he0] = vec_hesd[it_he0];
                }
#endif
                vid0 = vid1;
                it_he0 = it_he1;
            } while ( it_he0 != polygon.m_FirstHEID );
        }


        // Generate sub-faces from original faces and add their topology
        std::vector<editable_polygon_type> vec_old_polygons;
        std::swap( vec_old_polygons, m_addP ); //Implicitly removes all P
        std::vector<editable_half_edge_type> vec_old_he;
        std::swap( vec_old_he, m_addHE ); //Implicitly removes all HE

        // Add subd polygons and implicitly HE
        for( unsigned int it_p=0; it_p<vec_old_polygons.size(); it_p++ )
        {
            const editable_polygon_type& polygon = vec_old_polygons[it_p];
            //Gather 3 VID + 3 Split HE VID \note We INIT the following vecs because -O3 complains otherwise
            feature_index_type vec_vid[3] = { cInvalidFeatureIndex, cInvalidFeatureIndex, cInvalidFeatureIndex };
            feature_index_type vec_split_he_vid[3] = { cInvalidFeatureIndex, cInvalidFeatureIndex, cInvalidFeatureIndex };
            unsigned int it_he( polygon.m_FirstHEID );
            unsigned int it_eip(0);
            do
            {
                vec_vid[it_eip] = vec_old_he[it_he].m_OriginVID;
                vec_split_he_vid[it_eip] = vec_hesd[it_he];
                it_he = vec_old_he[it_he].m_NextHEID;
                it_eip++;
            } while ( it_he != polygon.m_FirstHEID );
            // Add 4 triangles
            AddPolygon3( vec_vid[0],          vec_split_he_vid[0], vec_split_he_vid[2] ); //0ac
            AddPolygon3( vec_split_he_vid[0], vec_vid[1],          vec_split_he_vid[1] ); //a1b
            AddPolygon3( vec_split_he_vid[2], vec_split_he_vid[1], vec_vid[2] );          //cb2
            AddPolygon3( vec_split_he_vid[0], vec_split_he_vid[1], vec_split_he_vid[2] ); //abc
        }
    }
    bool bResult = EndEdition();
    //GEO_LOG_WARNING("EditableMeshSolidShape2::Subdivide() END #P = %d, #V = %d, #HE = %d", GetNumP(), GetNumV(), GetNumHE() );

    //\todo for the LOOP_PUSH scheme, here we would perform a Smoothing pas where the effect of each T^k+1 on each V^k+1 would be accumulated

    return bResult;
}

/*Slow but easy:
  - Rebuild the whole EditableMeshSolidShape2 adding all existing
    polygons except the ones in vec_pid
  - Note that spurious vertices only adjacent to removed polygons may
    be added in the first pass, but will be removed in EndEdition()
    FixDegenerateVertices
  \todo Optimizations are possible but not required so far
*/
bool EditableMeshSolidShape2::RemovePolygons( uint32 num_pid, const feature_index_type* vec_pid )
{
    ClearEditData(); //we'll use baked data to re-add all edit data
    {
        // Re-add all Vertices
        for( unsigned int it_v=0; it_v<m_NumV; it_v++ )
            m_addV.push_back( editable_vertex_type( m_vecPoints[it_v], m_vecV[it_v].m_OutHEID, m_vecV[it_v].m_NumEdges ) );
        // Re-add non-removed Polygons
        for( unsigned int it_p=0; it_p<m_NumP; it_p++ )
        {
            bool bRemove(false);
            for( unsigned int it_remove_pid=0; !bRemove && it_remove_pid < num_pid; it_remove_pid++ )
                bRemove = it_p == vec_pid[it_remove_pid];
            if( !bRemove )
            {
                std::vector<feature_index_type> vec_vid;
                unsigned int num_vid(0);
                unsigned int it_he(m_vecP[it_p].m_FirstHEID);
                do
                {
                    vec_vid.push_back( m_vecHE[it_he].m_OriginVID );
                    it_he = m_vecHE[it_he].m_NextHEID;
                    num_vid++;
                } while ( it_he != m_vecP[it_p].m_FirstHEID );
                AddPolygonN( num_vid, &vec_vid[0] );
            }
        }
    }
    ClearBakedData(); //Clear baked data to re-bake it in EndEdition
    bool bResult = EndEdition();
    return bResult;
}

bool EditableMeshSolidShape2::AddDCR( const IShape2* p_shape, const Transform2& tr_s2mss )
{
    if( p_shape->GetType() == eShape_Polygonal2 )
    {
        if( m_pDCR ) delete m_pDCR;
        m_pDCR = 0;
        const PolygonalShape2& polygonal( *static_cast<const PolygonalShape2*>(p_shape) );
        SetBakedDCR_StrictlyNonshared_UglyHack( Create_DCR_MeshSolidShape2_From_PolygonalShape2( *this, Transform2::Identity(), GetVecDefaultSDOF(),
                                                                                                 polygonal, tr_s2mss, polygonal.GetVecDefaultSDOF() ) );
        return true;
    }
    else
        return false;
}

#ifdef __GEO_MSS_ENABLE_BVH
bool EditableMeshSolidShape2::AddBVH()
{
    if( m_pBVH ) delete m_pBVH;
    m_pBVH = new BVH_MeshSolidShape2;
    uint32 num_entries( m_pDCR ? m_pDCR->m_NumElements : m_NumP ); //TEMP: If there's a DCR, only use those P in the BVH
    m_pBVH->Rebuild_BottomUp( num_entries,
                              boost::bind<void>( &GEBV_MeshSolidShape2_E<BVH_MeshSolidShape2::entry_index_type,BVH_MeshSolidShape2::bv_type>,
                                                 this, transform_type::Identity(), GetVecDefaultSDOF(),
                                                 _1, _2 ) );
    return true;
}
#endif

//---- Internal methods
void EditableMeshSolidShape2::ClearEditData()
{
    m_addV.clear();
    m_addP.clear();
    m_addHE.clear();
    m_addBHE.clear();
    m_addBP.clear();
    m_addL.clear();
}

void EditableMeshSolidShape2::RebuildBakedData()
{
    // Bake Vertices
    Vec2 *vecPoints( new Vec2[m_addV.size()] );
    vertex_type *vecV( new vertex_type[m_addV.size()] );
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
    {
        vecPoints[it_v] = m_addV[it_v].m_Pos;
        vecV[it_v] = m_addV[it_v];
    }
    // Bake Polygons and BoundaryPolygons
    polygon_type *vecP( new polygon_type[m_addP.size() + m_addBP.size()] );
    for( unsigned int it_p=0; it_p<m_addP.size(); it_p++ )
        vecP[it_p] = m_addP[it_p];
    for( unsigned int it_bp=0; it_bp<m_addBP.size(); it_bp++ )
    {
        //GEO_LOG_WARNING( "Baking BP[%d]", it_bp );
        unsigned int it_p( m_addP.size() + it_bp );
        vecP[it_p] = m_addBP[it_bp];
        // Reindex Boundary features
        vecP[it_p].m_FirstHEID += m_addHE.size();
    }
    // Bake HalfEdges and BoundaryHalfEdges
    half_edge_type *vecHE( new half_edge_type[m_addHE.size() + m_addBHE.size()] );
    for( unsigned int it_he=0; it_he<m_addHE.size(); it_he++ )
        vecHE[it_he] = m_addHE[it_he];
    for( unsigned int it_bhe=0; it_bhe<m_addBHE.size(); it_bhe++ )
    {
        //GEO_LOG_WARNING( "Baking BHE[%d]", it_bhe );
        unsigned int it_he( m_addHE.size() + it_bhe );
        vecHE[it_he] = m_addBHE[it_bhe];
        // Reindex Boundary features
        vecHE[it_he].m_NextHEID += m_addHE.size();
        vecHE[it_he].m_LeftPID += m_addP.size();
#ifdef __ENABLE_LINK_BOUNDARY_SYMMETRIC_HE
        //vecHE[it_he].m_SymHEID += 0; \note Already normal index
        vecHE[ vecHE[it_he].m_SymHEID ].m_SymHEID += m_addHE.size();
#endif
    }

    polygon_layer_type* vecL( new polygon_layer_type[m_addL.size()] );
    for( unsigned int it_l=0; it_l<m_addL.size(); it_l++ )
        vecL[it_l] = m_addL[it_l];
    // Set baked data as non-shared
    SetBakedData( false,
                  m_addV.size(), m_addP.size(), m_addHE.size(),
                  m_addBP.size(), m_addBHE.size(), m_addL.size(),
                  vecPoints, vecV, vecP, vecHE, vecL );
    // Clear temporal stuff
    delete [] vecL;
    delete [] vecHE;
    delete [] vecP;
    delete [] vecV;
    delete [] vecPoints;
}

bool EditableMeshSolidShape2::FixDegeneracies()
{
    bool bFixDV = FixDegenerateVertices();
    bool bFixDE = FixDegenerateEdges();
    bool bFixDP = FixDegeneratePolygons();
    return bFixDV || bFixDE || bFixDP;
}

void EditableMeshSolidShape2::RebuildTopology()
{
    //\todo Use GSimpleSpatialHash, see TriSurfaceShape3 and TetSolidShape3

    // Init vid topology
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
    {
        m_addV[it_v].m_OutHEID = cInvalidFeatureIndex;
        m_addV[it_v].m_NumEdges = 0;
    }

    /*\note HE topology is already correct, built incrementally on
      addition
    */
    /* Link V to adjacent outgoing HE, ensuring that if V is on the
       perimeter, the first CCW outgoing HE will be selected, and
       iterator_polygons_around_vertex_ccw will successfully iterate all
       Polygons in CCW order starting from this HE.
    */
    for( unsigned int it_he=0; it_he < m_addHE.size(); it_he++ )
    {
        feature_index_type vid0 = m_addHE[it_he].m_OriginVID;
        if( m_addV[vid0].m_OutHEID == cInvalidFeatureIndex || m_addHE[it_he].m_SymHEID == cInvalidFeatureIndex )
            m_addV[vid0].m_OutHEID = it_he;
        // Increase directed vertex degree
        if( m_addHE[it_he].m_SymHEID == cInvalidFeatureIndex )
        {
            m_addV[vid0].m_NumEdges += 2;
            m_addV[m_addHE[m_addHE[it_he].m_NextHEID].m_OriginVID].m_NumEdges += 2;
        }
        else
        {
            m_addV[vid0].m_NumEdges++;
            m_addV[m_addHE[m_addHE[it_he].m_NextHEID].m_OriginVID].m_NumEdges++;
        }
    }
    // undirected-edge vertex degree = directed/=2
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        m_addV[it_v].m_NumEdges /= 2;

    /* Compute Boundary Polygons
       \note Bounday Polygons and HalfEdges are added in separate
       arrays, and any added feature is identified using
       boundary-local indexing starting at 0 for the first BPID or
       HEID, however, all referenced preexisting features are
       identified using normal-indexing. Boundary-local index
       variables have a "b" prefix (eg: bheid0, bpid...)
    */
    // Gather unassigned boundary edges
    std::vector<editable_half_edge_type> vec_UBHE;
    for( unsigned int it_he=0; it_he < m_addHE.size(); it_he++ )
        if( m_addHE[it_he].m_SymHEID == cInvalidFeatureIndex )
        {
            feature_index_type vid0 = m_addHE[it_he].m_OriginVID;
            feature_index_type vid1 = m_addHE[ m_addHE[it_he].m_NextHEID ].m_OriginVID;
            // Add BHE with Next as VID, not as HEID, and without a valid LeftPID
            //GEO_LOG_WARNING( "UBHE[%d] = (%d,%d)", vec_UBHE.size(), vid1, vid0 );
            vec_UBHE.push_back( editable_half_edge_type( vid1,
                                                         vid0,
                                                         it_he,
                                                         cInvalidFeatureIndex ) );
        }
    GEO_ASSERT( vec_UBHE.size() > 2 );

    // Connect boundary edges into polygons while there is unassigned BHE
    do
    {
        // Alloc new pid
        feature_index_type bpid( m_addBP.size() );
        feature_index_type bheid0( m_addBHE.size() );
        feature_index_type first_vid( vec_UBHE.front().m_OriginVID );
        feature_count_type num_edges(0);
        bool bClosed(false);
        // Add BHE
        do
        {
            // Use first UBHE
            const editable_half_edge_type& ubhe( vec_UBHE.front() );
            feature_index_type vid0( ubhe.m_OriginVID );
            feature_index_type vid1( ubhe.m_NextHEID ); //\note Here next is VID, not HEID!
#ifdef __ENABLE_LINK_BOUNDARY_SYMMETRIC_HE
            // Link symmetric HE to this new polygon and its HE
            feature_index_type sheid0( ubhe.m_SymHEID );
            if( sheid0 != cInvalidFeatureIndex ) m_addHE[sheid0].m_SymHEID = m_addBHE.size(); //boundary-local index
#endif
            // Add fully linked HE
            num_edges++;
            m_addBHE.push_back( editable_half_edge_type( vid0,
                                                         bheid0 + num_edges, // Wrong for last edge, we'll fix it when bClosed
#ifdef __ENABLE_LINK_BOUNDARY_SYMMETRIC_HE
                                                         sheid0,
#else
                                                         cInvalidFeatureIndex, //\todo ??sheid0,
#endif
                                                         bpid ) );

            // Remove first UBHE
            //GEO_LOG_WARNING( "Consumed UBHE (%d,%d)", vid0, vid1 );
            vec_UBHE.front() = vec_UBHE.back();
            vec_UBHE.pop_back();
            // Advance or Close
            if( vid1 == first_vid )
            {
                //GEO_LOG_WARNING( "Closing with UBHE (%d,%d)...", vid0, vid1 );
                bClosed = true;
            }
            else
            {
                //GEO_LOG_WARNING( "Searching UBHE (%d,?)...", vid1 );
                // Find Next UBHE, put it in the Front of vec_UBHE
                unsigned int next_ubhe_index(0);
                bool bFound(false);
                while( !bFound && next_ubhe_index < vec_UBHE.size() )
                {
                    if( vec_UBHE[next_ubhe_index].m_OriginVID == vid1 )
                    {
                        //GEO_LOG_WARNING( "...Next UBHE (%d,%d)", vec_UBHE[next_ubhe_index].m_NextHEID, vec_UBHE[next_ubhe_index].m_NextHEID );
                        std::swap( vec_UBHE.front(), vec_UBHE[next_ubhe_index] );
                        bFound = true;
                    }
                    else
                        next_ubhe_index++;
                }
                //\note Failure implies open or non-manifold boundary (a vid appears more than once in the boundary)
                GEO_ASSERT( bFound );
            }
        } while( !bClosed );
        // Fix last bhe m_NextHEID
        m_addBHE.back().m_NextHEID = bheid0;
        // Create polygon
        m_addBP.push_back( editable_polygon_type( bheid0, num_edges ) );
    } while( vec_UBHE.size() > 0 );
}

/* Rebuild layers from m_addX data, assuming Topology is up-to-date
 */
void EditableMeshSolidShape2::RebuildLayers()
{
    const unsigned int cMaxLayers = m_addP.size();
    //1) Compute per-vertex layer from boundary
    std::vector<unsigned int> vecVertexLayer( m_addV.size(), cMaxLayers );
    // Set boundary V layers to 0
    for( unsigned int it_bhe=0; it_bhe<m_addBHE.size(); it_bhe++ )
        vecVertexLayer[ m_addBHE[it_bhe].m_OriginVID ] = 0;
    // Propagate layer from boundary (\todo O(MaxLayers*E) iteration could be optimized by O(V+E) flooding)
    for( unsigned int i=0; i<cMaxLayers; i++ )
        for( unsigned int it_he=0; it_he<m_addHE.size(); it_he++ )
            vecVertexLayer[ m_addHE[it_he].m_OriginVID ] = mal::Min( vecVertexLayer[ m_addHE[it_he].m_OriginVID ],
                                                                     vecVertexLayer[ m_addHE[ m_addHE[it_he].m_NextHEID ].m_OriginVID ]+1 );
    //2) Compute per-polygon layer from min adjacent vertex layer
    std::vector<unsigned int> vecPolygonLayer( m_addP.size(), cMaxLayers );
    for( unsigned int it_p=0; it_p<m_addP.size(); it_p++ )
    {
        editable_polygon_type& polygon( m_addP[it_p] );
        unsigned int it_he( polygon.m_FirstHEID );
        do
        {
            vecPolygonLayer[it_p] = mal::Min( vecPolygonLayer[it_p], vecVertexLayer[ m_addHE[it_he].m_OriginVID ] );
            it_he = m_addHE[it_he].m_NextHEID;
        } while ( it_he != polygon.m_FirstHEID );
    }
    //3) Sort polygons per layer and generate array of layer descriptors [0..num_layers-1] \todo BOUNDARY P do not change!!
    unsigned int max_layer(0);
    for( unsigned int it_p=0; it_p<m_addP.size(); it_p++ )
        if( max_layer < vecPolygonLayer[it_p] )
            max_layer = vecPolygonLayer[it_p];
    unsigned int num_layers(max_layer+1);
    // Add Polygons to layer-lists
    std::vector< std::vector< feature_index_type > > vecLayers( num_layers );
    for( unsigned int it_p=0; it_p<m_addP.size(); it_p++ )
        vecLayers[ vecPolygonLayer[it_p] ].push_back( it_p );
    // Build contiguous layers and compute PID mapping
    m_addL.resize( num_layers ); //sets all layers to 0 elements
    std::vector< feature_index_type > vecOld2NewPID( m_addP.size() );
    std::vector< editable_polygon_type > addP_sorted;
    unsigned int last_pid(0);
    for( unsigned int it_l=0; it_l<vecLayers.size(); it_l++ )
    {
        //GEO_LOG_WARNING( "Layer[%d] F%d,S%ld", it_l, last_pid, vecLayers[it_l].size() );
        m_addL[ it_l ].m_FirstPID = last_pid;
        m_addL[ it_l ].m_NumPolygons = vecLayers[it_l].size();
        for( unsigned int it_pil=0; it_pil<vecLayers[it_l].size(); it_pil++ )
        {
            vecOld2NewPID[ vecLayers[it_l][it_pil] ] = last_pid++;
            addP_sorted.push_back( m_addP[ vecLayers[it_l][it_pil] ] );
            //GEO_LOG_WARNING( "%d => %d", vecLayers[it_l][it_pil], last_pid-1 );
        }
    }
    // Swap unsorted with sorted
    std::swap( m_addP, addP_sorted );
    //4) Remap HE.m_LeftPID
    for( unsigned int it_he=0; it_he<m_addHE.size(); it_he++ )
        m_addHE[it_he].m_LeftPID = vecOld2NewPID[ m_addHE[it_he].m_LeftPID ];

    /*\note We do NOT remap BHE.m_LeftPID, as BP have their own array
      m_addBP that does NOT change due to layer reordering. */

    //GEO_LOG_WARNING( "RebuildLayers() => %d layers", num_layers );

    /*OPTIONAL (to improve cache behaviour)
      5) Sort V to match polygon layer order and remap
      6) Sort HE to match polygon order and remap
      \todo This is NOT necessary for derived TriSolidShape2 that does
            not store edges explicitly, so better spend time writing
            TriSolidShape2 instead of losing it doing 5),6) that are
            complex and error-prone and not that useful once
            TriSolidShape2 is available.
    */
}

bool EditableMeshSolidShape2::FixDegenerateVertices()
{
    //\todo Weld coincident vertices, see TriSurfaceShape3

    // Compute V degree (topology is NOT available, only m_addP data is up to date)
    std::vector< int > vecVertexDegree( m_addV.size(), 0 );
    for( unsigned int it_p=0; it_p<m_addP.size(); it_p++ )
    {
        editable_polygon_type& polygon( m_addP[it_p] );
        unsigned int it_he( polygon.m_FirstHEID );
        do
        {
            vecVertexDegree[ m_addHE[it_he].m_OriginVID ]++;
            it_he = m_addHE[it_he].m_NextHEID;
        } while ( it_he != polygon.m_FirstHEID );
    }
    // Remove V with 0-degree and remap to their new index
    std::vector< feature_index_type > vecOld2NewVID( m_addV.size() );
    unsigned int num_removed(0);
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
    {
        if( vecVertexDegree[it_v] > 0 )
        {
            vecOld2NewVID[it_v] = it_v - num_removed;
            m_addV[it_v-num_removed] = m_addV[it_v];
        }
        else
        {
            vecOld2NewVID[it_v] = cInvalidFeatureIndex;
            num_removed++;
        }
    }
    m_addV.resize( m_addV.size() - num_removed ); //crop last num_removed V
    // Remap P vid => nothing to do
    // Remap HE vid
    for( unsigned int it_he=0; it_he < m_addHE.size(); it_he++ )
        m_addHE[it_he].m_OriginVID = vecOld2NewVID[ m_addHE[it_he].m_OriginVID ];

    // GEO_LOG_WARNING("%d 0-degree vertices removed",num_removed);
    return num_removed > 0;
}

bool EditableMeshSolidShape2::FixDegenerateEdges() { return false; }
bool EditableMeshSolidShape2::FixDegeneratePolygons() { return false; }

feature_index_type EditableMeshSolidShape2::FindV( const Vec2& pos, Real epsilon_sq ) const
{
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        if( mal::NormSq( m_addV[it_v].m_Pos - pos ) < epsilon_sq )
            return feature_index_type(it_v);
    return cInvalidFeatureIndex;
}

feature_index_type EditableMeshSolidShape2::FindHE( feature_index_type vid0, feature_index_type vid1 ) const
{
    for( unsigned int it_he=0; it_he<m_addHE.size(); it_he++ )
        if( m_addHE[it_he].m_OriginVID == vid0
            && m_addHE[ m_addHE[it_he].m_NextHEID ].m_OriginVID == vid1 )
            return feature_index_type(it_he);
    return cInvalidFeatureIndex;
}

//-----------------------------------------------------------------------------
//---- Make_MeshSolidShape2_XXX Implementation
//-----------------------------------------------------------------------------

void Make_MeshSolidShape2_Box( EditableMeshSolidShape2& emss, const Vec2& half_sizes )
{
    emss.Clear();
    emss.BeginEdition();
    {
        feature_index_type vid0 = emss.AddVertex( Vec2(-half_sizes.x(),-half_sizes.y()) );
        feature_index_type vid1 = emss.AddVertex( Vec2( half_sizes.x(),-half_sizes.y()) );
        feature_index_type vid2 = emss.AddVertex( Vec2( half_sizes.x(), half_sizes.y()) );
        feature_index_type vid3 = emss.AddVertex( Vec2(-half_sizes.x(), half_sizes.y()) );
        emss.AddPolygon3( vid0, vid1, vid2 );
        emss.AddPolygon3( vid2, vid3, vid0 );
    }
    emss.EndEdition();
}

void Make_MeshSolidShape2_BoxWithHole( EditableMeshSolidShape2& emss, const Vec2& half_sizes, const Vec2& hole_half_sizes )
{
    emss.Clear();
    emss.BeginEdition();
    {
        feature_index_type vid0 = emss.AddVertex( Vec2(-half_sizes.x(),-half_sizes.y()) );
        feature_index_type vid1 = emss.AddVertex( Vec2( half_sizes.x(),-half_sizes.y()) );
        feature_index_type vid2 = emss.AddVertex( Vec2( half_sizes.x(), half_sizes.y()) );
        feature_index_type vid3 = emss.AddVertex( Vec2(-half_sizes.x(), half_sizes.y()) );

        feature_index_type vid4 = emss.AddVertex( Vec2(-hole_half_sizes.x(),-hole_half_sizes.y()) );
        feature_index_type vid5 = emss.AddVertex( Vec2( hole_half_sizes.x(),-hole_half_sizes.y()) );
        feature_index_type vid6 = emss.AddVertex( Vec2( hole_half_sizes.x(), hole_half_sizes.y()) );
        feature_index_type vid7 = emss.AddVertex( Vec2(-hole_half_sizes.x(), hole_half_sizes.y()) );

        emss.AddPolygon3( vid0, vid1, vid4 );
        emss.AddPolygon3( vid1, vid5, vid4 );

        emss.AddPolygon3( vid1, vid2, vid5 );
        emss.AddPolygon3( vid2, vid6, vid5 );

        emss.AddPolygon3( vid2, vid3, vid6 );
        emss.AddPolygon3( vid3, vid7, vid6 );

        emss.AddPolygon3( vid3, vid0, vid7 );
        emss.AddPolygon3( vid0, vid4, vid7 );
    }
    emss.EndEdition();
}

void Make_MeshSolidShape2_Box( EditableMeshSolidShape2& emss, const Vec2& half_sizes, unsigned int num_x, unsigned int num_y )
{
    GEO_ASSERT( !mal::IsNaN(half_sizes) );
    GEO_ASSERT( num_x > 1 && num_y > 1 );
    emss.Clear();
    emss.BeginEdition();
    {
        Real step_x( 2 * half_sizes.x() / (num_x-1) );
        Real step_y( 2 * half_sizes.y() / (num_y-1) );
        for( unsigned int it_y=0; it_y < num_y; it_y++ )
            for( unsigned int it_x=0; it_x < num_x; it_x++ )
                emss.AddVertex( -half_sizes + Vec2(it_x*step_x,it_y*step_y) );
        for( unsigned int it_y=0; it_y < num_y-1; it_y++ )
            for( unsigned int it_x=0; it_x < num_x-1; it_x++ )
            {
                feature_index_type vid0( it_y*num_x     + it_x   );
                feature_index_type vid1( it_y*num_x     + it_x+1 );
                feature_index_type vid2( (it_y+1)*num_x + it_x+1 );
                feature_index_type vid3( (it_y+1)*num_x + it_x   );
                emss.AddPolygon3( vid0, vid1, vid2 );
                emss.AddPolygon3( vid2, vid3, vid0 );
            }
    }
    emss.EndEdition();
}

//----------------------------------------------------------------
// Create DCR \todo move elsewhere
//----------------------------------------------------------------

/* Build the DCR for a PolygonalShape2 embedded in a MeshSolidShape2
   \pre The polygonal MUST be closed
   \pre The polygonal MUST be completely enclosed in the MSS2 Layer[0]

   \todo In order to AVOID storing MSS element indices by making the
   m_vecE array parallel to MSS.layer[0] element array, we SHOULD BE
   ABLE TO REORDER MSS layer[0] here, as it has NO ORDERING GUARANTEES.
   => ALTERNATIVELY, we could store the element indices during
   element/patch construction, but reorder them at bake-time to MATCH
   the EXISTING element order in the MSS, which would remain unchanged.
   => THIS is A LOT BETTER, as it avoids a two-way dependencty
   between the DCR and the MSS
   - We'll need to ENSURE that all MSS.layer[0] elements have a
   DCR.Patch associated, and that no non-MSS.layer[0] element
   has it, otherwise DCR element array CANNOT be parallel to
   layer[0], thus, we need to CROP Ext elements before
   computing the DCR.
   => By now we relaxed this so that the DCR uses any number of
      elements with parallel indices to MSS, but some may have no
      associated patch, some may have > 1.
*/

// per-patch and per-element data
struct dcr2_patch_id_type
{
    feature_index_type m_EID; //Element id
    int m_PIE; //Patch-in-element
    dcr2_patch_id_type() : m_EID(cInvalidFeatureIndex), m_PIE(-1) {}
    dcr2_patch_id_type( feature_index_type eid, int pie ) : m_EID(eid), m_PIE(pie) {}
};
struct dcr2_patch_geometry_type
{
    dcr2_patch_geometry_type() {}
    std::vector< Vec2 > m_vecPoints;
    dcr2_patch_id_type m_Next, m_Prev;
};
struct dcr2_element_geometry_type
{
    dcr2_element_geometry_type() {}
    std::vector< dcr2_patch_geometry_type > m_vecPatches;
};

/* Check if point is inside an mss polygon.
   \pre Polygon must be a triangle by now, but the test would work for
   any convex polygon.
*/
bool IsPointInsideMeshSolidShape2_Polygon( const Vec2& point,
                                           const MeshSolidShape2& mesh, const Transform2& mesh_tr, const Vec2* vec_mesh_sdof,
                                           feature_index_type pid )
{
    GEO_ASSERT( mesh.P_NumEdges(pid) == 3 );
    unsigned int it_he( mesh.P_FirstHEID(pid) );
    Vec2 e0( mesh_tr * mesh.V_Pos( mesh.HE_OriginVID(it_he), vec_mesh_sdof ) );
    do
    {
        Vec2 e1( mesh_tr * mesh.V_Pos( mesh.HE_FinalVID(it_he), vec_mesh_sdof ) );
        Vec2 n( mal::PerpendicularCW(e1-e0) );
        if( mal::Dot( point-e0, n ) < 0 ) return false;
        it_he = mesh.HE_Next(it_he);
        e0 = e1;
    } while ( it_he != mesh.P_FirstHEID(pid) );
    return true;
}

DCR_MeshSolidShape2* Create_DCR_MeshSolidShape2_From_PolygonalShape2( const MeshSolidShape2& mesh, const Transform2& mesh_tr, const Vec2* vec_mesh_sdof,
                                                                      const PolygonalShape2& polygonal, const Transform2& polygonal_tr, const Vec2* vec_polygonal_sdof )
{
    GEO_ASSERT( polygonal.IsClosed() );

    DCR_MeshSolidShape2* pDCR = new DCR_MeshSolidShape2();

    Transform2 inv_mesh_tr(mal::Inverse(mesh_tr));

    // All elements are initially empty
    std::vector< dcr2_element_geometry_type > vecElements; //SAME order as MSS elements
    vecElements.resize( mesh.GetNumP() );

    // Find Element that contains first polygonal vertex
    Vec2 p0 = polygonal_tr * polygonal.V_Pos( 0, vec_polygonal_sdof );
    unsigned int first_eid(cInvalidFeatureIndex);
    for( unsigned int it_p = 0; it_p<mesh.GetNumP() && first_eid == cInvalidFeatureIndex; it_p++ )
    {
        if( IsPointInsideMeshSolidShape2_Polygon( p0, mesh, mesh_tr, vec_mesh_sdof, it_p ) )
            first_eid = it_p;
    }
    GEO_ASSERT(first_eid != cInvalidFeatureIndex);
    unsigned int max_eid( first_eid );

    // Create first DCR patch
    /*\note This patch may need to be extended "backwards", V[0] will
     rarely be exactly on an element boundary. Instead of searching
     backwards here, we'll merge the first and last patches if they
     are in the same element.
    */
    dcr2_patch_id_type first_patch_id( first_eid, 0 );
    dcr2_patch_id_type current_patch_id = first_patch_id;
    vecElements[current_patch_id.m_EID].m_vecPatches.push_back( dcr2_patch_geometry_type() );
    vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_vecPoints.push_back( inv_mesh_tr*p0 );
    pDCR->m_NumPatches++;
    pDCR->m_NumVertices++;

    // Iterate over Polygonal edges, clipping against MSS elements and tracking "current element", adding contiguous patches (Vertices,Segments)
    unsigned int it_polygonal_edge = 0;
    feature_index_type last_crossed_mesh_heid = cInvalidFeatureIndex;
    while( it_polygonal_edge < polygonal.GetNumVertices() ) //\pre polygonal.IsClosed()
    {
        // Current polygonal segment is (p0,p1)
        Vec2 p1 = polygonal_tr * polygonal.V_Pos( (it_polygonal_edge+1) % polygonal.GetNumVertices(), vec_polygonal_sdof );

        // Check mesh element boundary features for crossing with (p0,p1)
        float lambda_p(0),lambda_e(0);
        feature_index_type next_eid( cInvalidFeatureIndex );
        unsigned int it_he( mesh.P_FirstHEID(current_patch_id.m_EID) );
        Vec2 e0( mesh_tr * mesh.V_Pos( mesh.HE_OriginVID(it_he), vec_mesh_sdof ) );
        do
        {
            /*\todo As we add q to BOTH the current patch (p0,q) and
              the next candidate (q,p1), the Symm(E) of crossed edge
              E that yields q should NOT BE CONSIDERED for crossing
              in (q,p1), as it may intersect (q,p1) due to numerical
              prec
            */
            Vec2 e1( mesh_tr * mesh.V_Pos( mesh.HE_FinalVID(it_he), vec_mesh_sdof ) );
            if( mesh.HE_Sym( it_he ) != last_crossed_mesh_heid )
            {
                if( np::Intersection_Segment2_Segment2(p0,p1,e0,e1,lambda_p,lambda_e) )
                {
                    next_eid = mesh.HE_RightPID( it_he );
                    last_crossed_mesh_heid = it_he;
                }
            }
            it_he = mesh.HE_Next(it_he);
            e0 = e1;
        } while ( it_he != mesh.P_FirstHEID(current_patch_id.m_EID)
                  && next_eid == cInvalidFeatureIndex );

        // Move to next patch or remain in the current one
        if( next_eid != cInvalidFeatureIndex )
        {
            // Find clipped point
            Vec2 q( p0 + lambda_p*(p1-p0) );
            // Add q to current patch
            vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_vecPoints.push_back( inv_mesh_tr*q );
            // Create next DCR patch and add q as first point
            dcr2_patch_id_type next_patch_id = dcr2_patch_id_type( next_eid, vecElements[next_eid].m_vecPatches.size() );
            vecElements[next_patch_id.m_EID].m_vecPatches.push_back( dcr2_patch_geometry_type() );
            vecElements[next_patch_id.m_EID].m_vecPatches[next_patch_id.m_PIE].m_vecPoints.push_back( inv_mesh_tr*q );
            // Connect current/next patches
            vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_Next = next_patch_id;
            vecElements[next_patch_id.m_EID].m_vecPatches[next_patch_id.m_PIE].m_Prev = current_patch_id;
            // Move to next patch/element
            current_patch_id = next_patch_id;
            max_eid = mal::Max<uint32>( max_eid, next_eid );
            pDCR->m_NumPatches++;
            pDCR->m_NumSegments++;
            pDCR->m_NumVertices += 2;
            // Advance segment
            p0 = q;
        }
        else
        {
            // Add p1 to current patch
            vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_vecPoints.push_back( inv_mesh_tr*p1 );
            // Advance to next polygonal edge
            p0 = p1;
            it_polygonal_edge++;
            pDCR->m_NumSegments++;
            pDCR->m_NumVertices++;
            /*\note A patch CAN enter and exit an Element E through
             the same boundary feature, but NOT in a single segment,
             therefore, once we find an internal V, we reset the
             last_crossed_mesh_heid to avoid any restrictions on
             exiting E, otherwise some valid crossings are NOT
             detected!
            */
            last_crossed_mesh_heid = cInvalidFeatureIndex;
        }
    }

    // Check last/first patch connection
    if( first_patch_id.m_EID == current_patch_id.m_EID
        && first_patch_id.m_PIE != current_patch_id.m_PIE )
    {
        // Merge first and last patch if different but in the same element \pre IsClosed()
        // Relink first/last
        dcr2_patch_id_type last_patch_id = vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_Prev;
        vecElements[first_patch_id.m_EID].m_vecPatches[first_patch_id.m_PIE].m_Prev = last_patch_id;
        vecElements[last_patch_id.m_EID].m_vecPatches[last_patch_id.m_PIE].m_Next = first_patch_id;
        // Remove last vertex
        vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_vecPoints.pop_back();
        pDCR->m_NumVertices--;
        // Append last to first
        vecElements[first_patch_id.m_EID].m_vecPatches[first_patch_id.m_PIE].m_vecPoints.insert( vecElements[first_patch_id.m_EID].m_vecPatches[first_patch_id.m_PIE].m_vecPoints.begin(),
                                                                                                 vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_vecPoints.begin(),
                                                                                                 vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_vecPoints.end() );
        // Remove last patch, already subsumed in the first one
        vecElements[current_patch_id.m_EID].m_vecPatches.pop_back();
        pDCR->m_NumPatches--;
    }
    else
    {
        // Connect last added patch with the first one \pre IsClosed()
        vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_Next = first_patch_id;
        vecElements[first_patch_id.m_EID].m_vecPatches[first_patch_id.m_PIE].m_Prev = current_patch_id;
    }

    //\todo If the polygonal is closed, the sequence of Elements and
    //Patches will be too. This WILL BE REQUIRED to build cross-patch
    //topology in the baked representation. If NOT closed, the first
    //element and last elements will contain an extreme of the
    //polygonal.

    // \todo Ideally the mesh should be cropped so that only layer[0]
    // elements contain patches, but we'll tolerate otherwise by now
    pDCR->m_NumElements = max_eid+1;
    if( pDCR->m_NumElements != mesh.L_NumP(0) )
        GEO_LOG_WARNING( "Create_DCR_MeshSolidShape2_From_PolygonalShape2() DCR uses <= %d elements and MSS layer[0] has %d",
                         pDCR->m_NumElements, mesh.L_NumP(0) );

    // Alloc DCR arrays
    pDCR->m_vecE = new DCR_MeshSolidShape2::Element[ pDCR->m_NumElements ];
    pDCR->m_vecP = new DCR_MeshSolidShape2::Patch[ pDCR->m_NumPatches ];
    pDCR->m_vecS = new DCR_MeshSolidShape2::Segment[ pDCR->m_NumSegments ];
    pDCR->m_vecV = new Vec2[ pDCR->m_NumVertices ];

    /*Baking:
      for each element in mesh that contains patches
        for each patch
           add vertices to global array
           add segments to global array
           add patch
      for each patch in geometric order
        bake cross-patch topology into existing segments
    */
    unsigned int last_patch_index(0);
    unsigned int last_segment_index(0);
    unsigned int last_vertex_index(0);
    for( unsigned int it_element=0; it_element<pDCR->m_NumElements; it_element++ )
    {
        DCR_MeshSolidShape2::Element& ed( pDCR->m_vecE[it_element] );
        ed.m_FirstPID = last_patch_index;
        ed.m_NumPatches = 0;
        ed.m_FirstVID = last_vertex_index;
        ed.m_NumVertices = 0;
        for( int i=0; i<3; i++ ) ed.m_BDOP[i] = Interval::Empty();
        ed.m_BDOP_BestSlabIdx = 0;
        // Gather element nodes
        uint32 vec_nid[3];
        mesh.P_VecVID( it_element, vec_nid, 3 );
        Vec2 vec_element_nodes[3] = { mesh.V_Pos_0( vec_nid[0] ), mesh.V_Pos_0( vec_nid[1] ), mesh.V_Pos_0( vec_nid[2] ) };
        // Compute element barycentric coordinates matrix and its 3x2 inverse
        Mat3x3 Bs( 1, 1, 1,
                   vec_element_nodes[0].x(), vec_element_nodes[1].x(), vec_element_nodes[2].x(),
                   vec_element_nodes[0].y(), vec_element_nodes[1].y(), vec_element_nodes[2].y() );
        Mat3x3 invBs( mal::Inverse( Bs ) );
        for( unsigned int it_patch=0; it_patch<vecElements[it_element].m_vecPatches.size(); it_patch++ )
        {
            DCR_MeshSolidShape2::Patch& pd( pDCR->m_vecP[last_patch_index] );
            pd.m_EID = it_element;
            pd.m_FirstSID = last_segment_index;
            pd.m_NumSegments = 0;
            // Add patch vertices and segments
            for( unsigned int it_point=0; it_point<vecElements[it_element].m_vecPatches[it_patch].m_vecPoints.size(); it_point++ )
            {
                // Add vertex
                pDCR->m_vecV[last_vertex_index] = vecElements[it_element].m_vecPatches[it_patch].m_vecPoints[it_point];
                // Add segment
                if( it_point > 0 )
                {
                    DCR_MeshSolidShape2::Segment& pst( pDCR->m_vecS[last_segment_index] );
                    // vertices
                    pst.SetVID(0,last_vertex_index-1);
                    pst.SetVID(1,last_vertex_index);
                    // neighbours
                    if( it_point > 1 ) pst.SetNSID(0,last_segment_index-1); //internal segment
                    else pst.SetNSID(0,cInvalidFeatureIndex); //\note will be last SID in previous patch
                    if( it_point < vecElements[it_element].m_vecPatches[it_patch].m_vecPoints.size()-1 ) pst.SetNSID(1,last_segment_index+1); //internal segment
                    else pst.SetNSID(1,cInvalidFeatureIndex); //\note will be first SID in next patch
                    last_segment_index++;
                    pd.m_NumSegments++;
                }
                // WIP: Init or grow BDOP
                {
                    // Compute barycentric coords
                    Vec3 b = invBs * mal::Concat(1,pDCR->m_vecV[last_vertex_index]);
                    /*\note As we add b[i] to the BDOP, geometrically
                       the BDOP.Min() plane will be the *farthest*
                       from v[i], and BDOP.Max() will be the *closest*
                    */
                    if( ed.m_NumVertices == 0 )
                    {
                        //Init BDOP
                        ed.m_BDOP[0].Set( b[0] );
                        ed.m_BDOP[1].Set( b[1] );
                        ed.m_BDOP[2].Set( b[2] );
                    }
                    else
                    {
                        // Grow BDOP
                        ed.m_BDOP[0].Merge( b[0] );
                        ed.m_BDOP[1].Merge( b[1] );
                        ed.m_BDOP[2].Merge( b[2] );
                    }
                }
                last_vertex_index++;
                ed.m_NumVertices++;
            }
            last_patch_index++;
            ed.m_NumPatches++;
        }
        /*TEMP
        GEO_LOG( "ED[%d] BDOP = (%f, %f) (%f,%f) (%f,%f)",
                 it_element,
                 ed.m_BDOP[0].Min(), ed.m_BDOP[0].Max(),
                 ed.m_BDOP[1].Min(), ed.m_BDOP[1].Max(),
                 ed.m_BDOP[2].Min(), ed.m_BDOP[2].Max() );
        */
        if( !ed.m_BDOP[0].IsEmpty()
            && !ed.m_BDOP[1].IsEmpty()
            && !ed.m_BDOP[2].IsEmpty() )
        {
            // Select best BDOP => B-Slab with smallest volume
            geo::Vec2 bdop_pos[4];
            geo::Real area_bdop[3] = {mal::Infinity<geo::Real>(),mal::Infinity<geo::Real>(),mal::Infinity<geo::Real>()};
            ed.m_BDOP_BestSlabIdx = 0;
            for( int i=0; i<3; i++ )
            {
                bdop_pos[0] = ed.m_BDOP[i].Min() * vec_element_nodes[i] + (1-ed.m_BDOP[i].Min()) * vec_element_nodes[(i+1)%3];
                bdop_pos[1] = ed.m_BDOP[i].Min() * vec_element_nodes[i] + (1-ed.m_BDOP[i].Min()) * vec_element_nodes[(i+2)%3];
                bdop_pos[2] = ed.m_BDOP[i].Max() * vec_element_nodes[i] + (1-ed.m_BDOP[i].Max()) * vec_element_nodes[(i+1)%3];
                bdop_pos[3] = ed.m_BDOP[i].Max() * vec_element_nodes[i] + (1-ed.m_BDOP[i].Max()) * vec_element_nodes[(i+2)%3];
                // Trapezoid area = h*(min(l1,l2) + 1/2*abs(l1-l2)) \todo Consider simpler formula in TetSolidShape3 (Trapezoid area = Big + 0.5*(Big-Small))
                geo::Real h = mal::Abs( mal::Dot( bdop_pos[0]-bdop_pos[2],
                                                  mal::Normalized( mal::PerpendicularCCW(bdop_pos[0]-bdop_pos[1]) ) ) );
                area_bdop[i] = h * ( mal::Min( mal::Norm(bdop_pos[0]-bdop_pos[1]),
                                               mal::Norm(bdop_pos[2]-bdop_pos[3]) )
                                     + 0.5f*mal::Abs( mal::Norm(bdop_pos[0]-bdop_pos[1]) - mal::Norm(bdop_pos[2]-bdop_pos[3]) ) );
                if( area_bdop[i] < area_bdop[ed.m_BDOP_BestSlabIdx] )
                    ed.m_BDOP_BestSlabIdx = i;
            }
        }
        else //\todo This should never happen if all DCR.E contain >0 DCR.V, but MAY HAPPEN as DCR.E are NOT NECESSARILY layer[0] (should crop ext elements as we do in 3D)
        {
            GEO_LOG_WARNING("Create_DCR_MeshSolidShape2_From_PolygonalShape2() DCR.E[%u] has empty BDOP, setting to whole DCR.E", it_element );
            ed.m_BDOP[0] = Interval(0,1);
            ed.m_BDOP[1] = Interval(0,1);
            ed.m_BDOP[2] = Interval(0,1);
            ed.m_BDOP_BestSlabIdx = 0;
        }
    }
    GEO_ASSERT( last_patch_index == pDCR->m_NumPatches );
    GEO_ASSERT( last_segment_index == pDCR->m_NumSegments );
    GEO_ASSERT( last_vertex_index == pDCR->m_NumVertices );
    // Bake cross-patch topology
    // IMPORTANT At this point, SID are expected to be CCW-ordered inside their DCR.P
    current_patch_id = first_patch_id;
    do
    {
        dcr2_patch_id_type next_patch_id = vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_Next;
        // Get current and next patches
        DCR_MeshSolidShape2::Element& current_ed( pDCR->m_vecE[current_patch_id.m_EID] );
        DCR_MeshSolidShape2::Patch& current_pd( pDCR->m_vecP[current_ed.m_FirstPID + current_patch_id.m_PIE] );
        DCR_MeshSolidShape2::Element& next_ed( pDCR->m_vecE[next_patch_id.m_EID] );
        DCR_MeshSolidShape2::Patch& next_pd( pDCR->m_vecP[next_ed.m_FirstPID + next_patch_id.m_PIE] );
        // Connect last/first segments of current/next patches
        DCR_MeshSolidShape2::Segment& current_last_sd( pDCR->m_vecS[current_pd.m_FirstSID + current_pd.m_NumSegments - 1] );
        DCR_MeshSolidShape2::Segment& next_first_sd( pDCR->m_vecS[next_pd.m_FirstSID] );
        current_last_sd.SetNSID(1,next_pd.m_FirstSID);
        next_first_sd.SetNSID(0,current_pd.m_FirstSID + current_pd.m_NumSegments - 1);
        // Connect DCR.P NPID and SID
        dcr2_patch_id_type prev_patch_id = vecElements[current_patch_id.m_EID].m_vecPatches[current_patch_id.m_PIE].m_Prev;
        DCR_MeshSolidShape2::Element& prev_ed(pDCR->m_vecE[prev_patch_id.m_EID]);
        current_pd.m_vecNPID[0] = prev_ed.m_FirstPID + prev_patch_id.m_PIE;
        current_pd.m_vecNPID[1] = next_ed.m_FirstPID + next_patch_id.m_PIE;
        current_pd.m_vecSID[0] = current_pd.m_FirstSID;
        current_pd.m_vecSID[1] = current_pd.m_FirstSID + current_pd.m_NumSegments - 1;
        // Advance
        current_patch_id = next_patch_id;
    } while( current_patch_id.m_EID != first_patch_id.m_EID
             || current_patch_id.m_PIE != first_patch_id.m_PIE );


    /*TODO: NOO! WE EXPECT per-patch segments to be
     * CONSECUTIVE in MANY PLACES, but here we are shuffling
     * them at the DCR.E level, not the DCR.P, thus, we're
     * MIXING the patches...
     There is 2 solutions
     - Store per-patch BST and keep DCR.P SID consecutive
     - Store per-element BST and BREAK consecutive DCR.P SID, which FREES the DCR.P.m_FirstSID variable but requires complex iteration over DCR.P.Segments EVERY-FUCKING-WHERE... this is a no go, specially in 3D
     => SO, per-patch BST
    */

#ifdef __ENABLE_MSS_DCR_PATCH_BDT
    // Bake per-patch BDT
    std::vector<uint32> mapSID( pDCR->m_NumSegments, 0xFFFFFFFF );
    uint32 processed_sid(0);
    for( unsigned int it_e=0; it_e<pDCR->m_NumElements; it_e++ )
    {
        const DCR_MeshSolidShape2::Element& ed( pDCR->m_vecE[it_e] );
        // Get inverse barycentric transform
        uint32 vec_nid[3];
        mesh.P_VecVID( it_e, vec_nid, 3 );
        Vec2 vec_pos[3] = { mesh.V_Pos_0( vec_nid[0] ), mesh.V_Pos_0( vec_nid[1] ), mesh.V_Pos_0( vec_nid[2] ) };
        Mat3x3 invBm( mal::Inverse( Mat3x3( 1, 1, 1,
                                            vec_pos[0].x(), vec_pos[1].x(), vec_pos[2].x(),
                                            vec_pos[0].y(), vec_pos[1].y(), vec_pos[2].y() ) ) ); //\todo invBm could be precomputed in DCR.ED!!
        // Compute per-patch BDT
        for( unsigned int it_pie=0; it_pie<ed.m_NumPatches; it_pie++ )
        {
            uint32 pid( ed.m_FirstPID + it_pie );
            const DCR_MeshSolidShape2::Patch& pd( pDCR->m_vecP[pid] );

            typedef std::pair<DCR_MeshSolidShape2::Segment,uint32> patch_segment_type;
            std::vector< patch_segment_type > vecPS;

            // Gather Patch-Segment sub-array with their current SID
            for( unsigned int it_sip=0; it_sip<pd.m_NumSegments; it_sip++ )
            {
                uint32 sid( pd.m_FirstSID + it_sip );
                vecPS.push_back( std::make_pair(pDCR->m_vecS[sid],sid) );
                processed_sid++;
            }
            GEO_ASSERT(!vecPS.empty());

            // Generate BDT
            // GEO_LOG( "**** Creating DCR.P[%u].BDT with %u DCR.S", pid, pd.m_NumSegments );
            std::vector<DCR_MeshSolidShape2::Patch::bdt_node_sip_range> stackBDTN;
            stackBDTN.push_back( DCR_MeshSolidShape2::Patch::bdt_node_sip_range(0,vecPS.size()) ); //root-range
            unsigned int num_bdtn(0);
            while( !stackBDTN.empty() )
            {
                // Pop PS subarray
                DCR_MeshSolidShape2::Patch::bdt_node_sip_range bdtn_sr( stackBDTN.back() );
                stackBDTN.pop_back();

                // Compute BDOP for current PS subarray
                bv::BDOP3 bdop;
                for( unsigned int it_sin=bdtn_sr.first; it_sin<bdtn_sr.second; it_sin++ )
                {
                    //Get the original SID for a PS in the sub-array and use it to retrieve its geometry
                    uint32 sid( vecPS[it_sin].second );
                    const DCR_MeshSolidShape2::Segment& sd( pDCR->m_vecS[sid] );
                    Vec3 vec_b[2] = { invBm * mal::Concat(1,pDCR->m_vecV[sd.GetVID(0)]),
                                      invBm * mal::Concat(1,pDCR->m_vecV[sd.GetVID(1)]) };
                    for( int i=0; i<2; i++ ) //vertices
                        for( int j=0; j<3; j++ ) //bdop-axis
                            bdop[j].Merge( mal::Clamp01(vec_b[i][j]) );
                }
                // Select longest axis \todo Consider better heuristics
                unsigned int best_axis(0);
                for( unsigned int i=1; i<3; i++ )
                    if( bdop[i].Length() > bdop[best_axis].Length() )
                        best_axis = i;

                // Sort PS subarray along chosen axis
                std::sort( vecPS.begin() + bdtn_sr.first, vecPS.begin() + bdtn_sr.second,
                           [pDCR, invBm, best_axis](const patch_segment_type& ps1, const patch_segment_type& ps2)
                           {
                               Vec3 ps1_b( invBm * mal::Concat(1,
                                                               Real(0.5) * (pDCR->m_vecV[pDCR->m_vecS[ps1.second].GetVID(0)]
                                                                            + pDCR->m_vecV[pDCR->m_vecS[ps1.second].GetVID(1)]) ) );
                               Vec3 ps2_b( invBm * mal::Concat(1,
                                                               Real(0.5) * (pDCR->m_vecV[pDCR->m_vecS[ps2.second].GetVID(0)]
                                                                            + pDCR->m_vecV[pDCR->m_vecS[ps2.second].GetVID(1)]) ) );
                               return ps1_b[best_axis] < ps2_b[best_axis];
                               // return ps1.second > ps2.second;
                           }
                    );

                // Create BDTN for current PS subarray at the first PS
                DCR_MeshSolidShape2::Segment& bdtn( vecPS[bdtn_sr.first].first );
                bdtn.BDTN_Init( bdop );
                // TEMP: DEBUG
                // {
                //     //TEMP: Log node
                //     BDOP3 bdop_r( bdtn.BDTN_BDOPq() );
                //     DCR_MeshSolidShape2::Segment::BDOP3q bdop_q( bdtn.BDTN_BDOPq() );
                //     GEO_LOG("BDTN[%u] best_axis %u, range [%u,%u]", num_bdtn, best_axis, bdtn_sr.first, bdtn_sr.second );
                //     GEO_LOG("- BDOP [%f,%f]x[%f,%f]x[%f,%f]", bdop[0].Min(), bdop[0].Max(), bdop[1].Min(), bdop[1].Max(), bdop[2].Min(), bdop[2].Max() );
                //     GEO_LOG("- BDOPr [%f,%f]x[%f,%f]x[%f,%f]", bdop_r[0].Min(), bdop_r[0].Max(), bdop_r[1].Min(), bdop_r[1].Max(), bdop_r[2].Min(), bdop_r[2].Max() );
                //     GEO_LOG("- BDOPq [%d,%d]x[%d,%d]x[%d,%d]", bdop_q[0].Min(), bdop_q[0].Max(), bdop_q[1].Min(), bdop_q[1].Max(), bdop_q[2].Min(), bdop_q[2].Max() );
                //     //TEMP: Check conservative quantization
                //     for( int i=0; i<3; i++ )
                //         GEO_ASSERT( bdop_r[i].Min() <= bdop[i].Min() && bdop_r[i].Max() >= bdop[i].Max() );
                // }
                num_bdtn++;

                // Skip the BDTN PS and split the remaining PS subarray at the median
                int length_sr( bdtn_sr.second - bdtn_sr.first );
                int remaining_sr( length_sr - 1 );
                if( remaining_sr > 1 )
                {
                    // GEO_LOG("Stacking 2 L/R");
                    stackBDTN.push_back( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( bdtn_sr.first+1, bdtn_sr.first+1+remaining_sr/2 ) );
                    stackBDTN.push_back( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( bdtn_sr.first+1+remaining_sr/2, bdtn_sr.second ) );
                }
                else if( remaining_sr == 1 )
                {
                    // GEO_LOG("Stacking 1 L");
                    stackBDTN.push_back( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( bdtn_sr.first+1, bdtn_sr.first+2 ) );
                }
            }

            //TEMP: TEST ordering backwards to check correct remapping, works
            // std::sort( vecPS.begin(), vecPS.end(),
            //            [](const patch_segment_type& a, const patch_segment_type& b){ return a.second > b.second; } );

            // Bake Element Segment sub-array into global segment sub-array and fill mapSID
            for( unsigned int it_sip=0; it_sip<vecPS.size(); it_sip++ )
            {
                uint32 old_sid( vecPS[it_sip].second );
                uint32 new_sid( pd.m_FirstSID + it_sip );
                pDCR->m_vecS[new_sid] = vecPS[it_sip].first;
                mapSID[old_sid] = new_sid;
                processed_sid++;
            }
        }
    }

    /*TEMP: Check correct remapping
    GEO_ASSERT( processed_sid == 2*pDCR->m_NumSegments );
    uint32 num_unmapped_sid(0);
    for( unsigned int it_s=0; it_s<pDCR->m_NumSegments; it_s++ )
        if( mapSID[it_s] == 0xFFFFFFFF )
            num_unmapped_sid++;
    GEO_LOG_ERROR_IF( num_unmapped_sid > 0, "NumUnmappedSID = %u", num_unmapped_sid );
    GEO_ASSERT( num_unmapped_sid == 0 );
    */

    // Remap SID in DCR.P.m_vecSID and DCR.S.m_vecNSID
    for( unsigned int it_p=0; it_p<pDCR->m_NumPatches; it_p++ )
    {
        pDCR->m_vecP[it_p].m_vecSID[0] = mapSID[ pDCR->m_vecP[it_p].m_vecSID[0] ];
        pDCR->m_vecP[it_p].m_vecSID[1] = mapSID[ pDCR->m_vecP[it_p].m_vecSID[1] ];
    }
    for( unsigned int it_s=0; it_s<pDCR->m_NumSegments; it_s++ )
    {
        pDCR->m_vecS[it_s].SetNSID( 0, mapSID[ pDCR->m_vecS[it_s].GetNSID(0) ] );
        pDCR->m_vecS[it_s].SetNSID( 1, mapSID[ pDCR->m_vecS[it_s].GetNSID(1) ] );
    }

    /*TEMP: DEBUG LOG
    for( unsigned int it_p=0; it_p < pDCR->m_NumPatches; it_p++ )
    {
        uint32 num_visited_nodes(0);
        const geo::DCR_MeshSolidShape2::Patch& pd( pDCR->m_vecP[it_p] );
        typedef std::pair< DCR_MeshSolidShape2::Patch::bdt_node_sip_range,uint32> stack_entry_type;
        std::vector< stack_entry_type > stackBDTN;
        stackBDTN.push_back( stack_entry_type(DCR_MeshSolidShape2::Patch::bdt_node_sip_range(0,pd.m_NumSegments),0) );
        GEO_LOG("DCR.P[%u].BDT",it_p);
        while( !stackBDTN.empty() )
        {
            num_visited_nodes++;
            // Pop PS subarray
            stack_entry_type se( stackBDTN.back() );
            stackBDTN.pop_back();
            DCR_MeshSolidShape2::Patch::bdt_node_sip_range node_sr( se.first );
            // Get node
            const DCR_MeshSolidShape2::Segment& bdtn( pDCR->m_vecS[pd.m_FirstSID+node_sr.first] );
            BDOP3 bdop( bdtn.BDTN_BDOPq() );
            GEO_LOG("BDTN level %u, range [%u,%u], volume %f", se.second, node_sr.first, node_sr.second, bdop[0].Length()*bdop[1].Length()*bdop[2].Length() );
            GEO_LOG("- BDOP [%f,%f]x[%f,%f]x[%f,%f]", bdop[0].Min(), bdop[0].Max(), bdop[1].Min(), bdop[1].Max(), bdop[2].Min(), bdop[2].Max() );
            // Recurse
            int length_sr( node_sr.second - node_sr.first );
            int remaining_sr( length_sr - 1 );
            if( remaining_sr > 1 )
            {
                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1, node_sr.first+1+remaining_sr/2 ), se.second+1 ) );
                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1+remaining_sr/2, node_sr.second ), se.second+1) );
            }
            else if( remaining_sr == 1 )
                stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1, node_sr.first+2 ), se.second+1) );
        }
        GEO_LOG_ERROR_IF( num_visited_nodes != pd.m_NumSegments, "Recursive visit skips some nodes %u < %u!!", num_visited_nodes, pd.m_NumSegments );
    }
    */

#endif //__ENABLE_MSS_DCR_PATCH_BDT

    return pDCR;
}

/* Given a DCR (level=0), refine it by subdivision up to max_level
  For each element T
    For each level k
      compute T^k and V^k
      classify DCR.V into T^k
      compute subsets of TT^k and VV^k strictly required to embed DCR.V and bake them
      - Baking should compute:
        - The fastest representation that allows computing TT^k and VV^k (=> See Loop subd support, ideally FIXED and SMALL)
        - The relationship DCR.V in TT^k
          - 1 TT per DCR.V, we can use array of ids, probably uint8 = 256 T
            - 2D: 1Tri -> 4tri subd, 256 = 4^4, supporting k<=4, L0=1tri, L1=4, L2=16, L3=64, L4=256
            - 3D: 1Tet -> 4tet+1oct, 1oct->6tet+8oct=14 => L0=1tet, L1=4tet+1oct=5, L2=4*(4tet+1oct) + 1*(6tet+8oct)=22tet+12oct=34, L3=22*(4+1)+12*(6+8)=160tet+118oct=278 >256

  \note For each level k all data in k-1 is assumed available, both
  precomputed and at runtime

  \note For element T, its subd DEPENDS on its 1-ring neighbours
  {T_i}, which may or may not be DCR T (can be internal, Layer[1]
  elements). In order to keep the Subd-DCR decoupled from the MSS, we
  should save 1-ring neighbourhood for DCR.T, which has variable
  cardinality and no specific ordering (in 2D could be CCW, but in 3D
  there's no unique sequential order)

  \note The BV(Element|BDOP|BSlab) of each element T depends on its
  geometry and on its 1-ring neighbours. For each geometry
  Tri/Tet/BDOP/BSlab, try to find tight bounds as in AABB|KDOP(BDOP)
  paper, using restricted convex combinations, BUT using T nodes and
  its 1-ring as input data.

  \note Final application of subd embedding should use barycentric a
  matrix-palette for the TT^k (mult by pose matrices, precomputed) and
  each DCR.V would only need to be multiplied by its containing TT
  matrix.

  Alternatively:
  - Following 2004_SmoothSubdivisionOfTetrahedralMeshes.pdf,
    precompute the list of influences and weights for DCR.V^\infinity,
    which due to local support should be relatively short (1-ring
    size).
  - The resulting baked data would probably be more compact/direct,
    and even faster to access/operate at runtime, BUT we'd lose the
    nested descriptions of the DCR-embedding geometry, which we MAY
    need in the future to allow adaptive subdivision of the simulation
    meshes.
*/
void Make_DCR_Subd()// Basic_DCR, uint32 max_level, Subd_DCR )
{
//     for( int i=0; i<100; i-- )
//     {
// //        int i=0;
//     }
}

#ifdef __GEO_MSS_ENABLE_BVH
/* Specialization of BV(BDOP) with tight AABB(BDOP) computation based
   on the papers:
   - [2007] Adaptive Deformations with Fast Tight Bounds
   - [2008] Tight and efficient surface bounds in meshless animation
*/
template <>
void GEBV_MeshSolidShape2_BDOP( const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                feature_index_type eid, bv::AABB2& bv )
{
    GEO_ASSERT( 0 != p_dcr && eid >= 0 && eid < p_dcr->m_NumElements );
    const DCR_MeshSolidShape2::Element& ed( p_dcr->m_vecE[eid] );
    // compute global node pos
    uint32 vec_nid[3];
    p_mesh->P_VecVID( eid, vec_nid, 3 );
    Vec2 node_pos[3] = { tr*vec_node_pos[ vec_nid[0] ], tr*vec_node_pos[ vec_nid[1] ], tr*vec_node_pos[ vec_nid[2] ] };
    // sort node indices along X,Y
    int snie_x[3] = { 0, 1, 2 };
    if( node_pos[snie_x[0]].x() > node_pos[snie_x[1]].x() ) std::swap( snie_x[0], snie_x[1] );
    if( node_pos[snie_x[1]].x() > node_pos[snie_x[2]].x() ) std::swap( snie_x[1], snie_x[2] );
    if( node_pos[snie_x[0]].x() > node_pos[snie_x[1]].x() ) std::swap( snie_x[0], snie_x[1] );
    int snie_y[3] = { 0, 1, 2 };
    if( node_pos[snie_y[0]].y() > node_pos[snie_y[1]].y() ) std::swap( snie_y[0], snie_y[1] );
    if( node_pos[snie_y[1]].y() > node_pos[snie_y[2]].y() ) std::swap( snie_y[1], snie_y[2] );
    if( node_pos[snie_y[0]].y() > node_pos[snie_y[1]].y() ) std::swap( snie_y[0], snie_y[1] );
    // compute min/max \todo Consider manual optimization
    bv.SetMinMax( Vec2( ed.m_BDOP[snie_x[0]].Max()*node_pos[snie_x[0]].x() + ed.m_BDOP[snie_x[2]].Min()*node_pos[snie_x[2]].x() + (1-ed.m_BDOP[snie_x[0]].Max()-ed.m_BDOP[snie_x[2]].Min()) * node_pos[snie_x[1]].x(),
                        ed.m_BDOP[snie_y[0]].Max()*node_pos[snie_y[0]].y() + ed.m_BDOP[snie_y[2]].Min()*node_pos[snie_y[2]].y() + (1-ed.m_BDOP[snie_y[0]].Max()-ed.m_BDOP[snie_y[2]].Min()) * node_pos[snie_y[1]].y() ),
                  Vec2( ed.m_BDOP[snie_x[2]].Max()*node_pos[snie_x[2]].x() + ed.m_BDOP[snie_x[0]].Min()*node_pos[snie_x[0]].x() + (1-ed.m_BDOP[snie_x[2]].Max()-ed.m_BDOP[snie_x[0]].Min()) * node_pos[snie_x[1]].x(),
                        ed.m_BDOP[snie_y[2]].Max()*node_pos[snie_y[2]].y() + ed.m_BDOP[snie_y[0]].Min()*node_pos[snie_y[0]].y() + (1-ed.m_BDOP[snie_y[2]].Max()-ed.m_BDOP[snie_y[0]].Min()) * node_pos[snie_y[1]].y() ) );
}

/* Specialization of BV(BDOP) with tight DOP2_K8(BDOP) computation based
   on the papers:
   - [2007] Adaptive Deformations with Fast Tight Bounds
   - [2008] Tight and efficient surface bounds in meshless animation
*/
template <>
void GEBV_MeshSolidShape2_BDOP( const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                feature_index_type eid, bv::DOP2_K8& bv )
{
    GEO_ASSERT( 0 != p_dcr && eid >= 0 && eid < p_dcr->m_NumElements );
    const DCR_MeshSolidShape2::Element& ed( p_dcr->m_vecE[eid] );

    // compute global node pos
    uint32 vec_nid[3];
    p_mesh->P_VecVID( eid, vec_nid, 3 );
    Vec2 node_pos[3] = { tr*vec_node_pos[ vec_nid[0] ], tr*vec_node_pos[ vec_nid[1] ], tr*vec_node_pos[ vec_nid[2] ] };
    // sort node indices along X,Y
    int snie_x[3] = { 0, 1, 2 };
    if( node_pos[snie_x[0]].x() > node_pos[snie_x[1]].x() ) std::swap( snie_x[0], snie_x[1] );
    if( node_pos[snie_x[1]].x() > node_pos[snie_x[2]].x() ) std::swap( snie_x[1], snie_x[2] );
    if( node_pos[snie_x[0]].x() > node_pos[snie_x[1]].x() ) std::swap( snie_x[0], snie_x[1] );
    int snie_y[3] = { 0, 1, 2 };
    if( node_pos[snie_y[0]].y() > node_pos[snie_y[1]].y() ) std::swap( snie_y[0], snie_y[1] );
    if( node_pos[snie_y[1]].y() > node_pos[snie_y[2]].y() ) std::swap( snie_y[1], snie_y[2] );
    if( node_pos[snie_y[0]].y() > node_pos[snie_y[1]].y() ) std::swap( snie_y[0], snie_y[1] );
    // compute axis 0,1 = X,Y
    bv.SetInterval<0>(
        Interval( ed.m_BDOP[snie_x[0]].Max()*node_pos[snie_x[0]].x() + ed.m_BDOP[snie_x[2]].Min()*node_pos[snie_x[2]].x() + (1-ed.m_BDOP[snie_x[0]].Max()-ed.m_BDOP[snie_x[2]].Min()) * node_pos[snie_x[1]].x(),
                  ed.m_BDOP[snie_x[2]].Max()*node_pos[snie_x[2]].x() + ed.m_BDOP[snie_x[0]].Min()*node_pos[snie_x[0]].x() + (1-ed.m_BDOP[snie_x[2]].Max()-ed.m_BDOP[snie_x[0]].Min()) * node_pos[snie_x[1]].x() ) );
    bv.SetInterval<1>(
        Interval( ed.m_BDOP[snie_y[0]].Max()*node_pos[snie_y[0]].y() + ed.m_BDOP[snie_y[2]].Min()*node_pos[snie_y[2]].y() + (1-ed.m_BDOP[snie_y[0]].Max()-ed.m_BDOP[snie_y[2]].Min()) * node_pos[snie_y[1]].y(),
                  ed.m_BDOP[snie_y[2]].Max()*node_pos[snie_y[2]].y() + ed.m_BDOP[snie_y[0]].Min()*node_pos[snie_y[0]].y() + (1-ed.m_BDOP[snie_y[2]].Max()-ed.m_BDOP[snie_y[0]].Min()) * node_pos[snie_y[1]].y() ) );
    // Compute remaining axis 2,3 = (1,1),(1,-1) //\todo Could use a for() loop instead, see DOP3_K14
    Real node_proj_u[3];
    Real node_proj_v[3];
    for( int i=0; i<3; i++ )
    {
        node_proj_u[i] = bv::GKDOP_Dot<2,8,2>( node_pos[i] );
        node_proj_v[i] = bv::GKDOP_Dot<2,8,3>( node_pos[i] );
    }
    // sort node indices along axis 2,3
    int snie_u[3] = { 0, 1, 2 };
    if( node_proj_u[snie_u[0]] > node_proj_u[snie_u[1]] ) std::swap( snie_u[0], snie_u[1] );
    if( node_proj_u[snie_u[1]] > node_proj_u[snie_u[2]] ) std::swap( snie_u[1], snie_u[2] );
    if( node_proj_u[snie_u[0]] > node_proj_u[snie_u[1]] ) std::swap( snie_u[0], snie_u[1] );
    int snie_v[3] = { 0, 1, 2 };
    if( node_proj_v[snie_v[0]] > node_proj_v[snie_v[1]] ) std::swap( snie_v[0], snie_v[1] );
    if( node_proj_v[snie_v[1]] > node_proj_v[snie_v[2]] ) std::swap( snie_v[1], snie_v[2] );
    if( node_proj_v[snie_v[0]] > node_proj_v[snie_v[1]] ) std::swap( snie_v[0], snie_v[1] );
    // compute min/max \todo Consider manual optimization
    bv.SetInterval<2>(
        Interval( ed.m_BDOP[snie_u[0]].Max()*node_proj_u[snie_u[0]] + ed.m_BDOP[snie_u[2]].Min()*node_proj_u[snie_u[2]] + (1-ed.m_BDOP[snie_u[0]].Max()-ed.m_BDOP[snie_u[2]].Min()) * node_proj_u[snie_u[1]],
                  ed.m_BDOP[snie_u[2]].Max()*node_proj_u[snie_u[2]] + ed.m_BDOP[snie_u[0]].Min()*node_proj_u[snie_u[0]] + (1-ed.m_BDOP[snie_u[2]].Max()-ed.m_BDOP[snie_u[0]].Min()) * node_proj_u[snie_u[1]] ) );
    bv.SetInterval<3>(
        Interval( ed.m_BDOP[snie_v[0]].Max()*node_proj_v[snie_v[0]] + ed.m_BDOP[snie_v[2]].Min()*node_proj_v[snie_v[2]] + (1-ed.m_BDOP[snie_v[0]].Max()-ed.m_BDOP[snie_v[2]].Min()) * node_proj_v[snie_v[1]],
                  ed.m_BDOP[snie_v[2]].Max()*node_proj_v[snie_v[2]] + ed.m_BDOP[snie_v[0]].Min()*node_proj_v[snie_v[0]] + (1-ed.m_BDOP[snie_v[2]].Max()-ed.m_BDOP[snie_v[0]].Min()) * node_proj_v[snie_v[1]] ) );
}

/* Specialization of BV(BDOP) for Sphere(BDOP)
   \todo Fallback to BV(E)
*/
template <>
void GEBV_MeshSolidShape2_BDOP( const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                feature_index_type eid, bv::Sphere2& bv )
{
    GEBV_MeshSolidShape2_E<feature_index_type,bv::Sphere2>( p_mesh, tr, vec_node_pos, eid, bv );
}

/* Specialization of BV(BSlab) for AABB2(BSlab)
   \todo Most code could be shared with DOP2_K8(BSlab)!!
*/
template <>
void GEBV_MeshSolidShape2_BSlab( const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                 feature_index_type eid, bv::AABB2& bv )
{
    GEO_ASSERT( 0 != p_dcr && eid >= 0 && eid < p_dcr->m_NumElements );
    const DCR_MeshSolidShape2::Element& ed( p_dcr->m_vecE[eid] );
    // Compute global node pos
    uint32 vec_nid[3];
    p_mesh->P_VecVID( eid, vec_nid, 3 );
    Vec2 node_pos[3] = { tr*vec_node_pos[ vec_nid[0] ], tr*vec_node_pos[ vec_nid[1] ], tr*vec_node_pos[ vec_nid[2] ] };
    // Compute B-Slab points
    Vec2 bslab_pos[4];
    int i( ed.m_BDOP_BestSlabIdx );
    Vec2 node_pos_i_times_min( ed.m_BDOP[i].Min() * node_pos[i] );
    Real one_minus_min( 1 - ed.m_BDOP[i].Min() );
    Vec2 node_pos_i_times_max( ed.m_BDOP[i].Max() * node_pos[i] );
    Real one_minus_max( 1 - ed.m_BDOP[i].Max() );
    bslab_pos[0] = node_pos_i_times_min + one_minus_min * node_pos[(i+1)%3];
    bslab_pos[1] = node_pos_i_times_min + one_minus_min * node_pos[(i+2)%3];
    bslab_pos[2] = node_pos_i_times_max + one_minus_max * node_pos[(i+1)%3];
    bslab_pos[3] = node_pos_i_times_max + one_minus_max * node_pos[(i+2)%3];
    // Build AABB \todo Optimize manually
    bv = bv::AABB2(bslab_pos[0]).Merge(bslab_pos[1]).Merge(bslab_pos[2]).Merge(bslab_pos[3]);
}

/* Specialization of BV(BSlab) for DOP2_K8(BSlab)
   \todo Most code could be shared with AABB2(BSlab)!!
*/
template <>
void GEBV_MeshSolidShape2_BSlab( const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                 feature_index_type eid, bv::DOP2_K8& bv )
{
    GEO_ASSERT( 0 != p_dcr && eid >= 0 && eid < p_dcr->m_NumElements );
    const DCR_MeshSolidShape2::Element& ed( p_dcr->m_vecE[eid] );
    // compute global node pos
    uint32 vec_nid[3];
    p_mesh->P_VecVID( eid, vec_nid, 3 );
    Vec2 node_pos[3] = { tr*vec_node_pos[ vec_nid[0] ], tr*vec_node_pos[ vec_nid[1] ], tr*vec_node_pos[ vec_nid[2] ] };
    // Compute B-Slab points
    Vec2 bslab_pos[4];
    int i( ed.m_BDOP_BestSlabIdx );
    /*TEMP: unoptimized code
    bslab_pos[0] = ed.m_BDOP[i].Min() * node_pos[i] + (1-ed.m_BDOP[i].Min()) * node_pos[(i+1)%3];
    bslab_pos[1] = ed.m_BDOP[i].Min() * node_pos[i] + (1-ed.m_BDOP[i].Min()) * node_pos[(i+2)%3];
    bslab_pos[2] = ed.m_BDOP[i].Max() * node_pos[i] + (1-ed.m_BDOP[i].Max()) * node_pos[(i+1)%3];
    bslab_pos[3] = ed.m_BDOP[i].Max() * node_pos[i] + (1-ed.m_BDOP[i].Max()) * node_pos[(i+2)%3];
    */
    Vec2 node_pos_i_times_min( ed.m_BDOP[i].Min() * node_pos[i] );
    Real one_minus_min( 1 - ed.m_BDOP[i].Min() );
    Vec2 node_pos_i_times_max( ed.m_BDOP[i].Max() * node_pos[i] );
    Real one_minus_max( 1 - ed.m_BDOP[i].Max() );
    bslab_pos[0] = node_pos_i_times_min + one_minus_min * node_pos[(i+1)%3]; //\todo %3 should be optimized too...
    bslab_pos[1] = node_pos_i_times_min + one_minus_min * node_pos[(i+2)%3];
    bslab_pos[2] = node_pos_i_times_max + one_minus_max * node_pos[(i+1)%3];
    bslab_pos[3] = node_pos_i_times_max + one_minus_max * node_pos[(i+2)%3];
    // Build DOP2_K8 \todo Optimize manually
    bv = bv::DOP2_K8(bslab_pos[0]).Merge(bslab_pos[1]).Merge(bslab_pos[2]).Merge(bslab_pos[3]);
}

/* Specialization of BV(BSlab) for Sphere(BSlab)
   \todo Fallback to BV(E)
*/
template <>
void GEBV_MeshSolidShape2_BSlab( const MeshSolidShape2* p_mesh, const DCR_MeshSolidShape2* p_dcr, const Transform2& tr, const MeshSolidShape2::sdof_type* vec_node_pos,
                                 feature_index_type eid, bv::Sphere2& bv )
{
    GEBV_MeshSolidShape2_E<feature_index_type,bv::Sphere2>( p_mesh, tr, vec_node_pos, eid, bv );
}

#endif //__GEO_MSS_ENABLE_BVH

/* IMPORTANT!! I expect cost KDOP(BSlab) << KDOP(BDOP) with similar volume, for BSlab to be worth it....
       Cost analysis:
       - 2007: \todo BEST SORT(3) and SORT(4) alg??
         - 2D AABB
           = {x,y} * ( {a,b} * {2+ + 2- + 3*}  + SORT(3) )
             => BUBBLE_SORT(3) = 3*2 iter * {1< + 1SWAP} = 3*2*4 = 24
             => SWAP = 3=
           = 2*( 2*6 + 24) )
           = 72 OPS \todo AMb SORT(3) especific SEGUR que es molt mes rapid!!
         - 2D DOP8
           = 4*( 2*6 + 24 )
           = 144 OPS \todo NO! falta cost de projectar nodes en els axis no trivials!!
         - 3D AABB
           = {x,y,z} * ( {a,b} * {3+ + 3- + 4*}  + SORT(4) )
             => BUBBLE_SORT(4) = 4*3 iter * {1< + 1SWAP} = 4*3*4 = 48
             => MERGE_SORT(4) = 2< + 2SWAP + 1< + 1= + 1< + 1= + 2= ==> 2 + 2*3 + 6 = 14
           = 3*( 2*10 + 48 )
           = 198 OPS
           ==> 3*(2*10+14) = 102 OPS amb MERGE_SORT(4)
         - 3D DOP18
           = 9*(2*10 + 14) = 306 OPS \todo NO! falta cost de projectar nodes en els axis no trivials!!

         GENERAL:
         ==> K/2 * ( 2 * (D*sum + D*sub + (D+1)*mul) + Project(D+1) + SORT(D+1) )

       - best_slab:
         - 2D AABB
           = {p0,p1,p2,p3} * {x,y} * {1+ + 1- + 2*} + {x,y} * { MIN(4) + MAX(4) }
             => MIN(4) == MAX(4) = 1= + 3 iter * ( 1< + 1= ) => 1+3*2 = 7
             => MinMax(4) = Cas pitjor == MIN(4) + MAX(4) = 14
           = (4*2*4 + 2*2*7)
           = 60 OPS
         - 2D DOP8
           = (4*2*4 + 4*7)
           = 90 OPS \todo NO! falta cost de projectar p0..p3 en els axis no trivials!!
         - 3D AABB
           = {p0,p1,p2,p3,p4,p5} * {x,y,z} * {1+ + 1- + 2*} + {x,y,z} * { MIN(6) + MAX(6) }
             => MIN(6) == MAX(6) = 1= + 5 iter * ( 1< + 1= ) => 1+5*2 = 11
             => MinMax(6) = Cas pitjor == Min + Max = 22
           = 6*3*4 + 3*2*11
           = 105 OPS
         - 3D DOP18
           => 6*3*4 + 9*2*11 = 270 OPS  \todo NO! falta cost de projectar p0..p5 en els axis no trivials!!

         GENERAL:
         ==> 2*D * D * (sum + sub + mul) + K/2 * ( Project(2*D) + MinMax(2*D) )
         ===> Usant temporals per guardar els termes repetits (1 punt obtingut amb Dmul + 1 float obtingut amb 1sub) per TOTS ELS PUNTS de cada slab
         ====> 2*D * ( D * (sum + mul) ) + 2*(sub + D*mul) + K/2 * ( Project(2*D) + MinMax(2*D) )
         => Project() costa 0 per als primers D eixos, 1sum+1mul per la resta d'eixos en un DOP2K8 i DOP2K18 DOP2K18 (eixos (+-1|0,+-1|0,+-1|0) i i 2sum+1mul per als DOP2K14(eixos (+-1,+-1,+-1))

         MinMax4(a,b,c,d)
         {
           min = a; max = a;
           if( b < min ) min = b;
           else if( b > max ) max = b;
           if( c < min ) min = c;
           else if( c > max ) max = c;
           if( d < min ) min = d;
           else if( d > max ) max = d;
         }
    */

} //namespace geo
