#include "TriangleElementSet2.h"
#include "Rotation.h"
#include <Mal/GSolvePolynomialEq.h>
#include <memory.h>

//#define __S2_MS_ENABLE_TRACE
#ifdef __S2_MS_ENABLE_TRACE
#include <util/SimpleTracer.h>
#include <string>
#include <stdio.h>
#define S2_MS_TRACE_INIT() UTIL_TRACE_INIT()
#define S2_MS_TRACE_SHUTDOWN() UTIL_TRACE_SHUTDOWN()
#define S2_MS_TRACE_BEGIN_SLICE(x) UTIL_TRACE_BEGIN_SLICE(x)
#define S2_MS_TRACE_END_SLICE(b_flush) UTIL_TRACE_END_SLICE(b_flush)
#define S2_MS_TRACE_BEGIN_SCOPE(name) UTIL_TRACE_BEGIN_SCOPE(name)
#define S2_MS_TRACE_END_SCOPE() UTIL_TRACE_END_SCOPE()
#define S2_MS_TRACE_LOCAL(name,value) UTIL_TRACE_LOCAL(name,value)
#define S2_MS_TRACE_LOCAL_ARRAY(name,vec_value,num_values) UTIL_TRACE_LOCAL_ARRAY(name,vec_value,num_values)
#define S2_MS_TRACE_GLOBAL(name,value) UTIL_TRACE_GLOBAL(name,value)
#define S2_MS_TRACE_GLOBAL_ARRAY(name,vec_value,num_values) UTIL_TRACE_GLOBAL_ARRAY(name,vec_value,num_values)
#else
#define S2_MS_TRACE_INIT()
#define S2_MS_TRACE_SHUTDOWN()
#define S2_MS_TRACE_BEGIN_SLICE(x)
#define S2_MS_TRACE_END_SLICE(b_flush)
#define S2_MS_TRACE_BEGIN_SCOPE(name)
#define S2_MS_TRACE_END_SCOPE()
#define S2_MS_TRACE_LOCAL(name,value)
#define S2_MS_TRACE_LOCAL_ARRAY(name,vec_value,num_values)
#define S2_MS_TRACE_GLOBAL(name,value)
#define S2_MS_TRACE_GLOBAL_ARRAY(name,vec_value,num_values)
#endif

#define __ENABLE_PLASTIC_STRAIN_DEVIATION //From Obrien02, seems to improve plasticity by removing volume-change strain
#define REGULARIZE_CCM_THRESHOLD 1e-6f

namespace S2 {
namespace ms {
namespace fem {

const Real cEpsilon_dR_Numerical( 0.0001 ); //10e-3 works better in transition areas, 10e-2 works better in SVD-critical-point, which is not present in DAPD

//! ElementCache, constant during a subset of element lifetime
struct ElementCache
{
    // Degenerate
    int32 m_NoC; // local node-of-collapse, constant during degeneration
    Real m_DetF; // constant during CG iter, implicit in det(F) = det(Ds)/det(Dm)
    // CLFEM, PD_CCM (constant during CG iter)
    Mat2x2 m_R;
    // PD_CCM (constant during CG iter)
    //Mat2x2 m_S; //implicit in F = R*S
    // PK (constant during CG iter)
    Mat2x2 m_Ds;
    Mat2x2 m_F;
    //Mat2x2 m_U, m_DiagF, m_Vt; ??
    Mat2x2 m_InverseTransposeF; // used in implicit integration
#ifdef __ENABLE_PLASTICITY
    // Plasticity
    Vec3 m_Ep; //Plastic strain
#endif
    bool m_IsActive; //TEMP: Active constraint
};

class Cache: public IForce2::ICache
{
public:
    Cache( unsigned int num_elements, const Real* vec_inv_mass )
    : m_vecInvMass(vec_inv_mass), m_vecEC(0), m_vecPrevX(0), m_vecX(0), m_vecV(0), m_Time(0), m_TimeStep(0)
        {
            m_vecEC = new ElementCache[num_elements];
        }
    ~Cache()
        {
            if( m_vecEC ) delete[] m_vecEC;
            m_vecEC = 0;
        }
public:
    const Real *m_vecInvMass;
    ElementCache *m_vecEC;
    //\name Begin/EndEvaluation stuff
    //@{
    const Vec2 *m_vecPrevX;
    const Vec2 *m_vecX;
    const Vec2 *m_vecV;
    Real m_Time;
    Real m_TimeStep;
    //@}
#ifdef __S2_MS_ENABLE_STATS
public:
    struct Stats
    {
        uint32 m_Num_Degenerate;
        float32 m_Sum_Degenerate_DetF;
        Stats() : m_Num_Degenerate(0), m_Sum_Degenerate_DetF(0) {}
        void Begin() { m_Num_Degenerate = 0; m_Sum_Degenerate_DetF = 0; }
        void End() {}
    };
    Stats m_Stats;
#endif
};

int ComputeNoC( const Vec2 &a0, const Vec2 &a1, const Vec2 &a2,
                const Vec2 &b0, const Vec2 &b1, const Vec2 &b2,
                Real area, Real degenerate_threshold_det_F );

inline void LameParameters_From_YoungAndPoisson( Real young_modulus, Real poisson_ratio,
                                                 Real &lame_mu, Real &lame_lambda )
{
    lame_mu = Real(0.5) * young_modulus / (1+poisson_ratio);
    lame_lambda = young_modulus * poisson_ratio / ( (Real(1)+poisson_ratio) * (Real(1)-Real(2)*poisson_ratio) );
}

TriangleElementSet2::TriangleElementSet2()
: m_NumNodes(0)
, m_NumElements(0)
, m_TotalArea(0)
, m_vecRefPos(0)
, m_vecED(0)
, m_vecLED(0)
{}

TriangleElementSet2::~TriangleElementSet2()
{
    ClearBakedData();
}

void TriangleElementSet2::SetBakedData( bool b_shared,
                                              uint32 num_nodes, uint32 num_elements,
                                              const Vec2 *vec_ref_pos, const ElementData *vec_ed, const LinearED *vec_led )
{
    // Save baked data refs
    m_NumNodes = num_nodes;
    m_NumElements = num_elements;
    m_vecRefPos = vec_ref_pos;
    m_vecED = vec_ed;
    m_vecLED = vec_led;
}

void TriangleElementSet2::ClearBakedData()
{
    //\todo Following MeshSolidShape2 we should do this: if( m_pBuffer ) delete [] m_pBuffer;
    // Clear non-shared baked data
    //m_pBuffer = 0;
    m_vecRefPos = 0;
    m_vecED = 0;
    m_vecLED = 0;
}

//----------------------------------------------------------------
// Specific cache creation/Destruction
//----------------------------------------------------------------
IForce2::ICache* TriangleElementSet2::CreateCache( const Real* vec_inv_mass )
{
    // Alloc cached data
    Cache* p_cache = new Cache( m_NumElements, vec_inv_mass );
    // Init cached data
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache& ec( p_cache->m_vecEC[it_e] );
        ec.m_NoC = cInvalidNoC;
        ec.m_R = Mat2x2::Identity();
        ec.m_F = Mat2x2::Identity();
        ec.m_DetF = Real(1);
#ifdef __ENABLE_PLASTICITY
        ec.m_Ep = Vec3::Zero();
#endif
        //\todo Ds, InverseTransposeF...
        //\todo Consider full init, including F, R, etc...
        ec.m_IsActive = false;
    }
    return p_cache;
}

void TriangleElementSet2::DestroyCache( IForce2::ICache* p_cache )
{
    if( p_cache ) delete p_cache;
}

//----------------------------------------------------------------
// IForce virtual API
//----------------------------------------------------------------
void TriangleElementSet2::BeginEvaluation( IForce2::ICache* p_cache, const Vec2 *vec_prev_x, const Vec2 *vec_x, const Vec2 *vec_v, Real t ) const
{
    S2_MS_TRACE_BEGIN_SCOPE( "TriangleElementSet2::BeginEvaluation" );

    MS_ASSERT( 0 != p_cache );
    Cache &cache( *static_cast<Cache*>( p_cache ) );

#ifdef __S2_MS_ENABLE_STATS
    cache.m_Stats.Begin();
#endif
    // Save evaluation-constant params
    cache.m_vecPrevX = vec_prev_x;
    cache.m_vecX = vec_x;
    cache.m_vecV = vec_v;
    cache.m_TimeStep = t - cache.m_Time;
    cache.m_Time = t;
    // Compute common evaluation-constant EC stuff
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        ElementCache& ec( cache.m_vecEC[it_e] );
        Vec2 x0( vec_x[ ed.m_vecNID[0] ] );
        Vec2 x1( vec_x[ ed.m_vecNID[1] ] );
        Vec2 x2( vec_x[ ed.m_vecNID[2] ] );
        // Ds
        TriangleElement2::Compute_D( x0, x1, x2, ec.m_Ds );
        // F = Ds * Dm^-1
        ec.m_F = ec.m_Ds * ed.m_InvDm;
        ec.m_DetF = mal::Det( ec.m_F );
        /*
        MS_LOG_WARNING( "Element %d det(F) = %f", it_e, ec.m_DetF );
        //TEMP!
        {
            if( mal::Abs(ec.m_DetF) > 10 )
            {
                MS_LOG_WARNING( "Element %d |det(F)| > 10 (%f)", it_e, ec.m_DetF );
            }
        }
        */
        ec.m_InverseTransposeF = (mal::Abs(ec.m_DetF) > 0.000001f ) //\todo This is the same epsilon inside Inverse()... but we should find a better way than copying the constant...
                                 ? mal::Inverse( mal::Transposed( ec.m_F ) )
                                 : Mat2x2::Identity();
        // Check NaN
        MS_ASSERT( !mal::IsNaN( ec.m_F ) );
        MS_ASSERT( !mal::IsNaN( ec.m_DetF ) );
        MS_ASSERT( !mal::IsNaN( ec.m_InverseTransposeF ) );

        S2_MS_STAT_ADD( cache.m_Stats.m_Num_Degenerate, ((ec.m_DetF <= 0) ? 1 : 0) );
        S2_MS_STAT_ADD( cache.m_Stats.m_Sum_Degenerate_DetF, ((ec.m_DetF <= 0) ? ec.m_DetF : 0) );

#ifdef __S2_MS_ENABLE_TRACE //TEMP
        {
            char element_id[128];
            sprintf( element_id, "E[%u]", it_e );
            S2_MS_TRACE_GLOBAL( (std::string(element_id) + ".m_DetF").c_str(), ec.m_DetF );
        }
#endif //__S2_MS_ENABLE_TRACE
    }
    // Update inversion stuff: ec.m_NoC
    UpdateDegeneration( p_cache );

    // Update plastic flow, creep, whatever... \note Does nothing if m_TimeStep <= 0
    UpdatePlasticity( p_cache );

    // Update rotations
    switch( m_Params.m_RM )
    {
    case Params::eRM_Id: ComputeRotations_Id(p_cache); break;
    case Params::eRM_QR: ComputeRotations_QR(p_cache); break;
    case Params::eRM_MSVD: ComputeRotations_MSVD(p_cache); break;
        /*
    case Params::eRM_PDP: ComputeRotations_PD_Project(p_cache); break;
    case Params::eRM_PDR: ComputeRotations_PD_Reflect(p_cache); break;
        */
    case Params::eRM_DAPD: ComputeRotations_DAPD(p_cache); break;
    default: break;
    };

    S2_MS_TRACE_END_SCOPE();
}

// \note Forces are ACCUMULATED on vec_f, NOT SET into vec_f
void TriangleElementSet2::f( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    // Elastic + Friction + Plastic forces
    f_e( p_cache, vec_f );
    //f_d( vec_f ); \todo Needs vec_node_inv_mass for Rayleigh M term...
    f_p( p_cache, vec_f );
}

// \note Forces are ACCUMULATED on vec_f, NOT SET into vec_f
void TriangleElementSet2::f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    // Accumulate specific Elastic forces
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L: L_f_e( p_cache, vec_f ); break;
        // Corotational
    case Params::eMM_C_WRP: C_f_e( p_cache, vec_f ); break;
    case Params::eMM_C_LCM: C_LCM_f_e( p_cache, vec_f ); break;
    case Params::eMM_C_LCMH: C_LCM_f_e( p_cache, vec_f ); break;
    case Params::eMM_C_CCM: C_CCM_f_e( p_cache, vec_f ); break;
        // Hyperelastic
    case Params::eMM_H_LCM: H_LCM_f_e( p_cache, vec_f ); break;
    case Params::eMM_H_CCM: H_CCM_f_e( p_cache, vec_f ); break;
    case Params::eMM_H_NH_C0: H_NHC0_f_e( p_cache, vec_f ); break;
    case Params::eMM_H_NH_C1: H_NHC1_f_e( p_cache, vec_f ); break;
    default: break;
    };
}

// \note Forces are ACCUMULATED on vec_f, NOT SET into vec_f
void TriangleElementSet2::f_d( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    // Accumulate specific Elastic forces
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L: L_f_d( p_cache, vec_f ); break;
        // Corotational-Warped
    case Params::eMM_C_WRP: C_f_d( p_cache, vec_f ); break;
    case Params::eMM_C_LCM:
    case Params::eMM_C_LCMH:
    case Params::eMM_C_CCM:
        C_f_d( p_cache, vec_f ); //\todo HACK damping uses WRP, R is available, forces but may be inaccurate
        break;
        // Hyperelastic \todo HACK, we use C_f_d by now, because R IS available after H_XXX_f() evaluations... ugly but works
    case Params::eMM_H_LCM:
    case Params::eMM_H_CCM:
    case Params::eMM_H_NH_C0:
    case Params::eMM_H_NH_C1:
        C_f_d( p_cache, vec_f );
        break;
    default: break;
    };
}

// \note Forces are ACCUMULATED on vec_f, NOT SET into vec_f
void TriangleElementSet2::f_p( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
#ifdef __ENABLE_PLASTICITY
    // Accumulate specific Plastic forces
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L: L_f_p( p_cache, vec_f ); break;
        // Corotational
    case Params::eMM_C_WRP: C_f_p( p_cache, vec_f ); break;
        // Hyperelastic-Corotational \todo HACK, we use C_f_p by now, because R IS available after C_CCM_f() evaluations... ugly but works
    case Params::eMM_C_LCM:
    case Params::eMM_C_LCMH:
    case Params::eMM_C_CCM:
        C_f_p( p_cache, vec_f );
        break;
        // Hyperelastic \todo HACK, we use C_f_p by now, because R IS available after H_XXX_f() evaluations... ugly but works
    case Params::eMM_H_LCM:
    case Params::eMM_H_CCM:
    case Params::eMM_H_NH_C0:
    case Params::eMM_H_NH_C1:
        C_f_p( p_cache, vec_f );
        break;
    default: break;
    };
#endif
}

void TriangleElementSet2::df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const
{
    // Accumulate specific Elastic forces
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L: L_df_x( p_cache, vec_dx, vec_df ); break;
        // Corotational
    case Params::eMM_C_WRP: C_df_x( p_cache, vec_dx, vec_df ); break;
        // Hyperelastic-Corotational
    case Params::eMM_C_LCM: C_LCM_df_x( p_cache, vec_dx, vec_df ); break;
    case Params::eMM_C_LCMH: C_LCMH_df_x( p_cache, vec_dx, vec_df ); break;
    case Params::eMM_C_CCM: C_CCM_df_x( p_cache, vec_dx, vec_df ); break;
        // Hyperelastic
    case Params::eMM_H_LCM: H_LCM_df_x( p_cache, vec_dx, vec_df ); break;
    case Params::eMM_H_CCM: H_CCM_df_x( p_cache, vec_dx, vec_df ); break;
        // Hyperelastic \todo HACK, we use C_df_x by now
    case Params::eMM_H_NH_C0:
    case Params::eMM_H_NH_C1:
        C_df_x( p_cache, vec_dx, vec_df );
        break;
    default: break;
    };
}

void TriangleElementSet2::df_v( IForce2::ICache* p_cache, const Vec2 *vec_dv, Vec2 *vec_df ) const
{
}

Real TriangleElementSet2::V( IForce2::ICache* p_cache ) const
{
    // Accumulate specific Elastic forces
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L: return L_V( p_cache ); break;
        // Corotational
    case Params::eMM_C_WRP: return C_V( p_cache ); break;
        // Hyperelastic-Corotational
    case Params::eMM_C_LCM: return C_LCM_V( p_cache ); break;
    case Params::eMM_C_LCMH: return C_LCM_V( p_cache ); break;
    case Params::eMM_C_CCM: return C_CCM_V( p_cache ); break;
        // Hyperelastic
    case Params::eMM_H_LCM: return H_LCM_V( p_cache ); break;
    case Params::eMM_H_CCM: return H_CCM_V( p_cache ); break;
    case Params::eMM_H_NH_C0: return H_NHC0_V( p_cache ); break;
    case Params::eMM_H_NH_C1: return H_NHC1_V( p_cache ); break;
    default: return Real(0); break;
    };
}
void TriangleElementSet2::EndEvaluation( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache &cache( *static_cast<Cache*>( p_cache ) );

    // Unset evaluation-constant params
    cache.m_vecPrevX = 0;
    cache.m_vecX = 0;
    cache.m_vecV = 0;
#ifdef __S2_MS_ENABLE_STATS
    cache.m_Stats.End();
#endif
}

void TriangleElementSet2::ResetPlasticity( IForce2::ICache* p_cache )
{
    MS_ASSERT( 0 != p_cache );
    Cache &cache( *static_cast<Cache*>( p_cache ) );
#ifdef __ENABLE_PLASTICITY
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache &ec( cache.m_vecEC[it_e] );
        ec.m_Ep = Vec3::Zero();
    }
#endif
}

void TriangleElementSet2::GetNID( unsigned int eid,
                                        unsigned int &nid0, unsigned int &nid1, unsigned int &nid2 ) const
{
    nid0 = m_vecED[eid].m_vecNID[0];
    nid1 = m_vecED[eid].m_vecNID[1];
    nid2 = m_vecED[eid].m_vecNID[2];
}

Transform2 TriangleElementSet2::GetTransformE2W( const IForce2::ICache* p_cache, unsigned int eid, const Vec2* vec_pos ) const
{
    MS_ASSERT( 0 != p_cache );
    const Cache& cache( *static_cast<const Cache*>( p_cache ) );

    Transform2 Te2w;
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L:
        Te2w.m_Rot = Mat2x2::Identity();
        break;
        // Corotational
    case Params::eMM_C_WRP:
        Te2w.m_Rot = cache.m_vecEC[eid].m_R;
        break;
        // Hyperelastic-Corotational
    case Params::eMM_C_LCM:
    case Params::eMM_C_LCMH:
    case Params::eMM_C_CCM:
        Te2w.m_Rot = cache.m_vecEC[eid].m_R;
        break;
        // Hyperelastic: TEMP R=U*V^T is available, by now
    case Params::eMM_H_LCM:
    case Params::eMM_H_CCM:
    case Params::eMM_H_NH_C0:
    case Params::eMM_H_NH_C1:
        Te2w.m_Rot = cache.m_vecEC[eid].m_R;
        break;
    default:
        Te2w.m_Rot = Mat2x2::Identity();
        break;
    }
    const ElementData &ed( m_vecED[eid] );
    Te2w.m_Pos = mal::Rcp( Real(3) )*( vec_pos[ed.m_vecNID[0]] + vec_pos[ed.m_vecNID[1]] + vec_pos[ed.m_vecNID[2]] );
    return Te2w;
}

unsigned int TriangleElementSet2::GetNoC( const IForce2::ICache* p_cache, unsigned int eid ) const
{
    MS_ASSERT( 0 != p_cache );
    const Cache& cache( *static_cast<const Cache*>( p_cache ) );
    return cache.m_vecEC[eid].m_NoC;
}

void TriangleElementSet2::UpdateConstraints_ActiveSet( IForce2::ICache* p_cache, const Vec2* vec_vel, const Vec2* vec_acc )
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache& ec( cache.m_vecEC[it_e] );
        ec.m_IsActive = cInvalidNoC != ec.m_NoC;
        if( ec.m_IsActive )
        {
            const ElementData& ed( m_vecED[it_e] );
            //\todo We're computing the element normal TOO MANY TIMES, consider computing it wrt NoC ONCE and reusing it in DAPD and Constraints
            Vec2 normal(0,0);
            uint32 vid0( ed.m_vecNID[ec.m_NoC] );
            uint32 vid1( ed.m_vecNID[(ec.m_NoC + 1) % 3] );
            uint32 vid2( ed.m_vecNID[(ec.m_NoC + 2) % 3] );
            {
                Vec2 x12( cache.m_vecX[vid2] - cache.m_vecX[vid1] );
                Real dist12_sq( mal::NormSq(x12) );
                if( dist12_sq > 1e-6 )
                    normal = mal::PerpendicularCW( x12 ) / mal::Sqrt(dist12_sq);
            }
            if( mal::NormSq(normal) > 0 )
            {
                // Active if approaching relative vel or acc \todo Consider sampling vid1,vid2 vel/acc at their barycenter
                Real v0n_b12( mal::Dot(vec_vel[vid0] - Real(0.5)*(vec_vel[vid1]+vec_vel[vid2]),normal) );
                Real a0n_b12( mal::Dot(vec_acc[vid0] - Real(0.5)*(vec_acc[vid1]+vec_acc[vid2]),normal) );
                ec.m_IsActive = /*v0n_b12 < 0 || */ a0n_b12 < 0;
            }
        }
    }
}

void TriangleElementSet2::ApplyConstraints_NoC_Reaction( const IForce2::ICache* p_cache, Vec2* vec_v ) const
{
    MS_ASSERT( 0 != p_cache );
    const Cache& cache( *static_cast<const Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementCache& ec( cache.m_vecEC[it_e] );
        if( ec.m_IsActive ) //cInvalidNoC != ec.m_NoC ) //Degenerate and has proper NoC \todo ADD NORMAL FORCE activation CRITERIA AS IN CONTACTS
        {
            MS_ASSERT(cInvalidNoC != ec.m_NoC);
            const ElementData& ed( m_vecED[it_e] );
            /*\todo
            Vec2 normal = Compute_Normal_DAPD_Degenerate( ec.m_NoC,
                                                          cache.m_vecX[ ed.m_vecNID[0] ], cache.m_vecX[ ed.m_vecNID[1] ], cache.m_vecX[ ed.m_vecNID[2] ] );
            */
            Vec2 normal(0,0);
            uint32 vid0( ed.m_vecNID[ec.m_NoC] );
            uint32 vid1( ed.m_vecNID[(ec.m_NoC + 1) % 3] );
            uint32 vid2( ed.m_vecNID[(ec.m_NoC + 2) % 3] );
            {
                Vec2 x12( cache.m_vecX[vid2] - cache.m_vecX[vid1] );
                Real dist12_sq( mal::NormSq(x12) );
                if( dist12_sq > 1e-6 )
                    normal = mal::PerpendicularCW( x12 ) / mal::Sqrt(dist12_sq);
            }
            if( mal::NormSq(normal) > 0 )
            {
                // Apply constraint to all nodes (ad-hoc weights may be wrong...)
                Real v0n_b12( mal::Dot(vec_v[vid0] - Real(0.5)*(vec_v[vid1]+vec_v[vid2]),normal) );
                //\todo Alternate activation crashes CG, we MUST use a constant active set per CG solve as with contacts
                // if( v0n_b12 < 0 )
                {
                    Real cWeight(1);
                    vec_v[vid0] -= cWeight*v0n_b12*normal;
                    vec_v[vid1] += 0.5f*cWeight*v0n_b12*normal;
                    vec_v[vid2] += 0.5f*cWeight*v0n_b12*normal;
                }
            }
        }
    }
}

//----------------------------------------------------------------
// Internal methods
//----------------------------------------------------------------

/*!  IMPORTANT: The original version of UpdateDegeneration() version
   relied blindly on the cached value of ec.m_NoC to determine if the
   element *was* degenerate or not at m_vecPrevX. This works as long
   as m_vecPrevX is EXACTLY the state that was passed to
   BeginEvaluation() as vec_x at the *previous* timestep and used by
   UpdateDegeneration() to test for degeneration.

     However, this condition may not be fulfilled for several reasons:
   Reinitialization, Teleport, Newton-Raphson iterative vec_x
   approximations with the same vec_prev_x, etc...

     Therefore, we added the check (*), which is rather expensive, but
   covers all cases. If this becomes an issue, a faster specialized
   UpdateDegeneration_RelyOnCacheNoC() version could be added and used
   where safe.

   \todo This still fails if: prev and current positions are the SAME
   and the element is degenerate. This happens during mouse-drag, not
   sure why, but results in discont. If it happens, the element is
   degenerate but NoC is invalid, and we keep trying to compute the
   crossing despite the previous pos being degenerate, which does not
   make sense. __NEW_CROSSING_DETECTION should fix it, try if for a
   while, specially with FIE, and accept it if no problems appear.
*/
void TriangleElementSet2::UpdateDegeneration( IForce2::ICache* p_cache ) const
{
#define __NEW_CROSSING_DETECTION //\todo should be default
#ifdef __NEW_CROSSING_DETECTION
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Check all elements for degeneration according to det(F)
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache &ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //remains/becomes undegenerate
            ec.m_NoC = cInvalidNoC;
        else //remains/becomes degenerate
        {
            /*\note This code does NOT assume that ec.m_NoC is up-to
              date, and therefore re-checks if prevX was undegenerate,
              computing NoC if so. This is *required* for
              Newton-Raphson methods, where several
              Begin/EndEvaluation() are performed for configuration
              pairs (prevX,X_k), where prevX is fixed but X_k
              changes.
            */
            const ElementData &ed( m_vecED[it_e] );
            Vec2 a0( cache.m_vecPrevX[ ed.m_vecNID[0] ] );
            Vec2 a1( cache.m_vecPrevX[ ed.m_vecNID[1] ] );
            Vec2 a2( cache.m_vecPrevX[ ed.m_vecNID[2] ] );
            Mat2x2 prevDs;
            TriangleElement2::Compute_D( a0, a1, a2, prevDs );
            if( mal::Det(prevDs*ed.m_InvDm) > m_Params.m_DegenerateThresholdDetF ) // was undegenerate (*)
            {
                //compute NoC
                Vec2 b0( cache.m_vecX[ ed.m_vecNID[0] ] );
                Vec2 b1( cache.m_vecX[ ed.m_vecNID[1] ] );
                Vec2 b2( cache.m_vecX[ ed.m_vecNID[2] ] );
                ec.m_NoC = ComputeNoC( a0, a1, a2,
                                       b0, b1, b2,
                                       ed.m_Area,
                                       m_Params.m_DegenerateThresholdDetF );
                if( cInvalidNoC == ec.m_NoC )
                    MS_LOG_ERROR( "New Degenerate element %u CANNOT SOLVE Eq2", it_e );
            }
            else if( cInvalidNoC == ec.m_NoC )
            {
                /*\todo The element *was* degenerate and remains so,
                  but has no valid NoC, therefore, no continuous-time
                  crossing can be found, we should handle it with
                  discrete-time heuristic.

                  \todo Consider flagging this elements with NoC =
                  cDiscreteNoC and use such flag to select different
                  method.
                */
            }
            //else: remains degenerate with the same NoC
        }
    }

#ifdef __DISABLED_CODE //\todo THIS version would be cheaper, but may fail if (prevX,X) configurations are not sequential, because it assumes ec.m_NoC to be always correct. Fails for NR.
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Check all elements for degeneration according to det(F)
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache &ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //remains/becomes undegenerate
            ec.m_NoC = cInvalidNoC;
        else if( cInvalidNoC == ec.m_NoC ) //becomes degenerate or remains degenerate but discrete-time
        {
            const ElementData &ed( m_vecED[it_e] );
            Vec2 a0( cache.m_vecPrevX[ ed.m_vecNID[0] ] );
            Vec2 a1( cache.m_vecPrevX[ ed.m_vecNID[1] ] );
            Vec2 a2( cache.m_vecPrevX[ ed.m_vecNID[2] ] );
            Mat2x2 prevDs;
            TriangleElement2::Compute_D( a0, a1, a2, prevDs );
            if( mal::Det(prevDs*ed.m_InvDm) > m_Params.m_DegenerateThresholdDetF ) // was undegenerate (*)
            {
                //compute NoC
                Vec2 b0( cache.m_vecX[ ed.m_vecNID[0] ] );
                Vec2 b1( cache.m_vecX[ ed.m_vecNID[1] ] );
                Vec2 b2( cache.m_vecX[ ed.m_vecNID[2] ] );
                ec.m_NoC = ComputeNoC( a0, a1, a2,
                                       b0, b1, b2,
                                       ed.m_Area,
                                       m_Params.m_DegenerateThresholdDetF );
                if( cInvalidNoC == ec.m_NoC )
                    MS_LOG_ERROR( "New Degenerate element %u CANNOT SOLVE Eq2", it_e );
            }
            else // remains degenerate but discrete
            {
                /*\todo The element *was* degenerate and remains so,
                  but has no valid NoC, therefore, no continuous-time
                  crossing can be found, we should handle it with
                  discrete-time heuristic.

                  \todo Consider flagging this elements with NoC =
                  cDiscreteNoC and use such flag to select different
                  method.
                */
            }
        }
        //else: remains degenerate with the same NoC \note ASSUMES
        //ec.m_NoC is up to date and NOT persistent from an unrelated
        //BeginEvaluation(), which is WRONG FOR Newton-Raphson

    }
#  endif

#else //TEMPORAL: Old code, remove it when new sufficiently tested
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Check all elements for degeneration according to det(F)
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache &ec( cache.m_vecEC[it_e] );
        const ElementData &ed( m_vecED[it_e] );
        Vec2 a0( cache.m_vecPrevX[ ed.m_vecNID[0] ] );
        Vec2 a1( cache.m_vecPrevX[ ed.m_vecNID[1] ] );
        Vec2 a2( cache.m_vecPrevX[ ed.m_vecNID[2] ] );
        Mat2x2 prevDs;
        TriangleElement2::Compute_D( a0, a1, a2, prevDs );
        if( cInvalidNoC == ec.m_NoC // No cached degeneration
            ||
            mal::Det(prevDs*ed.m_InvDm) > m_Params.m_DegenerateThresholdDetF ) // was undegenerate (*)
        {
            if( ec.m_DetF <= m_Params.m_DegenerateThresholdDetF ) //becomes degenerate
            {
                //compute NoC
                Vec2 b0( cache.m_vecX[ ed.m_vecNID[0] ] );
                Vec2 b1( cache.m_vecX[ ed.m_vecNID[1] ] );
                Vec2 b2( cache.m_vecX[ ed.m_vecNID[2] ] );
                MS_ASSERT( mal::Det(prevDs*ed.m_InvDm) > m_Params.m_DegenerateThresholdDetF && ec.m_DetF <= m_Params.m_DegenerateThresholdDetF ); //TEMPORAL: this is what ComputeNoC expects...
                ec.m_NoC = ComputeNoC( a0, a1, a2,
                                       b0, b1, b2,
                                       ed.m_Area,
                                       m_Params.m_DegenerateThresholdDetF );
                if( cInvalidNoC == ec.m_NoC )
                    MS_LOG_ERROR( "New Degenerate element %u CANNOT SOLVE Eq2", it_e );
            }
            //else remains undegenerate
        }
        else //was degenerate AND has a valid NoC
        {
            if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //becomes undegenerate
            {
                ec.m_NoC = cInvalidNoC;
            }
            //else remains degenerate with same NoC
        }
    }
#endif
}

/*! Compute Node-of-Collapse
  \note Internally computes ToC
  \todo Should compute Direction-of-Collapse (F1,F2) in 3D
*/
int ComputeNoC( const Vec2 &a0, const Vec2 &a1, const Vec2 &a2,
                const Vec2 &b0, const Vec2 &b1, const Vec2 &b2,
                Real area, Real degenerate_threshold_det_F )
{
    int noc( TriangleElementSet2::cInvalidNoC );
    // Gather node positions and displacements
    Vec2 d0( b0 - a0 );
    Vec2 d1( b1 - a1 );
    Vec2 d2( b2 - a2 );
    // Compute ToC, when det(F) = m_DegenerateThresholdDetF => det(Q) = m_DegenerateThresholdDetF * det(P) => det(Q) = m_DegenerateThresholdDetF * 2 * Area(P)
    Real K( degenerate_threshold_det_F * 2 * area );
    Real A(   d1.x() * d2.y()
              + d2.x() * d0.y()
              + d0.x() * d1.y()
              - d1.x() * d0.y()
              - d0.x() * d2.y()
              - d2.x() * d1.y() );
    Real B( (   d1.x() * a2.y() + a1.x() * d2.y() )
            + ( d2.x() * a0.y() + a2.x() * d0.y() )
            + ( d0.x() * a1.y() + a0.x() * d1.y() )
            - ( d1.x() * a0.y() + a1.x() * d0.y() )
            - ( d0.x() * a2.y() + a0.x() * d2.y() )
            - ( d2.x() * a1.y() + a2.x() * d1.y() ) );
    Real C(   a1.x() * a2.y()
              + a2.x() * a0.y()
              + a0.x() * a1.y()
              - a1.x() * a0.y()
              - a0.x() * a2.y()
              - a2.x() * a1.y()
              - K );
    Real toc1(0), toc2(1);
    int num_roots( mal::GSolvePolynomialEq2<Real>( A, B, C, toc1, toc2 ) );
    if( num_roots > 0 )
    {
        //MS_ASSERT( num_roots == 1 );
        //MS_LOG_WARNING( "Potentially degenerated %d with toc = %f, %f", it_e, toc1, toc2 );
        MS_ASSERT( toc1 <= toc2 );
        //\todo We knot that det(F(0)) > DTDF, and that det(F(1)) <= DTDF, so there MUST BE STRICTLY 1 crossing in toc = [0,1]
        Real toc( toc1 > 0 ? toc1 : toc2 );
        if( toc > 1 )
        {
            MS_LOG_ERROR( "toc %f > 1", toc );
        }
        /* analyze inversion interval:
           0 <= toc1 <= toc2 < 1 => Collapse->Invert->Uninvert->Uncollapse
           toc1 < 0 <= toc2 < 1 => Collapse->Invert
           toc1 <= toc2 < 0 => Extrapolated collapse outside [0,1], no actual collapse
           1 <= toc1 <= toc2 < 1 => Extrapolated collapse outside [0,1], no actual collapse
        */
        /*TEMP
          if( !(toc2 < 0 || toc1 > 1) //Inside [0,1]
          &&
          (toc1 < 0 || toc2 > 1) //2 roots inside [0,1] //\todo MAY UNINVERT ACROSS AXIS DIFERENT than collapsed one!!
          )
        */
        {
            // Compute DoC
            Vec2 q0( a0 + toc*d0 );
            Vec2 q1( a1 + toc*d1 );
            Vec2 q2( a2 + toc*d2 );
            Real sql01( mal::NormSq( q1-q0 ) );
            Real sql12( mal::NormSq( q2-q1 ) );
            Real sql20( mal::NormSq( q0-q2 ) );
            // Choose NoC, vertex opposite to longest edge at ToC
            noc = (sql01 >= sql12)
                  ? (sql01 >= sql20) ? 2 : 1
                  : (sql12 >= sql20) ? 0 : 1;
            /*TEMP: Debug
              std::cout << "From "
              << a0 << ","
              << a1 << ","
              << a2 << " to "
              << m_vecX[ ed.m_vecNID[0] ] << ","
              << m_vecX[ ed.m_vecNID[1] ] << ","
              << m_vecX[ ed.m_vecNID[2] ]
              << " ToC " << toc1 << "," << toc2
              << " NoC " << noc
              << " PoC " << ( (sql01 >= sql12) ? (sql01 >= sql20) ? q2 : q1 : (sql12 >= sql20) ? q0 : q1 )
              << " Det(Q(ToC)) " << mal::Det( mal::GMat2x2_From_Columns(q1-q0,q2-q0) ) << " ?= D = " << D << std::endl;
            */
            //MS_LOG_WARNING( "Actually degenerate %d with noc = %d, toc = %f", it_e, noc, toc );
            //toc = toc1;
        }
        /*TEMP
          else
          {
          MS_LOG_WARNING( "No 0 or 2 degenerations in range [0,1]" );
          }
        */
    }
    else
    {
        MS_LOG_ERROR( "Cannot solve A*t^2 + B*t + C = 0 for A = %e, B = %e, C = %e", A, B, C );
    }
    return noc;
}

/* Update plastic strain from current total and elastic strain.
   \todo Formulas are wrong due to strain matrix frobenius norm != strain vector norm
   \todo If we end up saving per-element strain or computing it from total_strain G = 0.5(F + F^T) - I, it may be easier to perform all computations with G instead of the strain vector
*/
void TriangleElementSet2::UpdatePlasticity( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );

#ifdef __ENABLE_PLASTICITY
    if( m_Params.m_PlasticYield <= 0
        || m_Params.m_PlasticMax <= 0
        || m_Params.m_PlasticCreepPerSecond <= 0
        || cache.m_TimeStep <= 0 ) return;
    Real plastic_creep_factor( mal::Min( 1.0f, m_Params.m_PlasticCreepPerSecond * cache.m_TimeStep ) );

    // Update plastic strain for all elements
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L:
        {
            for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
            {
                ElementCache &ec( cache.m_vecEC[it_e] );
                const ElementData& ed( m_vecED[it_e] );
                const LinearED& led( m_vecLED[it_e] );
#ifdef __COMPUTE_STRAIN_FROM_B
                node_index_type nid0( ed.m_vecNID[0] );
                node_index_type nid1( ed.m_vecNID[1] );
                node_index_type nid2( ed.m_vecNID[2] );
                Vec2 u0( cache.m_vecX[nid0] - m_vecRefPos[nid0] );
                Vec2 u1( cache.m_vecX[nid1] - m_vecRefPos[nid1] );
                Vec2 u2( cache.m_vecX[nid2] - m_vecRefPos[nid2] );
                // Compute total strain = B_e * u_e
                Mat3x6 B;
                TriangleElement2::Compute_B( m_vecRefPos[nid0], m_vecRefPos[nid1], m_vecRefPos[nid2], B );
                Vec6 u( u0[0], u0[1], u1[0], u1[1], u2[0], u2[1] );
                Vec3 Et = B * u;
#else //From F
                // Compute total strain == 1/2(F+F^T) - I
                Mat2x2 strain_matrix = Real(0.5) * ( (ec.m_F + mal::Transposed(ec.m_F)) ) - Mat2x2::Identity();
                Vec3 Et = Vec3( strain_matrix(0,0), strain_matrix(1,1), 2*strain_matrix(0,1) );
#endif
                // Compute elastic strain
                Vec3 Ee = Et - ec.m_Ep;
#ifdef __ENABLE_PLASTIC_STRAIN_DEVIATION
                // Remove volume-change (=> Ee' = Ee - I * Trace(Ee)/2 but with Ee in vector form)
                Real Ee_vol = Real(0.5)*(Ee[0]+Ee[1]);
                Ee[0] -= Ee_vol;
                Ee[1] -= Ee_vol;
#endif

                // Update plastic strain
                Real sq_norm2_Ee = mal::Sq(Ee[0]) + mal::Sq(Ee[1]) + Real(0.5)*mal::Sq(Ee[2]); // Frobenius norm of the 2x2 Ee tensor //Last term is 2*Sq(Ee(1,1))  = 2*Sq( 0.5*Ee[2] ) != mal::NormSq( Ee ); //\todo THIS should be the Frobenius norm of the 2x2 strain tensor [ g_xx, g_xy, g_yx, g_yy ], which is NOT the same as tehe norm of the vector [ g_xx, g_yy, 2*g_xy ] (there's a 2 factor in g_xy term that gets squared and shouldn't)
                // If yield threshold exceeded...
                if( sq_norm2_Ee > mal::Sq( m_Params.m_PlasticYield ) )
                {
                    // Update Ep according to creep
                    Real norm2_Ee = mal::Sqrt( sq_norm2_Ee );
                    ec.m_Ep += ( plastic_creep_factor * ( norm2_Ee - m_Params.m_PlasticYield ) / norm2_Ee ) * Ee;
                    // Clamp Ep according to max yield
                    Real sq_norm2_Ep = mal::Sq( ec.m_Ep[0] ) + mal::Sq( ec.m_Ep[1] ) + Real(0.5)*mal::Sq( ec.m_Ep[2] ); //Frobenius norm of the 2x2 Ep tensor //mal::NormSq( ec.m_Ep );
                    if( sq_norm2_Ep > mal::Sq( m_Params.m_PlasticMax ) )
                        ec.m_Ep *= m_Params.m_PlasticMax / mal::Sqrt( sq_norm2_Ep );
                    //MS_LOG_WARNING("Ep = %f", mal::Norm( ec.m_Ep ) );
                }
            }
        }
        break;
        // Corotational
    case Params::eMM_C_WRP:
        {
            for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
            {
                ElementCache &ec( cache.m_vecEC[it_e] );
                const ElementData& ed( m_vecED[it_e] );
                const LinearED& led( m_vecLED[it_e] );
#ifdef __COMPUTE_STRAIN_FROM_B
                node_index_type nid0( ed.m_vecNID[0] );
                node_index_type nid1( ed.m_vecNID[1] );
                node_index_type nid2( ed.m_vecNID[2] );
                //\todo u_i computation is the ONLY difference with
                //Linear version, consider merging both someway,
                //probably computing per-element total_strain in a
                //separate (specific) pass and using it to update
                //plasticity in a common pass.
                Mat2x2 Rt( ec.m_R.Transposed() );
                Vec2 u0( Rt * cache.m_vecX[nid0] - m_vecRefPos[nid0] );
                Vec2 u1( Rt * cache.m_vecX[nid1] - m_vecRefPos[nid1] );
                Vec2 u2( Rt * cache.m_vecX[nid2] - m_vecRefPos[nid2] );
                // Compute total strain = B_e * u_e
                Mat3x6 B; //\todo Could be precomputed
                TriangleElement2::Compute_B( m_vecRefPos[nid0], m_vecRefPos[nid1], m_vecRefPos[nid2], B );
                Vec6 u( u0[0], u0[1], u1[0], u1[1], u2[0], u2[1] );
                Vec3 Et = B * u;
#else //From F => S = R^T * F
                // Compute total strain == 1/2(S+S^T) - I
                Mat2x2 S = mal::Transposed(ec.m_R) * ec.m_F;
                Mat2x2 strain_matrix = Real(0.5) * ( (S + mal::Transposed(S)) ) - Mat2x2::Identity();
                Vec3 Et = Vec3( strain_matrix(0,0), strain_matrix(1,1), 2*strain_matrix(0,1) );
                //MS_LOG_WARNING( "Diff = %f", mal::NormSq( Et - Et2 ) ); => OK
#endif

                // Compute elastic strain
                Vec3 Ee = Et - ec.m_Ep;
#ifdef __ENABLE_PLASTIC_STRAIN_DEVIATION
                // Remove volume-change (=> Ee' = Ee - I * Trace(Ee)/2 but with Ee in vector form)
                Real Ee_vol = Real(0.5)*(Ee[0]+Ee[1]);
                Ee[0] -= Ee_vol;
                Ee[1] -= Ee_vol;
#endif
                // Update plastic strain
                Real sq_norm2_Ee = mal::Sq(Ee[0]) + mal::Sq(Ee[1]) + Real(0.5)*mal::Sq(Ee[2]); // Frobenius norm of the 2x2 Ee tensor //Last term is 2*Sq(Ee(1,1))  = 2*Sq( 0.5*Ee[2] ) != mal::NormSq( Ee ); //\todo THIS should be the Frobenius norm of the 2x2 strain tensor [ g_xx, g_xy, g_yx, g_yy ], which is NOT the same as tehe norm of the vector [ g_xx, g_yy, 2*g_xy ] (there's a 2 factor in g_xy term that gets squared and shouldn't)
                // If yield threshold exceeded...
                if( sq_norm2_Ee > mal::Sq( m_Params.m_PlasticYield ) )
                {
                    // Update Ep according to creep
                    Real norm2_Ee = mal::Sqrt( sq_norm2_Ee );
                    ec.m_Ep += ( plastic_creep_factor * ( norm2_Ee - m_Params.m_PlasticYield ) / norm2_Ee ) * Ee;
                    // Clamp Ep according to max yield
                    Real sq_norm2_Ep = mal::Sq( ec.m_Ep[0] ) + mal::Sq( ec.m_Ep[1] ) + Real(0.5)*mal::Sq( ec.m_Ep[2] ); //Frobenius norm of the 2x2 Ep tensor //mal::NormSq( ec.m_Ep );
                    if( sq_norm2_Ep > mal::Sq( m_Params.m_PlasticMax ) )
                        ec.m_Ep *= m_Params.m_PlasticMax / mal::Sqrt( sq_norm2_Ep );
                    //MS_LOG_WARNING("Ep = %f", mal::Norm( ec.m_Ep ) );
                }
            }
        }
        break;
        // Hyperelastic-Corotational
    case Params::eMM_C_LCM:
    case Params::eMM_C_LCMH:
    case Params::eMM_C_CCM:
        //\todo
        break;
        // Hyperelastic
    case Params::eMM_H_LCM: break;
    case Params::eMM_H_CCM: break;
    case Params::eMM_H_NH_C0: break;
    case Params::eMM_H_NH_C1: break;
        //\todo
    default: break;
    };
#endif
}

//----------------------------------------------------------------

void TriangleElementSet2::ComputeRotations_Id( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache& ec( cache.m_vecEC[it_e] );
        ec.m_R = Mat2x2::Identity();
    }
}

void TriangleElementSet2::ComputeRotations_QR( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache& ec( cache.m_vecEC[it_e] );
        ec.m_R = Compute_R_QR( ec.m_F );
        MS_ASSERT( !mal::IsNaN( ec.m_R ) );
    }
}

void TriangleElementSet2::ComputeRotations_MSVD( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache& ec( cache.m_vecEC[it_e] );
        ec.m_R = Compute_R_SVD( ec.m_F, ec.m_DetF, m_Params.m_DegenerateThresholdDetF );
        MS_ASSERT( !mal::IsNaN( ec.m_R ) );
    }
}

void TriangleElementSet2::ComputeRotations_PD_Project( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache& ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //Undegenerate
        {
            ec.m_R = Compute_R_PD( ec.m_F, ec.m_DetF );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else if( cInvalidNoC != ec.m_NoC ) //Degenerate and has proper NoC
        {
            const ElementData& ed( m_vecED[it_e] );
            ec.m_R = Compute_R_PDP_Degenerate( ec.m_F, ec.m_DetF, m_Params.m_DegenerateThresholdDetF,
                                               ec.m_NoC,
                                               cache.m_vecX[ ed.m_vecNID[0] ], cache.m_vecX[ ed.m_vecNID[1] ], cache.m_vecX[ ed.m_vecNID[2] ],
                                               ed.m_Area, ed.m_InvDm );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else //Degenerate but no proper NoC
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element %u has invalid NoC. Using R=Id", it_e );
            ec.m_R = Mat2x2::Identity();
        }
    }
}

void TriangleElementSet2::ComputeRotations_PD_Reflect( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache& ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //Undegenerate
        {
            ec.m_R = Compute_R_PD( ec.m_F, ec.m_DetF );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else if( cInvalidNoC != ec.m_NoC ) //Degenerate and has proper NoC
        {
            const ElementData& ed( m_vecED[it_e] );
            ec.m_R = Compute_R_PDR_Degenerate( ec.m_F, ec.m_DetF, m_Params.m_DegenerateThresholdDetF,
                                               ec.m_NoC,
                                               cache.m_vecX[ ed.m_vecNID[0] ], cache.m_vecX[ ed.m_vecNID[1] ], cache.m_vecX[ ed.m_vecNID[2] ] );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else //Degenerate but no proper NoC
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element %u has invalid NoC. Using R=Id", it_e );
            ec.m_R = Mat2x2::Identity();
        }
    }
}

void TriangleElementSet2::ComputeRotations_DAPD( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache& ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //Undegenerate
        {
            ec.m_R = Compute_R_PD( ec.m_F, ec.m_DetF );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else if( cInvalidNoC != ec.m_NoC ) //Degenerate and has proper NoC
        {
            const ElementData& ed( m_vecED[it_e] );
            ec.m_R = Compute_R_DAPD_Degenerate( ec.m_F, ec.m_DetF, m_Params.m_DegenerateThresholdDetF,
                                                ec.m_NoC,
                                                cache.m_vecX[ ed.m_vecNID[0] ], cache.m_vecX[ ed.m_vecNID[1] ], cache.m_vecX[ ed.m_vecNID[2] ],
                                                ed.m_Area, ed.m_InvDm,
                                                m_Params.m_DAPD_L_Factor, m_Params.m_DAPD_NL_Factor, m_Params.m_DAPD_NL_Exponent );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else //Degenerate but no proper NoC
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element %u has invalid NoC. Using R=Id", it_e );
            ec.m_R = Mat2x2::Identity();
        }
    }
}

static bool Compute_dR_From_RM( const Params& params,
                                const Cache& cache,
                                const TriangleElementSet2::ElementData& ed, const ElementCache& ec,
                                const Mat2x2& Rt, const Vec2* vec_dx,
                                Mat2x2& dR )
{
    bool bCorrect( false );
    if( params.m_DM != Params::eDM_Truncated )
    {
        TriangleElementSet2::node_index_type nid0( ed.m_vecNID[0] );
        TriangleElementSet2::node_index_type nid1( ed.m_vecNID[1] );
        TriangleElementSet2::node_index_type nid2( ed.m_vecNID[2] );
        if( params.m_DM == Params::eDM_Exact
            || params.m_DM == Params::eDM_Exact_Inv )
        {
            if( params.m_RM == Params::eRM_DAPD )
            {
                if( params.m_DM == Params::eDM_Exact_Inv ) //We use Inv as it should work regardless of element degeneration state
                    bCorrect = Compute_dR_DAPD( params.m_DegenerateThresholdDetF, params.m_DAPD_L_Factor, params.m_DAPD_NL_Factor, params.m_DAPD_NL_Exponent,
                                                ec.m_NoC,
                                                cache.m_vecX[nid0], cache.m_vecX[nid1], cache.m_vecX[nid2],
                                                vec_dx[nid0], vec_dx[nid1], vec_dx[nid2],
                                                ed.m_InvDm, ec.m_F, ec.m_R, Rt, ed.m_Area,
                                                dR );
                else if( ec.m_DetF > params.m_DegenerateThresholdDetF )
                    bCorrect = Compute_dR( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2],
                                           ed.m_InvDm, ec.m_F, ec.m_R, Rt,
                                           dR );
            }
            else if( params.m_RM == Params::eRM_MSVD ) //McAdams dR is only correct for MSVD inversion handling
                bCorrect = Compute_dR( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2],
                                       ed.m_InvDm, ec.m_F, ec.m_R, Rt,
                                       dR );
            else if( params.m_RM == Params::eRM_QR )
                bCorrect = Compute_dR_QR_YX( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2],
                                             ed.m_InvDm, ec.m_F,
                                             dR );
        }
        else if( params.m_DM == Params::eDM_Numerical )
        {
            if( params.m_RM == Params::eRM_DAPD ) //Accepts Inv
                bCorrect =  Compute_dR_DAPD_Numerical( params.m_DegenerateThresholdDetF, params.m_DAPD_L_Factor, params.m_DAPD_NL_Factor, params.m_DAPD_NL_Exponent,
                                                       ec.m_NoC,
                                                       cache.m_vecX[nid0], cache.m_vecX[nid1], cache.m_vecX[nid2],
                                                       vec_dx[nid0], vec_dx[nid1], vec_dx[nid2],
                                                       ed.m_InvDm, ed.m_Area,
                                                       cEpsilon_dR_Numerical,
                                                       dR );
            else if( params.m_RM == Params::eRM_MSVD
                     || ec.m_DetF > params.m_DegenerateThresholdDetF )
                bCorrect = Compute_dR_Numerical( cache.m_vecX[nid0], cache.m_vecX[nid1], cache.m_vecX[nid2],
                                                 vec_dx[nid0], vec_dx[nid1], vec_dx[nid2],
                                                 ed.m_InvDm, params.m_DegenerateThresholdDetF,
                                                 cEpsilon_dR_Numerical,
                                                 dR );
        }
    }
    return bCorrect;
}

//----------------------------------------------------------------
// Linear
void TriangleElementSet2::L_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Vec2 u0( cache.m_vecX[nid0] - m_vecRefPos[nid0] ); //\todo Consider precomputing vec_u ONCE at BeginEvaluation() to avoid repeating x_i-r_i
        Vec2 u1( cache.m_vecX[nid1] - m_vecRefPos[nid1] );
        Vec2 u2( cache.m_vecX[nid2] - m_vecRefPos[nid2] );
        vec_f[ nid0 ] -= led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2;
        vec_f[ nid1 ] -= led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2;
        vec_f[ nid2 ] -= led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2;
    }
}

/*! Linear Rayleigh Damping:
  C_e = m_RayleighCoeff[0] * M_e + m_RayleighCoeff[1] * K_e;
  fd_e = -C_e * V_e
*/
void TriangleElementSet2::L_f_d( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Vec2 v0( cache.m_vecV[nid0] );
        Vec2 v1( cache.m_vecV[nid1] );
        Vec2 v2( cache.m_vecV[nid2] );
        vec_f[ nid0 ] -= ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid0]) ) * v0
                         + m_Params.m_RayleighCoeff[1] * (led.m_K00 * v0 + led.m_K01 * v1 + led.m_K02 * v2);
        vec_f[ nid1 ] -= ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid1]) ) * v1
                         + m_Params.m_RayleighCoeff[1] * (led.m_K10 * v0 + led.m_K11 * v1 + led.m_K12 * v2);
        vec_f[ nid2 ] -= ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid2]) ) * v2
                         + m_Params.m_RayleighCoeff[1] * (led.m_K20 * v0 + led.m_K21 * v1 + led.m_K22 * v2);
    }
}

/*! Linear Plasticity:
  P_e = precomputed
  fp_e = +P_e * plastic_strain;
  TEMP: In Muller04 scheme, f_p SUBSTRACTS the excess of elastic force
  added when using strain_total in f_e computation, instead of
  strain_elastic = strain_total - strain_plastic. For this reason, the
  sign of plastic forces (+=) is the opposite of elastic ones (-=)
*/
void TriangleElementSet2::L_f_p( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
#ifdef __ENABLE_PLASTICITY
    if( m_Params.m_PlasticYield <= 0
        || m_Params.m_PlasticMax <= 0
        || m_Params.m_PlasticCreepPerSecond <= 0 ) return;
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        vec_f[ nid0 ] += led.m_P0 * ec.m_Ep;
        vec_f[ nid1 ] += led.m_P1 * ec.m_Ep;
        vec_f[ nid2 ] += led.m_P2 * ec.m_Ep;
    }
#endif
}

//----------------------------------------------------------------
// Corotational
void TriangleElementSet2::C_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        // Rotate world frame X to reference frame (no need to translate)
        Mat2x2 Rt( ec.m_R.Transposed() );
        Vec2 u0( Rt * cache.m_vecX[nid0] - m_vecRefPos[nid0] ); //\todo NO POINT IN precomputing vec_u ONCE at BeginEvaluation(), because depends on per-element Rt
        Vec2 u1( Rt * cache.m_vecX[nid1] - m_vecRefPos[nid1] );
        Vec2 u2( Rt * cache.m_vecX[nid2] - m_vecRefPos[nid2] );
        // Compute elastic forces and rotate to world frame
        vec_f[ nid0 ] -= ec.m_R * ( led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 );
        vec_f[ nid1 ] -= ec.m_R * ( led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 );
        vec_f[ nid2 ] -= ec.m_R * ( led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 );

        //TEMP
        MS_ASSERT( !mal::IsNaN( vec_f[ nid0 ] ) );
        MS_ASSERT( !mal::IsNaN( vec_f[ nid1 ] ) );
        MS_ASSERT( !mal::IsNaN( vec_f[ nid2 ] ) );
    }
}

/*! Corotational Rayleigh Damping:
  C_e = m_RayleighCoeff[0] * M_e + m_RayleighCoeff[1] * K_e;
  fd_e = - R_e * ( C_e * R_e^T * V_e )
*/
void TriangleElementSet2::C_f_d( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Mat2x2 Rt( ec.m_R.Transposed() );
        Vec2 v0( Rt * cache.m_vecV[nid0] );
        Vec2 v1( Rt * cache.m_vecV[nid1] );
        Vec2 v2( Rt * cache.m_vecV[nid2] );
        vec_f[ nid0 ] -= ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid0]) ) * v0
                                    + m_Params.m_RayleighCoeff[1] * (led.m_K00 * v0 + led.m_K01 * v1 + led.m_K02 * v2) );
        vec_f[ nid1 ] -= ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid1]) ) * v1
                                    + m_Params.m_RayleighCoeff[1] * (led.m_K10 * v0 + led.m_K11 * v1 + led.m_K12 * v2) );
        vec_f[ nid2 ] -= ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid2]) ) * v2
                                    + m_Params.m_RayleighCoeff[1] * (led.m_K20 * v0 + led.m_K21 * v1 + led.m_K22 * v2) );

        //TEMP
        MS_ASSERT( !mal::IsNaN( vec_f[ nid0 ] ) );
        MS_ASSERT( !mal::IsNaN( vec_f[ nid1 ] ) );
        MS_ASSERT( !mal::IsNaN( vec_f[ nid2 ] ) );
    }
}

/*! Corotational Plasticity:
  P_e = precomputed
  fp_e = +P_e * plastic_strain;
  TEMP: In Muller04 scheme, f_p SUBSTRACTS the excess of elastic force
  added when using strain_total in f_e computation, instead of
  strain_elastic = strain_total - strain_plastic. For this reason, the
  sign of plastic forces (+=) is the opposite of elastic ones (-=)
*/
void TriangleElementSet2::C_f_p( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
#ifdef __ENABLE_PLASTICITY
    if( m_Params.m_PlasticYield <= 0
        || m_Params.m_PlasticMax <= 0
        || m_Params.m_PlasticCreepPerSecond <= 0 ) return;
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        vec_f[ nid0 ] += ec.m_R * (led.m_P0 * ec.m_Ep);
        vec_f[ nid1 ] += ec.m_R * (led.m_P1 * ec.m_Ep);
        vec_f[ nid2 ] += ec.m_R * (led.m_P2 * ec.m_Ep);
        //TEMP
        MS_ASSERT( !mal::IsNaN( vec_f[ nid0 ] ) );
        MS_ASSERT( !mal::IsNaN( vec_f[ nid1 ] ) );
        MS_ASSERT( !mal::IsNaN( vec_f[ nid2 ] ) );
    }
#endif
}

void TriangleElementSet2::C_LCM_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute P
        Mat2x2 S( ec.m_R.Transposed() * ec.m_F ); //F = R*S
        Mat2x2 S_minus_Id( S-Mat2x2::Identity() );
        Mat2x2 P( ec.m_R * ( 2 * m_LameMu * S_minus_Id
                             + m_LameLambda * mal::Trace(S_minus_Id) * Mat2x2::Identity() ) );

        // Compute PK force
        Mat2x2 H( - ed.m_Area * P * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
}

void TriangleElementSet2::C_CCM_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute \bar S from \bar R and F
        Mat2x2 S( ec.m_R.Transposed() * ec.m_F ); //F = R*S
        Real det_S( mal::Det(S) ); //\note == det(F) \todo Ensure and subst det_S by ec.m_DetF if correct

        // Compute P
        Mat2x2 P( ec.m_R * ( 2*m_LameMu*(S-Mat2x2::Identity())
                             + m_LameLambda * (det_S-1) * Mat2x2( S(1,1), -S(1,0),
                                                                  -S(0,1), S(0,0) ) ) );

        // Compute PK force
        Mat2x2 H( - ed.m_Area * P * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
}

//----------------------------------------------------------------
// P-K
inline void ComputeSVD_Invertible( const Mat2x2 &F, Mat2x2 &U, Vec2 &vec_diag_F, Mat2x2 &Vt )
{
    mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
    // Force pure rotation U,Vt by fixing potential inversion/reflection
    if( mal::Det(U) < 0 ) { mal::GSetColumn<1>( U, -mal::GColumn<1>(U) ); vec_diag_F[1] = -vec_diag_F[1]; }
    if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); vec_diag_F[1] = -vec_diag_F[1]; }
}

//---- ComputeDiagonalP_XXX
inline Vec2 ComputeDiagonalP_LCM( const Vec2 &vec_diag_F, Real lame_mu, Real lame_lambda )
{
    /* TEMP: Literal implementation from ITF, slower
       Mat2x2r diag_F( vec_diag_F[0], 0,
       0, vec_diag_F[1] );
       Mat2x2r E( diag_F - Mat2x2r::Identity() );
       Mat2x2r diag_P( 2*m_LameMu * E
       + m_LameLambda * mal::Trace( E ) * Mat2x2r::Identity() );
    */
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Real trE( vec_diag_F[0] + vec_diag_F[1] - 2 );
    return Vec2( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * trE,
                 2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * trE );
}

inline Vec2 ComputeDiagonalP_CCM( const Vec2 &vec_diag_F, Real lame_mu, Real lame_lambda )
{
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Real J( det_F );
    return Vec2( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * (J-1) * vec_diag_F[1],
                 2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * (J-1) * vec_diag_F[0] );
}

inline Vec2 ComputeDiagonalP_NHC0( const Vec2 &vec_diag_F, Real lame_mu, Real lame_lambda, Real ecie_e )
{
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Vec2 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                             mal::Max( ecie_e, vec_diag_F[1] ) );
    Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
    return Vec2( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ),
                 lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
}

inline Vec2 ComputeDiagonalP_NHC1( const Vec2 &vec_diag_F, Real lame_mu, Real lame_lambda, Real ecie_e, Real ecie_k )
{
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    bool bIsDegenerate0( vec_diag_F[0] < ecie_e );
    bool bIsDegenerate1( vec_diag_F[1] < ecie_e );
    Vec2 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                             mal::Max( ecie_e, vec_diag_F[1] ) );
    Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
    Vec2 vec_diag_P( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ),
                     lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
    Real h01( lame_lambda / clamped_J );
    //\todo Cases can be generalized if delta0,delta1 are 0 for non-clamped axis
    if( bIsDegenerate0 && bIsDegenerate1 )
    {
        Real delta0( vec_diag_F[0] - ecie_e );
        Real delta1( vec_diag_F[1] - ecie_e );
        vec_diag_P[0] += h01*delta1 + 2*ecie_k*delta0;
        vec_diag_P[1] += h01*delta0 + 2*ecie_k*delta1;
    }
    else if( bIsDegenerate0 && !bIsDegenerate1 )
    {
        Real delta0( vec_diag_F[0] - ecie_e );
        vec_diag_P[0] += 2*ecie_k*delta0;
        vec_diag_P[1] += h01*delta0;
    }
    else if( !bIsDegenerate0 && bIsDegenerate1 )
    {
        Real delta1( vec_diag_F[1] - ecie_e );
        vec_diag_P[0] += h01*delta1;
        vec_diag_P[1] += 2*ecie_k*delta1;
    }
    return vec_diag_P;
}

// P-K Forces
void TriangleElementSet2::H_LCM_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        ec.m_R = U*Vt; //TEMP: for debug viz AND for some hacked methods that use co-rotational code path instead...

        // Compute diag_P
        Vec2 vec_diag_P( ComputeDiagonalP_LCM( vec_diag_F, m_LameMu, m_LameLambda ) );
        Mat2x2 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );

        // Compute PK force
        Mat2x2 H( - ed.m_Area * U * diag_P * Vt * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
}

void TriangleElementSet2::H_CCM_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        ec.m_R = U*Vt; //TEMP: for debug viz AND for some hacked methods that use co-rotational code path instead...

        // Compute diag_P
        Vec2 vec_diag_P( ComputeDiagonalP_CCM( vec_diag_F, m_LameMu, m_LameLambda ) );
        Mat2x2 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );

        // Compute PK force
        Mat2x2 H( - ed.m_Area * U * diag_P * Vt * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
}

void TriangleElementSet2::H_NHC0_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        ec.m_R = U*Vt; //TEMP: for debug viz AND for some hacked methods that use co-rotational code path instead...

        // Compute diag_P
        Vec2 vec_diag_P( ComputeDiagonalP_NHC0( vec_diag_F, m_LameMu, m_LameLambda, m_Params.m_ECIE_e_threshold ) );
        Mat2x2 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );

        // Compute PK force
        Mat2x2 H( - ed.m_Area * U * diag_P * Vt * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
}

void TriangleElementSet2::H_NHC1_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        ec.m_R = U*Vt; //TEMP: for debug viz AND for some hacked methods that use co-rotational code path instead...

        // Compute diag_P
        Vec2 vec_diag_P( ComputeDiagonalP_NHC1( vec_diag_F,
                                                m_LameMu, m_LameLambda,
                                                m_Params.m_ECIE_e_threshold, m_Params.m_ECIE_k_factor * m_Params.m_YoungModulus ) );
        Mat2x2 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );

        // Compute PK force
        Mat2x2 H( - ed.m_Area * U * diag_P * Vt * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
}

Real TriangleElementSet2::L_V( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic energy
    Real acc_twice_V(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Vec2 u0( cache.m_vecX[nid0] - m_vecRefPos[nid0] );
        Vec2 u1( cache.m_vecX[nid1] - m_vecRefPos[nid1] );
        Vec2 u2( cache.m_vecX[nid2] - m_vecRefPos[nid2] );
        // Ve = 1/2 * Ue^T * Ke * Ue
        acc_twice_V += mal::Dot( u0, led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 );
        acc_twice_V += mal::Dot( u1, led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 );
        acc_twice_V += mal::Dot( u2, led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 );
    }
    return Real(0.5) * acc_twice_V;
}
Real TriangleElementSet2::C_V( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic energy
    Real acc_twice_V(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        // Rotate world frame X to reference frame (no need to translate)
        Mat2x2 Rt( ec.m_R.Transposed() );
        Vec2 u0( Rt * cache.m_vecX[nid0] - m_vecRefPos[nid0] );
        Vec2 u1( Rt * cache.m_vecX[nid1] - m_vecRefPos[nid1] );
        Vec2 u2( Rt * cache.m_vecX[nid2] - m_vecRefPos[nid2] );
        // Ve = 1/2 * Ue^T * Ke * Ue, reference-system invariant
        acc_twice_V += mal::Dot( u0, led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 );
        acc_twice_V += mal::Dot( u1, led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 );
        acc_twice_V += mal::Dot( u2, led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 );
    }
    return Real(0.5) * acc_twice_V;
}
Real TriangleElementSet2::C_LCM_V( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic energy
    Real acc_V(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        // Compute S from R
        Mat2x2 S( ec.m_R.Transposed() * ec.m_F ); //F = R*S
        Mat2x2 S_minus_Id( S-Mat2x2::Identity() );
        acc_V += ed.m_Area * ( m_LameMu * ( S_minus_Id ).NormSqF()  //|S-I|^2_Frobenius
                               + Real(0.5)*m_LameLambda*mal::Sq( mal::Trace(S_minus_Id) ) );
    }
    return acc_V;
}
Real TriangleElementSet2::C_CCM_V( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic energy
    Real acc_V(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute S from R
        Mat2x2 S( ec.m_R.Transposed() * ec.m_F ); //F = R*S
        Real det_S( mal::Det(S) ); //\note == det(F) \todo Ensure and subst det_S by ec.m_DetF if correct

        /*\todo Elastic energy is tricky here... first term is
          taken from the Siggraph2012 course notes on LinearFEM
          but using "S", the second one is taken from CCM "J"
          term using det_S instead. I think it's correct, but the
          resulting energy is not conserved, and CCM one is.

          The results seem to be quite depending on the timestep,
          which means it's highly nonlinear... for dt=0.0001
          trajectories are much different than dt=0.001
        */
        acc_V += ed.m_Area * ( m_LameMu * ( S-Mat2x2::Identity() ).NormSqF()  //|S-I|^2_Frobenius
                               + Real(0.5)*m_LameLambda*mal::Sq( det_S - Real(1) ) );
    }
    return acc_V;
}
Real TriangleElementSet2::H_LCM_V( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic energy
    Real acc_V(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        // \todo May be faster to use F instead of computing SVD again
        acc_V += ed.m_Area * ( m_LameMu * ( mal::Sq(vec_diag_F[0]-Real(1)) + mal::Sq(vec_diag_F[1]-Real(1)) )
                               + Real(0.5)*m_LameLambda*mal::Sq( (vec_diag_F[0]-Real(1)) + (vec_diag_F[1]-Real(1)) ) );
    }
    return acc_V;
}
Real TriangleElementSet2::H_CCM_V( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic energy
    Real acc_V(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        // \todo May be faster to use F instead of computing SVD again
        acc_V += ed.m_Area * ( m_LameMu * ( mal::Sq(vec_diag_F[0]-Real(1)) + mal::Sq(vec_diag_F[1]-Real(1)) )
                               + Real(0.5)*m_LameLambda*mal::Sq( ec.m_DetF - Real(1) ) );
    }
    return acc_V;
}
Real TriangleElementSet2::H_NHC0_V( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic energy
    Real acc_V(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        //
        Real ecie_e( m_Params.m_ECIE_e_threshold );
        Vec2 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                                 mal::Max( ecie_e, vec_diag_F[1] ) );
        Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
        // Clamped energy
        Real acc_density_Ve(0);
        acc_density_Ve += Real(0.5)*m_LameMu * ( mal::Sq(vec_clamped_diag_F[0]) + mal::Sq(vec_clamped_diag_F[1]) - Real(2) )
                          - m_LameMu * mal::Log(clamped_J)
                          + Real(0.5)*m_LameLambda*mal::Sq( mal::Log(clamped_J) );
        // Extrapolated energy
        /*TEMP: Trying to extract proper energy from NHC0 we see
          that it's NOT consistent... we cannot extrapolate it, if
          forces are clamped, so is energy? and the result is that
          initial TotalEnergy is not preserved...
        */
        Real g0( m_LameMu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( m_LameLambda * mal::Log(clamped_J) - m_LameMu ) );
        Real g1( m_LameMu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( m_LameLambda * mal::Log(clamped_J) - m_LameMu ) );
        Real h01( m_LameLambda / clamped_J );
        Real delta0( mal::Min<Real>( vec_diag_F[0] - ecie_e, 0 ) );
        Real delta1( mal::Min<Real>( vec_diag_F[1] - ecie_e, 0 ) );
        acc_density_Ve += g0 * delta0
                          + g1 * delta1
                          + h01 * delta0 * delta1;
        acc_V += ed.m_Area * acc_density_Ve;
    }
    return acc_V;
}
Real TriangleElementSet2::H_NHC1_V( IForce2::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic energy
    Real acc_V(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        //
        Real ecie_e( m_Params.m_ECIE_e_threshold );
        Real ecie_k( m_Params.m_ECIE_k_factor * m_Params.m_YoungModulus );
        bool bIsDegenerate0( vec_diag_F[0] < ecie_e );
        bool bIsDegenerate1( vec_diag_F[1] < ecie_e );
        Vec2 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                                 mal::Max( ecie_e, vec_diag_F[1] ) );
        Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
        // Clamped energy
        Real acc_density_Ve(0);
        acc_density_Ve += Real(0.5)*m_LameMu * ( mal::Sq(vec_clamped_diag_F[0]) + mal::Sq(vec_clamped_diag_F[1]) - Real(2) )
                          - m_LameMu * mal::Log(clamped_J)
                          + Real(0.5)*m_LameLambda*mal::Sq( mal::Log(clamped_J) );
        // Extrapolated energy
        Real g0( m_LameMu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( m_LameLambda * mal::Log(clamped_J) - m_LameMu ) );
        Real g1( m_LameMu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( m_LameLambda * mal::Log(clamped_J) - m_LameMu ) );
        Real h01( m_LameLambda / clamped_J );
        if( bIsDegenerate0 && bIsDegenerate1 )
        {
            Real delta0( vec_diag_F[0] - ecie_e );
            Real delta1( vec_diag_F[1] - ecie_e );
            acc_density_Ve += g0 * delta0
                              + g1 * delta1
                              + h01 * delta0 * delta1
                              + ecie_k * ( mal::Sq(delta0) + mal::Sq(delta1) );
        }
        else if( bIsDegenerate0 && !bIsDegenerate1 )
        {
            Real delta0( vec_diag_F[0] - ecie_e );
            acc_density_Ve += g0 * delta0 + ecie_k*mal::Sq(delta0);
        }
        else if( !bIsDegenerate0 && bIsDegenerate1 )
        {
            Real delta1( vec_diag_F[1] - ecie_e );
            acc_density_Ve += g1 * delta1 + ecie_k*mal::Sq(delta1);
        }
        acc_V += ed.m_Area * acc_density_Ve;
    }
    return acc_V;
}

//----------------------------------------------------------------
// df_x()
//----------------------------------------------------------------
void TriangleElementSet2::L_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Vec2 dx0( vec_dx[nid0] );
        Vec2 dx1( vec_dx[nid1] );
        Vec2 dx2( vec_dx[nid2] );
        vec_df[ nid0 ] -= led.m_K00 * dx0 + led.m_K01 * dx1 + led.m_K02 * dx2;
        vec_df[ nid1 ] -= led.m_K10 * dx0 + led.m_K11 * dx1 + led.m_K12 * dx2;
        vec_df[ nid2 ] -= led.m_K20 * dx0 + led.m_K21 * dx1 + led.m_K22 * dx2;
    }
}

void TriangleElementSet2::C_df_x( IForce2::ICache* p_cache, const Vec2* vec_dx, Vec2* vec_df ) const
{
    S2_MS_TRACE_BEGIN_SCOPE( "TriangleElementSet2::C_df_x" );

    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );

    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        // Rotate world frame X to reference frame (no need to translate)
        Mat2x2 Rt( ec.m_R.Transposed() );
        Vec2 ldx0( Rt * vec_dx[nid0] );
        Vec2 ldx1( Rt * vec_dx[nid1] );
        Vec2 ldx2( Rt * vec_dx[nid2] );
        //\todo THIS CODE IS SO LONG THAT IT MAY BE FASTER to use the Piola-Kirchhoff version from McAdams2011 instead of precomputed K_ij
        // Add truncated dR main term: df = R_e * K_e * R_e^T * dx
        vec_df[ nid0 ] -= ec.m_R * ( led.m_K00 * ldx0 + led.m_K01 * ldx1 + led.m_K02 * ldx2 );
        vec_df[ nid1 ] -= ec.m_R * ( led.m_K10 * ldx0 + led.m_K11 * ldx1 + led.m_K12 * ldx2 );
        vec_df[ nid2 ] -= ec.m_R * ( led.m_K20 * ldx0 + led.m_K21 * ldx1 + led.m_K22 * ldx2 );

        // Compute df terms that depend on dR
        Mat2x2 dR;
        if( Compute_dR_From_RM(m_Params,cache,ed,ec,Rt,vec_dx,dR) )
        {
            // Add exact dR force terms
            // Term dR_e * K_e * (R_e^T x_e - r_e) = dR_e * K_e * u_e
            Vec2 u0( Rt * cache.m_vecX[nid0] - m_vecRefPos[nid0] );
            Vec2 u1( Rt * cache.m_vecX[nid1] - m_vecRefPos[nid1] );
            Vec2 u2( Rt * cache.m_vecX[nid2] - m_vecRefPos[nid2] );
            vec_df[ nid0 ] -= dR * ( led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 );
            vec_df[ nid1 ] -= dR * ( led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 );
            vec_df[ nid2 ] -= dR * ( led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 );

            // Term R_e * K_e * dR_e^T * x_e
            Mat2x2 dRt( dR.Transposed() );
            Vec2 dRt_times_x0( dRt * cache.m_vecX[nid0] );
            Vec2 dRt_times_x1( dRt * cache.m_vecX[nid1] );
            Vec2 dRt_times_x2( dRt * cache.m_vecX[nid2] );
            vec_df[ nid0 ] -= ec.m_R * ( led.m_K00 * dRt_times_x0 + led.m_K01 * dRt_times_x1 + led.m_K02 * dRt_times_x2 );
            vec_df[ nid1 ] -= ec.m_R * ( led.m_K10 * dRt_times_x0 + led.m_K11 * dRt_times_x1 + led.m_K12 * dRt_times_x2 );
            vec_df[ nid2 ] -= ec.m_R * ( led.m_K20 * dRt_times_x0 + led.m_K21 * dRt_times_x1 + led.m_K22 * dRt_times_x2 );
        }

        MS_ASSERT( !mal::IsNaN( vec_df[ nid0 ] ) );
        MS_ASSERT( !mal::IsNaN( vec_df[ nid1 ] ) );
        MS_ASSERT( !mal::IsNaN( vec_df[ nid2 ] ) );
    }

    S2_MS_TRACE_END_SCOPE();
}

/* SVD-free LCM implementation from McAdams2011 & tech doc, using R,S and dR from DAPD
   \todo SIMPLER and probably FASTER than C_df_x DAPD with dR terms...
*/
void TriangleElementSet2::C_LCM_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Mat2x2 Rt( mal::Transposed( ec.m_R ) );
        Mat2x2 S( Rt * ec.m_F );
        Mat2x2 dDs;
        TriangleElement2::Compute_D( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], dDs );
        Mat2x2 dF( dDs * ed.m_InvDm );
        // Compute dP terms independent of dR, from McAdams2011
        Mat2x2 dP( 2*m_LameMu * dF
                   + m_LameLambda*mal::Trace(Rt*dF) * ec.m_R );
        // Compute dP terms that depend on dR, from McAdams2011
        Mat2x2 dR;
        if( Compute_dR_From_RM(m_Params,cache,ed,ec,Rt,vec_dx,dR) )
            dP += ( m_LameLambda*mal::Trace(S-Mat2x2::Identity()) - 2*m_LameMu ) * dR;

        // Compute PK force
        Mat2x2 dH( - ed.m_Area * dP * ed.m_InvDm.Transposed() );
        vec_df[ nid0 ] += - mal::GColumn<0>(dH) - mal::GColumn<1>(dH); //dh0 = -dh1-dh2;
        vec_df[ nid1 ] += mal::GColumn<0>(dH);
        vec_df[ nid2 ] += mal::GColumn<1>(dH);
    }
}

/* Hybrid LCM from DCNLFEM */
void TriangleElementSet2::C_LCMH_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Mat2x2 Rt( mal::Transposed( ec.m_R ) );
        Mat2x2 S( Rt * ec.m_F );
        Mat2x2 dDs;
        TriangleElement2::Compute_D( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], dDs );
        Mat2x2 dF( dDs * ed.m_InvDm );

        //TEMPORAL Using param m_ECIE_k_factor, consider removing it or defining specific param here
        Mat2x2 dP;
        if( m_Params.m_DM == Params::eDM_Truncated
            &&
            m_Params.m_ECIE_k_factor > 0.1 ) //Hybrid LCM
        {
            const Real cDimension = 2;
            const Real cRadiusWRP = 1; //\todo Just use 1... m_rParams.m_ECIE_k_factor;
            Real weight_LCM = mal::Clamp01( mal::Rcp(cRadiusWRP) * mal::Trace( mal::Abs( S - Mat2x2::Identity() ) )/cDimension );
            dP = 2*m_LameMu * (weight_LCM*dF + (1-weight_LCM)*ec.m_R*mal::GMatNxN_Symmetric_Part(Rt*dF) )
                 + m_LameLambda*mal::Trace(Rt*dF) * ec.m_R;
        }
        else if( m_Params.m_DM == Params::eDM_Truncated
                 &&
                 m_Params.m_ECIE_k_factor > 10 ) //WRP
        {
            dP = 2*m_LameMu * ec.m_R*mal::GMatNxN_Symmetric_Part(Rt*dF) + m_LameLambda*mal::Trace(Rt*dF) * ec.m_R;
        }
        else //LCM
        {
            dP = 2*m_LameMu * dF + m_LameLambda*mal::Trace(Rt*dF) * ec.m_R;
        }

        // Compute dP terms that depend on dR, from McAdams2011
        Mat2x2 dR;
        if( Compute_dR_From_RM(m_Params,cache,ed,ec,Rt,vec_dx,dR) )
            dP += ( m_LameLambda*mal::Trace(S-Mat2x2::Identity()) - 2*m_LameMu ) * dR;

        // Compute PK force
        Mat2x2 dH( - ed.m_Area * dP * ed.m_InvDm.Transposed() );
        vec_df[ nid0 ] += - mal::GColumn<0>(dH) - mal::GColumn<1>(dH); //dh0 = -dh1-dh2;
        vec_df[ nid1 ] += mal::GColumn<0>(dH);
        vec_df[ nid2 ] += mal::GColumn<1>(dH);
    }
}

/* SVD-free CCM implementation \note Requires regularization of 1/0 terms */
void TriangleElementSet2::C_CCM_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache& cache( *static_cast<Cache*>( p_cache ) );
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Mat2x2 Rt( mal::Transposed( ec.m_R ) ); //\note This R comes from SVD1 U*Vt
        Mat2x2 dDs;
        TriangleElement2::Compute_D( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], dDs );
        Mat2x2 dF( dDs * ed.m_InvDm );
        // Compute dP terms independent of dR, from McAdams2011
        Mat2x2 dP( 2*m_LameMu * dF ); //Robust dP term
        // Compute invF-dependent dP terms \todo There MUST BE a better regularization, see testFE, SPLIT seems the most reasonable approach
        //\todo See adj() in http://en.wikipedia.org/wiki/Adjugate_matrix#Cayley.E2.80.93Hamilton_formula ... I think there MUST BE a non-singular expression for dP
        Mat2x2 invF( mal::Inverse( ec.m_F, ec.m_DetF + mal::Sign(ec.m_DetF)*REGULARIZE_CCM_THRESHOLD ) ); //\note Regularized for det(F) = 0
        Mat2x2 invFt( mal::Transposed( invF ) );
        dP += m_LameLambda * ( (2*mal::Sq(ec.m_DetF) - ec.m_DetF) * mal::Trace(invF*dF) * invFt
                               - (mal::Sq(ec.m_DetF) - ec.m_DetF) * invFt * mal::Transposed(dF) * invFt );
        // Compute dP terms that depend on dR
        Mat2x2 dR;
        if( Compute_dR_From_RM(m_Params,cache,ed,ec,Rt,vec_dx,dR) )
            dP += - 2*m_LameMu * dR;

        // Compute PK force
        Mat2x2 dH( - ed.m_Area * dP * ed.m_InvDm.Transposed() );
        vec_df[ nid0 ] += - mal::GColumn<0>(dH) - mal::GColumn<1>(dH); //dh0 = -dh1-dh2;
        vec_df[ nid1 ] += mal::GColumn<0>(dH);
        vec_df[ nid2 ] += mal::GColumn<1>(dH);
    }
}


/* \todo Use SVD differentials */
void TriangleElementSet2::H_LCM_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const
{
    //\todo By now use corotational version, SVD df is more complex
    return C_LCM_df_x( p_cache, vec_dx, vec_df );

#ifdef __STILL_DEVELOPING //\todo This would be the generic implementation based in SVD, not in dR, from ECIE and SiggCourse12
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Mat2x2 dDs;
        TriangleElement2::Compute_D( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], dDs );
        Mat2x2 dF( dDs * ed.m_InvDm );

        //--- Compute dP
        //-- (a) Directly, using dP(F;dF) matrix formula \todo BUT THIS DOES NOT FIX DEGENERATE F!!
        //\todo For comparison...

        //-- (b) From diagonalization, using dP(...) = U * delta_diag_P( d_diag_F ) * V^t
        //\sa "Robust Quasistatic Finite Elements and Flesh Simulation" and "Energetically Consistent Invertible Elasticity"
        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat2x2 U, Vt;
        Vec2 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        /* Compute diag_P
        Vec2 vec_diag_P( ComputeDiagonalP_LCM( vec_diag_F, m_LameMu, m_LameLambda ) );
        Mat2x2 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );
        */
        Mat2x2 d_diag_F( mal::Transposed(U) * dF * mal::Transposed(Vt) ); // U^T * \delta F * V

        /*\todo If I get it right, here we need to compute
          \frac{\partial P}{\partial \hat F}, which is a 2x2x2x2
          tensor (converted into a 4x4 matrix), and multiply it by
          \delta \hat F 2x2 matrix (converted into a 4x1 vector) to
          get a new 2x2 matrix (from a 4x1 vector). The 4x4 matrix
          must be "projected" to ensure positive-definiteness using
          eigenvalue clamping, according to both RQFEFS and ECIE.
          Moreover, all uses of singular values must refer to the
          corrected ones (potentially negative), and the entries of
          the 4x4 matrix CAN be problematic for nearly repeated
          ones. ECIE supplementary doc explains how to compute terms
          robustly, but is hell...

          All in all, this part seems REALLY SLOW (in 3D it involves a
          9x9 matrix and 9x1 vectors, a 3x3 eigenanalysis... in
          addition to the 3D SVD...)

          The non-diagonalized P(F) version, found in Sigcourse2012
          does NOT seem to be in danger of not being positive
          definite, even though it should be equivalent for
          non-inverted configurations.

          If this is as SLOW as it seems to be, the PD-P-CCM version
          could be considerably faster, but maybe less accurate. In a
          non-linear solver, it may be worth running more inexact but
          cheap iterations instead of fewer but more expensive exact
          iterations.
        */
        Mat2x2 d_diag_P( ComputeDeltaDiagonalP_LCM( vec_diag_F, vec_delta_diag_F, m_LameMu, m_LameLambda ) );
        Mat2x2 dP( U * d_diag_P * Vt );

        // Compute PK force
        Mat2x2 dH( - ed.m_Area * dP * ed.m_InvDm.Transposed() );
        vec_df[ nid0 ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //dh0 = -dh1-dh2;
        vec_df[ nid1 ] += mal::GColumn<0>(H);
        vec_df[ nid2 ] += mal::GColumn<1>(H);
    }
#endif
}

/* \todo Use SVD differentials */
void TriangleElementSet2::H_CCM_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const
{
    //\todo By now use corotational version, SVD df is more complex
    return C_CCM_df_x( p_cache, vec_dx, vec_df );
}

//----------------------------------------------------------------
// EditableTriangleElementSet2
//----------------------------------------------------------------

EditableTriangleElementSet2::EditableTriangleElementSet2()
: m_LastElement(0)
, m_pAllocRefPos(0)
, m_pAllocED(0)
, m_pAllocLED(0)
{}

EditableTriangleElementSet2::~EditableTriangleElementSet2()
{
    ClearEditData();
}

void EditableTriangleElementSet2::ClearEditData()
{
    if( m_pAllocRefPos ) delete[] m_pAllocRefPos;
    if( m_pAllocED ) delete[] m_pAllocED;
    if( m_pAllocLED ) delete[] m_pAllocLED;
    m_NumNodes = 0;
    m_LastElement = 0;
    m_NumElements = 0;
    m_pAllocRefPos = 0;
    m_pAllocED = 0;
    m_pAllocLED = 0;
}

//\todo if method is given per-element, params.m_Method == eMethod_PerElement and alloc a per-element method descriptor
void EditableTriangleElementSet2::BeginEdition( unsigned int num_nodes, const Vec2 *vec_r, unsigned int num_elements, const Params &params )
{
    ClearEditData();
    m_NumNodes = num_nodes;
    m_NumElements = num_elements;
    //\todo here we could copy existing baked data as edit data to edit incrementally
    // Set params
    m_Params = params;
    // Alloc common stuff
    m_pAllocRefPos = new Vec2[num_nodes];
    memcpy( &m_pAllocRefPos[0], &vec_r[0], num_nodes * sizeof(Vec2) );
    m_pAllocED = new ElementData[num_elements];
    // Alloc specific stuff
    m_pAllocLED = new LinearED[num_elements];
    // Clear baked data
    ClearBakedData();
}

void EditableTriangleElementSet2::AddElement( unsigned int i, unsigned int j, unsigned int k )
{
    MS_ASSERT( node_index_type(i) == i && node_index_type(j) == j && node_index_type(k) == k );
    MS_ASSERT( m_LastElement < m_NumElements );
    ElementData& ed( m_pAllocED[m_LastElement++] );
    ed.m_vecNID[0] = node_index_type(i);
    ed.m_vecNID[1] = node_index_type(j);
    ed.m_vecNID[2] = node_index_type(k);
}

bool EditableTriangleElementSet2::EndEdition()
{
    MS_ASSERT( m_LastElement == m_NumElements );
    LameParameters_From_YoungAndPoisson( m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                         m_LameMu, m_LameLambda );
    for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
    {
        ElementData& ed( m_pAllocED[it_e] );
        Vec2 r0( m_pAllocRefPos[ ed.m_vecNID[0] ] );
        Vec2 r1( m_pAllocRefPos[ ed.m_vecNID[1] ] );
        Vec2 r2( m_pAllocRefPos[ ed.m_vecNID[2] ] );
        // Common element init
        ed.m_Area = TriangleElement2::Compute_Area( r0, r1, r2 );
        m_TotalArea += ed.m_Area;
        // Specific element init
        // Linear
        Mat6x6 K;
        TriangleElement2::Compute_K( r0, r1, r2,
                                     m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                     K );
        // K0X
        m_pAllocLED[it_e].m_K00 = mal::GRange<0,0,1,1>( K );
        m_pAllocLED[it_e].m_K01 = mal::GRange<0,2,1,3>( K );
        m_pAllocLED[it_e].m_K02 = mal::GRange<0,4,1,5>( K );
        // K1X
        m_pAllocLED[it_e].m_K10 = mal::GRange<2,0,3,1>( K );
        m_pAllocLED[it_e].m_K11 = mal::GRange<2,2,3,3>( K );
        m_pAllocLED[it_e].m_K12 = mal::GRange<2,4,3,5>( K );
        // K2X
        m_pAllocLED[it_e].m_K20 = mal::GRange<4,0,5,1>( K );
        m_pAllocLED[it_e].m_K21 = mal::GRange<4,2,5,3>( K );
        m_pAllocLED[it_e].m_K22 = mal::GRange<4,4,5,5>( K );

#ifdef __ENABLE_PLASTICITY
        Mat6x3 P;
        TriangleElement2::Compute_P( r0, r1, r2,
                                     m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                     //m_Params.m_PlasticYield, m_Params.m_PlasticMax, m_Params.m_PlasticCreep,
                                     P );
        m_pAllocLED[it_e].m_P0 = mal::GRange<0,0,1,2>( P ); //0,0 => (1,2)
        m_pAllocLED[it_e].m_P1 = mal::GRange<2,0,3,2>( P ); //2,0 => (3,2)
        m_pAllocLED[it_e].m_P2 = mal::GRange<4,0,5,2>( P ); //4,0 => (5,2)
#endif

        // C-H
        Mat2x2 Dm;
        TriangleElement2::Compute_D( r0, r1, r2, Dm );
        ed.m_InvDm = mal::Inverse( Dm );
    }
    // Init baked data
    TriangleElementSet2::SetBakedData( false, m_NumNodes, m_NumElements, m_pAllocRefPos, m_pAllocED, m_pAllocLED );
    return true;
}

void EditableTriangleElementSet2::SetParams( const Params &params )
{
    MS_ASSERT( IsValid() );
    m_Params = params;
    LameParameters_From_YoungAndPoisson( m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                         m_LameMu, m_LameLambda );
    // Rebuild material-dependent stuff directly
    for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
    {
        ElementData& ed( m_pAllocED[it_e] );
        Vec2 r0( m_pAllocRefPos[ ed.m_vecNID[0] ] );
        Vec2 r1( m_pAllocRefPos[ ed.m_vecNID[1] ] );
        Vec2 r2( m_pAllocRefPos[ ed.m_vecNID[2] ] );
        // Specific element init
        // Linear
        Mat6x6 K;
        TriangleElement2::Compute_K( r0, r1, r2,
                                     m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                     K );
        // K0X
        m_pAllocLED[it_e].m_K00 = mal::GRange<0,0,1,1>( K );
        m_pAllocLED[it_e].m_K01 = mal::GRange<0,2,1,3>( K );
        m_pAllocLED[it_e].m_K02 = mal::GRange<0,4,1,5>( K );
        // K1X
        m_pAllocLED[it_e].m_K10 = mal::GRange<2,0,3,1>( K );
        m_pAllocLED[it_e].m_K11 = mal::GRange<2,2,3,3>( K );
        m_pAllocLED[it_e].m_K12 = mal::GRange<2,4,3,5>( K );
        // K2X
        m_pAllocLED[it_e].m_K20 = mal::GRange<4,0,5,1>( K );
        m_pAllocLED[it_e].m_K21 = mal::GRange<4,2,5,3>( K );
        m_pAllocLED[it_e].m_K22 = mal::GRange<4,4,5,5>( K );
    }
}

}}} //namespace S2::ms::fem
