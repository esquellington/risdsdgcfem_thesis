#include "TetrahedronElementSet3.h"
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

#define __USE_EXACT_FORCE_DIFFERENTIAL
#define REGULARIZE_CCM_PD_U_THRESHOLD 1e-6f

#define __USE_GSL
#ifdef __USE_GSL
#  include <gsl/gsl_poly.h>
#endif

#include <iostream>//TEMP
#include <Mal/GSerialization.h>

//#define __USE_OMP

namespace mal
{
//Solve ax^3 + bx^2 + cx + d = 0 returns x1 <= x2 <= x3
template<typename T> int GSolvePolynomialEq3( T a, T b, T c, T d, T& x1, T& x2, T& x3 )
{
#ifdef __USE_GSL
    //\see http://www.gnu.org/software/gsl/manual/html_node/Cubic-Equations.html#Cubic-Equations
    if( a != T(0) ) //TEMP: Should use some epsilon....
    {
        double r1,r2,r3;
        int num_roots = gsl_poly_solve_cubic( b/a, c/a, d/a, &r1, &r2, &r3 );
        x1 = r1;
        x2 = r2;
        x3 = r3;
        return num_roots;
    }
    else
    {
        x3 = 0;
        return GSolvePolynomialEq2( b, c, d, x1, x2 );
    }
#else
    return 0;
#endif
}
} //namespace mal

namespace S2 {
namespace ms {
namespace fem {

const Real cEpsilon_dR_Numerical( 0.0001 ); //10e-3 works better in transition areas, 10e-2 works better in SVD-critical-point, which is not present in PD_U

//! ElementCache, constant during a subset of element lifetime
struct ElementCache3
{
    // Degenerate
    TetrahedronElement3::DoC m_DoC; // local direction-of-collapse, either V-F or E-E
    Real m_DetF; // constant during CG iter, implicit in det(F) = det(Ds)/det(Dm)
    // CLFEM, PD_CCM (constant during CG iter)
    Mat3x3 m_R;
    // PD_CCM (constant during CG iter)
    //Mat3x3 m_S; //implicit in F = R*S
    // PK (constant during CG iter)
    Mat3x3 m_Ds;
    Mat3x3 m_F;
    //Mat2x2 m_U, m_DiagF, m_Vt; ??
    Mat3x3 m_InverseTransposeF; // used in implicit integration
#ifdef __ENABLE_PLASTICITY
    // Plasticity
    Vec3 m_Ep; //Plastic strain
#endif
};

class Cache3: public IForce3::ICache
{
public:
    Cache3( unsigned int num_elements, const Real* vec_inv_mass )
    : m_vecInvMass(vec_inv_mass), m_vecEC(0), m_vecPrevX(0), m_vecX(0), m_vecV(0), m_Time(0), m_TimeStep(0)
        {
            m_vecEC = new ElementCache3[num_elements];
        }
    ~Cache3()
        {
            if( m_vecEC ) delete[] m_vecEC;
            m_vecEC = 0;
        }
public:
    const Real *m_vecInvMass;
    ElementCache3 *m_vecEC;
    //\name Begin/EndEvaluation stuff
    //@{
    const Vec3 *m_vecPrevX;
    const Vec3 *m_vecX;
    const Vec3 *m_vecV;
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

TetrahedronElement3::DoC ComputeDoC( const Vec3 &a0, const Vec3 &a1, const Vec3 &a2, const Vec3 &a3,
                                     const Vec3 &b0, const Vec3 &b1, const Vec3 &b2, const Vec3 &b3,
                                     Real volume, Real degenerate_threshold_det_F );

inline void LameParameters_From_YoungAndPoisson( Real young_modulus, Real poisson_ratio,
                                                 Real &lame_mu, Real &lame_lambda )
{
    lame_mu = Real(0.5) * young_modulus / (1+poisson_ratio);
    lame_lambda = young_modulus * poisson_ratio / ( (Real(1)+poisson_ratio) * (Real(1)-Real(2)*poisson_ratio) );
}

TetrahedronElementSet3::TetrahedronElementSet3()
: m_NumNodes(0)
, m_NumElements(0)
, m_TotalVolume(0)
, m_vecRefPos(0)
, m_vecED(0)
, m_vecLED(0)
{
}

TetrahedronElementSet3::~TetrahedronElementSet3()
{
    ClearBakedData();
}

void TetrahedronElementSet3::SetBakedData( bool b_shared,
                                           uint32 num_nodes, uint32 num_elements,
                                           const Vec3 *vec_ref_pos, const ElementData *vec_ed, const LinearED *vec_led )
{
    // Save baked data refs
    m_NumNodes = num_nodes;
    m_NumElements = num_elements;
    m_vecRefPos = vec_ref_pos;
    m_vecED = vec_ed;
    m_vecLED = vec_led;
}

void TetrahedronElementSet3::ClearBakedData()
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
IForce3::ICache* TetrahedronElementSet3::CreateCache( const Real* vec_inv_mass )
{
    // Alloc cached data
    Cache3* p_cache = new Cache3( m_NumElements, vec_inv_mass );
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    // Init cached data
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( p_cache->m_vecEC[it_e] );
        ec.m_DoC = TetrahedronElement3::DoC();
        ec.m_R = Mat3x3::Identity();
        ec.m_F = Mat3x3::Identity();
        ec.m_DetF = Real(1);
        //\todo Ds, InverseTransposeF...
        //\todo Consider full init, including F, R, etc...
    }
    return p_cache;
}

void TetrahedronElementSet3::DestroyCache( IForce3::ICache* p_cache )
{
    if( p_cache ) delete p_cache;
}

//----------------------------------------------------------------
// IForce virtual API
//----------------------------------------------------------------
void TetrahedronElementSet3::BeginEvaluation( IForce3::ICache* p_cache, const Vec3 *vec_prev_x, const Vec3 *vec_x, const Vec3 *vec_v, Real t ) const
{
    S2_MS_TRACE_BEGIN_SCOPE( "TetrahedronElementSet3::BeginEvaluation" );

    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

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
        ElementCache3& ec( cache.m_vecEC[it_e] );
        Vec3 x0( vec_x[ ed.m_vecNID[0] ] );
        Vec3 x1( vec_x[ ed.m_vecNID[1] ] );
        Vec3 x2( vec_x[ ed.m_vecNID[2] ] );
        Vec3 x3( vec_x[ ed.m_vecNID[3] ] );
        // Ds
        TetrahedronElement3::Compute_D( x0, x1, x2, x3, ec.m_Ds );
        // F = Ds * Dm^-1
        ec.m_F = ec.m_Ds * ed.m_InvDm;
        ec.m_DetF = mal::Det( ec.m_F );
        ec.m_InverseTransposeF = (mal::Abs(ec.m_DetF) > 0.000001f ) //\todo This is the same epsilon inside Inverse()... but we should find a better way than copying the constant...
                                 ? mal::Inverse( mal::Transposed( ec.m_F ) )
                                 : Mat3x3::Identity();
        // Check NaN
        MS_ASSERT( !mal::IsNaN( ec.m_F ) );
        MS_ASSERT( !mal::IsNaN( ec.m_DetF ) );
        MS_ASSERT( !mal::IsNaN( ec.m_InverseTransposeF ) );

        S2_MS_STAT_ADD( cache.m_Stats.m_Num_Degenerate, ((ec.m_DetF <= 0) ? 1 : 0) );
        S2_MS_STAT_ADD( cache.m_Stats.m_Sum_Degenerate_DetF, ((ec.m_DetF <= 0) ? ec.m_DetF : 0) );
    }
    // Update inversion stuff: ec.m_DoC
    UpdateDegeneration( p_cache );

    // Update plastic flow, creep, whatever... \note Does nothing if m_TimeStep <= 0
    //\todo UpdatePlasticity( p_cache );

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

#ifdef __S2_MS_ENABLE_STATS
void TetrahedronElementSet3::ComputeStats( const IForce3::ICache* p_cache, TetrahedronElementSet3::Stats& stats ) const
{
    MS_ASSERT( 0 != p_cache );
    const Cache3& cache( *static_cast<const Cache3*>( p_cache ) );
    stats.m_Num_Degenerate = cache.m_Stats.m_Num_Degenerate;
    stats.m_Sum_Degenerate_DetF = cache.m_Stats.m_Sum_Degenerate_DetF;
}
#endif

// \note Forces are ACCUMULATED on vec_f, NOT SET into vec_f
void TetrahedronElementSet3::f( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    // Elastic + Friction + Plastic forces
    f_e( p_cache, vec_f );
    //f_d( p_cache, vec_f ); \todo Needs vec_node_inv_mass for Rayleigh M term...
    f_p( p_cache, vec_f );
}

// \note Forces are ACCUMULATED on vec_f, NOT SET into vec_f
void TetrahedronElementSet3::f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
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
void TetrahedronElementSet3::f_d( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    // Accumulate specific Elastic forces
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L: L_f_d( p_cache, vec_f ); break;
        // Corotational
    case Params::eMM_C_WRP: C_f_d( p_cache, vec_f ); break;
    case Params::eMM_C_LCM:
    case Params::eMM_C_LCMH:
    case Params::eMM_C_CCM:
        C_f_d( p_cache, vec_f );
        break; //\todo HACK, we use C_f_d by now, because R IS available after H_CCM_PD_f() evaluations... ugly but works
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

void TetrahedronElementSet3::f_p( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    //\todo
}

void TetrahedronElementSet3::df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const
{
    // Accumulate specific Elastic forces
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L: L_df_x( p_cache, vec_dx, vec_df ); break;
        // Corotational
    case Params::eMM_C_WRP: C_df_x( p_cache, vec_dx, vec_df ); break;
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

void TetrahedronElementSet3::df_v( IForce3::ICache* p_cache, const Vec3 *vec_dv, Vec3 *vec_df ) const
{
}
Real TetrahedronElementSet3::V( IForce3::ICache* p_cache ) const
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
void TetrahedronElementSet3::EndEvaluation( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Unset evaluation-constant params
    cache.m_vecPrevX = 0;
    cache.m_vecX = 0;
    cache.m_vecV = 0;
#ifdef __S2_MS_ENABLE_STATS
    cache.m_Stats.End();
#endif
}

void TetrahedronElementSet3::GetNID( unsigned int eid,
                                     unsigned int &nid0, unsigned int &nid1, unsigned int &nid2, unsigned int &nid3 ) const
{
    nid0 = m_vecED[eid].m_vecNID[0];
    nid1 = m_vecED[eid].m_vecNID[1];
    nid2 = m_vecED[eid].m_vecNID[2];
    nid3 = m_vecED[eid].m_vecNID[3];
}

Transform3 TetrahedronElementSet3::GetTransformE2W( const IForce3::ICache* p_cache, unsigned int eid, const Vec3* vec_pos ) const
{
    MS_ASSERT( 0 != p_cache );
    const Cache3& cache( *static_cast<const Cache3*>( p_cache ) );

    Transform3 Te2w;
    switch( m_Params.m_MM )
    {
        // Linear
    case Params::eMM_L:
        Te2w.m_Rot = Mat3x3::Identity();
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
        Te2w.m_Rot = Mat3x3::Identity();
        break;
    }
    const ElementData &ed( m_vecED[eid] );
    Te2w.m_Pos = mal::Rcp( Real(4) )*( vec_pos[ed.m_vecNID[0]] + vec_pos[ed.m_vecNID[1]] + vec_pos[ed.m_vecNID[2]] + vec_pos[ed.m_vecNID[3]] );

    return Te2w;
}

TetrahedronElement3::DoC TetrahedronElementSet3::GetDoC( const IForce3::ICache* p_cache, unsigned int eid ) const
{
    MS_ASSERT( 0 != p_cache );
    const Cache3& cache( *static_cast<const Cache3*>( p_cache ) );
    return cache.m_vecEC[eid].m_DoC;
}

//----------------------------------------------------------------
// Internal methods
//----------------------------------------------------------------

/*IMPORTANT: see TriangleElementSet2 for important info on this method */
void TetrahedronElementSet3::UpdateDegeneration( IForce3::ICache* p_cache ) const
{
#define __NEW_CROSSING_DETECTION //\todo should be default, see TriangleElementSet2
#ifdef __NEW_CROSSING_DETECTION

    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __USE_OMP
#  pragma omp parallel for
#endif

    // Check all elements for degeneration according to det(F)
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //remains/becomes undegenerate
            ec.m_DoC = TetrahedronElement3::DoC(); //invalid
        else //remains/becomes degenerate
        {
            const ElementData &ed( m_vecED[it_e] );
            Vec3 a0( cache.m_vecPrevX[ ed.m_vecNID[0] ] );
            Vec3 a1( cache.m_vecPrevX[ ed.m_vecNID[1] ] );
            Vec3 a2( cache.m_vecPrevX[ ed.m_vecNID[2] ] );
            Vec3 a3( cache.m_vecPrevX[ ed.m_vecNID[3] ] );
            Mat3x3 prevDs;
            TetrahedronElement3::Compute_D( a0, a1, a2, a3, prevDs );
            if( mal::Det(prevDs*ed.m_InvDm) > m_Params.m_DegenerateThresholdDetF ) // was undegenerate (*)
            {
                //compute DoC
                Vec3 b0( cache.m_vecX[ ed.m_vecNID[0] ] );
                Vec3 b1( cache.m_vecX[ ed.m_vecNID[1] ] );
                Vec3 b2( cache.m_vecX[ ed.m_vecNID[2] ] );
                Vec3 b3( cache.m_vecX[ ed.m_vecNID[3] ] );
                ec.m_DoC = ComputeDoC( a0, a1, a2, a3,
                                       b0, b1, b2, b3,
                                       ed.m_Volume,
                                       m_Params.m_DegenerateThresholdDetF );
                if( ec.m_DoC.IsInvalid() )
                    MS_LOG_ERROR( "New Degenerate element %u CANNOT SOLVE Eq3", it_e );
            }
            else if( ec.m_DoC.IsInvalid() ) // remains degenerate but discrete
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
            //else: remains degenerate with the same DoC
        }
    }

#else

    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __USE_OMP
#  pragma omp parallel for
#endif

    // Check all elements for degeneration according to det(F)
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( cache.m_vecEC[it_e] );
        const ElementData &ed( m_vecED[it_e] );
        Vec3 a0( cache.m_vecPrevX[ ed.m_vecNID[0] ] );
        Vec3 a1( cache.m_vecPrevX[ ed.m_vecNID[1] ] );
        Vec3 a2( cache.m_vecPrevX[ ed.m_vecNID[2] ] );
        Vec3 a3( cache.m_vecPrevX[ ed.m_vecNID[3] ] );
        Mat3x3 prevDs;
        TetrahedronElement3::Compute_D( a0, a1, a2, a3, prevDs );
        if( ec.m_DoC.IsInvalid() // No cached degeneration
            ||
            mal::Det(prevDs*ed.m_InvDm) > m_Params.m_DegenerateThresholdDetF ) // was undegenerate (*)
        {
            if( ec.m_DetF <= m_Params.m_DegenerateThresholdDetF ) //becomes degenerate
            {
                //compute DoC
                Vec3 b0( cache.m_vecX[ ed.m_vecNID[0] ] );
                Vec3 b1( cache.m_vecX[ ed.m_vecNID[1] ] );
                Vec3 b2( cache.m_vecX[ ed.m_vecNID[2] ] );
                Vec3 b3( cache.m_vecX[ ed.m_vecNID[3] ] );
                ec.m_DoC = ComputeDoC( a0, a1, a2, a3,
                                       b0, b1, b2, b3,
                                       ed.m_Volume,
                                       m_Params.m_DegenerateThresholdDetF );
                if( ec.m_DoC.IsInvalid() )
                    MS_LOG_ERROR( "New Degenerate element %d CANNOT SOLVE Eq3", it_e );
            }
            //else remains undegenerate
        }
        else //was degenerate AND has a valid NoC
        {
            if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //becomes undegenerate
            {
                ec.m_DoC = TetrahedronElement3::DoC(); //invalid
            }
            //else remains degenerate with same DoC
        }
    }
#endif
}

/*! Compute Direction-of-Collapse
  \note Internally computes ToC
*/
TetrahedronElement3::DoC ComputeDoC( const Vec3 &a0, const Vec3 &a1, const Vec3 &a2, const Vec3 &a3,
                                     const Vec3 &b0, const Vec3 &b1, const Vec3 &b2, const Vec3 &b3,
                                     Real volume, Real degenerate_threshold_det_F )
{
    TetrahedronElement3::DoC doc; //Invalid by default
    // Gather node displacements
    Vec3 d0( b0 - a0 );
    Vec3 d1( b1 - a1 );
    Vec3 d2( b2 - a2 );
    Vec3 d3( b3 - a3 );
    // Compute displacement differences
    Vec3 d10( d1 - d0 );
    Vec3 d20( d2 - d0 );
    Vec3 d30( d3 - d0 );
    // Compute position differences
    Vec3 a10( a1 - a0 );
    Vec3 a20( a2 - a0 );
    Vec3 a30( a3 - a0 );
    // Compute ToC, when det(F) = m_DegenerateThresholdDetF => det(Ds) = m_DegenerateThresholdDetF * det(Dm) => det(Ds) = m_DegenerateThresholdDetF * 6 * Volume(Dm)
    Real K( degenerate_threshold_det_F * Real(6) * volume );
    Real A( + d10.x() * d20.y() * d30.z()
            + d20.x() * d30.y() * d10.z()
            + d30.x() * d10.y() * d20.z()
            - d30.x() * d20.y() * d10.z()
            - d20.x() * d10.y() * d30.z()
            - d10.x() * d30.y() * d20.z() );

    Real B( + ( a10.x() * d20.y() * d30.z() + d10.x() * a20.y() * d30.z() + d10.x() * d20.y() * a30.z() )
            + ( a20.x() * d30.y() * d10.z() + d20.x() * a30.y() * d10.z() + d20.x() * d30.y() * a10.z() )
            + ( a30.x() * d10.y() * d20.z() + d30.x() * a10.y() * d20.z() + d30.x() * d10.y() * a20.z() )
            - ( a30.x() * d20.y() * d10.z() + d30.x() * a20.y() * d10.z() + d30.x() * d20.y() * a10.z() )
            - ( a20.x() * d10.y() * d30.z() + d20.x() * a10.y() * d30.z() + d20.x() * d10.y() * a30.z() )
            - ( a10.x() * d30.y() * d20.z() + d10.x() * a30.y() * d20.z() + d10.x() * d30.y() * a20.z() ) );

    Real C( + ( a10.x() * a20.y() * d30.z() + a10.x() * d20.y() * a30.z() + d10.x() * a20.y() * a30.z() )
            + ( a20.x() * a30.y() * d10.z() + a20.x() * d30.y() * a10.z() + d20.x() * a30.y() * a10.z() )
            + ( a30.x() * a10.y() * d20.z() + a30.x() * d10.y() * a20.z() + d30.x() * a10.y() * a20.z() )
            - ( a30.x() * a20.y() * d10.z() + a30.x() * d20.y() * a10.z() + d30.x() * a20.y() * a10.z() )
            - ( a20.x() * a10.y() * d30.z() + a20.x() * d10.y() * a30.z() + d20.x() * a10.y() * a30.z() )
            - ( a10.x() * a30.y() * d20.z() + a10.x() * d30.y() * a20.z() + d10.x() * a30.y() * a20.z() ) );

    Real D( + a10.x() * a20.y() * a30.z()
            + a20.x() * a30.y() * a10.z()
            + a30.x() * a10.y() * a20.z()
            - a30.x() * a20.y() * a10.z()
            - a20.x() * a10.y() * a30.z()
            - a10.x() * a30.y() * a20.z()
            - K );

    Real vec_toc[3];
    int num_roots( mal::GSolvePolynomialEq3<Real>( A, B, C, D, vec_toc[0], vec_toc[1], vec_toc[2] ) );
    if( num_roots > 0 )
    {
        Real toc(2); //no toc > 1 is accepted
        int num_roots01(0);
        for( int i=0; i<num_roots; i++ )
        {
            if( vec_toc[i] > Real(0) && vec_toc[i] <= Real(1) && vec_toc[i] < toc )
            {
                toc = vec_toc[i];
                num_roots01++;
            }
        }
        /*
        MS_LOG( "Potentially degenerated with toc = %f, %f, %f, #roots = %d, #roots01 = %d, min_toc = %f",
                vec_toc[0], vec_toc[1], vec_toc[2], num_roots, num_roots01, toc );
        */
        if( num_roots01 != 1 && num_roots01 != 3 ) return doc; //TEMPORAL TO TEST
        MS_ASSERT( num_roots01 == 1 || num_roots01 == 3 );
        MS_ASSERT( Real(0) < toc && toc <= Real(1) );
        //\todo We knot that det(F(0)) > DTDF, and that det(F(1)) <= DTDF, so there MUST BE STRICTLY 1 or 3 crossing in toc = [0,1]
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
            Vec3 q0( a0 + toc*d0 );
            Vec3 q1( a1 + toc*d1 );
            Vec3 q2( a2 + toc*d2 );
            Vec3 q3( a3 + toc*d3 );
            // Compute collapse normal using the smallest singular-value from SVD
            /*\todo THIS SEEMS TO WORK FINE, but I'm not sure if it's
              justified... the matrix is not exaclty the
              covariance... */
            Mat3x3 Q( mal::GMat3x3_From_Columns( q1-q0, q2-q0, q3-q0 ) );
            Mat3x3 U,Vt;
            Vec3 S;
            mal::GSingularValueDecomposition_USVt( Q, U, S, Vt );
            Vec3 N = mal::GRow<2>( Vt ); //Singular Vector from the smallest Singular Value
            /* \todo non-unitary normal computed by adding "some"
               non-unitary face normals SHOULD also work, and would
               avoid expensive SVD...
            for( )
            {
                if( mal::Dot(edge_ij,N) < 0 )
                    N -= cross;
                else
                    N += cross;
            }
            */
            // Compute face areas (actually, per-face unit-height prism volumes, which is the same but avoids sqrt()
            Real vec_area[4]; //index is actually the VertexInTetrahedron of the vertex opposite to the actual face
            vec_area[0] = mal::Abs( mal::Dot( N, mal::Cross( q2-q1, q3-q1 ) ) );
            vec_area[1] = mal::Abs( mal::Dot( N, mal::Cross( q2-q0, q3-q0 ) ) );
            vec_area[2] = mal::Abs( mal::Dot( N, mal::Cross( q1-q0, q3-q0 ) ) );
            vec_area[3] = mal::Abs( mal::Dot( N, mal::Cross( q1-q0, q2-q0 ) ) );
            Real sum_area( vec_area[0] + vec_area[1] + vec_area[2] + vec_area[3] );
            // Check V-F case by searching for a face with Area_i = 1/2 \sum_j Area_j
            //\todo Precision may be critical here, consider doing FIRST E-E tests and selecting the V-F case with F_i CLOSEST TO 1/2 \sum Area_i otherwise
            for( int it_fit=0; it_fit<4; it_fit++ )
                if( mal::ApproxEq( vec_area[it_fit], Real(0.5)*sum_area ) )
                    doc = TetrahedronElement3::DoC( it_fit );
            // If not found, we have an E-E case
            if( doc.IsInvalid() )
            {
                // Determine which outgoing edge from V_0 crosses another one on the opposite face F_0
                if( mal::Dot( mal::Cross(q1-q0, q2-q0), mal::Cross(q1-q0, q3-q0) ) < Real(0) ) //i,j,k,l = 0,1,2,3
                    doc = TetrahedronElement3::DoC(0,1);
                else if( mal::Dot( mal::Cross(q2-q0, q1-q0), mal::Cross(q2-q0, q3-q0) ) < Real(0) ) //i,j,k,l = 0,2,1,3
                    doc = TetrahedronElement3::DoC(0,2);
                else if( mal::Dot( mal::Cross(q3-q0, q1-q0), mal::Cross(q3-q0, q2-q0) ) < Real(0) ) //i,j,k,l = 0,3,1,2
                    doc = TetrahedronElement3::DoC(0,3);
                else
                {
                    MS_LOG_ERROR("MISSED POTENTIAL E-E crossing, forcing DoC = V-F 0");
                    doc = TetrahedronElement3::DoC(0);
                }
            }
            else
            {
                //MS_LOG("V-F crossing DoC = %d", doc.m_VIT0 );
            }
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
        MS_LOG_ERROR( "Cannot solve A*t^3 + B*t^2 + C*t + D = 0 for A = %f, B = %f, C = %f, D = %f, #roots = %d", A, B, C, D, num_roots );
    }
    return doc;
}

//----------------------------------------------------------------
// COROTATIONAL

void TetrahedronElementSet3::ComputeRotations_Id( IForce3::ICache *p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( cache.m_vecEC[it_e] );
        ec.m_R = Mat3x3::Identity();
    }
}

void TetrahedronElementSet3::ComputeRotations_QR( IForce3::ICache *p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( cache.m_vecEC[it_e] );
        ec.m_R = Compute_R_QR( ec.m_F );
        MS_ASSERT( !mal::IsNaN( ec.m_R ) );
    }
}

void TetrahedronElementSet3::ComputeRotations_MSVD( IForce3::ICache *p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( cache.m_vecEC[it_e] );
        ec.m_R = Compute_R_SVD( ec.m_F, ec.m_DetF, m_Params.m_DegenerateThresholdDetF );
        MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        MS_ASSERT( mal::Det( ec.m_R ) >= 0 );
    }
}

void TetrahedronElementSet3::ComputeRotations_PD_Project( IForce3::ICache *p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //Undegenerate
        {
            ec.m_R = Compute_R_PD( ec.m_F, ec.m_DetF );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else if( ec.m_DoC.IsValid() ) //Degenerate and has proper DoC
        {
            const ElementData& ed( m_vecED[it_e] );
            ec.m_R = Compute_R_PDP_Degenerate( ec.m_F, ec.m_DetF, m_Params.m_DegenerateThresholdDetF,
                                               ec.m_DoC,
                                               cache.m_vecX[ ed.m_vecNID[0] ], cache.m_vecX[ ed.m_vecNID[1] ], cache.m_vecX[ ed.m_vecNID[2] ], cache.m_vecX[ ed.m_vecNID[3] ],
                                               ed.m_Volume, ed.m_InvDm );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else //Degenerate but no proper DoC
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element %u has invalid DoC. Using R=Id", it_e );
            ec.m_R = Mat3x3::Identity();
        }
    }
}

void TetrahedronElementSet3::ComputeRotations_PD_Reflect( IForce3::ICache *p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //Undegenerate
        {
            ec.m_R = mal::GRotation3x3_PolarDecomposition( ec.m_F, ec.m_DetF );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else if( ec.m_DoC.IsValid() ) //Degenerate and has proper DoC
        {
            const ElementData& ed( m_vecED[it_e] );
            ec.m_R = Compute_R_PDR_Degenerate( ec.m_F, ec.m_DetF, m_Params.m_DegenerateThresholdDetF,
                                               ec.m_DoC,
                                               cache.m_vecX[ ed.m_vecNID[0] ], cache.m_vecX[ ed.m_vecNID[1] ], cache.m_vecX[ ed.m_vecNID[2] ], cache.m_vecX[ ed.m_vecNID[3] ] );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else //Degenerate but no proper DoC
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element %u has invalid DoC. Using R=Id", it_e );
            ec.m_R = Mat3x3::Identity();
        }
    }
}


void TetrahedronElementSet3::ComputeRotations_DAPD( IForce3::ICache *p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        ElementCache3& ec( cache.m_vecEC[it_e] );
        if( ec.m_DetF > m_Params.m_DegenerateThresholdDetF ) //Undegenerate
        {
            ec.m_R = mal::GRotation3x3_PolarDecomposition( ec.m_F, ec.m_DetF );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else if( ec.m_DoC.IsValid() ) //Degenerate and has proper DoC
        {
            const ElementData& ed( m_vecED[it_e] );
            ec.m_R = Compute_R_DAPD_Degenerate( ec.m_F, ec.m_DetF, m_Params.m_DegenerateThresholdDetF,
                                                ec.m_DoC,
                                                cache.m_vecX[ ed.m_vecNID[0] ], cache.m_vecX[ ed.m_vecNID[1] ], cache.m_vecX[ ed.m_vecNID[2] ], cache.m_vecX[ ed.m_vecNID[3] ],
                                                ed.m_Volume, ed.m_InvDm,
                                                m_Params.m_DAPD_L_Factor, m_Params.m_DAPD_NL_Factor, m_Params.m_DAPD_NL_Exponent );
            MS_ASSERT( !mal::IsNaN( ec.m_R ) );
        }
        else //Degenerate but no proper DoC
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element %u has invalid DoC. Using R=Id", it_e );
            ec.m_R = Mat3x3::Identity();
        }
    }
}

static bool Compute_dR_From_RM( const Params& params,
                                const Cache3& cache,
                                const TetrahedronElementSet3::ElementData& ed, const ElementCache3& ec,
                                const Mat3x3& Rt, const Vec3* vec_dx,
                                Mat3x3& dR )
{
    bool bCorrect( false );
    if( params.m_DM != Params::eDM_Truncated )
    {
        TetrahedronElementSet3::node_index_type nid0( ed.m_vecNID[0] );
        TetrahedronElementSet3::node_index_type nid1( ed.m_vecNID[1] );
        TetrahedronElementSet3::node_index_type nid2( ed.m_vecNID[2] );
        TetrahedronElementSet3::node_index_type nid3( ed.m_vecNID[3] );
        if( params.m_DM == Params::eDM_Exact
            || params.m_DM == Params::eDM_Exact_Inv )
        {
            if( params.m_RM == Params::eRM_DAPD )
            {
                if( params.m_DM == Params::eDM_Exact_Inv ) //We use Inv as it should work regardless of element degeneration state
                    bCorrect = Compute_dR_DAPD( params.m_DegenerateThresholdDetF, params.m_DAPD_L_Factor, params.m_DAPD_NL_Factor, params.m_DAPD_NL_Exponent,
                                                ec.m_DoC,
                                                cache.m_vecX[nid0], cache.m_vecX[nid1], cache.m_vecX[nid2], cache.m_vecX[nid3],
                                                vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid3],
                                                ed.m_InvDm, ec.m_F, ec.m_R, Rt, ed.m_Volume,
                                                dR );
                else if( ec.m_DetF > params.m_DegenerateThresholdDetF )
                    bCorrect = Compute_dR( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid3],
                                           ed.m_InvDm, ec.m_F, ec.m_R, Rt,
                                           dR );
            }
            else if( params.m_RM == Params::eRM_MSVD ) //McAdams dR is only correct for MSVD inversion handling
                bCorrect = Compute_dR( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid3],
                                       ed.m_InvDm, ec.m_F, ec.m_R, Rt,
                                       dR );
            else if( params.m_RM == Params::eRM_QR )
                bCorrect = false; /*Compute_dR_QR_YX( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid3],
                                             ed.m_InvDm, ec.m_F,
                                             dR );*/
        }
        else if( params.m_DM == Params::eDM_Numerical )
        {
            if( params.m_RM == Params::eRM_DAPD ) //Accepts Inv
                bCorrect = Compute_dR_DAPD_Numerical( params.m_DegenerateThresholdDetF, params.m_DAPD_L_Factor, params.m_DAPD_NL_Factor, params.m_DAPD_NL_Exponent,
                                                      ec.m_DoC,
                                                      cache.m_vecX[nid0], cache.m_vecX[nid1], cache.m_vecX[nid2], cache.m_vecX[nid3],
                                                      vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid3],
                                                      ed.m_InvDm, ed.m_Volume,
                                                      cEpsilon_dR_Numerical,
                                                      dR );
            else if( params.m_RM == Params::eRM_MSVD
                     || ec.m_DetF > params.m_DegenerateThresholdDetF )
                bCorrect = Compute_dR_Numerical( cache.m_vecX[nid0], cache.m_vecX[nid1], cache.m_vecX[nid2], cache.m_vecX[nid3],
                                                 vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid3],
                                                 ed.m_InvDm, params.m_DegenerateThresholdDetF,
                                                 cEpsilon_dR_Numerical,
                                                 dR );
        }
    }
    return bCorrect;
}

//----------------------------------------------------------------
// Linear
void TetrahedronElementSet3::L_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

//\todo Accumulation may fail if 2 threads write on the same vec_f[i]
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        Vec3 u0( cache.m_vecX[nid0] - m_vecRefPos[nid0] ); //\todo Consider precomputing vec_u to avoid repeating x_i-r_i
        Vec3 u1( cache.m_vecX[nid1] - m_vecRefPos[nid1] );
        Vec3 u2( cache.m_vecX[nid2] - m_vecRefPos[nid2] );
        Vec3 u3( cache.m_vecX[nid3] - m_vecRefPos[nid3] );
#ifdef __USE_OMP
        Vec3 f0( led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 + led.m_K03 * u3 );
        Vec3 f1( led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 + led.m_K13 * u3 );
        Vec3 f2( led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 + led.m_K23 * u3 );
        Vec3 f3( led.m_K30 * u0 + led.m_K31 * u1 + led.m_K32 * u2 + led.m_K33 * u3 );
#  pragma omp critical
        vec_f[ nid0 ] -= f0;
        vec_f[ nid1 ] -= f1;
        vec_f[ nid2 ] -= f2;
        vec_f[ nid3 ] -= f3;
#else
        vec_f[ nid0 ] -= led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 + led.m_K03 * u3;
        vec_f[ nid1 ] -= led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 + led.m_K13 * u3;
        vec_f[ nid2 ] -= led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 + led.m_K23 * u3;
        vec_f[ nid3 ] -= led.m_K30 * u0 + led.m_K31 * u1 + led.m_K32 * u2 + led.m_K33 * u3;
#endif
    }
}

/*! Linear Rayleigh Damping:
  C_e = m_RayleighCoeff[0] * M_e + m_RayleighCoeff[1] * K_e;
  fd_e = -C_e * V_e
*/
void TetrahedronElementSet3::L_f_d( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

//\todo Accumulation may fail if 2 threads write on the same vec_f[i]
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    // Acc per-element damping forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        Vec3 v0( cache.m_vecV[nid0] );
        Vec3 v1( cache.m_vecV[nid1] );
        Vec3 v2( cache.m_vecV[nid2] );
        Vec3 v3( cache.m_vecV[nid3] );
#ifdef __USE_OMP
        Vec3 f0( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid0]) ) * v0
                 + m_Params.m_RayleighCoeff[1] * (led.m_K00 * v0 + led.m_K01 * v1 + led.m_K02 * v2 + led.m_K03 * v3) );
        Vec3 f1( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid1]) ) * v1
                 + m_Params.m_RayleighCoeff[1] * (led.m_K10 * v0 + led.m_K11 * v1 + led.m_K12 * v2 + led.m_K13 * v3) );
        Vec3 f2( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid2]) ) * v2
                 + m_Params.m_RayleighCoeff[1] * (led.m_K20 * v0 + led.m_K21 * v1 + led.m_K22 * v2 + led.m_K23 * v3) );
        Vec3 f3( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid3]) ) * v3
                 + m_Params.m_RayleighCoeff[1] * (led.m_K30 * v0 + led.m_K31 * v1 + led.m_K32 * v2 + led.m_K33 * v3) );
#  pragma omp critical
        vec_f[ nid0 ] -= f0;
        vec_f[ nid1 ] -= f1;
        vec_f[ nid2 ] -= f2;
        vec_f[ nid3 ] -= f3;
#else
        vec_f[ nid0 ] -= ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid0]) ) * v0
                         + m_Params.m_RayleighCoeff[1] * (led.m_K00 * v0 + led.m_K01 * v1 + led.m_K02 * v2 + led.m_K03 * v3);
        vec_f[ nid1 ] -= ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid1]) ) * v1
                         + m_Params.m_RayleighCoeff[1] * (led.m_K10 * v0 + led.m_K11 * v1 + led.m_K12 * v2 + led.m_K13 * v3);
        vec_f[ nid2 ] -= ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid2]) ) * v2
                         + m_Params.m_RayleighCoeff[1] * (led.m_K20 * v0 + led.m_K21 * v1 + led.m_K22 * v2 + led.m_K23 * v3);
        vec_f[ nid3 ] -= ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid3]) ) * v3
                         + m_Params.m_RayleighCoeff[1] * (led.m_K30 * v0 + led.m_K31 * v1 + led.m_K32 * v2 + led.m_K33 * v3);
#endif
    }
}

//----------------------------------------------------------------
// Corotational
void TetrahedronElementSet3::C_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

//\todo Accumulation may fail if 2 threads write on the same vec_f[i]
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        // Rotate world frame X to reference frame (no need to translate)
        Mat3x3 Rt( ec.m_R.Transposed() );
        Vec3 u0( Rt * cache.m_vecX[nid0] - m_vecRefPos[nid0] );
        Vec3 u1( Rt * cache.m_vecX[nid1] - m_vecRefPos[nid1] );
        Vec3 u2( Rt * cache.m_vecX[nid2] - m_vecRefPos[nid2] );
        Vec3 u3( Rt * cache.m_vecX[nid3] - m_vecRefPos[nid3] );
#ifdef __USE_OMP
        Vec3 f0( ec.m_R * ( led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 + led.m_K03 * u3 ) );
        Vec3 f1( ec.m_R * ( led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 + led.m_K13 * u3 ) );
        Vec3 f2( ec.m_R * ( led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 + led.m_K23 * u3 ) );
        Vec3 f3( ec.m_R * ( led.m_K30 * u0 + led.m_K31 * u1 + led.m_K32 * u2 + led.m_K33 * u3 ) );
#  pragma omp critical
        vec_f[ nid0 ] -= f0;
        vec_f[ nid1 ] -= f1;
        vec_f[ nid2 ] -= f2;
        vec_f[ nid3 ] -= f3;
#else
        // Compute elastic forces and rotate to world frame
        vec_f[ nid0 ] -= ec.m_R * ( led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 + led.m_K03 * u3 );
        vec_f[ nid1 ] -= ec.m_R * ( led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 + led.m_K13 * u3 );
        vec_f[ nid2 ] -= ec.m_R * ( led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 + led.m_K23 * u3 );
        vec_f[ nid3 ] -= ec.m_R * ( led.m_K30 * u0 + led.m_K31 * u1 + led.m_K32 * u2 + led.m_K33 * u3 );
#endif
    }
}

/*! Corotational Rayleigh Damping:
  C_e = m_RayleighCoeff[0] * M_e + m_RayleighCoeff[1] * K_e;
  fd_e = - R_e * ( C_e * R_e^T * V_e )
*/
void TetrahedronElementSet3::C_f_d( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

//\todo Accumulation may fail if 2 threads write on the same vec_f[i]
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    // Acc per-element damping forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        Mat3x3 Rt( ec.m_R.Transposed() );
        Vec3 v0( Rt * cache.m_vecV[nid0] );
        Vec3 v1( Rt * cache.m_vecV[nid1] );
        Vec3 v2( Rt * cache.m_vecV[nid2] );
        Vec3 v3( Rt * cache.m_vecV[nid3] );
#ifdef __USE_OMP
        Vec3 f0( ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid0]) ) * v0
                            + m_Params.m_RayleighCoeff[1] * (led.m_K00 * v0 + led.m_K01 * v1 + led.m_K02 * v2 + led.m_K03 * v3) ) );
        Vec3 f1( ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid1]) ) * v1
                            + m_Params.m_RayleighCoeff[1] * (led.m_K10 * v0 + led.m_K11 * v1 + led.m_K12 * v2 + led.m_K13 * v3) ) );
        Vec3 f2( ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid2]) ) * v2
                            + m_Params.m_RayleighCoeff[1] * (led.m_K20 * v0 + led.m_K21 * v1 + led.m_K22 * v2 + led.m_K23 * v3) ) );
        Vec3 f3( ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid3]) ) * v3
                            + m_Params.m_RayleighCoeff[1] * (led.m_K30 * v0 + led.m_K31 * v1 + led.m_K32 * v2 + led.m_K33 * v3) ) );
#  pragma omp critical
        vec_f[ nid0 ] -= f0;
        vec_f[ nid1 ] -= f1;
        vec_f[ nid2 ] -= f2;
        vec_f[ nid3 ] -= f3;
#else
        vec_f[ nid0 ] -= ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid0]) ) * v0
                                    + m_Params.m_RayleighCoeff[1] * (led.m_K00 * v0 + led.m_K01 * v1 + led.m_K02 * v2 + led.m_K03 * v3) );
        vec_f[ nid1 ] -= ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid1]) ) * v1
                                    + m_Params.m_RayleighCoeff[1] * (led.m_K10 * v0 + led.m_K11 * v1 + led.m_K12 * v2 + led.m_K13 * v3) );
        vec_f[ nid2 ] -= ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid2]) ) * v2
                                    + m_Params.m_RayleighCoeff[1] * (led.m_K20 * v0 + led.m_K21 * v1 + led.m_K22 * v2 + led.m_K23 * v3) );
        vec_f[ nid3 ] -= ec.m_R * ( ( m_Params.m_RayleighCoeff[0] * mal::Rcp(cache.m_vecInvMass[nid3]) ) * v3
                                    + m_Params.m_RayleighCoeff[1] * (led.m_K30 * v0 + led.m_K31 * v1 + led.m_K32 * v2 + led.m_K33 * v3) );
#endif
    }
}

// P-K Forces
void TetrahedronElementSet3::C_LCM_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute P
        Mat3x3 S( ec.m_R.Transposed() * ec.m_F ); //F = R*S
        Mat3x3 S_minus_Id( S-Mat3x3::Identity() ); //\todo THIS IS NOT CORRECT FOR DAPD, as S is NOT SYMMETRIC!
        Mat3x3 P( ec.m_R * ( 2 * m_LameMu * S_minus_Id
                             + m_LameLambda * mal::Trace(S_minus_Id) * Mat3x3::Identity() ) );

        // Compute PK force
        Mat3x3 H( - ed.m_Volume * P * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H) - mal::GColumn<2>(H); //h0 = -h1-h2-h3;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
        vec_f[ ed.m_vecNID[3] ] += mal::GColumn<2>(H);
    }
}

void TetrahedronElementSet3::C_CCM_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute \bar S from \bar R and F
        Mat3x3 S( ec.m_R.Transposed() * ec.m_F ); //F = R*S
        Real det_S( mal::Det(S) ); //\note == det(F) \todo Ensure and subst det_S by ec.m_DetF if correct

        // Compute P
        Mat3x3 P( ec.m_R * ( 2*m_LameMu*(S-Mat3x3::Identity())
                             + m_LameLambda * (det_S-1) * mal::Transposed( mal::Adjugate( S ) ) ) );

        // Compute PK force
        Mat3x3 H( - ed.m_Volume * P * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H) - mal::GColumn<2>(H); //h0 = -h1-h2-h3;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
        vec_f[ ed.m_vecNID[3] ] += mal::GColumn<2>(H);
    }
}

//----------------------------------------------------------------
// P-K
inline void ComputeSVD_Invertible( const Mat3x3 &F, Mat3x3 &U, Vec3 &vec_diag_F, Mat3x3 &Vt )
{
    mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt ); // Automagically forces pure rotation U,Vt by fixing potential inversion/reflection
}

//---- ComputeDiagonalP_XXX
inline Vec3 ComputeDiagonalP_LCM( const Vec3 &vec_diag_F, Real lame_mu, Real lame_lambda )
{
#ifdef __TO_PORT_TO_3D
    /* TEMP: Literal implementation from ITF, slower
       Mat3x3r diag_F( vec_diag_F[0], 0,
       0, vec_diag_F[1] );
       Mat3x3r E( diag_F - Mat3x3r::Identity() );
       Mat3x3r diag_P( 2*m_LameMu * E
       + m_LameLambda * mal::Trace( E ) * Mat3x3r::Identity() );
    */
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Real trE( vec_diag_F[0] + vec_diag_F[1] - 2 );
    return Vec3( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * trE,
                 2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * trE );
#else
    return Vec3::Zero();
#endif
}

inline Vec3 ComputeDiagonalP_CCM( const Vec3 &vec_diag_F, Real lame_mu, Real lame_lambda )
{
#ifdef __TO_PORT_TO_3D
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Real J( det_F );
    return Vec3( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * (J-1) * vec_diag_F[1],
                 2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * (J-1) * vec_diag_F[0] );
#else
    return Vec3::Zero();
#endif
}

inline Vec3 ComputeDiagonalP_NHC0( const Vec3 &vec_diag_F, Real lame_mu, Real lame_lambda, Real ecie_e )
{
#ifdef __TO_PORT_TO_3D
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Vec3 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                             mal::Max( ecie_e, vec_diag_F[1] ) );
    Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
    return Vec3( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ),
                 lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
#else
    return Vec3::Zero();
#endif
}

inline Vec3 ComputeDiagonalP_NHC1( const Vec3 &vec_diag_F, Real lame_mu, Real lame_lambda, Real ecie_e, Real ecie_k )
{
#ifdef __TO_PORT_TO_3D
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    bool bIsDegenerate0( vec_diag_F[0] < ecie_e );
    bool bIsDegenerate1( vec_diag_F[1] < ecie_e );
    Vec3 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                             mal::Max( ecie_e, vec_diag_F[1] ) );
    Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
    Vec3 vec_diag_P( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ),
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
#else
    return Vec3::Zero();
#endif
}

// P-K Forces
void TetrahedronElementSet3::H_LCM_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __TO_PORT_TO_3D
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        ec.m_R = U*Vt; //TEMP: for debug viz only

        // Compute diag_P
        Vec3 vec_diag_P( ComputeDiagonalP_LCM( vec_diag_F, m_LameMu, m_LameLambda ) );
        Mat3x3 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );

        // Compute PK force
        Mat3x3 H( - ed.m_Volume * U * diag_P * Vt * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
#endif
}

void TetrahedronElementSet3::H_CCM_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __TO_PORT_TO_3D
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        ec.m_R = U*Vt; //TEMP: for debug viz only

        // Compute diag_P
        Vec3 vec_diag_P( ComputeDiagonalP_CCM( vec_diag_F, m_LameMu, m_LameLambda ) );
        Mat3x3 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );

        // Compute PK force
        Mat3x3 H( - ed.m_Volume * U * diag_P * Vt * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
#endif
}

void TetrahedronElementSet3::H_NHC0_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __TO_PORT_TO_3D
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        ec.m_R = U*Vt; //TEMP: for debug viz only

        // Compute diag_P
        Vec3 vec_diag_P( ComputeDiagonalP_NHC0( vec_diag_F, m_LameMu, m_LameLambda, m_Params.m_ECIE_e_threshold ) );
        Mat3x3 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );

        // Compute PK force
        Mat3x3 H( - ed.m_Volume * U * diag_P * Vt * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
#endif
}

void TetrahedronElementSet3::H_NHC1_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __TO_PORT_TO_3D
    // Acc per-element forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        ec.m_R = U*Vt; //TEMP: for debug viz only

        // Compute diag_P
        Vec3 vec_diag_P( ComputeDiagonalP_NHC1( vec_diag_F,
                                                m_LameMu, m_LameLambda,
                                                m_Params.m_ECIE_e_threshold, m_Params.m_ECIE_k_factor * m_Params.m_YoungModulus ) );
        Mat3x3 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );

        // Compute PK force
        Mat3x3 H( - ed.m_Volume * U * diag_P * Vt * ed.m_InvDm.Transposed() );
        vec_f[ ed.m_vecNID[0] ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //h0 = -h1-h2;
        vec_f[ ed.m_vecNID[1] ] += mal::GColumn<0>(H);
        vec_f[ ed.m_vecNID[2] ] += mal::GColumn<1>(H);
    }
#endif
}

Real TetrahedronElementSet3::L_V( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element elastic energy
    Real acc_twice_V(0);
#ifdef __USE_OMP
#  pragma omp parallel for reduction(+:acc_twice_V)
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        Vec3 u0( cache.m_vecX[nid0] - m_vecRefPos[nid0] ); //\todo Consider precomputing vec_u to avoid repeating x_i-r_i
        Vec3 u1( cache.m_vecX[nid1] - m_vecRefPos[nid1] );
        Vec3 u2( cache.m_vecX[nid2] - m_vecRefPos[nid2] );
        Vec3 u3( cache.m_vecX[nid3] - m_vecRefPos[nid3] );
        // Ve = 1/2 * Ue^T * Ke * Ue
        acc_twice_V += mal::Dot( u0, led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 + led.m_K03 * u3 );
        acc_twice_V += mal::Dot( u1, led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 + led.m_K13 * u3 );
        acc_twice_V += mal::Dot( u2, led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 + led.m_K23 * u3 );
        acc_twice_V += mal::Dot( u3, led.m_K30 * u0 + led.m_K31 * u1 + led.m_K32 * u2 + led.m_K33 * u3 );
    }
    return Real(0.5) * acc_twice_V;
}
Real TetrahedronElementSet3::C_V( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element elastic energy
    Real acc_twice_V(0);
#ifdef __USE_OMP
#  pragma omp parallel for reduction(+:acc_twice_V)
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        // Rotate world frame X to reference frame (no need to translate)
        Mat3x3 Rt( ec.m_R.Transposed() );
        Vec3 u0( Rt * cache.m_vecX[nid0] - m_vecRefPos[nid0] );
        Vec3 u1( Rt * cache.m_vecX[nid1] - m_vecRefPos[nid1] );
        Vec3 u2( Rt * cache.m_vecX[nid2] - m_vecRefPos[nid2] );
        Vec3 u3( Rt * cache.m_vecX[nid3] - m_vecRefPos[nid3] );
        // Ve = 1/2 * Ue^T * Ke * Ue, reference-system invariant
        acc_twice_V += mal::Dot( u0, led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 + led.m_K03 * u3 );
        acc_twice_V += mal::Dot( u1, led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 + led.m_K13 * u3 );
        acc_twice_V += mal::Dot( u2, led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 + led.m_K23 * u3 );
        acc_twice_V += mal::Dot( u3, led.m_K30 * u0 + led.m_K31 * u1 + led.m_K32 * u2 + led.m_K33 * u3 );
    }
    return Real(0.5) * acc_twice_V;
}
Real TetrahedronElementSet3::C_LCM_V( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );
    Real acc_V(0);
#ifdef __USE_OMP
#  pragma omp parallel for reduction(+:acc_V)
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        // Compute S from R
        Mat3x3 S( ec.m_R.Transposed() * ec.m_F ); //F = R*S
        Mat3x3 S_minus_Id( S-Mat3x3::Identity() );
        acc_V += ed.m_Volume * ( m_LameMu * ( S_minus_Id ).NormSqF()  //|S-I|^2_Frobenius
                                 + Real(0.5)*m_LameLambda*mal::Sq( mal::Trace(S_minus_Id) ) );
    }
    return acc_V;
}
Real TetrahedronElementSet3::C_CCM_V( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __TO_PORT_TO_3D
    // Acc per-element elastic energy
    Real acc_V(0);
#ifdef __USE_OMP
#  pragma omp parallel for reduction(+:acc_V)
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute S from R
        Mat3x3 S( ec.m_R.Transposed() * ec.m_F ); //F = R*S
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
        acc_V += ed.m_Volume * ( m_LameMu * ( S-Mat3x3::Identity() ).NormSqF()  //|S-I|^2_Frobenius
                               + Real(0.5)*m_LameLambda*mal::Sq( det_S - Real(1) ) );
    }
    return acc_V;
#else
    return Real(0);
#endif
}
Real TetrahedronElementSet3::H_LCM_V( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element elastic energy
    Real acc_V(0);
#ifdef __USE_OMP
#  pragma omp parallel for reduction(+:acc_V)
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        // \todo May be faster to use F instead of computing SVD again
        acc_V += ed.m_Volume * ( m_LameMu * ( mal::Sq(vec_diag_F[0]-Real(1)) + mal::Sq(vec_diag_F[1]-Real(1)) + mal::Sq(vec_diag_F[2]-Real(1)) )
                               + Real(0.5)*m_LameLambda*mal::Sq( (vec_diag_F[0]-Real(1)) + (vec_diag_F[1]-Real(1)) + (vec_diag_F[2]-Real(1)) ) );
    }
    return acc_V;
}
Real TetrahedronElementSet3::H_CCM_V( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __TO_PORT_TO_3D
    // Acc per-element elastic energy
    Real acc_V(0);
#ifdef __USE_OMP
#  pragma omp parallel for reduction(+:acc_V)
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        // \todo May be faster to use F instead of computing SVD again
        acc_V += ed.m_Volume * ( m_LameMu * ( mal::Sq(vec_diag_F[0]-Real(1)) + mal::Sq(vec_diag_F[1]-Real(1)) )
                               + Real(0.5)*m_LameLambda*mal::Sq( ec.m_DetF - Real(1) ) );
    }
    return acc_V;
#else
    return Real(0);
#endif
}
Real TetrahedronElementSet3::H_NHC0_V( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __TO_PORT_TO_3D
    // Acc per-element elastic energy
    Real acc_V(0);
#ifdef __USE_OMP
#  pragma omp parallel for reduction(+:acc_V)
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        //
        Real ecie_e( m_Params.m_ECIE_e_threshold );
        Vec3 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
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
        Real delta0( mal::Min( vec_diag_F[0] - ecie_e, 0 ) );
        Real delta1( mal::Min( vec_diag_F[1] - ecie_e, 0 ) );
        acc_density_Ve += g0 * delta0
                          + g1 * delta1
                          + h01 * delta0 * delta1;
        acc_V += ed.m_Volume * acc_density_Ve;
    }
    return acc_V;
#else
    return Real(0);
#endif
}

Real TetrahedronElementSet3::H_NHC1_V( IForce3::ICache* p_cache ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

#ifdef __TO_PORT_TO_3D
    // Acc per-element elastic energy
    Real acc_V(0);
#ifdef __USE_OMP
#  pragma omp parallel for reduction(+:acc_V)
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );

        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        //
        Real ecie_e( m_Params.m_ECIE_e_threshold );
        Real ecie_k( m_Params.m_ECIE_k_factor * m_Params.m_YoungModulus );
        bool bIsDegenerate0( vec_diag_F[0] < ecie_e );
        bool bIsDegenerate1( vec_diag_F[1] < ecie_e );
        Vec3 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
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
        acc_V += ed.m_Volume * acc_density_Ve;
    }
    return acc_V;
#else
    return Real(0);
#endif
}

//----------------------------------------------------------------
// df_x()
//----------------------------------------------------------------
void TetrahedronElementSet3::L_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element elastic forces
//\todo Accumulation may fail if 2 threads write on the same vec_f[i]
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        Vec3 dx0( vec_dx[nid0] );
        Vec3 dx1( vec_dx[nid1] );
        Vec3 dx2( vec_dx[nid2] );
        Vec3 dx3( vec_dx[nid3] );
#ifdef __USE_OMP
        Vec3 df0( led.m_K00 * dx0 + led.m_K01 * dx1 + led.m_K02 * dx2 + led.m_K03 * dx3 );
        Vec3 df1( led.m_K10 * dx0 + led.m_K11 * dx1 + led.m_K12 * dx2 + led.m_K13 * dx3 );
        Vec3 df2( led.m_K20 * dx0 + led.m_K21 * dx1 + led.m_K22 * dx2 + led.m_K23 * dx3 );
        Vec3 df3( led.m_K30 * dx0 + led.m_K31 * dx1 + led.m_K32 * dx2 + led.m_K33 * dx3 );
#  pragma omp critical
        vec_df[ nid0 ] -= df0;
        vec_df[ nid1 ] -= df1;
        vec_df[ nid2 ] -= df2;
        vec_df[ nid3 ] -= df3;
#else
        vec_df[ nid0 ] -= led.m_K00 * dx0 + led.m_K01 * dx1 + led.m_K02 * dx2 + led.m_K03 * dx3;
        vec_df[ nid1 ] -= led.m_K10 * dx0 + led.m_K11 * dx1 + led.m_K12 * dx2 + led.m_K13 * dx3;
        vec_df[ nid2 ] -= led.m_K20 * dx0 + led.m_K21 * dx1 + led.m_K22 * dx2 + led.m_K23 * dx3;
        vec_df[ nid3 ] -= led.m_K30 * dx0 + led.m_K31 * dx1 + led.m_K32 * dx2 + led.m_K33 * dx3;
#endif
    }
}

void TetrahedronElementSet3::C_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element elastic forces
//\todo Accumulation may fail if 2 threads write on the same vec_f[i]
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const LinearED& led( m_vecLED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        // Rotate world frame X to reference frame (no need to translate)
        Mat3x3 Rt( ec.m_R.Transposed() );
        Vec3 ldx0( Rt * vec_dx[nid0] );
        Vec3 ldx1( Rt * vec_dx[nid1] );
        Vec3 ldx2( Rt * vec_dx[nid2] );
        Vec3 ldx3( Rt * vec_dx[nid3] );
#ifdef __USE_OMP
        Vec3 df0( ec.m_R * ( led.m_K00 * ldx0 + led.m_K01 * ldx1 + led.m_K02 * ldx2 + led.m_K03 * ldx3 ) );
        Vec3 df1( ec.m_R * ( led.m_K10 * ldx0 + led.m_K11 * ldx1 + led.m_K12 * ldx2 + led.m_K13 * ldx3 ) );
        Vec3 df2( ec.m_R * ( led.m_K20 * ldx0 + led.m_K21 * ldx1 + led.m_K22 * ldx2 + led.m_K23 * ldx3 ) );
        Vec3 df3( ec.m_R * ( led.m_K30 * ldx0 + led.m_K31 * ldx1 + led.m_K32 * ldx2 + led.m_K33 * ldx3 ) );
#  pragma omp critical
        vec_df[ nid0 ] -= df0;
        vec_df[ nid1 ] -= df1;
        vec_df[ nid2 ] -= df2;
        vec_df[ nid3 ] -= df3;
#else
        // Compute elastic forces and rotate to world frame
        vec_df[ nid0 ] -= ec.m_R * ( led.m_K00 * ldx0 + led.m_K01 * ldx1 + led.m_K02 * ldx2 + led.m_K03 * ldx3 );
        vec_df[ nid1 ] -= ec.m_R * ( led.m_K10 * ldx0 + led.m_K11 * ldx1 + led.m_K12 * ldx2 + led.m_K13 * ldx3 );
        vec_df[ nid2 ] -= ec.m_R * ( led.m_K20 * ldx0 + led.m_K21 * ldx1 + led.m_K22 * ldx2 + led.m_K23 * ldx3 );
        vec_df[ nid3 ] -= ec.m_R * ( led.m_K30 * ldx0 + led.m_K31 * ldx1 + led.m_K32 * ldx2 + led.m_K33 * ldx3 );
#endif

        // The next part is not ready for OMP, yet...
#ifdef __USE_OMP
#  pragma omp critical
#endif
        {
            // Compute df terms that depend on dR
            Mat3x3 dR;
            if( Compute_dR_From_RM(m_Params,cache,ed,ec,Rt,vec_dx,dR) )
            {
                // Add dR force terms
                // Term dR_e * K_e * (R_eP^T x_e - r_e) = dR_e * K_e * u_e
                Vec3 u0( Rt * cache.m_vecX[nid0] - m_vecRefPos[nid0] );
                Vec3 u1( Rt * cache.m_vecX[nid1] - m_vecRefPos[nid1] );
                Vec3 u2( Rt * cache.m_vecX[nid2] - m_vecRefPos[nid2] );
                Vec3 u3( Rt * cache.m_vecX[nid3] - m_vecRefPos[nid3] );
                vec_df[ nid0 ] -= dR * ( led.m_K00 * u0 + led.m_K01 * u1 + led.m_K02 * u2 + led.m_K03 * u3 );
                vec_df[ nid1 ] -= dR * ( led.m_K10 * u0 + led.m_K11 * u1 + led.m_K12 * u2 + led.m_K13 * u3 );
                vec_df[ nid2 ] -= dR * ( led.m_K20 * u0 + led.m_K21 * u1 + led.m_K22 * u2 + led.m_K23 * u3 );
                vec_df[ nid3 ] -= dR * ( led.m_K30 * u0 + led.m_K31 * u1 + led.m_K32 * u2 + led.m_K33 * u3 );
                // Term R_e * K_e * dR_e^T * x_e
                Mat3x3 dRt( dR.Transposed() );
                Vec3 dRt_times_x0( dRt * cache.m_vecX[nid0] );
                Vec3 dRt_times_x1( dRt * cache.m_vecX[nid1] );
                Vec3 dRt_times_x2( dRt * cache.m_vecX[nid2] );
                Vec3 dRt_times_x3( dRt * cache.m_vecX[nid3] );
                vec_df[ nid0 ] -= ec.m_R * ( led.m_K00 * dRt_times_x0 + led.m_K01 * dRt_times_x1 + led.m_K02 * dRt_times_x2 + led.m_K03 * dRt_times_x3 );
                vec_df[ nid1 ] -= ec.m_R * ( led.m_K10 * dRt_times_x0 + led.m_K11 * dRt_times_x1 + led.m_K12 * dRt_times_x2 + led.m_K13 * dRt_times_x3 );
                vec_df[ nid2 ] -= ec.m_R * ( led.m_K20 * dRt_times_x0 + led.m_K21 * dRt_times_x1 + led.m_K22 * dRt_times_x2 + led.m_K23 * dRt_times_x3 );
                vec_df[ nid3 ] -= ec.m_R * ( led.m_K30 * dRt_times_x0 + led.m_K31 * dRt_times_x1 + led.m_K32 * dRt_times_x2 + led.m_K33 * dRt_times_x3 );
            }
        }
    }
}

void TetrahedronElementSet3::C_LCM_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        Mat3x3 Rt( mal::Transposed( ec.m_R ) );
        Mat3x3 S( Rt * ec.m_F );
        Mat3x3 dDs;
        TetrahedronElement3::Compute_D( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid3], dDs );
        Mat3x3 dF( dDs * ed.m_InvDm );
        // Compute dP terms independent of dR, from McAdams2011
        Mat3x3 dP( 2*m_LameMu * dF
                   + m_LameLambda*mal::Trace(Rt*dF) * ec.m_R );

        // Compute df terms that depend on dR, from McAdams2011
        Mat3x3 dR;
        if( Compute_dR_From_RM(m_Params,cache,ed,ec,Rt,vec_dx,dR) )
            dP += ( m_LameLambda*mal::Trace(S-Mat3x3::Identity()) - 2*m_LameMu ) * dR;

        // Compute PK force
        Mat3x3 dH( - ed.m_Volume * dP * ed.m_InvDm.Transposed() );
        vec_df[ nid0 ] += - mal::GColumn<0>(dH) - mal::GColumn<1>(dH) - mal::GColumn<2>(dH); //dh0 = -dh1-dh2-dh3;
        vec_df[ nid1 ] += mal::GColumn<0>(dH);
        vec_df[ nid2 ] += mal::GColumn<1>(dH);
        vec_df[ nid3 ] += mal::GColumn<2>(dH);
    }
}

void TetrahedronElementSet3::C_LCMH_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        Mat3x3 Rt( mal::Transposed( ec.m_R ) );
        Mat3x3 S( Rt * ec.m_F );
        Mat3x3 dDs;
        TetrahedronElement3::Compute_D( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid3], dDs );
        Mat3x3 dF( dDs * ed.m_InvDm );

        //TEMPORAL Using param m_ECIE_k_factor, consider removing it or defining specific param here
        Mat3x3 dP;
        if( m_Params.m_DM == Params::eDM_Truncated
            &&
            m_Params.m_ECIE_k_factor > 0.1 ) //Hybrid LCM
        {
            const Real cDimension = 3;
            const Real cRadiusWRP = m_Params.m_ECIE_k_factor;
            Real weight_LCM = mal::Clamp01( mal::Rcp(cRadiusWRP) * mal::Trace( mal::Abs( S - Mat3x3::Identity() ) )/cDimension );
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

        // Compute df terms that depend on dR, from McAdams2011
        Mat3x3 dR;
        if( Compute_dR_From_RM(m_Params,cache,ed,ec,Rt,vec_dx,dR) )
            dP += ( m_LameLambda*mal::Trace(S-Mat3x3::Identity()) - 2*m_LameMu ) * dR;

        // Compute PK force
        Mat3x3 dH( - ed.m_Volume * dP * ed.m_InvDm.Transposed() );
        vec_df[ nid0 ] += - mal::GColumn<0>(dH) - mal::GColumn<1>(dH) - mal::GColumn<2>(dH); //dh0 = -dh1-dh2-dh3;
        vec_df[ nid1 ] += mal::GColumn<0>(dH);
        vec_df[ nid2 ] += mal::GColumn<1>(dH);
        vec_df[ nid3 ] += mal::GColumn<2>(dH);
    }
}

void TetrahedronElementSet3::C_CCM_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const
{
    MS_ASSERT( 0 != p_cache );
    Cache3& cache( *static_cast<Cache3*>( p_cache ) );

    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        node_index_type nid3( ed.m_vecNID[3] );
        Mat3x3 Rt( mal::Transposed( ec.m_R ) ); //\note This R comes from SVD1 U*Vt
        Mat3x3 dDs;
        TetrahedronElement3::Compute_D( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], vec_dx[nid2], dDs );
        Mat3x3 dF( dDs * ed.m_InvDm );
        // Compute dP terms independent of dR
        Mat3x3 dP( 2*m_LameMu * dF ); //Robust dP term
        // Compute invF-dependent dP terms \todo There MUST BE a better regularization, see testFE, SPLIT seems the most reasonable approach
        Mat3x3 invF( mal::Inverse( ec.m_F, ec.m_DetF + mal::Sign(ec.m_DetF)*REGULARIZE_CCM_PD_U_THRESHOLD ) ); //\note Regularized for det(F) = 0
        Mat3x3 invFt( mal::Transposed( invF ) );
        dP += m_LameLambda * ( (2*mal::Sq(ec.m_DetF) - ec.m_DetF) * mal::Trace(invF*dF) * invFt
                               - (mal::Sq(ec.m_DetF) - ec.m_DetF) * invFt * mal::Transposed(dF) * invFt );

        // Compute df terms that depend on dR
        Mat3x3 dR;
        if( Compute_dR_From_RM(m_Params,cache,ed,ec,Rt,vec_dx,dR) )
            dP += - 2*m_LameMu * dR;

        // Compute PK force
        Mat3x3 dH( - ed.m_Volume * dP * ed.m_InvDm.Transposed() );
        vec_df[ nid0 ] += - mal::GColumn<0>(dH) - mal::GColumn<1>(dH) - mal::GColumn<2>(dH); //dh0 = -dh1-dh2-dh3;
        vec_df[ nid1 ] += mal::GColumn<0>(dH);
        vec_df[ nid2 ] += mal::GColumn<1>(dH);
        vec_df[ nid3 ] += mal::GColumn<2>(dH);
    }
}

/* \todo Use SVD differentials */
void TetrahedronElementSet3::H_LCM_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const
{
    //\todo By now use corotational version, SVD df is more complex
    return C_LCM_df_x( p_cache, vec_dx, vec_df );

#ifdef __STILL_DEVELOPING
    // Acc per-element elastic forces
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ElementData& ed( m_vecED[it_e] );
        const ElementCache3& ec( cache.m_vecEC[it_e] );
        node_index_type nid0( ed.m_vecNID[0] );
        node_index_type nid1( ed.m_vecNID[1] );
        node_index_type nid2( ed.m_vecNID[2] );
        Mat3x3 dDs;
        TetrahedronElement2::Compute_D( vec_dx[nid0], vec_dx[nid1], vec_dx[nid2], dDs );
        Mat3x3 dF( dDs * ed.m_InvDm );

        //--- Compute dP
        //-- (a) Directly, using dP(F;dF) matrix formula \todo BUT THIS DOES NOT FIX DEGENERATE F!!
        //\todo For comparison...

        //-- (b) From diagonalization, using dP(...) = U * delta_diag_P( d_diag_F ) * V^t
        //\sa "Robust Quasistatic Finite Elements and Flesh Simulation" and "Energetically Consistent Invertible Elasticity"
        // Compute SVD \todo Constant during CG iter, consider doing it in BeginEvaluation() and caching
        Mat3x3 U, Vt;
        Vec3 vec_diag_F;
        ComputeSVD_Invertible( ec.m_F, U, vec_diag_F, Vt );
        /* Compute diag_P
        Vec3 vec_diag_P( ComputeDiagonalP_LCM( vec_diag_F, m_LameMu, m_LameLambda ) );
        Mat3x3 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );
        */
        Mat3x3 d_diag_F( mal::Transposed(U) * dF * mal::Transposed(Vt) ); // U^T * \delta F * V

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
        Mat3x3 d_diag_P( ComputeDeltaDiagonalP_LCM( vec_diag_F, vec_delta_diag_F, m_LameMu, m_LameLambda ) );
        Mat3x3 dP( U * d_diag_P * Vt );

        // Compute PK force
        Mat3x3 dH( - ed.m_Volume * dP * ed.m_InvDm.Transposed() );
        vec_df[ nid0 ] += - mal::GColumn<0>(H) - mal::GColumn<1>(H); //dh0 = -dh1-dh2;
        vec_df[ nid1 ] += mal::GColumn<0>(H);
        vec_df[ nid2 ] += mal::GColumn<1>(H);
    }
#endif
}

/* \todo Use SVD differentials */
void TetrahedronElementSet3::H_CCM_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const
{
    //\todo By now use corotational version, SVD df is more complex
    return C_CCM_df_x( p_cache, vec_dx, vec_df );
}

//----------------------------------------------------------------
// EditableTetrahedronElementSet3
//----------------------------------------------------------------

EditableTetrahedronElementSet3::EditableTetrahedronElementSet3()
: m_LastElement(0)
, m_pAllocRefPos(0)
, m_pAllocED(0)
, m_pAllocLED(0)
{}

EditableTetrahedronElementSet3::~EditableTetrahedronElementSet3()
{
    ClearEditData();
}

void EditableTetrahedronElementSet3::ClearEditData()
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

//\todo if method is given per-element, params.m_MM == eMM_MaterialPerElement and alloc a per-element method descriptor
void EditableTetrahedronElementSet3::BeginEdition( unsigned int num_nodes, const Vec3 *vec_r, unsigned int num_elements, const Params &params )
{
    ClearEditData();
    m_NumNodes = num_nodes;
    m_NumElements = num_elements;
    //\todo here we could copy existing baked data as edit data to edit incrementally
    // Set params
    m_Params = params;
    // Alloc common stuff
    m_pAllocRefPos = new Vec3[num_nodes];
    memcpy( &m_pAllocRefPos[0], &vec_r[0], num_nodes * sizeof(Vec3) );
    m_pAllocED = new ElementData[num_elements];
    // Alloc specific stuff
    m_pAllocLED = new LinearED[num_elements];
    // Clear baked data
    ClearBakedData();
}

void EditableTetrahedronElementSet3::AddElement( unsigned int i, unsigned int j, unsigned int k, unsigned int l )
{
    MS_ASSERT( node_index_type(i) == i && node_index_type(j) == j && node_index_type(k) == k && node_index_type(l) == l );
    MS_ASSERT( m_LastElement < m_NumElements );
    ElementData& ed( m_pAllocED[m_LastElement++] );
    ed.m_vecNID[0] = node_index_type(i);
    ed.m_vecNID[1] = node_index_type(j);
    ed.m_vecNID[2] = node_index_type(k);
    ed.m_vecNID[3] = node_index_type(l);
}

bool EditableTetrahedronElementSet3::EndEdition()
{
    MS_ASSERT( m_LastElement == m_NumElements );
    LameParameters_From_YoungAndPoisson( m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                         m_LameMu, m_LameLambda );
    Real acc_total_volume(0);
#pragma omp parallel for reduction(+:acc_total_volume)
    for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
    {
        ElementData& ed( m_pAllocED[it_e] );
        Vec3 r0( m_pAllocRefPos[ ed.m_vecNID[0] ] );
        Vec3 r1( m_pAllocRefPos[ ed.m_vecNID[1] ] );
        Vec3 r2( m_pAllocRefPos[ ed.m_vecNID[2] ] );
        Vec3 r3( m_pAllocRefPos[ ed.m_vecNID[3] ] );
        // Common element init
        ed.m_Volume = TetrahedronElement3::Compute_Volume( r0, r1, r2, r3 );
        acc_total_volume += ed.m_Volume;
        // Specific element init
        // Linear
        Mat12x12 K;
        TetrahedronElement3::Compute_K( r0, r1, r2, r3,
                                        m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                        K );
        // K0X
        m_pAllocLED[it_e].m_K00 = mal::GRange<0,0,2,2>( K );
        m_pAllocLED[it_e].m_K01 = mal::GRange<0,3,2,5>( K );
        m_pAllocLED[it_e].m_K02 = mal::GRange<0,6,2,8>( K );
        m_pAllocLED[it_e].m_K03 = mal::GRange<0,9,2,11>( K );
        // K1X
        m_pAllocLED[it_e].m_K10 = mal::GRange<3,0,5,2>( K );
        m_pAllocLED[it_e].m_K11 = mal::GRange<3,3,5,5>( K );
        m_pAllocLED[it_e].m_K12 = mal::GRange<3,6,5,8>( K );
        m_pAllocLED[it_e].m_K13 = mal::GRange<3,9,5,11>( K );
        // K2X
        m_pAllocLED[it_e].m_K20 = mal::GRange<6,0,8,2>( K );
        m_pAllocLED[it_e].m_K21 = mal::GRange<6,3,8,5>( K );
        m_pAllocLED[it_e].m_K22 = mal::GRange<6,6,8,8>( K );
        m_pAllocLED[it_e].m_K23 = mal::GRange<6,9,8,11>( K );
        // K3X
        m_pAllocLED[it_e].m_K30 = mal::GRange<9,0,11,2>( K );
        m_pAllocLED[it_e].m_K31 = mal::GRange<9,3,11,5>( K );
        m_pAllocLED[it_e].m_K32 = mal::GRange<9,6,11,8>( K );
        m_pAllocLED[it_e].m_K33 = mal::GRange<9,9,11,11>( K );

        // C-H
        Mat3x3 Dm;
        TetrahedronElement3::Compute_D( r0, r1, r2, r3, Dm );
        ed.m_InvDm = mal::Inverse( Dm );
    }
    m_TotalVolume = acc_total_volume;
    // Init baked data
    TetrahedronElementSet3::SetBakedData( false, m_NumNodes, m_NumElements, m_pAllocRefPos, m_pAllocED, m_pAllocLED );
    return true;
}

void EditableTetrahedronElementSet3::SetParams( const Params &params )
{
    MS_ASSERT( IsValid() );
    m_Params = params;
    LameParameters_From_YoungAndPoisson( m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                         m_LameMu, m_LameLambda );
    // Rebuild material-dependent stuff directly
#ifdef __USE_OMP
#  pragma omp parallel for
#endif
    for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
    {
        ElementData& ed( m_pAllocED[it_e] );
        Vec3 r0( m_pAllocRefPos[ ed.m_vecNID[0] ] );
        Vec3 r1( m_pAllocRefPos[ ed.m_vecNID[1] ] );
        Vec3 r2( m_pAllocRefPos[ ed.m_vecNID[2] ] );
        Vec3 r3( m_pAllocRefPos[ ed.m_vecNID[3] ] );
        // Specific element init
        // Linear
        Mat12x12 K;
        TetrahedronElement3::Compute_K( r0, r1, r2, r3,
                                        m_Params.m_YoungModulus, m_Params.m_PoissonRatio,
                                        K );
        // K0X
        m_pAllocLED[it_e].m_K00 = mal::GRange<0,0,2,2>( K );
        m_pAllocLED[it_e].m_K01 = mal::GRange<0,3,2,5>( K );
        m_pAllocLED[it_e].m_K02 = mal::GRange<0,6,2,8>( K );
        m_pAllocLED[it_e].m_K03 = mal::GRange<0,9,2,11>( K );
        // K1X
        m_pAllocLED[it_e].m_K10 = mal::GRange<3,0,5,2>( K );
        m_pAllocLED[it_e].m_K11 = mal::GRange<3,3,5,5>( K );
        m_pAllocLED[it_e].m_K12 = mal::GRange<3,6,5,8>( K );
        m_pAllocLED[it_e].m_K13 = mal::GRange<3,9,5,11>( K );
        // K2X
        m_pAllocLED[it_e].m_K20 = mal::GRange<6,0,8,2>( K );
        m_pAllocLED[it_e].m_K21 = mal::GRange<6,3,8,5>( K );
        m_pAllocLED[it_e].m_K22 = mal::GRange<6,6,8,8>( K );
        m_pAllocLED[it_e].m_K23 = mal::GRange<6,9,8,11>( K );
        // K3X
        m_pAllocLED[it_e].m_K30 = mal::GRange<9,0,11,2>( K );
        m_pAllocLED[it_e].m_K31 = mal::GRange<9,3,11,5>( K );
        m_pAllocLED[it_e].m_K32 = mal::GRange<9,6,11,8>( K );
        m_pAllocLED[it_e].m_K33 = mal::GRange<9,9,11,11>( K );
    }
}

}}} //namespace S2::ms::fem
