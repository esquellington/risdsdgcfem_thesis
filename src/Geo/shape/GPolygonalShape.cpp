#include "GPolygonalShape.h"
#include <Mal/GRandom.h> //For Distort()

#ifdef __GEO_ENABLE_IMPLICIT_RBF
#  include "MatrixD.h"
#  include <iostream> //TEMPORAL!!
#endif

namespace geo {

//-----------------------------------------------------------------------------
//---- PolygonalShape2 Implementation
//-----------------------------------------------------------------------------

/*TEMP:
void PolygonalShape2::ComputeBV( bv::IBoundingVolume &bv,
                                 const GPolygonalShape<2>::transform_type &transform,
                                 const GPolygonalShape<2>::sdof_type *vec_sdof ) const
{
    const GPolygonalShape<2>::sdof_type *actual_sdof( ( 0 != vec_sdof ) ? vec_sdof : m_vecPointsVDL );
    switch( bv.GetType() )
    {
    case bv::eBV_Sphere2:
        GEO_ASSERT( false ); //Not yet implemented
        bv.As<bv::Sphere2>().SetPosRadius( transform.m_Pos, m_Radius );
        break;
    case bv::eBV_AABB2:
        {
            bv::AABB2 aabb( transform*actual_sdof[0], Vec2(m_Radius) );
            if( m_Radius > 0 )
                for( unsigned int i=1; i < m_NumVertices; i++ )
                    aabb.Merge( transform*actual_sdof[i], m_Radius );
            else
                for( unsigned int i=1; i < m_NumVertices; i++ )
                    aabb.Merge( transform*actual_sdof[i] );
            bv.As<bv::AABB2>() = aabb;
        }
        break;
    case bv::eBV_LSS2:
        GEO_ASSERT( false ); //Not yet implemented
        bv.As<bv::LSS2>().SetPosRadius( transform.m_Pos, transform.m_Pos, m_Radius );
        break;
    case bv::eBV_Void: break;
    case bv::eBV_Infinite: break;
    default:
        GEO_ASSERT( false ); //wrong type or dimension
        break;
    }
}
*/

#ifdef __GEO_ENABLE_IMPLICIT_RBF

#define __ENABLE_REGULARIZATION
//#define __ENABLE_NORMALS_OFFSET_DOUBLE  //two-sided offset
//#define __ENABLE_NORMALS_OFFSET_SINGLE  //single-sided offset
#define __ENABLE_NORMALS_OFFSET_CONFIGURABLE //\todo THIS OPTION SUBSUMES ALL THE OTHER, SIMPLIFY CODE ASAP!

class ImplicitRBF2: public IImplicit2
{
public:
    enum ENormalFlags {
        eNF_None    = 0,
        eNF_Inside  = 1,
        eNF_Outside = 2,
        eNF_Both    = (eNF_Inside | eNF_Outside)
    };

public:
    ImplicitRBF2( const PolygonalShape2 *p_ps2,
                  Real regularization_coeff = 0.0f,
                  Flags32 nf = eNF_Both,
                  Real normal_offset = 0.1f,
                  Real normal_tension = 0.1f,
                  int32 max_iter = 10,
                  Real epsilon = 0 )
#ifdef __ENABLE_NORMALS_OFFSET_DOUBLE
    : m_K(3*p_ps2->GetNumVertices()+1)
#elif defined __ENABLE_NORMALS_OFFSET_SINGLE
    : m_K(2*p_ps2->GetNumVertices()+1)
#elif defined __ENABLE_NORMALS_OFFSET_CONFIGURABLE
    : m_K(0) //\todo size depends on flags!!
#else
    : m_K(p_ps2->GetNumVertices()+1)
#endif
    , m_vecCoeffsRBF(m_K+3) //\todo depends on m_K, we resize it later if required
    , m_vecRBF(0)
        {
#ifdef __ENABLE_NORMALS_OFFSET_DOUBLE
            GEO_ASSERT( p_ps2->HasVDL(eVDL_Normals) );
            m_vecRBF = new Vec2[m_K];
            Vec2 *p_rbf(m_vecRBF);
            for( int i=0; i<p_ps2->GetNumVertices(); i++ )
            {
                *p_rbf++ = p_ps2->GetPoint(i);
                *p_rbf++ = p_ps2->GetPoint(i) - normal_offset*p_ps2->GetNormal(i);
                *p_rbf++ = p_ps2->GetPoint(i) + normal_offset*p_ps2->GetNormal(i);
            }
#elif defined(__ENABLE_NORMALS_OFFSET_SINGLE)
            GEO_ASSERT( p_ps2->HasVDL(eVDL_Normals) );
            m_vecRBF = new Vec2[m_K];
            Vec2 *p_rbf(m_vecRBF);
            for( int i=0; i<p_ps2->GetNumVertices(); i++ )
            {
                *p_rbf++ = p_ps2->GetPoint(i);
                *p_rbf++ = p_ps2->GetPoint(i) - normal_offset*p_ps2->GetNormal(i);
            }
#elif defined(__ENABLE_NORMALS_OFFSET_CONFIGURABLE)
            if( nf == eNF_Both )
            {
                GEO_ASSERT( p_ps2->HasVDL(eVDL_Normals) );
                m_K = 3*p_ps2->GetNumVertices()+1;
                m_vecRBF = new Vec2[m_K];
                Vec2 *p_rbf(m_vecRBF);
                for( int i=0; i<p_ps2->GetNumVertices(); i++ )
                {
                    *p_rbf++ = p_ps2->GetPoint(i);
                    *p_rbf++ = p_ps2->GetPoint(i) - normal_offset*p_ps2->GetNormal(i);
                    *p_rbf++ = p_ps2->GetPoint(i) + normal_offset*p_ps2->GetNormal(i);
                }
            }
            else if( nf.Test( eNF_Inside ) )
            {
                GEO_ASSERT( p_ps2->HasVDL(eVDL_Normals) );
                m_K = 2*p_ps2->GetNumVertices()+1;
                m_vecRBF = new Vec2[m_K];
                Vec2 *p_rbf(m_vecRBF);
                for( int i=0; i<p_ps2->GetNumVertices(); i++ )
                {
                    *p_rbf++ = p_ps2->GetPoint(i);
                    *p_rbf++ = p_ps2->GetPoint(i) - normal_offset*p_ps2->GetNormal(i);
                }
            }
            else if( nf.Test( eNF_Outside ) )
            {
                GEO_ASSERT( p_ps2->HasVDL(eVDL_Normals) );
                m_K = 2*p_ps2->GetNumVertices()+1;
                m_vecRBF = new Vec2[m_K];
                Vec2 *p_rbf(m_vecRBF);
                for( int i=0; i<p_ps2->GetNumVertices(); i++ )
                {
                    *p_rbf++ = p_ps2->GetPoint(i);
                    *p_rbf++ = p_ps2->GetPoint(i) + normal_offset*p_ps2->GetNormal(i);
                }
            }
            else
            {
                m_K = p_ps2->GetNumVertices()+1;
                m_vecRBF = new Vec2[m_K];
                Vec2 *p_rbf(m_vecRBF);
                for( int i=0; i<p_ps2->GetNumVertices(); i++ )
                    *p_rbf++ = p_ps2->GetPoint(i);
            }
            m_vecCoeffsRBF.Resize(m_K+3); //\todo Must resize it becase m_K has just been computed...
#else
            m_vecRBF = new Vec2[m_K];
            for( int i=0; i<p_ps2->GetNumVertices(); i++ ) m_vecRBF[i] = p_ps2->GetPoint(i);
#endif
            // TEMPORAL + constraint at barycenter
            m_vecRBF[m_K-1] = p_ps2->GetBarycenter();
            int rows( m_K + 3 ); //K Points + Linear(1 + D dimension)
            int cols( rows );
            S2::LA::MatrixD A(rows,cols);
            S2::LA::VectorD b(rows);
            //--- Fill A and b
            // RBF coeffs, LT + D
            for( int i=0; i<m_K; i++ )
            {
                for( int j=0; j<i; j++ )
                {
                    Real dsq = (m_vecRBF[i]-m_vecRBF[j]).NormSq();
                    A(i,j) = dsq*mal::Log(dsq)*Real(0.5); //d^2*log(d) == d^2*log(d^2)*0.5b
                    //A(i,j) = dsq*dsq; //alternative kernel
                    //A(i,j) = mal::Rcp(dsq*dsq+0.01); //alternative kernel
                }
                A(i,i) = 0; //Otherwise Log(0) = -infinity!!
                b[i] = 0; //value h_i = F(p_i)
            }
            // TEMPORAL +1 constraint at barycenter
            b[m_K-1] = 1;

#ifdef __ENABLE_NORMALS_OFFSET_DOUBLE
            // Int/Ext normal offsets have non-zero values
            for( int i=0; i<p_ps2->GetNumVertices(); i++ )
            {
                b[3*i+1] =  normal_tension;
                b[3*i+2] = -normal_tension;
            }
#elif defined(__ENABLE_NORMALS_OFFSET_SINGLE)
            // Int normal offsets have non-zero values
            for( int i=0; i<p_ps2->GetNumVertices(); i++ )
                b[2*i+1] =  normal_tension;
#elif defined(__ENABLE_NORMALS_OFFSET_CONFIGURABLE)
            if( nf == eNF_Both )
            {
                for( int i=0; i<p_ps2->GetNumVertices(); i++ )
                {
                    b[3*i+1] =  normal_tension;
                    b[3*i+2] = -normal_tension;
                }
            }
            else if( nf.Test( eNF_Inside ) )
            {
                for( int i=0; i<p_ps2->GetNumVertices(); i++ )
                    b[2*i+1] = normal_tension;
            }
            else if( nf.Test( eNF_Outside ) )
            {
                for( int i=0; i<p_ps2->GetNumVertices(); i++ )
                    b[2*i+1] = -normal_tension;
            }
#endif

            // Linear part, LT + D
            for( int j=0; j < m_K; j++ ) A(m_K,j) = 1;
            for( int j=0; j < m_K; j++ ) A(m_K+1,j) = m_vecRBF[j].x();
            for( int j=0; j < m_K; j++ ) A(m_K+2,j) = m_vecRBF[j].y();
            for( int i=m_K; i<rows; i++ ) b[i] = 0; //Linear part values
            // Zeros part, LT + D
            for( int i=m_K; i<rows; i++ )
                for( int j=m_K; j<=i; j++ )
                    A(i,j) = 0;
            // Symmetrize LT to UT
            for( int i=0; i<rows; i++ )
                for( int j=0; j<i; j++ )
                    A(j,i) = A(i,j);

#ifdef __ENABLE_REGULARIZATION
            // \todo Regularize means this? see regtutorial.pdf or HermiteRBF.pdf
            for( int i=0; i<rows; i++ ) if( A(i,i) == 0 ) A(i,i) = regularization_coeff;
#endif

//#define __ENABLE_TRACE_SOLVE_SLE
#ifdef __ENABLE_TRACE_SOLVE_SLE
            //DEBUG
            std::cout << "A = " << std::endl;
            for( int i=0; i<rows; i++ )
            {
                std::cout << "[ ";
                for( int j=0; j<cols; j++ ) std::cout << A(i,j) << " ";
                std::cout << "]" << std::endl;
            }
            std::cout << "b = ";
            std::cout << "[ ";
            for( int i=0; i<rows; i++ ) std::cout << b[i] << " ";
            std::cout << "]" << std::endl;
#endif

#define SLE_SOLVER_LU 0
#define SLE_SOLVER_CG 1
#define SLE_SOLVER_GS 2
#define SLE_SOLVER_JC 3
#define SLE_SOLVER SLE_SOLVER_LU

#if SLE_SOLVER == SLE_SOLVER_LU
            // Solve with LU
            SolveLU(A,m_vecCoeffsRBF,b);
            Real eps( S2::LA::VectorD(A*m_vecCoeffsRBF - b).Norm2() );
            if( eps > epsilon )
                std::cout << "WARNING!! SolveLU Epsilon = " << eps << std::endl;
#elif SLE_SOLVER == SLE_SOLVER_CG
            // Solve with dense CG
            Real eps(epsilon);
            int iter(max_iter);
            //\todo m_vecCoeffsRBF could be an approximate solution?
            SolveCG(A,m_vecCoeffsRBF,b,eps,iter);
            Real eps2( S2::LA::VectorD(A*m_vecCoeffsRBF - b).Norm2() );
            if( eps > epsilon || eps2 > epsilon )
                std::cout << "WARNING!! SolveCG Epsilon = " << eps << ", " << eps2 << ", Iter = " << iter << std::endl;
#elif SLE_SOLVER == SLE_SOLVER_GS
            // Solve with iterative GS
            Real eps(epsilon);
            int iter(max_iter);
            //\todo m_vecCoeffsRBF should be an approximate solution
            SolveGS(A,m_vecCoeffsRBF,b,eps,iter);
            eps = S2::LA::VectorD(A*m_vecCoeffsRBF - b).Norm2();
            if( eps > epsilon )
                std::cout << "WARNING!! SolveGS Epsilon = " << eps << ", Iter = " << iter << std::endl;
#elif SLE_SOLVER == SLE_SOLVER_JC
            // Solve with iterative GS
            Real eps(epsilon);
            int iter(max_iter);
            //\todo m_vecCoeffsRBF should be an approximate solution
            SolveJacobi(A,m_vecCoeffsRBF,b,eps,iter);
            eps = S2::LA::VectorD(A*m_vecCoeffsRBF - b).Norm2();
            if( eps > epsilon )
                std::cout << "WARNING!! SolveJacobi Epsilon = " << eps << ", Iter = " << iter << std::endl;
#endif

#ifdef __ENABLE_TRACE_SOLVE_SLE
            // Debug
            std::cout << "x = ";
            std::cout << "[ ";
            for( int i=0; i<rows; i++ ) std::cout << m_vecCoeffsRBF[i] << " ";
            std::cout << "]" << std::endl;
#endif

//#define __ENABLE_DISTORT_RBF
#ifdef __ENABLE_DISTORT_RBF
            for( int it_point=0; it_point<m_K; it_point++ )
                //m_vecRBF[it_point] += mal::RandomRadialVec<2,Real>( -0.025, 0.025 );
                //m_vecRBF[it_point] += Vec2(0.2,0.2); //! SEMBLA QUE fitting RBF és invariant a translacions!!
                m_vecRBF[it_point] += mal::Random<Real>(-0.01f,0.01f)*(m_vecRBF[it_point]-m_vecRBF[0]);
#endif

        }
    ~ImplicitRBF2()
        {
            if( m_vecRBF ) delete m_vecRBF;
        }
    Real F( const Vec2 &point ) const
        {
            Real f(0);
            for( int i=0; i < m_K; i++ )
            {
                Real dsq( (point-m_vecRBF[i]).NormSq() );
                f += (dsq>0) ? m_vecCoeffsRBF[i]*dsq*mal::Log(dsq)*Real(0.5) : 0; //Log(0) = -infinity!!
                //f += m_vecCoeffsRBF[i]*dsq*dsq; //alternative kernel
                //f += m_vecCoeffsRBF[i]*mal::Rcp(dsq*dsq+0.01); //alternative kernel
            }
            return (f + m_vecCoeffsRBF[m_K]
                    + m_vecCoeffsRBF[m_K+1]*point.x()
                    + m_vecCoeffsRBF[m_K+2]*point.y());
        }
    Vec2 dF_dp( const Vec2 &point ) const { return Vec2(0,0); }

private:
    int m_K;
    Vec2 *m_vecRBF;
    S2::LA::VectorD m_vecCoeffsRBF;
};

const IImplicit2 *PolygonalShape2::GetImplicit2() const
{
    if( 0 == m_pImplicit2 ) m_pImplicit2 = new ImplicitRBF2(this);//\note default params = 0.1f, ImplicitRBF2::eNF_Both, 0.1f, 0.1f
    return m_pImplicit2;
}

#endif //__GEO_ENABLE_IMPLICIT_RBF

//-----------------------------------------------------------------------------
//---- EditablePolygonalShape2 Implementation
//-----------------------------------------------------------------------------

EditablePolygonalShape2::EditablePolygonalShape2()
: PolygonalShape2()
, m_MaxVertices(0)
, m_ExplicitVDL(0)
, m_bIsEditing(false)
#ifdef __GEO_ENABLE_IMPLICIT_RBF
// Tweaks
, m_RegularizationCoeff(0.1f)
, m_NormalFlags(ImplicitRBF2::eNF_Both)
, m_NormalOffset(0.1f)
, m_NormalTension(0.1f)
, m_MaxIter(10)
, m_Epsilon(0)
#endif
{}

EditablePolygonalShape2::~EditablePolygonalShape2()
{
    Dealloc();
}

void EditablePolygonalShape2::Dealloc()
{
    GEO_ASSERT( !m_bIsEditing );
    m_MaxVertices = 0;
    m_NumVertices = 0;
    m_ExplicitVDL.Set(0);
    m_AvailableVDL.Set(0);
    if( 0 != m_vecPointsVDL ) { delete [] m_vecPointsVDL; m_vecPointsVDL = 0; }
    if( 0 != m_vecNormalsVDL ) { delete [] m_vecNormalsVDL; m_vecNormalsVDL = 0; }
    if( 0 != m_vecTangentsVDL ) { delete [] m_vecTangentsVDL; m_vecTangentsVDL = 0; }
    if( 0 != m_vecLambdasVDL ) { delete [] m_vecLambdasVDL; m_vecLambdasVDL = 0; }
    if( 0 != m_vecFlagsVDL ) { delete [] m_vecFlagsVDL; m_vecFlagsVDL = 0; }
}

void EditablePolygonalShape2::Clear()
{
    GEO_ASSERT( !m_bIsEditing );
    m_NumVertices = 0;
    m_bIsClosed = false;
#ifdef __GEO_ENABLE_IMPLICIT_RBF

    if( 0 != m_pImplicit2 )
    {
        delete m_pImplicit2;
        m_pImplicit2 = 0;
    }
#endif
}

void EditablePolygonalShape2::Alloc( unsigned int max_vertices,
                                     Flags32 explicit_vdl, Flags32 available_vdl )
{
    GEO_ASSERT( !m_bIsEditing );
    GEO_ASSERT( m_MaxVertices == 0 && m_NumVertices == 0 );
    GEO_ASSERT( max_vertices > 1 );
    GEO_ASSERT( explicit_vdl.Test( eVDL_Points ) );
    m_MaxVertices = max_vertices;
    m_ExplicitVDL = explicit_vdl;
    m_AvailableVDL = available_vdl;
    m_AvailableVDL.Enable( eVDL_Points ); //!< Points must exist, at least
    // Reserve requided VDL
    m_vecPointsVDL = new Vec2[max_vertices];
    if( m_ExplicitVDL.Test(eVDL_Normals) || m_AvailableVDL.Test(eVDL_Normals) ) m_vecNormalsVDL = new Vec2[max_vertices];
    if( m_ExplicitVDL.Test(eVDL_Tangents) || m_AvailableVDL.Test(eVDL_Tangents) ) m_vecTangentsVDL = new Vec2[max_vertices];
    if( m_ExplicitVDL.Test(eVDL_Lambdas) || m_AvailableVDL.Test(eVDL_Lambdas) ) m_vecLambdasVDL = new Real[max_vertices];
    m_vecFlagsVDL = new Flags32[max_vertices];
    // Reset VDL \todo OPTIMIZE using MemZero!!
    for( unsigned int it=0; it<m_MaxVertices; it++ )
    {
        m_vecPointsVDL[it] = Vec2(0,0);
        if( 0 != m_vecNormalsVDL ) m_vecNormalsVDL[it] = Vec2(0,0);
        if( 0 != m_vecTangentsVDL ) m_vecTangentsVDL[it] = Vec2(0,0);
        if( 0 != m_vecLambdasVDL ) m_vecLambdasVDL[it] = Real(0);
        m_vecFlagsVDL[it].Reset();
    }
}

void EditablePolygonalShape2::BeginEdition()
{
    GEO_ASSERT( m_MaxVertices > 0 && !m_bIsEditing );
    m_bIsEditing = true;
}

bool EditablePolygonalShape2::EndEdition()
{
    GEO_ASSERT( m_bIsEditing );
    GEO_ASSERT( (!m_bIsClosed && m_NumVertices > 1)
                || (m_bIsClosed && m_NumVertices > 2) );
    m_bIsEditing = false;
    if( Validate() )
    {
        RecomputeFlags();
        if( m_AvailableVDL.Test(eVDL_Lambdas) ) RecomputeLambdas();
        if( m_AvailableVDL.Test(eVDL_Tangents) ) RecomputeTangents();
        if( m_AvailableVDL.Test(eVDL_Normals) ) RecomputeNormals();
#ifdef __GEO_ENABLE_IMPLICIT_RBF
        // Recompute implicit if existed
        if( 0 != m_pImplicit2 )
        {
            delete m_pImplicit2;
            m_pImplicit2 = new ImplicitRBF2( this,
                                             m_RegularizationCoeff,
                                             m_NormalFlags,
                                             m_NormalOffset,
                                             m_NormalTension,
                                             m_MaxIter,
                                             m_Epsilon );
        }
#endif
        return true;
    }
    return false;
}

int EditablePolygonalShape2::AddPoint( const Vec2 &point )
{
    GEO_ASSERT( m_bIsEditing );
    GEO_ASSERT( m_ExplicitVDL.Test(eVDL_Points) );
    GEO_ASSERT( m_NumVertices < m_MaxVertices );
    m_vecPointsVDL[m_NumVertices++] = point;
    return m_NumVertices-1;
}

void EditablePolygonalShape2::SetPoint( int vid, const Vec2 &point )
{
    GEO_ASSERT( m_bIsEditing );
    GEO_ASSERT( m_ExplicitVDL.Test(eVDL_Points) && vid >= 0 && vid < (int)m_NumVertices );
    m_vecPointsVDL[vid] = point;
}
void EditablePolygonalShape2::SetNormal( int vid, const Vec2 &normal )
{
    GEO_ASSERT( m_bIsEditing );
    GEO_ASSERT( m_ExplicitVDL.Test(eVDL_Normals) && vid >= 0 && vid < (int)m_NumVertices );
    m_vecNormalsVDL[vid] = normal;
    m_vecFlagsVDL[vid].Enable(eVF_HasNormal);
}
void EditablePolygonalShape2::SetTangent( int vid, const Vec2 &tangent )
{
    GEO_ASSERT( m_bIsEditing );
    GEO_ASSERT( m_ExplicitVDL.Test(eVDL_Tangents) && vid >= 0 && vid < (int)m_NumVertices );
    m_vecTangentsVDL[vid] = tangent;
    m_vecFlagsVDL[vid].Enable(eVF_HasTangent);
}
void EditablePolygonalShape2::SetLambda( int vid, Real lambda )
{
    GEO_ASSERT( m_bIsEditing );
    GEO_ASSERT( m_ExplicitVDL.Test(eVDL_Lambdas) && vid >= 0 && vid < (int)m_NumVertices );
    //\todo This may be too restrictive, consider allowing arbitrary lambdas and scaling them to [0,1] aftwerards
    GEO_ASSERT( lambda >= Real(0) && lambda <= Real(1) );
    m_vecLambdasVDL[vid] = lambda;
    m_vecFlagsVDL[vid].Enable(eVF_HasLambda);
}

int EditablePolygonalShape2::Refine( Real lambda )
{
    GEO_ASSERT( lambda > Real(0) && lambda < Real(1) ); //Strict (0,1), cannot refine an existing vertex!
    GEO_ASSERT( m_bIsEditing );
    GEO_ASSERT( m_NumVertices < m_MaxVertices );
    GEO_ASSERT( m_AvailableVDL.Test(eVDL_Lambdas) && m_vecLambdasVDL != 0 );

    // Look for segment that contains lambda
    unsigned int vid(1);
    while( vid < m_NumVertices && lambda > m_vecLambdasVDL[vid] ) vid++;
    // Either found or not found but closed, thus, vid0 is m_NumVertices-1
    GEO_ASSERT( m_bIsClosed || vid < m_NumVertices );
    // If lambda is on an existing vertex, we don't refine it
    if( m_vecLambdasVDL[vid] == lambda ) return -1;
    // Shift all following vertices
    for( unsigned int i=m_NumVertices; i > vid; i-- )
    {
        m_vecPointsVDL[i] = m_vecPointsVDL[i-1];
        if( 0 != m_vecNormalsVDL ) m_vecNormalsVDL[i] = m_vecNormalsVDL[i-1];
        if( 0 != m_vecTangentsVDL ) m_vecTangentsVDL[i] = m_vecTangentsVDL[i-1];
        m_vecLambdasVDL[i] = m_vecLambdasVDL[i-1];
        m_vecFlagsVDL[i] = m_vecFlagsVDL[i-1];
        //\todo Copy other explicit data (available non-explicit data will be re-generated)
    }

    // Insert refined vertex at lambda in (lambda[vid0],lambda[vid1]),
    // wrapping vid1 if it's the IsClosed() and it's the last vtx
    unsigned int vid0(vid-1);
    unsigned int vid1(vid);
    Real local_lambda(0);
    if( vid < m_NumVertices )
        local_lambda = (lambda - m_vecLambdasVDL[vid0]) / (m_vecLambdasVDL[vid1]-m_vecLambdasVDL[vid0]);
    else
    {
        vid1 = vid1 % m_NumVertices;
        local_lambda = (lambda - m_vecLambdasVDL[vid0]) / ( Real(1) - m_vecLambdasVDL[vid0] ); //Lambda is 1 at vid=0 when wrapping
    }
    GEO_ASSERT( Real(0) < local_lambda && local_lambda < Real(1) );
    m_vecPointsVDL[vid] = (Real(1)-local_lambda)*m_vecPointsVDL[vid0] + local_lambda*m_vecPointsVDL[vid1];
    m_vecLambdasVDL[vid] = lambda;
    m_NumVertices++;
    return vid;
}

void EditablePolygonalShape2::SetClosed( bool b_closed )
{
    GEO_ASSERT( m_bIsEditing );
    if( m_NumVertices > 2 ) m_bIsClosed = b_closed;
    else
    {
        m_bIsClosed = false;
        if( b_closed ) GEO_LOG_WARNING( "SetClosed(false) on an EditablePolygonalShape2 with < 3 vertices!" );
    }
}

void EditablePolygonalShape2::SetRadius( Real radius )
{
    GEO_ASSERT( m_bIsEditing );
    m_Radius = radius;
}

/* Refine where needed in order to obtain an approximate uniform
   sampling
*/
void EditablePolygonalShape2::Uniformize()
{
    GEO_ASSERT( m_AvailableVDL.Test(eVDL_Lambdas) && m_vecLambdasVDL != 0 );

    // Look for smallest segment
    Real min_length(1); //max lambda
    for( unsigned int vid=1; vid < m_NumVertices; vid++ )
    {
        Real length( m_vecLambdasVDL[vid] - m_vecLambdasVDL[vid-1] );
        GEO_ASSERT( length > 0 );
        if( length < min_length ) min_length = length;
    }
    if( m_bIsClosed )
    {
        Real length( Real(1) - m_vecLambdasVDL[m_MaxVertices-1] );
        if( length < min_length ) min_length = length;
    }
    //GEO_LOG_WARNING( "MinLength=%f", min_length );

    BeginEdition();
    {
        // Re-sample whole polygonal to smallest segment
        // length. m_NumVertices may grow during this loop, no problem
        for( unsigned int vid=1; vid < m_NumVertices && m_NumVertices < m_MaxVertices; vid++ )
        {
            Real length( m_vecLambdasVDL[vid] - m_vecLambdasVDL[vid-1] );
            if( length > min_length )
            {
                //GEO_LOG_WARNING( "Refining vid=%d", vid );
                //refine as floor(times)
                int num_refines( int(mal::Floor( length / min_length )) - 1 );
                for( int i=0; i<num_refines; i++ )
                    Refine( m_vecLambdasVDL[vid-1] + (i+1)*length/(num_refines+1) );
            }
        }
        // Refine last segment if closed
        if( m_bIsClosed && m_NumVertices < m_MaxVertices )
        {
            int vid( m_NumVertices );
            Real length( Real(1) - m_vecLambdasVDL[vid-1] );
            if( length > min_length )
            {
                //GEO_LOG_WARNING( "Refining vid=%d", vid );
                //refine as floor(times)
                int num_refines( int(mal::Floor( length / min_length )) - 1 );
                for( int i=0; i<num_refines; i++ )
                    Refine( m_vecLambdasVDL[vid-1] + (i+1)*length/(num_refines+1) );
            }
        }
    }
    EndEdition();
}

/*\todo Support different sets of Explicit VDL!!
  copy/interpolate other explicitly specified VDL

  \todo If num_vertices * 2^num_levels <= max_vertices, there's no need to realloc!!
*/
void EditablePolygonalShape2::Subdivide( int num_levels )
{
    for( int i=0; i<num_levels; i++ )
    {
        BeginEdition();

        Vec2 *vec_points( m_vecPointsVDL );
        Vec2 *vec_normals( m_vecNormalsVDL );
        Vec2 *vec_tangents( m_vecTangentsVDL );
        Real *vec_lambdas( m_vecLambdasVDL );
        Flags32 *vec_flags( m_vecFlagsVDL );

        int new_num_vertices( (m_bIsClosed) ? 2*m_NumVertices : 2*m_NumVertices-1 );
        m_vecPointsVDL = new Vec2[new_num_vertices];
        // Copy
        for( unsigned int it_point=0; it_point<m_NumVertices; it_point++ )
            m_vecPointsVDL[2*it_point] = vec_points[it_point];
        // Interpolate
        for( unsigned int it_point=0; it_point<m_NumVertices-1; it_point++ )
            m_vecPointsVDL[2*it_point+1] = Real(0.5)*(vec_points[it_point]+vec_points[it_point+1]);
                                           //+ Real(0.1)*m_vecNormalsVDL[it_point];
        if( m_bIsClosed )
            m_vecPointsVDL[new_num_vertices-1] = Real(0.5)*(vec_points[m_NumVertices-1]+vec_points[0]);

        // Delete/Realloc stuff
        delete vec_points;
        if( vec_normals )
        {
            delete vec_normals;
            m_vecNormalsVDL = new Vec2[new_num_vertices];
        }
        if( vec_tangents )
        {
            delete vec_tangents;
            m_vecTangentsVDL = new Vec2[new_num_vertices];
        }
        if( vec_lambdas )
        {
            delete vec_lambdas;
            m_vecLambdasVDL = new Real[new_num_vertices];
        }
        if( vec_flags )
        {
            delete vec_flags;
            m_vecFlagsVDL = new Flags32[new_num_vertices];
        }

        m_NumVertices = new_num_vertices;

        // Recompute everything necessary...
        EndEdition();
    }
}

void EditablePolygonalShape2::Simplify( Real threshold )
{
}

//\todo Handle other explicit VDL apart from Points
void EditablePolygonalShape2::Distort( Real magnitude )
{
    BeginEdition();
    for( unsigned int it_point=0; it_point<m_NumVertices; it_point++ )
        m_vecPointsVDL[it_point] += mal::RandomRadialVec<2,Real>( -magnitude, magnitude );
    EndEdition();
}

void EditablePolygonalShape2::Distort( Real magnitude, Real frequency )
{
    BeginEdition();
    std::vector<Vec2> vec_points(m_NumVertices);
    for( unsigned int i=0; i<m_NumVertices; i++ ) vec_points[i] = m_vecPointsVDL[i];
    for( unsigned int it_point=0; it_point<m_NumVertices; it_point++ )
    {
        Vec2 p( vec_points[it_point] );
        Vec2 n( mal::Normalized( mal::PerpendicularCW( vec_points[(it_point+1)%m_NumVertices]
                                                       - vec_points[((int)it_point-1+m_NumVertices)%m_NumVertices] ) ) );
        Real d( magnitude * mal::Sin(frequency*p.x() ) * mal::Sin(frequency*p.y()) );
        m_vecPointsVDL[it_point] += d * n;
    }
    EndEdition();
}

bool EditablePolygonalShape2::Validate()
{
    // Validate:
    // - P_i+2 != P_i+1 != P_i (degenerate edge(s))
    // - Lambda_i+1 > Lambda_i
    // - \todo N_i == 1 if edited
    // - \todo T_i != 0 if edited
    //...
    return true;
}

void EditablePolygonalShape2::RecomputeNormals()
{
    if( m_ExplicitVDL.Test(eVDL_Normals) )
    {
        GEO_ASSERT(false);
    }
    else
    {
        // No explicit tangents, compute them from Points
        for( unsigned int it=0; it<m_NumVertices; it++ )
        {
            int it_prev( m_bIsClosed ? (it-1+m_NumVertices) % m_NumVertices : mal::Max<int>(0,it-1) );
            int it_next( m_bIsClosed ? (it+1) % m_NumVertices : mal::Min<int>(it+1,m_NumVertices-1) );
            Vec2 tangent( m_vecPointsVDL[it_next] - m_vecPointsVDL[it_prev]);
            m_vecNormalsVDL[it] = Vec2( tangent.y(), -tangent.x() ).Normalized();
        }
    }
}
void EditablePolygonalShape2::RecomputeTangents()
{
    if( m_ExplicitVDL.Test(eVDL_Tangents) )
    {
        GEO_ASSERT(false);
    }
    else
    {
        // No explicit tangents, compute them from Points
        for( unsigned int it=0; it<m_NumVertices; it++ )
        {
            int it_prev( m_bIsClosed ? (it-1+m_NumVertices) % m_NumVertices : mal::Max<int>(0,it-1) );
            int it_next( m_bIsClosed ? (it+1) % m_NumVertices : mal::Min<int>(it+1,m_NumVertices-1) );
            m_vecTangentsVDL[it] = Real(0.5f) * (m_vecPointsVDL[it_next] - m_vecPointsVDL[it_prev]);
        }
    }
}
void EditablePolygonalShape2::RecomputeLambdas()
{
    if( m_ExplicitVDL.Test(eVDL_Lambdas) )
    {
        GEO_ASSERT(false);
    }
    else
    {
        // No explicit lambdas, compute them from arc-length
        Real total_length(0);
        m_vecLambdasVDL[0] = 0;
        for( unsigned int it=1; it<m_NumVertices; it++ )
        {
            total_length += (m_vecPointsVDL[it]-m_vecPointsVDL[it-1]).Norm();
            m_vecLambdasVDL[it] = total_length;
        }
        if( m_bIsClosed ) total_length += (m_vecPointsVDL[0]-m_vecPointsVDL[m_NumVertices-1]).Norm();
        // Compute normalized lambdas using arc-lengths
        for( unsigned int it=0; it<m_NumVertices; it++ ) m_vecLambdasVDL[it] /= total_length;
    }
}
void EditablePolygonalShape2::RecomputeFlags()
{
    // Flags always exist, recompute the curvature part from existing points
    // Assume closed polygonal...
    for( unsigned int it=0; it<m_NumVertices; it++ )
    {
        int it_prev( (it - 1 + m_NumVertices) % m_NumVertices );
        int it_next( (it + 1) % m_NumVertices );
        Vec2 v0( m_vecPointsVDL[it] - m_vecPointsVDL[it_prev] );
        Vec2 perp0( v0.y(), -v0.x() );
        Vec2 v1( m_vecPointsVDL[it_next] - m_vecPointsVDL[it] );
        Real dot( mal::Dot(perp0,v1) );
       /*\todo This comparison does not represent cos(a) because
        vectors are not unitary!! Use a PROPER angular limit or,
        better, a coplanarity test as in krm::CEditableTriMesh */
        if( dot < -Real(0.001f) ) m_vecFlagsVDL[it].Enable(eVF_Convex);
        else if ( dot > Real(0.001f) ) m_vecFlagsVDL[it].Enable(eVF_Concave);
    }
    //...overwrite terminal vertex curvature with eVF_Planar if the polygonal is open
    if( !m_bIsClosed )
    {
        m_vecFlagsVDL[0].Disable(eVF_Convex | eVF_Concave);
        m_vecFlagsVDL[m_NumVertices-1].Disable(eVF_Convex | eVF_Concave);
    }
}

#ifdef __GEO_ENABLE_IMPLICIT_RBF
void EditablePolygonalShape2::GetImplicitRBF_Params( Real &reg_coeff, int32 &nf, Real &normal_offset, Real &normal_tension,
                                                     int32 &max_iter, Real &epsilon ) const
{
    reg_coeff = m_RegularizationCoeff;
    nf = m_NormalFlags;
    normal_offset = m_NormalOffset;
    normal_tension = m_NormalTension;
    max_iter = m_MaxIter;
    epsilon = m_Epsilon;
}

void EditablePolygonalShape2::SetImplicitRBF_Params( Real reg_coeff, int32 nf, Real normal_offset, Real normal_tension,
                                                     int32 max_iter, Real epsilon )
{
    m_RegularizationCoeff = reg_coeff;
    m_NormalFlags = nf;
    m_NormalOffset = normal_offset;
    m_NormalTension = normal_tension;
    m_MaxIter = max_iter;
    m_Epsilon = epsilon;
    // Recompute implicit if existed
    if( 0 != m_pImplicit2 )
    {
        //GEO_LOG( "Updating ImplicitRBF with reg=%f", m_RegularizationCoeff );
        delete m_pImplicit2;
        m_pImplicit2 = new ImplicitRBF2( this, m_RegularizationCoeff, m_NormalFlags, m_NormalOffset, m_NormalTension,
                                         m_MaxIter, m_Epsilon );
    }
}
#endif //#__GEO_ENABLE_IMPLICIT_RBF

//-----------------------------------------------------------------------------
//---- PolygonalShape3 Implementation
//-----------------------------------------------------------------------------

/*TEMP:
void PolygonalShape3::ComputeBV( bv::IBoundingVolume &bv,
                                 const GPolygonalShape<3>::transform_type &transform,
                                 const GPolygonalShape<3>::sdof_type *vec_sdof ) const
{
    const GPolygonalShape<3>::sdof_type *actual_sdof( ( 0 != vec_sdof ) ? vec_sdof : m_vecPointsVDL );
    switch( bv.GetType() )
    {
    case bv::eBV_Sphere3:
        GEO_ASSERT( false ); //Not yet implemented
        bv.As<bv::Sphere3>().SetPosRadius( transform.m_Pos, m_Radius );
        break;
    case bv::eBV_AABB3:
        {
            bv::AABB3 aabb( transform*actual_sdof[0], Vec3(m_Radius) );
            if( m_Radius > 0 )
                for( unsigned int i=1; i < m_NumVertices; i++ )
                    aabb.Merge( transform*actual_sdof[i], m_Radius );
            else
                for( unsigned int i=1; i < m_NumVertices; i++ )
                    aabb.Merge( transform*actual_sdof[i] );
            bv.As<bv::AABB3>() = aabb;
        }
        break;
    case bv::eBV_LSS3:
        GEO_ASSERT( false ); //Not yet implemented
        bv.As<bv::LSS3>().SetPosRadius( transform.m_Pos, transform.m_Pos, m_Radius );
        break;
    case bv::eBV_Void: break;
    case bv::eBV_Infinite: break;
    default:
        GEO_ASSERT( false ); //wrong type or dimension
        break;
    }
}
*/

} //namespace geo
