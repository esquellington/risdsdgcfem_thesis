#ifndef GEO_SHAPE_GPOLYGONALSHAPE_H
#define GEO_SHAPE_GPOLYGONALSHAPE_H

#include <Geo/shape/IShape.h>
#include <Geo/bv/bv.h> // BV types for ComputeBV()
#include <memory.h> //TEMPORAL req by memcpy()

//#define __GEO_ENABLE_IMPLICIT_RBF // Deprecated experimental functionality... works but is not practical

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/*\todo IMPORTANT: We keep PolygonalShape2 and PolygonalShape3
  SEPARATED because in 2d we have a draft implementation of
  __GEO_ENABLE_IMPLICIT_RBF that is NOT available in 3d yet... when
  Implicit is an "optional" and dimension-independent API, 2/3d
  polygonal may be merged as we did with Sphere, Capsule, etc...
*/

namespace geo {

enum EVertexDataLayer {
    // Standard VDL
    eVDL_None     = 0,
    eVDL_Points   = (1<<0), //Vec
    eVDL_Normals  = (1<<1), //Vec
    eVDL_Tangents = (1<<2), //Vec
    eVDL_Lambdas  = (1<<3), //Real param [0..1]
    eVDL_Flags    = (1<<4), //Flags32
    eVDL_All      = ( eVDL_Points | eVDL_Normals | eVDL_Tangents
                      | eVDL_Lambdas | eVDL_Flags),
    //\todo User-defined VDL
    eVDL_User0    = (1<<5),
    eVDL_User1    = (1<<6)
    //...
};

enum EVertexFlags {
    // Curvature
    eVF_Planar     = 0,
    eVF_Concave    = (1<<0),
    eVF_Convex     = (1<<1),
    eVF_Saddle     = (eVF_Convex | eVF_Concave), //only possible in 3d
    // Optional VDL
    eVF_HasNormal  = (1<<2),
    eVF_HasTangent = (1<<3),
    eVF_HasLambda  = (1<<4)
};

/*! 2/3d polygonal shape.
  - May be open or closed with counter-clockwise orientation.
  - Topology is implicit in the vertex order. No vertex is repeated if the polygonal is closed.
  - All per-vertex magnitudes are stored in local coords.
    - \todo Vertices may be re-transformed in different
      methods... this is shitty, maybe should be cached in IObject
  - TEMPORAL: The shape is NOT deformable, DOF are not used by now.
  - PolygonalShape2 can work as a 2D polygon if it's closed.

  \todo DOF are NOT used in this Polygonal shape by now... they will
  be required when extended to support deformation:
  - DefautSDOF will store the same as now, the "reference pose" for the Polygonal
  - Actual SDOF will store instance SDOF, which, at least should be
    the PointsVDL, but MAY include Normals and Tangents... and it's
    not clear if these and other data (lambdas, flags) should be
    recomputed in every single method that depends on instantaneous
    SDOF if up to date values are required... They could depend on the
    SDOF, but have no allocated memory where to be stored once and
    used many times...

  \todo Accept user-defined VDL, which must be described at creation with:
  - element type => stride
  - merge function: Any vertex fusion process will require "merging"
    the attributes of merged vertices, usually:
    - max/min
    - average
    - copy
  - Vertex duplication usually implies copy of vtx attributes, but
    could get it's own per-VDL user-defined function too.
*/
template <unsigned D>
class GPolygonalShape: public GIShapeD<D>
{
public:
    enum EConstants
    {
        cVersion = 1 //Initial version
    };
public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = D;
    typedef mal::GTransform<Real,D> transform_type;
    typedef mal::GVec<Real,D> sdof_type;
    typedef mal::GVec<Real,D> vec_type;
    //@}

    typedef GPolygonalShape<D> self_type;

public:
    GPolygonalShape( unsigned int count, bool b_closed, Real radius )
    : m_NumVertices(count), m_bIsClosed(b_closed), m_Radius(radius)
    , m_vecPointsVDL(0), m_vecNormalsVDL(0), m_vecTangentsVDL(0), m_vecLambdasVDL(0), m_vecFlagsVDL(0)
        {}
    virtual ~GPolygonalShape()
        {
            if( 0 != m_vecPointsVDL ) delete [] m_vecPointsVDL;
            if( 0 != m_vecNormalsVDL ) delete [] m_vecNormalsVDL;
            if( 0 != m_vecTangentsVDL ) delete [] m_vecTangentsVDL;
            if( 0 != m_vecLambdasVDL ) delete [] m_vecLambdasVDL;
            if( 0 != m_vecFlagsVDL ) delete [] m_vecFlagsVDL;
        }

    //!\name Shape params
    //@{
    void Init( unsigned int num_vertices, bool b_closed, Real radius,
               const vec_type *vec_points, const vec_type *vec_normals, const vec_type *vec_tangents,
               const Real *vec_lambdas, const Flags32 *vec_flags )
        {
            GEO_ASSERT( 0 == m_NumVertices && num_vertices > 1  && 0 != vec_points );
            m_NumVertices = num_vertices;
            m_bIsClosed = b_closed;
            m_Radius = radius;
            if( vec_points != 0 )
            {
                m_AvailableVDL.Enable( eVDL_Points );
                m_vecPointsVDL = new vec_type[m_NumVertices];
                memcpy( &m_vecPointsVDL[0], &vec_points[0], m_NumVertices * sizeof(vec_type) );
            }
            if( vec_normals != 0 )
            {
                m_AvailableVDL.Enable( eVDL_Normals );
                m_vecNormalsVDL = new vec_type[m_NumVertices];
                memcpy( &m_vecNormalsVDL[0], &vec_normals[0], m_NumVertices * sizeof(vec_type) );
            }
            if( vec_tangents != 0 )
            {
                m_AvailableVDL.Enable( eVDL_Tangents );
                m_vecTangentsVDL = new vec_type[m_NumVertices];
                memcpy( &m_vecTangentsVDL[0], &vec_tangents[0], m_NumVertices * sizeof(vec_type) );
            }
            if( vec_lambdas != 0 )
            {
                m_AvailableVDL.Enable( eVDL_Lambdas );
                m_vecLambdasVDL = new Real[m_NumVertices];
                memcpy( &m_vecLambdasVDL[0], &vec_lambdas[0], m_NumVertices * sizeof(Real) );
            }
            if( vec_flags != 0 )
            {
                m_AvailableVDL.Enable( eVDL_Flags );
                m_vecFlagsVDL = new Flags32[m_NumVertices];
                memcpy( &m_vecFlagsVDL[0], &vec_flags[0], m_NumVertices * sizeof(Flags32) );
            }
        }

    finline unsigned int GetNumVertices() const { return m_NumVertices; }
    finline bool IsClosed() const { return m_bIsClosed; }
    finline Real GetRadius() const { return m_Radius; }
    const vec_type *GetPointsVDL() const { return m_vecPointsVDL; }
    const vec_type *GetNormalsVDL() const { return m_vecNormalsVDL; }
    const vec_type *GetTangentsVDL() const { return m_vecTangentsVDL; }
    const Real *GetLambdasVDL() const { return m_vecLambdasVDL; }
    const Flags32 *GetFlagsVDL() const { return m_vecFlagsVDL; }
    //@}

    //!\name Dynamic IShape implementation
    //@{
    //EShapeType GetType() const { return ShapeTypeOf<self_type>(); }
    unsigned int GetNumSDOF() const { return m_NumVertices; }
    const sdof_type *GetVecDefaultSDOF() const { return m_vecPointsVDL; }
    const Real *GetVecDefaultDOF() const { return reinterpret_cast<const Real*>( m_vecPointsVDL ); }
    //TEMP: void ComputeBV( bv::IBoundingVolume &bv, const transform_type &transform, const sdof_type *vec_sdof ) const { }
    void ComputeBVD( bv::GBoundingVolumeD<D> &bv, const transform_type &transform, const sdof_type *vec_sdof ) const
        {
            const GPolygonalShape<D>::sdof_type *actual_sdof( ( 0 != vec_sdof ) ? vec_sdof : m_vecPointsVDL );
            switch( bv.GetType() )
            {
            case bv::eBV_Sphere2:
            case bv::eBV_Sphere3:
                GEO_ASSERT( false ); //Not yet implemented
                //bv.As<bv::Sphere2>().SetPosRadius( transform.m_Pos, m_Radius );
                break;
            case bv::eBV_AABB2:
            case bv::eBV_AABB3:
                {
                    bv::GAABB<D> aabb( transform*actual_sdof[0], vec_type(m_Radius) );
                    if( m_Radius > 0 )
                        for( unsigned int i=1; i < m_NumVertices; i++ )
                            aabb.Merge( transform*actual_sdof[i], m_Radius );
                    else
                        for( unsigned int i=1; i < m_NumVertices; i++ )
                            aabb.Merge( transform*actual_sdof[i] );
                    bv.template As< bv::GAABB<D> >() = aabb;
                }
                break;
            case bv::eBV_LSS2:
            case bv::eBV_LSS3:
                GEO_ASSERT( false ); //Not yet implemented
                //bv.As<bv::LSS2>().SetPosRadius( transform.m_Pos, transform.m_Pos, m_Radius );
                break;
            case bv::eBV_Void: break;
            case bv::eBV_Infinite: break;
            default:
                GEO_ASSERT( false ); //wrong type or dimension
                break;
            }
        }

    GIDomainSamplerD<D> *CreateDomainSampler() const { return 0; }
    //@}

    //\name Extended Polygonal interface
    //@{
    finline bool HasVDL( Flags32 flags_vdl ) const { return m_AvailableVDL.Includes(flags_vdl); }
    finline const vec_type &GetPoint( int vid ) const { GEO_ASSERT(HasVDL(eVDL_Points)); return m_vecPointsVDL[vid]; }

    finline vec_type V_Pos( uint32 vid, const sdof_type *vec_sdof ) const { GEO_ASSERT(HasVDL(eVDL_Points)); return vec_sdof[vid]; }
    finline vec_type V_Pos_0( uint32 vid ) const { GEO_ASSERT(HasVDL(eVDL_Points)); return m_vecPointsVDL[vid]; }

    finline vec_type GetPointAtLambda( Real lambda ) const { GEO_ASSERT(false); }
    //Optional per-vertex data
    finline vec_type GetNormal( int vid ) const { GEO_ASSERT(HasVDL(eVDL_Normals)); return m_vecNormalsVDL[vid]; }
    finline vec_type GetTangent( int vid ) const { GEO_ASSERT(HasVDL(eVDL_Tangents)); return m_vecTangentsVDL[vid]; }
    finline Real GetLambda( int vid ) const { GEO_ASSERT(HasVDL(eVDL_Lambdas)); return m_vecLambdasVDL[vid]; }
    finline Flags32 GetVertexFlags( int vid ) const { GEO_ASSERT(HasVDL(eVDL_Flags)); return Flags32(m_vecFlagsVDL[vid]); }
    vec_type GetBarycenter() const
        {
            // \todo Could be pre-computed and only updated when Points actually change...
            GEO_ASSERT( m_NumVertices > 0 );
            vec_type acc_points( vec_type::Zero() );
            for( unsigned int it_vtx=0; it_vtx<m_NumVertices; it_vtx++ ) acc_points += m_vecPointsVDL[it_vtx];
            return acc_points / (Real)m_NumVertices;
        }
    //@}

protected:
    unsigned int m_NumVertices;
    bool m_bIsClosed;
    Real m_Radius;

    Flags32 m_AvailableVDL;
    vec_type *m_vecPointsVDL;
    vec_type *m_vecNormalsVDL;
    vec_type *m_vecTangentsVDL;
    Real *m_vecLambdasVDL;
    Flags32 *m_vecFlagsVDL;
};

#ifdef __GEO_ENABLE_IMPLICIT_RBF
class IImplicit2
{
public:
    IImplicit2() {}
    virtual ~IImplicit2() {}
    virtual Real F( const Vec2 &point ) const = 0;
    virtual Vec2 dF_dp( const Vec2 &point ) const = 0;
    //virtual void Domain( Vec2 &min, Vec2 &max ); \todo Actually domain may be infinite... but interesting domain is not...
};
#endif //__GEO_ENABLE_IMPLICIT_RBF

//! 2D
class PolygonalShape2: public GPolygonalShape<2>
{
public:
    PolygonalShape2( unsigned int count = 0, bool b_closed = false, Real radius = 0 )
    : GPolygonalShape<2>(count,b_closed,radius)
#ifdef __GEO_ENABLE_IMPLICIT_RBF
    , m_pImplicit2(0)
#endif
        {}
    ~PolygonalShape2()
        {
#ifdef __GEO_ENABLE_IMPLICIT_RBF
            if(m_pImplicit2) delete m_pImplicit2;
#endif
        }

    EShapeType GetType() const { return eShape_Polygonal2; }
    /*TEMP:
    void ComputeBV( bv::IBoundingVolume &bv,
                    const GPolygonalShape<2>::transform_type &transform,
                    const GPolygonalShape<2>::sdof_type *vec_sdof ) const;
    */

    //\todo 2d only, work as polygon
    bool IsConvex() const { return false; }
    bool IsSelfIntersecting() const { return false; }
    Real ComputeDistance( const vec_type &point ) { return Real(0); }

#ifdef __GEO_ENABLE_IMPLICIT_RBF
    const IImplicit2 *GetImplicit2() const;
#endif

protected:
#ifdef __GEO_ENABLE_IMPLICIT_RBF
    mutable IImplicit2 *m_pImplicit2;
#endif
};

//! Editable 2D

/*! Allows edition of a polygonal from Points, as well as optional and
  per-vertex edition of Normals, Tangents and Lambdas.

  Output PolygonalShape2 will have the required available_vdl, which
  may be provided explicitly or internally computed from explicit_vdl
  data.

  \todo It should work the SAME WAY as EditableMeshSolidShape2, with
  Begin/End, temporary dynamic arrays and a SetBakedData() method to
  init it with shared or exclusive data. Now it's ALWAYS exclusive,
  all passed ptr are copied.
*/
class EditablePolygonalShape2: public PolygonalShape2
{
public:
    EditablePolygonalShape2();
    ~EditablePolygonalShape2();

    void Alloc( unsigned int max_vertices,
                Flags32 explicit_vdl = eVDL_Points,
                Flags32 available_vdl = eVDL_All );
    void Dealloc();
    inline bool IsAlloc() const { return m_MaxVertices > 0; }
    void Clear();

    //! \name Begin/End edition protocol
    //{@
    void BeginEdition(); //\todo NOT INCREMENTAL! should work as TetSolidShape3, etc...
    bool EndEdition();

    int AddPoint( const Vec2 &point );
    void SetPoint( int vid, const Vec2 &point );
    void SetNormal( int vid, const Vec2 &normal );
    void SetTangent( int vid, const Vec2 &tangent );
    void SetLambda( int vid, Real lambda );

    //!< Insert a vtx at specified lambda, returns new vid that can be SetXXX immediately, or -1 if error
    int Refine( Real lambda );

    void SetClosed( bool b_closed );
    void SetRadius( Real radius );
    //@}

    void Uniformize(); //!< Refine where needed in order to obtain an approximate uniform sampling
    void Subdivide( int num_levels ); //!< subdivide edges num_levels times
    void Simplify( Real threshold );  //!< remove vertices below error threshold
    void Distort( Real magnitude );   //!< random-dir per-vertex distortion
    void Distort( Real magnitude, Real frequency ); //!< Sinusoidal per-vertex distortion

#ifdef __GEO_ENABLE_IMPLICIT_RBF
    //\name TEMPORAL: Parameters for tweaking
    //@{
    void GetImplicitRBF_Params( Real &reg_coeff, int32 &nf, Real &normal_offset, Real &normal_tension,
                                int32 &max_iter, Real &epsilon ) const;
    void SetImplicitRBF_Params( Real reg_coeff, int32 nf, Real normal_offset, Real normal_tension,
                                int32 max_iter, Real epsilon );
#endif
    //@}

private:
    bool Validate();
    void RecomputeNormals();
    void RecomputeTangents();
    void RecomputeLambdas();
    void RecomputeFlags();

private:
    unsigned int m_MaxVertices;
    Flags32 m_ExplicitVDL;
    bool m_bIsEditing;

#ifdef __GEO_ENABLE_IMPLICIT_RBF
    //-- Tweaks
    //Common
    Real m_RegularizationCoeff;
    int32 m_NormalFlags;
    Real m_NormalOffset;
    Real m_NormalTension;
    //Iterative
    int32 m_MaxIter;
    Real m_Epsilon;
#endif
};

//! 3D
class PolygonalShape3: public GPolygonalShape<3>
{
public:
    PolygonalShape3( unsigned int count = 0, bool b_closed = false, Real radius = 0 )
    : GPolygonalShape<3>(count,b_closed,radius) {}
    EShapeType GetType() const { return eShape_Polygonal3; }
    /*TEMP:
    void ComputeBV( bv::IBoundingVolume &bv,
                    const GPolygonalShape<3>::transform_type &transform,
                    const GPolygonalShape<3>::sdof_type *vec_sdof ) const;
    */
};

} //namespace geo

#endif // GEO_SHAPE_GPOLYGONALSHAPE_H
