#ifndef GEO_SHAPE_PATHSHAPE2_H
#define GEO_SHAPE_PATHSHAPE2_H

#include <Geo/shape/IShape.h>
#include <Geo/bv/bv.h> // BV types for ComputeBV()
#include <vector>

namespace geo
{

//\todo CONSIDER MOVING IT SOMEWHERE ELSE...
inline Vec2 Eval_Bezier3( const Vec2 &p0, const Vec2 &p1, const Vec2 &bezier3_a, const Vec2 &bezier3_b, Real lambda )
{
    //return p0 + mal::Clamp01(lambda) * (p1-p0) + Vec2(0,mal::Sin(lambda));
    Real l( mal::Clamp01(lambda) );
    Real sq_l( l*l );
    Real one_minus_l( Real(1) - l );
    Real sq_one_minus_l( one_minus_l * one_minus_l );
    return (sq_one_minus_l*one_minus_l) * p0
        + (Real(3) * sq_one_minus_l * l) * bezier3_a
        + (Real(3) * one_minus_l * sq_l) * bezier3_b
        + (sq_l * l) * p1;
}

/* Path Shape2 2d with Line and Bezier edges
   Based in MeshSolidShape2 API
   \todo Consider adding other GPolygonalShape functionality such as
   Lambdas, Normals, etc... or even unifying both shapes BUT
   preserving current PathShape2 API, as GPolygonalShape API is
   DEPRECATED at the time of writing.
*/
class PathShape2: public IShape2
{
public:
    enum EConstants
    {
        cVersion = 1 //Initial version
    };
public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = 2;
    typedef mal::GTransform<Real,2> transform_type;
    typedef mal::GVec<Real,2> sdof_type;
    typedef mal::GVec<Real,2> vec_type;
    //@}

public:
    //struct vertex_type {}; //\todo unnecessary
    struct edge_type
    {
        enum ECurveType { eCT_Line, eCT_Bezier3 } m_CurveType; //\todo Bezier2? Hermite3? Arc? \todo CONSIDER using NURBS2/3 instead of Arc!!
        Vec2 m_ParamA, m_ParamB; // Bezier points, H3 tangents..., Arc (radius, ccw_angle0, ccw_angle1)
    };

public:
    PathShape2();
    ~PathShape2();

    //!\name Dynamic IShape implementation
    //@{
    EShapeType GetType() const { return eShape_Path2; }
    unsigned int GetNumSDOF() const { return m_NumV; }
    const sdof_type *GetVecDefaultSDOF() const { return m_vecPoints; }
    const Real *GetVecDefaultDOF() const { return reinterpret_cast<const Real*>(m_vecPoints); }
    void ComputeBVD( bv::BoundingVolume2 &bv, const transform_type &transform, const sdof_type *vec_sdof ) const;

    IDomainSampler2 *CreateDomainSampler() const { return 0; }
    //@}

    //! Init from external arrays
    void SetBakedData( bool b_shared,
                       uint32 num_v, uint32 num_e, bool b_closed,
                       const Vec2 *vec_points, const edge_type *vec_edges );

    inline uint32 GetNumV() const { return m_NumV; }
    inline uint32 GetNumE() const { return m_NumE; }

    inline uint32 GetNumAllocV() const { return m_NumAllocV; }
    inline uint32 GetNumAllocE() const { return m_NumAllocE; }

    inline const Vec2 *GetVecPoints() const { return m_vecPoints; }
//    inline const vertex_type *GetVecV() const { return m_vecV; }
    inline const edge_type *GetVecE() const { return m_vecE; }

    // \todo
    inline bool IsClosed() const { return m_bClosed; }
    inline bool IsSimplyConnected() const { return true; } //!< \todo Check! false if there are holes
    inline bool IsConnected() const { return true; } //!< Must be, by construction
    inline bool IsConvex() const { return false; } //!< \todo

    // Vertex access
    inline Vec2 V_Pos_0( uint32 vid ) const { return m_vecPoints[vid]; }
    inline Vec2 V_Pos( uint32 vid, const sdof_type *vec_sdof ) const { return vec_sdof[vid]; }
    // Edge access
    inline uint32 E_OriginVID( uint32 eid ) const { return eid; }
    inline uint32 E_FinalVID( uint32 eid ) const { return eid+1; }
    inline const edge_type &E_Data( uint32 eid ) const { return m_vecE[eid]; }
    inline edge_type::ECurveType E_CurveType( uint32 eid ) const { return m_vecE[eid].m_CurveType; }
    // Barycenter
    Vec2 Barycenter_0() const;
    Vec2 Barycenter( const sdof_type *vec_sdof ) const;

protected:
    void ClearBakedData();

protected:

    uint32 m_NumV;
    uint32 m_NumE;
    uint32 m_NumAllocV;  //= m_NumV
    uint32 m_NumAllocE;

    bool m_bClosed;

    const Vec2 *m_vecPoints;
    //const vertex_type *m_vecV;
    const edge_type *m_vecE;

    uint32 *m_pBuffer; //!< If exists, the shape is NOT SHARED
};

/* Editable Path Shape 2D
   \note Based in EditableMeshSolidShape2 API
*/
class EditablePathShape2: public PathShape2
{
public:
    EditablePathShape2();
    ~EditablePathShape2();

    void Clear();
    void Set( const PathShape2 &ps2 );
    void BeginEdition(); //\note Incremental, preserves BakedData if any, call Clear() to reset completely
    bool EndEdition();

    feature_index_type SetFirstPoint( const Vec2 &point );
    feature_index_type AddLineTo( const Vec2 &point );
    feature_index_type AddBezier3To( const Vec2 &point, const Vec2 &bezier3_param_a, const Vec2 &bezier3_param_b );
    //feature_index_type AddBezier2To( const Vec2 &point, const Vec2 &bezier2_param );
    //feature_index_type AddHermite3To( const Vec2 &point, const Vec2 &hermite3_param_a, const Vec2 &hermite3_param_b ); //\todo STORE as bezier3 to simplify BV?!?!?!
    void Close(); //\note Does nothing if last-first vertex are coincident, does AddLineTo between them otherwise
    void Transform( const Transform2 &tr );
    void Scale( Real scale );

private:

    void ClearEditData();
    void ClearBakedData();
    void RebuildBakedData();
    bool FixDegeneracies();
    bool FixDegenerateEdges(); //!< Collapse too-short and aligned edges
    feature_index_type FindV( const Vec2& pos, Real epsilon_sq ) const;

private:

    struct editable_vertex_type//: public vertex_type
    {
        editable_vertex_type( const Vec2 &pos ) : m_Pos(pos) {}
        Vec2 m_Pos;
    };
    struct editable_edge_type: public edge_type
    {
        editable_edge_type( edge_type::ECurveType ct, const Vec2 &param_a, const Vec2 &param_b )
            { m_CurveType = ct; m_ParamA = param_a; m_ParamB = param_b; }
    };

    std::vector<editable_vertex_type> m_addV;
    std::vector<editable_edge_type> m_addE;
};

//!\name Make methods
//@{
void Make_PathShape2_Box( EditablePathShape2 &epss, const Vec2 &half_sizes );
void Make_PathShape2_Mushroom( EditablePathShape2 &epss, const Vec2 &half_sizes );
bool Make_PathShape2_SVG( EditablePathShape2 &eps, const char *svg_string );
//@}

} //namespace geo

#endif // GEO_SHAPE_PATHSHAPE2_H
