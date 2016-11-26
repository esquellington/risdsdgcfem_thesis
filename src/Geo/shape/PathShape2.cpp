#include "PathShape2.h"
#include <memory.h> //TEMPORAL req by memcpy()
#include <stdio.h> //req by sscanf()

namespace geo
{

//-----------------------------------------------------------------------------
//---- PathShape2 Implementation
//-----------------------------------------------------------------------------

const Real cLengthEpsilonSq(0.00001f);

PathShape2::PathShape2()
: m_NumV(0), m_NumE(0)
, m_NumAllocV(0), m_NumAllocE(0)
, m_bClosed(false)
, m_vecPoints(0), m_vecE(0)
, m_pBuffer(0)
{}

PathShape2::~PathShape2()
{
    ClearBakedData();
}

void PathShape2::ClearBakedData()
{
    if( m_pBuffer ) delete [] m_pBuffer;
    m_pBuffer = 0;
    m_NumV = 0; m_vecPoints = 0;
    m_NumE = 0; m_vecE = 0;
    m_NumAllocV = 0; m_NumAllocE = 0;
}

void PathShape2::ComputeBVD( bv::BoundingVolume2 &bv, const transform_type &transform, const sdof_type *vec_sdof ) const
{
    const sdof_type *actual_sdof( ( 0 != vec_sdof ) ? vec_sdof : m_vecPoints );
    switch( bv.GetType() )
    {
    case bv::eBV_Sphere2:
        GEO_ASSERT( false ); //Not yet implemented
        //bv.As<bv::Sphere2>().SetPosRadius( transform.m_Pos, m_Radius );
        break;
    case bv::eBV_AABB2:
        {
            // Vertices
            bv::AABB2 aabb( transform*actual_sdof[0], vec_type(0) );
            for( unsigned int it_v=1; it_v < m_NumV; it_v++ )
                aabb.Merge( transform*actual_sdof[it_v] );
            // Other points (Bezier, etc...)
            for( unsigned int it_e=0; it_e < m_NumE; it_e++ )
                switch( m_vecE[it_e].m_CurveType )
                {
                case edge_type::eCT_Bezier3:
                    aabb.Merge( transform*m_vecE[it_e].m_ParamA );
                    aabb.Merge( transform*m_vecE[it_e].m_ParamB );
                    break;
                case edge_type::eCT_Line:
                default:
                    break;
                }
            bv.As<bv::AABB2>() = aabb;
        }
        break;
    case bv::eBV_LSS2:
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

/* Init from external arrays

   If shared, arrays are NOT allocated and copied, only
   referenced. Otherwise all data allocated in m_pBuffer and copied,
   and the shape is self-contained.
*/
void PathShape2::SetBakedData( bool b_shared,
                               uint32 num_v, uint32 num_e, bool b_closed,
                               const Vec2 *vec_points,
                               const edge_type *vec_e )
{
    GEO_ASSERT( m_NumV == 0 && m_vecPoints == 0 && num_v > 0
                && m_NumE == 0 && m_vecE == 0 && num_e > 0
                && m_pBuffer == 0 );
    m_NumV = num_v;
    m_NumE = num_e;
    m_NumAllocV = num_v;
    m_NumAllocE = num_e;
    m_bClosed = b_closed;
    if( b_shared )
    {
        m_vecPoints = vec_points;
        m_vecE = vec_e;
    }
    else
    {
        // Alloc and fill single buffer
        size_t size_points_4aligned = 4*((sizeof(Vec2)*m_NumAllocV+3)/4);
        size_t size_e_4aligned = 4*((sizeof(edge_type)*m_NumAllocE+3)/4);
        size_t total_size_4aligned = size_points_4aligned + size_e_4aligned;
        m_pBuffer = new uint32[ total_size_4aligned ];
        Vec2 *p_buffer_points = reinterpret_cast<Vec2*>( &m_pBuffer[0] );
        edge_type *p_buffer_e = reinterpret_cast<edge_type*>( &m_pBuffer[ size_points_4aligned ] );
        memcpy( p_buffer_points, vec_points, sizeof(Vec2)*m_NumAllocV );
        memcpy( p_buffer_e, vec_e, sizeof(edge_type)*m_NumAllocE );
        // Save const pointers
        m_vecPoints = p_buffer_points;
        m_vecE = p_buffer_e;
    }

    /*TEMPORAL
    GEO_LOG_WARNING( "SetBakedData( %d, %d, %d; %d, %d )", num_v, num_p, num_he, num_boundary_p, num_boundary_he );
    for( unsigned int it_p=GetNumP(); it_p < GetNumP() + GetNumBoundaryP(); it_p++ )
    {
        GEO_LOG_WARNING( "BP[%d] = (%d,%d)", it_p, P_FirstHEID(it_p), P_NumEdges(it_p) );
    }
    */
}

// Barycenter
Vec2 PathShape2::Barycenter_0() const
{
    return Barycenter( GetVecDefaultSDOF() );
}
Vec2 PathShape2::Barycenter( const sdof_type *vec_sdof ) const
{
    //\todo Hackish version, ignores curves, ignores IsClosed()... just computes average vtx pos
    Vec2 acc_pos(0,0);
    for( unsigned int i=0; i<m_NumV; i++ ) acc_pos += vec_sdof[i];
    return acc_pos * mal::Rcp<Real>(m_NumV);
}

//-----------------------------------------------------------------------------
//---- EditablePathShape2 Implementation
//-----------------------------------------------------------------------------
EditablePathShape2::EditablePathShape2()
{
}

EditablePathShape2::~EditablePathShape2()
{
    Clear();
}

void EditablePathShape2::Set( const PathShape2 &ps2 )
{
    Clear();
    // Copy baked data
    for( unsigned int i=0; i<ps2.GetNumV(); i++ )
        m_addV.push_back( editable_vertex_type( ps2.GetVecPoints()[i] ) );
    for( unsigned int i=0; i<ps2.GetNumE(); i++ )
        m_addE.push_back( editable_edge_type( ps2.GetVecE()[i].m_CurveType,
                                              ps2.GetVecE()[i].m_ParamA,
                                              ps2.GetVecE()[i].m_ParamB ) );
    EndEdition();
}

void EditablePathShape2::Clear()
{
    ClearEditData();
    PathShape2::ClearBakedData();
}

void EditablePathShape2::BeginEdition()
{
    // Clear anything added outside Begin/End
    ClearEditData();
    // Copy baked data
    for( unsigned int i=0; i<m_NumV; i++ )
        m_addV.push_back( editable_vertex_type( m_vecPoints[i] ) );
    for( unsigned int i=0; i<m_NumE; i++ )
        m_addE.push_back( editable_edge_type( m_vecE[i].m_CurveType,
                                              m_vecE[i].m_ParamA,
                                              m_vecE[i].m_ParamB ) );
    // Delete baked data
    PathShape2::ClearBakedData();
}

bool EditablePathShape2::EndEdition()
{
    FixDegeneracies();
    //RebuildTopology();
    RebuildBakedData(); //\todo Copy all dyn data to const arrays
    ClearEditData();
    return true;
}

feature_index_type EditablePathShape2::SetFirstPoint( const Vec2 &point )
{
    GEO_ASSERT( m_addV.size() == 0 );
    m_addV.push_back( editable_vertex_type( point ) );
    return feature_index_type(m_addV.size()-1);
}

feature_index_type EditablePathShape2::AddLineTo( const Vec2 &point )
{
    //\todo Consider looking for existing V
    GEO_ASSERT( m_addV.size() < cInvalidFeatureIndex-1 );
    if( m_addV.size() == 0 )
    {
        GEO_LOG_WARNING("EditablePathShape2::AddLineTo() with 0 points, interpreting as SetFirstPoint()");
        return SetFirstPoint(point);
    }
    else
    {
        m_addV.push_back( editable_vertex_type( point ) );
        m_addE.push_back( editable_edge_type( edge_type::eCT_Line, Vec2::Zero(), Vec2::Zero() ) );
        return feature_index_type(m_addV.size()-1);
    }
}

feature_index_type EditablePathShape2::AddBezier3To( const Vec2 &point, const Vec2 &bezier3_param_a, const Vec2 &bezier3_param_b )
{
    //\todo Consider looking for existing V
    GEO_ASSERT( m_addV.size() < cInvalidFeatureIndex-1 );
    m_addV.push_back( editable_vertex_type( point ) );
    m_addE.push_back( editable_edge_type( edge_type::eCT_Bezier3, bezier3_param_a, bezier3_param_b ) );
    return feature_index_type(m_addV.size()-1);
}

void EditablePathShape2::Close()
{
    GEO_ASSERT( m_addE.size() > 1 );
    // Enforce strict equality if first/last are close enough, AddLineTo between them otherwise
    if( mal::NormSq( m_addV.front().m_Pos - m_addV.back().m_Pos ) <= cLengthEpsilonSq )
        m_addV.back().m_Pos = m_addV.front().m_Pos;
    else
        AddLineTo( m_addV.front().m_Pos );
}

void EditablePathShape2::Transform( const Transform2 &tr )
{
    // Transform vertices
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        m_addV[it_v].m_Pos = tr * m_addV[it_v].m_Pos;
    // Other points (Bezier, etc...)
    for( unsigned int it_e=0; it_e < m_addE.size(); it_e++ )
        switch( m_addE[it_e].m_CurveType )
        {
        case edge_type::eCT_Bezier3:
            m_addE[it_e].m_ParamA = tr*m_addE[it_e].m_ParamA;
            m_addE[it_e].m_ParamB = tr*m_addE[it_e].m_ParamB;
            break;
        case edge_type::eCT_Line:
        default:
            break;
        }
}

void EditablePathShape2::Scale( Real scale )
{
    // Transform vertices
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        m_addV[it_v].m_Pos *= scale;
    // Other points (Bezier, etc...)
    for( unsigned int it_e=0; it_e < m_addE.size(); it_e++ )
        switch( m_addE[it_e].m_CurveType )
        {
        case edge_type::eCT_Bezier3:
            m_addE[it_e].m_ParamA *= scale;
            m_addE[it_e].m_ParamB *= scale;
            break;
        case edge_type::eCT_Line:
        default:
            break;
        }
}

//---- Internal methods
void EditablePathShape2::ClearEditData()
{
    m_addV.clear();
    m_addE.clear();
}

void EditablePathShape2::RebuildBakedData()
{
    // Bake Vertices
    Vec2 *vecPoints( new Vec2[m_addV.size()] );
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        vecPoints[it_v] = m_addV[it_v].m_Pos;
    // Bake Edges
    edge_type *vecE( new edge_type[m_addE.size()] );
    for( unsigned int it_e=0; it_e<m_addE.size(); it_e++ )
        vecE[it_e] = m_addE[it_e];
    // Set baked data as non-shared
    bool b_closed( m_addE.size() > 2
                   &&
                   m_addV.front().m_Pos == m_addV.back().m_Pos );
                   //mal::NormSq( m_addV.front().m_Pos - m_addV.back().m_Pos ) <= cLengthEpsilonSq );
    SetBakedData( false, m_addV.size(), m_addE.size(), b_closed, vecPoints, vecE );
    // Clear temporal stuff
    delete [] vecPoints;
    delete [] vecE;
}

bool EditablePathShape2::FixDegeneracies()
{
    bool bFixDE = FixDegenerateEdges();
    return bFixDE;
};
bool EditablePathShape2::FixDegenerateEdges() { return false; }

feature_index_type EditablePathShape2::FindV( const Vec2& pos, Real epsilon_sq ) const
{
    for( unsigned int it_v=0; it_v<m_addV.size(); it_v++ )
        if( mal::NormSq( m_addV[it_v].m_Pos - pos ) <= epsilon_sq )
            return feature_index_type(it_v);
    return cInvalidFeatureIndex;
}

//-----------------------------------------------------------------------------
//---- Make_PathShape2_XXX Implementation
//-----------------------------------------------------------------------------

void Make_PathShape2_Box( EditablePathShape2 &epss, const Vec2 &half_sizes )
{
    epss.Clear();
    epss.BeginEdition();
    {
        feature_index_type vid0 = epss.SetFirstPoint( Vec2(-half_sizes.x(),-half_sizes.y()) );
        feature_index_type vid1 = epss.AddLineTo( Vec2( half_sizes.x(),-half_sizes.y()) );
        feature_index_type vid2 = epss.AddLineTo( Vec2( half_sizes.x(), half_sizes.y()) );
        feature_index_type vid3 = epss.AddLineTo( Vec2(-half_sizes.x(), half_sizes.y()) );
        epss.Close();
    }
    epss.EndEdition();
}

void Make_PathShape2_Mushroom( EditablePathShape2 &epss, const Vec2 &half_sizes )
{
    epss.Clear();
    epss.BeginEdition();
    {
        feature_index_type vid0 = epss.SetFirstPoint( Vec2(-half_sizes.x(),-half_sizes.y()) );
        feature_index_type vid1 = epss.AddLineTo( Vec2( half_sizes.x(),-half_sizes.y()) );
        feature_index_type vid2 = epss.AddLineTo( Vec2( half_sizes.x(), half_sizes.y()) );
        feature_index_type vid3 = epss.AddBezier3To( Vec2(-half_sizes.x(), half_sizes.y()),
                                                     Vec2( 8*half_sizes.x(), 3*half_sizes.y()),
                                                     Vec2(-8*half_sizes.x(), 3*half_sizes.y()) );
        epss.Close();
    }
    epss.EndEdition();
}

/* Makes a Path from an SVG path element definition
   \note Assumes svg_string begins at the first command in a SVG path definition d="command params..."
   \note Code adapted from libsvgtiny svgtiny_parse_path() function, see http://www.cplusplus.com/reference/clibrary/cstdio/scanf/
   \note sscanf() format from svgtiny_parse_path() extended to support ',' and ' ' separation between float components.
   \todo Adapt it to use Vec2f and make it look prettier.
*/
bool Make_PathShape2_SVG( EditablePathShape2 &eps, const char *svg_string )
{
    eps.Clear();
    eps.BeginEdition();

    float last_x = 0, last_y = 0;
    float last_cubic_x = 0, last_cubic_y = 0;
    float last_quad_x = 0, last_quad_y = 0;
    unsigned int num_commands(0);

    // parse d and build path
    const char *s( svg_string );
    while( *s )
    {
        char command[2];
        int plot_command;
        float x, y, x1, y1, x2, y2, rx, ry, rotation, large_arc, sweep;
        int n;
        if( sscanf(s, " %1[Mm] %f%*[, ]%f %n", command, &x, &y, &n) == 3 ) // moveto (M, m) (2 arguments)
        {
            /*LOG(("moveto"));*/
            do
            {
                if( *command == 'm' ) //lower case => relative coords
                {
                    x += last_x;
                    y += last_y;
                }
                last_cubic_x = last_quad_x = last_x = x;
                last_cubic_y = last_quad_y = last_y = y;
                s += n;
                // Sometimes Move is used as LineTo (Inkscape), so we support it...
                if( num_commands == 0 ) eps.SetFirstPoint( Vec2(x,y) );
                else eps.AddLineTo( Vec2(x,y) );
                num_commands++;
            } while( sscanf(s, "%f%*[, ]%f %n", &x, &y, &n) == 2 );
        }
        else if( sscanf(s, " %1[Ll] %f%*[, ]%f %n", command, &x, &y, &n) == 3 ) //lineto (L, l) (2 arguments)
        {
            /*LOG(("moveto or lineto"));*/
            GEO_ASSERT( num_commands > 0 );
            do
            {
                if( *command == 'l' ) //lower case => relative coords
                {
                    x += last_x;
                    y += last_y;
                }
                last_cubic_x = last_quad_x = last_x = x;
                last_cubic_y = last_quad_y = last_y = y;
                s += n;
                eps.AddLineTo( Vec2(x,y) );
                num_commands++;
            } while( sscanf(s, "%f%*[, ]%f %n", &x, &y, &n) == 2 );
        }
        else if( sscanf(s, " %1[Zz] %n", command, &n) == 1 ) //closepath (Z, z) (no arguments)
        {
            /*LOG(("closepath"));*/
            GEO_ASSERT( num_commands > 0 );
            s += n;
            eps.Close();
        }
        else if (sscanf(s, " %1[Hh] %f %n", command, &x, &n) == 2) // horizontal lineto (H, h) (1 argument)
        {
            /*LOG(("horizontal lineto"));*/
            GEO_ASSERT( num_commands > 0 );
            do
            {
                if( *command == 'h' ) x += last_x;
                last_cubic_x = last_quad_x = last_x = x;
                last_cubic_y = last_quad_y = last_y;
                s += n;
                eps.AddLineTo( Vec2(x,last_y) );
                num_commands++;
            } while( sscanf(s, "%f %n", &x, &n) == 1 );
        }
        else if( sscanf(s, " %1[Vv] %f %n", command, &y, &n) == 2 ) // vertical lineto (V, v) (1 argument)
        {
            /*LOG(("vertical lineto"));*/
            GEO_ASSERT( num_commands > 0 );
            do
            {
                if( *command == 'v' ) y += last_y;
                last_cubic_x = last_quad_x = last_x;
                last_cubic_y = last_quad_y = last_y = y;
                s += n;
                eps.AddLineTo( Vec2(last_x,y) );
                num_commands++;
            } while( sscanf(s, "%f %n", &x, &n) == 1 );
        }
        else if( sscanf(s, " %1[Cc] %f%*[, ]%f%*[, ]%f%*[, ]%f%*[, ]%f%*[, ]%f %n", command,
                        &x1, &y1, &x2, &y2, &x, &y, &n) == 7 ) // curveto (C, c) (6 arguments)
        {
            /*LOG(("curveto"));*/
            GEO_ASSERT( num_commands > 0 );
            do
            {
                if( *command == 'c' )
                {
                    x1 += last_x;
                    y1 += last_y;
                    x2 += last_x;
                    y2 += last_y;
                    x += last_x;
                    y += last_y;
                }
                x1;
                y1;
                last_cubic_x = x2;
                last_cubic_y = y2;
                last_quad_x = last_x = x;
                last_quad_y = last_y = y;
                s += n;
                eps.AddBezier3To( Vec2(x,y), Vec2(x1,y1), Vec2(x2,y2) );
                num_commands++;
            } while( sscanf(s, "%f%*[, ]%f%*[, ]%f%*[, ]%f%*[, ]%f%*[, ]%f %n",
                            &x1, &y1, &x2, &y2, &x, &y, &n) == 6 );
        }
#ifdef __DISABLE_UNSUPPORTED_COMMANDS //TEMPORAL: Consider adding them...
        else if (sscanf(s, " %1[Ss] %f %f %f %f %n", command,
                        &x2, &y2, &x, &y, &n) == 5) // shorthand/smooth curveto (S, s) (4 arguments)
        {
            /*LOG(("shorthand/smooth curveto"));*/
            do {
                p[i++] = svgtiny_PATH_BEZIER;
                x1 = last_x + (last_x - last_cubic_x);
                y1 = last_y + (last_y - last_cubic_y);
                if (*command == 's') {
                    x2 += last_x;
                    y2 += last_y;
                    x += last_x;
                    y += last_y;
                }
                p[i++] = x1;
                p[i++] = y1;
                p[i++] = last_cubic_x = x2;
                p[i++] = last_cubic_y = y2;
                p[i++] = last_quad_x = last_x = x;
                p[i++] = last_quad_y = last_y = y;
                s += n;
            } while (sscanf(s, "%f %f %f %f %n",
                            &x2, &y2, &x, &y, &n) == 4);

            /* quadratic Bezier curveto (Q, q) (4 arguments) */
        } else if (sscanf(s, " %1[Qq] %f %f %f %f %n", command,
                          &x1, &y1, &x, &y, &n) == 5) {
            /*LOG(("quadratic Bezier curveto"));*/
            do {
                p[i++] = svgtiny_PATH_BEZIER;
                last_quad_x = x1;
                last_quad_y = y1;
                if (*command == 'q') {
                    x1 += last_x;
                    y1 += last_y;
                    x += last_x;
                    y += last_y;
                }
                p[i++] = 1./3 * last_x + 2./3 * x1;
                p[i++] = 1./3 * last_y + 2./3 * y1;
                p[i++] = 2./3 * x1 + 1./3 * x;
                p[i++] = 2./3 * y1 + 1./3 * y;
                p[i++] = last_cubic_x = last_x = x;
                p[i++] = last_cubic_y = last_y = y;
                s += n;
            } while (sscanf(s, "%f %f %f %f %n",
                            &x1, &y1, &x, &y, &n) == 4);

            /* shorthand/smooth quadratic Bezier curveto (T, t)
               (2 arguments) */
        } else if (sscanf(s, " %1[Tt] %f %f %n", command,
                          &x, &y, &n) == 3) {
            /*LOG(("shorthand/smooth quadratic Bezier curveto"));*/
            do {
                p[i++] = svgtiny_PATH_BEZIER;
                x1 = last_x + (last_x - last_quad_x);
                y1 = last_y + (last_y - last_quad_y);
                last_quad_x = x1;
                last_quad_y = y1;
                if (*command == 't') {
                    x1 += last_x;
                    y1 += last_y;
                    x += last_x;
                    y += last_y;
                }
                p[i++] = 1./3 * last_x + 2./3 * x1;
                p[i++] = 1./3 * last_y + 2./3 * y1;
                p[i++] = 2./3 * x1 + 1./3 * x;
                p[i++] = 2./3 * y1 + 1./3 * y;
                p[i++] = last_cubic_x = last_x = x;
                p[i++] = last_cubic_y = last_y = y;
                s += n;
            } while (sscanf(s, "%f %f %n",
                            &x, &y, &n) == 2);

            /* elliptical arc (A, a) (7 arguments) */
        } else if (sscanf(s, " %1[Aa] %f %f %f %f %f %f %f %n", command,
                          &rx, &ry, &rotation, &large_arc, &sweep,
                          &x, &y, &n) == 8) {
            do {
                p[i++] = svgtiny_PATH_LINE;
                if (*command == 'a') {
                    x += last_x;
                    y += last_y;
                }
                p[i++] = last_cubic_x = last_quad_x = last_x
                       = x;
                p[i++] = last_cubic_y = last_quad_y = last_y
                       = y;
                s += n;
            } while (sscanf(s, "%f %f %f %f %f %f %f %n",
                            &rx, &ry, &rotation, &large_arc, &sweep,
                            &x, &y, &n) == 7);

        }
#endif //__DISABLE_UNSUPPORTED_COMMANDS
        else
        {
            GEO_LOG_ERROR( "Load_PathShape2_SVG: parse failed at '%s' ", s);
            return false;
            break;
        }
    }
    return eps.EndEdition();
}

} //namespace geo
