#include "Clip.h"
#include <Geo/np/RayCast.h>
#include <Geo/np/Overlap.h>

namespace geo {
namespace np {

/* Clip a Triangle3 against a Plane3 (Halfspace)
   - Computes a convex polygon resulting from the clipping, with the same orientation as the input triangle.
   - Returns the clipped polygon vertex count [3..4], or 0 if no intersection
   - The whole Triangle returned if wholly inside the halfspace
   \pre vec_clipped_points MUST have, at least, 4 valid slots
*/
unsigned int Clip_Triangle3_Plane3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                    const Vec3& plane_normal, Real plane_coeff_d,
                                    Vec3* vec_polygon_clipped_pos,
                                    const Context* p_context )
{
    // Trivial triangle outside/inside plane cases
    Real vec_dist[3] = { mal::Dot(tri0,plane_normal)+plane_coeff_d,
                         mal::Dot(tri1,plane_normal)+plane_coeff_d,
                         mal::Dot(tri2,plane_normal)+plane_coeff_d };
    if( vec_dist[0] > 0 && vec_dist[1] > 0 && vec_dist[2] > 0 ) return 0;
    else if( vec_dist[0] <= 0 && vec_dist[1] <= 0 && vec_dist[2] <= 0 )
    {
        vec_polygon_clipped_pos[0] = tri0;
        vec_polygon_clipped_pos[1] = tri1;
        vec_polygon_clipped_pos[2] = tri2;
        return 3;
    }
    // Handle triangle-pierces-plane case
    Vec3 vec_segment_pos0[3] = { tri0, tri1, tri2 };
    Vec3 vec_segment_dir[3] = { tri1-tri0, tri2-tri1, tri0-tri2 };
    // Clip each segment against plane and, if intersecting, add begin/endpoints unless repeated
    unsigned int num_cp(0);
    for( unsigned int it_segment=0; it_segment < 3; it_segment++ )
    {
        RayHit3 rh;
        const Vec3& pos0( vec_segment_pos0[it_segment] );
        const Vec3& dir( vec_segment_dir[it_segment] );
        if( GRayCast_HalfSpace<3>( pos0, dir, Interval(0,1), plane_normal, plane_coeff_d, rh ) )
        {
            // Add beginpoint unconditionally
            vec_polygon_clipped_pos[num_cp++] = pos0 + rh.m_Interval.Min()*dir;
            // Add endpoint if unique \todo Consider instead comparing actual beginpoint with last added endpoint, if any, instead of using interval with epsilon
            if( rh.m_Interval.Max() < 1-1e-6 ) vec_polygon_clipped_pos[num_cp++] = pos0 + rh.m_Interval.Max()*dir;
            //
            /*TEMP: This approach seemed safer a priori, but fails much more than interval epsilon...
            // Add beginpoint if unique
            Vec3 s0( pos0 + rh.m_Interval.Min()*dir );
            if( num_cp == 0 || mal::NormSq(s0-vec_polygon_clipped_pos[num_cp]) > 1e-6 ) vec_polygon_clipped_pos[num_cp++] = s0;
            // Add endpoint unconditionally
            vec_polygon_clipped_pos[num_cp++] = pos0 + rh.m_Interval.Max()*dir;
            */
        }
    }
    if( num_cp < 3 || num_cp > 4 )
    {
        GEO_LOG("Clip_Triangle3_Plane3() UNEXPECTED #CP = %u", num_cp );
        return 0;
    }
    GEO_ASSERT( num_cp == 0 || (num_cp >= 3 && num_cp <= 4) );
    return num_cp;
}

/* Clip a Triangle3 against a Tetrahedron3:
   - Computes a convex polygon resulting from the clipping, with the same orientation as the input triangle.
   - Returns the clipped polygon vertex count [3..6], or 0 if no intersection
   \pre vec_clipped_points MUST have, at least, 6 valid positions
   \todo Consider generalization to Clip_Polygon3_Tetrahedron3 with Tri3->Poly3? or even Clip_Polygon3_Polyedron3?
*/
unsigned int Clip_Triangle3_Tetrahedron3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                          const Vec3& tet0, const Vec3& tet1, const Vec3& tet2, const Vec3& tet3,
                                          Vec3* vec_polygon_clipped_pos,
                                          const Context* p_context )
{
    /*TEMP SOMETIMES the clipped poly is WAY LARGER than the original triangle, assert and debug... BUT IT DOES NOT HAPPEN!?!?!?!?
    Vec3 tet_barycenter( mal::Rcp<Real>(4.0)*(tet0+tet1+tet2+tet3) );
    Real tet_extent( mal::Max( mal::Max( mal::Norm(tet0-tet_barycenter),
                                         mal::Norm(tet1-tet_barycenter) ),
                               mal::Max( mal::Norm(tet2-tet_barycenter),
                                         mal::Norm(tet3-tet_barycenter) ) ) );
    tet_extent *= 1.01;
    */

    std::vector<Vec3> vec_tmp_pos; //\todo Subst std::vector with std::array
    geo::np::RayHit3 rh;
    //-- Test 6 tetrahedron edges against the triangle \todo Could optimize, I guess...
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tet0, tet1-tet0, Interval(0,1), tri0, tri1, tri2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tet0, tet2-tet0, Interval(0,1), tri0, tri1, tri2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tet0, tet3-tet0, Interval(0,1), tri0, tri1, tri2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tet1, tet2-tet1, Interval(0,1), tri0, tri1, tri2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tet2, tet3-tet2, Interval(0,1), tri0, tri1, tri2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tet3, tet1-tet3, Interval(0,1), tri0, tri1, tri2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    //-- Test 3 triangle points inside tetrahedron \todo Could optimize as Overlap_Point3_Tetrahedron3 recomputes inv_B on each call!
    if( geo::np::Overlap_Point3_Tetrahedron3( tri0, tet0, tet1, tet2, tet3 ) ) { vec_tmp_pos.push_back( tri0 ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::Overlap_Point3_Tetrahedron3( tri1, tet0, tet1, tet2, tet3 ) ) { vec_tmp_pos.push_back( tri1 ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::Overlap_Point3_Tetrahedron3( tri2, tet0, tet1, tet2, tet3 ) ) { vec_tmp_pos.push_back( tri2 ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    //-- Test 3 triangle edges against 4 tetrahedron faces \todo Could optimize, I guess...
    // tri.e0
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri0, tri1-tri0, Interval(0,1), tet0, tet1, tet2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri0, tri1-tri0, Interval(0,1), tet0, tet2, tet3, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri0, tri1-tri0, Interval(0,1), tet0, tet3, tet1, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri0, tri1-tri0, Interval(0,1), tet1, tet2, tet3, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    // tri.e1
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri1, tri2-tri1, Interval(0,1), tet0, tet1, tet2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri1, tri2-tri1, Interval(0,1), tet0, tet2, tet3, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri1, tri2-tri1, Interval(0,1), tet0, tet3, tet1, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri1, tri2-tri1, Interval(0,1), tet1, tet2, tet3, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    // tri.e2
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri2, tri0-tri2, Interval(0,1), tet0, tet1, tet2, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri2, tri0-tri2, Interval(0,1), tet0, tet2, tet3, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri2, tri0-tri2, Interval(0,1), tet0, tet3, tet1, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    if( geo::np::RayCast_Triangle3_DoubleSided_COHERENT( tri2, tri0-tri2, Interval(0,1), tet1, tet2, tet3, rh ) ) { vec_tmp_pos.push_back( rh.m_Point ); }//GEO_ASSERT(mal::Norm(vec_tmp_pos.back()-tet_barycenter)<=tet_extent); }
    //---- Compute convex hull of clipped points
    unsigned int num_clipped(0);
    if( !vec_tmp_pos.empty() )
    {
        // Compute barycenter \todo May be sligtly off if there's coincidences, but not too much I hope... OTHERWISE, remove duplicates BEFORE GrahamScan, and ENFORCE <=6 cp
        Vec3 barycenter(0,0,0);
        for( unsigned int it_v=0; it_v<vec_tmp_pos.size(); it_v++ )
            barycenter += vec_tmp_pos[it_v];
        barycenter = barycenter / vec_tmp_pos.size();
        // Using barycenter and clipped[0], perform graham scan to build the CCW convex hull
        std::sort( vec_tmp_pos.begin(), vec_tmp_pos.end(),
                   // p<q lambda grijander
                   [ b = barycenter,
                     d0 = mal::SafeNormalized(vec_tmp_pos[0]-barycenter),
                     n = mal::SafeNormalized( mal::Cross(tri1-tri0,tri2-tri0) ) ] //\todo SafeNormalized() added because escalunya crashed sometimes, but consider previous non-degeneracy test/treatment
                   ( const Vec3& p, const Vec3& q )
                   {
                       //\todo This looks EXPENSIVE, could be omtimized I guess...
                       Real cos_p_unnormalized( mal::Dot(p-b,d0) );
                       Real cos_q_unnormalized( mal::Dot(q-b,d0) );
                       Real sin_p_unnormalized( mal::Dot(n,mal::Cross(d0,p-b)) );
                       Real sin_q_unnormalized( mal::Dot(n,mal::Cross(d0,q-b)) );
                       // ATan2 works with unnormalized sin/cos, only ratio is required
                       return mal::ATan2(sin_p_unnormalized,cos_p_unnormalized)
                           < mal::ATan2(sin_q_unnormalized,cos_q_unnormalized);

                       /*TEMP: This version failed in some cases, and required length computation
                       Real dist_pb( mal::Norm(p-b) );
                       Real dist_qb( mal::Norm(q-b) );
                       Real cos_p( mal::Dot(p-b,d0)/dist_pb );
                       Real cos_q( mal::Dot(q-b,d0)/dist_qb );
                       Real sin_p( mal::Dot(n,mal::Cross(d0,p-b))/dist_pb );
                       Real sin_q( mal::Dot(n,mal::Cross(d0,q-b))/dist_qb );
                       if( sin_p * sin_q < 0 ) return sin_p >= sin_q; //different sides from d0, angle_p < angle_q iff sin_p > sin_q
                       else if( sin_p > 0 || sin_q > 0 ) return cos_p >= cos_q; //left side, largest cos(angle) yields smallest angle wrt d0
                       else return cos_p <= cos_q; ////right side, smallest cos(angle) yields smallest angle wrt d0
                       */
                   } );
        // Remove near-coincident points leveraging sorted array
        if( vec_tmp_pos.size()>6 ) GEO_LOG_WARNING("Clip_Triangle3_Tetrahedron3 with %ld > 6 points, merging coincident", vec_tmp_pos.size() );
        Real epsilon_length_sq( p_context->m_Epsilon_LengthSq );
        unsigned int num_iter(0);
        unsigned int num_coincidences(0);
        do
        {
            num_clipped = 1;
            while( num_clipped < vec_tmp_pos.size() )
            {
                if( mal::NormSq(vec_tmp_pos[num_clipped]-vec_tmp_pos[num_clipped-1]) < epsilon_length_sq )
                {
                    num_coincidences++;
                    for( unsigned int i=num_clipped; i<vec_tmp_pos.size()-1; i++ )
                        vec_tmp_pos[i] = vec_tmp_pos[i+1];
                    vec_tmp_pos.pop_back();
                }
                else
                    num_clipped++;
            }
            // Remove last if coincident with first
            if( mal::NormSq(vec_tmp_pos.front()-vec_tmp_pos.back()) < p_context->m_Epsilon_LengthSq )
            {
                num_coincidences++;
                vec_tmp_pos.pop_back();
                num_clipped--;
            }
            // 2xEpsilon for next iter
            epsilon_length_sq *= 2;
            num_iter++;
        } while( num_clipped > 6 && num_iter < 10 );
        if( num_iter > 1 ) GEO_LOG_WARNING("Clip_Triangle3_Tetrahedron3 required %d > 1 merge coincident iter, final eps = %f", num_iter, mal::Sqrt(epsilon_length_sq) );

        /* TEMP: debug large edges....
        for( auto it_pos : vec_tmp_pos )
        {
            GEO_ASSERT(mal::Norm(it_pos-tet_barycenter)<=tet_extent);
            // GEO_ASSERT(mal::Norm(it_pos-tet_barycenter)<=0.4);
        }
        for( unsigned int it_v=0; it_v<vec_tmp_pos.size(); it_v++ )
        {
            GEO_ASSERT(mal::Norm(vec_tmp_pos[it_v]-vec_tmp_pos[(it_v+1)%vec_tmp_pos.size()])<=0.25);
        }
        */

        // Build final polygon
        num_clipped = vec_tmp_pos.size();
        //\todo This MAY HAPPEN if eps is too small, the ONLY way to ensure num_clipped<=6 is to enforce it by merging pairs of closest points until num_clipped<=6
        GEO_LOG_ASSERT( num_clipped <= 6, "Clip_Triangle3_Tetrahedron3 with %d > 6 points AFTER merging %d coincidences with eps^2=%f",
                        num_clipped, num_coincidences, epsilon_length_sq );
        for( unsigned int it_v=0; it_v<num_clipped; it_v++ )
            vec_polygon_clipped_pos[it_v] = vec_tmp_pos[it_v];
    }
    return num_clipped;
}

#ifdef __DISABLED_WHILE_DEVELOPING
/* Simple HalfSpace class
   - stores (normal,coeff_d)
   - normal is NOT NECESSARILY unitary
   \todo Put somewhere else globally accessible
*/
template<unsigned D>
class GHalfSpace
{
public:
    enum EConstants { cDimension = D };
    typedef mal::GVec<Real,D> vec_type;
public:
    inline GHalfSpace() : m_Normal(0), m_CoeffD(0) {}
    inline GHalfSpace( const vec_type& non_unit_normal, Real coeff_d ) : m_Normal( mal::Normalized(non_unit_normal) ), m_CoeffD(coeff_d) {}
    inline GHalfSpace( const vec_type& non_unit_normal, const vec_type& point ) : m_Normal(mal::Normalized(non_unit_normal)), m_CoeffD( -mal::Dot(m_Normal,point) ) {}
    inline const vec_type& GetNormal() const { return m_Normal; }
    inline Real GetCoeffD() const { return m_CoeffD; }
private:
    vec_type m_Normal;
    Real m_CoeffD;
};
typedef GHalfSpace<2> HalfSpace2;
typedef GHalfSpace<3> HalfSpace3;

//\todo MOVE TO Overlap.h when HalfSpace3 is global
bool Overlap_Point3_Polyhedron3( const Vec3& point,
                                 const HalfSpace3* vec_polyhedron_hs, unsigned int num_hs,
                                 const Context* p_context )
{
    for( unsigned int it_hs=0; it_hs<num_hs; it_hs++ )
        if( mal::Dot(point,vec_polyhedron_hs[it_hs].GetNormal()) + vec_polyhedron_hs[it_hs].GetCoeffD() > Real(0) )
            return false;
    return true;
}

/* RayCast against a Polyhedron3
   - If ray_pos is INSIDE the Polyhedron3, rh.m_Point == ray_pos, and
     rh.m_Interval.Min() == interval.Min()
   - Otherwise, rh.m_Point is the ray-polyhedron entry/piercing point,
     and rh.m_Interval.Min() is the entry/piercing lambda
*/
bool RayCast_Polyhedron3( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                          const HalfSpace3* vec_polyhedron_hs, unsigned int num_hs,
                          RayHit3& rh, RayCache3* p_rc = 0 )
{
    GEO_LOG_ASSERT(false, "THIS methid MAY BE BROKEN, review it if ever needed");
    // Init rh interval to original ray interval
    rh.m_Interval = interval;
    for( unsigned int it_hs=0; it_hs<num_hs && !rh.m_Interval.IsEmpty(); it_hs++ )
    {
        // Clip original ray against polyhedron faces halfspaces, with tmp_rh containing the interval overlapping negative-halfspace in case of intersection
        RayHit3 tmp_rh;
        if( GRayCast_HalfSpace<3>( ray_pos, ray_dir, interval,
                                   vec_polyhedron_hs[it_hs].GetNormal(), vec_polyhedron_hs[it_hs].GetCoeffD(),
                                   tmp_rh ) )
        {
            // Clip rh interval against this face tmp_rh interval
            if( tmp_rh.m_Interval.Min() > rh.m_Interval.Min() )
                rh.m_Normal = vec_polyhedron_hs[it_hs].GetNormal();
            rh.m_Interval.Intersect(tmp_rh.m_Interval);
        }
        else //No overlap means ray completely outside hs
            return false;
    }
    if( !rh.m_Interval.IsEmpty() )
    {
        rh.m_Point = ray_pos + ray_dir*rh.m_Interval.Min();
        rh.m_FeatureId = feature_id();
        return true;
    }
    else
        return false;
}

/* CONVENTION: Polyhedron halfspace normals point outwards
  The brute-force approach would be:
  - Find all polygon-edge vs polyhedron-face intersection points
  - Find all polyhedron-edge vs polygon-face intersection points
  - Compute convex hull of all intersection points, using Graham Scan
    in the given polugon orientation (CCW)

  \todo DUE TO NUMERIC ERRORS, THERE MAY BE more clipped points than
  expected... consider using std::vector<> to gather them and remove
  duplicates before storing in vec_polygon_clipped_pos.
*/
unsigned int Clip_Polygon3_Polyhedron3( const Vec3* vec_polygon_pos, unsigned int num_v,
                                        const HalfSpace3* vec_polyhedron_hs, unsigned int num_hs,
                                        const Vec3* vec_polyhedron_edge_p0, const Vec3* vec_polyhedron_edge_p1, unsigned int num_e,
                                        Vec3* vec_polygon_clipped_pos,
                                        const Context* p_context )
{
    GEO_LOG_ASSERT(false, "THIS methid is broken due to polygon edges being inside the polyhedron");
    unsigned int num_clipped(0);
    //---- Gather clipped points from all edge/face pairs
    // Clip polygon edges against polyhedron
    geo::np::RayHit3 rh;
    for( unsigned int it_e=0; it_e < num_v; it_e++ )
        if( RayCast_Polyhedron3( vec_polygon_pos[it_e], vec_polygon_pos[(it_e+1)%num_v] - vec_polygon_pos[it_e], Interval(0,1),
                                 vec_polyhedron_hs, num_hs,
                                 rh ) )
        {
            // \todo BROKEN... this FAILS for poly vertices inside the tet... ray-polyhedron piercing point is Min() if ray_pos outside polyhedron or Max() otherwise
            if( rh.m_Interval.Min() > 0 ) vec_polygon_clipped_pos[num_clipped++] = rh.m_Point;
            else vec_polygon_clipped_pos[num_clipped++] = vec_polygon_pos[it_e] + rh.m_Interval.Max()*(vec_polygon_pos[it_e]-vec_polygon_pos[it_e]);
        }
    /* Clip polygon vertices against polyhedron
       \todo NO! this is redundant, as clipped edges already consider interior vertices!
    for( unsigned int it_v=0; it_v<num_v; it_v++ )
        if( Overlap_Point3_Polyhedron3( vec_polygon_pos[it_v], vec_polyhedron_hs, num_hs, p_context ) )
            vec_polygon_clipped_pos[num_clipped++] = vec_polygon_pos[it_v];
    */
    // Clip polyhedron edges against polygon face
    GEO_ASSERT(num_v == 3); //\todo This could be generalized to actual polygons...
    for( unsigned int it_e=0; it_e < num_e; it_e++ )
        if( geo::np::RayCast_Triangle3_DoubleSided( vec_polyhedron_edge_p0[it_e], vec_polyhedron_edge_p1[it_e]-vec_polyhedron_edge_p0[it_e],
                                                    Interval(0,1),
                                                    vec_polygon_pos[0], vec_polygon_pos[1], vec_polygon_pos[2],
                                                    rh ) )
            vec_polygon_clipped_pos[num_clipped++] = rh.m_Point;//vec_polyhedron_edge_p0[it_e] + rh.m_Interval.Min()*(vec_polyhedron_edge_p1[it_e]-vec_polyhedron_edge_p0[it_e]);
    //---- Compute convex hull of clipped points
    if( num_clipped > 0 )
    {
        // Compute barycenter
        Vec3 barycenter(0,0,0);
        for( unsigned int it_v=0; it_v<num_clipped; it_v++ )
            barycenter += vec_polygon_clipped_pos[it_v];
        barycenter = barycenter / num_clipped;
        // Using barycenter and clipped[0], perform graham scan to build the CCW convex hull
        std::vector<Vec3> vec_tmp_pos;
        for( unsigned int it_v=0; it_v<num_clipped; it_v++ )
            vec_tmp_pos.push_back(vec_polygon_clipped_pos[it_v]);
        std::sort( vec_tmp_pos.begin(), vec_tmp_pos.end(),
                   // p<q lambda grijander
                   [ b = barycenter,
                     d0 = vec_polygon_pos[0]-barycenter,
                     n = mal::Normalized( mal::Cross(vec_polygon_pos[1]-vec_polygon_pos[0],vec_polygon_pos[2]-vec_polygon_pos[0]) ) ]
                   ( const Vec3& p, const Vec3& q )
                   {
                       Real cos_p( mal::Dot(p-b,d0) );
                       Real cos_q( mal::Dot(q-b,d0) );
                       Real sin_p( mal::Dot(n,mal::Cross(d0,p-b)) );
                       Real sin_q( mal::Dot(n,mal::Cross(d0,q-b)) );
                       //\todo THIs IS WRONG USE ATan2 as in specific ClipTri2Tet!
                       if( sin_p * sin_q < 0 ) return sin_p >= sin_q; //different sides from d0, angle_p < angle_q iff sin_p > sin_q
                       else if( sin_p > 0 || sin_q > 0 ) return cos_p >= cos_q; //left side, largest cos(angle) yields smallest angle wrt d0
                       else return cos_p <= cos_q; ////right side, smallest cos(angle) yields smallest angle wrt d0
                   } );
        // \todo Remove near-coincident points leveraging sorted array
        // Build final polygon
        for( unsigned int it_v=0; it_v<num_clipped; it_v++ )
            vec_polygon_clipped_pos[it_v] = vec_tmp_pos[it_v];
    }
    return num_clipped;
}

/* Clip a Triangle3 against a Tetrahedron3:
   - Computes a convex polygon resulting from the clipping, with the same orientation as the input triangle.
   - Returns the clipped polygon vertex count [3..6], or 0 if no intersection
   \pre vec_clipped_points MUST have, at least, 6 valid positions
   \todo Consider generalization to Clip_Polygon3_Tetrahedron3 with Tri3->Poly3? or even Clip_Polygon3_Polyedron3?
*/
unsigned int Clip_Triangle3_Tetrahedron3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                          const Vec3& tet0, const Vec3& tet1, const Vec3& tet2, const Vec3& tet3,
                                          Vec3* vec_polygon_clipped_pos,
                                          const Context* p_context )
{
    /*\todo As a Tet is also a Convex Polyhedron, consider clipping
      triangle edges sequentially against arbitrary planes (tet
      faces). This could be useful in other cases, such as Hexaedral
      grids, MeshSurfaces with arbitrary polygons, voronoi cells,
      etc...
      \todo Reconsider it if ever required, but by now Clip_Polygon3_Polyedron3() is BROKEN and won't be fixed.
    */

    Vec3 vec_polygon_pos[3] = { tri0, tri1, tri2 };
    HalfSpace3 vec_polyhedron_hs[4] = { HalfSpace3(mal::Cross(tet2-tet1,tet3-tet1),tet1),
                                        HalfSpace3(mal::Cross(tet0-tet2,tet3-tet2),tet2),
                                        HalfSpace3(mal::Cross(tet0-tet3,tet1-tet3),tet3),
                                        HalfSpace3(mal::Cross(tet2-tet0,tet1-tet0),tet0) };
    Vec3 vec_polyhedron_edge_p0[6] = { tet0, tet0, tet0, tet1, tet2, tet3 };
    Vec3 vec_polyhedron_edge_p1[6] = { tet1, tet2, tet3, tet2, tet3, tet1 };
    unsigned int num_clipped = Clip_Polygon3_Polyhedron3( vec_polygon_pos, 3,
                                                          vec_polyhedron_hs, 4,
                                                          vec_polyhedron_edge_p0, vec_polyhedron_edge_p1, 6,
                                                          vec_polygon_clipped_pos, p_context );
    GEO_ASSERT(num_clipped < 7);
    return num_clipped;
}

#endif //__DISABLED_WHILE_DEVELOPING

}} //namespace geo::np
