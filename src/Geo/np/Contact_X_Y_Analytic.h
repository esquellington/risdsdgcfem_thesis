#ifndef GEO_NP_CONTACT_X_Y_ANALYTIC_H
#define GEO_NP_CONTACT_X_Y_ANALYTIC_H

#include "ContactData.h"
#include "Context.h"

namespace geo {
namespace np {

/*! Sphere Vs Sphere
  - Single contact point
  - Normal lies along the line connecting the centers
  - If coincident, normal = (1,0,0) is returned
*/
template< unsigned D >
inline bool GContact_Sphere_Sphere_Analytic( const mal::GVec<Real,D> &pos1, Real radius1,
                                             const mal::GVec<Real,D> &pos2, Real radius2,
                                             GContactData<D> &cd, GContactCache<D> *p_cc,
                                             const Context *p_context = g_pDefaultContext )
{
    cd.Begin();
    mal::GVec<Real,D> d( pos1 - pos2 );
    Real dist_sq( d.NormSq() );
    if( dist_sq < mal::Sq(radius1+radius2) )
    {
        if( dist_sq > p_context->m_Epsilon_LengthSq )
        {
            Real dist( mal::Sqrt(dist_sq) );
            cd.m_AvgNormal = d / dist;
            cd.m_AvgDepth = radius1 + radius2 - dist;
        }
        else
        {
            cd.m_AvgNormal = mal::GVec<Real,D>::UnitX();
            cd.m_AvgDepth = radius1 + radius2;
        }
        cd.AddCP( pos1 - radius1*cd.m_AvgNormal,
                  pos2 + radius2*cd.m_AvgNormal,
                  cd.m_AvgNormal, cd.m_AvgDepth );
        return cd.End();
    }
    return false;
}

/*! Sphere Vs Plane
  - Single contact point
  - Normal is the plane normal
*/
template< unsigned D >
inline bool GContact_Sphere_Plane_Analytic( const mal::GVec<Real,D> &sphere_pos, Real sphere_radius,
                                            const mal::GVec<Real,D> &plane_normal, Real plane_coeff_d,
                                            GContactData<D> &cd, GContactCache<D> *p_cc,
                                            const Context *p_context = g_pDefaultContext )
{
    cd.Begin();
    Real dist( mal::Dot(sphere_pos,plane_normal) + plane_coeff_d );
    if( dist < sphere_radius )
    {
        cd.m_AvgNormal = plane_normal;
        cd.m_AvgDepth = sphere_radius - dist;
        // Single point
        cd.AddCP( sphere_pos - sphere_radius*cd.m_AvgNormal, //on sphere
                  sphere_pos - (sphere_radius-cd.m_AvgDepth)*cd.m_AvgNormal, //on plane
                  cd.m_AvgNormal, cd.m_AvgDepth );
        return cd.End();
    }
    return false;
}

}} //namespace geo::np

#endif // GEO_NP_CONTACT_X_Y_ANALYTIC_H
