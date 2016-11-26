#ifndef S2_MS_FEM_LINEAR_TRIANGLE_2_H
#define S2_MS_FEM_LINEAR_TRIANGLE_2_H

#include "../Config.h"

//#define __ENABLE_TRACE_FEM
#ifdef __ENABLE_TRACE_FEM
#  include <iostream>
#  include <Mal/GSerialization.h>
#endif

namespace S2 {

typedef mal::GMat<Real,3,6> Mat3x6;
typedef mal::GMat<Real,3,3> Mat3x3;
typedef mal::GMat<Real,6,6> Mat6x6;
typedef mal::GVec<Real,6> Vec6;

namespace ms {
namespace fem {

/*Linear 2d Triangle element
  \todo Optimize:
  - K = B^T E B performs lots of multiplications with 0 due to B and E structure, consider computing K_ij explicitly.
*/
class LinearTriangle2
{
public:
    inline LinearTriangle2() : m_K(Real(0)) {}
    inline void Init( uint16 nid0, uint16 nid1, uint16 nid2,
                      const Vec2 &r0, const Vec2 &r1, const Vec2 &r2,
                      Real young_modulus, Real poisson_ratio )
        {
            // Save element params
            m_vec_nid[0] = nid0;
            m_vec_nid[1] = nid1;
            m_vec_nid[2] = nid2;
            m_vec_r[0] = r0;
            m_vec_r[1] = r1;
            m_vec_r[2] = r2;
            m_Barycenter0 = mal::Rcp( Real(3) )*(r0+r1+r2);
            m_YoungModulus = young_modulus;
            m_PoissonRatio = poisson_ratio;
            // Compute element area and stiffness matrix
            m_Area = Compute_K( r0, r1, r2, young_modulus, poisson_ratio, m_K );
            // Compute D_m^-1 for P-K forces
            Compute_InvDm( r0, r1, r2, m_InvDm );
        }
    inline void SetMaterialParams( Real young_modulus, Real poisson_ratio )
        {
            Init( m_vec_nid[0], m_vec_nid[1], m_vec_nid[2],
                  m_vec_r[0], m_vec_r[1], m_vec_r[2],
                  young_modulus, poisson_ratio ); //TEMP: Inefficient, but avoids repeated code...
        }

    inline const uint16 &nid( unsigned int i ) const { return m_vec_nid[i]; }
    inline const Vec2 &r( unsigned int i ) const { return m_vec_r[i]; }
    inline const Vec2 &Barycenter0() const { return m_Barycenter0; }
    inline const Mat6x6 &K() const { return m_K; }
    inline Real Area() const { return m_Area; }
    inline const Mat2x2 &InvDm() const { return m_InvDm; }
    inline void ComputeF( const Vec2 &x0, const Vec2 &x1, const Vec2 &x2, Mat2x2 &F ) const
        {
            // \see FEM.tex , were F = D * Dm^-1 appears as F = Q*P^-1
            Compute_F( x0, x1, x2, m_InvDm, F );
        }

    //inline Transform2 Te( const Vec2 &r1, const Vec2 &r2, const Vec2 &r3 );

public:
    static Real Compute_Area( const Vec2 &r0, const Vec2 &r1, const Vec2 &r2 )
        {
            /*
            return Real(0.5)*( r1.x()*r2.y() + r2.x()*r0.y() + r0.x()*r1.y()
                               - r1.x()*r0.y() - r0.x()*r2.y() - r2.x()*r1.y() );
            */
            /*\note This is the same as 1/2 det Q, where Q = [ r1-r0 r2-r0 ],
              which is FASTER to compute (4 SUB + 2 MUL
              instead of 2 ADD + 3 SUB + 6 MUL)
            */
            return( Real(0.5) *( ( r1.x() - r0.x() ) * ( r2.y() - r0.y() )
                                 - ( r2.x() - r0.x() ) * ( r1.y() - r0.y() ) ) );
        }

    static void Compute_InvDm( const Vec2 &r0, const Vec2 &r1, const Vec2 &r2, Mat2x2 &inv_Dm )
        {
            Mat2x2 Dm;
            Dm(0,0) = r1.x() - r0.x(); Dm(0,1) = r2.x() - r0.x();
            Dm(1,0) = r1.y() - r0.y(); Dm(1,1) = r2.y() - r0.y();
            inv_Dm = mal::Inverse(Dm);
        }

    static void Compute_F( const Vec2 &r0, const Vec2 &r1, const Vec2 &r2, const Mat2x2 &inv_Dm, Mat2x2 &F )
        {
            Mat2x2 Q;
            Q(0,0) = r1.x() - r0.x(); Q(0,1) = r2.x() - r0.x();
            Q(1,0) = r1.y() - r0.y(); Q(1,1) = r2.y() - r0.y();
            F = Q * inv_Dm;
        }

    static Real Compute_B( const Vec2 &r0, const Vec2 &r1, const Vec2 &r2,
                           Mat3x6 &B )
        {
            Real area( Compute_Area(r0,r1,r2) );
            Real rcp_two_area = mal::Rcp( Real(2) * area );
            Real x02( r0.x() - r2.x() ); x02 *= rcp_two_area;
            Real x10( r1.x() - r0.x() ); x10 *= rcp_two_area;
            Real x21( r2.x() - r1.x() ); x21 *= rcp_two_area;
            Real y01( r0.y() - r1.y() ); y01 *= rcp_two_area;
            Real y12( r1.y() - r2.y() ); y12 *= rcp_two_area;
            Real y20( r2.y() - r0.y() ); y20 *= rcp_two_area;
            B(0,0) = y12; B(0,1) =   0; B(0,2) = y20; B(0,3) =   0; B(0,4) = y01; B(0,5) =   0;
            B(1,0) =   0; B(1,1) = x21; B(1,2) =   0; B(1,3) = x02; B(1,4) =   0; B(1,5) = x10;
            B(2,0) = x21; B(2,1) = y12; B(2,2) = x02; B(2,3) = y20; B(2,4) = x10; B(2,5) = y01;
            return area;
        }
    static void Compute_E( Real young_modulus, Real poisson_ratio,
                           Mat3x3 &E )
        {
            Real a( young_modulus / ((1+poisson_ratio)*(1-2*poisson_ratio)) );
            Real a_times_nu( a * poisson_ratio );
            E(0,0) = a - a_times_nu; E(0,1) = a_times_nu;     E(0,2) = Real(0);
            E(1,0) = a_times_nu;     E(1,1) = a - a_times_nu; E(1,2) = Real(0);
            E(2,0) = 0;              E(2,1) = 0;              E(2,2) = Real(0.5)*a - a_times_nu;
        }
    static Real Compute_K( const Vec2 &r0, const Vec2 &r1, const Vec2 &r2,
                           Real young_modulus, Real poisson_ratio,
                           Mat6x6 &K )
        {
            Mat3x6 B;
            Real area = Compute_B( r0, r1, r2, B );
            Mat3x3 E;
            Compute_E(young_modulus, poisson_ratio, E);
            K = area * ( B.Transposed() * (E * B) );
#ifdef __ENABLE_TRACE_FEM
            std::cout << "Area = " << area << std::endl;
            std::cout << "B = " << B << std::endl;
            std::cout << "E = " << E << std::endl;
            std::cout << "K = " << K << std::endl;
#endif
            return area;
        }
private:
    uint16 m_vec_nid[3];
    Vec2 m_vec_r[3];
    Mat6x6 m_K;
    Mat2x2 m_InvDm; //D_m^-1
    Real m_Area;
    Vec2 m_Barycenter0;
    Real m_YoungModulus, m_PoissonRatio;
};

}}} //namespace S2::ms::fem

#endif //S2_MS_FEM_LINEAR_TRIANGLE_2_H
