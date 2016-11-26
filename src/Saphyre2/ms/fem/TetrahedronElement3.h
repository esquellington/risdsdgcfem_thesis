#ifndef S2_MS_FEM_TETRAHEDRON_ELEMENT_3_H
#define S2_MS_FEM_TETRAHEDRON_ELEMENT_3_H

#include "../Config.h"

//#define __ENABLE_TRACE_FEM
#ifdef __ENABLE_TRACE_FEM
#  include <iostream>
#  include <Mal/GSerialization.h>
#endif

namespace S2 {

typedef mal::GMat<Real,6,12> Mat6x12;
typedef mal::GMat<Real,6,6> Mat6x6;
typedef mal::GMat<Real,12,12> Mat12x12;

namespace ms {
namespace fem {

/*Linear 3d Tetrahedron element
  \todo Optimize:
  - K = B^T E B performs lots of multiplications with 0 due to B and E structure, consider computing K_ij explicitly.
*/
class TetrahedronElement3
{
public:
    struct DoC //Possible cases are VF => vit = 00,11,22,33 (vit0==vit1) and EE => vit = 01,02,03,12,13,23 (vit0<vit1)
    {
        inline DoC() : m_VIT0(3), m_VIT1(0) {} //Invalid
        inline explicit DoC( int vit ) : m_VIT0(vit), m_VIT1(vit) {} //V-F
        inline DoC( int vit0, int vit1 ) : m_VIT0(vit0), m_VIT1(vit1) {} //E-E
        inline bool IsValid() const { return m_VIT0 <= m_VIT1; }
        inline bool IsInvalid() const { return m_VIT0 > m_VIT1; }
        inline bool IsVF() const { return m_VIT0 == m_VIT1; }
        inline bool IsEE() const { return m_VIT0 < m_VIT1; }
        uint8 m_VIT0, m_VIT1; //\todo Actually, 2 bits per vit would be enough (0..3)
    };

public:
    inline TetrahedronElement3() : m_K(Real(0)) {}
    inline void Init( uint16 nid0, uint16 nid1, uint16 nid2, uint16 nid3,
                      const Vec3 &r0, const Vec3 &r1, const Vec3 &r2, const Vec3 &r3,
                      Real young_modulus, Real poisson_ratio )
        {
            // Save element params
            m_vec_nid[0] = nid0;
            m_vec_nid[1] = nid1;
            m_vec_nid[2] = nid2;
            m_vec_nid[3] = nid3;
            m_vec_r[0] = r0;
            m_vec_r[1] = r1;
            m_vec_r[2] = r2;
            m_vec_r[3] = r3;
            m_Barycenter0 = mal::Rcp( Real(4) )*(r0+r1+r2+r3);
            m_YoungModulus = young_modulus;
            m_PoissonRatio = poisson_ratio;
            // Compute element area and stiffness matrix
            m_Volume = Compute_K( r0, r1, r2, r3, young_modulus, poisson_ratio, m_K );
            // Compute D_m^-1 for P-K forces
            Mat3x3 Dm;
            Compute_D( r0, r1, r2, r3, Dm );
            m_InvDm = mal::Inverse( Dm );
        }
    inline void SetMaterialParams( Real young_modulus, Real poisson_ratio )
        {
            Init( m_vec_nid[0], m_vec_nid[1], m_vec_nid[2], m_vec_nid[3],
                  m_vec_r[0], m_vec_r[1], m_vec_r[2], m_vec_r[3],
                  young_modulus, poisson_ratio ); //TEMP: Inefficient, but avoids repeated code...
        }

    inline const uint16 &nid( unsigned int i ) const { return m_vec_nid[i]; }
    inline const Vec3 &r( unsigned int i ) const { return m_vec_r[i]; }
    inline const Vec3 &Barycenter0() const { return m_Barycenter0; }
    inline const Mat12x12 &K() const { return m_K; }
    inline Real Volume() const { return m_Volume; }
    inline const Mat3x3 &InvDm() const { return m_InvDm; }
    inline void ComputeF( const Vec3 &x0, const Vec3 &x1, const Vec3 &x2, const Vec3 &x3, Mat3x3 &F ) const
        {
            // \see FEM.tex , were F = D * Dm^-1 appears as F = Q*P^-1
            Compute_F( x0, x1, x2, x3, m_InvDm, F );
        }

public:
    static Real Compute_Volume( const Vec3 &r0, const Vec3 &r1, const Vec3 &r2, const Vec3 &r3 )
        {
            Mat3x3 D;
            Compute_D(r0,r1,r2,r3,D);
            return mal::Det(D) / Real(6);
        }

    static void Compute_D( const Vec3 &r0, const Vec3 &r1, const Vec3 &r2, const Vec3 &r3,
                           Mat3x3 &D )
        {
            D = mal::GMat3x3_From_Columns( r1-r0, r2-r0, r3-r0 );
        }

    static void Compute_F( const Vec3 &r0, const Vec3 &r1, const Vec3 &r2, const Vec3 &r3, const Mat3x3 &inv_Dm,
                           Mat3x3 &F )
        {
            Mat3x3 Dm;
            Compute_D( r0, r1, r2, r3, Dm );
            F = Dm * inv_Dm;
        }

    static Real Compute_B( const Vec3 &r0, const Vec3 &r1, const Vec3 &r2, const Vec3 &r3,
                           Mat6x12 &B )
        {
            Real volume = Compute_Volume( r0, r1, r2, r3 );
            MS_ASSERT( volume > Real(0) );
            //invJ from AFEM-16 p.7
            mal::GMat<Real,4,4> invJ;
            Real x1(r0.x()),x2(r1.x()),x3(r2.x()),x4(r3.x());
            Real y1(r0.y()),y2(r1.y()),y3(r2.y()),y4(r3.y());
            Real z1(r0.z()),z2(r1.z()),z3(r2.z()),z4(r3.z());
            invJ(0,0) = x2*(y3*z4 - y4*z3) + x3*(y4*z2-y2*z4) + x4*(y2*z3 - y3*z2);
            invJ(1,0) = x1*(y4*z3 - y3*z4) + x3*(y1*z4-y4*z1) + x4*(y3*z1 - y1*z3);
            invJ(2,0) = x1*(y2*z4 - y4*z2) + x2*(y4*z1-y1*z4) + x4*(y1*z2 - y2*z1);
            invJ(3,0) = x1*(y3*z2 - y2*z3) + x2*(y1*z3-y3*z1) + x3*(y2*z1 - y1*z2);
            invJ(0,1) = (y4-y2)*(z3-z2) - (y3-y2)*(z4-z2); invJ(0,2) = (x3-x2)*(z4-z2) - (x4-x2)*(z3-z2); invJ(0,3) = (x4-x2)*(y3-y2) - (x3-x2)*(y4-y2);
            invJ(1,1) = (y3-y1)*(z4-z3) - (y3-y4)*(z1-z3); invJ(1,2) = (x4-x3)*(z3-z1) - (x1-x3)*(z3-z4); invJ(1,3) = (x3-x1)*(y4-y3) - (x3-x4)*(y1-y3);
            invJ(2,1) = (y2-y4)*(z1-z4) - (y1-y4)*(z2-z4); invJ(2,2) = (x1-x4)*(z2-z4) - (x2-x4)*(z1-z4); invJ(2,3) = (x2-x4)*(y1-y4) - (x1-x4)*(y2-y4);
            invJ(3,1) = (y1-y3)*(z2-z1) - (y1-y2)*(z3-z1); invJ(3,2) = (x2-x1)*(z1-z3) - (x3-x1)*(z1-z2); invJ(3,3) = (x1-x3)*(y2-y1) - (x1-x2)*(y3-y1);
            invJ = mal::Rcp( Real(6) * volume ) * invJ;

            Real a1(invJ(0,1)), b1(invJ(0,2)), c1(invJ(0,3));
            Real a2(invJ(1,1)), b2(invJ(1,2)), c2(invJ(1,3));
            Real a3(invJ(2,1)), b3(invJ(2,2)), c3(invJ(2,3));
            Real a4(invJ(3,1)), b4(invJ(3,2)), c4(invJ(3,3));

            //B from AFEM-16 p.13
            B(0,0) = a1; B(0,1) =  0; B(0,2) =  0; B(0,3) = a2; B(0,4) =  0; B(0,5) =  0; B(0,6) = a3; B(0,7) =  0; B(0,8) =  0; B(0,9) = a4; B(0,10) =  0; B(0,11) =  0;
            B(1,0) =  0; B(1,1) = b1; B(1,2) =  0; B(1,3) =  0; B(1,4) = b2; B(1,5) =  0; B(1,6) =  0; B(1,7) = b3; B(1,8) =  0; B(1,9) =  0; B(1,10) = b4; B(1,11) =  0;
            B(2,0) =  0; B(2,1) =  0; B(2,2) = c1; B(2,3) =  0; B(2,4) =  0; B(2,5) = c2; B(2,6) =  0; B(2,7) =  0; B(2,8) = c3; B(2,9) =  0; B(2,10) =  0; B(2,11) = c4;
            B(3,0) = b1; B(3,1) = a1; B(3,2) =  0; B(3,3) = b2; B(3,4) = a2; B(3,5) =  0; B(3,6) = b3; B(3,7) = a3; B(3,8) =  0; B(3,9) = b4; B(3,10) = a4; B(3,11) =  0;
            B(4,0) =  0; B(4,1) = c1; B(4,2) = b1; B(4,3) =  0; B(4,4) = c2; B(4,5) = b2; B(4,6) =  0; B(4,7) = c3; B(4,8) = b3; B(4,9) =  0; B(4,10) = c4; B(4,11) = b4;
            B(5,0) = c1; B(5,1) =  0; B(5,2) = a1; B(5,3) = c2; B(5,4) =  0; B(5,5) = a2; B(5,6) = c3; B(5,7) =  0; B(5,8) = a3; B(5,9) = c4; B(5,10) =  0; B(5,11) = a4;

            return volume;
        }
    static void Compute_E( Real young_modulus, Real poisson_ratio,
                           Mat6x6 &E )
        {
            //E from AFEM-16 p.14
            E = Mat6x6::Zero();
            Real a( young_modulus / ((Real(1)+poisson_ratio)*(Real(1)-Real(2)*poisson_ratio)) );
            Real a_nu( a*poisson_ratio );
            Real a_k1( a*(Real(1) - poisson_ratio) );
            Real a_k2( a*(Real(0.5) - poisson_ratio) );
            E(0,0) = a_k1; E(0,1) = a_nu; E(0,2) = a_nu;
            E(1,0) = a_nu; E(1,1) = a_k1; E(1,2) = a_nu;
            E(2,0) = a_nu; E(2,1) = a_nu; E(2,2) = a_k1;
            E(3,3) = a_k2;
            E(4,4) = a_k2;
            E(5,5) = a_k2;
        }
    static Real Compute_K( const Vec3 &r0, const Vec3 &r1, const Vec3 &r2, const Vec3 &r3,
                           Real young_modulus, Real poisson_ratio,
                           Mat12x12 &K )
        {
            Mat6x12 B;
            Real volume = Compute_B( r0, r1, r2, r3, B );
            Mat6x6 E;
            Compute_E (young_modulus, poisson_ratio, E );
            K = volume * ( B.Transposed() * (E * B) );
#ifdef __ENABLE_TRACE_FEM
            std::cout << "Volume = " << area << std::endl;
            std::cout << "B = " << B << std::endl;
            std::cout << "E = " << E << std::endl;
            std::cout << "K = " << K << std::endl;
#endif
            return volume;
        }
private:
    uint16 m_vec_nid[4];
    Vec3 m_vec_r[4];
    Mat12x12 m_K;
    Mat3x3 m_InvDm; //D_m^-1
    Real m_Volume;
    Vec3 m_Barycenter0;
    Real m_YoungModulus, m_PoissonRatio;
};

}}} //namespace S2::ms::fem

#endif //S2_MS_FEM_TETRAHEDRON_ELEMENT_3_H
