/* This file contains an implementation of SVD for 3x3 matrix adapted from PhysBAM.

  The pieces required by the original function
Fast_Singular_Value_Decomposition() were scattered accross several
files and classes (MATRIX_3X3, SYMMETRIC_MATRIX_3X3,
SYMMETRIC_MATRIX_2X2, VECTOR_3D), and it was quite impractical to add
them as a dependency for Mal. I claim no ownership on the original or adapted code.

  The original code is part of PhysBAM_Tools and can be found at: http://physbam.stanford.edu/links/getcode.html

  What follows is the Copyright notice from: http://physbam.stanford.edu/links/backhistdisclaimcopy.html

Copyright 1999-2010 Andrew Selle, Andy Lutimirski, Avi Robinson-Mosher, Bridget Vuong, Christopher Allocco, Craig Schroeder, Don Hatch, Douglas Enright, Duc Nguyen, Eftychios Sifakis, Eilene Hao, Elliot English, Eran Guendelman, Fen Zhao, Frank Losasso, Frederic Gibou, Geoffrey Irving, Huamin Wang, Igor Neverov, Jared Go, Jeffrey Hellrung, Jeong-Mo Hong, Jerry Talton, Jiayi Chong, Jonathan Su, Jon Gretarsson, Joseph Teran, Joyce Pan, Justin Solomon, Kevin Der, Mark A. Wicks, Michael Lentine, Michael Turitzin, Mike Rodgers, Neil Molino, Nick Rasmussen, Nipun Kwatra, Paul, James White, Rachel Weinstein, Ranjitha Kumar, Robert Bridson, Robert Travis, Ron Fedkiw, Ryan Kautzman, Sergey Koltakov, Sergey Levine, Silvia Salinas-Blemker, Tamar Shinar, Unnur Gretarsdottir, Wen Zheng, Zhaosheng Bao. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE PHYSBAM PROJECT ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PHYSBAM PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef MAL_GMAT_DECOMPOSITION_SVD3x3_PHYSBAM_H
#define MAL_GMAT_DECOMPOSITION_SVD3x3_PHYSBAM_H

#ifndef MAL_GMAT_DECOMPOSITION_H
#error "You must only include <Mal/GMatDecomposition.h> instead of this implementation file"
#endif

namespace mal
{

namespace FROM_PhysBAM
{

template <typename T>
GVec<T,3> Orthogonal_Vector( const GVec<T,3>& v )
{
    T abs_x( Abs(v.x()) );
    T abs_y( Abs(v.y()) );
    T abs_z( Abs(v.z()) );
    if( abs_x < abs_y )
        return (abs_x < abs_z)
            ? GVec<T,3>(0,v.z(),-v.y())
            : GVec<T,3>(v.y(),-v.x(),0);
    else
        return (abs_y < abs_z)
            ? GVec<T,3>(-v.z(),0,v.x())
            : GVec<T,3>(v.y(),-v.x(),0);
}

template <typename T>
class GSymmetricMat2x2
{
public:
    T x11,x21,x22;

    /* Column-wise layout
    0 2
    1 3
    */
    GSymmetricMat2x2() {}
    GSymmetricMat2x2( const GSymmetricMat2x2 &other )
    : x11(other.x11), x21(other.x21), x22(other.x22)
        {}

    GSymmetricMat2x2( T y11, T y21, T y22 )
    : x11(y11), x21(y21), x22(y22)
        {}

    GSymmetricMat2x2 Cofactor_Matrix() const
        {
            return GSymmetricMat2x2( x22, -x21, x11 );
        }

    GVec<T,2> Largest_Column_Normalized() const
        {
            T sqr11=Sq(x11),sqr12=Sq(x21),sqr22=Sq(x22);
            T scale1=sqr11+sqr12,scale2=sqr12+sqr22;
            if(scale1>scale2) return GVec<T,2>(x11,x21)/Sqrt(scale1);
            else if(scale2>0) return GVec<T,2>(x21,x22)/Sqrt(scale2);
            else return GVec<T,2>(1,0);
        }
};

template <typename T>
GSymmetricMat2x2<T> Transpose_Times_With_Symmetric_Result( const GMat<T,3,2>& A, const GMat<T,3,2>& B ) // A^t*B and assume symmetric result, 9 mults, 6 adds
{
    /* layout
       0 3 6
       1 4 7
       2 5 8
    */
    return GSymmetricMat2x2<T>( A(0,0)*B(0,0) + A(1,0)*B(1,0) + A(2,0)*B(2,0),   //A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2],
                                A(0,1)*B(0,0) + A(1,1)*B(1,0) + A(2,1)*B(2,0),   //A.x[3]*B.x[0]+A.x[4]*B.x[1]+A.x[5]*B.x[2],
                                A(0,1)*B(0,1) + A(1,1)*B(1,1) + A(2,1)*B(2,1) ); //A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5] );
}

template <typename T> class GSymmetricMat3x3;

template <typename T>
GSymmetricMat2x2<T> Conjugate_With_Transpose( const GMat<T,3,2>& A, const GSymmetricMat3x3<T>& B ) // 21 mults, 12 adds
{
    GMat<T,3,3> full_B( B.x11, B.x21, B.x31,
                        B.x21, B.x22, B.x32,
                        B.x31, B.x32, B.x33 );
    return Transpose_Times_With_Symmetric_Result(full_B*A,A);
}

template <typename T>
class GSymmetricMat3x3
{
public:
    T x11,x21,x31,x22,x32,x33;

    /* Column-wise layout
    0 3 6
    1 4 7
    2 5 8
    */

    GSymmetricMat3x3() {}
    GSymmetricMat3x3( const GSymmetricMat3x3 &other )
    : x11(other.x11), x21(other.x21), x31(other.x31), x22(other.x22), x32(other.x32), x33(other.x33)
        {}
    GSymmetricMat3x3( T y11, T y21, T y31, T y22, T y32, T y33 )
    : x11(y11), x21(y21), x31(y31), x22(y22), x32(y32), x33(y33)
        {}

    void From_Normal_Equations( const GMat<T,3,3> &F )
        {
            x11 = Sq(F(0,0)) + Sq(F(1,0)) + Sq(F(2,0));
            x21 = F(0,1)*F(0,0) + F(1,1)*F(1,0) + F(2,1)*F(2,0); //x[3]*x[0]+x[4]*x[1]+x[5]*x[2];
            x31 = F(0,2)*F(0,0) + F(1,2)*F(1,0) + F(2,2)*F(2,0); //x[6]*x[0]+x[7]*x[1]+x[8]*x[2];
            x22 = Sq(F(0,1)) + Sq(F(1,1)) + Sq(F(2,1));
            x32 = F(0,2)*F(0,1) + F(1,2)*F(1,1) + F(2,2)*F(2,1); //x[6]*x[3]+x[7]*x[4]+x[8]*x[5];
            x33 = Sq(F(0,2)) + Sq(F(1,2)) + Sq(F(2,2));
        }

    GSymmetricMat3x3 Cofactor_Matrix() const
        {
            return GSymmetricMat3x3( x22*x33-x32*x32,
                                     x32*x31-x21*x33,
                                     x21*x32-x22*x31,
                                     x11*x33-x31*x31,
                                     x21*x31-x11*x32,
                                     x11*x22-x21*x21 );
        }

    GVec<T,3> Largest_Column_Normalized() const
        {
            T sqr11=Sq(x11),sqr12=Sq(x21),sqr13=Sq(x31),sqr22=Sq(x22),sqr23=Sq(x32),sqr33=Sq(x33);
            T scale1=sqr11+sqr12+sqr13,scale2=sqr12+sqr22+sqr23,scale3=sqr13+sqr23+sqr33;
            if(scale1>scale2){if(scale1>scale3) return GVec<T,3>(x11,x21,x31)/Sqrt(scale1);}
            else if(scale2>scale3) return GVec<T,3>(x21,x22,x32)/Sqrt(scale2);
            if(scale3>0) return GVec<T,3>(x31,x32,x33)/Sqrt(scale3);else return GVec<T,3>(1,0,0);
        }

    /* lambda_x > lambda_y > lambda_z
       reference: Smith, O. "Eigenvalues of a symmetric 3 x 3 matrix". Commun. ACM 4 (4), p. 168, 1961 (thanks, Gene)
    */
    GVec<T,3> Fast_Eigenvalues() const
        {
            const T cOneThird( T(1) / T(3) );
            const T cOneSixth( T(1) / T(6) );
            const T cRootThree( Sqrt(3) );
            T m=(T)cOneThird*(x11+x22+x33);
            T a11=x11-m,a22=x22-m,a33=x33-m,a12_sqr=x21*x21,a13_sqr=x31*x31,a23_sqr=x32*x32;
            T p=(T)cOneSixth*(a11*a11+a22*a22+a33*a33+2*(a12_sqr+a13_sqr+a23_sqr));
            T q=(T).5*(a11*(a22*a33-a23_sqr)-a22*a13_sqr-a33*a12_sqr)+x21*x31*x32;
            T sqrt_p=Sqrt(p),disc=p*p*p-q*q;
            T phi=(T)cOneThird*ATan2(Sqrt(Max((T)0,disc)),q),c=Cos(phi),s=Sin(phi);
            T sqrt_p_cos=sqrt_p*c,root_three_sqrt_p_sin=(T)cRootThree*sqrt_p*s;
            GVec<T,3> lambda( m+2*sqrt_p_cos,
                              m-sqrt_p_cos-root_three_sqrt_p_sin,
                              m-sqrt_p_cos+root_three_sqrt_p_sin );
            // SORT DECREASINGLY, x33 must be the smallest eigenvalue
            if( lambda[0] < lambda[1] )
            {
                if( lambda[1] < lambda[2] ) //2>1>0
                    return GVec<T,3>( lambda[2], lambda[1], lambda[0] );
                else if( lambda[0] < lambda[2] ) //1>2>0
                    return GVec<T,3>( lambda[1], lambda[2], lambda[0] );
                else  //if( lambda[0] >= lambda[2] ) //1>0>2
                    return GVec<T,3>( lambda[1], lambda[0], lambda[2] );
            }
            else //lambda[0] >= lambda[1]
            {
                if( lambda[1] > lambda[2] ) //0>1>2
                    return GVec<T,3>( lambda[0], lambda[1], lambda[2] );
                else if( lambda[0] > lambda[2] ) //0>2>1
                    return GVec<T,3>( lambda[0], lambda[2], lambda[1] );
                else //if( lambda[0] <= lambda[2] ) //2>0>1
                    return GVec<T,3>( lambda[2], lambda[0], lambda[1] );
            }
            //exchange_sort(lambda.x33,lambda.x22,lambda.x11);
            //return lambda;
        }

    GMat<T,3,3> Fast_Eigenvectors( const GSymmetricMat3x3<T>& A, const GVec<T,3> & lambda ) const
    {
        // flip if necessary so that first eigenvalue is the most different
        bool flipped( false );
        GVec<T,3> lambda_flip( lambda );
        if( lambda[0] - lambda[1] < lambda[1] - lambda[2] )
        {
            lambda_flip[0] = lambda[2];
            lambda_flip[2] = lambda[0]; //exchange(lambda_flip.x11,lambda_flip.x33);
            flipped = true;
        }

        // get first eigenvector
        GSymmetricMat3x3<T> A_minus_lambda_flip_0( A );
        A_minus_lambda_flip_0.x11 -= lambda_flip[0];
        A_minus_lambda_flip_0.x22 -= lambda_flip[0];
        A_minus_lambda_flip_0.x33 -= lambda_flip[0];
        GVec<T,3> v1 = A_minus_lambda_flip_0.Cofactor_Matrix().Largest_Column_Normalized(); //VECTOR<T,3> v1=(A-lambda_flip.x11).Cofactor_Matrix().Largest_Column_Normalized(); // 3a + 12m+6a + 9m+6a+1d+1s = 21m+15a+1d+1s

        // form basis for orthogonal complement to v1, and reduce A to this space
        GVec<T,3> v1_orthogonal( Normalized( Orthogonal_Vector(v1) ) ); // VECTOR<T,3> v1_orthogonal=v1.Unit_Orthogonal_Vector(); // 6m+2a+1d+1s (tweak: 5m+1a+1d+1s)
        GMat<T,3,2> other_v;
        GSetColumn<0>( other_v, v1_orthogonal );
        GSetColumn<1>( other_v, Cross( v1, v1_orthogonal ) ); //MATRIX<T,3,2> other_v(v1_orthogonal,VECTOR<T,3>::Cross_Product(v1,v1_orthogonal)); // 6m+3a (tweak: 4m+1a)
        GSymmetricMat2x2<T> A_reduced( Conjugate_With_Transpose(other_v,A) ); //SYMMETRIC_MATRIX<T,2> A_reduced=SYMMETRIC_MATRIX<T,2>::Conjugate_With_Transpose(other_v,A); // 21m+12a (tweak: 18m+9a)

        // find third eigenvector from A_reduced, and fill in second via cross product
        GSymmetricMat2x2<T> A_reduced_minus_lambda_flip_2( A_reduced );
        A_reduced_minus_lambda_flip_2.x11 -= lambda_flip[2];
        A_reduced_minus_lambda_flip_2.x22 -= lambda_flip[2];
        GVec<T,3> v3( other_v * A_reduced_minus_lambda_flip_2.Cofactor_Matrix().Largest_Column_Normalized() ); //VECTOR<T,3> v3=other_v*(A_reduced-lambda_flip.x33).Cofactor_Matrix().Largest_Column_Normalized(); // 6m+3a + 2a + 5m+2a+1d+1s = 11m+7a+1d+1s (tweak: 10m+6a+1d+1s)
        GVec<T,3> v2( Cross(v3,v1) ); //VECTOR<T,3> v2=VECTOR<T,3>::Cross_Product(v3,v1); // 6m+3a

        // finish
        return (flipped)
            ? GMat3x3_From_Columns( v3, v2, -v1 )
            : GMat3x3_From_Columns( v1, v2,  v3 ); //flipped?MATRIX<T,3>(v3,v2,-v1):MATRIX<T,3>(v1,v2,v3);
    }

    void Fast_Solve_Eigenproblem( GVec<T,3> & eigenvalues, GMat<T,3,3>& eigenvectors ) const
    {
        eigenvalues = Fast_Eigenvalues();
        eigenvectors = Fast_Eigenvectors( *this, eigenvalues );
    }

};

template <typename T>
inline void GFastSVD3x3( const GMat<T,3,3> &F,
                         GMat<T,3,3> &U, GVec<T,3> &diag_F, GMat<T,3,3> &Vt )
{
    // decompose normal equations
    GVec<T,3> lambda;
    GSymmetricMat3x3<T> normal_equation_matrix;
    normal_equation_matrix.From_Normal_Equations( F );
    GMat<T,3,3> V;
    normal_equation_matrix.Fast_Solve_Eigenproblem( lambda, V );
    Vt = Transposed(V);

    // compute singular values
    if( lambda[2] < 0 )
    {
        lambda[0] = Max( lambda[0], T(0) );
        lambda[1] = Max( lambda[1], T(0) );
        lambda[2] = Max( lambda[2], T(0) );
    }
    diag_F[0] = Sqrt( lambda[0] );
    diag_F[1] = Sqrt( lambda[1] );
    diag_F[2] = Sqrt( lambda[2] );

    // Negate smallest singular value if F inverted
    if( Det(F) < 0 ) diag_F[2] = -diag_F[2];

    // compute singular vectors
    GSetColumn<0>( U, Normalized(F * GColumn<0>(V)) ); //U.Column(1)=(*this*V.Column(1)).Normalized(); // 15m+8a+1d+1s
    GVec<T,3> v1_orthogonal = Normalized( Orthogonal_Vector( GColumn<0>(U) ) ); //VECTOR<T,3> v1_orthogonal=U.Column(1).Unit_Orthogonal_Vector(); // 6m+2a+1d+1s
    GMat<T,3,2> other_v;
    GSetColumn<0>( other_v, v1_orthogonal );
    GSetColumn<1>( other_v, Cross( GColumn<0>(U), v1_orthogonal ) ); //MATRIX<T,3,2> other_v(v1_orthogonal,VECTOR<T,3>::Cross_Product(U.Column(1),v1_orthogonal)); // 6m+3a
    GSetColumn<1>( U, other_v * ( Transposed(other_v) * Normalized(F * GColumn<1>(V)) ) ); //U.Column(2)=other_v*(other_v.Transpose_Times(*this*V.Column(2))).Normalized(); // 6m+3a + 6m+4a + 9m+6a + 6m+2a+1d+1s = 27m+15a+1d+1s
    GSetColumn<2>( U, Cross( GColumn<0>(U), GColumn<1>(U) ) ); //  U.Column(3)=VECTOR<T,3>::Cross_Product(U.Column(1),U.Column(2)); // 6m+3a
}

} //namespace FROM_PhysBAM

} //namespace mal

#endif //MAL_GMAT_DECOMPOSITION_SVD3x3_PHYSBAM_H

//----------------------------------------------------------------
// What follows is the original FastSVD 3x3 code from PhysBAM,
// scattered accross the classes:
// - MATRIX_3X3
// - SYMMETRIC_MATRIX_3X3
// - SYMMETRIC_MATRIX_2X2
// - VECTOR_3D
//----------------------------------------------------------------

//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Silvia Salinas-Blemker, Craig Schroeder,
//     Andrew Selle, Tamar Shinar, Jerry Talton, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifdef __DISABLED_CODE_PHYSBAM
template<class T>
class SYMMETRIC_MATRIX<T,3>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=3,n=3};

    T x11,x21,x31,x22,x32,x33;
};

SYMMETRIC_MATRIX<T,3> Normal_Equations_Matrix() const // 18 mults, 12 adds
{return SYMMETRIC_MATRIX<T,3>(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],x[3]*x[0]+x[4]*x[1]+x[5]*x[2],x[6]*x[0]+x[7]*x[1]+x[8]*x[2],
                              x[3]*x[3]+x[4]*x[4]+x[5]*x[5],x[6]*x[3]+x[7]*x[4]+x[8]*x[5],x[6]*x[6]+x[7]*x[7]+x[8]*x[8]);}

SYMMETRIC_MATRIX Cofactor_Matrix() // 12 mults, 6 adds
{return SYMMETRIC_MATRIX(x22*x33-x32*x32,x32*x31-x21*x33,x21*x32-x22*x31,x11*x33-x31*x31,x21*x31-x11*x32,x11*x22-x21*x21);}

//#####################################################################
// Function Fast_Eigenvalues
//#####################################################################
// lambda_x > lambda_y > lambda_z
// reference: Smith, O. "Eigenvalues of a symmetric 3 x 3 matrix". Commun. ACM 4 (4), p. 168, 1961 (thanks, Gene)
template<class T> DIAGONAL_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Fast_Eigenvalues() const // 24 mults, 20 adds, 1 atan2, 1 sincos, 2 sqrts
{
    if(!IS_SAME<T,double>::value) return DIAGONAL_MATRIX<T,3>(SYMMETRIC_MATRIX<double,3>(*this).Fast_Eigenvalues());
    // now T is double
    T m=(T)one_third*(x11+x22+x33);
    T a11=x11-m,a22=x22-m,a33=x33-m,a12_sqr=x21*x21,a13_sqr=x31*x31,a23_sqr=x32*x32;
    T p=(T)one_sixth*(a11*a11+a22*a22+a33*a33+2*(a12_sqr+a13_sqr+a23_sqr));
    T q=(T).5*(a11*(a22*a33-a23_sqr)-a22*a13_sqr-a33*a12_sqr)+x21*x31*x32;
    T sqrt_p=sqrt(p),disc=p*p*p-q*q;
    T phi=(T)one_third*atan2(sqrt(max((T)0,disc)),q),c=cos(phi),s=sin(phi);
    T sqrt_p_cos=sqrt_p*c,root_three_sqrt_p_sin=(T)root_three*sqrt_p*s;
    DIAGONAL_MATRIX<T,3> lambda(m+2*sqrt_p_cos,m-sqrt_p_cos-root_three_sqrt_p_sin,m-sqrt_p_cos+root_three_sqrt_p_sin);
    exchange_sort(lambda.x33,lambda.x22,lambda.x11);return lambda;
}
//#####################################################################
// Function Fast_Eigenvectors
//#####################################################################
namespace{
template<class T> MATRIX<T,3>
Fast_Eigenvectors(const SYMMETRIC_MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& lambda) // 71 mults, 44 adds, 3 divs, 3 sqrts
{
    if(!IS_SAME<T,double>::value) PHYSBAM_FATAL_ERROR();
    // T is now always double

    // flip if necessary so that first eigenvalue is the most different
    bool flipped=false;
    DIAGONAL_MATRIX<T,3> lambda_flip(lambda);
    if(lambda.x11-lambda.x22<lambda.x22-lambda.x33){ // 2a
        exchange(lambda_flip.x11,lambda_flip.x33);
        flipped=true;}

    // get first eigenvector
    VECTOR<T,3> v1=(A-lambda_flip.x11).Cofactor_Matrix().Largest_Column_Normalized(); // 3a + 12m+6a + 9m+6a+1d+1s = 21m+15a+1d+1s  //TEMP; A - lambda_flip = SYMM - DIAG => SYMM

    // form basis for orthogonal complement to v1, and reduce A to this space
    VECTOR<T,3> v1_orthogonal=v1.Unit_Orthogonal_Vector(); // 6m+2a+1d+1s (tweak: 5m+1a+1d+1s) //\todo Unit_Orthogonal_Vector
    MATRIX<T,3,2> other_v(v1_orthogonal,VECTOR<T,3>::Cross_Product(v1,v1_orthogonal)); // 6m+3a (tweak: 4m+1a)
    SYMMETRIC_MATRIX<T,2> A_reduced=SYMMETRIC_MATRIX<T,2>::Conjugate_With_Transpose(other_v,A); // 21m+12a (tweak: 18m+9a) //\todo Conjugate_With_Transpose

    // find third eigenvector from A_reduced, and fill in second via cross product
    VECTOR<T,3> v3=other_v*(A_reduced-lambda_flip.x33).Cofactor_Matrix().Largest_Column_Normalized(); // 6m+3a + 2a + 5m+2a+1d+1s = 11m+7a+1d+1s (tweak: 10m+6a+1d+1s) //\todo Largest_Column_Normalized
    VECTOR<T,3> v2=VECTOR<T,3>::Cross_Product(v3,v1); // 6m+3a

    // finish
    return flipped?MATRIX<T,3>(v3,v2,-v1):MATRIX<T,3>(v1,v2,v3);
}
}
//#####################################################################
// Function Fast_Solve_Eigenproblem
//#####################################################################
template<class T> void SYMMETRIC_MATRIX<T,3>::
Fast_Solve_Eigenproblem(DIAGONAL_MATRIX<T,3>& eigenvalues,MATRIX<T,3>& eigenvectors) const // roughly 95 mults, 64 adds, 3 divs, 5 sqrts, 1 atan2, 1 sincos
{
    if(!IS_SAME<T,double>::value){
        DIAGONAL_MATRIX<double,3> eigenvalues_double;MATRIX<double,3> eigenvectors_double;
        SYMMETRIC_MATRIX<double,3>(*this).Fast_Solve_Eigenproblem(eigenvalues_double,eigenvectors_double);
        eigenvalues=DIAGONAL_MATRIX<T,3>(eigenvalues_double);eigenvectors=MATRIX<T,3>(eigenvectors_double);return;}
    // now T is double
    eigenvalues=Fast_Eigenvalues();
    eigenvectors=Fast_Eigenvectors(*this,eigenvalues);
}


template<class T>
class MATRIX<T,3>:public MATRIX_BASE<T,MATRIX<T,3> >
{
public:
    typedef T SCALAR;typedef MATRIX_BASE<T,MATRIX<T,3> > BASE;
    enum WORKAROUND1 {m=3,n=3};
    using BASE::operator*;using BASE::Times_Transpose;using BASE::Transpose_Times;

    T x[9];
};

//#####################################################################
// Function Fast_Singular_Value_Decomposition
//#####################################################################
// U and V rotations, smallest singular value possibly negative
template<class T> void MATRIX<T,3>::
Fast_Singular_Value_Decomposition(MATRIX<T,3>& U,DIAGONAL_MATRIX<T,3>& singular_values,MATRIX<T,3>& V) const // 182 mults, 112 adds, 6 divs, 11 sqrts, 1 atan2, 1 sincos
{
    if(!IS_SAME<T,double>::value){
        MATRIX<double,3> U_double,V_double;DIAGONAL_MATRIX<double,3> singular_values_double;
        MATRIX<double,3>(*this).Fast_Singular_Value_Decomposition(U_double,singular_values_double,V_double);
        U=MATRIX<T,3>(U_double);singular_values=DIAGONAL_MATRIX<T,3>(singular_values_double);V=MATRIX<T,3>(V_double);return;}
    // now T is double

    // decompose normal equations
    DIAGONAL_MATRIX<T,3> lambda;
    Normal_Equations_Matrix().Fast_Solve_Eigenproblem(lambda,V); // 18m+12a + 95m+64a+3d+5s+1atan2+1sincos

    // compute singular values
    if(lambda.x33<0) lambda=lambda.Clamp_Min(0);
    singular_values=lambda.Sqrt(); // 3s
    if(Determinant()<0) singular_values.x33=-singular_values.x33; // 9m+5a

    // compute singular vectors
    U.Column(1)=(*this*V.Column(1)).Normalized(); // 15m+8a+1d+1s
    VECTOR<T,3> v1_orthogonal=U.Column(1).Unit_Orthogonal_Vector(); // 6m+2a+1d+1s
    MATRIX<T,3,2> other_v(v1_orthogonal,VECTOR<T,3>::Cross_Product(U.Column(1),v1_orthogonal)); // 6m+3a
    U.Column(2)=other_v*(other_v.Transpose_Times(*this*V.Column(2))).Normalized(); // 6m+3a + 6m+4a + 9m+6a + 6m+2a+1d+1s = 27m+15a+1d+1s
    U.Column(3)=VECTOR<T,3>::Cross_Product(U.Column(1),U.Column(2)); // 6m+3a
}

#endif //__DISABLED_CODE_PHYSBAM
