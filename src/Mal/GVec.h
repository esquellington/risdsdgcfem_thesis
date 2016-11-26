#ifndef MAL_GVEC_H
#define MAL_GVEC_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>

namespace mal
{

template<typename T, int N>
class GVec
{
public:
    typedef T real_type;
    static const unsigned dimension = N;
    static const unsigned size_in_reals = N;

    finline static GVec Zero() { return GVec(T(0)); }
    finline static GVec UnitX() { GVec res(T(0)); res[0] = T(1); return res; }
    finline static GVec UnitY() { MAL_STATIC_ASSERT(N>1); GVec res(T(0)); res[1] = T(1); return res; }
    finline static GVec UnitZ() { MAL_STATIC_ASSERT(N>2); GVec res(T(0)); res[2] = T(1); return res; }

public:
    //!\name Construction
    //@{
    finline GVec() {}

    template <typename S>
    finline explicit GVec( S val )
    { for(int i=0; i<N; i++) data[i] = T(val); } //!< fill with single value

    template <typename S0, typename S1>
    finline explicit GVec( S0 val0, S1 val1 )
    { MAL_STATIC_ASSERT(N==2); data[0] = T(val0); data[1] = T(val1); }

    template <typename S0, typename S1, typename S2>
    finline explicit GVec( S0 val0, S1 val1, S2 val2 )
    { MAL_STATIC_ASSERT(N==3); data[0] = T(val0); data[1] = T(val1); data[2] = T(val2); }

    template <typename S0, typename S1, typename S2, typename S3>
    finline explicit GVec( S0 val0, S1 val1, S2 val2, S3 val3 )
    { MAL_STATIC_ASSERT(N==4); data[0] = T(val0); data[1] = T(val1); data[2] = T(val2); data[3] = T(val3); }

    template <typename S0, typename S1, typename S2, typename S3, typename S4>
    finline explicit GVec( S0 val0, S1 val1, S2 val2, S3 val3, S4 val4 )
    { MAL_STATIC_ASSERT(N==5); data[0] = T(val0); data[1] = T(val1); data[2] = T(val2); data[3] = T(val3); data[4] = T(val4); }

    template <typename S0, typename S1, typename S2, typename S3, typename S4, typename S5>
    finline explicit GVec( S0 val0, S1 val1, S2 val2, S3 val3, S4 val4 , S5 val5 )
    { MAL_STATIC_ASSERT(N==6); data[0] = T(val0); data[1] = T(val1); data[2] = T(val2); data[3] = T(val3); data[4] = T(val4); data[5] = T(val5); }

    template <typename S>
    finline GVec( const GVec<S,N> &v ) { for(int i=0; i<N; i++) data[i] = T(v[i]); }

    template <typename S>
    finline explicit GVec( const S* scalar_array ) { for(int i=0; i<N; i++) data[i] = T(scalar_array[i]); }
    //@}

    //!\name Access
    //@{
    finline const T& operator[]( int i ) const { MAL_ASSERT( IsInRangeCO(i,0,N) ); return data[i]; }
    finline T& operator[]( int i ) { MAL_ASSERT( IsInRangeCO(i,0,N) ); return data[i]; }
    finline const T& x() const { return data[0]; }
    finline const T& y() const { MAL_STATIC_ASSERT(N>1); return data[1]; }
    finline const T& z() const { MAL_STATIC_ASSERT(N>2); return data[2]; }
    finline T& x() { return data[0]; }
    finline T& y() { MAL_STATIC_ASSERT(N>1); return data[1]; }
    finline T& z() { MAL_STATIC_ASSERT(N>2); return data[2]; }
    //@}

    //!\name Add/Sub
    //@{
    finline GVec operator+( const GVec &v ) const { GVec sum; for(int i=0; i<N; i++) sum[i] = data[i] + v[i]; return sum; }
    finline GVec operator-( const GVec &v ) const { GVec sub; for(int i=0; i<N; i++) sub[i] = data[i] - v[i]; return sub; }
    finline GVec operator+=( const GVec &v ) { for(int i=0; i<N; i++) data[i] += v[i]; return *this; }
    finline GVec operator-=( const GVec &v ) { for(int i=0; i<N; i++) data[i] -= v[i]; return *this; }
    //@}

    //!\name Products
    //@{
    finline GVec operator*( T val ) const { GVec prod; for(int i=0; i<N; i++) prod[i] = val*data[i]; return prod; }
    finline GVec operator*=( T val ) { for(int i=0; i<N; i++) data[i] *= val; return *this; }
    friend finline GVec operator*( T val, const GVec &v ) { return v*val; }

    finline GVec operator/( T val ) const { T inv_val = T(1)/val; return (*this) * inv_val; }
    finline GVec operator/=( T val ) { T inv_val = T(1)/val; return (*this)*= inv_val; }

    finline T operator*( const GVec& v ) const { T acc(0); for(int i=0; i<N; i++) acc += data[i] * v[i]; return acc; }
    finline GVec operator%( const GVec& v ) const { MAL_STATIC_ASSERT(N==3); return GVec( data[1]*v[2] - data[2]*v[1],
                                                                                          data[2]*v[0] - data[0]*v[2],
                                                                                          data[0]*v[1] - data[1]*v[0] ); }
    //@}

    //!\name Unary
    //@{
    finline GVec operator-() const { GVec neg; for(int i=0; i<N; i++) neg[i] = -data[i]; return neg; }
    finline GVec Normalized() const { MAL_ASSERT( 0 != Norm() ); return (*this)/Norm(); }
    finline GVec Normalize() { T norm = Norm(); MAL_ASSERT( 0 != norm ); (*this) /= norm; return *this; }
    finline T NormAndNormalize() { T norm = Norm(); MAL_ASSERT( 0 != norm ); (*this) /= norm; return norm; }
    //@}

    //!\name Norms
    //@{
    finline T NormSq() const { return ( (*this) * (*this) ); }
    finline T Norm() const { return Sqrt( NormSq() ); }
    //@}

    //!\name Utility
    //@{
    finline T* AsArray() { return &data[0]; }
    finline const T* AsArray() const { return &data[0]; }
    finline void ToArray( T* real_array ) const { for(int i=0; i<N; i++) real_array[i] = data[i]; }
    finline void FromArray( const T* real_array ) { for(int i=0; i<N; i++) data[i] = real_array[i]; }
    finline GVec Lerp( const GVec& p, T lambda ) const { return (lambda*(*this) + (T(1)-lambda)*p);}
    //@}

    //\name Comparison
    //@{
    finline bool operator==( const GVec& v ) const { for(int i=0; i<N; i++)
                                                         if( data[i] != v[i] )
                                                             return false;
                                                     return true; }
    finline bool operator!=( const GVec& v ) const { return !(*this == v ); }
    //@}

private:
    T data[N];
};

//! Dot product
template < typename T, int N >
finline T Dot( const GVec<T,N> &vec0, const GVec<T,N> &vec1 ) { return vec0 * vec1; }

//! Normalization
template < typename T, int N >
finline GVec<T,N> Normalized( const GVec<T,N> &vec ) { return vec.Normalized(); }

//! Safe Normalization
template < typename T, int N >
finline GVec<T,N> SafeNormalized( const GVec<T,N> &vec, T epsilon = Epsilon<T>() )
{
    if( vec.NormSq() > Sq(1e-6) ) return vec.Normalized();
    else if( sizeof(T) < sizeof(double) ) //retry with higher-precision Real...
    {
        GVec<double,N> vec_d( vec );
        if( vec_d.NormSq() > Sq(1e-9) ) return GVec<T,N>(vec_d.Normalized());
        else return GVec<T,N>(0);
    }
    else return GVec<T,N>(0);
}

//! Norm
template < typename T, int N >
finline T Norm( const GVec<T,N> &vec ) { return vec.Norm(); }

//! NormSq
template < typename T, int N >
finline T NormSq( const GVec<T,N> &vec ) { return vec.NormSq(); }

//---- Misc
template < typename T, int N >
finline GVec<T,N> Clamp( const GVec<T,N> &vec, const GVec<T,N> &vec_min, const GVec<T,N> &vec_max )
{
    GVec<T,N> res;
    for( unsigned int i=0; i<N; i++ ) res[i] = Clamp( vec[i], vec_min[i], vec_max[i] );
    return res;
}

template < typename T, int N >
finline T Min( const GVec<T,N> &vec )
{
    T res( vec[0] );
    for( unsigned int i=1; i<N; i++ ) if( res > vec[i] ) res = vec[i];
    return res;
}

template < typename T, int N >
finline T Max( const GVec<T,N> &vec )
{
    T res( vec[0] );
    for( unsigned int i=1; i<N; i++ ) if( res < vec[i] ) res = vec[i];
    return res;
}

template < typename T, int N >
finline T Sum( const GVec<T,N> &vec )
{
    T res(0);
    for( unsigned int i=0; i<N; i++ ) res += vec[i];
    return res;
}

template < typename T, int N >
finline bool IsZero( const GVec<T,N> &vec, T epsilon = Epsilon<T>() )
{
    T acc(0);
    for( unsigned int i=0; i<N; i++ )
        acc += Sq(vec[i]);
    return acc < Sq(epsilon);
}

template < typename T, int N >
finline bool IsNaN( const GVec<T,N> &vec )
{
    for( unsigned int i=1; i<N; i++ )
        if( IsNaN(vec[i]) )
            return true;
    return false;
}

//---- Size-specific vector ops

//! Cross product
template < typename T >
finline GVec<T,3> Cross( const GVec<T,3> &vec0, const GVec<T,3> &vec1 ) { return vec0 % vec1; }

//! Scalar triple product v0*(v1xv2)
template < typename T >
finline T ScalarTriple( const GVec<T,3> &vec0, const GVec<T,3> &vec1, const GVec<T,3> &vec2 ) { return Dot(vec0,Cross(vec1,vec2)); }

//! 2d perpendicular (CW orientation)
template < typename T >
finline GVec<T,2> PerpendicularCW( const GVec<T,2> &vec )
{
    return GVec<T,2>( -vec[1], vec[0] );
}
//! 2d perpendicular (CCW orientation)
template < typename T >
finline GVec<T,2> PerpendicularCCW( const GVec<T,2> &vec )
{
    return GVec<T,2>( vec[1], -vec[0] );
}

//----------------------------------------------------------------
//\todo SHOULD BE IN GVecUtils.h
template < int I0, int I1, typename T, int N >
inline GVec<T,I1-I0+1> GRange( const GVec<T,N> &v )
{
    MAL_STATIC_ASSERT( I0 >= 0 && I0 < I1 && I1 < N );
    GVec<T,I1-I0+1> r;
    for(int i=I0; i<=I1; i++) r[i-I0] = v[i];
    return r;
}
template < int I0, int I1, typename T, int N >
inline void GSetRange( GVec<T,N> &dst, const GVec<T,I1-I0+1> &src )
{
    MAL_STATIC_ASSERT( I0 >= 0 && I0 < I1 && I1 < N );
    for(int i=I0; i<=I1; i++) dst[i] = src[i-I0];
}
template < typename T, int N1, int N2 >
inline GVec<T,N1+N2> Concat( const GVec<T,N1>& v1, const GVec<T,N1> v2 )
{
    GVec<T,N1+N2> r;
    for(int i=0; i<N1; i++) r[i] = v1[i];
    for(int i=0; i<N2; i++) r[N1+i] = v2[i];
    return r;
}
template < typename S, typename T, int N >
inline GVec<T,N+1> Concat( const S& s, const GVec<T,N> v )
{
    GVec<T,N+1> r;
    r[0] = T(s);
    for(int i=0; i<N; i++) r[1+i] = v[i];
    return r;
}
template < typename T, int N, typename S >
inline GVec<T,N+1> Concat( const GVec<T,N> v, const S& s )
{
    GVec<T,N+1> r;
    for(int i=0; i<N; i++) r[i] = v[i];
    r[N] = T(s);
    return r;
}
//\todo SHOULD BE IN GVecUtils.h
//----------------------------------------------------------------

} //namespace mal

#endif
