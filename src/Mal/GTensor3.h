#ifndef MAL_G_TENSOR3_H
#define MAL_G_TENSOR3_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>

namespace mal
{

template <typename T, int Dim1, int Dim2, int Dim3>
class GTensor3
{    
public:

    typedef T real_type;
    const static int cDimension1 = Dim1;
    const static int cDimension2 = Dim2;
    const static int cDimension3 = Dim3;
    const static int size_in_reals = Dim1*Dim2*Dim3;
    
public:    
    //!\name Construction
    //@{
    finline GTensor3() {}
    //@}
    
    //! \name Accessors
    //@{
    finline const T &operator()( int i, int j, int k ) const { MAL_ASSERT( IsInRangeCO(i,0,Dim1)
                                                                           && IsInRangeCO(j,0,Dim2)
                                                                           && IsInRangeCO(k,0,Dim3));
                                                               return data[i][j][k]; }
    finline T &operator()( int i, int j, int k ) { MAL_ASSERT( IsInRangeCO(i,0,Dim1)
                                                               && IsInRangeCO(j,0,Dim2)
                                                               && IsInRangeCO(k,0,Dim3));
                                                   return data[i][j][k]; }
    //@}

private:
    T data[Dim1][Dim2][Dim3];
};

//---- Misc
template <typename T, int Dim1, int Dim2, int Dim3>
finline bool IsNaN( const GTensor3<T,Dim1,Dim2,Dim3> &t3 )
{
    for( int i=0; i<Dim1; i++ )
        for(int j=0; j<Dim2; j++)
             for(int k=0; k<Dim3; k++)
                 if( IsNaN( t3(i,j,k) ) )
                     return true;
    return false;
}

} //namespace Mal

#endif // MAL_G_TENSOR3_H
