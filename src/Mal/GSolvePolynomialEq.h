#ifndef MAL_GSOLVEPOLYNOMIALEQ_H
#define MAL_GSOLVEPOLYNOMIALEQ_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>

namespace mal
{

//Solve ax + b = 0 returns x
template<typename T> int GSolvePolynomialEq1( T a, T b, T& x )
{
    if( IsZero(a) ) return 0;
    x = -b / a;
    return 1;
}

/* Solve ax^2 + bx + c = 0 returns x1 <= x2
   x1 = ( -b - sqrt(b^2 - 4ac) ) / 2a
   x2 = ( -b + sqrt(b^2 - 4ac) ) / 2a
*/
template<typename T> int GSolvePolynomialEq2( T a, T b, T c, T& x1, T& x2 )
{
    if( IsZero(a) )
    {
        int num_roots = GSolvePolynomialEq1(b,c,x1);
        x2 = x1;
        return num_roots;
    }
    T discriminant( Sq(b) - 4*a*c );
    if( discriminant < 0 ) return 0;
    T root( Sqrt(discriminant) );
    T inv_divisor( Rcp( T(2)*a ) );
    x1 = (-b - root) * inv_divisor;
    x2 = (-b + root) * inv_divisor;
    if( x1 > x2 ) { T tmp(x1); x1 = x2; x2 = tmp; } //CAN happen if inv_divisor < 0
    return (x1 == x2) ? 1 : 2;
}

/*Solve ax^3 + bx^2 + cx + d = 0 returns x1 <= x2 <= x3
// \todo http://en.wikipedia.org/wiki/Cubic_function#General_formula_of_roots
template<typename T> bool GSolvePolynomialEq3( T a, T b, T c, T d, T& x1, T& x2, T& x3 )
{
    return false;
}
*/

} //namespace mal

#endif
