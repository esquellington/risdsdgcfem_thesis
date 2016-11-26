#ifndef S2_MS_SPH_KERNELS_H
#define S2_MS_SPH_KERNELS_H

#include "Config.h"

namespace S2 { namespace ms {


//---- DCEW96 Kernels -------------------------------------------

/*! 2D Poly3 kernel
  
  \param r = \| r1-r0 \|
  \param pow_h = [h^0,h^1,...,h^9]
  
  \pre r >= 0 and h >= 0

  if r<h
     W(r,h) = k_{norm} * (1-r/h)^3
  else
     W(r,h) = 0

  where k_{norm} is a normalization factor that ensures that the
  kernel adds up to 1 inside a disk with radius h.  
     k_{norm} = \frac{1}{ int_{disk(h)}{W(r,h)} }
     
  integrating analitically we obtain:
     k_{norm} = 10/(Pi*h^2)
     
  \sa http://integrals.wolfram.com/index.jsp?expr=%281-x%2Fh%29%5E3+*+2+*+Pi+*+x
*/
inline Real Kernel_DCEW96_2D( Real r, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return 0;
    Real tmp = Real(1.0f) - (r/pow_h[1]);
    return (10.0f/(MAL_CONSTANT_PI*pow_h[2])) * tmp*tmp*tmp;
}


/*! dW/dr

  \param r = \| r1-r0 \|
  \param pow_h = [h^0,h^1,...,h^9]
  
  \pre r >= 0 and h >= 0

  if r<h
     W(r,h) = k (1-r/h)^2
  else
     W(r,h) = 0

  where k gathers all constants, including the kernel normalization
  constant k_{normal}
     k = -30/(Pi*h^3)
*/
inline Real dWdr_Kernel_DCEW96_2D( Real r, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return 0;
    Real tmp = Real(1.0f) - (r/pow_h[1]);
    return (-30.0f/(MAL_CONSTANT_PI*pow_h[3])) * tmp*tmp;
}


//---- Keiser Kernels -------------------------------------------

/*! 2D Poly6 kernel
  
  \param r = \| r1-r0 \|
  \param pow_h = [h^0,h^1,...,h^9]
  
  \pre r >= 0 and h >= 0

  if r<h
     W(r,h) = k_{norm} * (1-r^2/h^2)^3
  else
     W(r,h) = 0

  where
     k_{norm} = ??
*/
inline Real Kernel_Poly6_2D( Real r, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return 0;
    Real tmp = Real(1.0f) - (r*r/pow_h[2]);
    return (4.0f/(MAL_CONSTANT_PI*pow_h[2])) * tmp*tmp*tmp;
}

/*! Gradient of Poly6 Kernel in 2D

  \param r = \| r1-r0 \|
  \param u = (r1-r0)/r
  \param pow_h = [h^0,h^1,...,h^9]

  \pre r >= 0 and h >= 0
  
  if r<h
     \nabla W(r1-r0,h) = \frac{-24*r}{\pi h^4} * (1-r^2/h^2)^2 * u
  else
     \nabla W(r,h) = 0
  
  \sa "Thesis Keiser"
*/
inline Vec2 GradientKernel_Poly6_2D( Real r, const Vec2 &u, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return Vec2(0,0);
    Real tmp =  1 - r*r/pow_h[2];
    return ( (-24.0f*r/(MAL_CONSTANT_PI*pow_h[4])) * tmp*tmp )* u;
}

/*! Spiky Kernel
  
  \param r = \| r1-r0 \|
  \param pow_h = [h^0,h^1,...,h^9]
  
  \pre r >= 0 and h >= 0
  
  if r<h
     W(r,h) = \frac{10}{\pi h^5} (h-r)^3
  else
     W(r,h) = 0

  \sa "Thesis Keiser"
*/
inline Real Kernel_Spiky_2D( Real r, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return 0;
    Real tmp = pow_h[1] - r;
    return (10.0f/(MAL_CONSTANT_PI*pow_h[5])) * tmp*tmp*tmp;
}

/*! Gradient of Spiky Kernel

  \param r = \| r1-r0 \|
  \param u = (r1-r0)/r
  \param pow_h = [h^0,h^1,...,h^9]

  \pre r >= 0 and h >= 0
  
  if r<h
     \nabla W(r1-r0,h) = \frac{-30}{\pi h^5} * (h-r)^2 * u
  else
     \nabla W(r,h) = 0
  
  \sa "Thesis Keiser"
*/
inline Vec2 GradientKernel_Spiky_2D( Real r, const Vec2 &u, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return Vec2(0,0);
    Real tmp =  pow_h[1] - r;
    return ( (-30.0f/(MAL_CONSTANT_PI*pow_h[5])) * tmp*tmp )* u;
}


/*! Viscosity Kernel (positive laplacian everywhere...)
  
  \param r = \| r1-r0 \|
  \param pow_h = [h^0,h^1,...,h^9]
  
  \pre r >= 0 and h >= 0
    
  if r<h
     W(r,h) = \frac{-30}{11\pi h^2}
     ( -\frac{r^3}{2h^3) + \frac{r^2}{h^2} + \frac{r}{2h} - 1)
  else
     W(r,h) = 0
     
  \sa "Thesis Keiser"
*/
inline Real Kernel_Viscosity_Keiser_2D( Real r, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return 0;
    return (-30.0f/(11.0f*MAL_CONSTANT_PI*pow_h[2]))
        * ( -r*r*r/(2*pow_h[3]) + r*r/pow_h[2] + r/(2*pow_h[1]) - 1 );
}


/*! Laplacian of the Viscosity Kernel (positive everywhere...)
  
  \param r = \| r1-r0 \|
  \param pow_h = [h^0,h^1,...,h^9]
  
  \pre r >= 0 and h >= 0
    
  if r<h
     \nabla^2 W(r,h) = \frac{-30}{11\pi h^2} ( -\frac{3r}{h^3) + \frac{2}{h^2} )
  else
     \nabla^2 W(r,h) = 0
     
  \sa "Thesis Keiser"
*/
inline Real LaplacianKernel_Viscosity_Keiser_2D( Real r, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return 0;
    return (-30.0f/(11.0f*MAL_CONSTANT_PI*pow_h[2])) * ( -3.0f*r/pow_h[3] + 2/pow_h[2] );
}


//---- Clavet Kernels -------------------------------------------

/*! 2D Poly2 kernel
  
  \param r = \| r1-r0 \|
  \param pow_h = [h^0,h^1,...,h^9]
  
  \pre r >= 0 and h >= 0

  if r<h
     W(r,h) = k_{norm} * (1-r/h)^2
  else
     W(r,h) = 0

  where k_{norm} is a normalization factor that ensures that the
  kernel adds up to 1 inside a disk with radius h.  
     k_{norm} = \frac{1}{ int_{disk(h)}{W(r,h)} }
     
  integrating analitically we obtain:
     k_{norm} = 6/(Pi*h^2)
     
  \sa http://integrals.wolfram.com/index.jsp?expr=%281-x%2Fh%29%5E2+*+%282+*+Pi+*+x%29
*/
inline Real Kernel_Clavet_Poly2( Real r, const Real *pow_h )
{
    MS_ASSERT( r >= 0 );
    if( r > pow_h[1] ) return 0;
    Real tmp = Real(1.0f) - (r/pow_h[1]);
    return (6.0f/(MAL_CONSTANT_PI*pow_h[2])) * tmp*tmp;
}

} } // namespace S2::ms

#endif // S2_MS_SPH_KERNELS_H
