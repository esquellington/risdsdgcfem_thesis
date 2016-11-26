#ifndef MAL_GFXP_H
#define MAL_GFXP_H

#include <Mal/Config.h>

namespace mal
{

template<int RF>
class GFxp
{
public:
    static const unsigned cRadix = RF;
    
public:
    //!\name Construction 
    //@{
    finline GFxp() {}
    finline GFxp( int val ) { m_Bits = (val << RF); }
    finline GFxp( unsigned int val ) { m_Bits = (val << RF); }
    finline GFxp( float val ) { m_Bits = int( val * (1<<RF) ); }
    finline GFxp( double val ) { m_Bits = int( val * (1<<RF) ); }
    //@}

    //!\name Access
    //@{
    finline int GetBits() { return m_Bits; }
    finline void SetBits( int bits ) { m_Bits = bits; }
    //@}

    //!\name Operators
    //@{
    finline GFxp operator+( GFxp val ) const { GFxp sum; sum.m_Bits = (m_Bits + val.m_Bits); return sum; }
    finline GFxp operator-( GFxp val ) const { GFxp sub; sub.m_Bits = (m_Bits - val.m_Bits); return sub; }
    finline GFxp operator*( GFxp val ) const { GFxp prod; prod.m_Bits = (m_Bits*val.m_Bits) >> RF; return prod; }
    finline GFxp operator/( GFxp val ) const { return( GetFloat() / val.GetFloat() ); return *this; } //!\todo SLOW
    
    finline GFxp operator+=( GFxp val ) { m_Bits += val.m_Bits; return *this; }
    finline GFxp operator-=( GFxp val ) { m_Bits -= val.m_Bits; return *this; }        
    finline GFxp operator*=( GFxp val ) { m_Bits = (m_Bits*val.m_Bits) >> RF; return *this; }
    finline GFxp operator/=( GFxp val ) { *this = GFxp( GetFloat() / val.GetFloat() ); return *this; }

    finline GFxp operator-() const { GFxp neg; neg.m_Bits = -m_Bits; return neg; }

    finline bool operator==( GFxp val ) const { return m_Bits == val.m_Bits; }
    finline bool operator!=( GFxp val ) const { return m_Bits != val.m_Bits; }
    //@}

    //!\name Utility
    //@{
    finline float GetFloat() const { return float(m_Bits)/(1<<RF); }
    finline double GetDouble() const { return double(m_Bits)/(1<<RF); }
    //@}
    
private:
    int m_Bits;
};

// ---- Function Templates specialized from RealUtils.h
template<int RF>
inline GFxp<RF> Sqrt( GFxp<RF> value ) { return GFxp<RF>( Sqrt( value.GetFloat() ) ); }

} //namespace mal

#endif
