#ifndef MAL_GSRV_H
#define MAL_GSRV_H

#include <pla_types.h>  //!< Required for uint8,uint32 and pla_type_id
#include <Mal/Config.h>
#include <Mal/RealUtils.h>

//#define __USE_TR1_ALIGNMENT_OF
#ifdef __USE_TR1_ALIGNMENT_OF
#  include <tr1/type_traits>
#endif

namespace mal
{

/*! Generic Structured Real Vector Header

  Describes an SRV in just 32b.

  SRV can be used as a lightweight representation for DOF vectors and
  related magnitudes (derivatives, jacobian matrices, ...)
  - reflection: type and structure
  - flat-real vector interpretation/operation (useful for generic algs and SIMD ops)
  - no virtual method overhead
  - dynamic SRT type-checking

  Statically allocated data follows the header, maybe padded using m_Align4.

  \note Support interleaved types? eg: pos,quat ?
  \todo Maybe we don't need a Stride, but we DO need to obtain
        srt::size_in_reals in order to know the flat Real vector
        size!!!
  \todo GSRM: same with a matrix...
  \todo Add arithmetic ops that work on flat real arrays? support SIMD?

  \todo Keep it at 32b, but maybe redistribute bits to allow count > 255?

  We could reduce bits for
  - TypeID may only be: Scalar, Vec2, Vec3, Vec4, Angle, Quat, Vec2_Angle, Vec3_Quat => 8 values, 4b are enough
  - Align4 for maximum of 16B-align can only be 0..3 => 2b are enough
  => Therefore there's room for 32-6 / 2 = 13b for #rows,#cols, 2^13 = 8*1024 > 8000

  Alternatively, we could pad the header size to 16B and enforce
  16B-alignment on it, so that the data is also 16B-aligned, this way
  we'd have "unlimited" 32b counts... and save the align4 field.
*/
template <typename RT>
class GSRV
{
public:
    typedef RT real_type;

public:
    //!\name Attribs. Should fill 32b
    //@{
    uint8 m_TypeID;  //!< element type
    uint8 m_Align4;  //!< 4B alignment from SRV to m_Data[0], required by SIMD
    uint8 m_Count;   //!< #elements \todo In GSRM would be #Rows
    uint8 m_UNUSED;  //32b padding \todo In GSRM would be #Columns
    //@}

    finline GSRV() : m_Count(0), m_TypeID(0), m_Align4(0) {}
    finline GSRV( uint32 type_id, uint32 align4, uint32 count )
    : m_TypeID(type_id), m_Align4(align4), m_Count(count)
        {
            MAL_ASSERT( type_id < 256 && align4 < 256 && count < 256 );
        }
    finline void Init( uint32 type_id, uint32 align4, uint32 count )
        {
            MAL_ASSERT( type_id < 256 && align4 < 256 && count < 256 );
            m_TypeID = type_id;
            m_Align4 = align4;
            m_Count = count;
        }

    //!\name Retrieval of Raw, Real and Structured Real data
    //@{
    // Raw Array
    finline const uint8 *VecRaw() const { return reinterpret_cast<const uint8*>(this) + (m_Align4<<2); }
    finline uint8 *VecRaw() { return reinterpret_cast<uint8*>(this) + (m_Align4<<2); }
    // Real Array
    finline const real_type *VecReal() const { return reinterpret_cast<const real_type*>(VecRaw()); }
    finline real_type *VecReal() { return reinterpret_cast<real_type*>(VecRaw()); }
    // Real Elements
    finline real_type R( int i ) const { return VecReal()[i]; }
    finline real_type &R( int i ) { return VecReal()[i]; }
    // SRT Array
    template <typename SRT>
    finline const SRT *VecSRT() const { MAL_ASSERT( pla_type_id<SRT>::value == m_TypeID ); return reinterpret_cast<const SRT*>( VecRaw() ); }
    template <typename SRT>
    finline SRT *VecSRT() { MAL_ASSERT( pla_type_id<SRT>::value == m_TypeID ); return reinterpret_cast<SRT*>( VecRaw() ); }
    // SRT Elements
    template <typename SRT>
    finline const SRT &SR( int i ) const { return VecSRT<SRT>()[i]; }
    template <typename SRT>
    finline SRT &SR( int i ) { return VecSRT<SRT>()[i]; }
    // Single SRT elements access
    template <typename SRT>
    finline const SRT &SR() const { MAL_ASSERT( m_Count == 1 ); return VecSRT<SRT>()[0]; }
    template <typename SRT>
    finline SRT &SR() { MAL_ASSERT( m_Count == 1 ); return VecSRT<SRT>()[0]; }
    //@}
};

/*! Generic Structured Real Vector Allocation

  Statically allocates a fully typed SRV and fills its GSRV header
*/

template <typename SRT, unsigned D>
class GSRVA: public GSRV<typename SRT::real_type>
{
public:
    typedef typename SRT::real_type real_type;
    typedef GSRV<real_type> base_type;

public:
    finline GSRVA()
    : base_type( pla_type_id<SRT>::value,
                 //sizeof(SRT)/4,
                 (reinterpret_cast<uint8*>(m_Data) - reinterpret_cast<uint8*>(this))/4, //std::tr1::alignment_of<SRT>::value/4,
                 D ) {}

    //\todo could add explicit SRT retrieval methods here... operator[] / ()... etc...

private:
    SRT m_Data[D];
};

/*\todo Dynamically allocated SRV are a bit trickier, because they
  need to ensure CONTIGUOUS allocation of the header and the actual
  data... The only way I can think of is reserving N+1 SRT elements and using the first as the header.
  \todo OR, MAYBE, if the header was at the END of the SRV, no extra element would be required. => NO, this is not possible, the HEADER must be at the start of the object to KNOW the element count!
  \todo ALSO, alignment MUST be taken into account, specially for SIMD (16B-aligment in SSE)
*/
template <typename SRT>
GSRV<typename SRT::real_type> *AllocSRV( unsigned int count )
{
    MAL_ASSERT( sizeof(SRT) >= 4 );
    SRT *pSRT( new SRT[count+1] );
    GSRV<typename SRT::real_type> *pSRV( reinterpret_cast< GSRV<typename SRT::real_type>* >(pSRV) );
    pSRV->Init( pla_type_id<SRT>::value,
                 (reinterpret_cast<uint8*>(&pSRT[1]) - reinterpret_cast<uint8*>(&pSRT[0]))/4,
                 count );
    return pSRV;
}

template <typename SRT>
void DeallocSRV( GSRV<typename SRT::real_type> *p_srv )
{
    SRT *pSRT( reinterpret_cast<SRT*>(p_srv) );
    delete[] pSRT; //count+1 elements
}

} //namespace mal

#endif //MAL_GSRV_H
