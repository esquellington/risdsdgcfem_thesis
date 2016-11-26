#ifndef UTIL_ID32_H
#define UTIL_ID32_H

#include <util/Config.h>
#include <string>
#include <string.h>

namespace util
{

/*! 32-bit Id that store a 6-char identifier with an optional ending digit.
- Alphabet = 6 lower-case letters in [a-z] + 5 symbols { _&[]. }
- Optional digit = '#' + [0-3]
\note: Id("foo") == Id("foo#0")
\note: App must call Id32::Init() to init internal tables.
*/
class Id32
{
private:
    static bool sIsInitialized;
    static uint8 sC2B[255];  //!< Contains B5 representation for chars, and B2 rep for digits [0-3]
    static char sB2C[32];    //!< Inverse map from B5 to char
    static const unsigned int cMaxStrLength = 8; //6 chars + '#' + 1 digit [0,1,2,3]
    static const unsigned int cInvalidId32 = 0;
    static const unsigned int cBitsInvalid = 0;            //5 lower bits = 0
    static const unsigned int cBitsSeparator = 0x000000FF; //8 lower bits = 1
    static const unsigned int cMaskLower5 = 0x0000001F;    //5 lower bits = 1
    static const unsigned int cMaskLower2 = 0x00000003;    //2 lower bits = 1
    
public:
    static void Init(); //!< Must be called to init tables
    
    finline Id32() : m_Bits(cInvalidId32) {}
    finline Id32( const Id32 &id ) { m_Bits = id.m_Bits; }
    finline explicit Id32( const char *str ) { Set(str); }
    finline explicit Id32( const std::string &str ) { Set(str.c_str()); }

    void Set( const char *str, unsigned int length = 0 )
    {
        UTIL_ASSERT( sIsInitialized );
        
        m_Bits = 0;
        unsigned int str_length = (length>0) ? length : strlen(str);
        UTIL_ASSERT( str_length > 0 && str_length < cMaxStrLength );
        for( unsigned int i=0; i<str_length; i++ )
        {
            uint8 b5 = sC2B[ (int)str[i] ];
            switch(b5)
            {
            case cBitsInvalid:   //unknown character, terminate
                m_Bits = cInvalidId32;
                break;
            case cBitsSeparator: //ending digits separator, skip it, set upper bits from digit [0-3] and force exit
                {                    
                    UTIL_ASSERT( i == str_length-2 );
                    uint8 b2 = sC2B[ (int)str[i+1] ];
                    UTIL_ASSERT( b2 < 4 );
                    m_Bits |= (b2 << 30);
                    i+=2;
                }
                break;
            default:  //known character, set bits and continue
                m_Bits |= (b5 << (5*i));
                break;
            }
        }
    }

    finline bool IsValid() const { return cInvalidId32 != m_Bits; }

    void ToStr( char *str ) const
    {
        UTIL_ASSERT( sIsInitialized ) ;
        
        if( IsValid() )
        {
            // print up to 5 valid
            int i=0;
            uint32 bits = m_Bits;
            uint8 b5 = bits & cMaskLower5;
            while( b5 != 0 && i < 6 )
            {
                str[i++] = sB2C[b5];
                bits = bits >> 5;
                b5 = bits & cMaskLower5;
            }
            // print digits if != 0 and ending mark after
            uint8 b2 = (m_Bits >> 30) & cMaskLower2;
            if( 0 == b2 )
                str[i] = 0;
            else
            {
                str[i] = '#';
                str[i+1] = '0' + b2;
                str[i+2] = 0;
            }
        }
        else
        {
            str[0]='E';str[1]='R';str[2]='R';str[3]='O';str[4]='R';str[5]=0;
        }
    }
    
private:
    uint32 m_Bits;
};

} //namespace util

#endif //UTIL_ID32_H
