#include "Id32.h"

bool util::Id32::sIsInitialized(false);
uint8 util::Id32::sC2B[255];
char util::Id32::sB2C[32];

namespace util
{

void Id32::Init()
{
    UTIL_ASSERT( !sIsInitialized ) ;
    
    // Most chars are invalid...
    for( int i=0; i<255; i++ ) sC2B[i] = cBitsInvalid;
    // 26 chars [a-z] mapped to [1,26]
    sC2B['a'] = 1; sC2B['b'] = 2; sC2B['c'] = 3; sC2B['d'] = 4; sC2B['e'] = 5;
    sC2B['f'] = 6; sC2B['g'] = 7; sC2B['h'] = 8; sC2B['i'] = 9; sC2B['j'] =10;
    sC2B['k'] =11; sC2B['l'] =12; sC2B['m'] =13; sC2B['n'] =14; sC2B['o'] =15;
    sC2B['p'] =16; sC2B['q'] =17; sC2B['r'] =18; sC2B['s'] =19; sC2B['t'] =20;
    sC2B['u'] =21; sC2B['v'] =22; sC2B['w'] =23; sC2B['x'] =24; sC2B['y'] =25; sC2B['z'] =26;
    // 5 simbols { _&[]. } mapped to [27-31]
    sC2B['_'] = 27; sC2B['&'] = 28; sC2B['['] = 29; sC2B[']'] = 30; sC2B['.'] = 31;
    // 1 mark '#' mapped to cBitsSeparator
    sC2B['#'] = cBitsSeparator;
    // Digits [0-3] mapped to their binary value (in 2 bits)
    sC2B['0'] = 0; sC2B['1'] = 1; sC2B['2'] = 2; sC2B['3'] = 3;

    // Inverse map, set to 'F' by default
    for( int i=0; i<32; i++ ) sB2C[i] = 'F';
    sB2C[0] = 'F';
    
    sB2C[1] = 'a'; sB2C[2] = 'b'; sB2C[3] = 'c'; sB2C[4] = 'd'; sB2C[5] = 'e';
    sB2C[6] = 'f'; sB2C[7] = 'g'; sB2C[8] = 'h'; sB2C[9] = 'i'; sB2C[10] ='j';
    sB2C[11] ='k'; sB2C[12] ='l'; sB2C[13] ='m'; sB2C[14] ='n'; sB2C[15] ='o';
    sB2C[16] ='p'; sB2C[17] ='q'; sB2C[18] ='r'; sB2C[19] ='s'; sB2C[20] ='t';
    sB2C[21] ='u'; sB2C[22] ='v'; sB2C[23] ='w'; sB2C[24] ='x'; sB2C[25] ='y'; sB2C[26] ='z';
    
    sB2C[27] = '_'; sB2C[28] = '&'; sB2C[29] = '['; sB2C[30] = ']'; sB2C[31] = '.';

    // Separator and digits are NEVER inverted, so no need to store them

    sIsInitialized = true;
}

} //namespace util
