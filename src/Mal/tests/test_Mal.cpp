//Test whole Mal

#include "tests.h"
#include <iostream>

int main()
{
    int failed = 0;
    failed += mal::TestAllHeaders() ? 0:1;
    failed += mal::TestTransform()  ? 0:1;
    failed += mal::TestTemplateInstantiation() ? 0:1;
    failed += mal::TestQuatPrecInFxp() ? 0:1;
    failed += mal::TestInterval() ? 0:1;
    if( 0 == failed )
        std::cout << "testMal OK" << std::endl;
    else
        std::cout << "testMal ERRORS = " << failed << std::endl;
    return (failed == 0);
}
