//Test whole Util

#include "tests.h"
#include <iostream>

int main()
{
    int failed = 0;
    failed += (util::TestItemStream())? 0:1;
    failed += (util::TestChronos())? 0:1;
    failed += (util::TestTracer())? 0:1;
    if( 0 == failed )
        std::cout << "testUtil OK" << std::endl;
    else
        std::cout << "testUtil ERRORS = " << failed << std::endl;
    return (failed == 0);
}
