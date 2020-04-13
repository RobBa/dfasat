// tails.cpp

// main() provided by Catch in file tests.cpp file.

#include "source/utility/catch.hpp"
#include "source/inputdata.h"

TEST_CASE( "2: Empty tail does not have a future", "[multi-file:2]" ) {

	tail *t = NULL;
    	t = new tail(0, 0, NULL);
	REQUIRE( t->future() == NULL);

}


