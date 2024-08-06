#include <gtest/gtest.h>

#include "config_sbu.h"
#include "obdetailstest.cpp"
#include "invectortest.cpp"

int main(int argc, char** argv) {
#ifdef _WIN32
	_putenv_s("BABEL_DATADIR", LOCAL_OB_DATADIR);
	_putenv_s("BABEL_LIBDIR", LOCAL_OB_LIBDIR);
#else
	setenv("BABEL_DATADIR", LOCAL_OB_DATADIR, 1);
	setenv("BABEL_LIBDIR", LOCAL_OB_LIBDIR, 1);
#endif
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
