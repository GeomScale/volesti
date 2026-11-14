#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "misc/matrix_market_reader.h"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef HPolytope<Point> Hpolytope;

TEST_CASE("matrix_market_reader_degen2") {
    std::cout << "\n--- Testing Matrix Market reader on degen2" << std::endl;

    Hpolytope P = matrix_market_to_hpolytope<Point>(
        "../test/netlib/degen2.mm",
        "../test/netlib/degen2_bounds.mm"
    );

    CHECK(P.dimension() == 758);
    CHECK(P.num_of_hyperplanes() > 444);

    std::cout << "Dimension: " << P.dimension() << std::endl;
    std::cout << "Constraints: " << P.num_of_hyperplanes() << std::endl;
    std::cout << "Matrix Market reader test PASSED" << std::endl;
}
