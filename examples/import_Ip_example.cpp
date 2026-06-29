#include <iostream>
#include "volesti/io/lp_importer.h"
#include "volesti/hpolytope.h"

int main()
{
    auto data = volesti::io::import_lp("data/sample.lp");

    HPolytope<double> P(data.A, data.b);

    std::cout << "Polytope dimension: "
              << P.dimension() << std::endl;

    return 0;
}