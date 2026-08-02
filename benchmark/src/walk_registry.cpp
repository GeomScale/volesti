#include "walk_registry.hpp"
#include "walk_adapters.hpp" 

using namespace std;

// Existing registry functions
WalkRegistry& get_walk_registry() {
    static WalkRegistry registry;
    return registry;
}

void register_walk(const string& name, RunFunction fn) {
    get_walk_registry()[name] = move(fn);
}

void initialize_all_walks() {
    // Uniform Walks
    register_walk("BallWalk", execute_walk<BallWalkType>);
    register_walk("BilliardWalk", execute_walk<BilliardWalkType>);
    register_walk("AcceleratedBilliardWalk", execute_walk<AcceleratedBilliardWalkType>);
    register_walk("SparseBilliardWalk", execute_walk<SparseBilliardWalkType>);
    register_walk("CDHRWalk", execute_walk<CDHRWalkType>);
    register_walk("RDHRWalk", execute_walk<RDHRWalkType>);
    register_walk("DikinWalk", execute_walk<DikinWalkType>);
    register_walk("JohnWalk", execute_walk<JohnWalkType>);
    register_walk("VaidyaWalk", execute_walk<VaidyaWalkType>);

    // Gaussian Walks
    register_walk("GaussianBallWalk", execute_walk<GaussianBallWalkType>);
    register_walk("GaussianCDHRWalk", execute_walk<GaussianCDHRWalkType>);

    // Shake-and-Bake Walks
    register_walk("ShakeAndBakeWalk", execute_walk<ShakeAndBakeWalkType>);
    register_walk("BilliardShakeAndBakeWalk", execute_walk<BilliardSBWalkType>);

    // Boundary Walks
    register_walk("BCDHRWalk", execute_walk<BCDHRWalkType>);
    register_walk("BRDHRWalk", execute_walk<BRDHRWalkType>);

    //Riemannian Walk
    register_walk("CRHMCWalk", execute_walk<CRHMCWalk>);
}