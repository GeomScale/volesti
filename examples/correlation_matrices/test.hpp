template <typename NT, typename WalkType, typename RNGType>
void naive_uniform_test(int n, int num_points, int walkL, int nreflex){
    typedef Cartesian<NT>                                               Kernel;
    typedef typename Kernel::Point                                      Point;
    typedef std::vector<Point>                                          PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>             MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 
    typedef typename WalkType::template Walk<Spectrahedron<Point>, RNGType>   RandomWalk;
    typedef RandomPointGenerator<RandomWalk>                            Generator;

    std::chrono::steady_clock::time_point start, end;
    double time;

    PointList randPoints;
    std::cout << "Direct implementation : " << n << std::endl;
    
    start = std::chrono::steady_clock::now();

    direct_uniform_sampling<NT, WalkType, RNGType, Point>(n, num_points, walkL, randPoints, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << "(ms)" << std::endl;

    write_to_file("sampling.txt", randPoints);
    check_output<VT,MT,PointList>(randPoints, num_points, n);
}

template <typename NT, typename WalkType, typename RNGType>
void new_test(unsigned int n, unsigned int const num_points, unsigned int walkL, unsigned int nreflex){

    std::cout << "Improved implementation : " << std::endl;

    typedef Cartesian<NT>           Kernel;
    typedef typename Kernel::Point  Point;
    typedef std::vector<Point> PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT; 
    
    typedef CorreSpectra<Point>     CorreSpectraType;
    typedef RandomPointGenerator<WalkType> Generator;

    PointList randPoints;
    CorreSpectraType P(n);

    const unsigned int d = P.dimension();
    Point startingPoint(d);
    RNGType rng(d);
    
    auto start = std::chrono::steady_clock::now();

    // uniform_correlation_sampling<WalkType>(P, randPoints, rng, walkL, num_points, startingPoint, 0);
    VT p(3);
    p[0] = 1./2;
    p[1] = 1./3;
    p[2] = 1./8;

    auto end = std::chrono::steady_clock::now();

    double time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Improved time : " << time << std::endl;

    // write_to_file("sampling_new.txt", randPoints);
    // check_output<VT,MT,PointList>(randPoints, num_points, n);
    
    // MT samples(d,num_points);
    // int j = 0;
    // for (typename std::vector<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); ++rpit, ++j)
    //     samples.col(j) = (*rpit).getCoefficients();

    // VT score = univariate_psrf<NT, VT, MT>(samples);
    // std::cout << "psrf = " << score.maxCoeff() << std::endl;

    // CHECK(score.maxCoeff() < 2.2);
    // return samples;
}