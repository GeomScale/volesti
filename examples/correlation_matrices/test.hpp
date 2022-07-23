template <typename NT, typename WalkType, typename RNGType>
void naive_uniform_test(int n, int num_points, int walkL){
    typedef Cartesian<NT>                                               Kernel;
    typedef typename Kernel::Point                                      Point;
    typedef std::vector<Point>                                          PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>             MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 
    // typedef typename WalkType::template Walk<Spectrahedron<Point>, RNGType>   RandomWalk;
    // typedef RandomPointGenerator<RandomWalk>                            Generator;

    std::chrono::steady_clock::time_point start, end;
    double time;

    PointList randPoints;
    std::cout << "Direct implementation : " << n << " - time : ";
    
    start = std::chrono::steady_clock::now();

    direct_uniform_sampling<NT, WalkType, RNGType, Point>(n, num_points, walkL, randPoints, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << " (ms)" << std::endl;

    // write_to_file("sampling.txt", randPoints);
    // check_output<VT,MT,PointList>(randPoints, num_points, n);
}

template <typename NT, typename WalkType, typename RNGType>
void new_test(unsigned int n, unsigned int const num_points, unsigned int walkL){

    std::cout << "Improved implementation : " << n << " - time : ";

    typedef Cartesian<NT>           Kernel;
    typedef typename Kernel::Point  Point;
    typedef std::vector<Point> PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT; 
    
    typedef CorreSpectra<Point>     CorreSpectraType;

    std::chrono::steady_clock::time_point start, end;
    double time;

    PointList randPoints;
    CorreSpectraType P(n);

    const unsigned int d = P.dimension();
    Point startingPoint(d);
    RNGType rng(d);
    
    start = std::chrono::steady_clock::now();
    
    uniform_correlation_sampling<WalkType>(P, randPoints, rng, walkL, num_points, startingPoint, 0);
    
    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << " (ms)" << std::endl;

    // write_to_file("sampling_new.txt", randPoints);
    // check_output<VT,MT,PointList>(randPoints, num_points, n);
}
