
Rcpp::NumericMatrix sample_R(Rcpp::NumericMatrix A, int walk_len, int rnum,  Rcpp::NumericVector internal_point, bool gaussian,
                             double variance, bool ball_walk, double delta, bool Vpoly, bool coord, bool verbose){

    bool rand_only=false, NN=false, birk=false;

    // define polytopes
    HPolytope<NT> P;
    VPolytope<NT> VP;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    int n=A.cols()-1, m = A.rows()-1, i, j=0;
    if (verbose) {
        if(!Vpoly) {
            std::cout<<"dimension = "<<n<<"\nnumber of facets = "<<m<<"\nnumber of points to sample = "<<rnum<<std::endl;
        } else {
            std::cout<<"dimension = "<<n<<"\nnumber of vertices = "<<m<<"\nnumber of points to sample = "<<rnum<<std::endl;
        }
    }
    std::vector<std::vector<NT> > Pin(m+1, std::vector<NT>(n+1));
    NT a = 1.0 / (2.0 * variance);
    for (i=0; i<m+1; i++){
        for(j=0; j<n+1; j++){
            Pin[i][j]=A(i,j);
        }
    }
    // construct polytope
    if (!Vpoly) {
        P.init(Pin);
    } else {
        VP.init(Pin);
    }

    // compute chebychev ball for H-polytope or an internal ball for V-polytope
    // we are going to need radius even an internal point is given as input
    std::pair<Point, NT> CheBall;
    if (!Vpoly) {
        CheBall = P.chebyshev_center();
    } else {
        CheBall = VP.chebyshev_center();
    }
    // print chebychev ball in verbose mode
    if (verbose) {
        if (!Vpoly) {
            std::cout << "\nChebychev center = " << std::endl;
            for (i = 0; i < n; i++) {
                std::cout << CheBall.first[i] << " ";
            }
            std::cout << "\nradius of chebychev ball = " << CheBall.second << std::endl;
        } else {
            std::cout << "\ncenter of internal ball  = " << std::endl;
            for (i = 0; i < n; i++) {
                std::cout << CheBall.first[i] << " ";
            }
            std::cout << "\nradius of internal ball = " << CheBall.second << std::endl;
        }
    }
    if (internal_point.size() == n){// if it is given as an input
        std::vector<NT> temp_p;
        std::cout<<"internal point given as an input"<<std::endl;
        for (int j=0; j<n; j++){
            temp_p.push_back(internal_point[j]);
        }
        CheBall.first = Point( n , temp_p.begin() , temp_p.end() );
    }
    Point p = CheBall.first;

    std::list<Point> randPoints;
    vars var(rnum,n,walk_len,1,0.0,0.0,0,0.0,0,CheBall.second,rng,urdist,urdist1,
             delta,verbose,rand_only,false,NN,birk,ball_walk,coord);
    if (ball_walk) {
        if (delta<0.0){
            if(!gaussian) {
                var.delta = 4.0 * CheBall.second / std::sqrt(NT(n));
            } else {
                var.delta = 4.0 * CheBall.second / std::sqrt(std::max(NT(1.0), a) * NT(n));
            }
        }
    }
    if(verbose) std::cout<<"\ncomputing first random point..."<<std::endl;
    Point q = get_point_on_Dsphere(n, CheBall.second);
    p=p+q;
    if (!Vpoly) {
        rand_point_generator(P, p, 1, 50 * n, randPoints, var);
    } else {
        rand_point_generator(VP, p, 1, 50 * n, randPoints, var);
    }
    if(verbose) std::cout<<"\nfirst random point computed!"<<std::endl;
    if(verbose) std::cout<<"p = ";
    if(verbose) p.print();

    randPoints.clear();
    if(verbose) std::cout<<"\nsampling points..."<<std::endl;
    if (!gaussian){
        if (!Vpoly) {
            rand_point_generator(P, p, rnum, walk_len, randPoints, var);
        } else {
            rand_point_generator(VP, p, rnum, walk_len, randPoints, var);
        }
    } else {
        vars_g var2(n, walk_len, 0, 0, 1, 0, CheBall.second, rng, 0, 0, 0, delta, false, verbose,
                    rand_only, false, NN, birk, ball_walk, coord);
        if (!Vpoly) {
            rand_gaussian_point_generator(P, p, rnum, walk_len, randPoints, a, var2);
        } else {
            rand_gaussian_point_generator(VP, p, rnum, walk_len, randPoints, a, var2);
        }
    }
    if(verbose) std::cout<<"\nsampling completed!\n"<<std::endl;

    Rcpp::NumericMatrix PointSet(n,rnum);

    typename std::list<Point>::iterator rpit=randPoints.begin();
    typename std::vector<NT>::iterator qit;
    j = 0;
    for ( ; rpit!=randPoints.end(); rpit++, j++) {
        qit = (*rpit).iter_begin(); i=0;
        for ( ; qit!=(*rpit).iter_end(); qit++, i++){
            PointSet(i,j)=*qit;
        }
    }

    return PointSet;

}
