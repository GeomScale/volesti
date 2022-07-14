
template<typename NT, typename ConvexBodyType, typename Point>
std::pair<NT, int> intersection(ConvexBodyType &P, const Point &x, const Point &v, const unsigned int k){
    double tau, tmp;
    int j = 0;
    if(v[0] > 0){
        tau = (1-x[0])/v[0];   
    }else{
        tau = -(1 + x[0])/v[0];
    }
    for(int i = 1; i < k; ++i){
        if(v[i] > 0){
            tmp = (1 - x[i])/v[i];
        }else{
            tmp = -(1 + x[i])/v[i];
        }
        if(tau > tmp){
            tau = tmp;
            j = i;
        }
    }
    tmp = P.positiveLinearIntersection(x.getCoefficients(), v.getCoefficients());
    std::cout << tau << " , " << tmp << std::endl;
    if(tau > tmp){
        tau = tmp;
        j = -1;
    }
    std::pair<double, int> res(tau,j);
    return res;
}

template<typename ConvexBodyType, typename Point>
void reflection(ConvexBodyType &P, Point &p, Point &v, const int flag){
    if(flag != -1){
        v.set_coord(flag, - v.getCoefficients()(flag));
        return;
    }
    P.compute_reflection(v, p, flag);
}

// This function is taken and simplified from uniform_billiard_walk.hpp

template<typename NT, typename ConvexBodyType, typename Point, typename RNGType>
Point BilliardWalkSpectra(ConvexBodyType &P, Point& q, unsigned int const& walk_length, RNGType &rng, NT const _Len){
    unsigned int k = P.dimension();
    int nreflex = 50*k;
    NT L, tau;
    Point p = q, v;
    std::pair<double, int> pbpair;
    for (unsigned int j=0; j<walk_length; ++j){
        L = rng.sample_urdist() * _Len;
        v = GetDirection<Point>::apply(k, rng);
        Point p0 = p;
        int it = 0;
        while (it < nreflex)
        {
            pbpair = intersection<NT>(P, p, v, k);
            tau = pbpair.first;
            if (L <= tau) {
                p += (L * v);
                break;
            }
            tau = 0.995 * tau; // 0.995: to approximate boundary points?
            p += tau * v; // A point (almost) on the boundary
            L -= tau;
            reflection(P, p, v, pbpair.second);
            it++;
        }
        if (it == nreflex){
            p = p0;
        }
    }
    return p;
}

template <typename NT, typename RNGType, typename Point, typename PointList>
void direct_sampling2(unsigned int n, unsigned int num_points, unsigned int walkL, PointList &randPoints, int nburns){
    
    Spectrahedron<Point> spectra = prepare_input<NT, Point>(n);
    int d = spectra.dimension();
    Point p(d);
    RNGType rng(d);
    
    std::pair<Point, NT> inner_ball = spectra.ComputeInnerBall();
    NT diameter = d * NT(6) * inner_ball.second;
    
    for (unsigned int i = 0; i < num_points; ++i){
        p = BilliardWalkSpectra(spectra, p, walkL, rng, diameter);
        randPoints.push_back(p);
    }
}