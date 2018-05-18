







template <typename K>
class stdHPolytope{
private:
    typedef std::vector<K>        stdCoeffs;
    typedef std::vector<stdCoeffs>  stdMatrix;
    int d;
    int m;
    stdMatrix A;
    stdCoeffs b;
    stdCoeffs HypSum;
    //EXPERIMENTAL
    //typedef std::vector<flann::Index<flann::L2<double> > >  Flann_trees;

public:
    stdHPolytope() {}

    // constructor: cube(d)
    stdHPolytope(int dim){
        d=dim;
        for(int i=0; i<d; ++i){
            stdCoeffs coeffs;
            coeffs.push_back(K(1));
            for(int j=0; j<d; ++j){
                if(i==j)
                    coeffs.push_back(K(1));
                else coeffs.push_back(K(0));
            }
            A.push_back(coeffs);
        }
        for(int i=0; i<d; ++i){
            stdCoeffs coeffs;
            coeffs.push_back(K(1));
            for(int j=0; j<d; ++j){
                if(i==j)
                    coeffs.push_back(K(-1));
                else coeffs.push_back(K(0));
            }
            A.push_back(coeffs);
        }
    }

    int dimension(){
        return d;
    }

    int num_of_hyperplanes(){
        return m;
    }

    K get_coeff(int i, int j){
        return A[i][j];
    }

    void put_coeff(int i, int j, K value){
        A[i][j] = value;
    }

    // default initialize: cube(d)

    int init(int dim, stdMatrix Pin, stdCoeffs bin){
        d = dim;
        m=bin.size();
        HypSum.resize(m);
        typename stdMatrix::iterator pit=Pin.begin();
        typename stdCoeffs::iterator rit;
        for( ; pit<Pin.end(); ++pit){
            rit=(*pit).begin();
            for ( ; rit<(*pit).begin(); ++rit){
                HypSum[]+=//develop
            }
            A.push_back(*pit);
        }
        
        
        int i;
        for (i=0; i<m; i++){
            b[i]=bin[i];
        }
        //double tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        //std::random_shuffle (_A.begin(), _A.end());
        //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //std::shuffle (_A.begin(), _A.end(), std::default_random_engine(seed));
        //std::shuffle (_A.begin(), _A.end(), std::default_random_engine(seed));
        //boost::random::shuffle_order_engine<stdMatrix>(_A.begin(), _A.end());
        //double tstop = (double)clock()/(double)CLOCKS_PER_SEC;
        //std::cout << "Shuffle time = " << tstop - tstart << std::endl;
        return 0;
    }

    // print polytope in input format
    int print() {
        std::cout<<" "<<A.size()<<" "<<d+1<<" float"<<std::endl;
        for(typename stdMatrix::iterator mit=A.begin(); mit<A.end(); ++mit){
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit)
                std::cout<<*lit<<" ";
            std::cout<<std::endl;
        }
        return 0;
    }
    
    

    int is_in(Point p) {
        //std::cout << "Running is in" << std::endl;
        //exit(1);
        K sum;
        for (int i=0; i<m i++){
            sum=0.0;
            for (int j=0; j<d; j++){
                sum+=A[i][j]*p[j];
            }
            if (sum>b[i]){
                return 0;
            }
        }
        
        return -1;
    }
    
    /*
    int chebyshev_center(Point& center, double& radius){
        typedef CGAL::Linear_program_from_iterators
                <K**,                             // for A
                K*,                              // for b
                CGAL::Const_oneset_iterator<CGAL::Comparison_result>,  // for r
                bool*,                           // for fl
                K*,                              // for l
                bool*,                           // for fu
                K*,                              // for u
                K*>                              // for c
                Program;
        typedef CGAL::Quadratic_program_solution<ET> Solution;

        //std::cout<<"Cheb"<<std::endl;
        K* b = new K[_A.size()];
        K** A_col = new K*[_d+1];
        for(size_t i = 0; i < _d+1; ++i)
            A_col[i] = new K[_A.size()];

        stdMatrix B(_A);
        std::random_shuffle (B.begin(), B.end());
        for(size_t i=0; i<B.size(); ++i){
            K sum_a2 = 0;
            b[i] = B[i][0];
            for(size_t j=0; j<_d; ++j){
                A_col[j][i] = B[i][j+1];
                sum_a2 += std::pow(B[i][j+1],2);
            }
            A_col[_d][i] = std::sqrt(sum_a2);
        }

        CGAL::Const_oneset_iterator<CGAL::Comparison_result>
                r(CGAL::SMALLER);

        bool* fl = new bool[_d+1]();
        bool* fu = new bool[_d+1]();

        for(size_t i=0; i<_d+1; ++i)
            fl[i]=false;
        for(size_t i=0; i<_d+1; ++i)
            fu[i]=false;
        fl[_d]=true;

        K* l = new K[_d+1]();
        K* u = new K[_d+1]();
        K* c = new K[_d+1]();

        for(size_t i=0; i<_d+1; ++i)
            l[i]=K(0);
        for(size_t i=0; i<_d+1; ++i)
            u[i]=K(0);
        for(size_t i=0; i<_d+1; ++i)
            c[i]=K(0);
        c[_d]=K(-1);

        Program lp (_d+1, int(_A.size()), A_col, b, r,
                    fl, l, fu, u, c, 0);

        //CGAL::Quadratic_program_options options;
        //options.set_verbosity(1);                         // verbose mode
        //options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule
        //options.set_auto_validation(true);                // automatic self-check
        //Solution s = CGAL::solve_linear_program(lp, ET(), options);

        Solution s = CGAL::solve_linear_program(lp, ET());
        // output solution
        //std::cout << s;
        if (s.is_infeasible()){
            std::cout << "The polytope P is unbounded and Vol(P)=0\n";
            exit(-1);
        }
        else {
            assert (s.is_optimal());
            Solution::Variable_value_iterator it = s.variable_values_begin();
            std::vector<double> vecp;
            for(; it!=s.variable_values_end()-1; ++it){
                //std::cout<<CGAL::to_double(*it)<<" ";
                vecp.push_back(CGAL::to_double(*it));
            }
            center = Point(_d,vecp.begin(),vecp.end());
            //std::cout << center;
            radius = CGAL::to_double(*it);
            //std::cout << radius << std::endl;
        }
        // deallocate memory
        delete [] l;
        delete [] u;
        delete [] c;
        delete [] fl;
        delete [] fu;
        for(size_t i = 0; i < _d+1; ++i)
            delete [] A_col[i];
        delete [] b;
        return 0;
    } */

    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by _A
    std::pair<Point,Point> line_intersect(Point p,
                                          Point v){
        //std::cout<<"line-polytope-intersection"<<std::endl;
        K lambda_min=-1.0 * numeric_limits<double>::max();
        K lambda_max= numeric_limits<double>::max();
        
        for (int i=0; i<m; i++){
            for (int j=0; j<d; j++){
                
            }
        }
        
        return std::pair<Point,Point> (r+(lambda_min*v),r+(lambda_max*v));
    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          int rand_coord){
        //std::cout<<"line-polytope-intersection"<<std::endl;
        K lamda=0;
        //std::vector<NT> new_lamdas(_A.size());
        //std::vector<NT> new_lamdas;
        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        //std::vector<NT>::iterator lamdait = lamdas.begin();
        for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
            typename stdCoeffs::iterator cit;
            Point::Cartesian_const_iterator rit;
            rit=r.cartesian_begin();
            //Point::Cartesian_const_iterator vit;
            //vit=v.cartesian_begin();
            cit=ait->begin();
            K sum_nom=(*cit);
            ++cit;
            //here we just need the "rand_coord" coordinate of c
            //std::cout<<*(cit+rand_coord)<<"c= "<<std::endl;
            //for(typename stdCoeffs::iterator cit2=ait->begin() ; cit2 < ait->end() ; ++cit2){
            //  std::cout<<*cit2<<" ";
            //}
            //std::cout<<std::endl;
            K sum_denom= *(cit+rand_coord);
            //std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
            //         v.cartesian_begin()-v.cartesian_end()<<std::endl;
            for( ; cit < ait->end() ; ++cit, ++rit){
                //std::cout << sum_nom << " " << sum_denom <<std::endl;
                //std::cout << int(rit-r.cartesian_begin()) << " " << int(vit-v.cartesian_begin()) <<std::endl;
                sum_nom -= *cit * (*rit);
                //sum_denom += *cit * (*vit);
            }
            //std::cout << sum_nom << " / "<< sum_denom<<std::endl;
            if(sum_denom==K(0)){
                //std::cout<<"div0"<<std::endl;
                ;
            }
            else{
                lamda = sum_nom*(1/sum_denom);
                //lamdas[ait-_A.begin()] = lamda;

                if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
                if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
                if(lamda<min_plus && lamda>0) min_plus=lamda;
                if(lamda>max_minus && lamda<0) max_minus=lamda;
            }
            //std::cout<<r+(lamda*v)<<"\n"<<lamda<<std::endl;
        }
        return std::pair<NT,NT> (min_plus,max_minus);
    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<NT> &lamdas,
                                          bool init){
        //std::cout<<"line-polytope-intersection"<<std::endl;
        K lamda=0;
        std::vector<NT>::iterator lamdait = lamdas.begin();

        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        int mini, maxi;

        if(init){ //first time compute the innerprod cit*rit
            for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
                typename stdCoeffs::iterator cit;
                Point::Cartesian_const_iterator rit;
                rit=r.cartesian_begin();
                cit=ait->begin();
                K sum_nom=(*cit);
                ++cit;
                K sum_denom= *(cit+rand_coord);
                //std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
                //         std::endl;
                for( ; cit < ait->end() ; ++cit, ++rit){
                    sum_nom -= *cit * (*rit);
                }
                lamdas[ait-_A.begin()] = sum_nom;
                if(sum_denom==K(0)){
                    //std::cout<<"div0"<<std::endl;
                    ;
                }
                else{
                    lamda = sum_nom*(1/sum_denom);

                    if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
                    if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
                    if(lamda<min_plus && lamda>0) min_plus=lamda;
                    if(lamda>max_minus && lamda<0) max_minus=lamda;
                    //TEST
                    //if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;mini=ait-_A.begin();}
                    //if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;
                    //	 maxi=ait-_A.begin();}
                    //if(lamda<min_plus && lamda>0) {min_plus=lamda;mini=ait-_A.begin();}
                    //if(lamda>max_minus && lamda<0) {max_minus=lamda;maxi=ait-_A.begin();}
                }
            }
        } else {//only a few opers no innerprod
            for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
                typename stdCoeffs::iterator cit;
                cit=ait->begin();
                ++cit;

                NT c_rand_coord = *(cit+rand_coord);
                NT c_rand_coord_prev = *(cit+rand_coord_prev);

                *lamdait = *lamdait
                        + c_rand_coord_prev * (r_prev[rand_coord_prev] - r[rand_coord_prev]);

                if(c_rand_coord==K(0)){
                    //std::cout<<"div0"<<std::endl;
                    ;
                } else {
                    lamda = (*lamdait) / c_rand_coord;

                    if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
                    if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
                    if(lamda<min_plus && lamda>0) min_plus=lamda;
                    if(lamda>max_minus && lamda<0) max_minus=lamda;
                    //TEST
                    //if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;mini=ait-_A.begin();}
                    //if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;
                    //	 maxi=ait-_A.begin();}
                    //if(lamda<min_plus && lamda>0) {min_plus=lamda;mini=ait-_A.begin();}
                    //if(lamda>max_minus && lamda<0) {max_minus=lamda;maxi=ait-_A.begin();}
                }
                ++lamdait;
            }
        }
        //std::cout<<"Oresult: "<<mini<<" "<<maxi<<std::endl;
        return std::pair<NT,NT> (min_plus,max_minus);
    }

    //void rotate(){
    //std::cout<<_A<<std::endl;
    //  exit(1);
    //}

private:
    int            _d; //dimension
    stdMatrix      _A; //inequalities
    //EXPERIMENTAL
    //Flann_trees    flann_trees; //the (functional) duals of A lifted to answer NN queries
    //defined for every d coordinate
};
