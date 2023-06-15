// VolEsti (volume computation and sampling library)

// Copyright (c) 2018-2020 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

#ifndef SAMPLERS_SIMPLEX_HPP
#define SAMPLERS_SIMPLEX_HPP

template <typename NT, typename RNGType, typename Point>
void Sam_Unit(unsigned int dim,
              unsigned int num,
              std::list<Point> &points,
              double seed = std::numeric_limits<double>::signaling_NaN())
{

    unsigned int j,i,x_rand,M=2147483647,pr,divisors,pointer;  // M is the largest possible integer
    std::vector<unsigned int> x_vec;
    std::vector<NT> y;

    boost::random::uniform_int_distribution<> uidist(1,M);
    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!std::isnan(seed)) {
        unsigned rng_seed2 = seed;
        rng.seed(rng_seed2);
    }
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);

    if (dim<=60){

        x_vec.assign(dim+1,0);

        for (i=0; i<num; i++){

            // Set all the point's coordinates equal to zero and clear the integers' list
            y.assign(dim,0); pointer=0;

            // Generate d distinct integers
            while ( pointer<dim ){

                x_rand = uidist(rng);
                // Check if this integer is selected first time
                if ( std::find(x_vec.begin(), x_vec.begin()+pointer, x_rand)
                     == x_vec.begin()+pointer )
                {
                    pointer++;
                    x_vec[pointer]=x_rand;
                }

            }

            // Sort the integers' list
            std::sort( x_vec.begin(), x_vec.end() );

            // Construct the point's coordinates
            for(j=0; j<dim; j++){
                y[j]=((NT)(x_vec[j+1]-x_vec[j])) / ((NT)M);
            }

            // Define the new point and add it to the point-list
            points.push_back(Point(dim,y.begin(),y.end()));

        }
    }else if(dim<=80){

        x_vec.assign(dim+1,0);

        std::vector<bool> filter;
        bool t1=true,t2=true;
        pr=3*dim+1;
        if(pr%2==0) pr+=1;

        while(t1){
            t2=true;
            divisors=(int)floor(sqrt((NT)pr))+1;
            for (i=3; i<divisors+1; i+=2){
                if (pr%i==0){
                    t2=false;
                    break;
                }
            }
            if(t2) break;
            pr+=2;
        }

        // Generate the number of points requested
        for (i=0; i<num; i++){

            // Set all the point's coordinates equal to zero and the filter equal to true
            y.assign(dim,0); filter.assign(pr,true);
            pointer=0;

            // Generate d distinct integers
            while ( pointer<dim ){
                x_rand = uidist(rng);

                // Check if this integer is the first that is mapped to the specific filter's position
                if ( filter[x_rand%pr] ){
                    filter[x_rand%pr]=false; pointer++;
                    x_vec[pointer]=x_rand;
                }
            }

            // Sort the integers' list
            std::sort( x_vec.begin(), x_vec.end() );

            // Construct the point's coordinates
            for(j=0; j<dim; j++){
                y[j]=((NT)(x_vec[j+1]-x_vec[j])) / ((NT)M);
            }

            // Define the new point and add it to the point-vector
            points.push_back(Point(dim,y.begin(),y.end()));

        }
    }else{

        std::vector<NT> x_vec2;
        NT Ti,sum;

        x_vec2.assign(dim+1,0.0);

        // Generate the number of points requested
        for (i=0; i<num; i++){

            // Set all the point's coordinates equal to zero and clear the integers' list
            pointer=0; sum=0.0;

            while ( pointer<dim+1 ){

                x_rand = uidist(rng);
                Ti=-log(((NT)x_rand)/((NT)M));
                sum+=Ti;
                x_vec2[pointer]= Ti; pointer++;

            }

            for (j=0; j<dim; j++){
                x_vec2[j]/=sum;
            }

            points.push_back(Point(dim,x_vec2.begin(),x_vec2.end()-1));

        }

    }

    return;

}

template <typename NT, typename RNGType, typename Point>
void Sam_Canon_Unit(unsigned int dim,
                    unsigned int num,
                    std::list<Point> &points,
                    double seed = std::numeric_limits<double>::signaling_NaN())
{

    unsigned int j,i,x_rand,M=2147483647,pointer;  // M is the largest possible integer
    std::vector<NT> y;
    dim--;
    boost::random::uniform_int_distribution<> uidist(1,M);

    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!std::isnan(seed)) {
        unsigned rng_seed2 = seed;
        rng.seed(rng_seed2);
    }

    std::vector<NT> x_vec2;
    NT Ti,sum;

    x_vec2.assign(dim+1,0.0);

    // Generate the number of points requested
    for (i=0; i<num; i++){

        // Set all the point's coordinates equal to zero and clear the integers' list
        pointer=0; sum=0.0;

        while ( pointer<dim+1 ) {

            x_rand = uidist(rng);
            Ti = -log(((NT) x_rand) / ((NT) M));
            sum += Ti;
            x_vec2[pointer] = Ti;
            pointer++;

        }

        for (j=0; j<dim+1; j++) {
            x_vec2[j] /= sum;
        }

        points.push_back(Point(dim+1, x_vec2.begin(), x_vec2.end()));

    }

    return;

}


//Owen mapping for sample from an arbitrary simplex given in V-represantation
template <typename Vpolytope, typename PointList>
void Sam_arb_simplex(const Vpolytope &P, unsigned int num, PointList &points){

    typedef typename Vpolytope::MT MT;
    typedef typename Vpolytope::NT NT;
    typedef typename Vpolytope::rngtype RNGType;
    typedef typename Vpolytope::PolytopePoint Point;

    MT V = P.get_mat();
    std::vector<Point> vec_point;

    unsigned int n=V.rows(),j,i,k,x_rand,M=2147483647,pr,divisors,pointer;  // M is the largest possible integer
    unsigned int dim = V.cols();
    std::vector<NT> temp_p(dim, 0.0);
    std::vector<unsigned int> x_vec;
    std::vector<NT> y;

    for (int i = 0; i < V.rows(); ++i) {
        for (int j = 0; j < V.cols(); ++j) {
            temp_p[j] = V(i,j);
        }
        vec_point.push_back(Point(dim, temp_p.begin(), temp_p.end()));
    }
    typename std::vector<Point>::iterator it_beg = vec_point.begin();

    NT Xj;
    Point p0=*it_beg;

    boost::random::uniform_int_distribution<> uidist(1,M);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    if (dim<=60){

        x_vec.assign(dim+2,0);  x_vec[dim+1]=M;

        for (i=0; i<num; i++){

            // Set all the point's coordinates equal to zero and clear the integers' list
            y.assign(dim,0); pointer=0;

            // Generate d distinct integers
            while ( pointer<dim ){

                x_rand = uidist(rng);
                // Check if this integer is selected first time
                if ( std::find(x_vec.begin(), x_vec.begin()+pointer, x_rand)
                     == x_vec.begin()+pointer )
                {
                    pointer++;
                    x_vec[pointer]=x_rand;
                }

            }

            // Sort the integers' list
            std::sort( x_vec.begin(), x_vec.end() );

            // Construct the point's coordinates
            for(j=0; j<dim+1; j++){
                Point pk=*(it_beg+j);
                Xj=((NT)(x_vec[j+1]-x_vec[j])) / ((NT)M);
                for (k=0; k<dim; k++){
                    y[k]+=Xj*pk[k];
                }
            }

            // Define the new point and add it to the point-list
            points.push_back(Point(dim,y.begin(),y.end()));

        }
    }else if(dim<=80){

        std::vector<bool> filter;
        x_vec.assign(dim+2,0);  x_vec[dim+1]=M;

        bool t1=true,t2=true;
        pr=3*dim+1;
        if(pr%2==0) pr+=1;

        while(t1){
            t2=true;
            divisors=(int)floor(sqrt((NT)pr))+1;
            for (i=3; i<divisors+1; i+=2){
                if (pr%i==0){
                    t2=false;
                    break;
                }
            }
            if(t2) break;
            pr+=2;
        }

        // Generate the number of points requested
        for (i=0; i<num; i++){

            // Set all the point's coordinates equal to zero and the filter equal to true
            y.assign(dim,0); filter.assign(pr,true);
            pointer=0;

            // Generate d distinct integers
            while ( pointer<dim ){
                x_rand = uidist(rng);

                // Check if this integer is the first that is mapped to the specific filter's position
                if ( filter[x_rand%pr] ){
                    filter[x_rand%pr]=false; pointer++;
                    x_vec[pointer]=x_rand;
                }
            }

            // Sort the integers' list
            std::sort( x_vec.begin(), x_vec.end() );

            // Construct the point's coordinates
            for(j=0; j<dim+1; j++){
                Point pk=*(it_beg+j);
                Xj=((NT)(x_vec[j+1]-x_vec[j])) / ((NT)M);
                for (k=0; k<dim; k++){
                    y[k]+=Xj*pk[k];
                }
            }

            // Define the new point and add it to the point-list
            points.push_back(Point(dim,y.begin(),y.end()));

        }
    }else{
        std::vector<NT> x_vec2;
        NT Ti,sum;
        x_vec2.assign(dim+1,0.0);

        // Generate the number of points requested
        for (i=0; i<num; i++){

            // Set all the point's coordinates equal to zero and clear the integers' list
            pointer=0; sum=0.0; y.assign(dim,0);
            while ( pointer<dim+1 ){

                x_rand = uidist(rng);
                Ti=-log(((NT)x_rand)/((NT)M));
                sum+=Ti;
                x_vec2[pointer]= Ti; pointer++;

            }

            for (j=0; j<dim+1; j++){
                x_vec2[j]/=sum;
                Point pk=*(it_beg+j);
                Xj=((NT)(x_vec[j+1]-x_vec[j])) / ((NT)M);
                for (k=0; k<dim; k++){
                    y[k]+=x_vec2[j]*pk[k];
                }
            }

            points.push_back(Point(dim,y.begin(),y.end()));

        }

    }

    return;

}


#endif
