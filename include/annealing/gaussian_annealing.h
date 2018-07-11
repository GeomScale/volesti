// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

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


#ifndef GAUSSIAN_ANNEALING_H
#define GAUSSIAN_ANNEALING_H

#include <complex>

//Functions below from http://www.cplusplus.com/forum/beginner/130733/
//Function for average
NT mean ( std::vector<NT> &v )
{
    NT return_value = 0.0;
    int n = v.size();

    for ( int i=0; i < n; i++)
    {
        return_value += v[i];
    }

    return ( return_value / (NT(n)));
}
//****************End of average funtion****************


//Function for variance
NT variance ( std::vector<NT> &v , NT m )
{
    NT sum = 0.0;
    NT temp;
    int n=v.size();

    for ( int j =0; j < n; j++)
    {
        //temp = pow(v[j] - m,2.0);
        temp=(v[j]-m)*(v[j]-m);
        sum += temp;
    }

    return sum/(NT(n -1));
}
//****************End of variance funtion****************



template <class T1>
int get_first_gaussian(T1 &K, NT radius, NT &error, std::vector<NT> &a_vals, NT frac, vars var){

    int m=K.num_of_hyperplanes(), dim=var.n, its=0;
    NT sum, lower=0.0, upper=1.0, sigma_sqd, t, mid;
    std::vector<NT> dists(m,0);
    for(int i=0; i<m; i++){
        sum=0.0;
        for(int j=1; j<dim+1; j++){
            sum+=K.get_coeff(i,j)*K.get_coeff(i,j);
        }
        dists[i]=K.get_coeff(i,0)/std::sqrt(sum);
    }

    NT d=std::pow(10.0,10.0);

    while(its<10000){
        its+=1;
        sum=0.0;
        for(int i=0; i<dists.size(); i++){
            sum+=std::exp(-upper*std::pow(dists[i],2.0))/(2.0*dists[i]*std::sqrt(M_PI*upper));
        }

        sigma_sqd = 1/(2.0*upper);

        //t = (d-sigma_sqd*(NT(dim)))/(sigma_sqd*std::sqrt(((double)dim)));
        //sum+=std::exp(std::pow(-t,2.0)/8.0);
        if(sum>frac*error){//} || t<=1){
            upper=upper*10;
        }else{
            break;
        }
    }

    if(its==10000){
        std::cout<<"Cannot obtain sharp enough starting Gaussian"<<std::endl;
        exit(-1);
    }

    //get a_0 with binary search
    while(upper-lower>0.0000001){
        mid = (upper+lower)/2.0;
        sum=0.0;
        for(int i=0; i<dists.size(); i++){
            sum += std::exp(-mid*std::pow(dists[i],2)) / (2*dists[i]*std::sqrt(M_PI*mid));
        }

        sigma_sqd = 1.0/(2.0*mid);

        //t = (d-sigma_sqd)/(sigma_sqd*std::sqrt(((double)dim)));
        //sum += std::exp(std::pow(-t,2.0)/8.0);
        if(sum<frac*error){//} && t>1){
            upper=mid;
        }else{
            lower=mid;
        }
    }

    a_vals.push_back((upper+lower)/2.0);
    error = (1.0-frac)*error;

    return 1;
}


template <class T1>
int get_next_gaussian(T1 K,std::vector<NT> &a_vals, NT a, int N, NT ratio, NT C, Point &p, vars var){
    //its++;
    NT last_a = a, last_ratio=0.1, avg,average;
    bool done=false, print=var.verbose;
    int k=1, i;
    std::vector<NT> fn(N,0);

    //sample N points using hit and run
    std::list<Point> randPoints;
    vars var2=var;
    var2.coordinate=false;
    rand_gaussian_point_generator(K, p, N, 1, randPoints, last_a, var);
    if(print) std::cout<<"size of randPoints = "<<randPoints.size()<<" N = "<<N<<std::endl;

    while(!done){
        a = last_a*std::pow(ratio,(NT(k)));
        if(print) std::cout<<"a = "<<a<<std::endl;

        i=0;
        for(std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit){
            fn[i] = eval_exp(*pit,a)/eval_exp(*pit, last_a);
            //if(print) std::cout<<"fn[i] = "<<fn[i]<<std::endl;
            i++;
        }
        avg=mean(fn);
        //average = std::accumulate( fn.begin(), fn.end(), 0.0)/fn.size();
        //if(print) std::cout<<"avg = "<<avg<<std::endl;
        //if(print) std::cout<<"average = "<<average<<std::endl;
        //if(print) std::cout<<"var = "<<variance(fn,avg)<<std::endl;
        //if(print) std::cout<<"var/avg^2 = "<<variance(fn,avg)/std::pow(avg,2.0)<<std::endl;
        //if(print) std::cout<<"avg/last_ratio = "<<avg/last_ratio<<std::endl;
        //if(print) std::cout<<"C = "<<C<<std::endl;
        //if(print) std::cout<<"k = "<<k<<"\n"<<std::endl;

        if(variance(fn,avg)/std::pow(avg,2.0)>=C || avg/last_ratio<1.00001){
            if(k!=1){
                k=k/2;
            }
            done=true;
        }else{
            k=2*k;
        }
        last_ratio = avg;
    }
    if(print) std::cout<<"k = "<<k<<std::endl;
    a_vals.push_back(last_a*std::pow(ratio, (NT(k)) ) );
    if(print) std::cout<<"Next a found = "<<last_a*std::pow(ratio, (NT(k)) ) <<std::endl;

    return 1;
}


template <class T1>
int get_annealing_schedule(T1 K, std::vector<NT> &a_vals, NT &error, NT radius, NT ratio, NT C, vars var){
    bool print=var.verbose;
    get_first_gaussian(K, radius, error, a_vals, 0.1, var);
    if(print) std::cout<<"first gaussian computed"<<std::endl;
    if(print) std::cout<<"first gaussian a_0 = "<<a_vals[0]<<std::endl;
    if(print) std::cout<<"error is = "<<error<<"\n"<<std::endl;
    NT a_stop=0.0, curr_fn=2.0, curr_its=1.0;
    int it=0, dim=K.dimension(), steps=((int)150/error)+1;
    std::list<Point> randPoints;

    if(a_vals[0]<a_stop){
        a_vals[0]=a_stop;
    }

    Point p(K.dimension());

    vars var2=var;
    var2.coordinate=false;
    if(print) std::cout<<"Computing the sequence of gaussians.."<<std::endl;
    if(print) std::cout<<"N = "<<500 * ((int) C) + ((int) (dim * dim / 2))<<"\n"<<std::endl;
    if(print) std::cout<<"ratio = "<<ratio<<std::endl;
    while(curr_fn/curr_its>1.001 && a_vals[it]>=a_stop) {
        get_next_gaussian(K, a_vals, a_vals[it], 500 * ((int) C) + ((int) (dim * dim / 2)), ratio, C, p, var);
        it++;
        if(print) std::cout<<"gaussian a_"<<it<<" = "<<a_vals[it]<<std::endl;

        curr_fn = 0;
        curr_its = 0;

        for (int j = 0; j < steps; j++) {
            //if(print) std::cout<<"first time\n";
            //p.print();
            rand_gaussian_point_generator(K, p, 1, 1, randPoints, a_vals[it-1], var);
            //p.print();
            //if(!K.is_in(p)){
                //std::cout<<"point not in K";
                //exit(-1);
            //}
            curr_its += 1.0;
            curr_fn += eval_exp(p, a_vals[it]) / eval_exp(p, a_vals[it - 1]);
        }
    }
    if (a_vals[it]>a_stop){
        a_vals.pop_back();
        a_vals[it-1]=a_stop;
    }else{
        a_vals[it] = a_stop;
    }

    return 1;
}


#endif
