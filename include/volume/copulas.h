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

#ifndef COPULAS_H
#define COPULAS_H


template <typename Point, typename RNGType, typename NT>
std::vector<std::vector<NT> > twoParHypFam(const int dim,
                                           const int num,
                                           const int num_slices,
                                           const std::vector<NT> &pl1, const std::vector<NT> &pl2,
                                           double seed = std::numeric_limits<double>::signaling_NaN())
{

    int i,j,col,row;
    std::vector<NT> vec1,vec2,Zs1,Zs2;
    NT sum1,sum2;
    std::list<Point> points;
    typename std::list<Point>::iterator rpit;
    std::pair< std::vector<NT>,std::vector<NT> > result;
    Point p;

    Sam_Canon_Unit<NT, RNGType> (dim, num, points, seed);

    std::vector<std::vector<int> > Matrix(num_slices);
    std::vector<std::vector<NT> > pos_Matrix(num_slices);
    for (i=0; i<num_slices; i++){
        Matrix[i].resize(num_slices);
        pos_Matrix[i].resize(num_slices);
    }
    for (i=0; i<num_slices; i++){
        for (j=0; j<num_slices; j++){
            Matrix[i][j]=0;
        }
    }

    rpit = points.begin();
    //for (i=0; i<num; i++){
    for( ; rpit!=points.end(); ++rpit){
        //p=points[i];
        p = (*rpit);
        //std::cout<<p<<std::endl;
        sum1=0.0; sum2=0.0;
        for (j=0; j<dim; j++){
            sum1+=p[j]*pl1[j];
            sum2+=p[j]*pl2[j];
            //sum+=p[j];

        }
        //std::cout<<sum<<std::endl;
        vec1.push_back(sum1);
        vec2.push_back(sum2);
    }
    std::sort( vec1.begin(), vec1.end() );
    std::sort( vec2.begin(), vec2.end() );

    NT pos = 1.0 / NT(num_slices);

    for (i=1; i<num_slices; i++){
        //std::cout<<((int)std::floor(i*(0.01)*((double)num) ))<<std::endl;
        //std::cout<<((int)std::floor(i*(0.01)*((double)num) ))<<std::endl;
        Zs1.push_back(vec1[((int)std::floor(i*(pos)*((NT)num) ))]);
        Zs2.push_back(vec2[((int)std::floor(i*(pos)*((NT)num) ))]);
    }



    //std::cout<<"hello2"<<std::endl;
    rpit = points.begin();
    //for (i=0; i<num; i++){
    for( ; rpit!=points.end(); ++rpit){
        //p=points[i];
        p = (*rpit);
        //std::cout<<"dimension is: "<<p.dimension()<<std::endl;
        sum1=0.0; sum2=0.0;
        col=-1; row=-1;
        for (j=0; j<dim; j++){
            sum1+=p[j]*pl1[j];
            sum2+=p[j]*pl2[j];
        }
        //std::cout<<"hello3"<<std::endl;
        for (unsigned int j=0; j<Zs1.size(); j++){
            if (sum1<Zs1[j]){
                col=j;
                break;
            }
        }
        for (unsigned int j=0; j<Zs2.size(); j++){
            if (sum2<Zs2[j]){
                row=j;
                break;
            }
        }
        if (col==-1){
            col=num_slices - 1;
        }

        if (row==-1){
            row=num_slices - 1;
        }
        //std::cout<<"hello4"<<std::endl;
        Matrix[row][col]++;
    }

    for (i=0; i<num_slices; i++){
        for (j=0; j<num_slices; j++){
            pos_Matrix[i][j]=((NT)Matrix[i][j])/((NT)num);
        }
    }

    return pos_Matrix;

}


template <typename Point, typename RNGType, typename ellipsoid, typename NT>
std::vector<std::vector<NT> > hypfam_ellfam(int dim,
                                            int num,
                                            int num_slices,
                                            std::vector<NT>  pl,
                                            ellipsoid G,
                                            double seed = std::numeric_limits<double>::signaling_NaN())
{

    int i,j,col,row;
    std::vector<NT> vec1,vec2,Zs1,Cs;
    NT sum1,sum2;
    std::list<Point> points;
    typename std::list<Point>::iterator rpit;
    std::pair< std::vector<NT>,std::vector<NT> > result;
    Point p;

    Sam_Canon_Unit<NT, RNGType> (dim, num, points, seed);

    std::vector<std::vector<int> > Matrix(num_slices);
    std::vector<std::vector<NT> > pos_Matrix(num_slices);
    for (i=0; i<num_slices; i++){
        Matrix[i].resize(num_slices);
        pos_Matrix[i].resize(num_slices);
    }
    for (i=0; i<num_slices; i++){
        for (j=0; j<num_slices; j++){
            Matrix[i][j]=0;
        }
    }

    rpit = points.begin();
    //for (i=0; i<num; i++){
    for( ; rpit!=points.end(); ++rpit){
        //p=points[i];
        p = (*rpit);
        //std::cout<<p<<std::endl;
        sum1=0.0;
        sum2=G.mat_mult(p);
        for (j=0; j<dim; j++){
            sum1+=p[j]*pl[j];
            //sum2+=p[j]*pl2[j];
            //sum+=p[j];

        }
        //std::cout<<sum<<std::endl;
        vec1.push_back(sum1);
        vec2.push_back(sum2);
    }
    std::sort( vec1.begin(), vec1.end() );
    std::sort( vec2.begin(), vec2.end() );


    NT pos = 1.0 / NT(num_slices);
    for (i=1; i<num_slices; i++){
        //std::cout<<((int)std::floor(i*(0.01)*((double)num) ))<<std::endl;
        //std::cout<<((int)std::floor(i*(0.01)*((double)num) ))<<std::endl;
        Zs1.push_back(vec1[((int)std::floor(i*(pos)*((NT)num) ))]);
        Cs.push_back(vec2[((int)std::floor(i*(pos)*((NT)num) ))]);
    }



    rpit = points.begin();
    //for (i=0; i<num; i++){
    for( ; rpit!=points.end(); ++rpit){
        //p=points[i];
        p = (*rpit);
        //std::cout<<"dimension is: "<<p.dimension()<<std::endl;
        sum1=0.0; sum2=0.0;
        col=-1; row=-1;
        sum2=G.mat_mult(p);
        for (j=0; j<dim; j++){
            sum1+=p[j]*pl[j];
            //sum2+=p[j]*pl2[j];
        }
        //std::cout<<"hello3"<<std::endl;
        for (unsigned int j=0; j<Zs1.size(); j++){
            if (sum1<Zs1[j]){
                col=j;
                break;
            }
        }
        for (unsigned int j=0; j<Cs.size(); j++){
            if (sum2<Cs[j]){
                row=j;
                break;
            }
        }
        if (col==-1){
            col=num_slices - 1;
        }

        if (row==-1){
            row=num_slices - 1;
        }
        //std::cout<<"hello4"<<std::endl;
        Matrix[row][col]++;
    }

    for (i=0; i<num_slices; i++){
        for (j=0; j<num_slices; j++){
            pos_Matrix[i][j]=((NT)Matrix[i][j])/((NT)num);
        }
    }

    return pos_Matrix;

}

#endif
