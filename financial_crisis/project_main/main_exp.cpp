//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_d.h>
#include <vector>
#include "boost/random.hpp"
#include <boost/math/constants/constants.hpp>
#include <../project_include/Eigen/Dense>
#include <math.h>
#include <chrono>
#include <ctime>

#include <algorithm>
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/Gmpq.h>
#include <CGAL/point_generators_d.h>

//--------------------//
#include <../project_include/approx_utils.h>
#include <../project_include/random_samplers_proj.h>

#include <../project_include/exact_utils.h>
#include <../project_include/exact_vol.h>
#include <../project_include/approximate_vol.h>





typedef boost::mt19937 RNGType;
typedef CGAL::Gmpq                  NT3;
//typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel2;
//typedef Kernel2::RT					RT;

//typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
//typedef CGAL::Gmpq 			EXACT_NT;
//typedef double NT;
typedef CGAL::Cartesian_d<double> 	      Kernel; 
typedef Kernel::Point_d Point_d;
typedef CGAL::Timer				  timer;
typedef Kernel::Hyperplane_d 		Plane_d;

using namespace Eigen;


int main(int argc, char* argv[]){
	
	
	//int choice=atoi(argv[1]);
	
	
    if(!strcmp(argv[1],"-h")||!strcmp(argv[1],"-help")){
        std::cerr<<
                    "Usage:\n"<<
                    "-r 1: Run rejection for arbitrary simplex and two arbitrary hyperplanes compared to rejection from unit simplex\n"<<
                    "-r 2: Run rejection for 2 hyperplanes intersecting the unit simplex\n"<<
                    "-r 3: Run rejection for two families of parallel hyperplanes intersecting the facet simplex\n"<<
                    "-r 4: Run rejection for one family of parallel hyperplanes and one family of cocentric, at the origin, ellipsoids intersecting the facet simplex\n"<<
                    "-ve 1 1: Run VolEsti for an ellipsoid intersecting a unit simplex using random inscribed ball\n"<<
                    "-ve 1 2: Run VolEsti for an ellipsoid intersecting a unit simplex using inscribed ball solving conical optimization problem\n"<<
                    "-ve 2: Run VolEsti for an ellipsoid and two parallel hyperplanes intersecting the unit simplex\n"<<
                    "-ve 3: Run VolEsti for two parallel hyperplanes and two cocentric ellipsoids intersecting the unit simplex [non convex body]\n"<< 
                    "-varsi: Run Varsi formula \n"<<
                    "-Lawn 1: Run multiprecision Lawrence formula for two parallel hyperplanes intersecting the unit simplex\n"<<
                    "-Lawn 2: Run double precision Lawrence formula for two parallel hyperplanes intersecting the unit simplex\n"<<
                    std::endl;
        exit(-1);
    }
      
	if ( !strcmp(argv[1],"-r") ){
		int choice=atoi(argv[2]);
		
		
		if (choice==1){
			
			timer tim;
			int dim=atoi(argv[3]), num=atoi(argv[4]),j,i;
			double time1,time2,size=100.0,vol_unit,vol_arb;
            std::vector<Point_d> v;
            std::vector<double> pl1,pl2;
            int M=2147483647,diair=21474836,x_rand;
            RNGType rng;
            boost::random::uniform_int_distribution<> uidist(-M,M);
            
			
			CGAL::Random_points_on_sphere_d<Point_d> gen (dim, size); v.clear();
			for (j=0; j<(dim+1); j++){
				v.push_back (*gen++);
				x_rand = uidist(rng);
				pl1.push_back( ((double)x_rand)/( (double)diair ) );
				x_rand = uidist(rng);
				//x_rand+=13*j;
				pl2.push_back( ((double)(x_rand))/( (double)diair ) );
			}
			Plane_d plane1(dim,pl1.begin(),pl1.end());
			Plane_d plane2(dim,pl2.begin(),pl2.end());
			
			tim.reset();
			tim.start();
			vol_arb=sampleToArb(v.begin(),v.end(),plane1,plane2,num);
			time1=tim.time();
			
			tim.reset();
			tim.start();
			vol_unit=sampleToUnit(v.begin(), v.end(), plane1, plane2, num);
			time2=tim.time();
			
			std::cout<<"\n";
			std::cout<<"Rejection in  arbitrary simplex, volume: "<<vol_arb<<std::endl;
			std::cout<<"Rejection in  arbitrary simplex, time: "<<time1<<std::endl;
			std::cout<<"\n";
			std::cout<<"Using rejection in  unit simplex, volume: "<<vol_unit<<std::endl;
			std::cout<<"Using rejection in  unit simplex, time: "<<time2<<std::endl;
            
          /*  std::ofstream outputFile3;
            outputFile3.open("hyperplane1.txt");
            std::ofstream outputFile4;
            outputFile4.open("hyperplane2.txt");
            for (i=0; i<dim+1; i++){
                if (i<dim){
                    outputFile3<<plane1[i]<<" ";
                    outputFile4<<plane2[i]<<" ";
                    
                }else{
                    outputFile3<<-plane1[i]<<"\n";
                    outputFile4<<-plane1[i]<<"\n";
                }
            }*/
            
			
			
		}else if (choice==2){
			
			int dim=atoi(argv[3]), num=atoi(argv[4]),i,j;
			std::vector<double> planes,pl1,pl2;
            double time,vol_unit; timer tim;
			
			std::ifstream inputFile;
			std::string inputFileName = argv[5]; //reading the input txt file
			inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
			std::istream_iterator< double >  input_begin( inputFile );
			std::istream_iterator< double >  input_end;
			for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
				planes.push_back(*it);
			}
			
			
            for (i=0; i<dim+1; i++){
                pl1.push_back(planes[i]);
                pl2.push_back(planes[i+dim+1]);
            }
            Plane_d plane1(dim,pl1.begin(),pl1.end());
			Plane_d plane2(dim,pl2.begin(),pl2.end());
            
            tim.reset();
			tim.start();
			vol_unit=sample_2hyp(dim, num, plane1, plane2);
			time=tim.time();
			
			std::cout<<"\n";
            std::cout<<"Rejection in  unit simplex, volume: "<<vol_unit<<std::endl;
			std::cout<<"Rejection in  unit simplex, time: "<<time<<std::endl;
			
		}else if(choice==3){
         
            int dim=atoi(argv[3]), num=atoi(argv[4]),i,j;
			std::vector<double> planes,pl1,pl2;
            double time, z1, z2;
            timer tim;
            std::vector<std::vector<double> > pos_Matrix;
            
            std::ifstream inputFile;
			std::string inputFileName = argv[5]; //reading the input txt file
			inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
			std::istream_iterator< double >  input_begin( inputFile );
			std::istream_iterator< double >  input_end;
			for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
				planes.push_back(*it);
			}
			
			for (i=0; i<dim; i++){
                pl1.push_back(planes[i]);
                pl2.push_back(planes[i+dim]);
            }
            tim.reset();
			tim.start();
            pos_Matrix=get_par_4hyp_vols(dim, num, pl1,pl2);
            time=tim.time();
            
            std::cout<<"\n";
            std::cout<<"Excecutional time: "<<time<<std::endl;
            std::cout<<"The results are in comp_planes_res.txt"<<std::endl;
            
            std::ofstream outputFile3;
            outputFile3.open("comp_planes_res.txt");
            for (i=0; i<100; i++){
                for (j=0; j<100; j++){
                    if (j<99){
                        outputFile3<<pos_Matrix[i][j]<<" ";
                    }else{
                        outputFile3<<pos_Matrix[i][j]<<"\n";
                    }
                }
            }
            
            
        }else if(choice==4){
            
            int dim=atoi(argv[3]), num=atoi(argv[4]),i,j;
            //std::cout<<"hello"<<std::endl;
			std::vector<double> planes,pl1,ellips;
            double time, z1, z2;
            timer tim;
            std::vector<std::vector<double> > pos_Matrix;
            
            std::ifstream inputFile;
			std::string inputFileName = argv[5]; //reading the input txt file
			inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
			std::istream_iterator< double >  input_begin( inputFile );
			std::istream_iterator< double >  input_end;
			for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
				planes.push_back(*it);
			}
			
			for (i=0; i<dim; i++){
                if (i<dim){
                    pl1.push_back(planes[i]);
                }
            }
            
            std::ifstream inputFile2;
			std::string inputFileName2 = argv[6]; //reading the input txt file
			inputFile2.open(inputFileName2.c_str()); //converting the input to a c-style string
			
			std::istream_iterator< double >  input_begin2( inputFile2 );
			std::istream_iterator< double >  input_end2;
			for(std::istream_iterator< double > it2 = input_begin2; it2 != input_end2; ++it2){
				ellips.push_back(*it2);
			}
            //std::cout<<"hello"<<std::endl;
            MatrixXd C(dim,dim);
            VectorXd center(dim);
            
            for (i=0; i<dim; i++){
                center(i)=0;
                for (j=0; j<dim; j++){
                    C(i,j)=ellips[j+i*dim];
                }
            }
            
            ellipsoids G(C,center,0,dim);          
            tim.reset();
			tim.start();
            pos_Matrix=get_par_hypellips_vols(dim, num,pl1,G);
            time=tim.time();
            
            std::cout<<"\n";
            std::cout<<"Excecutional time: "<<time<<std::endl;
            std::cout<<"The results are in comp_plane_ellips_res.txt"<<std::endl;
            
            std::ofstream outputFile3;
            outputFile3.open("comp_plane_ellips_res.txt");
            for (i=0; i<100; i++){
                for (j=0; j<100; j++){
                    if (j<99){
                        outputFile3<<pos_Matrix[i][j]<<" ";
                    }else{
                        outputFile3<<pos_Matrix[i][j]<<"\n";
                    }
                }
            }
            
        }
		
	}else if( !strcmp(argv[1],"-ve")) {
        int choice=atoi(argv[2]);
		
		if (choice==1){
            int choice2=atoi(argv[3]);
            
            if (choice2==1){
                
                int dim=atoi(argv[4]), num=atoi(argv[5]), walk_len=atoi(argv[6]), i,j;
                std::vector<double> vec;
                double time, time2, c0, e=atof(argv[7]),vol, vol_sam;
                timer tim;
                
               // std::cout<<"hello"<<std::endl;
                std::ifstream inputFile2;
                std::string inputFileName2 = argv[8]; //reading the input txt file
                inputFile2.open(inputFileName2.c_str()); //converting the input to a c-style string
			
                std::istream_iterator< double >  input_begin2( inputFile2 );
                std::istream_iterator< double >  input_end2;
                for(std::istream_iterator< double > it2 = input_begin2; it2 != input_end2; ++it2){
                    vec.push_back(*it2);
                }
            
                MatrixXd C(dim,dim);
                VectorXd center(dim);
              //  std::cout<<"hello"<<std::endl;
                for (i=0; i<dim+1; i++){
                   // std::cout<<"hello"<<std::endl;
                    if (i<dim){
                        for (j=0; j<dim; j++){
                            C(i,j)=vec[j+i*dim];
                          //  std::cout<<C(i,j)<<std::endl;
                        }
                    }else{
                        for (j=0; j<dim; j++){
                            center(j)=vec[j+dim*dim];
                          //  std::cout<<"center "<<center(i)<<std::endl;
                        }
                        c0=vec[dim+dim*dim];
                       // std::cout<<c0<<std::endl;
                    }
                }
               // std::cout<<center<<std::endl;
               // std::cout<<c0<<std::endl;
                ellipsoids G(C, center, c0, dim);
                int rnum = std::pow(e,-2) * 400 * dim * std::log(dim);
                
                tim.reset();
                tim.start();
                vol=VolEsti_ellips2(G, walk_len, rnum);
                time=tim.time();
                
                tim.reset();
                tim.start();
                vol_sam=sample_cut_ellipsoid(dim, num, G);
                time2=tim.time();
                
                std::cout<<"\n";
                std::cout<<"VolEsti Volume: "<<vol<<std::endl;
                std::cout<<"VolEsti Excecutional time: "<<time<<std::endl;
                
                std::cout<<"Rejection Volume: "<<vol_sam<<std::endl;
                std::cout<<"Rejection Excecutional time: "<<time2<<std::endl;
                
                
            }else if (choice2==2){
                
                int dim=atoi(argv[4]), num=atoi(argv[5]), walk_len=atoi(argv[6]), i,j;
                std::vector<double> vec,cntr;
                double time, time2, c0, RR, e=atof(argv[7]),vol, vol_sam;
                timer tim;
                
                
                std::ifstream inputFile;
                std::string inputFileName = argv[9]; //reading the input txt file
                inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
                std::istream_iterator< double >  input_begin( inputFile );
                std::istream_iterator< double >  input_end;
                for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
                    vec.push_back(*it);
                }
			
                for (i=0; i<dim; i++){
                    if (i<dim){
                        cntr.push_back(vec[i]);
                        //std::vector<<cntr[i]<<"
                    }else{
                        RR=vec[i];
                    }
                }
                RR=vec[dim];
                vec.clear();
                std::ifstream inputFile2;
                std::string inputFileName2 = argv[8]; //reading the input txt file
                inputFile2.open(inputFileName2.c_str()); //converting the input to a c-style string
			
                std::istream_iterator< double >  input_begin2( inputFile2 );
                std::istream_iterator< double >  input_end2;
                for(std::istream_iterator< double > it2 = input_begin2; it2 != input_end2; ++it2){
                    vec.push_back(*it2);
                }
            
                MatrixXd C(dim,dim);
                VectorXd center(dim);
                
                for (i=0; i<dim+1; i++){
                    if (i<dim){
                        for (j=0; j<dim; j++){
                            C(i,j)=vec[j+i*dim];
                        }
                    }else{
                        for (j=0; j<dim; j++){
                            center(j)=vec[j+dim*dim];
                        }
                        c0=vec[dim+dim*dim];
                    }
                }
                
                ellipsoids G(C, center, c0, dim);
                int rnum = std::pow(e,-2) * 400 * dim * std::log(dim);
                
                tim.reset();
                tim.start();
                vol=VolEsti_ellips4(G, walk_len, rnum, cntr, RR);
                time=tim.time();
                
                tim.reset();
                tim.start();
                vol_sam=sample_cut_ellipsoid(dim, num, G);
                time2=tim.time();
                
                std::cout<<"\n";
                std::cout<<"VolEsti Volume: "<<vol<<std::endl;
                std::cout<<"VolEsti Excecutional time: "<<time<<std::endl;
                std::cout<<"\n";
                std::cout<<"Rejection Volume: "<<vol_sam<<std::endl;
                std::cout<<"Rejection Excecutional time: "<<time2<<std::endl;
                
            }
            
            
        }else if (choice==2){
            
            int dim=atoi(argv[3]), num=atoi(argv[4]), walk_len=atoi(argv[5]), i,j;
            std::vector<double> vec,pl;
            double time, c0, e=atof(argv[6]),vol,z1,z2,vol_sam,time2;
            timer tim;
                
                
            std::ifstream inputFile2;
            std::string inputFileName2 = argv[7]; //reading the input txt file
            inputFile2.open(inputFileName2.c_str()); //converting the input to a c-style string
			
            std::istream_iterator< double >  input_begin2( inputFile2 );
            std::istream_iterator< double >  input_end2;
            for(std::istream_iterator< double > it2 = input_begin2; it2 != input_end2; ++it2){
                vec.push_back(*it2);
            }
           
            MatrixXd C(dim,dim);
            VectorXd center(dim);
              
            for (i=0; i<dim+1; i++){
                if (i<dim){
                    for (j=0; j<dim; j++){
                        C(i,j)=vec[j+i*dim];
                    }
                }else{
                    for (j=0; j<dim; j++){
                        center(j)=vec[j+dim*dim];
                    }
                    c0=vec[dim+dim*dim];
                }
            }
                
            ellipsoids G(C, center, c0, dim);
            int rnum = std::pow(e,-2) * 400 * dim * std::log(dim);
                
            vec.clear();
            std::ifstream inputFile;
            std::string inputFileName = argv[8]; //reading the input txt file
            inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
            std::istream_iterator< double >  input_begin( inputFile );
            std::istream_iterator< double >  input_end;
            for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
                vec.push_back(*it);
            }
			
            for (i=0; i<dim+2; i++){
                if (i<dim){
                    pl.push_back(vec[i]);
                }
            }
            z1=vec[dim]; z2=vec[dim+1];
                
            tim.reset();
            tim.start();
            vol=VolEsti_ellips5(G, walk_len, rnum, pl, z1, z2);
            time=tim.time();
            
            tim.reset();
            tim.start();
            vol_sam=sample_cut_ellipsoid_and_hyps(dim, num, G, pl, z1, z2);
            time2=tim.time();
              
            std::cout<<"\n";
            std::cout<<"VolEsti Volume: "<<vol<<std::endl;
            std::cout<<"VolEsti Excecutional time: "<<time<<std::endl;
            std::cout<<"\n";
            std::cout<<"Rejection Volume: "<<vol_sam<<std::endl;
            std::cout<<"Rejection Excecutional time: "<<time2<<std::endl;
            
        }else if (choice==3){
            
            int dim=atoi(argv[3]), num=atoi(argv[4]), walk_len=atoi(argv[5]), i,j;
            std::vector<double> vec,pl;
            double time, time2, c1, c2, e=atof(argv[6]),vol,z1,z2,vol_sam;
            timer tim;
                
                
            std::ifstream inputFile2;
            std::string inputFileName2 = argv[7]; //reading the input txt file
            inputFile2.open(inputFileName2.c_str()); //converting the input to a c-style string
			
            std::istream_iterator< double >  input_begin2( inputFile2 );
            std::istream_iterator< double >  input_end2;
            for(std::istream_iterator< double > it2 = input_begin2; it2 != input_end2; ++it2){
                vec.push_back(*it2);
            }
           
            MatrixXd C(dim,dim);
            VectorXd center(dim);
           //   std::cout<<"hello"<<std::endl;
            for (i=0; i<dim+1; i++){
                if (i<dim){
                    for (j=0; j<dim; j++){
                        C(i,j)=vec[j+i*dim];
                    }
                }else{
                    for (j=0; j<dim; j++){
                        center(j)=vec[j+dim*dim];
                    }
                    c1=vec[dim+i*dim];
                    c2=vec[dim+i*dim+1];
                }
            }
          //  std::cout<<"hello"<<std::endl;
            ellipsoids G1(C, center, c1, dim);
            ellipsoids G2(C, center, c2, dim);
            int rnum = std::pow(e,-2) * 400 * dim * std::log(dim);
                
            vec.clear();
            std::ifstream inputFile;
            std::string inputFileName = argv[8]; //reading the input txt file
            inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
            std::istream_iterator< double >  input_begin( inputFile );
            std::istream_iterator< double >  input_end;
            for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
                vec.push_back(*it);
            }
			
            for (i=0; i<dim; i++){
                if (i<dim){
                    pl.push_back(vec[i]);
                }
            }
            z1=vec[dim]; z2=vec[dim+1];
          //      std::cout<<"hello"<<std::endl;
            tim.reset();
            tim.start();
            vol=VolEsti_ellips_nonConvex(G1, G2, walk_len, rnum, pl, z1, z2);
            time=tim.time();
            
            tim.reset();
            tim.start();
            vol_sam=sample_cut_ellipsoid_and_hyps_nonConv(dim, num, G1, G2, pl, z1, z2);
            time2=tim.time();
            
            std::cout<<"\n";
            std::cout<<"Volume: "<<vol<<std::endl;
            std::cout<<"Excecutional time: "<<time<<std::endl;
            std::cout<<"\n";
            std::cout<<"Rejection Volume: "<<vol_sam<<std::endl;
            std::cout<<"Rejection Excecutional time: "<<time2<<std::endl;
            
        }
        
        
    }else if( !strcmp(argv[1],"-varsi")){
        
        int dim=atoi(argv[2]),i;
        std::vector<double> vec,pl;
        double time,vol,z1;
        timer tim;
            
        std::ifstream inputFile;
        std::string inputFileName = argv[3]; //reading the input txt file
        inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
        std::istream_iterator< double >  input_begin( inputFile );
        std::istream_iterator< double >  input_end;
        for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
            vec.push_back(*it);
        }
			
        for (i=0; i<dim+2; i++){
            if (i<dim){
                pl.push_back(vec[i]);
            }
        }
        z1=vec[dim];
        
        tim.reset();
        tim.start();
        vol=vol_Ali(pl,-z1,dim);
        for (i=0; i<dim; i++){
            vol=vol/((double)(i+1));
        }
        time=tim.time();
        
        std::cout<<"\n";
        std::cout<<"Volume: "<<vol<<std::endl;
        std::cout<<"Excecutional time: "<<time<<std::endl;
        
    }else if( !strcmp(argv[1],"-lawn")){
        int choice2=atoi(argv[2]);
        
        if (choice2==1){
            int dim=atoi(argv[3]),i;
            std::vector<double> vec,pl;
            double time,z1,z2,vol1,vol2;
            timer tim;
            
            std::ifstream inputFile;
            std::string inputFileName = argv[4]; //reading the input txt file
            inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
            std::istream_iterator< double >  input_begin( inputFile );
            std::istream_iterator< double >  input_end;
            for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
                vec.push_back(*it);
            }
			
            for (i=0; i<dim+2; i++){
                if (i<dim){
                    pl.push_back(vec[i]);
                }
            }
            z1=vec[dim]; z2=vec[dim+1];
           // std::cout<<"z1: "<<z1<<" z2: "<<z2<<std::endl;
            tim.reset();
            tim.start();
            NT3 vol=Lawn1(pl, dim, z1, z2);
            for (i=0; i<dim; i++){
                vol=vol/NT3(i+1);
            }
            time=tim.time();
        
            vol1=vol_Ali(pl,-z1,dim);
            vol2=vol_Ali(pl,-z2,dim);
          //  std::cout<<"vol1 Varsi: "<<vol1<<" vol2 Varsi: "<<vol2<<std::endl;
            vol2=vol2-vol1;
            for (i=0; i<dim; i++){
                vol2=vol2/((double)(i+1));
            }
            
            std::cout<<"\n";
            std::cout<<"Lawrence Volume: "<<CGAL::to_double(vol)<<std::endl;
            std::cout<<"Lawrence excecutional time: "<<time<<std::endl;
            std::cout<<"\n";
            std::cout<<"Exact Volume with Varsi: "<<vol2<<std::endl;
            
        }else if (choice2==2){
            
            int dim=atoi(argv[3]),i;
            std::vector<double> vec,pl;
            double time,vol,z1,z2,vol1,vol2;
            timer tim;
            
            std::ifstream inputFile;
            std::string inputFileName = argv[4]; //reading the input txt file
            inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
			
            std::istream_iterator< double >  input_begin( inputFile );
            std::istream_iterator< double >  input_end;
            for(std::istream_iterator< double > it = input_begin; it != input_end; ++it){
                vec.push_back(*it);
            }
			
            for (i=0; i<dim+2; i++){
                if (i<dim){
                    pl.push_back(vec[i]);
                }
            }
            z1=vec[dim]; z2=vec[dim+1];
        
            tim.reset();
            tim.start();
            vol=Lawn2(pl, dim, z1, z2);
            for (i=0; i<dim; i++){
                vol=vol/((double)(i+1));
            }
            time=tim.time();
        
            vol1=vol_Ali(pl,-z1,dim);
            vol2=vol_Ali(pl,-z2,dim);
            vol2=vol2-vol1;
            for (i=0; i<dim; i++){
                vol2=vol2/((double)(i+1));
            }
			
			std::cout<<"\n";
            std::cout<<"Lawrence volume: "<<vol<<std::endl;
            std::cout<<"Lawrence excecutional time: "<<time<<std::endl;
            std::cout<<"\n";
            std::cout<<"Exact Volume with Varsi: "<<vol2<<std::endl;
            
        }
        
        
    }
	
	
	
	
	
	
	
	
	
	return 0;
	
}
