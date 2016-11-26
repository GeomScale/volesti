/** \example main_exp_par_auto.cpp
 * This is an example for experimenting.
 */

#include <string>
#include "../source/Auto_random_kd_forest.h"

void input_cmd(int argc, char *argv[], size_t& N, size_t& D, std::string& datafile,
					 size_t& Q, std::string& queryfile, int& k, double& epsilon,
					int& file_option, std::string& matchfile, std::string& datatype);
void print_cmd(char *argv[]);

//./rkd_sam 100000 960 ~/parallel/rkd_forest/Datasets/mygist/mygist_base.txt 100 ~/parallel/rkd_forest/Datasets/mygist/mygist_query.txt 1 0.0 0 ~/parallel/rkd_forest/Datasets/mygist/mygist_match.txt float

// ./rkd_sam 10000 128 ../../nicolo_paper/data/siftsmall_data.txt 100 ../../nicolo_paper/data/siftsmall_query.txt 2 0.0 0 ../../nicolo_paper/data/siftsmall_match.txt int
int main(int argc, char *argv[]) {
  size_t N, D, Q;
  int k, file_option;
  double epsilon;
  std::string datafile, queryfile, matchfile, datatype;
  std::vector<std::vector<std::pair<float, int> > > res;
	input_cmd(argc, argv, N, D, datafile, Q, queryfile, k, epsilon, file_option, matchfile, datatype);
	print_cmd(argv);

	if(datatype == "int") {
		Auto_random_kd_forest<int> RKDf(N, D, datafile, Q, queryfile, k, epsilon, res, NULL, file_option);
		RKDf.check_miss(matchfile, Q, k, res);
	} else {
		Auto_random_kd_forest<float> RKDf(N, D, datafile, Q, queryfile, k, epsilon, res, NULL, file_option);
		RKDf.check_miss(matchfile, Q, k, res);
	}
		
	std::cout << "Done" << std::endl;
	return 0;
}

void input_cmd(int argc, char *argv[], size_t& N, size_t& D, std::string& datafile,
					 size_t& Q, std::string& queryfile, int& k, double& epsilon,
					int& file_option, std::string& matchfile, std::string& datatype) {
	if(argc != 11) {
    std::cout << "Wrong number of cmd arguments!\n";
    std::cout << "Usage: " << argv[0] << " N D "
        "datafile Q queryfile k epsilon file_option matchfile datatype\n";
    std::cout << "Exiting...\n";
    exit(1);
  }				 
	
	N           = atoi(argv[1]);
	D						= atoi(argv[2]);
	datafile 		= argv[3];
	Q			    	= atoi(argv[4]);
	queryfile 	= argv[5];
	k						= atoi(argv[6]);
	epsilon			= atof(argv[7]);
	file_option = atoi(argv[8]);
	matchfile 	= argv[9];
	datatype		= argv[10];
}

void print_cmd(char *argv[]) {
  std::cout << "N = " << argv[1] << "\n";
  std::cout << "D = " << argv[2] << "\n";
  std::cout << "data_file = " << argv[3] << "\n";
  std::cout << "Q = " << argv[4] << "\n";
  std::cout << "query_file = " << argv[5] << "\n";
  std::cout << "k = " << argv[6] << "\n";
  std::cout << "epsilon = " << argv[7] << "\n";
  std::cout << "file_option = " << argv[8] << "\n";
  std::cout << "match_file = " << argv[9] << "\n";
  std::cout << "data_type = " << argv[10] << "\n";
}
