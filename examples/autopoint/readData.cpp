#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define  MAXBUFSIZE ((int)1e6)

template <typename NT>
Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> readMatrix(const char * filename) {
    int cols=0, rows=0;
    NT buff[MAXBUFSIZE];

    ifstream infile;
    infile.open(filename);
    while(!infile.eof()) {
        string line;
        getline(infile,line);
        int temp_cols=0;
        stringstream stream(line);
        while(!stream.eof()) stream>>buff[cols*rows+temp_cols++];
        if(temp_cols==0) continue;
        if (cols==0) cols=temp_cols;
        rows++;
    }
    infile.close();
    Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> result(rows,cols);
    for(int i =0;i<rows;i++) {
        for (int j=0;j<cols;j++)
        {result(i,j)=(buff[cols*i+j]);
        }
    }
    return result;
    
}
