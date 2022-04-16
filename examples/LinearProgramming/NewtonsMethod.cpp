#include "InteriorPointMethod.h"

double threshold= 0.00001;
double e=0.000001;
void NewtonsMethod(mat& A,vec& b,vec& c,vec& p){
int max_iter=100;
int iter=0;
vec Grad=vec(A.n_cols);
Grad.fill(1);
mat H=mat(A.n_cols,A.n_cols);

while(iter<max_iter && norm(Grad)>threshold ){
vec temp= A*p;
//Compute the Gradient
//x_inverse.transform([](double val) {return (1.0/val);});

Grad=c;
for(int i=0;i< A.n_cols;i++){
for(int k=0;k< A.n_rows;k++){
//Grad(i)-=A(k,i)/max(b(k)-temp(k),threshold);
Grad(i)-=A(k,i)/max(b(k)-temp(k),e);
}
}
//cout<<"Grad=" <<Grad;
//Compute the Hessian matrix
H.fill(0);

for(int i=0;i<A.n_cols;i++){
for(int j=0;j<A.n_cols;j++){
for(int k=0;k<A.n_rows;k++){
//max((b(k)-temp(k))*(b(k)-temp(k)),threshold)
H(i,j)+=A(k,j)*A(k,i)/(max(b(k)-temp(k),e)*max(b(k)-temp(k),e));
//H(i,j)+=A(k,j)*A(k,i)/max((b(k)-temp(k))*(b(k)-temp(k)),e);

}
}
}
//cout<<"Grad= " <<Grad;
//cout<<"p= "<<p;
//update p
//p=p-Grad;
vec dx=-solve(H,Grad,solve_opts::fast);
double lambda=-as_scalar(Grad.t()*dx);
if(lambda/2<threshold){
  break;
}
//cout<<p.t();
p=p+dx;
iter++;
}

}
