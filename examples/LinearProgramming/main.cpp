#include "InteriorPointMethod.h"
/*
We implement a simple Barrier+Newton Interior point method for linear programming.
Solving the min tc^T*x-sum_i(bi-Axi) proble for increasing values of t we can solve the original.
Using Newton's method allows multiplicative update of the parameter t
with the guarantee that we converge epsilon close to the solution.
*/
int main() {
  // Load Data
  // Polytope Ax<= b
  int m;
  cin >> m;
  int n;
  cin >> n;
  mat A = mat(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      cin >> A(i, j);
    }
  }
  vec b = vec(m);
  for (int i = 0; i < m; i++) {
    cin >> b(i);
  }
  vec c = vec(n);
  for (int i = 0; i < n; i++) {
    cin >> c(i);
  }
  vec point = vec(n);
  for (int i = 0; i < n; i++) {
    cin >> point(i);
  }
  // Apply Interior point method Barrier + Newton Method
  //IPMprimaldual(A, b, c, point);
  cout << "The optimum is: " << IPMprimaldual(A, b, c, point) << "\n";
  cout << "At the point: " << point.t();
  return 0;
}
