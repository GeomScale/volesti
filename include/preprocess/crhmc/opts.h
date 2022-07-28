#ifndef OPTS_H
#define OPTS_H
template <typename Type>  class opts {
public:
  int ipmMaxIter=200;
  Type ipmDistanceTol=1e-8;
  Type ipmDualTol=1e-12;
  int maxNZ=30;

};
#endif
