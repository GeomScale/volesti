#ifndef OPTS_H
#define OPTS_H
template <typename Type> class opts {
public:
  const int ipmMaxIter = 200;
  const Type ipmDistanceTol = 1e-8;
  const Type ipmDualTol = 1e-12;
  int maxNZ = 30;
};
#endif
