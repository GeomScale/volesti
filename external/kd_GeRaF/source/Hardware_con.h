/**
 @file Hardware_con.h
 Provides hardware info for parallel execution,
 if the compiler is old.
 */

#ifndef HARDWARE_CON_H
#define HARDWARE_CON_H

#include <sys/sysinfo.h>

  /**
   * \brief Returns number of concurrent threads
   *  supported by the implementation. The value
   *  should be considered only a hint. 
   *
   * @return - number of concurrent threads supported.
   * If the value is not well defined or not computable,
   * returns 0. 
   */
  unsigned int my_hardware_concurrency()
  {
    #if defined(PTW32_VERSION) || defined(__hpux)
      return pthread_num_processors_np();
    #elif defined(__APPLE__) || defined(__FreeBSD__)
      int count;
      size_t size=sizeof(count);
      return sysctlbyname("hw.ncpu",&count,&size,NULL,0)?0:count;
    #elif defined(BOOST_HAS_UNISTD_H) && defined(_SC_NPROCESSORS_ONLN)
      int const count=sysconf(_SC_NPROCESSORS_ONLN);
      return (count>0)?count:0;
    #elif defined(_GNU_SOURCE)
      return get_nprocs();
    #endif
    return 0;
  }
  
  /*
  windows version
  unsigned thread::hardware_concurrency()
  {
    SYSTEM_INFO info={{0}};
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors;
  }
  */

#endif /* HARDWARE_CON_H */
