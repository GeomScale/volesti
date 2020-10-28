// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// This file is converted from BNMin1 (https://www.mrao.cam.ac.uk/~bn204/oof/bnmin1.html) by Apostolos Chalkis

// Original copyright notice:

/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2008

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2.

   \file bnmin_main.cxx

*/
#ifndef BNMIN_MAIN_H
#define BNMIN_MAIN_H

#include <string>
#include <stdexcept>

#include <boost/format.hpp>

//#include "bnmin_main1.h"
//#include "config.h"

//namespace Minim {

  inline const char * version(void)
  {
    //return PACKAGE_VERSION;
    return "11";
  }

    class BaseErr:
            public std::runtime_error
    {
    public:
        BaseErr(const std::string &s):
                std::runtime_error(s)
        {
        }

    };

    class NParsErr:
            public BaseErr
    {
    public:
        NParsErr(const std::string &fname,
                 size_t expected,
                 size_t received):
                BaseErr( (boost::format("In function %s expected %i but received %i pars ")
                          % fname
                          % expected
                          % received).str())
        {
        }


    };

  /*BaseErr::BaseErr(const std::string &s):
    std::runtime_error(s)
  {
  }

  NParsErr::NParsErr(const std::string &fname,
		     size_t expected,
		     size_t received):
    BaseErr( (boost::format("In function %s expected %i but received %i pars ") 
	      % fname
	      % expected
	      % received).str())
  {
  }*/
    

#endif

//}


