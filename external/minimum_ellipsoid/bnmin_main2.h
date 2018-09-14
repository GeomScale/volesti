/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2008

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2.

   \file bnmin_main.cxx

*/
#include <boost/format.hpp>

#include "bnmin_main1.h"
//#include "config.h"

namespace Minim {

  const char * version(void)
  {
    //return PACKAGE_VERSION;
    return "11";
  }

  BaseErr::BaseErr(const std::string &s):
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
  }
    



}


