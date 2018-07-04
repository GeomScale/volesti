/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2008

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2.

   \file bnmin_main.hxx

   The main include file for BNMin1 Library
   
   \mainpage A simple minimisation / inference library

*/

#include <string>
#include <stdexcept>

#ifndef __BNMIN_BNMIN_MAIN_HPP__
#define __BNMIN_BNMIN_MAIN_HPP__
namespace Minim {


  const char * version(void);

  /** \brief Base class for run-time errors within the library 
   */
  class BaseErr:
    public std::runtime_error
  {
  public:
    BaseErr(const std::string &s);

  };

  class NParsErr:
    public BaseErr
  {
  public:
    NParsErr(const std::string &fname,
	     size_t expected,
	     size_t received);
    

  };


}
#endif

