// Copyright 2012 National and Kapodistrian University of Athens, Greece.
//
// This file is part of HeaDaCHe.
//
// HeaDaCHe is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// HeaDaCHe is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDaCHe,
// see <http://www.gnu.org/licenses/>.

#ifndef PRINT_FUNCTIONS_H
#define PRINT_FUNCTIONS_H

#include <iostream>
#include <vector>

template <class T>
std::ostream& operator<<(std::ostream& ost,const std::vector<T> &V){
  if(!V.size())
    return ost;
  typename std::vector<T>::const_iterator it=V.begin();
  std::cout<<(*it++);
  for(;it!=V.end();it++)
    ost<<","<<(*it);
  return ost;
}

#endif //PRINT_FUNCTIONS_H
