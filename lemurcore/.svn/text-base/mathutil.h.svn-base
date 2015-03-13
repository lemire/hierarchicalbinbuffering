// Lemur OLAP library (c) 2003 National Research Council of Canada by Daniel Lemire, and Owen Kaser
 /**
 *  This program is free software; you can
 *  redistribute it and/or modify it under the terms of the GNU General Public
 *  License as published by the Free Software Foundation (version 2). This
 *  program is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 *  details. You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef MATHUTIL_H
#define MATHUTIL_H
#include "common.h"
/**
  *@author Daniel Lemire
  */

class MathUtil {
public:

/*  static bool increment( vector<int>& index, const vector<int>& bounds) {
    return increment(index,bounds,0u);
  }*/
  
  static int log(int Basis, int Number) {
    int answer = 1;
    while (Number > 1) {
      Number /= Basis;
      ++answer;
      if (Number == 1) return answer;
    }
    cerr << "WARNING: The value "<< Number<<" should be a power of "<< Basis<< ". "<<endl;
    return answer;
  }
  
  static int power(int Basis, int Number) { 
    int answer = 1;
    while (Number > 1) {
      answer *= Basis;
      --Number;
    }
    return answer;
  }
  
  // very inefficient?
  static bool increment( vector<int>& index, const vector<int>& start, const vector<int>& bounds, uint dim = 0) {
    if(++index[dim] >= bounds[dim]) {
      if(dim == index.size() - 1) return false;
      index[dim] = start[dim];
      return increment(index, start, bounds, ++dim);
    }
    return true;
  }

  // amount to calling increment offset times. Hopefully. Should be much faster. DL
  static bool add( const uint64 offset, vector<int>& index, 
      const vector<int>& start, const vector<int>& bounds, uint dim = 0) {
    index[dim] += offset;
    const int overflow = bounds[dim] - 1 - index[dim];
    if(overflow > 0) {
      if(dim == index.size() - 1) return false;
      index[dim] = start[dim];
      return add(overflow,index, start, bounds, ++dim);
    }
    return true;
  }
  
  // very inefficient?
  static bool increment( vector<map<int,float>::const_iterator>& index, 
      const vector<map<int,float>::const_iterator>& start, 
      const vector<map<int,float>::const_iterator>& bounds, uint dim = 0) {
    if(++index[dim] == bounds[dim]) {
      if(dim == index.size() - 1) return false;
      index[dim] = start[dim];
      return increment(index, start, bounds, ++dim);
    }
    return true;
  }


  static double L2_dist( vector<double> v1, vector<double> v2 ) {
    assert(v1.size() == v2.size());
    double acc = 0;
    for (uint i = 0; i < v1.size(); ++i)
      acc += (v1[i]-v2[i])*(v1[i]-v2[i]);
    return sqrt(acc);
  }


  MathUtil();
  ~MathUtil();
};


#endif

