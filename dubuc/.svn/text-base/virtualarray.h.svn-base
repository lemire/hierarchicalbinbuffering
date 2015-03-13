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

#ifndef VIRTUALARRAY_H
#define VIRTUALARRAY_H

#include <cassert>
#include <cmath>
#include <functional>
//#include "../lemurcore/common.h"
//#include "../lemurcore/fileutil.h"

using namespace std;

typedef unsigned long long uint64;



/*
 * This class can serve as a virtual array where the value of the components 
 * are given by some function.
 * This is useful to simulate megabytes worth of data without actually using
 * that much...
 * A safe copy constructor has not been implemented.
 * 
 */
template <class DataType, class Functor>
class VirtualArray {
  public:
    VirtualArray(uint64 size) : mArraySize(size), mF(){
    }

    virtual ~VirtualArray(){
    }
    const DataType  operator[](const uint64& pos) const {return mF(pos);}
    virtual uint64 size() const { return mArraySize; }

  protected:
    uint64 mArraySize; // array size
    Functor mF;
};


#endif
