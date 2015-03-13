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
#ifndef COMMON_H
#define COMMON_H


using namespace std;


// STL
//#include <stl.h>
#include <iterator>
#include <vector>
#include <map>
#include <deque>

//Standard
#include <sstream>
#include <iostream>
#include <cstdarg>
#include <fstream>
#include <exception>
#include <string>
#include <cstdio>
#include <cassert>
#include <cmath>

typedef unsigned long long uint64;

// own stuff
#include "molasseexception.h"
#include "optimizationoptions.h"
#include "mathutil.h"
#include "fileutil.h"

#ifndef null
#define null 0
#endif



#endif
