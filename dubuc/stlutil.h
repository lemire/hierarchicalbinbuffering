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

#ifndef STLUTIL_H
#define STLUTIL_H

#include <vector>
#include <map>
#include "dubucbuffer.h"

using namespace std;

class STLUtil {
  public:
    /*inline static vector<DubucBuffer> createDubucBuffers(const vector<int>& basis) {
      vector<DubucBuffer> answer(basis.size());
      for(uint k = 0; k < basis.size(); ++k)  
        answer[k] = DubucBuffer(basis[k]);		
      return answer;
    }*/
   /* 
    inline static vector<LazyTransform> createLazyTransforms(const vector<int>& basis, const int N) {
      vector<LazyTransform> answer(basis.size(),N);
      for(uint k = 0; k < basis.size(); ++k) answer[k] = LazyTransform(basis[k],N);			
      return answer;
    }*/
    static void print (const map<int,float>& m) {
      cout << " map<int,float> (size = " << m.size() << ") " ;
      for(map<int,float>::const_iterator i = m.begin(); i != m.end(); ++i)
        cout << "["<< i->first << "] = " << i->second<< " ";
      cout << endl;
    }
    static void print (const vector<int>& m) {
      cout << "vector<int> (size = " << m.size() << ") ";
      for(uint k = 0; k < m.size(); ++k)
        cout << m[k] << " " ;
      cout << endl;
    }
};

#endif

