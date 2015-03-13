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


#include "virtualarray.h"
#include "olabuffer.h"
#include "counted_ptr.h"

template< class Container>
float longRangeSum(Container & data, int begin, int end) {
  float sum = 0.0f;
  for(int k = begin; k < end; ++k)
    sum += data[k];
  return sum;    
}

int main() {
  vector<float> data(9);
  for(uint k = 0; k < data.size(); ++k) data[k] =k;
  OlaBuffer<float> ob(2,1);
  counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
  int begin = 0, end =5;
  RangedCubicPolynomial rcp (1,0,0,0,begin,end);
  for(uint k = 0; k < buffer->size(); ++k)
    cout << (*buffer)[k] << " ";
  cout << endl;
  float fastway = ob.query(rcp,data,*buffer);
  float slowway = longRangeSum(data,begin,end);
  cout << " fastway = " << fastway << " longway = " << slowway << endl;
}
  



