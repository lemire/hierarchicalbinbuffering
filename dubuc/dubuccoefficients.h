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

#ifndef DUBUCCOEFFICIENTS_H
#define DUBUCCOEFFICIENTS_H

//#include <functional>
#include <vector>
using namespace std;





  
  struct function : public unary_function<int,float> {
      virtual float operator() (const int l) const = 0;
      virtual ~function() {}
  };
  
  struct top: function {
      top(const int b, const int r) : mB(b), mR(r){}
      inline float operator() (const int l) const {return mR / (float) mB - l;}
      int mB, mR;
      virtual ~top() {}			
  };	

  struct bottom: function {
      bottom(const int m) : mM(m) {}
      inline float operator() (const int l) const {return mM - l;}
      int mM;
      virtual ~bottom() {}
  };





class DubucCoefficients {
  public:
    DubucCoefficients (const int b , const int N ):	
      mCoefficients(vector<vector<float> >(2 * N ,vector<float>(b,0.0f))), mN(N), mB(b) {
        for(int r = 0; r < b; ++r)  
          for(int m = 0; m < 2 * N; ++m ) {
            mCoefficients[m][r] = DDCoefs(m - N + 1,r);
          }
    }
    
    virtual ~DubucCoefficients() {}
    
    inline float coefficients(int m, int r ) const {
        assert(m +  mN - 1 >= 0);
        assert(m +  mN - 1 < (int) mCoefficients.size());
        assert(r >= 0);
        assert(r < mB);
        return mCoefficients[m + mN - 1][r];
    }    

    inline float leftCoefficients(int m, int r) const {
        return leftCoefs(m , r);
    }
    
  //  inline float rightCoefficients(int m, int r) {
   //     return leftCoefs(2 * mN - 1 - m, r);
   // }
     
    //const vector<vector<float> > & getCoefficients() const {return mCoefficients;}
    //const vector<vector<float> > & getLeftCoefficients() const {return mLeftCoefficients;}
    
  protected:
    vector<vector<float> > mCoefficients;
    vector<vector<float> > mLeftCoefficients;



    static inline float product(const function& f, 
        const int start, const int end, const int exclu) {
      float total = 1.0f;
      for(int i = start; i < exclu; ++i) 
        total *= f(i);
      for(int i = exclu + 1; i <= end; ++i) 
        total *= f(i);
      return total;
    }

    inline float leftCoefs( const int m, const int r) const {
      const float upper = product(top(mB,r),0,2*mN-1,m);
      const float lower = product(bottom(m),0,2*mN-1,m);
      return upper / lower;
    }

    inline float DDCoefs(const int m, const int r) const {
      const float upper = product(top(mB,r),1-mN,mN,m);
      const float lower = product(bottom(m),1-mN,mN,m);
      return upper / lower;
    }

    int mN, mB;


};


#endif


