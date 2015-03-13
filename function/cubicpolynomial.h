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
#ifndef CUBICPOLYNOMIAL_H
#define CUBICPOLYNOMIAL_H

#include <functional>
#include <cassert>

using namespace std;






/*
 * This is just an object-function representing the polynomial
 * a0 + a1 * x + a2 * x**2 + a3 x **3
 */
class CubicPolynomial : public unary_function<int,float>  {
    public:
      
      CubicPolynomial(float a0, float a1, float a2, float a3):	mA0(a0),mA1(a1),mA2(a2),mA3(a3){}
      
      CubicPolynomial(const CubicPolynomial& CP): mA0(CP.mA0),mA1(CP.mA1),mA2(CP.mA2),mA3(CP.mA3)		{}
      
      virtual ~CubicPolynomial() {}
      
      static CubicPolynomial monome(const int degree) {// convenience method!!!
        assert(degree >= 0);
        assert(degree < 4);
        if(degree == 0) return CubicPolynomial(1,0,0,0);
        if(degree == 1) return CubicPolynomial(0,1,0,0);
        if(degree == 2) return CubicPolynomial(0,0,1,0);
        if(degree == 3) return CubicPolynomial(0,0,0,1);
        return CubicPolynomial(0,0,0,0);
      }

      CubicPolynomial& operator=(const CubicPolynomial& rhs) {
        this->mA0 = rhs.mA0;
        this->mA1 = rhs.mA1;
        this->mA2 = rhs.mA2;
        this->mA3 = rhs.mA3;
        return *this;
      }
      
      inline float operator()(const int& x) const { 
        // somewhat optimized
        float answer = mA0;
        float poly = x;
        answer += mA1 * poly;
        poly *= x;// = x*x
        answer += mA2 * poly;
        poly *= x;// = x*x*x
        answer += mA3 * poly;
        return answer;
      }
      
      float mA0,mA1, mA2, mA3;
};


class RangedFunction : public  unary_function<int,float> {
  public:
    RangedFunction(int Start, int End) : mStart(Start), mEnd(End) {}
    virtual float operator()(const int& x) const = 0;
    virtual ~RangedFunction() {}
    int mStart, mEnd;			
};

/*
 * This is an object-function representing the polynomial
 * a0 + a1 * x + a2 * x**2 + a3 x **3
 * over the range Start <= x < End, and zero elsewhere (x < Start, x >= End).
 */
class RangedCubicPolynomial : public CubicPolynomial, public RangedFunction  {
  public:
      RangedCubicPolynomial(float a0, float a1, float a2, float a3, int Start, int End): 
        CubicPolynomial(a0,a1,a2,a3), RangedFunction(Start,End) {}

      RangedCubicPolynomial(const CubicPolynomial& CP, int Start, int End):
        CubicPolynomial(CP), RangedFunction(Start,End) {}

      RangedCubicPolynomial(const RangedCubicPolynomial& CP):
        CubicPolynomial(CP), RangedFunction(CP.mStart, CP.mEnd) {}

      virtual ~RangedCubicPolynomial() {}
      
      inline float operator()(const int& x) const {
        if((x < mStart) || (x >= mEnd)) return 0.0f;	
        return CubicPolynomial::operator()(x);
      }
      static RangedCubicPolynomial monome(const int degree,const int start, const int end) {
        // convenience method!!!
        return RangedCubicPolynomial(CubicPolynomial::monome(degree),start,end);
      }

      RangedCubicPolynomial& operator=(const RangedCubicPolynomial& rhs) {
        (CubicPolynomial)(*this) = (const CubicPolynomial&) rhs;
        mStart = rhs.mStart;
        mEnd = rhs.mEnd;
        return *this;
      }
      
};

#endif

