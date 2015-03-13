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

#ifndef OLABUFFER_H
#define OLABUFFER_H

#include <vector>
#include <cassert>
#include "counted_ptr.h"
#include "dubuccoefficients.h"
#include "cubicpolynomial.h"
#include <iostream>

typedef long long int64;
using namespace std;

/*
 * Ola: Online Computation of local moments
 * by Daniel Lemire, Ph.D. National Research Council of Canada, December 2001-2002.
 * Please foward questions to lemire@ondelette.com.
 * 
 *  This class can be used to compute really fast local moments. The computation has complexity O(log_b n+b)
 *  and the storage is only n/b. You can choose b = sqrt(n-1) for really cheap storage and O(sqrt(n)
 *  query time or a constant b for expensive storage (n/b) but logarithmic query time.
 *
 *  Example of how to use this class
 *
 *  
 *  int b = 2, N = 1 ; 
 *  OlaBuffer< float > ob(b,N);
 *  counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
 *  RangedCubicPolynomial rcp(0,1,0,0,begin,end);
 *  // we defined f(x) = x over the interval [begin, end)
 *  float answer = ob.query(rcp , data , * buffer); 
 *  // we computed really fast the scalar product <data,rcp>
 *
 *  There is one trick to worry about: the size of your data source cannot be arbitrary. To alleviate
 *  this problem, you need to pad you data (either really so, or virtually through wrapping).
 *  A method called recommendedPaddedLength will suggest to you a total size you can use.
 *  
 *  
 */
 
template < class DataType>
class OlaBuffer {

  public:

    enum{
      verboseTransformOnce = false, 
      verboseTransform = false, 
      verbose = false, 
      verboseInterpolate = false,
      validateRange = false,
      verboseUpdate = false
 
    };
    
    /*
     * The OlaBuffer object is used to precompute buffers and execute range queries (1D).
     * The constructor is cheap, only precompute Lagrange coefficients. You can use
     * the same OlaBuffer object for several different data sources, but you have to keep
     * track of your memory buffers.
     *
     * Typically way to use this object is to first compute the buffer (see method below)
     * and then call queries on it.
     * 
     * b is the basis (must be chosen carefully in light of the data source),
     * See the recommendedPaddedLength method
     *
     * N is the number of moments to buffer (N = 1, 2,...)
     *
     */
    OlaBuffer(int b, int N) : mB(b), mN(N), mDC(b,N) { assert(N>0); assert(b>1); }
  
    // this is thrown when a data stream is too small to be buffered, should never be thrown?
    class TooSmallException{
      public: TooSmallException() {}
    };
    // this is thrown if there is a mismatch between data size and the chosen basis
    class InvalidBasisVsDataSizeException {
      public: InvalidBasisVsDataSizeException() {}
    };
    // gets thrown only if validateRange is true and there is an error in computing the range
    class InvalidRangeException {
      public: InvalidRangeException() {}
    };



    

    /*
     * Using a buffer precomputed using the computeBuffer method and some data source,
     * compute the query given by the RangedFunction (the query is just the scalar product
     * of the Ranged Function with the data set). Note that for this to work, the 
     * RangedFunction must be a polynomial of degree at most mN.
     *
     * This should have log_b (data.size()) in complexity.
     */
    template <class Container>  
    float query(RangedFunction& f, const Container& data, vector<DataType>& buffer) 
       const throw(InvalidBasisVsDataSizeException){
   //   cout << " query with " << f.mStart << " to " << f.mEnd << endl;
      assert(f.mStart >=0);
      assert(f.mEnd >= f.mStart);
      assert((uint) f.mEnd <=  data.size());
      float sum = 0.0f;
      int scale = 1;
      pair<int64,int64> begin = imperfectRange(f.mStart,scale,data.size());
      pair<int64,int64> end = imperfectRange(f.mEnd,scale,data.size());
      if(begin.second > end.first) end.first = begin.second; // overlap
      for (int64 index = begin.first; index < begin.second; ++index) { 
        if(verbose) {
           cout << "(first) sum += (f("<<index
             <<") - interpolate("<<index<<","<<scale<<",f,"<<data.size()<<")) * data["<<index<<"];"<<endl;
             cout << " interpolate("<<index<<","<<scale<<",f,"<<data.size()<<") = " 
             <<  interpolate(index,scale,f,data.size())<< endl;
            cout << " data["<<index<<"] = " << data[index] << endl;
            cout << endl;
        }
        sum += (f(index) - interpolate(index,scale,f,data.size())) * data[index];
        if(verbose ) cout << " sum = " << sum << endl;
      }
      for (int64 index = end.first; index < end.second; ++index) {
        if(verbose) {
           cout << "(second) sum += (f("<<index
             <<") - interpolate("<<index<<","<<scale<<",f,"<<data.size()<<")) * data["<<index<<"];"<<endl;
             cout << " interpolate("<<index<<","<<scale<<",f,"<<data.size()<<") = " 
             <<  interpolate(index,scale,f,data.size())<< endl;
            cout << " data["<<index<<"] = " << data[index] << endl;
            cout << endl;
        }
         sum += (f(index) - interpolate(index,scale,f,data.size())) * data[index];
        if(verbose ) cout << " sum = " << sum << endl;
      }
      if(validateRange) {
        for(uint index = 0;  index < (uint64) data.size(); ++index) {
          if((index >= (uint) begin.first) && (index < (uint) begin.second)) continue;
          if((index >= (uint) end.first) && (index < (uint) end.second)) continue;
          if(abs(f(index) - interpolate(index,scale,f,data.size())) > 0.0001) {
            testImperfectRange(f.mStart,scale,data.size());
            testImperfectRange(f.mEnd,scale,data.size());
             cout << " bad range prediction!" <<endl;
            cout << " index = " << index << endl;
            cout << " interpolate("<<index<<","<<scale<<",f,"<<data.size()<<") = " 
             <<  interpolate(index,scale,f,data.size())<< endl;
            cout << " f( " << index << " ) = " << f(index) << endl;
            cout << " begin.first = " << begin.first << ", begin.second = " << begin.second << endl;
            cout << " end.first = " << end.first << ", end.second = " << end.second << endl;
            throw InvalidRangeException();
          }
        }
      }
      for (scale = mB; (mB*scale > 0) &&
         ( (uint64) data.size() / ( mB * scale) + 1 >=  (uint) 2 * mN ) ; scale *= mB) {
        if(verbose) cout << "*************8 intermediate scale = " << scale << " mB = " << mB << endl;
        if(data.size() % scale  != 1) throw InvalidBasisVsDataSizeException(); 
        pair<int64,int64> begin = imperfectRange(f.mStart,scale,data.size());
        pair<int64,int64> end = imperfectRange(f.mEnd,scale,data.size());
        if(begin.second > end.first) end.first = begin.second; // overlap
        for (int64 index = begin.first; index < begin.second; index+=scale) { 
          if(verbose) {
             cout << "(first) sum += (f("<<index
             <<") - interpolate("<<index<<","<<scale<<",f,"<<data.size()<<
             ")) * buffer["<<index/mB<<"];"<<endl;
             cout << " interpolate("<<index<<","<<scale<<",f,"<<data.size()<<") = " 
             <<  interpolate(index,scale,f,data.size())<< endl;
             cout << "f("<<index<<") = " << f(index) << endl;
             cout << endl;
          }
          sum += (f(index) - interpolate(index,scale,f,data.size())) * buffer[index/mB];
          if(verbose ) cout << " sum = " << sum << endl;
        }
        for (int64 index = end.first; index < end.second; index+=scale) {
          if(verbose) {
             cout << "(second) sum += (f("<<index
             <<") - interpolate("<<index<<","<<scale<<",f,"<<data.size()<<
             ")) * buffer["<<index/mB<<"];"<<endl;
             cout << " interpolate("<<index<<","<<scale<<",f,"<<data.size()<<") = " 
             <<  interpolate(index,scale,f,data.size())<< endl;
             cout << "f("<<index<<") = " << f(index) << endl;
             cout << endl;
          }
          sum += (f(index) - interpolate(index,scale,f,data.size())) * buffer[index/mB];
          if(verbose ) cout << " sum = " << sum << endl;
        }
        if(validateRange) {
          for(int64 index = 0; (uint64)index < (uint64) data.size(); index+=scale) {
            if((index >= begin.first) && (index < begin.second)) continue;
            if((index >= end.first) && (index < end.second)) continue;
            if(abs(f(index) - interpolate(index,scale,f,data.size())) > 0.0001)
              throw InvalidRangeException();
          }
        }
      }
      int64 blockbegin = f.mStart / scale * scale + (f.mStart % scale != 0 ? scale : 0);
      int64 blockend = f.mEnd / scale  * scale + 1;
      if(verbose) 
        cout << " scale = "<< scale << " blockbegin = " << blockbegin << " blockend = " << blockend << endl;
      for(int64 index = blockbegin; index < blockend ; index += scale) { 
        if(verbose) {
          cout << " index = " << index << endl;
          cout << " b[" << index<< " / "<< mB<< " ] = " << buffer[index/mB] << endl;
          cout << " f(index) = " << f(index) << endl;
        }
        sum += f(index) * buffer[index/mB];
        if(verbose ) cout << " sum = " << sum << endl;
        if(verbose)  cout << endl;
      }
      return sum;
    }
    
    /*
     * Compute the buffer which can be used by the query method.
     * You'd do that only once as it is expensive (linear complexity).
     * Make sure the size() reported by
     *  your container matches the recommendedPaddedLength 
     *  (see below for convenience method).
     */
    template <class Container>  
    counted_ptr<vector<DataType> >  computeBuffer (Container& data ) throw ( TooSmallException ) {
      if(verboseTransform) {
        for(int x = 1; ((uint64) data.size() / (mB * x) + 1 >= 2); x*=mB) {
          cout << " x = " << x << endl;
          cout << "residual = "<< ((uint64) data.size() / (mB * x) + 1) << endl; 
        }
      }
      int scale = 1;
      counted_ptr<vector<DataType> > buffer = transformOnce(data,scale);
      if(verboseTransform) {
        for (uint k = 0; k < buffer->size();++k) cout << " 1buf["<<k<<"] = "<< (* buffer)[k]<< " ";
        cout << endl;
      }
      for (scale = mB ; 
          (mB*scale > 0 ) && ((uint64) data.size() / (mB * scale) + 1 >= (uint) 2 * mN); scale *= mB) {
        if(verboseTransform) cout << " data.size() = " << data.size() 
          << " scale = " << scale << " mN = "<< mN << endl;
        counted_ptr<vector<DataType> > newbuffer  = transformOnce(*buffer, scale / mB);
        if(verboseTransform) {
          for (uint k = 0; k < newbuffer->size();++k) 
            cout << " new buf["<<k<<"] = "<< (* newbuffer)[k]<< " ";
          cout << endl;
        }
        for (uint k = 0; k < newbuffer->size();++k) (* buffer)[k * scale ] = (* newbuffer)[k];
        if(verboseTransform) {
          for (uint k = 0; k < buffer->size();++k) cout << " buf["<<k<<"] = "<< (* buffer)[k]<< " ";
          cout << endl;
        }
      }
      if(verboseTransform) cout << "scale = " << scale << endl;
      return buffer;
    }

    void updateBuffer(vector<DataType>& buffer, const int64 pos, const DataType change) {
      //cout << " updating! N= "<< mN << " b = "<< mB  << endl;
      map<int64, DataType> deltas;
      typename map<int64, DataType>::iterator iter;
      deltas[pos] = change;
      for (int64 scale = 1; (mB*scale <= 0 ) ||
        ( (( (int64) buffer.size() - 1 )*mB+1) / (mB * scale) + 1 >= 2 * mN) ;
        scale *= mB) {
          //cout << " scale = " << scale << endl;
          map<int64,DataType> olddeltas(deltas); // make a temp copy for the keys
          for(iter = olddeltas.begin(); iter != olddeltas.end(); iter++) {
              int64 index = iter->first;
              DataType value = iter->second;
              if( (index/scale)/mB * mB == index/scale ) continue;
              propagate(index, value, scale, deltas,buffer.size());              
              if( index /mB * mB == index) buffer[index/mB] += value;
              deltas.erase(index);
          }
          // print it out
          //cout << " printing content" << endl;
          //for(iter = deltas.begin(); iter != deltas.end(); iter++) {
           // cout << iter->first << ": " << iter->second  << endl;
          //}

      }
      //cout << " leftover ... " << endl;
      //for(iter = deltas.begin(); iter != deltas.end(); iter++) {
       //     cout << iter->first << ": " << iter->second  << endl;
      //}
      for(iter = deltas.begin(); iter != deltas.end(); iter++) {
        int64 index = iter->first;
        DataType value = iter->second;
        if( index /mB * mB == index) buffer[index/mB] += value;
        //buffer.erase(index);
      }
    }




    inline void propagate(int64 pos, DataType change, int64 scale, map<int64,DataType>& deltas, int buffer_size) {
        const int64 i = pos / scale;
        const int64 k = i / mB;
        const int r = i % mB;
        const int buffersize = (buffer_size - 1 ) / scale + 1;
        if ( k - mN + 1 < 0 ) { // left
            const int min = 0 , max = 2 * mN ;
            for(int m = min ; m < max ; ++m) {
              deltas[ m * scale * mB ] += mDC.leftCoefficients(m  ,i ) * change;
            }
        } else if (k + mN >= buffersize ) { // right
          const int min = 0 , max = 2 * mN;
          for(int m = min ; m < max ; ++m) {
            const int reversedi = (mB*(buffer_size-1))/scale  - i;
            deltas[(buffersize - 2*mN + m )  * scale *mB ] += 
                mDC.leftCoefficients(2 * mN - 1 - m, reversedi) * change;
            }
        } else { // middle
          const int min =  - mN + 1, max =  mN + 1 ;
          for(int m = min ; m < max ; ++m) {
              deltas[(k + m) * scale * mB] += mDC.coefficients(m, r ) * change;
            }
          }
       
    }



     /*
     * obselete (?) updateBuffer function
     */
    bool ___updateBuffer(vector<DataType> & buffer, const int64 pos, const DataType change, const int64 scale = 1) {
        if(
            (mB*scale <= 0 ) 
            ||
            ( (( (int) buffer.size() - 1 )*mB+1) / (mB * scale) + 1 < 2 * mN)
        ) return false; // we stop here
        assert(pos % scale == 0);
        const int buffersize = (buffer.size() - 1 ) / scale + 1;
        const int64 i = pos / scale;
        const int64 k = i / mB;
        assert(k >= 0); 
        assert(k < buffersize);
        const int r = i % mB;
        if( r == 0 ){
          if(! updateBuffer(buffer, k*mB*scale,change,scale*mB) ) {
            if(verboseUpdate)
              cout << "%%%%%  r = "<<r<< " scale = "<<scale
                <<" buffer[" << k*scale << "] += " << change << endl;
            buffer[k*scale] += change ;
          }
          return true;
        }
        if ( k - mN + 1 < 0 ) { // left
          if(verboseUpdate) cout << " left update " << endl;
          const int min = 0 , max = 2 * mN ;
          for(int m = min ; m < max ; ++m) {
            assert(m >= 0);
            assert(m * scale < (int) buffer.size());
            if(!updateBuffer(buffer, m*scale*mB,mDC.leftCoefficients(m  ,i ) * change,scale*mB)) {
              if(verboseUpdate) 
                cout << "(left)%buffer[" << m*scale << "] += " << 
                  mDC.leftCoefficients(m  ,i ) * change << endl;
              buffer[ m * scale ] += mDC.leftCoefficients(m  ,i ) * change;
            }
         }
         return false;
        } else if (k + mN >= buffersize ) { // right
          const int min = 0 , max = 2 * mN;
          for(int m = min ; m < max ; ++m) {
            assert(buffersize - 2*mN + m >= 0);
            assert((buffersize - 2*mN + m )  * scale < (int) buffer.size());
            const int reversedi = (mB*(buffer.size()-1))/scale  - i;
            if(! updateBuffer(buffer, (buffersize - 2*mN + m )  * scale * mB,
                mDC.leftCoefficients(2 * mN - 1 - m, reversedi) * change, scale*mB )) {
              if(verboseUpdate) cout << "(right)%% scale = "<<scale
              <<" buffer[" << (buffersize - 2*mN + m )  * scale 
              << "] += " 
              << mDC.leftCoefficients(2 * mN - 1 - m, reversedi) * change 
              << endl;
              buffer[(buffersize - 2*mN + m )  * scale  ] += 
                mDC.leftCoefficients(2 * mN - 1 - m, reversedi) * change;
            }
         }
         return false;
        } else { // middle
          const int min =  - mN + 1, max =  mN + 1 ;
          for(int m = min ; m < max ; ++m) {
            assert(k + m >= 0);
            assert((k+ m )*scale< (int) buffer.size());
            if( ! updateBuffer(buffer, (k + m) * scale*mB,mDC.coefficients(m, r ) * change,scale*mB) ) {
              if(verboseUpdate) cout << "(middle)%%%% scale = "<<scale<<" buffer[" << (k+m)*scale << "] += " 
              << mDC.coefficients(m, r ) * change << endl;
              buffer[(k + m) * scale] += mDC.coefficients(m, r ) * change;
            }
          }
          return false;
        }
    }

    
    /*
     * Using one out of every "scale" value, interpolate at current index.
     * Note that we need index/scale+1 == Length.
     *
     * You would typically not call this method, except for debugging purposes maybe or
     * if you want to solve an interpolation problem.
     */
    inline float interpolate(int64 index, int scale, RangedFunction& f, int64 Length) const {
    if(verboseInterpolate)  cout << "*********************interpolate " << scale << endl;
    assert(index >= 0);
    assert(index < Length);
    assert(scale > 0); 
    assert(scale < Length);
    assert((Length - 1) / scale * scale  == Length - 1); 
    int axis = index / (scale*mB) * (scale*mB) ;
    int r = ((index - axis) / scale) % mB;
    if(verboseInterpolate)
      cout << " index = " << index << " scale = " << scale << " axis = " << axis << " r = " << r << endl;
    assert ( r < mB);
    assert ( r >= 0);
    if( r == 0) return f(index);
    float answer = 0.0f;
    if(axis < ( mN  - 1 ) * mB * scale) { // left side
      for(int m = 0; m < 2 * mN;++m)  {
        int64 i = m * scale * mB; 
         if(verboseInterpolate) cout << " left side " << i <<  " mN = " 
          << mN << " axis = " << axis << " m = "
          << m << " Length = "<< Length << "index = " << index 
          << " f( "<<i <<") = " << f(i) 
          << "mDC.leftCoefficients("<<m << ","<<index/scale<<")="
          << mDC.leftCoefficients(m,index/scale) << endl;
        assert(i>= 0);
        assert( i< Length);
        answer += f( i) * mDC.leftCoefficients(m,index / scale);
      }
      return answer;
    }
    if(Length - index <= scale * mB * (mN - 1 )) {// right side
      for(int m = 0; m < 2 * mN;++m) {
        int64 i = m * scale * mB + Length - 1 - (2 * mN -1 ) * scale * mB ;
        if(verboseInterpolate)
          cout << " right side " << i <<  " mN = " << mN << " axis = " << axis << " m = " 
          << m << " Length = "<< Length << "index = " << index 
          << " f( "<<i <<") = " << f(i) 
          << "mDC.leftCoefficients("<< 2*mN - 1 - m<< ","<<(Length  - 1 - index) / scale <<")="<< 
          mDC.leftCoefficients(2*mN - 1 - m,(Length  - 1 - index) / scale )<< endl;
        assert(i >= 0);
        assert(i < Length);
        answer += f( i) * mDC.leftCoefficients(2*mN - 1 - m,(Length  - 1 - index) / scale /*mB  - r*/);
      }
      return answer;
    }
    // need to find neigbours
    for(int m =0 ; m < 2 * mN; ++m) {
        int64 i = (m-mN+1) * scale * mB + axis;  
        if(verboseInterpolate) cout << " middle " <<i << " Length = " << Length << 
          " f( " << i << ") = " << f(i) << " m = " << m << " mN = " << mN << " r = " << r <<
            endl;
        assert( i >= 0);
        assert( i < Length);
        answer +=  f(i) * mDC.coefficients(m-mN+1,r);
    }
    return answer ;
  }

 
    /*
     * How many "levels" can you expect in the hierarchical
     * structure. For example, with mN=1 and Length = 2, the answer is 1.
     *
     * The height of the tree for an array of length Length is defined
     * as the largest integer level such that
     *
     * Length/b^level +1 >= 2N 
     *
     * Convenience method.
     * 
     */
    int levels(const int64 Length) const {
        int levels = 0;
        int64 scale = mB;
        while ((Length / scale + 1 >= 2 *mN) && (Length / scale > 0) ){
          //cout << " scale = "<< scale << " Length = " << Length << " mN = " << mN << endl;
          scale *= mB; levels++; 
          //cout << " scale = "<< scale << " Length = " << Length << " mN = " << mN << endl;
          //cout << Length / scale + 1 << endl;
          if(scale <=0) break;
        }
        return levels;
    }
    
     /*
     * Sometimes, data sets come in a given length not of your choosing.
     * This method will suggest a "padded length" for you to use. Simply
     * use some wrapper to simulate a larger array, you can just padd with
     * whatever you like (such as zeroes).
     *
     * Convenience method.
     * 
     */
    int64 computeRecommendedPaddedLength(const int64 Length) const {
      // For padding purposes, returns the smallest allowed length 
      //equal or larger than the given argument
      const int steps = levels(Length);
      const int64 TransformRatio = power(steps);
      int leftover = (Length - 1) / TransformRatio;
      if ( (Length - 1) % TransformRatio > 0) ++leftover;
      return leftover * TransformRatio + 1;
    }

     // this determines the lower range and top of the lazy transform
    inline pair<int64,int64> testImperfectRange(const int64 x, const int scale, const int64 n) const {
      cout << " testing range code...******** x="<< x << " n = " << n << endl;
      assert(scale > 0);
      // first, some special cases...
      if(x == 0) return pair<int64,int64>(0,0); // no error 
      if(x == n) return pair<int64,int64>(n,n); // no error
      // next, we do the general case
      cout << " x/(scale *mB) = " << x/(scale*mB) << endl;
      cout << " x/(scale *mB) -mN + 1 = " << x/(scale*mB) - mN + 1 << endl;
      int64 lower = (x /( scale * mB )- mN ) * mB*scale + scale;
      int64 higher = ( x / ( scale * mB) + mN ) * mB * scale ;
      cout << " lower = " << lower << " higher = " << higher << endl;
      // next, if near left hand side, need to start from begining
      if( x  <= scale * mB * (2 * mN - 1) ) lower = scale;
      if ( n - x  < scale * mB * 2 * mN ) higher = n  - scale;
      // that should do it
      cout << "======================== predicted range : " << lower << " , " << higher << endl;
       assert(lower >= 0);
      assert(higher <= n);
      return pair<int64, int64>( 0, n);
//      return pair<int64,int64>(lower,higher);
    }

   
   
  protected:
    inline pair<int64,int64> imperfectRange(const int64 x, const int scale, const int64 n) const {
      assert(scale > 0);
      // first, some special cases...
      if(x == 0) return pair<int64,int64>(0,0); // no error 
      if(x == n) return pair<int64,int64>(n,n); // no error
      // next, we do the general case
      int64 lower = (x /( scale * mB )- mN ) * mB*scale + scale;
      int64 higher = ( x / ( scale * mB) + mN ) * mB * scale ;
      // next, if near left hand side, need to start from begining
      if( x  <= scale * mB * (2 * mN - 1) ) lower = scale;
      if ( n - x  < scale * mB * 2 * mN ) higher = n  - scale;
      // that should do it
       assert(lower >= 0);
      assert(higher <= n);
      return pair<int64,int64>(lower,higher);
    }
  
   
    // used to compute the transform
    template<class GenericContainer>
    counted_ptr<vector<DataType> > 
    transformOnce (GenericContainer& data, uint scale ) throw ( TooSmallException ) {
      int64 k;
      int min, max, r;
      const int buffersize = data.size() /( mB * scale) + 1;
      if(verboseTransformOnce) { 
        cout << " scale = " << scale << " data.size() = " << data.size() << " buffersize = " << buffersize
        << " mB = " << mB << endl;
        for(uint p = 0; p < data.size(); ++p) cout << " data["<<p<<"]= " << data[p] ;
        cout << endl;
      }
      if (buffersize < 2 * mN) throw TooSmallException();
      counted_ptr<vector<DataType> >  buffer ( new vector<DataType>(buffersize, 0));
      for (uint i = 0; i*scale < data.size() ; ++i) {
        if(verboseTransformOnce) 
          cout << " using data point data["<< i*scale << " ] = " << data[i*scale]<<endl;
        assert(i >=0); assert(i * scale < data.size());
        if(data[i * scale] == 0) continue;
        k = i / mB;
        assert(k >= 0); 
        assert((uint)k < buffer->size());
        r = i % mB;
        if( r == 0 ){
          (*buffer)[k] += data[i * scale];
          continue;
        }
        if ( k - mN + 1 < 0 ) { // left
          min = 0 ; max = 2 * mN ;
          for(int m = min ; m < max ; ++m) {
            (*buffer)[ m ] += mDC.leftCoefficients(m  ,i ) * data[ i * scale];
          }
        } else if (k + mN >= buffersize ) { // right
          min = 0 ; max = 2 * mN;
          for(int m = min ; m < max ; ++m) {
            (*buffer)[buffersize - 2*mN + m  ] += 
              mDC.leftCoefficients(2 * mN - 1 - m, data.size()/scale - 1 - i) * data[ i * scale];
          }
        } else { // middle
          min =  - mN + 1; max =  mN + 1 ; 
          for(int m = min ; m < max ; ++m) {
            assert(k + m >= 0);
            assert(k+ m < (int64) buffer->size());
            (*buffer)[k + m ] += mDC.coefficients(m, r ) * data[ i * scale];
          }
        }
      }
      return buffer;
    }

 
   inline int64 power(const int p) const {
      int64 answer = 1;
      for( int k = 0; k < p; ++k)
        answer *= mB;
      return answer;
    }
 

    int mB, mN;
    DubucCoefficients mDC;

};


#endif





