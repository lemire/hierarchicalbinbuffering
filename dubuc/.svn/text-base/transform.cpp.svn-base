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


#include "externalarray.h"
#include "olabuffer.h"
#include "counted_ptr.h"


/*
 * This file contains a bunch of regression tests
 * that must be run each time olabuffer.h is modified.
 * It basically insures, within reason, that the Ola
 * algorithm is properly implemented.
 */





/*
 * This class means that one of the regression tests failed.
 */
class TestFailedException {
  public:
    TestFailedException() {}
    TestFailedException(float f) {cout << "Exception : "<< f<< endl;}
};

void transformDeltas(int b, int N, int64 size, bool verbose = false) {
  if(verbose) cout << " Testing deltas b = "<< b << " N = " << N << " size = " << size << endl;
  OlaBuffer< float > ob(b,N);
  for(int64 k = 0; k < size; ++k) {
    vector<float> data(size,0);
    data[k] = 1;
    if(verbose) cout << " k = " << k << " size = " << size << endl;
    counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
    if(verbose) {
      for(uint j = 0; j < buffer->size(); ++j) cout << "b["<<j<< "]=" << (*buffer)[j]<< " ";
      cout << endl;
    }
    float total = 0.0;
    for (uint k = 0; k < buffer->size(); ++k)  {
      total += (*buffer)[k];
    }
    if( abs(total - 1.0f) > 0.00001f )  throw new TestFailedException(total);
  }
  if(verbose) cout << "    *Test succesful* " << endl; 
}



void checkErrorRange(int b, int N, int64 size, bool verbose = true) {
  OlaBuffer< float > ob(b,N);
  int64 begin = 0;
  for(int64 end = begin + 1 ; end <= size;++end) {
     int64 lower = size; int64 higher = 0;
     cout << " testing range from " << begin << " to " << end << endl;
     RangedCubicPolynomial rcp(1.2,0.9,0,0,begin,end);
     cout << " rcp " << endl;
     for(int64 k = 0; k < size; ++k) cout << rcp(k) << " " ;
     cout << endl;
     cout << " interpolate " << endl;
     for(int64 k = 0; k < size; ++k) cout << ob.interpolate(k,1,rcp,size) << " " ;
     cout << endl;
      for(int64 k = 0; k < size ; ++k) {
       if(abs(ob.interpolate(k,1,rcp,size) - rcp(k)) > 0.0001) {
         if(k < lower) lower = k;
         if(k > higher) higher = k;
       }
     }  
     ++higher ;
     cout << " range = ["<<lower << " , " << higher << " ] ( b = "<< b <<", N = "<< N 
       << ", size = "<< size << " )" << endl;
     cout << endl;
     ob.testImperfectRange(end,1,size);
  }  
}

void checkConstantInterpolation(int b, int N, int64 size, bool verbose = false) {
   OlaBuffer< float > ob(b,N);
   RangedCubicPolynomial rcp(1,0,0,0,0,size);
   for(int64 k = 0; k < size ; ++k) {
     if(abs(ob.interpolate(k,1,rcp,size)-1)> 0.00001) { 
       cout << "constant b = " << b << " N = "<< N << " size = " << size << endl;
       cout << "ob.interpolate("<< k<< ",b,rcp,size) = " << ob.interpolate(k,1,rcp,size) << endl;
       throw new TestFailedException(k);    
     }
   }
}


void checkLinearInterpolation(int b, int N, int64 size, bool verbose = false) {
   OlaBuffer< float > ob(b,N);
   RangedCubicPolynomial rcp(0,1,0,0,0,size);
   for(int64 k = 0; k < size ; ++k) {
     if(abs(ob.interpolate(k,1,rcp,size)-k)> 0.00001) {
       cout << "linear b = " << b << " N = "<< N << " size = " << size << endl;
       cout << "ob.interpolate("<< k<< ",b,rcp,size) = " << ob.interpolate(k,1,rcp,size) << endl;
       throw new TestFailedException(k);    
     }
   }
  
}


void rangeSums(int b, int N, int64 size, bool verbose = false) {
  if(verbose) cout << " Testing range sums b = "<< b << " N = " << N << " size = " << size << endl;
  OlaBuffer< float > ob(b,N);
  if (ob.computeRecommendedPaddedLength(size)  != size) 
    cout << " Failure in recommendedPaddedLength size = "<< size 
      << " recommended = "<< ob.computeRecommendedPaddedLength(size)  << endl;
  for(int64 k = 0; k < size; ++k) {
    vector<float> data(size,0);
    data[k] = 1;// **
    if(verbose) cout << "k = " << k << endl;
    counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
    if(verbose) {
      cout << "------------------" << endl;
      cout << "------------------" << endl;
    }
    if(verbose) {
      for(uint j = 0; j < buffer->size(); ++j) cout << "b["<<j<< "]=" << (*buffer)[j]<< " ";
      cout << endl;
    }
    for(int64 begin = 0; begin < size; ++begin) {
      for (int64 end = begin + 1 ; end <= size; ++end) {
        if(verbose) cout << endl << "======================="<< begin << " --> " << end << "(k="<< k 
          << " b="<< b<< " N="<< N<<", size = "<< size<<" )" << endl;
        RangedCubicPolynomial rcp(1,0,0,0,begin,end);
        if(verbose) {
          for(int64 j = begin; j < end; ++j) cout << " rcp(" << j<<") = " << rcp(j) ;
          cout << endl;
        }
        assert(rcp(begin) == 1.0); 
        float answer = ob.query(rcp , data , * buffer);
        if((begin <= k) && (k < end)) {
          if( abs(answer - 1.0f) > 0.00001f ) {
            cout << "k = "<< k << " is in range f("<<k << " ) = "<< rcp(k)<<endl;
            for(uint j = 0; j < buffer->size(); ++j) cout << "b["<<j<< "]=" << (*buffer)[j]<< " ";
            cout << endl;
            throw new TestFailedException(answer);    
          }
        } else {
          if( abs(answer ) > 0.00001f ) {
            cout << "k = "<< k << " is NOT in range f("<<k << " ) = "<< rcp(k)<<endl;
            for(uint j = 0; j < buffer->size(); ++j) cout << "b["<<j<< "]=" << (*buffer)[j]<< " ";
            cout << endl;
            throw new TestFailedException(answer);    
          }
        }
      }
    }
  } 
  if(verbose) cout << "    *Test succesful* " << endl; 
}

void rangeFirstMoments(int b, int N, int64 size, bool verbose = false) {
  if(verbose) cout << " Testing range first moments b = "<< b << " N = " << N << " size = " << size << endl;
  OlaBuffer< float > ob(b,N);
  for(int64 k = 0; k < size; ++k) {
    vector<float> data(size,0);
    data[k] = 1;// **
    if(verbose) cout << "k = " << k << endl;
    counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
    if(verbose) {
      cout << "------------------" << endl;
      cout << "------------------" << endl;
    }
    if(verbose) {
      for(uint j = 0; j < buffer->size(); ++j) cout << "b["<<j<< "]=" << (*buffer)[j]<< " ";
      cout << endl;
    }
    for(int64 begin = 0; begin < size; ++begin) {
      for (int64 end = begin + 1 ; end <= size; ++end) {
        if(verbose) cout << endl << "======================="<< begin << " --> " << end << "(k="<< k 
          << " b  ="<< b<< " N = "<< N<<", size = "<< size<<" )" << endl;
        RangedCubicPolynomial rcp(0,1,0,0,begin,end);
        if(verbose) {
          for(int64 j = begin; j < end; ++j) cout << " rcp(" << j<<") = " << rcp(j) ;
          cout << endl;
        }
        //assert(rcp(begin) == 1.0); 
        float answer = ob.query(rcp , data , * buffer);
        if((begin <= k) && (k < end)) {
          if( abs(answer - k) > 0.00001f ) {
            cout << "k = "<< k << " is in range f("<<k << " ) = "<< rcp(k)<<endl;
            for(uint j = 0; j < buffer->size(); ++j) cout << "b["<<j<< "]=" << (*buffer)[j]<< " ";
            cout << endl;
            throw new TestFailedException(answer);    
          }
        } else {
          if( abs(answer ) > 0.00001f ) {
            cout << "k = "<< k << " is NOT in range f("<<k << " ) = "<< rcp(k)<<endl;
            for(uint j = 0; j < buffer->size(); ++j) cout << "b["<<j<< "]=" << (*buffer)[j]<< " ";
            cout << endl;
            throw new TestFailedException(answer);    
          }
        }
      }
    }
  } 
  if(verbose) cout << "    *Test succesful* " << endl; 
}


void checkUpdate(int b, int N, int64 size, bool verbose = false) {
  if(verbose) cout << " Testing updates b = "<< b << " N = " << N << " size = " << size << endl;
  OlaBuffer< float > ob(b,N);
  cout << " beta (number of levels) = " << ob.levels(size) << endl;
  for(int64 k = 0 ; k < size ; ++k) {
    if(verbose) cout << " =======================" << endl;
    if(verbose) cout << " testing update of delta at k = " << k << " N = " << N << " size = " << size << endl;
    vector<float> data(size,0);
    data[k] = 1;// **
    if(verbose) cout << "k = " << k << endl;
    counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
    ob.updateBuffer(*buffer,k,-1);
    for(int i = (int) buffer->size() -1 ; i >= 0 ; --i) {
      if(abs((*buffer)[i]) > 0.00001) {
        cout << " buffer should have gone to zero " << endl;
        for (int l = 0; l < (int) buffer->size() ; ++l) cout << (*buffer)[l] << " ";
        cout << endl;
        throw TestFailedException((*buffer)[i]);
      }
    }
    counted_ptr<vector<float> > newbuffer = ob.computeBuffer(data);
    ob.updateBuffer(*buffer,k,1.0f);
    for(int i = 0; i < (int) buffer->size(); ++i) {
       if(abs((*buffer)[i] - (*newbuffer)[i]) > 0.00001) { 
          cout << " updating should be the same as tranforming " << endl;
          for (int l = 0; l < (int) buffer->size() ; ++l) cout << (*buffer)[l] << " ";
          cout << endl;
          for (int l = 0; l < (int) newbuffer->size() ; ++l) cout << (*newbuffer)[l] << " ";
          cout << endl;
          throw TestFailedException((*buffer)[i] - (*newbuffer)[i]);
       }
    }
   }
  if(verbose) cout << "    *Test succesful* " << endl; 
  } 



int main() {
  bool verbose = false;
  checkUpdate(2,1,5,verbose);
  checkUpdate(2,1,9,verbose);
  checkUpdate(2,1,17,verbose);
  checkUpdate(2,2,21,verbose);
  checkUpdate(2,3,21,verbose);
  checkUpdate(2,2,13,verbose);
  cout << "updates b = 2 ok " << endl;
  checkUpdate(4,1,5);
  checkUpdate(4,1,9);
  checkUpdate(4,2,13);
  cout << "updates b = 4 ok" << endl;
  rangeSums(2,1,5,verbose);
  rangeSums(2,1,9,verbose);
  rangeSums(2,2,13,verbose);
  rangeSums(2,3,21,verbose);
  cout << "range sums b = 2 ok " << endl;
  rangeFirstMoments(2,1,5);
  rangeFirstMoments(2,1,9);
  rangeFirstMoments(2,2,13);
  rangeFirstMoments(2,3,21);
  cout << "range first moment b = 2 ok" << endl;
  checkConstantInterpolation(4,1,9);
  checkConstantInterpolation(4,2,13);
  checkConstantInterpolation(4,3,21);
  cout << "constant ok" << endl;
  checkLinearInterpolation(4,1,9);
  checkLinearInterpolation(4,2,13);
  checkLinearInterpolation(4,3,21);
  cout << "linear ok" << endl;
  rangeFirstMoments(4,1,5);
  rangeFirstMoments(4,1,9);
  rangeFirstMoments(4,2,13);
  cout << "range first moment b = 4 ok" << endl;
  rangeSums(4,1,5);
  rangeSums(4,1,9);
  rangeSums(4,2,13);
  rangeSums(4,3,21);
  cout << "range sums b = 4 ok " << endl;
  transformDeltas(4,1,9);
  transformDeltas(4,2,13);
  cout << "deltas ok " << endl;
  cout << "If you made it that far, the code should be mostly bug free." << endl;
}





