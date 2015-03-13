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
#include "externalarray.h"
#include "olabuffer.h"
#include "counted_ptr.h"

#include <climits>
#include <ctime>
#include <sys/resource.h>
#include <sys/time.h>

#ifdef DO_PAPI
extern "C" {
#include "papi.h"
}
#endif

using namespace std;
/*
 * Thrown when something goes wrong with Ola algorithm.
 */
class OlaBenchmarkException {
  public:
    OlaBenchmarkException() {}

};


/*
 * This file contains the benchmarking code
 * to determine how good ola is.
 */

template< class Container>
float longRangeSum(Container & data, int64 begin, int64 end) {
  float sum = 0.0f;
  for(int64 k = begin; k < end; ++k)
    sum += data[k];
  return sum;    
}
template< class Container>
float longFirstMoment(Container & data, int64 begin, int64 end) {
  float sum = 0.0f;
  for(int64 k = begin; k < end; ++k)
    sum += k* data[k];
  return sum;    
}
inline double diff(rusage start,rusage end) {
   return end.ru_utime.tv_sec+end.ru_utime.tv_usec*1e-6 -
                  (start.ru_utime.tv_sec+start.ru_utime.tv_usec*1e-6);
}

inline double diff( struct timeval start, struct timeval end) {
   return end.tv_sec+end.tv_usec*1e-6 -
                  start.tv_sec+start.tv_usec*1e-6;
}
   

inline double diff(clock_t start,clock_t end) {
  cerr << " do not use me, I wrap around at 72 minutes" << endl;
  return  (double)(end - start) / CLOCKS_PER_SEC;
}


 /*
 * Example of a function that can be used with VirtualArray
 */
template<class DataType>
class Sine: public unary_function<int64,DataType> {
  public:
    Sine() {}
    virtual ~Sine(){}
    DataType  operator()(const int64 x) const {return sin((double)x);}
};


/*
 * This is just in case someone had doubts about the validity of the Ola algorithm.
 */
void testRangeSums(int b, int N, int64 size, int MAXTRIALS=50000, bool verbose = false) throw(OlaBenchmarkException){
  VirtualArray< float, Sine<float> > data(size);
  OlaBuffer< float > ob(b,N);
  counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
  int effectivesize = size > INT_MAX ? INT_MAX : size;
  int number = 0;
  for(int begin = 0; begin < effectivesize; ++begin) {
    for (int end = begin + 1 ; end <= effectivesize; ++end) {
      if(verbose) cout << " testing range sum " << begin << " , " << end<< endl;
      RangedCubicPolynomial rcp(1,0,0,0,begin,end);
      float answer = ob.query(rcp , data , * buffer);
      float longanswer = longRangeSum(data,begin,end);
      ++number;
      if(abs(answer - longanswer) > 0.001) {
        throw OlaBenchmarkException();
      }
    }
    if(number > MAXTRIALS) break;
  }
  if(verbose) cout << " done testing range sums " << endl;
} 

vector<pair<int64,int64> > ranges(int MAXTRIALS, int64 size) {
    srand(423432434);
    vector<pair<int64,int64> > container;
    while(container.size() < (uint) MAXTRIALS) {
      int64 x1 = (int64)(rand()/((double)RAND_MAX)*size) ;
      int64 x2 = (int64)(rand()/((double)RAND_MAX)*size) ;
      if( x1 > x2) container.push_back(pair<int64,int64>(x2,x1));
      else container.push_back(pair<int64,int64>(x1,x2));
    }
    return container;    
}

// Daniel starts to weep: macros/conditionals
// do not despair: nobody has to maintain this....

#ifdef USECLK
  #define TIMER(x)    x = clock()
  #define DECL_TIMER(x) clock_t x
#else
  #ifdef USE_EXTERNAL
     #define TIMER(x) gettimeofday(&x, 0)
     #define DECL_TIMER(x) struct timeval x
  #else
     #define TIMER(x) getrusage(RUSAGE_SELF, & x);
     #define DECL_TIMER(x) struct rusage x
  #endif
#endif


pair<double,double> updates(int b, int N, int64 size, int MAXTRIALS=50000 , bool verbose = false) {
  if(verbose) 
    cout << " == Updates === virtual array of size "<< size <<" N = " << N << " b = " << b << endl;
  VirtualArray< float, Sine<float> > data(size);
  DECL_TIMER(start); DECL_TIMER(end); TIMER(start);
  OlaBuffer< float > ob(b,N);
  counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
  TIMER(end);
  double Init =  diff (start,end);
  if(verbose) 
    cout << " It took " << Init<< " s to build a buffer of size " << buffer->size() << endl;
  int effectivesize = size > INT_MAX ? INT_MAX : size;
  if(verbose) cout << " beta (number of levels) = " << ob.levels(size) << endl;
  srand(432512); // fix seed
  TIMER(start);
  for(int k = 0 ; k < MAXTRIALS; ++k ) {
    int x = (int)(rand()/((double)RAND_MAX)*effectivesize) ;
    float change = 1.0;//(rand()- RAND_MAX/2.0f)/((float)RAND_MAX); // doesn't matter
    ob.updateBuffer(*buffer,x, change);
  }
  TIMER(end);
  double NombreDeSecondes =  diff(start,end);
  if(verbose) { 
    cout << " [update] Computations took " << NombreDeSecondes << endl;
    cout << " For " << MAXTRIALS << " range sums " << endl;
  }
  return pair<double,double>(Init,NombreDeSecondes);

}


#ifdef DO_PAPI

float rstart, pstart, flips;
long_long fp_issues_start;

#endif

void maybe_start_timing(void) {
#ifdef DO_PAPI
  PAPI_flips(&rstart,&pstart,&fp_issues_start, &flips);
#endif
}


void maybe_stop_timing(void) {
#ifdef DO_PAPI

  float rend, pend, flips1;
  long_long fp_issues_end;

  PAPI_flips(&rend,&pend,&fp_issues_end, &flips1);
  if (true) 
    cout << (fp_issues_end - fp_issues_start) << "\t";
  cout << (pend - pstart) << endl;

#endif
}

/*
 * Compute and benchmark Ola-based range sums
 */
pair<double,double> fastRangeSums(int b, int N, int64 size, int MAXTRIALS=50000 , bool verbose = false) {
  if(verbose) 
    cout << " == Fast Range Sums === virtual array of size "<< size <<" N = " << N << " b = " << b << endl;
#ifdef USE_EXTERNAL
  // let's not use virtual base classes or templstes; this is a temp hack
  ExternalArray<float> data(size,"OwenFile.deleteme");  // preinited data[i]=sin(i)
#else
  VirtualArray< float, Sine<float> > data(size);
#endif
  if(verbose) cout <<" b = " << b << " N = " << N << endl;
  if(verbose) cout <<" data.size() = " << data.size() << endl;

  DECL_TIMER(start); DECL_TIMER(end); TIMER(start);
  OlaBuffer< float > ob(b,N);
  counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
  TIMER(end);
  double Init =  diff(start,end);
  if(verbose) 
    cout << " It took " << Init<< " s to build a buffer of size " << buffer->size() << endl;
  int effectivesize = size > INT_MAX ? INT_MAX : size;
  vector<pair<int64,int64> > container = ranges(MAXTRIALS, effectivesize);
  TIMER(start);
  float average = 0.0;
#ifdef DO_PAPI
  cout << " rangelen  time " << endl << "---------------------" << endl;
#endif 
  for(vector<pair<int64,int64> >::iterator iter = container.begin();
      iter != container.end(); ++iter) {
      int begin = iter->first;
      int end = iter->second;
#ifdef DO_PAPI
      cout << iter->second - iter->first << " ";
#endif      
      maybe_start_timing();
      RangedCubicPolynomial rcp(1,0,0,0,begin,end);
      float answer = ob.query(rcp , data , * buffer);
      maybe_stop_timing();
      average += answer;
  }
  TIMER(end);
  double NombreDeSecondes =  diff(start,end);
  if(verbose) { 
    cout << " [frs] Computations took " << NombreDeSecondes << endl;
    cout << " For " << MAXTRIALS << " range sums " << endl;
    cout << " average was " << average / MAXTRIALS << endl; 
  }
  return pair<double,double>(Init,NombreDeSecondes);
}

/*
 * First moments the fast way....
 */
pair<double,double> fastFirstMoments(int b, int N, int64 size, int MAXTRIALS=50000, bool verbose=false ) {
  if(verbose)
      cout << " == Fast First Moments === virtual array of size "
      << size <<" N = " << N << " b = " << b << endl;
  VirtualArray< float, Sine<float> > data(size);
  DECL_TIMER(start); DECL_TIMER(end); TIMER(start);
  OlaBuffer< float > ob(b,N);
  counted_ptr<vector<float> > buffer = ob.computeBuffer(data);
  TIMER(end);
  double Init =  diff(start,end);
  if(verbose) 
    cout << " It took " << Init<< " s to build a buffer of size " << buffer->size() << endl;
  int effectivesize = size > INT_MAX ? INT_MAX : size;
  vector<pair<int64,int64> > container = ranges(MAXTRIALS, effectivesize);
  TIMER(start);
  float average = 0.0;
#ifdef DO_PAPI
  cout << " rangelen  time " << endl << "---------------------" << endl;
#endif 
  for(vector<pair<int64,int64> >::iterator iter = container.begin();
      iter != container.end(); ++iter) {
      int begin = iter->first;
      int end = iter->second;
#ifdef DO_PAPI
      cout << iter->second - iter->first << " ";
#endif      
      maybe_start_timing();
      RangedCubicPolynomial rcp(0,1,0,0,begin,end);
      float answer = ob.query(rcp , data , * buffer);
      maybe_stop_timing();
      average += answer;
  }
  TIMER(end);
  double NombreDeSecondes = diff(start,end);
  if(verbose) { 
    cout << " [ffm] Computations took " << NombreDeSecondes << endl;
    cout << " For " << MAXTRIALS << " range sums " << endl;
    cout << " average was " << average / MAXTRIALS << endl; 
  }
  return pair<double,double>(Init,NombreDeSecondes);
}
/*
 * Doing range sums the loooooonnnnng way.
 */
pair<double,double> slowRangeSums(int64 size, int MAXTRIALS=50000, bool verbose =false ) {
  if(verbose) cout << " == Slow Range Sums === virtual array of size "<< size << endl;

#ifdef USE_EXTERNAL
  // let's not use virtual base classes or templstes; this is a temp hack
  ExternalArray<float> data(size,"OwenFile.deleteme");  // preinited data[i]=sin(i)
#else
  VirtualArray< float, Sine<float> > data(size);
#endif
  //int effectivesize = size > INT_MAX ? INT_MAX : size;
  int64 effectivesize = size;
  vector<pair<int64,int64> > container = ranges(MAXTRIALS, effectivesize);
  int64 totaltime = 0;
  float average = 0.0;
  for(vector<pair<int64,int64> >::iterator iter = container.begin();
      iter != container.end(); ++iter) {
      int64 beginr = iter->first;
      int64 endr = iter->second;
      DECL_TIMER(start); DECL_TIMER(end); TIMER(start);
      float answer =  longRangeSum(data,beginr,endr);
      TIMER(end);
      totaltime += (int64) diff(start,end);
      average += answer;
  }
  double NombreDeSecondes = ((double)totaltime) / CLOCKS_PER_SEC;
  if(verbose) { 
    cout << " [srs] Computations took " << NombreDeSecondes << endl;
    cout << " For " << MAXTRIALS << " range sums " << endl;
    cout << " average was " << average / MAXTRIALS << endl; 
  }
  return pair<double,double>(0,NombreDeSecondes);
}/*
 * First moments the long way...
 */
pair<double,double> slowFirstMoments(int64 size, int MAXTRIALS=50000 , bool verbose=false) {
  if(verbose) cout << " == Slow First Moments === virtual array of size "<< size << endl;
  VirtualArray< float, Sine<float> > data(size);
  //int effectivesize = size > INT_MAX ? INT_MAX : size;
  int64 effectivesize = size;
  vector<pair<int64,int64> > container = ranges(MAXTRIALS, effectivesize);
  int64 totaltime = 0;
  float average = 0.0;
  for(vector<pair<int64,int64> >::iterator iter = container.begin();
      iter != container.end(); ++iter) {
      int64 beginr = iter->first;
      int64 endr = iter->second;
      DECL_TIMER(start); DECL_TIMER(end); TIMER(start);
      float answer =  longFirstMoment(data,beginr,endr);
      TIMER(end);
      totaltime += (int64) diff(start,end);
      average += answer;
  }
  double NombreDeSecondes = ((double)totaltime) / CLOCKS_PER_SEC;
  if(verbose) { 
    cout << " [sfm] Computations took " << NombreDeSecondes << endl;
    cout << " For " << MAXTRIALS << " range sums " << endl;
    cout << " average was " << average / MAXTRIALS << endl; 
  } 
  return pair<double,double>(0.0,NombreDeSecondes);
}


typedef vector<pair<double,double> > row;

void print(row  currentrow) {
      for(row::iterator celliter = currentrow.begin(); celliter != currentrow.end(); ++celliter) {
        cout << celliter->first << "," << celliter->second ;
        if(celliter + 1 != currentrow.end()) cout << " , ";
      }
      cout << endl;
      cout.flush();
}

void reportParams(int N, int b, int64 size, OlaBuffer<float> &ob)
{
  cout << " size = " << size << " b="<< b<<",beta=" << ob.levels(size) << endl;
  cout << " external array is : " << 4*size/ (1024*1024.0*1024.0)<< " GB" <<endl;
  cout << " bin arrays are : " << 4*b/ (1024.0)<< " KB" <<endl;
  cout << "N (aka M and M' in paper) is " << N <<endl;
}

void run(int N, int b, bool firstMomentsToo=false, int64 size =(1LL << 30)+1) {
  N /= 2;   // paper M and M' match code's N
  OlaBuffer< float > ob(b,N);
  int small = 200;  // ofk hack
  reportParams(N,b,size,ob);
  
  if (N*N*b > 100000) small = small * 100000 / (N*N*b);
  row currentrow;
  cout <<" query time (sums) / buffer construction " << endl;
  currentrow.push_back(fastRangeSums(b,N,size,small, true)); 
  if (firstMomentsToo)
    currentrow.push_back(fastFirstMoments(b,N,size,small, true)); 
  print(currentrow); 
}

void runUpdates(int N, int b, int64 size =(1LL << 30)+1) {
  N /= 2;   // paper
  OlaBuffer< float > ob(b,N);
  reportParams(N,b,size,ob);
  row currentrow;
  if (false) {
    int small = 2000;
    cout << "update time for  " << small << " updates" << endl;
    currentrow.push_back(updates(b,N,size,small, true));
    cout.flush();
  }
  int big = 200000;
  cout << "update time for  " << big << " updates" << endl;
  currentrow.push_back(updates(b,N,size,big, true)); 
  print(currentrow); 
}

 
int main() {

  // suitable b's for n= 4B
  int bValues [] = {1<<5, 1<<7, 1<<10, 1<<15, 1<<20, -1};
  int bTypical = bValues[1];
  int bSmallish = bValues[0];
  int bBiggish = bValues[3];

  // suitable Ns
  int nValues [] = {2,4,8,16,-1};
  int nTypical = 4;
 

  cout << " Benchmarking tool for the Ola algorithm " << endl;
  cout << " Build: " << __DATE__ << " " << __TIME__ << endl;
#ifdef USE_EXTERNAL
  cout << " Times are WALL CLOCK" << endl;
#else
  cout << " Times are user; ensure no significant system times" << endl;
#endif
  cout << " Benchmarks can take a long time to compute "<< endl;
  cout << " Values are always given as pairs init time + total operation time " << endl;

  /* sections of tests to disable/enable */

  bool doRangeTests = false,
    doVaryBonTypicalN =false,
    doVaryNonTypicalB = false,
    doVaryNonSmallB = false,
    doVaryNonBigB = false,
    doPapiTest1 = false,//ignore warnings about this being unused
    doPapiTest2 = false,//ignore warnings about this being unused
    doTestTypicalUpdate = false,
    doSmallerExternalTest = false,
    doSmallerNaiveTest = false,
    doUpdatesVsb = false,
    doUpdatesVsN = false,
    doNaiveSumTest = false,
    doRepeatedb128Small = false;

#ifdef USE_EXTERNAL
   doSmallerExternalTest = true;
   doSmallerNaiveTest = false;
  // doRepeatedb128Small = true;
#endif

#ifdef DO_PAPI
   doPapiTest1 = true;
   doPapiTest2 = true;
#endif

  if(doRangeTests) {// not necessary
    testRangeSums(4,1,9,true);
    testRangeSums(4,1,1025,true);
  }
  if(true) {
    cout << "Starting benchmarking from hellifax "<<endl;

    if (doVaryBonTypicalN) {
     cout << "testing effects of varying b" << endl;
    for(int bidx=0; bValues[bidx] != -1; ++bidx)
      run(nTypical,bValues[bidx]);
    }

    if (doVaryNonTypicalB) {
    cout << "testing effects of varying N" << endl;
    for(int Nidx=0; nValues[Nidx] != -1; ++Nidx)
      run(nValues[Nidx],bTypical);
    }

    if (doVaryNonSmallB) {
    cout << "testing effects of varying N on small b" << endl;
    for(int Nidx=0; nValues[Nidx] != -1; ++Nidx) 
      run(nValues[Nidx],bSmallish);
    }

    if (doVaryNonBigB) {
    cout << "testing effects of varying N on biggish b" << endl;
    for(int Nidx=0; nValues[Nidx] != -1; ++Nidx)
      run(nValues[Nidx],bBiggish);
    }

#ifdef DO_PAPI
    if (doPapiTest1) {
      cout << "testing dependence of Ola time on range length (multi-level)" << endl;
      run(nTypical,bTypical,true);
    }

    if (doPapiTest2) {
      cout << "testing dependence of Ola time on range length (single-level)" << endl;
      run(nTypical,bBiggish,true);
    }
#endif //DO_PAPI

    if (doTestTypicalUpdate) {
      cout << "Testing incremental updates with typical N and b" << endl;
      runUpdates(nTypical,bTypical);
    }

    if (doSmallerExternalTest) {
      cout << "Testing multilevel Ola on external array (multi-level)" << endl;
      cout << "Trying various b values" << endl;
      int64 smaller_n = (1LL<<28)+1;
      for(int bidx=0; bValues[bidx] != -1; ++bidx)
	run(nTypical,bValues[bidx], false, smaller_n);
    }


    if (doRepeatedb128Small) {
      cout << "Testing multilevel Ola on external array b=128" << endl;
      cout << "Trying various b a few times" << endl;
      int64 smaller_n = (1LL<<28)+1;
      for (int ii=0; ii < 10; ++ii)
      for(int bidx=0; bValues[bidx] != -1; ++bidx)
	run(nTypical,bValues[bidx], false, smaller_n);

    }

    if (doSmallerNaiveTest) {
      cout << "Testing obvious no-precomputation algorithm, smaller data" << 
	endl;
      int small = 20;
      row currentrow;
      int64 smaller_n = (1LL<<28)+1;
      cout << "time for  " << small << " naive rangesums" << endl;
      cout << "data size is "<< smaller_n << endl;
      currentrow.push_back(slowRangeSums(smaller_n, small, true ));
      print(currentrow);
    } 

    if (doNaiveSumTest) {
      cout << "Testing obvious no-precomputation algorithm, full data" << 
	endl;
      int small = 100;
      row currentrow;
      int64  n = (1LL<<30)+1;
      cout << "time for  " << small << " naive rangesums" << endl;
      cout << "data size is "<< n << endl;
      currentrow.push_back(slowRangeSums(n, small, true ));
      print(currentrow);
    } 





    if (doUpdatesVsb) {
      cout << "Testing incremental updates vs b, N=2" << endl;
      for(int bidx=0; bValues[bidx] != -1; ++bidx)
	runUpdates(nTypical,bValues[bidx]);
    }


    if (doUpdatesVsN) {
      cout << "Testing incremental updates vs N, b=128" << endl;
      for(int nidx=0; nValues[nidx] != -1; ++nidx)
        runUpdates(nValues[nidx],bTypical);
    }







    cout << "Done with benchmarking from hellifax."<<endl;
    exit(0);

  }

  // heritage stuff

  for(int N = 1; N <4; ++N) {
    int small = 2000;
    uint maxx = 4096;//512;
    int large = 20000;
    // test single scale ola
    //table data;
    cout << " Slow computations = " << small << " N = " << N << endl;
    cout << " x / size=x^2+1 / frs / sfm " << endl;
     for(int64 x = 64; x <= maxx; x*=8) {
      row currentrow;
      currentrow.push_back(pair<double,double>(x,x*x+1));
      currentrow.push_back(slowRangeSums(x*x+1,small)); 
      currentrow.push_back(slowFirstMoments(x*x+1,small)); 
      print(currentrow);
//      data.push_back(currentrow);
    }
     for(int64 x = 64; x <= maxx; x*=8) {
      int64 size = x*x+1;
      assert(size >0);
      OlaBuffer< float > ob8(8,N);
      OlaBuffer< float > ob64(64,N);
      OlaBuffer< float > obx(x,N);
      cout << " size = " << size << " b=8,beta=" << ob8.levels(size) 
        << " b=64,beta=" << ob64.levels(size) << " b=" <<x << ",beta=" << obx.levels(size) <<endl;
    }
    cout <<  endl;
    cout << " Updates for Ola  N = " << N <<  endl;
    cout << " x/size=x^2+1 / updates 8 / upt 64,  sampling = " << large << endl;
    for(int64 x = 64; x <= maxx; x*=8) {
      int64 size = x*x+1;
      row currentrow;
      currentrow.push_back(pair<double,double>(x,size));
      currentrow.push_back(updates(8,N,size,large));
      currentrow.push_back(updates(64,N,size,large));
      print(currentrow);
     // data.push_back(currentrow);
    }
   //printTable(data);
    //data.clear();
    cout << " Multi-Scale Ola size = " << small << " N = " << N << endl;
    cout << " x/size=x^2+1 / frs8/frs64/ffm8/ffm64  " << endl;
    for(int64 x = 64; x<=maxx; x*=8) {
      int64 size = x*x+1;
      row currentrow;
      currentrow.push_back(pair<double,double>(x,size));
      currentrow.push_back(fastRangeSums(8,N,size,small));
      currentrow.push_back(fastRangeSums(64,N,size, small));
      currentrow.push_back(fastFirstMoments(8,N,size,small));
      currentrow.push_back(fastFirstMoments(64,N,size,small));
      print(currentrow);
      //data.push_back(currentrow);
    }
   //printTable(data);
    //data.clear();
    cout << " Single Scale Ola size = " << small << " N = " << N << endl;
    cout << " x / size=x^2+1 / frs / ffm " << endl;
    for(int64 x = 64; x<=maxx; x*=8) {
      row currentrow;
      currentrow.push_back(pair<double,double>(x,x*x+1));
      currentrow.push_back(fastRangeSums(x,N,x*x+1,small));
      currentrow.push_back(fastFirstMoments(x,N,x*x+1,small));
      print(currentrow);
    //  data.push_back(currentrow);
    }
    // building table
   //printTable(data);
    //data.clear();
   // building table
 //  printTable(data);
    
    //table done
    // next, we run our experiments on data sets of the order of the gigabyte ...
     /*
     * First, we do some toy examples...
     */  
    if(false) {
      int big = 500;
      for(int64 size = 1024; size < 100000; size *= 8) {
        fastRangeSums(2,N,size+1,big);
        fastRangeSums(4,N,size+1,big);
        fastRangeSums(8,N,size+1,big);
        slowRangeSums(size+1,big); 
      }
      for(int64 size = 1024; size < 100000; size *= 8) {
        fastFirstMoments(2,N,size+1,big);
        fastFirstMoments(4,N,size+1,big);
        fastFirstMoments(8,N,size+1,big);
        slowFirstMoments(size+1,big); 
      }
    }
    cout << " Done. "<< endl;
 }

  
}





