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

#ifndef DATACUBE_H
#define DATACUBE_H

#include "common.h"
#include "cubicpolynomial.h"

/**
  * REFERENCE IMPLEMENTATION OF HOW A DATACUBE MUST BEHAVE.
  *
  * This is meant to be the base class and serves as
  * the external API. Notice here that a data cube
  * is taken to be a very simple thing: a big
  * multidimensional array. A more realistic
  * interface should be developed later.
  *
  * None of this is meant to be thread safe or
  * any such thing.
  *
  *@author Daniel Lemire
  */

class NoMoreValuesException {public: NoMoreValuesException(){}}; 

template <class DataType, class LongDataType>
class DataCube {
public: 
  DataCube(vector<int> Shape);
  DataCube(const DataCube<DataType,LongDataType>& DC);
  DataCube<DataType,LongDataType>& operator=(const DataCube<DataType,LongDataType>& DC);
  virtual ~DataCube();
  virtual void open(const char *FileName = null);
  virtual void close();
  virtual DataType get(const int i1, ...) ;
  virtual void put (const DataType, const int i1, ...);
  virtual DataType get(const uint64 indices)  ;
  virtual DataType get(const vector<int>&) ;
  virtual void put (const DataType value, uint64 indices);
  virtual void put (const DataType value, const vector<int>& );
  virtual void fillWithZeroes();
  virtual void fillWith(const DataType value);
  virtual vector<int> getShape() const;
  virtual vector<int> maxIndex() const;
  virtual uint64 getVolume() const;
  virtual double getApproxDenseStorageSize() const;
  // non trivial queries
  // Note that all the prefixSum queries are now obselete. Better
  // call polynomialSum which is the most general case.
  virtual LongDataType prefixSum(const int i1, ...); /* to be removed */
  virtual uint64 prefixSumCost(const int i1, ...); /*to be removed */
  virtual LongDataType prefixSum(const vector<int>&); /* to be removed*/
  virtual uint64 prefixSumCost(const vector<int>&); /* to be removed */
  // this should cover a wide range of range sums, all but max, min, percentiles.
  // you can do all the count, average, standard deviation and stuff using
  // polynomial range sums
  virtual float rangeQuery(const vector<RangedCubicPolynomial>& polys);
  // the next function returns the next non-zero values within index and bounds.
  // if you keep calling this function, you are garanteed to visit all non-zero
  // values. This should be faster than get.
  virtual DataType next(vector<int>& index, 
      const vector<int>& start, const vector<int>& bound) throw(NoMoreValuesException);
  enum{
    NEXT_BUFFER=16//size of the IO buffer used to make next faster
  };	
protected:
  virtual uint64 computeOffset (const vector<int>&) const ;
  inline vector<int> parseIndices (const int i1, va_list& arguments) const;
  inline vector<int> parseIndices (const uint64 indices) const;
  vector<int> mShape;
  vector<uint64> mOffsetBase;//use in computeOffset for fast computations
  fstream mFileStream;// used to be a pointer, should be faster on stack
  const char * mFileName;
private:


};

  
/*
* Implementation follows
*/

template <class DataType, class LongDataType>
DataCube<DataType,LongDataType>::DataCube(const DataCube<DataType,LongDataType>& DC) :
  mShape(DC.mShape), mOffsetBase(DC.mOffsetBase), mFileStream(), mFileName(DC.mFileName) { 
  if(mFileName != NULL) open(mFileName);// that might be unsafe to open a file twice...?
}
template <class DataType, class LongDataType>
DataCube<DataType,LongDataType>& DataCube<DataType,LongDataType>::operator=(
    const DataCube<DataType,LongDataType>& DC) {
  mShape = DC.mShape;
  mOffsetBase = DC.mOffsetBase;
  if(DC.mFileName != NULL) {
    close();
    open(DC.mFileName);
  } else close();
  return *this;
}

template <class DataType, class LongDataType>
LongDataType DataCube<DataType,LongDataType>::prefixSum(const int i1, ...){
va_list arguments;
va_start(arguments,i1);
return prefixSum( parseIndices(i1, arguments) );
}

template <class DataType, class LongDataType>
uint64 DataCube<DataType,LongDataType>::prefixSumCost(const int i1, ...){
va_list arguments;
va_start(arguments,i1);
return prefixSumCost( parseIndices(i1, arguments) );
}



template <class DataType, class LongDataType>
uint64 DataCube<DataType,LongDataType>::prefixSumCost(const vector<int>& idx){
vector<int> Bounds (idx);
for(uint i = 0; i < Bounds.size() ; ++i) Bounds[i]++;
vector<int> Indices(Bounds.size(),0);
vector<int> start(Bounds.size(),0);
uint64 number = 0;
while(MathUtil::increment( Indices, start, Bounds)) {
  ++number;
}
return number;
}

template <class DataType, class LongDataType>
LongDataType DataCube<DataType,LongDataType>::prefixSum(const vector<int>& idx){
// what was I smoking when I wrote this code? Oh well... it works!
vector<int> Bounds = idx;
for(uint i = 0; i < Bounds.size() ; ++i) Bounds[i]++;
vector<int> Indices(Bounds.size(),0);
vector<int> start(Bounds.size(),0);
LongDataType sum (get(Indices));
while(MathUtil::increment( Indices, start, Bounds)) {
  sum += get(Indices);
}
return sum;
}

template <class DataType, class LongDataType>
float DataCube<DataType,LongDataType>::rangeQuery(const vector<RangedCubicPolynomial>& polys) {
vector<int> start(polys.size(),0);
vector<int> Bounds(polys.size(),0);
for(uint dim = 0 ; dim < polys.size() ; ++ dim) {
  start[dim] = polys[dim].mStart;
  Bounds[dim] = polys[dim].mEnd;
  if(start[dim] == Bounds[dim]) return 0.0f; // nothing to do
}
vector<int> Indices = start;
float sum = 0.0f;
do {
  float value = 1.0f;
  for (uint dim = 0; dim < polys.size(); ++ dim) value *= polys[dim](Indices[dim]);
  sum += get(Indices) * value;
} while(MathUtil::increment( Indices, start, Bounds)); 
return sum;
}


template <class DataType, class LongDataType>
vector<int> DataCube<DataType,LongDataType>::getShape() const {
return mShape;
}

template <class DataType, class LongDataType>   // seems widely useful, albeit trivial
vector<int> DataCube<DataType,LongDataType>::maxIndex() const {
  vector<int> answer = getShape();
  for (uint i=0; i < answer.size(); ++i) 
    --answer[i];
  return answer;
}




template <class DataType, class LongDataType>
uint64 DataCube<DataType,LongDataType>::getVolume() const {
  uint64 prod = 1;
  for(uint i = 0; i < this->mShape.size(); ++i )
    prod *= this->mShape[i];
  return prod;
}


template <class DataType, class LongDataType>
double DataCube<DataType,LongDataType>::getApproxDenseStorageSize() const {
  double prod = 1.0;
  for(uint i = 0; i < mShape.size(); ++i )
    prod *= mShape[i];
  return  prod * sizeof( DataType);
}

template <class DataType, class LongDataType>
DataCube<DataType,LongDataType>::DataCube(vector<int> Shape) : mShape(Shape), mOffsetBase(mShape.size(),1) ,
mFileStream(), mFileName(NULL) {
if(OptimizationOptions::verbose) cout << "DataCube::DataCube " << endl;
for(uint i = 1; i < mOffsetBase.size(); ++i ) mOffsetBase[i] = mOffsetBase[i-1] * mShape[i - 1];
}

template <class DataType, class LongDataType>
DataCube<DataType,LongDataType>::~DataCube(){
if(OptimizationOptions::verbose) cout << "DataCube::~DataCube " << endl;
//close(); was a bad idea, better let the user close it, this simplifies copy constructor business
}

template <class DataType, class LongDataType>
void DataCube<DataType,LongDataType>::open(const char *FileName) {
  if(OptimizationOptions::verbose) cout << "DataCube::open " << endl;
  close();
  if(FileName == null) mFileName = (const char *) tempnam(null,"molasse"); else mFileName = FileName;
  if(OptimizationOptions::verboseIO) cout << "opening file "<< mFileName << endl;
  mFileStream.open( mFileName,std::ios::binary | std::ios::in | std::ios::out ) ;//| std::ios::trun
  if(! mFileStream.good()) {  // ofk: is_open() didn't cut it...
    if (mFileStream.is_open()) mFileStream.close();
    mFileStream.clear();
    mFileStream.open( mFileName,std::ios::binary | std::ios::in | std::ios::out| std::ios::trunc);
  }
  _assert(mFileStream.is_open());
  _assert(!mFileStream.fail());
  if(FileUtil::getFileSize(mFileName) < getVolume()*sizeof(DataType)) fillWithZeroes(); //automagical!
  //mFileStream.rdbuf()->pubsetbuf(new char[1024*1024],1024*1024);
}

template <class DataType, class LongDataType>
void DataCube<DataType,LongDataType>::close() {
  if(OptimizationOptions::verbose) cout << "DataCube::close " << endl;
  mFileStream.close();
  mFileName = NULL;
}


template <class DataType, class LongDataType>
uint64 DataCube<DataType,LongDataType>::computeOffset (const vector<int>& idx) const {
  uint64 offset = 0;
  for(uint i = 0; i < mShape.size(); ++i ) offset += idx[i] * mOffsetBase[i];
  return offset;
}


template <class DataType, class LongDataType>
vector<int> DataCube<DataType,LongDataType>::parseIndices (const int i1, va_list& arguments) const {
  // this type of highly modular code might not get me a job with Oracle
  vector<int> indices(mShape.size(),0);
  indices[0] = i1;
  for(uint i = 1; i < mShape.size(); ++i) indices[i] = va_arg(arguments, int);
  return indices;
}

template <class DataType, class LongDataType>
vector<int> DataCube<DataType,LongDataType>::parseIndices (const uint64 indices) const {
  vector<int> idx(mShape.size(),0);
  idx[0] = indices % mShape[0];
  uint64 base = mShape[0];  
  for(uint i = 1; i < mShape.size(); ++i) {
    idx[i] = ( indices / base ) % mShape[i];
    base *= mShape[i];
  }
  return idx;
}

template <class DataType, class LongDataType>
DataType DataCube<DataType,LongDataType>::get (const int i1, ...) {
  if(OptimizationOptions::verbosePutGet) cout << "getting" << endl;
  va_list arguments;
  va_start(arguments,i1);
  return get(parseIndices(i1, arguments));
}


template <class DataType, class LongDataType>
DataType DataCube<DataType,LongDataType>::get (uint64 indices)  {
  return get(parseIndices(indices));
}


template <class DataType, class LongDataType>
DataType DataCube<DataType,LongDataType>::get(const vector<int>& idx)  {
  uint64 offset = computeOffset(idx);
  _assert(mFileStream != null);
  if(OptimizationOptions::verboseIO) cout << "going to position "<< offset << endl;
  mFileStream.seekg( offset * sizeof(DataType) );
  _assert(! mFileStream.fail());
  DataType answer;
  mFileStream.read((char *)& answer,sizeof(DataType));
  _assert(! mFileStream.fail());
  return answer;
}
template <class DataType, class LongDataType>
void DataCube<DataType,LongDataType>::put (DataType value, uint64 indices)  {
  put(value, parseIndices(indices));
}

template <class DataType, class LongDataType>
void DataCube<DataType,LongDataType>::put(DataType value, const int i1,...) {
  if(OptimizationOptions::verbosePutGet) cout << "putting" << endl;
  va_list arguments;
  va_start(arguments,i1);
  put(value, parseIndices(i1,arguments));
}

template <class DataType, class LongDataType>
void DataCube<DataType,LongDataType>::put(DataType value,const vector<int>& idx) {
  uint64 offset = computeOffset(idx);
  _assert(mFileStream != null);
  if(OptimizationOptions::verboseIO) cout << "going to position "<< offset << endl;
  mFileStream.seekp( offset * sizeof(DataType) );
  _assert(! mFileStream.fail());
  mFileStream.write((char *) & value, sizeof(DataType));
  _assert(! mFileStream.fail());
}

template <class DataType, class LongDataType>
void DataCube<DataType,LongDataType>::fillWithZeroes() {
  if(OptimizationOptions::verbose) cout << "DataCube::fillWithZeroes " << endl;
  DataType * BigArray = new DataType[mShape[0]];
  memset(BigArray,0,mShape[0] * sizeof(DataType));
  int NumberOfTimes = 1;
  for(uint i = 1 ; i < mShape.size() ; ++i ) NumberOfTimes *= mShape[i];
  _assert(mFileStream != null);
  mFileStream.seekp(0);
  for(int times = 0; times < NumberOfTimes; ++times) mFileStream.write((char *) BigArray,mShape[0] * sizeof(DataType));
  _assert(! mFileStream.fail());
  delete[] BigArray;
  if(OptimizationOptions::verbose) cout << "DataCube::fillWithZeroes ok" << endl;
}


template <class DataType, class LongDataType>
void DataCube<DataType,LongDataType>::fillWith(const DataType value) {
  vector<int> indices(mShape.size(),0);
  put(value, indices);
  vector<int> start(mShape.size(),0);
  while(MathUtil::increment( indices,start, mShape)) {
    put(value,indices);
  }
}

template <class DataType, class LongDataType>
DataType DataCube<DataType,LongDataType>::next(vector<int>& index, 
    const vector<int> & start, const vector<int>& bound) throw(NoMoreValuesException) {
  if(! MathUtil::increment(index,start,bound) ) throw NoMoreValuesException();
  uint64 offset = computeOffset(index);
  _assert(mFileStream != null);
  mFileStream.seekg( offset * sizeof(DataType) );
  _assert(! mFileStream.fail());
  DataType answer[NEXT_BUFFER];// read 16 values at a time instead
  mFileStream.read((char *)& answer,sizeof(DataType)*NEXT_BUFFER);
  uint NumberRead = mFileStream.gcount() / sizeof(DataType);
  uint64 ZeroesFound = 0;
  while(NumberRead > 0) {// this loop is suboptimal in the sense that we can keep reading way past "bound".
    for(uint i = 0; i < NumberRead; ++i) {
      if(answer[i] !=0) {
        ZeroesFound += i;
        if(NumberRead < NEXT_BUFFER) mFileStream.clear();// return the stream to a good state 
        if( ! MathUtil::add(ZeroesFound,index,start,bound) )
          throw NoMoreValuesException();// found a non-zero value past "bound": we read too far
        return answer[i];
      }
    }
    ZeroesFound += NumberRead;
    mFileStream.read((char *)& answer,sizeof(DataType)*NEXT_BUFFER);
    NumberRead = mFileStream.gcount() / sizeof(DataType); 
  }
  mFileStream.clear();
  throw NoMoreValuesException();
}
#endif


