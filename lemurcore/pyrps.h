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

#ifndef PYRPS_H
#define PYRPS_H

#include "common.h"
#include "datacube.h"

/**
  * The idea of PyRPS is that it allows
  * to add a secondary data cube to an
  * existing data cube so that prefix
  * sums can be done very quickly.
  *
  * And it can be somewhat confusing as it
  * both extends DataCube and also wrap
  * a datacube, but that's because there
  * are really two data cubes, the original
  * one (PyRPS would store a pointer to it)
  * and the PyRPS cube which can be used for
  * fast prefixSums. Because there are really
  * *two* cubes, there must really be two
  * objects (unfortunately).
  *
  * Usage: create a data cube, feed it into
  * a PyRPS object. Next, do not update directly
  * your data cube, but go throuhg PyRPS. Your
  * puts are now slower, but you've got fast
  * prefixSums.
  *
  *@author Daniel Lemire
  */
template <class DataType, class LongDataType>
class PyRPS : public DataCube<LongDataType, LongDataType> {
public:
//	ChunkedDataCube(vector<int> Shape, vector<int> ShapeOfChunks, int NumberOfChunksInMemory = 10);
  PyRPS(DataCube<DataType,LongDataType> * dc, vector<int> Basis);
  PyRPS(const PyRPS<DataType,LongDataType> & DC);
  virtual ~PyRPS();
  virtual void transformBuffer(); // need to be called after constructor
  virtual void fillWithZeroes();
  virtual void fillWith(const DataType value);// a somewhat optimized version
  virtual uint64 prefixSumCost(const vector<int>&)  ;
  virtual LongDataType prefixSum(const vector<int>&)  ;//use the buffer to answer the query fast!
  virtual void put (const DataType, const vector<int>&);// need an update to the buffer
  virtual DataType get(const vector<int>&) ;
  using  DataCube<LongDataType, LongDataType>::get;
  using DataCube<LongDataType, LongDataType>::put;
  using DataCube<LongDataType, LongDataType>::prefixSum;
  using DataCube<LongDataType, LongDataType>::prefixSumCost;
  virtual float rangeQuery(const vector<RangedCubicPolynomial>& polys);
  virtual DataType next(vector<int>& index, 
      const vector<int>& start, const vector<int>& bound) throw(NoMoreValuesException);

protected:
  // this next line should probably never be used publicly
  virtual void copyCubeToBuffer();

  // this is to read directly the transformed buffer, not very useful?
  virtual LongDataType getBuffer(const int i1, ...)  ;
  // this would probably never be used?, please note that it is static!
  static void transform(vector<LongDataType> & data, const int basis);
  //
  DataCube<DataType, LongDataType> * mDataCube;
  const vector<int> mBasis;
  uint64 dimensionProductExceptFor(const uint ExceptionDimension) const;
  uint64 dimensionProductUpTo(const uint UpToDimension) const;
};

/*template <class DataType, class LongDataType>
DataType PyRPS<DataType,LongDataType>::prefixSum(const int i1, ...){
    va_list arguments;
    va_start(arguments,i1);
    parseIndices(i1, arguments);
    return prefixSum();
}*/


template <class DataType, class LongDataType>
LongDataType PyRPS<DataType,LongDataType>::getBuffer (const int i1, ...) {
  if(OptimizationOptions::verbosePutGet) cout << "getting" << endl;
  va_list arguments;
  va_start(arguments,i1);
  return DataCube<DataType,LongDataType>::get( parseIndices(i1, arguments) );
}


template <class DataType, class LongDataType>
PyRPS<DataType,LongDataType>::PyRPS(DataCube<DataType,LongDataType> * dc, vector<int> Basis) : DataCube<LongDataType, LongDataType>(dc->getShape()), mDataCube(dc), mBasis(Basis) {
  for (uint k = 0; k < mBasis.size() ; ++k) _assert(mBasis[k] > 0);
  _assert(mBasis.size() == this->mShape.size());
  if(OptimizationOptions::verbosePyRPS) 
    for(uint k = 0; k < mBasis.size() ; ++k) 
      cout <<"mBasis["<<k<<"] " <<mBasis[k];
}

template <class DataType, class LongDataType>
PyRPS<DataType,LongDataType>::PyRPS(const PyRPS<DataType,LongDataType> & DC): DataCube<LongDataType, LongDataType>(DC), mDataCube(DC.mDataCube), mBasis(DC.mBasis)  {}

template <class DataType, class LongDataType>
PyRPS<DataType,LongDataType>::~PyRPS() {}

template <class DataType, class LongDataType>
void PyRPS<DataType,LongDataType>::copyCubeToBuffer() {
  if(OptimizationOptions::verbosePyRPS) cout <<" copyCubeToBuffer "<<endl;
  // this is very slow and stupid
  vector<int> Indices(this->mShape.size(),0);
  vector<int> start(this->mShape.size(),0);
  DataCube<LongDataType,LongDataType>::put(mDataCube->get(Indices),Indices);
  while(MathUtil::increment( Indices, start, this->mShape)) 
    DataCube<LongDataType,LongDataType>::put(mDataCube->get(Indices),Indices);
  if(OptimizationOptions::verbosePyRPS) cout <<" copyCubeToBuffer ok"<<endl;
}

template <class DataType, class LongDataType>
void PyRPS<DataType,LongDataType>::transformBuffer(){
  if(OptimizationOptions::verbosePyRPS) cout << "transforming buffer"<<endl;
  // first I've got to init the file
  copyCubeToBuffer();
  // this is truly expensive!
  for(uint Dim = 0; Dim < this->mShape.size() ; ++Dim) {
    if(OptimizationOptions::verbosePyRPS) cout << "transforming buffer Dim = " << Dim << endl;
    const uint64 NumberOfTransformsToDo = dimensionProductExceptFor(Dim);
    if(OptimizationOptions::verbosePyRPS) cout << "transforming buffer NumberOfTransformsToDo = " << NumberOfTransformsToDo << endl;
    const uint64 Sampling = dimensionProductUpTo(Dim);
    if(OptimizationOptions::verbosePyRPS) cout << "transforming buffer Sampling = "<< Sampling << endl;
    uint64 StartingPosition = 0;
    int index = 0;
    vector<DataType> Buffer(this->mShape[Dim],0);// temporary buffer
    for(uint64 TransformNumber = 0; TransformNumber < NumberOfTransformsToDo; ++TransformNumber) {
       StartingPosition = (TransformNumber % Sampling) + (TransformNumber / Sampling)* Sampling * this->mShape[Dim];
        // we recover the data and buffer it
        for (index = 0 ; index <this->mShape[Dim]; ++index) {
            const vector<int>& position = parseIndices(StartingPosition + Sampling * index);	
            Buffer[index] = DataCube<LongDataType,LongDataType>::get(position);
        }
        transform(Buffer, mBasis[Dim]);
        //finally, we store it inside the new datacube
        for (index = 0 ; index < this->mShape[Dim]; ++index){ 
          const vector<int>& position = parseIndices(StartingPosition + Sampling * index) ;
          DataCube<LongDataType,LongDataType>::put(Buffer[index],position);
        }
    }
  }
}
template <class DataType, class LongDataType>
void PyRPS<DataType,LongDataType>::transform(vector<LongDataType>& data, const int basis) {
    uint i = 0; int b = 0;
    // the most expensive part is here
    for(i = 0; i  < data.size(); i+= basis) 
      for(b = 1; b < basis; ++b) 
        data[i + b] +=  data[i + b - 1] ;
    // higher level steps are much less expensive, so we can put them in a loop
    for(uint sampling = basis; sampling < data.size(); sampling *= basis) { 
      for(i = 0; i < data.size(); i+= basis*sampling) 
        for(b = 1; b < basis; ++b) 
          data[i + (b+1)*sampling - 1] += data[i + b*sampling - 1] ;
    }

}

template <class DataType, class LongDataType>
uint64 PyRPS<DataType,LongDataType>::dimensionProductExceptFor(const uint ExceptionDimension) const {
    uint64 product = 1;
    _assert(ExceptionDimension < this->mShape.size());
    for(uint dim = 0; dim < this->mShape.size() ; ++dim) {
        if(ExceptionDimension != dim) product *= this->mShape[dim];
    }
    return product;
}

template <class DataType, class LongDataType>
uint64 PyRPS<DataType,LongDataType>::dimensionProductUpTo(const uint UpToDimension) const {
    uint64 product = 1;
    // next line probably too complicated
    _assert(UpToDimension < this->mShape.size());
    for(uint dim = 0; dim < UpToDimension ; ++dim) {
        product *= this->mShape[dim];
    }
    return product;
}

/*
* This is a test method. Delete eventually
*/
template <class DataType, class LongDataType>
uint64 PyRPS<DataType,LongDataType>::prefixSumCost(const vector<int>& Indices)  {
  if(OptimizationOptions::verbosePyRPS) cout << "prefix sum" <<endl;
  vector<deque<int> > BufferOfIndices(this->mShape.size());
  for (uint dim = 0; dim < this->mShape.size(); ++dim) {
    BufferOfIndices[dim].push_front(Indices[dim]);
    if(OptimizationOptions::verbosePyRPS) cout << "starting with "<<Indices[dim] <<" size = "<<BufferOfIndices[dim].size() <<endl;
    for(int power = mBasis[dim]; power < this->mShape[dim] ; power *= mBasis[dim]) {
            int NewIndex = power * ( (Indices[dim] + 1) / power ) - 1;
            if( NewIndex == -1 ) break;
            if( NewIndex != BufferOfIndices[dim].front() ) BufferOfIndices[dim].push_front(NewIndex);
            if(OptimizationOptions::verbosePyRPS) cout << "adding " << NewIndex <<" size = " << BufferOfIndices[dim].size()<<endl;
    }
  }
  if(OptimizationOptions::verbosePyRPS) cout << "prefix sum, buffer of indices done" <<endl;
  vector<uint64> CumulativeNumberOfIndices(this->mShape.size() + 1,0);
  CumulativeNumberOfIndices[0] = 1;
  for(uint dim = 1; dim < this->mShape.size() + 1; ++dim)
    CumulativeNumberOfIndices[dim] = CumulativeNumberOfIndices[dim - 1] * BufferOfIndices[dim - 1].size();
  const uint64 NumberToVisit = CumulativeNumberOfIndices[this->mShape.size()];
  return NumberToVisit;
}

template <class DataType, class LongDataType>
LongDataType PyRPS<DataType,LongDataType>::prefixSum(const vector<int>& Indices) {
  if(OptimizationOptions::verbosePyRPS) cout << "prefix sum" <<endl;
  vector<deque<int> > BufferOfIndices(this->mShape.size());
  for (uint dim = 0; dim < this->mShape.size(); ++dim) {
    BufferOfIndices[dim].push_front(Indices[dim]);
    if(OptimizationOptions::verbosePyRPS) cout << "starting with "<<Indices[dim] <<" size = "<<BufferOfIndices[dim].size() <<endl;
    for(int power = mBasis[dim]; power < this->mShape[dim] ; power *= mBasis[dim]) {
            int NewIndex = power * ( (Indices[dim] + 1) / power ) - 1;
            if( NewIndex == -1 ) break;
            if( NewIndex != BufferOfIndices[dim].front() ) BufferOfIndices[dim].push_front(NewIndex);
            if(OptimizationOptions::verbosePyRPS) cout << "adding " << NewIndex <<" size = " << BufferOfIndices[dim].size()<<endl;
    }
  }
  if(OptimizationOptions::verbosePyRPS) cout << "prefix sum, buffer of indices done" <<endl;
  vector<uint64> CumulativeNumberOfIndices(this->mShape.size() + 1,0);
  CumulativeNumberOfIndices[0] = 1;
  for(uint dim = 1; dim < this->mShape.size() + 1; ++dim)
    CumulativeNumberOfIndices[dim] = CumulativeNumberOfIndices[dim - 1] * BufferOfIndices[dim - 1].size();
  const uint64 NumberToVisit = CumulativeNumberOfIndices[this->mShape.size()];
  LongDataType sum (0);
  if(OptimizationOptions::verbosePyRPS) cout << "prefix sum, cumul indices done, must visit "<< NumberToVisit <<endl;
  //
  vector<int> MovingIndices(Indices);
  for(uint64 index = 0; index < NumberToVisit; ++index) {
    for(uint dim = 0;  dim < this->mShape.size(); ++dim) {
        MovingIndices[dim] = BufferOfIndices[dim][(int)( index / CumulativeNumberOfIndices[dim]) % BufferOfIndices[dim].size()];
    }
    sum += DataCube<LongDataType, LongDataType>::get(MovingIndices);
    if(OptimizationOptions::verbosePyRPS) cout << " sum = " << sum << endl;
  }
  return sum;
}


template <class DataType, class LongDataType>
void PyRPS<DataType,LongDataType>::put(const DataType value, const vector<int>& Indices) {
  const DataType difference = value - mDataCube->get(Indices);
  if(difference == DataType(0)) return; // no need to do anything!
  // next we update the data cube
  mDataCube->put(value,Indices);
  // must update tree
  vector<deque<int> > BufferOfIndices(this->mShape.size());
  for(uint dim = 0; dim < this->mShape.size(); ++dim) {
    int end = Indices[dim];
    for (int power = mBasis[dim]; power < this->mShape[dim]; power *= mBasis[dim]) {
      int start = power * (Indices[dim]/power) + power - 1;
      const int nextpower = power * mBasis[dim];
      end = nextpower * (Indices[dim]/nextpower) + nextpower - 1;
      for (int index = start ; index < end ; index += power)
        BufferOfIndices[dim].push_front(index);
    }
    BufferOfIndices[dim].push_front(end);
  }
  vector<uint64> CumulativeNumberOfIndices(this->mShape.size() + 1,0);
  CumulativeNumberOfIndices[0] = 1 ;
  for(uint dim = 1; dim < this->mShape.size() + 1; ++dim)
    CumulativeNumberOfIndices[dim] = CumulativeNumberOfIndices[dim  - 1] * BufferOfIndices[dim - 1].size();
  const uint64 NumberToVisit = CumulativeNumberOfIndices[this->mShape.size()];
  vector<int> MovingIndices(Indices);
  for(uint64 index = 0; index < NumberToVisit; ++index) {
    for(uint dim = 0;  dim < this->mShape.size(); ++dim)
        MovingIndices[dim] = BufferOfIndices[dim][(int)( index / CumulativeNumberOfIndices[dim]) % BufferOfIndices[dim].size()];
        DataCube<LongDataType,LongDataType>::put(
            DataCube<LongDataType,LongDataType>::get(MovingIndices) + difference,
            MovingIndices
        );
  }
}

template <class DataType, class LongDataType>
DataType PyRPS<DataType,LongDataType>::get (const vector<int>& Indices)  {
   return mDataCube->get(Indices);  
}

template <class DataType, class LongDataType>
void PyRPS<DataType,LongDataType>::fillWithZeroes() {
  // clean data cube
  mDataCube->fillWithZeroes();
  // clean buffer
  DataCube<LongDataType,LongDataType>::fillWithZeroes();
}

template <class DataType, class LongDataType>
void PyRPS<DataType,LongDataType>::fillWith(const DataType value) {
  mDataCube->fillWith(value);
  transformBuffer();
}

template <class DataType, class LongDataType>
float PyRPS<DataType,LongDataType>::rangeQuery(const vector<RangedCubicPolynomial>& polys) {
  cout << "[WARNING] PyRPS should check whether this is a range sum and if so, do it fast. TODO!!!" <<endl;
  return mDataCube->rangeQuery(polys);
}

template <class DataType, class LongDataType> 
DataType PyRPS<DataType,LongDataType>::next(vector<int>& index, 
    const vector<int>& start,  const vector<int>& bound) 
  throw(NoMoreValuesException) {
  return mDataCube->next(index,start,bound);
}
#endif


  
