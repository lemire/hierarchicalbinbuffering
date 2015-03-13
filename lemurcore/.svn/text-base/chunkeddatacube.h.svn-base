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
#ifndef CHUNKEDDATACUBE_H
#define CHUNKEDDATACUBE_H

#include "common.h"
#include "datacube.h"

/**
  * Untested code.
  *
  * This is an implementation of the data cube
  * with a small chunk in memory, rest on disk model.
  * This uses a smarter (?) file format than the
  * stupid linear structure used for the DataCube
  * class and it should be a lot closer to what
  * people do in the industry. For small data cubes,
  * expect this class to underperform however.
  *
  *@author Daniel Lemire
  */


template <class DataType, class LongDataType>
class ChunkedDataCube : public DataCube<DataType, LongDataType> {
public:
  static vector<int> sqrtChunkShapes( vector<int> Shape) {
    vector<int> answer(Shape.size());
    for (uint i=0; i < Shape.size(); ++i)
      answer[i] = int(ceil(sqrt(double(Shape[i]))));  
    return answer;
  }
  ChunkedDataCube(const vector<int>& Shape, const vector<int>& ShapeOfChunks);
  ChunkedDataCube(const ChunkedDataCube<DataType,LongDataType> & DC);
  virtual ~ChunkedDataCube(){};
  virtual DataType next(vector<int>& index,
      const vector<int>& start,  const vector<int>& bound) throw(NoMoreValuesException);
  virtual uint64 getChunkIndex(const vector<int>&) const;
  virtual int getIndexInChunk(const vector<int>&) const;

protected:
  enum {verbose = false};
  virtual uint64 computeOffset (const vector<int>&) const ;
  virtual uint64 sizeOfSlices (const vector<int>& indices, const int i) const ;
  const vector<int> mShapeOfChunks;

private:
  // nothing

};


template <class DataType, class LongDataType>
uint64 ChunkedDataCube<DataType,LongDataType>::computeOffset (const vector<int>& idx) const {
  const uint64 offset = getChunkIndex(idx)  + getIndexInChunk(idx);
  return offset;
}

template <class DataType, class LongDataType>
ChunkedDataCube<DataType,LongDataType>::ChunkedDataCube(const vector<int>& Shape,
    const vector<int>& ShapeOfChunks) :
    DataCube<DataType,LongDataType> (Shape),
    mShapeOfChunks(ShapeOfChunks) {
  _assert(mShapeOfChunks.size() == this->mShape.size());
}
template <class DataType, class LongDataType>
ChunkedDataCube<DataType,LongDataType>::ChunkedDataCube(const ChunkedDataCube<DataType,LongDataType> & DC):
DataCube<DataType,LongDataType>(DC), mShapeOfChunks(DC.mShapeOfChunks)
{}

/**
 *  This is not nearly as easy as I thought it would.
 *  This is related to something called Morton's indexing
 *  which is related to Peano curves and Hilbert SFC (Space Filling Curves).
 *
 *  See for example Tensor Product Formulation for
 *  Hilbert Space-Filling Curves by Lin et al.
 *
 * The nasty thing about what we do is that we still favor one dimension
 * over the others, just less so. What would really be needed is Hilbert
 * Space Filling Indexing. If we could implement it efficiently for
 * data cubes, it would be progress.
 *
 * I finally got it right for matrices I think, I'm not sure about
 * higher dimensional cases.
 *
 */
template <class DataType, class LongDataType>
uint64 ChunkedDataCube<DataType,LongDataType>::getChunkIndex(const vector<int>& indices) const {
  uint64 offset = 0;
  vector<int> chunk_indexes(indices);
  for(uint i = 0; i < chunk_indexes.size(); ++i)
    chunk_indexes[i]/= mShapeOfChunks[i];
    for(uint i = 0; i < this->mShape.size(); ++i ) {
    offset += chunk_indexes[i] * sizeOfSlices(chunk_indexes,i);
  }
    return offset;
}

/**
 * This method is used by the getChunkIndex method.
 * (Note, this is actually buffered now as mSizeOfSlices for speed.
 *
 * Return the size of a i dimensional slice in dimension 0,...,i-1 in chunk space.
 *
 * For i = 0, returns the size (in cells) of one block
 * For i = 1, returns the size (in cels) of "column" of blocks
 * And so on... i is a dimension here.
 *
 * The indices parameter should contain an index to the location of
 * the current chunk in chunk space. This means we can't buffer this computation very
 * efficiently. I tried hard, but alas, I failed. I think we could
 * buffer this a bit... but would it be worth it? Not sure.
 *
 */
template <class DataType, class LongDataType>
uint64 ChunkedDataCube<DataType,LongDataType>::sizeOfSlices(const vector<int>& chunk_indexes, const int i) const {
    uint64 SizeOfSlices = mShapeOfChunks[i];
    for(uint dim = (uint) i+1; dim < this->mShape.size(); ++dim) {
      SizeOfSlices *= min(mShapeOfChunks[dim],
        this->mShape[dim]-
        chunk_indexes[dim] *
        this->mShapeOfChunks[dim]);
    }
    for(int dim = 0; dim < i; ++dim) {
      SizeOfSlices  *=  this->mShape[dim] ;
    }
    return SizeOfSlices ;
}

template <class DataType, class LongDataType>
int ChunkedDataCube<DataType,LongDataType>::getIndexInChunk(const vector<int>& indices) const {
  int index = 0, base = 1;
  for(uint i = 0; i < this->mShape.size(); ++i ) {
    index += ( indices[i] % mShapeOfChunks[i] ) * base;
    const int boundary_effect =
      indices[i] / mShapeOfChunks[i] * mShapeOfChunks[i] +
      mShapeOfChunks[i] - this->mShape[i];
    if(boundary_effect > 0)
      base *=  mShapeOfChunks[i] - boundary_effect;
    else base *=  mShapeOfChunks[i];
  }
  return index;
}

template <class DataType, class LongDataType>
DataType ChunkedDataCube<DataType,LongDataType>::next(vector<int>& index,
    const vector<int>& start,  const vector<int>& bound)
  throw(NoMoreValuesException) {
    cout <<"Warning!!! This is not implemented at all." <<endl;
    _assert(false);
    return -9999999;
}

#endif
