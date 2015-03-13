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
#ifndef RAMCUBE_H
#define RAMCUBE_H

#include "common.h"
#include "datacube.h"
#include "fileutil.h"

/**
  * This builds a cube entirely in virtual memory, as a big md-array (unchunked)
  * To resemble Richard Liu's thesis, the array should be chunked and
  * and memory mapped to a backing file
  *
  * Currently, no files are updated....
  *
  * I am just getting tired of waiting for seeks on files that aren't that big... -ofk
  *
  */

template <class DataType, class LongDataType>
class RAMCube : public DataCube<DataType,LongDataType> {
public: 
  RAMCube(const vector<int>& Shape);
  RAMCube(const RAMCube<DataType,LongDataType>& DC);
  virtual ~RAMCube();
  virtual void open(const char *FileName = null);
  virtual void close();
  virtual void fillWithZeroes();
  virtual DataType get(const vector<int>&) ;
  virtual void put (DataType,const vector<int>&);
  using  DataCube<DataType, LongDataType>::get;
  using DataCube<DataType, LongDataType>::put;
  virtual DataType next(vector<int>& index, 
      const vector<int>& start,  const vector<int>& bound) throw(NoMoreValuesException);

  
private:
  DataType *bigMemBuf;
  bool mUsingTemp;// if there is no file specified, then we never use the disk!
};

template <class DataType, class LongDataType>
RAMCube<DataType,LongDataType>::RAMCube(const vector<int>& Shape) : DataCube<DataType,LongDataType>(Shape) , bigMemBuf(), mUsingTemp(0){
  if(OptimizationOptions::verbose) cout << "RAMCube::RAMCube " << endl;
  // Couldn't understand what you were after so I will assume this is a mistake.
  // We most certainly do not want the volume of a cube to be a floating point
  // number as the volume is the number of cells. The number of cells in an 
  // integer. Now, if we want to get fancy about it, we just need to introduce
  // a new method later. No? DL

  // Not a mistake!  
  // The difficulty is I want to use getVolume (today!) on some outrageously BIG cubes whose
  // volumes exceed 64 bits.  I promise to abort on these cube :)
  // To know whether to abort, I need to calculate volume *w/o risk of overflow*.
  // Approximate volume would be okay for that. 
  //    So I added a new method to datacube to give approx storage as a double.
  //    I'm still pretty sure that we could've used doubles throughout, as long as
  //    the true volume requires fewer bits than the mantissa of a double
  //    contains.  This is cryptic but would reduce the datacube API by one method. -OFK
  //
  //  Yes, I buy your safety argument... alas, this could be confusing for someone looking
  //  at it for the first time... this is especially confusing since we can using uint64. Now,
  //  if you beat uint64 in storage requirement... hmmm... let's see... well... that's pretty
  //  much the storage of all hard drives over the planet or something crazy like that... should
  //  be good enough for another century or so...
  //  
  // bigMemBuf = new DataType[ getVolume() ];   moved to open.
}

template <class DataType, class LongDataType>
RAMCube<DataType,LongDataType>::RAMCube(const RAMCube<DataType,LongDataType>& DC) : 
DataCube<DataType,LongDataType>(DC) , bigMemBuf(DC.bigMemBuf), mUsingTemp(DC.mUsingTemp) {
  if(OptimizationOptions::verbose) cout << "RAMCube::RAMCube " << endl;
}


template <class DataType, class LongDataType>
RAMCube<DataType,LongDataType>::~RAMCube(){
  if(OptimizationOptions::verbose) cout << "RAMCube::~RAMCube " << endl;
  //delete [] bigMemBuf; simpler if you let close handle it
}

template <class DataType, class LongDataType>
void RAMCube<DataType,LongDataType>::open(const char *FileName) {
  if(OptimizationOptions::verbose) cout << "RAMCube::open " << endl;
  bigMemBuf = new DataType[ getVolume() ];
  // do the memory allocation now

  if(FileName == null) {
    this->mFileName = (const char *) tempnam(null,"molasse");
    mUsingTemp = true;
 }	else {
    this->mFileName = FileName;
    mUsingTemp = false;// need to save to disk
  }
  if(mUsingTemp) return; // we are done!
  if(!FileUtil::fileExists(this->mFileName) || (FileUtil::getFileSize(this->mFileName) < getVolume()*sizeof(DataType)))
  {// if the file doesn't exist or is smaller than needed, then fill it will zeroes
    this->mFileStream.open( this->mFileName,std::ios::binary | std::ios::in | std::ios::out| std::ios::trunc);
    _assert(this->mFileStream.is_open());
      _assert(!this->mFileStream.fail());

    DataCube<DataType, LongDataType>::fillWithZeroes();// thanks! sean!

    this->mFileStream.close();
  } else {// I think we only want to read the content of the file in memory when the file already existed and is proper size.
      this->mFileStream.open( this->mFileName,std::ios::binary | std::ios::in | std::ios::out);
      _assert(this->mFileStream.is_open());
      _assert(!this->mFileStream.fail());

      this->mFileStream.seekg(0);
      _assert(! this->mFileStream.fail());

      this->mFileStream.read((char *)bigMemBuf,sizeof(DataType)*getVolume());
      _assert(! this->mFileStream.fail());
      this->mFileStream.close();
  }
}

template <class DataType, class LongDataType>
void RAMCube<DataType,LongDataType>::close() {
  if(OptimizationOptions::verbose) cout << "DataCube::close " << endl;
//	delete [] bigMemBuf;// totally unsafe!!! Careful!
  if(! mUsingTemp) {
    try {
      this->mFileStream.open( this->mFileName,std::ios::binary | std::ios::in | std::ios::out);
      _assert(this->mFileStream.is_open());
      _assert(!this->mFileStream.fail());

      this->mFileStream.seekg(0);
      _assert(! this->mFileStream.fail());

      this->mFileStream.write((char *)bigMemBuf,sizeof(DataType)*getVolume());
      _assert(! this->mFileStream.fail());

      this->mFileStream.close();
    } catch (...) {// no matter what is thrown... we can't help it now!
          cerr<< "Couldn't save to file in ramcube!" << endl;
    }
  }
  delete[] bigMemBuf;
}

template <class DataType, class LongDataType>
DataType RAMCube<DataType,LongDataType>::get (const vector<int>& idx) {
  uint64 offset = computeOffset(idx);
  return bigMemBuf[offset];
}

template <class DataType, class LongDataType>
void RAMCube<DataType,LongDataType>::put(DataType value,const vector<int>& idx) {
  uint64 offset = computeOffset(idx);
  bigMemBuf[offset] = value;
}

template <class DataType, class LongDataType>
void RAMCube<DataType,LongDataType>::fillWithZeroes() {
  if(OptimizationOptions::verbose) cout << "RAMCube::fillWithZeroes " << endl;
  memset(bigMemBuf,0, sizeof(DataType) * getVolume() );  
  if(OptimizationOptions::verbose) cout << "DataCube::fillWithZeroes ok" << endl;
}

template <class DataType, class LongDataType>
DataType RAMCube<DataType,LongDataType>::next(vector<int>& index, 
    const vector<int>& start, const vector<int>& bound) throw(NoMoreValuesException) {
  if(! MathUtil::increment(index,start,bound) ) throw NoMoreValuesException(); 
  uint64 offset = computeOffset(index); 
  uint64 volume = getVolume();
  for(uint64 i = offset; i < volume ; ++i) {// shouldn't be getVolume, should rely on bound
    if(bigMemBuf[i] != 0) {
      uint64 indexoffset = i - offset;
      if( ! MathUtil::add(indexoffset, index, start, bound) ) throw NoMoreValuesException();  
      return bigMemBuf[i];	
    }		
  }
  throw NoMoreValuesException();  
}

#endif


