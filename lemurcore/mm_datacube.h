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
#ifndef MM_DATACUBE_H
#define MM_DATACUBE_H

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/mman.h>
#include "common.h"
#include "datacube.h"
#include "fileutil.h"


template <class DataType, class LongDataType>
class mmDataCube : public DataCube<DataType, LongDataType>
{
public:
    mmDataCube(const vector<int>& Shape);
  mmDataCube(const mmDataCube<DataType,LongDataType> & DC);
  virtual ~mmDataCube();

  virtual void put (DataType, const vector<int>&);
  virtual DataType get(const vector<int>&) ;
  using  DataCube<DataType, LongDataType>::get;
  using DataCube<DataType, LongDataType>::put;
  virtual void fillWithZeroes();
  virtual void open(const char *FileName = null);
  virtual void close();
 
private:
  int fd; //file descriptor
  DataType *data; //file data

};



template <class DataType, class LongDataType>
mmDataCube<DataType,LongDataType>::mmDataCube(const vector<int>& Shape) : DataCube<DataType,LongDataType> (Shape), fd(-1), data(NULL)
{
    if(OptimizationOptions::verbose)
    cout << "mmDataCube::mmDataCube " << endl;

}


template <class DataType, class LongDataType>
mmDataCube<DataType,LongDataType>::mmDataCube(const mmDataCube<DataType,LongDataType> & DC): DataCube<DataType,LongDataType>(DC), fd(-1), data(NULL)
{}


template <class DataType, class LongDataType>
mmDataCube<DataType,LongDataType>::~mmDataCube()
{
  if(OptimizationOptions::verbose)
    cout << "mmDataCube::~mmDataCube " << endl;

}


template <class DataType, class LongDataType>
void mmDataCube<DataType,LongDataType>::open(const char *FileName)
{
  if(OptimizationOptions::verbose)
    cout << "DataCube::open " << endl;

  close();

  if(FileName == null)
    mFileName = (const char *) tempnam(null,"molasse");
  else
    mFileName = FileName;

    if(OptimizationOptions::verboseIO)
    cout << "opening file "<< mFileName << endl;

  //if the file doesn't exist, or if it exists and the file size is zero....make it the right size
  if(!FileUtil::fileExists(mFileName) || (FileUtil::getFileSize(mFileName) < getVolume()*sizeof(DataType)))
  {
    mFileStream.open( mFileName,std::ios::binary | std::ios::in | std::ios::out| std::ios::trunc);
    _assert(mFileStream.is_open());
      _assert(!mFileStream.fail());

    DataCube<DataType, LongDataType>::fillWithZeroes();

    mFileStream.close();
  }

  fd = ::open(mFileName, O_RDWR | O_CREAT);
  _assert(fd != -1);
  
  data = (DataType *)mmap((caddr_t)0, sizeof(DataType)*getVolume(), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  _assert(data != MAP_FAILED);

}

template <class DataType, class LongDataType>
void mmDataCube<DataType,LongDataType>::close()
{
    if(OptimizationOptions::verbose)
    cout << "DataCube::close " << endl;
  
  if (data != NULL)
  {
    munmap(data, sizeof(DataType)*getVolume());
    ::close(fd);
    data = NULL;
  }

  mFileName = NULL;
}



template <class DataType, class LongDataType>
DataType mmDataCube<DataType,LongDataType>::get(const vector<int>& idx)
{
    uint64 offset = computeOffset(idx);

    if(OptimizationOptions::verboseIO)
    cout << "going to position "<< data[offset] << endl;

  DataType answer;
  answer = data[offset];
  return answer;
}

template <class DataType, class LongDataType>
void mmDataCube<DataType,LongDataType>::put(DataType value,const vector<int>& idx)
{
    uint64 offset = computeOffset(idx);

  if(OptimizationOptions::verboseIO)
    cout << "offset " << offset*sizeof(DataType) << endl;

  *(data + offset) = value;

  msync((data + offset), sizeof(DataType), MS_SYNC);
}

template <class DataType, class LongDataType>
void mmDataCube<DataType,LongDataType>::fillWithZeroes()
{
    if(OptimizationOptions::verbose)
    cout << "DataCube::fillWithZeroes " << endl;

  //fill memory with zeroes and sync file and memory.
  for (uint64 index = 0; index < getVolume(); index++)
    *(data + index) = 0;
  
  msync(data, getVolume(), MS_SYNC);

  if(OptimizationOptions::verbose)
    cout << "DataCube::fillWithZeroes ok" << endl;
}

#endif

