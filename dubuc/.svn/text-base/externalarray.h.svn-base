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

#ifndef EXTERNALARRAY_H
#define EXTERNALARRAY_H

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/mman.h>

#include <cassert>
//#include "../lemurcore/common.h"
//#include "../lemurcore/fileutil.h"

using namespace std;
typedef unsigned long long uint64;
#include "../lemurcore/fileutil.h"
/*
 * This class can serve as a somewhat unsafe drop in replacement for something like
 * a vector<int> assuming that the file has the proper content
 * and size.
 *
 * Use like this....
 *
 * ExternalArray<double> array(10000000);
 * array[10] = 2;
 * ...
 *
 * The array is cleaned up when it gets garbage collected.
 *
 * A safe copy constructor has not been implemented.
 * 
 */
template <class DataType>
class ExternalArray {
  public:
    ExternalArray(uint64 size, char * FileName = null) : mArraySize(size){
      if(FileName == null) mFileName = (char *) tempnam(null,"molasse");
      else mFileName = FileName;
      //if the file doesn't exist, or if it exists and the file size is zero....make it the right size
      if(!FileUtil::fileExists(mFileName) || 
          (FileUtil::getFileSize(mFileName) < mArraySize * sizeof(DataType)))  {
    
        fstream FileStream( mFileName,std::ios::binary | std::ios::in | std::ios::out| std::ios::trunc);
        assert(FileStream.is_open());
        assert(!FileStream.fail());
        DataType BigArray[1024];
        memset(BigArray,0,1024*sizeof(DataType));
        for(uint times = 0; times * 1024 < mArraySize + 1023;++times)
          FileStream.write((char *) BigArray, sizeof(BigArray));
        FileStream.close();
      }
      mFD = ::open(mFileName, O_RDWR | O_CREAT);
      assert(mFD != -1);
      mData = (DataType *)
        mmap((caddr_t)0, sizeof(DataType)*mArraySize, PROT_READ | PROT_WRITE, MAP_SHARED, mFD, 0);
      assert(mData != MAP_FAILED);
    }

    virtual ~ExternalArray(){
      if (mData != NULL)  {
	//	cerr << "the destructor for external array runs" << endl; // temp
        munmap(mData, sizeof(DataType)*mArraySize);
        ::close(mFD);
        mData = NULL;
      }
      mFileName = NULL;
    }
    const DataType & operator[](uint64 pos) const {return mData[pos];}
    DataType & operator[](uint64 pos) {return mData[pos];}
    virtual uint64 size() const { return mArraySize; }

  protected:
    int mFD; //file descriptor
    DataType *mData; //file data
    char * mFileName; // might be useful to know the name of the file
    uint64 mArraySize; // array size
};

#endif
