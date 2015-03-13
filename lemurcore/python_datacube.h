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
#include <boost/python.hpp>
#include <python2.2/Python.h>
#include <iostream>

#include "datacube.h"
#include "mm_datacube.h"
#include "chunkeddatacube.h"
#include "ramcube.h"
#include "cubicpolynomial.h"

using namespace boost::python;
#include <boost/python/list.hpp>

using namespace std;
#include <vector>

static vector<int> convert(list l);


template <class DataCubeType, class DataType, class LongDataType> 
class python_datacube
{
  public:
    python_datacube() {}
    python_datacube(list &l);
    python_datacube(list &l, list &chunksize);

    void open(const char *filename = null);
    void close();
    
    DataType get(list &l);
    void put(const DataType value, list &l);
    
    void fillWithZeroes();
    void fillWith(const DataType value);
    
    list getShape();
    uint64 getVolume();
    
    LongDataType prefixSum(list &l);
    uint64 prefixSumCost(list &l);
    
    float rangeQuery(list &l);
            
  private:
    list lShape;
    vector<int> vShape;
    DataCubeType * cube;
};


static vector<int> convert(list l)
{
  int list_length = PyList_Size(l.ptr());
  vector<int> shape;
  
  //copy list to vector
  for (int i = 0; i < list_length; i++)
  {
    shape.push_back(PyInt_AsLong(PyList_GetItem(l.ptr(), i)));
  }
  
  return shape;
}



//********************************
//python_datacube Member functions
//********************************



template <class DataCubeType, class DataType, class LongDataType>
python_datacube<DataCubeType, DataType, LongDataType>::python_datacube(list &l) : lShape(l),
cube(new DataCubeType(convert(l)))
{	
  
}


template <class DataCubeType, class DataType, class LongDataType>
python_datacube<DataCubeType, DataType, LongDataType>::python_datacube(list &l, list &chunksize) : lShape(l),
cube(new DataCubeType(convert(l), convert(chunksize)))
{	
  
}


template <class DataCubeType, class DataType, class LongDataType>
void python_datacube<DataCubeType, DataType, LongDataType>::open(const char *filename)
{
  cube->open(filename);
}



template <class DataCubeType, class DataType, class LongDataType>
void python_datacube<DataCubeType, DataType, LongDataType>::close()
{
  cube->close();
}




template <class DataCubeType, class DataType, class LongDataType>
DataType python_datacube<DataCubeType, DataType, LongDataType>::get(list &l)
{
  return cube->get(convert(l));
}



template <class DataCubeType, class DataType, class LongDataType>
void python_datacube<DataCubeType, DataType, LongDataType>::put(const DataType value, list &l)
{
  cube->put(value, convert(l));
}



template <class DataCubeType, class DataType, class LongDataType>
void python_datacube<DataCubeType, DataType, LongDataType>::fillWithZeroes()
{
  cube->fillWithZeroes();
}



template <class DataCubeType, class DataType, class LongDataType>
void python_datacube<DataCubeType, DataType, LongDataType>::fillWith(const DataType value)
{
  cube->fillWith(value);	
}



template <class DataCubeType, class DataType, class LongDataType>
list python_datacube<DataCubeType, DataType, LongDataType>::getShape()
{
  return lShape;
}




template <class DataCubeType, class DataType, class LongDataType>
uint64 python_datacube<DataCubeType, DataType, LongDataType>::getVolume()
{	
  return cube->getVolume();
}




template <class DataCubeType, class DataType, class LongDataType>
LongDataType python_datacube<DataCubeType, DataType, LongDataType>::prefixSum(list &l)
{
  return cube->prefixSum(convert(l));
}


template <class DataCubeType, class DataType, class LongDataType>
uint64 python_datacube<DataCubeType, DataType, LongDataType>::prefixSumCost(list &l)
{
  return cube->prefixSumCost(convert(l));
}


template <class DataCubeType, class DataType, class LongDataType>
float python_datacube<DataCubeType, DataType, LongDataType>::rangeQuery(list &l)
{
  int list_length = PyList_Size(l.ptr());
  float f1, f2, f3, f4;
  int start, end;
  
  vector<RangedCubicPolynomial> polys;
  RangedCubicPolynomial *tempPoly;
  PyObject *templist;
  
  //copy list to vector
  for (int i = 0; i < list_length; i++)
  {	
    templist = PyList_GetItem(l.ptr(), i);
    
    //get first 4 float values
    f1 = PyFloat_AsDouble(PyList_GetItem(templist, 0));
    f2 = PyFloat_AsDouble(PyList_GetItem(templist, 1));
    f3 = PyFloat_AsDouble(PyList_GetItem(templist, 2));
    f4 = PyFloat_AsDouble(PyList_GetItem(templist, 3));
    
    //get start and end
    start = PyInt_AsLong(PyList_GetItem(templist, 4));
    end = PyInt_AsLong(PyList_GetItem(templist, 5));
    
    tempPoly = new RangedCubicPolynomial(f1, f2, f3,f4, start, end);
    
    polys.push_back(*tempPoly);
  }	
  
  return cube->rangeQuery(polys);
}
