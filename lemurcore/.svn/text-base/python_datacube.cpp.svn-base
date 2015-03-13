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
    
#include "python_datacube.h"

    using namespace boost::python;
   
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(open_overloads, open, 0, 1)
    
    BOOST_PYTHON_MODULE(python_datacube)
    {	    
        class_< python_datacube<DataCube<short, long long>, short, long long> >("python_datacube")
      .def(init< list& >())
            .def("open", &python_datacube<DataCube<short, long long>, short, long long>::open, open_overloads())
      .def("close", &python_datacube<DataCube<short, long long>, short, long long>::close)
      .def("get", &python_datacube<DataCube<short, long long>, short, long long>::get)
      .def("put", &python_datacube<DataCube<short, long long>, short, long long>::put)
      .def("fillWithZeroes", &python_datacube<DataCube<short, long long>, short, long long>::fillWithZeroes)
      .def("fillWith", &python_datacube<DataCube<short, long long>, short, long long>::fillWith)
      .def("getShape", &python_datacube<DataCube<short, long long>, short, long long>::getShape)
      .def("getVolume", &python_datacube<DataCube<short, long long>, short, long long>::getVolume)
      .def("prefixSum", &python_datacube<DataCube<short, long long>, short, long long>::prefixSum)
      .def("prefixSumCost", &python_datacube<DataCube<short, long long>, short, long long>::prefixSumCost)
      .def("rangeQuery", &python_datacube<DataCube<short, long long>, short, long long>::rangeQuery)
            ;
      
  
  class_< python_datacube<mmDataCube<short, long long>, short, long long> >("mm_datacube")
      .def(init< list& >())
            .def("open", &python_datacube<mmDataCube<short, long long>, short, long long>::open, open_overloads())
      .def("close", &python_datacube<mmDataCube<short, long long>, short, long long>::close)
      .def("get", &python_datacube<mmDataCube<short, long long>, short, long long>::get)
      .def("put", &python_datacube<mmDataCube<short, long long>, short, long long>::put)
      .def("fillWithZeroes", &python_datacube<mmDataCube<short, long long>, short, long long>::fillWithZeroes)
      .def("fillWith", &python_datacube<mmDataCube<short, long long>, short, long long>::fillWith)
      .def("getShape", &python_datacube<mmDataCube<short, long long>, short, long long>::getShape)
      .def("getVolume", &python_datacube<mmDataCube<short, long long>, short, long long>::getVolume)
      .def("prefixSum", &python_datacube<mmDataCube<short, long long>, short, long long>::prefixSum)
      .def("prefixSumCost", &python_datacube<mmDataCube<short, long long>, short, long long>::prefixSumCost)
      .def("rangeQuery", &python_datacube<mmDataCube<short, long long>, short, long long>::rangeQuery)
            ;


  /*class_< python_datacube<DataCube<float, double>, float, double> >("python_datacube_f")
      .def(init< list& >())
            .def("open", &python_datacube<DataCube<float, double>, float, double>::open, open_overloads())
      .def("close", &python_datacube<DataCube<float, double>, float, double>::close)
      .def("get", &python_datacube<DataCube<float, double>, float, double>::get)
      .def("put", &python_datacube<DataCube<float, double>, float, double>::put)
      .def("fillWithZeroes", &python_datacube<DataCube<float, double>, float, double>::fillWithZeroes)
      .def("fillWith", &python_datacube<DataCube<float, double>, float, double>::fillWith)
      .def("getShape", &python_datacube<DataCube<float, double>, float, double>::getShape)
      .def("getVolume", &python_datacube<DataCube<float, double>, float, double>::getVolume)
      .def("prefixSum", &python_datacube<DataCube<float, double>, float, double>::prefixSum)
      .def("prefixSumCost", &python_datacube<DataCube<float, double>, float, double>::prefixSumCost)   
            .def("rangeQuery", &python_datacube<DataCube<float, double>, float, double>::rangeQuery)
      ;	*/
      
  class_< python_datacube<ChunkedDataCube<short, long long>, short, long long> >("chunked_datacube")
      .def(init< list& , list& >())	
            .def("open", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::open, open_overloads())
      .def("close", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::close)
      .def("get", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::get)
      .def("put", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::put)
      .def("fillWithZeroes", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::fillWithZeroes)
      .def("fillWith", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::fillWith)
      .def("getShape", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::getShape)
      .def("getVolume", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::getVolume)
      .def("prefixSum", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::prefixSum)
      .def("prefixSumCost", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::prefixSumCost)   
      .def("rangeQuery", &python_datacube<ChunkedDataCube<short, long long>, short, long long>::rangeQuery)   
            ;
      
      
  class_< python_datacube<RAMCube<short, long long>, short, long long> >("ram_datacube")
      .def(init< list& >())
            .def("open", &python_datacube<RAMCube<short, long long>, short, long long>::open, open_overloads())
      .def("close", &python_datacube<RAMCube<short, long long>, short, long long>::close)
      .def("get", &python_datacube<RAMCube<short, long long>, short, long long>::get)
      .def("put", &python_datacube<RAMCube<short, long long>, short, long long>::put)
      .def("fillWithZeroes", &python_datacube<RAMCube<short, long long>, short, long long>::fillWithZeroes)
      .def("fillWith", &python_datacube<RAMCube<short, long long>, short, long long>::fillWith)
      .def("getShape", &python_datacube<RAMCube<short, long long>, short, long long>::getShape)
      .def("getVolume", &python_datacube<RAMCube<short, long long>, short, long long>::getVolume)
      .def("prefixSum", &python_datacube<RAMCube<short, long long>, short, long long>::prefixSum)
      .def("prefixSumCost", &python_datacube<RAMCube<short, long long>, short, long long>::prefixSumCost)
      .def("rangeQuery", &python_datacube<RAMCube<short, long long>, short, long long>::rangeQuery)
            ;
  
    }
