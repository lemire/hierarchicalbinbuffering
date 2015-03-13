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

#include "molasseexception.h"
#include <sstream>
MolasseException::MolasseException(const char *Expression, const char* FileName, unsigned int LineNumber , Code ErrorCode){
  stringstream os;
  mErrorCode = ErrorCode;
  os << "Error code: " << mErrorCode << std::endl;
  os << "Broken condition: " << Expression << " in file " << FileName	<< "  at line " << LineNumber;
  mMessage = os.str();
  cout << mMessage << endl;// we do this to make sure the exception is noticed
}

MolasseException::~MolasseException() throw() {
}

void MolasseException::throwIt(const char * Expression, const char * FileName,	unsigned int LineNumber, Code ErrorCode) {
    throw MolasseException( Expression,  FileName, LineNumber, ErrorCode);;
}

const char* MolasseException::what() const throw() {
  return mMessage.c_str();
}

MolasseException::Code MolasseException::getErrorCode() {
  return mErrorCode;
}

