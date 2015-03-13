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
#ifndef MOLASSEEXCEPTION_H
#define MOLASSEEXCEPTION_H

#include "common.h"

/**
  *@author Daniel Lemire
  */

class MolasseException : public std::exception {

public: 
  typedef enum {GENERAL_ERROR = 0} Code;
  MolasseException(const char *Expression, const char* FileName, unsigned int LineNumber , Code ErrorCode);
  ~MolasseException() throw();
  const char* what() const throw();
  Code getErrorCode();

  static void throwIt( const char*Expression, const char* FileName,	unsigned int LineNumber, Code ErrorCode = GENERAL_ERROR);



private:
  std::string mMessage;
  Code mErrorCode;

};

#define _assert(exp) if( ! (exp) ) { \
          MolasseException::throwIt( #exp, __FILE__, __LINE__); }

#define _assertCode(exp, errCode) if( ! (exp) ) { \
          MolasseException::throwIt( #exp, __FILE__, __LINE__, errCode ); }


#endif
