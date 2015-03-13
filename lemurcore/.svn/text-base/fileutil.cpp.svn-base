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
#include "fileutil.h"

//simple function to determine if a file exists
bool FileUtil::fileExists(const char* fileName)
{
   struct stat my_stat;
   return (stat(fileName, &my_stat) == 0);
}

//returns the size of an existing file in bytes
// WARNING: 0 if there is an error or if file doesn't exist
uint64 FileUtil::getFileSize(const char* fileName)
{
  if(fileExists	(fileName)) {
    struct stat my_stat;
    //assumes file exists (probably not the best strategy)// now ok?
    stat(fileName, &my_stat);
    if ((my_stat.st_size) > 0)
      return my_stat.st_size;
    else return 0;
  } else {
          return 0 ;// that's evil, but user probably doesn't use this method to check whether the file exists or not 
  }
}


