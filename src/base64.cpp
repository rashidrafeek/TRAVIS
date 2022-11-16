/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm.

    ---------------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/


// This must always be the first include directive
#include "config.h"

#include "base64.h"
#include "string.h"


const char *GetRevisionInfo_base64(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_base64() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


static const char  BASE64_table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
//static const int   BASE64_INPUT_SIZE = 57;


bool isbase64(char c)
{
   return c && strchr(BASE64_table, c) != NULL;
}

inline char BASE64_value(char c)
{
   const char *p = strchr(BASE64_table, c);
   if(p) {
      return (char)(p-BASE64_table);
   } else {
      return 0;
   }
}

int UnBase64(unsigned char *dest, const unsigned char *src, int srclen)
{
   *dest = 0;
   if(*src == 0) 
   {
      return 0;
   }
   unsigned char *p = dest;
   do
   {

      char a = BASE64_value(src[0]);
      char b = BASE64_value(src[1]);
      char c = BASE64_value(src[2]);
      char d = BASE64_value(src[3]);
      *p++ = (a << 2) | (b >> 4);
      *p++ = (b << 4) | (c >> 2);
      *p++ = (c << 6) | d;
      if(!isbase64(src[1])) 
      {
         p -= 2;
         break;
      } 
      else if(!isbase64(src[2])) 
      {
         p -= 2;
         break;
      } 
      else if(!isbase64(src[3])) 
      {
         p--;
         break;
      }
      src += 4;
      while(*src && (*src == 13 || *src == 10)) src++;
   }
   while(srclen-= 4);
   *p = 0;
   return (int)(p-dest);
}


