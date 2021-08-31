/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2021.

    Please cite:  M. Brehm, M. Thomas: "An Efficient Lossless Compression Algorithm
                  for Trajectories of Atom Positions and Volumetric Data",
                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.

    --------------------------------------------------------------------------------

    LibBQB is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************************/



//=======================================================================
// Copyright (C) 1998-2013 William Hallahan
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//=======================================================================

//**********************************************************************
//  Class Implementation File: BQBLargeInteger.cpp
//  Author: Bill Hallahan
//  Date: March 11, 1998
//
//  Abstract:
//
//    This file contains the implementation for class BQBLargeInteger.
//    An instance of BQBLargeInteger can be used to store and calculate
//    integer values with a huge number of digits. The internal data
//    format is a both a boolean value named m_negative_flag that that
//    stores the sign of the number, and an array of integers that
//    contains a positive binary value. The maximum number of integers
//    that can be contained in the internal array is limited to the
//    maximum signed value which can be stored in an integer, The
//    maximum number of digits that can be used for octal and hexadecimal
//    input, and output, is equal to the the maximum value that can be
//    stored in a signed integer. The input and output conversion
//    methods for decimal numbers supports up to 19,346,699 digits or
//    almost 20 million decimal digits. This implementation only works
//    on little-endian machines.
//
//**********************************************************************


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"

#if defined __GNUC__
#include <stdio.h>
#include <string.h>
#endif

#include "bqb_largeinteger.h"


const char *GetRevisionInfo_bqb_largeinteger(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_largeinteger() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


#define INT_BIT_LENGTH_IS_POWER_OF_TWO


std::string g_sBQBLIString;
std::string g_sBQBLIString2;
std::string g_sBQBLIString3;


//======================================================================
//  Constants.
//======================================================================

namespace
{
    const unsigned int INTEGER_LOW_HALF_MASK = 0xFFFF;
    const unsigned int INTEGER_BIT_COUNT = 32;
    const unsigned int INTEGER_HALF_BIT_COUNT = 16;
    const int INTEGER_MAXIMUM_BASE = 90;
    const double D_MAX_UNSIGNED_VALUE =  4294967296.0;
    
    const char DIGIT_ARRAY[INTEGER_MAXIMUM_BASE] =
    {
        '0', '1', '2', '3', '4', '5', '6', '7',
        '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
        'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
        'O', 'P', 'Q', 'R', 'S', 'T', 'W', 'X',
        'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
        'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
        'o', 'p', 'q', 'r', 's', 't', 'w', 'x',
        'y', 'z', '.', ',', '?', '/', '\\', '|',
        '<', '>', ':', ';', '"', '\'', '{', '}',
        '[', ']', '`', '~', '!', '@', '#', '$',
        '%', '^', '&', '*', '(', ')', '-', '_',
        '+', '='
    };

    const unsigned int BIT_MASK_ARRAY[5] =
    {
        0xAAAAAAAA,
        0xCCCCCCCC,
        0xF0F0F0F0,
        0xFF00FF00,
        0xFFFF0000
    };

    const unsigned int BIT_MASK_TABLE_LENGTH =
        sizeof(BIT_MASK_ARRAY) / sizeof(unsigned int);

#ifdef LARGE_INT_ALTERNATE_POP_COUNT
    const unsigned int MASK_1 = 011111111111;
    const unsigned int MASK_3 = 033333333333;
    const unsigned int MASK_7 = 030707070707;
#endif
}

//======================================================================
//  istream operator for input.
//======================================================================

std::istream & operator >>(std::istream & is, BQBLargeInteger & value)
{
    //------------------------------------------------------------------
    //  Get the input string.
    //------------------------------------------------------------------

    const unsigned int max_string_length = 2048;
    char * number_ptr = new char [max_string_length];

    if (number_ptr != 0)
    {
        is.get(number_ptr, max_string_length - 1);
        int length = (int)is.gcount();
        number_ptr[length] = '\0';

        unsigned int base;

        switch (is.flags() & std::ios_base::basefield)
        {
        case std::ios_base::oct:
            base = 8;
            break;
        case std::ios_base::dec:
            base = 10;
            break;
        case std::ios_base::hex:
            base = 16;
            break;
        default:
            base = 10;
            break;
        }

        value.SetValue(number_ptr, base);

        delete [] number_ptr;
    }

    return is;
}

//======================================================================
//  ostream operator for output.
//======================================================================

std::ostream & operator <<(std::ostream & os, const BQBLargeInteger & value)
{
    unsigned int base;

    switch (os.flags() & std::ios_base::basefield)
    {
    case std::ios_base::oct:
        base = 8;
        break;
    case std::ios_base::dec:
        base = 10;
        break;
    case std::ios_base::hex:
        base = 16;
        break;
    default:
        base = 10;
        break;
    }

    std::string number_string;
    value.GetNumberString(number_string, base);
    os << number_string.c_str();

    return os;
}

//======================================================================
//  Constructor: BQBLargeInteger::BQBLargeInteger
//======================================================================

BQBLargeInteger::BQBLargeInteger()
  : m_array_length(0)
  , m_integer_length(0)
  , m_default_base(10)
  , m_array_ptr(0)
  , m_negative_flag(false)
{
    SetIntegerLength(1);
    m_array_ptr[0] = 0;
}

//======================================================================
//  Copy constructor
//======================================================================

BQBLargeInteger::BQBLargeInteger(const BQBLargeInteger & value)
  : m_array_length(0)
  , m_integer_length(0)
  , m_default_base(10)
  , m_array_ptr(0)
  , m_negative_flag(false)
{
    Copy(value);
}

//======================================================================
//  Conversion constructor for signed long.
//======================================================================

BQBLargeInteger::BQBLargeInteger(long value)
  : m_array_length(0)
  , m_integer_length(0)
  , m_default_base(10)
  , m_array_ptr(0)
  , m_negative_flag(false)
{
    SetIntegerLength(1);
    m_negative_flag = value < 0;
    m_array_ptr[0] =
        m_negative_flag ? (unsigned int)(-value) : (unsigned int)(value);
}

//======================================================================
//  Conversion constructor for unsigned long.
//======================================================================

BQBLargeInteger::BQBLargeInteger(unsigned long value)
  : m_array_length(0)
  , m_integer_length(0)
  , m_default_base(10)
  , m_array_ptr(0)
  , m_negative_flag(false)
{
    SetIntegerLength(1);
    m_array_ptr[0] = value;
}

//======================================================================
//  Conversion constructor for signed int.
//======================================================================

BQBLargeInteger::BQBLargeInteger(int value)
  : m_array_length(0)
  , m_integer_length(0)
  , m_default_base(10)
  , m_array_ptr(0)
  , m_negative_flag(false)
{
    SetIntegerLength(1);
    m_negative_flag = value < 0;
    m_array_ptr[0] =
        m_negative_flag ? (unsigned int)(-value) : (unsigned int)(value);
}

//======================================================================
//  Conversion constructor for unsigned int.
//======================================================================

BQBLargeInteger::BQBLargeInteger(unsigned int value)
  : m_array_length(0)
  , m_integer_length(0)
  , m_default_base(10)
  , m_array_ptr(0)
  , m_negative_flag(false)
{
    SetIntegerLength(1);
    m_array_ptr[0] = value;
}

//======================================================================
//  Constructor: BQBLargeInteger::BQBLargeInteger
//
//  The number is interpreted using the default base that is set
//  using the SetDefaultBase method.
//======================================================================

BQBLargeInteger::BQBLargeInteger(const char * psznumber_ptr)
: m_array_length(0)
, m_integer_length(0)
, m_default_base(10)
, m_array_ptr(0)
, m_negative_flag(false)
{
    SetValue(psznumber_ptr, m_default_base);
}

//======================================================================
//  Special constructor to allow initializing using an array of
//  binary data.
//======================================================================

BQBLargeInteger::BQBLargeInteger(const char * binary_data_ptr,
                           unsigned int length)
  : m_array_length(0)
  , m_integer_length(0)
  , m_default_base(10)
  , m_array_ptr(0)
  , m_negative_flag(false)
{
    SetBinaryValue(binary_data_ptr, length);
}

//======================================================================
//  Destructor
//======================================================================

BQBLargeInteger::~BQBLargeInteger()
{
    delete [] m_array_ptr;
}

//======================================================================
// Method: BQBLargeInteger::SetDefaultBase
//
// The default base is the value used for the constructor and operator
// equals methods that take only a string pointer for a number.
//======================================================================

void BQBLargeInteger::SetDefaultBase(unsigned int default_base)
{
    if ((default_base == 0) || (default_base > ((unsigned int)INTEGER_MAXIMUM_BASE)))
    {
        //throw BQBLargeIntegerException("Illegal base");
		printf("BQBLargeInteger::SetDefaultBase(): Error: Illegal base %u.\n",default_base);
		abort();
    }

    m_default_base = default_base;
}

//======================================================================
//  This method allows getting this large integer's binary value.
//  The binary value is a positive number in little endian format.
//  The buffer that is passed must be large enough to contain the
//  number. This method can be called passing a null pointer for
//  the buffer to just obtain the buffer length in bytes.
//======================================================================

unsigned int BQBLargeInteger::GetBinaryValue(unsigned char * binary_data_ptr) const
{
    //------------------------------------------------------------------
    //  Get the length of the output binary array.
    //------------------------------------------------------------------

    unsigned int length = sizeof(unsigned int) * m_integer_length;

    //------------------------------------------------------------------
    //  Discard leading zero bytes.
    //------------------------------------------------------------------

    char * source_ptr = reinterpret_cast<char *>(&m_array_ptr[0]);
    unsigned int leading_zero_count = 0;

    for (int i = length - 1; i >= 0; --i)
    {
        if (source_ptr[i] == 0)
        {
            ++leading_zero_count;
        }
        else
        {
            break;
        }
    }

    length -= leading_zero_count;

    //------------------------------------------------------------------
    //  If the passed array pointer is not zero then copy the data
    //  to the output array.
    //------------------------------------------------------------------

    if (binary_data_ptr != 0)
    {
        for (unsigned int j = 0; j < length; ++j)
        {
            binary_data_ptr[j] = source_ptr[j];
        }
    }

    return length;
}

//======================================================================
//  This method allows setting this large integer's binary value.
//  The binary value is a positive number in little endian format.
//  The buffer that is passed must be large enough to contain the
//  number. This method can be called passing a null pointer for
//  the buffer to just obtain the buffer length in bytes.
//======================================================================

void BQBLargeInteger::SetBinaryValue(const char * binary_data_ptr,
                                  unsigned int length)
{
    //------------------------------------------------------------------
    //  Set the sign for a positive number.
    //------------------------------------------------------------------

    m_negative_flag = false;

    //------------------------------------------------------------------
    //  Set the length of internal array of unsigned integers
    //  that will store the byte array.
    //------------------------------------------------------------------

    unsigned int array_length = (length + 3) >> 2;
    SetIntegerLength(array_length);

    //------------------------------------------------------------------
    //  Copy the bytes to the internal unsigned integer array.
    //------------------------------------------------------------------

    char * destination_ptr = reinterpret_cast<char *>(m_array_ptr);
    unsigned int i = 0;

    for (i = 0; i < length; ++i)
    {
        destination_ptr[i] = binary_data_ptr[i];
    }

    //------------------------------------------------------------------
    //  Determine if there are any extra bytes in the array that
    //  need to be zero padded.
    //------------------------------------------------------------------

    unsigned int index_one_past_final_byte = m_integer_length << 2;

    if (length < index_one_past_final_byte)
    {
        for (i = length; i < index_one_past_final_byte; ++i)
        {
            destination_ptr[i] = 0;
        }
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    Normalize();

    return;
}

//==================================================================
//  Get this large integer value represented in an stl string.
//==================================================================

void BQBLargeInteger::GetNumberString(std::string & number_string,
                                   unsigned int base) const
{
    //------------------------------------------------------------------
    //  Convert this large integer to a character string representation
    //  in the desired base.
    //------------------------------------------------------------------

    number_string.erase();

    //------------------------------------------------------------------
    //  Use special code for bases 8, 10, and 16.
    //------------------------------------------------------------------

    switch (base)
    {
    case 8:

        GetBase8NumberString(number_string);
        break;

    case 10:

        GetBase10NumberString(number_string);
        break;

    case 16:

        GetBase16NumberString(number_string);
        break;

    default:

        //------------------------------------------------------------------
        //  Generic code for converting to any base up to base 90.
        //------------------------------------------------------------------

        BQBLargeInteger x_temp = *this;

        //--------------------------------------------------------------
        //  The characters are converted in reverse order from the
        //  order they are displayed.
        //--------------------------------------------------------------

        do
        {
            BQBLargeInteger x_mult = x_temp / base;
            BQBLargeInteger x_part = x_mult * base;
            unsigned int digit_index = (unsigned int)(x_temp - x_part);
            number_string += DIGIT_ARRAY[digit_index];
            x_temp = x_mult;
        }
        while(!x_temp.IsZero());

        //------------------------------------------------------------------
        //  Reverse the order of the digits so the number can be
        //  properly displayed.
        //------------------------------------------------------------------

        char * number_ptr = const_cast<char *>(number_string.data());
        unsigned int length = (unsigned int)number_string.length();

        for (unsigned int i = 0; i < (length >> 1); ++i)
        {
            char temp = number_ptr[i];
            number_ptr[i] = number_ptr[(length - i) - 1];
            number_ptr[(length - i) - 1] = temp;
        }

        break;
    }

    //------------------------------------------------------------------
    //  If the number is negative then make the first character a
    //  minus sign.
    //------------------------------------------------------------------

    if (IsNegative())
    {
        number_string.insert(0, "-");
    }

    return;
}

//==================================================================
//  Set this large integer value using the number in a character
//  string.
//==================================================================

bool BQBLargeInteger::SetValue(const char * number_ptr, unsigned int base)
{
    //------------------------------------------------------------------
    //  Set this instance to the value zero.
    //------------------------------------------------------------------

    SetToZero();

    //------------------------------------------------------------------
    //  If the passed string pointer is equal to zero then exit.
    //------------------------------------------------------------------

    bool success_flag = number_ptr != 0;

    if (success_flag)
    {
        size_t length = ::strlen(number_ptr);

        if (length > 0)
        {
            //----------------------------------------------------------
            //  Skip any leading white space.
            //----------------------------------------------------------

            size_t index = 0;

            for (index = 0; index < length; ++index)
            {
                if (!::isspace((int)(number_ptr[index])))
                {
                    break;
                }
            }

            //----------------------------------------------------------
            //  Test for + or - character.
            //----------------------------------------------------------

            m_negative_flag = number_ptr[index] == '-';

            if (m_negative_flag)
            {
                ++index;
            }

            if (number_ptr[index] == '+')
            {
                ++index;
            }

            //----------------------------------------------------------
            //  Skip white space.
            //----------------------------------------------------------

            for (;index < length; ++index)
            {
                if (!::isspace((int)(number_ptr[index])))
                {
                    break;
                }
            }

            //----------------------------------------------------------
            //  Use special code for bases 8, 10, and 16.
            //----------------------------------------------------------

            switch (base)
            {
            case 8:

                success_flag = SetValueWithBase8String(number_ptr + index);
                break;

            case 10:

                success_flag = SetValueWithBase10String(number_ptr + index);
                break;

            case 16:

                success_flag = SetValueWithBase16String(number_ptr + index);
                break;

            default:
                {
                    //--------------------------------------------------
                    //  Loop and add each digit to the value.
                    //--------------------------------------------------

                    success_flag = true;

                    for (size_t i = index; i < length; ++i)
                    {
                        //----------------------------------------------
                        //  Get the numeric value of a digit.
                        //----------------------------------------------

                        int digit = ::tolower((int)(number_ptr[i]));

                        if (::isxdigit(digit))
                        {
                            if ((digit >= 'a') && (digit <= 'f'))
                            {
                                digit = digit + 10 - (int)('a');
                            }
                            else
                            {
                                digit = digit - (int)('0');
                            }
                            
                            if ((unsigned int)(digit) < base)
                            {
                                operator *=(base);
                                operator +=(digit);
                            }
                            else
                            {
                                success_flag = false;
                                break;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                break;
            }
        }
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    Normalize();

    return success_flag;    
}

//======================================================================
//  This method returns the position of the leading bit in this
//  large integer. If the value of this large integer is zero then
//  this function will return the value zero.
//======================================================================

unsigned int BQBLargeInteger::LeadingBitPosition() const
{
    unsigned int leading_bit_position = 0;

    //------------------------------------------------------------------
    //  The number has been normalized, so there must be a bit set
    //  in either the most significant term or the next-to most
    //  significant term.
    //------------------------------------------------------------------

    int leading_value_index = m_integer_length - 1;

    if (m_array_ptr[m_integer_length - 1] == 0)
    {
        leading_value_index--;
    }

    if (leading_value_index > 0)
    {
        unsigned int leading_value = m_array_ptr[leading_value_index];

        //------------------------------------------------------------------
        //  Find the position of the leading bit in the leading value.
        //------------------------------------------------------------------

        unsigned int bit_position_addend = INTEGER_HALF_BIT_COUNT;

        for (unsigned int i = BIT_MASK_TABLE_LENGTH - 1; int(i) >= 0; --i)
        {
            unsigned int bit_mask = BIT_MASK_ARRAY[i];

            if ((leading_value & bit_mask) != 0)
            {
                leading_bit_position += bit_position_addend;
                leading_value = leading_value & bit_mask;
                bit_position_addend = bit_position_addend >> 1;
            }
        }

        leading_bit_position = (INTEGER_BIT_COUNT * leading_value_index) + leading_bit_position;
    }


    return leading_bit_position;
}

//======================================================================
//  Return the count of the number of bits set in this large
//  integer value.
//======================================================================

unsigned int BQBLargeInteger::PopulationCount() const
{
    unsigned int bit_count = 0;

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        bit_count += PopulationCount32Bit(m_array_ptr[i]);
    }

    return bit_count;
}

//======================================================================
//  operator = for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator =(const BQBLargeInteger & value)
{
    Copy(value);
    return *this;
}

BQBLargeInteger BQBLargeInteger::operator =(long value)
{
    Copy(BQBLargeInteger(value));
    return *this;
}

BQBLargeInteger BQBLargeInteger::operator =(unsigned long value)
{
    Copy(BQBLargeInteger(value));
    return *this;
}

BQBLargeInteger BQBLargeInteger::operator =(int value)
{
    Copy(BQBLargeInteger(value));
    return *this;
}

BQBLargeInteger BQBLargeInteger::operator =(unsigned int value)
{
    Copy(BQBLargeInteger(value));
    return *this;
}

BQBLargeInteger BQBLargeInteger::operator =(const char * psznumber_ptr)
{
    SetValue(psznumber_ptr, m_default_base);
    return *this;
}

//======================================================================
//  Unary operator + for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator +=(const BQBLargeInteger & addend)
{
    //------------------------------------------------------------------
    //  If the signs of this large integer and the sign of the addend
    //  are the same then add the positive array values. If the signs
    //  are NOT the same then subtract the positive array for the
    //  addend from the positive array for this large integer.
    //------------------------------------------------------------------

    bool negative_flag_0 = IsNegative();
    bool negative_flag_1 = addend.IsNegative();

    if (negative_flag_0 ^ negative_flag_1)
    {
        m_negative_flag = negative_flag_0 ^ SubtractPositiveArray(addend);
    }
    else
    {
        AddPositiveArray(addend);
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    Normalize();

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator +=(long addend)
{
    return operator +=(BQBLargeInteger(addend));
}

BQBLargeInteger BQBLargeInteger::operator +=(unsigned long addend)
{
    return operator +=(BQBLargeInteger(addend));
}

BQBLargeInteger BQBLargeInteger::operator +=(int addend)
{
    return operator +=(BQBLargeInteger(addend));
}

BQBLargeInteger BQBLargeInteger::operator +=(unsigned int addend)
{
    return operator +=(BQBLargeInteger(addend));
}

//======================================================================
//  Unary operator - for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator -=(const BQBLargeInteger & subtrahend)
{
    //------------------------------------------------------------------
    //  If the signs of this large integer and the sign of the
    //  subtrahend are the same then subtract the positive array
    //  values of the subtrahend from the positive array values
    //  for this large integer. If the signs are NOT the same
    //  then add the positive array values.
    //------------------------------------------------------------------

    bool negative_flag_0 = IsNegative();
    bool negative_flag_1 = subtrahend.IsNegative();

    if (negative_flag_0 ^ negative_flag_1)
    {
        AddPositiveArray(subtrahend);
    }
    else
    {
        m_negative_flag = negative_flag_0 ^ SubtractPositiveArray(subtrahend);
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    Normalize();

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator -=(long subtrahend)
{
    return operator -=(BQBLargeInteger(subtrahend));
}

BQBLargeInteger BQBLargeInteger::operator -=(unsigned long subtrahend)
{
    return operator -=(BQBLargeInteger(subtrahend));
}

BQBLargeInteger BQBLargeInteger::operator -=(int subtrahend)
{
    return operator -=(BQBLargeInteger(subtrahend));
}

BQBLargeInteger BQBLargeInteger::operator -=(unsigned int subtrahend)
{
    return operator -=(BQBLargeInteger(subtrahend));
}

//======================================================================
//  Unary operator * for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator *=(const BQBLargeInteger & x_multiplier)
{
    //------------------------------------------------------------------
    //  Create a temporary variable to contain the product and set
    //  the sign of the product.
    //------------------------------------------------------------------

    BQBLargeInteger product = BQBLargeInteger(0);

    //------------------------------------------------------------------
    //  Multiply using base 65536 digits.
    //------------------------------------------------------------------

    product.SetIntegerLength(m_integer_length + x_multiplier.m_integer_length);

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        //--------------------------------------------------------------
        //  Get the base 65536 digits for the first multiplier.
        //--------------------------------------------------------------

        unsigned int a0 = m_array_ptr[i] & INTEGER_LOW_HALF_MASK;
        unsigned int a1 = m_array_ptr[i] >> INTEGER_HALF_BIT_COUNT;

        for (unsigned int j = 0; j < x_multiplier.m_integer_length; ++j)
        {
            //----------------------------------------------------------
            //  Get the base 65536 digits for the second multiplier.
            //----------------------------------------------------------

            unsigned int b0 = x_multiplier.m_array_ptr[j] & INTEGER_LOW_HALF_MASK;
            unsigned int b1 = x_multiplier.m_array_ptr[j] >> INTEGER_HALF_BIT_COUNT;

            //----------------------------------------------------------
            //  Calculate the product as shown. All multiplier
            //  variables contain INTEGER_HALF_BIT_COUNT
            //  bit quantities. The partial products are the
            //  length of an unsigned integer.
            //
            //
            //
            //                          b1    b0
            //                     X    a1    a0
            //          -------------------------
            //                       a0b0H a0b0L 
            //                 a0b1H a0b1L
            //                 a1b0H a1b0L
            //           a1b1H a1b1L
            //          -------------------------
            //             p1H   p1L   p0H   p0L
            //
            //----------------------------------------------------------

            //------------------------------------------------------------------
            //  Calculate the 64-bit partial products.
            //------------------------------------------------------------------

            unsigned int a0b0 = a0 * b0;
            unsigned int a0b1 = a0 * b1;
            unsigned int a1b0 = a1 * b0;
            unsigned int a1b1 = a1 * b1;

            //------------------------------------------------------------------
            //  Add the partial products to obtain the low dword.
            //------------------------------------------------------------------

            unsigned int a0b1L = a0b1 << INTEGER_HALF_BIT_COUNT;
            unsigned int p0 = a0b0 + a0b1L;
            unsigned int low_dword_carry = (p0 < a0b1L) ? 1 : 0;
            unsigned int a1b0L = a1b0 << INTEGER_HALF_BIT_COUNT;
            unsigned int temp = p0 + a1b0L;
            low_dword_carry += (temp < p0) ? 1 : 0;
            p0 = temp;

            //------------------------------------------------------------------
            //  Add the partial products to obtain the high dword.
            //  The high dword cannot carry because the largest 32-bit
            //  product is:
            //
            //  0xFFFFFFFE00000001 = 0xFFFFFFFF * 0xFFFFFFFF
            //------------------------------------------------------------------

            unsigned int a0b1H = a0b1 >> INTEGER_HALF_BIT_COUNT;
            unsigned int p1 = a1b1 + a0b1H;
            unsigned int a1b0H = a1b0 >> INTEGER_HALF_BIT_COUNT;
            p1 = p1 + a1b0H + low_dword_carry;

            //----------------------------------------------------------
            //  The 64 bit partial product has been calculated.
            //  Accumulate the product into the appropriate
            //  locations in the output integer array.
            //----------------------------------------------------------

            //----------------------------------------------------------
            //  Accumulate the the low 32 bit sum.
            //----------------------------------------------------------

            unsigned int index = i + j;
            AccumulateWithCarry(product, index, p0);

            //----------------------------------------------------------
            //  Accumulate the the high 32 bit sum.
            //----------------------------------------------------------

            AccumulateWithCarry(product, index + 1, p1);
        }
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    product.Normalize();

    //------------------------------------------------------------------
    //  Set the sign of the product.
    //------------------------------------------------------------------

    product.m_negative_flag = IsNegative() ^ x_multiplier.IsNegative();

    //------------------------------------------------------------------
    //  Copy the final product to this integer.
    //------------------------------------------------------------------

    Copy(product);

    return *this;
}

//======================================================================
//  Method BQBLargeInteger::AccumulateWithCarry()
//  This protected method is called from the
//  unary operator *(const BQBLargeInteger &) method.
//======================================================================

void BQBLargeInteger::AccumulateWithCarry(BQBLargeInteger &product,
                                       int index,
                                       unsigned int value)
{
    bool sum_overflows = true;

    while ((sum_overflows) && (index < (int)(product.m_integer_length)))
    {
        unsigned int temp = product.m_array_ptr[index];
        temp = temp + value;
        sum_overflows = temp < value;
        value = 1;
        product.m_array_ptr[index] = temp;
        ++index;
    }
}

BQBLargeInteger BQBLargeInteger::operator *=(long multiplier)
{
    return operator *=(BQBLargeInteger(multiplier));
}

BQBLargeInteger BQBLargeInteger::operator *=(unsigned long multiplier)
{
    return operator *=(BQBLargeInteger(multiplier));
}

BQBLargeInteger BQBLargeInteger::operator *=(int multiplier)
{
    return operator *=(BQBLargeInteger(multiplier));
}

BQBLargeInteger BQBLargeInteger::operator *=(unsigned int multiplier)
{
    return operator *=(BQBLargeInteger(multiplier));
}

//======================================================================
//  Unary operator / for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator /=(const BQBLargeInteger & divisor)
{

	printf("BQBLargeInteger::operator /=: Warning: Division seems to yield wrong results!\n");
	abort();

    //------------------------------------------------------------------
    //  If the divisor is equal to zero then cause a divide by zero
    //  exception.
    //------------------------------------------------------------------

    if (divisor.IsZero())
    {
        //throw BQBLargeIntegerException("Divide by zero");
        // Comment out the line above and uncomment the line below
        // to cause a regular integer divide by zero exception.
        //unsigned int temp = 1 / divisor.m_array_ptr[0];
		printf("BQBLargeInteger::operator /=: Error: Divide by zero.\n");
		abort();
    }

    //--------------------------------------------------------------
    //  Determine the sign of the quotient.
    //--------------------------------------------------------------

    bool quotient_sign = IsNegative() ^ divisor.IsNegative();

    //------------------------------------------------------------------
    //  Set the numerator variable to this object's value so that
    //  this object can return the quotient.
    //------------------------------------------------------------------

    BQBLargeInteger numerator = *this;

    // Make this numerator positive.
    numerator.m_negative_flag = false;

    //------------------------------------------------------------------
    //  Set the initial value of the quotient to zero.
    //------------------------------------------------------------------

    SetToZero();

    //------------------------------------------------------------------
    //  If the numerator is zero or the the denominator is greater than
    //  the numerator then the quotient is zero.
    //------------------------------------------------------------------

    if (!numerator.IsZero())
    {
        //--------------------------------------------------------------
        //  Shift the denominator to the left until it becomes the
        //  greatest shifted value that is less than or equal to the
        //  numerator.
        //--------------------------------------------------------------

        int scale_shift = 0;
        BQBLargeInteger denom = divisor;

        while (denom < numerator)
        {
            denom.ShiftLeft(1);
            ++scale_shift;
        }

        //--------------------------------------------------------------
        //  Perform long division.
        //--------------------------------------------------------------

        while (numerator >= denom)
        {
            BQBLargeInteger difference = numerator - denom;

            denom.ShiftRight(1);
            ShiftLeft(1);

            if (!difference.IsNegative())
            {
                ++m_array_ptr[0];
                numerator = difference;
            }
        }

        //--------------------------------------------------------------
        //  Adjust the numerator for the denominator scaling.
        //--------------------------------------------------------------

        *this <<= scale_shift;

        //--------------------------------------------------------------
        //  Normalize the result.
        //--------------------------------------------------------------

        Normalize();

        m_negative_flag = quotient_sign;
    }

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator /=(long divisor)
{
    return operator /=(BQBLargeInteger(divisor));
}

BQBLargeInteger BQBLargeInteger::operator /=(unsigned long divisor)
{
    return operator /=(BQBLargeInteger(divisor));
}

BQBLargeInteger BQBLargeInteger::operator /=(int divisor)
{
    return operator /=(BQBLargeInteger(divisor));
}

BQBLargeInteger BQBLargeInteger::operator /=(unsigned int divisor)
{
    return operator /=(BQBLargeInteger(divisor));
}

//======================================================================
//  operator <<= for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator <<=(const BQBLargeInteger & shift_count)
{
    if (shift_count.FitsIn32Bits())
    {
        if (shift_count.IsNegative())
        {
            ShiftRight(shift_count.m_array_ptr[0]);
        }
        else
        {
            ShiftLeft(shift_count.m_array_ptr[0]);
        }
    }
    else
    {
        SetToZero();
    }

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator <<=(long shift_count)
{
    return operator <<=(BQBLargeInteger(shift_count));
}

BQBLargeInteger BQBLargeInteger::operator <<=(unsigned long shift_count)
{
    return operator <<=(BQBLargeInteger(shift_count));
}

BQBLargeInteger BQBLargeInteger::operator <<=(int shift_count)
{
    return operator <<=(BQBLargeInteger(shift_count));
}

BQBLargeInteger BQBLargeInteger::operator <<=(unsigned int shift_count)
{
    return operator <<=(BQBLargeInteger(shift_count));
}

//======================================================================
//  operator >>= for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator >>=(const BQBLargeInteger & shift_count)
{
    if (shift_count.FitsIn32Bits())
    {
        if (shift_count.IsNegative())
        {
            ShiftLeft(shift_count.m_array_ptr[0]);
        }
        else
        {
            ShiftRight(shift_count.m_array_ptr[0]);
        }
    }
    else
    {
        SetToZero();
    }

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator >>=(long shift_count)
{
    return operator >>=(BQBLargeInteger(shift_count));
}

BQBLargeInteger BQBLargeInteger::operator >>=(unsigned long shift_count)
{
    return operator >>=(BQBLargeInteger(shift_count));
}

BQBLargeInteger BQBLargeInteger::operator >>=(int shift_count)
{
    return operator >>=(BQBLargeInteger(shift_count));
}

BQBLargeInteger BQBLargeInteger::operator >>=(unsigned int shift_count)
{
    return operator >>=(BQBLargeInteger(shift_count));
}

//======================================================================
//  Unary operator %= for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator %=(const BQBLargeInteger & divisor)
{
    BQBLargeInteger x_temp(*this);
    operator /=(divisor);
    operator *=(divisor);
    return operator -=(x_temp);
}

BQBLargeInteger BQBLargeInteger::operator %=(long divisor)
{
    return operator %=(BQBLargeInteger(divisor));
}

BQBLargeInteger BQBLargeInteger::operator %=(unsigned long divisor)
{
    return operator %=(BQBLargeInteger(divisor));
}

BQBLargeInteger BQBLargeInteger::operator %=(int divisor)
{
    return operator %=(BQBLargeInteger(divisor));
}

BQBLargeInteger BQBLargeInteger::operator %=(unsigned int divisor)
{
    return operator %=(BQBLargeInteger(divisor));
}

//======================================================================
//  operator ^= for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator ^=(const BQBLargeInteger & value)
{
    unsigned int largest_integer_length =
        m_integer_length < value.m_integer_length ? value.m_integer_length : m_integer_length;

    SetIntegerLength(largest_integer_length);

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        m_array_ptr[i] ^= value.GetSafeArrayValue(i);
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    Normalize();

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator ^=(long value)
{
    return operator %=(BQBLargeInteger(value));
}

BQBLargeInteger BQBLargeInteger::operator ^=(unsigned long value)
{
    return operator %=(BQBLargeInteger(value));
}

BQBLargeInteger BQBLargeInteger::operator ^=(int value)
{
    return operator %=(BQBLargeInteger(value));
}

BQBLargeInteger BQBLargeInteger::operator ^=(unsigned int value)
{
    return operator %=(BQBLargeInteger(value));
}

//======================================================================
//  operator &= for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator &=(const BQBLargeInteger & value)
{
    unsigned int largest_integer_length =
        m_integer_length < value.m_integer_length ? value.m_integer_length : m_integer_length;

    SetIntegerLength(largest_integer_length);

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        m_array_ptr[i] &= value.GetSafeArrayValue(i);
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    Normalize();

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator &=(unsigned long value)
{
    return operator &=(BQBLargeInteger(value));
}

BQBLargeInteger BQBLargeInteger::operator &=(int value)
{
    return operator &=(BQBLargeInteger(value));
}

BQBLargeInteger BQBLargeInteger::operator &=(unsigned int value)
{
    return operator &=(BQBLargeInteger(value));
}

//======================================================================
//  operator |= for integral types.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator |=(const BQBLargeInteger & value)
{
    unsigned int largest_integer_length =
        m_integer_length < value.m_integer_length ? value.m_integer_length : m_integer_length;

    SetIntegerLength(largest_integer_length);

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        m_array_ptr[i] |= value.GetSafeArrayValue(i);
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    Normalize();

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator |=(long value)
{
    return operator |=(BQBLargeInteger(value));
}

BQBLargeInteger BQBLargeInteger::operator |=(unsigned long value)
{
    return operator |=(BQBLargeInteger(value));
}

BQBLargeInteger BQBLargeInteger::operator |=(int value)
{
    return operator |=(BQBLargeInteger(value));
}

BQBLargeInteger BQBLargeInteger::operator |=(unsigned int value)
{
    return operator |=(BQBLargeInteger(value));
}

//======================================================================
//  operators with no arguments. The * operator and the & operator
//  with no arguments do not need to be overloaded.
//======================================================================

BQBLargeInteger BQBLargeInteger::operator !()
{
    if (IsZero())
    {
        m_array_ptr[0] = (unsigned int)(-1);
        m_integer_length = 1;
    }
    else
    {
        SetToZero();
    }
    
    return *this;
}

BQBLargeInteger BQBLargeInteger::operator ~()
{
    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        m_array_ptr[i] = ~m_array_ptr[i];
    }

    //------------------------------------------------------------------
    //  Normalize the result.
    //------------------------------------------------------------------

    Normalize();

    return *this;
}

BQBLargeInteger BQBLargeInteger::operator +()
{
    return *this;
}

BQBLargeInteger BQBLargeInteger::operator -()
{
    m_negative_flag = ! m_negative_flag;
    return *this;
}

//======================================================================
//  The prefix form of the increment and decrement operators.
//======================================================================

const BQBLargeInteger BQBLargeInteger::operator ++()
{
    operator +=(1);
    return *this;
}

const BQBLargeInteger BQBLargeInteger::operator --()
{
    operator -=(1);
    return *this;
}

//======================================================================
//  The postfix form of the increment and decrement operators.
//======================================================================

const BQBLargeInteger BQBLargeInteger::operator ++(int)
{
    BQBLargeInteger x_temp = *this;
    operator +=(1);
    return x_temp;
}

const BQBLargeInteger BQBLargeInteger::operator --(int)
{
    BQBLargeInteger x_temp = *this;
    operator -=(1);
    return x_temp;
}

//======================================================================
//  Unary operator == for integral types.
//======================================================================

bool BQBLargeInteger::operator ==(const BQBLargeInteger & value) const
{
    unsigned int largest_integer_length =
        m_integer_length < value.m_integer_length ? value.m_integer_length : m_integer_length;

    bool bEqual = true;

    for (unsigned int i = largest_integer_length - 1; (int)(i) >= 0; --i)
    {
        if (GetSafeArrayValue(i) != value.GetSafeArrayValue(i))
        {
            bEqual = false;
            break;
        }
    }

    return bEqual;
}

bool BQBLargeInteger::operator ==(long value) const
{
    return operator ==(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator ==(unsigned long value) const
{
    return operator ==(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator ==(int value) const
{
    return operator ==(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator ==(unsigned int value) const
{
    return operator ==(BQBLargeInteger(value));
}

//======================================================================
//  Unary operator != for integral types.
//======================================================================

bool BQBLargeInteger::operator !=(const BQBLargeInteger & value) const
{
    return ! operator ==(value);
}

bool BQBLargeInteger::operator !=(long value) const
{
    return operator !=(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator !=(unsigned long value) const
{
    return operator !=(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator !=(int value) const
{
    return operator !=(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator !=(unsigned int value) const
{
    return operator !=(BQBLargeInteger(value));
}

//======================================================================
//  Unary operator < for integral types.
//======================================================================

bool BQBLargeInteger::operator <(const BQBLargeInteger & value) const
{
    //------------------------------------------------------------------
    //  If this value is negative and value is positive then
    //  return true.
    //------------------------------------------------------------------

    bool is_less_than_flag = false;

    if (IsNegative())
    {
        if (!value.IsNegative())
        {
            is_less_than_flag = true;
        }
        else
        {
            is_less_than_flag = LessThanPositiveArrayCompare(value, *this);
        }
    }
    else
    {
        is_less_than_flag = LessThanPositiveArrayCompare(*this, value);
    }

    return is_less_than_flag;
}

bool BQBLargeInteger::operator <(long value) const
{
    return operator <(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator <(unsigned long value) const
{
    return operator <(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator <(int value) const
{
    return operator <(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator <(unsigned int value) const
{
    return operator <(BQBLargeInteger(value));
}

//======================================================================
//  Unary operator > for integral types.
//======================================================================

bool BQBLargeInteger::operator >(const BQBLargeInteger & value) const
{
    return ((!operator <(value)) && (!operator ==(value)));
}

bool BQBLargeInteger::operator >(long value) const
{
    return operator >(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator >(unsigned long value) const
{
    return operator >(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator >(int value) const
{
    return operator >(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator >(unsigned int value) const
{
    return operator >(BQBLargeInteger(value));
}

//======================================================================
//  Unary operator <= for integral types.
//======================================================================

bool BQBLargeInteger::operator <=(const BQBLargeInteger & value) const
{
    return ! operator >(value);
}

bool BQBLargeInteger::operator <=(long value) const
{
    return operator <=(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator <=(unsigned long value) const
{
    return operator <=(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator <=(int value) const
{
    return operator <=(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator <=(unsigned int value) const
{
    return operator <=(BQBLargeInteger(value));
}

//======================================================================
//  Unary operator >= for integral types.
//======================================================================

bool BQBLargeInteger::operator >=(const BQBLargeInteger & value) const
{
    return ! operator <(value);
}

bool BQBLargeInteger::operator >=(long value) const
{
    return operator <(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator >=(unsigned long value) const
{
    return operator <(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator >=(int value) const
{
    return operator <(BQBLargeInteger(value));
}

bool BQBLargeInteger::operator >=(unsigned int value) const
{
    return operator <(BQBLargeInteger(value));
}

//======================================================================
//  Cast conversion operators.
//======================================================================

BQBLargeInteger::operator long() const
{
    return (long)(m_array_ptr[0]);
}

BQBLargeInteger::operator unsigned long() const
{
    return (unsigned long)(m_array_ptr[0]);
}

BQBLargeInteger::operator int() const
{
    return (int)(m_array_ptr[0]);
}

BQBLargeInteger::operator unsigned int() const
{
    return (m_array_ptr[0]);
}

BQBLargeInteger::operator short() const
{
    return (short)(m_array_ptr[0]);
}

BQBLargeInteger::operator unsigned short() const
{
    return (unsigned short)(m_array_ptr[0]);
}

BQBLargeInteger::operator char() const
{
    return (char)(m_array_ptr[0]);
}

BQBLargeInteger::operator unsigned char() const
{
    return (unsigned char)(m_array_ptr[0]);
}

BQBLargeInteger::operator bool() const
{
    unsigned int temp = 0;

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        temp = temp | m_array_ptr[i];
    }

    return temp != 0;
}

BQBLargeInteger::operator float() const
{
#if defined __GNUC__
    double Temp = double(*this);
    return float(Temp);
#else
    long double Temp = long double(*this);
    return float(Temp);
#endif
}

BQBLargeInteger::operator double() const
{
    double temp = 0;

    for (int i = m_integer_length - 1; i >= 0; i--)
    {
        temp = D_MAX_UNSIGNED_VALUE * temp;
        temp = temp + m_array_ptr[i];
    }

    return temp;
/*
#if defined __GNUC__
    double Temp = double(*this);
    return Temp;
#else
    long double Temp = long double(*this);
    return double(Temp);
#endif
*/
}

BQBLargeInteger::operator long double() const
{
    //------------------------------------------------------------------
    //  If this value is negative then calculate the absolute value.
    //------------------------------------------------------------------

    BQBLargeInteger x_temp;

    //------------------------------------------------------------------
    //  Calculate the value as a long double.
    //------------------------------------------------------------------

    long double temp = 0;

//    for (unsigned int i = m_integer_length - 1; (int)(i) >= 0; ++i)
    for (int i = m_integer_length - 1; i >= 0; i--)
    {
        temp = D_MAX_UNSIGNED_VALUE * temp;
        temp = temp + m_array_ptr[i];
    }

    return temp;
}

#ifdef _DEBUG

//======================================================================
//  Dump the internal format of this large integer.
//======================================================================

void BQBLargeInteger::DebugDump(const char * pszText)
{
    if (pszText != 0)
    {
        std::cout << pszText << std::endl;
    }

    std::cout << "Length = " << m_integer_length << std::endl;
    
    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        std::cout << i << "    " << m_array_ptr[i] << std::endl;
    }

    return;
}

#endif

//======================================================================
//  Function to copy an instance of a large integer.
//======================================================================

void BQBLargeInteger::Copy(const BQBLargeInteger & that)
{
    SetIntegerLength(that.m_integer_length);

    m_negative_flag = that.m_negative_flag;

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        m_array_ptr[i] = that.m_array_ptr[i];
    }

    return;
}

//======================================================================
//  Add to this large integer.
//======================================================================

void BQBLargeInteger::AddPositiveArray(const BQBLargeInteger & addend)
{
    //------------------------------------------------------------------
    //  Make sure this large integer's array is as least as large
    //  as the addends array.
    //------------------------------------------------------------------

    if (m_integer_length < addend.m_integer_length)
    {
        SetIntegerLength(addend.m_integer_length, true);
    }

    //------------------------------------------------------------------
    //  Add the addend to this large integer.
    //------------------------------------------------------------------

    unsigned int carry = 0;

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        unsigned int addend_array_value = addend.GetSafeArrayValue(i);
        m_array_ptr[i] = m_array_ptr[i] + addend_array_value + carry;
        carry = (m_array_ptr[i] < addend_array_value) ? 1 : 0;
    }

    //------------------------------------------------------------------
    //  If a carry occurred in the most significant word then increase
    //  the size of this large integer.
    //------------------------------------------------------------------

    if (carry == 1)
    {
        SetIntegerLength(m_integer_length + 1, true);
        m_array_ptr[m_integer_length - 1] = 1;
    }

    return;
}

//======================================================================
//  Subtract from this large integer.
//======================================================================

bool BQBLargeInteger::SubtractPositiveArray(const BQBLargeInteger & subtrahend)
{
    //------------------------------------------------------------------
    //  Make sure this large integer's array is the same size
    //  as the subtrahend's array.
    //------------------------------------------------------------------

    if (m_integer_length < subtrahend.m_integer_length)
    {
        SetIntegerLength(subtrahend.m_integer_length, true);
    }

    //------------------------------------------------------------------
    //  Save the most significant array value. If this increases
    //  after the subtraction then the minuend is less than the
    //  subtrahend.
    //------------------------------------------------------------------

    unsigned int most_significant_array_value = m_array_ptr[m_integer_length - 1];

    //------------------------------------------------------------------
    //  Subtract the subtrahend from this large integer.
    //------------------------------------------------------------------

    bool borrow = false;

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        unsigned int minuend_array_value = m_array_ptr[i];
        unsigned int subtrahend_array_value = subtrahend.GetSafeArrayValue(i);
        unsigned int difference = minuend_array_value - subtrahend_array_value;
        m_array_ptr[i] = difference - (borrow ? 1 : 0);

        borrow = borrow & (difference == 0);

        if (!borrow)
        {
            borrow = minuend_array_value < subtrahend_array_value;
        }
    }

    //------------------------------------------------------------------
    //  If the result taken as a two's complement integer is negative
    //  then the magnitude of the minuend was less than the magnitude
    //  of the subtrahend. In this case, take the two's complement of
    //  the result to get a positive result and return the value true
    //  from this method to signal the the sign of the result should
    //  be toggled.
    //------------------------------------------------------------------

    bool reverse_sign =
        m_array_ptr[m_integer_length - 1] > most_significant_array_value;

    if (reverse_sign)
    {
        //-------------------------------------------------------------
        //  Add one to the one's complement of the result to obtain
        //  the two's complement value.
        //-------------------------------------------------------------

        operator ~();
        BQBLargeInteger xOne(1);
        AddPositiveArray(xOne);
    }

    return reverse_sign;
}

//======================================================================
//  Function to shift this large integer to the left.
//======================================================================

void BQBLargeInteger::ShiftLeft(unsigned int shift_count)
{
    if (shift_count != 0)
    {
        //--------------------------------------------------------------
        //  The array could grow as a result of this shift. Determine
        //  the new size of the array. Find the first non-zero value
        //  in the array starting at the most significant array value.
        //--------------------------------------------------------------

        bool bNonZero = false;

        unsigned int nonzero_array_index = 0;

        for (nonzero_array_index = m_integer_length - 1;
             (int)(nonzero_array_index) >= 0;
             --nonzero_array_index)
        {
            if (m_array_ptr[nonzero_array_index] != 0)
            {
                bNonZero = true;
                break;
            }
        }

        //--------------------------------------------------------------
        //  If the value is zero then it is not necessary to shift
        //  this value.
        //--------------------------------------------------------------

        if (bNonZero)
        {
            //----------------------------------------------------------
            //  Determine the position of the bit in the most
            //  significant word.
            //----------------------------------------------------------

            unsigned int most_significant_value = m_array_ptr[nonzero_array_index];
            unsigned int bit_position = 0;

            for (bit_position = 0;
                 bit_position < INTEGER_BIT_COUNT;
                 ++bit_position)
            {
                //------------------------------------------------------
                //  Test the most significant bit.
                //------------------------------------------------------

                if ((int)(most_significant_value) < 0)
                {
                    break;
                }

                most_significant_value <<= 1;
            }

            bit_position = INTEGER_BIT_COUNT - 1 - bit_position;

            //----------------------------------------------------------
            //  Calculate the position of the most significant bit
            //  in the integer array.
            //----------------------------------------------------------

            unsigned int array_bit_position =
                INTEGER_BIT_COUNT * nonzero_array_index + bit_position;

            //----------------------------------------------------------
            //  Calculate the position of the most significant bit
            //  after shifting this number to the left.
            //----------------------------------------------------------

            unsigned int final_shift_position = array_bit_position + shift_count;

            //----------------------------------------------------------
            //  Make the number array large enough for the left shift.
            //----------------------------------------------------------

            unsigned int integer_length =
                (final_shift_position + INTEGER_BIT_COUNT) / INTEGER_BIT_COUNT;

            SetIntegerLength(integer_length, true);

            //----------------------------------------------------------
            //  Do all shifts that are a multiple of 32.
            //----------------------------------------------------------

            unsigned int array_shift_value = shift_count / INTEGER_BIT_COUNT;

            if (array_shift_value != 0)
            {
                for (unsigned int i = m_integer_length - 1;
                     (int)(i) >= 0;
                     --i)
                {
                    m_array_ptr[i] = GetSafeArrayValue(i - array_shift_value);
                }
            }

            //----------------------------------------------------------
            //  Do the remaining shifts.
            //----------------------------------------------------------

            unsigned int remaining_shift_count = shift_count
                - (array_shift_value * INTEGER_BIT_COUNT);
    
            for (unsigned int i = m_integer_length - 1;
                  (int)(i) >= 0;
                  --i)
            {
                m_array_ptr[i] = (m_array_ptr[i] << remaining_shift_count)
                    | (GetSafeArrayValue(i - 1)
                        >> (INTEGER_BIT_COUNT - remaining_shift_count));
            }

            //----------------------------------------------------------
            //  Normalize the result.
            //----------------------------------------------------------

            Normalize();
        }
    }
    
    return;
}

//======================================================================
//  Function to shift this large integer to the right.
//======================================================================

void BQBLargeInteger::ShiftRight(unsigned int shift_count)
{
    if (shift_count != 0)
    {
        //--------------------------------------------------------------
        //  Do all shifts that are a multiple of 32.
        //--------------------------------------------------------------

        if (shift_count < m_integer_length * INTEGER_BIT_COUNT)
        {
            unsigned int array_shift_value = shift_count / INTEGER_BIT_COUNT;

            if (array_shift_value != 0)
            {
                for (unsigned int i = 0; i < m_integer_length - array_shift_value; ++i)
                {
                    m_array_ptr[i] = m_array_ptr[i + array_shift_value];
                }
            }

            //----------------------------------------------------------
            //  Do the remaining shifts.
            //----------------------------------------------------------

            unsigned int remaining_shift_count = shift_count
                - (array_shift_value * INTEGER_BIT_COUNT);

            for (unsigned int i = 0; i < m_integer_length; ++i)
            {
                m_array_ptr[i] = (m_array_ptr[i] >> remaining_shift_count)
                    | (GetSafeArrayValue(i + 1)
                        << (INTEGER_BIT_COUNT - remaining_shift_count));
            }

            //----------------------------------------------------------
            //  Normalize the result.
            //----------------------------------------------------------

            Normalize();
        }
        else
        {
            SetToZero();
        }
    }

    return;
}

//======================================================================
//  Function to test for this large integer equal to zero.
//======================================================================

bool BQBLargeInteger::IsZero() const
{
    unsigned int temp = 0;

    for (unsigned int i = 0; i < m_integer_length; ++i)
    {
        temp |= m_array_ptr[i];
    }

    return temp == 0;
}

//======================================================================
//  Function to test for this large integer equal to zero.
//======================================================================

bool BQBLargeInteger::IsNegative() const
{
    return m_negative_flag;
}

//======================================================================
//  Function to test if the value fits in 32 bits.
//======================================================================

bool BQBLargeInteger::FitsIn32Bits() const
{
    return m_integer_length == 1;
}

//======================================================================
//  Member Function: BQBLargeInteger::SetIntegerLength
//  Author: Bill Hallahan
//  Date: March 30, 1998
//
//  Abstract:
//
//    This function is called to set the large integer length.
//    If the new integer length is greater than the current array
//    length OR the new buffer length is less than one fourth the
//    current array length then the array is reallocated. If the
//    buffer is resized to a be smaller array OR if the passed
//    boolean value 'copy_data_flag' is the value true then the
//    data in the old array is copied to the newly allocated
//    array before the old array memory is freed.
//
//
//  Input:
//
//    integer_length      The new integer length.
//
//    copy_data_flag         If this flag is true then the
//                           data in the array before the length
//                           is modified is copied into the new
//                           buffer. The value is sign extended
//                           if necessary..
//
//
//  Output:
//
//    This function has no return value.
//
//======================================================================

void BQBLargeInteger::SetIntegerLength(unsigned int integer_length,
                                    bool copy_data_flag)
{
    //------------------------------------------------------------------
    //  If the new data length is greater than the current buffer
    //  length then allocate a larger buffer for the data.
    //  Also if the new data length is less than one fourth the
    //  current buffer length then allocate a new array.
    //------------------------------------------------------------------

    bool reallocate_array_memory = integer_length > m_array_length;

    //------------------------------------------------------------------
    //  Also if the new data length is less than one fourth the
    //  current buffer length then allocate a new array.
    //------------------------------------------------------------------

    if (!reallocate_array_memory)
    {
        reallocate_array_memory = integer_length < (m_array_length >> 2);
        copy_data_flag = true;
    }

    //------------------------------------------------------------------
    //  Conditionally reallocate the large integer array.
    //------------------------------------------------------------------

    if (reallocate_array_memory)
    {
        //--------------------------------------------------------------
        //  Allocate a new buffer.
        //--------------------------------------------------------------

        unsigned int * array_ptr = new unsigned int [integer_length];

        bool memory_allocated = array_ptr != 0;

        if (memory_allocated)
        {
            if (copy_data_flag)
            {
                //------------------------------------------------------
                //  Copy the data from the old buffer to the new one.
                //------------------------------------------------------

                unsigned int copy_length = integer_length;

                if (m_integer_length < integer_length)
                {
                    copy_length = m_integer_length;
                }

                for (unsigned int i = 0; i < copy_length; ++i)
                {
                    array_ptr[i] = m_array_ptr[i];
                }

                //------------------------------------------------------
                // If the copied data is shorter than the new array
                // length then zero-fill the rest of the array.
                //------------------------------------------------------

                if (copy_length < integer_length)
                {
                    for (unsigned int j = copy_length; j < integer_length; ++j)
                    {
                        array_ptr[j] = 0;
                    }
                }

            }
            else
            {
                //------------------------------------------------------
                //  Fill the array with zeros.
                //------------------------------------------------------

                for (unsigned int i = 0; i < integer_length; ++i)
                {
                    array_ptr[i] = 0;
                }
            }

            //----------------------------------------------------------
            //  Delete the old buffer.
            //----------------------------------------------------------

            delete [] m_array_ptr;

            //----------------------------------------------------------
            //  Save the new buffer pointer.
            //----------------------------------------------------------

            m_array_ptr = array_ptr;

            //----------------------------------------------------------
            //  Set the buffer length and the integer length.
            //----------------------------------------------------------

            m_array_length = integer_length;
        }
        else
        {
            //----------------------------------------------------------
            //  Memory allocation failed.  Throw an exception.
            //----------------------------------------------------------

            //throw BQBLargeIntegerException("Memory allocation failed");
			printf("BQBLargeInteger::SetIntegerLength(): Error: Memory allocation failed.\n");
			abort();
        }
    }

    //------------------------------------------------------------------
    //  Set the length of the integer.
    //------------------------------------------------------------------

    m_integer_length = integer_length;

    return;
}

//======================================================================
//  Set this large integer value to zero.
//======================================================================

void BQBLargeInteger::SetToZero()
{
    SetIntegerLength(1);
    m_array_ptr[0] = 0;
    return;
}

//======================================================================
//  Get the array value. If the index is out of range then return
//  the value zero.
//======================================================================

unsigned int BQBLargeInteger::GetSafeArrayValue(unsigned int index) const
{
    unsigned int array_value;

    if (index < m_integer_length)
    {
        array_value = m_array_ptr[index];
    }
    else
    {
        array_value = 0;
    }

    return array_value;
}

//======================================================================
//  Set the integer length to the smallest length for this value.
//======================================================================

void BQBLargeInteger::Normalize()
{
    //------------------------------------------------------------------
    //  Discard all leading zeros.
    //------------------------------------------------------------------

    unsigned int integer_length = 1;

    for (unsigned int i = m_integer_length - 1; i > 0; --i)
    {
        if (m_array_ptr[i] != 0)
        {
            integer_length = i + 1;
            break;
        }
    }

    m_integer_length = integer_length;

    return;
}

//======================================================================
//  Private Method to compare the positive arrays of two large integers.
//======================================================================

bool BQBLargeInteger::LessThanPositiveArrayCompare(const BQBLargeInteger & value_0,
                                                const BQBLargeInteger & value_1) const
{
    unsigned int largest_integer_length =
        value_0.m_integer_length > value_1.m_integer_length ?
        value_0.m_integer_length : value_1.m_integer_length;

    bool is_less_than_flag = false;

    for (unsigned int i = largest_integer_length - 1; (int)(i) >= 0; --i)
    {
        unsigned int array_value_0 = value_0.GetSafeArrayValue(i);
        unsigned int array_value_1 = value_1.GetSafeArrayValue(i);

        if (array_value_0 < array_value_1)
        {
            is_less_than_flag = true;
        }
        else if (array_value_0 > array_value_1)
        {
            break;
        }
    }

    return is_less_than_flag;
}

//======================================================================
//  Calculate the population count for a 32 bit integer.
//======================================================================

unsigned int BQBLargeInteger::PopulationCount32Bit(unsigned int value) const
{
#ifndef LARGE_INT_ALTERNATE_POP_COUNT

    unsigned int shift_value = 1;

    for (unsigned int i = 0; i < BIT_MASK_TABLE_LENGTH; ++i)
    {
        unsigned int bit_mask = BIT_MASK_ARRAY[i];
        value = ((value & bit_mask) >> shift_value) + (value & ~bit_mask);
        shift_value = shift_value << 1;
    }

    return value;

#else

    unsigned int temp =
        value - ((value >> 1) & MASK_3)
        - ((value >> 2) & MASK_1);
    temp = (temp + (temp >> 3)) & MASK_7;
    return temp % 63;

#endif
}

//======================================================================
//  Get this large integer value in an stl string in base 8.
//======================================================================

void BQBLargeInteger::GetBase8NumberString(std::string & number_string) const
{
    //--------------------------------------------------------------
    //  Get the base 8 digits.
    //--------------------------------------------------------------

    std::string base8_digits_string;
    GetBase8Digits(base8_digits_string);

    //--------------------------------------------------------------
    //  Convert the base 8 digits to ascii characters.
    //--------------------------------------------------------------

    number_string.erase();

    for (unsigned int i = 0; i < base8_digits_string.length(); ++i)
    {
        number_string += base8_digits_string[i] + '0';
    }

    return;
}

//======================================================================
//  Get the base 8 digits in an STL string.
//======================================================================

void BQBLargeInteger::GetBase8Digits(std::string & base8_digits_string) const
{
    base8_digits_string.erase();

    //--------------------------------------------------------------
    //  Calculate the base 8 digits. Convert 3 bits at a time
    //  to an octal digit starting at the most significant bits.
    //--------------------------------------------------------------

    unsigned int current_bit_position_plus_one =
        INTEGER_BIT_COUNT * m_integer_length;

    unsigned int extra_bits = current_bit_position_plus_one % 3;
    unsigned int index = m_integer_length - 1;

    while (current_bit_position_plus_one != 0)
    {
        char base8_digit;

        switch (extra_bits)
        {
        case 0:

            base8_digits_string += (char)((m_array_ptr[index] >> 29) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 26) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 23) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 20) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 17) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 14) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 11) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 8) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 5) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 2) & 7);
            break;

        case 1:

            base8_digit = (char)(((GetSafeArrayValue(index + 1) & 3) << 1)
                    | ((m_array_ptr[index] >> 31) & 1));
            base8_digits_string +=  base8_digit;

            base8_digits_string += (char)((m_array_ptr[index] >> 28) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 25) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 22) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 19) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 16) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 13) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 10) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 7) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 4) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 1) & 7);
            break;

        case 2:

            base8_digit = (char)(((GetSafeArrayValue(index + 1) & 1) << 2)
                    | ((m_array_ptr[index] >> 30) & 3));
            base8_digits_string +=  base8_digit;

            base8_digits_string += (char)((m_array_ptr[index] >> 27) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 24) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 21) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 18) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 15) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 12) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 9) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 6) & 7);
            base8_digits_string += (char)((m_array_ptr[index] >> 3) & 7);
            base8_digits_string += (char)(m_array_ptr[index] & 7);
            break;

        default:

            break;
        }

        --index;
        extra_bits = (extra_bits + 1) % 3;
        current_bit_position_plus_one -= INTEGER_BIT_COUNT;
    }

    //------------------------------------------------------------------
    //  Strip any extra sign extension characters.
    //------------------------------------------------------------------

    StripLeadingZeroDigits(base8_digits_string, 0);

    return;
}

//======================================================================
//  Get this large integer value in an stl string in base 10.
//======================================================================

void BQBLargeInteger::GetBase10NumberString(std::string & number_string) const
{
    //------------------------------------------------------------------
    //  Determine whether the product is a negative value.
    //------------------------------------------------------------------

    BQBLargeInteger x_temp = *this;

    //------------------------------------------------------------------
    //  Get the base 8 digits.
    //------------------------------------------------------------------

    std::string digits_string;
    x_temp.GetBase8Digits(digits_string);

    //------------------------------------------------------------------
    //  Convert the base 8 string to a base 10 string using the
    //  algorithm from "Semi-Numerical Methods", by Knuth.
    //  Double the K leading octal digits using decimal arithmetic
    //  and subtract them from the K + 1 leading digits using
    //  decimal arithmetic.
    //------------------------------------------------------------------

    char * digit_array_ptr = const_cast<char *>(digits_string.data());
    unsigned int digit_length = (unsigned int)(digits_string.length());

    char * subtrahend_ptr = new char [digit_length];

    if (subtrahend_ptr == 0)
    {
        //throw BQBLargeIntegerException("Memory allocation failed");
		printf("BQBLargeInteger::GetBase10NumberString(): Error: Memory allocation failed.\n");
		abort();
    }

    for (unsigned int k = 0; k < digit_length - 1; ++k)
    {
        //--------------------------------------------------------------
        //  Double the K leading octal digits using base 10 arithmetic
        //  and copy these digits into the K + 1'st location in the
        //  subtrahend array.
        //--------------------------------------------------------------

        unsigned int j = 0;

        for (j = k + 1; (int)(j) >= 0; --j)
        {
            subtrahend_ptr[j] = 0;
        }

        for (j = k; (int)(j) >= 0; --j)
        {
            char doubled_digit = digit_array_ptr[j] << 1;

            if (doubled_digit > 9)
            {
                subtrahend_ptr[j + 1] += doubled_digit - 10;
                subtrahend_ptr[j] += 1;
            }
            else
            {
                subtrahend_ptr[j + 1] += doubled_digit;
            }
        }

        //--------------------------------------------------------------
        //  Subtract the doubled digits from the original number
        //  using decimal arithmetic.
        //--------------------------------------------------------------

        for (unsigned int m = k + 1; (int)(m) >= 0; --m)
        {
            char difference = digit_array_ptr[m] - subtrahend_ptr[m];

            if (difference < 0)
            {
                digit_array_ptr[m] =  difference + 10;

                if ((int)(m - 1) >= 0)
                {
                    digit_array_ptr[m - 1] -= 1;
                }
            }
            else
            {
                digit_array_ptr[m] = difference;
            }
        }
    }

    delete [] subtrahend_ptr;

    //------------------------------------------------------------------
    //  Convert the digits to characters. First skip all leading zeros.
    //------------------------------------------------------------------

    unsigned int first_nonzero_position = 0;

    for (first_nonzero_position = 0;
          first_nonzero_position < digits_string.length();
          ++first_nonzero_position)
    {
        if (digits_string[first_nonzero_position] != 0)
        {
            break;
        }
    }

    number_string.erase();

    for (unsigned int j = first_nonzero_position;
          j < digits_string.length();
          ++j)
    {
        number_string += digits_string[j] + '0';
    }

    return;
}

//======================================================================
//  Get this large integer value in an stl string in base 16.
//======================================================================

void BQBLargeInteger::GetBase16NumberString(std::string & number_string) const
{
    number_string.erase();

    //------------------------------------------------------------------
    //  Convert the digits.
    //------------------------------------------------------------------

    for (unsigned int i = m_integer_length - 1; (int)(i) >= 0; --i)
    {
        unsigned int hex_digits = m_array_ptr[i];

        for (int iShift = 28; iShift >= 0; iShift -= 4)
        {
            number_string += DIGIT_ARRAY[((hex_digits >> iShift) & 0x0F)];
        }
    }

    //------------------------------------------------------------------
    //  Strip any extra sign extension characters.
    //------------------------------------------------------------------

    StripLeadingZeroDigits(number_string, '0');

    return;
}

//======================================================================
//  Strip any extra zero digits.
//======================================================================

void BQBLargeInteger::StripLeadingZeroDigits(std::string & number_string,
                                          char zero_digit) const
{
    //------------------------------------------------------------------
    //  Find the position past any extra zero digit characters.
    //------------------------------------------------------------------

    int start_position = 0;
    int digit_length = (int)number_string.length();

    for (int index = 0; index < digit_length - 1; ++index)
    {
        if (number_string[index] == zero_digit)
        {
            ++start_position;
        }
        else
        {
            break;
        }
    }

    //------------------------------------------------------------------
    //  If necessary then strip leading zero digit characters.
    //------------------------------------------------------------------

    if (start_position != 0)
    {
        if (start_position > 0)
        {
            number_string = number_string.substr(start_position);
        }
        else
        {
            number_string = number_string.substr(digit_length - 1);
        }
    }

    return;
}

//======================================================================
//  Set the value of this large integer value using a base 8 string.
//  Leading sign characters and space characters must be removed
//  before calling this method.
//======================================================================

bool BQBLargeInteger::SetValueWithBase8String(const char * base8_digits_ptr)
{
    //------------------------------------------------------------------
    //  Get the length of the input string.
    //------------------------------------------------------------------

    size_t digit_length = ::strlen(base8_digits_ptr);

    //------------------------------------------------------------------
    //  Convert the octal digit characters to octal values.
    //------------------------------------------------------------------

    std::string digits_string = base8_digits_ptr;

    char * digit_array_ptr = const_cast<char *>(digits_string.data());
    digit_length = (unsigned int)(digits_string.length());

    //------------------------------------------------------------------
    //  Convert the digit characters to digit values.
    //------------------------------------------------------------------

    bool success_flag = false;

    for (int i = 0; i < (int)(digit_length); ++i)
    {
        char digit = digit_array_ptr[i];

        success_flag = ((digit >= '0') && (digit < '8'));

        if (success_flag)
        {
            digit_array_ptr[i] = digit - '0';
        }
        else
        {
            break;
        }
    }

    if (success_flag)
    {
        //--------------------------------------------------------------
        //  Strip any extra sign extension characters.
        //--------------------------------------------------------------

        StripLeadingZeroDigits(digits_string, 0);

        //--------------------------------------------------------------
        //  Set the value of this large integer using the base 8
        //  digit values.
        //--------------------------------------------------------------

        SetValueWithBase8DigitValues(digits_string);
    }

    return success_flag;
}

//======================================================================
//  Set the value of this large integer using the base 8 digit values.
//======================================================================

void BQBLargeInteger::SetValueWithBase8DigitValues(const std::string & base8_digit_string)
{
    //------------------------------------------------------------------
    //  Set the length of the array.
    //------------------------------------------------------------------

    unsigned int number_of_digits = (unsigned int)base8_digit_string.length();
    unsigned int current_bit_position_plus_one = 3 * number_of_digits;
    unsigned int integer_length =
        (current_bit_position_plus_one + (INTEGER_BIT_COUNT - 1)) / INTEGER_BIT_COUNT;

    SetIntegerLength(integer_length);

    //------------------------------------------------------------------
    //  Add each digit to the array.
    //------------------------------------------------------------------

    current_bit_position_plus_one -= 3;
    unsigned int array_index = integer_length - 1;

    for (unsigned int octal_digit_index = 0;
          octal_digit_index < number_of_digits;
          ++octal_digit_index)
    {
        // The bit length is currently 32-bits.
#ifdef INT_BIT_LENGTH_IS_POWER_OF_TWO
        unsigned int shift_value =
            current_bit_position_plus_one & (INTEGER_BIT_COUNT - 1);
#else
        unsigned int shift_value =
            current_bit_position_plus_one % INTEGER_BIT_COUNT;
#endif
        current_bit_position_plus_one -= 3;

        unsigned int octal_digit = 0;

        switch (shift_value)
        {
        case 0:
        case 1:
        case 2:
            {
                octal_digit = (unsigned int)(base8_digit_string[octal_digit_index]);
                m_array_ptr[array_index] |= octal_digit << shift_value;
            }

            --array_index;
            break;

        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
        case 17:
        case 18:
        case 19:
        case 20:
        case 21:
        case 22:
        case 23:
        case 24:
        case 25:
        case 26:
        case 27:
        case 28:
        case 29:
            {
                octal_digit = (unsigned int)(base8_digit_string[octal_digit_index]);
                m_array_ptr[array_index] |= octal_digit << shift_value;
            }

            break;

        case 30:
            {
                octal_digit = (unsigned int)(base8_digit_string[octal_digit_index]);
                m_array_ptr[array_index] |= (octal_digit & 3) << 30;
                
                if (array_index + 1 < integer_length)
                {
                    m_array_ptr[array_index + 1] |= (octal_digit & 4) >> 2;
                }
            }

            break;

        case 31:
            {
                octal_digit = (unsigned int)(base8_digit_string[octal_digit_index]);
                m_array_ptr[array_index] |= (octal_digit & 1) << 31;
                
                if (array_index + 1 < integer_length)
                {
                    m_array_ptr[array_index + 1] |= ((octal_digit & 6) >> 1);
                }
            }

            break;

        default:
            
            break;
        }
    }

    return;
}

//======================================================================
//  Set the value of this large integer value using a base 10 string.
//  Leading sign characters and space characters must be removed
//  before calling this method.
//======================================================================

bool BQBLargeInteger::SetValueWithBase10String(const char * digits_ptr)
{
    //------------------------------------------------------------------
    //  Convert the decimal string to octal digit values.
    //------------------------------------------------------------------

    std::string base8_digit_string;

    bool success_flag = ConvertDecimalStringToOctalDigits(digits_ptr,
                                                          base8_digit_string);

    //------------------------------------------------------------------
    //  Set the value of this large integer using the base 8
    //  digit values.
    //------------------------------------------------------------------

    if (success_flag)
    {
        SetValueWithBase8DigitValues(base8_digit_string);
    }

    return success_flag;
}

//======================================================================
//  Convert a base 10 string of characters to base 8 values.
//  The passed string contains decimal character digits representing
//  a positive value in base 10. Leading sign characters and space
//  characters must be removed before calling this method.
//======================================================================

bool BQBLargeInteger::ConvertDecimalStringToOctalDigits(const char * decimal_digits_ptr,
                                                     std::string & digits_string) const
{
    unsigned int i;

    //------------------------------------------------------------------
    //  Get the length of the input string.
    //------------------------------------------------------------------

    size_t digit_length = ::strlen(decimal_digits_ptr);

    //--------------------------------------------------------------
    //  Convert the decimal digit characters to decimal values.
    //--------------------------------------------------------------

    digits_string = decimal_digits_ptr;

    //--------------------------------------------------------------
    //  Insert leading zeros to allow for number growth.
    //
    //  The algorithm used to calculate the number of digits
    //  of growth is not exact, but it provides a conservative
    //  value that guarantees that there will be enough extra
    //  digits for growth.
    //--------------------------------------------------------------

    unsigned int padding = (unsigned int)(((200 * digit_length) / 1000) + 1);

    char * padding_ptr = new char [padding];

    if (padding_ptr == 0)
    {
        //throw BQBLargeIntegerException("Memory allocation failed");
		printf("BQBLargeInteger::ConvertDecimalStringToOctalDigits(): Error: Memory allocation failed.\n");
		abort();
    }

    for (i = 0; i < padding; ++i)
    {
        padding_ptr[i] = 0;
    }

    digits_string.insert(0, padding_ptr, padding);

    delete [] padding_ptr;

    //--------------------------------------------------------------
    //  Get the new digit length and the buffer pointer.
    //--------------------------------------------------------------

    digit_length = (unsigned int)(digits_string.length());
    char * digit_array_ptr = const_cast<char *>(digits_string.data());

    //--------------------------------------------------------------
    //  Convert the digit characters to digit values.
    //--------------------------------------------------------------

    bool success_flag = false;

    for (i = padding; i < digit_length; ++i)
    {
        char digit = digit_array_ptr[i];

        success_flag = ::isdigit(digit) != 0;

        if (success_flag)
        {
            digit_array_ptr[i] = digit - '0';
        }
        else
        {
            break;
        }
    }

    //--------------------------------------------------------------
    //  Convert the base 10 digits to base 8 digits.
    //--------------------------------------------------------------

    if (success_flag)
    {
        //----------------------------------------------------------
        //  Convert the base 8 string to base 10 digits using the
        //  algorithm from "Semi-Numerical Methods", by Knuth.
        //  Double the K leading octal digits using octal arithmetic
        //  and add them from the K + 1 leading digits using octal
        //  arithmetic.
        //----------------------------------------------------------

        char * addend_ptr = new char [digit_length + 1];

        if (addend_ptr == 0)
        {
            //throw BQBLargeIntegerException("Memory allocation failed");
			printf("BQBLargeInteger::ConvertDecimalStringToOctalDigits(): Error: Memory allocation failed.\n");
			abort();
		}

        for (unsigned int k = 0; k < digit_length - 1; ++k)
        {
            //------------------------------------------------------
            //  Double the K leading decimal digits using base 8
            //  arithmetic and copy these digits into the K + 1'st
            //  location in the addend array.
            //------------------------------------------------------

            unsigned int j = 0;

            for (j = k + 1; (int)(j) >= 0; --j)
            {
                addend_ptr[j] = 0;
            }

            for (j = k; (int)(j) >= 0; --j)
            {
                char doubled_digit = digit_array_ptr[j] << 1;
                
                if (doubled_digit > 7)
                {
                    addend_ptr[j + 1] += doubled_digit - 8;
                    addend_ptr[j] += 1;
                }
                else
                {
                    addend_ptr[j + 1] += doubled_digit;
                }
            }

            //----------------------------------------------------------
            //  Add the doubled digits from the original number
            //  using octal arithmetic.
            //----------------------------------------------------------

            for (unsigned int m = k + 1; (int)(m) >= 0; --m)
            {
                char sum = digit_array_ptr[m] + addend_ptr[m];
                
                if (sum > 7)
                {
                    digit_array_ptr[m] = sum - 8;
                    
                    if ((int)(m - 1) >= 0)
                    {
                        digit_array_ptr[m - 1] += 1;
                    }
                }
                else
                {
                    digit_array_ptr[m] = sum;
                }
            }
        }

        delete [] addend_ptr;

        //--------------------------------------------------------------
        //  Remove all leading zeros.
        //--------------------------------------------------------------

        unsigned int start_position = 0;

        for (i = 0; i < digits_string.length(); ++i)
        {
            if ((digits_string[i] != 0) || (i == digits_string.length() - 1))
            {
                start_position = i;
                break;
            }
        }

        if (start_position != 0)
        {
            digits_string = digits_string.substr(start_position);
        }
    }

    return success_flag;
}

//======================================================================
//  Set the value of this large integer value using a base 16 string.
//  Leading sign characters and space characters must be removed
//  before calling this method.
//======================================================================

bool BQBLargeInteger::SetValueWithBase16String(const char * digits_ptr)
{
    //------------------------------------------------------------------
    //  Convert the hexadecimal string to hexadecimal digit values.
    //------------------------------------------------------------------

    bool success_flag = true;
    std::string hex_number_string = digits_ptr;

    unsigned int number_of_digits = (unsigned int)(hex_number_string.length());
    char * hex_digits_ptr = const_cast<char *>(hex_number_string.data());

    for (unsigned int i = 0; i < number_of_digits; ++i)
    {
        int digit = (int)(hex_digits_ptr[i]);

        success_flag = ::isxdigit(digit) != 0;

        if (success_flag)
        {
            if (::isdigit(digit))
            {
                hex_digits_ptr[i] = (char)(digit - '0');
            }
            else
            {
                hex_digits_ptr[i] = (char)(::toupper(digit)) + 10 - 'A';
            }
        }
        else
        {
            break;
        }
    }

    //------------------------------------------------------------------
    //  Set the value of this large integer using the base 8
    //  digit values.
    //------------------------------------------------------------------

    if (success_flag)
    {
        //--------------------------------------------------------------
        //  Strip any extra sign extension characters.
        //--------------------------------------------------------------

        StripLeadingZeroDigits(hex_number_string, 0);

        //--------------------------------------------------------------
        //  Set the value of this large integer using the base 16
        //  digit values.
        //--------------------------------------------------------------

        SetValueWithBase16DigitValues(hex_number_string);
    }

    return success_flag;
}

//======================================================================
//  Set the value of this large integer value using a base 16 string.
//  Leading space characters must be removed before calling this method.
//======================================================================

void BQBLargeInteger::SetValueWithBase16DigitValues(const std::string & base16_digit_values_string)
{
    //------------------------------------------------------------------
    //  Set the length of the array.
    //------------------------------------------------------------------

    unsigned int number_of_digits = (unsigned int)base16_digit_values_string.length();
    unsigned int integer_length = (number_of_digits + 7) >> 3;

    SetIntegerLength(integer_length);

    //------------------------------------------------------------------
    //  Add each digit to the array.
    //  Handle extra digits that do not fill a complete word.
    //------------------------------------------------------------------

    unsigned int array_index = integer_length - 1;
    unsigned int current_digit_index = 0;
    unsigned int extra_digits =
        number_of_digits - ((integer_length - 1) << 3);
    unsigned int hex_digit;

    switch (extra_digits)
    {
    case 7:
        {
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 24;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 20;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 16;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 12;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 8;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 4;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index--] |= hex_digit;
        }

        break;

    case 6:
        {
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 20;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 16;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 12;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 8;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 4;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index--] |= hex_digit;
        }

        break;

    case 5:
        {
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 16;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 12;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 8;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 4;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index--] |= hex_digit;
        }

        break;

    case 4:
        {
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 12;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 8;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 4;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index--] |= hex_digit;
        }

        break;

    case 3:
        {
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 8;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 4;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index--] |= hex_digit;
        }

        break;

    case 2:
        {
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index] |= hex_digit << 4;
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index--] |= hex_digit;
        }

        break;

    case 1:
        {
            hex_digit = base16_digit_values_string[current_digit_index++];
            m_array_ptr[array_index--] |= hex_digit;
        }

        break;

    case 0:
    default:

        break;
    }

    if (number_of_digits > 7)
    {
        //-------------------------------------------------------------
        //  The number of digits remaining is a multiple of 8.
        //  Copy these digits into the term array.
        //-------------------------------------------------------------

        for (unsigned int hex_digitIndex = current_digit_index;
              hex_digitIndex < number_of_digits;
              hex_digitIndex += 8)
        {
            unsigned int hex_digit_0 = base16_digit_values_string[hex_digitIndex];
            unsigned int hex_digit_1 = base16_digit_values_string[hex_digitIndex + 1];
            unsigned int hex_digit_2 = base16_digit_values_string[hex_digitIndex + 2];
            unsigned int hex_digit_3 = base16_digit_values_string[hex_digitIndex + 3];
            unsigned int hex_digit_4 = base16_digit_values_string[hex_digitIndex + 4];
            unsigned int hex_digit_5 = base16_digit_values_string[hex_digitIndex + 5];
            unsigned int hex_digit_6 = base16_digit_values_string[hex_digitIndex + 6];
            unsigned int hex_digit_7 = base16_digit_values_string[hex_digitIndex + 7];
                m_array_ptr[array_index--] =
                    (hex_digit_0 << 28) | (hex_digit_1 << 24) |
                    (hex_digit_2 << 20) | (hex_digit_3 << 16) |
                    (hex_digit_4 << 12) | (hex_digit_5 << 8) |
                    (hex_digit_6 << 4) | hex_digit_7;
        }
    }

    return;
}

//======================================================================
//    Global function declarations.
//======================================================================

//======================================================================
//  operator + for integral types.
//======================================================================

//======================================================================
//  Addition of two instances of this class.
//======================================================================

BQBLargeInteger operator +(const BQBLargeInteger & addend_a,
                        const BQBLargeInteger & addend_b)
{
    return BQBLargeInteger(addend_a) += addend_b;
}

//======================================================================
//  Addition of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator +(const BQBLargeInteger & x_addend,
                        int addend)
{
    return x_addend + BQBLargeInteger(addend);
}

BQBLargeInteger operator +(int addend,
                        const BQBLargeInteger & x_addend)
{
    return x_addend + addend;
}

//======================================================================
//  Addition of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator +(const BQBLargeInteger & x_addend,
                        unsigned int addend)
{
    return x_addend + BQBLargeInteger(addend);
}

BQBLargeInteger operator +(unsigned int addend,
                        const BQBLargeInteger & x_addend)
{
    return x_addend + addend;
}

//======================================================================
//  operator - for integral types.
//======================================================================

//======================================================================
//  Subtraction of two instances of this class.
//======================================================================

BQBLargeInteger operator -(const BQBLargeInteger & minuend,
                        const BQBLargeInteger & subtrahend)
{
    return BQBLargeInteger(minuend) -= subtrahend;
}

//======================================================================
//  Subtraction with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator -(const BQBLargeInteger & minuend,
                        int subtrahend)
{
    return minuend - BQBLargeInteger(subtrahend);
}

BQBLargeInteger operator -(int minuend,
                        const BQBLargeInteger & subtrahend)
{
    return BQBLargeInteger(minuend) - subtrahend;
}

//======================================================================
//  Subtraction with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator -(const BQBLargeInteger & minuend,
                        unsigned int subtrahend)
{
    return minuend - BQBLargeInteger(subtrahend);
}

BQBLargeInteger operator -(unsigned int minuend,
                        const BQBLargeInteger & subtrahend)
{
    return BQBLargeInteger(minuend) - subtrahend;
}

//======================================================================
//  operator * for integral types.
//======================================================================

//======================================================================
//  Multiplication of two instances of this class.
//======================================================================

BQBLargeInteger operator *(const BQBLargeInteger & multiplier_a,
                        const BQBLargeInteger & multiplier_b)
{
    return BQBLargeInteger(multiplier_a) *= multiplier_b;
}

//======================================================================
//  Multiplication of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator *(const BQBLargeInteger & x_multiplier,
                        int multiplier)
{
    return x_multiplier * BQBLargeInteger(multiplier);
}

BQBLargeInteger operator *(int multiplier,
                        const BQBLargeInteger & x_multiplier)
{
    return x_multiplier * multiplier;
}

//======================================================================
//  Multiplication of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator *(const BQBLargeInteger & x_multiplier,
                        unsigned int multiplier)
{
    return x_multiplier * BQBLargeInteger(multiplier);
}

BQBLargeInteger operator *(unsigned int multiplier,
                        const BQBLargeInteger & x_multiplier)
{
    return x_multiplier * multiplier;
}

//======================================================================
//  operator / for integral types.
//======================================================================

//======================================================================
//  Division of two instances of this class.
//======================================================================

BQBLargeInteger operator /(const BQBLargeInteger & x_dividend,
                        const BQBLargeInteger & x_denominator)
{
    return BQBLargeInteger(x_dividend) /= x_denominator;
}

//======================================================================
//  Division with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator /(const BQBLargeInteger & x_dividend,
                        int divisor)
{
    return x_dividend / BQBLargeInteger(divisor);
}

BQBLargeInteger operator /(int dividend,
                        const BQBLargeInteger & x_denominator)
{
    return BQBLargeInteger(dividend) / x_denominator;
}

//======================================================================
//  Division with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator /(const BQBLargeInteger & x_dividend,
                        unsigned int divisor)
{
    return x_dividend / BQBLargeInteger(divisor);
}

BQBLargeInteger operator /(unsigned int dividend,
                        const BQBLargeInteger & x_denominator)
{
    return BQBLargeInteger(dividend) / x_denominator;
}

//======================================================================
//  operator << for integral types.
//======================================================================

//======================================================================
//  Left shift of two instances of this class.
//======================================================================

BQBLargeInteger operator <<(const BQBLargeInteger & value,
                         const BQBLargeInteger & x_shift_count)
{
    return BQBLargeInteger(value) <<= x_shift_count;
}

//======================================================================
//  Left shift with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator <<(const BQBLargeInteger & value,
                         int shift_count)
{
    return value << BQBLargeInteger(shift_count);
}

BQBLargeInteger operator <<(int value,
                         const BQBLargeInteger & x_shift_count)
{
    return BQBLargeInteger(value) << x_shift_count;
}

//======================================================================
//  Left shift with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator <<(const BQBLargeInteger & value,
                         unsigned int shift_count)
{
    return value << BQBLargeInteger(shift_count);
}

BQBLargeInteger operator <<(unsigned int value,
                         const BQBLargeInteger & x_shift_count)
{
    return BQBLargeInteger(value) << x_shift_count;
}

//======================================================================
//  operator >> for integral types.
//======================================================================

//======================================================================
//  Right shift of two instances of this class.
//======================================================================

BQBLargeInteger operator >>(const BQBLargeInteger & value,
                         const BQBLargeInteger & x_shift_count)
{
    return BQBLargeInteger(value) >>= x_shift_count;
}

//======================================================================
//  Right shift with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator >>(const BQBLargeInteger & value,
                         int shift_count)
{
    return value >> BQBLargeInteger(shift_count);
}

BQBLargeInteger operator >>(int value,
                         const BQBLargeInteger & x_shift_count)
{
    return BQBLargeInteger(value) >> x_shift_count;
}

//======================================================================
//  Right shift with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator >>(const BQBLargeInteger & value,
                         unsigned int shift_count)
{
    return value >> BQBLargeInteger(shift_count);
}

BQBLargeInteger operator >>(unsigned int value,
                         const BQBLargeInteger & x_shift_count)
{
    return BQBLargeInteger(value) >> x_shift_count;
}

//======================================================================
//  operator % for integral types.
//======================================================================

//======================================================================
//  Modulus operator with two instances of this class.
//======================================================================

BQBLargeInteger operator %(const BQBLargeInteger & x_dividend,
                        const BQBLargeInteger & x_denominator)
{
    return BQBLargeInteger(x_dividend) %= x_denominator;
}

//======================================================================
//  Modulus operator with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator %(const BQBLargeInteger & x_dividend,
                        int divisor)
{
    return x_dividend % BQBLargeInteger(divisor);
}

BQBLargeInteger operator %(int dividend,
                        const BQBLargeInteger & x_denominator)
{
    return BQBLargeInteger(dividend) % x_denominator;
}

//======================================================================
//  Modulus operator with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator %(const BQBLargeInteger & x_dividend,
                        unsigned int divisor)
{
    return x_dividend % BQBLargeInteger(divisor);
}

BQBLargeInteger operator %(unsigned int dividend,
                        const BQBLargeInteger & x_denominator)
{
    return BQBLargeInteger(dividend) % x_denominator;
}

//======================================================================
//  operator ^ for integral types.
//======================================================================

//======================================================================
//  Bitwise Exclusive OR of two instances of this class.
//======================================================================

BQBLargeInteger operator ^(const BQBLargeInteger & value_a,
                        const BQBLargeInteger & value_b)
{
    return BQBLargeInteger(value_a) ^= value_b;
}

//======================================================================
//  Bitwise Exclusive OR of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator ^(const BQBLargeInteger & x_value,
                        int value)
{
    return x_value % BQBLargeInteger(value);
}

BQBLargeInteger operator ^(int value,
                        const BQBLargeInteger & x_value)
{
    return x_value % value;
}

//======================================================================
//  Bitwise Exclusive OR of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator ^(const BQBLargeInteger & x_value,
                        unsigned int value)
{
    return x_value % BQBLargeInteger(value);
}

BQBLargeInteger operator ^(unsigned int value,
                        const BQBLargeInteger & x_value)
{
    return x_value % value;
}

//======================================================================
//  operator & for integral types.
//======================================================================

//======================================================================
//  Bitwise AND of two instances of this class.
//======================================================================

BQBLargeInteger operator &(const BQBLargeInteger & value_a,
                        const BQBLargeInteger & value_b)
{
    return BQBLargeInteger(value_a) &= value_b;
}

//======================================================================
//  Bitwise AND of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator &(const BQBLargeInteger & x_value,
                        int value)
{
    return x_value & BQBLargeInteger(value);
}

BQBLargeInteger operator &(int value,
                        const BQBLargeInteger & x_value)
{
    return x_value & value;
}

//======================================================================
//  Bitwise AND of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator &(const BQBLargeInteger & x_value,
                        unsigned int value)
{
    return x_value & BQBLargeInteger(value);
}

BQBLargeInteger operator &(unsigned int value,
                        const BQBLargeInteger & x_value)
{
    return x_value & value;
}

//======================================================================
//  operator | for integral types.
//======================================================================

//======================================================================
//  Bitwise OR of two instances of this class.
//======================================================================

BQBLargeInteger operator |(const BQBLargeInteger & value_a,
                        const BQBLargeInteger & value_b)
{
    return BQBLargeInteger(value_a) |= value_b;
}

//======================================================================
//  Bitwise OR of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator |(const BQBLargeInteger & x_value,
                        int value)
{
    return x_value | BQBLargeInteger(value);
}

BQBLargeInteger operator |(int value,
                        const BQBLargeInteger & x_value)
{
    return x_value | value;
}

//======================================================================
//  Bitwise OR of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator |(const BQBLargeInteger & x_value,
                        unsigned int value)
{
    return x_value | BQBLargeInteger(value);
}

BQBLargeInteger operator |(unsigned int value,
                        const BQBLargeInteger & x_value)
{
    return x_value | value;
}


void BQBLargeInteger::SelfTest() {
    // This hex number contains 1024 'f' digits.
    const char * pszNumber = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
    BQBLargeInteger xMultiplier0;
    xMultiplier0.SetDefaultBase(16);
    xMultiplier0 = pszNumber;
    // Another way to set the value is:
    //xMultiplier0.SetValue(pszNumber, 16);
    BQBLargeInteger xMultiplier1 = xMultiplier0;

    // Time raising the 1024 digit hex number to the 11th power
    // by doing ten multiplications.  This is not an efficient way
    // to raise a number to a power, but this is to time the
    // multiplication code.
    //DWORD dwStartTime = GetTickCount();

    BQBLargeInteger xProduct;
    BQBLargeInteger xBitTest;

    int i = 0;
    for (i = 0; i < 10; ++i)
    {
        xProduct = xMultiplier0 * xMultiplier1;
    }

    //DWORD dwStopTime = GetTickCount();
    //DWORD dwElapsedTime = dwStopTime - dwStartTime;

    //std::cout << "Time = " << dwElapsedTime << " milliseconds. " << std::endl;

    std::cout << "Product = " << std::hex << xProduct << std::endl;

    //------------------------------------------------------------------
    //  Test the large unsigned integer.
    //------------------------------------------------------------------

    xMultiplier0 = 0xFFFFFFFF;
    xMultiplier1 = 0xFFFFFFFF;
    xProduct = xMultiplier0 * xMultiplier1;
    std::cout << "Product = " << std::hex << xProduct << std::endl;

    // Negative number tests
    xMultiplier0 = -2;
    xMultiplier1 = 2;
    xProduct = xMultiplier0 * xMultiplier1;
    std::cout << "Product = " << std::hex << xProduct << std::endl;

    BQBLargeInteger xDividend = BQBLargeInteger(-16);
    BQBLargeInteger xDivisor = BQBLargeInteger(4);
    BQBLargeInteger xQuotient = xDividend / xDivisor;

    std::cout << "Quotient = " << std::hex << xQuotient << std::endl;

    xQuotient <<= 1;
    std::cout << "Quotient <<= 1 = " << std::hex << xQuotient << std::endl;

    xQuotient <<= 1;
    std::cout << "Quotient <<= 1 = " << std::hex << xQuotient << std::endl;

    xQuotient <<= 1;
    std::cout << "Quotient <<= 1 = " << std::hex << xQuotient << std::endl;

    xQuotient <<= 2;
    std::cout << "Quotient <<= 2 = " << std::hex << xQuotient << std::endl;

    xQuotient >>= 1;
    std::cout << "Quotient >>= 1 = " << std::hex << xQuotient << std::endl;

    xQuotient >>= 1;
    std::cout << "Quotient >>= 1 = " << std::hex << xQuotient << std::endl;

    xQuotient >>= 1;
    std::cout << "Quotient >>= 1 = " << std::hex << xQuotient << std::endl;

    xQuotient >>= 1;
    std::cout << "Quotient >>= 1 = " << std::hex << xQuotient << std::endl;

    xQuotient >>= 1;
    std::cout << "Quotient >>= 1 = " << std::hex << xQuotient << std::endl;

    xQuotient >>= 1;
    std::cout << "Quotient >>= 1 " << std::hex << xQuotient << std::endl;

    // A trivial comparison test.
    BQBLargeInteger xA = BQBLargeInteger(6);
    BQBLargeInteger xB = BQBLargeInteger(10);

    if (xA <= xB)
    {
        std::cout << "xA <= xB " << std::dec << xA << " <= " << xB << std::endl;
    }

    // Shift tests
    BQBLargeInteger x_value = BQBLargeInteger(1);

    for (i = 0; i < 256; i++)
    {
        x_value <<= 1;
        std::cout << std::dec << x_value << std::endl;
        std::cout << std::hex << x_value << std::endl;
    }

    for (i = 0; i < 256; i++)
    {
        x_value >>= 1;
        std::cout << std::dec << x_value << std::endl;
        std::cout << std::hex << x_value << std::endl;
    }

    // Another comparison test.
    if (x_value == 0)
    {
        std::cout << x_value << " == 0 Failure" << std::endl;
    }
    else
    {
        std::cout << x_value << " != 0 Success" << std::endl;        
    }

    // More number tests.
    xBitTest.SetValue("FFFFFFFF", 16);
    unsigned int uiLeadingBitPosition = xBitTest.LeadingBitPosition();
    std::cout << "Leading Bit Position = " << uiLeadingBitPosition << std::endl;

    xBitTest.SetValue("FFFFFFFF", 16);
    std::cout << "Population Count = " << xBitTest.PopulationCount() << std::endl;

    xBitTest.SetValue("FFFFFFEF", 16);
    std::cout << "Population Count = " << xBitTest.PopulationCount() << std::endl;
    
    xBitTest.SetValue("FFFFFFCF", 16);
    std::cout << "Population Count = " << xBitTest.PopulationCount() << std::endl;
    
    xBitTest.SetValue("FFFFFF8F", 16);
    std::cout << "Population Count = " << xBitTest.PopulationCount() << std::endl;
    
    xBitTest.SetValue("0FFFFFFF", 16);
    std::cout << "Population Count = " << xBitTest.PopulationCount() << std::endl;
    
    xBitTest.SetValue("F0EFFFFF", 16);
    std::cout << "Population Count = " << xBitTest.PopulationCount() << std::endl;
    
    xBitTest.SetValue("1", 16);
    std::cout << "Population Count = " << xBitTest.PopulationCount() << std::endl;

    xBitTest.SetValue("49249249", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;

    xBitTest.SetValue("49249249", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;

    xBitTest.SetValue("24924924", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;

    xBitTest.SetValue("24924924", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;

    xBitTest.SetValue("1234567012", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;

    xBitTest.SetValue("77777777777777777777777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;

    xBitTest.SetValue("200", 8);
    xBitTest.SetValue("4294967296", 10);
    xBitTest.SetValue("777777777777777777777777777777", 8);
    xBitTest.SetValue("77", 8);
    xBitTest.SetValue("1234567890123456789", 10);

    xBitTest.SetValue("12345678901234567890123456789012345678901234567890123456789012345678901234567890987654321", 10);
    std::cout << "Decimal Value = " << std::dec << xBitTest << std::endl;
    xBitTest.SetValue("7", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("77", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("7777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("77777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("7777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("77777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("7777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("77777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("7777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("77777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("777777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("7777777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("77777777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("777777777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("7777777777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("77777777777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("777777777777777777777", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;
    xBitTest.SetValue("12345670123456701234567654321", 8);
    std::cout << "Octal Value = " << std::oct << xBitTest << std::endl;

    xBitTest.SetValue("F", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("FFFFFFFFFFFFFFFFFFFFF", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    xBitTest.SetValue("123456789ABCDEF0123456789abcdef0123456789abcdefedcba987654321", 16);
    std::cout << "Hex Value = " << std::hex << xBitTest << std::endl;
    std::cout << "Decimal Value = " << std::dec << xBitTest << std::endl;

    xBitTest.SetValue("1", 10);
    xBitTest = xBitTest << 65535;
    xBitTest--;
    std::cout << "((1 << 65536) - 1) in base 10 = " << std::dec << xBitTest << std::endl;
    std::cout << "((1 << 65536) - 1) in base 8 = " << std::oct << xBitTest << std::endl;
    std::cout << "((1 << 65536) - 1) in base 16 = " << std::hex << xBitTest << std::endl;

    BQBLargeInteger x_value_0;
    BQBLargeInteger x_value_1;

    x_value_0.SetValue("1", 10);
    x_value_1.SetValue("1", 10);
    BQBLargeInteger x_result = x_value_0 + x_value_1;
    std::cout << x_value_0 << " + " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("1", 10);
    x_value_1.SetValue("-2", 10);
    x_result = x_value_0 + x_value_1;
    std::cout << x_value_0 << " + " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("-2", 10);
    x_value_1.SetValue("1", 10);
    x_result = x_value_0 + x_value_1;
    std::cout << x_value_0 << " + " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("-2", 10);
    x_value_1.SetValue("-1", 10);
    x_result = x_value_0 + x_value_1;
    std::cout << x_value_0 << " + " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("-2", 10);
    x_value_1.SetValue("-1", 10);
    x_result = x_value_0 + x_value_1;
    std::cout << x_value_0 << " + " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("1", 10);
    x_value_1.SetValue("2", 10);
    x_result = x_value_0 - x_value_1;
    std::cout << x_value_0 << " - " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("1", 10);
    x_value_1.SetValue("-2", 10);
    x_result = x_value_0 - x_value_1;
    std::cout << x_value_0 << " - " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("-1", 10);
    x_value_1.SetValue("2", 10);
    x_result = x_value_0 - x_value_1;
    std::cout << x_value_0 << " - " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("-1", 10);
    x_value_1.SetValue("-2", 10);
    x_result = x_value_0 - x_value_1;
    std::cout << x_value_0 << " - " << x_value_1 << " = " << x_result << std::endl;

    x_value_0.SetValue("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", 10);
    x_result = x_value_0 + 1;
    std::cout << x_value_0 << " + 1 = " << x_result << std::endl;

    x_result = x_value_0 - 1;
    std::cout << x_value_0 << " - 1 = " << x_result << std::endl;

    x_value_1 = x_value_0 - x_result;
    std::cout << x_value_0 << " - " << x_result << " = " << x_value_1 << std::endl;

    x_value_0.SetValue("10", 10);
    x_value_1 = x_value_0;
    x_result = x_value_0 * x_value_1;
    std::cout << x_result << " = " << x_value_0 << " * " << x_value_1 << std::endl;

    x_value_0.SetValue("1000000000000000", 10);
    std::cout << "x_value_0 = " << x_value_0 << std::endl;

    x_value_1 = x_value_0;
    x_result = x_value_0 * x_value_1;
    std::cout << x_result << " = " << x_value_0 << " * " << x_value_1 << std::endl;

    x_value_0.SetValue("1000000", 10);
    std::cout << x_value_0 << " = " << x_value_0 << std::endl;

    x_value_0.SetValue("100000000000000000000000", 10);
    x_value_1 = x_value_0;
    x_result = x_value_0 * x_value_1;
    std::cout << x_result << " = " << x_value_0 << " * " << x_value_1 << std::endl;

    xDividend.SetValue("1048576", 10);
    xDivisor.SetValue("1024", 10);
    xQuotient = xDividend / xDivisor;
    std::cout << xDividend << " / " << xDivisor << " = " << xQuotient << std::endl;
}


void BQBLargeInteger::mprintf(unsigned int base) {

	std::string s;

	GetNumberString(s,base);

	::printf("%s",s.c_str());
}


void BQBLargeInteger::mfprintf(FILE *a, unsigned int base) {

	std::string s;

	GetNumberString(s,base);

	::fprintf(a,"%s",s.c_str());
}


const char* BQBLargeInteger::string() const {

	if (*this != 0)
		GetNumberString(g_sBQBLIString,10);
	else
		g_sBQBLIString = "0";
	return g_sBQBLIString.c_str();
}


const char* BQBLargeInteger::string2() const {

	if (*this != 0)
		GetNumberString(g_sBQBLIString2,10);
	else
		g_sBQBLIString2 = "0";
	return g_sBQBLIString2.c_str();
}


const char* BQBLargeInteger::string3() const {

	if (*this != 0)
		GetNumberString(g_sBQBLIString3,10);
	else
		g_sBQBLIString3 = "0";
	return g_sBQBLIString3.c_str();
}



