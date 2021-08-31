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

//======================================================================
//  Class Definition File: BQBLargeInteger.h
//  Author: Bill Hallahan
//  Date: March 11, 1998
//
//  Abstract:
//
//    This file contains the definition for class BQBLargeInteger.
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
//          How To Use This Class
//
//    Except for the the caveats given below, 'BQBLargeInteger' can be
//    used as if were a built-in data type.
//
//    There are certain unavoidable differences which are documented
//    here:
//
//    - This module should be compiled with overflow checking disabled.
//      This class grows the array that stores the binary number so
//      that overflows cannot occur.
//
//    - The compiler will not allow large integer types to be the
//      integer argument in a 'switch' statement. The large integer
//      must be cast to an built-in type.
//
//    - Instances of this class cannot be used to index simple arrays.
//      If an instance of this class must be used as an index for an
//      array then it should be cast to a built-in type.
//
//    - Shift values for the left and right shift operators are limited
//      to the maximum value that fits in an unsigned integer. If the
//      shift value is larger than this then the shifted large integer
//      is set to the value zero. Note that the internal format does
//      not use two's complement arithmetic so that shifted values
//      maintain the same sign as the original value.
//
//    - Literals with unlimited bit-length are not supported by any
//      compiler. To initialize a large integer use a constructor,
//      an equals operator, SetBinaryValue() method, or one of the
//      overloaded SetValue() methods.
//
//======================================================================


#ifndef BQB_LARGEINTEGER_H
#define BQB_LARGEINTEGER_H


//======================================================================
//  Include files.
//======================================================================


// This must always be the first include directive
#include "bqb_config.h"

#include <iostream>
#include <stdio.h>
//#include <exception>


extern std::string g_sBQBLIString;
extern std::string g_sBQBLIString2;
extern std::string g_sBQBLIString3;


//======================================================================
// BQBLargeIntegerException class for errors.
//======================================================================
/*
class BQBLargeIntegerException : public std::exception
{
public:

    BQBLargeIntegerException(const char * exception_text)
        : std::exception(exception_text)
    {
    }

    ~BQBLargeIntegerException() throw()
    {
    }
};*/

//======================================================================
//  Class definition for class BQBLargeInteger.
//======================================================================

class BQBLargeInteger
{
private:

    //==================================================================
    //  Data members.
    //==================================================================

    unsigned int m_array_length;
    unsigned int m_integer_length;
    unsigned int m_default_base;
    unsigned int * m_array_ptr;
    bool m_negative_flag;

public:

	static void SelfTest();

	void mprintf(unsigned int base);

	void mfprintf(FILE *a, unsigned int base);

    //==================================================================
    //  istream operator for input.
    //==================================================================

    friend std::istream & operator >>(std::istream & is,
                                      BQBLargeInteger & value);

    //==================================================================
    //  ostream operator for output.
    //==================================================================

    friend std::ostream & operator <<(std::ostream & os,
                                      const BQBLargeInteger & value);

    //==================================================================
    //  Constructors
    //==================================================================

    BQBLargeInteger();

    //==================================================================
    //  Copy constructor
    //==================================================================

    BQBLargeInteger(const BQBLargeInteger & value);

    //==================================================================
    //  Conversion constructors.
    //==================================================================

    explicit BQBLargeInteger(long value);

    explicit BQBLargeInteger(unsigned long value);

    explicit BQBLargeInteger(int value);

    explicit BQBLargeInteger(unsigned int value);

    explicit BQBLargeInteger(const char * psznumber_ptr);

    //==================================================================
    //  Special constructor to allow initializing using a binary array.
    //==================================================================

    BQBLargeInteger(const char * binary_data_ptr,
                 unsigned int length);

    //==================================================================
    //  Destructor
    //==================================================================

    virtual ~BQBLargeInteger();

    //======================================================================
    // Method: BQBLargeInteger::SetDefaultBase
    //
    // The default base is the value used for the constructor and operator
    // equals methods that take only a string pointer for a number.
    //======================================================================

    void SetDefaultBase(unsigned int default_base);

    //==================================================================
    //  This method allows getting this large integer's binary value.
    //  The binary value is a positive number in little endian format.
    //  The buffer that is passed must be large enough to contain the
    //  number. This method can be called passing a null pointer for
    //  the buffer to just obtain the buffer length in bytes.
    //==================================================================

    unsigned int GetBinaryValue(unsigned char * binary_data_ptr = 0) const;

    //==================================================================
    //  This method allows setting this large integer's binary value.
    //  The binary value is a positive number in little endian format.
    //  The buffer that is passed must be large enough to contain the
    //  number. This method can be called passing a null pointer for
    //  the buffer to just obtain the buffer length in bytes.
    //==================================================================

    void SetBinaryValue(const char * binary_data_ptr,
                        unsigned int length);

    //==================================================================
    //  Get this large integer value represented in an stl string.
    //==================================================================

    void GetNumberString(std::string & number_string,
                         unsigned int base) const;

    //==================================================================
    //  Set this large integer value using the number in a character
    //  string.
    //==================================================================

    bool SetValue(const char * number_ptr, unsigned int base);

    //==================================================================
    //  This method returns the position of the leading bit in this
    //  large integer.
    //==================================================================

    unsigned int LeadingBitPosition() const;

    //==================================================================
    //  Return the count of the number of bits set in this large
    //  integer.
    //==================================================================

    unsigned int PopulationCount() const;

    //==================================================================
    //  operator = for integral types.
    //==================================================================

    BQBLargeInteger operator =(const BQBLargeInteger & value);

    BQBLargeInteger operator =(long value);

    BQBLargeInteger operator =(unsigned long value);

    BQBLargeInteger operator =(int value);

    BQBLargeInteger operator =(unsigned int value);

    BQBLargeInteger operator =(const char * psznumber_ptr);

    //==================================================================
    //  Unary operator + for integral types.
    //==================================================================

    BQBLargeInteger operator +=(const BQBLargeInteger & addend);

    BQBLargeInteger operator +=(long addend);

    BQBLargeInteger operator +=(unsigned long addend);

    BQBLargeInteger operator +=(int addend);

    BQBLargeInteger operator +=(unsigned int addend);

    //==================================================================
    //  Unary operator - for integral types.
    //==================================================================

    BQBLargeInteger operator -=(const BQBLargeInteger & subtrahend);

    BQBLargeInteger operator -=(long subtrahend);

    BQBLargeInteger operator -=(unsigned long subtrahend);

    BQBLargeInteger operator -=(int subtrahend);

    BQBLargeInteger operator -=(unsigned int subtrahend);

    //==================================================================
    //  Unary operator * for integral types.
    //==================================================================

    BQBLargeInteger operator *=(const BQBLargeInteger & multiplier);

protected:

    void AccumulateWithCarry(BQBLargeInteger & product,
                             int index,
                             unsigned int value);

public:

    BQBLargeInteger operator *=(long multiplier);

    BQBLargeInteger operator *=(unsigned long multiplier);

    BQBLargeInteger operator *=(int multiplier);

    BQBLargeInteger operator *=(unsigned int multiplier);

    //==================================================================
    //  Unary operator / for integral types.
    //==================================================================

    BQBLargeInteger operator /=(const BQBLargeInteger & divisor);

    BQBLargeInteger operator /=(long divisor);

    BQBLargeInteger operator /=(unsigned long divisor);

    BQBLargeInteger operator /=(int divisor);

    BQBLargeInteger operator /=(unsigned int divisor);

    //==================================================================
    //  operator <<= for integral types.
    //==================================================================

    BQBLargeInteger operator <<=(const BQBLargeInteger & shift_count);

    BQBLargeInteger operator <<=(long shift_count);

    BQBLargeInteger operator <<=(unsigned long shift_count);

    BQBLargeInteger operator <<=(int shift_count);

    BQBLargeInteger operator <<=(unsigned int shift_count);

    //==================================================================
    //  operator >>= for integral types.
    //==================================================================

    BQBLargeInteger operator >>=(const BQBLargeInteger & shift_count);

    BQBLargeInteger operator >>=(long shift_count);

    BQBLargeInteger operator >>=(unsigned long shift_count);

    BQBLargeInteger operator >>=(int shift_count);

    BQBLargeInteger operator >>=(unsigned int shift_count);

    //==================================================================
    //  Unary operator %= for integral types.
    //==================================================================

    BQBLargeInteger operator %=(const BQBLargeInteger & divisor);

    BQBLargeInteger operator %=(long divisor);

    BQBLargeInteger operator %=(unsigned long divisor);

    BQBLargeInteger operator %=(int divisor);

    BQBLargeInteger operator %=(unsigned int divisor);

    //==================================================================
    //  operator ^= for integral types.
    //==================================================================

    BQBLargeInteger operator ^=(const BQBLargeInteger & value);

    BQBLargeInteger operator ^=(long value);

    BQBLargeInteger operator ^=(unsigned long value);

    BQBLargeInteger operator ^=(int value);

    BQBLargeInteger operator ^=(unsigned int value);

    //==================================================================
    //  operator &= for integral types.
    //==================================================================

    BQBLargeInteger operator &=(const BQBLargeInteger & value);

    BQBLargeInteger operator &=(unsigned long value);

    BQBLargeInteger operator &=(int value);

    BQBLargeInteger operator &=(unsigned int value);

    //==================================================================
    //  operator |= for integral types.
    //==================================================================

    BQBLargeInteger operator |=(const BQBLargeInteger & value);

    BQBLargeInteger operator |=(long value);

    BQBLargeInteger operator |=(unsigned long value);

    BQBLargeInteger operator |=(int value);

    BQBLargeInteger operator |=(unsigned int value);

    //==================================================================
    //  Unary operators. The unary * operator and the unary & operator
    //  do not need to be overloaded.
    //==================================================================

    BQBLargeInteger operator !();

    BQBLargeInteger operator ~();

    BQBLargeInteger operator +();

    BQBLargeInteger operator -();

    //==================================================================
    //  The prefix form of the increment and decrement operators.
    //==================================================================

    const BQBLargeInteger operator ++();

    const BQBLargeInteger operator --();

    //==================================================================
    //  The postfix form of the increment and decrement operators.
    //==================================================================

    const BQBLargeInteger operator ++(int);

    const BQBLargeInteger operator --(int);

    //==================================================================
    //  Unary operator == for integral types.
    //==================================================================

    bool operator ==(const BQBLargeInteger & value) const;

    bool operator ==(long value) const;

    bool operator ==(unsigned long value) const;

    bool operator ==(int value) const;

    bool operator ==(unsigned int value) const;

    //==================================================================
    //  Unary operator != for integral types.
    //==================================================================

    bool operator !=(const BQBLargeInteger & value) const;

    bool operator !=(long value) const;

    bool operator !=(unsigned long value) const;

    bool operator !=(int value) const;

    bool operator !=(unsigned int value) const;

    //==================================================================
    //  Unary operator < for integral types.
    //==================================================================

    bool operator <(const BQBLargeInteger & value) const;

    bool operator <(long value) const;

    bool operator <(unsigned long value) const;

    bool operator <(int value) const;

    bool operator <(unsigned int value) const;

    //==================================================================
    //  Unary operator > for integral types.
    //==================================================================

    bool operator >(const BQBLargeInteger & value) const;

    bool operator >(long value) const;

    bool operator >(unsigned long value) const;

    bool operator >(int value) const;

    bool operator >(unsigned int value) const;

    //==================================================================
    //  Unary operator <= for integral types.
    //==================================================================

    bool operator <=(const BQBLargeInteger & value) const;

    bool operator <=(long value) const;

    bool operator <=(unsigned long value) const;

    bool operator <=(int value) const;

    bool operator <=(unsigned int value) const;

    //==================================================================
    //  Unary operator >= for integral types.
    //==================================================================

    bool operator >=(const BQBLargeInteger & value) const;

    bool operator >=(long value) const;

    bool operator >=(unsigned long value) const;

    bool operator >=(int value) const;

    bool operator >=(unsigned int value) const;

    //==================================================================
    //  Cast conversion operators.
    //==================================================================

    operator long() const;

    operator unsigned long() const;

    operator int() const;

    operator unsigned int() const;

    operator short() const;

    operator unsigned short() const;

    operator char() const;

    operator unsigned char() const;

    operator bool() const;

    operator float() const;

    operator double() const;

    operator long double() const;
    
#ifdef _DEBUG

    void DebugDump(const char * pszText = 0);

#endif

	const char* string() const;
	const char* string2() const;
	const char* string3() const;


private:

    //==================================================================
    //  The following operators are not declared. The default
    //  behavior of these operators is the desired behavior.
    //
    //  operator ()
    //  operator []
    //  operator new
    //  operator new[]
    //  operator delete
    //  operator delete[]
    //  operator ,
    //  operator ->*() const;
    //  operator &&()
    //  operator ||()
    //
    //  The following unary operators are not declared. The default
    //  behavior of these operators is the desired behavior.
    //
    //  operator *() 
    //  operator &()
    //
    //==================================================================

    //==================================================================
    //  Private member functions.
    //==================================================================

    void Copy(const BQBLargeInteger & that);

    void AddPositiveArray(const BQBLargeInteger & addend);

    bool SubtractPositiveArray(const BQBLargeInteger & subtrahend);

    void ShiftLeft(unsigned int shift_count);

    void ShiftRight(unsigned int shift_count);

    bool IsZero() const;

    bool IsNegative() const;

    bool FitsIn32Bits() const;

    void SetIntegerLength(unsigned int integer_length,
                          bool copy_data_flag = false);

    void SetToZero();

    unsigned int GetSafeArrayValue(unsigned int iIndex) const;

    void Normalize();

    bool LessThanPositiveArrayCompare(const BQBLargeInteger & value_0,
                                      const BQBLargeInteger & value_1) const;

    unsigned int PopulationCount32Bit(unsigned int value) const;

    void GetBase8NumberString(std::string & number_string) const;

    void GetBase8Digits(std::string & base8_digits_string) const;

    void GetBase10NumberString(std::string & number_string) const;

    void GetBase16NumberString(std::string & number_string) const;

    void StripLeadingZeroDigits(std::string & number_string,
                                char zero_digit) const;

    bool SetValueWithBase8String(const char * digits_ptr);

    void SetValueWithBase8DigitValues(const std::string & base8_digit_string);

    bool SetValueWithBase10String(const char * digits_ptr);

    bool ConvertDecimalStringToOctalDigits(const char * decimal_digits_ptr,
                                           std::string & digits_string) const;

    bool SetValueWithBase16String(const char * digits_ptr);

    void SetValueWithBase16DigitValues(const std::string & base16_digit_values_string);
};

//======================================================================
//  Global function declarations.
//======================================================================

//======================================================================
//  operator + for integral types.
//======================================================================

//======================================================================
//  Addition of two instances of this class.
//======================================================================

BQBLargeInteger operator +(const BQBLargeInteger & addend_a,
                        const BQBLargeInteger & addend_b);

//======================================================================
//  Addition of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator +(const BQBLargeInteger & x_addend,
                        int addend);

BQBLargeInteger operator +(int addend,
                        const BQBLargeInteger & x_addend);

//======================================================================
//  Addition of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator +(const BQBLargeInteger & x_addend,
                        unsigned int addend);

BQBLargeInteger operator +(unsigned int addend,
                        const BQBLargeInteger & x_addend);

//======================================================================
//  operator - for integral types.
//======================================================================

//======================================================================
//  Subtraction of two instances of this class.
//======================================================================

BQBLargeInteger operator -(const BQBLargeInteger & minuend,
                        const BQBLargeInteger & subtrahend);

//======================================================================
//  Subtraction with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator -(const BQBLargeInteger & minuend,
                        int subtrahend);

BQBLargeInteger operator -(int minuend,
                        const BQBLargeInteger & subtrahend);

//======================================================================
//  Subtraction with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator -(const BQBLargeInteger & minuend,
                        unsigned int subtrahend);

BQBLargeInteger operator -(unsigned int minuend,
                        const BQBLargeInteger & subtrahend);

//======================================================================
//  operator * for integral types.
//======================================================================

//======================================================================
//  Multiplication of two instances of this class.
//======================================================================

BQBLargeInteger operator *(const BQBLargeInteger & multiplier_a,
                        const BQBLargeInteger & multiplier_b);

//======================================================================
//  Multiplication of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator *(const BQBLargeInteger & x_multiplier,
                        int multiplier);

BQBLargeInteger operator *(int multiplier,
                        const BQBLargeInteger & x_multiplier);

//======================================================================
//  Multiplication of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator *(const BQBLargeInteger & x_multiplier,
                        unsigned int multiplier);

BQBLargeInteger operator *(unsigned int multiplier,
                        const BQBLargeInteger & x_multiplier);

//======================================================================
//  operator / for integral types.
//======================================================================

//======================================================================
//  Division of two instances of this class.
//======================================================================

BQBLargeInteger operator /(const BQBLargeInteger & x_dividend,
                        const BQBLargeInteger & x_denominator);

//======================================================================
//  Division with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator /(const BQBLargeInteger & x_dividend,
                        int divisor);

BQBLargeInteger operator /(int dividend,
                        const BQBLargeInteger & x_denominator);

//======================================================================
//  Division with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator /(const BQBLargeInteger & x_dividend,
                        unsigned int divisor);

BQBLargeInteger operator /(unsigned int dividend,
                        const BQBLargeInteger & x_denominator);

//======================================================================
//  operator << for integral types.
//======================================================================

//======================================================================
//  Left shift of two instances of this class.
//======================================================================

BQBLargeInteger operator <<(const BQBLargeInteger & value,
                         const BQBLargeInteger & x_shift_count);

//======================================================================
//  Left shift with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator <<(const BQBLargeInteger & value,
                         int shift_count);

BQBLargeInteger operator <<(int value,
                         const BQBLargeInteger & x_shift_count);

//======================================================================
//  Left shift with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator <<(const BQBLargeInteger & value,
                         unsigned int shift_count);

BQBLargeInteger operator <<(unsigned int value,
                         const BQBLargeInteger & x_shift_count);

//======================================================================
//  operator >> for integral types.
//======================================================================

//======================================================================
//  Right shift of two instances of this class.
//======================================================================

BQBLargeInteger operator >>(const BQBLargeInteger & value,
                         const BQBLargeInteger & x_shift_count);

//======================================================================
//  Right shift with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator >>(const BQBLargeInteger & value,
                         int shift_count);

BQBLargeInteger operator >>(int value,
                         const BQBLargeInteger & x_shift_count);

//======================================================================
//  Right shift with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator >>(const BQBLargeInteger & value,
                         unsigned int shift_count);

BQBLargeInteger operator >>(unsigned int value,
                         const BQBLargeInteger & x_shift_count);

//======================================================================
//  operator % for integral types.
//======================================================================

//======================================================================
//  Modulus operator with two instances of this class.
//======================================================================

BQBLargeInteger operator %(const BQBLargeInteger & x_dividend,
                        const BQBLargeInteger & x_denominator);

//======================================================================
//  Modulus operator with an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator %(const BQBLargeInteger & x_dividend,
                        int divisor);

BQBLargeInteger operator %(int dividend,
                        const BQBLargeInteger & x_denominator);

//======================================================================
//  Modulus operator with an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator %(const BQBLargeInteger & x_dividend,
                        unsigned int divisor);

BQBLargeInteger operator %(unsigned int dividend,
                        const BQBLargeInteger & x_denominator);

//======================================================================
//  operator ^ for integral types.
//======================================================================

//======================================================================
//  Bitwise Exclusive OR of two instances of this class.
//======================================================================

BQBLargeInteger operator ^(const BQBLargeInteger & value_a,
                        const BQBLargeInteger & value_b);

//======================================================================
//  Bitwise Exclusive OR of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator ^(const BQBLargeInteger & x_value,
                        int value);

BQBLargeInteger operator ^(int value,
                        const BQBLargeInteger & x_value);

//======================================================================
//  Bitwise Exclusive OR of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator ^(const BQBLargeInteger & x_value,
                        unsigned int value);

BQBLargeInteger operator ^(unsigned int value,
                        const BQBLargeInteger & x_value);

//======================================================================
//  operator & for integral types.
//======================================================================

//======================================================================
//  Bitwise AND of two instances of this class.
//======================================================================

BQBLargeInteger operator &(const BQBLargeInteger & value_a,
                        const BQBLargeInteger & value_b);

//======================================================================
//  Bitwise AND of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator &(const BQBLargeInteger & x_value,
                        int value);

BQBLargeInteger operator &(int value,
                        const BQBLargeInteger & x_value);

//======================================================================
//  Bitwise AND of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator &(const BQBLargeInteger & x_value,
                        unsigned int value);

BQBLargeInteger operator &(unsigned int value,
                        const BQBLargeInteger & x_value);

//======================================================================
//  operator | for integral types.
//======================================================================

//======================================================================
//  Bitwise OR of two instances of this class.
//======================================================================

BQBLargeInteger operator |(const BQBLargeInteger & value_a,
                        const BQBLargeInteger & value_b);

//======================================================================
//  Bitwise OR of an instance of the BQBLargeInteger class and a signed integer.
//======================================================================

BQBLargeInteger operator |(const BQBLargeInteger & x_value,
                        int value);

BQBLargeInteger operator |(int value,
                        const BQBLargeInteger & x_value);

//======================================================================
//  Bitwise OR of an instance of the BQBLargeInteger class and an unsigned integer.
//======================================================================

BQBLargeInteger operator |(const BQBLargeInteger & x_value,
                        unsigned int value);

BQBLargeInteger operator |(unsigned int value,
                        const BQBLargeInteger & x_value);

#endif
