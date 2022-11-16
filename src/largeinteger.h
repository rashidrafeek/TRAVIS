/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

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
//  Class Definition File: LargeInteger.h
//  Author: Bill Hallahan
//  Date: March 11, 1998
//
//  Abstract:
//
//    This file contains the definition for class LargeInteger.
//    An instance of LargeInteger can be used to store and calculate
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
//    Except for the the caveats given below, 'LargeInteger' can be
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


#ifndef LARGEINTEGER_H
#define LARGEINTEGER_H


//======================================================================
//  Include files.
//======================================================================


// This must always be the first include directive
#include "config.h"

#include <iostream>
#include <stdio.h>
//#include <exception>


extern std::string g_sLIString;
extern std::string g_sLIString2;
extern std::string g_sLIString3;


//======================================================================
// LargeIntegerException class for errors.
//======================================================================
/*
class LargeIntegerException : public std::exception
{
public:

    LargeIntegerException(const char * exception_text)
        : std::exception(exception_text)
    {
    }

    ~LargeIntegerException() throw()
    {
    }
};*/

//======================================================================
//  Class definition for class LargeInteger.
//======================================================================

class LargeInteger
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
                                      LargeInteger & value);

    //==================================================================
    //  ostream operator for output.
    //==================================================================

    friend std::ostream & operator <<(std::ostream & os,
                                      const LargeInteger & value);

    //==================================================================
    //  Constructors
    //==================================================================

    LargeInteger();

    //==================================================================
    //  Copy constructor
    //==================================================================

    LargeInteger(const LargeInteger & value);

    //==================================================================
    //  Conversion constructors.
    //==================================================================

    explicit LargeInteger(long value);

    explicit LargeInteger(unsigned long value);

    explicit LargeInteger(int value);

    explicit LargeInteger(unsigned int value);

    explicit LargeInteger(const char * psznumber_ptr);

    //==================================================================
    //  Special constructor to allow initializing using a binary array.
    //==================================================================

    LargeInteger(const char * binary_data_ptr,
                 unsigned int length);

    //==================================================================
    //  Destructor
    //==================================================================

    virtual ~LargeInteger();

    //======================================================================
    // Method: LargeInteger::SetDefaultBase
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

    LargeInteger operator =(const LargeInteger & value);

    LargeInteger operator =(long value);

    LargeInteger operator =(unsigned long value);

    LargeInteger operator =(int value);

    LargeInteger operator =(unsigned int value);

    LargeInteger operator =(const char * psznumber_ptr);

    //==================================================================
    //  Unary operator + for integral types.
    //==================================================================

    LargeInteger operator +=(const LargeInteger & addend);

    LargeInteger operator +=(long addend);

    LargeInteger operator +=(unsigned long addend);

    LargeInteger operator +=(int addend);

    LargeInteger operator +=(unsigned int addend);

    //==================================================================
    //  Unary operator - for integral types.
    //==================================================================

    LargeInteger operator -=(const LargeInteger & subtrahend);

    LargeInteger operator -=(long subtrahend);

    LargeInteger operator -=(unsigned long subtrahend);

    LargeInteger operator -=(int subtrahend);

    LargeInteger operator -=(unsigned int subtrahend);

    //==================================================================
    //  Unary operator * for integral types.
    //==================================================================

    LargeInteger operator *=(const LargeInteger & multiplier);

protected:

    void AccumulateWithCarry(LargeInteger & product,
                             int index,
                             unsigned int value);

public:

    LargeInteger operator *=(long multiplier);

    LargeInteger operator *=(unsigned long multiplier);

    LargeInteger operator *=(int multiplier);

    LargeInteger operator *=(unsigned int multiplier);

    //==================================================================
    //  Unary operator / for integral types.
    //==================================================================

    LargeInteger operator /=(const LargeInteger & divisor);

    LargeInteger operator /=(long divisor);

    LargeInteger operator /=(unsigned long divisor);

    LargeInteger operator /=(int divisor);

    LargeInteger operator /=(unsigned int divisor);

    //==================================================================
    //  operator <<= for integral types.
    //==================================================================

    LargeInteger operator <<=(const LargeInteger & shift_count);

    LargeInteger operator <<=(long shift_count);

    LargeInteger operator <<=(unsigned long shift_count);

    LargeInteger operator <<=(int shift_count);

    LargeInteger operator <<=(unsigned int shift_count);

    //==================================================================
    //  operator >>= for integral types.
    //==================================================================

    LargeInteger operator >>=(const LargeInteger & shift_count);

    LargeInteger operator >>=(long shift_count);

    LargeInteger operator >>=(unsigned long shift_count);

    LargeInteger operator >>=(int shift_count);

    LargeInteger operator >>=(unsigned int shift_count);

    //==================================================================
    //  Unary operator %= for integral types.
    //==================================================================

    LargeInteger operator %=(const LargeInteger & divisor);

    LargeInteger operator %=(long divisor);

    LargeInteger operator %=(unsigned long divisor);

    LargeInteger operator %=(int divisor);

    LargeInteger operator %=(unsigned int divisor);

    //==================================================================
    //  operator ^= for integral types.
    //==================================================================

    LargeInteger operator ^=(const LargeInteger & value);

    LargeInteger operator ^=(long value);

    LargeInteger operator ^=(unsigned long value);

    LargeInteger operator ^=(int value);

    LargeInteger operator ^=(unsigned int value);

    //==================================================================
    //  operator &= for integral types.
    //==================================================================

    LargeInteger operator &=(const LargeInteger & value);

    LargeInteger operator &=(unsigned long value);

    LargeInteger operator &=(int value);

    LargeInteger operator &=(unsigned int value);

    //==================================================================
    //  operator |= for integral types.
    //==================================================================

    LargeInteger operator |=(const LargeInteger & value);

    LargeInteger operator |=(long value);

    LargeInteger operator |=(unsigned long value);

    LargeInteger operator |=(int value);

    LargeInteger operator |=(unsigned int value);

    //==================================================================
    //  Unary operators. The unary * operator and the unary & operator
    //  do not need to be overloaded.
    //==================================================================

    LargeInteger operator !();

    LargeInteger operator ~();

    LargeInteger operator +();

    LargeInteger operator -();

    //==================================================================
    //  The prefix form of the increment and decrement operators.
    //==================================================================

    const LargeInteger operator ++();

    const LargeInteger operator --();

    //==================================================================
    //  The postfix form of the increment and decrement operators.
    //==================================================================

    const LargeInteger operator ++(int);

    const LargeInteger operator --(int);

    //==================================================================
    //  Unary operator == for integral types.
    //==================================================================

    bool operator ==(const LargeInteger & value) const;

    bool operator ==(long value) const;

    bool operator ==(unsigned long value) const;

    bool operator ==(int value) const;

    bool operator ==(unsigned int value) const;

    //==================================================================
    //  Unary operator != for integral types.
    //==================================================================

    bool operator !=(const LargeInteger & value) const;

    bool operator !=(long value) const;

    bool operator !=(unsigned long value) const;

    bool operator !=(int value) const;

    bool operator !=(unsigned int value) const;

    //==================================================================
    //  Unary operator < for integral types.
    //==================================================================

    bool operator <(const LargeInteger & value) const;

    bool operator <(long value) const;

    bool operator <(unsigned long value) const;

    bool operator <(int value) const;

    bool operator <(unsigned int value) const;

    //==================================================================
    //  Unary operator > for integral types.
    //==================================================================

    bool operator >(const LargeInteger & value) const;

    bool operator >(long value) const;

    bool operator >(unsigned long value) const;

    bool operator >(int value) const;

    bool operator >(unsigned int value) const;

    //==================================================================
    //  Unary operator <= for integral types.
    //==================================================================

    bool operator <=(const LargeInteger & value) const;

    bool operator <=(long value) const;

    bool operator <=(unsigned long value) const;

    bool operator <=(int value) const;

    bool operator <=(unsigned int value) const;

    //==================================================================
    //  Unary operator >= for integral types.
    //==================================================================

    bool operator >=(const LargeInteger & value) const;

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

    void Copy(const LargeInteger & that);

    void AddPositiveArray(const LargeInteger & addend);

    bool SubtractPositiveArray(const LargeInteger & subtrahend);

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

    bool LessThanPositiveArrayCompare(const LargeInteger & value_0,
                                      const LargeInteger & value_1) const;

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

LargeInteger operator +(const LargeInteger & addend_a,
                        const LargeInteger & addend_b);

//======================================================================
//  Addition of an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator +(const LargeInteger & x_addend,
                        int addend);

LargeInteger operator +(int addend,
                        const LargeInteger & x_addend);

//======================================================================
//  Addition of an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator +(const LargeInteger & x_addend,
                        unsigned int addend);

LargeInteger operator +(unsigned int addend,
                        const LargeInteger & x_addend);

//======================================================================
//  operator - for integral types.
//======================================================================

//======================================================================
//  Subtraction of two instances of this class.
//======================================================================

LargeInteger operator -(const LargeInteger & minuend,
                        const LargeInteger & subtrahend);

//======================================================================
//  Subtraction with an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator -(const LargeInteger & minuend,
                        int subtrahend);

LargeInteger operator -(int minuend,
                        const LargeInteger & subtrahend);

//======================================================================
//  Subtraction with an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator -(const LargeInteger & minuend,
                        unsigned int subtrahend);

LargeInteger operator -(unsigned int minuend,
                        const LargeInteger & subtrahend);

//======================================================================
//  operator * for integral types.
//======================================================================

//======================================================================
//  Multiplication of two instances of this class.
//======================================================================

LargeInteger operator *(const LargeInteger & multiplier_a,
                        const LargeInteger & multiplier_b);

//======================================================================
//  Multiplication of an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator *(const LargeInteger & x_multiplier,
                        int multiplier);

LargeInteger operator *(int multiplier,
                        const LargeInteger & x_multiplier);

//======================================================================
//  Multiplication of an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator *(const LargeInteger & x_multiplier,
                        unsigned int multiplier);

LargeInteger operator *(unsigned int multiplier,
                        const LargeInteger & x_multiplier);

//======================================================================
//  operator / for integral types.
//======================================================================

//======================================================================
//  Division of two instances of this class.
//======================================================================

LargeInteger operator /(const LargeInteger & x_dividend,
                        const LargeInteger & x_denominator);

//======================================================================
//  Division with an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator /(const LargeInteger & x_dividend,
                        int divisor);

LargeInteger operator /(int dividend,
                        const LargeInteger & x_denominator);

//======================================================================
//  Division with an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator /(const LargeInteger & x_dividend,
                        unsigned int divisor);

LargeInteger operator /(unsigned int dividend,
                        const LargeInteger & x_denominator);

//======================================================================
//  operator << for integral types.
//======================================================================

//======================================================================
//  Left shift of two instances of this class.
//======================================================================

LargeInteger operator <<(const LargeInteger & value,
                         const LargeInteger & x_shift_count);

//======================================================================
//  Left shift with an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator <<(const LargeInteger & value,
                         int shift_count);

LargeInteger operator <<(int value,
                         const LargeInteger & x_shift_count);

//======================================================================
//  Left shift with an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator <<(const LargeInteger & value,
                         unsigned int shift_count);

LargeInteger operator <<(unsigned int value,
                         const LargeInteger & x_shift_count);

//======================================================================
//  operator >> for integral types.
//======================================================================

//======================================================================
//  Right shift of two instances of this class.
//======================================================================

LargeInteger operator >>(const LargeInteger & value,
                         const LargeInteger & x_shift_count);

//======================================================================
//  Right shift with an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator >>(const LargeInteger & value,
                         int shift_count);

LargeInteger operator >>(int value,
                         const LargeInteger & x_shift_count);

//======================================================================
//  Right shift with an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator >>(const LargeInteger & value,
                         unsigned int shift_count);

LargeInteger operator >>(unsigned int value,
                         const LargeInteger & x_shift_count);

//======================================================================
//  operator % for integral types.
//======================================================================

//======================================================================
//  Modulus operator with two instances of this class.
//======================================================================

LargeInteger operator %(const LargeInteger & x_dividend,
                        const LargeInteger & x_denominator);

//======================================================================
//  Modulus operator with an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator %(const LargeInteger & x_dividend,
                        int divisor);

LargeInteger operator %(int dividend,
                        const LargeInteger & x_denominator);

//======================================================================
//  Modulus operator with an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator %(const LargeInteger & x_dividend,
                        unsigned int divisor);

LargeInteger operator %(unsigned int dividend,
                        const LargeInteger & x_denominator);

//======================================================================
//  operator ^ for integral types.
//======================================================================

//======================================================================
//  Bitwise Exclusive OR of two instances of this class.
//======================================================================

LargeInteger operator ^(const LargeInteger & value_a,
                        const LargeInteger & value_b);

//======================================================================
//  Bitwise Exclusive OR of an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator ^(const LargeInteger & x_value,
                        int value);

LargeInteger operator ^(int value,
                        const LargeInteger & x_value);

//======================================================================
//  Bitwise Exclusive OR of an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator ^(const LargeInteger & x_value,
                        unsigned int value);

LargeInteger operator ^(unsigned int value,
                        const LargeInteger & x_value);

//======================================================================
//  operator & for integral types.
//======================================================================

//======================================================================
//  Bitwise AND of two instances of this class.
//======================================================================

LargeInteger operator &(const LargeInteger & value_a,
                        const LargeInteger & value_b);

//======================================================================
//  Bitwise AND of an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator &(const LargeInteger & x_value,
                        int value);

LargeInteger operator &(int value,
                        const LargeInteger & x_value);

//======================================================================
//  Bitwise AND of an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator &(const LargeInteger & x_value,
                        unsigned int value);

LargeInteger operator &(unsigned int value,
                        const LargeInteger & x_value);

//======================================================================
//  operator | for integral types.
//======================================================================

//======================================================================
//  Bitwise OR of two instances of this class.
//======================================================================

LargeInteger operator |(const LargeInteger & value_a,
                        const LargeInteger & value_b);

//======================================================================
//  Bitwise OR of an instance of the LargeInteger class and a signed integer.
//======================================================================

LargeInteger operator |(const LargeInteger & x_value,
                        int value);

LargeInteger operator |(int value,
                        const LargeInteger & x_value);

//======================================================================
//  Bitwise OR of an instance of the LargeInteger class and an unsigned integer.
//======================================================================

LargeInteger operator |(const LargeInteger & x_value,
                        unsigned int value);

LargeInteger operator |(unsigned int value,
                        const LargeInteger & x_value);

#endif
