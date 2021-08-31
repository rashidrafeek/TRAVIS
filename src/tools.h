/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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


#ifndef TOOLS_H
#define TOOLS_H


// This must always be the first include directive
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <stdarg.h>
#include <time.h>
#include <wchar.h>
#include <ctype.h>
#include <setjmp.h>
#include <vector>
#include <float.h>
//#include <assert.h>


//#include <sys/time.h>
//#include <termios.h>

class CxString;

/*************************************************************/


#define UNUSED(a) ((void)a)


#define Pi (3.1415926535897932385)

#ifndef M_PI
	#define M_PI (3.1415926535897932385)
#endif

#define MIN(a,b) (((a)<=(b)) ? (a) : (b))
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))

#define MIN3(a,b,c) MIN(MIN(a,b),c)
#define MED3(a,b,c) MAX( MIN(a,b), MIN( MAX(a,b), c ) )
#define MAX3(a,b,c) MAX(MAX(a,b),c)

#define MIN4(a,b,c,d) MIN(MIN(a,b),MIN(c,d))
#define MAX4(a,b,c,d) MAX(MAX(a,b),MAX(c,d))


#define SQR(x)   ((x)*(x))


#define RUNTYPE_ANALYZE  1
#define RUNTYPE_SCANMOL  2


#define LAMMPS_FIELD_ID       1
#define LAMMPS_FIELD_ELEMENT  2
#define LAMMPS_FIELD_XU       3
#define LAMMPS_FIELD_YU       4
#define LAMMPS_FIELD_ZU       5
#define LAMMPS_FIELD_X        6
#define LAMMPS_FIELD_Y        7
#define LAMMPS_FIELD_Z        8
#define LAMMPS_FIELD_VX       9
#define LAMMPS_FIELD_VY      10
#define LAMMPS_FIELD_VZ      11
#define LAMMPS_FIELD_FX      12
#define LAMMPS_FIELD_FY      13
#define LAMMPS_FIELD_FZ      14
#define LAMMPS_FIELD_Q       15
#define LAMMPS_FIELD_MAX     15


#ifdef TARGET_WINDOWS

	inline double roundf(double f) { return floor(f+0.5); }

	inline double nextafter(double from, double to) { return _nextafter(from,to); }

#endif



#ifdef TARGET_WINDOWS

	#define myisnan(X) (_isnan(X))

	#define GREY 8
	#define BLUE 9
	#define GREEN 10
	#define CYAN 11
	#define RED 12
	#define PINK 13
	#define YELLOW 14
	#define WHITE 15

	#define UML_AE ((unsigned char)142)
	#define UML_ae ((unsigned char)132)
	#define UML_OE ((unsigned char)153)
	#define UML_oe ((unsigned char)148)
	#define UML_UE ((unsigned char)154)
	#define UML_ue ((unsigned char)129)
	#define UML_ss ((unsigned char)225)

	#define mypow(X,Y) (pow(X,Y))

	#define iabs(X) (abs(X))

#else // Linux or generic

	#define myisnan(X) (std::isnan(X))

	#define GREY 0
	#define BLUE 4
	#define GREEN 2
	#define CYAN 6
	#define RED 1
	#define PINK 5
	#define YELLOW 3
	#define WHITE 7

	#define UML_AE ((unsigned char)196)
	#define UML_ae ((unsigned char)228)
	#define UML_OE ((unsigned char)214)
	#define UML_oe ((unsigned char)246)
	#define UML_UE ((unsigned char)220)
	#define UML_ue ((unsigned char)252)
	#define UML_ss ((unsigned char)223)

	#define mypow(X,Y) (std::pow(X,Y))

	#define iabs(X) (std::abs(X))

#endif


#include <errno.h>
#include <sys/types.h>


class CxWordArray;
class CxIntArray;

//#define max(a,b)  (((a) > (b)) ? (a) : (b))

//#define SAVEPOS	mprintf(GREEN,"You may jump back to here by entering \"$\" <---\n\n"); if (setjmp(g_JumpBuf)!=0)	mprintf(GREEN,"\n<--- Going Back\n\n");


inline double pow2(double x) { return x*x; }
inline double pow3(double x) { return x*x*x; }
inline double pow4(double x) { return (x*x)*(x*x); }
inline double pow5(double x) { return (x*x*x)*(x*x); }


extern FILE *g_pLogFile;
extern jmp_buf g_JumpBuf;

char* fgets_bin(char *buf, int n, FILE *a);

void SavePosition();
void LoadPosition();

//void AskString(const char *s, char *buf, const char *def, ...);
//void AskString_ND(const char *s, char *buf, ...);


bool IsValidInteger(const char *s);
bool IsValidUnsigned(const char *s);
bool IsValidFloat(const char *s);



#ifdef __GNUG__ // Variadic Argument Type Checking of GCC

	void AskString(const char *s, CxString *buf, const char *def, ...) __attribute__ ((format (printf, 1, 4)));
	void AskString_ND(const char *s, CxString *buf, ...) __attribute__ ((format (printf, 1, 3)));

	bool AskYesNo(const char *s, int def, ...) __attribute__ ((format (printf, 1, 3)));
	int AskUnsignedInteger(const char *s, int def, ...) __attribute__ ((format (printf, 1, 3)));
	double AskFloat(const char *s, double def, ...) __attribute__ ((format (printf, 1, 3)));
	bool AskYesNo_ND(const char *s, ...) __attribute__ ((format (printf, 1, 2)));
	int AskUnsignedInteger_ND(const char *s, ...) __attribute__ ((format (printf, 1, 2)));
	double AskFloat_ND(const char *s, ...) __attribute__ ((format (printf, 1, 2)));
	int AskRangeInteger(const char *s, int mini, int maxi, int def, ...) __attribute__ ((format (printf, 1, 5)));
	int AskRangeInteger_ND(const char *s, int mini, int maxi, ...) __attribute__ ((format (printf, 1, 4)));
	double AskRangeFloat(const char *s, double mini, double maxi, double def, ...) __attribute__ ((format (printf, 1, 5)));
	int AskInteger(const char *s, int def, ...) __attribute__ ((format (printf, 1, 3)));

#else

	void AskString(const char *s, CxString *buf, const char *def, ...);
	void AskString_ND(const char *s, CxString *buf, ...);

	bool AskYesNo(const char *s, int def, ...);
	int AskUnsignedInteger(const char *s, int def, ...);
	double AskFloat(const char *s, double def, ...);
	bool AskYesNo_ND(const char *s, ...);
	int AskUnsignedInteger_ND(const char *s, ...);
	double AskFloat_ND(const char *s, ...);
	int AskRangeInteger(const char *s, int mini, int maxi, int def, ...);
	int AskRangeInteger_ND(const char *s, int mini, int maxi, ...);
	double AskRangeFloat(const char *s, double mini, double maxi, double def, ...);
	int AskInteger(const char *s, int def, ...);

#endif


int AskMolecule(const char *s);

void InitColor();
void TextColor(int fg);
void TextNormal();

int HalfBox();
double HalfBox_Exact();
int HalfBoxSq3();
bool FileExist(const char *s);
void FreeFileName(const char *s);
//void FreeFileName(char *pre, char *s, char *post);
void saxonize(char *buf);


#ifdef __GNUG__  // Variadic Argument Type Checking of GCC

	void mprintf(const char *s, ...) __attribute__ ((format (printf, 1, 2)));
	void mprintf_nos(const char *s, ...) __attribute__ ((format (printf, 1, 2)));
	void mfprintf(FILE *a, const char *s, ...) __attribute__ ((format (printf, 2, 3)));
	void mprintf(int color, const char *s, ...) __attribute__ ((format (printf, 2, 3)));
	void bprintf(const char *s, ...) __attribute__ ((format (printf, 1, 2)));
	void inpprintf(const char *s, ...) __attribute__ ((format (printf, 1, 2)));
	void eprintf(const char *s, ...) __attribute__ ((format (printf, 1, 2)));

#else

	void mprintf(const char *s, ...);
	void mprintf_nos(const char *s, ...);
	void mfprintf(FILE *a, const char *s, ...);
	void mprintf(int color, const char *s, ...);
	void bprintf(const char *s, ...);
	void inpprintf(const char *s, ...);
	void eprintf(const char *s, ...);

#endif


//void myget(char *s);
void myget(CxString *st);
int mystricmp(const char *s1, const char *s2);
void RemoveDoubleBackslash(char *buf);
const char* TreeElement(const char *s);
bool ParseIntList(const char *s, CxIntArray *la);
bool ParseIntList(const char *s, std::vector<int> &la);
bool ParseIntList(const char *s, CxWordArray *wa);
void SortSingleMolAtomTypes();
void xAddAtom(const char *s, bool virt);
bool isinteger(const char *s);

bool ContainsDigit(const char *s);
bool ContainsDigitOrSpecial(const char *s);
//void ReplaceDigits(char *s);
void ReplaceDigits(CxString *s);
const char* RemovePath(const char *s);

FILE *OpenFileWrite(const char *buf, bool text, CxString *out = NULL);
bool IsTTY(FILE *f);
double ZeroDivide(double a, double b);
bool isdigit(char c);
bool isnumeric(char c);
char *FormatBytes(double i);

double dec(double a, double dig);
void decomp(double &a, double &b, double &c);
double maxbound(double a, double r);
double minbound(double a, double r);
double majorticks(double lower, double upper);

int GetSignificantDigit(double d, int sig);
void CreateTicks(double mi, double ma, int &major, int &minor, bool printinfo=true);

double MaxDiff_DoubleArray(double *p, int n);

//void ProtectCharacters(char *dest, const char *src, const char *rep, const char *prot);
void ProtectCharacters(CxString *dest, const char *src, const char *rep, const char *prot);

double NormalDistIntegral(double x);

double GammaLn(double xx);
double FactorialLn(int n);
double BinomialCoeff(int n, int k);

void HSV2RGB(double h, double s, double v, double &r, double &g, double &b);

void RGB2HSV(double r, double g, double b, double &h, double &s, double &v);


inline int floori(int a, int b) {
	if (0 < (a ^ b)) {
		return a / b;
	} else {
		if (a % b != 0)
			return a / b - 1;
		else
			return a / b;
	}
}


inline int PeriodicImage1D( double val, double mi, double ma ) {

	return (int)floor( (val-mi) / (ma-mi) );
}


double Chop(double d, double t);




#endif


