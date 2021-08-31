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



// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <stdarg.h>
#include "bqb_interface.h"
#include "bqb_format.h"


const char *GetRevisionInfo_bqb_tools(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_tools() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



int bqb_strcmp_nocase(const char *s1, const char *s2) {

    char *buf, *buf2;
    unsigned long z;
    int i;


    buf = new char[strlen(s1)+1];
    buf2 = new char[strlen(s2)+1];
    for (z=0;z<strlen(s1);z++)
        buf[z] = (char)toupper(s1[z]);
    buf[strlen(s1)] = 0;
    for (z=0;z<strlen(s2);z++)
        buf2[z] = (char)toupper(s2[z]);
    buf2[strlen(s2)] = 0;
    i = strcmp(buf,buf2);
    delete[] buf;
    delete[] buf2;
    return i;
}



void PushNumeration(std::vector<int> &ia, int n, int s1, int s2) {

/*	int t, i, z, a;

//	printf("Push %d ...\n",n);
	t = (int)(floor(mylog2(n+1))+0.5);
//	printf("t=%d\n",t);
	i = n;
	for (z=0;z<t;z++) {
		a = i - 2*(((int)(ceil(i/2.0)+0.5)) - 1);
//		printf("%d\n",a);
		if (a == 1)
			ia.push_back(s1);
		else
			ia.push_back(s2);
		i = ((int)(ceil(i/2.0)+0.5)) - 1;
	}*/

	int z2;

	z2 = n-1;
	while (true) {
		if (z2 & 1)
			ia.push_back(s2);
		else
			ia.push_back(s1);
		if (z2 < 2)
			break;
		z2 = (z2 - 2) / 2;
	}
}



int GetAtomOrd(const char *s) {

	// First the one-letter labels ...
	if (bqb_strcmp_nocase(s, "H"  ) == 0) return   1;
	if (bqb_strcmp_nocase(s, "D"  ) == 0) return   1;
	if (bqb_strcmp_nocase(s, "B"  ) == 0) return   5;
	if (bqb_strcmp_nocase(s, "C"  ) == 0) return   6;
	if (bqb_strcmp_nocase(s, "N"  ) == 0) return   7;
	if (bqb_strcmp_nocase(s, "O"  ) == 0) return   8;
	if (bqb_strcmp_nocase(s, "F"  ) == 0) return   9;
	if (bqb_strcmp_nocase(s, "P"  ) == 0) return  15;
	if (bqb_strcmp_nocase(s, "S"  ) == 0) return  16;
	if (bqb_strcmp_nocase(s, "K"  ) == 0) return  19;
	if (bqb_strcmp_nocase(s, "V"  ) == 0) return  23;
	if (bqb_strcmp_nocase(s, "Y"  ) == 0) return  39;
	if (bqb_strcmp_nocase(s, "I"  ) == 0) return  53;
	if (bqb_strcmp_nocase(s, "W"  ) == 0) return  74;
	if (bqb_strcmp_nocase(s, "U"  ) == 0) return  92;

	// Then the two-letter labels.
	if (bqb_strcmp_nocase(s, "He" ) == 0) return   2;
	if (bqb_strcmp_nocase(s, "Li" ) == 0) return   3;
	if (bqb_strcmp_nocase(s, "Be" ) == 0) return   4;
	if (bqb_strcmp_nocase(s, "Ne" ) == 0) return  10;
	if (bqb_strcmp_nocase(s, "Na" ) == 0) return  11;
	if (bqb_strcmp_nocase(s, "Mg" ) == 0) return  12;
	if (bqb_strcmp_nocase(s, "Al" ) == 0) return  13;
	if (bqb_strcmp_nocase(s, "Si" ) == 0) return  14;
	if (bqb_strcmp_nocase(s, "Cl" ) == 0) return  17;
	if (bqb_strcmp_nocase(s, "Ar" ) == 0) return  18;
	if (bqb_strcmp_nocase(s, "Ca" ) == 0) return  20;
	if (bqb_strcmp_nocase(s, "Sc" ) == 0) return  21;
	if (bqb_strcmp_nocase(s, "Ti" ) == 0) return  22;
	if (bqb_strcmp_nocase(s, "Cr" ) == 0) return  24;
	if (bqb_strcmp_nocase(s, "Mn" ) == 0) return  25;
	if (bqb_strcmp_nocase(s, "Fe" ) == 0) return  26;
	if (bqb_strcmp_nocase(s, "Co" ) == 0) return  27;
	if (bqb_strcmp_nocase(s, "Ni" ) == 0) return  28;
	if (bqb_strcmp_nocase(s, "Cu" ) == 0) return  29;
	if (bqb_strcmp_nocase(s, "Zn" ) == 0) return  30;
	if (bqb_strcmp_nocase(s, "Ga" ) == 0) return  31;
	if (bqb_strcmp_nocase(s, "Ge" ) == 0) return  32;
	if (bqb_strcmp_nocase(s, "As" ) == 0) return  33;
	if (bqb_strcmp_nocase(s, "Se" ) == 0) return  34;
	if (bqb_strcmp_nocase(s, "Br" ) == 0) return  35;
	if (bqb_strcmp_nocase(s, "Kr" ) == 0) return  36;
	if (bqb_strcmp_nocase(s, "Rb" ) == 0) return  37;
	if (bqb_strcmp_nocase(s, "Sr" ) == 0) return  38;
	if (bqb_strcmp_nocase(s, "Zr" ) == 0) return  40;
	if (bqb_strcmp_nocase(s, "Nb" ) == 0) return  41;
	if (bqb_strcmp_nocase(s, "Mo" ) == 0) return  42;
	if (bqb_strcmp_nocase(s, "Tc" ) == 0) return  43;
	if (bqb_strcmp_nocase(s, "Ru" ) == 0) return  44;
	if (bqb_strcmp_nocase(s, "Rh" ) == 0) return  45;
	if (bqb_strcmp_nocase(s, "Pd" ) == 0) return  46;
	if (bqb_strcmp_nocase(s, "Ag" ) == 0) return  47;
	if (bqb_strcmp_nocase(s, "Cd" ) == 0) return  48;
	if (bqb_strcmp_nocase(s, "In" ) == 0) return  49;
	if (bqb_strcmp_nocase(s, "Sn" ) == 0) return  50;
	if (bqb_strcmp_nocase(s, "Sb" ) == 0) return  51;
	if (bqb_strcmp_nocase(s, "Te" ) == 0) return  52;
	if (bqb_strcmp_nocase(s, "Xe" ) == 0) return  54;
	if (bqb_strcmp_nocase(s, "Cs" ) == 0) return  55;
	if (bqb_strcmp_nocase(s, "Ba" ) == 0) return  56;
	if (bqb_strcmp_nocase(s, "La" ) == 0) return  57;
	if (bqb_strcmp_nocase(s, "Ce" ) == 0) return  58;
	if (bqb_strcmp_nocase(s, "Pr" ) == 0) return  59;
	if (bqb_strcmp_nocase(s, "Nd" ) == 0) return  60;
	if (bqb_strcmp_nocase(s, "Pm" ) == 0) return  61;
	if (bqb_strcmp_nocase(s, "Sm" ) == 0) return  62;
	if (bqb_strcmp_nocase(s, "Eu" ) == 0) return  63;
	if (bqb_strcmp_nocase(s, "Gd" ) == 0) return  64;
	if (bqb_strcmp_nocase(s, "Tb" ) == 0) return  65;
	if (bqb_strcmp_nocase(s, "Dy" ) == 0) return  66;
	if (bqb_strcmp_nocase(s, "Ho" ) == 0) return  67;
	if (bqb_strcmp_nocase(s, "Er" ) == 0) return  68;
	if (bqb_strcmp_nocase(s, "Tm" ) == 0) return  69;
	if (bqb_strcmp_nocase(s, "Yb" ) == 0) return  70;
	if (bqb_strcmp_nocase(s, "Lu" ) == 0) return  71;
	if (bqb_strcmp_nocase(s, "Hf" ) == 0) return  72;
	if (bqb_strcmp_nocase(s, "Ta" ) == 0) return  73;
	if (bqb_strcmp_nocase(s, "Re" ) == 0) return  75;
	if (bqb_strcmp_nocase(s, "Os" ) == 0) return  76;
	if (bqb_strcmp_nocase(s, "Ir" ) == 0) return  77;
	if (bqb_strcmp_nocase(s, "Pt" ) == 0) return  78;
	if (bqb_strcmp_nocase(s, "Au" ) == 0) return  79;
	if (bqb_strcmp_nocase(s, "Hg" ) == 0) return  80;
	if (bqb_strcmp_nocase(s, "Tl" ) == 0) return  81;
	if (bqb_strcmp_nocase(s, "Pb" ) == 0) return  82;
	if (bqb_strcmp_nocase(s, "Bi" ) == 0) return  83;
	if (bqb_strcmp_nocase(s, "Po" ) == 0) return  84;
	if (bqb_strcmp_nocase(s, "At" ) == 0) return  85;
	if (bqb_strcmp_nocase(s, "Rn" ) == 0) return  86;
	if (bqb_strcmp_nocase(s, "Fr" ) == 0) return  87;
	if (bqb_strcmp_nocase(s, "Ra" ) == 0) return  88;
	if (bqb_strcmp_nocase(s, "Ac" ) == 0) return  89;
	if (bqb_strcmp_nocase(s, "Th" ) == 0) return  90;
	if (bqb_strcmp_nocase(s, "Pa" ) == 0) return  91;
	if (bqb_strcmp_nocase(s, "Np" ) == 0) return  93;
	if (bqb_strcmp_nocase(s, "Pu" ) == 0) return  94;
	if (bqb_strcmp_nocase(s, "Am" ) == 0) return  95;
	if (bqb_strcmp_nocase(s, "Cm" ) == 0) return  96;
	if (bqb_strcmp_nocase(s, "Bk" ) == 0) return  97;
	if (bqb_strcmp_nocase(s, "Cf" ) == 0) return  98;
	if (bqb_strcmp_nocase(s, "Es" ) == 0) return  99;
	if (bqb_strcmp_nocase(s, "Fm" ) == 0) return 100;
	if (bqb_strcmp_nocase(s, "Md" ) == 0) return 101;
	if (bqb_strcmp_nocase(s, "No" ) == 0) return 102;
	if (bqb_strcmp_nocase(s, "Lr" ) == 0) return 103;

	return 200;
}



const char* GetAtomOrdLabel( int ord ) {

	switch( ord ) {
		case   1: return "H" ;
		case   2: return "He";
		case   3: return "Li";
		case   4: return "Be";
		case   5: return "B" ;
		case   6: return "C" ;
		case   7: return "N" ;
		case   8: return "O" ;
		case   9: return "F" ;
		case  10: return "Ne";
		case  11: return "Na";
		case  12: return "Mg";
		case  13: return "Al";
		case  14: return "Si";
		case  15: return "P" ;
		case  16: return "S" ;
		case  17: return "Cl";
		case  18: return "Ar";
		case  19: return "K" ;
		case  20: return "Ca";
		case  21: return "Sc";
		case  22: return "Ti";
		case  23: return "V" ;
		case  24: return "Cr";
		case  25: return "Mn";
		case  26: return "Fe";
		case  27: return "Co";
		case  28: return "Ni";
		case  29: return "Cu";
		case  30: return "Zn";
		case  31: return "Ga";
		case  32: return "Ge";
		case  33: return "As";
		case  34: return "Se";
		case  35: return "Br";
		case  36: return "Kr";
		case  37: return "Rb";
		case  38: return "Sr";
		case  39: return "Y" ;
		case  40: return "Zr";
		case  41: return "Nb";
		case  42: return "Mo";
		case  43: return "Tc";
		case  44: return "Ru";
		case  45: return "Rh";
		case  46: return "Pd";
		case  47: return "Ag";
		case  48: return "Cd";
		case  49: return "In";
		case  50: return "Sn";
		case  51: return "Sb";
		case  52: return "Te";
		case  53: return "I" ;
		case  54: return "Xe";
		case  55: return "Cs";
		case  56: return "Ba";
		case  57: return "La";
		case  58: return "Ce";
		case  59: return "Pr";
		case  60: return "Nd";
		case  61: return "Pm";
		case  62: return "Sm";
		case  63: return "Eu";
		case  64: return "Gd";
		case  65: return "Tb";
		case  66: return "Dy";
		case  67: return "Ho";
		case  68: return "Er";
		case  69: return "Tm";
		case  70: return "Yb";
		case  71: return "Lu";
		case  72: return "Hf";
		case  73: return "Ta";
		case  74: return "W" ;
		case  75: return "Re";
		case  76: return "Os";
		case  77: return "Ir";
		case  78: return "Pt";
		case  79: return "Au";
		case  80: return "Hg";
		case  81: return "Tl";
		case  82: return "Pb";
		case  83: return "Bi";
		case  84: return "Po";
		case  85: return "At";
		case  86: return "Rn";
		case  87: return "Fr";
		case  88: return "Ra";
		case  89: return "Ac";
		case  90: return "Th";
		case  91: return "Pa";
		case  92: return "U" ;
		case  93: return "Np";
		case  94: return "Pu";
		case  95: return "Am";
		case  96: return "Cm";
		case  97: return "Bk";
		case  98: return "Cf";
		case  99: return "Es";
		case 100: return "Fm";
		case 101: return "Md";
		case 102: return "No";
		case 103: return "Lr";
	}
	return "X";
}



bool BQBFileExist(const char *s) {

	FILE *a;

	a = fopen(s,"rb");
	if (a == NULL)
		return false;
	fclose(a);
	return true;
}



FILE* CBQBTools::BQBOpenFileWrite(const char *buf, bool text) {

	FILE *a;


	if (text)
		a = fopen(buf,"wt");
	else
		a = fopen(buf,"wb");
	if (a == NULL) {
		m_IF.eprintf("Could not open file \"%s\" for writing.\n",buf);
		m_IF.eprintf("Be sure to have writing permission for this directory.\n");
		m_IF.eprintf("\n");
		m_IF.eprintf("Aborting execution.\n");
		abort();
	}

	return a;
}



long CBQBTools::FindIndexPosition(FILE *a) {

	unsigned int i;
	unsigned char uc;
	long l;
	char buf[3];

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf(">>> FindIndexPosition >>>\n");

	if (a == NULL) {
		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("    File pointer == NULL.\n");
			m_IF.printf("<<< FindIndexPosition <<<\n");
			m_IF.printf("\n");
		}
		return 0;
	}

	fseek(a,-3,SEEK_END);
	(void)!fread(buf,3,1,a);

//	if (memcmp(buf,"IDX",3) != 0) {
	if ((buf[0] != 'I') || (buf[1] != 'D') || (buf[2] != 'X')) { // To prevent Valgrind "uninitialized value" warning
		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("    Last three characters differ from \"IDX\".\n");
			m_IF.printf("<<< FindIndexPosition <<<\n");
			m_IF.printf("\n");
		}
		return 0;
	} else if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("    Found \"IDX\".\n");

	fseek(a,-7,SEEK_END);

	(void)!fread(&uc,1,1,a);
	i = ((unsigned int)uc) << 24;
	(void)!fread(&uc,1,1,a);
	i += ((unsigned int)uc) << 16;
	(void)!fread(&uc,1,1,a);
	i += ((unsigned int)uc) << 8;
	(void)!fread(&uc,1,1,a);
	i += uc;

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("    Found index frame length %u.\n",i);

	fseek(a,-((int)i),SEEK_END);

	l = ftell(a);

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("    Found index frame offset %ld.\n",l);
		m_IF.printf("<<< FindIndexPosition <<<\n");
		m_IF.printf("\n");
	}

	return l;
}



void BQBFormatTime(unsigned long eta, char *buf) {

	char tbuf[256], tbuf2[256];


	if (buf == NULL)
		return;

	if ((eta/60) > 0)
		sprintf(tbuf,"%02lus",eta%60);
	else
		sprintf(tbuf,"%2lus",eta%60);

	eta /= 60;
	if (eta > 0) {

		strcpy(tbuf2,tbuf);

		if ((eta/60) > 0)
			sprintf(tbuf,"%02lum",eta%60);
		else
			sprintf(tbuf,"%2lum",eta%60);

		strcat(tbuf,tbuf2);

		eta /= 60;
		if (eta > 0) {

			strcpy(tbuf2,tbuf);

			if ((eta/60) > 0)
				sprintf(tbuf,"%02luh",eta%24);
			else
				sprintf(tbuf,"%2luh",eta%24);

			strcat(tbuf,tbuf2);

			eta /= 24;
		}
		if (eta > 0) {

			strcpy(tbuf2,tbuf);

			if ((eta/60) > 0)
				sprintf(tbuf,"%02lud",eta);
			else
				sprintf(tbuf,"%2lud",eta);

			strcat(tbuf,tbuf2);
		}
	}

	strcpy(buf,tbuf);
}



bool bqbisinteger(const char *s) {
	int z;
	z = 0;
	if (strlen(s) > 0)
		if (s[0] == '-')
			z++;
	for (;z<(int)strlen(s);z++)
		if ((s[z] < '0') || (s[z] > '9'))
			return false;
	return true;
}



void CBQBStatistics::PushStatistics() {

	m_oaStatStack.push_back(m_oStat);
}



void CBQBStatistics::PopStatistics() {

	if (m_oaStatStack.size() == 0) {
		m_IF.eprintf("CBQBStatistics::PopStatistics(): Error: Statistics stack is empty.\n");
		abort();
	}

	m_oStat.m_lSize = m_oaStatStack.back().m_lSize;
	m_oStat.m_lOverhead = m_oaStatStack.back().m_lOverhead;
	m_oStat.m_lAlphabet = m_oaStatStack.back().m_lAlphabet;
	m_oStat.m_lHuffmanTables = m_oaStatStack.back().m_lHuffmanTables;
	m_oStat.m_lTableSwitch = m_oaStatStack.back().m_lTableSwitch;
	m_oStat.m_lHuffmanData = m_oaStatStack.back().m_lHuffmanData;

	m_oaStatStack.pop_back();
}



void CBQBStatistics::PopDiffStatistics() {

	if (m_oaStatStack.size() < 2) {
		m_IF.eprintf("CBQBStatistics::PopDiffStatistics(): Error: Statistics stack has size of %lu < 2.\n",(unsigned long)m_oaStatStack.size());
		abort();
	}

	m_oStat.m_lSize += m_oaStatStack[m_oaStatStack.size()-2].m_lSize - m_oaStatStack.back().m_lSize;
	m_oStat.m_lOverhead += m_oaStatStack[m_oaStatStack.size()-2].m_lOverhead - m_oaStatStack.back().m_lOverhead;
	m_oStat.m_lAlphabet += m_oaStatStack[m_oaStatStack.size()-2].m_lAlphabet - m_oaStatStack.back().m_lAlphabet;
	m_oStat.m_lHuffmanTables += m_oaStatStack[m_oaStatStack.size()-2].m_lHuffmanTables - m_oaStatStack.back().m_lHuffmanTables;
	m_oStat.m_lTableSwitch += m_oaStatStack[m_oaStatStack.size()-2].m_lTableSwitch - m_oaStatStack.back().m_lTableSwitch;
	m_oStat.m_lHuffmanData += m_oaStatStack[m_oaStatStack.size()-2].m_lHuffmanData - m_oaStatStack.back().m_lHuffmanData;

	m_oaStatStack.pop_back();
}



void CBQBStatistics::ResetStatistics() {

	m_oaStatStack.clear();
	m_oStat.reset();
}



void CBQBStatistics::PopIgnoreStatistics() {

	if (m_oaStatStack.size() == 0) {
		m_IF.eprintf("CBQBStatistics::PopStatistics(): Error: Statistics stack is empty.\n");
		abort();
	}

	m_oaStatStack.pop_back();
}


