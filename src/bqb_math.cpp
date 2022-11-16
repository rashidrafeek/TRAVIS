/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2022.

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

#include "bqb_math.h"
#include "bqb_interface.h"


const char *GetRevisionInfo_bqb_math(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_math() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}




void CBQBMath::DumpMatrix(const CBQBDMatrix3 &m) const {

	int z/*, z2*/;
	for (z=0;z<3;z++)
	{
/*		for (z2=0;z2<3;z2++)
			if ((GetAt(z2,z) < 0) && (GetAt(z2,z) > -0.000000001))
        			GetAt(z2,z) = 0;*/
		m_IF.printf("( %6.3f %6.3f %6.3f )\n",m.GetAt(0,z),m.GetAt(1,z),m.GetAt(2,z));
	}
}


void CBQBMath::DumpMatrixSmall(const CBQBDMatrix3 &m) const {

	int z/*, z2*/;
	for (z=0;z<3;z++)
	{
/*		for (z2=0;z2<3;z2++)
			if ((GetAt(z2,z) < 0) && (GetAt(z2,z) > -0.000000001))
        			GetAt(z2,z) = 0;*/
		m_IF.printf("( %12.9f %12.9f %12.9f )\n",m.GetAt(0,z),m.GetAt(1,z),m.GetAt(2,z));
	}
}


void CBQBMath::DumpMatrix(const CBQBDMatrixMN &m) const {

	int z, z2;
	for (z=0;z<m.GetRows();z++)
	{
		for (z2=0;z2<m.GetCols();z2++)
			if (fabs(m.GetAt(z,z2)) < 0.000005)
				m_IF.printf(" %11.5f",0.0);
			else
				m_IF.printf(" %11.5f",m.GetAt(z,z2));
		m_IF.printf("\n");
	}
}


void CBQBMath::DumpMatrixSmall(const CBQBDMatrixMN &m) const {

	int z, z2;
	for (z=0;z<m.GetRows();z++)
	{
		for (z2=0;z2<m.GetCols();z2++)
			if (fabs(m.GetAt(z,z2)) < 0.00000005)
				m_IF.printf(" %10.7f",0.0);
			else
				m_IF.printf(" %10.7f",m.GetAt(z,z2));
		m_IF.printf("\n");
	}
}


void CBQBMath::InvertMatrix(CBQBDMatrix3 &m) {

	#define _inva11 m.GetAt(0,0)
	#define _inva12 m.GetAt(0,1)
	#define _inva13 m.GetAt(0,2)
	#define _inva21 m.GetAt(1,0)
	#define _inva22 m.GetAt(1,1)
	#define _inva23 m.GetAt(1,2)
	#define _inva31 m.GetAt(2,0)
	#define _inva32 m.GetAt(2,1)
	#define _inva33 m.GetAt(2,2)

	CBQBDMatrix3 t;
	double det;

	det = _inva11 * (_inva33*_inva22 - _inva32*_inva23)
	    - _inva21 * (_inva33*_inva12 - _inva32*_inva13)
	    + _inva31 * (_inva23*_inva12 - _inva22*_inva13);

	if (det == 0) {
		m_IF.eprintf("CBQBMath::InvertMatrix(): Error: Matrix is singular.\n");
		return;
	}

	t.GetAt(0,0) =  (_inva33*_inva22 - _inva32*_inva23);
	t.GetAt(0,1) = -(_inva33*_inva12 - _inva32*_inva13);
	t.GetAt(0,2) =  (_inva23*_inva12 - _inva22*_inva13);
	t.GetAt(1,0) = -(_inva33*_inva21 - _inva31*_inva23);
	t.GetAt(1,1) =  (_inva33*_inva11 - _inva31*_inva13);
	t.GetAt(1,2) = -(_inva23*_inva11 - _inva21*_inva13);
	t.GetAt(2,0) =  (_inva32*_inva21 - _inva31*_inva22);
	t.GetAt(2,1) = -(_inva32*_inva11 - _inva31*_inva12);
	t.GetAt(2,2) =  (_inva22*_inva11 - _inva21*_inva12);

	t *= 1.0/det;

	m = t;
}


double CBQBMath::Angle(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2) const {
	double t;
	if ((vec1.GetLength() == 0) || (vec2.GetLength() == 0)) {
		m_IF.eprintf("CBQBMath::Angle(): Warning: Indeterminate angle between ( %f | %f | %f ) and ( %f | %f | %f ).\n",vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2]);
		return 0;
	}
	t = ::DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;
	return acos(t);
}


double CBQBMath::Angle_Deg(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2) const {
	double t;
	if ((vec1.GetLength() == 0) || (vec2.GetLength() == 0)) {
		m_IF.eprintf("CBQBMath::Angle_Deg(): Warning: Indeterminate angle between ( %f | %f | %f ) and ( %f | %f | %f ).\n",vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2]);
		return 0;
	}
	t = ::DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;
	t = acos(t);
	return fabs(t*180.0 / Pi);
}


double CBQBMath::DotP(const CBQBDVectorN &vec1, const CBQBDVectorN &vec2) const {
	if (vec1.GetDim() != vec2.GetDim()) {
		m_IF.eprintf("double CBQBMath::DotP(const CBQBDVectorN &, const CBQBDVectorN &): Dimension mismatch (%d vs %d).\n",vec1.GetDim(),vec2.GetDim());
		abort();
	}
	int z;
	double f=0;
	for (z=0;z<vec1.GetDim();z++)
		f += vec1[z]*vec2[z];
	return f;
}


double CBQBMath::Angle(const CBQBDVectorN &vec1, const CBQBDVectorN &vec2) const {
	double t;
	if ((vec1.GetLength() == 0) || (vec2.GetLength() == 0)) {
		m_IF.eprintf("CBQBMath::Angle_Deg(): Warning: Indeterminate angle.\n");
		return 0;
	}
	t = DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;
	return acos(t);
}


double CBQBMath::Angle_Deg(const CBQBDVectorN &vec1, const CBQBDVectorN &vec2) const {
	double t;
	if ((vec1.GetLength() == 0) || (vec2.GetLength() == 0)) {
		m_IF.eprintf("CBQBMath::Angle_Deg(): Warning: Indeterminate angle.\n");
		return 0;
	}
	t = DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;
	t = acos(t);
	return fabs(t*180.0 / Pi);
}




#ifndef BQB_INSIDE_TRAVIS




void CBQBDMatrixMN::Dump(FILE *a) const
{
	int z, z2;
	for (z=0;z<m_iRows;z++)
	{
		for (z2=0;z2<m_iCols;z2++)
			if (fabs(GetAt(z,z2)) < 0.000005)
				fprintf(a," %11.5f",0.0);
			else
				fprintf(a," %11.5f",GetAt(z,z2));
		fprintf(a,"\n");
	}
}


void CBQBDMatrixMN::DumpSmall(FILE *a) const
{
	int z, z2;
	for (z=0;z<m_iRows;z++)
	{
		for (z2=0;z2<m_iCols;z2++)
			if (fabs(GetAt(z,z2)) < 0.00000005)
				fprintf(a," %10.7f",0.0);
			else
				fprintf(a," %10.7f",GetAt(z,z2));
		fprintf(a,"\n");
	}
}



#endif



