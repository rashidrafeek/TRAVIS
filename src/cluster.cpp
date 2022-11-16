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

#include "cluster.h"
#include "grace.h"
#include "maintools.h"
#include "globalvar.h"


const char *GetRevisionInfo_cluster(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_cluster() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


CClusterNode::CClusterNode()
{
	m_vaAtoms.SetName("CClusterNode::m_vaAtoms");
	m_iaAtoms.SetName("CClusterNode::m_iaAtoms");
}


CClusterNode::~CClusterNode()
{
}


CClusterAnalysis::CClusterAnalysis()
{
	m_poaClusterPoly = NULL;
	m_pDistCache = NULL;
	m_pClusterDistanceDF = NULL;
	m_pPolymerDF = NULL;
	m_pClusterCountDF = NULL;
	m_pClusterDistributionDF = NULL;
	m_pClusterSizeDF = NULL;
	m_pClusterDistance2D = NULL;
	m_pClusterDistribution2D = NULL;
	m_pClusterSize2D = NULL;
	m_pClusterCount2D = NULL;
	m_pRandom = NULL;

	m_oaAtomGroups.SetName("CClusterAnalysis::m_oaAtomGroups");
	m_iaAtomList.SetName("CClusterAnalysis::m_iaAtomList");
	m_oaClusters.SetName("CClusterAnalysis::m_oaClusters");
	m_oaBaseClusters.SetName("CClusterAnalysis::m_oaBaseClusters");
	m_oaTopClusters.SetName("CClusterAnalysis::m_oaTopClusters");
	m_faMinDist.SetName("CClusterAnalysis::m_faMinDist");
	m_faMaxDist.SetName("CClusterAnalysis::m_faMaxDist");
	m_faHetNumber.SetName("CClusterAnalysis::m_faHetNumber");
	m_faCenter.SetName("CClusterAnalysis::m_faCenter");
	m_faIntHetNumber.SetName("CClusterAnalysis::m_faIntHetNumber");
}


CClusterAnalysis::~CClusterAnalysis()
{
	if (m_poaClusterPoly != NULL)
		delete[] m_poaClusterPoly;
	if (m_pDistCache != NULL)
		delete[] m_pDistCache;
}


void CClusterAnalysis::AddCluster(int i)
{
	CClusterNode *n;

	try { n = new CClusterNode(); } catch(...) { n = NULL; }
	if (n == NULL) NewException((double)sizeof(CClusterNode),__FILE__,__LINE__,__PRETTY_FUNCTION__);

//	n->m_iParent = -1;
	n->m_iIndex = m_oaClusters.GetSize();
	n->m_fDistance = 1.0E30;
	n->m_iMonomers = 1;
	n->m_pParent = NULL;
	if (m_bDistCache)
		n->m_iaAtoms.SetMaxSize(i);
			else n->m_vaAtoms.SetMaxSize(i);
	n->m_iChildren[0] = -1;

	m_oaClusters.Add(n);
	m_oaBaseClusters.Add(n);
	m_poaClusterPoly[0].Add(n);
}


void CClusterAnalysis::AddParticle(double x, double y, double z)
{
	CClusterNode *n;

	n = (CClusterNode*)m_oaClusters[m_oaClusters.GetSize()-1];
	n->m_vaAtoms.Add(CxDVector3(x,y,z));
}


void CClusterAnalysis::AddParticle(int i)
{
	CClusterNode *n;

	n = (CClusterNode*)m_oaClusters[m_oaClusters.GetSize()-1];
	n->m_iaAtoms.Add(i);
}


void CClusterAnalysis::BuildTree()
{
	int z, i, j, st/*, m*/;
	double d;
	CClusterNode *n, *c1, *c2;

	n = NULL;

//	mprintf("BuildTree online.\n");
	for (z=0;z<m_oaClusters.GetSize();z++)
		m_oaTopClusters.Add(m_oaClusters[z]);

	for (z=0;z<m_oaClusters.GetSize();z++)
	{
		n = (CClusterNode*)m_oaClusters[z];
		FindNearestNeighbor(n);
	}

//	mprintf("Initialization done.\n");

	st = 0;
	while (m_oaTopClusters.GetSize() > 1)
	{
		st++;
//		mprintf("Step %d, %d Top-Clusters remaining.\n",st,m_oaTopClusters.GetSize());
		d = 1.0E30;
		i = -1;
		for (z=0;z<m_oaTopClusters.GetSize();z++)
		{
			if (((CClusterNode*)m_oaTopClusters[z])->m_fNNDistance < d)
			{
				d = ((CClusterNode*)m_oaTopClusters[z])->m_fNNDistance;
				i = z;
			}
		}
		if (i == -1) {
			eprintf("CClusterAnalysis::BuildTree(): Internal error.\n");
			abort();
		}
		j = ((CClusterNode*)m_oaTopClusters[i])->m_iNextNeighbor;

		c1 = (CClusterNode*)m_oaTopClusters[i];
		c2 = (CClusterNode*)m_oaClusters[j];

		try { n = new CClusterNode(); } catch(...) { n = NULL; }
		if (n == NULL) NewException((double)sizeof(CClusterNode),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if (m_bDistCache)
			n->m_iaAtoms.SetMaxSize(c1->m_iaAtoms.GetSize()+c2->m_iaAtoms.GetSize());
				else n->m_vaAtoms.SetMaxSize(c1->m_vaAtoms.GetSize()+c2->m_vaAtoms.GetSize());

		n->m_pParent = NULL;
//		n->m_iParent = -1;
		n->m_fDistance = d;
//		n->m_laChildren.Add(((CClusterNode*)m_oaTopClusters[i])->m_iIndex);
//		n->m_laChildren.Add(j);
		n->m_iChildren[0] = c1->m_iIndex;
		n->m_iChildren[1] = j;
		n->m_iMonomers = c1->m_iMonomers + c2->m_iMonomers;
//		m = n->m_iMonomers - 1;
		m_poaClusterPoly[n->m_iMonomers-1].Add(n);
//		mprintf("  Combining %d and %d, d=%f...\n",((CClusterNode*)m_oaTopClusters[i])->m_iIndex+1,j+1,d);
		n->m_iIndex = m_oaClusters.GetSize();

/*		m_pClusterSizeDF->AddToBinInt_Fast(m);

		m_pClusterDistribution2DF->AddToBinInt_Fast(m,d-((CClusterNode*)m_oaTopClusters[i])->m_fDistance);

		if (((CClusterNode*)m_oaTopClusters[i])->m_iMonomers == 1)
			m_pClusterDistribution2DF->AddToBinInt_Fast(0,d);
				else m_pClusterDistribution2DF->AddToBinInt_Fast(((CClusterNode*)m_oaTopClusters[i])->m_iMonomers-1,d-((CClusterNode*)m_oaTopClusters[i])->m_fDistance);

		if (((CClusterNode*)m_oaClusters[j])->m_iMonomers == 1)
			m_pClusterDistribution2DF->AddToBinInt_Fast(0,d);
				else m_pClusterDistribution2DF->AddToBinInt_Fast(((CClusterNode*)m_oaClusters[j])->m_iMonomers-1,d-((CClusterNode*)m_oaClusters[j])->m_fDistance);
*/
//		((CClusterNode*)m_oaTopClusters[i])->m_iParent = n->m_iIndex;
		((CClusterNode*)m_oaTopClusters[i])->m_pParent = n;

//		((CClusterNode*)m_oaClusters[j])->m_iParent = n->m_iIndex;
		((CClusterNode*)m_oaClusters[j])->m_pParent = n;

		if (m_bDistCache)
		{
			for (z=0;z<c1->m_iaAtoms.GetSize();z++)
				n->m_iaAtoms.Add(c1->m_iaAtoms[z]);

			for (z=0;z<c2->m_iaAtoms.GetSize();z++)
				n->m_iaAtoms.Add(c2->m_iaAtoms[z]);
		} else
		{
			for (z=0;z<c1->m_vaAtoms.GetSize();z++)
				n->m_vaAtoms.Add(c1->m_vaAtoms[z]);

			for (z=0;z<c2->m_vaAtoms.GetSize();z++)
				n->m_vaAtoms.Add(c2->m_vaAtoms[z]);
		}

		m_oaClusters.Add(n);
		m_oaTopClusters.Add(n);
		m_oaTopClusters.RemoveAt_NoShrink(i,1);
		for (z=0;z<m_oaTopClusters.GetSize();z++)
		{
			if (((CClusterNode*)m_oaTopClusters[z])->m_iIndex == j)
			{
				m_oaTopClusters.RemoveAt_NoShrink(z,1);
				break;
			}
		}
//		for (z=0;z<m_oaTopClusters.GetSize();z++)
//			FindNearestNeighbor((CClusterNode*)m_oaTopClusters[z]);

		FindNearestNeighbor(n);

		TraceNeighbors();
	}
	m_pTop = n;
//	m_pClusterDistribution2DF->AddToBinInt_Fast(m_iMonomers-1,m_fMaxDist-n->m_fDistance);
//	mprintf("BuildTree done.\n");
}


void CClusterAnalysis::FindNearestNeighbor(CClusterNode *n)
{
	int z2, z3, z4;
	double t;
	CClusterNode *n2;

	n->m_fNNDistance = 1.0E30;
	n->m_iNextNeighbor = -1;

	if (m_bDistCache)
	{
		for (z2=0;z2<m_oaTopClusters.GetSize();z2++)
		{
			n2 = (CClusterNode*)m_oaTopClusters[z2];
			if (n->m_iIndex == n2->m_iIndex)
				continue;
			for (z3=0;z3<n->m_iaAtoms.GetSize();z3++)
			{
				for (z4=0;z4<n2->m_iaAtoms.GetSize();z4++)
				{
	//				mprintf("  z2=%d, z3=%d, z4=%d, dist=%f.\n",z2,z3,z4,n->m_fNNDistance);
					t = Distance_Cache(n->m_iaAtoms[z3],n2->m_iaAtoms[z4]);
					if (t < n->m_fNNDistance)
					{
						n->m_fNNDistance = t;
						n->m_iNextNeighbor = n2->m_iIndex;
					}
				}
			}
		}
	} else
	{
		for (z2=0;z2<m_oaTopClusters.GetSize();z2++)
		{
			n2 = (CClusterNode*)m_oaTopClusters[z2];
			if (n->m_iIndex == n2->m_iIndex)
				continue;
			for (z3=0;z3<n->m_vaAtoms.GetSize();z3++)
			{
				for (z4=0;z4<n2->m_vaAtoms.GetSize();z4++)
				{
	//				mprintf("  z2=%d, z3=%d, z4=%d, dist=%f.\n",z2,z3,z4,n->m_fNNDistance);
					t = Distance(&n->m_vaAtoms[z3],&n2->m_vaAtoms[z4]);
					if (t < n->m_fNNDistance)
					{
						n->m_fNNDistance = t;
						n->m_iNextNeighbor = n2->m_iIndex;
					}
				}
			}
		}
	}
}


void CClusterAnalysis::TraceNeighbors()
{
	int z;
	CClusterNode *n, *n2;

	for (z=0;z<m_oaTopClusters.GetSize();z++)
	{
		n = (CClusterNode*)m_oaTopClusters[z];

		if (n->m_iNextNeighbor == -1)
			continue;

		n2 = (CClusterNode*)m_oaClusters[n->m_iNextNeighbor];

//		while (n2->m_iParent != -1)
//			n2 = (CClusterNode*)m_oaClusters[n2->m_iParent];

		while (n2->m_pParent != NULL)
			n2 = n2->m_pParent;

		n->m_iNextNeighbor = n2->m_iIndex;
	}
}


void CClusterAnalysis::DumpDot(const char *s)
{
	FILE *a;
	int z;
	CClusterNode *n;

	a = OpenFileWrite(s,true);

	mfprintf(a,"digraph test {\n");
	for (z=0;z<m_oaClusters.GetSize();z++)
	{
		n = (CClusterNode*)m_oaClusters[z];
//		if (n->m_laChildren.GetSize() == 0)
		if (n->m_iChildren[0] == -1)
			mfprintf(a,"  n%d [label=\"%d X%f\"];\n",z+1,n->m_iIndex+1,n->m_fPosX);
				else mfprintf(a,"  n%d [label=\"X%f Y%.2f\"];\n",z+1,n->m_fPosX,n->m_fDistance);
	}

	mfprintf(a,"\n");

	REC_DumpDot(a,m_pTop);

	mfprintf(a,"}\n");

	fclose(a);
}


void CClusterAnalysis::REC_DumpDot(FILE *a, CClusterNode *n)
{
/*	int z;

	for (z=0;z<n->m_laChildren.GetSize();z++)
	{
		mfprintf(a,"  n%d -> n%d;\n",n->m_iIndex+1,(int)(n->m_laChildren[z]+1));
		REC_DumpDot(a,(CClusterNode*)m_oaClusters[n->m_laChildren[z]]);
	}*/

	if (n->m_iChildren[0] != -1)
	{
		mfprintf(a,"  n%d -> n%d;\n",n->m_iIndex+1,(int)(n->m_iChildren[0]+1));
		REC_DumpDot(a,(CClusterNode*)m_oaClusters[n->m_iChildren[0]]);
		mfprintf(a,"  n%d -> n%d;\n",n->m_iIndex+1,(int)(n->m_iChildren[1]+1));
		REC_DumpDot(a,(CClusterNode*)m_oaClusters[n->m_iChildren[1]]);
	}
}


int CA_compareX(const void *arg1, const void *arg2 )
{
	if ((*((CClusterNode**)arg1))->m_fPosX > (*((CClusterNode**)arg2))->m_fPosX)
		return 1;
			else return -1;
}


void CClusterAnalysis::DumpAgr(const char *s)
{
	int z;
//	FILE *a;
	CxDoubleArray **fa;
//	CClusterNode *n, *n2;
	int *p;
	double /*d, */d2/*, d2l*/;
	bool b, c;
	int i;
	CGrace *g;
//	char buf[256];
	CxString buf;

	REC_SortX(1,0,m_pTop);

	qsort(&m_oaBaseClusters[0],m_oaBaseClusters.GetSize(),sizeof(CxObject*),CA_compareX);

	for (z=0;z<m_oaBaseClusters.GetSize();z++)
	{
//		mprintf("%d: %f\n",z+1,((CClusterNode*)m_oaBaseClusters[z])->m_fPosX);
		((CClusterNode*)m_oaBaseClusters[z])->m_fPosX = z+1;
	}

	REC_AvgX(m_pTop);

	try { g = new CGrace(); } catch(...) { g = NULL; }
	if (g == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	// Revised Version, better and clearer
	if (m_bPOVTrajectory && m_bPOVDiagram)
	{
		g->CurrentGraph()->m_bShowXAxis = false;
		g->CurrentGraph()->m_bShowFrame = false;
		g->CurrentGraph()->m_fViewportX1 = 0.10;
		g->CurrentGraph()->m_fViewportX2 = 1.30;
		g->CurrentGraph()->m_fViewportY1 = 0.6;
		g->CurrentGraph()->m_fViewportY2 = 0.95;
		g->CurrentGraph()->m_fYAxisBarWidth = 2.0;
		g->SetRangeX(-5.0,m_oaBaseClusters.GetSize()+1);
		g->SetRangeY(250.0,300.0);

		REC_DumpAgr(m_pTop,g,-1,0);
	} else
	{
		try { fa = new CxDoubleArray*[m_oaBaseClusters.GetSize()]; } catch(...) { fa = NULL; }
		if (fa == NULL) NewException((double)m_oaBaseClusters.GetSize()*sizeof(CxDoubleArray*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		try { p = new int[m_oaBaseClusters.GetSize()]; } catch(...) { p = NULL; }
		if (p == NULL) NewException((double)m_oaBaseClusters.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<m_oaBaseClusters.GetSize();z++)
		{
			try { fa[z] = new CxDoubleArray("CClusterAnalysis::DumpAgr():fa[z]"); } catch(...) { fa[z] = NULL; }
			if (fa[z] == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			p[z] = 0;
		}
	//	d = 0;

	//	DumpDot("C:\\Software\\test.dot");

	//	a = OpenFileWrite(s,true);

		
		for (z=0;z<m_oaBaseClusters.GetSize();z++)
		{
			g->AddDataset();
			g->SetSetLineColor(z,0,0,0);
//			sprintf(buf,"Molecule %d",((CClusterNode*)m_oaBaseClusters[z])->m_iIndex+1);
			buf.sprintf("Molecule %d",((CClusterNode*)m_oaBaseClusters[z])->m_iIndex+1);
			g->SetDatasetName(z,buf);
		}

	//	mfprintf(a,"0");
		for (z=0;z<m_oaBaseClusters.GetSize();z++)
		{
			fa[z]->Add(((CClusterNode*)m_oaBaseClusters[z])->m_fPosX);
			fa[z]->Add(0);
			REC_DumpPoints(fa[z],(CClusterNode*)m_oaBaseClusters[z]);
	//		mfprintf(a,";  %f",((CClusterNode*)m_oaBaseClusters[z])->m_fPosX);
			g->AddXYTupel(z,((CClusterNode*)m_oaBaseClusters[z])->m_fPosX,0);
		}
	//	mfprintf(a,"\n");

	/*	for (z=0;z<fa[0]->GetSize()/2;z++)
			mfprintf(a,"%f; %f\n",fa[0]->GetAt(z*2),fa[0]->GetAt(z*2+1));
		fclose(a);
		return;*/

	//	d2 = 0; // ?????

		i = 0;
		while (true)
		{
			i++;
	//		mprintf("Stage %d starting.\n",i);
	//		d2l = d2;
			d2 = 1.0e35;
			for (z=0;z<m_oaBaseClusters.GetSize();z++)
			{
	//			mprintf("    Cluster %d: %d Pairs, Pos=%d.\n",z+1,fa[z]->GetSize()/2,p[z]);
				if ((p[z]+1)*2+1 < fa[z]->GetSize())
				{
	//				mprintf("    %f < %f?\n",fa[z]->GetAt((p[z]+1)*2+1),d2);
					if (fa[z]->GetAt((p[z]+1)*2+1) < d2)
						d2 = fa[z]->GetAt((p[z]+1)*2+1);
				}
			}
	//		mprintf("d2=%f.\n",d2);
			if (d2 > 1.0e30)
				break;
	_next:
			b = false;
			c = false;
			for (z=0;z<m_oaBaseClusters.GetSize();z++)
			{
				if ((p[z]+1)*2+1 < fa[z]->GetSize())
				{
					if (fa[z]->GetAt((p[z]+1)*2+1) <= d2)
					{
						b = true;
						p[z]++;
					}
					if (!c)
					{
						c = true;
	//					mfprintf(a,"%f",d2);
					}
	//				mfprintf(a,";  %f",fa[z]->GetAt(p[z]*2));
					g->AddXYTupel(z,fa[z]->GetAt(p[z]*2),d2);
				}
			}
			if (!c)
			{
	//			mfprintf(a,"%f",d2+1.0);
				for (z=0;z<m_oaBaseClusters.GetSize();z++)
				{
					g->AddXYTupel(z,fa[0]->GetAt(p[0]*2),d2+100.0);
	//				mfprintf(a,";  %f",fa[0]->GetAt(p[0]*2));
				}
	//			mfprintf(a,"\n");
				break;
			}
	//		mfprintf(a,"\n");
			if (b)
				goto _next;
		}

		g->SetRangeX(0,m_oaBaseClusters.GetSize()+1);
		g->SetRangeY(0,m_fMaxDist);
	}

	g->MakeTicks();
	g->SetLabelX("Molecule");


	g->WriteAgr(s,true);

//	fclose(a);
}


void CClusterAnalysis::REC_DumpAgr(CClusterNode *n, CGrace *g, int ec, unsigned long color)
{
	int z;
	CClusterNode *n2;

	if (n->m_iMonomers == 1)
		return;

	if (ec == -1)
	{
		for (z=0;z<m_oaPOVExtendedCluster.GetSize();z++)
		{
			if (((CExtendedCluster*)m_oaPOVExtendedCluster[z])->m_iClusterIndex == n->m_iIndex)
			{
				ec = z;
				break;
			}
		}
	}

	if (ec != -1)
	{
		for (z=0;z<m_iaPOVSMCluster.GetSize();z++)
		{
			if (m_iaPOVSMCluster[z] == ec)
			{
//				mprintf("z=%d, ec=%d, color=%d.\n",z,ec,m_iaPOVSMColor[z]);
				color = m_iaPOVSMColor[z];
				break;
			}
		}
	} else color = 0x808080;

	if (color == 0xFFFF00)
		color = 0xD0D000;

	g->AddDataset();
	g->SetSetLineColor((unsigned char)((color&0xFF0000)>>16), (unsigned char)((color&0xFF00)>>8), (unsigned char)(color&0xFF));

	if (ec == -1)
		g->SetSetLineWidth(g->CurrentGraph()->m_oaDatasets.GetSize()-1,2.0);
			else g->SetSetLineWidth(g->CurrentGraph()->m_oaDatasets.GetSize()-1,2.5);

//	sprintf(buf,"Cluster %d",z+1);
//	g->SetDatasetName(buf);

	n2 = (CClusterNode*)m_oaClusters[n->m_iChildren[0]];
	if (n2->m_iMonomers == 1)
		g->AddXYTupel(n2->m_fPosX,0);
			else g->AddXYTupel(n2->m_fPosX,n2->m_fDistance);
	g->AddXYTupel(n2->m_fPosX,n->m_fDistance);

	n2 = (CClusterNode*)m_oaClusters[n->m_iChildren[1]];
	g->AddXYTupel(n2->m_fPosX,n->m_fDistance);
	if (n2->m_iMonomers == 1)
		g->AddXYTupel(n2->m_fPosX,0);
			else g->AddXYTupel(n2->m_fPosX,n2->m_fDistance);

	if (((CClusterNode*)m_oaClusters[n->m_iChildren[0]])->m_iMonomers > 1)
		REC_DumpAgr((CClusterNode*)m_oaClusters[n->m_iChildren[0]],g,ec,color);

	if (((CClusterNode*)m_oaClusters[n->m_iChildren[1]])->m_iMonomers > 1)
		REC_DumpAgr((CClusterNode*)m_oaClusters[n->m_iChildren[1]],g,ec,color);
}


void CClusterAnalysis::REC_SortX(double pos, int depth, CClusterNode *n)
{
	n->m_fPosX = pos;

	if (n->m_iChildren[0] == -1)
		return;

	if (((CClusterNode*)m_oaClusters[n->m_iChildren[0]])->m_fDistance > ((CClusterNode*)m_oaClusters[n->m_iChildren[1]])->m_fDistance)
	{
		REC_SortX(pos+mypow(0.5,depth+1),depth+1,(CClusterNode*)m_oaClusters[n->m_iChildren[0]]);
		REC_SortX(pos-mypow(0.5,depth+1),depth+1,(CClusterNode*)m_oaClusters[n->m_iChildren[1]]);
	} else
	{
		REC_SortX(pos-mypow(0.5,depth+1),depth+1,(CClusterNode*)m_oaClusters[n->m_iChildren[0]]);
		REC_SortX(pos+mypow(0.5,depth+1),depth+1,(CClusterNode*)m_oaClusters[n->m_iChildren[1]]);
	}
}


void CClusterAnalysis::REC_AvgX(CClusterNode *n)
{
	double d;
//	int z;

	if (n->m_iChildren[0] == -1)
		return;

	d = 0;

/*	for (z=0;z<n->m_laChildren.GetSize();z++)
	{
		REC_AvgX((CClusterNode*)m_oaClusters[n->m_laChildren[z]]);
		d += ((CClusterNode*)m_oaClusters[n->m_laChildren[z]])->m_fPosX;
	}*/

	REC_AvgX((CClusterNode*)m_oaClusters[n->m_iChildren[0]]);
	d += ((CClusterNode*)m_oaClusters[n->m_iChildren[0]])->m_fPosX;
	REC_AvgX((CClusterNode*)m_oaClusters[n->m_iChildren[1]]);
	d += ((CClusterNode*)m_oaClusters[n->m_iChildren[1]])->m_fPosX;

	d /= 2.0; //n->m_laChildren.GetSize();
//	mprintf("PosX=%f.\n",d);
	n->m_fPosX = d;
}


void CClusterAnalysis::REC_DumpPoints(CxDoubleArray *fa, CClusterNode *n)
{
	CClusterNode *n2;

//	if (n->m_iParent != -1)
	if (n->m_pParent != NULL)
	{
//		n2 = (CClusterNode*)m_oaClusters[n->m_iParent];
		n2 = n->m_pParent;
		fa->Add(n->m_fPosX);
		fa->Add(n2->m_fDistance);
		fa->Add(n2->m_fPosX);
		fa->Add(n2->m_fDistance);
		REC_DumpPoints(fa,n2);
	}
}


void CClusterAnalysis::BinDistances(CTimeStep *ts)
{
	int z, m, c, z2, i, z3, z0;
	double tf, mi, ma, ce, tm;
	CClusterNode *n;
	CExtendedCluster *ec;
	bool tb;
//	char buf[256];
	CxString buf;
	FILE *a;

	tf = (double)m_iCounter * g_iStride * g_fTimestepLength / 1000.0;

	ce = 0;
	c = 0;
	mi = 1e30;
	ma = 0;
	i = 0;

	if (m_bPOVTrajectory)
	{
		for (z=0;z<g_oaSingleMolecules.GetSize();z++)
		{
			m_iaPOVSMCluster[z] = -1;
			m_iaPOVSMColor[z] = -1;
		}
	}

	for (z=0;z<m_oaClusters.GetSize();z++)
	{
		n = (CClusterNode*)m_oaClusters[z];

		if (n->m_fDistance > 1.0E10)
		{
			m_pClusterSizeDF->AddToBinInt_Fast(0);
			m_pClusterDistributionDF->AddToBinInt_Fast(0,n->m_pParent->m_fDistance);
			if (m_b2DPlots)
			{
				m_pClusterSize2D->AddToBin(tf,1.0);
				m_pClusterDistribution2D->AddToBin(tf,0.0,n->m_pParent->m_fDistance);
			}
			continue;
		}

		m = n->m_iMonomers-1;
		m_pClusterSizeDF->AddToBinInt_Fast(m);

		if (n->m_pParent != NULL)
			m_pClusterDistributionDF->AddToBinInt_Fast(m,n->m_pParent->m_fDistance-n->m_fDistance);
				else m_pClusterDistributionDF->AddToBinInt_Fast(m,m_fMaxDist-n->m_fDistance);

		m_pClusterDistanceDF->AddToBin_Fast(n->m_fDistance);

		if (m_bHetMeasure)
		{
			if (n->m_fDistance < mi)
				mi = n->m_fDistance;
			if (n->m_fDistance > ma)
				ma = n->m_fDistance;

			ce += n->m_fDistance;
			c++;
		}

		if (m_b2DPlots)
		{
			m_pClusterSize2D->AddToBin(tf,sqrt(m+1)*sqrt(m_iMonomers));
			if (n->m_pParent != NULL)
				m_pClusterDistribution2D->AddToBin(tf,sqrt(m+1)*sqrt(m_iMonomers),n->m_pParent->m_fDistance-n->m_fDistance);
					else m_pClusterDistribution2D->AddToBin(tf,sqrt(m+1)*sqrt(m_iMonomers),m_fMaxDist-n->m_fDistance);
			m_pClusterDistance2D->AddToBin(tf,n->m_fDistance);
		}

		if (m_bDiffPlot)
		{
			if (n->m_pParent != NULL)
				m_pClusterDiff2D->AddToBin_IntX(n->m_iMonomers-1,n->m_pParent->m_fDistance - n->m_fDistance,1.0);
		}

		if (m_bClusterTopo)
		{
			if ((n->m_iMonomers > 1) && (n->m_iMonomers <= m_iClusterTopoMax))
			{
				ec = new CExtendedCluster();
				ec->m_iCount = 1;
				ec->m_iCountUndir = 1;
				ec->m_fSignificance = n->m_pParent->m_fDistance - n->m_fDistance;
				ec->m_fSignificanceUndir =  n->m_pParent->m_fDistance - n->m_fDistance;
				ec->CreateFrom(n,ts,this,true);

				tb = false;

				if (m_bClusterTopoDrawDirected || m_bClusterTopoDrawAtom)
				{
					for (z2=0;z2<ec->m_oaBonds.GetSize();z2++)
						if (((CxIntArray*)ec->m_oaBonds[z2])->GetSize() == 0)
							goto _skip1;
					ec->BuildAtomCodes();
					if (AddToClusterTopo(ec,&i))
						tb = true;

					if (m_bClusterTopoTrajectories)
					{
#ifdef TARGET_WINDOWS
//						sprintf(buf,"ClusterTopo\\xtraj_directed_%02d_%04d.xyz",n->m_iMonomers,i+1);
						buf.sprintf("ClusterTopo\\xtraj_directed_%02d_%04d.xyz",n->m_iMonomers,i+1);
#else
//						sprintf(buf,"ClusterTopo/xtraj_directed_%02d_%04d.xyz",n->m_iMonomers,i+1);
						buf.sprintf("ClusterTopo/xtraj_directed_%02d_%04d.xyz",n->m_iMonomers,i+1);
#endif
						a = fopen(buf,"a+t");
						WriteClusterXYZ(n,ts,a);
						fclose(a);
					}
_skip1:;
				}

				if (m_bClusterTopoDrawUndirected)
				{
					for (z2=0;z2<ec->m_oaRelBondsUndir.GetSize();z2++)
					{
						if (((CxIntArray*)ec->m_oaRelBondsUndir[z2])->GetSize() == 0)
						{
//								mprintf("Skip A.\n");
							goto _skip2;
						}
					}
					ec->BuildAtomCodesUndir();
					if (AddToClusterTopoUndir(ec,&i))
						tb = true;

					if (m_bClusterTopoTrajectories)
					{
#ifdef TARGET_WINDOWS
//						sprintf(buf,"ClusterTopo\\xtraj_undir_%02d_%04d.xyz",n->m_iMonomers,i+1);
						buf.sprintf("ClusterTopo\\xtraj_undir_%02d_%04d.xyz",n->m_iMonomers,i+1);
#else
//						sprintf(buf,"ClusterTopo/xtraj_undir_%02d_%04d.xyz",n->m_iMonomers,i+1);
						buf.sprintf("ClusterTopo/xtraj_undir_%02d_%04d.xyz",n->m_iMonomers,i+1);
#endif
						a = fopen(buf,"a+t");
						WriteClusterXYZ(n,ts,a);
						fclose(a);
					}
_skip2:;
				}

				if (!tb)
					delete ec;
			}
		}

		if (m_bCutClusters)
		{
			if (n->m_iMonomers <= m_iCutClusterMaxSize)
			{
				if (n->m_pParent != NULL)
				{
					if (n->m_pParent->m_fDistance - n->m_fDistance >= m_pCutClusterDifference[n->m_iMonomers-2])
					{
						m_pCutClusterCounter[n->m_iMonomers-2]++;

						WriteClusterXYZ(n,ts,m_pCutClusterFiles[n->m_iMonomers-2]);
					}
				}
			}
		}
	}

	if (m_bPOVTrajectory)
	{
		for (z0=m_iPOVMonomersMax;z0>=m_iPOVMonomersMin;z0--)
		{
			for (z=0;z<m_oaClusters.GetSize();z++)
			{
_beg:
				n = (CClusterNode*)m_oaClusters[z];

				if (n->m_fDistance > 1.0E10)
					continue;

				if (n->m_iMonomers == z0)
				{
					ec = new CExtendedCluster();
					ec->CreateFrom(n,ts,this,false);

					for (z2=0;z2<ec->m_iaMonomerSM.GetSize();z2++)
					{
						if (m_iaPOVSMCluster[ec->m_iaMonomerSM[z2]] != -1)
						{
							z++;
							goto _beg;
						}
					}

					if (IsConnected(ec))
					{
/*						int colstat[10];

						for (z2=0;z2<10;z2++)
							colstat[z2] = 0;

						for (z2=0;z2<ec->m_iaMonomerSM.GetSize();z2++)
						{
							for (z3=1;z3<100;z3++)
							{
								if (((CxIntArray*)m_oaPOVSMColorHistory[ec->m_iaMonomerSM[z2]])->GetAt(z3) != -1)
									colstat[((CxIntArray*)m_oaPOVSMColorHistory[ec->m_iaMonomerSM[z2]])->GetAt(z3)]++;
							}
						}

						z3 = 0;
						i = -1;

						for (z2=0;z2<10;z2++)
						{
							if (colstat[z2] > z3)
							{
								z3 = colstat[z2];
								i = z2;
							}
						}*/

						i = -1;
						for (z3=1;z3<100;z3++)
						{
							if (((CxIntArray*)m_oaPOVSMColorHistory[ec->m_iaMonomerSM[0]])->GetAt(z3) == -1)
								continue;

							for (z2=1;z2<ec->m_iaMonomerSM.GetSize();z2++)
								if (((CxIntArray*)m_oaPOVSMColorHistory[ec->m_iaMonomerSM[z2]])->GetAt(z3) != ((CxIntArray*)m_oaPOVSMColorHistory[ec->m_iaMonomerSM[0]])->GetAt(z3))
									goto _not;
							i = ((CxIntArray*)m_oaPOVSMColorHistory[ec->m_iaMonomerSM[0]])->GetAt(z3);
							break;
_not:;
						}

						if (i == -1)
						{
							i = m_iPOVColCount;
							m_iPOVColCount++;
							if (m_iPOVColCount > 9)
								m_iPOVColCount = 0;
						}

						for (z2=0;z2<ec->m_iaMonomerSM.GetSize();z2++)
							((CxIntArray*)m_oaPOVSMColorHistory[ec->m_iaMonomerSM[z2]])->SetAt(0,i);

						for (z2=0;z2<n->m_iaAtoms.GetSize();z2++)
						{
							m_iaPOVSMCluster[g_laAtomSMIndex[m_iaAtomList[n->m_iaAtoms[z2]]]] = m_oaPOVExtendedCluster.GetSize();
							m_iaPOVSMColor[g_laAtomSMIndex[m_iaAtomList[n->m_iaAtoms[z2]]]] = POVColorFunction(i);
						}

						m_oaPOVExtendedCluster.Add(ec);

//						mprintf("Put");
					} else delete ec;
				}
			}
		}
	}


	if (m_bHetMeasure)
	{
		ce /= c;

		tm = 0;
		for (z=0;z<m_oaClusters.GetSize();z++)
		{
			n = (CClusterNode*)m_oaClusters[z];

			if (n->m_fDistance > 1.0E10)
				continue;

			tm += pow2(n->m_fDistance - ce);
		}
		tm /= c;
		tm = sqrt(tm) / ce;

		m_faMinDist.Add(mi);
		m_faMaxDist.Add(ma);
		m_faHetNumber.Add((ma-mi)/(ma+mi));
		m_faCenter.Add(ce);
		m_faIntHetNumber.Add(tm);

		tm = 0;
		for (z=0;z<m_oaClusters.GetSize();z++)
		{
			n = (CClusterNode*)m_oaClusters[z];

			if (n->m_fDistance > 1.0E10)
				continue;

			tm += fabs(n->m_fDistance - ce);
		}
		tm /= c;
		tm = tm / ce;

		m_faIntHetNumber2.Add(tm);
	}
}


void CClusterAnalysis::Parse()
{
	int ti, ti2, z, z2, z3, z4, z5;
//	char buf[256], buf2[256];
	CxString buf, buf2;
	CAtomGroup *ag;
	CSingleMolecule *sm;
	CMolecule *m;

	mprintf(WHITE,"\n");
	mprintf(WHITE,"    #######################################################################\n");
	mprintf(WHITE,"    ####                 Hierarchical Cluster Analysis                 ####\n");
	mprintf(WHITE,"    #######################################################################\n");
	mprintf(WHITE,"    ####   Original implementation:                                    ####\n");
	mprintf(WHITE,"    ####     (c) Martin Brehm, 2013 (see PhD thesis)                   ####\n");
	mprintf(WHITE,"    ####   Present version:                                            ####\n");
	mprintf(WHITE,"    ####     (c) Martin Brehm, Tom Froembgen, Barbara Kirchner, 2022   ####\n");
	mprintf(WHITE,"    #######################################################################\n");
	mprintf(WHITE,"\n");

_canextmol:
	mprintf("\n");
//	sprintf(buf,"    Which molecule type to use for cluster analysis (");
	buf.sprintf("    Which molecule type to use for cluster analysis (");
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
//		sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//		strcat(buf,buf2);
		buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
		buf.strcat(buf2);

		if (z < g_oaMolecules.GetSize()-1)
//			strcat(buf,", ");
			buf.strcat(", ");
	}
//	strcat(buf,")? ");
	buf.strcat(")? ");

	ti = AskRangeInteger_ND("%s",1,g_oaMolecules.GetSize(),(const char*)buf) - 1;

	for (z=0;z<m_oaAtomGroups.GetSize();z++)
	{
		if (((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_iIndex == ti)
		{
			eprintf("This molecule type is already used for cluster analysis.\n");
			goto _canextmol;
		}
	}

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
_caagain:
	AskString("    Which atoms to take into account from %s (*=all)? [*] ",&buf,"*",((CMolecule*)g_oaMolecules[ti])->m_sName);
	if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[ti],buf))
		goto _caagain;

	m_oaAtomGroups.Add(ag);

	mprintf("\n");
	if (AskYesNo("    Add another molecule type to the cluster analysis (y/n)? [no] ",false))
		goto _canextmol;

	for (z=0;z<m_oaAtomGroups.GetSize();z++)
	{
		ag = (CAtomGroup*)m_oaAtomGroups[z];
		m = ag->m_pMolecule;
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			for (z3=0;z3<ag->m_oaAtoms.GetSize();z3++)
			{
				for (z4=0;z4<((CxIntArray*)ag->m_oaAtoms[z3])->GetSize();z4++)
				{
					m_iaAtomList.Add(((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z3]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z3])->GetAt(z4)));
				}
			}
		}
	}
	m_iAtomListSize = m_iaAtomList.GetSize();

	if (g_bAdvanced2) {
		mprintf("\n");
		m_bCOMTrick = AskYesNo("    Use \"CoM Trick\" (y/n)? [no] ",false);
	} else
		m_bCOMTrick = false;


	mprintf("\n");
	mprintf("    In the following, you can define which intermolecular distances are taken into account\n");
	mprintf("    for the cluster recognition. For example, if you want to investigate a hydrogen bond\n");
	mprintf("    network, answer \"y\" and enter the hydrogen bond donor and acceptor sites. If you answer\n");
	mprintf("    \"n\", the closest distance of any pair of atoms from two molecules will define the distance\n");
	mprintf("    between these molecules.\n");
	mprintf("\n");

	m_bSelectBonds = AskYesNo("    Consider only selected cluster bonds (y) or all bonds (n)? [no] ",false);

	if (m_bSelectBonds)
	{
		m_iaAtomBondFrom.SetSize(m_iAtomListSize);
		m_iaAtomBondTo.SetSize(m_iAtomListSize);
		for (z=0;z<m_iAtomListSize;z++)
		{
			m_iaAtomBondFrom[z] = 0;
			m_iaAtomBondTo[z] = 0;
		}

		mprintf("\n    * Definition of bond donor atoms\n");
		for (z=0;z<m_oaAtomGroups.GetSize();z++)
		{
_againx1:
			AskString("      Which atoms from %s to allow as bond donor (e.g. C1,C5-7, *=all)? [no atoms] ",&buf,"",((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_sName);
			if (strlen(buf) == 0)
				continue;
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			if (!ag->ParseAtoms(((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule,buf))
				goto _againx1;
			for (z2=0;z2<ag->m_baAtomType.GetSize();z2++)
			{
//				mprintf("z2=%d\n",z2);
				for (z3=0;z3<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z3++)
				{
//					mprintf("z3=%d\n",z3);
					for (z4=0;z4<((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_laSingleMolIndex.GetSize();z4++)
					{
//						mprintf("z4=%d\n",z4);
						ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_laSingleMolIndex[z4]])->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z3));
//						mprintf("ti=%d\n",ti);
						for (z5=0;z5<m_iAtomListSize;z5++)
						{
//							mprintf("z5=%d: %d vs %d\n",z5,m_iaAtomList[z5],ti);
							if (m_iaAtomList[z5] == ti)
							{
								ti = z5;
								goto _found1;
							}
						}
						eprintf("CClusterAnalysis::Parse(): Internal Error A.\n");
_found1:
						m_iaAtomBondFrom[ti] = 1;
					}
				}
			}
			delete ag;
		}

		mprintf("\n    * Definition of bond acceptor atoms\n");
		for (z=0;z<m_oaAtomGroups.GetSize();z++)
		{
_againx2:
			AskString("      Which atoms from %s to allow as bond acceptor (e.g. C1,C5-7, *=all)? [no atoms] ",&buf,"",((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_sName);
			if (strlen(buf) == 0)
				continue;
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			if (!ag->ParseAtoms(((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule,buf))
				goto _againx2;
			for (z2=0;z2<ag->m_baAtomType.GetSize();z2++)
			{
				for (z3=0;z3<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z3++)
				{
					for (z4=0;z4<((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_laSingleMolIndex.GetSize();z4++)
					{
						ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_laSingleMolIndex[z4]])->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z3));
						for (z5=0;z5<m_iAtomListSize;z5++)
						{
							if (m_iaAtomList[z5] == ti)
							{
								ti = z5;
								goto _found2;
							}
						}
						eprintf("CClusterAnalysis::Parse(): Internal Error B.\n");
_found2:
						m_iaAtomBondTo[ti] = 1;
					}
				}
			}
			delete ag;
		}

/*		mprintf("\nBond donors:\n");
		for (z=0;z<m_iAtomListSize;z++)
			mprintf("    %d: %d\n",z,m_iaAtomBondFrom[z]);
		mprintf("\nBond acceptors:\n");
		for (z=0;z<m_iAtomListSize;z++)
			mprintf("    %d: %d\n",z,m_iaAtomBondTo[z]);*/

		mprintf("\n");
	}

	if (g_bAdvanced2)
		m_bIdealGas = AskYesNo("    Use \"ideal gas mode\" (replace all coordinates by random numbers) (y/n)? [no] ",false);
	else
		m_bIdealGas = false;

	if (m_bIdealGas)
		m_pRandom = new CRandom();

	m_fMaxDist = AskFloat("    Enter max. distance for cluster analysis distribution (in pm): [%d] ",(double)HalfBoxSq3(),HalfBoxSq3());
	m_iRes = AskUnsignedInteger("    Enter resolution of the cluster distance distribution: [%.0f] ",(int)floor(m_fMaxDist),floor(m_fMaxDist));

	if (g_bAdvanced2) {
		mprintf("\n");
		m_fCSDFThresh = AskFloat("    Which relative threshold to use for the modified cluster significance DF? [0.01] ",0.01);
		mprintf("\n");
	} else {
		mprintf("\n");
		mprintf("    Using a relative threshold of 0.01 for the modified cluster significance DF.\n");
		mprintf("    Switch on the advanced mode to modify this parameter.\n\n");
		m_fCSDFThresh = 0.01;
	}

	mprintf("    The following \"standard\" results from the cluster analyis are always written:\n");
	mprintf("      (*) Cluster Distance Distribution Function\n");
	mprintf("      (*) Cluster Size Distribution Function\n");
	mprintf("      (*) Cluster Count Distribution Function\n");
	mprintf("      (*) Polymer Distribution Function (cumulative and non-cumulative)\n");
	mprintf("      (*) Cluster Significance Distribution Function (both original and modified)\n");

	mprintf("\n");

	mprintf("    In addition, a dendrogram (\"cluster_diagram.agr\") is written for the first processed frame.\n");

	mprintf("\n");

	if (!g_bAdvanced2) {

		mprintf("    Some additional functions are available in the advanced mode.\n");

		m_bHetMeasure = false;
		m_b2DPlots = false;
		m_bDiffPlot = false;
		m_bCutClusters = false;
		m_bAnim = false;
		m_bClusterTopo = false;
		m_bPOVTrajectory = false;

	} else {

		mprintf("    The following questions refer to additional analyses. If you are unsure, leave them switched off.\n\n");

		m_bHetMeasure = AskYesNo("    Write out temporal development of heterogeneity measures (y/n)? [no] ",false);
		if (m_bHetMeasure)
		{
			if (g_fTimestepLength == 0)
			{
				g_fTimestepLength = AskFloat("    Enter the physical time distance between successive trajectory frames in fs: [0.5] ",0.5);
				mprintf("\n");
			}
		}

		m_b2DPlots = AskYesNo("    Create 2D plots with temporal development (y/n)? [no] ",false);
		if (m_b2DPlots)
		{
			if (g_fTimestepLength == 0)
				g_fTimestepLength = AskFloat("      Enter the physical time distance between successive trajectory frames in fs: [0.5] ",0.5);
			m_i2DResX = AskUnsignedInteger("      Which resolution to use for the temporal axis of the 2D plots? [150] ",150);
			m_i2DResY = AskUnsignedInteger("      Which resolution to use for the ordinate axis of the 2D plots? [150] ",150);
			m_i2DStride = AskUnsignedInteger("      Use each n-th time step for temporal axis (stride)? [1] ",1);
			mprintf("\n");
		}

		m_bDiffPlot = AskYesNo("    Create Cluster Distance Difference plot (y/n)? [no] ",false);
		if (m_bDiffPlot)
		{
			m_fDiffInterval = AskFloat("      Enter max. value of difference inverval (in pm): [200] ",200.0);
			mprintf("\n");
		}

		m_bCutClusters = AskYesNo("    Export cluster structures from cluster analysis (y/n)? [no] ",false);
		if (m_bCutClusters)
		{
			m_iCutClusterMaxSize = AskUnsignedInteger("      Enter the largest cluster size to export: [3] ",3);

			try { m_pCutClusterDifference = new double[m_iCutClusterMaxSize-1]; } catch(...) { m_pCutClusterDifference = NULL; }
			if (m_pCutClusterDifference == NULL) NewException((double)(m_iCutClusterMaxSize-1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			try { m_pCutClusterFiles = new FILE*[m_iCutClusterMaxSize-1]; } catch(...) { m_pCutClusterFiles = NULL; }
			if (m_pCutClusterFiles == NULL) NewException((double)(m_iCutClusterMaxSize-1)*sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			try { m_pCutClusterCounter = new int[m_iCutClusterMaxSize-1]; } catch(...) { m_pCutClusterCounter = NULL; }
			if (m_pCutClusterCounter == NULL) NewException((double)(m_iCutClusterMaxSize-1)*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z=0;z<m_iCutClusterMaxSize-1;z++)
			{
				buf.sprintf("cluster_%dmer.xyz",z+2);
				m_pCutClusterFiles[z] = OpenFileWrite(buf,true);
				m_pCutClusterCounter[z] = 0;
			}

			for (z=0;z<m_iCutClusterMaxSize-1;z++)
				m_pCutClusterDifference[z] = AskFloat("      Enter required distance difference for %d-mer (in pm): [50] ",50.0,z+2);

			mprintf("\n");
		}

		m_bAnim = AskYesNo("    Create Cluster Diagram temporal animation (y/n)? [no] ",false);
		if (m_bAnim)
		{
			m_iResX = AskUnsignedInteger("    Enter width (in pixel) of the animation images: [640] ",640);
			m_iResY = AskUnsignedInteger("    Enter height (in pixel) of the animation images: [480] ",480);
			mprintf("\n");
		}

		m_bClusterTopo = AskYesNo("    Perform Cluster Topology Analysis (y/n)? [no] ",false);

		if (m_bClusterTopo)
		{
			mprintf("\n");
			m_iClusterTopoMax = AskUnsignedInteger("    Analyze Topology of clusters up to how many monomers? [8] ",8);


			m_bClusterTopoAllBonds = AskYesNo("    Consider only hydrogen bonds (n), or allow also dispersive bonds (y)? [no] ",false);

			m_bUseBondCutoffs = AskYesNo("    Use bond cutoff values (allow for ring clusters, etc.) (y/n)? [yes] ",true);

			if (m_bUseBondCutoffs)
			{
				m_fClusterTopoHBThres = AskFloat("    Enter hydrogen bond threshold distance (pm): [200] ",200.0);

				if (m_bClusterTopoAllBonds)
					m_fClusterTopoDispThres = AskFloat("    Enter dispersive bond threshold distance (pm): [200] ",200.0);
			} else
			{
				m_fClusterTopoHBThres = 0;
				m_fClusterTopoDispThres = 0;
			}

			m_bClusterTopoDrawAtom = AskYesNo("    Draw atom-level graphs (y/n)? [no] ",false);
			m_bClusterTopoDrawDirected = AskYesNo("    Draw directed graphs (y/n)? [yes] ",true);
			m_bClusterTopoDrawUndirected = AskYesNo("    Draw undirected graphs (y/n)? [yes] ",true);

			if (m_oaAtomGroups.GetSize() > 1)
			{
				mprintf("\n    More than 1 type of molecule used for cluster analysis.\n\n");
				for (z=0;z<m_oaAtomGroups.GetSize();z++)
				{
					mprintf("    * Color for molecule %s\n",((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_sName);

					if (m_oaAtomGroups.GetSize() > 1)
					{
						switch(z)
						{
							case 0: ti2 = 255; break; // Red
							case 1: ti2 = 128; break; // Blue
							case 2: ti2 = 128; break; // Green
							case 3: ti2 = 255; break; // Yellow
							default: ti2 = 211;
						}
					} else ti2 = 211;
					ti = AskUnsignedInteger("      Enter red component (0-255): [%d] ",ti2,ti2);
					m_iaClusterTopoMolColor.Add(ti);

					if (m_oaAtomGroups.GetSize() > 1)
					{
						switch(z)
						{
							case 0: ti2 = 128; break; // Red
							case 1: ti2 = 128; break; // Blue
							case 2: ti2 = 255; break; // Green
							case 3: ti2 = 255; break; // Yellow
							default: ti2 = 211;
						}
					} else ti2 = 211;
					ti = AskUnsignedInteger("      Enter green component (0-255): [%d] ",ti2,ti2);
					m_iaClusterTopoMolColor.Add(ti);

					if (m_oaAtomGroups.GetSize() > 1)
					{
						switch(z)
						{
							case 0: ti2 = 128; break; // Red
							case 1: ti2 = 255; break; // Blue
							case 2: ti2 = 128; break; // Green
							case 3: ti2 = 0; break; // Yellow
							default: ti2 = 211;
						}
					} else ti2 = 211;
					ti = AskUnsignedInteger("      Enter blue component (0-255): [%d] ",ti2,ti2);
					m_iaClusterTopoMolColor.Add(ti);

				//	if (z < m_oaAtomGroups.GetSize() -1)
						mprintf("\n");
				}
			} else
			{
				m_iaClusterTopoMolColor.Add(211);
				m_iaClusterTopoMolColor.Add(211);
				m_iaClusterTopoMolColor.Add(211);
			}
			m_iClusterTopoTries = AskUnsignedInteger("    How many tries for graph drawing? [10] ",10);
			m_bClusterTopoTrajectories = AskYesNo("    Write trajectory for each cluster topology (y/n)? [no] ",false);
			mprintf("\n");
		}

		m_bPOVTrajectory = AskYesNo("    Render POV-Ray trajectory of whole box indicating clusters (y/n)? [no] ",false);
		if (m_bPOVTrajectory)
		{
			mprintf("\n");
			m_bClusterTopoDrawDirected = true;
			m_bClusterTopoDrawUndirected = true;
			m_bClusterTopoAllBonds = AskYesNo("    Consider only hydrogen bonds (n), or allow also dispersive bonds (y)? [no] ",false);
			m_fClusterTopoHBThres = AskFloat("    Enter hydrogen bond threshold distance (pm): [200] ",200.0);
			if (m_bClusterTopoAllBonds)
				m_fClusterTopoDispThres = AskFloat("    Enter dispersive bond threshold distance (pm): [200] ",200.0);
			m_fPOVAngleIncrement = AskFloat("    Rotate camera (degree per time step)? [0.3] ",0.3f) * Pi / 180.0;
			m_iPOVMonomersMin = AskUnsignedInteger("    Minimal monomer number for emphasized clusters? [5] ",5);
			m_iPOVMonomersMax = AskUnsignedInteger("    Maximal monomer number for emphasized clusters? [5] ",5);
			if (m_bAnim)
				m_bPOVDiagram = AskYesNo("    Colorize Cluster Diagram according to POV files (y/n)? [yes] ",true);
			mprintf("\n");
		}
	}

	mprintf("\n");
	if (g_bAdvanced2)
		m_bDistCache = AskYesNo("    Use distance caching (y/n)? [yes] ",true);
	else {
		mprintf("    Using distance caching. Switch on the advanced mode to modify this setting.\n");
		m_bDistCache = true;
	}

	mprintf(WHITE,"\n<<< End of Cluster Analysis <<<\n\n");
}


void CClusterAnalysis::Create()
{
	int z, z2;//, z3, z4;
//	char buf[64];
	CxString buf;
//	CAtomGroup *ag;
//	CSingleMolecule *sm;
//	CMolecule *m;

	if (m_bAnim)
	{
		(void)!system("mkdir cluster_agr");

#ifdef TARGET_WINDOWS
		m_fAnim = OpenFileWrite("cluster_agr\\render_cluster_anim",true);
#else
		m_fAnim = OpenFileWrite("cluster_agr/render_cluster_anim",true);
#endif
	}

	m_iCounter = 0;
	m_iMonomers = 0;

	for (z=0;z<m_oaAtomGroups.GetSize();z++)
		m_iMonomers += ((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_laSingleMolIndex.GetSize();

	try { m_poaClusterPoly = new CxObArray[m_iMonomers]; } catch(...) { m_poaClusterPoly = NULL; }
	if (m_poaClusterPoly == NULL) NewException((double)m_iMonomers*sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (m_bDistCache)
	{
		try { m_pAtomIndex = new int[g_iGesVirtAtomCount]; } catch(...) { m_pAtomIndex = NULL; }
		if (m_pAtomIndex == NULL) NewException((double)g_iGesVirtAtomCount*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		for (z=0;z<g_iGesVirtAtomCount;z++)
			m_pAtomIndex[z] = -1;

		for (z=0;z<m_iaAtomList.GetSize();z++)
			m_pAtomIndex[m_iaAtomList[z]] = z;

		try { m_pDistCache = new double[m_iaAtomList.GetSize()*m_iaAtomList.GetSize()]; } catch(...) { m_pDistCache = NULL; }
		if (m_pDistCache == NULL) NewException((double)m_iaAtomList.GetSize()*m_iaAtomList.GetSize()*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}	

	try { m_pClusterDistributionDF = new CDF(); } catch(...) { m_pClusterDistributionDF = NULL; }
	if (m_pClusterDistributionDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pClusterDistributionDF->m_iResolution = m_iMonomers;
	m_pClusterDistributionDF->m_fMinVal = 1;
	m_pClusterDistributionDF->m_fMaxVal = m_iMonomers+1.0;
	m_pClusterDistributionDF->SetLabelX("Cluster Size");
	m_pClusterDistributionDF->SetLabelY("Percentage");
	m_pClusterDistributionDF->m_bLeft = true;
	m_pClusterDistributionDF->Create();

	try { m_pClusterDistributionDF_Mod = new CDF(); } catch(...) { m_pClusterDistributionDF_Mod = NULL; }
	if (m_pClusterDistributionDF_Mod == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pClusterDistributionDF_Mod->m_iResolution = m_iMonomers;
	m_pClusterDistributionDF_Mod->m_fMinVal = 1;
	m_pClusterDistributionDF_Mod->m_fMaxVal = m_iMonomers+1.0;
	m_pClusterDistributionDF_Mod->SetLabelX("Cluster Size");
	m_pClusterDistributionDF_Mod->SetLabelY("Percentage");
	m_pClusterDistributionDF_Mod->m_bLeft = true;
	m_pClusterDistributionDF_Mod->Create();

	if (m_b2DPlots)
	{
		try { m_pClusterDistribution2D = new C2DF(); } catch(...) { m_pClusterDistribution2D = NULL; }
		if (m_pClusterDistribution2D == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_pClusterDistribution2D->m_iRes[0] = m_i2DResX;
		m_pClusterDistribution2D->m_iRes[1] = m_i2DResY;
		m_pClusterDistribution2D->m_fMinVal[0] = 0;
		m_pClusterDistribution2D->m_fMaxVal[0] = (double)m_i2DResX*m_i2DStride*g_iStride*g_fTimestepLength/1000.0;
		m_pClusterDistribution2D->m_fMinVal[1] = 0;
		m_pClusterDistribution2D->m_fMaxVal[1] = m_iMonomers;
		m_pClusterDistribution2D->SetLabelX("Simulation Time / ps");
		m_pClusterDistribution2D->SetLabelY("Cluster Size");
		m_pClusterDistribution2D->m_bContourLines = false;
		m_pClusterDistribution2D->m_fPlotExp = 0.1;
		m_pClusterDistribution2D->Create();
	}

	if (m_bDiffPlot)
	{
		try { m_pClusterDiff2D = new C2DF(); } catch(...) { m_pClusterDiff2D = NULL; }
		if (m_pClusterDiff2D == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_pClusterDiff2D->m_iRes[0] = m_iMonomers;
		m_pClusterDiff2D->m_iRes[1] = 150;
		m_pClusterDiff2D->m_fMinVal[0] = 0;
		m_pClusterDiff2D->m_fMaxVal[0] = m_iMonomers;
		m_pClusterDiff2D->m_fMinVal[1] = 0;
		m_pClusterDiff2D->m_fMaxVal[1] = m_fDiffInterval;
		m_pClusterDiff2D->SetLabelX("Cluster Size");
		m_pClusterDiff2D->SetLabelY("Cluster Distance Difference / pm");
		m_pClusterDiff2D->m_bContourLines = false;
		m_pClusterDiff2D->m_fPlotExp = 0.1;
		m_pClusterDiff2D->Create();
	}

/*	try { m_pClusterDistribution2DF = new CDF(); } catch(...) { m_pClusterDistribution2DF = NULL; }
	if (m_pClusterDistribution2DF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pClusterDistribution2DF->m_iResolution = m_iMonomers;
	m_pClusterDistribution2DF->m_fMinVal = 1;
	m_pClusterDistribution2DF->m_fMaxVal = m_iMonomers*((m_iMonomers+1.0)/m_iMonomers);
	m_pClusterDistribution2DF->SetLabelX("Cluster size");
	m_pClusterDistribution2DF->SetLabelY("Percentage");
	m_pClusterDistribution2DF->m_bLeft = true;
	m_pClusterDistribution2DF->Create();*/

	try { m_pClusterSizeDF = new CDF(); } catch(...) { m_pClusterSizeDF = NULL; }
	if (m_pClusterSizeDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pClusterSizeDF->m_iResolution = m_iMonomers;
	m_pClusterSizeDF->m_fMinVal = 1;
	m_pClusterSizeDF->m_fMaxVal = m_iMonomers+1.0;
	m_pClusterSizeDF->SetLabelX("Cluster Size");
	m_pClusterSizeDF->SetLabelY("Percentage");
	m_pClusterSizeDF->m_bLeft = true;
	m_pClusterSizeDF->Create();

	if (m_b2DPlots)
	{
		try { m_pClusterSize2D = new C2DF(); } catch(...) { m_pClusterSize2D = NULL; }
		if (m_pClusterSize2D == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_pClusterSize2D->m_iRes[0] = m_i2DResX;
		m_pClusterSize2D->m_iRes[1] = m_i2DResY;
		m_pClusterSize2D->m_fMinVal[0] = 0;
		m_pClusterSize2D->m_fMaxVal[0] = (double)m_i2DResX*m_i2DStride*g_iStride*g_fTimestepLength/1000.0;
		m_pClusterSize2D->m_fMinVal[1] = 0;
		m_pClusterSize2D->m_fMaxVal[1] = m_iMonomers;
		m_pClusterSize2D->SetLabelX("Simulation Time / ps");
		m_pClusterSize2D->SetLabelY("Cluster Size");
		m_pClusterSize2D->m_bContourLines = false;
		m_pClusterSize2D->m_fPlotExp = 0.3;
		m_pClusterSize2D->Create();
	}

	try { m_pClusterDistanceDF = new CDF(); } catch(...) { m_pClusterDistanceDF = NULL; }
	if (m_pClusterDistanceDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pClusterDistanceDF->m_iResolution = m_iRes;
	m_pClusterDistanceDF->m_fMinVal = 0;
	m_pClusterDistanceDF->m_fMaxVal = m_fMaxDist;
	m_pClusterDistanceDF->SetLabelX("Cutoff Distance / pm");
	m_pClusterDistanceDF->SetLabelY("Occurrence");
	m_pClusterDistanceDF->Create();

	if (m_b2DPlots)
	{
		try { m_pClusterDistance2D = new C2DF(); } catch(...) { m_pClusterDistance2D = NULL; }
		if (m_pClusterDistance2D == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_pClusterDistance2D->m_iRes[0] = m_i2DResX;
		m_pClusterDistance2D->m_iRes[1] = m_i2DResY;
		m_pClusterDistance2D->m_fMinVal[0] = 0;
		m_pClusterDistance2D->m_fMaxVal[0] = (double)m_i2DResX*m_i2DStride*g_iStride*g_fTimestepLength/1000.0;
		m_pClusterDistance2D->m_fMinVal[1] = 0;
		m_pClusterDistance2D->m_fMaxVal[1] = m_fMaxDist;
		m_pClusterDistance2D->SetLabelX("Simulation Time / ps");
		m_pClusterDistance2D->SetLabelY("Cutoff Distance / pm");
		m_pClusterDistance2D->m_bContourLines = false;
		m_pClusterDistance2D->m_fPlotExp = 0.15;
		m_pClusterDistance2D->Create();
	}

	try { m_pPolymerDF = new CDF(); } catch(...) { m_pPolymerDF = NULL; }
	if (m_pPolymerDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pPolymerDF->m_iResolution = m_iRes;
	m_pPolymerDF->m_fMinVal = 0;
	m_pPolymerDF->m_fMaxVal = m_fMaxDist;
	m_pPolymerDF->SetLabelX("Cutoff Distance / pm");
	m_pPolymerDF->SetLabelY("Occurrence");
	m_pPolymerDF->CreateMulti(m_iMonomers);
	m_pPolymerDF->SetLabelMulti(0,"Monomer");

	for (z=1;z<m_iMonomers;z++)
	{
//		sprintf(buf,"%d-mer",z+1);
		buf.sprintf("%d-mer",z+1);
		m_pPolymerDF->SetLabelMulti(z,buf);
	}

	try { m_pClusterCountDF = new CDF(); } catch(...) { m_pClusterCountDF = NULL; }
	if (m_pClusterCountDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pClusterCountDF->m_iResolution = m_iRes;
	m_pClusterCountDF->m_fMinVal = 0;
	m_pClusterCountDF->m_fMaxVal = m_fMaxDist;
	m_pClusterCountDF->SetLabelX("Cutoff Distance / pm");
	m_pClusterCountDF->SetLabelY("Number of Clusters");
	m_pClusterCountDF->Create();

	if (m_b2DPlots)
	{
		try { m_pClusterCount2D = new C2DF(); } catch(...) { m_pClusterCount2D = NULL; }
		if (m_pClusterCount2D == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_pClusterCount2D->m_iRes[0] = m_i2DResX;
		m_pClusterCount2D->m_iRes[1] = m_i2DResY;
		m_pClusterCount2D->m_fMinVal[0] = 0;
		m_pClusterCount2D->m_fMaxVal[0] = (double)m_i2DResX*m_i2DStride*g_iStride*g_fTimestepLength/1000.0;
		m_pClusterCount2D->m_fMinVal[1] = 0;
		m_pClusterCount2D->m_fMaxVal[1] = m_fMaxDist;
		m_pClusterCount2D->SetLabelX("Simulation Time / ps");
		m_pClusterCount2D->SetLabelY("Cutoff Distance / pm");
		m_pClusterCount2D->m_bContourLines = false;
		m_pClusterCount2D->m_fPlotExp = 0.5;
		m_pClusterCount2D->Create();
	}

	if (m_bClusterTopo)
	{
		m_iaClusterTopoSMTemp.SetSize(g_oaSingleMolecules.GetSize());

		if (m_bClusterTopoDrawDirected || m_bClusterTopoDrawAtom)
		{
			m_oaClusterTopo.SetSize(m_iMonomers);
			for (z=0;z<m_iMonomers;z++)
				m_oaClusterTopo[z] = new CxObArray();
		}
		if (m_bClusterTopoDrawUndirected)
		{
			m_oaClusterTopoUndir.SetSize(m_iMonomers);
			for (z=0;z<m_iMonomers;z++)
				m_oaClusterTopoUndir[z] = new CxObArray();
		}
	}

	if (m_bPOVTrajectory)
	{
		(void)!system("mkdir POV");
		m_iaPOVSMCluster.SetSize(g_oaSingleMolecules.GetSize());
		m_iaPOVSMColor.SetSize(g_oaSingleMolecules.GetSize());
		for (z=0;z<g_oaSingleMolecules.GetSize();z++)
		{
			m_oaPOVSMColorHistory.Add(new CxIntArray());
			((CxIntArray*)m_oaPOVSMColorHistory[z])->SetSize(100);
			for (z2=0;z2<100;z2++)
				((CxIntArray*)m_oaPOVSMColorHistory[z])->SetAt(z2,-1);
		}

#ifdef TARGET_WINDOWS
		m_fPOVScript = OpenFileWrite("POV\\pov_render_script.bat",true);
#else
		m_fPOVScript = OpenFileWrite("POV/pov_render_script.bat",true);
#endif

		m_iPOVFrameCounter = 1;
		m_iPOVColCount = 0;
		m_fPOVAngle = 0;
	}
}


void CClusterAnalysis::Process(CTimeStep *ts)
{
	int z0, z, z2, z3, ti;
	CMolecule *m;
	CSingleMolecule *sm;
	CAtomGroup *ag;
//	char buf[256];
	CxString buf;

	CleanUp();

	if (m_bDistCache)
		BuildDistCache(ts);

	if (m_bPOVTrajectory)
	{
		ts->FoldMolecules();
		for (z=0;z<g_oaSingleMolecules.GetSize();z++)
			for (z2=99;z2>0;z2--)
				((CxIntArray*)m_oaPOVSMColorHistory[z])->SetAt(z2,((CxIntArray*)m_oaPOVSMColorHistory[z])->GetAt(z2-1));
	}

	for (z0=0;z0<m_oaAtomGroups.GetSize();z0++)
	{
		ag = (CAtomGroup*)m_oaAtomGroups[z0];
		m = (CMolecule*)g_oaMolecules[ag->m_pMolecule->m_iIndex];
		for (z=0;z<m->m_laSingleMolIndex.GetSize();z++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z]];
			AddCluster(ag->m_iAtomGes);
			for (z2=0;z2<ag->m_baAtomType.GetSize();z2++)
			{
				for (z3=0;z3<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z3++)
				{
					ti = ((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z3));
					if (m_bDistCache)
						AddParticle(m_pAtomIndex[ti]);
							else AddParticle(ts->m_vaCoords[ti][0],ts->m_vaCoords[ti][1],ts->m_vaCoords[ti][2]);
				}
			}
		}
	}

	BuildTree();
	BinDistances(ts);
	BinPoly();

	if (m_iCounter == 0)
	{
		DumpAgr("cluster_diagram.agr");
	}

	if (m_bAnim)
	{

#ifdef TARGET_WINDOWS
//		sprintf(buf,"cluster_agr\\cluster_anim_%05lu.agr",g_iSteps/g_iStride);
		buf.sprintf("cluster_agr\\cluster_anim_%05lu.agr",g_iSteps/g_iStride);
#else
//		sprintf(buf,"cluster_agr/cluster_anim_%05lu.agr",g_iSteps/g_iStride);
		buf.sprintf("cluster_agr/cluster_anim_%05lu.agr",g_iSteps/g_iStride);
#endif

		DumpAgr(buf);
		mfprintf(m_fAnim,"echo 'Printing frame %lu' ; \n",g_iSteps/g_iStride);
//		sprintf(buf,"cluster_anim_%05lu.agr",g_iSteps/g_iStride);
		buf.sprintf("cluster_anim_%05lu.agr",g_iSteps/g_iStride);
		mfprintf(m_fAnim,"xmgrace %s -batch gracebatch -nosafe -hardcopy ; \n",(const char*)buf);
		mfprintf(m_fAnim,"mv output.png frame%05lu.png ; \n",g_iSteps/g_iStride);
	}

	if (m_bPOVTrajectory)
	{

#ifdef TARGET_WINDOWS
//		sprintf(buf,"POV\\frame_%06d.pov",(int)g_iSteps);
		buf.sprintf("POV\\frame_%06d.pov",(int)g_iSteps);
#else
//		sprintf(buf,"POV/frame_%06d.pov",(int)g_iSteps);
		buf.sprintf("POV/frame_%06d.pov",(int)g_iSteps);
#endif

		RenderStepPOV(ts,buf);

		for (z=0;z<m_oaPOVExtendedCluster.GetSize();z++)
			delete (CExtendedCluster*)m_oaPOVExtendedCluster[z];
		m_oaPOVExtendedCluster.RemoveAll_KeepSize();

#ifdef TARGET_WINDOWS
		mfprintf(m_fPOVScript,"pvengine64 /NR /EXIT +a0.3 +w600 +h600 +fn -j0.0 frame_%06d.pov\n",m_iPOVFrameCounter);
#else
		mfprintf(m_fPOVScript,"povray +a0.3 +w600 +h600 +fn -j0.0 frame_%06d.pov\n",m_iPOVFrameCounter);
#endif

		fflush(m_fPOVScript);

		m_iPOVFrameCounter++;
		m_fPOVAngle += m_fPOVAngleIncrement;
	}

	m_iCounter++;
}


void CClusterAnalysis::CleanUp()
{
	int z;

	for (z=0;z<m_oaClusters.GetSize();z++)
		delete (CClusterNode*)m_oaClusters[z];
	m_oaClusters.RemoveAll_KeepSize();

	m_oaTopClusters.RemoveAll_KeepSize();

	m_oaBaseClusters.RemoveAll_KeepSize();

	for (z=0;z<m_iMonomers;z++)
		m_poaClusterPoly[z].RemoveAll_KeepSize();
}


void CClusterAnalysis::BinPoly()
{
	int z, z2, z3, l, u;
	double tf;
	CClusterNode *n;

	tf = (double)m_iCounter * g_iStride * g_fTimestepLength / 1000.0;

	for (z=0;z<m_iMonomers;z++)
	{
		for (z2=0;z2<m_poaClusterPoly[z].GetSize();z2++)
		{
			n = (CClusterNode*)(m_poaClusterPoly[z])[z2];

			if (n->m_iChildren[0] == -1)
				l = 0;
					else l = (int)(n->m_fDistance / m_fMaxDist * m_iRes);

			if (l >= m_iRes)
				l = m_iRes-1;

//			if (n->m_iParent == -1)
//				u = m_iRes-1;
//					else u = ((CClusterNode*)m_oaClusters[n->m_iParent])->m_fDistance / m_fMaxDist * m_iRes - 1;

			if (n->m_pParent == NULL)
				u = m_iRes-1;
					else u = (int)(n->m_pParent->m_fDistance / m_fMaxDist * m_iRes - 1);

			if (u >= m_iRes)
				u = m_iRes-1;

			for (z3=l;z3<=u;z3++)
			{
				m_pPolymerDF->AddToBin_Multi_Int_Fast(z,z3,z+1.0);
				m_pClusterCountDF->AddToBinInt_Fast(z3);
			}

			if (m_b2DPlots)
			{
				if (n->m_iChildren[0] == -1)
					l = 0;
						else l = (int)(n->m_fDistance / m_fMaxDist * m_pClusterCount2D->m_iRes[1]);

				if (l >= m_pClusterCount2D->m_iRes[1])
					l = m_pClusterCount2D->m_iRes[1]-1;

				if (n->m_pParent == NULL)
					u = m_pClusterCount2D->m_iRes[1]-1;
						else u = (int)(n->m_pParent->m_fDistance / m_fMaxDist * m_pClusterCount2D->m_iRes[1] - 1);

				if (u >= m_pClusterCount2D->m_iRes[1])
					u = m_pClusterCount2D->m_iRes[1]-1;

				for (z3=l;z3<=u;z3++)
					m_pClusterCount2D->AddToBin(tf,z3*m_fMaxDist/m_pClusterCount2D->m_iRes[1]);
			}
		}
	}
}


void CClusterAnalysis::BuildDistCache(CTimeStep *ts)
{
	int z, z2;
	CxDVector3 *a, *b;
	double t;

	if (m_bIdealGas)
	{
		for (z=0;z<m_iAtomListSize;z++)
		{
			t = m_pRandom->RandomUniform();
			ts->m_vaCoords[m_iaAtomList[z]][0] = t*g_fBoxX;
			ts->m_vaCoords[m_iaAtomList[z]][1] = m_pRandom->RandomUniform()*g_fBoxY;
			ts->m_vaCoords[m_iaAtomList[z]][2] = m_pRandom->RandomUniform()*g_fBoxZ;
		}
	}

	if (m_bCOMTrick)
	{
		for (z=0;z<m_iAtomListSize;z++)
		{
			a = &ts->m_vaCoords[MassCenter(m_iaAtomList[z])];
			for (z2=z+1;z2<m_iAtomListSize;z2++)
			{
				if (m_bSelectBonds)
				{
					if (((m_iaAtomBondFrom[z] != 0) && (m_iaAtomBondTo[z2] != 0)) || ((m_iaAtomBondFrom[z2] != 0) && (m_iaAtomBondTo[z] != 0)))
					{
						b = &ts->m_vaCoords[m_iaAtomList[z2]];
						t = Distance(a,b);
					} else t = 10000000.0;
				} else
				{
					b = &ts->m_vaCoords[MassCenter(m_iaAtomList[z2])];
					t = Distance(a,b);
				}
				m_pDistCache[z*m_iAtomListSize+z2] = t;
				m_pDistCache[z2*m_iAtomListSize+z] = t;
			}
		}
	}
	else // Normal distances
	{
		for (z=0;z<m_iAtomListSize;z++)
		{
			a = &ts->m_vaCoords[m_iaAtomList[z]];
			for (z2=z+1;z2<m_iAtomListSize;z2++)
			{
				if (m_bSelectBonds)
				{
					if (((m_iaAtomBondFrom[z] != 0) && (m_iaAtomBondTo[z2] != 0)) || ((m_iaAtomBondFrom[z2] != 0) && (m_iaAtomBondTo[z] != 0)))
					{
						b = &ts->m_vaCoords[m_iaAtomList[z2]];
						t = Distance(a,b);
					} else t = 10000000.0;
				} else
				{
					b = &ts->m_vaCoords[m_iaAtomList[z2]];
					t = Distance(a,b);
				}
				m_pDistCache[z*m_iAtomListSize+z2] = t;
				m_pDistCache[z2*m_iAtomListSize+z] = t;
			}
		}
	}
}


void CClusterAnalysis::BuildClusterDistribution()
{
	int z, z2;

	for (z=0;z<m_iMonomers;z++)
	{
//		for (z2=0;z2<m_iRes;z2++)
//			m_pClusterDistributionDF->m_pBin[z] += m_pPolymerDF->m_pMultiBin[z][z2];

//		m_pClusterDistributionDF->m_pBin[z] *= 100.0 / m_iRes / m_iMonomers / m_iCounter;

		m_pClusterDistributionDF->m_pBin[z] *= 100.0 * (z+1.0) / m_iMonomers / m_iCounter / m_fMaxDist;

		m_pClusterSizeDF->m_pBin[z] *= 100.0 * (z+1.0) / m_iMonomers / m_iCounter;
	}

	if (m_b2DPlots)
	{
		for (z=0;z<m_i2DResY;z++)
		{
			for (z2=0;z2<m_i2DResX;z2++)
			{
				m_pClusterDistribution2D->m_pBin[z*m_i2DResY+z2] *= 100.0 * (z+1.0)/m_i2DResY / m_fMaxDist;
				m_pClusterSize2D->m_pBin[z*m_i2DResY+z2] *= 100.0 * (z+1.0)/m_i2DResY;
			}
		}
		m_pClusterCount2D->MultiplyBin(1.0/m_i2DStride);
	}
}


void CClusterAnalysis::WriteOutput(const char *multibuf)
{
	int z, z2, z3, i, j, timin, timax;
	double tf, tfm;
	CExtendedCluster *ec;
//	char buf[256], buf2[256];
	CxString buf, buf2;
	FILE *a;

	if (m_bPOVTrajectory)
	{
		fclose(m_fPOVScript);

#ifdef TARGET_WINDOWS
		a = OpenFileWrite("POV\\video.avs",true);
#else
		a = OpenFileWrite("POV/video.avs",true);
#endif

		mfprintf(a,"video = ImageSource(\"frame_%c06d.png\", start=1, end=%d, fps=15, pixel_type = \"rgb24\" )\n",'%',m_iPOVFrameCounter-1);
		mfprintf(a,"return video\n");
		fclose(a);
	}

	if (m_bClusterTopo)
	{
		mprintf(WHITE,"  * Cluster Topology Analysis\n");
		(void)!system("mkdir ClusterTopo");

		if (m_bClusterTopoDrawDirected || m_bClusterTopoDrawAtom)
		{
			mprintf("    * Directed Graphs\n");
			a = OpenFileWrite("cluster_topologies_directed.csv",true);
			mfprintf(a,"n-mer;  id;  count;  countperc;  avgsigni / pm;  signisum / pm\n");
			for (z=0;z<m_oaClusterTopo.GetSize();z++)
			{
				if (((CxObArray*)m_oaClusterTopo[z])->GetSize() == 0)
					continue;
				
				for (z2=0;z2<((CxObArray*)m_oaClusterTopo[z])->GetSize();z2++)
					((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iIndex = z2;

				// Stack Sort
				for (z2=0;z2<((CxObArray*)m_oaClusterTopo[z])->GetSize();z2++)
				{
					j = 0;
					i = -1;
					for (z3=z2;z3<((CxObArray*)m_oaClusterTopo[z])->GetSize();z3++)
					{
						if (((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z3))->m_iCount > j)
						{
							j = ((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z3))->m_iCount;
							i = z3;
						}
					}
					ec = (CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2);
					((CxObArray*)m_oaClusterTopo[z])->SetAt(z2,((CxObArray*)m_oaClusterTopo[z])->GetAt(i));
					((CxObArray*)m_oaClusterTopo[z])->SetAt(i,ec);
				}

				i = 0;
				for (z2=0;z2<((CxObArray*)m_oaClusterTopo[z])->GetSize();z2++)
					i += ((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iCount;

				if (m_bClusterTopoTrajectories)
				{
					for (z2=0;z2<((CxObArray*)m_oaClusterTopo[z])->GetSize();z2++)
					{

#ifdef TARGET_WINDOWS
//						sprintf(buf,"ClusterTopo\\xtraj_directed_%02d_%04d.xyz",z+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iIndex+1);
//						sprintf(buf2,"ClusterTopo\\traj_directed_%02d_%04d.xyz",z+1,z2+1);
						buf.sprintf("ClusterTopo\\xtraj_directed_%02d_%04d.xyz",z+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iIndex+1);
						buf2.sprintf("ClusterTopo\\traj_directed_%02d_%04d.xyz",z+1,z2+1);
#else
//						sprintf(buf,"ClusterTopo/xtraj_directed_%02d_%04d.xyz",z+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iIndex+1);
//						sprintf(buf2,"ClusterTopo/traj_directed_%02d_%04d.xyz",z+1,z2+1);
						buf.sprintf("ClusterTopo/xtraj_directed_%02d_%04d.xyz",z+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iIndex+1);
						buf2.sprintf("ClusterTopo/traj_directed_%02d_%04d.xyz",z+1,z2+1);
#endif

						if (rename(buf,buf2) != 0)
							eprintf("Could not rename %s to %s.\n",(const char*)buf,(const char*)buf2);
					}
				}

				mprintf("        %d %d-mer topologies found. Rendering:",((CxObArray*)m_oaClusterTopo[z])->GetSize(),z+1);
				for (z2=0;z2<((CxObArray*)m_oaClusterTopo[z])->GetSize();z2++)
				{
					mprintf(".");
					if (m_bClusterTopoDrawDirected)
					{

#ifdef TARGET_WINDOWS
//						sprintf(buf,"ClusterTopo\\cluster_directed_%02d_%04d",z+1,z2+1);
						buf.sprintf("ClusterTopo\\cluster_directed_%02d_%04d",z+1,z2+1);
#else
//						sprintf(buf,"ClusterTopo/cluster_directed_%02d_%04d",z+1,z2+1);
						buf.sprintf("ClusterTopo/cluster_directed_%02d_%04d",z+1,z2+1);
#endif

						DumpSchematicDirectedDOT((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2),buf);
					}
					if (m_bClusterTopoDrawAtom)
					{

#ifdef TARGET_WINDOWS
//						sprintf(buf,"ClusterTopo\\cluster_atom_%02d_%04d",z+1,z2+1);
						buf.sprintf("ClusterTopo\\cluster_atom_%02d_%04d",z+1,z2+1);
#else
//						sprintf(buf,"ClusterTopo/cluster_atom_%02d_%04d",z+1,z2+1);
						buf.sprintf("ClusterTopo/cluster_atom_%02d_%04d",z+1,z2+1);
#endif

						DumpAtomDOT((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2),buf);
					}
					mfprintf(a,"%d;  %d;  %d;  %f;  %f;  %f\n",z+1,z2+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iCount,100.0*((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iCount/i,((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_fSignificance/((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_iCount,((CExtendedCluster*)((CxObArray*)m_oaClusterTopo[z])->GetAt(z2))->m_fSignificance);
					fflush(a);
			//		sprintf(buf2,"del %s",buf);
			//		(void)!system(buf2);
				}
				mfprintf(a,"\n");
				mprintf("\n");
			}
			fclose(a);
		}

		if (m_bClusterTopoDrawUndirected)
		{
			mprintf("    * Unirected Graphs\n");
			a = OpenFileWrite("cluster_topologies_undirected.csv",true);
			mfprintf(a,"n-mer;  id;  count;  countperc;  avgsigni / pm;  signisum / pm\n");
			for (z=0;z<m_oaClusterTopoUndir.GetSize();z++)
			{
				if (((CxObArray*)m_oaClusterTopoUndir[z])->GetSize() == 0)
					continue;
				
				for (z2=0;z2<((CxObArray*)m_oaClusterTopoUndir[z])->GetSize();z2++)
					((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iIndex = z2;

				// Stack Sort
				for (z2=0;z2<((CxObArray*)m_oaClusterTopoUndir[z])->GetSize();z2++)
				{
					j = 0;
					i = -1;
					for (z3=z2;z3<((CxObArray*)m_oaClusterTopoUndir[z])->GetSize();z3++)
					{
						if (((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z3))->m_iCountUndir > j)
						{
							j = ((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z3))->m_iCountUndir;
							i = z3;
						}
					}
					ec = (CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2);
					((CxObArray*)m_oaClusterTopoUndir[z])->SetAt(z2,((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(i));
					((CxObArray*)m_oaClusterTopoUndir[z])->SetAt(i,ec);
				}

				i = 0;
				for (z2=0;z2<((CxObArray*)m_oaClusterTopoUndir[z])->GetSize();z2++)
					i += ((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iCountUndir;

				if (m_bClusterTopoTrajectories)
				{
					for (z2=0;z2<((CxObArray*)m_oaClusterTopoUndir[z])->GetSize();z2++)
					{
#ifdef TARGET_WINDOWS
//						sprintf(buf,"ClusterTopo\\xtraj_undir_%02d_%04d.xyz",z+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iIndex+1);
//						sprintf(buf2,"ClusterTopo\\traj_undir_%02d_%04d.xyz",z+1,z2+1);
						buf.sprintf("ClusterTopo\\xtraj_undir_%02d_%04d.xyz",z+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iIndex+1);
						buf2.sprintf("ClusterTopo\\traj_undir_%02d_%04d.xyz",z+1,z2+1);
#else
//						sprintf(buf,"ClusterTopo/xtraj_undir_%02d_%04d.xyz",z+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iIndex+1);
//						sprintf(buf2,"ClusterTopo/traj_undir_%02d_%04d.xyz",z+1,z2+1);
						buf.sprintf("ClusterTopo/xtraj_undir_%02d_%04d.xyz",z+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iIndex+1);
						buf2.sprintf("ClusterTopo/traj_undir_%02d_%04d.xyz",z+1,z2+1);
#endif

						if (rename(buf,buf2) != 0)
							eprintf("Could not rename %s to %s.\n",(const char*)buf,(const char*)buf2);
					}
				}

				mprintf("        %d %d-mer topologies found. Rendering:",((CxObArray*)m_oaClusterTopoUndir[z])->GetSize(),z+1);
				for (z2=0;z2<((CxObArray*)m_oaClusterTopoUndir[z])->GetSize();z2++)
				{
					mprintf(".");

#ifdef TARGET_WINDOWS
//					sprintf(buf,"ClusterTopo\\cluster_undir_%02d_%04d",z+1,z2+1);
					buf.sprintf("ClusterTopo\\cluster_undir_%02d_%04d",z+1,z2+1);
#else
//					sprintf(buf,"ClusterTopo/cluster_undir_%02d_%04d",z+1,z2+1);
					buf.sprintf("ClusterTopo/cluster_undir_%02d_%04d",z+1,z2+1);
#endif

					DumpSchematicUndirectedDOT((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2),buf);
					mfprintf(a,"%d;  %d;  %d;  %f;  %f;  %f\n",z+1,z2+1,((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iCountUndir,100.0*((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iCountUndir/i,((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_fSignificanceUndir/((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iCountUndir,((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_fSignificanceUndir);
//					mfprintf(a,"    %8.4f%c: %d\n",100.0*((CExtendedCluster*)((CxObArray*)m_oaClusterTopoUndir[z])->GetAt(z2))->m_iCountUndir/i,'%',z2+1);
					fflush(a);
			//		sprintf(buf2,"del %s",buf);
			//		(void)!system(buf2);
				}
				mfprintf(a,"\n");
				mprintf("\n");
			}
			fclose(a);
		}
	}

	mprintf(WHITE,"  * Cluster Distance Distribution\n");
//	mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",m_pClusterDistanceDF->m_fBinEntries,m_pClusterDistanceDF->m_fSkipEntries,ZeroDivide(m_pClusterDistanceDF->m_fSkipEntries,m_pClusterDistanceDF->m_fBinEntries+m_pClusterDistanceDF->m_fSkipEntries)*100.0);
//	m_pClusterDistanceDF->CalcMeanSD();
	m_pClusterDistanceDF->MultiplyBin(100.0/m_iCounter);
//	mprintf("    Mean value: %.10G pm    Standard deviation: %.10G pm\n",m_pClusterDistanceDF->m_fMean,m_pClusterDistanceDF->m_fSD);
//	mprintf("    Min. value: %.10G pm    Max. value:         %.10G pm\n",m_pClusterDistanceDF->m_fMinInput,m_pClusterDistanceDF->m_fMaxInput);
//	sprintf(buf,"cluster_distance_df%s.csv",multibuf);
	buf.sprintf("cluster_distance_df%s.csv",multibuf);
	mprintf("    Saving Cluster Distance distribution as %s ...\n",(const char*)buf);
	m_pClusterDistanceDF->Write("",buf,"",false);
//	sprintf(buf,"cluster_distance_df%s.agr",multibuf);
	buf.sprintf("cluster_distance_df%s.agr",multibuf);
	mprintf("    Saving Cluster Distance distribution AGR file as \"%s\"...\n",(const char*)buf);
	m_pClusterDistanceDF->WriteAgr("",buf,"","Cluster Distance distribution",false);

	mprintf(WHITE,"\n  * Cluster Significance Distribution (Original)\n");
//	sprintf(buf,"cluster_significance_df%s.csv",multibuf);
	buf.sprintf("cluster_significance_df_original%s.csv",multibuf);
	mprintf("    Saving Cluster Significance distribution (original) as %s ...\n",(const char*)buf);
	m_pClusterDistributionDF->WriteMulti("",buf,"");
//	sprintf(buf,"cluster_significance_df%s.agr",multibuf);
	buf.sprintf("cluster_significance_df_original%s.agr",multibuf);
	mprintf("    Saving Cluster Significance distribution (original) AGR file as \"%s\"...\n",(const char*)buf);
	m_pClusterDistributionDF->WriteMultiAgr("",buf,"","Cluster Significance distribution (original)",false);


/*	sprintf(buf,"cluster_sizeX_df%s.csv",multibuf);
	mprintf("    Saving Cluster size X distribution as %s ...\n",buf);
	m_pClusterDistribution2DF->WriteMulti("",buf,"");
	sprintf(buf,"cluster_sizeX_df%s.agr",multibuf);
	mprintf("    Saving Cluster size X distribution AGR file as \"%s\"...\n",buf);
	m_pClusterDistribution2DF->WriteMultiAgr("",buf,"","Cluster size X distribution",false);*/

	mprintf(WHITE,"\n  * Cluster Size Distribution\n");
//	sprintf(buf,"cluster_size_df%s.csv",multibuf);
	buf.sprintf("cluster_size_df%s.csv",multibuf);
	mprintf("    Saving Cluster Size distribution as %s ...\n",(const char*)buf);
	m_pClusterSizeDF->WriteMulti("",buf,"");
//	sprintf(buf,"cluster_size_df%s.agr",multibuf);
	buf.sprintf("cluster_size_df%s.agr",multibuf);
	mprintf("    Saving Cluster Size distribution AGR file as \"%s\"...\n",(const char*)buf);
	m_pClusterSizeDF->WriteMultiAgr("",buf,"","Cluster Size distribution",false);

	mprintf(WHITE,"\n  * Cluster Count Distribution\n");
//	mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",m_pClusterCountDF->m_fBinEntries,m_pClusterCountDF->m_fSkipEntries,ZeroDivide(m_pClusterCountDF->m_fSkipEntries,m_pClusterCountDF->m_fBinEntries+m_pClusterCountDF->m_fSkipEntries)*100.0);
//	m_pClusterCountDF->CalcMeanSD();
	m_pClusterCountDF->MultiplyBin(1.0/m_iCounter);
//	mprintf("    Mean value: %.10G pm    Standard deviation: %.10G pm\n",m_pClusterCountDF->m_fMean,m_pClusterCountDF->m_fSD);
//	mprintf("    Min. value: %.10G pm    Max. value:         %.10G pm\n",m_pClusterCountDF->m_fMinInput,m_pClusterCountDF->m_fMaxInput);
//	sprintf(buf,"cluster_count_df%s.csv",multibuf);
	buf.sprintf("cluster_count_df%s.csv",multibuf);
	mprintf("    Saving Cluster Count distribution as %s ...\n",(const char*)buf);
	m_pClusterCountDF->Write("",buf,"",false);
//	sprintf(buf,"cluster_count_df%s.agr",multibuf);
	buf.sprintf("cluster_count_df%s.agr",multibuf);
	mprintf("    Saving Cluster Count distribution AGR file as \"%s\"...\n",(const char*)buf);
	m_pClusterCountDF->WriteAgr("",buf,"","Cluster Count distribution",false);

	mprintf(WHITE,"\n  * Polymer Distribution\n");
//	mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",m_pPolymerDF->m_fBinEntries,m_pPolymerDF->m_fSkipEntries,ZeroDivide(m_pPolymerDF->m_fSkipEntries,m_pPolymerDF->m_fBinEntries+m_pPolymerDF->m_fSkipEntries)*100.0);
//	m_pPolymerDF->CalcMeanSD();
	m_pPolymerDF->MultiplyBin(1.0/m_iCounter);
//	mprintf("    Mean value: %.10G pm    Standard deviation: %.10G pm\n",m_pPolymerDF->m_fMean,m_pPolymerDF->m_fSD);
//	mprintf("    Min. value: %.10G pm    Max. value:         %.10G pm\n",m_pPolymerDF->m_fMinInput,m_pPolymerDF->m_fMaxInput);
//	sprintf(buf,"cluster_polymer_df%s.csv",multibuf);
	buf.sprintf("cluster_polymer_df%s.csv",multibuf);
	mprintf("    Saving Polymer distribution as %s ...\n",(const char*)buf);
	m_pPolymerDF->WriteMulti("",buf,"");
//	sprintf(buf,"cluster_polymer_df_cumulative%s.csv",multibuf);
	buf.sprintf("cluster_polymer_df_cumulative%s.csv",multibuf);
	mprintf("    Saving cumulative Polymer distribution as %s ...\n",(const char*)buf);
	m_pPolymerDF->WriteMulti_Cumulative("",buf,"");
//	sprintf(buf,"cluster_polymer_df%s.agr",multibuf);
	buf.sprintf("cluster_polymer_df%s.agr",multibuf);
	mprintf("    Saving Polymer distribution AGR file as \"%s\"...\n",(const char*)buf);
	m_pPolymerDF->WriteMultiAgr("",buf,"","Polymer distribution",false);
//	sprintf(buf,"cluster_polymer_df_cumulative%s.agr",multibuf);
	buf.sprintf("cluster_polymer_df_cumulative%s.agr",multibuf);
	mprintf("    Saving cumulative Polymer distribution AGR file as \"%s\"...\n",(const char*)buf);
	m_pPolymerDF->WriteMultiAgr_Cumulative("",buf,"","Cumulative Polymer distribution",false);


	mprintf(WHITE,"\n  * Cluster Significance Distribution (Modified)\n");
	mprintf("    Using a relative threshold of %.5f to determine integration borders.\n",m_fCSDFThresh);
	tfm = 0;
	for (z=0;z<m_pClusterDistanceDF->m_iResolution;z++)
		if (m_pClusterDistanceDF->m_pBin[z] > tfm)
			tfm = m_pClusterDistanceDF->m_pBin[z];
	timin = -1;
	timax = -1;
	for (z=0;z<m_pClusterDistanceDF->m_iResolution;z++)
		if (m_pClusterDistanceDF->m_pBin[z] >= m_fCSDFThresh*tfm) {
			timin = z;
			break;
		}
	for (z=m_pClusterDistanceDF->m_iResolution-1;z>=0;z--)
		if (m_pClusterDistanceDF->m_pBin[z] >= m_fCSDFThresh*tfm) {
			timax = z;
			break;
		}
	if ((timin == -1) || (timax == -1)) {
		eprintf("CClusterAnalysis::WriteOutput(): Internal error: Could not determine integration borders.\n");
		abort();
	}
	mprintf(
		"    The integration borders are %.3f ... %.3f pm (bin %d ... %d / %d).\n",
		m_pClusterDistanceDF->m_fMinVal+(timin+0.5)*(m_pClusterDistanceDF->m_fMaxVal-m_pClusterDistanceDF->m_fMinVal)/m_pClusterDistanceDF->m_iResolution,
		m_pClusterDistanceDF->m_fMinVal+(timax+0.5)*(m_pClusterDistanceDF->m_fMaxVal-m_pClusterDistanceDF->m_fMinVal)/m_pClusterDistanceDF->m_iResolution,
		timin+1,
		timax+1,
		m_pClusterDistanceDF->m_iResolution
	);
	tf = 0;
	for (z=0;z<m_pClusterDistributionDF_Mod->m_iResolution;z++) {
		m_pClusterDistributionDF_Mod->m_pBin[z] = 0;
		for (z2=timin;z2<=timax;z2++)
			m_pClusterDistributionDF_Mod->m_pBin[z] += m_pPolymerDF->m_pMultiBin[z][z2];
		tf += m_pClusterDistributionDF_Mod->m_pBin[z];
	}
	tf = 100.0 / tf;
	for (z=0;z<m_pClusterDistributionDF_Mod->m_iResolution;z++)
		m_pClusterDistributionDF_Mod->m_pBin[z] *= tf;
	buf.sprintf("cluster_significance_df_modified%s.csv",multibuf);
	mprintf("    Saving Cluster Significance distribution (modified) as %s ...\n",(const char*)buf);
	m_pClusterDistributionDF_Mod->WriteMulti("",buf,"");
	buf.sprintf("cluster_significance_df_modified%s.agr",multibuf);
	mprintf("    Saving Cluster Significance distribution (modified) AGR file as \"%s\"...\n",(const char*)buf);
	m_pClusterDistributionDF_Mod->WriteMultiAgr("",buf,"","Cluster Significance distribution (modified)",false);

	
	if (m_bHetMeasure)
	{
		mprintf(WHITE,"\n  * Heterogeneity Time Development\n");
//		sprintf(buf,"cluster_hetero_td%s.csv",multibuf);
		buf.sprintf("cluster_hetero_td%s.csv",multibuf);
		mprintf("    Saving Heterogeneity Time Development as %s ...\n",(const char*)buf);
		a = OpenFileWrite(buf,true);
		mfprintf(a,"# Time / ps;  Min. Distance;  Max. Distance;  Center / pm;  Het. Measure;  Int. Het. Measure;  Int. Het. Measure NoSqr\n");
		for (z=0;z<m_faMinDist.GetSize();z++)
			mfprintf(a,"%G;  %G;  %G;  %G;  %G;  %G;  %G\n",z*g_iStride*g_fTimestepLength/1000.0,m_faMinDist[z],m_faMaxDist[z],m_faCenter[z],m_faHetNumber[z],m_faIntHetNumber[z],m_faIntHetNumber2[z]);
		fclose(a);
	}

	if (m_b2DPlots)
	{
		mprintf(WHITE,"\n  * 2D Cluster Distance Distribution\n");
		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_TRIPLES"))
		{
//			sprintf(buf,"cluster_distance_df%s_triples.csv",multibuf);
			buf.sprintf("cluster_distance_df%s_triples.csv",multibuf);
			mprintf("      Saving 2D Cluster Distance Distribution triples as %s ...\n",(const char*)buf);
			m_pClusterDistance2D->Write("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATRIX"))
		{
//			sprintf(buf,"cluster_distance_df%s_matrix.csv",multibuf);
			buf.sprintf("cluster_distance_df%s_matrix.csv",multibuf);
			mprintf("      Saving 2D Cluster Distance Distribution matrix as %s ...\n",(const char*)buf);
			m_pClusterDistance2D->WriteCSV("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA"))
		{
//			sprintf(buf,"cluster_distance_df%s.nb",multibuf);
			buf.sprintf("cluster_distance_df%s.nb",multibuf);
			mprintf("      Saving 2D Cluster Distance Distribution Mathematica notebook as %s ...\n",(const char*)buf);
			m_pClusterDistance2D->WriteMathematicaNb("",buf,"",false);
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_GNUPLOT"))
		{
//			sprintf(buf,"cluster_distance_df%s",multibuf);
			buf.sprintf("cluster_distance_df%s",multibuf);
			mprintf("      Saving 2D Cluster Distance Distribution Gnuplot Input as %s.gp ...\n",(const char*)buf);
			m_pClusterDistance2D->WriteGnuplotInput("",buf,"",false);
		}

		mprintf(WHITE,"\n  * 2D Cluster Size Distribution\n");
		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_TRIPLES"))
		{
//			sprintf(buf,"cluster_size_df%s_triples.csv",multibuf);
			buf.sprintf("cluster_size_df%s_triples.csv",multibuf);
			mprintf("      Saving 2D Cluster Size Distribution triples as %s ...\n",(const char*)buf);
			m_pClusterSize2D->Write("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATRIX"))
		{
//			sprintf(buf,"cluster_size_df%s_matrix.csv",multibuf);
			buf.sprintf("cluster_size_df%s_matrix.csv",multibuf);
			mprintf("      Saving 2D Cluster Size Distribution matrix as %s ...\n",(const char*)buf);
			m_pClusterSize2D->WriteCSV("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA"))
		{
//			sprintf(buf,"cluster_size_df%s.nb",multibuf);
			buf.sprintf("cluster_size_df%s.nb",multibuf);
			mprintf("      Saving 2D Cluster Size Distribution Mathematica notebook as %s ...\n",(const char*)buf);
			m_pClusterSize2D->WriteMathematicaNb("",buf,"",false);
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_GNUPLOT"))
		{
//			sprintf(buf,"cluster_size_df%s",multibuf);
			buf.sprintf("cluster_size_df%s",multibuf);
			mprintf("      Saving 2D Cluster Size Distribution Gnuplot Input as %s.gp ...\n",(const char*)buf);
			m_pClusterSize2D->WriteGnuplotInput("",buf,"",false);
		}

		mprintf(WHITE,"\n  * 2D Cluster Significance Distribution\n");
		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_TRIPLES"))
		{
//			sprintf(buf,"cluster_significance_df%s_triples.csv",multibuf);
			buf.sprintf("cluster_significance_df%s_triples.csv",multibuf);
			mprintf("      Saving 2D Cluster Significance Distribution triples as %s ...\n",(const char*)buf);
			m_pClusterDistribution2D->Write("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATRIX"))
		{
//			sprintf(buf,"cluster_significance_df%s_matrix.csv",multibuf);
			buf.sprintf("cluster_significance_df%s_matrix.csv",multibuf);
			mprintf("      Saving 2D Cluster Significance Distribution matrix as %s ...\n",(const char*)buf);
			m_pClusterDistribution2D->WriteCSV("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA"))
		{
//			sprintf(buf,"cluster_significance_df%s.nb",multibuf);
			buf.sprintf("cluster_significance_df%s.nb",multibuf);
			mprintf("      Saving 2D Cluster Significance Distribution Mathematica notebook as %s ...\n",(const char*)buf);
			m_pClusterDistribution2D->WriteMathematicaNb("",buf,"",false);
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_GNUPLOT"))
		{
//			sprintf(buf,"cluster_significance_df%s",multibuf);
			buf.sprintf("cluster_significance_df%s",multibuf);
			mprintf("      Saving 2D Cluster Significance Distribution Gnuplot Input as %s.gp ...\n",(const char*)buf);
			m_pClusterDistribution2D->WriteGnuplotInput("",buf,"",false);
		}

		mprintf(WHITE,"\n  * 2D Cluster Count Distribution\n");
		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_TRIPLES"))
		{
//			sprintf(buf,"cluster_count_df%s_triples.csv",multibuf);
			buf.sprintf("cluster_count_df%s_triples.csv",multibuf);
			mprintf("      Saving 2D Cluster Count Distribution triples as %s ...\n",(const char*)buf);
			m_pClusterCount2D->Write("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATRIX"))
		{
//			sprintf(buf,"cluster_count_df%s_matrix.csv",multibuf);
			buf.sprintf("cluster_count_df%s_matrix.csv",multibuf);
			mprintf("      Saving 2D Cluster Count Distribution matrix as %s ...\n",(const char*)buf);
			m_pClusterCount2D->WriteCSV("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA"))
		{
//			sprintf(buf,"cluster_count_df%s.nb",multibuf);
			buf.sprintf("cluster_count_df%s.nb",multibuf);
			mprintf("      Saving 2D Cluster Count Distribution Mathematica notebook as %s ...\n",(const char*)buf);
			m_pClusterCount2D->WriteMathematicaNb("",buf,"",false);
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_GNUPLOT"))
		{
//			sprintf(buf,"cluster_count_df%s",multibuf);
			buf.sprintf("cluster_count_df%s",multibuf);
			mprintf("      Saving 2D Cluster Count Distribution Gnuplot Input as %s.gp ...\n",(const char*)buf);
			m_pClusterCount2D->WriteGnuplotInput("",buf,"",false);
		}
	}

	if (m_bDiffPlot)
	{
		mprintf(WHITE,"\n  * 2D Cluster Distance Difference Plot\n");

		m_pClusterDiff2D->NormalizeXCount();
		m_pClusterDiff2D->NormalizeBinIntegral(1000000.0);

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_TRIPLES"))
		{
//			sprintf(buf,"cluster_distdiff_df%s_triples.csv",multibuf);
			buf.sprintf("cluster_distdiff_df%s_triples.csv",multibuf);
			mprintf("      Saving 2D Cluster Distance Difference triples as %s ...\n",(const char*)buf);
			m_pClusterDiff2D->Write("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATRIX"))
		{
//			sprintf(buf,"cluster_distdiff_df%s_matrix.csv",multibuf);
			buf.sprintf("cluster_distdiff_df%s_matrix.csv",multibuf);
			mprintf("      Saving 2D Cluster Distance Difference matrix as %s ...\n",(const char*)buf);
			m_pClusterDiff2D->WriteCSV("",buf,"");
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA"))
		{
//			sprintf(buf,"cluster_distdiff_df%s.nb",multibuf);
			buf.sprintf("cluster_distdiff_df%s.nb",multibuf);
			mprintf("      Saving 2D Cluster Distance Difference Mathematica notebook as %s ...\n",(const char*)buf);
			m_pClusterDiff2D->WriteMathematicaNb("",buf,"",false);
		}

		if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_GNUPLOT"))
		{
//			sprintf(buf,"cluster_distdiff_df%s",multibuf);
			buf.sprintf("cluster_distdiff_df%s",multibuf);
			mprintf("      Saving 2D Cluster Distance Difference Gnuplot Input as %s.gp ...\n",(const char*)buf);
			m_pClusterDiff2D->WriteGnuplotInput("",buf,"",false);
		}
	}

	if (m_bCutClusters)
	{
		mprintf(WHITE,"\n  * 2D Cluster Distance Distribution\n\n");
		for (z=0;z<m_iCutClusterMaxSize-1;z++)
		{
			mprintf("      - %d %d-mer clusters written to file cluster_%dmer.xyz\n",m_pCutClusterCounter[z],z+2,z+2);
			fclose(m_pCutClusterFiles[z]);
		}
	}

	if (m_bAnim)
	{
		mprintf("\n");

#ifdef TARGET_WINDOWS
		a = OpenFileWrite("cluster_agr\\gracebatch",true);
#else
		a = OpenFileWrite("cluster_agr/gracebatch",true);
#endif

		mfprintf(a,"PRINT TO \"output.png\"\n");
		mfprintf(a,"HARDCOPY DEVICE \"PNG\"\n");
		mfprintf(a,"PAGE SIZE %d, %d\n",m_iResX,m_iResY);
		mfprintf(a,"DEVICE \"PNG\" FONT ANTIALIASING on\n");
		mfprintf(a,"DEVICE \"PNG\" OP \"compression:9\"\n");
		mfprintf(a,"PRINT\n");
		fclose(a);
		mprintf("    Saved batch script as \"render_cluster_anim\".\n");
		fclose(m_fAnim);
	}
	mprintf("\n");
}


CExtendedCluster::CExtendedCluster()
{
}


CExtendedCluster::~CExtendedCluster()
{
	int z;

	for (z=0;z<m_oaBonds.GetSize();z++)
		delete (CxIntArray*)m_oaBonds[z];

	for (z=0;z<m_oaRelBonds.GetSize();z++)
		delete (CxIntArray*)m_oaRelBonds[z];

	for (z=0;z<m_oaRelBondsUndir.GetSize();z++)
		delete (CxIntArray*)m_oaRelBondsUndir[z];
}


void CExtendedCluster::CreateFrom(CClusterNode *n, CTimeStep *ts, CClusterAnalysis *ca, bool wrap)
{
	int z, z2, z3, ti, ti2;
	CxIntArray *ia;
	double dist, hbthres, dispthres;
	CSingleMolecule *sm;

	m_iaAtoms.SetSize(n->m_iaAtoms.GetSize());
	for (z=0;z<n->m_iaAtoms.GetSize();z++)
		m_iaAtoms[z] = ca->m_iaAtomList[n->m_iaAtoms[z]];
//	m_iaAtoms.CopyFrom(&n->m_iaAtoms);
	m_oaBonds.SetSize(m_iaAtoms.GetSize());
	m_oaRelBonds.SetSize(m_iaAtoms.GetSize());
	m_iClusterIndex = n->m_iIndex;

	if (ca->m_bUseBondCutoffs)
	{
		hbthres = ca->m_fClusterTopoHBThres;
		if (!ca->m_bCOMTrick)
			if (hbthres < n->m_fDistance)
				hbthres = n->m_fDistance;
	} else hbthres = n->m_fDistance;

	if (ca->m_bUseBondCutoffs)
	{
		dispthres = ca->m_fClusterTopoDispThres;
		if (!ca->m_bCOMTrick)
			if (dispthres < n->m_fDistance)
				dispthres = n->m_fDistance;
	} else dispthres = n->m_fDistance;

//	mprintf("cluster = %f pm, hbthres = %f pm, dispthres = %f pm.\n",n->m_fDistance,hbthres,dispthres);

	m_iMonomers = n->m_iMonomers;
	m_iInterBonds = 0;

	for (z=0;z<m_iaAtoms.GetSize();z++)
	{
		ia = new CxIntArray();
		m_oaBonds[z] = ia;
		ia = new CxIntArray();
		m_oaRelBonds[z] = ia;
	}

	for (z=0;z<m_iaAtoms.GetSize();z++)
	{
		for (z2=z+1;z2<m_iaAtoms.GetSize();z2++)
		{
			{
				if (wrap)
					dist = FoldedLength(ts->m_vaCoords[m_iaAtoms[z]] - ts->m_vaCoords[m_iaAtoms[z2]]);
						else dist = (ts->m_vaCoords[m_iaAtoms[z]] - ts->m_vaCoords[m_iaAtoms[z2]]).GetLength();
			}

			if (g_laAtomSMIndex[m_iaAtoms[z]] == g_laAtomSMIndex[m_iaAtoms[z2]])
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[g_laAtomSMIndex[m_iaAtoms[z]]];
				for (z3=0;z3<sm->m_laBonds.GetSize()/2;z3++)
				{
					if (((sm->m_laBonds[z3*2] == m_iaAtoms[z]) && (sm->m_laBonds[z3*2+1] == m_iaAtoms[z2])) ||
						((sm->m_laBonds[z3*2] == m_iaAtoms[z2]) && (sm->m_laBonds[z3*2+1] == m_iaAtoms[z])))
					{
						((CxIntArray*)m_oaBonds[z])->Add(m_iaAtoms[z2]);
						((CxIntArray*)m_oaBonds[z2])->Add(m_iaAtoms[z]);

						((CxIntArray*)m_oaRelBonds[z])->Add(z2);
						((CxIntArray*)m_oaRelBonds[z2])->Add(z);
					}
				}
			} else
			{
				if (ca->m_bSelectBonds)
				{
	//				mprintf("%d -- %d, %d -- %d, (%d,%d), (%d,%d)\n",z,z2,n->m_iaAtoms[z],n->m_iaAtoms[z2],ca->m_iaAtomBondFrom[n->m_iaAtoms[z]],ca->m_iaAtomBondTo[n->m_iaAtoms[z2]],ca->m_iaAtomBondFrom[n->m_iaAtoms[z2]],ca->m_iaAtomBondTo[n->m_iaAtoms[z]]);
					if (((ca->m_iaAtomBondFrom[n->m_iaAtoms[z]] != 0) && (ca->m_iaAtomBondTo[n->m_iaAtoms[z2]] != 0)) || ((ca->m_iaAtomBondFrom[n->m_iaAtoms[z2]] != 0) && (ca->m_iaAtomBondTo[n->m_iaAtoms[z]] != 0)))
					{
	//					mprintf("  # accepted, dist = %f vs %f\n",dist,hbthres);
						if (dist <= hbthres+0.1)
						{
							m_iInterBonds++;

							((CxIntArray*)m_oaBonds[z])->Add(m_iaAtoms[z2]);
							((CxIntArray*)m_oaBonds[z2])->Add(m_iaAtoms[z]);

							((CxIntArray*)m_oaRelBonds[z])->Add(z2);
							((CxIntArray*)m_oaRelBonds[z2])->Add(z);
						}
					}
				} else
				{
					if (IsHydrogenBond(g_waAtomRealElement[m_iaAtoms[z]],g_waAtomRealElement[m_iaAtoms[z2]]))
					{
						if (dist <= hbthres+0.1)
						{
							m_iInterBonds++;

						//	mprintf("Adding bond A: %f\n",dist);

							((CxIntArray*)m_oaBonds[z])->Add(m_iaAtoms[z2]);
							((CxIntArray*)m_oaBonds[z2])->Add(m_iaAtoms[z]);

							((CxIntArray*)m_oaRelBonds[z])->Add(z2);
							((CxIntArray*)m_oaRelBonds[z2])->Add(z);
						}
					} else if (ca->m_bClusterTopoAllBonds)
					{
						if (dist <= dispthres+0.1)
						{
							m_iInterBonds++;

						//	mprintf("Adding bond B: %f\n",dist);

							((CxIntArray*)m_oaBonds[z])->Add(m_iaAtoms[z2]);
							((CxIntArray*)m_oaBonds[z2])->Add(m_iaAtoms[z]);

							((CxIntArray*)m_oaRelBonds[z])->Add(z2);
							((CxIntArray*)m_oaRelBonds[z2])->Add(z);
						}
					}
				}
			}
		}
	}

//	mprintf("%d Monomers, %d Bonds.\n",m_iMonomers,m_iInterBonds);

	if (ca->m_bClusterTopoDrawUndirected)
	{
		m_oaRelBondsUndir.SetSize(m_iMonomers);

		for (z=0;z<m_iMonomers;z++)
		{
			ia = new CxIntArray();
			m_oaRelBondsUndir[z] = ia;
		}

		for (z=0;z<m_iaAtoms.GetSize();z++)
		{
			for (z2=0;z2<m_iaMonomerSM.GetSize();z2++)
				if (g_laAtomSMIndex[m_iaAtoms[z]] == m_iaMonomerSM[z2])
					goto _found;
			m_iaMonomerSM.Add(g_laAtomSMIndex[m_iaAtoms[z]]);
_found:;
		}

//		mprintf("Building up %d-mer (%d), d=%fpm, %d InterBonds...\n",m_iMonomers,m_iaMonomerSM.GetSize(),n->m_fDistance,m_iInterBonds);

		for (z=0;z<m_iaAtoms.GetSize();z++)
		{
			ti = g_laAtomSMIndex[m_iaAtoms[z]];
			for (z2=0;z2<m_iaMonomerSM.GetSize();z2++)
			{
				if (ti == m_iaMonomerSM[z2])
				{
					ti = z2;
					goto _found2;
				}
			}
			eprintf("CExtendedCluster::CreateFrom(): Internal error 1.\n");
_found2:
			for (z2=0;z2<((CxIntArray*)m_oaRelBonds[z])->GetSize();z2++)
			{
		/*		if (ca->m_bSelectBonds)
				{
					if (((ca->m_iaAtomBondFrom[n->m_iaAtoms[z]] == 0) || (ca->m_iaAtomBondTo[n->m_iaAtoms[((CxIntArray*)m_oaRelBonds[z])->GetAt(z2)]] == 0)) && ((ca->m_iaAtomBondFrom[n->m_iaAtoms[((CxIntArray*)m_oaRelBonds[z])->GetAt(z2)]] == 0) || (ca->m_iaAtomBondTo[n->m_iaAtoms[z]] == 0)))
						goto _have;
				}*/
				ti2 = g_laAtomSMIndex[((CxIntArray*)m_oaBonds[z])->GetAt(z2)];
				if (ti2 != g_laAtomSMIndex[m_iaAtoms[z]])
				{
					for (z3=0;z3<m_iaMonomerSM.GetSize();z3++)
					{
						if (ti2 == m_iaMonomerSM[z3])
						{
							ti2 = z3;
							goto _found3;
						}
					}
					eprintf("CExtendedCluster::CreateFrom(): Internal error 2.\n");
_found3:
					for (z3=0;z3<((CxIntArray*)m_oaRelBondsUndir[ti])->GetSize();z3++)
					{
						if (((CxIntArray*)m_oaRelBondsUndir[ti])->GetAt(z3) == ti2)
							goto _have;
					}
					((CxIntArray*)m_oaRelBondsUndir[ti])->Add(ti2);
					((CxIntArray*)m_oaRelBondsUndir[ti2])->Add(ti);
//					mprintf("  Adding Bond %d - %d.\n",ti,ti2);
_have:;
				}
			}
		}
	}
}


void CClusterAnalysis::DumpAtomDOT(CExtendedCluster *ec, const char *s)
{
	int z, z2;
	double tf, tf2;
	FILE *a;
//	char buf[256];
	CxString buf;

//	sprintf(buf,"%s.dot",s);
	buf.sprintf("%s.dot",s);
	a = OpenFileWrite(buf,true);

	mfprintf(a,"graph molecule {\n");
	mfprintf(a,"  graph [pack=true,splines=true,overlap=false];\n");
	mfprintf(a,"  node [shape=none,fontsize=16,fontname=\"Arial\",margin=0,fixedsize=true,height=0.28];\n");
	mfprintf(a,"  edge [style=bold,len=0.70];\n");

	for (z=0;z<ec->m_iaAtoms.GetSize();z++)
		mfprintf(a,"  _%d [ label=\"%s\", width=%.2f ];\n",ec->m_iaAtoms[z]+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[ec->m_iaAtoms[z]]])->m_sName,strlen((const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[ec->m_iaAtoms[z]]])->m_sName)*0.17);

	for (z=0;z<ec->m_iaAtoms.GetSize();z++)
	{
		for (z2=0;z2<((CxIntArray*)ec->m_oaBonds[z])->GetSize();z2++)
		{
			if (((CxIntArray*)ec->m_oaRelBonds[z])->GetAt(z2) < z)
				continue;

			if (g_laAtomSMIndex[ec->m_iaAtoms[z]] == g_laAtomSMIndex[((CxIntArray*)ec->m_oaBonds[z])->GetAt(z2)])
			{
				tf = 2.5;
				tf2 = 0.4;
			} else
			{
				tf = 1.0;
				tf2 = 0.8;
			}

			mfprintf(a,"  _%d -- _%d [ len=%.3f, weight=%d, penwidth=%f ];\n",ec->m_iaAtoms[z]+1,((CxIntArray*)ec->m_oaBonds[z])->GetAt(z2)+1,tf2,100,tf);
		}
	}
	mfprintf(a,"}\n");

	fclose(a);

//	sprintf(buf,"dot %s -Tpng -o%s.png -Kneato -Gepsilon=0.000001 -Gns0limit=5000 -Gmclimit=5000 -v -Gstart=rand",s,s);

//	mprintf("Command: \"%s\".\n",buf);

	RenderFormula(s);

//	(void)!system(buf);
}


void CClusterAnalysis::DumpSchematicDirectedDOT(CExtendedCluster *ec, const char *s)
{
	int z, z2, cr, cg, cb;
	FILE *a;
//	char buf[256];
	CxString buf;

//	sprintf(buf,"%s.dot",s);
	buf.sprintf("%s.dot",s);
	a = OpenFileWrite(buf,true);

	mfprintf(a,"digraph molecule {\n");
    
	mfprintf(a,"  graph [ pack=true, overlap=false ];\n");
	mfprintf(a,"  node [ label = \"\", penwidth=3.0, shape = circle, height=0.3, style=\"filled\" ];\n");
	mfprintf(a,"  edge [ penwidth=3.0, arrowsize = 1.5];\n\n");

	for (z=0;z<m_iaClusterTopoSMTemp.GetSize();z++)
		m_iaClusterTopoSMTemp[z] = 0;

	for (z=0;z<ec->m_iaAtoms.GetSize();z++)
	{
		if (m_iaClusterTopoSMTemp[g_laAtomSMIndex[ec->m_iaAtoms[z]]] != 0)
			continue;
		m_iaClusterTopoSMTemp[g_laAtomSMIndex[ec->m_iaAtoms[z]]] = 1;
		for (z2=0;z2<m_oaAtomGroups.GetSize();z2++)
		{
			if (((CAtomGroup*)m_oaAtomGroups[z2])->m_pMolecule == (CMolecule*)g_oaMolecules[g_waAtomMolIndex[ec->m_iaAtoms[z]]])
			{
				cr = m_iaClusterTopoMolColor[z2*3];
				cg = m_iaClusterTopoMolColor[z2*3+1];
				cb = m_iaClusterTopoMolColor[z2*3+2];
				goto _found;
			}
		}
		cr = 211;
		cg = 211;
		cb = 211;
_found:
		mfprintf(a,"  _%d [ fillcolor=\"#%02X%02X%02X\" ];\n",g_laAtomSMIndex[ec->m_iaAtoms[z]],cr,cg,cb);
	}

	mfprintf(a,"\n");

	for (z=0;z<ec->m_iaAtoms.GetSize();z++)
	{
		for (z2=0;z2<((CxIntArray*)ec->m_oaBonds[z])->GetSize();z2++)
		{
			if (((CxIntArray*)ec->m_oaRelBonds[z])->GetAt(z2) < z)
				continue;

			if (g_laAtomSMIndex[ec->m_iaAtoms[z]] == g_laAtomSMIndex[((CxIntArray*)ec->m_oaBonds[z])->GetAt(z2)])
				continue;

			if (ec->IsHydrogenBond(g_waAtomRealElement[ec->m_iaAtoms[z]],g_waAtomRealElement[((CxIntArray*)ec->m_oaBonds[z])->GetAt(z2)]))
			{
				if (((CAtom*)g_oaAtoms[g_waAtomRealElement[ec->m_iaAtoms[z]]])->m_pElement->m_fMass > ((CAtom*)g_oaAtoms[g_waAtomRealElement[((CxIntArray*)ec->m_oaBonds[z])->GetAt(z2)]])->m_pElement->m_fMass)
					mfprintf(a,"  _%d -> _%d;\n",g_laAtomSMIndex[((CxIntArray*)ec->m_oaBonds[z])->GetAt(z2)],g_laAtomSMIndex[ec->m_iaAtoms[z]]);
						else mfprintf(a,"  _%d -> _%d;\n",g_laAtomSMIndex[ec->m_iaAtoms[z]],g_laAtomSMIndex[((CxIntArray*)ec->m_oaBonds[z])->GetAt(z2)]);
			} else mfprintf(a,"  _%d -> _%d [penwidth=1.5,arrowsize = 0];\n",g_laAtomSMIndex[((CxIntArray*)ec->m_oaBonds[z])->GetAt(z2)],g_laAtomSMIndex[ec->m_iaAtoms[z]]);
		}
	}
	mfprintf(a,"}\n");

	fclose(a);

	RenderFormula(s);

//	sprintf(buf,"dot %s -Tpng -o%s.png -Kneato -Gepsilon=0.000001 -Gns0limit=5000 -Gmclimit=5000 -v -Gstart=rand",s,s);

//	mprintf("Command: \"%s\".\n",buf);

//	(void)!system(buf);
}


void CClusterAnalysis::DumpSchematicUndirectedDOT(CExtendedCluster *ec, const char *s)
{
	int z, z2, cr, cg, cb;
	FILE *a;
//	char buf[256];
	CxString buf;

//	sprintf(buf,"%s.dot",s);
	buf.sprintf("%s.dot",s);
	a = OpenFileWrite(buf,true);

	mfprintf(a,"graph molecule {\n");
    
	mfprintf(a,"  graph [ pack=true, overlap=false ];\n");
	mfprintf(a,"  node [ label = \"\", penwidth=3.0, shape = circle, height=0.3, style=\"filled\" ];\n");
	mfprintf(a,"  edge [ penwidth=3.0 ];\n\n");

	for (z=0;z<ec->m_iaMonomerSM.GetSize();z++)
	{
		for (z2=0;z2<m_oaAtomGroups.GetSize();z2++)
		{
			if (((CAtomGroup*)m_oaAtomGroups[z2])->m_pMolecule == (CMolecule*)g_oaMolecules[((CSingleMolecule*)g_oaSingleMolecules[ec->m_iaMonomerSM[z]])->m_iMolType])
			{
				cr = m_iaClusterTopoMolColor[z2*3];
				cg = m_iaClusterTopoMolColor[z2*3+1];
				cb = m_iaClusterTopoMolColor[z2*3+2];
				goto _found;
			}
		}
		cr = 211;
		cg = 211;
		cb = 211;
_found:
		mfprintf(a,"  _%d [ fillcolor=\"#%02X%02X%02X\" ];\n",z,cr,cg,cb);
	}

	mfprintf(a,"\n");

	for (z=0;z<ec->m_iaMonomerSM.GetSize();z++)
	{
		for (z2=0;z2<((CxIntArray*)ec->m_oaRelBondsUndir[z])->GetSize();z2++)
		{
			if (((CxIntArray*)ec->m_oaRelBondsUndir[z])->GetAt(z2) < z)
				continue;

			mfprintf(a,"  _%d -- _%d;\n",z,((CxIntArray*)ec->m_oaRelBondsUndir[z])->GetAt(z2));
		}
	}
	mfprintf(a,"}\n");

	fclose(a);

	RenderFormula(s);

//	sprintf(buf,"dot %s -Tpng -o%s.png -Kneato -Gepsilon=0.000001 -Gns0limit=5000 -Gmclimit=5000 -v -Gstart=rand",s,s);

//	mprintf("Command: \"%s\".\n",buf);

//	(void)!system(buf);
}


void CExtendedCluster::BuildAtomCodes()
{
	int z, z2, c1, c2, i;
	double ac;

	m_faAtomCodes.SetSize(m_iaAtoms.GetSize());
	m_faTempAtomCodes.SetSize(m_iaAtoms.GetSize());
	m_faCountAC.SetSize(m_iaAtoms.GetSize());

	// Die Anfangswerte der AtomCodes: [Ordnungszahl] * 10.0 + [Zahl der Nicht-Wasserstoff-Bindungen]
	for (z=0;z<m_iaAtoms.GetSize();z++)
	{
		m_faAtomCodes[z] = 10.0 * ((CAtom*)g_oaAtoms[g_waAtomRealElement[m_iaAtoms[z]]])->m_pElement->m_fMass;

		i = 0;
		for (z2=0;z2<((CxIntArray*)m_oaBonds[z])->GetSize();z2++)
		{
			// Alle Wasserstoff-Atome ueberspringen
			if (((CAtom*)g_oaAtoms[g_waAtomRealElement[((CxIntArray*)m_oaBonds[z])->GetAt(z2)]])->m_pElement->m_fMass < 4.5)
				continue;
			m_faAtomCodes[z]++;
			i++;
		}
	}

	i = 0;
	do {
		c1 = CountDifferentAtomCodes();

//		if (g_bVerbose)
//			mprintf(WHITE,"\n  Cycle %d: %d different atom codes exist.\n\n",i+1,c1);

		for (z=0;z<m_iaAtoms.GetSize();z++)
		{
			m_faTempAtomCodes[z] = m_faAtomCodes[z] * 5.0;

			for (z2=0;z2<((CxIntArray*)m_oaRelBonds[z])->GetSize();z2++)
				m_faTempAtomCodes[z] += m_faAtomCodes[((CxIntArray*)m_oaRelBonds[z])->GetAt(z2)];
		}
		for (z=0;z<m_iaAtoms.GetSize();z++)
			m_faAtomCodes[z] = m_faTempAtomCodes[z];

		c2 = CountDifferentAtomCodes();
		i++;

	} while (c1 != c2);

	// Sortieren mittels StackSort
	for (z=0;z<m_iaAtoms.GetSize();z++)
	{
		ac = -1;
		i = -1;
		for (z2=z;z2<m_iaAtoms.GetSize();z2++)
		{
			if (m_faAtomCodes[z2] > ac)
			{
				ac = m_faAtomCodes[z2];
				i = z2;
			}
		}
		if (i != -1)
		{
			m_faAtomCodes[i] = m_faAtomCodes[z];
			m_faAtomCodes[z] = ac;
		} else
		{
			eprintf("CExtendedCluster::BuildAtomCodes(): Internal error.\n");
			return;
		}
	}

/*	mprintf("AtomCodes:\n");
	for (z=0;z<m_iaAtoms.GetSize();z++)
		mprintf("  - %.0f\n",m_faAtomCodes[z]);*/
}


int CExtendedCluster::CountDifferentAtomCodes()
{
	int z, z2, i;
	
	i = 0;
	for (z=0;z<m_faAtomCodes.GetSize();z++)
	{
		for (z2=0;z2<i;z2++)
			if (m_faCountAC[z2] == m_faAtomCodes[z])
				goto _next;
		m_faCountAC[i] = m_faAtomCodes[z];
		i++;
_next:;
	}
	return i;
}


void CExtendedCluster::BuildAtomCodesUndir()
{
	int z, z2, c1, c2, i;
	double ac;

	m_faAtomCodesUndir.SetSize(m_iaMonomerSM.GetSize());
	m_faTempAtomCodesUndir.SetSize(m_iaMonomerSM.GetSize());
	if (m_faCountAC.GetSize() < m_iaMonomerSM.GetSize())
		m_faCountAC.SetSize(m_iaMonomerSM.GetSize());

	// Die Anfangswerte der AtomCodes: [Ordnungszahl] * 10.0 + [Zahl der Nicht-Wasserstoff-Bindungen]
	for (z=0;z<m_iaMonomerSM.GetSize();z++)
	{
		m_faAtomCodesUndir[z] = 10.0 * ((CMolecule*)g_oaMolecules[((CSingleMolecule*)g_oaSingleMolecules[m_iaMonomerSM[z]])->m_iMolType])->m_fMass;

		i = 0;
		for (z2=0;z2<((CxIntArray*)m_oaRelBondsUndir[z])->GetSize();z2++)
		{
			m_faAtomCodesUndir[z]++;
			i++;
		}
	}

	i = 0;
	do {
		c1 = CountDifferentAtomCodesUndir();

//		if (g_bVerbose)
//			mprintf(WHITE,"\n  Cycle %d: %d different atom codes exist.\n\n",i+1,c1);

		for (z=0;z<m_iaMonomerSM.GetSize();z++)
		{
			m_faTempAtomCodesUndir[z] = m_faAtomCodesUndir[z] * 5.0;

			for (z2=0;z2<((CxIntArray*)m_oaRelBondsUndir[z])->GetSize();z2++)
				m_faTempAtomCodesUndir[z] += m_faAtomCodesUndir[((CxIntArray*)m_oaRelBondsUndir[z])->GetAt(z2)];
		}
		for (z=0;z<m_iaMonomerSM.GetSize();z++)
			m_faAtomCodesUndir[z] = m_faTempAtomCodesUndir[z];

		c2 = CountDifferentAtomCodesUndir();
		i++;

	} while (c1 != c2);

	// Sortieren mittels StackSort
	for (z=0;z<m_iaMonomerSM.GetSize();z++)
	{
		ac = -1;
		i = -1;
		for (z2=z;z2<m_iaMonomerSM.GetSize();z2++)
		{
			if (m_faAtomCodesUndir[z2] > ac)
			{
				ac = m_faAtomCodesUndir[z2];
				i = z2;
			}
		}
		if (i != -1)
		{
			m_faAtomCodesUndir[i] = m_faAtomCodesUndir[z];
			m_faAtomCodesUndir[z] = ac;
		} else
		{
			eprintf("CExtendedCluster::BuildAtomCodesUndir(): Internal error.\n");
			return;
		}
	}

/*	mprintf("AtomCodes:\n");
	for (z=0;z<m_iaAtoms.GetSize();z++)
		mprintf("  - %.0f\n",m_faAtomCodes[z]);*/
}


int CExtendedCluster::CountDifferentAtomCodesUndir()
{
	int z, z2, i;
	
	i = 0;
	for (z=0;z<m_faAtomCodesUndir.GetSize();z++)
	{
		for (z2=0;z2<i;z2++)
			if (m_faCountAC[z2] == m_faAtomCodesUndir[z])
				goto _next;
		m_faCountAC[i] = m_faAtomCodesUndir[z];
		i++;
_next:;
	}
	return i;
}


bool CClusterAnalysis::AddToClusterTopo(CExtendedCluster *ec, int *i)
{
	int z, z2;
	CxObArray *oa;
	CExtendedCluster *ec2;

	oa = ((CxObArray*)m_oaClusterTopo[ec->m_iMonomers-1]);

	for (z=0;z<oa->GetSize();z++)
	{
		ec2 = (CExtendedCluster*)oa->GetAt(z);
		for (z2=0;z2<ec->m_faAtomCodes.GetSize();z2++)
			if (ec->m_faAtomCodes[z2] != ec2->m_faAtomCodes[z2])
				goto _diff;
		ec2->m_iCount++;
		ec2->m_fSignificance += ec->m_fSignificance;
		ec2->m_fSignificanceUndir += ec->m_fSignificanceUndir;
		*i = z;
		goto _found;
_diff:;
	}
//	mprintf("Add\n");
	if (IsConnected(ec))
	{
		*i = oa->GetSize();
		oa->Add(ec);
		return true;
	}

_found:
	return false;
}


bool CClusterAnalysis::AddToClusterTopoUndir(CExtendedCluster *ec, int *i)
{
	int z, z2;
	CxObArray *oa;
	CExtendedCluster *ec2;

	oa = ((CxObArray*)m_oaClusterTopoUndir[ec->m_iMonomers-1]);

	for (z=0;z<oa->GetSize();z++)
	{
		ec2 = (CExtendedCluster*)oa->GetAt(z);
		for (z2=0;z2<ec->m_faAtomCodesUndir.GetSize();z2++)
			if (ec->m_faAtomCodesUndir[z2] != ec2->m_faAtomCodesUndir[z2])
				goto _diff;
		ec2->m_iCountUndir++;
		*i = z;
		goto _found;
_diff:;
	}
//	mprintf("Add\n");
	if (IsConnected(ec))
	{
		*i = oa->GetSize();
		oa->Add(ec);
//		mprintf("Take.\n");
		return true;
	}// else mprintf("Skip B.\n");

_found:
//	mprintf("Skip C.\n");
	return false;
}


void CClusterAnalysis::RenderFormula(const char *s)
{
	int z, zm, tries;
//	char buf[256], buf2[256];
	CxString buf, buf2;
	char *p, *q;
	double tf, mi, ma, av;
	FILE *a;
	
	tries = m_iClusterTopoTries;

//	mprintf("%d Tries:\n",tries);
	zm = -1;
	av = 0;
	ma = 0;
	mi = 1e30;
//	mprintf("      Command: dot %s.dot -Tpng -o%s.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=rand\n",s,s);
//	mprintf("      Optimizing (%d tries): ",tries);
//	mprintf(WHITE,"[");
	for (z=0;z<tries;z++)
	{
//		mprintf(WHITE,"#");

#ifdef TARGET_WINDOWS
//		sprintf(buf,"dot %s.dot -Tpng -o%s.%d.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=%d > ClusterTopo\\dot%d.log 2>&1",s,s,z,rand(),z);
		buf.sprintf("dot %s.dot -Tpng -o%s.%d.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=%d > ClusterTopo\\dot%d.log 2>&1",s,s,z,rand(),z);
#else
//		sprintf(buf,"dot %s.dot -Tpng -o%s.%d.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=%d > ClusterTopo/dot%d.log 2>&1",s,s,z,rand(),z);
		buf.sprintf("dot %s.dot -Tpng -o%s.%d.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=%d > ClusterTopo/dot%d.log 2>&1",s,s,z,rand(),z);
#endif

//		mprintf("    (%s)\n",buf);

		(void)!system(buf);

#ifdef TARGET_WINDOWS
//		sprintf(buf,"ClusterTopo\\dot%d.log",z);
		buf.sprintf("ClusterTopo\\dot%d.log",z);
#else
//		sprintf(buf,"ClusterTopo/dot%d.log",z);
		buf.sprintf("ClusterTopo/dot%d.log",z);
#endif

		a = fopen(buf,"rt");
		if (a == NULL)
		{
			eprintf("\nRenderFormula(): Error opening %s for reading.\n",(const char*)buf);
			eprintf("It seems that GraphViz is not working (try command \"dot\").\n\n");
			return;
		}
		while (!feof(a))
		{
//			(void)!fgets(buf,256,a);
			(void)!buf.fgets(256,a);
			if (strstr(buf,"final") != NULL)
			{
				p = buf.GetWritePointer();
	//			mprintf("    buf=\"%s\".\n",buf);
				while (((*p < '0') || (*p > '9')) && (*p != 0))
					p++;
	//			mprintf("    p=\"%s\".\n",p);
				if (*p == 0)
					continue;
				q = p;
//				mprintf("    *p = %d.\n",*p);
				while (((*p >= '0') && (*p <= '9')) || (*p == '.'))
				{
	//				mprintf("    *p = %d --> p++;\n",*p);
					p++;
				}
	//			mprintf("    p2=\"%s\".\n",p);
				*p = 0;
		//		mprintf("    q=\"%s\".\n",q);
				tf = atof(q);
//				mprintf("    tf = %g.\n",tf);
				if (tf > ma)
					ma = tf;
				if (tf < mi)
				{
					mi = tf;
					zm = z;
				}
				av += tf;
				goto _done;
			}
		}
		eprintf("\nError with GraphViz output. See dot%d.log.\n",z);
		eprintf("It seems like GraphViz is not working (try command \"dot\").\n\n");
		return;
_done:
		fclose(a);

#ifdef TARGET_WINDOWS
//		sprintf(buf,"ClusterTopo\\dot%d.log",z);
		buf.sprintf("ClusterTopo\\dot%d.log",z);
#else
//		sprintf(buf,"ClusterTopo/dot%d.log",z);
		buf.sprintf("ClusterTopo/dot%d.log",z);
#endif

		remove(buf);
	}
//	mprintf(WHITE,"]\n");
	av /= tries;
	for (z=0;z<tries;z++)
	{
		if (z == zm)
			continue;
//		sprintf(buf,"%s.%d.png",s,z);
		buf.sprintf("%s.%d.png",s,z);
		if (remove(buf) != 0)
			eprintf("      Error removing \"%s\".\n",(const char*)buf);
	}
//	sprintf(buf,"%s.%d.png",s,zm);
//	sprintf(buf2,"%s.png",s);
	buf.sprintf("%s.%d.png",s,zm);
	buf2.sprintf("%s.png",s);
	remove(buf2);
	if (rename(buf,buf2) != 0)
		eprintf("      Error renaming \"%s\" to \"%s\".\n",(const char*)buf,(const char*)buf2);
//	sprintf(buf,"%s.dot",s);
	buf.sprintf("%s.dot",s);
	remove(buf);
}


bool CExtendedCluster::IsHydrogenBond(int i1, int i2)
{
	CAtom *a1, *a2;

	a1 = (CAtom*)g_oaAtoms[i1];
	a2 = (CAtom*)g_oaAtoms[i2];

	if (mystricmp(a1->m_sName,"H") == 0)
	{
		if (mystricmp(a2->m_sName,"O") == 0)
			return true;
		if (mystricmp(a2->m_sName,"N") == 0)
			return true;
	}

	if (mystricmp(a2->m_sName,"H") == 0)
	{
		if (mystricmp(a1->m_sName,"O") == 0)
			return true;
		if (mystricmp(a1->m_sName,"N") == 0)
			return true;
	}

	return false;
}


bool CClusterAnalysis::IsConnected(CExtendedCluster *ec)
{
	int z;

	if (m_iaConnectedBuffer.GetSize() < ec->m_iMonomers)
		m_iaConnectedBuffer.SetSize(ec->m_iMonomers);

	for (z=0;z<ec->m_iMonomers;z++)
		m_iaConnectedBuffer[z] = 0;

	REC_IsConnected(ec,0);

	for (z=0;z<ec->m_iMonomers;z++)
		if (m_iaConnectedBuffer[z] == 0)
			return false;

	return true;
}


void CClusterAnalysis::REC_IsConnected(CExtendedCluster *ec, int i)
{
	int z;

	m_iaConnectedBuffer[i] = 1;

	for (z=0;z<((CxIntArray*)ec->m_oaRelBondsUndir[i])->GetSize();z++)
	{
		if (m_iaConnectedBuffer[((CxIntArray*)ec->m_oaRelBondsUndir[i])->GetAt(z)] != 0)
			continue;
		REC_IsConnected(ec,((CxIntArray*)ec->m_oaRelBondsUndir[i])->GetAt(z));
	}
}


void CClusterAnalysis::RenderStepPOV(CTimeStep *ts, const char *s)
{
	FILE *b;
	CMolecule *m;
	CSingleMolecule *sm;
	CElement *el, *el2;
	int z, z2, z3, z4, o, o2;
	CxDVector3 vec1, vec2, vec3, vec1b, vec2b, vec3b, vecA, vecB, vecC, vecD;
	CMolBond *mb;
	CxDVector3 cam;
	CExtendedCluster *ec;
	double lb, bleach;
//	double cr, cg, cb;

	b = OpenFileWrite(s,true);

	mfprintf(b,"// Written by TRAVIS\n");
	mfprintf(b,"// See http://www.travis-analyzer.de\n\n");
	mfprintf(b,"#version 3.6;\n");
	mfprintf(b,"\n");

	mfprintf(b,"/**** Atoms ****/\n");

	mfprintf(b,"#declare atom_draw          = true;\n");
	mfprintf(b,"#declare atom_r             = 0.65;\n");
	mfprintf(b,"#declare atom_specular      = 0.7;\n");
	mfprintf(b,"#declare atom_reflect       = 0;\n");
	mfprintf(b,"#declare atom_ambient       = 0.2;\n");
	mfprintf(b,"#declare atom_diffuse       = 0.7;\n");
	mfprintf(b,"//#declare atom_color         = < 1.0, 1.0, 1.0, 0, 0 >;\n");
	mfprintf(b,"//#declare atom_trans         = 0.7;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare atom_draw_halo1      = true;\n");
	mfprintf(b,"#declare atom_r_halo1         = 0.007;\n");
	mfprintf(b,"#declare atom_d_halo1         = 0.0125;\n");
	mfprintf(b,"#declare atom_color_halo1     = < 0, 0, 0, 0, 0 >;\n");
	mfprintf(b,"#declare atom_specular_halo1  = 0;\n");
	mfprintf(b,"#declare atom_reflect_halo1   = 0;\n");
	mfprintf(b,"#declare atom_ambient_halo1   = 1.0;\n");
	mfprintf(b,"#declare atom_diffuse_halo1   = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"/**** Bonds ****/\n");

	mfprintf(b,"#declare bond_draw        = true;\n");
	mfprintf(b,"#declare bond_r           = 0.015 * 0.8;\n");
	mfprintf(b,"#declare bond_cluster_r   = 0.015 * 0.3;\n");
	mfprintf(b,"#declare bond_specular    = 0.7;\n");
	mfprintf(b,"#declare bond_reflect     = 0;\n");
	mfprintf(b,"#declare bond_ambient     = 0.2;\n");
	mfprintf(b,"#declare bond_diffuse     = 0.7;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare bond_draw_halo1      = true;\n");
	mfprintf(b,"#declare bond_r_halo1         = 0.007;\n");
	mfprintf(b,"#declare bond_d_halo1         = 0;\n");
	mfprintf(b,"#declare bond_color_halo1     = < 0, 0, 0, 0, 0 >;\n");
	mfprintf(b,"#declare bond_specular_halo1  = 0;\n");
	mfprintf(b,"#declare bond_reflect_halo1   = 0;\n");
	mfprintf(b,"#declare bond_ambient_halo1   = 1.0;\n");
	mfprintf(b,"#declare bond_diffuse_halo1   = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare bond_draw_stub       = true;\n");
	mfprintf(b,"#declare bond_r_stub          = 0.004;\n");
	mfprintf(b,"#declare bond_l_stub          = 0.004;\n");
	mfprintf(b,"#declare bond_color_stub      = < 0, 0, 0, 0, 0 >;\n");
	mfprintf(b,"#declare bond_specular_stub   = 0;\n");
	mfprintf(b,"#declare bond_reflect_stub    = 0;\n");
	mfprintf(b,"#declare bond_ambient_stub    = 1.0;\n");
	mfprintf(b,"#declare bond_diffuse_stub    = 0;\n");

	mfprintf(b,"\n");

	bleach = 0.75;

	mfprintf(b,"/**** Element Colors ****/\n");
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		if (z == g_iVirtAtomType)
			continue;
		el = ((CAtom*)g_oaAtoms[z])->m_pElement;
		mfprintf(b,"#declare elem_%s_col   = < %f, %f, %f, 0, 0 >;\n",el->m_sLabel,el->m_iColorR/255.0*(1.0-bleach)+bleach,el->m_iColorG/255.0*(1.0-bleach)+bleach,el->m_iColorB/255.0*(1.0-bleach)+bleach);
	}

	mfprintf(b,"\n/**** Element Radii ****/\n");
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		if (z == g_iVirtAtomType)
			continue;
		el = ((CAtom*)g_oaAtoms[z])->m_pElement;
		mfprintf(b,"#declare elem_%s_r  = %f * 0.8;\n",el->m_sLabel,el->m_fRadius/1000.0);
	}
	mfprintf(b,"\n");

//	cam[0] = 3.0;
//	cam[1] = 2.0;
//	cam[2] = 8.0;

	lb = g_fBoxX;
	if (g_fBoxY > lb)
		lb = g_fBoxY;
	if (g_fBoxZ > lb)
		lb = g_fBoxZ;

	cam[0] = 8.544 * lb/1566.0 * sin(m_fPOVAngle);
	cam[1] = 2.0 * lb/1566.0;
	cam[2] = 8.544 * lb/1566.0 * cos(m_fPOVAngle);

	mfprintf(b,"global_settings {\n");
	mfprintf(b,"  assumed_gamma 1\n");
	mfprintf(b,"/*  radiosity {\n");
	mfprintf(b,"    pretrace_start 0.08\n");
	mfprintf(b,"    pretrace_end   0.04\n");
	mfprintf(b,"    count 100\n\n");
	mfprintf(b,"    nearest_count 5\n");
	mfprintf(b,"    error_bound 0.4\n");
	mfprintf(b,"    recursion_limit 1\n\n");
	mfprintf(b,"    low_error_factor .5\n");
	mfprintf(b,"    gray_threshold 0.0\n");
	mfprintf(b,"    minimum_reuse 0.015\n");
	mfprintf(b,"    brightness 1\n\n");
	mfprintf(b,"    adc_bailout 0.01/2\n");
	mfprintf(b,"  }*/\n");
	mfprintf(b,"}\n");

	mfprintf(b,"\ncamera {\n");
	mfprintf(b,"	location <%f, %f, %f>\n",cam[0],cam[1],cam[2]);
	mfprintf(b,"	sky y\n");
	mfprintf(b,"	right -0.3*x*image_width/image_height\n");
	mfprintf(b,"	up 0.3*y\n");
	mfprintf(b,"	look_at <0, 0, 0>\n");
	mfprintf(b,"}\n");
	mfprintf(b,"\n");

	mfprintf(b,"// Solid background\n");
	mfprintf(b,"background { rgb < 1.0, 1.0, 1.0 > }\n");
	mfprintf(b,"\n");

	mfprintf(b,"// Gradient background\n");
	mfprintf(b,"/*sky_sphere {\n");
	mfprintf(b,"  pigment {\n");
	mfprintf(b,"    gradient y\n");
	mfprintf(b,"    color_map {\n");
	mfprintf(b,"      [ 0 color rgb < 0.05, 0.05, 0.05 > ]\n");
	mfprintf(b,"      [ 1 color rgb < 0.20, 0.16, 0.50 > ]\n");
	mfprintf(b,"    }\n");
	mfprintf(b,"    scale 0.1\n");
	mfprintf(b,"    translate -0.05\n");
	mfprintf(b,"  }\n");
	mfprintf(b,"}*/\n\n");

	mfprintf(b,"/**** Invisible, only for Radiosity ****/\n");
	mfprintf(b,"sphere {\n");
	mfprintf(b,"  <0, 0, 0>, 1\n");
	mfprintf(b,"  texture {\n");
	mfprintf(b,"    pigment {color rgb < 1.0, 1.0, 1.0 > }\n");
	mfprintf(b,"    finish { diffuse 0 ambient 1 }\n");
	mfprintf(b,"  }\n");
	mfprintf(b,"  hollow on\n");
	mfprintf(b,"  no_shadow\n");
	mfprintf(b,"  no_image\n");
	mfprintf(b,"  scale 30000\n");
	mfprintf(b,"}\n\n");

//	mfprintf(b,"light_source { < -8, 20, 20 > color rgb 0.8 }\n");
	mfprintf(b,"light_source { < %f, 20, %f > color rgb 0.8 }\n",21.54*sin(m_fPOVAngle-0.3805),21.54*cos(m_fPOVAngle-0.3805));
	mfprintf(b,"//light_source { < 25, 12, 20 > color rgb 0.5 }\n\n");

	mfprintf(b,"union {\n");
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n",-g_fBoxX/2000.0,-g_fBoxY/2000.0,-g_fBoxZ/2000.0, g_fBoxX/2000.0,-g_fBoxY/2000.0,-g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n", g_fBoxX/2000.0,-g_fBoxY/2000.0,-g_fBoxZ/2000.0, g_fBoxX/2000.0, g_fBoxY/2000.0,-g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n", g_fBoxX/2000.0, g_fBoxY/2000.0,-g_fBoxZ/2000.0,-g_fBoxX/2000.0, g_fBoxY/2000.0,-g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n",-g_fBoxX/2000.0, g_fBoxY/2000.0,-g_fBoxZ/2000.0,-g_fBoxX/2000.0,-g_fBoxY/2000.0,-g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n",-g_fBoxX/2000.0,-g_fBoxY/2000.0, g_fBoxZ/2000.0, g_fBoxX/2000.0,-g_fBoxY/2000.0, g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n", g_fBoxX/2000.0,-g_fBoxY/2000.0, g_fBoxZ/2000.0, g_fBoxX/2000.0, g_fBoxY/2000.0, g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n", g_fBoxX/2000.0, g_fBoxY/2000.0, g_fBoxZ/2000.0,-g_fBoxX/2000.0, g_fBoxY/2000.0, g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n",-g_fBoxX/2000.0, g_fBoxY/2000.0, g_fBoxZ/2000.0,-g_fBoxX/2000.0,-g_fBoxY/2000.0, g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n",-g_fBoxX/2000.0,-g_fBoxY/2000.0,-g_fBoxZ/2000.0,-g_fBoxX/2000.0,-g_fBoxY/2000.0, g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n", g_fBoxX/2000.0,-g_fBoxY/2000.0,-g_fBoxZ/2000.0, g_fBoxX/2000.0,-g_fBoxY/2000.0, g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n", g_fBoxX/2000.0, g_fBoxY/2000.0,-g_fBoxZ/2000.0, g_fBoxX/2000.0, g_fBoxY/2000.0, g_fBoxZ/2000.0);
	mfprintf(b,"  cylinder { < %f, %f, %f >, < %f, %f, %f >, 0.01 open }\n",-g_fBoxX/2000.0, g_fBoxY/2000.0,-g_fBoxZ/2000.0,-g_fBoxX/2000.0, g_fBoxY/2000.0, g_fBoxZ/2000.0);
	mfprintf(b,"  no_shadow pigment { rgbft < 0, 0, 1, 0, 0> }\n  finish { reflection 0 specular 0.7 ambient 0.2 diffuse 0.7 }\n}\n\n");

	mfprintf(b,"#macro m_atom_color(col)\n");
	mfprintf(b,"  #ifdef(atom_color)\n");
	mfprintf(b,"    atom_color\n");
	mfprintf(b,"  #else #if (defined(atom_trans))\n");
	mfprintf(b,"    col + < 0, 0, 0, 0, atom_trans >\n");
	mfprintf(b,"  #else\n");
	mfprintf(b,"    col\n");
	mfprintf(b,"  #end #end\n");
	mfprintf(b,"#end\n");

	mfprintf(b,"\nunion {\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];

		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];

			for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
			{
				if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
					continue;

				el = ((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_pElement;

				for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
				{
					o = ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4);
					vec1 = ts->m_vaCoords[o];
					vec1 /= 1000.0;

		//			mfprintf(b,"#if (atom_draw)\n");
					mfprintf(b,"  sphere { <%g, %g, %g>, elem_%s_r*atom_r\n",vec1[0],vec1[1],vec1[2],el->m_sLabel);

					if (m_iaPOVSMCluster[m->m_laSingleMolIndex[z2]] != -1)
					{
						mfprintf(b,"    pigment { rgbft < %f, %f, %f, 0, 0 > } ",((m_iaPOVSMColor[m->m_laSingleMolIndex[z2]]&0x00FF0000)>>16)/255.0,((m_iaPOVSMColor[m->m_laSingleMolIndex[z2]]&0x0000FF00)>>8)/255.0,(m_iaPOVSMColor[m->m_laSingleMolIndex[z2]]&0x000000FF)/255.0);
					} else mfprintf(b,"    pigment { rgbft m_atom_color(elem_%s_col) } ",el->m_sLabel);

					mfprintf(b,"finish { reflection atom_reflect specular atom_specular ambient atom_ambient diffuse atom_diffuse } }\n");
		//			mfprintf(b,"#end\n");

					vec2 = cam - vec1;
					vec2.Normalize();

	//				if (shadow)
					{
		//				mfprintf(b,"#if (atom_draw_halo1)\n");
						mfprintf(b,"  disc { < %g - (%g) * atom_d_halo1, %g - (%g) * atom_d_halo1, %g - (%g) * atom_d_halo1 >,\n",vec1[0],vec2[0],vec1[1],vec2[1],vec1[2],vec2[2]);
					//	mfprintf(b,"    < %g, %g, %g >, (elem_%s_r * atom_r) + atom_r_halo1, elem_%s_r * atom_r\n",vec2[0],vec2[1],vec2[2],el->m_sLabel,el->m_sLabel);
						mfprintf(b,"    < %g, %g, %g >, (elem_%s_r * atom_r) + atom_r_halo1\n",vec2[0],vec2[1],vec2[2],el->m_sLabel);
						mfprintf(b,"    pigment { rgbft atom_color_halo1 } finish { reflection atom_reflect_halo1 specular atom_specular_halo1 ambient atom_ambient_halo1 diffuse atom_diffuse_halo1 } no_reflection no_radiosity }\n");
		//				mfprintf(b,"#end\n");
					}
				}
			}

			for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
			{
				mb = (CMolBond*)sm->m_oaBonds[z3];
				o = mb->m_iAtomOffset[0];
				o2 = mb->m_iAtomOffset[1];

				el = ((CAtom*)g_oaAtoms[g_waAtomRealElement[o]])->m_pElement;
				el2 = ((CAtom*)g_oaAtoms[g_waAtomRealElement[o2]])->m_pElement;

				vec1 = ts->m_vaCoords[o];
				vec1 /= 1000.0;

				vec2 = ts->m_vaCoords[o2];
				vec2 /= 1000.0;

				if ((vec1-vec2).GetLength() > 0.3)
					continue;

				vec3 = (vec1/el->m_fRadius + vec2/el2->m_fRadius) / (1.0/el->m_fRadius+1.0/el2->m_fRadius);

				vec3b = vec2 - vec1;
				vec3b.Normalize();

				mfprintf(b,"#local vec1c = < %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ),\n",vec1[0],vec3b[0],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"  %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ),\n",vec1[1],vec3b[1],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"  %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ) >;\n",vec1[2],vec3b[2],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"#local vec2c = < %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ),\n",vec2[0],vec3b[0],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);
				mfprintf(b,"  %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ),\n",vec2[1],vec3b[1],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);
				mfprintf(b,"  %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ) >;\n",vec2[2],vec3b[2],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);

			//	mfprintf(b,"#if (bond_draw)\n");
				if (m_iaPOVSMCluster[m->m_laSingleMolIndex[z2]] != -1)
				{
					mfprintf(b,"  cylinder { vec1c, vec2c, bond_r open\n");
					mfprintf(b,"    pigment { rgbft < %f, %f, %f, 0, 0 > } finish { reflection bond_reflect specular bond_specular ambient bond_ambient diffuse bond_diffuse } }\n",((m_iaPOVSMColor[m->m_laSingleMolIndex[z2]]&0x00FF0000)>>16)/255.0,((m_iaPOVSMColor[m->m_laSingleMolIndex[z2]]&0x0000FF00)>>8)/255.0,(m_iaPOVSMColor[m->m_laSingleMolIndex[z2]]&0x000000FF)/255.0);
				} else
				{
					mfprintf(b,"  cylinder { < %g, %g, %g >, vec1c, bond_r open\n",vec3[0],vec3[1],vec3[2]);
					mfprintf(b,"    pigment { rgbft m_atom_color(elem_%s_col) } finish { reflection bond_reflect specular bond_specular ambient bond_ambient diffuse bond_diffuse } }\n",el->m_sLabel);
					mfprintf(b,"  cylinder { < %g, %g, %g >, vec2c, bond_r open\n",vec3[0],vec3[1],vec3[2]);
					mfprintf(b,"    pigment { rgbft m_atom_color(elem_%s_col) } finish { reflection bond_reflect specular bond_specular ambient bond_ambient diffuse bond_diffuse } }\n",el2->m_sLabel);
				}
			//	mfprintf(b,"#end\n");

				vec3 = cam - (vec1 + vec2) / 2.0;
				vec3.Normalize();

				vec2b = vec2-vec1;
				vec1b = CrossP(vec3,vec2b);
				vec1b.Normalize();

				mfprintf(b,"#local vec3  = < %g, %g, %g >;\n",vec3[0],vec3[1],vec3[2]);
				mfprintf(b,"#local vec1b = < %g, %g, %g >;\n",vec1b[0],vec1b[1],vec1b[2]);

//				if (shadow)
				{
			//		mfprintf(b,"#if (bond_draw_halo1)\n");
					mfprintf(b,"  #local vecA = vec1c + vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecB = vec2c + vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecC = vec2c + vec1b * (bond_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecD = vec1c + vec1b * (bond_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  #local vecA = vec1c - vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecB = vec2c - vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecC = vec2c - vec1b * (bond_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecD = vec1c - vec1b * (bond_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
			//		mfprintf(b,"#end\n");
				}

//				if (shadow)
				{
					vec3 = vec2 - vec1;
					vec3.Normalize();

			//		mfprintf(b,"#if (bond_draw_stub)\n");
					mfprintf(b,"  #local vec3  = < %g, %g, %g >;\n",vec3[0],vec3[1],vec3[2]);
					mfprintf(b,"  difference {\n");
					mfprintf(b,"    cylinder { vec1c, vec1c + vec3 * bond_l_stub, bond_r + bond_r_stub }\n");
					mfprintf(b,"    cylinder { vec1c - vec3 * (bond_l_stub+0.001), vec1c + vec3 * (bond_l_stub+0.001), bond_r }\n");
					mfprintf(b,"    pigment { rgbft bond_color_stub } finish { reflection bond_reflect_stub specular bond_specular_stub ambient bond_ambient_stub diffuse bond_diffuse_stub } }\n");
					mfprintf(b,"  difference {\n");
					mfprintf(b,"    cylinder { vec2c, vec2c - vec3 * bond_l_stub, bond_r + bond_r_stub }\n");
					mfprintf(b,"    cylinder { vec2c + vec3 * (bond_l_stub+0.001), vec2c - vec3 * (bond_l_stub+0.001), bond_r }\n");
					mfprintf(b,"    pigment { rgbft bond_color_stub } finish { reflection bond_reflect_stub specular bond_specular_stub ambient bond_ambient_stub diffuse bond_diffuse_stub } }\n");
			//		mfprintf(b,"#end\n");
				}
			}

		}
	}

	for (z=0;z<m_oaPOVExtendedCluster.GetSize();z++)
	{
		ec = (CExtendedCluster*)m_oaPOVExtendedCluster[z];

		for (z2=0;z2<ec->m_iaAtoms.GetSize();z2++)
		{
			for (z3=0;z3<((CxIntArray*)ec->m_oaBonds[z2])->GetSize();z3++)
			{
				if (g_laAtomSMIndex[ec->m_iaAtoms[z2]] == g_laAtomSMIndex[((CxIntArray*)ec->m_oaBonds[z2])->GetAt(z3)])
					continue;

				if (ec->m_iaAtoms[z2] > ((CxIntArray*)ec->m_oaBonds[z2])->GetAt(z3))
					continue;

				o = ec->m_iaAtoms[z2];
				o2 = ((CxIntArray*)ec->m_oaBonds[z2])->GetAt(z3);

				el = ((CAtom*)g_oaAtoms[g_waAtomRealElement[o]])->m_pElement;
				el2 = ((CAtom*)g_oaAtoms[g_waAtomRealElement[o2]])->m_pElement;

				vec1 = ts->m_vaCoords[o];
				vec1 /= 1000.0;

				vec2 = ts->m_vaCoords[o2];
				vec2 /= 1000.0;

		//		if ((vec1-vec2).GetLength() > 0.3)
		//			continue;

				// Not used ??
				//vec3 = (vec1/el->m_fRadius + vec2/el2->m_fRadius) / (1.0/el->m_fRadius+1.0/el2->m_fRadius);

				vec3b = vec2 - vec1;
				vec3b.Normalize();

				mfprintf(b,"#local vec1c = < %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_cluster_r*bond_cluster_r ) ) ),\n",vec1[0],vec3b[0],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"  %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_cluster_r*bond_cluster_r ) ) ),\n",vec1[1],vec3b[1],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"  %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_cluster_r*bond_cluster_r ) ) ) >;\n",vec1[2],vec3b[2],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"#local vec2c = < %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_cluster_r*bond_cluster_r ) ) ),\n",vec2[0],vec3b[0],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);
				mfprintf(b,"  %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_cluster_r*bond_cluster_r ) ) ),\n",vec2[1],vec3b[1],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);
				mfprintf(b,"  %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_cluster_r*bond_cluster_r ) ) ) >;\n",vec2[2],vec3b[2],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);

		//		mfprintf(b,"#if (bond_draw)\n");
				mfprintf(b,"  cylinder { vec1c, vec2c, bond_cluster_r open\n");
				mfprintf(b,"    pigment { rgbft < %f, %f, %f, 0, 0 > } finish { reflection bond_reflect specular bond_specular ambient bond_ambient diffuse bond_diffuse } }\n",((m_iaPOVSMColor[g_laAtomSMIndex[ec->m_iaAtoms[z2]]]&0x00FF0000)>>16)/255.0,((m_iaPOVSMColor[g_laAtomSMIndex[ec->m_iaAtoms[z2]]]&0x0000FF00)>>8)/255.0,(m_iaPOVSMColor[g_laAtomSMIndex[ec->m_iaAtoms[z2]]]&0x000000FF)/255.0);
		//		mfprintf(b,"#end\n");

				vec3 = cam - (vec1 + vec2) / 2.0;
				vec3.Normalize();

				vec2b = vec2-vec1;
				vec1b = CrossP(vec3,vec2b);
				vec1b.Normalize();

				mfprintf(b,"#local vec3  = < %g, %g, %g >;\n",vec3[0],vec3[1],vec3[2]);
				mfprintf(b,"#local vec1b = < %g, %g, %g >;\n",vec1b[0],vec1b[1],vec1b[2]);

	//				if (shadow)
				{
			//		mfprintf(b,"#if (bond_draw_halo1)\n");
					mfprintf(b,"  #local vecA = vec1c + vec1b * (bond_cluster_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecB = vec2c + vec1b * (bond_cluster_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecC = vec2c + vec1b * (bond_cluster_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecD = vec1c + vec1b * (bond_cluster_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  #local vecA = vec1c - vec1b * (bond_cluster_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecB = vec2c - vec1b * (bond_cluster_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecC = vec2c - vec1b * (bond_cluster_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecD = vec1c - vec1b * (bond_cluster_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
			//		mfprintf(b,"#end\n");
				}

	//				if (shadow)
				{
					vec3 = vec2 - vec1;
					vec3.Normalize();

				//	mfprintf(b,"#if (bond_draw_stub)\n");
					mfprintf(b,"  #local vec3  = < %g, %g, %g >;\n",vec3[0],vec3[1],vec3[2]);
					mfprintf(b,"  difference {\n");
					mfprintf(b,"    cylinder { vec1c, vec1c + vec3 * bond_l_stub, bond_cluster_r + bond_r_stub }\n");
					mfprintf(b,"    cylinder { vec1c - vec3 * (bond_l_stub+0.001), vec1c + vec3 * (bond_l_stub+0.001), bond_cluster_r }\n");
					mfprintf(b,"    pigment { rgbft bond_color_stub } finish { reflection bond_reflect_stub specular bond_specular_stub ambient bond_ambient_stub diffuse bond_diffuse_stub } }\n");
					mfprintf(b,"  difference {\n");
					mfprintf(b,"    cylinder { vec2c, vec2c - vec3 * bond_l_stub, bond_cluster_r + bond_r_stub }\n");
					mfprintf(b,"    cylinder { vec2c + vec3 * (bond_l_stub+0.001), vec2c - vec3 * (bond_l_stub+0.001), bond_cluster_r }\n");
					mfprintf(b,"    pigment { rgbft bond_color_stub } finish { reflection bond_reflect_stub specular bond_specular_stub ambient bond_ambient_stub diffuse bond_diffuse_stub } }\n");
				//	mfprintf(b,"#end\n");
				}
			}
		}
	}

	mfprintf(b,"\n  no_shadow\n}\n\n");

	fclose(b);
}


void CClusterAnalysis::WriteClusterXYZ(CClusterNode *n, CTimeStep *ts, FILE *a)
{
	int z2;
	double avgx, avgy, avgz;

	mfprintf(a,"%d\nStep %lu\n",n->m_iaAtoms.GetSize(),g_iSteps);

	avgx = 0;
	avgy = 0;
	avgz = 0;
	for (z2=0;z2<n->m_iaAtoms.GetSize();z2++)
	{
		if (z2 != 0)
		{
//			mprintf("\nOrig: (%f|%f|%f), Avg: (%f|%f|%f)",ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][0],ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][1],ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][2],avgx/z2, avgy/z2, avgz/z2);
			while (ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][0]-(avgx/z2) > g_fBoxX/2)
				ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][0] -= g_fBoxX;
			while (ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][0]-(avgx/z2) < -g_fBoxX/2)
				ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][0] += g_fBoxX;
			while (ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][1]-(avgy/z2) > g_fBoxY/2)
				ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][1] -= g_fBoxY;
			while (ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][1]-(avgy/z2) < -g_fBoxY/2)
				ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][1] += g_fBoxY;
			while (ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][2]-(avgz/z2) > g_fBoxZ/2)
				ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][2] -= g_fBoxZ;
			while (ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][2]-(avgz/z2) < -g_fBoxZ/2)
				ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][2] += g_fBoxZ;
//			mprintf("\n  New: (%f|%f|%f)",ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][0],ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][1],ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][2]);
		}
		avgx += ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][0];
		avgy += ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][1];
		avgz += ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][2];
	}
	avgx /= n->m_iaAtoms.GetSize();
	avgy /= n->m_iaAtoms.GetSize();
	avgz /= n->m_iaAtoms.GetSize();

	for (z2=0;z2<n->m_iaAtoms.GetSize();z2++)
		mfprintf(a,"%s  %f  %f  %f\n",(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[m_iaAtomList[n->m_iaAtoms[z2]]]])->m_sName,(ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][0]-avgx)/100.0,(ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][1]-avgy)/100.0,(ts->m_vaCoords[m_iaAtomList[n->m_iaAtoms[z2]]][2]-avgz)/100.0);

	fflush(a);
}


unsigned int CClusterAnalysis::POVColorFunction(int i)
{
	switch(i % 10)
	{
		case 0: return 0x000000E0; // Dunkelblau
		case 1: return 0x00009000; // Dunkelgruen
		case 2: return 0x00C00000; // Dunkelrot
		case 3: return 0x00FFFF00; // Gelb
		case 4: return 0x00FF00FF; // Violett
		case 5: return 0x0000FFFF; // Cyan
		case 6: return 0x005050FF; // Hellblau
		case 7: return 0x0060FF60; // Hellgruen
		case 8: return 0x00FF3060; // Hellrot
		case 9: return 0x00E07000; // Braun
	}

	return 0;
}



