/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Sascha Gehrke.

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

#include "luzar.h"

#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_luzar(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_luzar() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


CLuzarCorrelation::CLuzarCorrelation()
{                      
    m_iInput = 0;  
    m_iDepth = 0;  
    m_pFFT = NULL; 
    m_pFFT2 = NULL;
} 


void CLuzarCorrelation::Init(int input, int depth)                                                                                                                 
{                                                                                                               
    m_iInput = input;                                                                                       
    m_iDepth = depth;                                                                                       
	m_iFFTSize = CalcSize(input);                                                          
	
	try { m_pFFT = new CFFT(); } catch(...) { m_pFFT = NULL; }                                     
	if (m_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);  
	
	m_pFFT->PrepareFFT_C2C(2*m_iFFTSize);                                                          
	
	try { m_pFFT2 = new CFFT(); } catch(...) { m_pFFT2 = NULL; }                                   
	if (m_pFFT2 == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__); 
	
	m_pFFT2->PrepareFFT_C2C(2*m_iFFTSize);                                                         
	
	try { m_pFFTback = new CFFT(); } catch(...) { m_pFFTback = NULL; }                             
	if (m_pFFTback == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__); 
	
	m_pFFTback->PrepareInverseFFT_C2C(2*m_iFFTSize);    
}	


int CLuzarCorrelation::CalcSize(int n)  // Calculate the next bigger  number which is factorable by twos, threes, and fives
{					                                                                    
	int m;

    while( n > 0 )	// Should be endless ... if not, you have a no-step trajectorie 
	{
        m = n;
        while ( (m%2) == 0 ) 
		m/=2;
        while ( (m%3) == 0 ) 	
		m/=3;
        while ( (m%5) == 0 ) 
		m/=5;
        if ( m == 1 )
            break; 
        n++;
	}
	return n;
}


void CLuzarCorrelation::Correlate(CxDoubleArray* input1, CxDoubleArray* input2, CxDoubleArray* output1, CxDoubleArray* output2)
{
	int i;	

	output1->SetSize(m_iDepth);
	output2->SetSize(m_iDepth);

	for ( i=0; i < m_iInput; i++ )	// Fill Input-Strings with values and zeros
	{
		m_pFFT->m_pInput[i*2] = (*input1)[i];
		m_pFFT->m_pInput[i*2+1] = 0;
        m_pFFT2->m_pInput[i*2] = (*input2)[i];
        m_pFFT2->m_pInput[i*2+1] = 0;
    }

	for ( i=m_iInput; i < 2*m_iFFTSize; i++ )	// Fill with zeros
	{
		m_pFFT->m_pInput[i*2] = 0;
		m_pFFT->m_pInput[i*2+1] = 0;
        m_pFFT2->m_pInput[i*2] = 0;
        m_pFFT2->m_pInput[i*2+1] = 0;
    }


	m_pFFT->DoFFT();		// Move into Fourier-space
	m_pFFT2->DoFFT();

	for ( i=0; i < 2*m_iFFTSize; i++ )	// Autocorrelate input1
	{
		m_pFFTback->m_pInput[i*2] = (double)((m_pFFT->m_pOutput[i*2]*m_pFFT->m_pOutput[i*2] + m_pFFT->m_pOutput[i*2+1]*m_pFFT->m_pOutput[i*2+1]));
        m_pFFTback->m_pInput[i*2+1] = 0;
	}

	m_pFFTback->DoFFT();		// Leaving Fourier-space 
	
	for ( i=0; i < m_iDepth; i++ )	// Fill autocorrelated function in output1
	{
		(*output1)[i] = (double)((double)m_pFFTback->m_pOutput[2*i] / ((double)m_iInput - i));    
	}

	for ( i=0; i < 2*m_iFFTSize; i++ )	// Crosscorrelate input1 and input2
	{
		m_pFFTback->m_pInput[i*2] = (double)((m_pFFT->m_pOutput[i*2]*m_pFFT2->m_pOutput[i*2]   + m_pFFT->m_pOutput[i*2+1]*m_pFFT2->m_pOutput[i*2+1]));	
		m_pFFTback->m_pInput[i*2+1] = (double)((m_pFFT->m_pOutput[i*2]*m_pFFT2->m_pOutput[i*2+1] - m_pFFT->m_pOutput[i*2+1]*m_pFFT2->m_pOutput[i*2]));
	}

	m_pFFTback->DoFFT();		// Leaving Fourier-space 

	for ( i=0; i < m_iDepth; i++ )	// Fill crosscorrelated function in output2
	{
		(*output2)[i] = (double)((double)m_pFFTback->m_pOutput[2*i] / ((double)m_iInput - i));   
	}
}


CLuzarCorrelation::~CLuzarCorrelation()
{
}


void CLuzarCorrelation::Fit(CxDoubleArray* c, CxDoubleArray* n, CxDoubleArray* k, CxDoubleArray* a, CxDoubleArray* b, int trunc)
{
	double RMS_alt = 1e25, conv = 1e-25, a_alt = 0.3, b_alt = 0.7, RMSa, RMSb, RMS;
	int i, j = 0;
	double multiplier, count, count_start = 1e7;
	bool repeat = true;
	double ps_corr = 1000 / (g_fTimestepLength * g_iStride);

	for ( i=1; i<c->GetSize(); i++ )
		k->Add( c->GetAt(i) - c->GetAt(i-1) );
 
	while ( repeat )
	{
		repeat = false;
		
		//mprintf("       k_f[1/ps]       k_b[1/ps]       tau_f[ps]       tau_b[ps]\n");
		mprintf("       k_f / ps^-1     k_b / ps^-1     tau_f / ps      tau_b / ps\n");
		mprintf("       ----------------------------------------------------------\n");

		a->Add(0.3);
		b->Add(0.7);

	        RMS = GetRMS((*a)[0], (*b)[0], c, n, k, trunc);
	        multiplier = 0.5;
			count = count_start;
	        while ( multiplier > 0.000000000000009 )
	        {
			while ( RMS_alt - RMS > conv )
			{
				a_alt = (*a)[j];
				b_alt = (*b)[j];
				RMS_alt = RMS;
				
				RMSa = GetRMS((*a)[j]*0.999999, (*b)[j], c, n, k, trunc);
				RMSb = GetRMS((*a)[j],(*b)[j]*0.999999, c, n, k, trunc);
				
				if ( RMS > RMSa )
				    (*a)[j] *= (1 - multiplier);
				else if ( RMS < RMSa )
				    (*a)[j] *= (1 + multiplier);
				
				if ( RMS > RMSb )
				    (*b)[j] *= (1 - multiplier);
				else if ( RMS < RMSb )
				    (*b)[j] *= (1 + multiplier);
				
				RMS = GetRMS((*a)[j], (*b)[j], c, n, k, trunc);
				
				if ( count == 0 )
				{
					multiplier = 0;
					break;
				}
				count--;
			}
			RMS_alt = 1e25;
			(*a)[j] = a_alt;
			(*b)[j] = b_alt;
			multiplier /= 5;
	        }
		mprintf("%15.10f %15.10f %15.5f %15.5f  ", (*a)[j]*ps_corr,(*b)[j]*ps_corr,1/(*a)[j]/ps_corr,1/(*b)[j]/ps_corr);
		if ( count > 0 )
			mprintf("\n");
		else
		{
			mprintf("not converged!!!\n");
			repeat = true;
		}

		mprintf("\n\n");

		if ( repeat )
		{
			count_start *= 100;
			repeat = AskYesNo("    Fitting is not converged. Retry with more steps (y/n)? [yes] ",true);
			mprintf("\n\n");
		}
	}
}


double CLuzarCorrelation::GetRMS(double a, double b, CxDoubleArray* c, CxDoubleArray* n, CxDoubleArray* k, int trunc)
{
    double RMS = 0.0;
    int i;

    for ( i=trunc; i<k->GetSize(); i++ )
    {
        RMS += pow2((a * c->GetAt(i) - b * n->GetAt(i) + k->GetAt(i)));
    }

    return sqrt( RMS / (k->GetSize()-trunc) );
}
