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

#include "bqb_interface.h"
#include "bqb_driver.h"
#include <stdarg.h>


const char *GetRevisionInfo_bqb_interface(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_interface() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



#define BQB_BUF_LEN 16384



/*******************************************************************************************
******      Global Functions      **********************************************************
*******************************************************************************************/



unsigned long g_iBQBLastError = BQB_OK;



CBQBInterface* BQBCreateInterface(unsigned long flags) {

	UNUSED(flags);

	CBQBInterface *p;

	p = new CBQBInterface();

	g_iBQBLastError = BQB_OK;

	return p;
}



bool BQBReleaseInterface(CBQBInterface *i) {

	if (!i->m_bActive) {
		printf("BQBReleaseInterface(): Error: Interface has already been released before.\n");
		printf("\n");
		abort();
	}

	i->m_bActive = false;

	delete i;

	g_iBQBLastError = BQB_OK;

	return true;
}



unsigned long BQBGetLastError() {
	return g_iBQBLastError;
}



/*******************************************************************************************
******      CBQBInterface      *************************************************************
*******************************************************************************************/



CBQBInterface::CBQBInterface() {

	m_pPrintCallback = NULL;
	m_pPrintCallbackVar = NULL;
	m_pEPrintCallback = NULL;
	m_pEPrintCallbackVar = NULL;
	m_pBPrintCallback = NULL;
	m_pBPrintCallbackVar = NULL;
	m_pPrintFile = NULL;
	m_pEPrintFile = NULL;
	m_pBPrintFile = NULL;

	m_iPrintLevel = BQB_PL_STANDARD;

	m_bActive = true;

	m_iLastError = 0;
}



CBQBInterface::~CBQBInterface() {

	unsigned long z;

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::~CBQBInterface(): Destructing instance ...\n");

	if (m_bActive) {
		this->eprintf("CBQBInterface::~CBQBInterface(): Error: Cannot manually destruct instance of this class.\n");
		this->printf("Use BQBReleaseInterface() for this purpose.\n");
		this->printf("\n");
		abort();
	}

	z = 0;
	std::set<CBQBFileWriter*>::iterator itfw;
	for ( itfw=m_oaListFileWriter.begin(); itfw!=m_oaListFileWriter.end(); ++itfw ) {
		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::~CBQBInterface(): Deleting CBQBFileWriter %lu/%lu ...\n",z+1,(unsigned long)m_oaListFileWriter.size());
		(*itfw)->m_bActive = false;
		delete (*itfw);
		z++;
	}

	z = 0;
	std::set<CBQBFileReader*>::iterator itfr;
	for ( itfr=m_oaListFileReader.begin(); itfr!=m_oaListFileReader.end(); ++itfr ) {
		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::~CBQBInterface(): Deleting CBQBFileReader %lu/%lu ...\n",z+1,(unsigned long)m_oaListFileReader.size());
		(*itfr)->m_bActive = false;
		delete (*itfr);
		z++;
	}

	z = 0;
	std::set<CBQBVolHeader*>::iterator itvh;
	for ( itvh=m_oaListVolHeader.begin(); itvh!=m_oaListVolHeader.end(); ++itvh ) {
		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::~CBQBInterface(): Deleting CBQBVolHeader %lu/%lu ...\n",z+1,(unsigned long)m_oaListVolHeader.size());
		(*itvh)->m_bActive = false;
		delete (*itvh);
		z++;
	}

	z = 0;
	std::set<CBQBAtomHeader*>::iterator itah;
	for ( itah=m_oaListAtomHeader.begin(); itah!=m_oaListAtomHeader.end(); ++itah ) {
		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::~CBQBInterface(): Deleting CBQBAtomHeader %lu/%lu ...\n",z+1,(unsigned long)m_oaListAtomHeader.size());
		(*itah)->m_bActive = false;
		delete (*itah);
		z++;
	}

	z = 0;
	std::set<CBQBCellInfo*>::iterator itci;
	for ( itci=m_oaListCellInfo.begin(); itci!=m_oaListCellInfo.end(); ++itci ) {
		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::~CBQBInterface(): Deleting CBQBCellInfo %lu/%lu ...\n",z+1,(unsigned long)m_oaListCellInfo.size());
		(*itci)->m_bActive = false;
		delete (*itci);
		z++;
	}

	z = 0;
	std::set<CBQBVolFrame*>::iterator itvf;
	for ( itvf=m_oaListVolFrame.begin(); itvf!=m_oaListVolFrame.end(); ++itvf ) {
		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::~CBQBInterface(): Deleting CBQBVolFrame %lu/%lu ...\n",z+1,(unsigned long)m_oaListVolFrame.size());
		(*itvf)->m_bActive = false;
		delete (*itvf);
		z++;
	}

	z = 0;
	std::set<CBQBAtomPosFrame*>::iterator itaf;
	for ( itaf=m_oaListAtomPosFrame.begin(); itaf!=m_oaListAtomPosFrame.end(); ++itaf ) {
		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::~CBQBInterface(): Deleting CBQBAtomPosFrame %lu/%lu ...\n",z+1,(unsigned long)m_oaListAtomPosFrame.size());
		(*itaf)->m_bActive = false;
		delete (*itaf);
		z++;
	}

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::~CBQBInterface(): Done.\n");
}



CBQBDriver* CBQBInterface::CreateDriver(unsigned long flags) {

	UNUSED(flags);

	CBQBDriver *p;

	p = new CBQBDriver(*this);

	return p;
}



bool CBQBInterface::DestroyDriver(CBQBDriver *d) {

	delete d;

	return true;
}



CBQBEngine* CBQBInterface::CreateEngine(unsigned long flags) {

	UNUSED(flags);

	CBQBEngine *p;

	p = new CBQBEngine(*this);

	return p;
}



bool CBQBInterface::DestroyEngine(CBQBEngine *e) {

	delete e;

	return true;
}



void CBQBInterface::printf(const char *s, ...) const {

	va_list args;
	static char buffer[BQB_BUF_LEN];


	if (strlen(s) >= BQB_BUF_LEN) {
		eprintf("CBQBInterface::printf(): Internal error: Buffer overflow (A).\n");
		abort();
	}

	buffer[BQB_BUF_LEN-1] = 0;

	va_start(args,s);
	vsprintf(buffer,s,args);
	va_end(args);

	if (buffer[BQB_BUF_LEN-1] != 0) {
		eprintf("CBQBInterface::printf(): Internal error: Buffer overflow (B).\n");
		abort();
	}

	if (m_pPrintCallbackVar != NULL)
		m_pPrintCallbackVar("%s",buffer);
	else if (m_pPrintCallback != NULL)
		m_pPrintCallback(buffer);
	else if (m_pPrintFile != NULL)
		::fprintf( m_pPrintFile, "%s", buffer );
	else
		::printf("%s",buffer);
}



void CBQBInterface::bprintf(const char *s, ...) const {

	va_list args;
	static char buffer[BQB_BUF_LEN];


	if (strlen(s) >= BQB_BUF_LEN) {
		eprintf("CBQBInterface::bprintf(): Internal error: Buffer overflow (A).\n");
		abort();
	}

	buffer[BQB_BUF_LEN-1] = 0;

	va_start(args,s);
	vsprintf(buffer,s,args);
	va_end(args);

	if (buffer[BQB_BUF_LEN-1] != 0) {
		eprintf("CBQBInterface::bprintf(): Internal error: Buffer overflow (B).\n");
		abort();
	}

	if (m_pBPrintCallbackVar != NULL)
		m_pBPrintCallbackVar("%s",buffer);
	else if (m_pBPrintCallback != NULL)
		m_pBPrintCallback(buffer);
	else if (m_pBPrintFile != NULL)
		::fprintf( m_pBPrintFile, "%s", buffer );
	else
		::printf("%s",buffer);
}



void CBQBInterface::eprintf(const char *s, ...) const {

	va_list args;
	static char buffer[BQB_BUF_LEN];


	if (strlen(s) >= BQB_BUF_LEN) {
		eprintf("CBQBInterface::eprintf(): Internal error: Buffer overflow (A).\n");
		abort();
	}

	buffer[BQB_BUF_LEN-1] = 0;

	va_start(args,s);
	vsprintf(buffer,s,args);
	va_end(args);

	if (buffer[BQB_BUF_LEN-1] != 0) {
		eprintf("CBQBInterface::eprintf(): Internal error: Buffer overflow (B).\n");
		abort();
	}

	if (m_pEPrintCallbackVar != NULL)
		m_pEPrintCallbackVar("%s",buffer);
	else if (m_pEPrintCallback != NULL)
		m_pEPrintCallback(buffer);
	else if (m_pEPrintFile != NULL)
		::fprintf( m_pEPrintFile, "%s", buffer );
	else
		::printf("%s",buffer);
}



void CBQBInterface::FlushLog() {

}



void CBQBInterface::SetPrintCallback( void (*fp)(const char *) ) {

	m_pPrintFile = NULL;
	m_pEPrintFile = NULL;
	m_pBPrintFile = NULL;
	m_pPrintCallbackVar = NULL;
	m_pEPrintCallbackVar = NULL;
	m_pBPrintCallbackVar = NULL;
	m_pPrintCallback = fp;
	m_pEPrintCallback = fp;
	m_pBPrintCallback = fp;
}



void CBQBInterface::SetPrintCallback( void (*fp)(const char *, ...) ) {

	m_pPrintFile = NULL;
	m_pEPrintFile = NULL;
	m_pBPrintFile = NULL;
	m_pPrintCallback = NULL;
	m_pEPrintCallback = NULL;
	m_pBPrintCallback = NULL;
	m_pPrintCallbackVar = fp;
	m_pEPrintCallbackVar = fp;
	m_pBPrintCallbackVar = fp;
}



void CBQBInterface::ResetPrint() {

	m_pPrintCallback = NULL;
	m_pEPrintCallback = NULL;
	m_pBPrintCallback = NULL;
	m_pPrintCallbackVar = NULL;
	m_pEPrintCallbackVar = NULL;
	m_pBPrintCallbackVar = NULL;
	m_pPrintFile = NULL;
	m_pEPrintFile = NULL;
	m_pBPrintFile = NULL;
}



void CBQBInterface::SetEPrintCallback( void (*fp)(const char *) ) {

	m_pEPrintFile = NULL;
	m_pEPrintCallbackVar = NULL;
	m_pEPrintCallback = fp;
}



void CBQBInterface::SetEPrintCallback( void (*fp)(const char *, ...) ) {

	m_pEPrintFile = NULL;
	m_pEPrintCallback = NULL;
	m_pEPrintCallbackVar = fp;
}



void CBQBInterface::ResetEPrint() {

	m_pEPrintCallback = NULL;
	m_pEPrintCallbackVar = NULL;
	m_pEPrintFile = NULL;
}



void CBQBInterface::SetBPrintCallback( void (*fp)(const char *) ) {

	m_pBPrintFile = NULL;
	m_pBPrintCallbackVar = NULL;
	m_pBPrintCallback = fp;
}



void CBQBInterface::SetBPrintCallback( void (*fp)(const char *, ...) ) {

	m_pBPrintFile = NULL;
	m_pBPrintCallback = NULL;
	m_pBPrintCallbackVar = fp;
}



void CBQBInterface::ResetBPrint() {

	m_pBPrintCallback = NULL;
	m_pBPrintCallbackVar = NULL;
	m_pBPrintFile = NULL;
}



void CBQBInterface::SetPrintFile( FILE *f ) {

	m_pPrintFile = f;
	m_pEPrintFile = f;
	m_pBPrintFile = f;
	m_pPrintCallback = NULL;
	m_pPrintCallbackVar = NULL;
	m_pEPrintCallback = NULL;
	m_pEPrintCallbackVar = NULL;
	m_pBPrintCallback = NULL;
	m_pBPrintCallbackVar = NULL;
}



void CBQBInterface::SetBPrintFile( FILE *f ) {

	m_pBPrintFile = f;
	m_pBPrintCallback = NULL;
	m_pBPrintCallbackVar = NULL;
}



void CBQBInterface::SetEPrintFile( FILE *f ) {
	
	m_pEPrintFile = f;
	m_pEPrintCallback = NULL;
	m_pEPrintCallbackVar = NULL;
}



void CBQBInterface::SetPrintLevel(int i) {

	m_iPrintLevel = i;
}



int CBQBInterface::GetPrintLevel() const {

	return m_iPrintLevel;
}



void CBQBInterface::SetLastError( unsigned long e ) {

	if (IsPL( BQB_PL_VERBOSE )) {
		if (e == BQB_OK)
			this->printf("CBQBInterface::SetLastError(): BQB_OK.\n");
		else
			this->printf("CBQBInterface::SetLastError(): %lu.\n",e);
	}

	m_iLastError = e;
}



CBQBFileWriter* CBQBInterface::CreateBQBFileWriter( const char *name, int history, unsigned long contents, unsigned long flags ) {

	CBQBFileWriter *fw;

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateBQBFileWriter(): Creating file writer (\"%s\", %d, %lu, %lu) ...\n",name,history,contents,flags);

	if (history < 1) {
		this->eprintf("CBQBInterface::CreateBQBFileWriter(): Error: History needs to be >= 1.\n");
		this->eprintf("\n");
		SetLastError( BQB_ERROR_PARAMETERS );
		return NULL;
	}

	fw = new CBQBFileWriter( *this );

	if (!fw->OpenFileWrite( name, history, contents, flags )) {
		if (IsPL( BQB_PL_STANDARD ))
			this->eprintf( "CBQBInterface::CreateBQBFileWriter(): Error: Could not open file \"%s\".\n", name );
		SetLastError( BQB_ERROR_FILE );
		return NULL;
	}

	m_oaListFileWriter.insert( fw );
	SetLastError( BQB_OK );

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateBQBFileWriter(): Done.\n");

	return fw;
}



CBQBFileReader* CBQBInterface::CreateBQBFileReader( const char *name, unsigned long flags ) {

	CBQBFileReader *fr;

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateBQBFileReader(): Creating file reader (\"%s\", %lu) ...\n",name,flags);

	fr = new CBQBFileReader( *this );
	if (!fr->OpenFileRead( name, flags )) {
		if (IsPL( BQB_PL_STANDARD ))
			this->eprintf( "CBQBInterface::CreateBQBFileReader(): Error: Could not open file \"%s\".\n", name );
		SetLastError( BQB_ERROR_FILE );
		return NULL;
	}

	m_oaListFileReader.insert( fr );
	SetLastError( BQB_OK );

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateBQBFileReader(): Done.\n");

	return fr;
}



CBQBVolHeader* CBQBInterface::CreateVolumetricHeader() {

	CBQBVolHeader *vh;

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateVolumetricHeader(): Creating volumetric header ...\n");

	vh = new CBQBVolHeader( *this );

	// Default Values
	vh->SetSigni( 5 );
	vh->SetEps( 12 );

	m_oaListVolHeader.insert( vh );

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateVolumetricHeader(): Done.\n");

	return vh;
}



CBQBAtomHeader* CBQBInterface::CreateAtomHeader() {

	CBQBAtomHeader *ah;

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateAtomHeader(): Creating atom header ...\n");

	ah = new CBQBAtomHeader( *this );

	// Default Values
	ah->SetPrecision( 6 );

	m_oaListAtomHeader.insert( ah );

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateAtomHeader(): Done.\n");

	return ah;
}



CBQBCellInfo* CBQBInterface::CreateCellInfo() {

	CBQBCellInfo *ci;

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateCellInfo(): Creating cell info ...\n");

	ci = new CBQBCellInfo( *this );
	m_oaListCellInfo.insert( ci );

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateCellInfo(): Done.\n");

	return ci;
}



CBQBVolFrame* CBQBInterface::CreateVolumetricFrame( const CBQBVolHeader *header, int id ) {

	CBQBVolFrame *vf;
	unsigned long z;


	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateVolumetricFrame(): Creating volumetric frame ...\n");

	if (&header->m_IF != this) {
		this->eprintf("CBQBInterface::CreateVolumetricFrame(): Error: Supplied CBQBVolHeader does not belong to this CBQBInterface.\n");
		this->eprintf("\n");
		abort();
	}

	if (!m_oaRecycledVolFrames.empty()) {

		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::CreateVolumetricFrame(): Recycling frame ...\n");

		if (m_oaRecycledVolFrames.back()->m_pVolHeader == header) {

			if (IsPL( BQB_PL_VERBOSE ))
				this->printf("CBQBInterface::CreateVolumetricFrame(): Using last frame (%lu/%lu) for recycling.\n",
					(unsigned long)m_oaRecycledVolFrames.size(), (unsigned long)m_oaRecycledVolFrames.size() );

			vf = m_oaRecycledVolFrames.back();
			m_oaRecycledVolFrames.pop_back();

			if (IsPL( BQB_PL_VERBOSE ))
				this->printf("CBQBInterface::CreateVolumetricFrame(): Done.\n");

			return vf;

		} else {

			for (z=0;z<m_oaRecycledVolFrames.size();z++) {

				if (m_oaRecycledVolFrames[z]->m_pVolHeader == header) {

					if (IsPL( BQB_PL_VERBOSE ))
						this->printf("CBQBInterface::CreateVolumetricFrame(): Using frame %lu/%lu for recycling.\n",
							z+1,(unsigned long)m_oaRecycledVolFrames.size());

					vf = m_oaRecycledVolFrames.back();
					m_oaRecycledVolFrames.erase( m_oaRecycledVolFrames.begin() + z );

					if (IsPL( BQB_PL_VERBOSE ))
						this->printf("CBQBInterface::CreateVolumetricFrame(): Done.\n");

					return vf;
				}
			}
		}
	}

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateVolumetricFrame(): Creating new frame ...\n");

	vf = new CBQBVolFrame( *this );
	vf->m_pVolHeader = header;
	vf->m_oCubeFrame.m_iRes[0] = header->GetResolutionX();
	vf->m_oCubeFrame.m_iRes[1] = header->GetResolutionY();
	vf->m_oCubeFrame.m_iRes[2] = header->GetResolutionZ();
	vf->m_oCubeFrame.m_iResXY = vf->m_oCubeFrame.m_iRes[0] * vf->m_oCubeFrame.m_iRes[1];
	vf->m_oCubeFrame.m_iResXYZ = vf->m_oCubeFrame.m_iRes[0] * vf->m_oCubeFrame.m_iRes[1] * vf->m_oCubeFrame.m_iRes[2];
	vf->m_oCubeFrame.m_faBin.resize( header->GetResolutionXYZ() );
	vf->m_oCubeFrame.m_iaExpo.resize( header->GetResolutionXYZ() );
	vf->m_oCubeFrame.m_iaMantis.resize( header->GetResolutionXYZ() );
	vf->m_oCubeFrame.m_iSigni = header->GetSigni();
	vf->m_oCubeFrame.m_iEps = header->GetEps();
	vf->m_oCubeFrame.m_iID = id;

	for (z=0;z<3;z++) {
		vf->m_oCubeFrame.m_fCenter[z] = header->m_faCenter[z];
		vf->m_oCubeFrame.m_fStrideA[z] = header->m_faStrideA[z];
		vf->m_oCubeFrame.m_fStrideB[z] = header->m_faStrideB[z];
		vf->m_oCubeFrame.m_fStrideC[z] = header->m_faStrideC[z];
		vf->m_oCubeFrame.m_iCenter[z] = header->m_iaCenter[z];
		vf->m_oCubeFrame.m_iStrideA[z] = header->m_iaStrideA[z];
		vf->m_oCubeFrame.m_iStrideB[z] = header->m_iaStrideB[z];
		vf->m_oCubeFrame.m_iStrideC[z] = header->m_iaStrideC[z];
	}

	m_oaListVolFrame.insert( vf );

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateVolumetricFrame(): Done.\n");

	return vf;
}



CBQBAtomPosFrame* CBQBInterface::CreateAtomPosFrame( const CBQBAtomHeader *header ) {

	CBQBAtomPosFrame *af;
	unsigned long z;


	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateAtomPosFrame(): Creating atom position frame ...\n");

	if (&header->m_IF != this) {
		this->eprintf("CBQBInterface::CreateAtomPosFrame(): Error: Supplied CBQBAtomHeader does not belong to this CBQBInterface.\n");
		this->eprintf("\n");
		abort();
	}

	if (!m_oaRecycledAtomPosFrames.empty()) {

		if (IsPL( BQB_PL_VERBOSE ))
			this->printf("CBQBInterface::CreateAtomPosFrame(): Recycling frame ...\n");

		if (m_oaRecycledAtomPosFrames.back()->m_pAtomHeader == header) {

			if (IsPL( BQB_PL_VERBOSE ))
				this->printf("CBQBInterface::CreateAtomPosFrame(): Using last frame (%lu/%lu) for recycling.\n",
					(unsigned long)m_oaRecycledAtomPosFrames.size(), (unsigned long)m_oaRecycledAtomPosFrames.size() );

			af = m_oaRecycledAtomPosFrames.back();
			m_oaRecycledAtomPosFrames.pop_back();

			if (IsPL( BQB_PL_VERBOSE ))
				this->printf("CBQBInterface::CreateAtomPosFrame(): Done.\n");

			return af;

		} else {

			for (z=0;z<m_oaRecycledAtomPosFrames.size();z++) {

				if (m_oaRecycledAtomPosFrames[z]->m_pAtomHeader == header) {

					if (IsPL( BQB_PL_VERBOSE ))
						this->printf("CBQBInterface::CreateAtomPosFrame(): Using frame %lu/%lu for recycling.\n",
							z+1,(unsigned long)m_oaRecycledAtomPosFrames.size());

					af = m_oaRecycledAtomPosFrames.back();
					m_oaRecycledAtomPosFrames.erase( m_oaRecycledAtomPosFrames.begin() + z );

					if (IsPL( BQB_PL_VERBOSE ))
						this->printf("CBQBInterface::CreateAtomPosFrame(): Done.\n");

					return af;
				}
			}
		}
	}

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateVolumetricFrame(): Creating new frame ...\n");

	af = new CBQBAtomPosFrame( *this );

	af->m_pAtomHeader = header;
	af->m_oAtomSet.m_oaAtoms.resize( header->m_iAtomCount );

	for (z=0;z<header->m_iAtomCount;z++) {
		af->m_oAtomSet.m_oaAtoms[ z ] = new CBQBAtom();
		af->m_oAtomSet.m_oaAtoms[ z ]->m_iOrd = header->m_iaOrd[z];
		af->m_oAtomSet.m_oaAtoms[ z ]->m_sLabel = header->m_saLabel[z];
	}
	af->m_oAtomSet.m_iSigni = header->GetPrecision();

	af->m_oAtomSet.m_bOrd = true;
	af->m_oAtomSet.m_bLabels = true;

	m_oaListAtomPosFrame.insert( af );

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::CreateAtomPosFrame(): Done.\n");

	return af;
}



bool CBQBInterface::WriteCubeFrame(
		FILE *outfile,
		const CBQBVolHeader *volhead,
		const CBQBVolFrame *volframe,
		const CBQBAtomHeader *atomhead,
		const CBQBAtomPosFrame *atomframe,
		const CBQBCellInfo *cellinfo,
		unsigned long flags
	) {

	UNUSED(outfile);
	UNUSED(volhead);
	UNUSED(volframe);
	UNUSED(atomhead);
	UNUSED(atomframe);
	UNUSED(cellinfo);
	UNUSED(flags);

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::WriteCubeFrame(): Writing CUBE frame ...\n");

	SetLastError( BQB_OK );

	if (IsPL( BQB_PL_VERBOSE ))
		this->printf("CBQBInterface::WriteCubeFrame(): Done.\n");

	return true;
}



/*******************************************************************************************
******      CBQBFileWriter      ************************************************************
*******************************************************************************************/



CBQBFileWriter::~CBQBFileWriter() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::~CBQBFileWriter(): Destructing instance.\n");

	if (m_bActive) {
		m_IF.eprintf("CBQBFileWriter::~CBQBFileWriter(): Error: Cannot manually destruct instance of this class.\n");
		m_IF.eprintf("This is done automatically when the CBQBInterface is released via BQBReleaseInterface().\n");
		m_IF.eprintf("\n");
		abort();
	}

	if (m_pParmPosVol != NULL) {
		delete m_pParmPosVol;
		m_pParmPosVol = NULL;
	}

	if (m_pEngine != NULL) {
		delete m_pEngine;
		m_pEngine = NULL;
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::~CBQBFileWriter(): Done.\n");
}



bool CBQBFileWriter::OpenFileWrite( const char *s, int history, unsigned long contents, unsigned long flags ) {

	UNUSED(flags);

	int z;


	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::OpenFileWrite: Opening BQB file for writing ...\n");

	if ((contents & BQB_CONTENTS_ATOMPOS) != 0)
		m_bContainsAtomPosition = true;

	if ((contents & BQB_CONTENTS_VOLUMETRIC) != 0)
		m_bContainsVolumetricData = true;

	if ((contents & BQB_CONTENTS_CELLINFO) != 0)
		m_bContainsCellInfo = true;

	m_pEngine = new CBQBEngine( m_IF );

	//m_pEngine->DisableDeleteOutFrames();

	m_iHistoryDepth = history;

/*	if (m_iHistoryDepth == 0) {
		if (m_IF.IsPL( BQB_PL_VERBOSE ))
			m_IF.printf("CBQBFileWriter::OpenFileWrite: History depth not set. Using 10 as default value.\n");
		m_iHistoryDepth = 10;
	}*/

	m_pParmPosVol = new CBQBParameterSet_PosAndVol( m_IF );

	// Set default values
	m_pParmPosVol->SetVolSigni(5);
	m_pParmPosVol->SetVolEps(12);
	m_pParmPosVol->SetVolOrder( MAX( MIN( m_iHistoryDepth-1, 8 ), 0 ) );
	m_pParmPosVol->SetVolOptOrder(true);
	m_pParmPosVol->SetVolHilbert(true);
	m_pParmPosVol->SetVolNbhFac(1.075);
	m_pParmPosVol->SetVolSplit(10);
	m_pParmPosVol->SetVolTableCount(6);
	m_pParmPosVol->SetVolOptTables(false);
	m_pParmPosVol->SetVolBlockLength(20);
	m_pParmPosVol->SetVolBW(false);
	m_pParmPosVol->SetVolMTF(false);
	m_pParmPosVol->SetVolMaxIter(10);
	m_pParmPosVol->SetVolRLE(true);
	m_pParmPosVol->SetVolMaxChunk(0);
	m_pParmPosVol->SetPosPrecision(6);
	m_pParmPosVol->SetPosOrder( MAX( MIN( m_iHistoryDepth-1, 8 ), 0 ) );
	m_pParmPosVol->SetPosOptOrder(true);
	m_pParmPosVol->SetPosSplit(14);
	m_pParmPosVol->SetPosTableCount(1);
	m_pParmPosVol->SetPosOptTables(false);
	m_pParmPosVol->SetPosBlockLength(40);
	m_pParmPosVol->SetPosBW(false);
	m_pParmPosVol->SetPosMTF(false);
	m_pParmPosVol->SetPosRLE(true);
	m_pParmPosVol->SetPosMaxIter(10);
	m_pParmPosVol->SetPosSortAtom(true);
	m_pParmPosVol->SetPosMaxChunk(0);

	m_pParmPosVol->SetVolUseExtra(true);
	m_pParmPosVol->SetVolExtraSRange(7);
	m_pParmPosVol->SetVolExtraTRange( MAX( MIN( m_iHistoryDepth-1, 6 ), 1 ) );
	m_pParmPosVol->SetVolExtraSOrder(3);
	m_pParmPosVol->SetVolExtraTOrder( MAX( MIN( m_iHistoryDepth-2, 2 ), 0 ) );
	m_pParmPosVol->SetVolExtraOffset(3);
	m_pParmPosVol->SetVolExtraCrossS(true);
	m_pParmPosVol->SetVolExtraCrossT(false);
	m_pParmPosVol->SetVolExtraWrap(true);
	m_pParmPosVol->SetVolExtraCrossRangeS(true);
	m_pParmPosVol->SetVolExtraCrossRangeT(false);
	m_pParmPosVol->SetVolExtraDistExpo(3.0);
	m_pParmPosVol->SetVolExtraTimeExpo(1.0);
	m_pParmPosVol->SetVolExtraPredCorr(true);

	m_pParmPosVol->SetPosUseExtra(true);
	m_pParmPosVol->SetPosExtraTRange( MAX( MIN( m_iHistoryDepth-1, 9 ), 1 ) );
	m_pParmPosVol->SetPosExtraTOrder( MAX( MIN( m_iHistoryDepth-2, 6 ), 0 ) );
	m_pParmPosVol->SetPosExtraTimeExpo(4.0);


	if (m_pParmPosVol->GetVolUseExtra())
		if (m_pParmPosVol->GetVolExtraTRange() > m_pParmPosVol->GetVolOrder())
			m_pParmPosVol->SetVolOrder( m_pParmPosVol->GetVolExtraTRange() );

	if (m_pParmPosVol->GetPosUseExtra())
		if (m_pParmPosVol->GetPosExtraTRange() > m_pParmPosVol->GetPosOrder())
			m_pParmPosVol->SetPosOrder( m_pParmPosVol->GetPosExtraTRange() );

	if (m_bContainsAtomPosition) {

		m_pEngine->m_pReadCacheXYZ = new CBQBReadXYZCache( m_IF );

		m_pEngine->m_pReadCacheXYZ->ActivateInterfaceMode();

		m_pEngine->m_pReadCacheXYZ->SetReadParameters( m_pParmPosVol->GetPosPrecision() );
	}

	if (m_bContainsVolumetricData) {

		m_pEngine->m_pReadCacheCube = new CBQBReadCubeCache( m_IF );

		m_pEngine->m_pReadCacheCube->ActivateInterfaceMode();

		m_pEngine->m_pReadCacheCube->SetReadParameters(
			m_pParmPosVol->GetVolEps(),
			m_pParmPosVol->GetVolSigni(),
			m_pParmPosVol->GetPosPrecision()
		);
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::OpenFileWrite: The history depth is %d.\n", m_iHistoryDepth );

	if (m_bContainsAtomPosition) {
		m_pEngine->m_pReadCacheXYZ->SetHistoryDepth( m_iHistoryDepth );
		m_pEngine->m_pReadCacheXYZ->Reset();
		m_pEngine->m_pReadCacheXYZ->SetStartLock();
	}

	if (m_bContainsVolumetricData) {
		m_pEngine->m_pReadCacheCube->SetHistoryDepth( m_iHistoryDepth );
		m_pEngine->m_pReadCacheCube->Reset();
		m_pEngine->m_pReadCacheCube->SetStartLock();
	}

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Opening compressed file \"%s\" ...\n",s);
	m_pBQBFile = new CBQBFile( m_IF );

	if (flags == BQB_OVERWRITE) {
		if (!m_pBQBFile->OpenWriteReplace( s )) {
			m_IF.eprintf("Error: Could not open file for writing.\n");
			return false;
		}
	} else { // APPEND is the default
		if (!m_pBQBFile->OpenWriteAppend( s )) {
			m_IF.eprintf("Error: Could not open file for writing.\n");
			return false;
		}
	}

	m_iReadPos = -1;

	if (m_bContainsVolumetricData) {
		m_oaVolumetricFrameHistory.resize( m_iHistoryDepth );
		m_iVolumetricWritePos = -1;
	}

	if (m_bContainsAtomPosition) {
		m_oaAtomPosFrameHistory.resize( m_iHistoryDepth );
		m_iAtomPosWritePos = -1;
	}

	m_iCommentLineWritePos = -1;

	m_saCommentLine1.resize( m_iHistoryDepth );
	m_saCommentLine2.resize( m_iHistoryDepth );

	m_bStartLock = true;

	if (m_bContainsVolumetricData) {
		m_pEngine->m_oaOutputCubeBuf.resize( m_pParmPosVol->GetVolOrder()+1 );
		for (z=0;z<m_pParmPosVol->GetVolOrder()+1;z++)
			m_pEngine->m_oaOutputCubeBuf[z] = NULL;
		m_pEngine->m_iOutputCubeBufPos = 0;
	}

	if (m_bContainsAtomPosition) {
		m_pEngine->m_oaOutputAtomBuf.resize( m_pParmPosVol->GetPosOrder()+1 );
		for (z=0;z<m_pParmPosVol->GetPosOrder()+1;z++)
			m_pEngine->m_oaOutputAtomBuf[z] = NULL;
		m_pEngine->m_iOutputAtomBufPos = 0;
	}

	m_IF.SetLastError( BQB_OK );

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::OpenFileWrite: Done.\n");

	return true;
}



void CBQBFileWriter::CloseFile() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::CloseFile(): Closing BQB file ...\n");

	if (!m_bActive) {
		m_IF.eprintf("CBQBFileWriter::CloseFile(): Error: File was already closed before.\n");
		m_IF.eprintf("\n");
		abort();
	}

	if (!m_pBQBFile->Close()) {
		m_IF.eprintf("CBQBFileWriter::CloseFile(): Error: Could not close BQB file.\n");
		m_IF.eprintf("\n");
		abort();
	}

	delete m_pBQBFile;
	m_pBQBFile = NULL;

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::CloseFile(): Destructing instances ...\n");

	if (m_pEngine != NULL) {
		delete m_pEngine;
		m_pEngine = NULL;
	}

	if (m_pParmPosVol != NULL) {
		delete m_pParmPosVol;
		m_pParmPosVol = NULL;
	}

	m_bActive = false;

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::CloseFile(): Done.\n");
}



void CBQBFileWriter::SetVolumetricHeader( CBQBVolHeader *header ) {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::SetVolumetricHeader(): Setting volumetric header ...\n");

	if (m_bAlreadyPushed) {
		m_IF.eprintf("CBQBFileWriter::SetVolumetricHeader(): Error: Cannot set volumetric header after call to PushFrame() functions.\n");
		m_IF.eprintf("\n");
		abort();
	}

	if (!m_bContainsVolumetricData) {
		m_IF.eprintf("CBQBFileWriter::SetVolumetricHeader(): Error: File does not contain volumetric data.\n");
		m_IF.eprintf("\n");
		abort();
	}

	if (&header->m_IF != &m_IF) {
		m_IF.eprintf("CBQBFileWriter::SetVolumetricHeader(): Error: Supplied CBQBVolHeader does not belong to this CBQBInterface.\n");
		m_IF.eprintf("\n");
		abort();
	}

	header->m_bAssigned = true;

	m_pVolHeader = header;

	if (m_pCellInfo != NULL) {
		m_pVolHeader->m_faStrideA[0] = m_pCellInfo->GetCellAX() / m_pVolHeader->GetResolutionX();
		m_pVolHeader->m_faStrideA[1] = m_pCellInfo->GetCellAY() / m_pVolHeader->GetResolutionX();
		m_pVolHeader->m_faStrideA[2] = m_pCellInfo->GetCellAZ() / m_pVolHeader->GetResolutionX();
		m_pVolHeader->m_faStrideB[0] = m_pCellInfo->GetCellBX() / m_pVolHeader->GetResolutionY();
		m_pVolHeader->m_faStrideB[1] = m_pCellInfo->GetCellBY() / m_pVolHeader->GetResolutionY();
		m_pVolHeader->m_faStrideB[2] = m_pCellInfo->GetCellBZ() / m_pVolHeader->GetResolutionY();
		m_pVolHeader->m_faStrideC[0] = m_pCellInfo->GetCellCX() / m_pVolHeader->GetResolutionZ();
		m_pVolHeader->m_faStrideC[1] = m_pCellInfo->GetCellCY() / m_pVolHeader->GetResolutionZ();
		m_pVolHeader->m_faStrideC[2] = m_pCellInfo->GetCellCZ() / m_pVolHeader->GetResolutionZ();
		m_pVolHeader->m_iaStrideA[0] = FloatToFixed( m_pVolHeader->m_faStrideA[0], 6 );
		m_pVolHeader->m_iaStrideA[1] = FloatToFixed( m_pVolHeader->m_faStrideA[1], 6 );
		m_pVolHeader->m_iaStrideA[2] = FloatToFixed( m_pVolHeader->m_faStrideA[2], 6 );
		m_pVolHeader->m_iaStrideB[0] = FloatToFixed( m_pVolHeader->m_faStrideB[0], 6 );
		m_pVolHeader->m_iaStrideB[1] = FloatToFixed( m_pVolHeader->m_faStrideB[1], 6 );
		m_pVolHeader->m_iaStrideB[2] = FloatToFixed( m_pVolHeader->m_faStrideB[2], 6 );
		m_pVolHeader->m_iaStrideC[0] = FloatToFixed( m_pVolHeader->m_faStrideC[0], 6 );
		m_pVolHeader->m_iaStrideC[1] = FloatToFixed( m_pVolHeader->m_faStrideC[1], 6 );
		m_pVolHeader->m_iaStrideC[2] = FloatToFixed( m_pVolHeader->m_faStrideC[2], 6 );
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::SetVolumetricHeader(): Done.\n");
}



void CBQBFileWriter::SetAtomHeader( CBQBAtomHeader *header ) {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::SetAtomHeader(): Setting atom header ...\n");

	if (m_bAlreadyPushed) {
		m_IF.eprintf("CBQBFileWriter::SetAtomHeader(): Error: Cannot set atom header after call to PushFrame() functions.\n");
		m_IF.eprintf("\n");
		abort();
	}

	if (!m_bContainsAtomPosition) {
		m_IF.eprintf("CBQBFileWriter::SetAtomHeader(): Error: File does not contain atom positions.\n");
		m_IF.eprintf("\n");
		abort();
	}

	if (&header->m_IF != &m_IF) {
		m_IF.eprintf("CBQBFileWriter::SetAtomHeader(): Error: Supplied CBQBAtomHeader does not belong to this CBQBInterface.\n");
		m_IF.eprintf("\n");
		abort();
	}

	header->m_bAssigned = true;

	m_pAtomHeader = header;

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::SetAtomHeader(): Done.\n");
}



void CBQBFileWriter::SetCellInfo( CBQBCellInfo *info ) {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::SetCellInfo(): Setting cell info ...\n");

	if (m_bAlreadyPushed) {
		m_IF.eprintf("CBQBFileWriter::SetCellInfo(): Error: Cannot set cell info after call to PushFrame() functions.\n");
		m_IF.eprintf("\n");
		abort();
	}

	if (!m_bContainsCellInfo) {
		m_IF.eprintf("CBQBFileWriter::SetCellInfo(): Error: File does not contain cell info.\n");
		m_IF.eprintf("\n");
		abort();
	}

	if (&info->m_IF != &m_IF) {
		m_IF.eprintf("CBQBFileWriter::SetCellInfo(): Error: Supplied CBQBCellInfo does not belong to this CBQBInterface.\n");
		m_IF.eprintf("\n");
		abort();
	}

	info->m_bAssigned = true;

	m_pCellInfo = info;

	if (m_pVolHeader != NULL) {
		m_pVolHeader->m_faStrideA[0] = m_pCellInfo->GetCellAX() / m_pVolHeader->GetResolutionX();
		m_pVolHeader->m_faStrideA[1] = m_pCellInfo->GetCellAY() / m_pVolHeader->GetResolutionX();
		m_pVolHeader->m_faStrideA[2] = m_pCellInfo->GetCellAZ() / m_pVolHeader->GetResolutionX();
		m_pVolHeader->m_faStrideB[0] = m_pCellInfo->GetCellBX() / m_pVolHeader->GetResolutionY();
		m_pVolHeader->m_faStrideB[1] = m_pCellInfo->GetCellBY() / m_pVolHeader->GetResolutionY();
		m_pVolHeader->m_faStrideB[2] = m_pCellInfo->GetCellBZ() / m_pVolHeader->GetResolutionY();
		m_pVolHeader->m_faStrideC[0] = m_pCellInfo->GetCellCX() / m_pVolHeader->GetResolutionZ();
		m_pVolHeader->m_faStrideC[1] = m_pCellInfo->GetCellCY() / m_pVolHeader->GetResolutionZ();
		m_pVolHeader->m_faStrideC[2] = m_pCellInfo->GetCellCZ() / m_pVolHeader->GetResolutionZ();
		m_pVolHeader->m_iaStrideA[0] = FloatToFixed( m_pVolHeader->m_faStrideA[0], 6 );
		m_pVolHeader->m_iaStrideA[1] = FloatToFixed( m_pVolHeader->m_faStrideA[1], 6 );
		m_pVolHeader->m_iaStrideA[2] = FloatToFixed( m_pVolHeader->m_faStrideA[2], 6 );
		m_pVolHeader->m_iaStrideB[0] = FloatToFixed( m_pVolHeader->m_faStrideB[0], 6 );
		m_pVolHeader->m_iaStrideB[1] = FloatToFixed( m_pVolHeader->m_faStrideB[1], 6 );
		m_pVolHeader->m_iaStrideB[2] = FloatToFixed( m_pVolHeader->m_faStrideB[2], 6 );
		m_pVolHeader->m_iaStrideC[0] = FloatToFixed( m_pVolHeader->m_faStrideC[0], 6 );
		m_pVolHeader->m_iaStrideC[1] = FloatToFixed( m_pVolHeader->m_faStrideC[1], 6 );
		m_pVolHeader->m_iaStrideC[2] = FloatToFixed( m_pVolHeader->m_faStrideC[2], 6 );
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::SetCellInfo(): Done.\n");
}



const CBQBVolHeader* CBQBFileWriter::GetVolumetricHeader() const {
	if (!m_bContainsVolumetricData) {
		m_IF.eprintf("CBQBFileWriter::GetVolumetricHeader(): Error: File does not contain volumetric data.\n");
		m_IF.eprintf("\n");
		abort();
	}
	return m_pVolHeader;
}



const CBQBAtomHeader* CBQBFileWriter::GetAtomHeader() const {
	if (!m_bContainsAtomPosition) {
		m_IF.eprintf("CBQBFileWriter::GetAtomHeader(): Error: File does not contain atom position data.\n");
		m_IF.eprintf("\n");
		abort();
	}
	return m_pAtomHeader;
}



const CBQBCellInfo* CBQBFileWriter::GetCellInfo() const {
	if (!m_bContainsCellInfo) {
		m_IF.eprintf("CBQBFileWriter::GetCellInfo(): Error: File does not contain cell info.\n");
		m_IF.eprintf("\n");
		abort();
	}
	return m_pCellInfo;
}



void CBQBFileWriter::PushVolumetricFrame( CBQBVolFrame *frame ) {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::PushVolumetricFrame(): Pushing volumetric frame ...\n");

	m_bAlreadyPushed = true;

	m_iVolumetricWritePos++;
	if (m_iVolumetricWritePos >= (int)m_oaVolumetricFrameHistory.size()) {
		if (m_bStartLock) {
			m_IF.eprintf("CBQBFileWriter::PushVolumetricFrame(): Error: Violation of start lock.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iVolumetricWritePos = 0;
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE )) {
		m_IF.printf("CBQBFileWriter::PushVolumetricFrame(): Current write position is now %d (history depth %d).\n",m_iVolumetricWritePos,m_iHistoryDepth);
		m_IF.printf("CBQBFileWriter::PushVolumetricFrame(): There are now %d pending volumetric frames.\n",GetPendingVolumetricFrameCount());
	}

	m_iVolumetricFramesPending++;

	if (m_iVolumetricFramesPending > m_iHistoryDepth) {
		m_IF.eprintf("CBQBFileWriter::PushVolumetricFrame(): Error: Tried to push more frames than history depth (%d).\n",m_iHistoryDepth);
		m_IF.eprintf("\n");
		abort();
	}

	frame->m_oCubeFrame.FloatToInt();

	if (m_oaVolumetricFrameHistory[ m_iVolumetricWritePos ] != NULL)
		m_IF.RecycleVolumetricFrame( m_oaVolumetricFrameHistory[ m_iVolumetricWritePos ] );

	m_oaVolumetricFrameHistory[ m_iVolumetricWritePos ] = frame;

	if (m_bContainsAtomPosition)
		if (m_oaAtomPosFrameHistory[ m_iVolumetricWritePos ] != NULL)
			m_oaVolumetricFrameHistory[ m_iVolumetricWritePos ]->m_oCubeFrame.m_pAtoms = &m_oaAtomPosFrameHistory[ m_iVolumetricWritePos ]->m_oAtomSet;

	m_pEngine->m_pReadCacheCube->PushFrame( &frame->m_oCubeFrame );

	m_IF.SetLastError( BQB_OK );

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::PushVolumetricFrame(): Done.\n");
}



void CBQBFileWriter::PushAtomPosFrame( CBQBAtomPosFrame *frame ) {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::PushAtomPosFrame(): Pushing atom position frame ...\n");

	m_bAlreadyPushed = true;

	m_iAtomPosWritePos++;
	if (m_iAtomPosWritePos >= (int)m_oaAtomPosFrameHistory.size()) {
		if (m_bStartLock) {
			m_IF.eprintf("CBQBFileWriter::PushAtomPosFrame(): Error: Violation of start lock.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iAtomPosWritePos = 0;
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE )) {
		m_IF.printf("CBQBFileWriter::PushAtomPosFrame(): Current write position is now %d (history depth %d).\n",m_iAtomPosWritePos,m_iHistoryDepth);
		m_IF.printf("CBQBFileWriter::PushAtomPosFrame(): There are now %d pending atom position frames.\n",GetPendingAtomPosFrameCount());
	}

	m_iAtomPosFramesPending++;

	if (m_iAtomPosFramesPending > m_iHistoryDepth) {
		m_IF.eprintf("CBQBFileWriter::PushAtomPosFrame(): Error: Tried to push more frames than history depth (%d).\n",m_iHistoryDepth);
		m_IF.eprintf("\n");
		abort();
	}

	frame->m_oAtomSet.FloatToInt();

	if (m_oaAtomPosFrameHistory[ m_iAtomPosWritePos ] != NULL)
		m_IF.RecycleAtomPositionFrame( m_oaAtomPosFrameHistory[ m_iAtomPosWritePos ] );

	m_oaAtomPosFrameHistory[ m_iAtomPosWritePos ] = frame;

	if (m_bContainsVolumetricData)
		if (m_oaVolumetricFrameHistory[ m_iAtomPosWritePos ] != NULL)
			m_oaVolumetricFrameHistory[ m_iAtomPosWritePos ]->m_oCubeFrame.m_pAtoms = &frame->m_oAtomSet;

	m_pEngine->m_pReadCacheXYZ->PushFrame( &frame->m_oAtomSet );

	m_IF.SetLastError( BQB_OK );

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::PushAtomPosFrame(): Done.\n");
}



void CBQBFileWriter::PushCommentLine( const char *line1, const char *line2 ) {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::PushCommentLine(): Pushing comment lines ...\n");

	m_bAlreadyPushed = true;

	m_iCommentLineWritePos++;
	if (m_iCommentLineWritePos >= (int)m_saCommentLine1.size()) {
		if (m_bStartLock) {
			m_IF.eprintf("CBQBFileWriter::PushCommentLine(): Error: Violation of start lock.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iCommentLineWritePos = 0;
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE )) {
		m_IF.printf("CBQBFileWriter::PushCommentLine(): Current write position is now %d (history depth %d).\n",m_iCommentLineWritePos,m_iHistoryDepth);
		m_IF.printf("CBQBFileWriter::PushCommentLine(): There are now %d pending comment lines.\n",GetPendingCommentLineCount());
	}

	m_iCommentLinesPending++;

	if (m_iCommentLinesPending > m_iHistoryDepth) {
		m_IF.eprintf("CBQBFileWriter::PushCommentLine(): Error: Tried to push more comment lines than history depth (%d).\n",m_iHistoryDepth);
		m_IF.eprintf("\n");
		abort();
	}

	if (line1 != NULL)
		m_saCommentLine1[ m_iCommentLineWritePos ] = line1;
	else
		m_saCommentLine1[ m_iCommentLineWritePos ] = "(Comment Line 1)";

	if (line2 != NULL)
		m_saCommentLine2[ m_iCommentLineWritePos ] = line2;
	else
		m_saCommentLine2[ m_iCommentLineWritePos ] = "(Comment Line 2)";

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::PushCommentLine(): Done.\n");
}



void CBQBFileWriter::Optimize( int mode ) {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::Optimize(): Optimizing parameters (mode=%d) ...\n",mode);

	if (!m_bStartLock) {
		m_IF.eprintf("CBQBFileWriter::Optimize(): Error: Can only call optimized before first call to WritePendingFrames().\n");
		m_IF.eprintf("\n");
		abort();
	}

	CBQBDriver driver( m_IF );

	// Bugfix MB 24.12.2020: Constructor already creates an engine which we do not need...
	if (driver.m_pEngine != NULL)
		delete driver.m_pEngine;

	driver.m_pEngine = m_pEngine;

	if (GetPendingVolumetricFrameCount() > 10)
		driver.m_iOptIncludeFirst = 0;
	else
		driver.m_iOptIncludeFirst = 1;

	if (m_bContainsAtomPosition) {

		driver.OptimizeXYZParameters( mode, m_pParmPosVol->GetPositionParameterSet(), true, true );

		m_pEngine->m_pReadCacheXYZ->RewindReadPos();

		if (m_pParmPosVol->GetPosUseExtra())
			if (m_pParmPosVol->GetPosExtraTRange() > m_pParmPosVol->GetPosOrder())
				m_pParmPosVol->SetPosOrder(m_pParmPosVol->GetPosExtraTRange());
	}

	if (m_bContainsVolumetricData) {

		driver.OptimizeCubeParameters( mode, m_pParmPosVol->GetVolumetricParameterSet() );

		m_pEngine->m_pReadCacheCube->RewindReadPos();

		if (m_pParmPosVol->GetVolUseExtra())
			if (m_pParmPosVol->GetVolExtraTRange() > m_pParmPosVol->GetVolOrder())
				m_pParmPosVol->SetVolOrder(m_pParmPosVol->GetVolExtraTRange());
	}

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		
		m_IF.printf("\n");
		m_IF.bprintf("    The optimal parameters for this trajectory are:\n");
		m_IF.printf("\n");

		m_IF.printf("%s",m_pParmPosVol->ToString(4).c_str());
		m_IF.printf("\n");

		m_IF.bprintf("    Parameter key:");
		m_IF.printf("  %s\n",m_pParmPosVol->ToKey().c_str());
		m_IF.printf("\n");
	}

	// To avoid deleting it...
	driver.m_pEngine = NULL;

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::Optimize(): Done.\n");
}



bool CBQBFileWriter::WritePendingFrames( bool check ) {

	CBQBBitSet bsat(m_IF), bscu(m_IF), bshe(m_IF);
	bool err, cl;
	int z, z2, histused, ao, co, ft, thu, id;
	const CBQBCubeFrame *cfr = NULL;
	CBQBCubeFrame *cfr2;
	char buf[512];


	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::WritePendingFrames(): Writing pending frames ...\n");

	m_bStartLock = false;

	if (m_iFramesWritten == 0) {

		if (m_IF.IsPL( BQB_PL_VERBOSE ))
			m_IF.printf("CBQBFileWriter::WritePendingFrames(): This is the first write call.\n");

		if (m_bContainsAtomPosition)
			m_pEngine->m_pReadCacheXYZ->LiftStartLock();

		if (m_bContainsVolumetricData)
			m_pEngine->m_pReadCacheCube->LiftStartLock();

		if (m_pParmPosVol->GetVolUseExtra()) {

			if (m_pParmPosVol->GetVolExtraPredCorr()) {

				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("      Initializing Volumetric Predictor Extrapolator (%d/%d)...\n",
						m_pParmPosVol->GetVolExtraTRange(),m_pParmPosVol->GetVolExtraTOrder());

				m_pEngine->m_pExtrapolator = new CBQBExtrapolator(m_IF);
				m_pEngine->m_pExtrapolator->Initialize(
					m_pVolHeader->GetResolutionX(), // resx,
					m_pVolHeader->GetResolutionY(), // resy,
					m_pVolHeader->GetResolutionZ(), // resz,
					1,              // srangex,
					1,              // srangey,
					1,              // srangez,
					m_pParmPosVol->GetVolExtraTRange(),   // trange,
					0,              // sorder,
					m_pParmPosVol->GetVolExtraTOrder(),   // torder,
					0,              // offsetx,
					0,              // offsety,
					0,              // offsetz,
					false,          // crosss,
					false,          // crosst,
					false,          // wrap,
					false,          // crossranges,
					false,          // crossranget
					0,              // distexpo
					m_pParmPosVol->GetVolExtraTimeExpo(), // timeexpo
					false
				);

				if (m_IF.IsPL(BQB_PL_STANDARD)) {
					if (m_pParmPosVol->IsVolExtraSRangeXYZEqual())
						m_IF.printf("      Initializing Volumetric Corrector Extrapolator (%d/%d)...\n",
							m_pParmPosVol->GetVolExtraSRangeX(),m_pParmPosVol->GetVolExtraSOrder());
					else
						m_IF.printf("      Initializing Volumetric Corrector Extrapolator (%d|%d|%d/%d)...\n",
							m_pParmPosVol->GetVolExtraSRangeX(),m_pParmPosVol->GetVolExtraSRangeY(),m_pParmPosVol->GetVolExtraSRangeZ(),m_pParmPosVol->GetVolExtraSOrder());
				}

				m_pEngine->m_pExtrapolatorCorr = new CBQBExtrapolator(m_IF);
				m_pEngine->m_pExtrapolatorCorr->Initialize(
					m_pVolHeader->GetResolutionX(), // resx,
					m_pVolHeader->GetResolutionY(), // resy,
					m_pVolHeader->GetResolutionZ(), // resz,
					m_pParmPosVol->GetVolExtraSRangeX(),      // srangex,
					m_pParmPosVol->GetVolExtraSRangeY(),      // srangey,
					m_pParmPosVol->GetVolExtraSRangeZ(),      // srangez,
					1,                  // trange,
					m_pParmPosVol->GetVolExtraSOrder(),       // sorder,
					0,                  // torder,
					m_pParmPosVol->GetVolExtraOffsetX(),      // offsetx,
					m_pParmPosVol->GetVolExtraOffsetY(),      // offsety,
					m_pParmPosVol->GetVolExtraOffsetZ(),      // offsetz,
					m_pParmPosVol->GetVolExtraCrossS(),       // crosss,
					false,              // crosst,
					m_pParmPosVol->GetVolExtraWrap(),         // wrap,
					m_pParmPosVol->GetVolExtraCrossRangeS(),  // crossranges,
					false,              // crossranget
					m_pParmPosVol->GetVolExtraDistExpo(),     // distexpo
					0,                  // timeexpo
					false
				);

			} else {

				if (m_IF.IsPL(BQB_PL_STANDARD)) {

					if (m_pParmPosVol->IsVolExtraSRangeXYZEqual())
						m_IF.printf("      Initializing Volumetric Extrapolator (%d/%d;%d/%d)...\n",
							m_pParmPosVol->GetVolExtraSRangeX(),
							m_pParmPosVol->GetVolExtraSOrder(),
							m_pParmPosVol->GetVolExtraTRange(),
							m_pParmPosVol->GetVolExtraTOrder()
						);
					else
						m_IF.printf("      Initializing Volumetric Extrapolator (%d|%d|%d/%d;%d|%d)...\n",
							m_pParmPosVol->GetVolExtraSRangeX(),
							m_pParmPosVol->GetVolExtraSRangeY(),
							m_pParmPosVol->GetVolExtraSRangeZ(),
							m_pParmPosVol->GetVolExtraSOrder(),
							m_pParmPosVol->GetVolExtraTRange(),
							m_pParmPosVol->GetVolExtraTOrder()
						);
				}

				m_pEngine->m_pExtrapolator = new CBQBExtrapolator(m_IF);
				m_pEngine->m_pExtrapolator->Initialize(
					m_pVolHeader->GetResolutionX(), // resx,
					m_pVolHeader->GetResolutionY(), // resy,
					m_pVolHeader->GetResolutionZ(), // resz,
					m_pParmPosVol->GetVolExtraSRangeX(),      // srangex,
					m_pParmPosVol->GetVolExtraSRangeY(),      // srangey,
					m_pParmPosVol->GetVolExtraSRangeZ(),      // srangez,
					m_pParmPosVol->GetVolExtraTRange(),       // trange,
					m_pParmPosVol->GetVolExtraSOrder(),       // sorder,
					m_pParmPosVol->GetVolExtraTOrder(),       // torder,
					m_pParmPosVol->GetVolExtraOffsetX(),      // offsetx,
					m_pParmPosVol->GetVolExtraOffsetY(),      // offsety,
					m_pParmPosVol->GetVolExtraOffsetZ(),      // offsetz,
					m_pParmPosVol->GetVolExtraCrossS(),       // crosss,
					m_pParmPosVol->GetVolExtraCrossT(),       // crosst,
					m_pParmPosVol->GetVolExtraWrap(),         // wrap,
					m_pParmPosVol->GetVolExtraCrossRangeS(),  // crossranges,
					m_pParmPosVol->GetVolExtraCrossRangeT(),  // crossranget
					m_pParmPosVol->GetVolExtraDistExpo(),     // distexpo
					m_pParmPosVol->GetVolExtraTimeExpo(),     // timeexpo
					false
				);
			}
		}

		if (m_pParmPosVol->GetPosUseExtra()) {

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Initializing Position Extrapolator (%d/%d)...\n",
					m_pParmPosVol->GetPosExtraTRange(),m_pParmPosVol->GetPosExtraTOrder());

			m_pEngine->m_pExtrapolatorXYZ = new CBQBExtrapolator(m_IF);
			m_pEngine->m_pExtrapolatorXYZ->InitializeXYZ(
				m_pParmPosVol->GetPosExtraTRange(),
				m_pParmPosVol->GetPosExtraTOrder(),
				m_pParmPosVol->GetPosExtraTimeExpo(),
				false
			);
		}

		if (m_pParmPosVol->GetVolUseExtra())
			m_iVolMaxOrder = m_pParmPosVol->GetVolExtraTRange();
		else
			m_iVolMaxOrder = m_pParmPosVol->GetVolOptOrder();

		if (m_pParmPosVol->GetPosUseExtra())
			m_iAtomPosMaxOrder = m_pParmPosVol->GetPosExtraTRange();
		else
			m_iAtomPosMaxOrder = m_pParmPosVol->GetPosOptOrder();

		m_iUseHistory = 0;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(
			"CBQBFileWriter::WritePendingFrames(): Have %d pending volumetric frames and %d pending atom frames.\n",
			GetPendingVolumetricFrameCount(),
			GetPendingAtomPosFrameCount()
		);

	while (true) {

		if (m_bContainsAtomPosition && !HavePendingAtomPosFrames())
			break;

		if (m_bContainsVolumetricData && !HavePendingVolumetricFrames())
			break;

		cl = HavePendingCommentLines();

		if (m_IF.IsPL( BQB_PL_VERBOSE ))
			m_IF.printf("CBQBFileWriter::WritePendingFrames(): Have %d pending volumetric frames and %d pending atom frames. Writing a frame ...\n",
				GetPendingVolumetricFrameCount(), GetPendingAtomPosFrameCount()
			);

		m_iReadPos++;

		if (m_iReadPos >= m_iHistoryDepth)
			m_iReadPos = 0;

		if (m_IF.IsPL( BQB_PL_VERBOSE ))
			m_IF.printf("CBQBFileWriter::WritePendingFrames(): Read position is now %d/%d.\n",
				m_iReadPos, m_iHistoryDepth
			);

		if (cl)
			sprintf( buf, "%s\n%s", m_saCommentLine1[ m_iReadPos ].c_str(), m_saCommentLine2[ m_iReadPos ].c_str() );
		else
			strcpy( buf, "(Comment Line 1)\n(Comment Line 2)" );

		if (m_oaAtomPosFrameHistory[ m_iReadPos ]->m_oAtomSet.m_sComment != NULL)
			delete[] m_oaAtomPosFrameHistory[ m_iReadPos ]->m_oAtomSet.m_sComment;
		m_oaAtomPosFrameHistory[ m_iReadPos ]->m_oAtomSet.m_sComment = new char[ strlen(buf)+1 ];
		strcpy( m_oaAtomPosFrameHistory[ m_iReadPos ]->m_oAtomSet.m_sComment, buf );

		m_iVolumetricFramesPending--;
		m_iAtomPosFramesPending--;
		m_iCommentLinesPending--;

		if (m_bContainsVolumetricData) {

			if (m_pEngine->m_pExtrapolator != NULL)
				m_pEngine->m_pExtrapolator->PushCubeFrame( &m_oaVolumetricFrameHistory[ m_iReadPos ]->m_oCubeFrame );

			cfr = m_pEngine->m_pReadCacheCube->GetNextFrame();
		}

		if (m_bContainsAtomPosition)
			if (m_pEngine->m_pExtrapolatorXYZ != NULL)
				m_pEngine->m_pExtrapolatorXYZ->PushAtomFrame( &m_oaAtomPosFrameHistory[ m_iReadPos ]->m_oAtomSet );

_again:
		if (m_iFramesWritten < m_iVolMaxOrder)
			co = m_iFramesWritten;
		else if (m_iUseHistory < m_iVolMaxOrder)
			co = m_iUseHistory;
		else
			co = m_iVolMaxOrder;

		if (m_iFramesWritten < m_iAtomPosMaxOrder)
			ao = m_iFramesWritten;
		else if (m_iUseHistory < m_iAtomPosMaxOrder)
			ao = m_iUseHistory;
		else
			ao = m_iAtomPosMaxOrder;

		if (m_IF.IsPL( BQB_PL_VERBOSE ))
			m_IF.printf("CBQBFileWriter::WritePendingFrames(): Using volumetric history of %d and atom position history of %d.\n", co, ao );

		bsat.Clear();

		histused = 0;

		m_pEngine->CompressAtomFrame(
			&bsat,
			ao,
			m_pParmPosVol->GetPositionParameterSet(),
			histused,
			(m_iFramesWritten==0), // Store static trajectory information?
			true, // Comment?
			NULL
		);

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (m_pEngine->m_pExtrapolatorXYZ == NULL)
				m_IF.printf("      Atoms: Order %d,     output size %9.3f KiB.\n",ao,bsat.GetByteLength()/1024.0);
			else
				m_IF.printf("      Atoms: History %2d,  output size %9.3f KiB.\n",histused,bsat.GetByteLength()/1024.0);
		}

		bscu.Clear();

		m_pEngine->CompressCubeFrame(
			&bscu,
			co,
			m_pParmPosVol->GetVolumetricParameterSet(),
			thu,
			(m_iFramesWritten==0), // Store static info?
			NULL
		);

		if (thu > histused)
			histused = thu;

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (m_pEngine->m_pExtrapolator == NULL)
				m_IF.printf("      Cube:  Order %d,     output size %9.3f KiB,  ratio %5.2f : 1.\n",co,bscu.GetByteLength()/1024.0,(cfr->m_iResXYZ*13.16666)/bscu.GetByteLength());
			else
				m_IF.printf("      Cube:  History %2d,  output size %9.3f KiB,  ratio %5.2f : 1.\n",thu,bscu.GetByteLength()/1024.0,(cfr->m_iResXYZ*13.16666)/bscu.GetByteLength());
		}

		// Compare

		if (check) {

			cfr2 = new CBQBCubeFrame(m_IF);
			m_pEngine->PushOutputCubeFrame(cfr2);
			cfr2->m_pAtoms = new CBQBAtomSet(m_IF);
			m_pEngine->PushOutputAtomFrame(cfr2->m_pAtoms);

			m_pEngine->DecompressAtomFrame(
				&bsat,
				BQB_FRAMETYPE_COMPTRAJ_VERSION,
				m_pParmPosVol->GetPositionParameterSet(),
				NULL, // histused
				false,
				true
			);

			m_pEngine->DecompressCubeFrame(
				&bscu,
				BQB_FRAMETYPE_COMPCUBE_VERSION,
				m_pParmPosVol->GetVolumetricParameterSet(),
				NULL, // histused
				false,
				true
			);

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Comparing input and output...\n");

			err = false;

			for (z=0;z<(int)cfr2->m_pAtoms->m_oaAtoms.size();z++)
				for (z2=0;z2<3;z2++)
					if (cfr->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2] != cfr2->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2]) {
						m_IF.eprintf("        Error in atom coordinate %d[%d]: %.6f (%ld) != %.6f (%ld)\n",
							z,z2,cfr->m_pAtoms->m_oaAtoms[z]->m_fCoord[z2],cfr->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2],
							cfr2->m_pAtoms->m_oaAtoms[z]->m_fCoord[z2],cfr2->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2]);
						err = true;
						goto _skerr;
					}

			for (z=0;z<cfr->m_iResXYZ;z++)
				if (!ExpMantisEqual(cfr->m_iaExpo[z],cfr->m_iaMantis[z],cfr2->m_iaExpo[z],cfr2->m_iaMantis[z])) {
					m_IF.eprintf("        Error in volumetric data element %7d: %.10G vs %.10G\n",
						z,cfr->m_faBin[z],cfr2->m_faBin[z]);
					err = true;
					goto _skerr;
				}
	_skerr:
			if (err) {
				if (m_iUseHistory != 0) {
					m_IF.eprintf("Errors occurred. Compressing frame again with history zero.\n");
					err = false;
					m_iUseHistory = 0;
					goto _again;
				} else {
					m_IF.eprintf("Errors occurred despite of zero history. Aborting.\n");
					m_IF.SetLastError( BQB_ERROR_COMPRESS_FAILED );
					return false;
				}
			}

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Done.\n");
		}

		// End Compare

		bshe.Clear();

		bshe.WriteBits(bsat.GetByteLength(),32);
		bshe.WriteBits(bscu.GetByteLength(),32);

		if (m_iFramesWritten == 0)
			ft = BQB_FRAMETYPE_COMPCUBESTART;
		else if (m_iUseHistory == 0)
			ft = BQB_FRAMETYPE_COMPCUBEKEY;
		else
			ft = BQB_FRAMETYPE_COMPCUBE;

		if (m_bContainsVolumetricData) {
			if (cfr->m_iID != -1)
				id = cfr->m_iID;
			else
				id = m_iFramesWritten+1;
		} else
			id = m_iFramesWritten+1;

		m_pBQBFile->CreateShortFrame(
			ft,
			BQB_FRAMETYPE_COMPCUBE_VERSION,
			id
		);

		m_pBQBFile->PushPayload( bshe.m_iaData );
		m_pBQBFile->PushPayload( bsat.m_iaData );
		m_pBQBFile->PushPayload( bscu.m_iaData );
		m_pBQBFile->FinalizeFrame( NULL );

		m_iFramesWritten++;
		m_iUseHistory++;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(
			"CBQBFileWriter::WritePendingFrames(): Finished writing. Now have %d pending volumetric frames and %d pending atom frames.\n",
			GetPendingVolumetricFrameCount(),
			GetPendingAtomPosFrameCount()
		);

	m_IF.SetLastError( BQB_OK );

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::WritePendingFrames(): Done.\n");

	return true;
}



bool CBQBFileWriter::HavePendingVolumetricFrames() const {

	if (m_iVolumetricFramesPending != 0)
		return true;
	else
		return false;
}



int CBQBFileWriter::GetPendingVolumetricFrameCount() const {

	return m_iVolumetricFramesPending;
}



bool CBQBFileWriter::HavePendingAtomPosFrames() const {

	if (m_iAtomPosFramesPending != 0)
		return true;
	else
		return false;
}



int CBQBFileWriter::GetPendingAtomPosFrameCount() const {

	return m_iAtomPosFramesPending;
}



bool CBQBFileWriter::HavePendingCommentLines() const {

	if (m_iCommentLinesPending != 0)
		return true;
	else
		return false;
}



int CBQBFileWriter::GetPendingCommentLineCount() const {

	return m_iCommentLinesPending;
}



bool CBQBFileWriter::WriteIndex() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::WriteIndex(): Writing index frame ...\n");
	else if (m_IF.IsPL( BQB_PL_STANDARD ))
		m_IF.printf("  Writing index frame...\n");

	if (!m_pBQBFile->WriteIndexFrame( true, NULL )) {
		m_IF.eprintf("CBQBFileWriter::WriteIndex(): Error: Could not write index frame.\n");
		return false;
	}

	m_IF.SetLastError( BQB_OK );

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileWriter::WriteIndex(): Done.\n");

	return true;
}


/*******************************************************************************************
******      CBQBFileReader      ************************************************************
*******************************************************************************************/



CBQBFileReader::~CBQBFileReader() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileReader::~CBQBFileReader(): Destructing instance.\n");

	if (m_bActive) {
		m_IF.eprintf("CBQBFileReader::~CBQBFileReader(): Error: Cannot manually destruct instance of this class.\n");
		m_IF.eprintf("This is done automatically when the CBQBInterface is released via BQBReleaseInterface().\n");
		m_IF.eprintf("\n");
		abort();
	}
}



void CBQBFileReader::CloseFile() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileReader::CloseFile(): Closing BQB file ...\n");

	if (!m_bActive) {
		m_IF.eprintf("CBQBFileReader::CloseFile(): Error: File was already closed before.\n");
		m_IF.eprintf("\n");
		abort();
	}
	m_bActive = false;

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileReader::CloseFile(): Done.\n");
}



bool CBQBFileReader::OpenFileRead( const char *s, unsigned long flags ) {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileReader::OpenFileRead(): Opening BQB file for reading (\"%s\", %lu) ...\n",s,flags);

	m_IF.SetLastError( BQB_OK );

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileReader::OpenFileRead(): Done.\n");

	return true;
}



bool CBQBFileReader::ReadFrame( unsigned long flags ) {

	UNUSED(flags);


	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileReader::ReadFrame(): Reading frame ...\n");

	m_IF.SetLastError( BQB_OK );

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBFileReader::ReadFrame(): Done.\n");

	return true;
}



const CBQBVolHeader* CBQBFileReader::GetVolumetricHeader() const {
	if (!m_bContainsVolumetricData) {
		m_IF.eprintf("CBQBFileReader::GetVolumetricHeader(): Error: File does not contain volumetric data.\n");
		m_IF.eprintf("\n");
		abort();
	}
	return m_pVolHeader;
}



const CBQBAtomHeader* CBQBFileReader::GetAtomHeader() const {
	if (!m_bContainsAtomPosition) {
		m_IF.eprintf("CBQBFileReader::GetAtomHeader(): Error: File does not contain atom position data.\n");
		m_IF.eprintf("\n");
		abort();
	}
	return m_pAtomHeader;
}



const CBQBCellInfo* CBQBFileReader::GetCellInfo() const {
	if (!m_bContainsCellInfo) {
		m_IF.eprintf("CBQBFileReader::GetCellInfo(): Error: File does not contain cell info.\n");
		m_IF.eprintf("\n");
		abort();
	}
	return m_pCellInfo;
}



const CBQBVolFrame* CBQBFileReader::GetCurrentVolumetricFrame() const {

	return NULL;
}



const CBQBAtomPosFrame* CBQBFileReader::GetCurrentAtomPosFrame() const {

	return NULL;
}



/*******************************************************************************************
******      CBQBVolHeader      *************************************************************
*******************************************************************************************/



CBQBVolHeader::~CBQBVolHeader() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBVolHeader::~CBQBVolHeader(): Destructing instance.\n");

	if (m_bActive) {
		m_IF.eprintf("CBQBVolHeader::~CBQBVolHeader(): Error: Cannot manually destruct instance of this class.\n");
		m_IF.eprintf("This is done automatically when the CBQBInterface is released via BQBReleaseInterface().\n");
		m_IF.eprintf("\n");
		abort();
	}
}



/*******************************************************************************************
******      CBQBAtomHeader      ************************************************************
*******************************************************************************************/



CBQBAtomHeader::~CBQBAtomHeader() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBAtomHeader::~CBQBAtomHeader(): Destructing instance.\n");

	if (m_bActive) {
		m_IF.eprintf("CBQBAtomHeader::~CBQBAtomHeader(): Error: Cannot manually destruct instance of this class.\n");
		m_IF.eprintf("This is done automatically when the CBQBInterface is released via BQBReleaseInterface().\n");
		m_IF.eprintf("\n");
		abort();
	}
}



/*******************************************************************************************
******      CBQBCellInfo      **************************************************************
*******************************************************************************************/



CBQBCellInfo::~CBQBCellInfo() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBCellInfo::~CBQBCellInfo(): Destructing instance.\n");

	if (m_bActive) {
		m_IF.eprintf("CBQBCellInfo::~CBQBCellInfo(): Error: Cannot manually destruct instance of this class.\n");
		m_IF.eprintf("This is done automatically when the CBQBInterface is released via BQBReleaseInterface().\n");
		m_IF.eprintf("\n");
		abort();
	}
}



/*******************************************************************************************
******      CBQBVolFrame      **************************************************************
*******************************************************************************************/



CBQBVolFrame::~CBQBVolFrame() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBVolFrame::~CBQBVolFrame(): Destructing instance.\n");

	if (m_bActive) {
		m_IF.eprintf("CBQBVolFrame::~CBQBVolFrame(): Error: Cannot manually destruct instance of this class.\n");
		m_IF.eprintf("This is done automatically when the CBQBInterface is released via BQBReleaseInterface().\n");
		m_IF.eprintf("\n");
		abort();
	}
}



void CBQBVolFrame::LoadFromArray( const double *dp, int mode ) {

	unsigned long x, y, z;

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBVolFrame::LoadFromArray(): Loading frame contents from array (mode=%d) ...\n",mode);

	if (mode == BQB_ORDER_ZYX) {

		for (x=0;x<m_pVolHeader->GetResolutionX();x++)
			for (y=0;y<m_pVolHeader->GetResolutionY();y++)
				for (z=0;z<m_pVolHeader->GetResolutionZ();z++)
					m_oCubeFrame.m_faBin[ x*m_oCubeFrame.m_iRes[1]*m_oCubeFrame.m_iRes[2] + y*m_oCubeFrame.m_iRes[2] + z ] =
						dp[ z*m_oCubeFrame.m_iResXY + y*m_oCubeFrame.m_iRes[0] + x ];

	} else if (mode == BQB_ORDER_XYZ) {

		for (x=0;x<m_pVolHeader->GetResolutionX();x++)
			for (y=0;y<m_pVolHeader->GetResolutionY();y++)
				for (z=0;z<m_pVolHeader->GetResolutionZ();z++)
					m_oCubeFrame.m_faBin[ x*m_oCubeFrame.m_iRes[1]*m_oCubeFrame.m_iRes[2] + y*m_oCubeFrame.m_iRes[2] + z ] =
						dp[ x*m_oCubeFrame.m_iRes[1]*m_oCubeFrame.m_iRes[2] + y*m_oCubeFrame.m_iRes[2] + z ];

	} else {

		m_IF.eprintf( "CBQBVolFrame::LoadFromArray(): Error: Unknown mode: %d.\n", mode );
		m_IF.eprintf("\n");
		abort();
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBVolFrame::LoadFromArray(): Done.\n");
}



void CBQBVolFrame::StoreToArray( double *dp, int mode ) const {

	unsigned long x, y, z;

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBVolFrame::StoreToArray(): Storing frame contents to array (mode=%d) ...\n",mode);

	if (mode == BQB_ORDER_ZYX) {

		for (x=0;x<m_pVolHeader->GetResolutionX();x++)
			for (y=0;y<m_pVolHeader->GetResolutionY();y++)
				for (z=0;z<m_pVolHeader->GetResolutionZ();z++)
					dp[ z*m_oCubeFrame.m_iResXY + y*m_oCubeFrame.m_iRes[0] + x ] = 
						m_oCubeFrame.m_faBin[ x*m_oCubeFrame.m_iRes[1]*m_oCubeFrame.m_iRes[2] + y*m_oCubeFrame.m_iRes[2] + z ];

	} else if (mode == BQB_ORDER_XYZ) {

		for (z=0;z<m_pVolHeader->GetResolutionZ();z++)
			for (y=0;y<m_pVolHeader->GetResolutionY();y++)
				for (x=0;x<m_pVolHeader->GetResolutionX();x++)
					dp[ x*m_oCubeFrame.m_iRes[1]*m_oCubeFrame.m_iRes[2] + y*m_oCubeFrame.m_iRes[2] + z ] = 
						m_oCubeFrame.m_faBin[ x*m_oCubeFrame.m_iRes[1]*m_oCubeFrame.m_iRes[2] + y*m_oCubeFrame.m_iRes[2] + z ];

	} else {

		m_IF.eprintf( "CBQBVolFrame::StoreToArray(): Error: Unknown mode: %d.\n", mode );
		m_IF.eprintf("\n");
		abort();
	}

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBVolFrame::StoreToArray(): Done.\n");
}



/*******************************************************************************************
******      CBQBAtomPosFrame      **********************************************************
*******************************************************************************************/



CBQBAtomPosFrame::~CBQBAtomPosFrame() {

	if (m_IF.IsPL( BQB_PL_VERBOSE ))
		m_IF.printf("CBQBAtomPosFrame::~CBQBAtomPosFrame(): Destructing instance.\n");

	if (m_bActive) {
		m_IF.eprintf("CBQBAtomPosFrame::~CBQBAtomPosFrame(): Error: Cannot manually destruct instance of this class.\n");
		m_IF.eprintf("This is done automatically when the CBQBInterface is released via BQBReleaseInterface().\n");
		m_IF.eprintf("\n");
		abort();
	}
}




