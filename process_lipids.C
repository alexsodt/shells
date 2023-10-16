#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int isAltSterol( const char * res )
{
if( !strcasecmp( res, "DHC7") ||
!strcasecmp( res, "ERG") ||
!strcasecmp( res, "CSUL") ||
!strcasecmp( res, "4BOH") ||
!strcasecmp( res, "25OH") ||
!strcasecmp( res, "24OH") ||
!strcasecmp( res, "CCHL") ||
!strcasecmp( res, "SMOOV") ||
!strcasecmp( res, "ZYM") ||
!strcasecmp( res, "DHC8") ||
!strcasecmp( res, "LATH") ||
!strcasecmp( res, "DESM") ||
!strcasecmp( res, "27OH") ||
!strcasecmp( res, "CHONE") ||
!strcasecmp( res, "7KETO") ||
!strcasecmp( res, "VD3") ||
!strcasecmp( res, "LANO") ||
!strcasecmp( res, "9THC") ||
!strcasecmp( res, "VITE") )
        return 1;
        return 0;
}

int nChains( char *resname )
{
	if( !strcasecmp( resname, "DNPC" ) ) return 2;
	if( !strcasecmp( resname, "DXPC" ) ) return 2;
	if( !strcasecmp( resname, "DYPC" ) ) return 2;
	if( !strcasecmp( resname, "DGPC" ) ) return 2;
	if( !strncasecmp( resname, "SAPI", 4 ) )
		return 2;
	if( !strcasecmp( resname, "PSM") ) return 2;
	if( !strncasecmp( resname, "CHL", 3 ) ) return 1;
	if( isAltSterol( resname ) ) return 1;

	if( !strncasecmp( resname, "CHOL", 4 ) ) return 1;
	if( !strncasecmp( resname, "CER1", 4 ) ) return 2;
	if( !strncasecmp( resname, "DHC7", 4) ) return 1;


	if( !strcasecmp( resname, "PPPE" ) ) return 2;
	if( !strcasecmp( resname, "PVPG" ) ) return 2;

	if( !strcasecmp( resname, "PVCL2" ) ) return 4;	
	if( !strcasecmp( resname, "PVCL" ) ) return 4;	

	if( strlen( resname ) == 4 )
	{
		// some reasonable standard for it probably being a lipid, override however you like.

		int isPC = !strncasecmp( resname + 2, "PC",2 );
		int isPS = !strncasecmp( resname + 2, "PS",2 );
		int isPE = !strncasecmp( resname + 2, "PE",2 );
		int isPI = !strncasecmp( resname + 2, "PI",2 );
		int isPG = !strncasecmp( resname + 2, "PG",2 );
	
		int isDi = resname[0] == 'D';

		int isOleoyl = resname[1] == 'O';

		int isPalm = resname[1] == 'P';

		if( (isPC || isPE || isPS || isPI || isPG ) && (isDi || isOleoyl || isPalm ) )
			return 2;

		if( (isPC || isPE || isPS || isPI ) && (!strncasecmp(resname, "SA",2) ||
					       !strncasecmp(resname, "PL",2) ||
					       !strncasecmp(resname, "PA",2) || !strncasecmp( resname, "SD", 2)  ) )
			return 2;

		if( !strncasecmp( resname, "MM", 2) ) 
			return 2;
		if( !strncasecmp( resname, "PP", 2) ) 
			return 2;
		if( !strncasecmp( resname, "SO", 2) ) 
			return 2;
	}
	return 0; // not a lipid.
}

int isPI( char *resname )
{
	if( strlen(resname) > 2 && !strncasecmp( resname+2, "PI", 2 ) )
		return 1;
	return 0;

}


int isChl( char *resname )
{
	if( !strncasecmp( resname, "CHL", 3 ) ) return 1;
	if( !strncasecmp( resname, "CHOL", 4 ) ) return 1;
	if( !strncasecmp( resname, "DHC7", 4) ) return 1;
	return 0;

}


int isSphingo( char *resname )
{
	if( !strcasecmp( resname, "SSM") ) return 1;
	if( !strcasecmp( resname, "PSM") ) return 1;
	if( !strncasecmp( resname, "CER1", 4) ) return 1;

	return 0;
}

int isPE( char *resname )
{
	if( !strcasecmp( resname, "PSM") ) return 1;
	if( !strncasecmp( resname, "CHL", 3 ) ) return 0;
	if( !strncasecmp( resname, "CHOL", 4 ) ) return 0;

	if( strlen( resname ) == 4 )
	{
		// some reasonable standard for it probably being a lipid, override however you like.

		int lipid_isPE = !strncasecmp( resname + 2, "PE",2 );

		return lipid_isPE;
	}
	return 0;
}

int isPS( char *resname )
{
	if( !strcasecmp( resname, "PSM") ) return 1;
	if( !strncasecmp( resname, "CHL", 3 ) ) return 0;
	if( !strncasecmp( resname, "CHOL", 4 ) ) return 0;

	if( strlen( resname ) == 4 )
	{
		// some reasonable standard for it probably being a lipid, override however you like.

		int lipid_isPS = !strncasecmp( resname + 2, "PS",2 );

		return lipid_isPS;
	}
	return 0; 
}

int isSaturated( char *resname )
{
	if( !strcasecmp( resname, "PSM") ) return 1;
	if( !strncasecmp( resname, "CHL", 3 ) ) return 0;
	if( !strncasecmp( resname, "CHOL", 4 ) ) return 0;

	if( strlen( resname ) == 4 )
	{
		// some reasonable standard for it probably being a lipid, override however you like.

		int isPC = !strncasecmp( resname + 2, "PC",2 );
		int isPE = !strncasecmp( resname + 2, "PE",2 );
		int isPS = !strncasecmp( resname + 2, "PS",2 );
	
		int isDi = resname[0] == 'D';

		int isOleoyl = resname[1] == 'O';

		int isPalm = resname[1] == 'P';

		if( isOleoyl ) return 0;

		if( (isPC || isPE || isPS ) && (isDi || isPalm ) )
			return 1;
	}
	return 0;
}

int whichChain (char *resName, char *atomName )
{
	if( nChains(resName) == 0 ) return -1;

	int sl = strlen(atomName)-1;

	if( !strcasecmp( resName, "PSM" ) || !strncasecmp( resName, "CER1",4 ) )
	{
		if( !strcasecmp( atomName, "C1F") ) return 1;
		if( !strcasecmp( atomName, "OF") ) return 1;
		if( !strcasecmp( atomName, "C2S") ) return 1;

		if( atomName[sl] == 'S' ) return 0;	
		if( atomName[sl] == 'F' ) return 1;	
	}
	else if( !strncasecmp( resName, "CHL", 3 ) || isAltSterol( resName ) )
	{
		return 0;
	}
	else if( !strncasecmp( resName, "CHOL", 4 ) )
	{
		return 0;
	}
	else if( !strcasecmp( resName, "PVCL2") || !strcasecmp( resName, "PVCL") )
	{
		if( atomName[0] != 'C' ) return -1;
			
		if( atomName[1] == 'A' ) return 0;
		if( atomName[1] == 'B' ) return 1;
		if( atomName[1] == 'C' ) return 2;
		if( atomName[1] == 'D' ) return 3;
	}
	else 
	{	
		if( atomName[0] != 'C' ) return -1;
	
		if( 
			!strcasecmp( atomName, "C1A" ) ||
			!strcasecmp( atomName, "C2A" ) ||
			!strcasecmp( atomName, "C3A" ) ||
			!strcasecmp( atomName, "C4A" ) ||
			!strcasecmp( atomName, "C5A" ) ||
			!strcasecmp( atomName, "D1A" ) ||
			!strcasecmp( atomName, "D2A" ) ||
			!strcasecmp( atomName, "D3A" ) ||
			!strcasecmp( atomName, "D4A" ) ||
			!strcasecmp( atomName, "D5A" ) ) 
			return 0;
		else if( 
			!strcasecmp( atomName, "C1B" ) ||
			!strcasecmp( atomName, "C2B" ) ||
			!strcasecmp( atomName, "C3B" ) ||
			!strcasecmp( atomName, "C4B" ) ||
			!strcasecmp( atomName, "C5B" ) ||
			!strcasecmp( atomName, "D1B" ) ||
			!strcasecmp( atomName, "D2B" ) ||
			!strcasecmp( atomName, "D3B" ) ||
			!strcasecmp( atomName, "D4B" ) ||
			!strcasecmp( atomName, "D5B" ) ) 
			return 1;
		

		if( sl > 1 )
		{
			if( atomName[1] == '2' ) return 0;
			if( atomName[1] == '3' ) return 1;
		}
	}
	return -1; 
}


int lipidType( const char * resname )
{
	if( !strcasecmp( resname, "DXPC") ) return 0;
	if( !strcasecmp( resname, "DYPC") ) return 0;
	if( !strcasecmp( resname, "DGPC") ) return 1;
	if( !strcasecmp( resname, "DNPC") ) return 2;
	if( !strcasecmp( resname, "PPPE") ) return 0; 
	if( !strcasecmp( resname, "PVPG") ) return 1; 
	if( !strcasecmp( resname, "PVCL") ) return 2; 
	if( !strcasecmp( resname, "PVCL2") ) return 2; 

	if( !strcasecmp( resname, "DSM") ) return 0;
	if( !strcasecmp( resname, "DPPC") ) return 0;
	if( !strcasecmp( resname, "PPPC") ) return 0;
	if( !strcasecmp( resname, "MMPC") ) return 0;
	if( !strcasecmp( resname, "SOPC") ) return 0;

	if( !strcasecmp( resname, "CHOL") ) return 0;
	if( !strcasecmp( resname, "DOPC") ) return 0;

	if( !strcasecmp( resname, "DOPE") ) return 0;
	if( !strcasecmp( resname, "DOPS") ) return 0;
	if( !strcasecmp( resname, "POPS") ) return 0;
	if( !strcasecmp( resname, "DOPG") ) return 0;
	if( !strcasecmp( resname, "CHL1") ) return 1;
	if( !strncasecmp( resname, "SAPI", 4 ) ) return 2;
	
	if( !strncasecmp(resname, "LSM", 3 ) ) return 0;
	if( !strncasecmp(resname, "NSM", 3 ) ) return 0;
	if( !strncasecmp(resname, "OAPE", 4 ) ) return 0;
	if( !strncasecmp(resname, "OAPS", 4 ) ) return 0;
	if( !strncasecmp(resname, "PAPC", 4 ) ) return 0;
	if( !strncasecmp(resname, "PAPS", 4 ) ) return 0;
	if( !strncasecmp(resname, "PDPE", 4 ) ) return 0;
	if( !strncasecmp(resname, "PLAO", 4 ) ) return 0;
	if( !strncasecmp(resname, "PLAS", 4 ) ) return 0;
	if( !strncasecmp(resname, "PLPC", 4 ) ) return 0;
	if( !strncasecmp(resname, "PLQS", 4 ) ) return 0;
	if( !strncasecmp(resname, "POPC", 4 ) ) return 0;
	if( !strncasecmp(resname, "POPE", 4 ) ) return 0;
	if( !strncasecmp(resname, "PSM", 3 ) ) return 0;
	if( !strncasecmp(resname, "SAPS", 4 ) ) return 0;
	if( !strncasecmp(resname, "SAPC", 4 ) ) return 0;
	if( !strncasecmp(resname, "SDPC", 4 ) ) return 0;
	if( !strncasecmp(resname, "SAPE", 4 ) ) return 2;
	if( !strncasecmp(resname, "SDPE", 4 ) ) return 2;
	return -1;
}




