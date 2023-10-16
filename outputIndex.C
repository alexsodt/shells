#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pdb.h"
#include "dcd.h"
#include "process_lipids.h"
#include <math.h>

int isProtRes( const char *resname )
{
	if( !strcasecmp( resname, "ALA" ) ) return 1;
	if( !strcasecmp( resname, "CYS" ) ) return 1;
	if( !strcasecmp( resname, "ASP" ) ) return 1;
	if( !strcasecmp( resname, "GLU" ) ) return 1;
	if( !strcasecmp( resname, "PHE" ) ) return 1;
	if( !strcasecmp( resname, "HIS" ) ) return 1;
	if( !strcasecmp( resname, "HSD" ) ) return 1;
	if( !strcasecmp( resname, "HSE" ) ) return 1;
	if( !strcasecmp( resname, "ILE" ) ) return 1;
	if( !strcasecmp( resname, "LYS" ) ) return 1;
	if( !strcasecmp( resname, "LEU" ) ) return 1;
	if( !strcasecmp( resname, "MET" ) ) return 1;
	if( !strcasecmp( resname, "ASN" ) ) return 1;
	if( !strcasecmp( resname, "PRO" ) ) return 1;
	if( !strcasecmp( resname, "GLN" ) ) return 1;
	if( !strcasecmp( resname, "ARG" ) ) return 1;
	if( !strcasecmp( resname, "SER" ) ) return 1;
	if( !strcasecmp( resname, "THR" ) ) return 1;
	if( !strcasecmp( resname, "VAL" ) ) return 1;
	if( !strcasecmp( resname, "TRP" ) ) return 1;
	if( !strcasecmp( resname, "TYR" ) ) return 1;

	return 0;
}

int main( int argc_in, char **argv_in )
{
	char buffer[4096];
	
	char *argv[argc_in];
	int argc = 0;
	int doChains = 0;

	for( int c = 0; c < argc_in; c++ )
	{
		if( !strncasecmp( argv_in[c], "--", 2) )
		{
			if( !strcasecmp( argv_in[c], "--chain") )
			{	
				doChains=1;
			}
			else
			{
				printf("Unknown flag '%s'.\n", argv_in[c] );
				exit(1);
			}
		}
		else
		{
			argv[argc] = argv_in[c];
			argc++;
		}
	}

	if( argc < 2 )
	{
		printf("Syntax: outputIndex [--chain] pdb\n");	
		return 0;
	}
	
	if( 1 >= argc )
	{
		printf("PDB required.\n");
		exit(1);
	}

	FILE *psfFile = fopen(argv[1], "r" );
	if( ! psfFile )
	{
		printf("Couldn't open PDB file '%s'.\n", argv[1] );
		return 0;
	}

	loadPSFfromPDB( psfFile );    

	rewind(psfFile);

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );

	loadPDB( psfFile, at );
			
	char pseg[256] = {'\0'};
	int pres = -1;

	int n_encoded = 0;

	for( int a = 0; a < curNAtoms(); a++ )
	{
		if( at[a].res == pres && (!at[a].segid || strlen(at[a].segid) == 1 || !strcasecmp( at[a].segid, pseg )) ) continue;


		pres = at[a].res;
		strcpy( pseg, at[a].segid );

		// get one unique atom name
		if( lipidType(at[a].resname) != -1 )
		{
			if(doChains )
			{
				for( int c = 0; c < nChains( at[a].resname ); c++ )
					printf("%s_%d\n", at[a].resname, c );
			}
			else
			{
				printf("%s\n", at[a].resname );
			}
	
			if( doChains )
				n_encoded += nChains( at[a].resname );
			else
				n_encoded += 1;
		}
	
		if( (isProtRes(at[a].resname) || (at[a].segid && !strncasecmp( at[a].segid, "PRO", 3)) || !strcasecmp( at[a].atname, "BB") )  && fabs( at[a].z) > 6 && fabs(at[a].z) < 18 )
		{
//			printf("%s\n", at[a].resname );
			n_encoded += 1;
		}
	}


	fclose(psfFile);	


}
