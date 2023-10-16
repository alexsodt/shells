#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include <math.h>

int main( int argc, char **argv )
{
	int nlip = 0;
	char **lipidNames = NULL;
 	int init = 1;

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	double *read = NULL;
	double *data  = NULL;
	double *data2 = NULL;
	int nphi = 0;
	int nphiSpace = 12;
	int nfiles = argc-1;
	for( int c = 1; c < argc; c++ )
	{
		FILE *theFile = fopen(argv[c],"r");

		if( !theFile )
		{
			printf("Couldn't open file '%s'.\n", argv[c] );
			exit(1);
		}

		getLine( theFile, buffer );


		char *t = buffer;	

		if( strncasecmp(t, "#Phi_bin", 8 ) )
		{
			printf("File '%s' has wrong format; expected '#Phi_bin' at first line.\n", argv[c] );
			exit(1);
		}

		t += strlen("#Phi_bin ");
	
		int llip = 0;
		int llipSpace = 10;

		char **lnames = (char **)malloc( sizeof(char *) * llipSpace );

		while( *t )
		{	
			if( llip == llipSpace )
			{
				llipSpace *= 2;
				lnames = (char **)realloc( lnames, sizeof(char *) * llipSpace ); 
			}
			lnames[llip] = (char *)malloc( sizeof(char) * (1+strlen(t)) );

			int s = 0;

			while( *t && *t != ' ' && *t != '\t' ) 
			{
				lnames[llip][s++] = *t;
				lnames[llip][s] = '\0';
				t+=1;
			}

			llip++;
		
			while( *t == ' ' || *t == '\t' ) t += 1;
		}

		int mapping[llip];

		for( int i = 0; i < llip; i++ )
			mapping[i] = i;

		if( init )
		{
			lipidNames = lnames;
			nlip = llip;
		} 
		else
		{
			if( nlip != llip )
			{
				printf("NLIP mismatch in file '%s'.\n", argv[c] );
				exit(1);
			}
		
			for( int i = 0; i < llip; i++ )
				mapping[i] = -1;

			for( int a = 0; a < nlip; a++ )
			for( int b = 0; b < llip; b++ )
			{
				if( !strcasecmp( lipidNames[a], lnames[b] ) )
					mapping[b] = a;
			}	
	
			for( int b = 0; b < llip; b++ )
			{
				if( mapping[b] == -1 )
				{
					printf("Couldn't map '%s' in file '%s'.\n", lnames[b], argv[c] );
					exit(1);
				}
			}
		}

		if( init )
		{
			data = (double *)malloc( sizeof(double) * nlip * nphiSpace );
			data2 = (double *)malloc( sizeof(double) * nlip * nphiSpace );
			memset( data, 0, sizeof(double) * nlip * nphiSpace ); 
			memset( data2, 0, sizeof(double) * nlip * nphiSpace ); 
			read = (double *)malloc( sizeof(double) * (1 + nlip ) );
		}

		int lphi = 0;

		while( !feof(theFile) )
		{	
			getLine( theFile, buffer );
			if( feof(theFile) ) break;


			int nr = readNDoubles( buffer, read, 1 + nlip );

			if( nr != 1 + nlip )
				break;

			if( lphi == nphiSpace )
			{
				nphiSpace *= 2;
				data = (double *)realloc( data, sizeof(double) * nlip * nphiSpace );
				data2 = (double *)realloc( data2, sizeof(double) * nlip * nphiSpace );
			}		

			if( lphi >= nphi )
			{
				for( int x = 0; x < nlip; x++ )
				{
					data[lphi*nlip+mapping[x]] = read[1+x];
					data2[lphi*nlip+mapping[x]] = read[1+x]*read[1+x];
				}
			}
			else
			{
				for( int x = 0; x < nlip; x++ )
				{
					data[lphi*nlip+mapping[x]] += read[1+x];
					data2[lphi*nlip+mapping[x]] += read[1+x]*read[1+x];
				}
			}

			lphi++;

			if( lphi >= nphi )
				nphi = lphi;
		}

		init = 0;

		fclose(theFile);
	}

	printf("#Phi_bin");
	for( int l = 0; l < nlip; l++ )
	{
		if( !strcasecmp( lipidNames[l] + strlen(lipidNames[l])-3, "err") )
			continue;
		printf(" %s", lipidNames[l] );
	}
	for( int l = 0; l < nlip; l++ )
	{
		if( !strcasecmp( lipidNames[l] + strlen(lipidNames[l])-3, "err") )
			continue;
		printf(" %s_err", lipidNames[l] );
	}
	printf("\n");		

	for( int p = 0; p < nphi; p++ )
	{
		printf("%d", p );

		for( int l = 0; l < nlip; l++ )	
		{
			if( !strcasecmp( lipidNames[l] + strlen(lipidNames[l])-3, "err") )
				continue;
			double av = data[p*nlip+l] / nfiles;
			double av2 = data2[p*nlip+l] / nfiles;
			double stde = sqrt(av2-av*av)/sqrt(nfiles-1);

			printf(" %lf", av );
		}
		
		for( int l = 0; l < nlip; l++ )	
		{
			if( !strcasecmp( lipidNames[l] + strlen(lipidNames[l])-3, "err") )
				continue;
			double av = data[p*nlip+l] / nfiles;
			double av2 = data2[p*nlip+l] / nfiles;
			double stde = sqrt(av2-av*av)/sqrt(nfiles-1);

			printf(" %lf",stde );
		}
		printf("\n");
	}
	


}
