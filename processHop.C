// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#define SUB_PRESS
#define LEN 100
#include "border_config.h"
#include "borders.h"

//#define MAX_BORDERS 25
//
//struct border_header
//{
//	int nencoded;
//	int origResID[MAX_ENCODED];
//	int origType[MAX_ENCODED];
//};
//
//struct border_data
//{
//	int nborders;
//	int bpartners[MAX_BORDERS];
//	float bvals[MAX_BORDERS];
//	int is_interacting_with_positive;
//	int which_chain;
//};

int main( int argc_in, char **argv_in )
{
	char *argv[argc_in];	
	int argc = 0;
	FILE *phiOut = NULL;

	for( int c = 0; c < argc_in; c++ )
	{
		if( !strncasecmp( argv_in[c], "--", 2) )
		{
			if( !strncasecmp( argv_in[c], "--phiout=", 9 ) )
				phiOut = fopen(argv_in[c]+9, "w");			
		}
		else
		{
			argv[argc] = argv_in[c];
			argc++;
		}
	}

	char buffer[4096];

	if( argc != 3 )
	{
		printf("Syntax: processHop border.bin thresh\n");	
		return 0;
	}

	FILE *borderBin = fopen( argv[1], "rb");
	double thresh = atof(argv[2]);

	border_header theHeader;

	fread( &theHeader, sizeof(border_header), 1, borderBin );
	int nframes = 0;
	border_data *frame_read = (border_data *)malloc( sizeof(border_data) * theHeader.nencoded );

	while( !feof(borderBin) )
	{
		int nr = fread( frame_read, sizeof(border_data), theHeader.nencoded, borderBin );
		if( nr == theHeader.nencoded )
			nframes++;
	}
	rewind(borderBin);
	fread( &theHeader, sizeof(border_header), 1, borderBin );

	free( frame_read );

	frame_read = (border_data *)malloc( sizeof(border_data) * theHeader.nencoded * nframes );

	int i = 0;

	while( !feof(borderBin) && i < nframes )
	{
		int nr = fread( frame_read + i * theHeader.nencoded, sizeof(border_data), theHeader.nencoded, borderBin );
		if( nr == theHeader.nencoded ) i++;
	}

	int n_encoded = theHeader.nencoded;

	for( int t_master = 0; t_master < nframes; t_master++ )
	{
/*		printf("res7 ");
		for( int ix = 0; ix < frame_read[t_master * n_encoded+14].nborders; ix++ )
		{
			int j  = frame_read[t_master*n_encoded+14].bpartners[ix];
			if( theHeader.origType[j] == 0  && j == 6)
				printf("%d ", j );
		}
		printf(" len ");
		for( int ix = 0; ix < frame_read[t_master * n_encoded+14].nborders; ix++ )
		{
			int j  = frame_read[t_master*n_encoded+14].bpartners[ix];
			if( theHeader.origType[j] == 0 && j == 6)
				printf("%lf ", frame_read[t_master*n_encoded+14].bvals[ix] );
		}
		printf("\n");
*/
		for( int i = 0; i < n_encoded; i++ )
		{
			for( int ix = 0; ix < frame_read[t_master*n_encoded+i].nborders; ix++ )
			{
				int j = frame_read[t_master*n_encoded+i].bpartners[ix];

				double val = frame_read[t_master*n_encoded+i].bvals[ix];

				if( val > thresh )
				{
					frame_read[t_master*n_encoded+i].bvals[ix] = -1;

					int done = 0;
					for( int t_backtrack = t_master-1; !done && t_backtrack >= 0; t_backtrack-- )
					{
						int got_j = 0;
						for( int ix2 = 0; !done && ix2 < frame_read[t_backtrack*n_encoded+i].nborders; ix2++ )
						{
							if( frame_read[t_backtrack*n_encoded+i].bpartners[ix2] == j )
							{
								got_j = 1;
								if( frame_read[t_backtrack*n_encoded+i].bvals[ix2] >= 0 )
									frame_read[t_backtrack*n_encoded+i].bvals[ix2] = -1;
								else 
									done = 1;
							}	
						}
						if( !got_j ) done = 1;
					}
				}
				else if( t_master > 0 )
				{
					int t_backtrack = t_master-1;
					for( int ix2 = 0; ix2 < frame_read[t_backtrack*n_encoded+i].nborders; ix2++ )
					{
						if( frame_read[t_backtrack*n_encoded+i].bpartners[ix2] == j )
						{
							if( frame_read[t_backtrack*n_encoded+i].bvals[ix2] < 0  )
								frame_read[t_master*n_encoded+i].bvals[ix] = -1;
						}	
					}
				} 
			}
		}
	}

	int *solvation_shell = (int*)malloc( sizeof(int) * n_encoded );
	int *borders = (int*)malloc( sizeof(int) * n_encoded * n_encoded ); 
	int *nborders = (int*)malloc( sizeof(int) * n_encoded );

	FILE *isWithPos = fopen("interacting_with_pos.txt","w");

	for( int t = 0; t < nframes; t++ )
	{
		int *at_code = theHeader.origType;
		int *orig_id = theHeader.origResID; 

		memset( nborders, 0, sizeof(int) * n_encoded );

		for( int i = 0; i < n_encoded; i++ )
		{
			for( int ix = 0; ix < frame_read[t*n_encoded+i].nborders; ix++ )
			{
				int j = frame_read[t*n_encoded+i].bpartners[ix];

				if( frame_read[t*n_encoded+i].bvals[ix] < 0 )
				{
					borders[i*n_encoded+nborders[i]] = j;
					nborders[i] += 1;
				}
			} 
		}

		for( int i = 0; i < n_encoded; i++ )
		{
			if( at_code[i] == 0 )
				solvation_shell[i] = 0;
			else
				solvation_shell[i] = -1;
		}			
	
		int done = 0;
		
		int xz = 0;
		while( !done )
		{
			done = 1;
			for( int i = 0; i < n_encoded; i++ )
			{
		
				int min_border = 10000;
				for( int p = 0; p < nborders[i]; p++ )
				{
					if( solvation_shell[borders[i*n_encoded+p]] != -1 )
					{
						if(  solvation_shell[borders[i*n_encoded+p]] < min_border )
							min_border = solvation_shell[borders[i*n_encoded+p]];
					}
				}
				if( min_border < 10000 )
				{
					if( solvation_shell[i] == -1 || min_border+1 < solvation_shell[i] ) 
					{
						solvation_shell[i] = min_border+1;	
						done=0;
					}
				}
			}
	
		}
		
		for( int i = 0; i < n_encoded; i++ )
		{
			int min_border = 10000;
		
			for( int p = 0; p < nborders[i]; p++ )
			{
				if( solvation_shell[borders[i*n_encoded+p]] < min_border )
					min_border = solvation_shell[borders[i*n_encoded+p]];
			}
			if( min_border >= solvation_shell[i] && solvation_shell[i] > 0 )
			{
				printf("Logical error.\n");
			}
		}

		for( int i = 0; i < n_encoded; i++ )
		{
			if( theHeader.origType[i] > 0 )
			{
				if( phiOut ) fprintf(phiOut, " %lf", frame_read[t*n_encoded+i].phi );
				printf("%d ", solvation_shell[i] ); 
				fprintf(isWithPos, "%d ", frame_read[t*n_encoded+i].is_interacting_with_positive ); 
			}
		}
		fprintf(isWithPos, "\n");
		printf("\n");
		if( phiOut ) fprintf(phiOut, "\n");
		fflush(stdout);
		fflush(isWithPos);
	}
	fclose(isWithPos);
}






