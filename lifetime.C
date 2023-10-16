#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

double pos_thresh = 0.5;

#define MAX_SHELLS 20

int main ( int argc, char ** argv )
{
	if( argc < 3 )
	{
		printf("Syntax lifetime index.file shells.file [interacting_with_pos.txt|null] [filter_width]\n");
		return -1;
	}

	FILE *index = fopen(argv[1],"r");
	FILE *shells = fopen(argv[2],"r");
	FILE *posFile = NULL;
	int start =0;
	int stop = -1;
	int filter_width = 1;

	if( argc > 3 )
	{
		if( strcasecmp( argv[3], "NULL" ) )
		{
			posFile = fopen( argv[3], "r" );

			if( !posFile )
			{
				printf("Couldn't open pos file '%s'.\n", argv[3] );
				exit(1);
			}
		}
	}

	if( argc > 3 ) 
		filter_width = atoi(argv[4]);
	
	if( !index )
	{
		printf("Couldn't open index file '%s'.\n", argv[1] );
		return -1;
	}
	
	if( !shells )
	{
		printf("Couldn't open index file '%s'.\n", argv[2] );
		return -1;
	}
	
	struct ltype
	{
		double *hist;
		char name[256];
		int *indices;
		int nind;
		int nindSpace;
		struct ltype *next;
	};

	struct ltype *lipids = NULL;

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	int ind = 0;
	while( !feof(index) )
	{
		getLine( index, buffer );
		if( feof(index) ) break;
	
		struct ltype *got = NULL;

		for( struct ltype *l = lipids; l; l = l->next )
		{
			if( !strcasecmp( l->name, buffer ) ) 
				got = l;
		}

		if( ! got )
		{
			got = (struct ltype *)malloc( sizeof(struct ltype)  );
			if( strlen(buffer)<256)
				strcpy( got->name, buffer );
			else
				strcpy( got->name, "ERROR");
			got->nind = 0;
			got->nindSpace = 10;			
			got->indices = (int *)malloc( sizeof(int) * got->nindSpace );
			got->hist = (double *)malloc( sizeof(double ) * MAX_SHELLS );
			got->next = lipids;
			lipids = got;
			memset( got->hist, 0, sizeof(double) * MAX_SHELLS );
		}			

		if( got->nind == got->nindSpace )
		{
			got->nindSpace *= 2;
			got->indices = (int *)realloc( got->indices, sizeof( int) * got->nindSpace );
		}

		got->indices[got->nind] = ind; 
		got->nind++;
		ind++;
	}

	int n_index = ind;

	int *shell_read = (int *)malloc( sizeof(int) * n_index );
	int *pos_read = (int *)malloc( sizeof(int) * n_index );
	memset( pos_read, 0, sizeof(int) * n_index );


	int max_time=0;
	int time_space = 10;

	char *shell_total = (char *)malloc( sizeof(char) * n_index * time_space ); 
	char *pos_total = (char *)malloc( sizeof(char) * n_index * time_space ); 

	while( !feof(shells) ) 
	{
		getLine( shells, buffer);
		if( feof(shells) ) break;

		int nr = readNInts( buffer, shell_read, n_index );

		if( nr != n_index )
		{
			printf("Failed to read %d shells (read %d).\n", n_index, nr );
			exit(1);
		}
		
		printf("%d", max_time );
		
		if( posFile )
		{		
			getLine( posFile, buffer);
			if( feof(posFile) ) break;

			if( max_time == 10000 )
			{
				printf("Debug this.\n");
			}
			int p_nr = readNInts( buffer, pos_read, n_index );
			if( p_nr != n_index )
			{	
				printf("Failed to read %d positive-res indicators (read %d).\n", n_index, p_nr );
			}
		}

		if( max_time == time_space )
		{
			time_space *= 2;
			shell_total = (char *)realloc( shell_total, sizeof(char) * n_index * time_space );
			pos_total = (char *)realloc( pos_total, sizeof(char) * n_index * time_space );
		}

		for( int i = 0; i < n_index; i++ )
		{
			shell_total[max_time*n_index+i] = (char) shell_read[i];
			pos_total[max_time*n_index+i] = (char) pos_read[i];
		}
	
		for( struct ltype *l = lipids; l; l = l->next )
		{
			if( strncasecmp( l->name, "SAPI",4 ) ) 
				continue;
			for( int ind = 0; ind < l->nind; ind++ )
			{
				printf(" %d", shell_total[max_time*n_index+l->indices[ind]]);
			}	
		}
		printf("\n");
		max_time++;
	}

	char *filtered = (char *)malloc( sizeof(char) * max_time );
	char *shell_total_filtered = (char *)malloc( sizeof(char) * n_index * max_time ); 

	int time_min = 0;
	int time_max = 3000;
	int nbins = 100;
	

	double lifetime_PIP2_in[nbins];
	double lifetime_PIP2_out[nbins];

	memset( lifetime_PIP2_in, 0, sizeof(double) * nbins );
	memset( lifetime_PIP2_out, 0, sizeof(double) * nbins );


	for( struct ltype *l = lipids; l; l = l->next )
	{
		for( int ind = 0; ind < l->nind; ind++ )
		{
			int x = l->indices[ind];

#ifdef LOW_PASS
			for( int t = 0; t < max_time; t++ )
			{
				int t_min = t - filter_width/2;
				int t_max = t + filter_width/2;


				if( t_min < 0 ) t_min = 0;
				if( t_max >= max_time ) t_max = max_time-1;

				double av_shell = 0;
				for( int t2 = t_min; t2 <= t_max; t2++ )
					av_shell += shell_total[t2*n_index+x];

				av_shell /= (t_max-t_min+1);
				long the_shell = lround(av_shell);

				filtered[t] = (char)the_shell;
			} 
#else
			for( int t = 0; t < max_time; t++ )
				filtered[t] = shell_total[t*n_index+x];
#endif			
			int done = 0;

			while( !done )
			{
				done = 1;

				int min_run_start = -1;
				int min_run_stop = -1;
				int min_run_len = max_time;

				char ppshell = -1;
				char pshell = filtered[0];
				int t_mark = 0;
				int switch_to = pshell;

				// remove runs shorter than max_time.
				for( int t = 1; t < max_time; t++ )
				{
					if( filtered[t] != pshell )
					{
						if( t -t_mark < min_run_len && t_mark > 0 && ppshell == filtered[t] )
						{
							min_run_len = t-t_mark;
							min_run_start = t_mark;
							min_run_stop = t;
							switch_to = filtered[t];
						}

						pshell = filtered[t];
						ppshell = filtered[t-1];
						t_mark = t;						
					}	
				} 
				if( min_run_len <filter_width )
				{
					done = 0;
					for( int t = min_run_start; t < min_run_stop; t++ )
						filtered[t] = switch_to;
				}
			}
			
		
			for( int t = 0; t < max_time; t++ )
			{
//				if( ind == 0 )  
//					printf("%d %d\n", t, filtered[t] );
				shell_total_filtered[x*max_time+t] = filtered[t];	
			}
//			printf("\n");
		} 
	}

	int max_shell_lifetime = 6;

	for( struct ltype *l = lipids; l; l = l->next )
	{
		double nav_lifetime_end[1+max_shell_lifetime];
		double av_lifetime_end[1+max_shell_lifetime];
		double nav_lifetime_in[1+max_shell_lifetime];
		double nav_lifetime_out[1+max_shell_lifetime];
		double av_lifetime_in[1+max_shell_lifetime];
		double av_lifetime_out[1+max_shell_lifetime];
		memset( av_lifetime_end, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( av_lifetime_in, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( av_lifetime_out, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( nav_lifetime_end, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( nav_lifetime_in, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( nav_lifetime_out, 0, sizeof(double) * (1+max_shell_lifetime) );
		
		double nav_lifetime_end_pos[1+max_shell_lifetime];
		double nav_lifetime_in_pos[1+max_shell_lifetime];
		double nav_lifetime_out_pos[1+max_shell_lifetime];
		double av_lifetime_end_pos[1+max_shell_lifetime];
		double av_lifetime_in_pos[1+max_shell_lifetime];
		double av_lifetime_out_pos[1+max_shell_lifetime];
		memset( av_lifetime_end_pos, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( av_lifetime_in_pos, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( av_lifetime_out_pos, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( nav_lifetime_end_pos, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( nav_lifetime_in_pos, 0, sizeof(double) * (1+max_shell_lifetime) );
		memset( nav_lifetime_out_pos, 0, sizeof(double) * (1+max_shell_lifetime) );


		for( int ind = 0; ind < l->nind; ind++ )	
		{
			int x = l->indices[ind];

			char pshell = shell_total_filtered[x*max_time+0];
			int t_mark = 0;
			int n_pos = 0;

			for( int t = 1; t <= max_time; t++ )
			{
				if( t == max_time )
				{
					int lifetime = t - t_mark;
					int from = pshell;
						
					av_lifetime_end[from] += lifetime;
					nav_lifetime_end[from] += 1;
					if( n_pos > pos_thresh * lifetime)
					{
						av_lifetime_end_pos[from] += lifetime;
						nav_lifetime_end_pos[from] += 1;
					}
				}
				else if( shell_total_filtered[x*max_time+t] != pshell )
				{
					int lifetime = t - t_mark;
					int from = pshell;
					int to = shell_total_filtered[x*max_time+t];

					if( from < to && from <= max_shell_lifetime )
					{
						int tbin = nbins*(lifetime-time_min)/(double)(time_max-time_min);
						if( tbin >= nbins ) tbin = nbins-1;

						if( !strncasecmp( l->name, "SAPI", 4 ) && from < 3 )
							lifetime_PIP2_out[tbin] += 1;	

						av_lifetime_out[from] += lifetime;
						nav_lifetime_out[from] += 1;
						if( n_pos > pos_thresh * lifetime)
						{
							av_lifetime_out_pos[from] += lifetime;
							nav_lifetime_out_pos[from] += 1;
						}
					}
					else if( from > to && from <= max_shell_lifetime )
					{
						int tbin = nbins*(lifetime-time_min)/(double)(time_max-time_min);
						if( tbin >= nbins ) tbin = nbins-1;

						if( !strncasecmp( l->name, "SAPI", 4 ) && to < 3 )
							lifetime_PIP2_in[tbin] += 1;	
						av_lifetime_in[from] += lifetime;
						nav_lifetime_in[from] += 1;

//						printf("n_pos: %d lifetime: %d\n", n_pos, lifetime );
						if( n_pos > pos_thresh * lifetime)
						{
							av_lifetime_in_pos[from] += lifetime;
							nav_lifetime_in_pos[from] += 1;
						}
					}
					t_mark = t;
					pshell = shell_total_filtered[x*max_time+t];
					n_pos = 0;
				}


				if( t < max_time && pos_total[t*n_index+x] ) 
					n_pos++; 				
			}

//			printf("n_pos at finish: %d\n", n_pos );
		}

		printf("Lipid %s\n", l->name );
		printf("Lifetime moving radially out:");
		for( int s = 0; s <= max_shell_lifetime; s++ )
			printf(" %lf", av_lifetime_out[s]/nav_lifetime_out[s] );
		printf("\n");
		printf("Lifetime moving radially in:");
		for( int s = 0; s <= max_shell_lifetime; s++ )
			printf(" %lf", av_lifetime_in[s]/nav_lifetime_in[s] );
		printf("\n");
		printf("Lifetime at end:");
		for( int s = 0; s <= max_shell_lifetime; s++ )
			printf(" %lf", av_lifetime_end[s]/nav_lifetime_end[s] );
		printf("\n");
		
		if( posFile )
		{
			printf("+ Lipid %s\n", l->name );
			printf("+ Lifetime moving radially out:");
			for( int s = 0; s <= max_shell_lifetime; s++ )
				printf(" %lf", av_lifetime_out_pos[s]/nav_lifetime_out_pos[s] );
			printf("\n");
			printf("+ Lifetime moving radially in:");
			for( int s = 0; s <= max_shell_lifetime; s++ )
				printf(" %lf", av_lifetime_in_pos[s]/nav_lifetime_in_pos[s] );
			printf("\n");
			printf("+ Lifetime at end:");
			for( int s = 0; s <= max_shell_lifetime; s++ )
				printf(" %lf", av_lifetime_end_pos[s]/nav_lifetime_end_pos[s] );
			printf("\n");
		}
	}
/*
	for( int b = 0; b < nbins; b++ )
		printf("%lf %lf %lf\n", time_min + (b+0.5)*(time_max-time_min)/nbins,
				lifetime_PIP2_in[b], lifetime_PIP2_out[b] );	 
	*/

}




