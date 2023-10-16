#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

#define MAX_SHELLS 20

int main ( int argc_in, char ** argv_in )
{
	int argc = 0;
	char *argv[argc_in];
	int print_traj = 0;
	FILE *leafletInput = NULL;
	FILE *phiInput = NULL;
	int n_error_bins = 10;
	int do_phi_hist = 0;
	int nphi = 12;
	int nshell_average = 2;
	int start =0;
	int stop = -1;
	int do_both = 1;
	int do_upper = 0;
	int do_lower = 0;
	int phi_traj_bin = -1;
	int do_counts = 0;
	int do_pts = 0;

	for( int c = 0; c < argc_in; c++ )
	{
		if( !strncasecmp( argv_in[c], "--", 2 ) )
		{
			if( !strcasecmp( argv_in[c], "--traj") ) 
				print_traj=1;
			else if( !strncasecmp( argv_in[c], "--nerr=",7) ) 
				n_error_bins = atoi(argv_in[c]+7);
			else if( !strncasecmp( argv_in[c], "--nphi=",7) ) 
				nphi = atoi(argv_in[c]+7);
			else if( !strncasecmp( argv_in[c], "--nshell=",9) ) 
				nshell_average = atoi(argv_in[c]+9);
			else if( !strncasecmp( argv_in[c], "--phitraj=",10) ) 
				phi_traj_bin = atoi(argv_in[c]+10);
			else if( !strncasecmp( argv_in[c], "--phihist", 9 ) )
			{
				do_phi_hist = 1;
			}
			else if( !strncasecmp( argv_in[c], "--phiinput=",11) ) 
				phiInput = fopen(argv_in[c]+11,"r");
			else if( !strncasecmp( argv_in[c], "--leaflet=",10) ) 
				leafletInput = fopen(argv_in[c]+10,"r");
			else if( !strncasecmp( argv_in[c], "--upper",7) ) 
			{	do_upper = 1; do_both = 0; }
			else if( !strncasecmp( argv_in[c], "--lower",7) ) 
			{	do_lower = 1; do_both = 0; }
			else if( !strncasecmp( argv_in[c], "--start=",8) ) 
				start = atoi(argv_in[c]+8);
			else if( !strncasecmp( argv_in[c], "--stop=",7) ) 
				stop = atoi(argv_in[c]+7);
			else if( !strncasecmp( argv_in[c], "--counts",8) ) 
				do_counts = 1;
			else if( !strncasecmp( argv_in[c], "--pts",5) ) 
				do_pts = 1;
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

	if( argc < 3 )
	{
		printf("Syntax histo index.file shells.file [interacting_with_pos.txt|null]\n");
		return -1;
	}

	FILE *index = fopen(argv[1],"r");
	FILE *shells = fopen(argv[2],"r");
	FILE *posFile = NULL;

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
		double *hist_err;
		double *inst_hist;
		double *pos_hist;
		double *phi_tot;
		double *phi_tot_err;
		FILE *ptsFile;
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
			got->ptsFile = NULL;
			got->indices = (int *)malloc( sizeof(int) * got->nindSpace );
			got->hist = (double *)malloc( sizeof(double ) * MAX_SHELLS );
			got->pos_hist = (double *)malloc( sizeof(double ) * MAX_SHELLS );
			got->inst_hist = (double *)malloc( sizeof(double ) * MAX_SHELLS );
			got->phi_tot = (double *)malloc( sizeof(double ) * nphi );
			memset( got->phi_tot, 0, sizeof(double) * nphi );
			if( n_error_bins > 0 )
			{
				got->hist_err = (double *)malloc( sizeof(double ) * MAX_SHELLS * n_error_bins );
				memset( got->hist_err, 0, sizeof(double) * MAX_SHELLS * n_error_bins);
				got->phi_tot_err = (double *)malloc( sizeof(double ) * nphi * n_error_bins );
				memset( got->phi_tot_err, 0, sizeof(double) * nphi * n_error_bins);
			}
			got->next = lipids;
			lipids = got;
			memset( got->hist, 0, sizeof(double) * MAX_SHELLS );
			memset( got->pos_hist, 0, sizeof(double) * MAX_SHELLS );
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
	
	int *leaflet = (int *)malloc( sizeof(int) * (n_index+1) );

	memset( leaflet, 0, sizeof(int) * n_index );

	if( leafletInput )
	{
		getLine( leafletInput, buffer );

		int nr = readNInts( buffer, leaflet, n_index );

		if( nr != n_index )
		{
			printf("Failed to read exactly %d integers (indicating leaflet) from file.\n");
			exit(1);
		}
	}

	double *phi_read = (double *)malloc( sizeof(double) * n_index );
	int *shell_read = (int *)malloc( sizeof(int) * n_index );
	int *pos_read = (int *)malloc( sizeof(int) * n_index );
	memset( pos_read, 0, sizeof(int) * n_index );

	int fr = 0;
	int nvals = 0;
	if( n_error_bins > 0 )
	{
		int fp = ftell( shells );

		while( !feof(shells) ) 
		{
			getLine( shells, buffer);
			if( feof(shells) ) break;
	
			int nr = readNInts( buffer, shell_read, n_index );
			if( nr == n_index )
				nvals++;
		}		
		fseek(shells, fp, SEEK_SET );
	}

	if( stop < nvals && stop >= 0 )
		nvals = stop;
 	int use_start = ( start < 0 ? 0 : start );

	int curp=0;
	while( !feof(shells) ) 
	{
		getLine( shells, buffer);
		if( feof(shells) ) break;

		int err_bin = 0;

//		if( nvals > 0 && do_phi_hist )
		{
			err_bin = curp * (double)n_error_bins / (double)(nvals-start);
		}
		int nr = readNInts( buffer, shell_read, n_index );
		
		if( nr != n_index )
		{
			printf("Failed to read %d shells (read %d).\n", n_index, nr );
			exit(1);
		}
	
		if( phiInput )
		{
			getLine( phiInput, buffer );
			nr = readNDoubles( buffer, phi_read, n_index );
		
			if( nr != n_index )
			{
				printf("Failed to read %d phi values (read %d).\n", n_index, nr );
				exit(1);			
			}
		}
	
		for( struct ltype *l = lipids; l; l = l->next )
			memset( l->inst_hist, 0, sizeof(double) * MAX_SHELLS );

		if( posFile )
		{		
			getLine( posFile, buffer);
			if( feof(posFile) ) break;
			int p_nr = readNInts( buffer, pos_read, n_index );
			if( p_nr != n_index )
			{	
				printf("Failed to read %d positive-res indicators (read %d).\n", n_index, p_nr );
			}
		}

		double ntot_ptraj=0;
		double nocc_ptraj=0;

		if( fr >= start && (fr < stop || stop < 0) )
		{
			for( struct ltype *l = lipids; l; l = l->next )
			{
				for( int i = 0; i < l->nind; i++ )
				{
					int leaflet_match =0;
					if( do_both || (do_upper && leaflet[l->indices[i]] == 1) || (do_lower && leaflet[l->indices[i]] == -1) )
						leaflet_match = 1;
					if( !leaflet_match ) continue;

					int s = shell_read[l->indices[i]];
					if( s < MAX_SHELLS )
					{
						l->hist[s] += 1;
						l->inst_hist[s] += 1;
						if( n_error_bins > 0 )
							l->hist_err[err_bin*MAX_SHELLS+s] += 1;

						if( do_phi_hist )
						{	
							int phi_bin = nphi * phi_read[l->indices[i]] / (2*M_PI);
							if( s <= nshell_average )
							{ 
								if(phi_bin == phi_traj_bin )
								{
									ntot_ptraj += 1;
									if( l == lipids )
										nocc_ptraj+=1;
								}

								l->phi_tot[phi_bin] += 1;
								if( n_error_bins > 0 )
									l->phi_tot_err[err_bin*nphi+phi_bin] += 1;
							}
						}

					}
					if( pos_read[l->indices[i]] )
						l->pos_hist[s] += 1;
				} 
			}
			curp++;
		}

		if( phi_traj_bin >= 0 )
			printf("%d %lf\n", fr, (double)nocc_ptraj/(double)ntot_ptraj );
		
		if( print_traj )
		{
			printf("%d", fr );
			for( struct ltype *l = lipids; l; l = l->next )
				printf(" %lf", l->inst_hist[1] );
			printf("\n");
		}

		fr++;

	}
	
	if( do_phi_hist )
	{
		printf("#Phi_bin");
		for( struct ltype *l = lipids; l; l = l->next )
			printf(" %s", l->name );
		if( n_error_bins > 0 )
		{
			for( struct ltype *l = lipids; l; l = l->next )
				printf(" %s_err", l->name );
		}
		printf("\n");

		for( int p = 0; p < nphi; p++ )
		{
			printf("%d", p );
			double ntot = 0;
			double ntot_err[n_error_bins];
			memset( ntot_err, 0, sizeof(double) * n_error_bins );

			if( n_error_bins > 0 )
			{
				for( struct ltype *l = lipids; l; l = l->next )
				{
					ntot += l->phi_tot[p];
					for( int b = 0; b < n_error_bins; b++ )
						ntot_err[b] += l->phi_tot_err[b*nphi+p];
				}
			}
			// the occupancy
			for( struct ltype *l = lipids; l; l = l->next )
			{
				double occ = l->phi_tot[p]/ntot;
				printf(" %lf", occ );
			}

			if( n_error_bins > 0 )
			{
				for( struct ltype *l = lipids; l; l = l->next )
				{
					double x=0,x2=0;
					for( int b = 0; b < n_error_bins; b++ )
					{
						double occ = l->phi_tot_err[b*nphi+p]/ntot_err[b];
						x += occ;
						x2 += occ*occ;
					}
					x /= n_error_bins;
					x2 /= n_error_bins;
					printf(" %lf", sqrt(x2-x*x)/sqrt(n_error_bins-1));
				}
			}
			printf("\n");
		}	
	}
	else
	{
		if( do_pts )
		{
			for( struct ltype *l = lipids; l; l = l->next )
			{	
				char fileName[256];
				if( strlen(l->name) < 250 ) 
				{
					sprintf(fileName, "%s.pts", l->name );
					l->ptsFile = fopen(fileName,"w");
				}
			}
		}	

		if( ! print_traj )
		{
		printf("Shell");
		for( struct ltype *l = lipids; l; l = l->next )
			printf(" %s", l->name );
		printf("\n");
	
		for( int s = 0; s < MAX_SHELLS; s++ )
		{	
			double gtot=0;
			for( struct ltype *l = lipids; l; l = l->next )
			{
				gtot += l->hist[s];
				if( l->ptsFile )
					fprintf(l->ptsFile, "%d ", s );
			}
			printf("%d", s );
	
			double tot = 0;
	
			for( struct ltype *l = lipids; l; l = l->next )
				printf(" %lf", (do_counts ? l->hist[s] : l->hist[s] / gtot) );

			if( n_error_bins > 0 && ! do_counts )
			{
				printf(" err ");
				for( struct ltype *l = lipids; l; l = l->next )
				{
					double x = 0, x2 = 0;
					for( int b = 0; b < n_error_bins; b++ )
					{
						double tot = 0;
						for( struct ltype *l2 = lipids; l2; l2 = l2->next )
							tot += l2->hist_err[b*MAX_SHELLS+s];	
						double inst = l->hist_err[b*MAX_SHELLS+s] / tot;
						x2 += inst*inst;
						x += inst;

						if( l->ptsFile )
							fprintf(l->ptsFile, " %lf", inst );
					}
					x2 /= n_error_bins;
					x /= n_error_bins;

					printf(" %lf", sqrt(x2-x*x)/sqrt(n_error_bins-1) );
				}
			}
						
			for( struct ltype *l = lipids; l; l = l->next )
			{
				if( l->ptsFile )
					fprintf(l->ptsFile, "\n" );
			}
			printf("\n");
		}



	
		if( posFile )
		{
			printf("+Shell");
			for( struct ltype *l = lipids; l; l = l->next )
				printf(" %s", l->name );
			printf("\n");
	
			for( int s = 0; s < MAX_SHELLS; s++ )
			{
				printf("%d", s );
		
				double tot = 0;
	
				for( struct ltype *l = lipids; l; l = l->next )
					printf(" %lf", l->pos_hist[s] );
				printf("\n");
			}
		
		}	
		}
	}
}




