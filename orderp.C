#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "alignSet.h"

// takes a hydrogen-less one, adds hydrogens if necessary.

struct fullh {
	struct atom_rec *ats;
	int nats;
	struct fullh *next;
	char *name;
};

double av_lipid_op( struct atom_rec *at, int nat, double *chain, double *nchain )
{
	static fullh *listh;
	struct atom_rec *full = NULL;
	int nat_full;
	static double opvec[16];
	static double nav[16];
	static int init = 1;
	static int cntr = 0;

	if( init )
	{
		memset( opvec, 0, sizeof(double) * 16 );
		memset( nav, 0, sizeof(double) * 16 );
		init = 0;
	}

	double op_sum = 0;
	double nop = 0;

	for( struct fullh *rec = listh; rec; rec = rec->next )
	{
		if( !strncasecmp(rec->name, at[0].resname, strlen(rec->name) ) )
		{
			full = rec->ats;
			nat_full = rec->nats;
			break;
		}
	}

	if( !full )
	{
		char filename[256];
		sprintf(filename, "%s.fullh.pdb", at[0].resname );
		FILE *theFile = fopen(filename, "r");
		if( !theFile )
		{
			printf("Error!!!\n");	
			exit(1);
		}

		char buffer[4096];
		int nat = 0;
		while( !feof(theFile) )
		{
			getLine( theFile, buffer );
			if( feof(theFile) ) break;
			if( !strncasecmp( buffer, "ATOM", 4 ) )
				nat++;
		}
		rewind(theFile);
		struct atom_rec *ats = (struct atom_rec *)malloc( sizeof(atom_rec) * nat );
		nat=0;

		while( !feof(theFile) )
		{
			getLine( theFile, buffer );

			if( feof(theFile) ) break;

			if( !strncasecmp( buffer, "ATOM", 4 ) )
			{
				readATOM( buffer, ats+nat );
				nat++;
			}
		}

		struct fullh * rec = (struct fullh *)malloc( sizeof(struct fullh) );
		rec->ats = ats;
		rec->nats = nat;
		rec->next = listh;
		listh = rec;
		rec->name = (char *)malloc( sizeof(char) * ( 1 + strlen(ats[0].resname) ) );
		strcpy( rec->name, ats[0].resname );

		full = ats;
		nat_full = nat;
	}

#ifdef DEBUG
	for( int a = 0; a < nat; a++ )
		printATOM( stdout,  1+a, 1, at+a );
#endif
	int curn_dbgprint = 1+nat;

	for( int a = 0; a < nat_full; a++ )
	{
		double at_bonded[3];
		char at_bonded_name[256];
		
		double at_off1[3];
		char at_off1_name[256];
		
		double at_off2[3];
		char at_off2_name[256];

		if( full[a].atname[0] == 'H') // a hydrogen, find who it is bonded to.
		{
			double mind = 1e10;

			int mind_list[10];
			int nmind = 0;
			for( int a2 = 0; a2 < nat_full; a2++ )
			{
				if( full[a2].atname[0] == 'H' ) continue;
				if( a2 == a ) continue;

				double dr[3] = { full[a2].x - full[a].x, full[a2].y - full[a].y, full[a2].z - full[a].z };
			
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				
				if( r < 1.3)
				{
					mind_list[nmind] = a2;
					nmind++;
				}		
			}

			if( nmind > 1 )
			{
				// hmm, not sure who this hydrogen is bonded to.
				printf("error 1.\n");
				exit(1);
			}	
			if( nmind == 0 )
			{
				// hmm, not sure who this hydrogen is bonded to.
				printf("error 2.\n");
				exit(1);
			
			}

			if( full[mind_list[0]].atname[0] != 'C' )
				continue;
			
			int a_bond = mind_list[0];

			strcpy( at_bonded_name, full[mind_list[0]].atname );

			at_bonded[0] = full[mind_list[0]].x;
			at_bonded[1] = full[mind_list[0]].y;
			at_bonded[2] = full[mind_list[0]].z;
	
			int nlist[10];
			int nn = 0;

			for( int a2 = 0; a2 < nat_full; a2++ )
			{
				if( full[a2].atname[0] == 'H' ) continue;
				if( a2 == a ) continue;
				if( a2 == a_bond ) continue;

				double dr[3] = { full[a2].x - full[a_bond].x, full[a2].y - full[a_bond].y, full[a2].z - full[a_bond].z };
			
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				
				if( r < 1.8 )
				{
					nlist[nn] = a2;
					nn++;
				}		
			}
		
			if( nn < 2 )
			{
				// this can happen if it's the last carbon of the chain.
				for( int a2 = 0; a2 < nat_full; a2++ )
				{
					if( full[a2].atname[0] == 'H' ) continue;
					if( a2 == nlist[0] ) continue;
					if( a2 == a ) continue;
					if( a2 == a_bond ) continue;
	
					double dr[3] = { full[a2].x - full[nlist[0]].x, full[a2].y - full[nlist[0]].y, full[a2].z - full[nlist[0]].z };
				
					double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
					
					if( r < 1.8 )
					{
						nlist[nn] = a2;
						nn++;
					}		
				}

			}

			double at_off1[3];
			char at_off1_name[256];
		
			double at_off2[3];
			char at_off2_name[256];
		
			strcpy( at_off1_name, full[nlist[0]].atname );
			at_off1[0] = full[nlist[0]].x;
			at_off1[1] = full[nlist[0]].y;
			at_off1[2] = full[nlist[0]].z;
			
			strcpy( at_off2_name, full[nlist[1]].atname );
			at_off2[0] = full[nlist[1]].x;
			at_off2[1] = full[nlist[1]].y;
			at_off2[2] = full[nlist[1]].z;

			// now populate the positions from the target structure.
		
			double s_at_bonded[3];
			double s_at_off1[3];
			double s_at_off2[3];
			int got[3]={0,0,0};

			for( int a2 = 0; a2 < nat; a2++ )
			{
				if( !strcasecmp( at[a2].atname, at_bonded_name ) ) 
				{
					s_at_bonded[0] = at[a2].x;
					s_at_bonded[1] = at[a2].y;
					s_at_bonded[2] = at[a2].z;
					got[0] = 1;
				}	
				if( !strcasecmp( at[a2].atname, at_off1_name ) ) 
				{
					s_at_off1[0] = at[a2].x;
					s_at_off1[1] = at[a2].y;
					s_at_off1[2] = at[a2].z;
					got[1] = 1;
				}	
				if( !strcasecmp( at[a2].atname, at_off2_name ) ) 
				{
					s_at_off2[0] = at[a2].x;
					s_at_off2[1] = at[a2].y;
					s_at_off2[2] = at[a2].z;
					got[2] = 1;
				}	
			}
			
			if( !got[0] || !got[1] || !got[2] )
			{
				printf("Error 4.\n");
				exit(1);
			}
	
			double coords1[12] = { at_bonded[0], at_bonded[1], at_bonded[2],
					       at_off1[0], at_off1[1], at_off1[2],
					       at_off2[0], at_off2[1], at_off2[2],
					       full[a].x, full[a].y, full[a].z };
			double coords2[12] = { s_at_bonded[0], s_at_bonded[1], s_at_bonded[2],
					       s_at_off1[0], s_at_off1[1], s_at_off1[2],
					       s_at_off2[0], s_at_off2[1], s_at_off2[2],
					       s_at_bonded[0], s_at_bonded[1], s_at_bonded[2] };

			int align[3] = {0,1,2};

			alignStructuresOnAtomSet( coords2, align, coords1, align, 3, 4 );	

			// now use coords1.

			double vec[3] = { coords1[9] - coords1[0], coords1[10] - coords1[1], coords1[11] - coords1[2] };
			double tp[3] = { full[a].x, full[a].y, full[a].z };
#ifdef DEBUG
			full[a].x = coords1[9];
			full[a].y = coords1[10];
			full[a].z = coords1[11];
			printATOM( stdout,  curn_dbgprint, 1, full+a );
			curn_dbgprint+=1;
			full[a].x = tp[9];
			full[a].y = tp[10];
			full[a].z = tp[11];
#endif		
			double rv = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
			vec[0] /= rv;
			vec[1] /= rv;
			vec[2] /= rv;

			
		

			
			char *name = full[a].atname;
			name += 1;	
			int cn = atoi( name );
			if( strlen(name) == 3 && (name[2] == 'A' || name[2] == 'B' || name[2] == 'C') ) 
			{
			}
			else
			{
				if( cn <= 18 && cn >= 1 )
				{
					op_sum += 0.5 * (3 * vec[2]*vec[2] - 1);
					nop += 1;		
					chain[cn-1] += 0.5 * (3 * vec[2]*vec[2] - 1); 
					nchain[cn-1] += 1; 
				}				
			}
		}
	}

#ifdef DEBUG
	exit(1);
#endif
	if(cntr % 1000 == 0 )	
	{
/*
		printf("op: ");
		for( int x = 0; x < 16; x++ )
			printf("%lf ", opvec[x] / nav[x] );
		printf("\n");
*/
	}
	cntr += 1;
	return op_sum / nop;	
}


double av_lipid_op_fullh( struct atom_rec *at, int nat, double *chain, double *nchain )
{
	int nat_full;
	static double opvec[16];
	static double nav[16];
	static int init = 1;
	static int cntr = 0;

	if( init )
	{
		memset( opvec, 0, sizeof(double) * 16 );
		memset( nav, 0, sizeof(double) * 16 );
		init = 0;
	}

	double op_sum = 0;
	double nop = 0;

	for( int a = 0; a < nat; a++ )
	{
		if( at[a].atname[0] == 'H') // a hydrogen, find who it is bonded to.
		{
			double mind = 1e10;

			int mind_list[10];
			int nmind = 0;
			for( int a2 = 0; a2 < nat; a2++ )
			{
				if( a2 == a ) continue;
				if( at[a2].atname[0] == 'H' ) continue;

				double dr[3] = { at[a2].x - at[a].x, at[a2].y - at[a].y, at[a2].z - at[a].z };
			
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				
				if( r < 1.3)
				{
					mind_list[nmind] = a2;
					nmind++;
				}		
			}

			if( nmind > 1 )
			{
				// hmm, not sure who this hydrogen is bonded to.
				printf("error 1.\n");
				exit(1);
			}
	
			if( nmind == 0 )
			{
				// hmm, not sure who this hydrogen is bonded to.
				printf("error 2.\n");
				exit(1);
			
			}

			if( at[mind_list[0]].atname[0] != 'C' )
				continue;
		
			double vec[3] = { at[a].x - at[mind_list[0]].x, at[a].y - at[mind_list[0]].y, at[a].z - at[mind_list[0]].z };
	
			
			double rv = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
			
			vec[0] /= rv;
			vec[1] /= rv;
			vec[2] /= rv;

			char *name = at[a].atname;
			name += 1;	
			int cn = atoi( name );

			if( !strcasecmp( at[a].resname, "SDPE" ) )
			{
				int t = 0;
				if( name[strlen(name)-1] == 'R' || name[strlen(name)-1] == 'S' || name[strlen(name)-1] == 'T' )
					t = 0;
				else if( name[strlen(name)-1] == 'X' || name[strlen(name)-1] == 'Y' || name[strlen(name)-1] == 'Z' )
					t = 1;
				else
					continue;
					
				op_sum += 0.5 * (3 * vec[2]*vec[2] - 1);
				nop += 1;		
				chain[t*22+cn-1] += 0.5 * (3 * vec[2]*vec[2] - 1); 
				nchain[t*22+cn-1] += 1; 
			}
			else if( strlen(name) == 3 && (name[2] == 'A' || name[2] == 'B' || name[2] == 'C') ) 
			{
			}
			else
			{
				if( cn <= 18 && cn >= 1 )
				{
					int t = 0;
					if( name[strlen(name)-1] == 'X' || name[strlen(name)-1] == 'Y' || name[strlen(name)-1] == 'Z' )
						t = 1;
					if( name[strlen(name)-1] == 'F' || name[strlen(name)-1] == 'G' || name[strlen(name)-1] == 'H' )
						t = 1;
 
					op_sum += 0.5 * (3 * vec[2]*vec[2] - 1);
					nop += 1;		
					chain[t*18+cn-1] += 0.5 * (3 * vec[2]*vec[2] - 1); 
					nchain[t*18+cn-1] += 1; 
				}				
			}
		}
	}

#ifdef DEBUG
	exit(1);
#endif
	if(cntr % 1000 == 0 )	
	{
/*
		printf("op: ");
		for( int x = 0; x < 16; x++ )
			printf("%lf ", opvec[x] / nav[x] );
		printf("\n");
*/
	}
	cntr += 1;
	return op_sum / nop;	
}

