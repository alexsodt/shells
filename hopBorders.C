// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "dcd.h"
#include "ctype.h"
#include "alignSet.h"
#include "border_config.h"
#include "process_lipids.h"
#include "borders.h"

#define SUB_PRESS
#define LEN 100

#define MAX_CHAINS 10
//#define COLOR_RHODOPSIN

static double lowp_cut = 6.0;
static double highp_cut = 18.0;

double POS_CUT = 3.5;

//#define OUTPUT_PS
#define REMOVE_COM

int N_TALLY = 6; 

int nChains( char *resname );
int whichChain (char *resName, char *atomName );
int isProtRes( const char *resname )
{
	if( !strcasecmp( resname, "DALA" ) ) return 1;
	if( !strcasecmp( resname, "DLEU" ) ) return 1;
	if( !strcasecmp( resname, "DVAL" ) ) return 1;

	if( !strcasecmp( resname, "TUBE") ) return 1;

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


int nhelices = 8;
int helix_start[8] = { 25, 70, 107, 148, 199, 248, 284, 177  };
int helix_stop[8] =  { 65, 102, 138, 177, 230, 279, 311, 198 };
double rhod_colors[8][3] =
{
	{ 0.0, 0.0, 1.0 },
	{ 0.0, 0.0, 0.85 },
	{ 0.0, 0.0, 0.72 },
	{ 0.0, 0.0, 0.59 },
	{ 0.0, 0.0, 0.46 },
	{ 0.0, 0.0, 0.33 },
	{ 0.0, 0.0, 0.2 },
	{ 0.5, 0.5, 0.5 }
/*	{ 0.0, 0.0, 1.0 },
	{ 1.0, 0.0, 1.0 },
	{ 0.5, 0.5, 1.0 },
	{ 1.0, 1.0, 0.0 },
	{ 0.4, 0.0, 0.0},
	{ 0.25, 0.25, 1.0 },
	{ 1.0, 0.5, 0.5 },
	{ 0.2, 0.5, 0.2 }*/
};
double writeFrame( double *atoms_in, int nat_in, double *colors, int seq, const char *unique, double Lx, double Ly, char *fileName, int *atCode, int doPS, int *borders, int *nborders, double *borderLens, int *solvationShell, double *area, int *local_shell=NULL  );
char code( int cnts[3] )
{
	static int done = 0;
	static int *array=  NULL;//[(N_TALLY+1)*(N_TALLY+1)];

	if( done == 0 )
	{
		array = (int *)malloc( sizeof(int) * (N_TALLY+1)*(N_TALLY+1) );

		char cur = 'A';

		for( int x = 0; x <= N_TALLY; x++ )
		for( int y = 0; y <= N_TALLY; y++ )
			array[x*(N_TALLY+1)+y] = 'a';

		for( int Neither = 0; Neither <= N_TALLY; Neither++ )
		for( int y = 0; y <= Neither; y++ )
		{
			int x = Neither - y;

			array[y*(N_TALLY+1)+x] = cur;

			if( cur == 'Z' )
				cur = 'a';
			else
				cur += 1;
		}

		done = 1;
	}	

	return array[cnts[0]*(N_TALLY+1)+cnts[1]];
}

struct elem
{
	int i, j,k;
	double PXX,PYY,PZZ;
};

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
//	int whichChain;
//};


int main( int argc_in, char **argv_in )
{
	char buffer[4096];
	FILE *alignPDB=NULL;
	int use_ref_align=0;
	int doCOM = 1;//!strncasecmp( argv[2], "COM", 3 );
	int doChains = 0;
	int multiProt = 0;
	int doPeripheral = 0;
	int nfrcRes = 0;
	int nfrcResSpace = 10;
	int *frc_res = (int *)malloc( sizeof(int) * nfrcResSpace );


	struct protein_align
	{
		char segid[256];
		double pcen[3];
		int npcen;
		double dphi;
		int nat;	
		int *align_at;
		int *align_at_ref;
		int nalign_use;
		struct protein_align *next;
	};
	
	protein_align *to_align = NULL;

	char *argv[argc_in];
	int argc = 0;

	for( int c = 0; c < argc_in; c++ )
	{
		if( !strncasecmp( argv_in[c], "--", 2) )
		{
			if( !strncasecmp( argv_in[c], "--align=",8) )
			{
				alignPDB = fopen(argv_in[c]+8,"r");
				use_ref_align = 1;
				printf("atan2(-1,0) + M_PI/2: %lf\n", atan2(-1,0) + M_PI/2 );
				if( !alignPDB )
				{
					printf("Couldn't open alignment reference pdb file '%s'.\n", argv_in[c]+8 );
					exit(1);
				}
			}
			else if( !strncasecmp( argv_in[c], "--frc=", 6 ) )
			{
				if( nfrcRes == nfrcResSpace )
				{
					nfrcResSpace *= 2;
					frc_res = (int*)realloc( frc_res, sizeof(int) * nfrcResSpace );
				}
				frc_res[nfrcRes] = atoi(argv_in[c]+6);
				nfrcRes++;
			}
			else if( !strcasecmp( argv_in[c], "--multiProt") )
				multiProt = 1;
			else if( !strcasecmp( argv_in[c], "--peripheral") )
				doPeripheral = 1;
			else if( !strcasecmp( argv_in[c], "--chain") )
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

	if( argc < 3 )
	{
		printf("Syntax: hopBorders psf dcd [dcd2 ... ]\n");	
		printf("Options:\n");
		printf("\t--chain ; Run in chain mode vs. whole-lipid.\n");
		printf("\t--align=XX.pdb ; report polar locations of lipids relative to reference structure\n");
		return 0;
	}

	int nat_align = 0;
	struct atom_rec *ref_align = NULL;
	if( use_ref_align )
	{
		loadPSFfromPDB( alignPDB );
		nat_align = curNAtoms();	
		ref_align = (struct atom_rec *)malloc( sizeof(struct atom_rec) * nat_align );
		rewind(alignPDB);
		loadPDB( alignPDB, ref_align, nat_align );
	}


	FILE *psfFile = fopen(argv[1], "r" );
	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[1] );
		return 0;
	}

	if( !strcasecmp( argv[1] + strlen(argv[1])-3, "pdb" ) ) 
		loadPSFfromPDB( psfFile );    
        else
		loadPSF( psfFile );
	fclose(psfFile);	
	
#define N_BORDER_BINS   100
#define BORDER_MAX     25.0

	double border_histogram16[N_BORDER_BINS];
	double border_histogram24[N_BORDER_BINS];

	memset( border_histogram16, 0, sizeof(double) * N_BORDER_BINS );	
	memset( border_histogram24, 0, sizeof(double) * N_BORDER_BINS );	
	

	double ll_border_histogram16[N_BORDER_BINS];
	double ll_border_histogram24[N_BORDER_BINS];

	memset( ll_border_histogram16, 0, sizeof(double) * N_BORDER_BINS );	
	memset( ll_border_histogram24, 0, sizeof(double) * N_BORDER_BINS );	

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );


	int np_link[100];
	int n_np = 0;

	int init_done = 0;

	int nat = curNAtoms();
	int type[nat];
	

	double *last_pos = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	memset( last_pos, 0, sizeof(double) * curNAtoms() * 3 );
	
	double *phi = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	double *cur_pos_local = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	memset( cur_pos_local, 0, sizeof(double ) * 3 * curNAtoms() );
	double *cur_pos = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	memset( cur_pos, 0, sizeof(double) * curNAtoms() * 3 );
	
	
	int *is_p = (int *)malloc( sizeof(int) * curNAtoms() );

	int did_init = 0;

	int bl[nat];

	int *res_id;
	int n_encoded = 0;
	int n_encoded_extended = 0;
	int *at_link = NULL;
	int *which_chain = NULL;
	int *at_code = NULL;
	int *at_link_top_bottom = NULL;


	int *frame_leaflet = NULL;
	double *prev_pos = NULL;

	double *colors = NULL;
	double *xy;
	double *saved_areas;
	double *xy_com;
	double *nxy_com;
	int nframes_done = 0;

	int frame = 0;
	int the_id = 0;
	int id_tb = 0;	

	int gf = 0;
	
	int n_real = 0;
			
	int is_set = 0;		
	int *old_borders = NULL;
	int *old_nborders = NULL;

	int nalign_use = 0;
	int *align_at_ref = NULL;
	int *align_at     = NULL;

	int *atom_to_lipid_res = (int *)malloc( sizeof(int) * curNAtoms() );
	int *lipid_res_to_site = (int *)malloc( sizeof(int) * curNAtoms() );	
	int *lipid_res_to_site_n = (int *)malloc( sizeof(int) * curNAtoms() );	
	int *pick_leaflet = (int *)malloc( sizeof(int) * curNAtoms() );
	for( int a = 0; a < curNAtoms(); a++ )
	{
		lipid_res_to_site[a] = -1;
		pick_leaflet[a] = -1;
	}

	FILE *theLog = fopen("hop.log","w");

	FILE *borderBin; 

	borderBin = fopen( "borders.bin", "wb");
	
	double *cur_set = (double *)malloc( sizeof(double) * 3 * nat );
	double *last_align = (double *)malloc( sizeof(double) * 3 * nat );
	int aligned = 0;

	int frame_freq = 1;

	// list of positive residues on the protein
	int *pos_list = (int *)malloc( sizeof(int) * nat );
	int *lipid_start = (int *)malloc( sizeof(int) * nat );
	int *lipid_stop  = (int *)malloc( sizeof(int) * nat );
	int npos = 0;

	for( int i = 0; i < nat; i++ )
	{
		lipid_start[i] = -1;
		lipid_stop[i] = -1;
	}

	for( int c = 2; c < argc; c++ )
	{
		FILE *dcdFile = fopen(argv[c], "r");

		if( ! dcdFile )
		{
			printf("Couldn't open dcd file '%s'.\n", argv[c] );
			return 0;
		}
			
		readDCDHeader(dcdFile);
		setSymmetric();
		int nframes = curNFrames();

		for( int f = 0; f < nframes; f++, gf++ )
		{
			fprintf(theLog, "Computing frame %d.\n", f );
			fflush(theLog);
			double cur_com[2] = {0,0};
	
			double La, Lb, Lc;
			double alpha,beta,gamma;
			
			loadFrame( dcdFile, at );
			if( !DCDsuccess() )
			{
				nframes = f;
				break;
			}
		
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );
			for( int x = 0; x < nat; x++ )
			{
				cur_set[x*3+0] = at[x].x;
				cur_set[x*3+1] = at[x].y;
				cur_set[x*3+2] = at[x].z;
			}
	
			if( aligned )
			{
				for( int p = 0; p < nat; p++ )
				{
					while( cur_set[p*3+0] - last_align[p*3+0] < -La/2 ) cur_set[3*p+0] += La; 
					while( cur_set[p*3+0] - last_align[p*3+0] >  La/2 ) cur_set[3*p+0] -= La; 
					while( cur_set[p*3+1] - last_align[p*3+1] < -Lb/2 ) cur_set[3*p+1] += Lb; 
					while( cur_set[p*3+1] - last_align[p*3+1] >  Lb/2 ) cur_set[3*p+1] -= Lb; 
					while( cur_set[p*3+2] - last_align[p*3+2] < -Lc/2 ) cur_set[3*p+2] += Lc; 
					while( cur_set[p*3+2] - last_align[p*3+2] >  Lc/2 ) cur_set[3*p+2] -= Lc; 
				}
			}

			for( int i = 0; i < nat; i++ )
			{
				at[i].x = cur_set[3*i+0];
				at[i].y = cur_set[3*i+1];
				at[i].z = cur_set[3*i+2];
			}

			memcpy( last_align, cur_set, sizeof(double) * 3 * nat );
			aligned = 1;

			double pcen[3] = { 0,0,0};
			int tnp =0;
			
			if( !init_done )
			{
				if( use_ref_align ) {

					if( multiProt )
					{
						char pseg[256];
						strcpy( pseg, at[0].segid );
						int init = 1;
						// first, find all the proteins to do:
						for( int a2 = 0; a2 < nat; a2++ )
						{
							if( (init || strcasecmp( at[a2].segid, pseg )) && !strncasecmp( at[a2].segid, "PRO", 3) )
							{
								protein_align * an_align = (protein_align *)malloc( sizeof(protein_align) );
								an_align->nalign_use = 0;
								an_align->align_at = (int *)malloc( sizeof(int) * nat_align );
								an_align->align_at_ref = (int *)malloc( sizeof(int) * nat_align );
								strcpy( an_align->segid, at[a2].segid );
								an_align->next = to_align;	
								to_align = an_align;
							}
							init=0;
							strcpy( pseg, at[a2].segid);
						}

						for( protein_align *an_align = to_align; an_align; an_align = an_align->next )
						{
							for( int a2 = 0; a2 < nat; a2++ )
							{
								if( strcasecmp( at[a2].segid, an_align->segid ) ) continue;

								for( int a1 = 0; a1 < nat_align; a1++ )
								{
									if( at[a2].res < ref_align[a1].res ) break;
		
									if( ref_align[a1].res == at[a2].res && !strcasecmp( ref_align[a1].atname, at[a2].atname) )
									{
										an_align->align_at_ref[an_align->nalign_use] = a1;
										an_align->align_at[an_align->nalign_use] = a2;
										an_align->nalign_use++;
										break;
									}
								}
							}
						} 
					}		
					else
					{
						align_at_ref = (int *)malloc( sizeof(int) * nat_align );
						align_at = (int *)malloc( sizeof(int) * nat_align );
						nalign_use = 0;

						for( int a1 = 0; a1 < nat_align; a1++ )
						{
							for( int a2 = 0; a2 < nat; a2++ )
							{
								if( strcasecmp( ref_align[a1].segid, at[a2].segid) ) continue;
	
								if( at[a2].res > ref_align[a1].res ) break;
	
								if( ref_align[a1].res == at[a2].res && !strcasecmp( ref_align[a1].atname, at[a2].atname) )
								{
									align_at_ref[nalign_use] = a1;
									align_at[nalign_use] = a2;
									nalign_use++;
									break;
								}
							}
						}
					}
				}
			}

			if( multiProt )
			{
				for( protein_align *an_align = to_align; an_align; an_align = an_align->next )
				{
					an_align->pcen[0] = 0;
					an_align->pcen[1] = 0;
					an_align->pcen[2] = 0;
					an_align->npcen = 0;		
				}

				for( int a = 0; a < nat; a++ )
				{
					if( at[a].segid && strlen(at[a].segid) > 1 && !strcasecmp( at[a].segid, "PROA") )
					{
						pcen[0] += at[a].x;
						pcen[1] += at[a].y;
						pcen[2] += at[a].z;
						tnp += 1;
					}

					if( strncasecmp( at[a].segid, "PRO", 3 ) ) continue;
		
					for( protein_align *an_align = to_align; an_align; an_align = an_align->next )
					{
						if( !strcasecmp( an_align->segid, at[a].segid ) )
						{
							an_align->pcen[0] += at[a].x;
							an_align->pcen[1] += at[a].y;
							an_align->pcen[2] += at[a].z;
							an_align->npcen+=1;
						}
					}
				}
				
				for( protein_align *an_align = to_align; an_align; an_align = an_align->next )
				{
					an_align->pcen[0] /= an_align->npcen;
					an_align->pcen[1] /= an_align->npcen;
					an_align->pcen[2] /= an_align->npcen;
				}
			}
			else
			{
					int first_wrap=-1;
				for( int a = 0; a < nat; a++ )
				{
					if( (at[a].segid && strlen(at[a].segid) > 1 && (!strncasecmp( at[a].segid, "PRO", 3) || !strcasecmp( at[a].atname, "BB"))) || isProtRes(at[a].resname) )
					{
						if( first_wrap == -1 ) 
							first_wrap = a;

						if( strcasecmp( at[a].segid, at[first_wrap].segid ) )
						{
							double dr[3] = { at[a].x - pcen[0]/tnp,
									 at[a].y - pcen[1]/tnp,
									 at[a].z - pcen[2]/tnp };

							while( dr[0]  < -La/2 ) dr[0] += La;
							while( dr[1]  < -La/2 ) dr[1] += Lb;
							while( dr[2]  < -La/2 ) dr[2] += Lc;
							while( dr[0]  > La/2 ) dr[0] -= La;
							while( dr[1]  > Lb/2 ) dr[1] -= Lb;
							while( dr[2]  > Lc/2 ) dr[2] -= Lc;

							pcen[0] += pcen[0]/tnp + dr[0];
							pcen[1] += pcen[1]/tnp + dr[1];
							pcen[2] += pcen[2]/tnp + dr[2];
							tnp+=1;
						}	
						else
						{
							pcen[0] += at[a].x;
							pcen[1] += at[a].y;
							pcen[2] += at[a].z;
							tnp += 1;
						}
					}
				}
			}

			if( tnp > 0 )
			{
				pcen[0] /= tnp;
				pcen[1] /= tnp;
				pcen[2] /= tnp;

				printf("pcen: %lf %lf %lf\n", pcen[0], pcen[1], pcen[2] );
			}

			// center in z around the protein.
			double align_com[3] = { pcen[0], pcen[1], pcen[2] };
			double cur_align[3] = { 0,0,0};
			char pseg[256] = {'\0'};
			int pres = -1;
			int pactive =0;
			for( int a = 0; a < curNAtoms(); a++ )
			{

				// if there's no segid info we have to ignore it.
				int match = 0;
				if( at[a].segid && strlen(at[a].segid) > 1 )	
				{
					if( !strcasecmp( at[a].segid, pseg) )
						match=1;
				}
				else
					match=1;

//				if( !strcasecmp( at[a].segid, pseg) || ((at[a].res != pres) && !( (!strncasecmp(at[a].segid, "PRO",3) || !strcasecmp( at[a].atname, "BB"))  && pactive))  ) 

				if( (pactive && (isProtRes(at[a].resname) || match )) || ( at[a].res != pres || !strcasecmp(at[a].resname, "TUBE") ) )
				{ // always align protein to pcen, lipids only if it's a new residue.
					cur_align[0] = align_com[0];
					cur_align[1] = align_com[1];
					cur_align[2] = align_com[2];
				}
	
				while( at[a].x - cur_align[0] < -La/2 ) at[a].x += La;
				while( at[a].y - cur_align[1] < -Lb/2 ) at[a].y += Lb;
				while( at[a].z - cur_align[2] < -Lc/2 ) at[a].z += Lc;
				while( at[a].x - cur_align[0] >  La/2 ) at[a].x -= La;
				while( at[a].y - cur_align[1] >  Lb/2 ) at[a].y -= Lb;
				while( at[a].z - cur_align[2] >  Lc/2 ) at[a].z -= Lc;
	
				if( (at[a].segid && !strncasecmp( at[a].segid, "PRO",3)) || !strcasecmp( at[a].atname, "BB") || isProtRes(at[a].resname) ) 
					pactive = 1;
				else
					pactive = 0;

				cur_align[0] = at[a].x;
				cur_align[1] = at[a].y;
				cur_align[2] = at[a].z;
				pres = at[a].res;
				if( at[a].segid )
					strcpy( pseg, at[a].segid );
			}

			int nl = 0;
			double lcom[3] = { 0,0,0};
			for( int a = 0; a < nat; a++ )
			{
				if( lipidType( at[a].resname ) != -1 )
				{	
					lcom[0] += at[a].x;
					lcom[1] += at[a].y;
					lcom[2] += at[a].z;
					nl++;
				}
			}

			lcom[0] /= nl;
			lcom[1] /= nl;
			lcom[2] /= nl;

			printf("lcom: %lf %lf %lf\n", lcom[0], lcom[1], lcom[2] );

			for( int a = 0; a < nat; a++ )
			{
				at[a].x -= lcom[0];
				at[a].y -= lcom[1];
				at[a].z -= lcom[2];
			}
	
			if( multiProt )
			{
				for( protein_align *an_align = to_align; an_align; an_align = an_align->next )
				{
					an_align->pcen[0] = 0;
					an_align->pcen[1] = 0;
					an_align->pcen[2] = 0;
					an_align->npcen = 0;		
				}

				for( int a = 0; a < nat; a++ )
				{
					if( strncasecmp( at[a].segid, "PRO", 3 ) ) continue;
		
					for( protein_align *an_align = to_align; an_align; an_align = an_align->next )
					{
						if( !strcasecmp( an_align->segid, at[a].segid ) )
						{
							an_align->pcen[0] += at[a].x;
							an_align->pcen[1] += at[a].y;
							an_align->pcen[2] += at[a].z;
							an_align->npcen+=1;
						}
					}
				}
				
				for( protein_align *an_align = to_align; an_align; an_align = an_align->next )
				{
					an_align->pcen[0] /= an_align->npcen;
					an_align->pcen[1] /= an_align->npcen;
					an_align->pcen[2] /= an_align->npcen;
				}
			}
	
			if( !init_done )
			{
				int pres = -1;
				char pseg[256] = { ' ' };

				for( int a = 0; a < curNAtoms(); a++ )
				{
					if( ( !strcasecmp( at[a].resname, "ARG" ) || !strcasecmp( at[a].resname, "LYS" )) && at[a].atname[0] == 'N' ) 
					{
						pos_list[npos] = a;
						npos++;
					}

					if( !strcasecmp( at[a].resname, "TUBE") ) 
					{
					}
					else
					{
					if( at[a].res == pres && (!at[a].segid || strlen(at[a].segid) == 1 || !strcasecmp( at[a].segid, pseg )) ) continue;
					}

					pres = at[a].res;
					strcpy( pseg, at[a].segid );

					// get one unique atom name
					if( lipidType(at[a].resname) != -1 )
					{
						if(doChains && nChains(at[a].resname) > 1)
							printf("Encoding %s %s [%d chains]\n", at[a].atname, at[a].resname, nChains(at[a].resname) );
						else if( doChains )
							printf("Encoding %s %s [%d chain]\n", at[a].atname, at[a].resname, nChains(at[a].resname) );
						else
							printf("Encoding %s %s\n", at[a].atname, at[a].resname );
						is_p[n_encoded] = 0;
	
						if( doChains )
							n_encoded += nChains( at[a].resname );
						else
							n_encoded += 1;
					}
	

					double dz = at[a].z;
		
					while( dz > Lc /2 ) dz -= Lc;
					while( dz < -Lc /2 ) dz += Lc;
					

					int frc_encode = 0;

					if( ((isProtRes(at[a].resname) || (at[a].segid && !strncasecmp( at[a].segid, "PRO", 3))) ) )
					{
						for( int x = 0; x < nfrcRes; x++ )
							if( at[a].res == frc_res[x] )
								frc_encode = 1;
					}
	
					
					if( frc_encode || ((isProtRes(at[a].resname) || (at[a].segid && !strncasecmp( at[a].segid, "PRO", 3)) || !strcasecmp( at[a].atname, "BB") )  && fabs(dz) > lowp_cut && fabs(dz) < highp_cut ))
					{
						printf("Encoding %s %s %s\n", at[a].segid, at[a].atname, at[a].resname );
						is_p[n_encoded] = 1;
						n_encoded += 1;
					}
				}

				printf("n_encoded: %d\n", n_encoded );
	
	
				frame_leaflet = (int *)malloc (sizeof(int) * n_encoded );
				at_link = (int *)malloc( sizeof( int ) * n_encoded );
				which_chain = (int *)malloc( sizeof( int ) * n_encoded );
				at_code = (int *)malloc( sizeof( int ) * n_encoded );
				at_link_top_bottom = (int *)malloc( sizeof(int) * n_encoded );
				prev_pos = (double *)malloc( sizeof(double) * 3 * n_encoded );
				colors = (double *)malloc( sizeof(double) * 3 * n_encoded );
				res_id = (int*)malloc( sizeof(int) * n_encoded );
				old_nborders = (int *)malloc( sizeof(int) * n_encoded );
				old_borders = (int *)malloc( sizeof(int) * n_encoded * n_encoded );
	
				int id = 0;
				pres = -1;
				int nchain_loop = 1;
				int cur_lipid[MAX_CHAINS];

				for( int c = 0; c < MAX_CHAINS; c++ )
					cur_lipid[c] = -1;

				for( int a = 0; a < curNAtoms(); a++ )
				{
					if( !strcasecmp( at[a].resname, "TUBE") ) 
					{
					}
					else
					{
						if( at[a].res == pres && (!at[a].segid || strlen(at[a].segid) == 1 || !strcasecmp( at[a].segid, pseg ) ) ) {
							if( cur_lipid[0] >= 0 && at[a].res == at[lipid_start[cur_lipid[0]]].res && !strcasecmp( at[a].resname, at[lipid_start[cur_lipid[0]]].resname) )
							{
								for( int c = 0; c < nchain_loop; c++ )
									lipid_stop[cur_lipid[c]] = a;
							}
							continue;
						}
					}

					if( lipidType( at[a].resname ) != -1 )	
					{

						nchain_loop = 1;

						if( doChains )
							nchain_loop = nChains( at[a].resname);
						
						lipid_res_to_site[at[a].res-1] = the_id;
						lipid_res_to_site_n[at[a].res-1] = nchain_loop;

						for( int c = 0; c < nchain_loop; c++ )
						{
							lipid_start[the_id] = a;
							lipid_stop[the_id] = a;
	
							at_link_top_bottom[id_tb] = a;
							id_tb += 1;

							cur_lipid[c] = the_id;
							which_chain[the_id] = c;	
							at_link[the_id] = a;
							at_code[the_id] = 1 + lipidType( at[a].resname );	
							switch( at_code[the_id] ) 
							{
								case 1:
									colors[3*the_id+0] = 1;
									colors[3*the_id+1] = 0;
									colors[3*the_id+2] = 0;
									break;
								case 2:
									colors[3*the_id+0] = 1;
									colors[3*the_id+1] = 1;
									colors[3*the_id+2] = 0;
									break;
								case 3:
									colors[3*the_id+0] = 1;
									colors[3*the_id+1] = 0;
									colors[3*the_id+2] = 1;
									break;
								default:
									colors[3*the_id+0] = 0.5;
									colors[3*the_id+1] = 0.5;
									colors[3*the_id+2] = 0.5;
									break;
							}

							if( c == 1 )
							{
								colors[3*the_id+0] = colors[3*the_id+0] *0.5 + 0.5;
								colors[3*the_id+1] = colors[3*the_id+1] *0.5 + 0.5;
								colors[3*the_id+2] = colors[3*the_id+2] *0.5 + 0.5;
							}

							res_id[the_id] = at[a].res;
	
				//			lipid_res_to_site[at[a].res-1] = the_id;
							the_id++;
						}
					}
					else
						cur_lipid[0] = -1;

					int frc_encode = 0;

					if( ((isProtRes(at[a].resname) || (at[a].segid && !strncasecmp( at[a].segid, "PRO", 3))) ) )
					{
						for( int x = 0; x < nfrcRes; x++ )
							if( at[a].res == frc_res[x] )
								frc_encode = 1;
					}
					
					if( frc_encode || (( isProtRes(at[a].resname) || (at[a].segid && !strncasecmp( at[a].segid, "PRO", 3)) || !strcasecmp( at[a].atname, "BB"))  && fabs( at[a].z) > lowp_cut && fabs(at[a].z) < highp_cut) )
					{

						at_link[the_id] = a;
						at_code[the_id] = 0;	
						which_chain[the_id] = 0;

#ifdef COLOR_RHODOPSIN
						colors[3*the_id+0] = 0.25;
						colors[3*the_id+1] = 0.25;
						colors[3*the_id+2] = 0.25;

						int gotit = 0;
						for( int h = 0; h < nhelices; h++ )
						{
							if( at[a].res >= helix_start[h] && at[a].res <= helix_stop[h] )
							{
								colors[3*the_id+0] = rhod_colors[h][0];
								colors[3*the_id+1] = rhod_colors[h][1];
								colors[3*the_id+2] = rhod_colors[h][2];
								gotit = 1;
							}
						}
						if( !gotit )
						{
							fprintf(theLog, "Couldn't assign color to %s %s %s %d\n", at[a].segid, at[a].resname, at[a].atname, at[a].res ); 
						}
#else
						colors[3*the_id+0] = 0.0;
						colors[3*the_id+1] = 0.0;
						colors[3*the_id+2] = 1.0;
#endif

						res_id[the_id] = at[a].res;
						the_id++;
					}

					pres = at[a].res;
					if( at[a].segid )
						strcpy(pseg,at[a].segid);
				}

				for( int a = 0; a < curNAtoms(); a++ )
				{
					atom_to_lipid_res[a] = -1;
					
					if( lipidType( at[a].resname ) != -1 )  
						atom_to_lipid_res[a] = at[a].res-1;
				}

				n_real = the_id;
				border_header theHeader;
				theHeader.nencoded = n_real;

				if( n_real > MAX_ENCODED )
				{
					printf("n_real %d MAX_ENCODED %d\n", n_real, MAX_ENCODED );
					exit(1);
				}

				for( int i = 0; i < n_real; i++ )
				{
					theHeader.origType[i] = at_code[i];
					theHeader.origResID[i] = at[at_link[i]].res;

					printf("HEADER %d %d\n", theHeader.origType[i], theHeader.origResID[i] );
				}
				fwrite( &theHeader, sizeof(border_header), 1, borderBin ); 
			}
		
			n_encoded_extended = n_encoded;
			n_encoded = n_real;

			if( !init_done )
			{
				xy = (double *)malloc( sizeof(double) * 3 * n_encoded );
				saved_areas = (double *)malloc( sizeof(double) * n_encoded );
				xy_com = (double *)malloc( sizeof(double) * 3 * n_encoded  );
				nxy_com = (double *)malloc( sizeof(double) * 3 * n_encoded);
			}
				
			init_done = 1;
						
			double avxy[2] = { 0, 0};
			int navxy = 0;

			for( int e = 0; e < n_encoded; e++ )
			{
				if(at_code[e] != 0 ) continue;

				if( navxy == 0 )
				{
					avxy[0] = at[at_link[e]].x;	
					avxy[1] = at[at_link[e]].y;	

					navxy += 1;
				}
				else
				{
					double ccen[2] = { avxy[0] / navxy, avxy[1] / navxy };
					double toadd[2] = { at[at_link[e]].x, at[at_link[e]].y };
					while( toadd[0] - ccen[0] < -La/2 ) toadd[0] += La;
					while( toadd[0] - ccen[0] >  La/2 ) toadd[0] -= La;
					while( toadd[1] - ccen[1] < -Lb/2 ) toadd[1] += Lb;
					while( toadd[1] - ccen[1] >  Lb/2 ) toadd[1] -= Lb;
	
					avxy[0] += toadd[0];
					avxy[1] += toadd[1];
					
					navxy += 1;
				}
			}

			if( navxy > 0 )
			{
				avxy[0] /= navxy;
				avxy[1] /= navxy;
			}

			for( int x = 0; x < curNAtoms(); x++ )
			{
				at[x].x -= avxy[0];
				at[x].y -= avxy[1];
	
				while( at[x].x < -La/2 ) at[x].x += La;
				while( at[x].y < -Lb/2 ) at[x].y += Lb;
				while( at[x].x > La/2 ) at[x].x -= La;
				while( at[x].y > Lb/2 ) at[x].y -= Lb;
			}
			
			for( int e = 0; e < n_encoded; e++ )
			{
				xy_com[3*e+0] = 0;
				xy_com[3*e+1] = 0;
				xy_com[3*e+2] = 0;
				nxy_com[e] = 0;
			}

			for( int a = 0; a < curNAtoms(); a++ )
			{
				if( atom_to_lipid_res[a] >= 0 && lipid_res_to_site[atom_to_lipid_res[a]] >= 0 )
				{
					int nc = lipid_res_to_site_n[atom_to_lipid_res[a]];

					if( pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]] == -1 )
					{
						for( int i = 0; i < nc; i++ )
							pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]+i] = a;
					}
					if( at[a].atname[0] == 'P' ) 
					{
						for( int i = 0; i < nc; i++ )
							pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]+i] = a;
					}
					if( !strcasecmp( at[a].resname, "CHL1" ) && at[a].atname[0] == 'O' )
					{
						for( int i = 0; i < nc; i++ )
							pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]+i] = a;
					}
					if( !strcasecmp( at[a].resname, "CHOL" ) && !strcasecmp( at[a].atname, "ROH") )
					{
						for( int i = 0; i < nc; i++ )
							pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]+i] = a;
					}

					for( int c = 0; c < lipid_res_to_site_n[atom_to_lipid_res[a]]; c++ )
					{
						if( whichChain( at[a].resname, at[a].atname ) != c ) continue;	

						int site = lipid_res_to_site[atom_to_lipid_res[a]]+c;
						int n = nxy_com[site];
						if(nxy_com[site] > 0 )
						{
							double cur[3] = { xy_com[3*(site)+0] / n,
									  xy_com[3*(site)+1] / n,
									  xy_com[3*(site)+2] / n };
	
							double toadd[3] = { at[a].x, at[a].y, at[a].z };
	
							while( toadd[0] - cur[0] < -La/2 ) toadd[0] += La;
							while( toadd[0] - cur[0] >  La/2 ) toadd[0] -= La;
							while( toadd[1] - cur[1] < -Lb/2 ) toadd[1] += Lb;
							while( toadd[1] - cur[1] >  Lb/2 ) toadd[1] -= Lb;
	
		
							xy_com[3*site+0] += toadd[0];					
							xy_com[3*site+1] += toadd[1];					
							xy_com[3*site+2] += toadd[2];					
						}
						else
						{
							xy_com[3*site+0] += at[a].x;
							xy_com[3*site+1] += at[a].y;
							xy_com[3*site+2] += at[a].z;
						}
						nxy_com[site] += 1;
					}
				}
			}


			double avz = 0;

			for( int e = 0; e < id_tb; e++ )
				avz += at[at_link_top_bottom[e]].z;

			avz /= id_tb;

			for( int e = 0; e < n_encoded; e++ )
			{
				xy[e*3+0] = at[at_link[e]].x;	
				xy[e*3+1] = at[at_link[e]].y;	
				xy[e*3+2] = at[at_link[e]].z;	
	
				if( doCOM && at_code[e] != 0 )
				{
					cur_pos[3*e+0] = La/2 + xy_com[e*3+0] / nxy_com[e];	
					cur_pos[3*e+1] = Lb/2 + xy_com[e*3+1] / nxy_com[e];	
					cur_pos[3*e+2] = xy_com[e*3+2] / nxy_com[e] - avz;	
				}
				else
				{
					cur_pos[3*e+0] = at[at_link[e]].x + La/2;	
					cur_pos[3*e+1] = at[at_link[e]].y + Lb/2;	
					cur_pos[3*e+2] = at[at_link[e]].z - avz;	
				}

				while( cur_pos[3*e+0] < -La/2 ) cur_pos[3*e+0] += La;
				while( cur_pos[3*e+1] < -Lb/2 ) cur_pos[3*e+1] += Lb;
				while( cur_pos[3*e+0] > La/2 ) cur_pos[3*e+0] -= La;
				while( cur_pos[3*e+1] > Lb/2 ) cur_pos[3*e+1] -= Lb;
			}
			
			double use_pcen[3] = {0,0,0};
			int nuse_p=0;
			for( int e = 0; e < n_encoded; e++ )
			{
				if( is_p[e] )
				{
					if( nuse_p == 0 )
					{
						use_pcen[0] = cur_pos[3*e+0];
						use_pcen[1] = cur_pos[3*e+1];
						use_pcen[2] = cur_pos[3*e+2];
					}
					else
					{
						double dr[3] = { cur_pos[3*e+0] - use_pcen[0]/nuse_p,
								 cur_pos[3*e+1] - use_pcen[1]/nuse_p,
								 cur_pos[3*e+2] - use_pcen[2]/nuse_p };
						while( dr[0] < -La/2 ) dr[0] += La;
						while( dr[1] < -Lb/2 ) dr[1] += Lb;
						while( dr[0] >  La/2 ) dr[0] -= La;
						while( dr[1] >  Lb/2 ) dr[1] -= Lb;

						use_pcen[0] += use_pcen[0]/nuse_p + dr[0];
						use_pcen[1] += use_pcen[1]/nuse_p + dr[1];
						use_pcen[2] += use_pcen[2]/nuse_p + dr[2];
					}
					nuse_p++;
				}
			}

			use_pcen[0]/=nuse_p;
			use_pcen[1]/=nuse_p;
			use_pcen[2]/=nuse_p;

			double base_del_phi = 0;

			if( multiProt )
			{
				if( use_ref_align )
				{	
					for( protein_align *an_align = to_align; an_align; an_align = an_align->next )
					{
						// estimate the rotation of the current system from the reference.
						double com1[3] = { 0,0,0};
						double com2[3] = { 0,0,0};
		
						for( int ax = 0; ax < an_align->nalign_use; ax++ )
						{
							com1[0] += at[an_align->align_at[ax]].x;
							com1[1] += at[an_align->align_at[ax]].y;
							com1[2] += at[an_align->align_at[ax]].z;
							
							com2[0] += ref_align[an_align->align_at_ref[ax]].x;
							com2[1] += ref_align[an_align->align_at_ref[ax]].y;
							com2[2] += ref_align[an_align->align_at_ref[ax]].z;
						}
		
						com1[0] /= an_align->nalign_use;
						com1[1] /= an_align->nalign_use;
						com1[2] /= an_align->nalign_use;
						
						com2[0] /= an_align->nalign_use;
						com2[1] /= an_align->nalign_use;
						com2[2] /= an_align->nalign_use;
						double del_phi = 0;		
						for( int ax = 0; ax < an_align->nalign_use; ax++ )
						{
							double dr1[2] = { (ref_align[an_align->align_at_ref[ax]].x - com2[0]),
									 (ref_align[an_align->align_at_ref[ax]].y - com2[1]) };
							double dr2[2] = {  (at[an_align->align_at[ax]].x - com1[0]),
									  (at[an_align->align_at[ax]].y - com1[1]) };
							double dphi = atan2( -dr2[1], dr2[0] ) - atan2( -dr1[1], dr1[0] );
							while( dphi > M_PI ) dphi -= 2*M_PI;
							while( dphi < -M_PI ) dphi += 2*M_PI;
		
							if( ax == 0 )
								del_phi = dphi;
							else
							{
								while( dphi - del_phi/ax  < -M_PI ) dphi += 2*M_PI;
								while( dphi - del_phi/ax  >  M_PI ) dphi -= 2*M_PI;
		
								del_phi += dphi;
							}
						}
		
						an_align->dphi = del_phi / an_align->nalign_use;
					}	
				}
	
				for( int e = 0; e < n_encoded; e++ )
				{
					
					protein_align *nearp = NULL;
					double near_r2 = 1e10;

					for( protein_align * an_align = to_align; an_align; an_align = an_align->next )
					{
						double dr[2] = { cur_pos[3*e+0] - (an_align->pcen[0]+La/2),
								 cur_pos[3*e+1] - (an_align->pcen[1]+Lb/2) };
	
						if( doPeripheral && ( cur_pos[3*e+2] < 0 && an_align->pcen[2] > 0 ) )
							continue;
						if( doPeripheral && ( cur_pos[3*e+2] > 0 && an_align->pcen[2] < 0 ) )
							continue;

						while( dr[0] >  La/2 ) dr[0] -= La;
						while( dr[1] >  Lb/2 ) dr[1] -= Lb;
						while( dr[0] < -La/2 ) dr[0] += La;
						while( dr[1] < -Lb/2 ) dr[1] += Lb;

						double r2 = dr[0]*dr[0]+dr[1]*dr[1];
						if( r2 < near_r2 )
						{
							near_r2 = r2;
							nearp = an_align;
						}
						
					}

					if( nearp )
					{
						double dr[2] = { cur_pos[3*e+0] - (nearp->pcen[0]+La/2),
								 cur_pos[3*e+1] - (nearp->pcen[1]+Lb/2) };
		
						while( dr[0] >  La/2 ) dr[0] -= La;
						while( dr[1] >  Lb/2 ) dr[1] -= Lb;
						while( dr[0] < -La/2 ) dr[0] += La;
						while( dr[1] < -Lb/2 ) dr[1] += Lb;
		
		
						phi[e] = -nearp->dphi + atan2( -dr[1], dr[0] ) + M_PI / 2;
		
						while( phi[e] < 0 ) phi[e] += 2*M_PI;
						while( phi[e] > 2*M_PI ) phi[e] -= 2*M_PI;

						//printf("nearp: %lf dphi: %lf\n", sqrt(dr[0]*dr[0]+dr[1]*dr[1]), phi[e] );
					}
					else
						phi[e] = 0;
				}
			}
			else
			{
				if( use_ref_align )
				{
	
					// estimate the rotation of the current system from the reference.
					double com1[3] = { 0,0,0};
					double com2[3] = { 0,0,0};
	
					for( int ax = 0; ax < nalign_use; ax++ )
					{
						com1[0] += at[align_at[ax]].x;
						com1[1] += at[align_at[ax]].y;
						com1[2] += at[align_at[ax]].z;
						
						com2[0] += ref_align[align_at_ref[ax]].x;
						com2[1] += ref_align[align_at_ref[ax]].y;
						com2[2] += ref_align[align_at_ref[ax]].z;
					}
	
					com1[0] /= nalign_use;
					com1[1] /= nalign_use;
					com1[2] /= nalign_use;
					
					com2[0] /= nalign_use;
					com2[1] /= nalign_use;
					com2[2] /= nalign_use;
					
					for( int ax = 0; ax < nalign_use; ax++ )
					{
						double dr1[2] = { (ref_align[align_at_ref[ax]].x - com2[0]),
								 (ref_align[align_at_ref[ax]].y - com2[1]) };
						double dr2[2] = {  (at[align_at[ax]].x - com1[0]),
								  (at[align_at[ax]].y - com1[1]) };
						double dphi = atan2( -dr2[1], dr2[0] ) - atan2( -dr1[1], dr1[0] );
						while( dphi > M_PI ) dphi -= 2*M_PI;
						while( dphi < -M_PI ) dphi += 2*M_PI;
	
						if( ax == 0 )
							base_del_phi = dphi;
						else
						{
							while( dphi - base_del_phi/ax  < -M_PI ) dphi += 2*M_PI;
							while( dphi - base_del_phi/ax  >  M_PI ) dphi -= 2*M_PI;
	
							base_del_phi += dphi;
						}
					}
	
					base_del_phi /= nalign_use;
		
					//printf("base_del_phi: %lf\n", base_del_phi );
				}
	
				for( int e = 0; e < n_encoded; e++ )
				{
					double dr[2] = { cur_pos[3*e+0] - use_pcen[0],
							 cur_pos[3*e+1] - use_pcen[1] };
	
					while( dr[0] >  La/2 ) dr[0] -= La;
					while( dr[1] >  Lb/2 ) dr[1] -= Lb;
					while( dr[0] < -La/2 ) dr[0] += La;
					while( dr[1] < -Lb/2 ) dr[1] += Lb;
	
	
					phi[e] = -base_del_phi + atan2( -dr[1], dr[0] ) + M_PI / 2;
	
					while( phi[e] < 0 ) phi[e] += 2*M_PI;
					while( phi[e] > 2*M_PI ) phi[e] -= 2*M_PI;
				}
			}


			border_data *theData = (border_data *)malloc( sizeof(border_data) * n_encoded );
				
			int *borders = (int *)malloc( sizeof(int) * n_encoded * n_encoded );
			int *nborders = (int *)malloc( sizeof(int) * n_encoded );
			double *borderLens = (double *)malloc( sizeof(double) * n_encoded * n_encoded );	
			int *solvation_shell = (int *)malloc( sizeof(int) * n_encoded );
			double *area = (double *)malloc( sizeof(double) * n_encoded );	
			for( int leaflet = 0; leaflet < 2; leaflet++ )
			{		
		

				int n_local = 0;

				int link[n_encoded];
				int at_code_local[n_encoded];
				double saved_areas_local[n_encoded];
				double local_colors[3*n_encoded];
				for( int i = 0; i < n_encoded; i++ )
				{
					double zpick = cur_pos[3*i+2];
					if( at_code[i] > 0 )
						zpick = at[pick_leaflet[i]].z;


					if( !( cur_pos[3*i+0] < 1 || cur_pos[3*i+0] > 0 ) )
					{
						//printf("Check.\n");
					} 
					if( zpick > 0 && leaflet == 0 )
					{
						local_colors[3*n_local+0] = colors[3*i+0];
						local_colors[3*n_local+1] = colors[3*i+1];
						local_colors[3*n_local+2] = colors[3*i+2];
						cur_pos_local[3*n_local+0] = cur_pos[3*i+0];
						cur_pos_local[3*n_local+1] = cur_pos[3*i+1];
						cur_pos_local[3*n_local+2] = cur_pos[3*i+2];
						link[n_local] = i;
						if( i == 18 + 37 || i == 18 + 38 )
						{
//							printf("%lf %lf %lf + Encoding global %d as %d %s:%d.\n", cur_pos[3*i+0], cur_pos[3*i+1], cur_pos[3*i+2], i, n_local, at[at_link[i]].resname, at[at_link[i]].res );
						}
						at_code_local[n_local] = at_code[i];
						n_local++;
					}
					else if( zpick <= 0 && leaflet == 1 )
					{	
						local_colors[3*n_local+0] = colors[3*i+0];
						local_colors[3*n_local+1] = colors[3*i+1];
						local_colors[3*n_local+2] = colors[3*i+2];
						cur_pos_local[3*n_local+0] = cur_pos[3*i+0];
						cur_pos_local[3*n_local+1] = cur_pos[3*i+1];
						cur_pos_local[3*n_local+2] = cur_pos[3*i+2];
						link[n_local] = i;
						if( i == 18 + 37 || i == 18 + 38 )
						{
//							printf(" - %lf %lf %lf Encoding global %d as %d %s:%d.\n", cur_pos[3*i+0], cur_pos[3*i+1], cur_pos[3*i+2], i, n_local, at[at_link[i]].resname, at[at_link[i]].res  );
						}
						at_code_local[n_local] = at_code[i];
						n_local++;
					}


				}
				fprintf(theLog, "n_encoded: %d leaflet %d n_local %d\n", n_encoded, leaflet, n_local );

				int doPS = 0;

				char fileName[256];
				sprintf(fileName, "junk.ps");
				if( 1 && f % frame_freq == 0 )
				{
					char name5[256];
					print5(f/frame_freq, name5 );
					sprintf(fileName, "frame%s.ps", name5  );
					if( leaflet == 0 ) doPS = 1;
					if( leaflet == 1 ) doPS = 0;
				}


				for( int i = 0; i < n_local; i++ )
				{
					if( is_p[link[i]] )
					{
						//printf("leaflet: %d seg: %s protein z: %lf\n", leaflet, at[at_link[link[i]]].segid, cur_pos_local[3*i+2] );
					}
				}

				writeFrame( cur_pos_local, n_local, local_colors, f, "hi", La, Lb, fileName, at_code_local, doPS, borders, nborders, borderLens, NULL, saved_areas_local );// + f * n_encoded); 

				for( int i = 0; i < n_local; i++ )
					saved_areas[link[i]] = saved_areas[i];

				for( int i = 0; i < n_local; i++ )
				{
					theData[link[i]].nborders = nborders[i];
					for( int x = 0; x < nborders[i]; x++ )
					{
						theData[link[i]].bpartners[x] = link[borders[i*n_local+x]];
						theData[link[i]].bvals[x] = borderLens[i*n_local+x];
					}
				}
			}
			
			for( int i = 0; i < n_encoded; i++ )
			{
				theData[i].phi = phi[i];
				theData[i].is_interacting_with_positive = 0;
				
				if( lipid_start[i] >= 0 && lipid_stop[i] >= 0 )
				{
					for( int a = lipid_start[i]; a <= lipid_stop[i]; a++ )
					{
						if( at[a].atname[0] != 'O') continue;

						for( int ip = 0; ip < npos; ip++ )
						{
							double dr[3] = { at[a].x - at[pos_list[ip]].x,
									 at[a].y - at[pos_list[ip]].y,
									 at[a].z - at[pos_list[ip]].z };
							while( dr[0] < -La/2) dr[0] += La;
							while( dr[1] < -Lb/2) dr[1] += Lb;
							while( dr[2] < -Lc/2) dr[2] += Lc;
							while( dr[0] >  La/2) dr[0] -= La;
							while( dr[1] >  Lb/2) dr[1] -= Lb;
							while( dr[2] >  Lc/2) dr[2] -= Lc;

							double l = sqrt( dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2] );
							if( l < POS_CUT )
								theData[i].is_interacting_with_positive = 1;
						}
					}
							
				}	
				else
				{
//					printf("Uh oh.\n");
				}
			}
	
			free(borderLens);
			free(nborders);
			free(borders);
			free(solvation_shell);
			free(area);

			fwrite( theData, sizeof(border_data), n_encoded, borderBin );

			free(theData);

			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
		}
		fclose(dcdFile); 
		nframes_done += nframes;
	}
	
	fclose(theLog);
	
	free(xy);
	free(saved_areas);
	free(xy_com);
	free(nxy_com);
				



	free(frame_leaflet);
	free(at_link);
	free(at_code);
	free(at_link_top_bottom);
	free(prev_pos);
	free(colors);
	free(res_id);
	free(old_nborders);
	free(old_borders);

	free(at);
	



	free(atom_to_lipid_res);
	free(lipid_res_to_site);
	free(pick_leaflet);
	free(cur_set);
	free(last_align);
	free(pos_list);
	free(lipid_start);
	free(lipid_stop);

	
	free(last_pos);
	free(cur_pos_local);
	free(cur_pos);
	free(is_p);
	
}


