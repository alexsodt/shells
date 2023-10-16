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
#define SUB_PRESS
#define LEN 100

//#define ASSIGN_RHODOPSIN
int nhelices = 8;
int helix_start[8] = { 25, 70, 107, 148, 199, 248, 284, 177  };
int helix_stop[8] =  { 65, 102, 138, 177, 230, 279, 311, 198 };
double av_lipid_op_fullh( struct atom_rec *at, int nat, double *chain, double *nchain );

//#define OUTPUT_PS
#define REMOVE_COM

int N_TALLY = 6; 

int lipidType( const char * resname )
{
	if( !strcasecmp( resname, "DPPC") ) return 0;
	if( !strcasecmp( resname, "SOPC") ) return 0;
	if( !strcasecmp( resname, "PPPC") ) return 0;
	if( !strcasecmp( resname, "MMPC") ) return 0;
	if( !strcasecmp( resname, "CHL1") ) return 1;
	if( !strcasecmp( resname, "CHOL") ) return 1;
	if( !strcasecmp( resname, "DOPC") ) return 2;
	if( !strcasecmp( resname, "POPC") ) return 0;
	return -1;
}

#define MAX_BORDER 12

struct LipidBorderData
{
        char type;  // DXPC or DNPC
        char nborders;
  	short int leaflet;
        short int shell;
        short int parentLipid;
        short int border_id[MAX_BORDER]; // who do we border?
        float border_len[MAX_BORDER]; // what length?
	float op;
};


double writeFrame( double *atoms_in, int nat_in, double *colors, int seq, const char *unique, double Lx, double Ly, char *fileName, int *atCode, int doPS, int *borders, int *nborders, double *borderLens, int *solvationShell, double *area, int *local_shell=NULL );
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

int main( int argc, char **argv )
{
	char buffer[4096];

	if( argc < 4 )
	{
		printf("Syntax: rehop shell.file psf dcd [dcd2 ... ]\n");	
		return 0;
	}

	int doCOM=1;
	FILE *shellFile = fopen(argv[1],"r");
	if( ! shellFile )
	{
		printf("Couldn't open shell file '%s'.\n", argv[1] );
		return 0;
	}

	FILE *psfFile = fopen(argv[2], "r" );
	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[2] );
		return 0;
	}

	if( !strcasecmp( argv[2] + strlen(argv[2])-3, "pdb" ) ) 
		loadPSFfromPDB( psfFile );    
        else
		loadPSF( psfFile );
	fclose(psfFile);	
	

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );


	int np_link[100];
	int n_np = 0;

	int init_done = 0;

	int nat = curNAtoms();
	int type[nat];
	

	double *last_pos = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	memset( last_pos, 0, sizeof(double) * curNAtoms() * 3 );
	
	double *cur_pos_local = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	memset( cur_pos_local, 0, sizeof(double ) * 3 * curNAtoms() );
	double *cur_pos = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	memset( cur_pos, 0, sizeof(double) * curNAtoms() * 3 );
	
	
	int *is_p = (int *)malloc( sizeof(int) * curNAtoms() );

	int did_init = 0;

	int bl[nat];

	int *sshell = NULL;
	int *res_id;
	int n_encoded = 0;
	int n_encoded_extended = 0;
	int *at_link = NULL;
	int *at_code = NULL;
	int *at_link_top_bottom = NULL;
	double *prev_pos = NULL;

	double *colors = NULL;
	double *xy;
	double *saved_areas;
	double *xy_com;
	double *nxy_com;
	double *max_border;
	double *sum_border;
	int nframes_done = 0;

	int frame = 0;
	int the_id = 0;
	int id_tb = 0;	

	int gf = 0;
	

	int n_real = 0;
			
	int is_set = 0;		
	int *old_borders = NULL;
	int *old_nborders = NULL;

	int *atom_to_lipid_res = (int *)malloc( sizeof(int) * curNAtoms() );
	int *lipid_res_to_site = (int *)malloc( sizeof(int) * curNAtoms() );	
	int *pick_leaflet = (int *)malloc( sizeof(int) * curNAtoms() );

	FILE *sevenFile = fopen("debug.7dat", "w");

	for( int a = 0; a < curNAtoms(); a++ )
	{
		lipid_res_to_site[a] = -1;
		pick_leaflet[a] = -1;
	}
	FILE *theLog = fopen("hop.log","w");

	double *shells; 
	int nprot = 0;

	double *cur_set = (double *)malloc( sizeof(double) * 3 * nat );
	double *last_align = (double *)malloc( sizeof(double) * 3 * nat );
	int aligned = 0;

	double nchol[2] = { 0,0};
	double nlip[2] = { 0,0};

	double cnt_chl[2][8][4];
	double cnt_lip[2][8][4];

	for( int l = 0; l < 2; l++ )
	for( int h = 0; h < 8; h++ )
	for( int sh = 0; sh < 4; sh++ )
	{
		cnt_chl[l][h][sh] = 0;
		cnt_lip[l][h][sh] = 0;
	}

	double first_shell_upper[3] = {0,0,0};
	double second_shell_upper[3] = { 0,0,0 };
	double first_shell_lower[3] = {0,0,0};
	double second_shell_lower[3] = { 0,0,0 };

	int prev_tp = -1;

	FILE *borderData = fopen("mc_borders.bin","wb");

	int *encode_start = (int*)malloc( sizeof(int) * curNAtoms() );
	int *encode_length = (int *)malloc( sizeof(int) * curNAtoms() );

	for( int c = 3; c < argc; c++ )
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
			
			double pcen[3]={0,0,0};
			int tnp=0;
			for( int a = 0; a < nat; a++ )
			{
				if( !strncasecmp( at[a].segid, "PRO", 3)  || !strcasecmp( at[a].atname, "BB")  )
				{
					pcen[0] += at[a].x;
					pcen[1] += at[a].y;
					pcen[2] += at[a].z;
					tnp += 1;
				}
			}

			pcen[0] /= tnp;
			pcen[1] /= tnp;
			pcen[2] /= tnp;

			// center in z around the protein.
			double align_com[3] = { pcen[0], pcen[1], pcen[2] };
			double cur_align[3] = { 0,0,0};
			char pseg[256];
			int pres = -1;
			int pactive =0;
			for( int a = 0; a < curNAtoms(); a++ )
			{
				if( strcasecmp( at[a].segid, pseg) || ((at[a].res != pres) && ( !( !strncasecmp(at[a].segid, "PRO",3) || !strcasecmp( at[a].atname, "BB") )  && pactive))  ) 
				{
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
	
				if( !strncasecmp( at[a].segid, "PRO",3)  || !strcasecmp( at[a].atname, "BB") ) 
					pactive = 1;
				else
					pactive = 0;

				cur_align[0] = at[a].x;
				cur_align[1] = at[a].y;
				cur_align[2] = at[a].z;
				pres = at[a].res;
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

			for( int a = 0; a < nat; a++ )
			{
				at[a].x -= lcom[0];
				at[a].y -= lcom[1];
				at[a].z -= lcom[2];
			}
			if( !init_done )
			{
				int pres = -1;
				char pseg[256];

				int cur_encode = 0;

				for( int a = 0; a < curNAtoms(); a++ )
				{
					if( at[a].res == pres && !strcasecmp( at[a].segid, pseg ) ) continue;

					if( cur_encode )
					{
						encode_length[n_encoded-1] = a - encode_start[n_encoded-1];
						cur_encode = 0;
					} 

					pres = at[a].res;
					strcpy( pseg, at[a].segid );

					// get one unique atom name
					if( lipidType(at[a].resname) != -1 )
					{
						encode_start[n_encoded] = a;
						encode_length[n_encoded] = 0;
						is_p[n_encoded] = 0;
						cur_encode = 1;
						n_encoded += 1;
					}
	
					if( (!strncasecmp( at[a].segid, "PRO", 3) || !strcasecmp( at[a].atname, "BB")   ) && fabs( at[a].z) > 6 && fabs(at[a].z) < 18 )
					{
						encode_start[n_encoded] = a;
						encode_length[n_encoded] = 0;
						cur_encode = 1;
						
						is_p[n_encoded] = 1;
						n_encoded += 1;
					}
				}
	
				if( cur_encode )
				{
					encode_length[n_encoded-1] = curNAtoms() - encode_start[n_encoded-1];
					cur_encode =0;
				}
	
				at_link = (int *)malloc( sizeof( int ) * n_encoded );
				at_code = (int *)malloc( sizeof( int ) * n_encoded );
				at_link_top_bottom = (int *)malloc( sizeof(int) * n_encoded );
				prev_pos = (double *)malloc( sizeof(double) * 3 * n_encoded );
				colors = (double *)malloc( sizeof(double) * 3 * n_encoded );
				res_id = (int*)malloc( sizeof(int) * n_encoded );
				old_nborders = (int *)malloc( sizeof(int) * n_encoded );
				old_borders = (int *)malloc( sizeof(int) * n_encoded * n_encoded );
	
				int id = 0;
				pres = -1;
				for( int a = 0; a < curNAtoms(); a++ )
				{
					if( at[a].res == pres && !strcasecmp( at[a].segid, pseg ) ) continue;

					if( lipidType( at[a].resname ) != -1 )	
					{
						at_link_top_bottom[id_tb] = a;
						id_tb += 1;

						at_link[the_id] = a;
						at_code[the_id] = 1 + lipidType( at[a].resname );	
						colors[3*the_id+0] = 1;
						colors[3*the_id+1] = 0;
						colors[3*the_id+2] = 0;
						res_id[the_id] = at[a].res;

						lipid_res_to_site[at[a].res-1] = the_id;
						the_id++;
					}
					
					if( (!strncasecmp( at[a].segid, "PRO", 3) || !strcasecmp( at[a].atname, "BB") ) && fabs( at[a].z) > 6 && fabs(at[a].z) < 18 )
					{
						at_link[the_id] = a;
						at_code[the_id] = 0;	
						colors[3*the_id+0] = 0;
						colors[3*the_id+1] = 0;
						colors[3*the_id+2] = 1;
						res_id[the_id] = at[a].res;
						the_id++;
					}

					pres = at[a].res;
					strcpy(pseg,at[a].segid);
				}

				for( int a = 0; a < curNAtoms(); a++ )
				{
					atom_to_lipid_res[a] = -1;
					
					if( lipidType( at[a].resname ) != -1 )  
						atom_to_lipid_res[a] = at[a].res-1;
				}

				n_real = the_id;
				init_done = 1;
			}
		
			n_encoded_extended = n_encoded;
			n_encoded = n_real;

			if( !sshell )
			{
				nframes_done = nframes;
				sshell = (int*)malloc( sizeof(int) * n_encoded * nframes );
				xy = (double *)malloc( sizeof(double) * 3 * n_encoded * nframes );
				saved_areas = (double *)malloc( sizeof(double) * n_encoded * nframes );
				xy_com = (double *)malloc( sizeof(double) * 3 * n_encoded * nframes );
				nxy_com = (double *)malloc( sizeof(double) * 3 * n_encoded * nframes );
				sum_border = (double *)malloc( sizeof(double) * n_encoded * nframes );
				max_border = (double *)malloc( sizeof(double) * n_encoded * nframes );
				shells = (double *)malloc( sizeof(double) * (n_encoded+10) );
			}
			else if( f == 0 )
			{
				nframes_done += nframes;
				sshell = (int*)realloc( sshell, sizeof(int) * n_encoded * nframes_done );
				xy = (double *)realloc( xy, sizeof(double) * 3 * n_encoded * nframes_done );
				saved_areas = (double *)realloc( saved_areas, sizeof(double) * n_encoded * nframes );
				xy_com = (double *)realloc( xy_com, sizeof(double) * 3 * n_encoded * nframes_done );
				nxy_com = (double *)realloc( nxy_com, sizeof(double) * 3 * n_encoded * nframes_done );
				sum_border = (double *)realloc( sum_border, sizeof(double) * n_encoded * nframes );
				max_border = (double *)realloc( max_border, sizeof(double) * n_encoded * nframes );
				shells = (double *)realloc( shells, sizeof(double) * (n_encoded+10) );
			}
			
			getLine( shellFile, buffer );
			int nr = readNDoubles( buffer, shells, n_encoded );
		
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
	
/*				while( at[x].x < -La/2 ) at[x].x += La;
				while( at[x].y < -Lb/2 ) at[x].y += Lb;
				while( at[x].x > La/2 ) at[x].x -= La;
				while( at[x].y > Lb/2 ) at[x].y -= Lb;*/
			}
			
			for( int e = 0; e < n_encoded; e++ )
			{
				xy_com[gf*3*n_encoded+3*e+0] = 0;
				xy_com[gf*3*n_encoded+3*e+1] = 0;
				xy_com[gf*3*n_encoded+3*e+2] = 0;
				nxy_com[gf*n_encoded+e] = 0;
			}

			for( int a = 0; a < curNAtoms(); a++ )
			{
				if( atom_to_lipid_res[a] >= 0 && lipid_res_to_site[atom_to_lipid_res[a]] >= 0 )
				{
					if( pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]] == -1 )
						pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]] = a;
					if( at[a].atname[0] == 'P' ) 
						pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]] = a;
					if( !strcasecmp( at[a].resname, "CHL1" ) && at[a].atname[0] == 'O' )
						pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]] = a;
					if( !strcasecmp( at[a].resname, "CHOL" ) && !strcasecmp( at[a].atname, "ROH") )
						pick_leaflet[lipid_res_to_site[atom_to_lipid_res[a]]] = a;
					int n = nxy_com[gf*n_encoded+lipid_res_to_site[atom_to_lipid_res[a]]];
					if(nxy_com[gf*n_encoded+lipid_res_to_site[atom_to_lipid_res[a]]] > 0 )
					{
						double cur[3] = { xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+0] / n,
								  xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+1] / n,
								  xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+2] / n };

						double toadd[3] = { at[a].x, at[a].y, at[a].z };

						while( toadd[0] - cur[0] < -La/2 ) toadd[0] += La;
						while( toadd[0] - cur[0] >  La/2 ) toadd[0] -= La;
						while( toadd[1] - cur[1] < -Lb/2 ) toadd[1] += Lb;
						while( toadd[1] - cur[1] >  Lb/2 ) toadd[1] -= Lb;

	
						xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+0] += toadd[0];					
						xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+1] += toadd[1];					
						xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+2] += toadd[2];					
					}
					else
					{
						xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+0] += at[a].x;
						xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+1] += at[a].y;
						xy_com[gf*3*n_encoded+3*lipid_res_to_site[atom_to_lipid_res[a]]+2] += at[a].z;
					}
					nxy_com[gf*n_encoded+lipid_res_to_site[atom_to_lipid_res[a]]] += 1;
				}
			}


			double avz = 0;

			for( int e = 0; e < id_tb; e++ )
				avz += at[at_link_top_bottom[e]].z;

			avz /= id_tb;

			for( int e = 0; e < n_encoded; e++ )
			{
				xy[gf*n_encoded*3+e*3+0] = at[at_link[e]].x;	
				xy[gf*n_encoded*3+e*3+1] = at[at_link[e]].y;	
				xy[gf*n_encoded*3+e*3+2] = at[at_link[e]].z;	
	
				if( doCOM && at_code[e] != 0 )
				{
					cur_pos[3*e+0] = La/2 + xy_com[gf*n_encoded*3+e*3+0] / nxy_com[gf*n_encoded+e];	
					cur_pos[3*e+1] = Lb/2 + xy_com[gf*n_encoded*3+e*3+1] / nxy_com[gf*n_encoded+e];	
					cur_pos[3*e+2] = xy_com[gf*n_encoded*3+e*3+2] / nxy_com[gf*n_encoded+e] - avz;	
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
	
	double first_shell_upper[3] = {0,0,0};
	double second_shell_upper[3] = { 0,0,0 };
	double first_shell_lower[3] = {0,0,0};
	double second_shell_lower[3] = { 0,0,0 };

#ifdef ASSIGN_RHODOPSIN
			int nhelices = 8;
			int helix_start[8] = { 25, 70, 107, 148, 199, 248, 284, 177  };
			int helix_stop[8] =  { 65, 102, 138, 177, 230, 279, 311, 198 };

			double hpos[2][8][3];
			double npos[2][8];

			for( int l = 0; l < 2; l++ )
			for( int h = 0; h < 8; h++ )
			{
				npos[l][h] = 0;
				hpos[l][h][0] = 0;
				hpos[l][h][1] = 0;
				hpos[l][h][2] = 0;
			}

			for( int e = 0; e < n_encoded; e++ )
			{
				if( at_code[e] == 0 )
				{
					int th = -1;
					int a = at_link[e]; 
					for( int h = 0; h < 8; h++ )
					{
						if( at[a].res >= helix_start[h] && at[a].res <= helix_stop[h] )
							th = h;
					}
					if( th >= 0 )
					{
						if( at[a].z >= 0 )
						{
							hpos[0][th][0] += cur_pos[3*e+0];
							hpos[0][th][1] += cur_pos[3*e+1];
							hpos[0][th][2] += cur_pos[3*e+2];
							npos[0][th] += 1;
						}	
						else
						{
							hpos[1][th][0] += cur_pos[3*e+0];
							hpos[1][th][1] += cur_pos[3*e+1];
							hpos[1][th][2] += cur_pos[3*e+2];
							npos[1][th] += 1;
						}	
					}
				}
			}

			for( int l = 0; l < 2; l++ )
			for( int h = 0; h < 8; h++ )
			{
				hpos[l][h][0] /= 1e-8 + npos[l][h];
				hpos[l][h][1] /= 1e-8 + npos[l][h];
				hpos[l][h][2] /= 1e-8 + npos[l][h];
			}
#endif
			char fileName[256];
	
			sprintf(fileName, "frame.%d.ps", f );
	
			int *borders = (int *)malloc( sizeof(int) * n_encoded * n_encoded );
			int *nborders = (int *)malloc( sizeof(int) * n_encoded );
			double *borderLens = (double *)malloc( sizeof(double) * n_encoded * n_encoded );	
			int *solvation_shell = (int *)malloc( sizeof(int) * n_encoded );

			int y = 0;
			for( int x = 0; x < n_encoded; x++ )
			{
				if( at_code[x] != 0 )
				{
					solvation_shell[x] = shells[y]; 
					y++;
				} 
			}
			
                        int olink[n_encoded];
                        int tp = 0;
                        for( int e = 0; e < n_encoded; e++ )
                        {
                                olink[e] = -1;
                                if( at_code[e] != 0 )
                                {
                                        olink[e] = tp;
                                        tp++;
                                }       
                        }

			if( prev_tp == -1 )
			{
				fwrite( &tp, sizeof(int), 1, borderData );
				fwrite( &nframes, sizeof(int), 1, borderData );
				prev_tp = tp;	
			}
			else
			{
				if( tp != prev_tp )
				{
					printf("Fatal error.\n");
					exit(1);
				}
			}

			double ind_op[n_encoded];

			for( int e = 0; e < n_encoded; e++ )
			{
			// need: encoded_link_forward	
			// reslen
 
				ind_op[e] = 0;

				if( 
					!strcasecmp( at[encode_start[e]].resname, "SOPC" ) ||
					!strcasecmp( at[encode_start[e]].resname, "POPC" ) ||
					!strcasecmp( at[encode_start[e]].resname, "PPPC" ) ||
					!strcasecmp( at[encode_start[e]].resname, "MMPC" ) ||
					!strcasecmp( at[encode_start[e]].resname, "DOPC" ) ||
					!strcasecmp( at[encode_start[e]].resname, "DPPC" ) )	
				{
					double l_op[64];
					double n_op[64];
					memset( l_op, 0, sizeof(double) * 64 );
					memset( n_op, 0, sizeof(double) * 64 );
					av_lipid_op_fullh( at+encode_start[e], encode_length[e], l_op, n_op );
					double num = 0, den = 0;
					for( int x = 18; x < 40; x++ )
					{
						num += l_op[x] / (n_op[x]+1e-10);
						den += 1;
					}

					fprintf(theLog, "frame %d ind %d op %lf\n", f, e, num/den );
					ind_op[e] = num / den;
				}
			}
	
					fflush(theLog);

			LipidBorderData *data = (LipidBorderData *)malloc( sizeof(LipidBorderData) * tp );

			for( int e = 0; e < tp; e++ )
				data[e].shell = shells[e];

			for( int e = 0; e < n_encoded; e++ )
			{
				if( olink[e] == -1 ) continue;

				data[olink[e]].type = at_code[e];
				data[olink[e]].parentLipid = at[at_link[e]].res;
				data[olink[e]].nborders = 0;
				data[olink[e]].op = ind_op[e];
			}

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

					if( zpick > 0 && leaflet == 0 )
					{
						local_colors[3*n_local+0] = colors[3*i+0];
						local_colors[3*n_local+1] = colors[3*i+1];
						local_colors[3*n_local+2] = colors[3*i+2];
						cur_pos_local[3*n_local+0] = cur_pos[3*i+0];
						cur_pos_local[3*n_local+1] = cur_pos[3*i+1];
						cur_pos_local[3*n_local+2] = cur_pos[3*i+2];
						link[n_local] = i;
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
						at_code_local[n_local] = at_code[i];
						n_local++;
					}
				}
				fprintf(theLog, "n_encoded: %d leaflet %d n_local %d\n", n_encoded, leaflet, n_local );

				int doPS = 0;
				writeFrame( cur_pos_local, n_local, local_colors, f, "hi", La, Lb, fileName, at_code_local, doPS, borders, nborders, borderLens, NULL, saved_areas_local );// + f * n_encoded); 

				for( int i = 0; i < n_local; i++ )
				{
					if( olink[link[i]] == -1 ) continue;
					int d = olink[link[i]];
					int lnb = 0;
					for( int bx = 0; bx < nborders[i]; bx++ )
					{
						if( olink[link[borders[i*n_local+bx]]] == -1 )
							continue;
						data[d].border_id[lnb] = olink[link[borders[i*n_local+bx]]];
						data[d].border_len[lnb] = borderLens[i*n_local+bx];
						lnb++;
					}
					data[d].leaflet = leaflet;
					data[d].nborders = lnb;
				}
			}

			fwrite(data, sizeof(LipidBorderData), tp, borderData );			


/*
			writeFrame( cur_pos, n_encoded, colors, gf, "hi", La, Lb, fileName, at_code, 0, borders, nborders, borderLens, NULL, NULL );

			
			for( int i = 0; i < n_encoded; i++ )
			{
				if( at_code[i] == 0 ) continue;
				int tr = res_id[i];
				int leaflet = (cur_pos[3*i+2] > 0 ? 0 : 1 );

				if( solvation_shell[i] == 1 )
				{
					if( leaflet == 0 )
						first_shell_upper[at_code[i]-1] += 1;
					else if( leaflet == 1 )
						first_shell_lower[at_code[i]-1] += 1;
				}
				else if( solvation_shell[i] == 2 )
				{
					if( leaflet == 0 )
						second_shell_upper[at_code[i]-1] += 1;
					else if( leaflet == 1 )
						second_shell_lower[at_code[i]-1] += 1;
				}
			}

			printf("f %d first shell upper %lf %lf %lf first shell lower %lf %lf %lf second shell upper %lf %lf %lf second shell lower %lf %lf %lf\n",
				f, 
				first_shell_upper[0], first_shell_upper[1], first_shell_upper[2],
				first_shell_lower[0], first_shell_lower[1], first_shell_lower[2],
				second_shell_upper[0], second_shell_upper[1], second_shell_upper[2],
				second_shell_lower[0], second_shell_lower[1], second_shell_lower[2] );
#ifdef ASSIGN_RHODOPSIN
			int l_cnt_chl[2][8][4];
			int l_cnt_lip[2][8][4];

			for( int l = 0; l < 2; l++ )
			for( int h = 0; h < 8; h++ )
			for( int sh = 0; sh < 4; sh++ )
			{
				l_cnt_chl[l][h][sh] = 0;
				l_cnt_lip[l][h][sh] = 0;
			}

			int l_nlip[2]={0,0};
			int l_nchl[2] = {0,0};
			for( int i = 0; i < n_encoded; i++ )
			{
				if( at_code[i] == 0 ) continue;
				int tr = res_id[i];
				int leaflet = (cur_pos[3*i+2] > 0 ? 0 : 1 );

				int near_h = -1;
				double near_r2 = 1e10;

				for( int h = 0; h < 8; h++ )
				{
					if( npos[leaflet][h] < 1 ) continue;

					double dr[3] = { 
						cur_pos[3*i+0] - hpos[leaflet][h][0],
						cur_pos[3*i+1] - hpos[leaflet][h][1],
						cur_pos[3*i+2] - hpos[leaflet][h][2] };
					dr[2] = 0;
					if( dr[0]*dr[0]+dr[1]*dr[1] < near_r2 )
					{
						near_r2 = dr[0]*dr[0]+dr[1]*dr[1];
						near_h = h;
					}
				} 
				if( near_h < 0 ) continue;

				int sshel = solvation_shell[i]-1;
				if( sshel < 0 ) sshel = 0;
				if( sshel > 3 ) sshel = 3;

				if( at_code[i] == 1 )
					cnt_lip[leaflet][near_h][sshel] += 1;
				else if( at_code[i] == 2 )
					cnt_chl[leaflet][near_h][sshel] += 1;
				
				if( at_code[i] == 1 )
					l_cnt_lip[leaflet][near_h][sshel] += 1;
				else if( at_code[i] == 2 )
					l_cnt_chl[leaflet][near_h][sshel] += 1;
		

#if 0
				if( solvation_shell[i] < 3 && at_code[i] == 1 )
				{	nlip[0] += 1; l_nlip[0] += 1; }
				else if( at_code[i] == 1 )
				{	nlip[1] += 1; l_nlip[1] += 1; }
				else if( solvation_shell[i] < 3 && at_code[i] == 2 )
				{	nchol[0] += 1; l_nchl[0] += 1; }
				else if( at_code[i] == 2 )
				{	nchol[1] += 1; l_nchl[1] += 1; }
#endif
			}
	
			for( int leaflet = 0; leaflet < 2; leaflet++ )
			for( int h = 0; h < 8; h++ )
			{
				printf(" h%d_%c", h, (leaflet == 0 ? 'u' : 'l') );
				for( int s = 0; s < 4; s++ )
					printf(" %d %d", l_cnt_lip[leaflet][h][s],
							 l_cnt_chl[leaflet][h][s] );
			}
			printf("\n");
#endif
*/
//			printf(" %d %d %d %d %d\n", f, l_nlip[0], l_nlip[1], l_nchl[0], l_nchl[1] ); 
			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
			free(data);
			free(borderLens);
			free(nborders);
			free(borders);
			free(solvation_shell);
		}
		fclose(dcdFile); 
	}
	printf("nlip %lf %lf nchl %lf %lf\n", nlip[0], nlip[1], nchol[0], nchol[1] );	
	fclose(theLog);
	
}






