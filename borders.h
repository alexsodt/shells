#ifndef __bordersh__
#define __bordersh__
#define MAX_BORDERS 25 

struct border_header
{
	int nencoded;
	int origResID[MAX_ENCODED];
	int origType[MAX_ENCODED];
};

struct border_data
{
	int nborders;
	int bpartners[MAX_BORDERS];
	float bvals[MAX_BORDERS];
	float phi; // location relative to an alignment of the protein
	int is_interacting_with_positive;
	int which_chain;
};
#endif
