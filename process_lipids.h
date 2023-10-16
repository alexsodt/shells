#ifndef __process_lipidsh__
#define __process_lipidsh__
int isAltSterol( const char * res );
int nChains( char *resname );
int isPI( char *resname );
int isChl( char *resname );
int isSphingo( char *resname );
int isPE( char *resname );
int isPS( char *resname );
int isSaturated( char *resname );
int whichChain (char *resName, char *atomName );
int lipidType( const char * resname );
#endif
