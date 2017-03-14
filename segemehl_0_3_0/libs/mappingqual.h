
/*
 *
 *	mappingqual.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 09.11.2015 09:51:14 EST  
 *
 */
#include "karlin.h"

double type3mappingset(mappingset_t *set, char *desc);
double type2mappingset(mappingset_t *set, char *desc, Uint n, karlin_t *stats);
double longestmatchqual(mappingset_t *set, Uint n, karlin_t* stats); 

