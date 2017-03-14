#ifndef MATCH_H
#define MATCH_H

/*
 *
 *	match.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/28/14 12:01:39 CEST  
 *
 */

#include "segemehl.h"
void
match(void *space, Suffixarray *s, fasta_t *reads, 
    segemehl_t *nfo);


#endif
