
/*
 *
 *	queryalign.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 29.07.2013 17:57:49 CEST  
 *
 */

#ifndef QUERYALIGN_H
#define QUERYALIGN_H

#include "mapfrag.h"
#include "multicharseq.h"

 mapseedlist_t*
bl_getGoodSeeds (matchstem_t **items, unsigned m, unsigned n, karlin_t *stats, 
    segemehl_t *nfo);

mappingset_t*
bl_seedAlign(Suffixarray *s, mappingset_t* set, MultiCharSeq *mseq, 
    char **seqs, char **qual, Uint len, char* qname, mapseedlist_t *l, 
    segemehl_t *nfo, Uint *enctab, bitvector* D, char ismate);


#endif

