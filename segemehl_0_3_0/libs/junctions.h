#ifndef JUNCTIONS_H
#define JUNCTIONS_H

/*
 *
 *	junctions.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10/30/16 21:05:10 CET  
 *
 */
#include "segemehl.h"
void bl_getSpliceJunctionsFromMappingSet(mappingset_t *set, MultiCharSeq *mseq, char *readname, segemehl_t *nfo);


#endif
