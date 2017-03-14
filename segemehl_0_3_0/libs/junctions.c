
/*
 *  junctions.c
 *  routines to report and store split read junctions
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 10/30/16 19:57:19 CET
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "biofiles.h"
#include "mapfrag.h"
#include "mathematics.h"
#include "segemehl.h"


void
bl_getSpliceJunctionsFromMappingSet(mappingset_t *set, MultiCharSeq *mseq, char *readname, segemehl_t *nfo) {

  Uint i,k;

  for(i=0; i < set->n; i++) {
    for(k=0; k < 2; k++) {

     
    //  bl_dumpMappingSet(nfo->multisplitdev, set);
      bl_dumpSpliceJunctions(&set->elem[i], k, mseq, nfo->readgroupid, readname, nfo);
    
       }
  }

  return;
}


