
/*
 *
 *	brendel.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/01/2016 03:25:21 PM CEST  
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include "basic-types.h"
#include "alignment.h"
#include "memory.h"
#include "mathematics.h"
#include "iupac.h"
#include "biofiles.h"

Alignment* splicedaligndp (char *read, unsigned int m, char *genome, unsigned int n, gene_t **model);
mapping_t* bl_dpsplicealign2map(Alignment *al, gene_t *model, MultiCharSeq *mseq, Uint vpos, Uint vlen, char strand, char *querydesc, char *query, char *qual, Uint ulen, char ismate) ;
char bl_checkSpliceAlign(mapping_t *m);
Alignment* splicedaligndpopt (char *read, unsigned int m, char *genome, unsigned int n, gene_t **model, Uint a, Uint b, Uint l, Uint r);

