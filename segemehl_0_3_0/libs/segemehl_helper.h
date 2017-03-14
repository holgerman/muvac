#ifndef SEGEMEHL_HELPER_H
#define SEGEMEHL_HELPER_H

/*
 *
 *	segemehl_helper.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 12/01/16 23:55:07 CET  
 *
 */
#include "segemehl.h"


typedef struct seseq_s {
  
  char *sequence;
  char *sequence_rc;
  char *quality;
  char *quality_re;
  char *convert;
  char *convert_rc;
  Uint len;
  char converttype;
  char converttype_rc;

} seseq_t;

void
se_getData(void *space, seseq_t *seq, char **seqs, char **quals, char bisulfite, char type);

void
se_segemehlSeqInit(void *space, seseq_t *seq, char *orig, char *qual, Uint len); 

void
getqualandseq(void *space, char *seq, char *qual, char **seqs, char **quals, Uint len);

void
convert(void *space, char **seqs, char *orig, Uint len, Uint phase, Uint type, segemehl_t *nfo);

void
se_segemehlSeqDestruct(void *space, seseq_t *seq);

#endif
