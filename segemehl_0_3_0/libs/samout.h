#ifndef SAMOUT_H
#define SAMOUT_H

#include "multicharseq.h"
#include "mapfrag.h"
#include "segemehl.h"
#include "biofiles.h"

/*
 *
 *	samout.h
 *  samout
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 21.03.2012 08:41:55 CET  
 *
 */

typedef struct samtag_s {
  char *tag;
  char *key;
  char *type;
  char *val;
} samtag_t;

typedef struct samrec_s {
  char *qname;
  unsigned flag;
  char *rname;
  uint64_t pos;
  uint8_t mapq;
  char *cigar;
  char *rnext;
  uint64_t pnext;
  int64_t tlen;
  char *seq;
  char *qual;

  Uint nooftags;
  samtag_t *tags;

} samrec_t;

typedef struct samlist_s {
  unsigned int noofrecs;
  samrec_t *recs;

} samlist_t;

typedef struct samheader_s{
  char *version;
  
  char **rnames;
  uint64_t *rlens;
  uint64_t nrnames;

  char **rgroups;
  char **rgroupsinfo;
  Uint nrgroups;
  
  char *cmd;

} samheader_t;
/*
typedef struct mapfragment_s {

  char *id;
  Uint subject;
  unsigned char rc;
  Uint i;
  Uint j;
  Uint p;
  Uint q;
  int scr; 
  int mat;
  int mis;
  int ins;
  int del;
  int edist;
  Alignment *al;
  double evalue;
  char *seq;
  char *qual;
  char *rname;

  char issplit;
  char usenextfieldforsplits;

  Uint fragno;
  Uint previdx;
  Uint prevpos;
  char prevflags;
  Uint nextidx;
  Uint nextpos;
  char nextflags;

  unsigned char skip;

} mapfrag_t;

typedef struct {
  Uint nooffrags;
  mapfrag_t *frags;

} mapfraglist_t;*/


samlist_t *sam_getSamList (fasta_t *reads, Uint id, mapping_t *l, MultiCharSeq *mseq, char ismultiple, Uint noofqueries, Uint noofmates, Uint queryid, Uint mateid, char unmappedtornext, segemehl_t *nfo); 
void sam_printSamrec (FILE *dev, samrec_t* r, char lf);
char* sam_printSamrec2Buffer (samrec_t* r, char lf);
void sam_printSamlist (samlist_t *l, mapping_t *f, segemehl_t*);
void sam_destructSamList (samlist_t *list);
mappingset_t* sam_mappingJoinFrags (mappingset_t *s, segemehl_t *nfo);
void sam_printEmptyAlign (char *desc, char* seq, char* qual, char hasPaired, char isQuery, char nomatemapped, char *nextchr, int64_t nextrpos, char nextrc, char ismultiple, char ischimeric, mapseed_t *seed, segemehl_t *nfo);
samheader_t* sam_getHeader (samheader_t *head, char *line);
samrec_t * sam_line2rec(char *line, Uint len, samheader_t *head);
samheader_t* sam_readFile(char* filename);
    //unsigned char gzip, struct access *index,)
samtag_t* sam_getTag (samrec_t *rec, char* key);
void sam_addReadGroup (samheader_t *head, char *id, char *info);
void sam_initHeader (samheader_t *head);
void sam_dumpHeader (samheader_t  *head);
void sam_destructHeader (samheader_t *head);
void sam_destruct(samrec_t *samrec);

#endif
