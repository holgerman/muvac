
/*
 *  samout.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 21.03.2012 08:51:37 CET
 *  
 */
#include <stdarg.h>
#include <string.h>
#include "samout.h"
#include "mapfrag.h"
#include "segemehl.h"
#include "stringutils.h"
#include "biofiles.h"
#include "mathematics.h"
#include "stringutils.h"
#include "info.h"
#include <inttypes.h>
#include "splitalign.h"
#include "filebuffer.h"
#include "debug.h"


/*---------------------------------- inits -----------------------------------
 *    
 * @brief inits
 * @author Steve Hoffmann 
 *   
 */
 
samrec_t* sam_initRec(samrec_t *samrec) {

  samrec->qname = NULL;
  samrec->flag = 0;
  samrec->rname = NULL;
  samrec->pos = 0;
  samrec->mapq = 0;
  samrec->cigar = NULL;
  samrec->rnext = NULL;
  samrec->pnext = 0;
  samrec->tlen = 0;
  samrec->seq = NULL;
  samrec->qual = NULL;
  samrec->nooftags = 0;
  samrec->tags = NULL;

  return samrec;
}

samrec_t* sam_init(samrec_t *samrec, char *qname) {

  sam_initRec(samrec);
  samrec->qname = my_strdup(qname);

  return samrec;
}


/*------------------------------- destructors --------------------------------
 *    
 * @brief destructors
 * @author Steve Hoffmann 
 *   
 */
 

void sam_destructTag(samtag_t *kv) {
//  if(kv->tag) FREEMEMORY(NULL, kv->tag);
  if(kv->tag) FREEMEMORY(NULL, kv->tag); 
  if(kv->val) FREEMEMORY(NULL, kv->val);
  if(kv->type) FREEMEMORY(NULL, kv->type);
  if(kv->key) FREEMEMORY(NULL, kv->key);
  return;
}

void sam_destruct(samrec_t *samrec) {
  Uint i;

  if(samrec->qname) FREEMEMORY(NULL, samrec->qname);
  if(samrec->rname) FREEMEMORY(NULL, samrec->rname);
  if(samrec->cigar) FREEMEMORY(NULL, samrec->cigar);
  if(samrec->rnext) FREEMEMORY(NULL, samrec->rnext);
  if(samrec->seq)   FREEMEMORY(NULL, samrec->seq);
  if(samrec->qual)  FREEMEMORY(NULL, samrec->qual);
  for(i=0; i < samrec->nooftags; i++) {
    sam_destructTag(&samrec->tags[i]);
  }
  if(samrec->tags) FREEMEMORY(NULL, samrec->tags);
  
 
  return;
}


/*--------------------------------- getters ----------------------------------
 *    
 * @brief getter functions
 * @author Steve Hoffmann 
 *   
 */
 

uint64_t sam_getPos(samrec_t *samrec) {
  return samrec->pos;
}

unsigned sam_getFlag(samrec_t *samrec) {
  return samrec->flag;
}

char *sam_getRname(samrec_t *samrec) {
  return samrec->rname;
}

uint8_t sam_getMapq(samrec_t *samrec) {
  return samrec->mapq;
}

char *sam_getCigar(samrec_t *samrec) {
  return samrec->cigar;
}

char *sam_getRnext(samrec_t *samrec) {
  return samrec->rnext;
}

uint64_t sam_getPnext(samrec_t *samrec) {
  return samrec->pnext;
}

uint64_t sam_getTlen(samrec_t *samrec) {
  return samrec->tlen;
}

char *sam_getSeq(samrec_t * samrec) {
  return samrec->seq;
}

char *sam_getQual(samrec_t *samrec) {
  return samrec->qual;
}

unsigned int sam_getNooftags (samrec_t *samrec) {
  return samrec->nooftags;
}

char sam_getRC(samrec_t *samrec){
  return ((samrec->flag & 0x10) > 0);
}

char sam_isMapped(samrec_t *samrec) {
  return !(samrec->flag & 0x4);
}

/*------------------------------- sam_getTag ---------------------------------
 *    
 * @brief get the value of a sam tag
 * @author Steve Hoffmann 
 *   
 */
 
samtag_t*
sam_getTag (samrec_t *rec, char* key)
{
  Uint i;
  samtag_t *tag = NULL;

  for(i=0; i < rec->nooftags; i++) {
    if(strncmp(rec->tags[i].key, key, 2) ==0) {
      tag = &rec->tags[i];
    }
  }
  
  return tag;
}

/*-------------------------- sam_getPrevSplitRefRc --------------------------
 *    
 * @brief get the reference strand of the prev trans split
 * @author Steve Hoffmann 
 *   
 */
 
char 
sam_getPrevSplitRefRc (samrec_t *samrec) {
  char rc = 0;
  samtag_t *tag;

  tag = sam_getTag(samrec, "XS");
  if(tag) {
    rc = (atoi(tag->val)) ? 0 : 1;
  }
  
  return rc;
}


/*-------------------------- sam_getPrevSplitRefPos --------------------------
 *    
 * @brief get the reference position of the prev trans split
 * @author Steve Hoffmann 
 *   
 */
 

uint64_t sam_getPrevSplitRefPos (samrec_t *samrec) {
  uint64_t pos = 0;
  samtag_t *tag;

  tag = sam_getTag(samrec, "XU");
  if(tag) {
    pos = strtoull(tag->val, NULL, 10);
  }
  
  return pos;
}

/*-------------------------- sam_getNextSplitRefRc --------------------------
 *    
 * @brief get the reference strand of the next trans split
 * @author Steve Hoffmann 
 *   
 */
 
char 
sam_getNextSplitRefRc (samrec_t *samrec) {
  char rc = 0;
  samtag_t *tag;

  tag = sam_getTag(samrec, "XT");
  if(tag) {
    rc = (atoi(tag->val)) ? 0 : 1;
  }
  
  return rc;
}


/*-------------------------- sam_getNextSplitRefPos --------------------------
 *    
 * @brief get the reference position of the next trans split
 * @author Steve Hoffmann 
 *   
 */
 

uint64_t sam_getNextSplitRefPos (samrec_t *samrec) {
  uint64_t pos = 0;
  samtag_t *tag;

  tag = sam_getTag(samrec, "XV");
  if(tag) {
    pos = strtoull(tag->val, NULL, 10);
  }
  
  return pos;
}

/*-------------------------- sam_getPrevSplitRefChr --------------------------
 *    
 * @brief get the reference chromosome of the prev trans split
 * @author Steve Hoffmann 
 *   
 */
 
char*
sam_getPrevSplitRefChr (samrec_t *samrec)
{
  samtag_t *tag;
  char *ref = NULL;

  tag = sam_getTag(samrec, "XP");
  if(tag) {
    ref = tag->val;
  }
  
  return ref;
}


/*-------------------------- sam_getNextSplitRefChr --------------------------
 *    
 * @brief get the reference chromosome of the next trans split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
sam_getNextSplit (samrec_t *samrec, samheader_t *head, uint64_t *pos, char *rc)
{
  samtag_t *tag;
  Uint i;
  stringset_t* set;
  char *ref = NULL;

  tag = sam_getTag(samrec, "XC");
  
  if(tag) {
    
    set = tokensToStringset(NULL, ",", tag->val, strlen(tag->val));
    assert(set->noofstrings == 7);  
    ref = set->strings[0].str;
    *pos = strtoull(set->strings[1].str, NULL, 10);
    *rc = (set->strings[2].str[0] == '-') ? 1 : 0;
  }
  
  if(!ref) {
    return -1;
  }

  for(i=0; i < head->nrnames; i++){
    if(!strcmp(ref, head->rnames[i])) break;
  }

  assert(i < head->nrnames);
  destructStringset(NULL, set);

  return i;
}


/*-------------------------- sam_getNextSplitRefIdx --------------------------
 *    
 * @brief get the index of the reference chromsome of the next trans split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
sam_getNextSplitRefIdx (samrec_t *samrec, samheader_t* head)
{
  char *ref = NULL;
  Uint i = 0;
  
  //ref = sam_getNextSplitRefChr(samrec);
  fprintf(stderr, "ref: '%s'\n", ref);

  if(!ref) {
    return -1;
  }

  for(i=0; i < head->nrnames; i++){
    if(!strcmp(ref, head->rnames[i])) break;
  }

  assert(i < head->nrnames);

  return i;
}
/*------------------------------ sam_getRefIdx -------------------------------
 *    
 * @brief get the index of the reference chromsome of the prev trans split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
sam_getRefIdx (samrec_t *samrec, samheader_t* head)
{
  char *ref;
  Uint i = 0;
  
  ref = sam_getRname(samrec);

  if(!ref) {
    return -1;
  }

  for(i=0; i < head->nrnames; i++){
    if(!strcmp(ref, head->rnames[i])) break;
  }

  assert(i < head->nrnames);

  return i;
}

/*-------------------------- sam_getPrevSplitRefIdx --------------------------
 *    
 * @brief get the index of the reference chromsome of the prev trans split
 * @author Steve Hoffmann 
 *   
 */
 
Uint
sam_getPrevSplitRefIdx (samrec_t *samrec, samheader_t* head)
{
  char *ref;
  Uint i = 0;
  
  ref = sam_getPrevSplitRefChr(samrec);

  if(!ref) {
    return -1;
  }

  for(i=0; i < head->nrnames; i++){
    if(!strcmp(ref, head->rnames[i])) break;
  }

  assert(i < head->nrnames);

  return i;
}

Uint
sam_getPrevSplit (samrec_t *samrec, samheader_t *head, uint64_t *pos, char *rc)
{
  samtag_t *tag;
  Uint i;
  stringset_t* set;
  char *ref = NULL;

  tag = sam_getTag(samrec, "XP");
  if(tag) {
    ref = tag->val;
  }

  
  if(tag) {
    
    set = tokensToStringset(NULL, ",", tag->val, strlen(tag->val));
    assert(set->noofstrings == 7);  
    ref = set->strings[0].str;
    *pos = strtoull(set->strings[1].str, NULL, 10);
    *rc = (set->strings[2].str[0] == '-') ? 1 : 0;
  }
  
  if(!ref) {
    return -1;
  }

  for(i=0; i < head->nrnames; i++){
    if(!strcmp(ref, head->rnames[i])) break;
  }

  assert(i < head->nrnames);

  destructStringset(NULL, set);
  return i;
}


/*--------------------------------- setters ----------------------------------
 *    
 * @brief setters
 * @author Steve Hoffmann 
 *   
 */
 

void sam_addTag(samrec_t *rec, const char *fmt, ...) {
  int len;
  Uint k = 0;
  char *xval, *ptr, *mycopy;
  char *saveptr;
  va_list ap;

  va_start(ap, fmt);
  len = vsnprintf(NULL, 0, fmt, ap);
  va_end(ap);

  xval = ALLOCMEMORY(space, NULL, char, len+1);
  
  va_start(ap, fmt);
  vsprintf(xval, fmt, ap); 
  va_end(ap);

  rec->tags = ALLOCMEMORY(NULL, rec->tags, samtag_t, 
      rec->nooftags+1);
  rec->tags[rec->nooftags].tag = xval;

  mycopy = my_strdup(xval);
  ptr = strtok_bl(mycopy, ":", &saveptr);

  while(ptr != NULL) {
    switch(k) {
      case 0:
        rec->tags[rec->nooftags].key = my_strdup(ptr);
        break;
      case 1:
        rec->tags[rec->nooftags].type = my_strdup(ptr);
        break;
      case 2:
        rec->tags[rec->nooftags].val = my_strdup(ptr);
        break;
      default:
        if(k < 2) { 
          fprintf(stderr, "samout error: malformed tag\n");
          exit(-1);
        } else {
          char *tmp = my_strdup(ptr);
          Uint tmplen = strlen(tmp);
          Uint len = strlen(rec->tags[rec->nooftags].val);
          rec->tags[rec->nooftags].val = 
            ALLOCMEMORY(NULL, rec->tags[rec->nooftags].val, char, len+tmplen+1);
          memmove(&rec->tags[rec->nooftags].val[len], tmp, tmplen);
          rec->tags[rec->nooftags].val[len+tmplen] = '\0';
          FREEMEMORY(NULL, tmp);
        }
    }
    ptr = strtok_bl(NULL, ":", &saveptr); 
    k++;
  }
 
  rec->nooftags++;
  FREEMEMORY(NULL, mycopy);

  return;
}

void sam_setMapqDbl(samrec_t *samrec, double mapq, double sigma, double maxqual, double fqual){
  sam_addTag(samrec, "QQ:i:%f", exp(mapq));
  sam_addTag(samrec, "QS:i:%f", sigma);
  sam_addTag(samrec, "QF:i:%f", fqual);
  sam_addTag(samrec, "QM:i:%f", maxqual);
  sam_addTag(samrec, "Q:i:%f", MIN(exp(fqual), exp(mapq)));

}


uint8_t sam_setMapq(samrec_t *samrec, double mapq) {
  uint8_t score;
  double phred = mapq * -4.34294481903252; 

  //PICARD WANTS A 0 QUAL
  if(mapq == 0.0)
    score = 0;
  else if(phred < 1) 
    score = 1;
  else if(phred >= 60)
    score = 60;
  else 
    score = (uint8_t) phred;

  samrec->mapq = score;

  return score;
}

void sam_setCigar(samrec_t *samrec, char *cigar) {
  samrec->cigar = cigar;

  return;
}

void sam_setNext(samrec_t *samrec, char *rnext, unsigned int pnext) {
 
  if(rnext) { 
    if(!strcmp(rnext, samrec->rname)) {
      samrec->rnext = my_strdup("=\0");
    } else { 
      samrec->rnext = my_strdup(rnext);
    }
  } else {
    samrec->rnext = my_strdup("*\0");
  }

  samrec->pnext = pnext+1; //0-based C representation to 1-based sam

  return;
}

void sam_setMultipleHits(samrec_t *samrec, Uint nh) {
  sam_addTag(samrec, "NH:i:%d", nh);

  return;
}

void sam_setEdist(samrec_t *samrec, Uint nm) {
  sam_addTag(samrec, "NM:i:%u\0", nm);

  return;
}

void sam_setMD(samrec_t *samrec, char *md) {
 
  sam_addTag(samrec, "MD:Z:%s\0", md);
  return;
}


/*-------------------------------- sam_setter --------------------------------
 *    
 * @brief setter functions for segemehl private tags
 * @author Steve Hoffmann 
 *   
 */
 
void sam_setReadGroup(samrec_t *samrec, char *readgroup) {
  sam_addTag(samrec, "RG:Z:%s", readgroup);
}


void sam_setMappingType(samrec_t *samrec, char mappingtype){ 
  sam_addTag(samrec, "YZ:Z:%d", mappingtype);
}


void sam_setSplitStart(samrec_t *samrec, Uint start) {
  sam_addTag(samrec, "XX:i:%d", start+1); //0-based C to 1-based sam

  return;
}

void sam_setSplitEnd(samrec_t *samrec, Uint end) {
  sam_addTag(samrec, "XY:i:%d", end+1); //0-based C to 1-based sam

  return;
}

void sam_setSplitNumber(samrec_t *samrec, Uint number, Uint noofparts, Uint totalparts) {
  sam_addTag(samrec, "XI:i:%d", number);  //XQ   starting split in this very alignment
  sam_addTag(samrec, "XH:i:%d", noofparts); //XL number of splits in this very alignment
  sam_addTag(samrec, "XJ:i:%d", totalparts);  //XQ total number of splits in all the alignms of this mapping
  return;
}

void sam_setPrevSplit(samrec_t *samrec, char *rprev, Uint pprev, Uint sprev, Uint prevu, Uint prevlen, Uint preverr, uint8_t prevqual) {
                                            //      chr , pos  , strand     , edist, qual, qstart
  sam_addTag(samrec, "XP:Z:%s,%"PRIu64",%c,%d,%d,%d,%d", rprev, pprev+1, (char)sprev, prevu+1, prevlen, preverr, prevqual);
  /*
  sam_addTag(samrec, "XP:Z:%s", rprev);
  sam_addTag(samrec, "XU:i:%"PRIu64"", pprev+1); //0-based C to 1-based sam
  sam_addTag(samrec, "XS:i:%d", sprev); 
  sam_addTag(samrec, "XR:i:%d", prevlen); 
  sam_addTag(samrec, "XE:i:%d", preverr); 
   */
  return;
}

void sam_setNextSplit(samrec_t *samrec, char *rnext, Uint pnext, Uint snext, Uint nextu, Uint nextlen, Uint nexterr, uint8_t nextqual) {

  sam_addTag(samrec, "XC:Z:%s,%"PRIu64",%c,%d,%d,%d,%d", rnext, pnext+1, (char)snext, nextu+1, nextlen, nexterr, nextqual);
/*  sam_addTag(samrec, "XC:Z:%s", rnext);
  sam_addTag(samrec, "XV:i:%"PRIu64"", pnext+1); //0-based C to 1-based sam
  sam_addTag(samrec, "XT:i:%d", snext);
  sam_addTag(samrec, "XG:i:%d", nextlen);
  sam_addTag(samrec, "XO:i:%d", nexterr);*/

  return;
}

void sam_setBisulfiteProtocol(samrec_t *samrec, segemehl_t *nfo) { 
  
  if (nfo->bisulfiterun == 1) {
    sam_addTag(samrec, "XB:Z:F%u/CT", nfo->bisulfiteprotocol);
  } else if (nfo->bisulfiterun == 2) {
    sam_addTag(samrec, "XB:Z:F%u/GA", nfo->bisulfiteprotocol);
  }

  return;
}

void sam_setBisulfiteMismatches(samrec_t *samrec, Uint mis, Uint misstrand, segemehl_t *nfo) {
  sam_addTag(samrec, "XD:i:%u", mis);
  sam_addTag(samrec, "XF:i:%u", misstrand);

  return;
}

void sam_setMatchId(samrec_t *samrec, Uint matchid) {

  sam_addTag(samrec, "HI:i:%u", matchid);

  return;
}

void sam_setSeedInformation(samrec_t *samrec, char maxevalue, char maxinterval, 
    uint64_t seedstart, uint64_t seedlen, uint64_t refidx, char *refname, 
    uint64_t refpos, char rc) {
  
  sam_addTag(samrec, "ZE:A:%d", maxevalue);
  sam_addTag(samrec, "ZI:A:%d", maxinterval);
  sam_addTag(samrec, "ZM:A:%d", rc);
  sam_addTag(samrec, "ZS:i:%u", seedstart);
  sam_addTag(samrec, "ZL:i:%u", seedlen);
  sam_addTag(samrec, "ZR:i:%u", refidx);
  sam_addTag(samrec, "ZP:i:%u", refpos);
  if(refname)
  sam_addTag(samrec, "ZZ:Z:%s", refname);

  return;
}



void sam_setPartialEdists(samrec_t *samrec, Uint *edist, Uint noofparts) { 
  Uint i, len, totallen=0;
  char *string = NULL;

  for(i=0; i < noofparts; i++){
    len = snprintf(NULL, 0, ",%d", edist[i]);

    string = ALLOCMEMORY(NULL, string, char, len+totallen+1);
//    if(i>0)
      totallen += snprintf(&string[totallen], len+1, ",%d", edist[i]);
//    else
//      totallen += snprintf(&string[totallen], len+1, "%d", edist[i]);
  }

  sam_addTag(samrec, "XM:B:I%s", string); //YN
  FREEMEMORY(NULL, string);
}

void sam_setPartialLength(samrec_t *samrec, Uint *length, Uint noofparts) { 
  Uint i, len, totallen=0;
  char *string = NULL;

  for(i=0; i < noofparts; i++){
    len = snprintf(NULL, 0, ",%d", length[i]);

    string = ALLOCMEMORY(NULL, string, char, len+totallen+1);
//    if(i > 0)
    totallen += snprintf(&string[totallen], len+1, ",%d", length[i]);
//    else
//    totallen += snprintf(&string[totallen], len+1, "%d", length[i]);
  }

  sam_addTag(samrec, "XL:B:I%s", string); //YL 
  FREEMEMORY(NULL, string);
}


/*---------------------------- sam_extractSplits -----------------------------
 *    
 * @brief extract splits
 * @author Steve Hoffmann 
 *   
 */

spliceevents_t*
sam_extractSplits(samrec_t *rec, samheader_t *head) {  
  uint64_t vdonpos, vaccpos, ulen=0, voff;
  Uint vaccidx, vdonidx, idx;
  char vaccrc=0, vdonrc=0, rc=0, *seq=NULL;
  char *cigar;
  spliceevents_t *events;
  splitalignment_t *aln;

  cigar = sam_getCigar(rec);
  voff = sam_getPos(rec);
  seq = sam_getSeq(rec);
  ulen = strlen(seq);
  rc = sam_getRC(rec);

  idx = sam_getRefIdx(rec, head);

  vdonidx = sam_getPrevSplit(rec, head, &vdonpos, &vdonrc);
  if(vdonidx == -1) {
    vdonpos = 0;
    vdonrc = rc;
  }

  vaccidx = sam_getNextSplit (rec, head, &vaccpos, &vaccrc);
  if(vaccidx==-1) {
    vaccpos = 0;
    vaccrc = rc;
  }

  //restoring 0-offset before passing to the SplitAlignment Routine!
  aln = bl_cigarGetSplitAlignment(seq, ulen, 0, voff-1, idx, rc, 
    vdonpos-1, vdonidx, vdonrc, 
    vaccpos-1, vaccidx, vaccrc, cigar); 
  //extract events and dump
  events = bl_splitAlignmentGetSpliceEvents (aln);

//#ifdef DEBUGSPLICEVENTS
  bl_dumpSpliceEvents (events);
  //destruct events and aln
  //bl_destructSpliceEvents (events);
//#endif

  bl_destructSplitAlignment (aln);
  FREEMEMORY(NULL, aln);

  return events;
}

/*-------------------------- sam_encodeSplitStrand ---------------------------
 *    
 * @brief simple encoder for strands of next/prev splits
 * @author Steve Hoffmann 
 *   
 */
 
Uint
sam_encodeSplitStrand (Uint in, char next)
{
  if(!in) return (Uint)'+';
  else return (Uint)'-';
  
/*  if(in == 0 && next) //plus strand next
    return 32;
  if(in == 0 && !next) //plus strand prev
    return 64;
	
  return 0; //otherwise minus strand */
}

/*--------------------------- sam_destructSamList ----------------------------
 *    
 * @brief destroy the sam list
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_destructSamList (samlist_t *list)
{
  Uint i;

  for(i=0; i < list->noofrecs; i++) { 
    sam_destruct(&list->recs[i]);
  }

  FREEMEMORY(NULL, list->recs);
  return ;
}



void sam_setPos(samrec_t *samrec, char *seq, char *qual,
    char *rname, uint64_t pos, Uint fragno, 
    Uint nfrags, char hasPaired,  char properlymapped, char segmentunmapped, char rc, 
    char nextrc, char isQuery, char nomatemapped, char ismultiple, char ischimeric) {

  samrec->seq  = my_strdup(seq);
  samrec->qual = my_strdup(qual);


#ifdef SORTEDUNMAPPED
  char sortedunmapped = 1;
#else
  char sortedunmapped =0;
#endif

  if(segmentunmapped) {
    samrec->flag |= 0x4;
    
    //expects that the unmapped read comes with the position of its mate 
    if(sortedunmapped && rname) {
      samrec->rname = my_strdup(rname);
      samrec->pos = pos; //assume that pos is already in 1-based sam 
    } else { 
      samrec->rname = my_strdup("*\0");
      samrec->pos = 0;
    }

    char *cigar = my_strdup("*\0");
    sam_setCigar(samrec, cigar);
  } else { 
    
    samrec->rname = my_strdup(rname);
    samrec->pos = pos+1; //0-based C representation to 1-based sam
  }

  //template has multiple segments IN SEQUENCING
  if (hasPaired) {
    samrec->flag |= 0x1;
  }

  //each segment properly alignned according to the aligner
  if (hasPaired && properlymapped) {
    samrec->flag |= 0x2;
  }

  //next segment in the template is unmapped
  if((fragno==0 && nomatemapped) || (fragno+1==nfrags && nomatemapped)) {
    samrec->flag |= 0x8;
  }

  //SEQ being reverse complemented
  if (rc) { 
    samrec->flag |= 0x10;
    samrec->qual=strrev(samrec->qual, strlen(samrec->qual));
  }

  //SEQ of next segment in the template in rc
  if (nextrc) {
      samrec->flag |= 0x20;
  }

  //first segment in template (query)
  if (hasPaired && fragno == 0 && isQuery) {
    samrec->flag |= 0x40;
  }

  //segment ist the last in template (mate)
  if(hasPaired && fragno == nfrags-1 && !isQuery) { 
    samrec->flag |= 0x80;
  }

  //multiple mappings
  if(ismultiple) {
    samrec->flag |= 0x100;
  }

  //chimeric alignment
  if(ischimeric) {
    samrec->flag |= 0x800;
  }


  return;
} 

/*----------------------------- bl_getClipCigar ------------------------------
 *    
 * @brief get a clip cigar string
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_getClipCigar (Uint cliplen)
{

  char *cigar;
  Uint len = snprintf(NULL, 0, "S%d", cliplen); 

  cigar = ALLOCMEMORY(NULL, NULL, cigar, len+1);
  snprintf(cigar, len, "S%d", cliplen); 
  
  return cigar;
}


/*--------------------------- sam_restoreClipping ----------------------------
 *    
 * @brief restore the clipping
 * @author Steve Hoffmann 
 *   
 */

  char*
sam_restoreClippingSeq (char **seq, char **qual, fasta_t *reads, Uint id, 
    mapping_t *l, Uint i, Uint *lclip, Uint *rclip, char hardclip, char hasqual)
{

  Uint nfrags = l->n;
  Uint ll=0, rr=0, myll=0, myrr=0;
  char *left = NULL, *right = NULL;
  char *leftqual = NULL, *rightqual = NULL;
  char *myqual = *qual;
  char *myseq = *seq;

  char rc = bl_getMapFragStrand(&l->f[i]);
  char mate = bl_getMapFragIsMate(&l->f[i]);

  if(!mate) { 
    bl_fastaGetClipPos(reads, id, &ll, &rr);
  } else { 
    bl_fastaGetMateClipPos(reads, id, &ll, &rr);
  }

  if(i==0 || bl_getMapFragIsMate(&l->f[i-1]) != mate) { //this is a first fragment of query or mate
  //  fprintf(stderr, "this is a first fragment\n");
    if(!rc && ll) { //if fragment is !rc attach all leftclips to the left
    //  fprintf(stderr, "attach to leftclip to left\n");
      myll = ll;
      if(!hardclip) { 
        if(!mate) {  //in this case of the query
          left = bl_fastaGetLeftClip(reads,id);
          if(hasqual) leftqual = bl_fastaGetLeftClipQual(reads,id);
        } else { //or of the mate
          left = bl_fastaGetMateLeftClip(reads,id);
          if(hasqual) leftqual = bl_fastaGetMateLeftClipQual(reads,id);
        }
      }
    } else if (rc && ll){ //if fragment is rc attach all leftclips to the right
   //   fprintf(stderr, "attach to leftclip to right\n");
      myrr = ll;
      if(!hardclip) { 
        if(!mate) {  //in this case of the query
          right = bl_fastaGetLeftClip(reads,id);
          if(hasqual) rightqual = bl_fastaGetLeftClipQual(reads,id);
        } else { //or of the mate
          right = bl_fastaGetMateLeftClip(reads,id);
          if(hasqual) rightqual = bl_fastaGetMateLeftClipQual(reads,id);
        }
        //and reverse it
        char *tmp = charIUPACcomplement(NULL, right, myrr); //reverse
        if(hasqual) rightqual = strrev(rightqual, myrr);
        FREEMEMORY(space, right);
        right = tmp;
      }
    }
  }


  if(i+1==nfrags || mate != bl_getMapFragIsMate(&l->f[i+1])) { //this is a last fragment

   // fprintf(stderr, "this is a last fragment\n");
    if(!rc && rr) { //if fragment is !rc attach all rightclips to right
     // fprintf(stderr, "attach to rightclip to right\n");
      myrr = rr;
      if(!hardclip) { 
        if(!mate) { //in this case of the query
          right = bl_fastaGetRightClip(reads,id);
          if(hasqual) rightqual = bl_fastaGetRightClipQual(reads,id);
        } else { //or of the mate
          right = bl_fastaGetMateRightClip(reads,id);
          if(hasqual) rightqual = bl_fastaGetMateRightClipQual(reads,id);
        }
      }
    } else if (rc && rr){ //if fragment is rc attach all rightclips to left
    //  fprintf(stderr, "attach to rightclip to left\n");
      myll = rr;
      if(!hardclip) { 
        if(!mate) { //in this case of the query
          left = bl_fastaGetRightClip(reads,id);
          if(hasqual) leftqual = bl_fastaGetRightClipQual(reads,id);
        } else { //or of the mate
          left = bl_fastaGetMateRightClip(reads,id);
          if(hasqual) leftqual = bl_fastaGetMateRightClipQual(reads,id);
        }

        char *tmp = charIUPACcomplement(NULL, left, myll); //reverse
        if(hasqual) leftqual = strrev(leftqual, myll);
        FREEMEMORY(space, left);
        left = tmp;
      }
    }
  } 

  if(myll && !hardclip) {
    char *tmp = concat(NULL, left, myseq, myll, strlen(myseq));
    FREEMEMORY(NULL, myseq);   
    myseq = tmp;
    if(hasqual) { 
      tmp = concat(NULL, leftqual, myqual, myll, strlen(myqual));
      FREEMEMORY(NULL, myqual);
      myqual = tmp;
    }
  }

  if(myrr && !hardclip) {
    myseq = concat(NULL, myseq, right, strlen(myseq), myrr);
    FREEMEMORY(space, right);
    if(hasqual) { 
      myqual = concat(NULL, myqual, rightqual, strlen(myqual), myrr);
      FREEMEMORY(space, rightqual);
    }
  }

  *lclip = myll;
  *rclip = myrr; 
  *seq = myseq;

  if(hasqual) *qual = myqual;

  return myseq;
}

/*------------------------------ sam_getSamList ------------------------------
 *    
 * @brief get the list of sam entries for a number of aligned frags
 * the fragment list is the container that holds the mapped fragments
 * for a single read (single end seq) or a pair (mate pair seq)
 * serialized:
 *

    //nextchr = my_strdup((CharSequence*)mseq->ref[l->frags[i].nextidx].ref)->description;    
    //sam_setNextSplit(&samlist->recs[i], nextchr, l->frags[i].nextpos, l->frags[i].strand);

 * @author Steve Hoffmann 
 *   
 */
 
samlist_t *sam_getSamList (fasta_t *reads, Uint id, mapping_t *l, 
    MultiCharSeq *mseq, char isMultiple, Uint noofqueries, Uint noofmates, 
    Uint queryid, Uint mateid, char unmappedtornext, segemehl_t *nfo) {

  Uint i, nfrags, edist, ll, rr;
  char *cigar, *md, *prevchr, *nextchr;
  //TODO take care of the following variables!
  char properlymapped = 1;
  char unmap = 0;
  char* lastchrom = NULL;
  uint64_t prevpos, nextpos;
  uint8_t prevqual, nextqual;
  Uint prevlen, nextlen, preverr, nexterr, prevstart, nextstart;
  int64_t minpos = 0, maxpos=0, tmp;
  Uint queryidx, mateidx;
  Uint bimis=0, bistrandmis=0;
  MultiCharSeqAlignment* al;
  Uint rightmost;
  Uint startidx = 0;
  Uint endidx = 0;
  char nomatemapped = 0;
 //TODO multiple hits
 //TODO trailing gaps
  samlist_t *samlist;


  nfrags = l->n;

  samlist = ALLOCMEMORY(NULL, NULL, samlist_t, 1);
  samlist->recs = ALLOCMEMORY(NULL, NULL, samrec_t, nfrags);
  samlist->noofrecs = nfrags;

  if(nfrags) {  
    lastchrom = bl_getMapFragRefDesc(&l->f[0]);
    minpos = bl_getMapFragP(&l->f[0]);
    maxpos = bl_getMapFragQ(&l->f[0]);
    rightmost = 0;
  }

  char hasPaired = bl_fastaHasMate(reads);
  char isPaired =  bl_isPairedMapping(l);
  char isChimeric = !l->consecutive;

  if(hasPaired && !isPaired) {
     nomatemapped = 1;
     properlymapped = 0;
  }

  queryidx = bl_getQueryStartIdx(l);
  mateidx = bl_getMateStartIdx(l);
  Uint totalqueryparts =0;
  Uint totalmateparts =0;

  for(i=0; i < nfrags; i++) {
    if(bl_getMapFragIsMate(&l->f[i]))
    totalmateparts += getPartialAlignNumber(l->f[i].mcsa->al);
    else
    totalqueryparts += getPartialAlignNumber(l->f[i].mcsa->al);
  }


  for(i=0; i < nfrags; i++) {
  
    char *seq  = bl_getMapFragQry(&l->f[i]);
    char *qual = bl_getMapFragQual(&l->f[i]);
    char isQuery = !bl_getMapFragIsMate(&l->f[i]);
    char rc = bl_getMapFragStrand(&l->f[i]);
    char mappingtype = bl_mappingGetType(l, !isQuery);
    char nextrc = 0;
    char *myqual;

    //all sequences should be already in IUPAC / non-seed
    if(nfo->bisulfite) {
      bl_reconvertBisulfite(seq, strlen(seq), nfo->bisulfite);
    }

    sam_init(&samlist->recs[i], bl_getMapFragQryDesc(&l->f[i]));
 


    if(!bl_getMapFragIsMate(&l->f[i])) {
      sam_setMatchId(&samlist->recs[i], queryid); 
      sam_setMultipleHits(&samlist->recs[i], noofqueries);
    } else {

      sam_setMatchId(&samlist->recs[i], mateid); 
      sam_setMultipleHits(&samlist->recs[i], noofmates);
    }
    
    if(properlymapped && lastchrom == bl_getMapFragRefDesc(&l->f[i]) ) {
      tmp = MIN(minpos, bl_getMapFragP(&l->f[i]));
      minpos = tmp;
      tmp = MAX(maxpos, bl_getMapFragQ(&l->f[i]));
      if(tmp != maxpos) rightmost = i;
      maxpos = tmp;
    } else {
      properlymapped = 0;
    }
      
    lastchrom = bl_getMapFragRefDesc(&l->f[i]);

    
#ifdef NOPICARD
    //NEXTRC: for the last READ this refers to the first READ
    if(nfrags > 1) { 
      if(i+1 < nfrags){ 
        nextrc = bl_getMapFragStrand(&l->f[i+1]);
      } else {
        nextrc = bl_getMapFragStrand(&l->f[0]);
      }
    }
#else
    if(hasPaired && isPaired) {
      if(isQuery) {
        nextrc = bl_getMapFragStrand(&l->f[mateidx]);
      } else {
        nextrc = bl_getMapFragStrand(&l->f[queryidx]);
      }
    }
#endif

    //check the clipping 
    ll = rr = 0;
    sam_restoreClippingSeq (&seq, &qual, reads, id, l, i, &ll, &rr, nfo->hardclip, (qual!=NULL));

    if(qual) {
      myqual = qual;
    } else {
      myqual = "*\0";
    }
    
    cigar = cigarstring(bl_getMapFragAlignment(&l->f[i]), ll, 
        rr, (nfo->hardclip) ? 'H':'S', 0, nfo->briefcigar);
   
    uint64_t mypos= bl_getMapFragP(&l->f[i]) - bl_getMapFragSubstart(&l->f[i]);

    sam_setPos(&samlist->recs[i], seq, myqual, bl_getMapFragRefDesc(&l->f[i]), 
        mypos, i, nfrags, hasPaired, properlymapped, unmap, 
        rc, nextrc, isQuery, nomatemapped, isMultiple, isChimeric);

    FREEMEMORY(space, qual);
    FREEMEMORY(space, seq);

   
    if(nfo->mappingqual == 0) { 
      sam_setMapq(&samlist->recs[i], MAX(l->mapqual, l->f[i].mapq_dbl));
    } else { 
      sam_setMapq(&samlist->recs[i], ((double)nfo->mappingqual)/-4.34294481903252);
    }
//    sam_setMapq(&samlist->recs[i], l->f[i].mapq_dbl);
//    sam_setMapqDbl(&samlist->recs[i], l->mapqual, l->sigma, l->maxqual, l->f[i].mapq_dbl);

    sam_setCigar(&samlist->recs[i], cigar);


    edist =  bl_getMapFragEdist(&l->f[i]);
    sam_setEdist(&samlist->recs[i], edist);

    md = mdstring(bl_getMapFragAlignment(&l->f[i]), 0);
    sam_setMD(&samlist->recs[i], md);
    FREEMEMORY(space, md);

    //next field for the query
    if(!bl_getMapFragIsMate(&l->f[i]) && mateidx != -1) {
      nextchr = bl_getMapFragRefDesc(&l->f[mateidx]);    
      nextpos = bl_getMapFragP(&l->f[mateidx]) - bl_getMapFragSubstart(&l->f[mateidx]);
      sam_setNext(&samlist->recs[i], nextchr, nextpos);
    }

    //next field for the mate
    if(bl_getMapFragIsMate(&l->f[i]) && queryidx != -1) {
       prevchr = bl_getMapFragRefDesc(&l->f[queryidx]);    
       prevpos = bl_getMapFragP(&l->f[queryidx]) - bl_getMapFragSubstart(&l->f[queryidx]);
       sam_setNext(&samlist->recs[i], prevchr, prevpos);
    }

    //set the RNEXT to CUR in case there is no mate
    if(!isPaired && unmappedtornext) {
      nextchr = bl_getMapFragRefDesc(&l->f[i]);
      nextpos = mypos;
      sam_setNext(&samlist->recs[i], nextchr, nextpos);
    }

    //once we have reached the mate, reset the counter
    if(i == mateidx) endidx = 0;
    
    startidx = endidx;
    Uint noofparts=0;
    Uint *edists = getSplitEdist(l->f[i].mcsa->al, &noofparts);
    Uint *lengths = getUPartialAlignlen(l->f[i].mcsa->al, &noofparts);
    endidx += noofparts;

    if(noofparts>1) { 
      sam_setPartialEdists(&samlist->recs[i], edists, noofparts); 
      sam_setPartialLength(&samlist->recs[i], lengths, noofparts); 
    }

    FREEMEMORY(NULL, edists);
    FREEMEMORY(NULL, lengths);

    //check next split
    if(i+1 < nfrags && l->f[i+1].issplit && l->f[i+1].mate == l->f[i].mate) {        
      nextchr = bl_getMapFragRefDesc(&l->f[i+1]);    
      
      //if(bl_getMapFragStrand(&l->f[i+1])) {
      //  nextpos = bl_getMapFragQ(&l->f[i+1]) - bl_getMapFragSubstart(&l->f[i+1]);
      //} else { 
        nextpos = bl_getMapFragP(&l->f[i+1]) - bl_getMapFragSubstart(&l->f[i+1]);
      //}

      nextlen = bl_getMapFragGetUlen(&l->f[i+1]);
      nexterr = bl_getMapFragEdist(&l->f[i+1]);
      nextstart = bl_getMapFragU(&l->f[i+1]);
     
      if(nfo->mappingqual == 0) { 
        nextqual = sam_setMapq(&samlist->recs[i+1], MAX(l->mapqual, l->f[i].mapq_dbl));
      } else { 
        nextqual = sam_setMapq(&samlist->recs[i+1], ((double)nfo->mappingqual)/-4.34294481903252);
      }

      sam_setNextSplit(&samlist->recs[i], nextchr, nextpos, 
            sam_encodeSplitStrand(bl_getMapFragStrand(&l->f[i+1]), 1), nextstart, nextlen, nexterr, nextqual);
    }
    
    //check prev split
    if(i > 0 && l->f[i-1].issplit && l->f[i-1].mate == l->f[i].mate) {
      prevchr = bl_getMapFragRefDesc(&l->f[i-1]);    

      //if(bl_getMapFragStrand(&l->f[i-1])) {
        prevpos = bl_getMapFragP(&l->f[i-1]) - bl_getMapFragSubstart(&l->f[i-1]);
      //} else { 
      //  prevpos = bl_getMapFragQ(&l->f[i-1]) - bl_getMapFragSubstart(&l->f[i-1]);
      //}

      prevlen = bl_getMapFragGetUlen(&l->f[i-1]);
      prevstart = bl_getMapFragU(&l->f[i-1]);
      preverr = bl_getMapFragEdist(&l->f[i-1]);
      prevqual = samlist->recs[i].mapq;

      sam_setPrevSplit(&samlist->recs[i], prevchr, prevpos, 
          sam_encodeSplitStrand(bl_getMapFragStrand(&l->f[i-1]), 0), prevstart, prevlen, preverr, prevqual);
    }

    //check cur split
    if(l->f[i].issplit) {
      sam_setSplitStart(&samlist->recs[i], bl_getMapFragU(&l->f[i]));
      sam_setSplitEnd(&samlist->recs[i], bl_getMapFragV(&l->f[i])); 
      Uint totalparts;
      if(bl_getMapFragIsMate(&l->f[i])) 
        totalparts = totalmateparts; 
      else
        totalparts = totalqueryparts;
      
      if(startidx + noofparts > totalparts) {
        fprintf(stdout, "wrong split idx\n");
      }

      sam_setSplitNumber(&samlist->recs[i], startidx, noofparts, totalparts); //i, nfrags
      //cave split might be mate segment!TODO!!!
    }

    //handle bisulfite tags
    if(nfo->bisulfite) {
      al = bl_getMapFragMCSA(&l->f[i]);
      bimis = getBisulfiteMismatches(al->al, nfo->bisulfite);
      bistrandmis = getWrongStrandBisulfiteMismatches(al->al, nfo->bisulfite);
      sam_setBisulfiteMismatches(&samlist->recs[i], bimis, bistrandmis, nfo);
      sam_setBisulfiteProtocol(&samlist->recs[i], nfo); 
    }
  
    sam_setReadGroup(&samlist->recs[i], nfo->readgroupid);
    sam_setMappingType(&samlist->recs[i], mappingtype);
  }


  if(nfrags > 1 && properlymapped) {
    for(i=0; i < nfrags; i++) {
      if(maxpos-minpos+1 > nfo->maxpairinsertsize) 
        ;//samlist->recs[i].flag &= 0xFF7F;
      samlist->recs[i].tlen = maxpos-minpos+1;
      if(i==rightmost) samlist->recs[i].tlen *= -1;
    }
  }

  if(!properlymapped) {
    for(i=0; i < nfrags; i++) {
      ;//samlist->recs[i].flag &= 0xFF7F;
    }
  }

  return samlist;
}

/*-------------------------- sam_printSamrec2Buffer --------------------------
 *    
 * @brief write the sam entry to a buffer
 * @author Steve Hoffmann 
 *   
 */
 
char*
sam_printSamrec2Buffer (samrec_t* r, char lf)
{
  char *tmp = NULL;
  Uint i;

  bl_bsprintf(&tmp, "%s\t%u\t%s\t%ju\t%u\t%s\t",  
      r->qname, r->flag, r->rname, r->pos, r->mapq, r->cigar);

  if(r->rnext) {
    bl_bsprintf(&tmp,"%s\t%ju\t%jd\t", r->rnext, r->pnext, r->tlen);
  } else {
    bl_bsprintf(&tmp,"*\t0\t0\t");
  }

  bl_bsprintf(&tmp, "%s\t%s\t", r->seq, r->qual);


  for(i=0; i < r->nooftags; i++) {
    bl_bsprintf(&tmp, "%s", r->tags[i].tag);
    if(i < r->nooftags-1) bl_bsprintf(&tmp,"\t");
  }

  bl_bsprintf(&tmp, "%c",lf);

  return tmp;
}

/*----------------------------- sam_printSamrec ------------------------------
 *    
 * @brief print one sam record
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_printSamrec (FILE *dev, samrec_t* r, char lf)
{
  unsigned int i;

   
  fprintf(dev, "%s\t%u\t%s\t%ju\t%u\t%s\t ",
      r->qname, r->flag, r->rname, r->pos, r->mapq, r->cigar);

  if(r->rnext) {
    fprintf(dev,"%s\t%ju\t%jd\t", r->rnext, r->pnext, r->tlen);
  } else {
    fprintf(dev,"*\t0\t0\t");
  }

  fprintf(dev, "%s\t%s\t", r->seq, r->qual);


  for(i=0; i < r->nooftags; i++) {
    fprintf(dev, "%s", r->tags[i].tag);
    if(i < r->nooftags-1) fprintf(dev,"\t");
  }

  fprintf(dev, "%c",lf);
  return ;
}

/*------------------------------ sam_getDevice -------------------------------
 *    
 * @brief get the approriate device for a samrec_t
 * @author Steve Hoffmann 
 *   
 */
 
FILE*
sam_getDevice (samrec_t *rec, bl_fileBin_t **fx, segemehl_t *nfo)
{
  FILE *dev;
  char *bisulfite;
  uint64_t pos = 0;
  char *chr = NULL;
  
  if(!nfo->bins) {
    if(rec || !nfo->samnomatchdev) { 
      dev = nfo->dev;
    } else { 
      dev = nfo->samnomatchdev;
    }
      
    if (nfo->threadno > 1) {
      pthread_mutex_lock(nfo->mtx3);
    }

  } else {
    
    if (rec) { 
      pos = sam_getPos(rec);
      chr = sam_getRname(rec);  
    } else {
      pos = 0;
      chr = NULL;
    }

    if (!nfo->bisulfitemerging){
      
      *fx = bl_fileBinsDomainGetBin(nfo->bins, chr, pos);
      bl_fileBinsLock(*fx);
      dev=bl_fileBinsOpen(NULL, *fx, "w");    

    } else {
      bisulfite = calloc(MAX_INT_LENGTH+1, 1);
      sprintf(bisulfite, "%u", (nfo->bisulfite + 1) % 2);
      *fx = bl_fileBinsDomainGetBin(nfo->bins, bisulfite, nfo->threadid);
      
      //DBG("info.bisulfite=%u\tbisulfite=%s\tthreadid=%u\tfilename=%s\tstart=%llu\tend=%llu\n",
      //    nfo->bisulfite, bisulfite, nfo->threadid, (*fx)->fname, (*fx)->id->start, (*fx)->id->end);
      
      bl_fileBinsLock(*fx);
      dev=bl_fileBinsOpen(NULL, *fx, "w");
      free(bisulfite);
    }
  } 

  return dev;
}


/*----------------------------- sam_closeDevice ------------------------------
 *    
 * @brief close the locked device
 * @author Steve Hoffmann 
 *   
 */

void
sam_closeDevice (bl_fileBin_t *bin, segemehl_t *nfo)
{
  if(nfo->bins){
    bl_fileBinsUnlock(bin);
  } else { 
    if (nfo->threadno > 1) {
      pthread_mutex_unlock(nfo->mtx3);
    }
  }

  return ;
}

/*----------------------------- sam_printSamlist -----------------------------
 *    
 * @brief print a list of sams
 * @author Steve Hoffmann 
 *   
 */

  void
sam_printSamlist (samlist_t *l, mapping_t* m, segemehl_t *nfo)

{
  char lf = '\n', alignlf = '\n';
  unsigned int i;
  FILE *dev = stdout;
  bl_fileBin_t *fx = NULL;
  MultiCharSeqAlignment *al;


  if((nfo->order || nfo->bisulfitemerging) && nfo->align){
    lf = 7;
    alignlf = 7;
  }

  for(i=0; i < l->noofrecs; i++) {
    if(!nfo->bufferedwrite || nfo->bins) { 
      dev = sam_getDevice(&l->recs[i], &fx, nfo);
      sam_printSamrec(dev, &l->recs[i], lf);
      
      if(nfo->align) {
        al = bl_getMapFragMCSA(&m->f[i]);  
        showAlignLF(al->al, dev, alignlf);
        fprintf(dev, "\n");
      }

      sam_closeDevice(fx, nfo);
      fx = NULL;
    } else { 
      char *line = sam_printSamrec2Buffer(&l->recs[i], lf);      
      bl_circBufferAddSave(&nfo->sambuffer[nfo->threadid], line, strlen(line));
      FREEMEMORY(NULL, line);
      if(nfo->align) {
        al = bl_getMapFragMCSA(&m->f[i]);  
        line = getAlignString(al->al, alignlf);
        bl_circBufferAddSave(&nfo->sambuffer[nfo->threadid], line, strlen(line));
        FREEMEMORY(NULL, line);
      }
    }
  }

  return ;
}



/*--------------------------- sam_printEmptyAlign ----------------------------
 *    
 * @brief print empty alignment in SAM format
 * @author Steve Hoffmann 
 *   
 */

void
sam_printEmptyAlign (char *desc, char* seq, char* qual, char hasPaired, char isQuery,
    char nomatemapped, char *nextchr, int64_t nextpos, char nextrc, char ismultiple, char ischimeric, mapseed_t *seed,  
    segemehl_t *nfo)
{

 // if(!nfo->samunmapped) return;

  samlist_t *samlist;
  char rc = 0;
  Uint nfrags = 1;
  Uint this = 0;
  char properlymapped = 0;
  char unmapped = 1; 
  bl_fileBin_t *fx = NULL;
  char lf = '\n';
  char *myqual;
  FILE *dev = NULL;

  samlist = ALLOCMEMORY(NULL, NULL, samlist_t, 1);
  samlist->recs = ALLOCMEMORY(NULL, NULL, samrec_t, 1);
  samlist->noofrecs = 1;

  if(qual) {
    myqual = qual;
  } else {
    myqual = "*\0";
  }


  sam_init(&samlist->recs[0], desc);
  
  //here the position of the unmapped mate is set to the pos of its mapped mate
  //it should be set to NULL if no mate has been aligned

  sam_setPos(&samlist->recs[0], seq, myqual, nextchr, nextpos, this, nfrags, hasPaired,
      properlymapped, unmapped, rc, nextrc, isQuery, nomatemapped, ismultiple, 
      ischimeric);

  if(!nomatemapped) { 
    sam_setNext(&samlist->recs[0], nextchr, nextpos);
  } else { 
    sam_setNext(&samlist->recs[0], nextchr, nextpos-1);
  }

  sam_setMapq(&samlist->recs[0], 0.0);
 
  if(seed) { 
    sam_setSeedInformation(&samlist->recs[0], seed->maxevalue, seed->maxinterval, 
        seed->u, seed->seedlen, seed->refidx, seed->refname, seed->refpos, seed->rc);
  } else {
     // set seed information when no seed is available?
  }
  
  sam_setReadGroup(&samlist->recs[0], nfo->readgroupid);
  sam_setMultipleHits(&samlist->recs[0], 1);
  sam_setMatchId(&samlist->recs[0], 0);

  if(!nfo->bufferedwrite || nfo->bins) {  
    dev = sam_getDevice(NULL, &fx, nfo);
    sam_printSamrec(dev, &samlist->recs[0], lf);
    sam_closeDevice(fx, nfo);
  } else {
    char *line = sam_printSamrec2Buffer(&samlist->recs[0], lf);
    bl_circBufferAddSave(&nfo->sambuffer[nfo->threadid], line, strlen(line));
    FREEMEMORY(NULL, line);
  }


  
    sam_destructSamList(samlist);
  FREEMEMORY(NULL, samlist);

  return ;
}

//#define JOINDEBUG

/*--------------------------- bl_mappingJoinFrags ----------------------------
 *    
 * @brief join split frags to obtain padded alignments
 * @author Steve Hoffmann 
 *   
 */

mappingset_t*
sam_mappingJoinFrags (mappingset_t *s, segemehl_t *nfo)
{
  unsigned int i, j, k=0, adjust =-1;
  Uint head, tail, next;
  MultiCharSeqAlignment **al = NULL;
  mappingset_t *newset;
  mapping_t *newmapping;

  newset = ALLOCMEMORY(NULL, NULL, mappingset_t, 1);
  bl_initMappingSet(newset);
  

  /*heads and tails are assigned to the first and the last sequence*/
  for(i=0; i < s->n; i++) {  
    next = 0;
    head = 0;
    tail = 0;

    newmapping = ALLOCMEMORY(NULL, NULL, mapping_t, 1);
    bl_initMapping(newmapping, s->elem[i].seq, s->elem[i].qual, 
        s->elem[i].lclip, s->elem[i].rclip);
 //   fprintf(stderr, "iterating over %d fragments\n", s->elem[i].n);
    //first fragment of alignment, we have to take care of the head
    for(j=1; j < s->elem[i].n ; j++) {
  //    fprintf(stderr, "this is fragment %d\n", j);
      if(s->elem[i].f[j].issplit == s->elem[i].f[j-1].issplit && 
          s->elem[i].f[j].mate == s->elem[i].f[j-1].mate &&
          bl_getMapFragStrand(&s->elem[i].f[j]) == bl_getMapFragStrand(&s->elem[i].f[j-1]) &&
          bl_getMapFragChrIdx(&s->elem[i].f[j]) == bl_getMapFragChrIdx(&s->elem[i].f[j-1]) && 

          ((bl_getMapFragP(&s->elem[i].f[j]) > bl_getMapFragQ(&s->elem[i].f[j-1]) && 
            !bl_getMapFragStrand(&s->elem[i].f[j])) || 
           (bl_getMapFragQ(&s->elem[i].f[j]) < bl_getMapFragP(&s->elem[i].f[j-1]) &&
            bl_getMapFragStrand(&s->elem[i].f[j])))) {

        //if all in order: add the fragments to a chain in the read order
        //clipping should not be a problem because it is masked since the sequence
        //was returned from fasta in libs/match.c
        if(k == 0) {
 
          if(j == 1 || s->elem[i].f[j-1].mate != s->elem[i].f[j-2].mate) { 
            head = 1;
          }
          if(j == s->elem[i].n-1 || s->elem[i].f[j].mate != s->elem[i].f[j+1].mate) { 
            tail = 1; 
          }

          al = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment*, k+2);
          al[0] = s->elem[i].f[j-1].mcsa;
          al[1] = s->elem[i].f[j].mcsa;
          k = 2;
          next = j+1;
        } else {   

          if(j == s->elem[i].n-1 || s->elem[i].f[j].mate != s->elem[i].f[j+1].mate) {             
            tail = 1;
          }

          al = ALLOCMEMORY(NULL, al, MultiCharSeqAlignment*, k+1);
          al[k] = s->elem[i].f[j].mcsa;
          next = j+1;
          k++;
        }

#ifdef JOINDEBUG
        char *cigar, *md;
        cigar = cigarstring(bl_getMapFragAlignment(&s->elem[i].f[j]), 
            s->elem[i].f[j].rclip, s->elem[i].f[j].lclip, 
            (nfo->hardclip) ? 'H':'S', 0, nfo->briefcigar);
        
        md = mdstring(s->elem[i].f[j].mcsa->al, 0);
        
        fprintf(stdout, "joining map frags %d - %d, cigar:%s, md:%s (head:%d,tail:%d)\n", 
            j-1, j, cigar, md, head, tail);
        
        FREEMEMORY(NULL, cigar);
        FREEMEMORY(NULL, md);
#endif

      } else {

        //otherwise check if there is already a chain that is in order and output it
        if(k > 0) { 
   
          Uint rev = bl_getMapFragStrand(&s->elem[i].f[j-1]);
          Uint lsize=0, rsize=0;
          MultiCharSeqAlignment **myal = al;
 
          if(head) { 
            if(!rev) {
              lsize = myal[0]->al->uoff;
            } else {
              rsize = myal[0]->al->ulen - (myal[0]->al->uoff + getUalignlen(myal[0]->al));
            }
          }

          if(tail) { 
            if(rev) {
              lsize = myal[k-1]->al->uoff;
            } else {
              rsize = myal[k-1]->al->ulen - (myal[k-1]->al->uoff + getUalignlen(myal[k-1]->al));
            }
          }
 
          MultiCharSeqAlignment* res = joinalignments(myal, k, 
              bl_getMapFragStrand(&s->elem[i].f[j-1]), 
              head, tail, lsize, rsize, nfo->cliptype);

          bl_addMapFrag(newmapping, res, NULL, bl_getMapFragIsMate(&s->elem[i].f[j-1]), 
              bl_getMapFragSplit(&s->elem[i].f[j-1]));
  

          FREEMEMORY(space, al);
          
          k = 0;
        } else {
          //if this is not the case start to report the alignments in between that are in trans 
          Uint lsize=0, rsize=0;
          Uint mystart = 0, myend = 0, nextstart = 0, prevend = 0;
          Uint rev, nextrev, prevrev;
          MultiCharSeqAlignment **myal = &s->elem[i].f[j-1].mcsa;            

          //unified coords of prev fragment (if available)
          if(j > 1) {
            prevrev = bl_getMapFragStrand(&s->elem[i].f[j-2]);
            prevend = bl_getMapFragV(&s->elem[i].f[j-2]) + 1;
          }
          //unified coords of this fragment
          rev = bl_getMapFragStrand(&s->elem[i].f[j-1]);           
          mystart = bl_getMapFragU(&s->elem[i].f[j-1]);
          myend = bl_getMapFragV(&s->elem[i].f[j-1]) + 1;
          //unified coords of next fragment
          nextrev = bl_getMapFragStrand(&s->elem[i].f[j]);
          nextstart = bl_getMapFragU(&s->elem[i].f[j]);
          
          if(j == 1 || s->elem[i].f[j-1].mate != s->elem[i].f[j-2].mate) { 
            head = 1;
            if(!rev) {
              //fprintf(stderr, "set lsize 1\n");
              lsize = myal[0]->al->uoff; //mystart
              assert(lsize == mystart);
            } else {
              //fprintf(stderr, "set rsize 2\n");
              rsize = myal[0]->al->ulen - (myal[0]->al->uoff + getUalignlen(myal[0]->al)); //mystart
              assert(rsize == mystart);
            }
           //fprintf(stderr, "this is head in rc:%d with rsize:%d and lsize:%d", rev, rsize, lsize);
          }

          //new mate follows
          if(s->elem[i].f[j-1].mate != s->elem[i].f[j].mate) { 
            //fprintf(stderr, "new mate follows up next\n");
            tail = 1; 
            if(rev) {
              lsize = myal[0]->al->uoff; //opposite of mystart
            } else {
              rsize = myal[0]->al->ulen - (myal[0]->al->uoff + getUalignlen(myal[0]->al)); //opposite of mystart
            }
          }


          //check dangling 5' end of trans alignments
          //there are 4 cases: 1: rev && !prevrev, 2: rev && prevrev, 3: !rev && prevrev, 4: !rev && !prevrev
          if(j >1 && s->elem[i].f[j-2].mate == s->elem[i].f[j-1].mate){ 
            //fprintf(stderr, "the trans alignment is inside a mate sequence 1\n");
            if((rev && !prevrev) || (rev && prevrev)) { 
              //fprintf(stderr, "strand switch (1st case)\n");
              if(adjust != j-2) { //radjust 
              rsize = mystart - prevend;
              }
            } else if(adjust != j-2){ //ladjust
              lsize = mystart - prevend;
            }
          }

          //check dangling 3' ends of trans alignments
          if(s->elem[i].f[j].mate == s->elem[i].f[j-1].mate) { 
            //fprintf(stderr, "the trans alignment is inside a mate sequence (2nd case) rev:%d nextrev:%d\n", rev, nextrev);
            if((rev && !nextrev) || (rev && nextrev)) { 
              // fprintf(stderr, "strand switch (2nd case)\n");
              lsize = nextstart-myend;
              //flag to indicate that dangling end has been taken
              //care of in j-1
              adjust = j-1; //ladjust
            } else {
              //fprintf(stderr, "all others (2nd case) nexstart:%d myend:%d\n", nextstart, myend);
              rsize = nextstart-myend;              
              //flag to indicate that dangling end has
              //been taken care of
              adjust = j-1; //radjust
            }
          }
          //fprintf(stderr, "alignment is shown now:\n");
#ifdef JOINDEBUG
          showAlign(s->elem[i].f[j-1].mcsa->al, stdout);
#endif          
          MultiCharSeqAlignment* res = joinalignments(myal, 1, 
              bl_getMapFragStrand(&s->elem[i].f[j-1]), 
              head, tail, lsize, rsize, nfo->cliptype);

          bl_addMapFrag(newmapping, res, NULL, bl_getMapFragIsMate(&s->elem[i].f[j-1]), 
              bl_getMapFragSplit(&s->elem[i].f[j-1]));
    

          //flag alignment as chimeric
          if(s->elem[i].n > 1 && j>1 && s->elem[i].f[j-2].mate == s->elem[i].f[j-1].mate) { 
            newmapping->consecutive = 0;
          }
        }
        head=0;
        tail=0;
      }
    }

    //finish the current mapping
    if(k > 0) { 
      Uint rev = bl_getMapFragStrand(&s->elem[i].f[j-1]);
      Uint lsize=0, rsize=0;
      MultiCharSeqAlignment **myal = al;

      if(head) { 
        if(!rev) {
          lsize = myal[0]->al->uoff;
        } else {
          rsize = myal[0]->al->ulen - (myal[0]->al->uoff + getUalignlen(myal[0]->al));
        }
      }

      if(tail) { 
        if(rev) {
          lsize = myal[k-1]->al->uoff;
        } else {
          rsize = myal[k-1]->al->ulen - (myal[k-1]->al->uoff + getUalignlen(myal[k-1]->al));
        }
      }
          
      //fprintf(stderr, "joinalignments is called in finish mapping\n");
      MultiCharSeqAlignment* res = joinalignments(myal, k, 
              bl_getMapFragStrand(&s->elem[i].f[j-1]), 
              head, tail, lsize, rsize, nfo->cliptype);

      bl_addMapFrag(newmapping, res, NULL, bl_getMapFragIsMate(&s->elem[i].f[j-1]), 
              bl_getMapFragSplit(&s->elem[i].f[j-1]));

      FREEMEMORY(space, al);
      k = 0;
    }

    if(next < s->elem[i].n) {

      Uint lsize=0, rsize=0;   
      MultiCharSeqAlignment **myal = &s->elem[i].f[j-1].mcsa;
      Uint mystart = 0, prevend = 0;
      Uint rev, prevrev;

      //unified coords of prev fragment (if available)
      if(j > 1) {
        prevrev = bl_getMapFragStrand(&s->elem[i].f[j-2]);
        prevend = bl_getMapFragV(&s->elem[i].f[j-2]) + 1;
      }
      //unified coords of this fragment
      rev = bl_getMapFragStrand(&s->elem[i].f[j-1]);           
      mystart = bl_getMapFragU(&s->elem[i].f[j-1]);


      //this might be the head
      if(j == 1 || s->elem[i].f[j-1].mate != s->elem[i].f[j-2].mate) { 
        head = 1;
        if(!rev) {
          lsize = myal[0]->al->uoff;
        } else {
          rsize = myal[0]->al->ulen - (myal[0]->al->uoff + getUalignlen(myal[0]->al));
        }
  //      fprintf(stderr, "setting head here.\n" );
      }
      //this is always the tail
      tail = 1; 
      if(rev) {
        lsize = myal[0]->al->uoff;
      } else {
        rsize = myal[0]->al->ulen - (myal[0]->al->uoff + getUalignlen(myal[0]->al));
      }

      //check dangling 3' end of trans alignments
      if(j >1 && s->elem[i].f[j-2].mate == s->elem[i].f[j-1].mate){ 
        if((rev && !prevrev) || (rev && prevrev)) { 
          if (adjust != j-2) { //radjust
          rsize = mystart - prevend;
          } 
        } else if(adjust != j-2){ //ladjust
          lsize = mystart - prevend;
        }
      }

#ifdef JOINDEBUG 
      showAlign(s->elem[i].f[j-1].mcsa->al, stdout);
#endif      
      MultiCharSeqAlignment* res = joinalignments(myal, 1, 
              bl_getMapFragStrand(&s->elem[i].f[j-1]), 
              head, tail, lsize, rsize, nfo->cliptype);

        bl_addMapFrag(newmapping, res, NULL, bl_getMapFragIsMate(&s->elem[i].f[j-1]), 
              bl_getMapFragSplit(&s->elem[i].f[j-1]));


 
        //flag alignment as chimeric if this was not the first segment in mapping
        if(s->elem[i].n > 1 && s->elem[i].f[j-2].mate == s->elem[i].f[j-1].mate) { 
          newmapping->consecutive = 0;
        }
    }
    bl_addMapping(newset, newmapping);
    bl_removeMapping(newmapping);
    FREEMEMORY(NULL, newmapping);
  }

  return newset;
}

/*------------------------------ sam_initHeader ------------------------------
 *    
 * @brief initalize header structure
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_initHeader (samheader_t *head)
{
  head->version = NULL;
  head->rnames = NULL;
  head->rlens = NULL;
  head->nrnames = 0;
  head->rgroups = NULL;
  head->rgroupsinfo = NULL;
  head->nrgroups = 0;
  head->cmd = NULL;

  return ;
}


/*---------------------------- sam_destructHeader ----------------------------
 *    
 * @brief destroy the header
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_destructHeader (samheader_t *head)
{
  Uint i;
  FREEMEMORY(NULL, head->version);
  FREEMEMORY(NULL, head->rlens);
  FREEMEMORY(NULL, head->cmd);
  for(i=0; i< head->nrnames; i++){
    FREEMEMORY(NULL, head->rnames[i]);
  }
 
  for(i=0; i< head->nrgroups; i++){
    FREEMEMORY(NULL, head->rgroups[i]);
    FREEMEMORY(NULL, head->rgroupsinfo[i]);
  }

  FREEMEMORY(NULL, head->rnames);
  FREEMEMORY(NULL, head->rgroups);
  FREEMEMORY(NULL, head->rgroupsinfo);
  
  return ;
}



/*----------------------------- sam_addReadGroup -----------------------------
 *    
 * @brief add a read group to a header
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_addReadGroup (samheader_t *head, char *id, char *info)
{

 head->rgroups = ALLOCMEMORY(NULL, head->rgroups, char**, head->nrgroups+1);
 head->rgroupsinfo = ALLOCMEMORY(NULL, head->rgroupsinfo, char**, head->nrgroups+1);
 head->rgroups[head->nrgroups] = id;
 head->rgroupsinfo[head->nrgroups] = info;
 head->nrgroups++;

 return ;
}

/*------------------------ sam_getReadGroupFromHeader ------------------------
 *    
 * @brief extract the reference sequence name and its length
 * @author Steve Hoffmann 
 *   
 */
 
char*
sam_getReadGroupFromHeader (char **xval, Uint ntags, char **rgname)
{
    
  Uint i, len = 0, totallen = 0;
  char checkname=0;
  char *info = NULL;

  for(i=0; i < ntags; i++) {
    if(!strncmp(xval[i],"ID:", 3)) {
      *rgname = my_strdup(&xval[i][3]);
      checkname++;
    } else {
      len = snprintf(NULL, 0, "\t%s", xval[i]);
      info = ALLOCMEMORY(NULL, info, char, totallen+len+1);
      totallen += snprintf(&info[totallen], len+1, "\t%s", xval[i]);
    }
  }


  if(checkname == 1);
	
  return info;
}


/*------------------- sam_getReferenceSequencesFromHeader --------------------
 *    
 * @brief extract the reference sequence name and its length
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_getReferenceSequencesFromHeader (char **xval, Uint ntags, char**rname, uint64_t *len)
{
    
  Uint i;
  char checkname=0, checklen=0;

  for(i=0; i < ntags; i++) {
    if(!strncmp(xval[i],"SN:", 3)) {
      *rname = my_strdup(&xval[i][3]);
      checkname = 1;
    }
    if(!strncmp(xval[i],"LN:", 3)) {
       *len = strtoull(&xval[i][3], NULL, 10);
       checklen = 1;
    }
  }

  assert(checkname && checklen);
	
  return ;
}



/*------------------------------ sam_getHeader -------------------------------
 *    
 * @brief read the header from a sam file
 * @author Steve Hoffmann 
 *   
 */
 
samheader_t*
sam_getHeader (samheader_t *head, char *line)
{

  char *ptr, *saveptr;
  char *mycopy = my_strdup(line); 
  char **xval=NULL;
  Uint ntags = 0, i=0;
  char *rname = NULL;
  char *rgroup = NULL;
  uint64_t rlen =0;


  if(line[0] != '@')
    return NULL;
  
  ptr = strtok_bl(mycopy, "\t", &saveptr); 
  
  while (ptr != NULL) {
    xval = ALLOCMEMORY(NULL, xval, char*, ntags+1); 
    xval[ntags] = my_strdup(ptr);
    ptr = strtok_bl(NULL, "\t", &saveptr);
    ntags++;
  }


  //header line
  if(!strncmp(xval[0], "@HD", 3)) {
    //don't care
  } else if (!strncmp(xval[0], "@SQ", 3)) { 
    sam_getReferenceSequencesFromHeader(&xval[1], ntags-1, &rname, &rlen);
    head->rnames = ALLOCMEMORY(NULL, head->rnames, char*, head->nrnames+1);
    head->rlens = ALLOCMEMORY(NULL, head->rlens, uint64_t, head->nrnames+1);
    head->rnames[head->nrnames] = rname;
    head->rlens[head->nrnames] = rlen;
    head->nrnames++;
  } else if (!strncmp(xval[0], "@RG", 3)) { 
    char *info = sam_getReadGroupFromHeader(&xval[1], ntags-1, &rgroup);
    sam_addReadGroup(head, rgroup, info);
  } else if (!strncmp(xval[0], "@PG", 3)) { 
    //don't care
  } else if (!strncmp(xval[0], "@CO", 3)) {
    //don't care
  }

  for(i=0; i < ntags; i++) {
    FREEMEMORY(NULL, xval[i]);
  }

  FREEMEMORY(NULL, xval);
  FREEMEMORY(NULL, mycopy);
  return head;
}


/*----------------------------- sam_dumpHeader -------------------------------
 *    
 * @brief dump the header
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_dumpHeader (samheader_t  *head)
{
  uint64_t i;

  for(i=0; i < head->nrnames; i++) {
    fprintf(stderr, "found rname %s (%" PRIu64 ")\n", head->rnames[i], head->rlens[i]);
  }
	  
  for(i=0; i < head->nrgroups; i++) {
    fprintf(stderr, "found read group %s\n", head->rgroups[i]);
    fprintf(stderr, "additional info %s\n", head->rgroupsinfo[i]);
  }

  return ;
}

/*------------------------------- sam_line2rec -------------------------------
 *    
 * @brief read a sam line and store it in a samrec
 * @author Steve Hoffmann 
 *   
 */
 
samrec_t *
sam_line2rec(char *line, Uint len, samheader_t *head) {
  Uint i;
  stringset_t *set;
  char *qname, *rname, *cigar, *rnext;
  char *seq, *qual;
  uint64_t pos=0, pnext=0, tlen=0;
  uint8_t mapq;
  int16_t flag=0;
  samrec_t *rec = NULL;
  spliceevents_t *events = NULL;

  
  //header line?
  if(line[0] == '@' || line[0] == '\n') 
    return rec;

  rec = ALLOCMEMORY(NULL, NULL, samrec_t, 1); 
  set = tokensToStringset(NULL, "\t", line, len);

  for(i=0; i < set->noofstrings; i++) {
    //fprintf(stderr, "%d: '%s'\n", i, set->strings[i].str);
    char * cur = set->strings[i].str;
    
    switch(i+1) {
      case 1:
        qname = cur;
        rec = sam_init(rec, qname);
        break;
      case 2: 
        flag = atoi(cur);
        rec->flag = flag;
        break;
      case 3:
        rname = my_strdup(cur);
        rec->rname = rname;
        break;
      case 4:
        pos = strtoull(cur, NULL, 10);
        rec->pos = pos;
        break;
      case 5:
        mapq = atoi(cur);
        sam_setMapq(rec, mapq);
        break;
      case 6:
        cigar = my_strdup(cur);
        sam_setCigar(rec, cigar);
        break;
      case 7:
        rnext = my_strdup(cur);
        rec->rnext = rnext;
        break;
      case 8:
        pnext = strtoull(cur, NULL, 10);
        rec->pnext = pnext;
        break;
      case 9:
        tlen = strtoull(cur, NULL, 10);
        rec->tlen = tlen;
        break;
      case 10:
        seq = my_strdup(cur);
        rec->seq = seq;
        break;
      case 11:
        qual = my_strdup(cur);
        rec->qual = qual;
        break;
      default:
        sam_addTag(rec, cur, NULL);
    }
  }

  if(sam_isMapped(rec)) {
    events=sam_extractSplits(rec, head);
  }

  if(events) { 
    //bl_dumpSpliceEvents(events);
    bl_destructSpliceEvents(events);
    FREEMEMORY(NULL, events);
  }

  destructStringset(NULL, set);

  return rec;
}


samheader_t*
sam_readFile(char* filename)
    //unsigned char gzip, struct access *index,)

{
  FILE *fp;
  off_t offset = 0;
  char ch;
  char *buffer;
//  char *descrbuffer = NULL;
//  char *seqbuffer = NULL;
//  char *qualbuffer = NULL;
//  char idchar=0;
  int ret=0;
//  unsigned char desc = 0;

//  unsigned char qualdesc = 0;
//  unsigned char qual = 0;
//  unsigned char seq = 0;
  unsigned char gzip = 0;
  struct gzidxfile *gzf = NULL;
  struct access * index = NULL;

//  Uint descrlength = 0; 
//  Uint seqlen = 0;
  Uint buffersize = MAXBUFFERSIZE;
//  Uint n = startseq;
  Uint len = 0;  

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
  
  if(gzip) {
    fp = fopen(filename, "rb");
    gzf = bl_initgzidxfile(fp, index, offset, MEDIUMCHUNK);
  } else {
    fprintf(stderr, "open normal.\n");
    fp = fopen(filename, "r");
  }

  if (fp == NULL) {
    NFO("Couldn't open file '%s': %d. Exit forced.\n", filename, errno);
    exit(-1);
  }
  
  if(offset > 0) {
    ret = fseeko(fp, offset, SEEK_SET);
    if (ret == -1) {
      NFO("fseeko failed for file %s. Exit forced.\n", filename);
      exit(-1);
    }
  }

  samheader_t *head;
  head = ALLOCMEMORY(NULL, NULL, samheader_t, 1);
  sam_initHeader(head);

  while((ch= (gzip) ? bl_getgzidxc(gzf) : getc(fp)) != EOF) {


    if(ch == '\n' && buffer){
      //finalize buffer
      buffer = ALLOCMEMORY(space, buffer, char, len+1);  
      buffer[len] = '\0'; 
      //process buffer 
      if(buffer[0] == '@') {
        sam_getHeader(head, buffer);
      } else { 
        sam_line2rec(buffer, len, head);
      }
      //free current buffer
      FREEMEMORY(NULL, buffer);
      //and allocate a new one (or erase the old);
      len = 0;
      buffer = ALLOCMEMORY(space, NULL, char, buffersize);
    }

    if(ch != '\n') {
      len++;
      if(len == buffersize-1) {
        buffersize = 2*buffersize+1;
        buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      }
      buffer[len-1] = ch;
    }
  }


  //sam_dumpHeader(head);

  fclose(fp);
  FREEMEMORY(NULL, buffer);
  return head;
}


