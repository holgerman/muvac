
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
  if(kv->val) FREEMEMORY(NULL, kv->val); 

  return;
}

void sam_destruct(samrec_t *samrec) {
  Uint i;

  if(samrec->qname) FREEMEMORY(NULL, samrec->qname);
  if(samrec->rname) FREEMEMORY(NULL, samrec->rname);
  if(samrec->cigar) FREEMEMORY(NULL, samrec->cigar);
  if(samrec->rnext) FREEMEMORY(NULL, samrec->rnext);
  if(samrec->seq) FREEMEMORY(NULL, samrec->seq);
  if(samrec->qual) FREEMEMORY(NULL, samrec->qual);
  for(i=0; i < samrec->nooftags; i++) {
    sam_destructTag(&samrec->tags[i]);
  }
  FREEMEMORY(NULL, samrec->tags);

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

samtag_t *sam_getSamtag(samrec_t *samrec, Uint i) {
  if(i >= sam_getNooftags(samrec)) {
    return NULL;
  }
  return &samrec->tags[i];
}


/*--------------------------------- setters ----------------------------------
 *    
 * @brief setters
 * @author Steve Hoffmann 
 *   
 */
 

void sam_addTag(samrec_t *rec, const char *fmt, ...) {
  int len;
  char *xval;
  va_list ap;

  va_start(ap, fmt);
  len = vsnprintf(NULL, 0, fmt, ap);
  va_end(ap);

  xval = ALLOCMEMORY(space, NULL, char, len+1);
  
  va_start(ap, fmt);
  vsprintf(xval, fmt, ap); 
  va_end(ap);

  rec->tags = ALLOCMEMORY(NULL, rec->tags, samtag_t, rec->nooftags+1);
  rec->tags[rec->nooftags].val = xval;
  rec->nooftags++;

  return;
}



void sam_setMapq(samrec_t *samrec, uint8_t mapq) {
  samrec->mapq = mapq;
}

void sam_setCigar(samrec_t *samrec, char *cigar) {
  samrec->cigar = cigar;

  return;
}

void sam_setNext(samrec_t *samrec, char *rnext, unsigned int pnext) {
 
  if(!strcmp(rnext, samrec->rname)) {
    FREEMEMORY(NULL, rnext);
    samrec->rnext = my_strdup("=\0");
  } else { 
    samrec->rnext = rnext;
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

void sam_setSplitStart(samrec_t *samrec, Uint start) {
  sam_addTag(samrec, "XX:i:%d", start+1); //0-based C to 1-based sam

  return;
}

void sam_setSplitEnd(samrec_t *samrec, Uint end) {
  sam_addTag(samrec, "XY:i:%d", end+1); //0-based C to 1-based sam

  return;
}

void sam_setSplitNumber(samrec_t *samrec, Uint number, Uint noofsplits) {
  sam_addTag(samrec, "XQ:i:%d", number);
  sam_addTag(samrec, "XL:i:%d", noofsplits);
  
  return;
}

void sam_setPrevSplit(samrec_t *samrec, char *rprev, Uint pprev, Uint sprev) {
  sam_addTag(samrec, "XP:Z:%s", rprev);
  sam_addTag(samrec, "XU:i:%d", pprev+1); //0-based C to 1-based sam
  sam_addTag(samrec, "XS:i:%d", sprev); 

  return;
}

void sam_setNextSplit(samrec_t *samrec, char *rnext, Uint pnext, Uint snext) {
  sam_addTag(samrec, "XC:Z:%s", rnext);
  sam_addTag(samrec, "XV:i:%d", pnext+1); //0-based C to 1-based sam
  sam_addTag(samrec, "XT:i:%d", snext);

  return;
}

void sam_setBisulfiteProtocol(samrec_t *samrec, char protocol, char run, segemehl_t *nfo) { 
  
  if (run == 1) {
    sam_addTag(samrec, "XB:Z:F%u/CT", nfo->bisulfiteprotocol);
  } else if (run == 2) {
    sam_addTag(samrec, "XB:Z:F%u/GA", nfo->bisulfiteprotocol);
  }

  return;
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
  if(in == 0 && next) //plus strand next
    return 32;
  if(in == 0 && !next) //plus strand prev
    return 64;
	
  return 0; //otherwise minus strand
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
    Uint nfrags, char properlymapped, char segmentunmapped, char rc, 
    char nextrc) {

  samrec->seq  = my_strdup(seq);
  samrec->qual = my_strdup(qual);


  if(segmentunmapped) {
    samrec->pos = 0;
    samrec->flag |= 0x4;
    samrec->rname = my_strdup("*\0");
    char *cigar = my_strdup("*\0");
    sam_setCigar(samrec, cigar);

    return;
  } 
  
  samrec->rname = my_strdup(rname);
  samrec->pos = pos+1; //0-based C representation to 1-based sam

  if (fragno == 0) {
    samrec->flag |= 0x40;
  }

  if (nfrags > 1) {
    samrec->flag |= 0x1;

    if(fragno == nfrags-1) { 
      samrec->flag |= 0x80;
    } else if (nextrc) {
      samrec->flag |= 0x20;
    }
  }

  if (properlymapped) {
    samrec->flag |= 0x2;
  }

  if (rc) { 
    samrec->flag |= 0x10;
  }

  return;
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
    MultiCharSeq *mseq, segemehl_t *nfo) {
  Uint i, nfrags, lclip=0, rclip=0, edist;
  char *cigar, *md, *prevchr, *nextchr;
  //TODO take care of the following variables!
  char properlymapped = 1;
  char unmap = 0;
  char* lastchrom = NULL;
  int64_t minpos = 0, maxpos=0, tmp;
  Uint rightmost;

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

  for(i=0; i < nfrags; i++) {

    if(properlymapped && lastchrom == bl_getMapFragRefDesc(&l->f[i])) {
      tmp = MIN(minpos, bl_getMapFragP(&l->f[i]));
      minpos = tmp;
      tmp = MAX(maxpos, bl_getMapFragQ(&l->f[i]));
      if(tmp != maxpos) rightmost = i;
      maxpos = tmp;
    } else {
      properlymapped = 0;
    }

    lastchrom = bl_getMapFragRefDesc(&l->f[i]);

    //if first frag or last frag of query and mate
    //check the clipping
    if(i+1 == nfrags){
      //TODO clipping unchecked!!! this is a stub!!
      bl_fastaGetClipPos(reads, id, &lclip, &rclip);
      bl_fastaGetMateClipPos(reads, id, &lclip, &rclip);
    }

    sam_init(&samlist->recs[i], bl_getMapFragQryDesc(&l->f[i]));
    if(nfrags > 1) {
      samlist->recs[i].flag |= 0x1;
    }

    char *seq  = bl_getMapFragQry(&l->f[i]);
    char *qual = bl_getMapFragQual(&l->f[i]);
    char rc = bl_getMapFragStrand(&l->f[i]);
    char nextrc = 0;
    if(i+1 < nfrags){ 
      nextrc = bl_getMapFragStrand(&l->f[i+1]);
    }
  
    sam_setPos(&samlist->recs[i], seq, qual, bl_getMapFragRefDesc(&l->f[i]), 
        bl_getMapFragP(&l->f[i]), i, nfrags, properlymapped, unmap, rc, nextrc);

    FREEMEMORY(space, qual);
    FREEMEMORY(space, seq);

    
    sam_setMapq(&samlist->recs[i], l->f[i].mapq);

    cigar = cigarstring(bl_getMapFragAlignment(&l->f[i]), l->f[i].rclip, 
        l->f[i].lclip, (nfo->hardclip) ? 'H':'S', 0);
    sam_setCigar(&samlist->recs[i], cigar);


    edist =  bl_getMapFragEdist(&l->f[i]);
    sam_setEdist(&samlist->recs[i], edist);

    md = mdstring(l->f[i].mcsa->al, 0);
    sam_setMD(&samlist->recs[i], md);
    FREEMEMORY(space, md);

    //next field for all but last fragment
    if(i+1 < nfrags) {
      nextchr = my_strdup(bl_getMapFragRefDesc(&l->f[i+1]));    
      sam_setNext(&samlist->recs[i], nextchr, bl_getMapFragP(&l->f[i+1]));
    }

    //next field for last fragment
    if(i > 0 && i == nfrags-1) {
       prevchr = my_strdup(bl_getMapFragRefDesc(&l->f[0]));    
       sam_setNext(&samlist->recs[i], prevchr, bl_getMapFragP(&l->f[0]));
    }

    //check next split
    if(i+1 < nfrags && l->f[i+1].issplit) {        
      nextchr = bl_getMapFragRefDesc(&l->f[i+1]);    
      sam_setNextSplit(&samlist->recs[i], nextchr, bl_getMapFragP(&l->f[i+1]), 
            sam_encodeSplitStrand(bl_getMapFragStrand(&l->f[i+1]), 1));
    }
    
    //check prev split
    if(i > 0 && l->f[i-1].issplit) {
      prevchr = bl_getMapFragRefDesc(&l->f[i-1]);    
      sam_setPrevSplit(&samlist->recs[i], prevchr, bl_getMapFragP(&l->f[i-1]), 
          sam_encodeSplitStrand(bl_getMapFragStrand(&l->f[i-1]), 0));
    }

    //check cur split
    if(l->f[i].issplit) {
      sam_setSplitStart(&samlist->recs[i], bl_getMapFragU(&l->f[i]));
      sam_setSplitEnd(&samlist->recs[i], bl_getMapFragV(&l->f[i]));
      sam_setSplitNumber(&samlist->recs[i], i, nfrags); 
      //cave split might be mate segment!TODO!!!
    }
  }


  if(nfrags > 1) {
    for(i=0; i < nfrags; i++) {
      samlist->recs[i].tlen = maxpos-minpos+1;
      if(i==rightmost) samlist->recs[i].tlen *= -1;
    }
  }

  return samlist;
}


/*----------------------------- sam_printSamrec ------------------------------
 *    
 * @brief print one sam record
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_printSamrec (FILE *dev, samrec_t* r)
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
    fprintf(dev, "%s", r->tags[i].val);
    if(i < r->nooftags-1) fprintf(dev,"\t");
  }

  fprintf(dev, "\n");
  return ;
}


/*----------------------------- sam_printSamlist -----------------------------
 *    
 * @brief print a list of sams
 * @author Steve Hoffmann 
 *   
 */
 
void
sam_printSamlist (FILE *dev, samlist_t *l)

{
  unsigned int i;

  for(i=0; i < l->noofrecs; i++) {
    sam_printSamrec(dev, &l->recs[i]);
  }
	
  return ;
}


/*--------------------------- bl_mappingJoinFrags ----------------------------
 *    
 * @brief join split frags to obtain padded alignments
 * @author Steve Hoffmann 
 *   
 */
 
mappingset_t*
sam_mappingJoinFrags (mappingset_t *s, segemehl_t *nfo)
{
  unsigned int i, j, k=0;
  char *cigar, *md, strand;
  Uint head, tail, next;
  MultiCharSeqAlignment **al;

  /*heads and tails are assigned to the first and the last sequence*/
  for(i=0; i < s->n; i++) {  
    next = 0;
    head = 0;
    tail = 0;
    //first fragment of alignment, we have to take care of the head
    for(j=1; j < s->elem[i].n ; j++) {

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

        strand =  bl_getMapFragStrand(&s->elem[i].f[j-1]);

        if(k == 0) {
          /*if(next < j) { 
            fprintf(stderr, "alignment skipped (%d - %d).\n", next, j);
            for(u=next; u < j; u++) 
            showAlign(s->elem[i].f[u].mcsa->al, stderr);
          }*/

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


        cigar = cigarstring(bl_getMapFragAlignment(&s->elem[i].f[j]), 
            s->elem[i].f[j].rclip, s->elem[i].f[j].lclip, 
            (nfo->hardclip) ? 'H':'S', 0);
    
        md = mdstring(s->elem[i].f[j].mcsa->al, 0);

  //      fprintf(stderr, "joining map frags %d - %d, cigar:%s, md:%s (head:%d,tail:%d)\n", j-1, j, cigar, md, head, tail);
        FREEMEMORY(NULL, cigar);
        FREEMEMORY(NULL, md);
      } else {

        //otherwise check if there is already a chain that is in order and output it
        if(k > 0) { 
          joinalignments(al, k, strand, head, tail);
          FREEMEMORY(NULL, al);
          k = 0;
        } else {
          //if this is not the case start to report the alignments in between that are in trans 
          Uint dstart=0,fstart=0;
          Uint dsize=0, fsize=0;
          Uint lsize=0, rsize=0;

          Uint rev = bl_getMapFragStrand(&s->elem[i].f[j-1]);
          MultiCharSeqAlignment **myal = &s->elem[i].f[j-1].mcsa;

          if(j == 1 || s->elem[i].f[j-1].mate != s->elem[i].f[j-2].mate) { 
            head = 1;
            if(!rev) {
              lsize = myal[0]->al->uoff;
            } else {
              rsize = myal[0]->al->ulen - (myal[0]->al->uoff + getUalignlen(myal[0]->al));
            }
          }
          
          if(s->elem[i].f[j-1].mate != s->elem[i].f[j].mate) { 
            tail = 1; 
            if(rev) {
              lsize = myal[0]->al->uoff;
            } else {
              rsize = myal[0]->al->ulen - (myal[0]->al->uoff + getUalignlen(myal[0]->al));
            }
          }

          Uint mystart = 0;
          Uint myend = 0;

          if(bl_getMapFragStrand(&s->elem[i].f[j-1])) {
            mystart = s->elem[i].f[j-1].mcsa->al->uoff;
            myend   = mystart + getUalignlen(s->elem[i].f[j-1].mcsa->al);
          } else {
            mystart = s->elem[i].f[j-1].mcsa->al->ulen - 
               (s->elem[i].f[j-1].mcsa->al->uoff+getUalignlen(s->elem[i].f[j-1].mcsa->al));
            myend  = mystart + getUalignlen(s->elem[i].f[j-1].mcsa->al);
          }

          fprintf(stderr, "checking j-1:%d vs. j-2:%d\n", j-1, j-2);
          //check dangling 3' end of trans alignments
          if(j >1 && s->elem[i].f[j-2].mate == s->elem[i].f[j-1].mate){ 
          if(!bl_getMapFragStrand(&s->elem[i].f[j-2])) { 
             dstart = s->elem[i].f[j-2].mcsa->al->ulen - 
               (s->elem[i].f[j-2].mcsa->al->uoff+getUalignlen(s->elem[i].f[j-2].mcsa->al));
             dsize = dstart - myend; //add dsize to the right side
             rsize = dsize;
          } else {
            dstart = s->elem[i].f[j-2].mcsa->al->uoff;
            dsize = dstart - myend; //add dsize to the left side
            lsize = dsize;
          }

          fprintf(stderr, "3' end: trans alignment encountered! %d (of %d) with gap of size %d (uoff[j-2]:%d, len[j-2]:%d, invstart:%d, mystart:%d, myend:%d)\n", 
              j-1,  s->elem[i].n, dsize, s->elem[i].f[j-2].mcsa->al->uoff, getUalignlen(s->elem[i].f[j-2].mcsa->al), 
              s->elem[i].f[j-2].mcsa->al->ulen - 
               (s->elem[i].f[j-2].mcsa->al->uoff+getUalignlen(s->elem[i].f[j-2].mcsa->al)),
              mystart, myend);
          }
          
          fprintf(stderr, "checking j-1:%d vs. j:%d\n", j-1, j);
          //check dangling 5' ends of trans alignments
          if(s->elem[i].f[j].mate == s->elem[i].f[j-1].mate) { 
          if( !bl_getMapFragStrand(&s->elem[i].f[j])) { 
             fprintf(stderr, "5 ' case 1\n");
             fstart = s->elem[i].f[j].mcsa->al->ulen- s->elem[i].f[j].mcsa->al->uoff;
             fsize  = mystart - fstart; //left side
             lsize = fsize;
          } else {
            fprintf(stderr, "5' case 2\n");
            fstart = s->elem[i].f[j].mcsa->al->uoff +getUalignlen(s->elem[i].f[j].mcsa->al);
            fsize = mystart - fstart; //right side
            rsize = fsize;
          }
          
          fprintf(stderr, "5' end: trans alignment encountered! %d (of %d) with gap of size %d (dstart:%d, fstart:%d, uoff[j]:%d)\n", 
              j-1,  s->elem[i].n, fsize, fstart, s->elem[i].f[j-1].mcsa->al->uoff, s->elem[i].f[j].mcsa->al->uoff);
          }



          showAlign(s->elem[i].f[j-1].mcsa->al, stderr);
          fprintf(stderr, "my checkpoint rev:%d\n", bl_getMapFragStrand(&s->elem[i].f[j-1]));
          joinalignments2(myal, 1, bl_getMapFragStrand(&s->elem[i].f[j-1]), head, tail, lsize, rsize);
          fprintf(stderr, "my checkpoint end\n");
        }
        head=0;
        tail=0;
      }
    }
  
  if(k > 0) { 
    fprintf(stderr, "start finalize\n");
    joinalignments(al, k, strand, head, tail);
    FREEMEMORY(NULL, al);
    k = 0;
    fprintf(stderr, "end finalize\n");
  }

  if(next < s->elem[i].n) {
    fprintf(stderr, "missing alignment\n");
  }

  fprintf(stderr, "last included %d, current elem %d\n", next, s->elem[i].n-1); 
  }



  return NULL;
}
