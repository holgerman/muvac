
/*
 *  matealign.c
 *  align the mates
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 12/15/2012 20:56:30 CET
 *  
 */

#include "mapfrag.h"
#include "mathematics.h"
#include "bitvectoralg.h"
#include "sufarray.h"
#include "kdseed.h"
#include "segemehl.h"
#include "sort.h"
#include <math.h>
#include <string.h>
#include <limits.h>

/*------------------------- bl_greedypairMappingSet --------------------------
 *    
 * @brief pair a mapping based on the minimum distace, take care of the mem
 * @author Steve Hoffmann 
 *   
 */
 
mappingset_t*
bl_greedypairMappingSets(mappingset_t *l, mappingset_t *r) {
  unsigned int i, j;
  long long int d1, d2=4000000000;
  unsigned char found = 0;
  mapping_t *m1, *m2, *m3, *m4;
  mappingset_t *s = NULL;

  for(i=0; i < l->n; i++) {
    for(j=0; j < r->n; j++) {
      m1 = &l->elem[i];
      m2 = &r->elem[j];
      d1 = bl_distMapping(m2, m1);

      if(d1 != LLONG_MAX && (!found || d1 <= d2)) {
        found = 1; 
        m3 = m1;
        m4 = m2;
        d2 = d1;
      }
    }
  }

  if(found) {
    s=ALLOCMEMORY(NULL, NULL, mappingset_t, 1);
    s->elem = bl_concatMapping(m3, m4);
    s->n = 1;
  }

  return s;
}


/*------------------------------ bl_myersMCSA --------------------------------
 *    
 * @brief do myers bitvector to fill MCSA
 * @author Steve Hoffmann 
 *   
 */
 

unsigned int
bl_myersMCSA(MultiCharSeqAlignment *mcsa, MultiCharSeq *mseq, 
    unsigned int *enctab, bitvector *peq, bitvector *D, unsigned int maxE) {

  PairSint mb;

  myersbitmatrix (NULL, mcsa->query, mcsa->qrylen, 
      mcsa->refseq, mcsa->reflen, mseq->map, 
      mseq->mapsize, enctab, mcsa->qrylen-maxE, peq, &mb, D, mcsa->reflen);


  if(mb.a != -1 && mb.b <= maxE && mb.a < mcsa->reflen) {
    bitvectorbacktrack(mcsa->al, D, mcsa->reflen, mcsa->qrylen, 
        mb.a, mcsa->refseq, enctab, peq);

    Uint newedist = getEdist(mcsa->al);
    if(newedist != mb.b) { 
   
#ifdef MYERSALIGNDEBUG
      DBG("myers align error: mb.a=%d, mb.b=%d, maxE=%d, qrylen=%d, reflen=%d", 
          mb.a, mb.b, maxE, mcsa->qrylen, mcsa->reflen);
      showAlign(mcsa->al, DBGDEVICE);
#endif
    }
    if(newedist > maxE) {
      return 0;
    }

    return 1;
  } else {
    return 0;
  }

}


/*------------------------------- bl_matealign -------------------------------
 *    
 * @brief this function assigns up to four mate alignments in upstream, 
 * downstream in fwd and rev an generates a new mapping set if a mate can not 
 * be found for a sequence it is added as a singleton to the new mapping set
 * @author Steve Hoffmann 
 *   
 */

mappingset_t *
bl_matealign(mappingset_t *l, MultiCharSeq *mseq, char **seqs, char **qual, 
    char *qname, unsigned int m, unsigned int *enctab, bitvector *D, 
    unsigned int insertsize, char ismate, Suffixarray *arr, segemehl_t *nfo){

  unsigned int i, j, k, p, q, loff, ret,len;
  MultiCharSeqAlignment mcsa, *copy;
  mappingset_t *s;
  bitvector *peq[2];
  Uint strandseq[2];
  //char matefound =0 ;
  char skip;
  PairUint res[2];
  Uint *spos[2], val;
  Uint npos[2];

  unsigned int maxedist = getMaximumAlignmentEdist(m, nfo->accuracy);
  //get extended length of the sequence (includes edist on both sides)
  unsigned int ext = getExtendedAlignmentLength(m, nfo->accuracy);
  //get the offset for each side
  //unsigned int off = floor(((double)ext-m)/2.0);
  //add the insertsize to the extended length
  ext += insertsize;

  s = ALLOCMEMORY(NULL, NULL, mappingset_t, 1);
  bl_initMappingSet(s);  

  for(j=0; j < 2; j++){ 
    npos[j] = 0;
    res[j] = searchSuffix(NULL, arr, seqs[j], m);
    if(res[j].a <= res[j].b && res[j].b - res[j].a < 20) { 
      npos[j] = res[j].b - res[j].a + 1;
      spos[j] = ALLOCMEMORY(NULL, NULL, Uint, npos[j]);
      for(k=0, i=res[j].a; i <= res[j].b; i++, k++) {
        spos[j][k] = arr->suftab[i];
      }
      qsort(spos[j], npos[j], sizeof(Uint), cmp_Uint_qsort); 
    }
  }

  peq[0] = getpeq(NULL, seqs[0], m, mseq->map, mseq->mapsize, enctab);
  peq[1] = getpeq(NULL, seqs[1], m, mseq->map, mseq->mapsize, enctab);

  //for each existing mapping in set
  for(i=0; i < l->n; i++) {
    //and only for unpaired mappings
    if(!bl_isPairedMapping(&l->elem[i])) {  

      //reset flag
      //matefound = 0;

      //align fwd and reverse mate seq
      if(bl_getMapFragStrand(&l->elem[i].f[0])) {
        strandseq[0] = 0;
        strandseq[1] = 1;
      } else { 
        strandseq[0] = 1;
        strandseq[1] = 0;
      } 

      for(skip=0, j=0; j < 2; j++) {
        
        //if this is a split alignment start and beginning need to be checked
        for(q=0; q < MIN(2,l->elem[i].n); q++) { 

          if(q==0) { 
            p = l->elem[i].P;
          } else {
            p = l->elem[i].Q;
          }

          loff = ext/2; //off
          len = ext;      
          ret = 0;

          //search the suffix array first for full match before the more 
          //expensive alignment
          
          if(npos[j]) { 

            for(val=0, k=0; k < npos[strandseq[j]]; k++) { 
              if(DIST(spos[strandseq[j]][k], p) < 3000) {
                val = spos[strandseq[j]][k];
              }
            }

            if(val) { 
              initMultiCharSeqAlignment(NULL, &mcsa, mseq, val, 0, 
                  m, strandseq[j], qname, seqs[strandseq[j]], qual[strandseq[j]], m);

              insertMeop(mcsa.al, Match, m);

              if(nfo->bisulfite) {
                //required to enforce C->T/G->A and exclude T->C/A->G
                reevalMultiCharSeqAlignment(&mcsa);
              }

              ret =1;
            }
          } 

          if (!ret){ 

            //align the mate
            initMultiCharSeqAlignment(NULL, &mcsa, mseq, p, loff, 
                len, strandseq[j], qname, seqs[strandseq[j]], qual[strandseq[j]], m);

            ret = bl_myersMCSA(&mcsa, mseq, enctab, peq[strandseq[j]], D, maxedist); 
          }

          //add the pair if suitable to a new set
          if(ret && !bl_mappingsetHasPos(s, mcsa.refstart+mcsa.al->voff )) {
            //copy the new fragment
            copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
            copy = bl_copyMCSA(copy, &mcsa);
            //rescue the old mapping
            bl_addMapping(s, &l->elem[i]);
            //add the new fragment
            bl_addMapFrag(&s->elem[s->n-1], copy, NULL, ismate, 0);
            //skip next strand calculation if score good enough
            if(getEdist(mcsa.al) <= (Uint)(.95*((float)len))) skip=1;
            //matefound = 1;
          }

          //wrap the mcsa, it has been copied before
          wrapMultiCharSeqAlignment(NULL, &mcsa);
          if(skip) break;
        }
      }
    }
    //add the singleton in any case
    bl_addMapping(s, &l->elem[i]);

  }


  //clean up
  for(j=0; j < 2; j++) {
    for(i=0; i < mseq->mapsize; i++) {
      FREEMEMORY(space, peq[j][i]);
    }  
    FREEMEMORY(space, peq[j]);
  }

  if(npos[0]) FREEMEMORY(NULL, spos[0]);
  if(npos[1]) FREEMEMORY(NULL, spos[1]);

  return s;
}


/*---------------------------- bl_pairMateMapping ----------------------------
 *    
 * @brief add all potential query mate pairs to a new mapping
 * @author Steve Hoffmann 
 *   
 */

  mappingset_t *
bl_pairMateMapping (mappingset_t *l, mappingset_t *r, unsigned int insertsize)
{

  unsigned int i, j, k;
  long long int d;
  mappingset_t *s;
  MultiCharSeqAlignment *copy;
  char *pairedR;//, pairedL=0;

  s = ALLOCMEMORY(NULL, NULL, mappingset_t, 1);
  bl_initMappingSet(s);

  pairedR = ALLOCMEMORY(NULL, NULL, char, r->n);
  memset(pairedR, 0, sizeof(char)*r->n);

  //iterate the first mapping set
  for(i=0; i < l->n; i++) {
    //iterate the second mapping set
    for(j=0; j < r->n; j++) {
          //l has a query mapping and r has a mate mapping
      if((bl_isQueryMapping(&l->elem[i]) && !bl_isMateMapping(&l->elem[i]) &&
            !bl_isQueryMapping(&r->elem[j]) && bl_isMateMapping(&r->elem[j])) ||
          //l has a mate mapping and r has query mapping
          (bl_isQueryMapping(&r->elem[j]) && !bl_isMateMapping(&r->elem[j]) &&
           !bl_isQueryMapping(&l->elem[i]) && bl_isMateMapping(&l->elem[i])) ){ 
        //get the distance
        d = bl_distMapping(&l->elem[i], &r->elem[j]);
        //check insert size and add pair to the new set
#ifdef PAIRMATEDEBUG
        DBG("checking insertsize btw elems: i=%d, j=%d, d=%lld\n", i, j, d);
#endif
        if(d < insertsize) {
          pairedR[j] = 1;
          //pairedL = 1;
          //add first mapping
          bl_addMapping(s, &l->elem[i]);
          //add the fragments of the second mapping to first
          for(k=0; k < r->elem[j].n; k++) {
            //copy the new fragment
            copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
            copy = bl_copyMCSA(copy, r->elem[j].f[k].mcsa);
            bl_addMapFrag(&s->elem[s->n-1], copy, 
                r->elem[j].f[k].seed, r->elem[j].f[k].mate, 
                r->elem[j].f[k].issplit);
          }
      
          assert(s->elem[s->n-1].matestatus == 3);
        }
      }
    }
    //if the l mapping was not paired, add singleton
    //if(!pairedL) {  
    bl_addMapping(s, &l->elem[i]);
    //}
    //pairedL = 0;
  }

  //add all r mappings that were not paired
  for(j=0; j < r->n; j++) {
   // if(!pairedR[j]) {
      bl_addMapping(s, &r->elem[j]);
   // }
  }

  FREEMEMORY(space, pairedR);
  return s;
}



