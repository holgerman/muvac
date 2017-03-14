
/*
 *  match.c
 *  replacement for kdmatch
 *  it gets the seeds from kdseed.c,
 *  it uses queryalign, matealign and
 *  split align. the last three return a mapping_t
 *  data structure
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 12/21/2012 00:25:18 CET
 *  
 */

#include "../segemehl.h"
#include "karlin.h"
#include "mapfrag.h"
#include "bitvectoralg.h"
#include "sufarray.h"
#include "kdseed.h"
#include "mathematics.h"
#include "iupac.h"
#include "queryalign.h"
#include "matealign.h"
#include "splitalign.h"
#include "manout.h"
#include "vtprogressbar.h"
#include "samout.h"
#include "mappingqual.h"
#include "splitalign.h"
#include "pigeon.h"

/*--------------------------------- se_clip ----------------------------------
 *    
 * @brief clipping sequences
 * @author Steve Hoffmann 
 *   
 */

  void
se_clip (void *space, fasta_t *reads, Uint elem, segemehl_t *nfo)
{

  if(nfo->hardclip3Prime || nfo->hardclip5Prime) {
    bl_fastaHardClip(space, reads, elem, nfo->hardclip5Prime, 
        nfo->hardclip3Prime);
    if(bl_fastaHasMate(reads)) {
      bl_fastaMateHardClip(space, reads, elem, nfo->hardclip5Prime, 
          nfo->hardclip3Prime);
    } 
  } 

  if(nfo->softclip3Prime || nfo->softclip5Prime) {
    bl_fastaSoftClip(space, reads, elem, 
        nfo->softclip5Prime, nfo->softclip5PrimeLen, nfo->minclipscr5,
        nfo->softclip3Prime, nfo->softclip3PrimeLen, nfo->clipacc, nfo->polyAlen);
    if(bl_fastaHasMate(reads)) {
      bl_fastaMateSoftClip(space, reads, elem, 
          nfo->softclip5Prime, nfo->softclip5PrimeLen, nfo->minclipscr5,
          nfo->softclip3Prime, nfo->softclip3PrimeLen, nfo->clipacc, nfo->polyAlen);
    }
  }

  return ;
}


/*-------------------------------- se_convert --------------------------------
 *    
 * @brief a small helper for bisulfite conversion
 * @author Steve Hoffmann 
 *   
 */

  char **
se_convert (char **out, char *in, unsigned int len, segemehl_t *nfo, char type)
{        

  void *space = NULL;

  if (nfo->bisulfite){
    FREEMEMORY(space, out[1]);
    out[0] = ALLOCMEMORY(space, out[0], char, len+1);
    memmove(out[0], in, len+1);
    out[1] = charIUPACcomplement(space, out[0], len);
    bl_convertBisulfite(out[0], len, nfo->bisulfite, type);
    bl_convertBisulfite(out[1], len, nfo->bisulfite, type);
  }

  return out;
}


/*--------------------------------- se_jump ----------------------------------
 *    
 * @brief a helper function for jump step calculation
 * @author Steve Hoffmann 
 *   
 */

  unsigned int
se_jump (unsigned int len, segemehl_t *nfo)
{
  unsigned int jump;

  //init jumps
  if (nfo->jump == 0) {
    jump = floor(len/75) * 2;
    jump = (jump > 0) ? jump : 1;
  } else {
    jump = nfo->jump;
  }

  return jump;
}


void
match(void *space, Suffixarray *s, fasta_t *reads, 
    segemehl_t *nfo) {

  char *seqs[2], *mateseqs[2], *quals[2], *matequals[2];
  unsigned int *enctab, dim, maxL, wordno;
  unsigned int k, i, u, len, matelen=0, jump; //, setmatches;
  Uint querysplitedist =0;
  Uint matesplitedist = 0;
  mappingset_t *tmp;
  mapseedlist_t *seeds=NULL, *mateseeds=NULL;
  karlin_t stats;
  bitvector *D, *Mv, *Mh;
  matchstem_t *stems[2] = {NULL, NULL}, *matestems[2] = {NULL, NULL},
              *b0[2], *mateb0[2]; 
  int oldqscore, oldmscore, newqscore, newmscore;
  int scores[]={1,-2};
  int indel = -3;

  mappingset_t *maps, *matemaps, *savemaps=NULL, *savemaps2=NULL;


  //initialize data for stats and bv
  enctab = encodetab(nfo->seq->map, nfo->seq->mapsize);
  dim = getExtendedAlignmentLength(reads->maxlen+1000, nfo->accuracy);
  maxL = nfo->maxinsertsize;

  //add the maximum insert size if we use the myers for looking up mate
  if(bl_fastaHasMate(reads)) {
    dim += nfo->maxinsertsize;
  }

  wordno = (reads->maxlen/BITVECTOR_WORDSIZE)+1; 

  D = ALLOCMEMORY(space, NULL, bitvector, 3*(dim+1));
  Mv = &D[dim+1];
  Mh = &D[2*(dim+1)];

  for(i=0; i <= dim; i++) {
    D[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
    Mv[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
    Mh[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
  }  

  karlinunitcostpp(space, &stats.lambda, &stats.H, &stats.K);

  //iter the reads
  for(k=0; k < reads->noofseqs; k++) {
    seeds = NULL;
    mateseeds = NULL;
    //progress bar   
    if(!nfo->mute) se_updateProgressBar(k, nfo);
    //clipping
    se_clip(space, reads, k, nfo);
    //get the read len and ...
    len = bl_fastaGetSequenceLength(reads, k);

    /*attempt to split align poorly aligned reads*/
    querysplitedist = (((double)len)/100.0 * 2.0);

    //check the minium length
    if(len >= nfo->minsize) {
      //get the query (and mate) and convert for bisulfite if necessary
      seqs[0] = bl_fastaGetSequence(reads, k);
      seqs[1] = charIUPACcomplement(space, seqs[0], len);
      //get the quals and reverse them
      if(bl_fastaHasQuality(reads)) { 
        quals[0] = bl_fastaGetQuality(reads, k);
        quals[1] = ALLOCMEMORY(NULL, NULL, char, len+1);
        memmove(quals[1], quals[0], len+1);
        quals[1] = strrev(quals[1], len);
      } else {
        quals[0] = NULL;
        quals[1] = NULL;
      }
      //convert the query if necessary
      //se_convert(seqs, bl_fastaGetSequence(reads,k), len, nfo, 1);

      if (nfo->bisulfite){    
        seqs[0] = ALLOCMEMORY(space, NULL, char, len+1);
        memmove(seqs[0], bl_fastaGetSequence(reads, k), len+1);
        bl_convertBisulfite(seqs[0], len, nfo->bisulfite, 1);
        bl_convertBisulfite(seqs[1], len, nfo->bisulfite, 1);
      }

      stems[0] = NULL; stems[1] = NULL; 
      b0[0] = NULL; b0[1] = NULL;

      //init jumps 
      jump = se_jump(len, nfo);

      /* restrict search to one strand */
      for (u = 0; u < 2; u++){
        /* nfo->strand == 1 : search only on plus strand
         * => init stems[1] as empty
         * nfo->strand == 2 : search only on minus strand
         * => init stems[0] as empty
         * Note: empty but initialized stems are ignored
         * in function kdbest
         */
        if (nfo->strand == 2 - u){	 
          stems[u] = ALLOCMEMORY(space, NULL, matchstem_t, len);
          for (i = 0; i < len; i++){
            stems[u][i].branches = NULL;
            stems[u][i].noofbranches = 0;
          }
        }
      }

      char pigeonsearch = 0;
      if (nfo->bestonly && countNonMatchingChars(seqs[0], len) <= 1){
        maps = pigeon(space, s, nfo->seq, seqs, quals, len, 2, bl_fastaGetDescription(reads, k),
            enctab, D, 0);
       /* if(bl_hasQueryMapping(maps)) {
          fprintf(stdout, "found pigeon alignment\n");
        } else { 
          fprintf(stdout, "did not find pigeon alignment\n");
        }
        */  
        pigeonsearch = 1;
      }

      if(!pigeonsearch || !bl_hasQueryMapping(maps)) { 

        if(pigeonsearch) { 
        //remove the old map, a new one will be generated
        bl_removeMappingSet(maps);
        FREEMEMORY(space, maps);
        }

        if (stems[0] == NULL){
          stems[0]=kdseeds(space, s, seqs[0], len, jump, nfo->s_ext, nfo->p_mis,
              nfo->Xoff, nfo->k_p, b0[0], nfo->nosuflinks);
        }

        if (stems[1] == NULL){
          stems[1]=kdseeds(space, s, seqs[1], len, jump, nfo->s_ext, nfo->p_mis,
              nfo->Xoff, nfo->k_p, b0[1], nfo->nosuflinks);
        }

        /* convert for alignment */
        //se_convert(seqs, bl_fastaGetSequence(reads,k), len, nfo, 0);
        /* convert for alignment */
        if (nfo->bisulfite){
          FREEMEMORY(space, seqs[1]);
          memmove(seqs[0], bl_fastaGetSequence(reads, k), len+1);
          seqs[1] = charIUPACcomplement(space, seqs[0], len);
          bl_convertBisulfite(seqs[0], len, nfo->bisulfite, 0); 
          bl_convertBisulfite(seqs[1], len, nfo->bisulfite, 0); 
        }


        //filter the seeds
        seeds = bl_getGoodSeeds(stems, len, s->numofsuffixes, &stats, nfo);
        //map them
        maps = bl_seedAlign(s, nfo->seq, seqs, quals, len, 
            bl_fastaGetDescription(reads, k), seeds, nfo, enctab, D, 0);
      }
    } else {
      //minimum length criterium failed
      maps = ALLOCMEMORY(space, NULL, mappingset_t, 1);
      bl_initMappingSet(maps);
    } 

    //prepare the mates
    if(bl_fastaHasMate(reads)){
      matelen = bl_fastaGetMateLength(reads, k);
      if(matelen >= nfo->minsize) { 
        matesplitedist = (((double)matelen)/100.0 * 2.0);

        //get the mate seq
        mateseqs[0] = bl_fastaGetMate(reads, k);
        mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen);
        //get the mate qual
        if(bl_fastaHasQuality(reads)) { 
          matequals[0] = bl_fastaGetMateQuality(reads, k);
          matequals[1] = ALLOCMEMORY(NULL, NULL, char, matelen+1);
          memmove(matequals[1], matequals[0], matelen+1);
          matequals[1] = strrev(matequals[1], matelen);
        } else {
          matequals[0] = NULL;
          matequals[1] = NULL;
        }
        //convert the mate for direct alignment
        //se_convert(mateseqs, bl_fastaGetMate(reads,k), matelen, nfo, 0);
        if (nfo->bisulfite){
          mateseqs[0] = ALLOCMEMORY(space, NULL, char, matelen+1);
          memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
          bl_convertBisulfite(mateseqs[0], matelen, nfo->bisulfite, 0);
          bl_convertBisulfite(mateseqs[1], matelen, nfo->bisulfite, 0);
        }   

      }

      matestems[0] = NULL; matestems[1] = NULL;
      mateb0[0] = NULL; mateb0[1] = NULL;
    }

    //now go check the mates
    if (bl_fastaHasMate(reads) && matelen >= nfo->minsize) {
      //fprintf(stderr, "now mapping the mate\n");
      tmp = bl_matealign(maps, nfo->seq, mateseqs, matequals, 
          bl_fastaGetMateDescription(reads, k),
          matelen, enctab, D, maxL, 1, s, nfo);
      //tmp = pigeonmate(space, maps, s, nfo->seq, mateseqs, matequals, 
      //    matelen, 3, bl_fastaGetMateDescription(reads, k), enctab, D, 1); 

      
      bl_removeMappingSet(maps);
      FREEMEMORY(NULL, maps);
      maps = tmp;
      //did the attempt fail?
      //relative fail


      if(nfo->split && bl_mappingsetHasPaired(maps) && ( 
            !bl_hasQueryMappingMaxEdist(maps, querysplitedist) ||
            !bl_hasMateMappingMaxEdist(maps, matesplitedist)
            )) {
        //          fprintf(stdout, "saving am map\n");
        savemaps = maps;
        maps = ALLOCMEMORY(space, NULL, mappingset_t, 1);
        bl_initMappingSet(maps);
      } else {
        //          fprintf(stdout, "no map saved\n");
        savemaps = NULL;
      }

      if(!bl_mappingsetHasPaired(maps)) {
        //fprintf(stdout, "mapping mate with matealign has failed\n");
        //set up the heuristics for jumping seeds
        jump = se_jump(matelen, nfo);
        /* convert for mate seed search */
        //se_convert(mateseqs, bl_fastaGetMate(reads,k), matelen, nfo, 1);
        if (nfo->bisulfite){
          FREEMEMORY(space, mateseqs[1]);
          memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
          mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen);
          bl_convertBisulfite(mateseqs[0], matelen, nfo->bisulfite, 1);
          bl_convertBisulfite(mateseqs[1], matelen, nfo->bisulfite, 1);
        } 

        matestems[0] = NULL; matestems[1] = NULL;
        mateb0[0] = NULL; mateb0[1] = NULL;

        /* restrict search to one strand (bisulfite only)*/
        for (u = 0; u < 2; u++){
          /* nfo->strand == 1 : search only on plus strand
           * => search for mate only on minus strand
           * => init stems[0] as empty
           * nfo->strand == 2 : search only on minus strand
           * => search for mate only on plus strand
           * => init stems[1] as empty
           * Note: empty but initialized stems are ignored
           * in function kdbest
           */
          if (nfo->strand == u + 1){
            matestems[u] = ALLOCMEMORY(space, NULL, matchstem_t, matelen);
            for (i = 0; i < matelen; i++){
              matestems[u][i].branches = NULL;
              matestems[u][i].noofbranches = 0;
            }
          }
        }

        if (nfo->bestonly && countNonMatchingChars(mateseqs[0], matelen) <= nfo->k_p){
            kdbest(space, s, mateseqs, matelen, nfo->s_ext, nfo->p_mis,
              nfo->Xoff, nfo->k_p, matestems, mateb0);

        }

        if (matestems[0] == NULL){
          matestems[0]=kdseeds(space, s, mateseqs[0], matelen, 
              jump, nfo->s_ext, nfo->p_mis,
              nfo->Xoff, nfo->k_p, mateb0[0], nfo->nosuflinks);
        }
        
        if (matestems[1] == NULL){
          matestems[1]=kdseeds(space, s, mateseqs[1], matelen, 
              jump, nfo->s_ext, nfo->p_mis,
              nfo->Xoff, nfo->k_p, mateb0[1], nfo->nosuflinks);
        }

        /* convert query for mate alignment */
        //se_convert(mateseqs, bl_fastaGetMate(reads,k), matelen, nfo, 0);
        /* convert for mate alignment */
        if (nfo->bisulfite){
          FREEMEMORY(space, mateseqs[1]);
          memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
          mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen);
          bl_convertBisulfite(mateseqs[0], matelen, nfo->bisulfite, 0);
          bl_convertBisulfite(mateseqs[1], matelen, nfo->bisulfite, 0);
        }   

        //filter the seeds
        mateseeds = bl_getGoodSeeds(matestems, matelen, s->numofsuffixes, &stats, nfo);
        //align mate to the seeds
        matemaps = bl_seedAlign(s, nfo->seq, mateseqs, matequals, matelen, 
            bl_fastaGetMateDescription(reads, k),
            mateseeds, nfo, enctab, D, 1);
        //now go check the query (bisulfite is already in alignment mode)
        if(len >= nfo->minsize) { 
          tmp = bl_matealign(matemaps, nfo->seq, seqs, quals, 
              bl_fastaGetDescription(reads, k),
              len, enctab, D, maxL, 0, s, nfo);
          bl_removeMappingSet(matemaps);
          FREEMEMORY(NULL, matemaps);
          matemaps = tmp;
        }

        //attempt failed again?
        if(!bl_mappingsetHasPaired(matemaps)) {
          //  fprintf(stdout, "pairing alignment\n");
          //pair mate maps or copy singletons to matemap
          tmp = bl_pairMateMapping(maps, matemaps, 3000000000);
          bl_removeMappingSet(maps);
          bl_removeMappingSet(matemaps);
          FREEMEMORY(NULL, maps);
          FREEMEMORY(NULL, matemaps);
          maps = tmp;

          //            bl_dumpMappingSet(stdout, maps);


        } else {

          //if the mate first attempt resulted in a mate match
          //accepted it if split is turned off
          if(!nfo->split || 
              ( bl_hasQueryMappingMaxEdist(matemaps, querysplitedist) 
                && bl_hasMateMappingMaxEdist(matemaps, matesplitedist))
            ) { 

            bl_removeMappingSet(maps);
            FREEMEMORY(NULL, maps);
            maps = matemaps;
          } else {
            //otherwise check if it is a better alignment
            //than the one stored in the query first attempt
            if(savemaps && 
                bl_getMappingMaxScore(savemaps, scores, indel) > 
                bl_getMappingMaxScore(matemaps, scores, indel)) {
              //if this is not the case remove the mate maps
              bl_removeMappingSet(matemaps);
              FREEMEMORY(NULL, matemaps);

            } else {
              //if this is the case remove savemaps
              //and store the matempas in savemaps
              if(savemaps) { 
                //                  fprintf(stdout, "overwriting the following savemap\n");
                //                  bl_dumpMappingSet(stdout, savemaps);
                bl_removeMappingSet(savemaps);
                FREEMEMORY(NULL, savemaps);
              }
              savemaps = matemaps;
            }
          }
        } 
      } 
    }

    bl_sortMappingSetByScore(maps, scores, indel);

    mappingset_t* querysplitmap;
    mappingset_t* matesplitmap;

    if(nfo->split && savemaps) {
      //        fprintf(stdout, "following alignments have been put away\n");
      //        bl_dumpMappingSet(stdout, savemaps);
      bl_sortMappingSetByScore(savemaps, scores, indel);
    }

    /*attempt to split align unmapped mate 1 reads or
     * previously removed potentially suboptimal alignments*/     
    if(nfo->split && stems[0] && len >= nfo->minsize && (!bl_hasQueryMapping(maps) ||
          !bl_hasQueryMappingMaxEdist(maps,querysplitedist)) ) {
      savemaps2 = NULL;

      //init a split map
      querysplitmap = ALLOCMEMORY(space, NULL, mappingset_t, 1);
      bl_initMappingSet(querysplitmap);
      //caluclate split map
      querysplitmap = bl_splitAlign(space, s,  nfo->seq, querysplitmap,
          bl_fastaGetDescription(reads, k), 
          stems, seqs, quals, len, &stats, enctab, D, 0, nfo);

      //       fprintf(stdout, "querysplitmap result\n");
      //        bl_dumpMappingSet(stdout, querysplitmap);

      //if the splitmap run was successful and
      //if there are query maps with more than querysplit edist
      //save the old set, remove all affected query maps
      //and replace them by the new split
      if(bl_hasQueryMapping(maps) && bl_hasQueryMapping(querysplitmap) 
          && bl_hasQueryMappingMaxEdist(maps, querysplitedist))  {

        savemaps2 = maps;

        //          fprintf(stdout, "merging 0\n");
        maps = bl_copyMappingSet(savemaps2);
        bl_removeBadMatesEdistQM(maps, querysplitedist, 0);
        tmp = bl_pairMateMapping(maps, querysplitmap, 3000000000);
        bl_removeMappingSet(maps);
        FREEMEMORY(NULL, maps);
        maps = tmp;
      } else {
        //          fprintf(stdout, "merging 1\n");
        tmp = bl_pairMateMapping(maps, querysplitmap, 3000000000);
        bl_removeMappingSet(maps);
        FREEMEMORY(NULL, maps);
        maps = tmp;
      }
      //        fprintf(stdout, "result 0\n");
      //        bl_dumpMappingSet(stdout, maps);

      bl_removeMappingSet(querysplitmap);
      FREEMEMORY(NULL, querysplitmap);

      bl_sortMappingSetByScore(maps, scores, indel);
      //     fprintf(stdout, "following alignment was obtained\n"); 
      //     bl_dumpMappingSet(stdout, maps);

      if(savemaps2) { 

        bl_getMappingScoreQM (&savemaps2->elem[0], scores, indel, &oldqscore, &oldmscore);
        if(maps->n > 0)
          bl_getMappingScoreQM (&maps->elem[0], scores, indel, &newqscore, &newmscore);

        if(maps->n == 0 || (oldmscore + oldqscore > newqscore + newmscore)) {
          //        fprintf(stdout, "mapping correction not accepted\n");
          bl_removeMappingSet(maps);
          FREEMEMORY(space, maps);
          maps = savemaps2;
        } else {
          //        fprintf(stdout, "mapping correction accepted\n");
          bl_removeMappingSet(savemaps2);
          FREEMEMORY(space, savemaps2);
        }
      }
    }

    if(nfo->split && len >= nfo->minsize) { 
      bl_fixSplitAlign(s, maps, nfo->seq, bl_fastaGetDescription(reads, k), 
          seqs, quals, len, 0);

      if(bl_fastaHasMate(reads) && matelen >= nfo->minsize) {   
        bl_crosscorrection(s, maps, nfo->seq, 
            bl_fastaGetMateDescription(reads, k), mateseqs, matequals, matelen, len, matesplitedist, 1);
      }
    }

    //matestems are not necessarily available if there was a suffix array search
    //if query alignment was successful.
    if(nfo->split && bl_fastaHasMate(reads) && matelen >= nfo->minsize 
        && matestems[0] && (!bl_hasMateMapping(maps) || 
          !bl_hasMateMappingMaxEdist(maps,matesplitedist) ||
          !bl_getMappingIsConsecutive(&maps->elem[0]) || bl_getMappingRange(&maps->elem[0]) > 10000)) {  
      //fprintf(stderr, "splitmap 2\n");
      savemaps2 = NULL;

      //init a split map
      matesplitmap = ALLOCMEMORY(space, NULL, mappingset_t, 1);
      bl_initMappingSet(matesplitmap);
      //caluclate split map
      matesplitmap = bl_splitAlign(space, s,  nfo->seq, matesplitmap,
          bl_fastaGetMateDescription(reads, k), 
          matestems, mateseqs, matequals, matelen, &stats, enctab, D, 1, nfo);

      //        fprintf(stdout, "matesplitmap result\n");
      //        bl_dumpMappingSet(stdout, matesplitmap);


      if(bl_hasMateMapping(maps) && bl_hasMateMapping(matesplitmap) 
          && bl_hasMateMappingMaxEdist(maps, matesplitedist))  {

        //          fprintf(stdout, "merging 2\n");
        savemaps2 = maps;
        maps = bl_copyMappingSet(savemaps2);
        //          fprintf(stdout, "before remove bad mates\n");
        //          bl_dumpMappingSet(stdout, maps);

        bl_removeBadMatesEdistQM(maps, matesplitedist, 1);

        //          fprintf(stdout, "after removing bad mates\n");
        //          bl_dumpMappingSet(stdout, maps);

        tmp = bl_pairMateMapping(maps, matesplitmap, 3000000000);
        bl_removeMappingSet(maps);
        FREEMEMORY(NULL, maps);
        maps = tmp;
      } else {
        //          fprintf(stdout, "merging 3\n");
        //          bl_dumpMappingSet(stdout, matesplitmap);
        //          bl_dumpMappingSet(stdout, maps);
        tmp = bl_pairMateMapping(maps, matesplitmap, 3000000000);
        bl_removeMappingSet(maps);
        FREEMEMORY(NULL, maps);
        maps = tmp;
      }

      //        fprintf(stdout, "result 1\n");
      //        bl_dumpMappingSet(stdout, maps);
      bl_removeMappingSet(matesplitmap);
      FREEMEMORY(NULL, matesplitmap);
      bl_sortMappingSetByScore(maps, scores, indel);
      //     fprintf(stdout, "following alignment was obtained\n"); 
      //     bl_dumpMappingSet(stdout, maps);

      if(savemaps2) { 

        bl_getMappingScoreQM (&savemaps2->elem[0], scores, indel, &oldqscore, &oldmscore);
        if(maps->n > 0)
          bl_getMappingScoreQM (&maps->elem[0], scores, indel, &newqscore, &newmscore);

        if(maps->n ==0 || (oldmscore + oldqscore > newqscore + newmscore)) {
          //        fprintf(stdout, "mapping correction not accepted\n");
          bl_removeMappingSet(maps);
          FREEMEMORY(space, maps);
          maps = savemaps2;
        } else {
          //        fprintf(stdout, "mapping correction accepted\n");
          bl_removeMappingSet(savemaps2);
          FREEMEMORY(space, savemaps2);
        }
      }
    }

    /*check alignment ends for splits and gaps between splits*/
    if(nfo->split && bl_fastaHasMate(reads) && matelen >= nfo->minsize) { 

      bl_fixSplitAlign(s, maps, nfo->seq, bl_fastaGetMateDescription(reads, k), 
          mateseqs, matequals, matelen, 1);
      if(len >= nfo->minsize) { 
        bl_crosscorrection(s, maps, nfo->seq, 
            bl_fastaGetDescription(reads, k), seqs, quals, len, matelen, querysplitedist, 0);
      }
    }



    /*restore old alignment if the split attempt was not successful*/
    if(nfo->split && savemaps)  {

      //      fprintf(stdout, "savemaps is\n");
      //      bl_dumpMappingSet(stdout, savemaps);

      bl_getMappingScoreQM (&savemaps->elem[0], scores, indel, &oldqscore, &oldmscore);
      if(maps->n > 0)
        bl_getMappingScoreQM (&maps->elem[0], scores, indel, &newqscore, &newmscore);

      if(maps->n == 0 || (oldmscore + oldqscore > newqscore + newmscore)) {
        bl_removeMappingSet(maps);
        FREEMEMORY(space, maps);
        maps = savemaps;
      } else {
        bl_removeMappingSet(savemaps);
        FREEMEMORY(space, savemaps);
      } 
    }

    //      fprintf(stdout, "no savemaps\n");
    savemaps = NULL;
    mappingset_t* newmaps = sam_mappingJoinFrags(maps, nfo);


    /*GET THE BEST SEEDS (FOR UNMAPPED READS)*/ 
    mapseed_t *bestseed = NULL;
    if(seeds && seeds->n) { 
      bestseed = bl_getMapSeedListBest(seeds);
      bl_getMapSeedAdditionalInformation (bestseed, s, nfo->seq);
    }

    mapseed_t *bestmateseed = NULL;
    if(bl_fastaHasMate(reads) && mateseeds && mateseeds->n) { 
      bestmateseed = bl_getMapSeedListBest(mateseeds);
      bl_getMapSeedAdditionalInformation (bestmateseed, s, nfo->seq);
    }

    /*OUTPUT*/
    se_output(newmaps, reads, k, bestseed, bestmateseed, s, &stats, nfo);


    /*LOOP CLEANUP BEGIN*/ 
    if(nfo->bisulfite) { 
      FREEMEMORY(space, seqs[0]);
      if(bl_fastaHasMate(reads)) { 
        FREEMEMORY(space, mateseqs[0]);
      }
    } 
    FREEMEMORY(space, seqs[1]);
    if(bl_fastaHasQuality(reads)) { 
      FREEMEMORY(space, quals[1]);
    }

    if(bl_fastaHasMate(reads)) {
      FREEMEMORY(space, mateseqs[1]);
      if(bl_fastaHasQuality(reads)) { 
        FREEMEMORY(space, matequals[1]);
      }
    }

    if(seeds) { 
      bl_wrapSeedList(seeds);
      FREEMEMORY(space, seeds);
    }

    if(mateseeds) {
      bl_wrapSeedList(mateseeds);
      FREEMEMORY(space, mateseeds);
    }

    if(stems[0]) { 
      bl_kdMatchstemDestruct(space, stems[0], len);
      bl_kdMatchstemDestruct(space, stems[1], len);
      stems[0] = NULL;
      stems[1] = NULL;
    }

    bl_removeMappingSet(maps);
    FREEMEMORY(space, maps);
    bl_removeMappingSet(newmaps);
    FREEMEMORY(space, newmaps);

    if(matestems[0]) {
      bl_kdMatchstemDestruct(space, matestems[0], matelen);
      bl_kdMatchstemDestruct(space, matestems[1], matelen);
      matestems[0] = NULL;
      matestems[1] = NULL;
    }
    /*LOOP CLEANUP END*/

    //    } else {
    //      fprintf(stdout, "minimum length\n");
    //    }

    //if there are no valid query alignments 
    //if(!bl_mappingSetHasQueryMatches(nmaps)) {
    //...do the split

    //if there are mates
    //if (bl_mappingSetHasMateMatches(maps)) {
    //attach them

    //}
    //}

    //if there are no valid mate alignments 
    //if(bl_fastaHasMate(reads) && !bl_mappingSetHasMateMatches(nmaps)) {
    //...do the split

    //if there are queries
    // if(bl_mappingHasQueryMatches(maps)) {
    //attach them 
    //}
    //}

}

wrapBitmatrix(space, D, 3*(dim+1));
FREEMEMORY(space, D);
FREEMEMORY(space, enctab);
}
