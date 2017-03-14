
/*
 *  deleteme.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/23/14 13:38:54 CEST
 *  
 */



/*--------------------------- se_kdAlignSplitChain ---------------------------
 *    
 * @brief align a chain of fragments using local multi spliced alignment
 * @author Steve Hoffmann 
 *   
 */

mapping_t*
se_kdAlignSplitChain (void *space, branchChain_t *chains, Uint noofchains,
    Suffixarray *arr, MultiCharSeq *seq, char *querydesc, matchstem_t **stems,
    char **seqs, Uint qrylen, int *scores, int indel, int transition, 
    spliceevents_t *events, segemehl_t *nfo) {

  Uint i, j, start, k, edist=0, floff = 0,
    ulen=0, vlen=0, maxedist, totalcover = 0, totaledist = 0,
     beststart=0;
  gmatchlist_t *list; 
  unsigned char purge = 0;
  int score = 0;
  branchChain_t *cur;
  MultiCharSeqAlignment *a;
  mapping_t *map; int **bd;

  maxedist = qrylen - ceil((nfo->accuracy * qrylen)/100);
  list = bl_gmatchlistInit(space, maxedist, 0);
  qsort(chains, noofchains, sizeof(branchChain_t), cmp_chainscores);
  
  /*DEBUG SHOW CHAINS*/
  //showChains(chains, noofchains, arr, stdout, seqs[1], qrylen); 
  cur = &chains[0];

  if (noofchains > 1) { 
    //fprintf(stdout, "nooffrags: %d; scr1:%d scr:2:%d\n", cur->nooffragments, chains[0].score, chains[1].score);
    if (((double)chains[0].score)*nfo->chainscorescale > chains[1].score) {
      //fprintf(stdout, "accept %f (%f)\n", ((double)chains[0].score)*nfo->chainscorescale, nfo->chainscorescale);
    } else {
      //fprintf(stdout, "reject %f (%f)\n", ((double)chains[0].score)*nfo->chainscorescale, nfo->chainscorescale);
      return list;
    }
  }

  if(cur->nooffragments <= 1) return list; 

  //minimize the distance of fragments on reference
  bd = bl_splitMinimizeElemDist (cur, &beststart, seq, arr);
  //init the multi char seq alignments for all fragments individually
  a = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, cur->nooffragments);

  for(i=0; i < cur->nooffragments; i++) { 
    if(cur->f[i]->strand) {
       floff = qrylen - cur->f[i]->end + maxedist;
    } else { 
       floff = cur->f[i]->start + maxedist;
    }

    start = MAX(0, (Lint) arr->suftab[bd[beststart][i]]-20);
    initMultiCharSeqAlignment(space, &a[i], seq, start,
        floff, qrylen+2*(maxedist+1+20), cur->f[i]->strand,
        querydesc, seqs[cur->f[i]->strand], qrylen);    
   // fprintf(stdout, "%s\n qrylen:%d, fragment:%d, start:%d, strand:%d, curstart:%d, curend:%d, maxedist:%d mapping to [%d,%d]\n", querydesc, qrylen, i, start, strands[i],  cur->f[i]->start, cur->f[i]->end, maxedist, start-floff, start-floff+qrylen+2*(maxedist+1+20));
  }

  //do the alignment
  a = bl_splitLocalMultiSpliceAlign(a, cur->noofragments, seqs, qrylen, seq, arr, bd, beststart);
  
  //enlist the fragments
  for(purge=0, i=0; i< cur->nooffragments; i++) {
   
    if(getUalignlen(a[i].al) >= nfo->minfragmentalignlen && 
       getValignlen(a[i].al) >= nfo->minfragmentalignlen &&
       getAlignScore(a[i].al, scores, indel) >= nfo->minfragmentalignscore) { 
   
      copy = ALLOCMEMORY(NULL, NULL, MultiCharSeqAlignment, 1);
      copy = bl_copyMCSA(&copy, &a[i]);
      map = bl_addMapFrag(map, copy, NULL, mate);
    } else {
      purge = 1;
    }
  }
  //calcuate total coverage of the query ...
  totalcover = bl_getMappingUlen(map)*100;
  totalcover /= qrylen;
  //and remove the mapping if thresholds are underscored
  if(bl_getMappingScore(map, scores, indel) < nfo->minsplicedalignscore ||
     totalcover < nfo->minsplicedaligncover && 
     map->n <= 1 || purge) {
    bl_removeMapping(map);
  }

  for (i=0; i < r-l+1; i++) {
    FREEMEMORY(space, bd[i]);
  }
  FREEMEMORY(space, a);  
  FREEMEMORY(space, bd);
 
  return map;
}

/*---------------------- bl_splitLocalMultiSpliceAlign -----------------------
 *    
 * @brief align a split chain
 * @author Steve Hoffmann 
 *   
 */
 

MultiCharSeqAlignment *
bl_splitLocalMultiSpliceAlign(MultiCharSeqAlignment *a, unsigned int n, char **seqs, 
    unsigned int qrylen, MultiCharSeq *seq, Suffixarray *arr, unsigned int **bd, int beststart) {


  //TODO: unify scoring
  int indel = -2;
  int transition = -10;
  int scores[]={1, -2};

  void *space = NULL;
  char **refseqs;
  unsigned int *strands, sub_start, sub_end, j;
  Alignment **aligns;
    unsigned int *reflens, i;
  int *M, **lmr, **lmv, **lmc;

  reflens = ALLOCMEMORY(space, NULL, Uint, n);
  strands = ALLOCMEMORY(space, NULL, Uint, n);
  starts = ALLOCMEMORY(space, NULL, Uint, n); 
  ends = ALLOCMEMORY(space, NULL, Uint, n);
  tstarts = ALLOCMEMORY(space, NULL, Uint, n); 
  tends = ALLOCMEMORY(space, NULL, Uint, n);
  lengths = ALLOCMEMORY(space, NULL, Uint, n);
  refseqs = ALLOCMEMORY(space, NULL, char*, n);
  aligns =  ALLOCMEMORY(space, NULL, Alignment*, n);
  diag = ALLOCMEMORY(space, NULL, PairUint, n);
  
  for(i=0; i < n; i++) {
    aligns[i]  = a[i].al;
    refseqs[i] = a[i].refseq;
    reflens[i] = a[i].reflen; 
    strands[i] = a[i].strand;
    lengths[i] = a[i].qrylen;
    tstarts[i] = a[i].qrystart;
    tends[i] = tstarts[i] + lengths[i]-1;

    if(strands[i]==0) { 
      starts[i] = a[i].qrystart;
      ends[i] = starts[i]+lengths[i]-1;
    } else {
      starts[i] = qrylen - (a[i].qrystart + lengths[i]);
      assert(qrylen >= a[i].qrystart+lengths[i]);
      ends[i] = starts[i]+lengths[i]-1;
      assert(ends[i] <=  qrylen);
    }

    if(strands[i]==0) { 
    diag[i].b = a[i].floff;
    } else {
       if(a[i].reflen >= MIN(maxmargin, maxedist+margin) + uroff + (cur->f[i]->end - cur->f[i]->start) - 1){ 
         diag[i].b = a[i].reflen - MIN(maxmargin, maxedist+margin) - uroff - (cur->f[i]->end - cur->f[i]->start) + 1;
       } else {
         diag[i].b = 0;
       }
    } 
    
    diag[i].a = cur->f[i]->start - a[i].qrystart;
  }

  //do the alignment

  M = localmultisplicedmatrixopt(space, seqs[0], seqs[1], qrylen, lengths,
      refseqs, reflens, strands, starts, ends, tstarts, tends, cur->nooffragments, indel, transition,
      constscr, scores, &lmv, &lmr, &lmc, &bestscr, &K, diag);


  if(M == NULL) {
    fprintf(stderr, "empty matrix returned for seqs: '%s'/'%s' (%d)\n", 
        seqs[0], seqs[1], qrylen);

    for(i=0; i < n; i++) {

      getMultiCharSeqIdxBounds(seq, a[i].subidx, &sub_start, &sub_end);
      fprintf(stderr, "fragment %d: %d in %d[%d,%d] '", 
          i, arr->suftab[bd[beststart][i]], a[i].subidx, sub_start, sub_end);
      for(j=0; j< qrylen; j++) fprintf(stderr, "%c", refseqs[i][j]);
      fprintf(stderr, "'(%d) strand:%d\n", reflens[i], strands[i]);
    }
    return NULL;
  }

#ifdef DEBUGKBAND 
  B = 
#endif 

  localmultisplicedtraceback(space, M, seqs[0], seqs[1], qrylen, 
      refseqs, reflens, strands, n, indel, transition,
      constscr, scores, aligns, lmv, lmr, lmc);
  
  for(i=0; i < n; i++) {
    FREEMEMORY(space, lmv[i]);
    FREEMEMORY(space, lmr[i]);
    FREEMEMORY(space, lmc[i]);
  }
  
  FREEMEMORY(space, reflens);
  FREEMEMORY(space, refseqs);
  FREEMEMORY(space, strands);
  FREEMEMORY(space, starts);
  FREEMEMORY(space, ends);
  FREEMEMORY(space, tstarts);
  FREEMEMORY(space, tends);
  FREEMEMORY(space, lengths);

  FREEMEMORY(space, M);
  FREEMEMORY(space, lmv);
  FREEMEMORY(space, lmr);
  FREEMEMORY(space, lmc);
  FREEMEMORY(space, aligns);

  return a;
}


