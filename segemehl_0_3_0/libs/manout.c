
/*
 *  manout.c
 *  the new output manager. simply takes a mapping set, determines the
 *  mapping policy and selects the right output functions
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 26.12.2012 02:20:31 CET
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include "mapfrag.h"
#include "samout.h"
#include "biofiles.h"
#include "segemehl.h"
#include "multicharseq.h"
#include "manoutformats.h"
#include "version.h"
#include <inttypes.h>
#include "mappingqual.h"
#include "info.h"
#include "junctions.h"

/*-------------------------- se_printUnmappedFastx ---------------------------
 *    
 * @brief print unmapped fastx
 * @author Steve Hoffmann 
 *   
 */

  void
se_printUnmappedFastx(FILE* dev, char *desc, char *seq, char *qual, mapseed_t *seed)
{

  if(qual) {
    if(seed)
      fprintf(dev, "@%s ef:%d;if:%d %" PRIu64 ":%" PRIu64 " %" PRIu64 ":%" PRIu64 ":%d\n%s\n+%s\n%s\n",		
          desc, 
          seed->maxevalue, seed->maxinterval,
          seed->u, seed->seedlen,
          seed->refidx, seed->refpos, seed->rc, 
          seq, desc, qual);
    else

      fprintf(dev, "@%s ef:0;if:0 0:0 0:0:0\n%s\n+%s\n%s\n",		
          desc, 
          seq, desc, qual);
  } else {
    if(seed)
      fprintf(dev, ">%s ef:%d;if:%d %" PRIu64 ":%" PRIu64 " %" PRIu64 ":%" PRIu64 ":%d\n%s\n", 
          desc, 
          seed->maxevalue, seed->maxinterval,
          seed->u, seed->seedlen,
          seed->refidx, seed->refpos, seed->rc,
          seq); 
    else

      fprintf(dev, ">%s ef:0;if:0 0:0 0:0:0\n%s\n",		
          desc, 
          seq);
  }

  return ;
}

/*-------------------------------- se_output ---------------------------------
 *    
 * @brief manage the output
 * @author Steve Hoffmann 
 *   
 */


void
se_output(mappingset_t *s, fasta_t * reads, unsigned int k, 
    mapseed_t *seed, mapseed_t *mateseed, Suffixarray *arr, karlin_t* stats, segemehl_t *nfo) {

  unsigned int i, noofqueries=0, noofmates=0;
  unsigned querycnt=0, matecnt=0;

  samlist_t *samlist;
  MultiCharSeq *mseq=nfo->seq;

  //remove reads with a bad accuracy
  bl_removeBadMatesAcc(s, bl_fastaGetSequenceLength(reads, k), 
      bl_fastaGetMateLength(reads, k), nfo->accuracy);

  //remove reads that are heavily softclipped
  bl_removeBadMatesCov(s, bl_fastaGetSequenceLength(reads, k),
      bl_fastaGetMateLength(reads, k), nfo->minsplicedaligncover);


  if(bl_fastaHasMate(reads)) {
    if(bl_mappingsetHasPaired(s)) {
      //proper pairs available
      if(nfo->bestonly) { 
        bl_removeSuboptimalMapping(s, nfo->scores, nfo->indel);
      }
      //still paired?
      if(bl_mappingsetHasPaired(s)) { 
        bl_removeUnpairedMapping(s);
      }
    } else { 
      //no proper pairs available
      if(nfo->bestonly) { 
        //eliminate mates and queries individually
        //fprintf(stdout, "removing mates individually\n");
        bl_removeSuboptimalMappingQM(s, nfo->scores, nfo->indel, 0);
        bl_removeSuboptimalMappingQM(s, nfo->scores, nfo->indel, 1);
      }

      //if there are only two mappings left and it is one query and one mate
      //we are allowed to merge them into one
      bl_mergeMappings(s);
    }
  } else {
    //remove suboptimal if best only is chosen
    if(nfo->bestonly) { 
      bl_removeSuboptimalMapping(s, nfo->scores, nfo->indel);
    }
  }

  bl_countMultipleMappings(s, &noofqueries, &noofmates);

  /*CALCULATE MAPPING QUAL IF NO VALUE IS SET*/
  if(!nfo->mappingqual) { 
    if(bl_fastaHasQuality(reads)) { 
      type3mappingset(s,bl_fastaGetDescription(reads, k));
    }
    longestmatchqual(s, arr->numofsuffixes, stats);
  }

  //Write BED files for splits
  if(nfo->split) { 
    bl_getSpliceJunctionsFromMappingSet(s, mseq, bl_fastaGetDescription(reads,k), nfo);
  }

  //in case there is one uniq alignment of one mate and the other is not mapped
  //RNEXT is filled with the artificial position of the unmapped mate
#ifdef SORTEDUNMAPPED
  char unmappedtornext = (bl_fastaHasMate(reads) && s->n==1 && (!bl_hasQueryMapping(s) || !bl_hasMateMapping(s)));
#else
  char unmappedtornext = 0 
#endif

  for(i=0; i < s->n; i++) {

    samlist = sam_getSamList(reads, k, &s->elem[i], mseq, (s->n > 1), noofqueries, noofmates, querycnt, matecnt, unmappedtornext, nfo);
    sam_printSamlist(samlist, &s->elem[i], nfo);

    querycnt += (bl_isQueryMapping(&s->elem[i])) ? 1 : 0;
    matecnt += (bl_isMateMapping(&s->elem[i])) ? 1 : 0;
   
    //use samlist information to store mate info in empty aligns - if there is a uniq hit
    if(s->n == 1) { 
      if(bl_fastaHasMate(reads)&& !bl_hasQueryMapping(s)) {
        char *seq  = bl_fastaGetSequenceNoClip(reads, k);
        char *qual = bl_fastaGetQualityNoClip(reads, k);
        char *desc = bl_fastaGetDescription(reads, k);
#ifdef CHRISTIAN1        
        char *nextrname = samlist->recs[0].rname;
        //should always be > 0, because the alignment exists uniqly
        uint64_t nextrpos = samlist->recs[0].pos;
        char nextrc = (samlist->recs[0].flag & 0x10);
        char unmappedmate = !bl_hasMateMapping(s);
#elif CHRISTIAN2
        char *nextrname = NULL;
        int64_t nextrpos = 0;
        char nextrc = 0;
        char unmappedmate = 1;
#else
        char *nextrname = samlist->recs[0].rname;
        //should always be > 0, because the alignment exists uniqly
        uint64_t nextrpos = samlist->recs[0].pos;
        char nextrc = 0;
        char unmappedmate = 1;
#endif
        sam_printEmptyAlign (desc, seq, qual, 1, 1, unmappedmate, nextrname, nextrpos, nextrc, 0, 0, seed, nfo);
      } 

      if(bl_fastaHasMate(reads) && !bl_hasMateMapping(s)) {
        char *seq  = bl_fastaGetMateNoClip(reads, k);
        char *qual = bl_fastaGetMateQualityNoClip(reads, k);
        char *desc = bl_fastaGetMateDescription(reads, k);
#ifdef CHRISTIAN1
        char *nextrname = samlist->recs[0].rname;
        //should always be > 0, because the alignment exists uniqly
        uint64_t nextrpos = samlist->recs[0].pos;
        char nextrc = (samlist->recs[0].flag & 0x10);
        char unmappedmate = !bl_hasQueryMappings(s);
#elif CHRISTIAN2
        char *nextrname = NULL;
        int64_t nextrpos = 0;
        char nextrc = 0;
        char unmappedmate = 1;
#else
        char *nextrname = samlist->recs[0].rname;
        //should always be > 0, because the alignment exists uniqly
        uint64_t nextrpos = samlist->recs[0].pos;
        char nextrc = 0;
        char unmappedmate = 1;
#endif
        sam_printEmptyAlign (desc, seq, qual, 1, 0, unmappedmate, nextrname, nextrpos, nextrc, 0, 0, mateseed, nfo);
      }
    }

    sam_destructSamList (samlist);
    FREEMEMORY(NULL, samlist);
  }

  //now dump the unaligned reads that have no (uniqely) aligned mate
  if(s->n != 1) {

    if(!bl_hasQueryMapping(s)) {
      char *seq  = bl_fastaGetSequenceNoClip(reads, k);
      char *qual = bl_fastaGetQualityNoClip(reads, k);
      char *desc = bl_fastaGetDescription(reads, k);
      char hasPaired = bl_fastaHasMate(reads);
      char *nextrname = NULL;
      int64_t nextrpos = 0;
      char nextrc = 0;

      sam_printEmptyAlign (desc, seq, qual, hasPaired, 1, 1, nextrname, nextrpos, nextrc, 0, 0, seed, nfo);
    }

    if(bl_fastaHasMate(reads) && !bl_hasMateMapping(s)) {
      char *seq  = bl_fastaGetMateNoClip(reads, k);
      char *qual = bl_fastaGetMateQualityNoClip(reads, k);
      char *desc = bl_fastaGetMateDescription(reads, k);
      char hasPaired = bl_fastaHasMate(reads);
      char *nextrname = NULL;
      int64_t nextrpos = 0;
      char nextrc = 0;

      sam_printEmptyAlign (desc, seq, qual, hasPaired, 0, 1, nextrname, nextrpos, nextrc, 0, 0, mateseed, nfo);
    }
  }

  //now dump to additional no-match-device if requested
  if(!bl_hasQueryMapping(s)) {
    char *seq  = bl_fastaGetSequenceNoClip(reads, k);
    char *qual = bl_fastaGetQualityNoClip(reads, k);
    char *desc = bl_fastaGetDescription(reads, k);

    if(nfo->nomatchdev){ 
      if (nfo->threadno > 1) pthread_mutex_lock(nfo->mtx2);
      se_printUnmappedFastx(nfo->nomatchdev, desc, seq, qual, seed);     
      fflush(nfo->nomatchdev);
      if (nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx2);
    }
  } 

  if(bl_fastaHasMate(reads) && !bl_hasMateMapping(s)) {
    char *seq  = bl_fastaGetMateNoClip(reads, k);
    char *qual = bl_fastaGetMateQualityNoClip(reads, k);
    char *desc = bl_fastaGetMateDescription(reads, k);
    if(nfo->nomatchdev){ 
      if (nfo->threadno > 1) pthread_mutex_lock(nfo->mtx2);
      se_printUnmappedFastx(nfo->nomatchdev, desc, seq, qual, mateseed);  
      fflush(nfo->nomatchdev);
      if (nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx2);
    }
  }



  if(nfo->threadno > 1) pthread_mutex_lock(nfo->mtx5);
  //STATS
  nfo->stats->total +=1;
  if(bl_fastaHasMate(reads)) {
    nfo->stats->total +=1;
  }

  if(bl_mappingsetHasPaired(s)) {
    //mapped in pair (+2 mapped reads)
    nfo->stats->mapped+=2;
    nfo->stats->paired+=1;

    if(bl_hasMultiplePairedMappings(s)) {
      //multiple pair (+2 uniquely mapped reads) 
      nfo->stats->multiplemapped+=2;
      nfo->stats->multiplepaired+=1;
    } else {
      //uniq pair (+2 multiple mapped reads)
      nfo->stats->uniquemapped+=2;
      nfo->stats->uniquepaired+=1;
    }
  } else {
    //unmapped in pair
    if(bl_hasQueryMapping(s)) {
      //query mapped (+1 mapped reads)
      nfo->stats->mapped+=1;
      nfo->stats->singlequerymapped+=1;
      if(bl_hasMultipleQueryMappings(s)) {
        nfo->stats->multiplemapped+=1;
      } else {
        nfo->stats->uniquemapped+=1;
      }
    } else {
      nfo->stats->unmapped+=1;
    }
    if(bl_hasMateMapping(s)) {
      //mate mapped (+1 mapped reads)
      nfo->stats->mapped+=1;
      nfo->stats->singlematemapped+=1;
      if(bl_hasMultipleMateMappings(s)) {
        nfo->stats->multiplemapped+=1;
      } else {
        nfo->stats->uniquemapped+=1;
      }
    } else if(bl_fastaHasMate(reads)) {
      nfo->stats->unmapped+=1;
    }
  }


  for(i=0; i < s->n; i++) {
    char isquerysplit = bl_isSplitMappingQM(&s->elem[i], 0);
    char ismatesplit = bl_isSplitMappingQM(&s->elem[i], 1);
    if(isquerysplit && ismatesplit) {
      //split pair
      nfo->stats->splitpair +=1;
    } else if(isquerysplit || ismatesplit) {
      //split single
      nfo->stats->singlesplit +=1;
    }
  }

  if(nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx5);

}

/*------------------------------- se_SAMHeader -------------------------------
 *    
 * @brief SAM header
 * @author Steve Hoffmann 
 *   
 */

  char*
se_SAMHeader (void *space, segemehl_t *info, Uint binno)
{

  Uint i,len=1000, curlen=0, size=0;
  Uint *seqlen;
  char sep, lf;
  char **seq;
  char *header;


  if(info->order){ 
    sep = 8;
    lf = 7;
  } else {
    sep = '\t';
    lf = '\n';
  }

  if(info->bins && binno != -1) {
    size = bl_fileBinsDomainsGetList(space, info->bins, &seq, &seqlen);
  } else { 
    size = info->fasta->noofseqs;
    seq = ALLOCMEMORY(space, NULL, char**, size);
    seqlen = ALLOCMEMORY(space, NULL, Uint*, size);

    for(i=0; i < size; i++) {
      seq[i] = bl_fastaGetDescription(info->fasta, i); 
      seqlen[i] = bl_fastaGetSequenceLength(info->fasta, i); 
    }
  }


  len += strlen(VERSION);
  if(info->cmdline)
    len += strlen(info->cmdline); 

  if(info->bins && binno != -1) {
    len += snprintf(NULL, 0, "@SQ%cSN:%s%cLN:%d%c", sep, seq[binno], sep, seqlen[binno], lf);
  } else { 
    for(i=0; i < size; i++) {
      len += snprintf(NULL, 0, "@SQ%cSN:%s%cLN:%d%c", sep, seq[i], sep, seqlen[i], lf);
    }
  }

  if(info->readgroupid) {
    len += snprintf(NULL, 0, "@RG%cID:%s%c", sep, info->readgroupid, lf);
  }

  if(info->readgroupinfo) { 
    len += snprintf(NULL, 0, "%s%c", info->readgroupinfo, lf);
  }

  header = calloc(len, sizeof(char));
  sprintf(header,"@HD%cVN:1.0",sep);
  curlen = strlen(header);

  if(info->order) sprintf(&header[curlen], "%cSO:coordinate", sep);

  curlen = strlen(header);
  sprintf(&header[curlen],"%c",lf);

  for(i=0; i < size; i++) {
    curlen = strlen(header);
    sprintf(&header[curlen],"@SQ%cSN:%s%cLN:%d%c", sep, seq[i], sep, seqlen[i], lf);
  }

  curlen = strlen(header);
  sprintf(&header[curlen],"@RG%cID:%s", sep, info->readgroupid);

  if(info->readgroupinfo) {
    curlen = strlen(header);
    sprintf(&header[curlen],"%s%c", info->readgroupinfo, lf);
  } else {
    curlen = strlen(header);
    sprintf(&header[curlen],"%c",lf);
  }

  curlen = strlen(header);
  sprintf(&header[curlen],"@PG%cID:segemehl", sep);

  curlen = strlen(header);
  sprintf(&header[curlen],"%cVN:%s", sep, VERSION);

  curlen = strlen(header);
  if(info->cmdline)
    sprintf(&header[curlen],"%cCL:%s", sep, info->cmdline);

  curlen = strlen(header);
  sprintf(&header[curlen],"%c",lf);

  if(info->order) {
    curlen = strlen(header);
    sprintf(&header[curlen-1], "%c", 29);
    sprintf(&header[curlen], "%c", '\n'); 
  } 

  FREEMEMORY(space, seq);
  FREEMEMORY(space, seqlen);

  return header;
}

/*----------------------------- se_initChromBins -----------------------------
 *    
 * @brief set up bins for chromosomes
 * @author Steve Hoffmann 
 *   
 */

  bl_fileBinDomains_t*
se_createChromDomains (void *space, fasta_t *f, Uint avgbins, Uint maxbins, 
    char *filetemplate, Uint tmplen)
{
  bl_fileBinDomains_t* domains;
  char **desc;
  Uint *size;
  Uint i, no, total=0;

  no = f->noofseqs;
  if(no > maxbins) return NULL;

  desc = ALLOCMEMORY(space, NULL, char*, no);
  size = ALLOCMEMORY(space, NULL, Uint, no);

  for(i=0; i < no; i++) {
    desc[i] = bl_fastaGetDescription(f, i);
    size[i] = bl_fastaGetSequenceLength(f, i);
    total += size[i];
  } 

  domains = bl_fileBinsDomainsInit(space, desc, size, no, total,  
      avgbins, maxbins, filetemplate, tmplen);

  FREEMEMORY(space, desc);
  FREEMEMORY(space, size);

  return domains;
}

/*----------------------------- se_initChromBins -----------------------------
 *    
 * @brief set up bins for chromosomes
 * @author Steve Hoffmann 
 *   
 */

  bl_fileBins_t*
se_createChromBins (void *space, fasta_t *f, int maxbins, char *template, 
    Uint tmplen)
{
  bl_fileBins_t* bins;
  char **desc;
  Uint i, no;

  no = f->noofseqs;
  if(no > maxbins) return NULL;

  bins = ALLOCMEMORY(space, NULL, bl_fileBins_t, 1);
  desc = ALLOCMEMORY(space, NULL, char*, no);
  bl_fileBinsInit(space, bins);

  for(i=0; i < no; i++) {
    desc[i] = bl_fastaGetDescription(f, i);
  }

  bl_fileBinsAdd(space, bins, no, bl_fileBinCClassAssign, desc, NULL, 
      template, tmplen);

  FREEMEMORY(space, desc);
  return bins;
}

/*-------------------------------- se_createBisulifteBins ---------------------
 *
 * @brief set up bin domains for matching runs and threads,
 *        domain names are simply 0...(noofdomains-1) as strings
 * @author Christian Otto
 *
 */

bl_fileBinDomains_t*
se_createBisulfiteBins (void *space, Uint noofdomains,
    Uint threadno, char *filetemplate, Uint tmplen){
  Uint i, j;
  bl_fileBinDomains_t *domains;
  bl_fileBinClass_t *class;

  domains = ALLOCMEMORY(space, NULL, bl_fileBinDomains_t, 1);
  domains->noofdomains = noofdomains;
  domains->exp = 0;
  domains->domain = ALLOCMEMORY(space, NULL, bl_fileBinDomain_t, noofdomains);

  for (i = 0; i < noofdomains; i++){    
    domains->domain[i].domainsize = threadno;    
    domains->domain[i].domainname = ALLOCMEMORY(space, NULL, char, log10(i+1) + 3);
    snprintf(domains->domain[i].domainname, log10(i+1)+2, "%d", i);

    bl_fileBinsInit(space, &domains->domain[i].bins);
    bl_fileBinsAdd(space, &domains->domain[i].bins, threadno, NULL, NULL, NULL,
        filetemplate, tmplen);

    domains->domain[i].bins.noofbins = threadno;

    for (j = 0; j < domains->domain[i].bins.noofbins; j++){
      class = ALLOCMEMORY(space, NULL, bl_fileBinClass_t, 1);
      class->start = j;
      class->end = j;
      class->classname = NULL;
      domains->domain[i].bins.b[j].id = class;
    }
  }
  /*
     DBG("domains: noofdomains=%u\texp=%u\n", domains->noofdomains, domains->exp);
     for (i = 0; i < domains->noofdomains; i++){
     DBG("domain %u: domainname=%s\tdomainsize=%u\tnoofbins=%u\n", i,
     domains->domain[i].domainname, domains->domain[i].domainsize,
     domains->domain[i].bins.noofbins);
     for (j = 0; j < domains->domain[i].bins.noofbins; j++){
     DBG("bin %u: filename=%s\tstart=%llu\tend=%llu\n", j, domains->domain[i].bins.b[j].fname,
     domains->domain[i].bins.b[j].id->start, domains->domain[i].bins.b[j].id->end);
     }
     }*/
  return domains;

}


/*--------------------------- registerOutputDevice ---------------------------
 *    
 * @brief select the output device for matches and print header
 * @author Steve Hoffmann 
 *   
 */

void
se_registerOutputDevice(void *space, segemehl_t *info) {
  char *buffer;

  if(info->outfile) {
    info->dev=fopen(info->outfile, "wb"); //w -> wb
    setvbuf(info->dev, NULL, _IOFBF, 524288);

    NFO("opening file '%s'.\n", info->outfile);

    if (info->dev == NULL) {
      fprintf(stderr, "Couldn't open file '%s'. Exit forced.\n", info->outfile);
      exit(EXIT_FAILURE);
    }
  }

  buffer = se_SAMHeader(space, info, -1);
  fprintf(info->dev, "%s", buffer);
  FREEMEMORY(space, buffer);
}



/*------------------------------- se_storeHeader ------------------------------
 *    
 * @brief read and store header (delimitted by '\n') from file
 * @author Steve Hoffmann 
 *   
 */
void
se_storeHeader(void *space, char *filename, char **header, Uint *headerlen){
  FILE *fp;
  int ret;

  fp = fopen(filename, "rb+");
  if(!fp) {
    fprintf(stderr,"Couldnt open file '%s'. Exit forced!\n", filename);
    exit(-1);
  }

  ret = bl_fgets(space, fp, header);


  fprintf(stderr, "header is '%s'.\n", *header);

  if (ret == EOF || ret < 0){
    fprintf(stderr, "Couldn't retrieve header information. End of file reached.\n");
    exit(-1);
  }
  fclose(fp);

  *headerlen = (Uint) ret;

  return;
}


/*-------------------------------- getDevice ---------------------------------
 *    
 * @brief get the device
 * @author Steve Hoffmann 
 *   
 */

  FILE*
getDevice (void *space, char *chr, Uint pos, bl_fileBin_t **fx, segemehl_t *nfo)
{
  FILE *dev;
  char *bisulfite;

  if(!nfo->bins) {
    dev = nfo->dev;
  } else {
    if (!nfo->bisulfitemerging){
      *fx = bl_fileBinsDomainGetBin(nfo->bins, chr, pos);
      bl_fileBinsLock(*fx);
      dev=bl_fileBinsOpen(space, *fx, "w");    
    } else {
      bisulfite = calloc(MAX_INT_LENGTH+1, 1);
      sprintf(bisulfite, "%u", (nfo->bisulfite + 1) % 2);
      *fx = bl_fileBinsDomainGetBin(nfo->bins, bisulfite, nfo->threadid);
      //DBG("info.bisulfite=%u\tbisulfite=%s\tthreadid=%u\tfilename=%s\tstart=%llu\tend=%llu\n",
      //    nfo->bisulfite, bisulfite, nfo->threadid, (*fx)->fname, (*fx)->id->start, (*fx)->id->end);
      bl_fileBinsLock(*fx);
      dev=bl_fileBinsOpen(space, *fx, "w");
      free(bisulfite);
    }
  } 

  return dev;
}

/*--------------------------- se_printMappingStats ---------------------------
 *    
 * @brief print the stats of a mapping run
 * @author Steve Hoffmann 
 *   
 */

void
se_printMappingStats(FILE *device, segemehl_t *nfo) {

  mappingstats_t *s = nfo->stats;

  fprintf(device, "total reads: %"PRIu64"\n", s->total);
  fprintf(device, "unmapped reads: %"PRIu64"\t%.1f\n", s->unmapped, ((double)s->unmapped/s->total)*100);
  fprintf(device, "mapped reads: %"PRIu64"\t%.1f\n", s->mapped, ((double)s->mapped/s->total)*100);
  fprintf(device, "uniquely mapped reads: %"PRIu64"\t%.1f\n", s->uniquemapped, ((double)s->uniquemapped/s->total)*100);
  fprintf(device, "multiply mapped reads: %"PRIu64"\t%.1f\n", s->multiplemapped, ((double)s->multiplemapped/s->total)*100);
  fprintf(device, "Pairs\n");
  fprintf(device, "mapped pairs: %"PRIu64"\t%.1f\n", s->paired, ((double)s->paired/((double)s->total/2))*100);
  fprintf(device, "unique pairs: %"PRIu64"\t%.1f\n", s->uniquepaired, ((double)s->uniquepaired/((double)s->total/2))*100);
  fprintf(device, "multiple pairs: %"PRIu64"\t%.1f\n", s->multiplepaired, ((double)s->multiplepaired/((double)s->total/2))*100);
  fprintf(device, "Singletons\n");
  fprintf(device, "single query: %"PRIu64"\t%.1f\n", s->singlequerymapped, ((double)s->singlequerymapped/((double)s->total))*100);
  fprintf(device, "single mate: %"PRIu64"\t%.1f\n", s->singlematemapped, ((double)s->singlematemapped/((double)s->total))*100);
  fprintf(device, "Splits\n");
  uint64_t totalsplits = s->splitpair*2+s->singlesplit;
  fprintf(device, "total reads split: %"PRIu64"\t%.1f\n", totalsplits, ((double)totalsplits/((double)s->total))*100);
  fprintf(device, "split pairs: %"PRIu64"\t%.1f\n", s->splitpair, ((double)s->splitpair/((double)s->total/2))*100);
  fprintf(device, "split single: %"PRIu64"\t%.1f\n", s->singlesplit, ((double)s->singlesplit/((double)s->total))*100);

}
