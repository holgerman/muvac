/*
 * merge.c
 * functions to merge matches
 *
 *  SVN
 *  Revision of last commit: $Rev: 348 $
 *  Author: $Author: steve $
 *  Date: $Date: 2012-08-24 08:46:52 -0400 (Fri, 24 Aug 2012) $
 *
 *  Id: $Id: merge.c 348 2012-08-24 12:46:52Z steve $
 *  Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/merge.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "info.h"
#include "basic-types.h"
#include "stringutils.h"
#include "mathematics.h"
#include "biofiles.h"
#include "fileBins.h"
#include "matchfilesfields.h"
#include "merge.h"



/*------------------------- bl_mergefilesInit ---------------------------------
 *    
 * @brief init container for multiple merge files
 * @author Christian Otto
 *   
 */
void bl_mergefilesInit(bl_mergefiles_t *files, Uint nooffiles){
    files->f = ALLOCMEMORY(NULL, NULL, bl_mergefile_t, nooffiles);
    files->nooffiles = nooffiles;
}


/*----------------------- bl_mergefilesDestruct -------------------------------
 *    
 * @brief destruct container for multiple merge files
 * @author Christian Otto
 *   
 */
void bl_mergefilesDestruct(bl_mergefiles_t *files){
  if (files->f != NULL){
    FREEMEMORY(NULL, files->f);
    files->f = NULL;
  }
  files->nooffiles = 0;
}

/*--------------------------- bl_mergefileInit --------------------------------
 *    
 * @brief init merge file container
 * @author Christian Otto
 *   
 */
void bl_mergefileInit(bl_mergefile_t *file, FILE *fp){
  file->fp = fp;
  file->eof = 0;
  file->buffer = NULL;
  file->complete = 0;
  file->entry = ALLOCMEMORY(NULL, NULL, bl_mergefilematch_t, 1);
  bl_mergefilematchInit(file->entry);
}


/*------------------------ bl_mergefileDestruct -------------------------------
 *    
 * @brief destruct merge file container
 * @author Christian Otto
 *   
 */
void bl_mergefileDestruct(bl_mergefile_t *file){
  file->fp = NULL;
  file->eof = 0;
  file->complete = 0;
  if (file->buffer != NULL){
    FREEMEMORY(NULL, file->buffer);
    file->buffer = NULL;
  }
  if (file->entry != NULL){
    bl_mergefilematchDestruct(file->entry);
    FREEMEMORY(NULL, file->entry);
    file->entry = NULL;
  }
}


/*------------------------ bl_mergefilematchInit ------------------------------
 *    
 * @brief init merge file entry
 * @author Christian Otto
 *   
 */
void bl_mergefilematchInit(bl_mergefilematch_t *entry){
  entry->qname = NULL;
  entry->matchid = -1;
  entry->read = NULL;
  entry->mate = NULL;
}


/*---------------------- bl_mergefilematchDestruct ----------------------------
 *    
 * @brief destruct merge file entry
 * @author Christian Otto
 *   
 */
void bl_mergefilematchDestruct(bl_mergefilematch_t *entry){
  entry->qname = NULL;
  entry->matchid = -1;
  if (entry->read != NULL){
    sam_destruct(entry->read);
    FREEMEMORY(NULL, entry->read);
    entry->read = NULL;
  }
  if (entry->mate != NULL){
    sam_destruct(entry->mate);
    FREEMEMORY(NULL, entry->mate);
    entry->mate = NULL;
  }
}


/*--------------------- bl_mergefilematchCompareUnmapSeed ---------------------
 *    
 * @brief compare two sam record entries regarding best found seed
 *        regarding maximum sum of seed length (stored in SAM tag ZL). If
 *        there is no unmapped alignment or ZL is not found, the first entry
 *        is deemed better. Returns -1 if first is better than or equal to
 *        second entry, 1 otherwise.
 *        Christian Otto
 *   
 */
int bl_mergeCompareUnmapped(samrec_t *i, samrec_t *j){
  int tmpi, tmpj;
  samtag_t *tag = NULL;
  
  if (i == NULL || j == NULL){
    return -1;
  }

  tag = sam_getTag(i, "ZL");
  if (!tag){
    return -1;
  }
  tmpi = atoi(tag->val);

  tag = sam_getTag(j, "ZL");
  if (!tag){
    return -1;
  }
  tmpj = atoi(tag->val);

  if (tmpi >= tmpj){
    return -1;
  }
  else {
    return 1;
  }
}


/*------------------ bl_mergefilematchRemoveRedundantUnmapped -----------------
 *    
 * @brief remove redundant unmapped alignment records in mergefilematch list.
 *        These are reads/mates for which either aligned records exist or
 *        "better" unmapped records (regarding the found seed length).
 *        Christian Otto
 *   
 */
void bl_mergefilematchRemoveRedundantUnmapped(bl_mergefilematch_t **list, Uint n){
  int i, cmp, tmpiread, tmpimate, bestread, bestmate, tmpbestread, tmpbestmate;

  bestread = -1, bestmate = -1;
  for (i = 0; i < n; i++){
    if (list[i]->read != NULL){
      /* compare against previous best read */
      if (bestread != -1){
	tmpbestread = (list[bestread]->read->flag & 0x4);
	tmpiread = (list[i]->read->flag & 0x4);
	/* previous was aligned */
	if (!tmpbestread){
	  /* remove current if unmapped */
	  if (tmpiread){
	    sam_destruct(list[i]->read);
	    FREEMEMORY(NULL, list[i]->read);
	    list[i]->read = NULL;
	  }
	}
	else {
	  /* remove best if current is aligned */
	  if (!tmpiread){
	    sam_destruct(list[bestread]->read);
	    FREEMEMORY(NULL, list[bestread]->read);
	    list[bestread]->read = NULL;
	    bestread = i;
	  }
	  /* otherwise: take best unmapped */
	  else {
	    cmp = bl_mergeCompareUnmapped(list[bestread]->read, list[i]->read);
	    if (cmp < 0){
	      sam_destruct(list[i]->read);
	      FREEMEMORY(NULL, list[i]->read);
	      list[i]->read = NULL;	    
	    }
	    else {
	      sam_destruct(list[bestread]->read);
	      FREEMEMORY(NULL, list[bestread]->read);
	      list[bestread]->read = NULL;
	      bestread = i;
	    }
	  }
	}
      }
      else {
	bestread = i;
      }
    }
    if (list[i]->mate != NULL){
      /* compare against previous best mate */
      if (bestmate != -1){
	tmpbestmate = (list[bestmate]->mate->flag & 0x4);
	tmpimate = (list[i]->mate->flag & 0x4);
	/* previous was aligned */
	if (!tmpbestmate){
	  /* remove current if unmapped */
	  if (tmpimate){
	    sam_destruct(list[i]->mate);
	    FREEMEMORY(NULL, list[i]->mate);
	    list[i]->mate = NULL;
	  }
	}
	else {
	  /* remove best if current is aligned */
	  if (!tmpimate){
	    sam_destruct(list[bestmate]->mate);
	    FREEMEMORY(NULL, list[bestmate]->mate);
	    list[bestread]->mate = NULL;
	    bestmate = i;
	  }
	  /* otherwise: take best unmapped */
	  else {
	    cmp = bl_mergeCompareUnmapped(list[bestmate]->mate, list[i]->mate);
	    if (cmp < 0){
	      sam_destruct(list[i]->mate);
	      FREEMEMORY(NULL, list[i]->mate);
	      list[i]->mate = NULL;	    
	    }
	    else {
	      sam_destruct(list[bestmate]->mate);
	      FREEMEMORY(NULL, list[bestmate]->mate);
	      list[bestmate]->mate = NULL;
	      bestmate = i;
	    }
	  }
	}
      }
      else {
	bestmate = i;
      }
    }
  }
  return;
}


/*------------------- bl_mergefilematchRemoveNextInformation ------------------
 *    
 * @brief remove/unset next information in terms of RNEXT, PNEXT, and TLEN.
 *        Christian Otto
 *   
 */
void bl_mergefilematchRemoveNextInformation(bl_mergefilematch_t **list, Uint n){
  int i;
  for (i = 0; i < n; i++){
    if (list[i]->read != NULL){
      FREEMEMORY(NULL, list[i]->read->rnext);
      list[i]->read->rnext = my_strdup("*\0");
      list[i]->read->pnext = 0;
      list[i]->read->tlen = 0;
    }
    if (list[i]->mate != NULL){
      FREEMEMORY(NULL, list[i]->mate->rnext);
      list[i]->mate->rnext = my_strdup("*\0");
      list[i]->mate->pnext = 0;
      list[i]->mate->tlen = 0;
    }
  }
  return;
}

/*------------------- bl_mergefilematchRemoveSuboptimalPairs ------------------
 *    
 * @brief remove suboptimal paired alignments in a mergefilematch list which
 *        are alignments that exceeded the minimal pair edit distance.
 *        Christian Otto
 *   
 */
void bl_mergefilematchRemoveSuboptimalPairs(bl_mergefilematch_t **list, Uint n){
  int i, edist, bestedist = -1;

  for (i = 0; i < n; i++){
    if (list[i]->read != NULL && list[i]->mate != NULL){
      edist = atoi(sam_getTag(list[i]->read, "NM")->val);
      edist += atoi(sam_getTag(list[i]->mate, "NM")->val);
      if (bestedist == -1 || edist < bestedist){
	bestedist = edist;
      }
    }
  }
  for (i = 0; i < n; i++){
    if (list[i]->read != NULL && list[i]->mate != NULL){
      edist = atoi(sam_getTag(list[i]->read, "NM")->val);
      edist += atoi(sam_getTag(list[i]->mate, "NM")->val);
      if (edist > bestedist) {      
	sam_destruct(list[i]->read);
	FREEMEMORY(NULL, list[i]->read);
	list[i]->read = NULL;
	sam_destruct(list[i]->mate);
	FREEMEMORY(NULL, list[i]->mate);
	list[i]->mate = NULL;
      }
    }
  }
  return;
}


/*----------------- bl_mergefilematchRemoveSuboptimalSingletons ---------------
 *    
 * @brief remove suboptimal singleton alignments in a mergefilematch list which
 *        are alignments that exceeded the minimal read/mate edit distance.
 *        Note that the minimal edit distance may be different for read and mate
 *        alignments.
 *        Christian Otto
 *   
 */
void bl_mergefilematchRemoveSuboptimalSingletons(bl_mergefilematch_t **list, Uint n){
  int i, edist, bestreadedist = -1, bestmateedist = -1;

  for (i = 0; i < n; i++){
    if (list[i]->read != NULL && !(list[i]->read->flag & 0x4)){
      edist = atoi(sam_getTag(list[i]->read, "NM")->val);
      if (bestreadedist == -1 || edist < bestreadedist){
	bestreadedist = edist;
      }
    }
    if (list[i]->mate != NULL && !(list[i]->mate->flag & 0x4)){
      edist = atoi(sam_getTag(list[i]->mate, "NM")->val);
      if (bestmateedist == -1 || edist < bestmateedist){
	bestmateedist = edist;
      }
    }
  }
  for (i = 0; i < n; i++){
    if (list[i]->read != NULL && !(list[i]->read->flag & 0x4)){
      edist = atoi(sam_getTag(list[i]->read, "NM")->val);
      if (edist > bestreadedist){
	sam_destruct(list[i]->read);
	FREEMEMORY(NULL, list[i]->read);
	list[i]->read = NULL;
      }
    }
    if (list[i]->mate != NULL && !(list[i]->mate->flag & 0x4)){
      edist = atoi(sam_getTag(list[i]->mate, "NM")->val);      
      if (edist > bestmateedist){
	sam_destruct(list[i]->mate);
	FREEMEMORY(NULL, list[i]->mate);
	list[i]->mate = NULL;
      }
    }
  }
  return;
}


/*--------------------- bl_mergefileComparePairingState -----------------------
 *    
 * @brief compare two merge file entries regarding SAM flag in case of
 *        paired-end reads (i.e. best pairing state), returns -1 if
 *        first is better, 1 if second is better, 0 otherwise
 *        @author Christian Otto
 *   
 */
int bl_mergefilematchComparePairingState(bl_mergefilematch_t *i, bl_mergefilematch_t *j,
					 Uint *pairingstate){
  int tmpi, tmpj, tmpiread, tmpimate, tmpjread, tmpjmate;

  tmpiread = (i->read != NULL) && !(i->read->flag & 0x4);
  tmpimate = (i->mate != NULL) && !(i->mate->flag & 0x4);
  tmpjread = (j->read != NULL) && !(j->read->flag & 0x4);
  tmpjmate = (j->mate != NULL) && !(j->mate->flag & 0x4);

  tmpi = tmpiread & tmpimate;
  tmpj = tmpjread & tmpjmate;

  /* check if either one is fully mapped and properly paired */
  if (!(tmpi | tmpj)){
    *pairingstate = NO_PAIR;
  }
  else {
    *pairingstate = PAIR;
    if (tmpi != tmpj){
      return tmpj - tmpi;
    }
    else {
      assert((i->read->flag & 0x2) == (i->mate->flag & 0x2));
      assert((j->read->flag & 0x2) == (j->mate->flag & 0x2));
      tmpi = (i->read->flag & 0x2);
      tmpj = (j->read->flag & 0x2);
      if (tmpi != tmpj){
	return tmpj - tmpi;
      }
    }
  }
  return 0;
}


/*------------------------- bl_mergefileFastaIDCompare -------------------------------
 *    
 * @brief compare two fasta descriptions if they contain the same fasta ID,
 *        in case of paired-end data, it disregards /1 or /2 at the end or
 *        any  differences after the first white space
 *        returns 1 if both descriptions contain the same ID,  0 otherwise
 * @author Christian Otto
 *   
 */
unsigned char
bl_mergefileFastaIDCompare(char *desc1, Uint desc1len, char *desc2, Uint desc2len) {

  char *id, *id2, *tok1, *tok2;
  unsigned char res;

  id = ALLOCMEMORY(NULL, NULL, char, desc1len+2); 
  id2 = ALLOCMEMORY(NULL, NULL, char, desc2len+2); 

  strcpy(id, desc1);
  strcpy(id2, desc2);

  tok1 = strtok(id, "/");
  tok2 = strtok(id2, "/");
  res = (strcmp(tok1, tok2)==0);

  if(!res) { 
    FREEMEMORY(NULL, id);
    FREEMEMORY(NULL, id2);

    id = ALLOCMEMORY(NULL, NULL, char, desc1len+2); 
    id2 = ALLOCMEMORY(NULL, NULL, char, desc2len+2); 

    strcpy(id, desc1);
    strcpy(id2, desc2);

    tok1 = strtok(id, " ");
    tok2 = strtok(id2, " ");
    res = (strcmp(tok1, tok2)==0);
  }

  FREEMEMORY(NULL, id);
  FREEMEMORY(NULL, id2);
  return res;
}


/*------------------------- bl_mergeParseLine ---------------------------------
 *    
 * @brief parses a SAM-formatted line (single or paired-end) and 
 *        inserts the data in the given container
 *        NOTE: split reads not supported up to now
 * @author Christian Otto
 *   
 */
unsigned char bl_mergeParseLine(samheader_t *head, bl_mergefilematch_t *entry, char *line, Uint *len){
  samrec_t *rec = NULL;
  samtag_t *tag = NULL;
  int matchid = -1;

  /* parse record */
  rec = sam_line2rec(line, *len, head);

  /* get matchid (in case of mapped) */
  tag = sam_getTag(rec, "HI");
  if (tag == NULL || strcmp(tag->type, "i") || atoi(tag->val) < 0){
    DBG("Error in reading HI tag for SAM entry: %sExit forced.\n",
	sam_printSamrec2Buffer(rec, '\n'));
    exit(-1);
  }
  matchid = atoi(tag->val);

  if (rec->flag & 0x800){
    DBG("Split reads not supported yet. Exit forced.\n", NULL);
    exit(-1);
  }
  
  /* abort if non-valid flags if paired (either first or second in pair) */
  if (rec->flag & 0x1){
    if (!((rec->flag & 0x40) ^ 
	  (rec->flag & 0x80))){
      DBG("Invalid SAM flag for entry: %sExit forced.\n",
	  sam_printSamrec2Buffer(rec, '\n'));
      exit(-1);
    }
  }
  else {
    if ((rec->flag & 0x2) ||	
	(rec->flag & 0x40) ||
	(rec->flag & 0x80)){
      DBG("Invalid SAM flag for entry: %sExit forced.\n",
	  sam_printSamrec2Buffer(rec, '\n'));
      exit(-1);
    }
  }
  
  if (entry->qname == NULL){
    entry->qname = rec->qname;
    entry->matchid = matchid;
  }
  else {
    if (! bl_mergefileFastaIDCompare(entry->qname, strlen(entry->qname),
				     rec->qname, strlen(rec->qname)) ||
	entry->matchid != matchid){
      sam_destruct(rec);
      FREEMEMORY(NULL, rec);
      return 1;
    }
  }
    
  /* assign as read if unpaired or mate 1 */
  if (!(rec->flag & 0x1) || (rec->flag & 0x40)){
    if (entry->read != NULL){
      DBG("Multiple alignments for read %s with same HI tag value found. Exit forced.\n",
	  entry->qname);
      exit(-1);
    }
    entry->read = rec;
  }
  /* assign as mate otherwise */
  else {
    if (entry->mate != NULL){
      DBG("Multiple alignments for read %s with same HI tag value found. Exit forced.\n",
	  entry->qname);
      exit(-1);
    }
    entry->mate = rec;
  }
  /* set input line as processed via len */
  *len = 0;
  
  if (rec->flag & 0x8){
    return 1;
  }
  else {
    return 0;
  }
}


/*-------------------------- bl_mergeReadNext ---------------------------------
 *    
 * @brief read next entry in mergefile
 * @author Christian Otto
 *   
 */
void bl_mergeReadNext(samheader_t *head, bl_mergefile_t *file){
  Uint buffersize = 1024, len = 0;
  char *buffer = NULL, ch;

  if (file->buffer != NULL){
    len = strlen(file->buffer);
    file->complete = bl_mergeParseLine(head, file->entry, file->buffer, &len);
    assert(len == 0);
    FREEMEMORY(NULL, file->buffer);
    file->buffer = NULL;
  }

  if (!file->complete && !file->eof){
    buffer = ALLOCMEMORY(NULL, NULL, char, buffersize);
    len = 0;
    while((ch = getc(file->fp)) != EOF) {
      /* extend buffer */
      if(len == buffersize-1) {
        buffersize = 2*buffersize+1;
        buffer = ALLOCMEMORY(NULL, buffer, char, buffersize);
      }
      /* process buffer */
      if(ch == '\n' && len > 0) {
        buffer[len] = '\0';

        file->complete = bl_mergeParseLine(head, file->entry, buffer, &len);

        if (file->complete){
	  if (len > 0){
	    file->buffer = buffer;
	  }
          break;
        } else {
          len = 0;
          continue;
        }
      } else {
        if(ch != '\n') buffer[len++] = ch;
      }
    }
    /* set end of file */
    if (!file->eof && ch == EOF){	
      file->eof = 1;
      if (file->entry->qname != NULL){
	file->complete = 1;
      }
    }
    if (file->buffer == NULL){
      FREEMEMORY(NULL, buffer);
    }
  }
}


/*------------------------- bl_mergeUpdateTag ---------------------------------
 *    
 * @brief update SAM tags NH and HI for given record
 * @author Christian Otto
 *   
 */
void bl_mergeUpdateTag(samrec_t *rec, Uint matchid, Uint noofmatches){
  samtag_t *tag;

  /* get and update HI tag (if necessary) */
  tag = sam_getTag(rec, "HI");
  if (tag == NULL || strcmp(tag->type, "i")){
    DBG("HI tag is missing or invalid in SAM entry: %s",  
	  sam_printSamrec2Buffer(rec, '\n'));
    exit(-1);
  }
  FREEMEMORY(NULL, tag->tag);
  FREEMEMORY(NULL, tag->val);
  bl_asprintf(&tag->tag, "HI:i:%d", matchid);
  bl_asprintf(&tag->val, "%d", matchid);

  /* get and update NH tag (if necessary) */
  tag = sam_getTag(rec, "NH");
  if (tag == NULL || strcmp(tag->type, "i")){
    DBG("NH tag is missing or invalid in SAM entry: %s",  
	  sam_printSamrec2Buffer(rec, '\n'));
    exit(-1);
  }
  FREEMEMORY(NULL, tag->tag);
  FREEMEMORY(NULL, tag->val);
  bl_asprintf(&tag->tag, "NH:i:%d", noofmatches);
  bl_asprintf(&tag->val, "%d", noofmatches);
}


/*-------------------- bl_mergeCountAlignedMappings ---------------------------
 *    
 * @brief count aligned mappings (i.e. excluding unmapped) for query and mate
 *         in list of matches
 * @author Christian Otto
 *   
 */
void bl_mergeCountAlignedMappings(bl_mergefilematch_t **list, int n,
				  Uint *nqueries, Uint *nmates){
  Uint i, noofmates = 0, noofqueries = 0;
  for (i=0; i < n; i++){
    if (list[i]->read != NULL && !(list[i]->read->flag & 0x4)){
      noofqueries++;
    }
    if (list[i]->mate != NULL && !(list[i]->mate->flag & 0x4)){
      noofmates++;
    }
  }
  *nqueries = noofqueries;
  *nmates = noofmates;

  return;
}


/*----------------------- bl_mergeCountMappings -------------------------------
 *    
 * @brief count mappings for query and mate in list of matches
 * @author Christian Otto
 *   
 */
void bl_mergeCountMappings(bl_mergefilematch_t **list, int n,
			   Uint *nqueries, Uint *nmates){
  Uint i, noofmates = 0, noofqueries = 0;
  for (i=0; i < n; i++){
    if (list[i]->read != NULL){
      noofqueries++;
    }
    if (list[i]->mate != NULL){
      noofmates++;
    }
  }
  *nqueries = noofqueries;
  *nmates = noofmates;

  return;
}


/*------------------------- bl_mergeBisulfiteBins -----------------------------
 *    
 * @brief  merging of bisulfite bins according to given read order between
 *         different bisulfite matching runs for each bin separately
 * @author Christian Otto
 *   
 */
void
se_mergeBisulfiteBins (bl_fileBinDomains_t *bsdomains, fasta_t *reads, samheader_t *head,
		       FILE *dev, bl_fileBinDomains_t *chrdomains, unsigned char remove,
                       Uint bestonly, mappingstats_t *stats){
  Uint i, j, k, l, curlen, noofbins, pairingstate,
    noofbest, noofqueries, noofmates,
    allocbest = 1000;    
  int cmp;
  char *curkey;
  FILE *fp;
  bl_mergefiles_t files;
  bl_mergefilematch_t **best;

  assert(bsdomains->noofdomains > 0);
  noofbins = bsdomains->domain[0].bins.noofbins;
  for (i = 0; i < bsdomains->noofdomains; i++){
    if (bsdomains->domain[i].bins.noofbins != bsdomains->domain[0].bins.noofbins){
      DBG("Inconsistent noofbins in domains. Exit forced.\n", NULL);
      exit(-1);
    }
  }

  best = ALLOCMEMORY(NULL, NULL, bl_mergefilematch_t **, allocbest);

  NFO("Merging bisulfite bins.\n", NULL);
  /* init and open files */
  bl_mergefilesInit(&files, noofbins * bsdomains->noofdomains);
  k = 0;
  for (i = 0; i < bsdomains->noofdomains; i++){
    for (j = 0; j < bsdomains->domain[i].bins.noofbins; j++){
      bl_mergefileInit(&files.f[k], bl_fileBinsOpen(NULL, &bsdomains->domain[i].bins.b[j], "r"));
      k++;
    }
  }

  /* perform merging */
  for (i = 0; i < reads->noofseqs; i++){
    noofbest = 0;
    
    /* get next key */
    curkey = bl_fastaGetDescription(reads, i);
    curlen = bl_fastaGetDescriptionLength(reads, i);
    //DBG("queryfile:id=%d\tkey=%s\n", i, curkey);
    /* 
     * find match(es) with current key,
     * best pairing state, minimal edist (qry edist + mate edist)
     */
    for (j = 0; j < files.nooffiles; j++){
      while (1){
	/* read next entry */
	if (!files.f[j].complete){
	  bl_mergeReadNext(head, &files.f[j]);
	}
	/* no match left */
	if (!files.f[j].complete){
	  break;
	}
	/*
	DBG("files.f[%d]: curkey=%s\n", j, files.f[j].entry->qname);
	if(files.f[j].entry->read != NULL){
	  sam_printSamrec(stderr, files.f[j].entry->read, '\n');
	}
	if(files.f[j].entry->mate != NULL){
	  sam_printSamrec(stderr, files.f[j].entry->mate, '\n');
	}
	*/
	/* 
	 * next match with different key
	 * Note: allow for partial matches (match one to the end)
	 * due to clipping of /1 or /2 in paired-end data
	 */
	if (! bl_mergefileFastaIDCompare(curkey, curlen, files.f[j].entry->qname,
					 strlen(files.f[j].entry->qname))){
	  break;
	}
	  
	/* compare current with previous best match based on pairing state (no pair, pair, proper pair) */
	if (noofbest > 0){
	  cmp = bl_mergefilematchComparePairingState(best[0], files.f[j].entry, &pairingstate);
	  if (cmp > 0){
	    /* new best match found => destruct previous ones */
	    for (k = 0; k < noofbest; k++){
	      bl_mergefilematchDestruct(best[k]);
	      FREEMEMORY(NULL, best[k]);
	    }
	    noofbest = 0;
	  }
	}
	else {
	  cmp = 0;
	}
	
	if (cmp >= 0){
	  /* extend best buffer */
	  if (noofbest == allocbest - 1){
	    allocbest *= 2;
	    best = ALLOCMEMORY(NULL, best, bl_mergefilematch_t **, allocbest);
	  }
	  /* append current match to best */
	  best[noofbest++] = files.f[j].entry;
	    
	  files.f[j].entry = ALLOCMEMORY(NULL, NULL, bl_mergefilematch_t, 1);
	  bl_mergefilematchInit(files.f[j].entry);
	}
	/* better match already found => clear data */
	else {
	  bl_mergefilematchDestruct(files.f[j].entry);
	}
	/*
	DBG("AFTER COMPARING PAIRING STATE\n", NULL);
	for (k = 0; k < noofbest; k++){
	  if (best[k]->read != NULL){
	    sam_printSamrec(stderr, best[k]->read, '\n');
	  }
	  if (best[k]->mate != NULL){
	    sam_printSamrec(stderr, best[k]->mate, '\n');
	  }
	}
	*/
	
	/* remove redundant unmapped alignments and next information */
	if (noofbest > 1 && pairingstate == NO_PAIR){
	  bl_mergefilematchRemoveRedundantUnmapped(best, noofbest);
	  bl_mergefilematchRemoveNextInformation(best, noofbest);
	}
	/*
	DBG("AFTER REMOVING UNMAPPED\n", NULL);
	for (k = 0; k < noofbest; k++){
	  if (best[k]->read != NULL){
	    sam_printSamrec(stderr, best[k]->read, '\n');
	  }
	  if (best[k]->mate != NULL){
	    sam_printSamrec(stderr, best[k]->mate, '\n');
	  }
	}
	*/

	/* apply best-only */
	if (noofbest > 1 && bestonly){
	  /* compare pair edit distance in case of PAIR or PROPER_PAIR */
	  if (pairingstate != NO_PAIR){
	    bl_mergefilematchRemoveSuboptimalPairs(best, noofbest);
	  }
	  /* otherwise: compare edist distance for reads and mates separately */
	  else {
	    bl_mergefilematchRemoveSuboptimalSingletons(best, noofbest);
	  }
	  /*
	  DBG("AFTER REMOVING SUBOPTIMAL\n", NULL);
	  for (k = 0; k < noofbest; k++){
	    if (best[k]->read != NULL){
	      sam_printSamrec(stderr, best[k]->read, '\n');
	    }
	    if (best[k]->mate != NULL){
	      sam_printSamrec(stderr, best[k]->mate, '\n');
	    }
	  }
	  */
	}
	files.f[j].complete = 0;
      }
    }
    

    /* STATS */
    bl_mergeCountAlignedMappings(best, noofbest, &noofqueries, &noofmates);
    stats->total +=1;
    if (bl_fastaHasMate(reads)){    
      stats->total +=1;
    }
    if(pairingstate == PAIR) {
      //mapped in pair (+2 mapped reads)
      stats->mapped+=2;
      stats->paired+=1;
      
      if (noofqueries > 1){
      //multiple pair (+2 multiple mapped reads) 
	stats->multiplemapped+=2;
	stats->multiplepaired+=1;
      } else {
	//unique pair (+2 uniquely mapped reads)
	stats->uniquemapped+=2;
	stats->uniquepaired+=1;
      }
    } else {
      //unmapped in pair
      if(noofqueries > 0){
	//query mapped (+1 mapped reads)
	stats->mapped+=1;
	if (bl_fastaHasMate(reads)){  
	  stats->singlequerymapped+=1;
	}
	if(noofqueries > 1) {
	  stats->multiplemapped+=1;
	} else {
	  stats->uniquemapped+=1;
	}
      } else {
	stats->unmapped+=1;
      }
      if(noofmates > 0){
	//mate mapped (+1 mapped reads)
	stats->mapped+=1;
	if (bl_fastaHasMate(reads)){ 
	  stats->singlematemapped+=1;
	}
	if(noofmates > 1) {
	  stats->multiplemapped+=1;
	} else {
	  stats->uniquemapped+=1;
	}
      } else if(bl_fastaHasMate(reads)) {
	stats->unmapped+=1;
      }
    }
    
    /* REPORTING */
    k = 0, l = 0;
    bl_mergeCountMappings(best, noofbest, &noofqueries, &noofmates);
    for (j = 0; j < noofbest; j++){
      /* updating SAM flag and tags */
      if (best[j]->read != NULL){
	if (noofqueries > 1){
	  best[j]->read->flag |= 0x100;
	}
	bl_mergeUpdateTag(best[j]->read, k++, noofqueries);
      }
      if (best[j]->mate != NULL){
	if (noofmates > 1){
	  best[j]->mate->flag |= 0x100;
	}
	bl_mergeUpdateTag(best[j]->mate, l++, noofmates);
      }

      /* report output */
      fp = dev;
	
      if(best[j]->read != NULL){
	/* select output device in case of chrdomains */
	if (chrdomains != NULL) {
	  fp = bl_fileBinsOpen(NULL, bl_fileBinsDomainGetBin(chrdomains, best[j]->read->rname,
							     best[j]->read->pos), "w");
	}
	sam_printSamrec(fp, best[j]->read, '\n');
      }
	
      if(best[j]->mate != NULL){
	/* select output device in case of chrdomains */
	if (chrdomains != NULL){
	  fp = bl_fileBinsOpen(NULL, bl_fileBinsDomainGetBin(chrdomains, best[j]->mate->rname,
							     best[j]->mate->pos), "w");
	}
	sam_printSamrec(fp, best[j]->mate, '\n');
      }

      /* clear data */
      bl_mergefilematchDestruct(best[j]);
      FREEMEMORY(NULL, best[j]);
    }
  }

  /* PREPARE FOR NEXT ITERATION */
  for (j = 0; j < files.nooffiles; j++){
    /* check whether match file is entirely processed */
    bl_mergeReadNext(head, &files.f[j]);
    if (!files.f[j].eof || files.f[j].complete){
      DBG("Files not yet entirely processed. Exit forced.\n", NULL);
      //DBG("files.f[%d]: key=%s\n", j, files.f[j].entry->qname);
      //if (files.f[j].entry->read != NULL){	  
      //  sam_printSamrec (stderr, files.f[j].entry->read, '\n');
      //}
      //if (files.f[j].entry->mate != NULL){	  
      //  sam_printSamrec (stderr, files.f[j].entry->mate, '\n');
      //}
      exit(-1);
    }
    /* destruct */
    bl_mergefileDestruct(&files.f[j]);
  }
  bl_mergefilesDestruct(&files);
  
  /* close file */
  for (i = 0; i < bsdomains->noofdomains; i++){
    for (j = 0; j < bsdomains->domain[i].bins.noofbins; j++){
      bl_fileBinsClose(&bsdomains->domain[i].bins.b[j]);
      if (remove){
	bl_rm(NULL, bsdomains->domain[i].bins.b[j].fname);
	bsdomains->domain[i].bins.b[j].unlinked = 1;
      }
    }
  }
  FREEMEMORY(NULL, best);  
}
