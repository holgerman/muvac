#ifndef MERGE_H
#define MERGE_H

/*
 * merge.h
 * functions to merge matches
 *
 *  SVN
 *  Revision of last commit: $Rev: 348 $
 *  Author: $Author: steve $
 *  Date: $Date: 2012-08-24 08:46:52 -0400 (Fri, 24 Aug 2012) $
 *
 *  Id: $Id: merge.h 348 2012-08-24 12:46:52Z steve $
 *  Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/merge.h $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fileBins.h"
#include "biofiles.h"
#include "samout.h"

#define NO_PAIR 0
#define PAIR 1

typedef struct {
  char* qname;
  int matchid;
  samrec_t *read;
  samrec_t *mate;
} bl_mergefilematch_t;

typedef struct {
  /* file pointer */
  FILE *fp;
  /* EOF read */
  unsigned char eof;
  /* current read alignment */
  bl_mergefilematch_t *entry;
  /* current entry complete? */
  unsigned char complete;
  char *buffer;
} bl_mergefile_t;

typedef struct {
  Uint nooffiles;
  bl_mergefile_t *f;
} bl_mergefiles_t;

void bl_mergefilesInit(bl_mergefiles_t *files, Uint nooffiles);
void bl_mergefilesDestruct(bl_mergefiles_t *files);
void bl_mergefileInit(bl_mergefile_t *file, FILE *fp);
void bl_mergefileDestruct(bl_mergefile_t *file);
void bl_mergefilematchInit(bl_mergefilematch_t *match);
int bl_mergeCompareUnmapped(samrec_t *i, samrec_t *j);
int bl_mergefilematchComparePairingState(bl_mergefilematch_t *i, bl_mergefilematch_t *j, Uint *pairingstate);
void bl_mergefilematchRemoveSuboptimalPairs(bl_mergefilematch_t **list, Uint n);
void bl_mergefilematchRemoveSuboptimalSingletons(bl_mergefilematch_t **list, Uint n);
unsigned char bl_mergefileFastaIDCompare(char *desc1, Uint desc1len, char *desc2, Uint desc2len);
void bl_mergefilematchDestruct(bl_mergefilematch_t *match);
unsigned char bl_mergeParseLine(samheader_t* head, bl_mergefilematch_t *match, char *line, Uint *len);
void bl_mergeReadNext(samheader_t *head, bl_mergefile_t *file);
void bl_mergeUpdateTag(samrec_t *rec, Uint matchid, Uint noofmatches);
void bl_mergeCountAlignedMappings(bl_mergefilematch_t **list, int n, Uint *nqueries, Uint *nmates);
void bl_mergeCountMappings(bl_mergefilematch_t **list, int n, Uint *nqueries, Uint *nmates);
void se_mergeBisulfiteBins (bl_fileBinDomains_t *bsdomains, fasta_t *reads, samheader_t *head,
			    FILE *dev, bl_fileBinDomains_t *chrdomains, unsigned char remove,
                            Uint bestonly, mappingstats_t *stats);
#endif /* MERGE_H */
