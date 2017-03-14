
/*
 *  maria.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 23.07.2012 17:40:37 CEST
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "alignment.h"
#include "debug.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "sort.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "zran.h"
#include "nw.h"
#include "matchfiles.h"
#include "evalmatchfiles.h"
#include "manout.h"
#include "matchfilesfields.h"
#include "matepairs.h"
#include "matchfiles.h"
#include "manopt.h"
#include "iupac.h"
#include "info.h"
#include "fqueue.h"
#include "list.h"

unsigned char mute = 0;
char *ntcode;



/*---------------------------- matchfileGetAligns ----------------------------
 *    
 * @brief get the alignments for matchloclinks
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileExpandSplitSites(void *space, matchfile_t *file, char *chromname, 
    Uint start, Uint end, unsigned char fields) {

  stringset_t *token;

  Uint buffersize=1024, startbin, //endbin, 
       len=0, k=0, 
       curstart=0, curend=0, 
       acceptorpos=0, donorpos=0, xstart=0, xend = 0, xno = 0, pos=0, 
       acceptorchridx, donorchridx, 
       acceptorflg = 0, donorflg = 0, curchromidx
         //,counter=0, pnext = 0, curcnt = 0
         ; 
  char *buffer = NULL, ch,
       //*curseq=NULL, *curqual=NULL, 
       //*curaln, 
       *filename, strand,
       *acceptorchr = NULL, *donorchr = NULL, 
       //*rnext, 
       *curchrom;
  unsigned char header = 1;
  int ret=0, 
      //readlen=0, u, 
      edist=0;

  matchfileindex_t *index;
  unsigned char gzip, fmt, curedist=0;
  FILE *fp = NULL;
  struct gzidxfile *gzf = NULL;
  
  

#ifdef DBGIDX
  Uint i=0, donorcnt=0, acceptorcnt=0;
#endif

  gzip = file->gzip;
  fmt = file->fmt;
  index = file->index;
  filename = file->filename;

  startbin = (start >> index->exp);
  k = bl_matchfileGetChromIndexNumber(index, chromname); 

  if(k >= index->noofchroms || startbin >= index->noofbins[k]) { 
    return 0;
  }

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  if (gzip) {  
    fp = fopen(filename, "rb");
    gzf = bl_initgzidxfile(fp, index->gzindex, index->bins[k][startbin].offset, CHUNK);
  } else {
    fp = fopen(filename, "r");
    ret = fseeko (fp, index->bins[k][startbin].offset, SEEK_SET);

    if(ret == -1) {
      DBGEXIT("fseeko failed for '%s' Exit forced!\n", filename);
    }
  }

  if(fp == NULL) {
    DBGEXIT("Couldn't open file %s. Exit forced!\n", filename);
  }

  curchromidx = k;
  curend = start;

   while((ch = (gzip) ? bl_getgzidxc(gzf) : getc(fp)) != EOF && 
      curstart <= end && curchromidx == k) {    

    if(len == buffersize-1) {
      buffersize = 2*buffersize+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }

    if(ch == '\n' && len > 0) {

      buffer = ALLOCMEMORY(space, buffer, char, len+1); 
      buffer[len] = '\0';

#ifdef DBGIDX
      DBG("buffer: %s\n", buffer);
#endif

      if(header) header = bl_matchfileIsHeader(buffer, len, fmt);

      if(!header) { 
        token = tokensToStringset(space, "\t", buffer, len);

        curchrom = bl_matchfileGetChrom(token, fmt);
        curstart = bl_matchfileGetStartPos(token, fmt);
        curend   = bl_matchfileGetEndPos(token, fmt);
        curchromidx =  bl_matchfileGetChromIndexNumber(index, curchrom); 

        //fprintf(stderr, "reading %s: %d - %d\n", curchrom, curstart, curend);
        /*last condition to avoid inclusion of 0-alignments in BAM files*/
        if (curstart != curend && 
            curend >= start && curstart <= end && curend+1 > 0 && curchromidx == k) {
          /*
             counter++;
             curseq   = bl_matchfileGetRead(token, fmt);
             readlen  = strlen(curseq);

             if(fields & MFREAD_QUAL) 
             curqual  = bl_matchfileGetQual(token, fmt);

             if(fields & MFREAD_MCNT)
             curcnt   = bl_matchfileGetMatchCnt(token, fmt);

             curaln   = bl_matchfileGetAln(token, fmt);
             rnext    = bl_matchfileGetRNext(token, fmt);
             pnext    = bl_matchfileGetPNext(token, fmt);
             */  
          edist    = bl_matchfileGetEdist(token, fmt);
          strand   = bl_matchfileGetStrand(token, fmt);

          if(fields & MFREAD_SPLITS) { 
            acceptorpos = bl_matchfileGetNextPos(token, fmt);
            acceptorchr = bl_matchfileGetNextChr(token, fmt);
            acceptorflg = bl_matchfileGetNextFlag(token, fmt);
            donorpos = bl_matchfileGetPrevPos(token, fmt);
            donorchr = bl_matchfileGetPrevChr(token, fmt);
            donorflg = bl_matchfileGetPrevFlag(token, fmt);
            xstart = bl_matchfileGetSplitStart(token, fmt);
            xend = bl_matchfileGetSplitEnd(token, fmt);
            xno = bl_matchfileGetSplitNumber(token, fmt);
          }

          if(edist > 255) {
            curedist = 255; 
          } else curedist = edist;

          donorchridx = -1;
          acceptorchridx = -1;
          if (acceptorchr || donorchr) { 

            if(donorchr) {

              donorchridx =  bl_matchfileGetChromIndexNumber(index, donorchr); 

              if(strand == '-') {
                if(donorflg & SPLIT_PREV_PLUS) {

                } else {

                }
                  pos = curend;               
              }

              if(strand == '+') {
                if(!(donorflg & SPLIT_PREV_PLUS)) {

                } else {

                }
                  pos = curstart;           
              }

              if(curchromidx < donorchridx || ( curchromidx == donorchridx && pos < donorpos)) {  
                fprintf(stdout, "%d\t%d\t%d\t%d\t%s\n",curchromidx, pos, donorchridx, donorpos, buffer);
              } else { 
                fprintf(stdout, "%d\t%d\t%d\t%d\t%s\n",donorchridx, donorpos, curchromidx, pos, buffer);
              }
            }


            if(acceptorchr) {

              acceptorchridx =  bl_matchfileGetChromIndexNumber(index, acceptorchr); 

              if(strand == '-') {

                if (acceptorflg & SPLIT_NEXT_PLUS) {
                              
                } else {

                }
                pos = curstart;                 
              }

              if(strand == '+') {
                if (!(acceptorflg & SPLIT_NEXT_PLUS)) { 
               
                } else {

                }
                     pos = curend;                 

              }

              if(curchromidx < donorchridx || (curchromidx == donorchridx && pos < donorpos)) {  
                fprintf(stdout, "%d\t%d\t%d\t%d\t%s\n",curchromidx, pos, acceptorchridx, acceptorpos, buffer);
              } else { 
                fprintf(stdout, "%d\t%d\t%d\t%d\t%s\n",acceptorchridx, acceptorpos, curchromidx, pos, buffer);
              }

            }
          } 
        }

        destructStringset(space, token);
      }

      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      len = 0;

    } else {
      if(ch != '\n') buffer[len++] = ch;
    }
  }

  FREEMEMORY(space, buffer);
  return 0;
}


#ifdef REALIGNTEST

int main(int argc, char **argv) {
  
  void *space = NULL;
  
  manopt_optionset optset;
  manopt_arg *unflagged; 
  manopt_arg *queries;
  manopt_arg *dbfilenames;

  manopt_arg *indices;
  matchfile_t **files = NULL;
  matchfile_t *file = NULL;
  fasta_t *fasta = NULL;
  unsigned char gzip = 0, test=0, saveindex=0, stats=0;//, splice=0
  char version[]="0.1", expand=0;
  int i;
  Uint j, nchr, mincover=0;
  Uint prefixlen=0;
  matchfileindex_t *idx2 = NULL;

   
  initIUPAC(1,1); 
  manopt_initoptionset(&optset, argv[0], NULL, 
      "Expand split sites of segemehl alignments\n",
      "SEGEMEHL is free software for non-commercial use \n  (C) 2008 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de"); 
  manopt(&optset, LISTOPT, 1, 'd', "database", 
      "list of path/filename(s) of database sequence(s)", "<file> [<file> ...]", 
      NULL, NULL);
  manopt(&optset, LISTOPT, 1, 'q', "query", 
      "path/filename of alignment file", "<file> [<file> ...]", NULL, NULL); 
  manopt(&optset, LISTOPT, 0, 'i', "index", 
      "path/filename of db index", "[<file> ... ]", NULL, NULL);



  unflagged = manopt_getopts(&optset, argc, argv);
  saveindex = manopt_isset(&optset, 'x', NULL);
  
  if(!(!manopt_isset(&optset, 'i', NULL) ^ !manopt_isset(&optset, 'x', NULL))) {
    manopt_help(&optset, "please give index filename using -i XOR -x option\n");
  } else if(unflagged->noofvalues > 1) { 
    manopt_help(&optset, "unknown argument(s)\n");
  }

  MSG("reading database sequences.\n"); 


  dbfilenames = manopt_getarg(&optset, 'd', "database");
  fasta = bl_fastxGetSet(space, dbfilenames->values, 
      dbfilenames->noofvalues, 1, 0, 0, 1);

  NFO("%d database sequences found.\n", fasta->noofseqs);
  MSG("reading query files.\n");

  queries = manopt_getarg(&optset, 'q', "query");
  if(queries->noofvalues > 30) {
    manopt_help(&optset, "currently no more than 30 query files allowed\n");
  }

  if(saveindex) {
    indices = manopt_getarg(&optset, 'x', "generate");
  } else {
    indices = manopt_getarg(&optset, 'i', "index");
  }

  if(indices->noofvalues != queries->noofvalues) {
    manopt_help(&optset, "please provide an index file name for each query file\n");
  }

  ntcode  = getNTcodekey(space);
  files = ALLOCMEMORY(space, NULL, matchfile_t*, queries->noofvalues);

  for(i=0; i < queries->noofvalues; i++) {

    files[i] = ALLOCMEMORY(space, NULL, matchfile_t, 1);  
    files[i]->fmt = 0;
    files[i]->index = NULL;
    files[i]->filename = queries->values[i];

    prefixlen = bl_fileprefixlen(files[i]->filename);

    if(strncmp(&files[i]->filename[prefixlen], ".gz", 3) == 0 || 
        strncmp(&files[i]->filename[prefixlen], ".gzip", 3) == 0) {
      gzip = 1;
    }

    files[i]->gzip = gzip;

    if(saveindex) {

      bl_matchfileIndex(space, files[i], fasta);
      bl_matchfileWriteIndex(files[i]->index, indices->values[i]);  
      if(stats) { 
        if(!test) { 
//          getstats(space, files[i], fasta, mincover, maxcover); 
        } else { 
//          teststats(space, files[i],fasta, mincover, maxcover);
        }
      }

      bl_matchfileWriteIndex(files[i]->index, indices->values[i]);  
      idx2 = bl_matchfileReadIndex(space, indices->values[i]);

      fprintf(stderr, "compare index (%p:%p):%d\n", 
          (void*)files[i]->index, (void *)idx2, 
          bl_compareIndices(files[i]->index, idx2));
      bl_matchfileDestructIndex(space, idx2); 
      FREEMEMORY(space, idx2);

    } else if(indices->values[i]) {
    
      MSG("reading index file\n");
      files[i]->index = bl_matchfileReadIndex(space, indices->values[i]);
    }
   
    /*
     *  typically stats should be present if stats
     *  failed however we are not allowed to set
     *  this one
     */
    if(files[i]->index->stats)
    files[i]->index->stats->mincover = mincover;

    /*
     * if stats are not present in any file
     * neither gev, dumpstats, nor call
     * can be called
     */
//    if (!files[i]->index->stats) {  
//      manopt_help(&optset, "please generate an index with stats (option -H) if the options -Z, -c, -S are used\n");
//    }
  }


    for(i=0; i < queries->noofvalues; i++) { 

      file = files[i];
      nchr = file->index->noofchroms;
      for(j=0; j < nchr; j++) { 

        bl_matchfileExpandSplitSites(space, file, file->index->chromnames[j], 
            file->index->matchstart[j], file->index->matchend[j], 255);
      }
    }


 
  bl_fastaDestruct(space, fasta);
  FREEMEMORY(space, fasta);

  if(files) {
    for(i=0; i < queries->noofvalues; i++) { 
      if (files[i]->index) {
        bl_matchfileDestructIndex(space, files[i]->index);
        FREEMEMORY(space, files[i]->index);
      }
      FREEMEMORY(space, files[i]);
    }
    FREEMEMORY(space, files);
  }
 
  manopt_destructoptionset(&optset);
  manopt_destructarg(unflagged);
  FREEMEMORY(space, unflagged);
    exit(EXIT_SUCCESS);
}


#endif

