
/*
 *  segemehl.c
 *  segemehl
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/10/2007 02:50:57 PM CEST
 *
 *  Revision of last commit: 
 *  $Rev: 103 $
 *  $Author: steve $
 *  $Date: 2008-12-10 15:18:18 +0100 (Wed, 10 Dec 2008) $
 *
 *
 *  $Id: segemehl.c 103 2008-12-10 14:18:18Z steve $
 *  $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/segemehl.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "memory.h"
#include "biofiles.h"
#include "fileio.h"
#include "filebuffer.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "mmchar.h"
#include "mathematics.h"
#include <sys/types.h>
#include <unistd.h>
#include <sys/times.h>
#include "vtprogressbar.h"
#include "manout.h"
#include <time.h>
#include "sort.h"
#include "list.h"
#include "biofiles.h"
#include "match.h"
//#include "kdmatch.h"
#include "debug.h"
#include "info.h"
#include <pthread.h>
#include "citation.h"
#include "kdseed.h"
#include "manopt.h"
#include "segemehl.h"
#include "manout.h"
#include "fileBins.h"
#include "seqclip.h"
#include "iupac.h"
#include "merge.h"
#include "samout.h"
#include "version.h"
#ifdef HASHING
#include "hash.h"
#endif

pthread_mutex_t updatemtx;
pthread_mutex_t fastamtx;
unsigned char mute=0;

void*
checkclock(void *args) {
  checkthreadinfo_t *t;

  sleep(2);
  cursorVisible();
  t = (checkthreadinfo_t*) args;
  initProgressBarVT();

  while (pthread_mutex_trylock(&updatemtx) != 0) {
    progressBarVT("reads matched.", t->noofseqs, (*t->counter), 25);
  }

  cursorVisible();
  fprintf(stderr, "\n");
  return NULL;
}


/*--------------------------- se_updateProgressBar ---------------------------
 *    
 * @brief keeping the user informed ... somehow
 * @author Steve Hoffmann 
 *   
 */

void
se_updateProgressBar(Uint k, segemehl_t *nfo) {

  if (!nfo->mute) {
    if (nfo->counter == NULL) {
      progressBarVT("reads matched.", nfo->reads->noofseqs, k, 25);
    } else {
      (*nfo->counter)++;
    }
  }
  return;
}


/*------------------------------ bl_fastaMaster ------------------------------
 *    
 * @brief a server of fasta chunks
 * @author Steve Hoffmann 
 *   
 */

  fasta_t*
se_fastaMaster (void *space, fasta_t *f, Uint size, segemehl_t *nfo)
{

  fasta_t *piece = NULL;

  //return full set if no threads are used
  if(nfo->threadno == 1) {
    if(nfo->nextfastaidx[0] ==0) { 
      nfo->nextfastaidx[0] = 1;
      return f;
    } else { 
      return NULL;
    }
  }

  //in case threads are used
  pthread_mutex_lock(nfo->mtx4);

  if(nfo->index) {
    if(nfo->nextfastaidx[0] < f->chunkindex->size) { 
      piece = bl_fastxCopyIndex (space, f, nfo->nextfastaidx[0], size);
      nfo->nextfastaidx[0] += piece->chunkindex->size; 
    }   
  } else { 
    if(nfo->nextfastaidx[0] < f->noofseqs) { 
      piece = bl_fastxCopy (space, f, nfo->nextfastaidx[0], size);
      nfo->nextfastaidx[0] += piece->noofseqs;
    }    
  }

  pthread_mutex_unlock(nfo->mtx4);

  return piece;
}


/*-------------------------------- matchSlave --------------------------------
 *    
 * @brief the slave for threaded matching
 * @author Steve Hoffmann 
 *   
 */

void*
matchSlave(void *args) {
  segemehl_t *t;

  t = (segemehl_t*) args;  
  match(t->space, t->sarray, t->reads, t);
  return NULL;
}


/*--------------------------- se_scanReadGroupFile ---------------------------
 *    
 * @brief get the id from read group file
 * @author Steve Hoffmann 
 *   
 */
/* 
   char*
   se_getIDfromReadGroupFile (char *rgfile)
   {
   stringset_t** set;
   Uint lines = 0;

   set = readcsv(NULL, rgfile, '\t', &lines);

   for(i=0; i < lines; i++) {
   for(j=0; j < set[i]->noofstrings; j++) {
   fprintf(stderr, "%d,%d: %s\n", i, j, set[i]->strings[j]);
   }
   }

   return ;
   }
   */
/*----------------------------------- main -----------------------------------
 *    
 * @brief the main function
 * @author Steve Hoffmann 
 *   
 */


int main(int argc, char** argv) {

  segemehl_t info, *th_info;
  manopt_arg *unflagged;
  manopt_optionset optset;
  manopt_intconstraint threadconstraint;
  manopt_intconstraint accuracyconstraint;
  manopt_intconstraint jumpconstraint;
  manopt_intconstraint bisulfiteconstraint;

  int *space = NULL;
  int i = 0, j, k, qfbaselen = 0, ch=0;
  Uint filebinbasenamelen=0, //splitfilebasenamelen=0, 
       clipseqlen3=0, dmno=0, desclen, headerlen;
  char oldch, newch, *qfbase, *splitfilebasename = NULL, *filebinbasename=NULL, 
    *version, *clipseq3=NULL, **header, *desc, *adapter=NULL, *buffer, *token;
  unsigned int *suflinktable,
               counter=0; 
  unsigned char index = 0,
                brief = 0;
  samheader_t *head = NULL;

  bl_fileBinDomains_t* bins;

  double difsuf,
         difmatch;
  time_t startsuf, endsuf;
  time_t startmatch, endmatch;

  pthread_t *threads;
  pthread_t clockthread;
  checkthreadinfo_t ch_info;
  manopt_arg *dbfilenames;
  threadconstraint.max = 3000;
  threadconstraint.min = 1;
  accuracyconstraint.max = 100;
  accuracyconstraint.min = 0;
  jumpconstraint.max = INT_MAX;
  jumpconstraint.min = 0;
  bisulfiteconstraint.min = 1;
  bisulfiteconstraint.max = 6;

  se_setdefault(&info);
  info.cmdline = malloc(strlen(argv[0])+2);
  strcpy(info.cmdline, argv[0]);

  for(i = 1; i < argc; i++) {    
    info.cmdline = realloc(info.cmdline, strlen(info.cmdline) + 
        strlen(argv[i]) + 10);
    strcat(info.cmdline," ");
    strcat(info.cmdline,argv[i]);
  }

  version = getNiceSVNVersion(VERSION);
  manopt_initoptionset(&optset, argv[0], NULL, 
      "Heuristic mapping of short sequences\n",
      "SEGEMEHL is free software for non-commercial use \n  (C) 2008 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de");

  manopt_blockseparator(&optset, "INPUT");
  manopt(&optset, LISTOPT, 1, 'd', "database", 
      "list of path/filename(s) of database sequence(s)", "<file> [<file>]", 
      NULL, NULL);
  manopt(&optset, STRINGOPT, 0, 'q', "query", 
      "path/filename of query sequences", "<file>", NULL, &info.queryfilename);
  manopt(&optset, STRINGOPT, 0, 'p', "mate", 
      "path/filename of mate pair sequences", "<file>", NULL, &info.matefilename);
  manopt(&optset, REQSTRINGOPT, 0, 'i', "index", 
      "path/filename of db index", "<file>", NULL, &info.idxfilename);
  manopt(&optset, REQSTRINGOPT, 0, 'j', "index2", 
      "path/filename of second db index", "<file>", NULL, &info.idx2filename);
  manopt(&optset, REQSTRINGOPT, 0, 'x', "generate", 
      "generate db index and store to disk", "<file>", NULL, &info.idxfilename); 
  manopt(&optset, REQSTRINGOPT, 0, 'y', "generate2", 
      "generate second db index and store to disk", "<file>", NULL, &info.idx2filename); 
  manopt(&optset, STRINGOPT, 0, 'B', "filebins",
      "file bins with basename <string> for easier data handling", 
      "<string>", NULL, &info.filebinbasename);
  /*manopt(&optset, REQUINTOPT, 0, 'F', "bisulfite", 
    "bisulfite mapping with methylC-seq/Lister et al. (=1) or bs-seq/Cokus et al. protocol (=2), PAR-CLIP with 4SU (=3) or with 6SG (=4)", "<n>",
    &bisulfiteconstraint, &info.bisulfiteprotocol);
    */
  manopt(&optset, REQUINTOPT, 0, 'F', "bisulfite", 
      "bisulfite mapping with methylC-seq/Lister et al. (=1) or bs-seq/Cokus et al. protocol (=2)", "<n>",
      &bisulfiteconstraint, &info.bisulfiteprotocol);

  manopt_blockseparator(&optset, "GENERAL");
  manopt(&optset, REQUINTOPT, 0, 'm', "minsize", 
      "minimum size of queries", "<n>", NULL, &info.minsize);
  manopt(&optset, FLAG, 0, 's', "progressbar", 
      "show a progress bar", NULL, NULL, &info.mute);
  manopt(&optset, FLAG, 0, 'b', "brief", 
      "brief output", NULL, NULL, &brief);
  manopt(&optset, FLAG, 0, 'c', "checkidx", 
      "check index", NULL, NULL, &info.check);
  manopt(&optset, FLAG, 0, 'e', "briefcigar", 
      "brief cigar string (M vs X and =)", NULL, NULL, &info.briefcigar); 
  manopt(&optset, REQUINTOPT, 0, 't', "threads", 
      "start <n> threads", "<n>", &threadconstraint, &info.threadno);
  manopt(&optset, REQSTRINGOPT, 0, 'o', "outfile", 
      "outputfile", "<string>", NULL, &info.outfile);
  manopt(&optset, REQSTRINGOPT, 0, 'u', "nomatchfilename", 
      "filename for unmatched reads", "<file>", NULL, &info.nomatchname); 
  manopt(&optset, REQSTRINGOPT, 0, 'G', "readgroupfile", 
      "filename to read @RG header", "<file>", NULL, &info.readgroupfile); 
  manopt(&optset, REQSTRINGOPT, 0, 'g', "readgroupid", 
      "read group id", "<string>", NULL, &info.readgroupid); 

  manopt_blockseparator(&optset, "SEEDPARAMS");
  manopt(&optset, REQUINTOPT, 0, 'D', "differences", 
      "search seeds initially with <n> differences", "<n>", NULL, &info.k_p);
  manopt(&optset, REQUINTOPT, 0, 'J', "jump", 
      "search seeds with jump size <n> (0=automatic)", "<n>", &jumpconstraint, &info.jump);
  manopt(&optset, FLAG, 0, 'z', "nosuflinks",
      "dont use suflinks (does not affect index construction, for short reads only, increases runtime!)", NULL, NULL, &info.nosuflinks);
  manopt(&optset, REQDBLOPT, 0, 'E', "evalue", 
      "max evalue", "<double>", NULL, &info.maxevalue); 
  manopt(&optset, REQDBLOPT, 0, 'w', "maxsplitevalue", 
      "max evalue for splits", "<double>", NULL, &info.maxsplitevalue);  
  manopt(&optset, REQUINTOPT, 0, 'M', "maxinterval", 
      "maximum width of a suffix array interval, i.e. a query seed will be omitted if it matches more than <n> times", "<n>", NULL, &info.M);  
  manopt(&optset, REQUINTOPT, 0, 'r', "maxout", 
      "maximum number of alignments that will be reported. If set to zero, all alignments will be reported", "<n>", NULL, &info.maxout);  
  manopt(&optset, STRINGOPT, 0, 'S', "splits", 
      "detect split/spliced reads", NULL, NULL, &info.splitfilebasename);
  /*  manopt(&optset, FLAG, 0, 'K', "SEGEMEHL",
      "output SEGEMEHL format (needs to be selected for brief)", NULL, NULL, &info.sam);
      */
  manopt(&optset, FLAG, 0, 'V', "MEOP",
      "output MEOP field for easier variance calling in SAM (XE:Z:)", NULL, NULL, &info.SAMmeop);

  manopt(&optset, FLAG, 0, 0, "nohead",
      "do not output header", NULL, NULL, &info.nohead);
  /*  manopt(&optset, FLAG, 0, 'Z', "PAIR",
      "output pairstatus flag XA:Z: field in SAM", NULL, NULL, &info.SAMpairstat);
      */
  manopt_blockseparator(&optset, "SEEDEXTENSIONPARAMS");
  manopt(&optset, REQUINTOPT, 0, 'n', "extensionpenalty", 
      "penalty for a mismatch during extension", "<n>", NULL, &info.p_mis);
  manopt(&optset, REQUINTOPT, 0, 'X', "dropoff", 
      "dropoff parameter for extension", "<n>", NULL, &info.Xoff);
  manopt_blockseparator(&optset, "ALIGNPARAMS");
  manopt(&optset, REQUINTOPT, 0, 'A', "accuracy", 
      "min percentage of matches per read in semi-global alignment", "<n>", &accuracyconstraint, &info.accuracy);
  manopt(&optset, REQUINTOPT, 0, 'W', "minsplicecover", 
      "min coverage for spliced transcripts", "<n>", &accuracyconstraint, &info.minsplicedaligncover);
  manopt(&optset, REQUINTOPT, 0, 'U', "minfragscore", 
      "min score of a spliced fragment", "<n>", NULL, &info.minfragmentalignscore);
  manopt(&optset, REQUINTOPT, 0, 'Z', "minfraglen", 
      "min length of a spliced fragment", "<n>", NULL, &info.minfragmentalignlen);
  manopt(&optset, REQDBLOPT, 0, 'l', "splicescorescale", 
      "report spliced alignment with score s only if <f>*s is larger than next best spliced alignment", "<f>", NULL, &info.chainscorescale);
  manopt(&optset, REQUINTOPT, 0, 'H', "hitstrategy", 
      "report only best scoring hits (=1) or all (=0)", NULL, NULL, &info.bestonly);  
  manopt(&optset, FLAG, 0, 0, "showalign", 
      "show alignments", NULL, NULL, &info.align);
  manopt(&optset, REQSTRINGOPT, 0, 'P', "prime5",
      "add 5' adapter", "<string>", NULL, &info.softclip5Prime);
  manopt(&optset, REQSTRINGOPT, 0, 'Q', "prime3",
      "add 3' adapter", "<string>", NULL, &info.softclip3Prime);
  manopt(&optset, REQINTOPT, 0, 'R', "clipacc",
      "clipping accuracy", "<n>", NULL, &info.clipacc);
  manopt(&optset, FLAG, 0, 'T', "polyA",
      "clip polyA tail", NULL, NULL, &info.polyA); 
  manopt(&optset, FLAG, 0, 'Y', "autoclip",
      "autoclip unknown 3prime adapter", NULL, NULL, &info.autoclip);
  manopt(&optset, FLAG, 0, 'C', "hardclip",
      "enable hard clipping", NULL, NULL, &info.hardclip);

  manopt(&optset, FLAG, 0, 'O', "order",
      "sorts the output by chromsome and position (might take a while!)", 
      "<n>", NULL, &info.order);
  manopt(&optset, REQINTOPT, 0, 'I', "maxpairinsertsize",
      "maximum size of the inserts (paired end) in case of multiple hits", "<n>", NULL, &info.maxpairinsertsize);
  unflagged = manopt_getopts(&optset, argc, argv);

  if (!manopt_isset(&optset, 'F', "bisulfite")){
    index = manopt_isset(&optset, 'x', NULL);

    if(!(!manopt_isset(&optset, 'i', NULL) ^ !manopt_isset(&optset, 'x', NULL))) {
      manopt_help(&optset, "please give index filename using -i XOR -x option\n");
    } else if(unflagged->noofvalues > 1) { 
      manopt_help(&optset, "unknown argument(s)\n");
    }
  }
  else {
    if (manopt_isset(&optset, 'i', NULL) && manopt_isset(&optset, 'x', NULL)){
      manopt_help(&optset, "please give C/T index filename using -i XOR -x option\n");
    }
    if (manopt_isset(&optset, 'j', NULL) && manopt_isset(&optset, 'y', NULL)){
      manopt_help(&optset, "please give G/A index filename using -j XOR -y option\n");

    }
    if (!manopt_isset(&optset, 'i', NULL) && !manopt_isset(&optset, 'x', NULL) &&
        !manopt_isset(&optset, 'j', NULL) && !manopt_isset(&optset, 'y', NULL)){
      manopt_help(&optset, "please give C/T and/or G/A index filename using (-i XOR -x) AND/OR (-j XOR -y) option\n");      
    } else if(unflagged->noofvalues > 1) { 
      manopt_help(&optset, "unknown argument(s)\n");
    }

    if (manopt_isset(&optset, 'q', "query")){
      if ((manopt_isset(&optset, 'i', NULL) ^ manopt_isset(&optset, 'x', NULL)) &&
          (manopt_isset(&optset, 'j', NULL) ^ manopt_isset(&optset, 'y', NULL))){
        info.bisulfitemerging = 1;
      }
      if (info.nomatchname != NULL){
        NFO("Warning: file with unmapped reads may contain reads that are mapped \
            in one but not in both matching runs.\n", NULL);
      }
      if (manopt_isset(&optset, 'S', "splits")){
        manopt_help(&optset, "split alignments not yet supported with bisulfite mapping\n");
      }
    }
  }

  if(manopt_isset(&optset, 'U', "minfragscore")){
    info.minsplicedalignscore = 2*info.minfragmentalignscore;
  }

  if(manopt_isset(&optset, 'b', "brief")) {
    info.rep_type = 5;
  }

  /*  if(manopt_isset(&optset, 'K', "SAM")) {
  //default
  } else {
  info.rep_type = 15; 
  }
  */

  if(manopt_isset(&optset, 'O', "order") && (
        !manopt_isset(&optset, 'o', "outfile") && 
        !manopt_isset(&optset, 'B', "filebins"))) {
    manopt_help(&optset, "option -O, --order may not be used when output is dumped to stdout!\n");
  }

  if(info.nomatchname != NULL) { 
    info.nomatchdev = fopen(info.nomatchname, "w");
    if(info.nomatchdev == NULL) {
      manopt_help(&optset, "could not open file for unmapped reads. Writing privileges set?\n");
    }
  }


  if(info.readgroupfile) {
    head = sam_readFile(info.readgroupfile);
  } else  {
    head = ALLOCMEMORY(NULL, NULL, samheader_t, 1);
    sam_initHeader(head);
    char *tempid;
    if(info.readgroupid) {
      tempid = ALLOCMEMORY(NULL, NULL, char, strlen(info.readgroupid)+1);
      memmove(tempid, info.readgroupid, strlen(info.readgroupid));
      tempid[strlen(info.readgroupid)] = '\0';
    } else { 
      tempid = ALLOCMEMORY(NULL, NULL, char, 3);
      memmove(tempid, "A1\0", 3);
    }
    sam_addReadGroup(head, tempid, NULL);
  }

  if(head->nrgroups != 1) {
    NFO("read group file contained %d read group IDs. Exactly one is required.\n", head->nrgroups);
    exit(-1);
  } else if(!head->rgroups[0]) {
    NFO("malformed read group id '%s' in file '%s'.\n", head->rgroups[0], info.readgroupfile);
    exit(-1);
  } else { 
    info.readgroupid = head->rgroups[0];
    info.readgroupinfo = head->rgroupsinfo[0];
  }

  NFO("reads assigned to read group '%s'\n", info.readgroupid);

  if(info.queryfilename) {
    NFO("reading queries in '%s'.\n", info.queryfilename);

#ifdef HASHING
    info.index = 1;
#endif

    if(info.index) {
      info.reads = bl_fastxGetSet(space, 
          &info.queryfilename, 1, 1, 0, 1, info.threadno); //info.threadno //1 triggers chunksize
    } else {

      info.reads = bl_fastxRead(space, 
          info.reads, info.queryfilename, 1, 0, 0, 0, 0, bl_fastxAdd);
    }
    NFO("%d query sequences found.\n", info.reads->noofseqs);

#ifdef HASHING
    if (!info.index){
      MSG("Hashing without fasta index\n");
      bl_fastxGetTags(space, info.reads);
    }
    else {
      MSG("Hashing with fasta index\n");
      bl_fastxGetTags3(space, info.reads);      
    }
    exit(-1);
#endif

    if (info.threadno > info.reads->noofseqs) {
      NFO("more threads than queries. Exit forced\n", NULL);
      exit(EXIT_FAILURE);
    }


    if(info.reads->noofseqs < 50 && info.autoclip) {
      NFO("A minimum of 50 queries is reccommended for autoclip.\n", info.reads->noofseqs);  
      MSG("Do you want to proceed with autoclip? (y/n): ");
      while((ch=getchar()) != 'n' && ch != 'y');
      if(ch == 'n') {
        MSG("Do you want to proceed without clipping? (y/n): ");
        while((ch=getchar()) != 'n' && ch != 'y');
        if(ch == 'n') exit(EXIT_SUCCESS);
        else info.autoclip = 0;
      } 
    }

    if(info.autoclip) {
      adapter = bl_seqclipFind3Prime(space, info.reads, 100000, 40, 10);
      NFO("found adapter sequence: '%s'\n", adapter);
      MSG("Do you want to clip? (y/n): ");
      while((ch=getchar()) != 'n' && ch != 'y');
      if(ch == 'n') {
        MSG("Do you want to proceed without clipping? (y/n): ");
        while((ch=getchar()) != 'n' && ch != 'y');
        if(ch == 'n') exit(EXIT_SUCCESS);
      } else {
        info.softclip3Prime = adapter;
      }
    }
  }

  if(info.softclip3Prime) {
    info.softclip3PrimeLen = strlen(info.softclip3Prime);
  }

  if(info.softclip5Prime) {
    info.softclip5PrimeLen = strlen(info.softclip5Prime);
  }

  if(info.queryfilename && info.matefilename) {
    NFO("reading mates in '%s'.\n", info.matefilename);

    if (info.index) {
      info.reads = bl_fastxGetMateSet(space, info.reads, 
          &info.matefilename, 1, 1, 0, 1, info.threadno);
    } else {
      info.reads = 
        bl_fastxRead(space, info.reads, 
            info.matefilename, 1, 0, 0, 0, 0, bl_fastxAddMate);
    }
    NFO("%d mate sequences found.\n", info.reads->noofseqs);
  }


  oldch = newch = ' ';
  for (k = 0; k < 2; k++){
    /* reset counter variables */
    info.counter = 0;
    counter = 0;
    info.totallength = 0;

    /* normal matching run */
    if(!info.bisulfiteprotocol){
      if (k == 1){
        break;
      }
      initIUPAC(1, 1);
      info.bisulfite = 0;
    } else {
      /* initialize bisulfite matching run */
      if (k == 0){
        initIUPAC(2, 1);

        if (manopt_isset(&optset, 'i', NULL) ^ manopt_isset(&optset, 'x', NULL)){
          info.bisulfiterun = 1;
        }
        else {
          info.bisulfiterun = 2;
        }

        /* bisulfite binning in case of two matching runs */
        if (info.bisulfitemerging){
          /* create domain for each matching run with bins for each thread */	  
          if (!info.filebinbasename) {
            qfbase = bl_basename(info.queryfilename);
            qfbaselen = strlen(qfbase);
            filebinbasename = ALLOCMEMORY(space, NULL, char, qfbaselen);
            memmove(filebinbasename, qfbase, bl_fileprefixlen(qfbase));
            filebinbasename[bl_fileprefixlen(qfbase)] = 0;
            info.filebinbasename = filebinbasename;
          }
          filebinbasenamelen = strlen(info.filebinbasename);

          NFO("creating bisulfite bins.\n", NULL);
          info.bins = se_createBisulfiteBins(space, 2, info.threadno, info.filebinbasename, filebinbasenamelen);

          if(info.bins == NULL) {
            NFO("Could not create bisulfite bins %s*! Exit forced.\n", 
                info.filebinbasename);
            exit(-1);
          }
        }
      } else {
        if (manopt_isset(&optset, 'i', NULL) ^ manopt_isset(&optset, 'x', NULL) &&
            manopt_isset(&optset, 'j', NULL) ^ manopt_isset(&optset, 'y', NULL)){
          info.bisulfiterun = 2;
          /* cleanup */
          destructMultiCharSeq(space, info.seq);
          bl_fastaDestruct(space, info.fasta);
          FREEMEMORY(space, info.fasta);
        }
        else {
          break;
        }
      }

      if (info.bisulfiterun == 1){
        oldch = 'C';
        newch = 'T';
        index = manopt_isset(&optset, 'x', NULL);
      } else if (info.bisulfiterun == 2){
        //reset fastaMaster pointer
        info.nextfastaidx[0] = 0;
        oldch = 'G';
        newch = 'A';
        info.idxfilename = info.idx2filename;
        index = manopt_isset(&optset, 'y', NULL);
      }
      /* 
       * set conversion accordingly to run
       * in bisulfite and PARCLIP with 4SG
       * info.bisulfite = 1 in run 1
       * info.bisulfite = 2 in run 2
       */
      info.bisulfite = info.bisulfiterun;
      /*
       * adjustment of conversion in PAR-CLIP with 4SU:
       * info.bisulfite = 3 in run 1
       * info.bisulfite = 4 in run 2
       */
      if (info.bisulfiteprotocol == 3){
        info.bisulfite = info.bisulfiterun + 2;
      }

      /* 
       * set strand accordingly in bisulfite with
       * Lister et al.'s protocol and PARCLIP with 4SU
       * info.strand == 1 (plus) in run 1
       * info.strand == 2 (minus) in run 2
       */
      if (info.bisulfiteprotocol == 1 ||
          info.bisulfiteprotocol == 3){
        info.strand = info.bisulfiterun;
      }
      /*
       * adjustment to PAR-CLIP with 4SG:
       * info.strand == 2 (minus) in run 1
       * info.strand == 1 (plus) in run 2
       */
      if (info.bisulfiteprotocol == 4){
        info.strand = 3 - info.bisulfiterun;
      }

//      NFO("info.bisulfiteprotocol=%d\tinfo.bisulfite=%d\tinfo.strand=%d\tseedconv:%c->%c\n",
//          info.bisulfiteprotocol, info.bisulfite, info.strand, oldch, newch);
//      NFO("bisulfite/parclip mapping run %d\n", info.bisulfiterun, oldch, newch);
    }

    MSG("reading database sequences.\n"); 

    dbfilenames = manopt_getarg(&optset, 'd', "database");
    info.fasta = bl_fastxGetSet(space, dbfilenames->values, 
        dbfilenames->noofvalues, 1, 0, 0, 1);

    NFO("%d database sequences found.\n", info.fasta->noofseqs);

    for(i=0; i < info.fasta->noofseqs; i++) {
      info.totallength += bl_fastaGetSequenceLength(info.fasta, i); 
    }

    for(i=0; i < info.fasta->noofseqs; i++) {
      desclen = bl_fastaGetDescriptionLength(info.fasta, i);
      desc = strclip(space, bl_fastaGetDescription(info.fasta, i), &desclen);
      FREEMEMORY(space, info.fasta->seqs[i]->description);
      info.fasta->seqs[i]->description = desc;
      info.fasta->seqs[i]->descrlen = desclen;
    }

    NFO("total length of db sequences: %u\n", info.totallength);
    //** INTERCEPT
    //  
    //  bl_fastaDestruct(space, info.fasta);
    //  FREEMEMORY(space, info.fasta);
    //  return EXIT_SUCCESS;
    //**

    if (info.bisulfiteprotocol){
      info.seq = concatCharSequences(space, info.fasta->seqs, info.fasta->noofseqs, (char)126, (char)127);

      /* character conversion */
      for (i=0; i < info.fasta->noofseqs; i++){
        strconvert(bl_fastaGetSequence(info.fasta, i), 
            bl_fastaGetSequenceLength(info.fasta, i), oldch, newch);
      }
    }

    if (!info.bisulfitemerging && !info.bins &&
        manopt_isset(&optset, 'B', "filebins") &&  
        manopt_isset(&optset, 'q', "query")) {

      if (!info.filebinbasename) {
        qfbase = bl_basename(info.queryfilename);
        qfbaselen = strlen(qfbase);
        filebinbasename = ALLOCMEMORY(space, NULL, char, qfbaselen);
        memmove(filebinbasename, qfbase, bl_fileprefixlen(qfbase));
        filebinbasename[bl_fileprefixlen(qfbase)] = 0;
        info.filebinbasename = filebinbasename;
      }

      filebinbasenamelen = strlen(info.filebinbasename);
      info.bins = se_createChromDomains(space, info.fasta, 500, 500, 
          info.filebinbasename, filebinbasenamelen); 

      if(info.bins == NULL) {
        NFO("Could not create bins %s*! Try w/o binning! Exit forced.\n", 
            info.filebinbasename);
        exit(-1);
      }
    }


    if(manopt_isset(&optset, 'S', "splits") && 
        manopt_isset(&optset, 'q', "query")) {

      if (info.splitfilebasename) {
        qfbase = info.splitfilebasename;
      } else if (info.outfile){ 
        qfbase = info.outfile;
      } else {
        qfbase = bl_basename(info.queryfilename); 
      }


      info.multisplitfilename = bl_changefilesuffix(qfbase, "mult.bed");
      info.singlesplitfilename = bl_changefilesuffix(qfbase, "sngl.bed");
      info.transsplitfilename = bl_changefilesuffix(qfbase, "trns.bed");

      NFO("writing multi splits to '%s'\n", info.multisplitfilename);
      NFO("writing single splits to '%s'\n", info.singlesplitfilename);
      NFO("writing trans splits to '%s'\n", info.transsplitfilename);
      //    info.multisplitdev = stdout;
      //    info.singlesplitdev = stdout;

      info.multisplitdev = fopen(info.multisplitfilename, "wb");
      setvbuf(info.multisplitdev, NULL, _IOFBF, 524288);
      info.singlesplitdev = fopen(info.singlesplitfilename, "wb");
      setvbuf(info.singlesplitdev, NULL, _IOFBF, 524288);
      info.transsplitdev = fopen(info.transsplitfilename, "wb");
      setvbuf(info.transsplitdev, NULL, _IOFBF, 524288);

      fprintf(info.multisplitdev, "track name=\"MultiSplit:%s\" description=\"segemehl multisplit predictions from %s\" visibility=2 itemRgb=\"On\"\n", 
          info.readgroupid, info.queryfilename);

      fprintf(info.singlesplitdev, "track name=\"SingleSplit:%s\" description=\"segemehl singlesplit predictions from %s\" visibility=2 itemRgb=\"On\"\n", 
          info.readgroupid, info.queryfilename);

      /*
         if(!manopt_isset(&optset, 'B', "filebins")) {
         info.splitdev = fopen(info.splitfilebasename, "w");
         } else { 
         splitfilebasenamelen = strlen(info.splitfilebasename);
         info.splitbins = se_createChromDomains(space, info.fasta, 150, 150, 
         info.splitfilebasename, splitfilebasenamelen); 

         if(info.splitbins == NULL) {
         NFO("Could not create splitbins %s*! Try w/o binning! Exit forced.\n", 
         info.splitfilebasename);
         exit(-1);
         }
         }
         */
      info.split = 1;
    }    

    if(index) {
      time (&startsuf);
      info.sarray = constructSufArr(space, info.fasta->seqs, 
          info.fasta->noofseqs, NULL, mute); 

      for(i=0; i < info.fasta->noofseqs; i++) {
        FREEMEMORY(space, info.fasta->seqs[i]->sequence); 
        info.fasta->seqs[i]->sequence = NULL;
      }

      if (info.check) {
        NFO("checking suffixarray %s\n", info.idxfilename);
        for(i=1; i < info.sarray->numofsuffixes-1; i++) {
          if(!mute) {
            progressBarVT("suffixes checked.", info.sarray->numofsuffixes, i, 25);
          }
          if (strcmp(&info.sarray->seq->sequences[info.sarray->suftab[i-1]],
                &info.sarray->seq->sequences[info.sarray->suftab[i]]) > 0) {
            NFO("suffix array '%s' corrupted! Exit forced.\n", info.idxfilename);
            exit(-1);
          }
        }
      }

      MSG("constructing lcp.\n");
      constructLcp(space, info.sarray);
      MSG("deleting inv_suftab\n");
      //     if(!info.check) {
      FREEMEMORY(space, info.sarray->inv_suftab);
      info.sarray->inv_suftab = NULL;

      MSG("constructing child tab.\n");
      constructchildtab(space, info.sarray);
      MSG("constructing suffix links.\n");
      MSG("constructing id.\n");
      computeId(space, info.sarray);
      MSG("constructing suflinks - bottom up.\n");
      suflinktable = getsufsucc(space, info.sarray);
      MSG("constructing suflinks - top down.\n");
      constructsuflinks(space, info.sarray, suflinktable);
      FREEMEMORY(space, suflinktable);
      time (&endsuf);
      difsuf = difftime (endsuf, startsuf);
      NFO("building  the suffix array has taken %f seconds.\n", difsuf);
      NFO("total length of suffix array was %u.\n", info.totallength);

    } else {

      time (&startsuf);
      NFO("reading suffix array '%s' from disk.\n", info.idxfilename);
      info.sarray=readSuffixarray(space, info.idxfilename, info.fasta->seqs, 
          info.fasta->noofseqs, mute);

      for(i=0; i < info.fasta->noofseqs; i++) {
        FREEMEMORY(space, info.fasta->seqs[i]->sequence); 
        info.fasta->seqs[i]->sequence = NULL;
      }

      time (&endsuf);
      difsuf = difftime (endsuf, startsuf);
      NFO("reading the suffix array has taken %f seconds.\n", difsuf);
    }

    if (info.check) {
      NFO("checking suffixarray %s\n", info.idxfilename);
      for(i=1; i < info.sarray->numofsuffixes-1; i++) {
        if(!mute) {
          progressBarVT("suffixes checked.", info.sarray->numofsuffixes, i, 25);
        }
        if (strcmp(&info.sarray->seq->sequences[info.sarray->suftab[i-1]],
              &info.sarray->seq->sequences[info.sarray->suftab[i]]) > 0) {
          NFO("suffix array '%s' corrupted! Exit forced.\n", info.idxfilename);
          exit(-1);
        }
      }
      checksuflinks(info.sarray, 0, info.sarray->numofsuffixes-1);
    }

    if(index && info.idxfilename) {
      NFO("writing suffix array '%s' to disk.\n", info.idxfilename);
      writeSuffixarray(info.sarray, info.idxfilename); 
    }

    if(info.queryfilename) {

      if (!info.bisulfiteprotocol)
        info.seq = info.sarray->seq;

      if(!info.bins && k == 0)
        se_registerOutputDevice(space, &info);


      if (info.polyA) {
        info.polyAlen = MIN(80, info.reads->maxlen);
        clipseq3 = ALLOCMEMORY(space, NULL, char, info.polyAlen+1);
        memset(&clipseq3[0], 'A', info.polyAlen);
        clipseqlen3 = info.polyAlen;
        clipseq3[info.polyAlen] = 0;
        info.minclipscr3 = 5;
      }

      if(info.softclip3Prime) {
        clipseqlen3 += info.softclip3PrimeLen;
        clipseq3 = ALLOCMEMORY(space, clipseq3, char, clipseqlen3 +1);
        memmove(&clipseq3[info.polyAlen], info.softclip3Prime, info.softclip3PrimeLen);
        clipseq3[clipseqlen3] = 0;
        //info.minclipscr3 = floor((((float)info.softclip3PrimeLen)*info.clipacc)/100.);
        info.minclipscr3 = 5;
      }

      // fprintf(stderr, "clipseq:%s\n", clipseq3);

      info.softclip3Prime = clipseq3;
      info.softclip3PrimeLen = clipseqlen3;


      if(info.softclip5Prime) {
        info.minclipscr5 = floor((((float)info.softclip5PrimeLen)*info.clipacc)/100.);
      }

#ifdef SEEDHASH
      NFO("creating seed hash with size.\n", info.hashsize);
      info.hash = suffixArrayCreateHash(space, info.sarray, info.hashsize);
#endif

      if(info.bufferedwrite) {
        info.sambuffer = bl_circBufferInitArray(info.threadno, 1000000, 
            info.dev, info.mtx3);
        info.snglbuffer = bl_circBufferInitArray(info.threadno, 1000000, 
            info.singlesplitdev, info.mtx6);
        info.multbuffer = bl_circBufferInitArray(info.threadno, 1000000, 
            info.multisplitdev, info.mtx7);
        info.trnsbuffer = bl_circBufferInitArray(info.threadno, 1000000, 
            info.transsplitdev, info.mtx8);
      }

      if (info.threadno > 1) {

        info.counter=&counter;
        NFO("starting %d threads.\n", info.threadno);


        th_info = ALLOCMEMORY(space, NULL, segemehl_t, info.threadno);
        threads = ALLOCMEMORY(space, NULL, pthread_t, info.threadno);
        ch_info.noofseqs = info.reads->noofseqs;
        ch_info.counter = &counter;

        if (!mute && !info.mute) {
          pthread_mutex_init(&updatemtx, NULL);
          pthread_mutex_lock(&updatemtx);
          pthread_create(&clockthread, NULL, checkclock, &ch_info);
        }


        time (&startmatch);

        for(i=0; i < info.threadno; i++) {
          memmove(&th_info[i], &info, sizeof(segemehl_t));
          th_info[i].reads = info.reads;
          th_info[i].threadid = i;
          pthread_create(&threads[i], NULL, matchSlave, &th_info[i]);
        }

        for(i=0; i < info.threadno; i++) {
          pthread_join(threads[i], NULL); 
        } 

        if(!mute && !info.mute) {
          /*notification via mutex - why use signals?*/
          pthread_mutex_unlock(&updatemtx);
          pthread_join(clockthread, NULL);
        }

        fflush(info.dev);
        time (&endmatch);
        difmatch = difftime (endmatch, startmatch);
        NFO("threaded matching w/ suffixarray has taken %f seconds.\n", difmatch);

        FREEMEMORY(space, th_info);
        FREEMEMORY(space, threads);

      } else {

        time (&startmatch);  
        initProgressBarVT();
        match(info.space, info.sarray, info.reads, &info); //match

        time (&endmatch);
        difmatch = difftime (endmatch, startmatch);
        NFO("matching w/ suffixarray has taken %f seconds.\n", difmatch); 
      }

      if(info.bufferedwrite) { 
        bl_circBufferEmptyArray(info.sambuffer, info.threadno); 
        bl_circBufferEmptyArray(info.snglbuffer, info.threadno); 
        bl_circBufferEmptyArray(info.multbuffer, info.threadno); 
        bl_circBufferEmptyArray(info.trnsbuffer, info.threadno); 

        bl_circBufferDestructArray(info.sambuffer, info.threadno); 
        bl_circBufferDestructArray(info.snglbuffer, info.threadno); 
        bl_circBufferDestructArray(info.multbuffer, info.threadno); 
        bl_circBufferDestructArray(info.trnsbuffer, info.threadno); 

        FREEMEMORY(space, info.sambuffer);
        FREEMEMORY(space, info.snglbuffer);
        FREEMEMORY(space, info.multbuffer);
        FREEMEMORY(space, info.trnsbuffer);
      }
    }


    destructSufArr(space, info.sarray);
    } /* END OF for (k = 0; k < 2; k++) */

    /* merge thread-bins */
    if (info.bisulfitemerging){

      bl_fileBinDomainsCloseAll(info.bins);
      bins = NULL;

      /* if no chromosomal binning --> write to output file */
      if (!manopt_isset(&optset, 'B', "filebins")){
        se_registerOutputDevice(space, &info);
      }
      /* otherwise initialize and write to chromosome bins */
      else {
        bins = se_createChromDomains(space, info.fasta, 500, 500, 
            info.filebinbasename, filebinbasenamelen); 
        if(bins == NULL) {
          NFO("Could not create bins %s*! Try w/o binning! Exit forced.\n", 
              info.filebinbasename);
          exit(-1);
        }
      }
      /* reset mapping stat */
      memset(info.stats, 0, sizeof(mappingstats_t));

      /* add reference seqs to sam header for merging  */
      buffer = se_SAMHeader(space, &info, -1);
      token = strtok(buffer, "\n");
      while(token != NULL){
	head = sam_getHeader(head, token);
	token = strtok(NULL, "\n");
      }
      /* do bisulfite merging and cleanup */
      se_mergeBisulfiteBins(info.bins, info.reads, head, info.dev, bins, 1, info.bestonly, info.stats);
      FREEMEMORY(space, buffer);

      /* destruct bins */
      bl_fileBinDomainsDestruct(space, info.bins);
      FREEMEMORY(space, info.bins);

      info.bins = bins;
    }

    if (info.outfile && !info.bins) {
      NFO("closing file '%s'.\n", info.outfile);

      fclose(info.dev);

      if(info.order) {
        NFO("Sorting file '%s'. Consider option '-B' for faster sorting!\n",
            info.outfile);
        /* read and store header information (until first '\n') */
        header = ALLOCMEMORY(space, NULL, char*, 1);
        se_storeHeader(space, info.outfile, header, &headerlen);
        /* 
         * replace header by repeated string to appear again 
         * in first line after sort 
         * NOTE: assumes ascending order in sort (not -r) and
         * sort order to consider not more than 20 columns
         */
        if (headerlen < 40){
          MSG("Warning: header may not be sorted at beginning of file.\n");
        }
        NFO("replacing tabs in '%s'\n", info.outfile);
        bl_freplacestr(info.outfile, "\000\t", 2, 29);

        /* sort file */
        bl_UnixSort(space, info.outfile, SORT[info.rep_type], SORTDELIM);

        /* write header back to file */
        NFO("re-writing header to '%s'\n", info.outfile);
        bl_freplacestr(info.outfile, *header, headerlen, '\n');
        FREEMEMORY(space, *header);
        FREEMEMORY(space, header);
      }
      if (info.align && (info.order || info.bisulfitemerging)){
        NFO("Expanding alignments in '%s'.\n", info.outfile);
        bl_freplacearr(info.outfile, "\007","\n", 1, EOF);
      }
      bl_freplacearr(info.outfile, "\007\010","\n\t", 2, 29);
    }

    if(info.nomatchname != NULL)
      fclose(info.nomatchdev);

    // if(info.splitdev != NULL) fclose(info.splitdev);

    if(info.bins) {
      bl_fileBinDomainsCloseAll(info.bins);

      if(info.order) {
        bl_fileBinDomainsUnixSort(space, info.bins, SORTBIN[info.rep_type], SORTDELIM);
      }

      if(info.align && (info.order || info.bisulfitemerging)) {
        MSG("Expanding alignments in all bins.\n");
        for(i=0; i < info.bins->noofdomains; i++) {
          for(j=0; j < info.bins->domain[i].bins.noofbins; j++) {
            bl_freplacearr(info.bins->domain[i].bins.b[j].fname, "\007",
                "\n", 1, EOF);
          }
        }
      }


      dmno = bl_fileBinsDomainsGetNoOfNames(space, info.bins); 
      header = ALLOCMEMORY(space, NULL, char**, dmno);
      for(i=0; i < dmno; i++) {
        header[i]= se_SAMHeader(space, &info, i);
      }


      bl_fileBinDomainsMerge(space, info.bins, info.filebinbasename, 
          filebinbasenamelen, "sam", 3, header, 1);

      for(i=0; i < dmno; i++) {
        FREEMEMORY(space, header[i]);
      }

      FREEMEMORY(space, header);

      bl_fileBinDomainsDestruct(space, info.bins);
      FREEMEMORY(space, info.bins);

    }

    if(info.splitbins) {
      /*
         bl_fileBinDomainsCloseAll(info.splitbins);

         if(info.order) {
         MSG("sorting split reads\n");
         bl_fileBinDomainsUnixSort(space, info.splitbins, "-k9,9n", SORTDELIM);
         }

         MSG("merging split reads\n");
         bl_fileBinDomainsMerge(space, info.splitbins, info.splitfilebasename, 
         splitfilebasenamelen, "spl", 3, NULL, 1);

         bl_fileBinDomainsDestruct(space, info.splitbins);
         */
      FREEMEMORY(space, info.splitbins);
    }

    MSG("\nStats:\n");
    se_printMappingStats(stderr, &info);
    
    se_destructInfo(space, &info);

    if(adapter) {
      FREEMEMORY(space, adapter);
    }
 
    if(filebinbasename) { 
      FREEMEMORY(space, filebinbasename);
    }
    
    if(splitfilebasename) { 
      FREEMEMORY(space, splitfilebasename);
    }

    manopt_destructoptionset(&optset);
    manopt_destructarg(unflagged);
    free(unflagged);
    sam_destructHeader(head);
    FREEMEMORY(space, head);
    FREEMEMORY(space, version);

    NFO("\nGoodbye.\n %s\n", citerand());  
    return EXIT_SUCCESS;

  }
