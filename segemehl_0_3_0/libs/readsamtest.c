
/*
 *  readsamtest.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 09.12.2016 16:47:47 CET
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


int main(int argc, char **argv) {
  char *filename;
  samheader_t *myhead;

  filename = argv[1];

  fprintf(stderr, "filename: %s\n", filename);
  myhead = sam_readFile(filename);

  return EXIT_SUCCESS;
}
