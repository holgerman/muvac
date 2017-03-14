
/*
 *  samtest.c
 *  test sam routines
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06/05/2015 07:14:35 AM EDT
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include "samout.h"
#include "manopt.h"
#include <pthread.h>

pthread_mutex_t updatemtx;
unsigned char mute=0;

int main(int argc, char **argv) 
{ 
  manopt_arg *unflagged;
  manopt_optionset optset;
  char *version = "test\0";
  char *qryfilename = NULL;

  manopt_initoptionset(&optset, argv[0], NULL, 
      "Heuristic mapping of short sequences\n",
      "SEGEMEHL is free software for non-commercial use \n  (C) 2008 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de");

  manopt(&optset, STRINGOPT, 0, 'q', "sam", 
      "path/filename to query (sam alignment) ", "<file>", NULL, &qryfilename);
   
  unflagged = manopt_getopts(&optset, argc, argv);

  if (!manopt_isset(&optset, 'q', "sam")){
    manopt_help(&optset, "no sam file given.\n");
  }


  if(qryfilename)
  sam_readFile(qryfilename);
    //unsigned char gzip, struct access *index,)


  manopt_destructoptionset(&optset);
  manopt_destructarg(unflagged);
  free(unflagged);


  return EXIT_SUCCESS;
}
