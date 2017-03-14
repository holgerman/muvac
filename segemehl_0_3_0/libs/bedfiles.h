#ifndef _BEDFILES_
#define _BEDFILES_

/*
 *
 *	bedfiles.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10/30/16 17:41:47 CET  
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringutils.h"
#include "basic-types.h"
#include "charsequence.h"
#include "zran.h"
#include "biofiles.h"

annotationtrack_t* bl_BEDread (void *space, char *filename);
void bl_BEDwrite (annotationtrack_t *track, FILE *dev);

#endif
