#ifndef VTPROGRESSBAR_H
#define VTPROGRESSBAR_H

/*
 * =====================================================================================
 * 
 *       Filename:  vtprogressbar.h
 * 
 *    Description:  header file for a simple vt100 progressbar
 * 
 *        Version:  1.0
 *        Created:  12/07/06 00:11:54 CET
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:  Steve Hoffmann (SH), shoffmann@zbh.uni-hamburg.de
 *        Company:  Center for Bioinformatics, Hamburg
 * 
 * =====================================================================================
 */
#include "basic-types.h"

void progressBarVT (char *message, Uint complete, Uint processed, Uint size);
void cursorInvisible();
void cursorVisible();
void initProgressBarVT();

#endif
