#ifndef FILEBUFFER_H
#define FILEBUFFER_H

/*
 *
 *	filebuffer.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 22.11.2016 16:46:12 CET  
 *
 */


#define _POSIX_SOURCE
#define _BSD_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pthread.h>


typedef struct {
    int beg; // Index of first element added to buffer.
    int end; // Index of most recent element added to buffer.
    size_t size; // Number of elements in circular buffer.
    char *buffer; 
    FILE *dev;
    pthread_mutex_t *mtx;
} circbuffer_t;


char  bl_circBufferRead(circbuffer_t *cb);
void  bl_circBufferAddByte(circbuffer_t *cb, char c);
void bl_circBufferSaveAddByte(circbuffer_t *cb, char c);
int bl_circBufferAddSave(circbuffer_t *cb, char *data, size_t len);
void bl_circBufferInit(circbuffer_t *cb, int size, FILE *dev, pthread_mutex_t *mtx);
int   bl_circBufferIsEmpty(circbuffer_t *cb);
void  bl_circBufferDestruct(circbuffer_t *cb);
int  bl_circBufferIsFull(circbuffer_t *cb);
circbuffer_t* bl_circBufferInitArray(int n, int size, FILE *dev, pthread_mutex_t *mtx);
void  bl_circBufferDestructArray(circbuffer_t *bufarr, int n);
void bl_circBufferEmptyArray(circbuffer_t *bufarr, int n);

#endif
