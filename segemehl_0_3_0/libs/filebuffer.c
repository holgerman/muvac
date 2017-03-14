
/*
 *  filebuffer.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 22.11.2016 16:45:13 CET
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
#include "filebuffer.h"
#include "stringutils.h"
#include "basic-types.h"
#include "memory.h"

char 
bl_circBufferRead(circbuffer_t *cb) {
  char mychar;

  if(cb->beg == cb->end) return 0;

  mychar = cb->buffer[cb->beg];
  cb->beg = (cb->beg + 1) % cb->size;

  return mychar;
}


void 
bl_circBufferAddByte(circbuffer_t *cb, char c) {
  cb->buffer[cb->end] = c;
  cb->end = ( ((cb->end) + 1) % cb->size);
  if (cb->end == cb->beg)
  {
    cb->beg = (cb->beg + 1) % cb->size; 
  }
}

void 
bl_circBufferSaveAddByte(circbuffer_t *cb, char c) {
  char mychar;

  if(bl_circBufferIsFull(cb)) {

    if(cb->mtx) pthread_mutex_lock(cb->mtx);
    
      while(!bl_circBufferIsEmpty(cb)) {
      mychar = bl_circBufferRead(cb);
      fprintf(cb->dev, "%c", mychar);
    }

    if(cb->mtx) pthread_mutex_unlock(cb->mtx); 
  }

  bl_circBufferAddByte(cb, c);
}

void
bl_circBufferEmpty(circbuffer_t *cb) {
   size_t right, left=0;

  if(cb->beg <= cb->end) { 
    right = cb->end - cb->beg;
  } else {
    //cb->beg > cb->end
    right = cb->size - cb->beg;
    left = cb->end;
  }

   if(cb->mtx) pthread_mutex_lock(cb->mtx);

    if(cb->beg <= cb->end) {
      fwrite(&cb->buffer[cb->beg], sizeof(char), right, cb->dev);
    } else {
      fwrite(&cb->buffer[cb->beg], sizeof(char), right, cb->dev);
      fwrite(cb->buffer, sizeof(char), left, cb->dev);
    }
    
    cb->beg = 0;
    cb->end = 0;

   if(cb->mtx) pthread_mutex_unlock(cb->mtx);
}

int
bl_circBufferAddSave(circbuffer_t *cb, char *data, size_t len) {

  size_t fspace=0, right = 0;


  //if data does not fit into buffer - empty the buffer and
  //write straight into file
  if(len > cb->size) {

    if(!bl_circBufferIsEmpty(cb)) {
      bl_circBufferEmpty(cb);
    }
    
    if(cb->mtx) pthread_mutex_lock(cb->mtx);
    fwrite(data, sizeof(char), len, cb->dev);
    if(cb->mtx) pthread_mutex_unlock(cb->mtx);

    return 0;
  } 

  if(cb->beg <= cb->end) { 
    fspace = cb->size - (cb->end - cb->beg);
    right = cb->size - cb->end;
  } else {
    fspace = cb->beg - cb->end;
  }

  //if remaining space does not suffice, dump buffer
  if (len > fspace) { 
    bl_circBufferEmpty(cb);
    right = cb->size;
  } 

  if(cb->beg <= cb->end) {
    if(len > right) {
      memmove(&cb->buffer[cb->end], data, right);
      memmove(cb->buffer, &data[right], len-right);
      cb->end = len-right;
    } else {
      memmove(&cb->buffer[cb->end], data, len);
      cb->end += len;
    }
  } else {
    memmove(&cb->buffer[cb->end], data, len);
    cb->end += len;
  }

  return 0;
}

void 
bl_circBufferInit(circbuffer_t *cb, int size, FILE *dev, pthread_mutex_t *mtx) {
    cb->buffer = calloc(size + 1, sizeof(char));
    cb->size  = size + 1;
    cb->beg = 0;
    cb->end = 0;
    cb->dev = dev;
    cb->mtx = mtx;
}

int 
bl_circBufferIsEmpty(circbuffer_t *cb) {
    return cb->end == cb->beg;
}

void 
bl_circBufferDestruct(circbuffer_t *cb) {
    free(cb->buffer); /* OK if null */ 
}

int 
bl_circBufferIsFull(circbuffer_t *cb) {
    return (cb->end + 1) % cb->size == cb->beg; 
}

circbuffer_t* 
bl_circBufferInitArray(int n, int size, FILE *dev, pthread_mutex_t *mtx) {
  unsigned int i;
  circbuffer_t *bufarr;

  bufarr = ALLOCMEMORY(NULL, NULL, circbuffer_t, n);
  for(i=0; i < n; i++) {
    bl_circBufferInit(&bufarr[i], size, dev, mtx);
  }

  return bufarr;
}

void 
bl_circBufferDestructArray(circbuffer_t *bufarr, int n) {
  unsigned int i;

  for(i=0; i < n; i++) {
    bl_circBufferDestruct(&bufarr[i]);
  }
}

void 
bl_circBufferEmptyArray(circbuffer_t *bufarr, int n) {
  unsigned int i;

  for(i=0; i < n; i++) {
    bl_circBufferEmpty(&bufarr[i]);
  }
}


#ifdef FILEBUFFERTEST

int 
main(int argc, char** argv) {
  unsigned int i;
  char c;
  circbuffer_t mybuffer;
  pthread_mutex_t mtx;
  FILE *fp = fopen("filebuffertest.txt", "wb");
  char mystring[ ] = "This is a test string with a tab\tand a line break\n";

  pthread_mutex_init(&mtx, NULL); 
  bl_circBufferInit(&mybuffer, 198, fp, &mtx);

  fprintf(mybuffer.dev, "now adding byte by byte\n");
  for(i=0; i < 50; i++) {
    bl_circBufferSaveAddByte(&mybuffer, mystring[(i%strlen(mystring))]);
  }
 
  fprintf(mybuffer.dev, "last char added '%c'\n",  mystring[((i-1)%strlen(mystring))]);

  for(i=0; i < 49; i++) {
    c= bl_circBufferRead(&mybuffer);
    if(c) fprintf(mybuffer.dev, "[%d,%d]='%c'\n", mybuffer.beg, mybuffer.end, c);
  }
    
  fprintf(mybuffer.dev, "now adding byte by byte\n");
  for(i=0; i < 25; i++) {
    bl_circBufferSaveAddByte(&mybuffer, mystring[(i%strlen(mystring))]);
  }

  bl_circBufferSaveAddByte(&mybuffer, '\n');
  
  fprintf(mybuffer.dev, "now adding data normally starting at [%d,%d]\n", mybuffer.beg, mybuffer.end);
  for(i=0; i < 3333; i++) {
    bl_circBufferAddSave(&mybuffer, mystring, strlen(mystring));
  }

  fprintf(mybuffer.dev, "empty buffer\n");
  bl_circBufferEmpty(&mybuffer);
  fclose(fp);

  bl_circBufferDestruct(&mybuffer);

}

#endif
