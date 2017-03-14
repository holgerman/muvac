
/*
 *  bitVector.c
 *  implementations
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/14/2007 04:15:14 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 93 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-07 16:58:47 +0100 (Sun, 07 Dec 2008) $
 *
 *  Id: $Id: bitArray.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/bitArray.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "memory.h"
#include "bitVector.h"
#include <inttypes.h>

inline bitvector
initbitvector(void *space, uint64_t len) {
  uint64_t n;

  bitvector a;

  n = len/BITVECTOR_WORDSIZE;
  n += (len % (BITVECTOR_WORDSIZE) > 0) ? 1 : 0;
  a = calloc(n, sizeof(bitvector_t));

  return a;
}

bitvector
resizebitvector(void *space, bitvector a, uint64_t len) {
  uint64_t n;

  n = len/BITVECTOR_WORDSIZE;
  n += (len % BITVECTOR_WORDSIZE > 0) ? 1 : 0;
  a = ALLOCMEMORY(space, a, bitvector_t, n);

  return a;
}

inline void
setbitvector(bitvector a, uint64_t len, unsigned char val) {
  uint64_t n;

  n = len/BITVECTOR_WORDSIZE;
  n += (len % BITVECTOR_WORDSIZE > 0) ? 1 : 0;
  memset(a,((val) ? 255 :0), n*(sizeof(bitvector_t)));
}

unsigned char
valbitvector(bitvector a, uint64_t len, unsigned char val) {
    uint64_t i;
    bitvector array;

    array = a;

    for(i=0; i < (len/BITVECTOR_WORDSIZE); i++) {  
        if (array[i] != (int) 255) 
          return 0;
    }
    for(i=0; i < (len%BITVECTOR_WORDSIZE); i++){
        if (bitvector_getbit(a, len-i-1)!= val) 
          return 0;
    }
    
    return 1;
}

void
dumpbitvector(bitvector a, uint64_t len) {
  uint64_t i;

  if(len ==0) return;

  for(i=len; i > 0; i--) {
    if(i % 32 == 0) printf("\n");
    printf("%d ", bitvector_getbit(a, i-1));    
  }
  printf("\n");
}

void
dumpbitvectorsep(bitvector a, uint64_t len, uint64_t sep) {
  uint64_t i;
  int cnt=len/sep;
  if(len ==0) return;

  for(i=len; i > 0; i--) {
    if(i % 32 == 0) printf("\n");
    else if(i % sep == 0) printf("| (%d) ", cnt);
    if(i % sep == 0) cnt--; 
    printf("%d ", bitvector_getbit(a, i-1));    
  }
  printf("\n");
}


void
bitvectorAND(bitvector dest, bitvector a, bitvector b, uint64_t len) {
    uint64_t i;
    uint64_t n = len/BITVECTOR_WORDSIZE;

    for(i=0; i <n ; i++) {
        dest[i] = a[i] & b[i];
    }
}

void
bitvectorOR(bitvector dest, bitvector a, bitvector b, uint64_t len) {
    uint64_t i;
    uint64_t n;

    n = len/BITVECTOR_WORDSIZE;

    for(i=0; i < n; i++) {
        dest[i] = a[i] | b[i];
    }
}

void
bitvectorNOT(bitvector dest, bitvector a, uint64_t len) {
    uint64_t i;
    uint64_t n = len/BITVECTOR_WORDSIZE;

    for(i=0; i < n; i++) {
        dest[i] = ~a[i];
    }
}

void
bitvectorXOR(bitvector dest, bitvector a, bitvector b, uint64_t len) {
    uint64_t i;
    uint64_t n = len/BITVECTOR_WORDSIZE;

    for(i=0; i < n; i++) {
        dest[i] = a[i] ^ b[i];
    }
}

unsigned char
bitvectorADD(bitvector dest, bitvector a, bitvector b, uint64_t len){
    uint64_t i,n = len/BITVECTOR_WORDSIZE;
    unsigned char carry = 0;

    for(i=0; i < n; i++) {
        dest[i] = a[i] + b[i] + carry;
        if(carry) {
          carry = ((dest[i] <= a[i]) || (dest[i] <= b[i]));
        } else { 
          carry = ((dest[i] < a[i]) || (dest[i] < b[i]));
        }
    }
    return carry;
}


void
bitvectorLSHIFT(bitvector dest, bitvector a, uint64_t len, uint64_t shift) {
    uint64_t i;

    uint64_t wordshift = shift/BITVECTOR_WORDSIZE;
    uint64_t offset = shift % BITVECTOR_WORDSIZE;
    uint64_t n = len/BITVECTOR_WORDSIZE;
    n += (len % BITVECTOR_WORDSIZE > 0) ? 1 : 0;
    uint64_t suboffset = BITVECTOR_WORDSIZE - offset;

    if (offset == 0) {
        for(i=n-1; i >= wordshift; --i) {
            dest[i] = a[i-wordshift];
        }
    } else {
        for(i=n-1; i > wordshift; --i) {
//            fprintf(stderr, "for %d [ws:%d]: %"PRIu64", %"PRIu64" -> %"PRIu64" %"PRIu64"\n", i, wordshift, a[i-wordshift], a[i-wordshift-1], a[i-wordshift] << offset, a[i-wordshift-1] >> suboffset);
            dest[i] = (a[i-wordshift] << offset) | (a[i-wordshift-1] >> suboffset);
        }
        dest[wordshift] = a[0] << offset; 
    }
}

void
bitvectorRSHIFT(bitvector dest, bitvector a, uint64_t len, uint64_t shift) {
    uint64_t i;

    uint64_t wordshift = shift/BITVECTOR_WORDSIZE;
    uint64_t offset = shift % BITVECTOR_WORDSIZE;
    uint64_t n = len/BITVECTOR_WORDSIZE;
    n += (len % BITVECTOR_WORDSIZE > 0) ? 1 : 0;
    uint64_t limit = n - wordshift -1;
    uint64_t suboffset = BITVECTOR_WORDSIZE - offset;

    if (offset == 0) {
        for(i=0; i <= limit; ++i) {
            dest[i] = a[i+wordshift];
        }
    } else {
        for(i=0; i < limit; ++i) {
            dest[i] = (a[i+wordshift] >> offset) | (a[i+wordshift+1] << suboffset);
        }
        dest[limit] = a[n-1] >> offset; 
    }
}



