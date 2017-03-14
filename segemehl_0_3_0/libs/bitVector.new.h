#ifndef BITVECTOR_H
#define BITVECTOR_H

/*
 *
 *	bitVector.h
 *  declarations for bit arrays
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/14/2007 04:15:27 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 93 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-07 16:58:47 +0100 (Sun, 07 Dec 2008) $
 *
 *  Id: $Id: bitArray.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/bitArray.h $
 */

#include "basic-types.h"

#define BITVECTOR_WORDSIZE  (sizeof(unsigned long long int)*8)


typedef unsigned long long int bitvector_t;
typedef bitvector_t* bitvector;

extern bitvector initbitvector(void *, uint64_t length);
void dumpbitvector(bitvector a, uint64_t len);
unsigned char valbitvector(bitvector a, uint64_t len, unsigned char val);
extern void setbitvector(bitvector a, uint64_t len, unsigned char val);
bitvector resizebitvector(void *space, bitvector, uint64_t len);
void wrapBitmatrix(void *space, bitvector *, uint64_t m);
void bitvectorLSHIFT(bitvector dest, bitvector a, uint64_t len, uint64_t shift);
void dumpbitvectorsep(bitvector a, uint64_t len, uint64_t sep);

static inline void
bitvector_setbit(bitvector a, uint64_t pos, unsigned char val) {
  uint64_t  byte,
       bits;
  bitvector_t mask=0;

  byte  = pos/BITVECTOR_WORDSIZE;
  bits = pos & (BITVECTOR_WORDSIZE - 1);
  mask = 1;
  mask <<= bits;

  a[byte] ^= ((bitvector_t)-val ^ a[byte]) & mask;
}

static inline void
bitvector_setword(bitvector a, uint64_t pos, bitvector_t val) {
  uint64_t  firstbyte, firstbits, carrybits;
  bitvector_t mask=0;

  firstbyte  = pos/BITVECTOR_WORDSIZE;
  firstbits = pos & (BITVECTOR_WORDSIZE - 1);
  mask = -1;
  mask <<= firstbits;

  //fprintf(stderr, "val: %d, firstbyte: %d, firstbits :%d, mask:%d\n", (int)val, firstbyte, firstbits, (int)mask);

  a[firstbyte] ^= ((val << firstbits) ^ a[firstbyte]) & mask;

  if(firstbits > 0) {
    carrybits = BITVECTOR_WORDSIZE - firstbits;
    mask = -1;
    mask >>= carrybits;
    a[firstbyte+1] ^= ((val >> carrybits) ^ a[firstbyte+1]) & mask;
  }

}

static inline unsigned char
bitvector_getbit(bitvector a, uint64_t pos) {
  uint64_t byte;
  uint64_t bits;
  bitvector_t mask=0;


  byte = pos/BITVECTOR_WORDSIZE;
  bits  = pos & (BITVECTOR_WORDSIZE - 1);
  mask = 1;
  mask <<= bits;

  return ((bitvector_t)a[byte] & mask)? 1 : 0;
}

static inline bitvector_t
bitvector_getword(bitvector a, uint64_t pos) {
  uint64_t byte;
  uint64_t bits;
  bitvector_t word=0;

  byte = pos/BITVECTOR_WORDSIZE;
  bits  = pos & (BITVECTOR_WORDSIZE - 1);

  word = (a[byte] >> bits);

  if(bits) {
    uint64_t mask = (1 << bits)-1;
    word |= ((bitvector_t)a[byte+1] & mask) << (BITVECTOR_WORDSIZE-bits);
  }

  return word;
}



#endif
