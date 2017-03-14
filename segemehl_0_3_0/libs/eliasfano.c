
/*
 *  eliasfano.c
 *  implementation of elias fano 
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 11/27/16 03:25:40 CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "basic-types.h"
#include "memory.h"
#include "bitVector.h"
#include "mathematics.h"
#include <time.h>

void
bl_biPAC(bitvector *bv, Uint n){
 Uint i;

   //read the chunks of 64bits into biPAC_word
   for(i=0; i < n/64; i++) {
     //bl_biPAC_word(&bv[0] + 64 * k, &out[0] + bit * k, bit)
   }
}

/*
 *
 *  first step: get average distance between the codes (ascending order)
 *  secnd step: count the leading zeros of average distance
 *
 *
 */

#define GETHIPOS(X)  (((X)/BITVECTOR_WORDSIZE) + 2)

bitvector
bl_buildEliasFano(uint64_t *in, uint64_t n, uint64_t* nbits, uint64_t* lowersize) {
  uint64_t delta, lomask, lastval;
  uint64_t msb, pos=0;
  uint64_t losize, hisize, i;
//  uint32_t sumoofones = 0;
  uint32_t blocksize = 1024;
  uint16_t blocks = ALLOCMEMORY(NULL, NULL, uint16_t, (n/blocksize)+1);
  memset(blocks, 0, sizeof(uint16_t)*((n/blocksize)+1));
  bitvector v, lo, hi;

  //get largest elem in array and get average delta
  delta = (in[n-1] - in[0]) / n;
  //get highest most significant bit by counting leading zeros -> asm bsr
  msb = ((sizeof(uint64_t) << 3) -1) - __builtin_clzll(delta);
  //get lowest bit mask
  lomask = ((uint64_t) 1 << msb)-1; 
  //get maximum size of compressed data
  losize = msb * n;
  hisize = (2+ceil(log2((in[n-1]-in[0])/n)))*n - losize;

  fprintf(stdout, "init a bitvector of size %"PRIu64" + %"PRIu64" = %"PRIu64" (%"PRIu64"MB)", losize, hisize, losize+hisize, (losize+hisize)/(8*1000*1000));

  //init bitvectors
  v = initbitvector(NULL, losize+hisize+(BITVECTOR_WORDSIZE*2));
  lo = &v[0];
  hi = &v[GETHIPOS(losize)];
 
  for(i=1, lastval=in[0]; i < n; i++) {
    bitvector_setword(lo, i*msb, (in[i]-lastval) & lomask); //in[i]-start 
    lastval = in[i];
  }

#ifdef DEBUGELIASFANO
  fprintf(stdout, "lo array follows with mask of size %"PRIu64" totalsize:%"PRIu64", losize:%"PRIu64"\n", msb, losize+hisize, losize);
  dumpbitvector(lo, losize);
#endif

  //store highest bits in unary encoding (for a number k its 0^k 1)
  //since only gaps are encoded, start offset does not need to be taken
  //into account
  for(i=1, lastval=in[0]; i < n; i++) {
    uint64_t val = (in[i]-lastval) >> msb;

#ifdef DEBUGELIASFANO    
   if(i >= 9474980) { 
    fprintf(stderr, "%"PRIu64": gap:%"PRIu64" -> val: %"PRIu64", pos:%"PRIu64"\n", 
        i, in[i]-lastval, val, pos);
   }
#endif

    bitvector_setbit(hi, pos, 1); 
    sumofones++;

    if(sumofones % blocksize == 0) {
      uint64_t blockid = sumofones/blocksize;
      blocks[blockid] = pos;
    }

    blocks[blockid] += 1;
    pos += val+1;
    lastval = in[i];
  }
    
  bitvector_setbit(hi, pos, 1);
#ifdef DEBUGELIASFANO 
  fprintf(stdout, "hi array follows with mask of size %"PRIu64" totalsize:%"PRIu64", hisize:%"PRIu64"\n", msb, losize+hisize, hisize);
  dumpbitvector(hi, hisize);
#endif

  *nbits = pos+1;
  *lowersize = losize;

  return v; 
}

int
unarydecode(bitvector code, uint64_t len, bitvector rest, uint64_t masksize, uint64_t offset, uint64_t *out) {
  uint64_t pos = 0;
  uint64_t val = 0, newval = 0;
  uint64_t mask =0;
  uint64_t noofwords = len/64;
  uint64_t lastval = offset;
  noofwords += (len % 64 > 0) ? 1 : 0;

  mask = (1 << masksize)-1;

  for(uint64_t k = 0; k < noofwords; k++) {
    uint64_t word = code[k];
    while (word != 0) {
      uint64_t t = word & -word;
      uint64_t r = __builtin_popcountl(t-1);//_mm_popcnt_u64 (t-1);
      newval = k * 64 +  r;
      uint64_t restval = (bitvector_getword(rest, pos*masksize));
      restval &= mask;
#ifdef DEBUGELIASFANO
      if(pos >= 9474980) { 
      fprintf(stderr, "%d: %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64"\n", 
          pos, (newval-val), masksize, ((newval-val) << masksize), restval, lastval);
      }
#endif
      out[pos] = (((newval-val) << masksize) | restval) + lastval;

      lastval = out[pos];
      pos++;
      val = newval + 1;
      //word &= word - 1;
      word ^= t;
    }
  }

  return pos; 
}



int main(int argc, char **argv) {

  uint64_t n = 9999999;
  uint64_t sum = 0;
  time_t epoch;
  epoch = time(NULL);
//  epoch = 1481992140;

  srand(epoch);  
  n = rand() % 4000000000;

  fprintf(stderr, "randomize %"PRIu64" numbers\n", n);
  
  uint64_t *numbers = ALLOCMEMORY(NULL, NULL, uint64_t, n);
    fprintf(stderr, "alloced number\n");
  uint64_t *out = ALLOCMEMORY(NULL, NULL, uint64_t, n);
    fprintf(stderr, "alloced outvec\n");


/*
  bitvector v;
  unsigned char val=4;
  fprintf(stderr, "wrodsize: %d\n", (int)BITVECTOR_WORDSIZE);
  fprintf(stderr, "%u\n", val);
  fprintf(stderr, "%u\n", -val);
  v = initbitvector(NULL, 200);
  bitvector_setbit(v, 40, 1);
  dumpbitvector(v, 200);
  bitvector_setbit(v, 40, 0);
  dumpbitvector(v, 200);
  bitvector_setword(v, 0, -1);
  dumpbitvector(v, 200;
  bitvector_setbit(v, 66, 1);
  dumpbitvector(v, 200);
  bitvector_setbit(v, 64, 1);
  dumpbitvector(v, 200);
  bitvector_setword(v, 1, -1);
  dumpbitvector(v, 200);
  bitvector_setword(v, 0, 0);
  dumpbitvector(v, 200);
  bitvector_setword(v, 16, -1);
  dumpbitvector(v, 200);
  
  for(int i=0; i < 17; i++) { 
    bitvectorLSHIFT(v, v, 200, 1);
    dumpbitvector(v, 200);
  }

  bitvector_setword(v, 33, 1431655765);
  dumpbitvector(v, 200);
*/


  for(uint64_t i=0; i < n; i++) {
    sum += rand() % 10000; 
    numbers[i] = sum;
  }

  fprintf(stderr, "randomized number\n");
  
  uint64_t nbits, losize;
  bitvector hi, lo;
  bitvector v = bl_buildEliasFano(numbers, n, &nbits, &losize);
  lo = &v[0];
  hi = &v[GETHIPOS(losize)];

  fprintf(stderr, "created EF index\n");

  //dumpbitvectorsep(lo, losize, 5);
  unarydecode(hi, nbits, lo, losize/n, numbers[0], out);

  fprintf(stderr, "decoded EF index\n");

  for(uint64_t i=0; i < n; i++) {
  //  fprintf(stderr, "%d: %"PRIu64" - %"PRIu64"\n", i, numbers[i], out[i]);
    if(numbers[i] != out[i]) {
      fprintf(stderr, "%"PRIu64" : %"PRIu64" %"PRIu64"\n", i, numbers[i], out[i]);
      fprintf(stderr, "used seed '%ju'\n", epoch);
    }
    assert(numbers[i]-out[i] == 0);
  }

  fprintf(stderr, "successfully compared in and out data\n");

  FREEMEMORY(NULL, v);
  FREEMEMORY(NULL, numbers);
  FREEMEMORY(NULL, out);

  return EXIT_SUCCESS;
}
