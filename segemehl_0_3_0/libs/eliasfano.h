
/*
 *
 *	eliasfano.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 12/12/16 01:26:06 CET  
 *
 */


void eliasFanoCreateDecodingTable(){
 
  int decodingTableDocIdNumber[256] = { 0 };
  int decodingTableHighBitsCarryover[256] = { 0 };
  int decodingTableHighBits[256] = { 0 };

  for (int i = 0; i < 256; i++)
  {
    int zeroCount = 0;
    //for all bits of the byte
    for (int j = 7; j >= 0; j--)
    {
      // in i there is a '1' at position j
      if ((i & (1 << j)) > 0)
      {
        // unary code of high bits of nth docid within this byte
        decodingTableHighBits[i, decodingTableDocIdNumber[i] ] = zeroCount; 
        fprintf(stderr, "writing: decodingTableHighBits[%d, decodingTableDocIdNumber[%d]=%d] = %d\n", i, i, decodingTableDocIdNumber[i], zeroCount); 

        // docIdNumber = number of docid = number of 1 within one byte
        decodingTableDocIdNumber[i]++;
        zeroCount = 0;
      }
      else
      {
        // count 0 since last 1 within i
        zeroCount++;
      }
    }
    // number of trailing zeros (zeros carryover), if whole byte=0 then unaryCodeLength+=8
    decodingTableHighBitsCarryover[i] = zeroCount;
  }
}
