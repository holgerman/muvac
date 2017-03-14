
/*
 *
 *	matealign.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 29.07.2013 18:07:27 CEST  
 *
 */

mappingset_t* bl_greedypairMappingSets(mappingset_t *l, mappingset_t *r);

unsigned int bl_myersMCSA(MultiCharSeqAlignment *mcsa, MultiCharSeq *mseq, 
    unsigned int *enctab, bitvector *peq, bitvector *D, unsigned int maxE);

mappingset_t *
bl_matealign(mappingset_t *l, MultiCharSeq *mseq, char **seqs, char **quals, char *qname, unsigned int m,
    unsigned int *enctab, bitvector *D, unsigned int maxL, char ismate, Suffixarray *s, segemehl_t *nfo);

mappingset_t *
bl_pairMateMapping (mappingset_t *l, mappingset_t *r, unsigned int insertsize);


