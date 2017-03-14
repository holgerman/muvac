  CC=gcc
  LD=${CC} 
  CFLAGS= -Wall -pedantic -std=c99 -O3 -g3 -DSORTEDUNMAPPED -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DDBGNFO -DSHOWALIGN -DDBGLEVEL=0 -DPROGNFO -Ilibs -Ilibs/sufarray
  LDFLAGS=-lm -lpthread -lz -lncurses 
  CTAGS=ctags > tags
  LIBS=-lob -lm -lpthread 
  
 
  SAMOUTOBJ =	libs/locus.o\
				libs/match.o\
				libs/mapfrag.o\
				libs/matealign.o\
				libs/queryalign.o\
				libs/splitalign.o\
				libs/samout.o\
                libs/debug.o\
				libs/info.o\
				libs/md5.o\
				libs/stringutils.o\
				libs/mathematics.o\
                libs/memory.o\
				libs/zran.o\
				libs/fileio.o\
				libs/biofiles.o\
				libs/bedfiles.o\
				libs/karlin.o\
				libs/mappingqual.o\
				libs/kdseed.o\
				libs/queue.o\
				libs/sufarray/sufarray.o\
				libs/sufarray/charsequence.o\
				libs/sufarray/multicharseq.o\
			    libs/radixsort.o\
			    libs/fileBins.o\
			 	libs/stack.o\
				libs/aluruSort.o\
				libs/sort.o\
				libs/vtprogressbar.o\
				libs/sufarray/mmchar.o\
                libs/bitArray.o\
				libs/bitVector.o\
				libs/bitvectoralg.o\
                libs/manout.o\
				libs/kdchain.o\
				libs/manopt.o\
				libs/container.o\
				libs/vstack.o\
				libs/vqueue.o\
				libs/alignment.o\
				libs/sw.o\
				libs/nw.o\
				libs/seqclip.o\
				libs/iupac.o\
				libs/hash.o\
				libs/matchfilesfields.o\
				libs/merge.o\
				libs/708.o\
				libs/pigeon.o\
				libs/brendel.o\
				libs/junctions.o\
				libs/filebuffer.o\
				libs/segemehl_helper.o\
				libs/segemehl.o

  SPLICETESTOBJ = libs/debug.o\
				libs/info.o\
                libs/memory.o\
				libs/stringutils.o\
				libs/sort.o\
				libs/vstack.o\
				libs/alignment.o\
				libs/mathematics.o\
				libs/iupac.o\
				libs/sufarray/charsequence.o\
				libs/bitArray.o\
				libs/bitVector.o\
                libs/708.o\
				libs/fileio.o\
				libs/sw.o\
				libs/seqclip.o\
				libs/zran.o\
				libs/biofiles.o\
				libs/splitalign_test.o

  SAMTESTOBJ =	libs/locus.o\
				libs/mapfrag.o\
				libs/matealign.o\
				libs/queryalign.o\
				libs/splitalign.o\
				libs/samout.o\
                libs/debug.o\
				libs/info.o\
				libs/md5.o\
				libs/stringutils.o\
				libs/mathematics.o\
                libs/memory.o\
				libs/zran.o\
				libs/fileio.o\
				libs/biofiles.o\
				libs/karlin.o\
				libs/kdseed.o\
				libs/queue.o\
				libs/sufarray/sufarray.o\
				libs/sufarray/charsequence.o\
				libs/sufarray/multicharseq.o\
			    libs/radixsort.o\
			    libs/fileBins.o\
			 	libs/stack.o\
				libs/aluruSort.o\
				libs/sort.o\
				libs/vtprogressbar.o\
				libs/sufarray/mmchar.o\
                libs/bitArray.o\
				libs/bitVector.o\
				libs/bitvectoralg.o\
                libs/manout.o\
				libs/kdchain.o\
				libs/manopt.o\
				libs/container.o\
				libs/vstack.o\
				libs/vqueue.o\
				libs/alignment.o\
				libs/sw.o\
				libs/seqclip.o\
				libs/iupac.o\
				libs/hash.o\
				libs/matchfilesfields.o\
				libs/merge.o\
				libs/708.o\
				libs/samtest.o



all: segemehl.x

segemehl.x: ${SAMOUTOBJ}
					gcc $(CFLAGS) ${SAMOUTOBJ} -o $@ $(LDFLAGS) -lform -lmenu

samout.mod.x: ${SAMOUTMOD}
					gcc $(CFLAGS) ${SAMOUTMOD} -o $@ $(LDFLAGS) -lform -lmenu

samtest.x: ${SAMTESTOBJ}
					gcc $(CFLAGS) ${SAMTESTOBJ} -o $@ $(LDFLAGS) -lform -lmenu

splicetest.x: ${SPLICETESTOBJ}
					gcc $(CFLAGS) ${SPLICETESTOBJ} -o $@ $(LDFLAGS) -lform -lmenu

clean: 
	rm -f ${SAMOUTOBJ} 
	rm -f ${SAMOUTMOD}
	rm -rf *~

