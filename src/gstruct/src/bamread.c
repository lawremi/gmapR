static char rcsid[] = "$Id: bamread.c 67539 2012-06-27 04:56:07Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bamread.h"
#include <stdlib.h>
#include "mem.h"
#include "complement.h"
#include "bool.h"
#include "list.h"
#include "samflags.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Bamread_next_bamline */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif



#ifdef HAVE_SAMTOOLS
#include <bam.h>
typedef uint8_t *BAM_Sequence_T;
typedef uint32_t *BAM_Cigar_T;
#endif

#define T Bamreader_T
struct T {

#ifdef HAVE_SAMTOOLS
  bamFile fp;
  bam_header_t *header;

  bam_index_t *idx;
  bam_iter_t iter;

  bam1_t *bam;
  bam1_core_t *core;
#endif

  int region_limited_p;
  int ndivs;
  char **divnames;
  unsigned int *divlengths;
};


#ifdef HAVE_SAMTOOLS
#ifndef bam_init_header_hash
/* Not declared in bam.h */
extern void bam_init_header_hash (bam_header_t *header);
#endif
#endif


void
Bamread_free (T *old) {

#ifdef HAVE_SAMTOOLS
  bam_index_destroy((*old)->idx);
  bam_destroy1((*old)->bam);
  if ((*old)->header != NULL) {
    bam_header_destroy((*old)->header);
  }
  bam_close((*old)->fp);
#endif

  FREE(*old);

  return;
}


T
Bamread_new (char *filename) {
  T new = (T) MALLOC(sizeof(*new));

#ifdef HAVE_SAMTOOLS
  if ((new->fp = bam_open(filename,"rb")) == NULL) {
    fprintf(stderr,"Cannot open BAM file %s\n",filename);
    return (T) NULL;
  }

  if ((new->header = bam_header_read(new->fp)) == NULL) {
    fprintf(stderr,"bam file has no SQ header lines\n");
  } else {
    new->ndivs = new->header->n_targets;
    new->divnames = new->header->target_name;
    new->divlengths = new->header->target_len;
  }

  bam_init_header_hash(new->header);

  if ((new->idx = bam_index_load(filename)) == NULL) {
    fprintf(stderr,"Warning: BAM file %s is not indexed\n",filename);
  }

  new->bam = bam_init1();
  new->core = &(new->bam->core);
#endif

  new->region_limited_p = 0;

  return new;
}

void
Bamread_write_header (T this) {
  int i;

  for (i = 0; i < this->ndivs; i++) {
    printf("@SQ\tSN:%s\tLN:%u\n",this->divnames[i],this->divlengths[i]);
  }

#if 0
  if (sam_read_group_id != NULL) {
    fprintf(fp,"@RG\tID:%s",sam_read_group_id);
    if (sam_read_group_platform != NULL) {
      fprintf(fp,"\tPL:%s",sam_read_group_platform);
    }
    if (sam_read_group_library != NULL) {
      fprintf(fp,"\tLB:%s",sam_read_group_library);
    }
    fprintf(fp,"\tSM:%s",sam_read_group_name);
    fprintf(fp,"\n");
  }
#endif

  return;
}


Genomicpos_T
Bamread_chrlength (T this, char *chr) {
  int32_t tid;

#ifdef HAVE_SAMTOOLS
  /* bam_parse_region(this->header,region,&tid,&chrstart,&chrend); */
  if ((tid = bam_get_tid(this->header,chr)) < 0) {
    fprintf(stderr,"chr %s is not in BAM file\n",chr);
    return 0U;
  } else {
    return this->divlengths[tid];
  }
#else
  return 0U;
#endif

}


#if 0
void
Bamread_reset (T this) {
#ifdef HAVE_SAMTOOLS
  bam_destroy1(this->bam);
  this->bam = bam_init1();
  this->core = &(this->bam->core);
#endif

  this->iter = NULL;
  this->region_limited_p = 0;

  return;
}
#endif


bool
Bamread_limit_region (T this, char *chr, Genomicpos_T chrstart, Genomicpos_T chrend) {
  int32_t tid;
  bool success = false;
#ifdef HAVE_SAMTOOLS
  /* bam_parse_region(this->header,region,&tid,&chrstart,&chrend); */
  if ((tid = bam_get_tid(this->header,chr)) < 0) {
    fprintf(stderr,"chr %s is not in BAM file\n",chr);
    success = false;
  } else {
    this->iter = bam_iter_query(this->idx,tid,(int) (chrstart - 1U),(int) chrend);
    this->region_limited_p = 1;
    success = true;
  }
#endif

  return success;
}


void
Bamread_unlimit_region (T this) {

#ifdef HAVE_SAMTOOLS
  if (this->region_limited_p == 1) {
    this->region_limited_p = 0;
    bam_iter_destroy(this->iter);
    this->iter = NULL;
  }
#endif

  return;
}


int
Bamread_nreads (int *npositions, T this, char *chr, Genomicpos_T chrpos1, Genomicpos_T chrpos2) {
  int nreads = 0;
  Genomicpos_T chrstart, chrend;
  int npositions_lowend, npositions_highend;

  if (chrpos1 < chrpos2) {
    chrstart = chrpos1;
    chrend = chrpos2;
  } else {
    chrstart = chrpos2;
    chrend = chrpos1;
  }
  
#ifdef HAVE_SAMTOOLS
  Genomicpos_T chrpos_low, chrpos_high;
  BAM_Cigar_T cigar;
  int type;
  int i;
  unsigned int length;
  bool lowend_p;
  int max_overhang1_lowend = 0, max_overhang2_lowend = 0,
    max_overhang1_highend = 0, max_overhang2_highend = 0;
  int overhang1, overhang2;
  int max_readlength_lowend = 0, max_readlength_highend = 0;
  int readlength;

  Bamread_limit_region(this,chr,chrstart,chrend);
  while (bam_iter_read(this->fp,this->iter,this->bam) >= 0) {
    readlength = this->core->l_qseq;
    chrpos_high = chrpos_low = this->core->pos + 1U;
    if (this->core->mtid != this->core->tid) {
      lowend_p = true;
    } else if (/*mate_chrpos_low*/this->core->mpos + 1U > chrpos_low) {
      lowend_p = true;
    } else {
      lowend_p = false;
    }

    cigar = bam1_cigar(this->bam);
    for (i = 0; i < this->core->n_cigar; i++) {
      type = (int) ("MIDNSHP"[cigar[i] & BAM_CIGAR_MASK]);
      length = cigar[i] >> BAM_CIGAR_SHIFT;

      if (type == 'S') {
	/* Ignore */

      } else if (type == 'H') {
	/* Ignore */

      } else if (type == 'M') {
	chrpos_high += length;

      } else if (type == 'N') {
	chrpos_high += length;

      } else if (type == 'P') {
	/* Do nothing */

      } else if (type == 'I') {
	/* Do nothing */

      } else if (type == 'D') {
	/* CHECK */
	chrpos_high += length;

      } else {
	fprintf(stderr,"Cannot parse type %c\n",type);
	exit(9);
      }
    }

    /* printf("read %d between %u and %u, at %u..%u\n",
       nreads,chrstart,chrend,chrpos_low,chrpos_high); */
    if (chrpos_low >= chrstart && chrpos_high <= chrend) {
      /* Internal to both bounds */
      nreads++;

    } else if (chrpos_high <= chrend) {
      /* Straddles low bound */
      if (chrpos_high > chrstart + readlength/2) {
	/* More than half of the read is in the intron */
	if (lowend_p == true) {
	  if ((overhang1 = readlength - (chrpos_high - chrstart)) > max_overhang1_lowend) {
	    max_overhang1_lowend = overhang1;
	  }
	  if (readlength > max_readlength_lowend) {
	    max_readlength_lowend = readlength;
	  }
	} else {
	  if ((overhang1 = readlength - (chrpos_high - chrstart)) > max_overhang1_highend) {
	    max_overhang1_highend = overhang1;
	  }
	  if (readlength > max_readlength_highend) {
	    max_readlength_highend = readlength;
	  }
	}
	nreads++;
      }

    } else if (chrpos_low >= chrstart) {
      /* Straddles high bound */
      if (chrpos_low < chrend - readlength/2) {
	/* More than half of the read is in the intron */
	if (lowend_p == true) {
	  if ((overhang2 = readlength - (chrend - chrpos_low)) > max_overhang2_lowend) {
	    max_overhang2_lowend = overhang2;
	  }
	  if (readlength > max_readlength_lowend) {
	    max_readlength_lowend = readlength;
	  }
	} else {
	  if ((overhang2 = readlength - (chrend - chrpos_low)) > max_overhang2_highend) {
	    max_overhang2_highend = overhang2;
	  }
	  if (readlength > max_readlength_highend) {
	    max_readlength_highend = readlength;
	  }
	}
	nreads++;
      }

    } else {
      /* Includes both bounds.  Don't want these reads.  */
    }
  }
  Bamread_unlimit_region(this);
#endif


  if (max_readlength_lowend == 0) {
    npositions_lowend = 0;
  } else if ((npositions_lowend = (chrend - chrstart + 1) + max_overhang1_lowend + max_overhang2_lowend - max_readlength_lowend) < 0) {
    npositions_lowend = 0;
  }

  if (max_readlength_highend == 0) {
    npositions_highend = 0;
  } else if ((npositions_highend = (chrend - chrstart + 1) + max_overhang1_highend + max_overhang2_highend - max_readlength_highend) < 0) {
    npositions_highend = 0;
  }
  
  if ((*npositions = npositions_lowend + npositions_highend) == 0) {
    *npositions = 1;
  }

  return nreads;
}



static void
parse_line (T this, char **acc, unsigned int *flag, int *mapq, char **chr, Genomicpos_T *chrpos,
	    char **mate_chr, Genomicpos_T *mate_chrpos, int *insert_length,
	    Intlist_T *cigartypes, Uintlist_T *cigarlengths, int *cigarlength,
	    int *readlength, char **read, char **quality_string) {
  int type;
  int i;
  unsigned int length;

#ifdef HAVE_SAMTOOLS
  BAM_Sequence_T seq;
  BAM_Cigar_T cigar;
  uint8_t *ptr;

  *acc = bam1_qname(this->bam);
  *flag = this->core->flag;
  *mapq = this->core->qual;

  if (this->core->tid < 0) {
    *chr = (char *) NULL;
  } else if (this->core->tid >= this->ndivs) {
    fprintf(stderr,"tid %d >= ndivs %d\n",this->core->tid,this->ndivs);
    exit(9);
  } else {
    *chr = this->divnames[this->core->tid];
  }
  *chrpos = this->core->pos + 1U;
  /* printf("%s:%u\n",*chr,*chrpos); */

  if (this->core->mtid < 0) {
    *mate_chr = (char *) NULL;
  } else if (this->core->mtid >= this->ndivs) {
    fprintf(stderr,"mtid %d >= ndivs %d\n",this->core->mtid,this->ndivs);
    exit(9);
  } else {
    *mate_chr = this->divnames[this->core->mtid];
  }
  *mate_chrpos = this->core->mpos + 1U;

  *insert_length = this->core->isize;

  *readlength = this->core->l_qseq;
  *read = (char *) CALLOC((*readlength)+1,sizeof(char));
  seq = bam1_seq(this->bam);
  for (i = 0; i < *readlength; i++) {
    switch (bam1_seqi(seq,i)) {
    case 1: (*read)[i] = 'A'; break;
    case 2: (*read)[i] = 'C'; break;
    case 4: (*read)[i] = 'G'; break;
    case 8: (*read)[i] = 'T'; break;
    case 15: (*read)[i] = 'N'; break;
    default: (*read)[i] = '?'; break;
    }
  }

  ptr = bam1_qual(this->bam);
  if (ptr[0] == 0xff) {
    *quality_string = (char *) NULL;
  } else {
    *quality_string = (char *) ptr;
  }

  /* Cigar */
  *cigarlength = 0;
  *cigartypes = (Intlist_T) NULL;
  *cigarlengths = (Uintlist_T) NULL;

  cigar = bam1_cigar(this->bam);
  for (i = 0; i < this->core->n_cigar; i++) {
    length = cigar[i] >> BAM_CIGAR_SHIFT;
    *cigarlengths = Uintlist_push(*cigarlengths,length);

    type = (int) ("MIDNSHP"[cigar[i] & BAM_CIGAR_MASK]);
    *cigartypes = Intlist_push(*cigartypes,type);

    if (type == 'S' || type == 'M' || type == 'I') {
      *cigarlength += (int) length;
    } else if (type == 'H') {
      *cigarlength += (int) length;
    } else if (type == 'D' || type == 'N') {
      /* Ignore */
    } else if (type == 'P') {
      /* Ignore */
    } else {
      fprintf(stderr,"bamread.c cannot cigar int of %d\n",cigar[i] & BAM_CIGAR_MASK);
      exit(9);
    }
  }

  *cigartypes = Intlist_reverse(*cigartypes);
  *cigarlengths = Uintlist_reverse(*cigarlengths);

#endif

  return;
}


int
Bamread_next_line (T this, char **acc, unsigned int *flag, int *mapq, char **chr, Genomicpos_T *chrpos,
		   char **mate_chr, Genomicpos_T *mate_chrpos,
		   Intlist_T *cigartypes, Uintlist_T *cigarlengths, int *cigarlength,
		   int *readlength, char **read, char **quality_string) {
  int insert_length;

#ifndef HAVE_SAMTOOLS
  return 0;
#else
  if (this->region_limited_p == 1) {
    if (bam_iter_read(this->fp,this->iter,this->bam) < 0) {
      return 0;
    } else {
      parse_line(this,&(*acc),&(*flag),&(*mapq),&(*chr),&(*chrpos),
		 &(*mate_chr),&(*mate_chrpos),&insert_length,
		 &(*cigartypes),&(*cigarlengths),&(*cigarlength),
		 &(*readlength),&(*read),&(*quality_string));
      return 1;
    }
  } else {
    if (bam_read1(this->fp,this->bam) < 0) {
      return 0;
    } else {
      parse_line(this,&(*acc),&(*flag),&(*mapq),&(*chr),&(*chrpos),
		 &(*mate_chr),&(*mate_chrpos),&insert_length,
		 &(*cigartypes),&(*cigarlengths),&(*cigarlength),
		 &(*readlength),&(*read),&(*quality_string));
      return 1;
    }
  }
#endif
}



struct Bamline_T {
  char *acc;
  unsigned int flag;
  int nhits;
  bool good_unique_p;			/* Good above good_unique_mapq.  Dependent on second_mapq. */
  int mapq;
  char splice_strand;
  char *chr;
  Genomicpos_T chrpos_low;
  char *mate_chr;
  Genomicpos_T mate_chrpos_low;
  int insert_length;
  Intlist_T cigar_types;
  Uintlist_T cigar_npositions;
  int cigar_querylength;
  int readlength;
  char *read;
  char *quality_string;
  unsigned char *aux_start;
  unsigned char *aux_end;
};


char *
Bamline_acc (Bamline_T this) {
  return this->acc;
}

unsigned int
Bamline_flag (Bamline_T this) {
  return this->flag;
}

static bool
concordantp (unsigned int flag) {
  if (flag & QUERY_UNMAPPED) {
    return false;
  } else if ((flag & PAIRED_READ) == 0U) {
    /* Not even a paired read, so consider it concordant */
    return true;
  } else if (flag & PAIRED_MAPPING) {
    return true;
  } else {
    return false;
  }
}


int
Bamline_concordantp (Bamline_T this) {
  if (this->flag & QUERY_UNMAPPED) {
    return false;
  } else if ((this->flag & PAIRED_READ) == 0U) {
    return false;
  } else if (this->flag & PAIRED_MAPPING) {
    return true;
  } else {
    return false;
  }
}

int
Bamline_lowend_p (Bamline_T this) {
  if (this->mate_chrpos_low > this->chrpos_low) {
    return true;
  } else {
    return false;
  }
}

int
Bamline_firstend_p (Bamline_T this) {
  if ((this->flag & PAIRED_READ) == 0) {
    return true;
  } else if (this->flag & FIRST_READ_P) {
    return true;
  } else if (this->flag & SECOND_READ_P) {
    return false;
  } else {
    fprintf(stderr,"Read is marked as paired (0x1), but neither first read nor second read bit is set\n");
    exit(9);
  }
}

int
Bamline_nhits (Bamline_T this) {
  return this->nhits;
}

bool
Bamline_good_unique_p (Bamline_T this) {
  return this->good_unique_p;
}

int
Bamline_mapq (Bamline_T this) {
  return this->mapq;
}

char *
Bamline_chr (Bamline_T this) {
  return this->chr;
}

Genomicpos_T
Bamline_chrpos_low (Bamline_T this) {
  return this->chrpos_low;
}

char *
Bamline_mate_chr (Bamline_T this) {
  return this->mate_chr;
}

Genomicpos_T
Bamline_mate_chrpos_low (Bamline_T this) {
  return this->mate_chrpos_low;
}

int
Bamline_insert_length (Bamline_T this) {
  return this->insert_length;
}

Intlist_T
Bamline_cigar_types (Bamline_T this) {
  return this->cigar_types;
}


Uintlist_T
Bamline_cigar_npositions (Bamline_T this) {
  return this->cigar_npositions;
}

void
Bamread_print_cigar (Bamline_T this) {
  Intlist_T p;
  Uintlist_T q;

  for (p = this->cigar_types, q = this->cigar_npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    printf("%u%c",Uintlist_head(q),Intlist_head(p));
  }
  return;
}

int
Bamline_cigar_querylength (Bamline_T this) {
  return this->cigar_querylength;
}

int
Bamline_readlength (Bamline_T this) {
  return this->readlength;
}

char *
Bamline_read (Bamline_T this) {
  return this->read;
}

char *
Bamline_quality_string (Bamline_T this) {
  return this->quality_string;
}


static void
aux_print (FILE *fp, unsigned char *s, unsigned char *aux_end) {
  unsigned char tag1, tag2;
  int type;

  while (s < aux_end) {
    tag1 = *s++;
    tag2 = *s++;
    type = *s++;
    fprintf(fp,"\t%c%c:",tag1,tag2);

    if (type == 'c') {
      fprintf(fp,"i:%d",* (int8_t *) s);
      s += 1;
    } else if (type == 'C') {
      fprintf(fp,"i:%u",* (uint8_t *) s);
      s += 1;
    } else if (type == 's') {
      fprintf(fp,"i:%d",* (int16_t *) s);
      s += 2;
    } else if (type == 'S') {
      fprintf(fp,"i:%u",* (uint16_t *) s);
      s += 2;
    } else if (type == 'i') {
      fprintf(fp,"i:%d",* (int32_t *) s);
      s += 4;
    } else if (type == 'I') {
      fprintf(fp,"i:%u",* (uint32_t *) s);
      s += 4;
    } else if (type == 'A') {
      fprintf(fp,"A:%c",* (char *) s);
      s += 1;
    } else if (type == 'f') {
      fprintf(fp,"f:%f",* (float *) s);
      s += 4;
    } else if (type == 'd') {
      fprintf(fp,"d:%f",* (double *) s);
      s += 8;
    } else if (type == 'Z' || type == 'H') {
      fprintf(fp,"Z:");
      while (*s) {
	fprintf(fp,"%c",*s++);
      }
      s++;
    } else {
      /* fprintf(stderr,"Unrecognized type %c\n",type); */
    }
  }

  return;
}


void
Bamline_print (FILE *fp, Bamline_T this, unsigned int newflag) {
  Intlist_T p;
  Uintlist_T q;

  fprintf(fp,"%s\t",this->acc);
  fprintf(fp,"%u\t",newflag);
  if (this->chr == NULL) {
    fprintf(fp,"*\t0\t");
  } else {
    fprintf(fp,"%s\t%u\t",this->chr,this->chrpos_low);
  }
  fprintf(fp,"%d\t",this->mapq);
  for (p = this->cigar_types, q = this->cigar_npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    fprintf(fp,"%u%c",Uintlist_head(q),Intlist_head(p));
  }
  fprintf(fp,"\t");
  if (this->mate_chr == NULL) {
    fprintf(fp,"*\t0\t");
  } else if (this->chr != NULL && strcmp(this->mate_chr,this->chr) == 0) {
    fprintf(fp,"=\t%u\t",this->mate_chrpos_low);
  } else {
    fprintf(fp,"%s\t%u\t",this->mate_chr,this->mate_chrpos_low);
  }
  fprintf(fp,"%d\t",this->insert_length);
  fprintf(fp,"%s\t",this->read);
  if (this->quality_string == NULL) {
    fprintf(fp,"*");
  } else {
    fprintf(fp,"%s",this->quality_string);
  }

  aux_print(fp,this->aux_start,this->aux_end);

  fprintf(fp,"\n");
  return;
}


    
static char
aux_splice_strand (T this) {
  char strand;

#ifndef HAVE_SAMTOOLS
  return ' ';
#else
  uint8_t *s;

  s = bam_aux_get(this->bam,"XS");
  if (s == NULL) {
    return ' ';
  } else if ((strand = bam_aux2A(s)) == '?') {
    return '?';
  } else {
    return strand;
  }
#endif
}

static int
aux_nhits (T this) {
#ifndef HAVE_SAMTOOLS
  return 1;
#else
  uint8_t *s;

  s = bam_aux_get(this->bam,"NH");
  if (s == NULL) {
    return 1;
  } else {
    return bam_aux2i(s);
  }
#endif
}


static bool
aux_good_unique_p (T this, int good_unique_mapq) {
#ifndef HAVE_SAMTOOLS
  return true;
#else
  uint8_t *s;

  s = bam_aux_get(this->bam,"X2");
  if (s != NULL) {
    if (bam_aux2i(s) < good_unique_mapq) {
      return true;
    } else {
      return false;
    }

  } else {
    s = bam_aux_get(this->bam,"NH");
    if (s == NULL) {
      return true;
    } else if (bam_aux2i(s) <= 1) {
      return true;
    } else {
      return false;
    }
  }
#endif
}


static char complCode[128] = COMPLEMENT_LC;


static char
find_strand (bool *canonicalp, char *donor1, char *donor2, char *acceptor1, char *acceptor2,
	     Genomicpos_T firstpos, Genomicpos_T secondpos, char *chr,
	     Genome_T genome, IIT_T chromosome_iit, Bamline_T this, bool trust_sam_p) {
  Chrnum_T chrnum;
  Genomicpos_T chroffset;
  char nt1, nt2, nt3, nt4;
  char truestrand;

  if ((truestrand = this->splice_strand) != ' ') {
    if (trust_sam_p == true) {
      *canonicalp = true;
      *donor1 = *donor2 = *acceptor1 = *acceptor2 = ' ';
      return truestrand;
    } else {
      chrnum = IIT_find_one(chromosome_iit,chr);
      chroffset = Interval_low(IIT_interval(chromosome_iit,chrnum)) - 1U;

      /* Look at genome inside of firstpos and secondpos to get dinucleotides */
      nt1 = Genome_get_char(genome,chroffset+firstpos+1);
      nt2 = Genome_get_char(genome,chroffset+firstpos+2);
      nt3 = Genome_get_char(genome,chroffset+secondpos-2);
      nt4 = Genome_get_char(genome,chroffset+secondpos-1);

      debug(printf("Got splice from %u to %u\n",firstpos,secondpos));
      debug(printf("Dinucleotides are %c%c to %c%c\n",nt1,nt2,nt3,nt4));

      if (truestrand == '+') {
	if (nt1 == 'G' && (nt2 == 'T' || nt2 == 'C') && nt3 == 'A' && nt4 == 'G') {
	  *canonicalp = true;
	} else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'A' && nt4 == 'C') {
	  *canonicalp = true;
	} else {
	  *canonicalp = false;
	}
	*donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;

      } else if (truestrand == '-') {
	if (nt1 == 'C' && nt2 == 'T' && (nt3 == 'A' || nt3 == 'G') && nt4 == 'C') {
	  *canonicalp = true;
	} else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'A' && nt4 == 'T') {
	  *canonicalp = true;
	} else {
	  *canonicalp = false;
	}
	*donor1 = complCode[(int) nt4]; *donor2 = complCode[(int) nt3]; *acceptor1 = complCode[(int) nt2]; *acceptor2 = complCode[(int) nt1];

      } else {
	fprintf(stderr,"Unrecognized truestrand %c\n",truestrand);
	abort();
      }

      return truestrand;
    }

  } else if (chromosome_iit == NULL) {
    fprintf(stderr,"Strand is not present in auxinfo\n");
    fprintf(stderr,"To determine strand, need to provide index file with -d flag\n");
    exit(9);

  } else {
    chrnum = IIT_find_one(chromosome_iit,chr);
    chroffset = Interval_low(IIT_interval(chromosome_iit,chrnum)) - 1U;

    /* Look at genome inside of firstpos and secondpos to determine truestrand */
    nt1 = Genome_get_char(genome,chroffset+firstpos+1);
    nt2 = Genome_get_char(genome,chroffset+firstpos+2);
    nt3 = Genome_get_char(genome,chroffset+secondpos-2);
    nt4 = Genome_get_char(genome,chroffset+secondpos-1);

    debug(printf("Got splice from %u to %u\n",firstpos,secondpos));
    debug(printf("Dinucleotides are %c%c to %c%c\n",nt1,nt2,nt3,nt4));

    if (nt1 == 'G' && (nt2 == 'T' || nt2 == 'C') && nt3 == 'A' && nt4 == 'G') {
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = true;
      return '+';
    } else if (nt1 == 'C' && nt2 == 'T' && (nt3 == 'A' || nt3 == 'G') && nt4 == 'C') {
      *donor1 = complCode[(int) nt4]; *donor2 = complCode[(int) nt3]; *acceptor1 = complCode[(int) nt2]; *acceptor2 = complCode[(int) nt1];
      *canonicalp = true;
      return '-';
    } else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'A' && nt4 == 'C') {
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = true;
      return '+';
    } else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'A' && nt4 == 'T') {
      *donor1 = complCode[(int) nt4]; *donor2 = complCode[(int) nt3]; *acceptor1 = complCode[(int) nt2]; *acceptor2 = complCode[(int) nt1];
      *canonicalp = true;
      return '-';
    } else {
      /* In GSNAP, will want to output sense information in SAM output. */
#if 0
      fprintf(stderr,"Splice %s:%u..%u is not (semi-)canonical: %c%c...%c%c.  Cannot determine sense.\n",
	      chr,firstpos,secondpos,nt1,nt2,nt3,nt4);
#endif
      *donor1 = nt1; *donor2 = nt2; *acceptor1 = nt3; *acceptor2 = nt4;
      *canonicalp = false;
      return ' ';
    }
  }
}



char
Bamline_strand (Bamline_T this, Genome_T genome, IIT_T chromosome_iit) {
  char strand = ' ';
  Genomicpos_T chrpos, firstpos, secondpos;
  Intlist_T p;
  Uintlist_T q;
  int type;
  bool canonicalp;
  char donor1, donor2, acceptor1, acceptor2;


  chrpos = this->chrpos_low;
  for (p = this->cigar_types, q = this->cigar_npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      chrpos += Uintlist_head(q);

    } else if (type == 'N') {
      firstpos = chrpos - 1U;
      chrpos += Uintlist_head(q);
      secondpos = chrpos;

      if (strand == ' ') {
	strand = find_strand(&canonicalp,&donor1,&donor2,&acceptor1,&acceptor2,
			     firstpos,secondpos,this->chr,genome,chromosome_iit,
			     this,/*trust_sam_p*/true);
	return strand;
      }

    } else if (type == 'P') {
      /* Do nothing */

    } else if (type == 'I') {
      /* Do nothing */

    } else if (type == 'D') {
      /* CHECK */
      chrpos += Uintlist_head(q);

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }
    debug(printf("  type = %c, chrpos = %u\n",type,chrpos));
  }

  return ' ';
}




Genomicpos_T
Bamline_chrpos_high (Bamline_T this) {
  Intlist_T p;
  Uintlist_T q;
  Genomicpos_T chrpos_high;
  int type;

  chrpos_high = this->chrpos_low;
  for (p = this->cigar_types, q = this->cigar_npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      chrpos_high += Uintlist_head(q);

    } else if (type == 'N') {
      chrpos_high += Uintlist_head(q);

    } else if (type == 'P') {
      /* Do nothing */

    } else if (type == 'I') {
      /* Do nothing */

    } else if (type == 'D') {
      /* CHECK */
      chrpos_high += Uintlist_head(q);

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }
    debug(printf("  type = %c, chrpos = %u\n",type,chrpos_high));
  }

  return chrpos_high - 1U;
}



static void
Bamline_regions (Uintlist_T *chrpos_lows, Uintlist_T *chrpos_highs, Bamline_T this) {
  Intlist_T p;
  Uintlist_T q;
  Genomicpos_T position;
  int type;

  position = this->chrpos_low;
  for (p = this->cigar_types, q = this->cigar_npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      *chrpos_lows = Uintlist_push(*chrpos_lows,position);
      position += Uintlist_head(q);
      *chrpos_highs = Uintlist_push(*chrpos_highs,position);

    } else if (type == 'N') {
      position += Uintlist_head(q);

    } else if (type == 'P') {
      /* Do nothing */

    } else if (type == 'I') {
      /* Do nothing */

    } else if (type == 'D') {
      /* CHECK */
      position += Uintlist_head(q);

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }
  }

  return;
}


static void
Bamline_splices (Uintlist_T *splice_lows, Uintlist_T *splice_highs, Intlist_T *splice_signs,
		 Bamline_T this) {
  Intlist_T p;
  Uintlist_T q;
  Genomicpos_T position;
  int type;

  position = this->chrpos_low;
  for (p = this->cigar_types, q = this->cigar_npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      position += Uintlist_head(q);

    } else if (type == 'N') {
      *splice_lows = Uintlist_push(*splice_lows,position);
      position += Uintlist_head(q);
      *splice_highs = Uintlist_push(*splice_highs,position);

      if (this->splice_strand == '+') {
	*splice_signs = Intlist_push(*splice_signs,+1);
      } else if (this->splice_strand == '-') {
	*splice_signs = Intlist_push(*splice_signs,-1);
      } else {
	*splice_signs = Intlist_push(*splice_signs,0);
      }

    } else if (type == 'P') {
      /* Do nothing */

    } else if (type == 'I') {
      /* Do nothing */

    } else if (type == 'D') {
      /* CHECK */
      position += Uintlist_head(q);

    } else {
      fprintf(stderr,"Cannot parse type %c\n",type);
      exit(9);
    }
  }

  return;
}




void
Bamline_free (Bamline_T *old) {
  if (*old) {
    FREE((*old)->acc);
    if ((*old)->chr != NULL) {
      FREE((*old)->chr);
    }
    if ((*old)->mate_chr != NULL) {
      FREE((*old)->mate_chr);
    }
    Intlist_free(&(*old)->cigar_types);
    Uintlist_free(&(*old)->cigar_npositions);
    FREE((*old)->read);
    if ((*old)->quality_string != NULL) {
      FREE((*old)->quality_string);
    }
    FREE(*old);
  }

  return;
}


static Bamline_T
Bamline_new (char *acc, unsigned int flag, int nhits, bool good_unique_p, int mapq,
	     char splice_strand, char *chr, Genomicpos_T chrpos_low,
	     char *mate_chr, Genomicpos_T mate_chrpos_low, int insert_length,
	     Intlist_T cigar_types, Uintlist_T cigar_npositions, int cigar_querylength, int readlength,
	     char *read, char *quality_string, unsigned char *aux_start, unsigned char *aux_end) {
  Bamline_T new = (Bamline_T) MALLOC(sizeof(*new));

  new->acc = (char *) CALLOC(strlen(acc)+1,sizeof(char));
  strcpy(new->acc,acc);

  new->flag = flag;
  new->nhits = nhits;
  new->good_unique_p = good_unique_p;

  new->mapq = mapq;
  
  new->splice_strand = splice_strand;
  if (chr == NULL) {
    new->chr = (char *) NULL;
  } else {
    new->chr = (char *) CALLOC(strlen(chr)+1,sizeof(char));
    strcpy(new->chr,chr);
  }
  new->chrpos_low = chrpos_low;

  if (mate_chr == NULL) {
    new->mate_chr = (char *) NULL;
  } else {
    new->mate_chr = (char *) CALLOC(strlen(mate_chr)+1,sizeof(char));
    strcpy(new->mate_chr,mate_chr);
  }
  new->mate_chrpos_low = mate_chrpos_low;

  new->insert_length = insert_length;

  new->cigar_types = cigar_types;		/* not copying */
  new->cigar_npositions = cigar_npositions;	/* not copying */

  new->cigar_querylength = cigar_querylength;
  new->readlength = readlength;

  /* No need to copy, since it is allocated already */
  new->read = read;

  if (quality_string == NULL) {
    new->quality_string = NULL;
  } else {
    new->quality_string = (char *) CALLOC(readlength+1,sizeof(char));
    strncpy(new->quality_string,quality_string,readlength);
  }

  new->aux_start = aux_start;
  new->aux_end = aux_end;

  return new;
}


Bamline_T
Bamread_next_bamline (T this, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		      bool need_unique_p, bool need_primary_p, bool need_concordant_p) {
  char *acc, *chr, *mate_chr, splice_strand;
  unsigned int flag;
  int nhits;
  bool good_unique_p;
  int mapq;
  Genomicpos_T chrpos_low, mate_chrpos_low;
  int insert_length;
  Intlist_T cigar_types;
  Uintlist_T cigarlengths;
  int cigar_querylength, readlength;
  char *read;
  char *quality_string;

#ifndef HAVE_SAMTOOLS
  return (Bamline_T) NULL;
#else
  if (this->region_limited_p == 1) {
    debug1(fprintf(stderr,"Region limited\n"));
    while (bam_iter_read(this->fp,this->iter,this->bam) >= 0) {
      debug1(fprintf(stderr,"Got line\n"));
      parse_line(this,&acc,&flag,&mapq,&chr,&chrpos_low,
		 &mate_chr,&mate_chrpos_low,&insert_length,
		 &cigar_types,&cigarlengths,&cigar_querylength,
		 &readlength,&read,&quality_string);
      if (mapq < minimum_mapq) {
	debug1(fprintf(stderr,"Fails because mapq %d < minimum %d\n",mapq,minimum_mapq));
	Intlist_free(&cigar_types);
	Uintlist_free(&cigarlengths);
	FREE(read);
      } else if ((nhits = aux_nhits(this)) > maximum_nhits) {
	debug1(fprintf(stderr,"Fails because nhits %d > maximum %d\n",nhits,maximum_nhits));
	Intlist_free(&cigar_types);
	Uintlist_free(&cigarlengths);
	FREE(read);
      } else if (need_unique_p == true && nhits > 1) {
	debug1(fprintf(stderr,"Fails because need unique and nhits %d\n",nhits));
	Intlist_free(&cigar_types);
	Uintlist_free(&cigarlengths);
	FREE(read);
      } else if (need_primary_p == true && (flag & NOT_PRIMARY) != 0) {
	debug1(fprintf(stderr,"Fails because need primary and flag is %u\n",flag));
	Intlist_free(&cigar_types);
	Uintlist_free(&cigarlengths);
	FREE(read);
      } else if (need_concordant_p == true && concordantp(flag) == false) {
	debug1(fprintf(stderr,"Fails because need concordant and flag is %u\n",flag));
	Intlist_free(&cigar_types);
	Uintlist_free(&cigarlengths);
	FREE(read);
      } else {
	debug1(fprintf(stderr,"Success\n"));
	splice_strand = aux_splice_strand(this);
	good_unique_p = aux_good_unique_p(this,good_unique_mapq);
	return Bamline_new(acc,flag,nhits,good_unique_p,mapq,splice_strand,chr,chrpos_low,
			   mate_chr,mate_chrpos_low,insert_length,
			   cigar_types,cigarlengths,cigar_querylength,
			   readlength,read,quality_string,
			   /*aux_start*/bam1_aux(this->bam),
			   /*aux_end*/this->bam->data + this->bam->data_len);
      }
    }
    return (Bamline_T) NULL;

  } else {
    while (bam_read1(this->fp,this->bam) >= 0) {
      parse_line(this,&acc,&flag,&mapq,&chr,&chrpos_low,
		 &mate_chr,&mate_chrpos_low,&insert_length,
		 &cigar_types,&cigarlengths,&cigar_querylength,
		 &readlength,&read,&quality_string);
      if (mapq >= minimum_mapq &&
	  (nhits = aux_nhits(this)) <= maximum_nhits &&
	  (need_unique_p == false || nhits == 1) &&
	  (need_primary_p == false || (flag & NOT_PRIMARY) == 0) &&
	  (need_concordant_p == false || concordantp(flag) == true)) {
	splice_strand = aux_splice_strand(this);
	good_unique_p = aux_good_unique_p(this,good_unique_mapq);
	return Bamline_new(acc,flag,nhits,good_unique_p,mapq,splice_strand,chr,chrpos_low,
			   mate_chr,mate_chrpos_low,insert_length,
			   cigar_types,cigarlengths,cigar_querylength,
			   readlength,read,quality_string,
			   /*aux_start*/bam1_aux(this->bam),
			   /*aux_end*/this->bam->data + this->bam->data_len);
      } else {
	Intlist_free(&cigar_types);
	Uintlist_free(&cigarlengths);
	FREE(read);
      }
    }
    return (Bamline_T) NULL;

  }
#endif

}


/************************************************************************
 *   Bamstore
 ************************************************************************/

struct Bamstore_T {
  Genomicpos_T chrpos;
  List_T bamlines;
};

void
Bamstore_free (Bamstore_T *old) {
  List_T p;
  Bamline_T bamline;

  for (p = (*old)->bamlines; p != NULL; p = List_next(p)) {
    bamline = (Bamline_T) List_head(p);
    Bamline_free(&bamline);
  }
  List_free(&(*old)->bamlines);

  FREE(*old);
  return;
}


Bamstore_T
Bamstore_new (Genomicpos_T chrpos) {
  Bamstore_T new = (Bamstore_T) MALLOC(sizeof(*new));

  new->chrpos = chrpos;
  new->bamlines = (List_T) NULL;

  return new;
}



Bamline_T
Bamstore_get (Table_T bamstore_chrtable, char *chr, Genomicpos_T low, char *acc,
	      Genomicpos_T mate_low) {
  List_T p, list = NULL;
  Bamline_T wanted = NULL, bamline;
  Bamstore_T bamstore;
  Uinttable_T bamstore_table;
  Chrom_T chrom;

  chrom = Chrom_from_string(chr,/*mitochondrial_string*/NULL,/*order*/0U);
  if ((bamstore_table = (Uinttable_T) Table_get(bamstore_chrtable,(void *) chrom)) == NULL) {
    fprintf(stderr,"Unexpected error.  No bamstore_table for chr %s\n",chr);
    Chrom_free(&chrom);
    return (Bamline_T) NULL;
  } else {
    Chrom_free(&chrom);
  }

  if ((bamstore = (Bamstore_T) Uinttable_get(bamstore_table,low)) == NULL) {
    /* May have been excluded for other reasons */
    return (Bamline_T) NULL;
  } else {
    for (p = bamstore->bamlines; p != NULL; p = List_next(p)) {
      bamline = (Bamline_T) List_head(p);
      if (strcmp(Bamline_acc(bamline),acc) == 0 && Bamline_mate_chrpos_low(bamline) == mate_low) {
	wanted = bamline;
      } else {
	list = List_push(list,(void *) bamline);
      }
    }

    List_free(&bamstore->bamlines);
    bamstore->bamlines = list;
    if (list == NULL) {
      Uinttable_remove(bamstore_table,low);
      Bamstore_free(&bamstore);
    }

    return wanted;
  }
}
  

void
Bamstore_add_at_low (Table_T bamstore_chrtable, char *chr, Genomicpos_T low,
		     Bamline_T bamline) {
  Bamstore_T bamstore;
  Uinttable_T bamstore_table;
  Chrom_T chrom;

  chrom = Chrom_from_string(chr,/*mitochondrial_string*/NULL,/*order*/0U);
  if ((bamstore_table = (Uinttable_T) Table_get(bamstore_chrtable,(void *) chrom)) == NULL) {
    bamstore_table = Uinttable_new(65522); /* estimate 65522 splice sites per chromosome */
    Table_put(bamstore_chrtable,(void *) chrom,(void *) bamstore_table);
  } else {
    Chrom_free(&chrom);
  }

  if ((bamstore = (Bamstore_T) Uinttable_get(bamstore_table,low)) == NULL) {
    bamstore = Bamstore_new(low);
    Uinttable_put(bamstore_table,low,(void *) bamstore);
  }

  bamstore->bamlines = List_push(bamstore->bamlines,(void *) bamline);
  return;
}


void
Bamstore_table_free (Uinttable_T *bamstore_table) {
  Genomicpos_T *keys;
  int n, i;
  Bamstore_T bamstore;

  if ((n = Uinttable_length(*bamstore_table)) > 0) {
    keys = (Genomicpos_T *) Uinttable_keys(*bamstore_table,/*sortp*/false);
    for (i = 0; i < n; i++) {
      bamstore = Uinttable_get(*bamstore_table,keys[i]);
      if (bamstore == NULL) {
	fprintf(stderr,"key is %u, value is NULL\n",keys[i]);
	abort();
      } else {
	Bamstore_free(&bamstore);
      }
    }
    FREE(keys);
  }

  return;
}


/************************************************************************
 *   Retrieving and assigning levels to all reads in a region
 ************************************************************************/

struct Bampair_T {
  Bamline_T bamline_low;
  Bamline_T bamline_high;
  Genomicpos_T chrpos_low;
  Genomicpos_T chrpos_high;
  int level;
};

Genomicpos_T
Bampair_chrpos_low (Bampair_T this) {
  return this->chrpos_low;
}

Genomicpos_T
Bampair_chrpos_high (Bampair_T this) {
  return this->chrpos_high;
}

int
Bampair_level (Bampair_T this) {
  return this->level;
}

bool
Bampair_good_unique_p (Bampair_T this) {
  if (this->bamline_low != NULL && this->bamline_low->good_unique_p == false) {
    return false;
  }
  if (this->bamline_high != NULL && this->bamline_high->good_unique_p == false) {
    return false;
  }
  return true;
}

bool
Bampair_uniquep (Bampair_T this) {
  if (this->bamline_low != NULL && this->bamline_low->nhits > 1) {
    return false;
  }
  if (this->bamline_high != NULL && this->bamline_high->nhits > 1) {
    return false;
  }
  return true;
}

bool
Bampair_primaryp (Bampair_T this) {
  if (this->bamline_low != NULL && (this->bamline_low->flag & NOT_PRIMARY) != 0) {
    return false;
  }
  if (this->bamline_high != NULL && (this->bamline_high->flag & NOT_PRIMARY) != 0) {
    return false;
  }
  return true;
}



static Bampair_T
Bampair_new (Bamline_T bamline_low, Bamline_T bamline_high) {
  Bampair_T new = (Bampair_T) MALLOC(sizeof(*new));

  new->bamline_low = bamline_low;
  new->bamline_high = bamline_high;
  if (bamline_low == NULL) {
    new->chrpos_low = bamline_high->chrpos_low;
    new->chrpos_high = Bamline_chrpos_high(bamline_high);
  } else if (bamline_high == NULL) {
    new->chrpos_low = bamline_low->chrpos_low;
    new->chrpos_high = Bamline_chrpos_high(bamline_low);
  } else {
    new->chrpos_low = bamline_low->chrpos_low;
    new->chrpos_high = Bamline_chrpos_high(bamline_high);
  }
  new->level = -1;
  return new;
}


void
Bampair_free (Bampair_T *old) {
  if (*old) {
    Bamline_free(&(*old)->bamline_low);
    Bamline_free(&(*old)->bamline_high);
    FREE(*old);
  }
  return;
}


void
Bampair_print (FILE *fp, Bampair_T this) {
  if (this->bamline_low != NULL) {
    Bamline_print(fp,this->bamline_low,this->bamline_low->flag);
  }
  if (this->bamline_high != NULL) {
    Bamline_print(fp,this->bamline_high,this->bamline_high->flag);
  }
  return;
}


void
Bampair_details (Uintlist_T *chrpos_lows, Uintlist_T *chrpos_highs,
		 Uintlist_T *splice_lows, Uintlist_T *splice_highs, Intlist_T *splice_signs,
		 Bampair_T this) {
  *chrpos_lows = (Uintlist_T) NULL;
  *chrpos_highs = (Uintlist_T) NULL;
  *splice_lows = (Uintlist_T) NULL;
  *splice_highs = (Uintlist_T) NULL;
  *splice_signs = (Intlist_T) NULL;

  if (this->bamline_low != NULL) {
    Bamline_regions(&(*chrpos_lows),&(*chrpos_highs),this->bamline_low);
    Bamline_splices(&(*splice_lows),&(*splice_highs),&(*splice_signs),this->bamline_low);
  }

  if (this->bamline_high != NULL) {
    Bamline_regions(&(*chrpos_lows),&(*chrpos_highs),this->bamline_high);
    Bamline_splices(&(*splice_lows),&(*splice_highs),&(*splice_signs),this->bamline_high);
  }

  *chrpos_lows = Uintlist_reverse(*chrpos_lows);
  *chrpos_highs = Uintlist_reverse(*chrpos_highs);
  *splice_lows = Uintlist_reverse(*splice_lows);
  *splice_highs = Uintlist_reverse(*splice_highs);
  *splice_signs = Intlist_reverse(*splice_signs);

  return;
}


List_T
Bamread_all_pairs (T bamreader, int minimum_mapq, int good_unique_mapq, int maximum_nhits,
		   bool need_unique_p, bool need_primary_p, bool need_concordant_p) {
  List_T lines = NULL;
  Bamline_T bamline_low, bamline;

  Table_T bamstore_chrtable;
  Uinttable_T bamstore_table;
  Chrom_T *chroms, chrom;
  int n, i;


  bamstore_chrtable = Table_new(100,Chrom_compare_table,Chrom_hash_table);

  while ((bamline = Bamread_next_bamline(bamreader,minimum_mapq,good_unique_mapq,maximum_nhits,
					 need_unique_p,need_primary_p,need_concordant_p)) != NULL) {
    if (Bamline_concordantp(bamline) == false) {
      /* Handle now */
      if (Bamline_firstend_p(bamline) == true) {
	lines = List_push(lines,Bampair_new(/*bamline_low*/bamline,/*bamline_high*/NULL));
      } else {
	lines = List_push(lines,Bampair_new(/*bamline_low*/NULL,/*bamline_high*/bamline));
      }

    } else if (Bamline_lowend_p(bamline) == true) {
      /* Wait for high end */
      Bamstore_add_at_low(bamstore_chrtable,Bamline_chr(bamline),Bamline_chrpos_low(bamline),
			  bamline);

    } else {
      /* This is the high end */
      bamline_low = Bamstore_get(bamstore_chrtable,Bamline_chr(bamline),Bamline_mate_chrpos_low(bamline),
				 Bamline_acc(bamline),Bamline_chrpos_low(bamline));
      if (bamline_low == NULL) {
	fprintf(stderr,"Hmm...low end not found for %s at %s:%u\n",
		Bamline_acc(bamline),Bamline_chr(bamline),Bamline_chrpos_low(bamline));
      } else {
	lines = List_push(lines,Bampair_new(bamline_low,/*bamline_high*/bamline));
      }
    }
  }

  if ((n = Table_length(bamstore_chrtable)) > 0) {
    chroms = (Chrom_T *) Table_keys(bamstore_chrtable,NULL);
    for (i = 0; i < n; i++) {
      chrom = chroms[i];
      bamstore_table = (Uinttable_T) Table_get(bamstore_chrtable,(void *) chrom);
      Bamstore_table_free(&bamstore_table);
      Uinttable_free(&bamstore_table);
    }
    for (i = 0; i < n; i++) {
      Chrom_free(&(chroms[i]));
    }
    FREE(chroms);
  }
  Table_free(&bamstore_chrtable);


  return List_reverse(lines);
}


static int
level_cmp (const void *x, const void *y) {
  Bampair_T a = * (Bampair_T *) x;
  Bampair_T b = * (Bampair_T *) y;

  if (a->chrpos_low < b->chrpos_low) {
    return -1;
  } else if (a->chrpos_low > b->chrpos_low) {
    return +1;
  } else if (a->chrpos_high < b->chrpos_high) {
    return -1;
  } else if (a->chrpos_high > b->chrpos_high) {
    return +1;
#if 0
  } else if (a->sign > b->sign) {
    return -1;
  } else if (a->sign < b->sign) {
    return +1;
#endif

  } else {
    return 0;
  }
}



int
Bampair_compute_levels (List_T bampairs, Genomicpos_T mincoord,
			Genomicpos_T maxcoord, int max_allowed_levels,
			double xfactor, Genomicpos_T min_pairlength) {
  int nbampairs, i;
  int maxlevel = -1, level;
  bool donep;
  Bampair_T *array, bampair;
  double *rightmost, xlow;

  if ((nbampairs = List_length(bampairs)) > 0) {
    array = (Bampair_T *) List_to_array(bampairs,NULL);
    qsort(array,nbampairs,sizeof(T),level_cmp);

    rightmost = (double *) CALLOC(max_allowed_levels,sizeof(double));
    for (i = 0; i < max_allowed_levels; i++) {
      rightmost[i] = 0.0;
    }

    for (i = 0; i < nbampairs; i++) {
      bampair = array[i];
      if (bampair->chrpos_high - bampair->chrpos_low < min_pairlength) {
	bampair->level = -1;
      } else {
	/* Find appropriate level */
	level = 0;
	donep = false;
	while (level < max_allowed_levels && !donep) {
	  xlow = xfactor * bampair->chrpos_low;
	  if (level > maxlevel) {
	    donep = true;
	    maxlevel = level;
	  } else if (rightmost[level] < xlow) {
	    donep = true;
	  } else {
	    level++;
	  }
	}

	if (level < max_allowed_levels) {
	  rightmost[level] = xfactor * (bampair->chrpos_high + 10);

#if 0
	  if (bampair->chrpos_high < mincoord || bampair->chrpos_low > maxcoord) {
	    /* Skip printing if both ends outside of region */
	  } else if (bampair->chrpos_low < mincoord || bampair->chrpos_high > maxcoord) {
	    /* Skip printing if either end outside of region */
	  } else {
	    bampair->level = level;
	  }
#else
	  bampair->level = level;
#endif
	}
      }
    }

    FREE(rightmost);
    FREE(array);
  }

  return maxlevel + 1;
}
