static char rcsid[] = "$Id: samread.c 219290 2019-05-21 01:14:10Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>

#include "samread.h"
#include "except.h"
#include "mem.h"
#include "assert.h"
#include "bool.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


char *
Samread_get_acc (unsigned int *flag, char *line) {
  char *acc, *p;
  int length;

  p = line;
  while (*p != '\0' && *p != '\t') p++;
  length = (p - line)/sizeof(char);
  acc = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(acc,line,length);
  acc[length] = '\0';

  if (*p == '\0') {
    fprintf(stderr,"Can't parse flag part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  *flag = strtoul(p,NULL,10);

  return acc;
}


/* ILLUMINA-A1CCE9_0004:1:1:1103:6310#0	0	20	33639850	255	55M21S	*	0	0	AAAAATTGTATACCGCAGATTCAGGCATGGATTCCGTGAAGGAACAACACCTAAANCCAAAGNTCGGAAGANCGGN	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDCCCCBCCBDBDCBCCDCDC@CCC&AAAA################	NM:i:2 */

char *
Samread_parse_line (char **acc, unsigned int *flag, int *mapq, char **chr, Genomicpos_T *chrpos, char **cigar,
		    char **mate_chr, Genomicpos_T *mate_chrpos_low, int *readlength, char **read, char **quality_string,
		    char **hardclip, char **hardclip_quality, char *line) {
  char *auxinfo, *p, *q, *r;
  char tag1, tag2;
  int length, i;

  debug(printf("Entering Samread_parse_line with %s\n",line));

  p = line;
  while (!isspace(*p)) p++;
  length = (p - line)/sizeof(char);
  *acc = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(*acc,line,length);
  (*acc)[length] = '\0';

  if (*p != '\0') {		/* Skip over tab */
    p++;
  }

  if (sscanf(p,"%u",&(*flag)) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",p);
    abort();
  } else {
    debug(printf("  flag = %u\n",*flag));
  }

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *chr = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(*chr,p,length);
  (*chr)[length] = '\0';

  debug(printf("  chr = %s\n",*chr));
  if (*q != '\0') {
    q++;
  }


  p = q;
  if (sscanf(p,"%u",&(*chrpos)) != 1) {
    fprintf(stderr,"Unable to find chrpos in %s\n",p);
    abort();
  } else {
    debug(printf("  chrpos = %u\n",*chrpos));
  }


  while (!isspace(*p)) p++;	/* Skip over chrpos */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chrpos part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }

  /* Read mapping quality */
  if (sscanf(p,"%d",&(*mapq)) != 1) {
    fprintf(stderr,"Unable to find mapq in %s\n",p);
    abort();
  } else {
    debug(printf("  mapq = %d\n",*mapq));
  }

  /* Skip past mapping quality */
  while (!isspace(*p)) p++;


  if (*p == '\0') {
    fprintf(stderr,"Can't parse cigar part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *cigar = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(*cigar,p,length);
  (*cigar)[length] = '\0';

  debug(printf("  cigar = %s\n",*cigar));
  

  /* mate chr */
  p = q;
  if (*p != '\0') {
    p++;			/* Should be a tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *mate_chr = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(*mate_chr,p,length);
  (*mate_chr)[length] = '\0';

  debug(printf("  mate_chr = %s\n",*mate_chr));
  if (*q == '\0') {
    fprintf(stderr,"Can't parse mate chr part of %s\n",line);
    abort();
  } else {
    q++;
  }

  /* mate chrpos low */
  p = q;
  if (sscanf(p,"%u",&(*mate_chrpos_low)) != 1) {
    fprintf(stderr,"Unable to find mate_chrpos_low in %s\n",p);
    abort();
  } else {
    debug(printf("  mate_chrpos_low = %u\n",*mate_chrpos_low));
  }

  while (!isspace(*p)) p++;	/* Skip over mate_chrpos */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse mate chrpos part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }


  /* Skip over insert size */
  while (!isspace(*p)) p++;
  if (*p == '\0') {
    fprintf(stderr,"Can't parse mate chrpos part of %s\n",line);
    abort();
  } else {
    p++;
  }


  q = p;
  while (!isspace(*q)) q++;
  *readlength = (q - p)/sizeof(char);
  if (*q == '\t') q++;
  debug(printf("  readlength = %d\n",*readlength));

  *read = (char *) MALLOC(((*readlength)+1)*sizeof(char));
  strncpy(*read,p,*readlength);
  (*read)[*readlength] = '\0';

  debug(printf("  read = %s\n",*read));

  p = q;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *quality_string = (char *) MALLOC(((*readlength)+1)*sizeof(char));
  if (length == *readlength) {
    strncpy(*quality_string,p,length);
    (*quality_string)[length] = '\0';

  } else {
    for (i = 0; i < *readlength; i++) {
      (*quality_string)[i] = ' ';
    }
    (*quality_string)[length] = '\0';
  }

  if (*q == '\t') q++;
  p = auxinfo = q;

  *hardclip = (char *) NULL;
  *hardclip_quality = (char *) NULL;
  while (*p != '\0' && *p != '\n') {
    tag1 = p[0];
    tag2 = p[1];

    if (tag1 == 'X' && tag2 == 'H') {
      debug(printf("Found tag XH\n"));
      /* XH:Z: */
      p += 5;

      r = p;
      while (!isspace(*r)) r++;
      length = (r - p)/sizeof(char);
      *hardclip = (char *) MALLOC((length+1) * sizeof(char));
      strncpy(*hardclip,p,length);
      (*hardclip)[length] = '\0';

      p = r;
      if (*p == '\t') {
	p++;
      }

    } else if (tag1 == 'X' && tag2 == 'I') {
      debug(printf("Found tag XI\n"));
      /* XI:Z: */
      p += 5;

      r = p;
      while (!isspace(*r)) r++;
      length = (r - p)/sizeof(char);
      *hardclip_quality = (char *) MALLOC((length+1) * sizeof(char));
      strncpy(*hardclip_quality,p,length);
      (*hardclip_quality)[length] = '\0';

      p = r;
      if (*p == '\t') {
	p++;
      }

    } else {
      while (*p != '\0' && *p != '\t') {
	p++;
      }
      if (*p == '\t') {
	p++;
      }
    }
  }

  return auxinfo;
}


char *
Samread_chr (char *line) {
  char *chr;
  unsigned int flag;
  int mapq;

  char *p, *q;
  int length;

  debug(printf("Entering Samread_chr with %s\n",line));

  p = line;
  while (!isspace(*p)) p++;
  length = (p - line)/sizeof(char);
#if 0
  *acc = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(*acc,line,length);
  (*acc)[length] = '\0';
#endif

  if (*p != '\0') {		/* Skip over tab */
    p++;
  }

  if (sscanf(p,"%u",&flag) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",p);
    abort();
  } else {
    debug(printf("  flag = %u\n",flag));
  }

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  chr = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(chr,p,length);
  chr[length] = '\0';

  debug(printf("  chr = %s\n",chr));
  if (*q != '\0') {
    q++;
  }

  return chr;
}


char *
Samread_chrinfo (Genomicpos_T *chrpos, char **cigar, char *line) {
  char *chr;
  unsigned int flag;
  int mapq;

  char *p, *q;
  int length;

  debug(printf("Entering Samread_chrinfo with %s\n",line));

  p = line;
  while (!isspace(*p)) p++;
  length = (p - line)/sizeof(char);
#if 0
  *acc = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(*acc,line,length);
  (*acc)[length] = '\0';
#endif

  if (*p != '\0') {		/* Skip over tab */
    p++;
  }

  if (sscanf(p,"%u",&flag) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",p);
    abort();
  } else {
    debug(printf("  flag = %u\n",flag));
  }

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  chr = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(chr,p,length);
  chr[length] = '\0';

  debug(printf("  chr = %s\n",chr));
  if (*q != '\0') {
    q++;
  }


  p = q;
  if (sscanf(p,"%u",&(*chrpos)) != 1) {
    fprintf(stderr,"Unable to find chrpos in %s\n",p);
    abort();
  } else {
    debug(printf("  chrpos = %u\n",*chrpos));
  }

  while (!isspace(*p)) p++;	/* Skip over chrpos */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chrpos part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }

  /* Read mapping quality */
  if (sscanf(p,"%d",&mapq) != 1) {
    fprintf(stderr,"Unable to find mapq in %s\n",p);
    abort();
  } else {
    debug(printf("  mapq = %d\n",mapq));
  }

  /* Skip past mapping quality */
  while (!isspace(*p)) p++;


  if (*p == '\0') {
    fprintf(stderr,"Can't parse cigar part of %s\n",line);
    abort();
  } else {
    p++;			/* Skip over tab */
  }
  q = p;
  while (!isspace(*q)) q++;
  length = (q - p)/sizeof(char);
  *cigar = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(*cigar,p,length);
  (*cigar)[length] = '\0';

  debug(printf("  cigar = %s\n",*cigar));

  return chr;
}


void
Samread_print_altered_single (FILE *fp, unsigned int newflag, int mapq, char *line) {
  char *p;
  int i;

  p = line;

  /* Print acc */
  while (!isspace(*p)) {
    putc(*p,fp);
    p++;
  }

  if (*p != '\0') {		/* Skip over tab */
    putc(*p,fp);
    p++;
  }

  fprintf(fp,"%u",newflag);

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);
    p++;			/* Skip over tab */
  }

  /* Print chr, chrpos */
  for (i = 0; i < 2; i++) {
    while (!isspace(*p)) {
      putc(*p,fp);			/* Print field */
      p++;
    }
    if (*p != '\0') {
      putc(*p,fp);
      p++;
    }
  }

  /* MAPQ */
  fprintf(fp,"%d",mapq);

  while (!isspace(*p)) p++;	/* Skip over mapq */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);
    p++;			/* Skip over tab */
  }


  /* print rest of line */
  fprintf(fp,"%s",p);

  return;
}

void
Samread_print_altered_paired (FILE *fp, unsigned int newflag, int mapq, char *mate_chr, Genomicpos_T mate_chrpos,
			      int insert_length, char *line) {
  char *p;
  int i;

  p = line;

  /* Print acc */
  while (!isspace(*p)) {
    putc(*p,fp);
    p++;
  }

  if (*p != '\0') {		/* Skip over tab */
    putc(*p,fp);
    p++;
  }

  /* New flag */
  fprintf(fp,"%u",newflag);

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);
    p++;			/* Skip over tab */
  }


  /* Print chr, chrpos */
  for (i = 0; i < 2; i++) {
    while (!isspace(*p)) {
      putc(*p,fp);			/* Print field */
      p++;
    }
    if (*p != '\0') {
      putc(*p,fp);
      p++;
    }
  }

  /* MAPQ */
  fprintf(fp,"%d",mapq);

  while (!isspace(*p)) p++;	/* Skip over mapq */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);
    p++;			/* Skip over tab */
  }


  /* Print cigar */
  while (!isspace(*p)) {
    putc(*p,fp);			/* Print field */
    p++;
  }
  if (*p != '\0') {
    putc(*p,fp);
    p++;
  }


  /* New mate chr */
  fprintf(fp,"%s",mate_chr);

  while (!isspace(*p)) p++;	/* Skip over mate_chr */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* New mate chrpos */
  fprintf(fp,"%u",mate_chrpos);

  while (!isspace(*p)) p++;	/* Skip over mate_chrpos */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* New insert length */
  fprintf(fp,"%d",insert_length);

  while (!isspace(*p)) p++;	/* Skip over insert_length */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* print rest of line */
  fprintf(fp,"%s",p);

  return;
}


void
Samread_print_altered_mate (FILE *fp, char *chr, Genomicpos_T chrpos, char *mate_chr, Genomicpos_T mate_chrpos,
			    int insert_length, char *line) {
  char *p, *q;
  int length, i;
  unsigned int flag;
  int mapq;

  p = line;

  /* Print acc */
  while (!isspace(*p)) {
    putc(*p,fp);
    p++;
  }

  if (*p != '\0') {		/* Skip over tab */
    putc(*p,fp);
    p++;
  }

  /* Flag */
  if (sscanf(p,"%u",&flag) != 1) {
    fprintf(stderr,"Unable to find flag in %s\n",p);
    abort();
  } else {
    fprintf(fp,"%u",flag);
  }

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse flag part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);
    p++;			/* Skip over tab */
  }


  /* Print chr, chrpos (could also use chr and chrpos given as parameters) */
  for (i = 0; i < 2; i++) {
    while (!isspace(*p)) {
      putc(*p,fp);			/* Print field */
      p++;
    }
    if (*p != '\0') {
      putc(*p,fp);
      p++;
    }
  }

  /* MAPQ */
  if (sscanf(p,"%d",&mapq) != 1) {
    fprintf(stderr,"Unable to find mapq in %s\n",p);
    abort();
  } else {
    fprintf(fp,"%d",mapq);
  }

  while (!isspace(*p)) p++;	/* Skip over mapq */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);
    p++;			/* Skip over tab */
  }


  /* Print cigar */
  while (!isspace(*p)) {
    putc(*p,fp);			/* Print field */
    p++;
  }
  if (*p != '\0') {
    putc(*p,fp);
    p++;
  }


  /* New mate chr */
  if (mate_chr == NULL) {
    fprintf(fp,"*");
  } else if (!strcmp(mate_chr,chr)) {
    fprintf(fp,"=");
  } else {
    fprintf(fp,"%s",mate_chr);
  }


  while (!isspace(*p)) p++;	/* Skip over mate_chr */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse mate_chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* New mate chrpos */
  if (mate_chr == NULL) {
    fprintf(fp,"0");
  } else {
    fprintf(fp,"%u",mate_chrpos);
  }

  while (!isspace(*p)) p++;	/* Skip over mate_chrpos */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse mate_chrpos part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* New insert length */
  assert(insert_length >= 0);
  if (chrpos < mate_chrpos) {
    fprintf(fp,"%d",insert_length);
  } else {
    fprintf(fp,"%d",-insert_length);
  }

  while (!isspace(*p)) p++;	/* Skip over insert_length */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse insert length part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* print rest of line */
  fprintf(fp,"%s",p);

  return;
}


void
Samread_print_altered (FILE *fp, unsigned int newflag, int mapq, int insert_length, char *line) {
  char *p;
  int i;

  p = line;

  /* Print acc */
  while (!isspace(*p)) {
    putc(*p,fp);
    p++;
  }

  if (*p != '\0') {		/* Skip over tab */
    putc(*p,fp);
    p++;
  }

  fprintf(fp,"%u",newflag);

  while (!isspace(*p)) p++;	/* Skip over flag */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* Skip chr and chrpos */
  for (i = 0; i < 2; i++) {
    while (!isspace(*p)) {
      putc(*p,fp);			/* Print field */
      p++;
    }
    if (*p != '\0') {
      putc(*p,fp);
      p++;
    }
  }

  fprintf(fp,"%d",mapq);

  while (!isspace(*p)) p++;	/* Skip over insert_length */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* Skip cigar, mate_chr, and mate_chrpos */
  for (i = 0; i < 3; i++) {
    while (!isspace(*p)) {
      putc(*p,fp);			/* Print field */
      p++;
    }
    if (*p != '\0') {
      putc(*p,fp);
      p++;
    }
  }


  fprintf(fp,"%d",insert_length);

  while (!isspace(*p)) p++;	/* Skip over insert_length */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    putc(*p,fp);		/* Print tab */
    p++;			/* Skip over tab */
  }

  /* print rest of line */
  fprintf(fp,"%s",p);

  return;
}


char
Samread_splice_strand (char *auxinfo) {
  char *p;
  char tag1, tag2;

  debug(printf("Entering Samread_splice_strand with %s\n",auxinfo));

  p = auxinfo;
  while (*p != '\0' && *p != '\n') {
    tag1 = p[0];
    tag2 = p[1];

    if (tag1 == 'X' && tag2 == 'S') {
      debug(printf("Found tag XS\n"));
      /* XS:A: */
      p += 5;

      if (*p == '+') {
	return '+';
      } else if (*p == '-') {
	return '-';
      } else if (*p == '?') {
	return '?';
      } else {
	fprintf(stderr,"Cannot parse strand %c after XS tag\n",*p);
	return ' ';
      }
    } else {
      while (*p != '\0' && *p != '\t') {
	p++;
      }
      if (*p == '\t') {
	p++;
      }
    }
  }

  return ' ';
}



Intlist_T
Samread_parse_cigar (Uintlist_T *npositions, int *readlength, int *softclip_length, char *cigar) {
  Intlist_T types = NULL;
  unsigned int npos;
  char *p, type;

  *npositions = (Uintlist_T) NULL;
  *readlength = 0;
  *softclip_length = 0;

  if (cigar[0] == '*') {
    return (Intlist_T) NULL;
  }

  p = cigar;
  while (*p != '\0') {
    if (sscanf(p,"%u",&npos) != 1) {
      fprintf(stderr,"Unable to parse cigar %s.  No number in %s\n",cigar,p);
      abort();
    } else {
      *npositions = Uintlist_push(*npositions,npos);
    }

    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
    if (*p == '\0') {
      fprintf(stderr,"Unable to parse cigar %s.  No letter after number %u\n",cigar,npos);
      exit(9);
    } else {
      type = *p++;
      types = Intlist_push(types,(int) type);
    }

    if (type == 'S') {
      if (npos > *softclip_length) {
	*softclip_length = npos;
      }
      *readlength += npos;
    } else if (type == 'M' || type == 'X' || type == 'I') {
      *readlength += npos;
    } else if (type == 'H') {
      *readlength += npos;
    } else if (type == 'D' || type == 'N') {
      /* Ignore */
    } else {
      fprintf(stderr,"Unable to parse cigar %s.  Do not recognize letter %c\n",cigar,type);
      exit(9);
    }
  }

  *npositions = Uintlist_reverse(*npositions);
  return Intlist_reverse(types);
}


void
Samread_print_cigar (Intlist_T types, Uintlist_T npositions) {
  Intlist_T p;
  Uintlist_T q;

  for (p = types, q = npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    printf("%u%c",Uintlist_head(q),Intlist_head(p));
  }
  return;
}


bool
Samread_cigar_hardclipped_p (char *cigar) {
  unsigned int npos;
  char *p, type;

  if (cigar[0] == '*') {
    return false;
  }

  p = cigar;

  while (*p != '\0') {
    if (sscanf(p,"%u",&npos) != 1) {
      fprintf(stderr,"Unable to parse cigar %s.  No number in %s\n",cigar,p);
      abort();
    }

    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
    if (*p == '\0') {
      fprintf(stderr,"Unable to parse cigar %s.  No letter after number %u\n",cigar,npos);
      exit(9);
    } else {
      type = *p++;
    }

    if (type == 'H') {
      return true;
    }
  }

  return false;
}



Genomicpos_T
Samread_chrpos_high (Intlist_T types, Uintlist_T npositions, Genomicpos_T chrpos_low) {
  Intlist_T p;
  Uintlist_T q;
  Genomicpos_T chrpos_high;
  int type;

  chrpos_high = chrpos_low;
  for (p = types, q = npositions; p != NULL; p = Intlist_next(p), q = Uintlist_next(q)) {
    if ((type = Intlist_head(p)) == 'S') {
      /* Ignore */

    } else if (type == 'H') {
      /* Ignore */

    } else if (type == 'M') {
      chrpos_high += Uintlist_head(q);

    } else if (type == 'X') {
      chrpos_high += Uintlist_head(q);

    } else if (type == 'N') {
      chrpos_high += Uintlist_head(q);

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


int
Samread_get_query_coordinates (int *query5, int *query3, Intlist_T types, Uintlist_T npositions,
			       int readlength, char *cigar) {
  int validlength;
  Intlist_T p;
  Uintlist_T q;
  int type;

  *query5 = 1;			/* 1-based */
  *query3 = readlength;
  validlength = 0;

  p = types;
  q = npositions;
  while (p != NULL) {
    if ((type = Intlist_head(p)) == 'S') {
      if (p == types) {
	*query5 = Uintlist_head(q) + 1; /* 1-based */
      } else if (Intlist_next(p) == NULL) {
	*query3 = readlength - Uintlist_head(q);
      } else {
	fprintf(stderr,"Cannot parse cigar %s.  Type S occurs in middle\n",cigar);
	exit(9);
      }
    } else if (type == 'H') {
      /* Do nothing */
    } else if (type == 'M') {
      validlength += Uintlist_head(q);
    } else if (type == 'X') {
      validlength += Uintlist_head(q);
    } else if (type == 'N') {
      /* Do nothing */
    } else if (type == 'I') {
      validlength += Uintlist_head(q);
    } else if (type == 'D') {
      /* Do nothing */
    }
    p = Intlist_next(p);
    q = Uintlist_next(q);
  }

  debug(printf("Got query %d to %d, with length %d\n",*query5,*query3,validlength));
  if (validlength != (*query3) - (*query5) + 1) {
    fprintf(stderr,"Validlength %d from cigar != %d - %d + 1\n",validlength,*query3,*query5);
    abort();
  }

  return validlength;
}



#if 0
Tableuint_T
Samread_get_chromosomes (char *nextchar, Uintlist_T *chrlengths, FILE *input) {
  List_T chromosomes = NULL;
  unsigned int chrlength;
  char *chr, *p, *q;
  int length;

  *chrlengths = (Uintlist_T) NULL;

  while (*nextchar == '@') {
    if (input == NULL || feof(input)) {
      *nextchar = EOF;
      *chrlengths = Uintlist_reverse(*chrlengths);
      return List_reverse(chromosomes);
    }

    if (fgets(&(Header[0]),HEADERLEN,input) == NULL) {
      /* File must terminate after @ */
      *nextchar = EOF;
      *chrlengths = Uintlist_reverse(*chrlengths);
      return List_reverse(chromosomes);
    } else {
      p = &(Header[0]);
      while (*p != '\0' && *p != ':') {
	p++;
      }
      p++;

      q = p;
      while (!isspace((int) *q)) {
	q++;
      }
      length = (q - p)/sizeof(char);
      chr = (char *) MALLOC((length+1)*sizeof(char));
      strncpy(chr,p,length);
      chr[length] = '\0';

      chromosomes = List_push(chromosomes,(void *) chr);

      p = q;
      while (*p != '\0' && *p != ':') {
	p++;
      }
      p++;
      chrlength = strtoul(p,NULL,10);
      *chrlengths = Uintlist_push(*chrlengths,chrlength);
    }
      
    *nextchar = fgetc(input);
  }

  *chrlengths = Uintlist_reverse(*chrlengths);
  return List_reverse(chromosomes);
}

static unsigned int
chromosome_length (List_T chromosomes, Uintlist_T chrlengths, char *desired_chr) {
  List_T p;
  Uintlist_T q;

  p = chromosomes;
  q = chrlengths;
  while (p != NULL) {
    if (!strcmp(desired_chr,(char *) List_head(p))) {
      return Uintlist_head(q);
    }
    p = List_next(p);
    q = Uintlist_next(q);
  }
  fprintf(stderr,"Could not find chromosome %s in SAM input\n",desired_chr);
  exit(9);
}

#endif



int
get_substrings (int *querylength, int **query_starts, Genomicpos_T **genomic_starts, Genomicpos_T **genomic_ends,
		char *cigar, Genomicpos_T chrpos_low) {
  int nsubstrings = 0;
  unsigned int npos;
  char *p, type;

  int querypos = 0;
  Genomicpos_T genomicpos = chrpos_low;
  Intlist_T query_starts_list = NULL;
  Uintlist_T genomic_starts_list = NULL, genomic_ends_list = NULL;

  if (cigar[0] == '*') {
    *querylength = 0;
    *query_starts = (int *) NULL;
    *genomic_starts = (Genomicpos_T *) NULL;
    *genomic_ends = (Genomicpos_T *) NULL;
    return 0;
  }

  query_starts_list = Intlist_push(NULL,querypos);
  genomic_starts_list = Uintlist_push(NULL,genomicpos);

  p = cigar;
  while (*p != '\0') {
    if (sscanf(p,"%u",&npos) != 1) {
      fprintf(stderr,"Unable to parse cigar %s in get_substrings.  No number in %s\n",cigar,p);
      abort();
    }

    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
    if (*p == '\0') {
      fprintf(stderr,"Unable to parse cigar %s.  No letter after number %u\n",cigar,npos);
      exit(9);
    } else {
      type = *p++;
    }

    if (type == 'S') {
      querypos += npos;

    } else if (type == 'M') {
      querypos += npos;
      genomicpos += npos;

    } else if (type == 'X') {
      querypos += npos;
      genomicpos += npos;

    } else if (type == 'I') {
      querypos += npos;

    } else if (type == 'H') {
      /* ? querypos += npos; */

    } else if (type == 'D') {
      genomicpos += npos;

    } else if (type == 'N') {
      genomic_ends_list = Uintlist_push(genomic_ends_list,genomicpos);
      /* nsubstrings++; */

      genomicpos += npos;

      query_starts_list = Intlist_push(query_starts_list,querypos);
      genomic_starts_list = Uintlist_push(genomic_starts_list,genomicpos);

    } else {
      fprintf(stderr,"Unable to parse cigar %s.  Do not recognize letter %c\n",cigar,type);
      exit(9);
    }
  }

  *querylength = querypos;
  genomic_ends_list = Uintlist_push(genomic_ends_list,genomicpos);
  /* nsubstrings++; */


  /* Convert lists to arrays */
  query_starts_list = Intlist_reverse(query_starts_list);
  *query_starts = Intlist_to_array(&nsubstrings,query_starts_list);
  Intlist_free(&query_starts_list);

  genomic_starts_list = Uintlist_reverse(genomic_starts_list);
  *genomic_starts = Uintlist_to_array(&nsubstrings,genomic_starts_list);
  Uintlist_free(&genomic_starts_list);

  genomic_ends_list = Uintlist_reverse(genomic_ends_list);
  *genomic_ends = Uintlist_to_array(&nsubstrings,genomic_ends_list);
  Uintlist_free(&genomic_ends_list);

  return nsubstrings;
}



int
Samread_compute_insert_length (int *querylength5, int *querylength3,
			       char *cigar5, Genomicpos_T chrpos_low_5, char *cigar3, Genomicpos_T chrpos_low_3) {
  int insert_length;
  int nsubstrings5, nsubstrings3, i, j;
  int *query_starts_5, *query_starts_3;
  Genomicpos_T *genomic_starts_5, *genomic_ends_5, *genomic_starts_3, *genomic_ends_3;
  Genomicpos_T pos5, pos3;

  if (cigar5[0] == '*' || cigar3[0] == '*') {
    return 0;
  }

  nsubstrings5 = get_substrings(&(*querylength5),&query_starts_5,&genomic_starts_5,&genomic_ends_5,cigar5,chrpos_low_5);
  nsubstrings3 = get_substrings(&(*querylength3),&query_starts_3,&genomic_starts_3,&genomic_ends_3,cigar3,chrpos_low_3);

  for (i = 0; i < nsubstrings5; i++) {
    for (j = 0; j < nsubstrings3; j++) {
      if (genomic_ends_5[i] < genomic_starts_3[j]) {
	/* No overlap */
      } else if (genomic_starts_5[i] > genomic_ends_3[j]) {
	/* No overlap */
      } else {
	pos5 = genomic_starts_5[i] - query_starts_5[i];
	pos3 = genomic_starts_3[j] - query_starts_3[j];

	FREE(query_starts_5);
	FREE(genomic_starts_5);
	FREE(genomic_ends_5);
	FREE(query_starts_3);
	FREE(genomic_starts_3);
	FREE(genomic_ends_3);

	if (pos5 > pos3) {
	  return (int) (pos5 - pos3);
	} else {
	  return (int) (pos3 - pos5);
	}
      }
    }
  }

  if (genomic_ends_5[nsubstrings5-1] < genomic_starts_3[0]) {
    insert_length = genomic_starts_3[0] - genomic_ends_5[nsubstrings5-1] + (*querylength5) + (*querylength3);
  } else if (genomic_ends_3[nsubstrings3-1] < genomic_starts_5[0]) {
    insert_length = genomic_starts_5[0] - genomic_ends_3[nsubstrings3-1] + (*querylength5) + (*querylength3);
  } else {
    insert_length = 0;
  }

  FREE(query_starts_5);
  FREE(genomic_starts_5);
  FREE(genomic_ends_5);

  FREE(query_starts_3);
  FREE(genomic_starts_3);
  FREE(genomic_ends_3);

  return insert_length;
}


