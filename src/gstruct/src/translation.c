static char rcsid[] = "$Id: translation.c 124200 2014-01-22 22:51:51Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "translation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For toupper */
#include "mem.h"
#include "complement.h"
#include "list.h"


#define IGNORE_MARGIN 6
#define HORIZON 99
#define MIN_NPAIRS 30


/* Finding coding sequence bounds */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


typedef enum {FRAME0, FRAME1, FRAME2, NOFRAME} Frame_T;
static char complCode[128] = COMPLEMENT_UC;
static char uppercaseCode[128] = UPPERCASE_U2T;

#define T Translation_T
struct T {
  int querypos;
  char aa;
  Frame_T frame;
};


static struct T *
Translation_array_new (int translationlen) {
  struct T *new;
  int i;

  new = (struct T *) CALLOC(translationlen,sizeof(struct T));
  for (i = 0; i < translationlen; i++) {
    new[i].querypos = i;
    new[i].aa = ' ';
    new[i].frame = NOFRAME;
  }

  return new;
}

#ifdef DEBUG
static void
Translation_dump (struct T *translation, int translationlen) {
  int i;

  for (i = 0; i < translationlen; i++) {
    printf("%d: ",i);
    switch (translation[i].frame) {
    case NOFRAME: printf("%c %c %c",' ',' ',' '); break;
    case FRAME0: printf("%c %c %c",translation[i].aa,' ',' '); break;
    case FRAME1: printf("%c %c %c",' ',translation[i].aa,' '); break;
    case FRAME2: printf("%c %c %c",' ',' ',translation[i].aa); break;
    }
    printf("\n");
  }
  return;
}
#endif


/************************************************************************/


char
Translation_get_codon (char a, char b, char c) {
  switch (b) {
  case 'T':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'F';
      case 'A': case 'G': case 'R': return 'L';
      default: 	return 'X';
      }
    case 'C': return 'L';
    case 'A': 
      switch (c) {
      case 'G': return 'M';
      case 'T': case 'A': case 'C': case 'H': return 'I';
      default: return 'X';
      }
    case 'G': return 'V';
    default: return 'X';
  }
  case 'C':
    switch (a) {
    case 'T': return 'S';
    case 'C': return 'P';
    case 'A': return 'T';
    case 'G': return 'A';
    default: return 'X';
    }
  case 'A':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'Y';
      case 'A': case 'G': case 'R': return '*';
      default: return 'X';
      }
    case 'C':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'H';
      case 'A': case 'G': case 'R': return 'Q';
      default: return 'X';
      }
    case 'A':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'N';
      case 'A': case 'G': case 'R': return 'K';
      default: return 'X';
      }
    case 'G':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'D';
      case 'A': case 'G': case 'R': return 'E';
      default: return 'X';
      }
    default: return 'X';
    }
  case 'G':
    switch (a) {
    case 'T':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'C';
      case 'A': return '*';
      case 'G': return 'W';
      default: return 'X';
      }
    case 'C': return 'R';
    case 'A':
      switch (c) {
      case 'T': case 'C': case 'Y': return 'S';
      case 'A': case 'G': case 'R': return 'R';
      default: return 'X';
      }
    case 'G': return 'G';
    default: return 'X';
    }
  default: return 'X';
  }
  return 'X';
}


static void
find_bounds_forward (Frame_T *translation_frame, int *translation_starti, 
		     int *translation_endi, int *translation_length,
		     bool *endstopp, struct T *translation, int startpos, int endpos,
		     bool fulllengthp) {
  int beststart0, beststart1, beststart2, bestend0, bestend1, bestend2;
  int bestorf0 = 0, bestorf1 = 0, bestorf2 = 0, orf0 = 0, orf1 = 0, orf2 = 0;
  int start0 = 0, start1 = 0, start2 = 0;
  bool needmet0p, needmet1p, needmet2p;
  bool endstop0p = false, endstop1p = false, endstop2p = false;
  char codon;
  int i, frame;

  if (fulllengthp == true) {
    needmet0p = needmet1p = needmet2p = true;
  } else {
    needmet0p = needmet1p = needmet2p = false;
  }

  for (i = startpos; i < endpos; i++) {
    debug(printf("%d %c: %d %d %d\n",i,translation[i].aa,orf0,orf1,orf2));
    frame = translation[i].frame;
    if ((codon = translation[i].aa) != ' ') {
      if (frame == FRAME0) {
	if (needmet0p) {
	  if (codon == 'M') {
	    orf0 = 1;
	    start0 = i;
	    needmet0p = false;
	  }
	} else if (codon == '*') {
	  orf0++;
	  if (orf0 > bestorf0) {
	    debug(printf("Frame 0: Best orf is %d\n",orf0));
	    bestorf0 = orf0;
	    beststart0 = start0;
	    bestend0 = i;
	    endstop0p = true;
	  }
	  needmet0p = true;
	} else {
	  debug(printf("Incrementing orf0\n"));
	  orf0++;
	}
      } else if (frame == FRAME1) {
	if (needmet1p) {
	  if (codon == 'M') {
	    orf1 = 1;
	    start1 = i;
	    needmet1p = false;
	  }
	} else if (codon == '*') {
	  orf1++;
	  if (orf1 > bestorf1) {
	    debug(printf("Frame 1: Best orf is %d\n",orf1));
	    bestorf1 = orf1;
	    beststart1 = start1;
	    bestend1 = i;
	    endstop1p = true;
	  }
	  needmet1p = true;
	} else {
	  debug(printf("Incrementing orf1\n"));
	  orf1++;
	}
      } else if (frame == FRAME2) {
	if (needmet2p) {
	  if (codon == 'M') {
	    orf2 = 1;
	    start2 = i;
	    needmet2p = false;
	  }
	} else if (codon == '*') {
	  orf2++;
	  if (orf2 > bestorf2) {
	    debug(printf("Frame 2: Best orf is %d\n",orf2));
	    bestorf2 = orf2;
	    beststart2 = start2;
	    bestend2 = i;
	    endstop2p = true;
	  }
	  needmet2p = true;
	} else {
	  debug(printf("Incrementing orf2\n"));
	  orf2++;
	}
      } else {
	fprintf(stderr,"No frame at %d\n",i);
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    beststart0 = start0;
    bestend0 = endpos-1;
    endstop0p = false;
  }
  if (orf1 > bestorf1) {
    debug(printf("Frame 1: Best orf is %d\n",orf1));
    bestorf1 = orf1;
    beststart1 = start1;
    bestend1 = endpos-1;
    endstop1p = false;
  }
  if (orf2 > bestorf2) {
    debug(printf("Frame 2: Best orf is %d\n",orf2));
    bestorf2 = orf2;
    beststart2 = start2;
    bestend2 = endpos-1;
    endstop2p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;
  if (bestorf1 > *translation_length) {
    *translation_length = bestorf1;
    *endstopp = endstop1p;
  }
  if (bestorf2 > *translation_length) {
    *translation_length = bestorf2;
    *endstopp = endstop2p;
  }

  if (bestorf2 == *translation_length) {
    debug(printf("Assigning frame 2\n"));
    *translation_frame = FRAME2;
    *translation_starti = beststart2;
    *translation_endi = bestend2;
  } else if (bestorf1 == *translation_length) {
    debug(printf("Assigning frame 1\n"));
    *translation_frame = FRAME1;
    *translation_starti = beststart1;
    *translation_endi = bestend1;
  } else if (bestorf0 == *translation_length) {
    debug(printf("Assigning frame 0\n"));
    *translation_frame = FRAME0;
    *translation_starti = beststart0;
    *translation_endi = bestend0;
  } else {
    abort();
  }

  debug(printf("Frame0: %d, Frame1: %d, Frame2: %d\n",bestorf0,bestorf1,bestorf2));
  debug(printf("Best orf is %d codons\n",*translation_length));
  debug(printf("Frame0: %d %d\n",beststart0,bestend0));
  debug(printf("Frame1: %d %d\n",beststart1,bestend1));
  debug(printf("Frame2: %d %d\n",beststart2,bestend2));
  debug(printf("Value of endstopp is %d\n",*endstopp));

  return;
}


static void
find_bounds_backward (Frame_T *translation_frame, int *translation_starti,
		      int *translation_endi, int *translation_length,
		      bool *endstopp, struct T *translation, int startpos, int endpos,
		      bool fulllengthp) {
  int beststart0, beststart1, beststart2, bestend0, bestend1, bestend2;
  int bestorf0 = 0, bestorf1 = 0, bestorf2 = 0, orf0 = 0, orf1 = 0, orf2 = 0;
  int start0 = endpos-1, start1 = endpos-1, start2 = endpos-1;
  bool needmet0p, needmet1p, needmet2p;
  bool endstop0p = false, endstop1p = false, endstop2p = false;
  char codon;
  int i, frame;

  if (fulllengthp == true) {
    needmet0p = needmet1p = needmet2p = true;
  } else {
    needmet0p = needmet1p = needmet2p = false;
  }

  for (i = endpos-1; i >= startpos; --i) {
    frame = translation[i].frame;
    if ((codon = translation[i].aa) != ' ') {
      if (frame == FRAME0) {
	if (needmet0p) {
	  if (codon == 'M') {
	    orf0 = 1;
	    start0 = i;
	    needmet0p = false;
	  }
	} else if (codon == '*') {
	  orf0++;
	  if (orf0 > bestorf0) {
	    debug(printf("Frame 0: Best orf is %d\n",orf0));
	    bestorf0 = orf0;
	    bestend0 = i;
	    beststart0 = start0;
	    endstop0p = true;
	  }
	  needmet0p = true;
	} else {
	  orf0++;
	}
      } else if (frame == FRAME1) {
	if (needmet1p) {
	  if (codon == 'M') {
	    orf1 = 1;
	    start1 = i;
	    needmet1p = false;
	  }
	} else if (codon == '*') {
	  orf1++;
	  if (orf1 > bestorf1) {
	    debug(printf("Frame 1: Best orf is %d\n",orf1));
	    bestorf1 = orf1;
	    bestend1 = i;
	    beststart1 = start1;
	    endstop1p = true;
	  }
	  needmet1p = true;
	} else {
	  orf1++;
	}
      } else if (frame == FRAME2) {
	if (needmet2p) {
	  if (codon == 'M') {
	    orf2 = 1;
	    start2 = i;
	    needmet2p = false;
	  }
	} else if (codon == '*') {
	  orf2++;
	  if (orf2 > bestorf2) {
	    debug(printf("Frame 2: Best orf is %d\n",orf2));
	    bestorf2 = orf2;
	    bestend2 = i;
	    beststart2 = start2;
	    endstop2p = true;
	  }
	  needmet2p = true;
	} else {
	  orf2++;
	}
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    bestend0 = 0;
    beststart0 = start0;
    endstop0p = false;
  }
  if (orf1 > bestorf1) {
    debug(printf("Frame 1: Best orf is %d\n",orf1));
    bestorf1 = orf1;
    bestend1 = 0;
    beststart1 = start1;
    endstop1p = false;
  }
  if (orf2 > bestorf2) {
    debug(printf("Frame 2: Best orf is %d\n",orf2));
    bestorf2 = orf2;
    bestend2 = 0;
    beststart2 = start2;
    endstop2p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;
  if (bestorf1 > *translation_length) {
    *translation_length = bestorf1;
    *endstopp = endstop1p;
  }
  if (bestorf2 > *translation_length) {
    *translation_length = bestorf2;
    *endstopp = endstop2p;
  }

  if (bestorf0 == *translation_length) {
    debug(printf("Assigning frame 0\n"));
    *translation_frame = FRAME0;
    *translation_starti = beststart0;
    *translation_endi = bestend0;
  } else if (bestorf1 == *translation_length) {
    debug(printf("Assigning frame 1\n"));
    *translation_frame = FRAME1;
    *translation_starti = beststart1;
    *translation_endi = bestend1;
  } else if (bestorf2 == *translation_length) {
    debug(printf("Assigning frame 2\n"));
    *translation_frame = FRAME2;
    *translation_starti = beststart2;
    *translation_endi = bestend2;
  } else {
    abort();
  }

  debug(printf("Frame0: %d, Frame1: %d, Frame2: %d\n",bestorf0,bestorf1,bestorf2));
  debug(printf("Best orf is %d codons\n",*translation_length));
  debug(printf("Frame0: %d %d\n",beststart0,bestend0));
  debug(printf("Frame1: %d %d\n",beststart1,bestend1));
  debug(printf("Frame2: %d %d\n",beststart2,bestend2));
  debug(printf("Value of endstopp is %d\n",*endstopp));

  return;
}


static void
find_bounds_forward_fromstart (Frame_T *translation_frame, int *translation_starti, 
			       int *translation_endi, int *translation_length,
			       bool *endstopp, struct T *translation, int startpos, int endpos,
			       int cds_startpos) {
  Frame_T frame_fromstart;
  int phase_fromstart;
  int beststart0, bestend0;
  int bestorf0 = 0, orf0 = 0;
  int start0 = 0;
  bool endstop0p = false;
  char codon;
  int i;

  phase_fromstart = (cds_startpos - 1) % 3;
  if (translation[0].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[0].frame;
  } else if (translation[1].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[1].frame;
  } else if (translation[2].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[2].frame;
  } else {
    fprintf(stderr,"Something wrong with Translation_T object\n");
    abort();
  }

  for (i = startpos; i < endpos && endstop0p == false; i++) {
    if (translation[i].querypos >= cds_startpos - 1 &&
	translation[i].frame == frame_fromstart &&
	(codon = translation[i].aa) != ' ') {
      debug(printf("%d %c\n",i,translation[i].aa));
      if (orf0 == 0) {
	start0 = i;
      }
      if (codon == '*') {
	orf0++;
	if (orf0 > bestorf0) {
	  debug(printf("Frame 0: Best orf is %d\n",orf0));
	  bestorf0 = orf0;
	  beststart0 = start0;
	  bestend0 = i;
	  endstop0p = true;
	}

      } else {
	debug(printf("Incrementing orf0\n"));
	orf0++;
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    beststart0 = start0;
    bestend0 = endpos-1;
    endstop0p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;

  debug(printf("Assigning frame %d\n",frame_fromstart));
  *translation_frame = frame_fromstart;
  *translation_starti = beststart0;
  *translation_endi = bestend0;

  return;
}


static void
find_bounds_backward_fromstart (Frame_T *translation_frame, int *translation_starti,
				int *translation_endi, int *translation_length,
				bool *endstopp, struct T *translation, int startpos, int endpos,
				int cds_startpos) {
  Frame_T frame_fromstart;
  int phase_fromstart;
  int beststart0, bestend0;
  int bestorf0 = 0, orf0 = 0;
  int start0 = endpos-1;
  bool endstop0p = false;
  char codon;
  int i;

  phase_fromstart = (endpos - cds_startpos) % 3;
  if (translation[endpos-1].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[endpos-1].frame;
  } else if (translation[endpos-2].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[endpos-2].frame;
  } else if (translation[endpos-3].querypos % 3 == phase_fromstart) {
    frame_fromstart = translation[endpos-3].frame;
  } else {
    fprintf(stderr,"Something wrong with Translation_T object\n");
    abort();
  }
    
  for (i = endpos-1; i >= startpos && endstop0p == false; --i) {
    if (translation[i].querypos >= cds_startpos - 1 &&
	translation[i].frame == frame_fromstart &&
	(codon = translation[i].aa) != ' ') {
      debug(printf("%d %c\n",i,translation[i].aa));
      if (orf0 == 0) {
	start0 = i;
      }
      if (codon == '*') {
	orf0++;
	if (orf0 > bestorf0) {
	  debug(printf("Frame 0: Best orf is %d\n",orf0));
	  bestorf0 = orf0;
	  bestend0 = i;
	  beststart0 = start0;
	  endstop0p = true;
	}
      } else {
	debug(printf("Incrementing orf0\n"));
	orf0++;
      }
    }
  }
  
  /* Handle last segments */
  if (orf0 > bestorf0) {
    debug(printf("Frame 0: Best orf is %d\n",orf0));
    bestorf0 = orf0;
    bestend0 = 0;
    beststart0 = start0;
    endstop0p = false;
  }

  /* Find overall best */
  *translation_length = bestorf0;
  *endstopp = endstop0p;

  debug(printf("Assigning frame %d\n",frame_fromstart));
  *translation_frame = frame_fromstart;
  *translation_starti = beststart0;
  *translation_endi = bestend0;

  return;
}



static struct T *
translate_pairs_forward (char *genome, int npairs, bool revcompp) {
  struct T *translation;
  int i, gpos = 0;
  char codon, nt2 = 'X', nt1 = 'X', nt0 = 'X';

  translation = Translation_array_new(npairs);

  /* Go backward so we put aa at beginning of codon */
  for (i = npairs-1; i >= 0; --i) {
    nt2 = nt1;
    nt1 = nt0;
    nt0 = revcompp ? complCode[(int) genome[i]] : uppercaseCode[(int) genome[i]];

    codon = Translation_get_codon(nt0,nt1,nt2);
    if (gpos < 2 && codon == 'X') {
      /* translation[i].aa = ' '; */
    } else {
      translation[i].aa = codon;
      switch (gpos % 3) {
      case 0: translation[i].frame = FRAME0; break;
      case 1: translation[i].frame = FRAME1; break;
      case 2: translation[i].frame = FRAME2; break;
      }
    }
    gpos++;
  }

  return translation;
}  

static struct T *
translate_pairs_backward (char *genome, int npairs, bool revcompp) {
  struct T *translation;
  int i, gpos = 0;
  char codon, nt2 = 'X', nt1 = 'X', nt0 = 'X';

  translation = Translation_array_new(npairs);

  /* Go forwards so we put aa at beginning of codon */
  for (i = 0; i < npairs; i++) {
    nt2 = nt1;
    nt1 = nt0;
    nt0 = revcompp ? complCode[(int) genome[i]] : uppercaseCode[(int) genome[i]];

    codon = Translation_get_codon(nt0,nt1,nt2);
    if (gpos < 2 && codon == 'X') {
      /* translation[i].aa = ' '; */
    } else {
      translation[i].aa = codon;
      switch (gpos % 3) {
      case 0: translation[i].frame = FRAME0; break;
      case 1: translation[i].frame = FRAME1; break;
      case 2: translation[i].frame = FRAME2; break;
      }
    }
    gpos++;
  }

  return translation;
}  


char *
Translation_via_genomic (int *translation_leftpos, int *translation_rightpos, int *translation_length,
			 char *genome, int startpos, int endpos, int npairs,
			 bool backwardp, bool revcompp, bool fulllengthp, int cds_startpos) {
  char *aa_g;
  char lastaa;
  struct T *translation;
  bool endstopp;
  int i, aapos = 0;
  Frame_T translation_frame;
  int translation_starti = 0, translation_endi = 0, phase;
  int minpos, maxpos;

  aa_g = (char *) CALLOC(npairs,sizeof(char));

  if (backwardp == false) {
    translation = translate_pairs_forward(genome,npairs,revcompp);
    if (cds_startpos > 0) {
      find_bounds_forward_fromstart(&translation_frame,&translation_starti,
				    &translation_endi,&(*translation_length),&endstopp,
				    translation,startpos,endpos,cds_startpos);
    } else {
      find_bounds_forward(&translation_frame,&translation_starti,
			  &translation_endi,&(*translation_length),&endstopp,
			  translation,startpos,endpos,fulllengthp);
      if (fulllengthp == true && *translation_length == 0) {
	/* fprintf(stderr,"No full length gene found.  Assuming partial length.\n"); */
	find_bounds_forward(&translation_frame,&translation_starti,
			    &translation_endi,&(*translation_length),&endstopp,
			    translation,startpos,endpos,/*fulllengthp*/false);
      }
    }

  } else {
    translation = translate_pairs_backward(genome,npairs,revcompp);
    if (cds_startpos > 0) {
      find_bounds_backward_fromstart(&translation_frame,&translation_starti,
				     &translation_endi,&(*translation_length),&endstopp,
				     translation,startpos,endpos,cds_startpos);
    } else {
      find_bounds_backward(&translation_frame,&translation_starti,
			   &translation_endi,&(*translation_length),&endstopp,
			   translation,startpos,endpos,fulllengthp);
      if (fulllengthp == true && *translation_length == 0) {
	/* fprintf(stderr,"No full length gene found.  Assuming partial length.\n"); */
	find_bounds_backward(&translation_frame,&translation_starti,
			     &translation_endi,&(*translation_length),&endstopp,
			     translation,startpos,endpos,/*fulllengthp*/false);
      }
    }
  }
  /* debug(printf("ref:\n")); */
  debug(Translation_dump(translation,/*translationlen*/npairs));

  for (i = 0; i < npairs; i++) {
    aa_g[i] = ' ';
  }

  if (translation_starti < 0 || translation_endi < 0) {
    *translation_leftpos = *translation_rightpos = -1;
  } else {
    minpos = endpos - 1;
    maxpos = startpos;
    if (backwardp == false) {
      /* forward */
      debug(printf("Translation is forward pairs %d..%d\n",translation_starti,translation_endi));
      for (i = translation_starti; i <= translation_endi; i++) {
	/* Avoid problems with genome position advancing prematurely */
	if (genome[i] != ' ') {
	  if (translation[i].frame == translation_frame) {
	    if ((aa_g[i] = translation[i].aa) != ' ') {
	      if (i < minpos) {
		minpos = i;
	      }
	      if (i > maxpos) {
		maxpos = i;
	      }
	      lastaa = aa_g[i];
	      aapos++;
	      /* aaphase_g[i] = 0; */
	    }
	  } else if (translation[i].frame != 3) {
	    if ((phase = translation_frame - translation[i].frame) < 0) {
	      phase += 3;
	    }
	    /* aaphase_g[i] = phase; */
	    switch (phase) {
	    case 0: abort();
	    case 1: aa_g[i] = '1'; break;
	    case 2: aa_g[i] = '2'; break;
	    default: abort();
	    }
	  }
	}
      }
      *translation_leftpos = minpos;
      if ((*translation_rightpos = maxpos + 2) >= endpos) {
	*translation_rightpos = endpos - 1;
      }
      if (lastaa == '*') {
	*translation_length -= 1;
      }
      
      for (i = 0; i < *translation_leftpos; i++) {
	aa_g[i] = '5';
      }
      for (i = *translation_rightpos + 1; i < npairs; i++) {
	aa_g[i] = '3';
      }

    } else {
      /* backward */
      debug(printf("Translation is backward pairs %d..%d\n",translation_starti,translation_endi));
      
      for (i = translation_starti; i >= translation_endi; --i) {
	/* Avoid problems with genome position advancing prematurely */
	if (genome[i] != ' ') {
	  if (translation[i].frame == translation_frame) {
	    if ((aa_g[i] = translation[i].aa) != ' ') {
	      if (i < minpos) {
		minpos = i;
	      }
	      if (i > maxpos) {
		maxpos = i;
	      }
	      lastaa = aa_g[i];
	      aapos++;
	      /* aaphase_g[i] = 0; */
	    }
	  } else if (translation[i].frame != 3) {
	    if ((phase = translation_frame - translation[i].frame) < 0) {
	      phase += 3;
	    }
	    /* aaphase_g[i] = phase; */
	    switch (phase) {
	    case 0: abort();
	    case 1: aa_g[i] = '1'; break;
	    case 2: aa_g[i] = '2'; break;
	    default: abort();
	    }
	  }
	}
      }
      if ((*translation_leftpos = minpos - 2) < startpos) {
	*translation_leftpos = startpos;
      }
      *translation_rightpos = maxpos;
      if (lastaa == '*') {
	*translation_length -= 1;
      }
      
      for (i = 0; i < *translation_leftpos; i++) {
	aa_g[i] = '5';
      }
      for (i = *translation_rightpos + 1; i < npairs; i++) {
	aa_g[i] = '3';
      }
    }
  }

  FREE(translation);

  return aa_g;
}

