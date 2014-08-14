static char rcsid[] = "$Id: dynprog_cigar.c 141011 2014-07-09 16:34:33Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "dynprog_cigar.h"
#include "mem.h"
#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Nucleotide comparisons */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


static char *
tokens_string (List_T tokens) {
  char *string, *token, *q;
  int stringlength = 0;
  List_T p;

  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    stringlength += strlen(token);
  }
  string = (char *) MALLOC((stringlength+1) * sizeof(char));

  q = string;
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    strcpy(q,token);
    q += strlen(token);
  }
  *q = '\0';
  return string;
}


static void
tokens_free (List_T *tokens) {
  List_T p;
  char *token;

  for (p = *tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FREE(token);
  }
  List_free(&(*tokens));

  return;
}

static List_T
push_token (bool *startp, List_T tokens, char *token) {
  char *copy;

  /* printf("Pushing token %s\n",token); */
  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);

  *startp = false;
  return List_push(tokens,(void *) copy);
}


static int
revise_md_string (List_T *md_tokens, int matchlength, int r, int c, int diaglength, int diaglength_postins,
		  char *rsequence, char *gsequence, char *gsequence_alt, bool revp) {
  char a, g;
  char token[10];
  bool startp;
  int i;

  debug(printf("Entered revise_md_string with matchlength %d, diaglength %d\n",matchlength,diaglength));
  r += diaglength - 1;
  c += diaglength - 1;

  for (i = diaglength_postins - 1; i >= 0; i--) {
    a = rsequence[r];
    g = gsequence[c];
    debug1(printf("Comparing %c with %c\n",a,g));
    if (a == g) {
      matchlength += 1;
    } else {
      if (matchlength > 0) {
	sprintf(token,"%d",matchlength);
	*md_tokens = push_token(&startp,*md_tokens,token);
	matchlength = 0;
      }
      sprintf(token,"%c",g);	/* Need to take plus strand version */
      *md_tokens = push_token(&startp,*md_tokens,token);
    }
    r--;
    c--;
  }

  debug(printf("Returning matchlength %d\n",matchlength));
  return matchlength;
}

static List_T
add_md_string_deletion_traceback (List_T md_tokens, int c, int dist,
				  char *gsequence, char *gsequence_alt, bool revp,
				  int matchlength) {
  char g;
  char token[10];
  bool startp;
  int i;

  if (matchlength > 0) {
    sprintf(token,"%d",matchlength);
    md_tokens = push_token(&startp,md_tokens,token);
  }

  c += dist - 1;
  for (i = dist - 1; i >= 0; i--) {
    g = gsequence[c];
    sprintf(token,"%c",g);	/* Need to take plus strand version */
    md_tokens = push_token(&startp,md_tokens,token);
    c--;
  }

  md_tokens = push_token(&startp,md_tokens,/*token*/"^");

  return md_tokens;
}


static List_T
add_md_string_deletion_candidate (List_T md_tokens, int matchlength, char *deletion_string) {
  char token[10];
  bool startp;

  if (matchlength > 0) {
    sprintf(token,"%d",matchlength);
    md_tokens = push_token(&startp,md_tokens,token);
  }

  if (deletion_string == NULL) {
    md_tokens = push_token(&startp,md_tokens,/*token*/"^???");
    abort();
  } else {
    md_tokens = push_token(&startp,md_tokens,/*token*/deletion_string);
    md_tokens = push_token(&startp,md_tokens,/*token*/"^");
  }

  return md_tokens;
}

static List_T
finish_md_string (List_T md_tokens, int matchlength) {
  char token[10];
  bool startp;

  if (matchlength > 0) {
    sprintf(token,"%d",matchlength);
    md_tokens = push_token(&startp,md_tokens,token);
  }

  return md_tokens;
}

static List_T
handle_outstanding_insertion (bool *startp, int *matchlength, int *init_softclip, List_T *md_tokens, List_T cigar_tokens,
			      int r, int c, int Mlength_postins, int *Ilength, int Mlength, int insertion_width,
			      char *rsequence, char *gsequence, char *gsequence_alt, bool revp,
			      bool finalp) {
  char token[10];

  debug(printf("startp is %d, finalp is %d, Mlength_postins %d, Ilength %d, Mlength %d\n",
	       *startp,finalp,Mlength_postins,*Ilength,Mlength));

  if (*Ilength == 0) {
    if (Mlength > 0) {
      sprintf(token,"%dM",Mlength);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength,Mlength,
				      rsequence,gsequence,gsequence_alt,revp);
    }

  } else if (*Ilength == insertion_width) {
    debug(printf("Case 0: Ilength is same as insertion_width %d\n",insertion_width));
    if (*startp == true || Mlength_postins > 0) {
      sprintf(token,"%dM",Mlength_postins);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength_postins + (*Ilength) + Mlength,Mlength_postins,
				      rsequence,gsequence,gsequence_alt,revp);
    }
    sprintf(token,"%dI",*Ilength);
    cigar_tokens = push_token(&(*startp),cigar_tokens,token);
    if (finalp == true || Mlength > 0) {
      sprintf(token,"%dM",Mlength);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength,Mlength,
				      rsequence,gsequence,gsequence_alt,revp);
    }

  } else if (*startp == true && finalp == true) {
    /* Not sure what to do in this case */
    debug(printf("Case 1: Ilength %d is different from insertion_width %d\n",*Ilength,insertion_width));
    if (Mlength_postins > Mlength) {
      sprintf(token,"%dM",Mlength_postins);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength_postins + (*Ilength) + Mlength,Mlength_postins,
				      rsequence,gsequence,gsequence_alt,revp);

      sprintf(token,"%dS",Mlength + (*Ilength));
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);

    } else {
      sprintf(token,"%dS",Mlength_postins + (*Ilength));
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *init_softclip = Mlength_postins;

      sprintf(token,"%dM",Mlength);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength,Mlength,
				      rsequence,gsequence,gsequence_alt,revp);
    }

  } else if (*startp == true && finalp == false) {
    debug(printf("Case 2: Ilength %d is different from insertion_width %d\n",*Ilength,insertion_width));
    sprintf(token,"%dS",Mlength_postins + (*Ilength));
    cigar_tokens = push_token(&(*startp),cigar_tokens,token);
    *init_softclip = Mlength_postins; /* Ilength wasn't factored into computation of finalc */

    if (Mlength > 0) {
      sprintf(token,"%dM",Mlength);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength,Mlength,
				      rsequence,gsequence,gsequence_alt,revp);
    }

  } else if (*startp == false && finalp == true) {
    debug(printf("Case 3: Ilength %d is different from insertion_width %d\n",*Ilength,insertion_width));
    if (Mlength_postins > 0) {
      sprintf(token,"%dM",Mlength_postins);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength_postins + (*Ilength) + Mlength,Mlength_postins,
				      rsequence,gsequence,gsequence_alt,revp);
    }

    sprintf(token,"%dS",(*Ilength) + Mlength);
    cigar_tokens = push_token(&(*startp),cigar_tokens,token);

  } else {
    debug(printf("Case 4: Ilength %d is different from insertion_width %d\n",*Ilength,insertion_width));
    /* *startp == false && finalp == false */
    if (Mlength_postins > 0) {
      sprintf(token,"%dM",Mlength_postins);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength_postins + (*Ilength) + Mlength,Mlength_postins,
				      rsequence,gsequence,gsequence_alt,revp);
    }

    sprintf(token,"%dI",*Ilength);
    cigar_tokens = push_token(&(*startp),cigar_tokens,token);

    if (Mlength > 0) {
      sprintf(token,"%dM",Mlength);
      cigar_tokens = push_token(&(*startp),cigar_tokens,token);
      *matchlength = revise_md_string(&(*md_tokens),*matchlength,r,c,Mlength,Mlength,
				      rsequence,gsequence,gsequence_alt,revp);
    }
  }

  *Ilength = 0;

  return cigar_tokens;
}


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
char *
Dynprog_cigar_8 (char **md_string, int *finalc,
		 Direction8_T **directions_nogap, Direction8_T **directions_Egap, Direction8_T **directions_Fgap,
		 int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
		 int *nindels, int queryoffset, int genomeoffset, bool revp,
		 Univcoord_T chroffset, Univcoord_T chrhigh) {
  char *cigar;
  List_T cigar_tokens = NULL, md_tokens = NULL;
  char token[10];
  int Mlength = 0, Ilength = 0, Mlength_postins = 0;
  int matchlength = 0;
  int insertion_width = 0;
  bool startp = true;
  int init_softclip = 0;

  int dist;
  Direction8_T dir;

  debug(printf("Starting Dynprog_cigar_8 at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  *finalc = c;
  while (c > 0 && r > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[c][r]) == DIAG) {
      /* printf("At r %d, c %d (%c), dir is DIAG, nindels %d\n",r,c,gsequence[c],nindels[c-1]); */

      if (nindels[c-1] < 0) {
	/* Genome modified with a deletion */
	cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						    r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						    rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
	Mlength = 1;

	sprintf(token,"%dD",-nindels[c-1]);
	cigar_tokens = push_token(&startp,cigar_tokens,token);

	md_tokens = add_md_string_deletion_candidate(md_tokens,matchlength,delstrings[c-1]);
	matchlength = 0;

      } else if (nindels[c-1] > 0) {
	insertion_width = nindels[c-1];
	if (Ilength == 0) {
	  Mlength_postins = Mlength; /* Cannot push M yet, since I could change to S */
	  Mlength = 0;
	}
	Ilength += 1;
	/* printf("Incrementing Ilength to be %d\n",Ilength); */

      } else if (Mlength == 0) {
	cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						    r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						    rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
	Mlength = 1;

      } else {
	Mlength += 1;
	/* printf("Incrementing Mlength to be %d\n",Mlength); */
      }
      r--; c--;

    } else if (dir == HORIZ) {
      /* printf("At r %d, c %d, dir is HORIZ\n",r,c); */
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
      Mlength = 0;

      dist = 1;
      while (c > 1 && directions_Egap[c][r] != DIAG) {
	dist++;
	c--;
      }
      c--;
      /* dir = directions_nogap[c][r]; */

      sprintf(token,"%dD",dist);
      cigar_tokens = push_token(&startp,cigar_tokens,token);

      md_tokens = add_md_string_deletion_traceback(md_tokens,c,dist,gsequence,gsequence_alt,revp,matchlength);
      matchlength = 0;

    } else {
      /* Must be VERT */
      /* printf("At r %d, c %d, dir is VERT\n",r,c); */
      if (1 || Ilength == 0) {
	if (Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  cigar_tokens = push_token(&startp,cigar_tokens,token);
	  matchlength = revise_md_string(&md_tokens,matchlength,r,c,Mlength,Mlength,
					 rsequence,gsequence,gsequence_alt,revp);
	  Mlength = 0;
	}
      }

      Ilength += 1;
      while (r > 1 && directions_Fgap[c][r] != DIAG) {
	Ilength += 1;
	r--;
      }
      r--;
      /* dir = directions_nogap[c][r]; */

#if 1
      if (Mlength_postins > 0) {
	sprintf(token,"%dM",Mlength_postins);
	cigar_tokens = push_token(&startp,cigar_tokens,token);
	Mlength_postins = 0;
      }
      sprintf(token,"%dI",Ilength);
      cigar_tokens = push_token(&startp,cigar_tokens,token);
      Ilength = 0;
#else
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
#endif

    }
  }

  cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
					      r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
					      rsequence,gsequence,gsequence_alt,revp,/*finalp*/true);

  *finalc = c + init_softclip;
  /* cigar_tokens = List_reverse(cigar_tokens); */
  cigar = tokens_string(cigar_tokens);
  tokens_free(&cigar_tokens);

  md_tokens = finish_md_string(md_tokens,matchlength);
  /* md_tokens = List_reverse(md_tokens); */
  *md_string = tokens_string(md_tokens);
  tokens_free(&md_tokens);

  debug(printf("Returning %s\n\n",cigar));

  return cigar;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* md_string does not appear to have substitutions */
char *
Dynprog_cigar_nogap_8 (int *nmismatches, char **md_string, int *finalc,
		       int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
		       int *nindels, int queryoffset, int genomeoffset, bool revp,
		       Univcoord_T chroffset, Univcoord_T chrhigh) {
  char *cigar;
  List_T cigar_tokens = NULL, md_tokens = NULL;
  char token[10];
  int Mlength = 0, Ilength = 0, Mlength_postins = 0;
  int matchlength = 0;
  int insertion_width = 0;
  bool startp = true;
  int init_softclip = 0;


  debug(printf("Starting Dynprog_cigar_8 at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  *nmismatches = 0;
  *finalc = c;
  while (c > 0 && r > 0) {  /* dir != STOP */
    /* printf("At r %d, c %d (%c), dir is DIAG, nindels %d\n",r,c,gsequence[c],nindels[c-1]); */
    if (rsequence[r] != gsequence[c]) {
      *nmismatches += 1;
    }

    if (nindels[c-1] < 0) {
      /* Genome modified with a deletion */
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
      Mlength = 1;

      sprintf(token,"%dD",-nindels[c-1]);
      cigar_tokens = push_token(&startp,cigar_tokens,token);
      
      md_tokens = add_md_string_deletion_candidate(md_tokens,matchlength,delstrings[c-1]);
      matchlength = 0;

    } else if (nindels[c-1] > 0) {
      insertion_width = nindels[c-1];
      if (Ilength == 0) {
	Mlength_postins = Mlength; /* Cannot push M yet, since I could change to S */
	Mlength = 0;
      }
      Ilength += 1;
      /* printf("Incrementing Ilength to be %d\n",Ilength); */

    } else if (Mlength == 0) {
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
      Mlength = 1;

    } else {
      Mlength += 1;
      /* printf("Incrementing Mlength to be %d\n",Mlength); */
    }
    r--; c--;
  }

  cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
					      r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
					      rsequence,gsequence,gsequence_alt,revp,/*finalp*/true);

  *finalc = c + init_softclip;
  /* cigar_tokens = List_reverse(cigar_tokens); */
  cigar = tokens_string(cigar_tokens);
  tokens_free(&cigar_tokens);

  md_tokens = finish_md_string(md_tokens,matchlength);
  /* md_tokens = List_reverse(md_tokens); */
  *md_string = tokens_string(md_tokens);
  tokens_free(&md_tokens);

  debug(printf("Returning %s\n\n",cigar));

  return cigar;
}
#endif



#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
char *
Dynprog_cigar_16 (char **md_string, int *finalc,
		  Direction16_T **directions_nogap, Direction16_T **directions_Egap, Direction16_T **directions_Fgap,
		  int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
		  int *nindels, int queryoffset, int genomeoffset, bool revp,
		  Univcoord_T chroffset, Univcoord_T chrhigh) {
  char *cigar;
  List_T cigar_tokens = NULL, md_tokens = NULL;
  char token[10];
  int Mlength = 0, Ilength = 0, Mlength_postins = 0;
  int matchlength = 0;
  int insertion_width = 0;
  bool startp = true;
  int init_softclip = 0;

  int dist;
  Direction16_T dir;

  debug(printf("Starting Dynprog_cigar_16 at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  *finalc = c;
  while (c > 0 && r > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[c][r]) == DIAG) {
      /* printf("At r %d, c %d (%c), dir is DIAG, nindels %d\n",r,c,gsequence[c],nindels[c-1]); */

      if (nindels[c-1] < 0) {
	/* Genome modified with a deletion */
	cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						    r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						    rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
	Mlength = 1;

	sprintf(token,"%dD",-nindels[c-1]);
	cigar_tokens = push_token(&startp,cigar_tokens,token);

	md_tokens = add_md_string_deletion_candidate(md_tokens,matchlength,delstrings[c-1]);
	matchlength = 0;

      } else if (nindels[c-1] > 0) {
	insertion_width = nindels[c-1];
	if (Ilength == 0) {
	  Mlength_postins = Mlength; /* Cannot push M yet, since I could change to S */
	  Mlength = 0;
	}
	Ilength += 1;
	/* printf("Incrementing Ilength to be %d\n",Ilength); */

      } else if (Mlength == 0) {
	cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						    r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						    rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
	Mlength = 1;

      } else {
	Mlength += 1;
	/* printf("Incrementing Mlength to be %d\n",Mlength); */
      }
      r--; c--;

    } else if (dir == HORIZ) {
      /* printf("At r %d, c %d, dir is HORIZ\n",r,c); */
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
      Mlength = 0;

      dist = 1;
      while (c > 1 && directions_Egap[c][r] != DIAG) {
	dist++;
	c--;
      }
      c--;
      /* dir = directions_nogap[c][r]; */

      sprintf(token,"%dD",dist);
      cigar_tokens = push_token(&startp,cigar_tokens,token);

      md_tokens = add_md_string_deletion_traceback(md_tokens,c,dist,gsequence,gsequence_alt,revp,matchlength);
      matchlength = 0;

    } else {
      /* Must be VERT */
      /* printf("At r %d, c %d, dir is VERT\n",r,c); */
      if (1 || Ilength == 0) {
	if (Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  cigar_tokens = push_token(&startp,cigar_tokens,token);
	  matchlength = revise_md_string(&md_tokens,matchlength,r,c,Mlength,Mlength,
					 rsequence,gsequence,gsequence_alt,revp);
	  Mlength = 0;
	}
      }

      Ilength += 1;
      while (r > 1 && directions_Fgap[c][r] != DIAG) {
	Ilength += 1;
	r--;
      }
      r--;
      /* dir = directions_nogap[c][r]; */

#if 1
      if (Mlength_postins > 0) {
	sprintf(token,"%dM",Mlength_postins);
	cigar_tokens = push_token(&startp,cigar_tokens,token);
	Mlength_postins = 0;
      }
      sprintf(token,"%dI",Ilength);
      cigar_tokens = push_token(&startp,cigar_tokens,token);
      Ilength = 0;
#else
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
#endif
    }
  }

  cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
					      r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
					      rsequence,gsequence,gsequence_alt,revp,/*finalp*/true);

  *finalc = c + init_softclip;
  /* cigar_tokens = List_reverse(cigar_tokens); */
  cigar = tokens_string(cigar_tokens);
  tokens_free(&cigar_tokens);

  md_tokens = finish_md_string(md_tokens,matchlength);
  /* md_tokens = List_reverse(md_tokens); */
  *md_string = tokens_string(md_tokens);
  tokens_free(&md_tokens);

  debug(printf("Final Mlength = %d\n",Mlength));
  debug(printf("Returning %s\n\n",cigar));

  return cigar;
}
#endif


#if defined(HAVE_SSE4_1) || defined(HAVE_SSE2)
/* md_string does not appear to have substitutions */
char *
Dynprog_cigar_nogap_16 (int *nmismatches, char **md_string, int *finalc,
			int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
			int *nindels, int queryoffset, int genomeoffset, bool revp,
			Univcoord_T chroffset, Univcoord_T chrhigh) {
  char *cigar;
  List_T cigar_tokens = NULL, md_tokens = NULL;
  char token[10];
  int Mlength = 0, Ilength = 0, Mlength_postins = 0;
  int matchlength = 0;
  int insertion_width = 0;
  bool startp = true;
  int init_softclip = 0;


  debug(printf("Starting Dynprog_cigar_16 at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  *nmismatches = 0;
  *finalc = c;
  while (c > 0 && r > 0) {  /* dir != STOP */
    /* printf("At r %d, c %d (%c), dir is DIAG, nindels %d\n",r,c,gsequence[c],nindels[c-1]); */
    if (rsequence[r] != gsequence[c]) {
      *nmismatches += 1;
    }

    if (nindels[c-1] < 0) {
      /* Genome modified with a deletion */
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
      Mlength = 1;

      sprintf(token,"%dD",-nindels[c-1]);
      cigar_tokens = push_token(&startp,cigar_tokens,token);

      md_tokens = add_md_string_deletion_candidate(md_tokens,matchlength,delstrings[c-1]);
      matchlength = 0;

    } else if (nindels[c-1] > 0) {
      insertion_width = nindels[c-1];
      if (Ilength == 0) {
	Mlength_postins = Mlength; /* Cannot push M yet, since I could change to S */
	Mlength = 0;
      }
      Ilength += 1;
      /* printf("Incrementing Ilength to be %d\n",Ilength); */

    } else if (Mlength == 0) {
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
      Mlength = 1;

    } else {
      Mlength += 1;
      /* printf("Incrementing Mlength to be %d\n",Mlength); */
    }
    r--; c--;
  }

  cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
					      r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
					      rsequence,gsequence,gsequence_alt,revp,/*finalp*/true);

  *finalc = c + init_softclip;
  /* cigar_tokens = List_reverse(cigar_tokens); */
  cigar = tokens_string(cigar_tokens);
  tokens_free(&cigar_tokens);

  md_tokens = finish_md_string(md_tokens,matchlength);
  /* md_tokens = List_reverse(md_tokens); */
  *md_string = tokens_string(md_tokens);
  tokens_free(&md_tokens);

  debug(printf("Returning %s\n\n",cigar));

  return cigar;
}
#endif



char *
Dynprog_cigar_std (char **md_string, int *finalc,
		   Direction32_T **directions_nogap, Direction32_T **directions_Egap, Direction32_T **directions_Fgap,
		   int r, int c, char *rsequence, char *gsequence, char *gsequence_alt, char **delstrings,
		   int *nindels, int queryoffset, int genomeoffset, bool revp,
		   Univcoord_T chroffset, Univcoord_T chrhigh) {
  char *cigar;
  List_T cigar_tokens = NULL, md_tokens = NULL;
  char token[10];
  int Mlength = 0, Ilength = 0, Mlength_postins = 0;
  int matchlength = 0;
  int insertion_width = 0;
  bool startp = true;
  int init_softclip = 0;

  int dist;
  Direction32_T dir;

  debug(printf("Starting Dynprog_cigar_std at r=%d,c=%d (roffset=%d, goffset=%d)\n",r,c,queryoffset,genomeoffset));

  *finalc = c;
  while (c > 0 && r > 0) {  /* dir != STOP */
    if ((dir = directions_nogap[c][r]) == DIAG) {
      /* printf("At r %d, c %d (%c), dir is DIAG, nindels %d\n",r,c,gsequence[c],nindels[c-1]); */

      if (nindels[c-1] < 0) {
	/* Genome modified with a deletion */
	cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						    r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						    rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
	Mlength = 1;

	sprintf(token,"%dD",-nindels[c-1]);
	cigar_tokens = push_token(&startp,cigar_tokens,token);

	md_tokens = add_md_string_deletion_candidate(md_tokens,matchlength,delstrings[c-1]);
	matchlength = 0;

      } else if (nindels[c-1] > 0) {
	insertion_width = nindels[c-1];
	if (Ilength == 0) {
	  Mlength_postins = Mlength; /* Cannot push M yet, since I could change to S */
	  Mlength = 0;
	}
	Ilength += 1;
	/* printf("Incrementing Ilength to be %d\n",Ilength); */

      } else if (Mlength == 0) {
	cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						    r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						    rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
	Mlength = 1;

      } else {
	Mlength += 1;
	/* printf("Incrementing Mlength to be %d\n",Mlength); */
      }
      r--; c--;

    } else if (dir == HORIZ) {
      /* printf("At r %d, c %d, dir is HORIZ\n",r,c); */
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
      Mlength = 0;

      dist = 1;
      while (c > 1 && directions_Egap[c][r] != DIAG) {
	dist++;
	c--;
      }
      c--;
      /* dir = directions_nogap[c][r]; */

      sprintf(token,"%dD",dist);
      cigar_tokens = push_token(&startp,cigar_tokens,token);

      md_tokens = add_md_string_deletion_traceback(md_tokens,c,dist,gsequence,gsequence_alt,revp,matchlength);
      matchlength = 0;

    } else {
      /* Must be VERT */
      /* printf("At r %d, c %d, dir is VERT\n",r,c); */
      if (Ilength == 0) {
	if (1 || Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  cigar_tokens = push_token(&startp,cigar_tokens,token);
	  matchlength = revise_md_string(&md_tokens,matchlength,r,c,Mlength,Mlength,
					 rsequence,gsequence,gsequence_alt,revp);
	  Mlength = 0;
	}
      }

      Ilength += 1;
      while (r > 1 && directions_Fgap[c][r] != DIAG) {
	Ilength += 1;
	r--;
      }
      r--;
      /* dir = directions_nogap[c][r]; */

#if 1
      if (Mlength_postins > 0) {
	sprintf(token,"%dM",Mlength_postins);
	cigar_tokens = push_token(&startp,cigar_tokens,token);
	Mlength_postins = 0;
      }
      sprintf(token,"%dI",Ilength);
      cigar_tokens = push_token(&startp,cigar_tokens,token);
      Ilength = 0;
#else
      cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
						  r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
						  rsequence,gsequence,gsequence_alt,revp,/*finalp*/false);
#endif
    }
  }

  cigar_tokens = handle_outstanding_insertion(&startp,&matchlength,&init_softclip,&md_tokens,cigar_tokens,
					      r,c,Mlength_postins,&Ilength,Mlength,insertion_width,
					      rsequence,gsequence,gsequence_alt,revp,/*finalp*/true);

  *finalc = c + init_softclip;
  /* cigar_tokens = List_reverse(cigar_tokens); */
  cigar = tokens_string(cigar_tokens);
  tokens_free(&cigar_tokens);

  md_tokens = finish_md_string(md_tokens,matchlength);
  /* md_tokens = List_reverse(md_tokens); */
  *md_string = tokens_string(md_tokens);
  tokens_free(&md_tokens);

  return cigar;
}
