 static char rcsid[] = "$Id: bam_print.c 138430 2014-06-06 21:26:50Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>		/* For qsort, rint */

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "chrom.h"
#include "genomicpos.h"
#include "intlist.h"
#include "uintlist.h"
#include "list.h"
#include "iit-read.h"
#include "interval.h"
#include "table.h"
#include "uinttable.h"
#include "maxent_hr.h"		/* For Maxent_hr_setup */

#ifdef BAM_INPUT
#include "bamread.h"
#endif

#include "samflags.h"
#include "samread.h"
#include "genome.h"
#include "complement.h"

#include "bamtally.h"
#include "gstruct.h"

#include "parserange.h"
#include "datadir.h"
#include "getopt.h"


#define E 2.71828182845

#define SPACER_YHEIGHT 20
#define AXIS_YHEIGHT 10
#define TALLY_YHEIGHT 50
#define EDGES_YHEIGHT 50
#define EXTENTS_YHEIGHT 50

#define SPLICES_YDELTA 15
#define RAWREADS_YDELTA 4

#define BLACK /*red*/0,/*green*/0,/*blue*/0
#define RED /*red*/1,/*green*/0,/*blue*/0
#define GRAY /*red*/0.8,/*green*/0.8,/*blue*/0.8
#define GREEN /*red*/0,/*green*/1,/*blue*/0
#define VERY_LIGHT_BLUE /*red*/0.87,/*green*/0.92,/*blue*/0.97
#define LIGHT_BLUE /*red*/0.62,/*green*/0.79,/*blue*/0.88
#define DARK_BLUE /*red*/0.42,/*green*/0.68,/*blue*/0.84
#define BEIGE /*red*/0.93,/*green*/0.91,/*blue*/0.67

/*
static int
COLORBLIND_RED[NGROUPS] = {77, 55, 228, 152, 255, 255, 166, 247, 153};
static int
COLORBLIND_GREEN[NGROUPS] = {175, 126, 26, 78, 127, 255, 86, 129, 153};
static int
COLORBLIND_BLUE[NGROUPS] = {74, 184, 28, 163, 0, 51, 40, 191, 153};
*/

#define BREWER_RED /*red*/0.894,/*green*/0.102,/*blue*/0.110
#define BREWER_GREEN /*red*/0.302,/*green*/0.686,/*blue*/0.290



#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Sequence differences */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Needed for Genome_T */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_mapdir = NULL;
static List_T map_iitfiles = NULL;

static char *arcfile = NULL;


/* Filters */
static int min_depth = 1;

/* For splicegene */
static int max_exonlength = 2000;

/* Output options */
static int xbase = 36;
static int pagewidth = 540;	/* 7.5 * 72 */
static int fontsize = 8;

static int ybase = 36;
static int pageheight = 700;	/* 10.5 * 72 - 20 for args at bottom */
static int argsheight = 12;

static List_T tracks = NULL;
static bool separate_tallies_p = false;
static bool separate_extents_p = false;
static bool print_edges_p = false;
static bool need_gstruct_p = false;
static bool print_seqdiffs_p = true;


static char *chromosome = NULL;
static int blocksize = 1000;

static bool need_canonical_p = false;

static char *desired_read_group = NULL;
static int minimum_mapq = 0;
static int good_unique_mapq = 35;
static int maximum_nhits = 1000000;
static bool need_unique_p = false;
static bool need_primary_p = false;
static bool ignore_duplicates_p = false;
static bool need_concordant_p = false;

static Genomicpos_T max_pairlength = 1000000;
static Genomicpos_T min_pairlength = 0;
static Genomicpos_T shortsplicedist = 400000;

static int mincount = 1;
static int mincount_end_alt = 2;
static int minsupport = 8;
static bool monochromep = false;


/* static int min_mlength = 0; */

/* For Illumina, subtract 64.  For Sanger, subtract 33 */
/* For BAM, quality is already adjusted down to 0 */
static int quality_score_adj = 0;

#if 0
static int quality_score_constant = -1;
#endif

static int alloclength = 400000;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Known genes */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfiles */

  /* Filters */
  {"depth", required_argument, 0, 0},	   /* min_depth */

  /* Output options */
  {"tracks", required_argument, 0, 't'}, /* tracks */
  {"fontsize", required_argument, 0, 'f'}, /* fontsize */
  {"sep-tallies", no_argument, 0, 0}, /* separate_tallies_p */
  {"sep-extents", no_argument, 0, 0}, /* separate_extents_p */
  {"no-seqdiffs", no_argument, 0, 0}, /* print_seqdiffs_p */

  {"read-group", required_argument, 0, 0},   /* desired_read_group */
  {"mapq", required_argument, 0, 'q'}, /* minimum_mapq */
  {"nhits", required_argument, 0, 'n'}, /* maximum_nhits */
  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* need_unique_p */
  {"primary", required_argument, 0, 'P'}, /* need_primary_p */
  {"ignore-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */
  {"allow-duplicates", no_argument, 0, 0}, /* ignore_duplicates_p */

  {"pairmax", required_argument, 0, 'p'}, /* alloclength */
  {"mincount", required_argument, 0, 's'}, /* mincount */
  {"monochrome", no_argument, 0, 0}, /* monochromep */


  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"BAM_PRINT\n");
  fprintf(stdout,"Part of GSTRUCT package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}



static void
print_program_usage ();


/* Coordinate */
static double
calc_xpos (unsigned int coord, unsigned int maxcoord, unsigned int mincoord,
	   double left, double right) {
  return (right-left)/((double) maxcoord-(double) mincoord)*((double) coord-(double) mincoord) + left;
}

/* Value */
static double
calc_ypos (int lineindex, double linesep, double top, double bottom) {
  return top - linesep*(lineindex+1);
}


static double prev_red = 0.0;
static double prev_green = 0.0;
static double prev_blue = 0.0;


static void
print_horizontal_line (double x1, double y1, double x2, double y2,
		       double left, double right,
		       double red, double green, double blue) {

  if (red != prev_red ||
      green != prev_green ||
      blue != prev_blue) {
    printf("%f %f %f setrgbcolor\n",red,green,blue);
    prev_red = red;
    prev_green = green;
    prev_blue = blue;
  }

  printf("newpath\n");
  if (x1 < left) {
    x1 = left;
  } else if (x1 > right) {
    x1 = right;
  }
  printf("%.2f %.2f moveto\n",x1,y1);

  if (x2 < left) {
    x2 = left;
  } else if (x2 > right) {
    x2 = right;
  }
  printf("%.2f %.2f lineto\n",x2,y2);
  printf("stroke\n");
  return;
}

static void
print_vertical_line (double x1, double y1, double x2, double y2,
		     double red, double green, double blue) {

  if (red != prev_red ||
      green != prev_green ||
      blue != prev_blue) {
    printf("%f %f %f setrgbcolor\n",red,green,blue);
    prev_red = red;
    prev_green = green;
    prev_blue = blue;
  }

  printf("newpath\n");
  printf("%.2f %.2f moveto\n",x1,y1);
  printf("%.2f %.2f lineto\n",x2,y2);
  printf("stroke\n");
  return;
}


static void
print_rect (double xpos1, double ypos1, double xpos2, double ypos2,
	    double red, double green, double blue) {
  if (red != prev_red ||
      green != prev_green ||
      blue != prev_blue) {
    printf("%f %f %f setrgbcolor\n",red,green,blue);
    prev_red = red;
    prev_green = green;
    prev_blue = blue;
  }

  printf("newpath\n");
  printf("%.2f %.2f moveto\n",xpos1,ypos1);
  printf("%.2f %.2f lineto\n",xpos1,ypos2);
  printf("%.2f %.2f lineto\n",xpos2,ypos2);
  printf("%.2f %.2f lineto\n",xpos2,ypos1);
  printf("closepath\n");
  printf("fill\n");
  return;
}


static void
print_arc_procedure () {
  printf("\
/ellipsedict 8 dict def\n\
ellipsedict /mtrx matrix put\n\
/ellipse\n\
  { ellipsedict begin\n\
    /endangle exch def\n\
    /startangle exch def\n\
    /yrad exch def\n\
    /xrad exch def\n\
    /y exch def\n\
    /x exch def\n\
    /savematrix mtrx currentmatrix def\n\
    x y translate\n\
    xrad yrad scale\n\
    0 0 1 startangle endangle arc\n\
    savematrix setmatrix\n\
    end\n\
  } def\n\
");
  return;
}


static void
print_arc (double xpos1, double xpos2, double ypos, double yheight,
	   double red, double green, double blue) {
  double x, xrad;

  if (red != prev_red ||
      green != prev_green ||
      blue != prev_blue) {
    printf("%f %f %f setrgbcolor\n",red,green,blue);
    prev_red = red;
    prev_green = green;
    prev_blue = blue;
  }

  x = (xpos1 + xpos2)/2.0;
  if (xpos2 >= xpos1) {
    xrad = (xpos2 - xpos1)/2.0;
  } else {
    xrad = (xpos1 - xpos2)/2.0;
  }

  printf("newpath\n");
  printf("%.2f %.2f %.2f %.2f 0 180 ellipse\n",x,ypos,xrad,yheight);
  printf("stroke\n");

  return;
}


static unsigned int
power (int base, int exponent) {
  unsigned int result = 1;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


static void
axis_print (Genomicpos_T mincoord, Genomicpos_T maxcoord,
	    double xorigin, double left, double right, double yorigin) {
  double proposed_interval;
  Genomicpos_T interval, big_interval, chrpos;
  double xpos;
  char *string;

  proposed_interval = (double) (maxcoord - mincoord + 1)/(double) 30;
  if (proposed_interval < 1.0) {
    interval = 1U;
    big_interval = 10;
  } else {
    interval = power(10,rint(log10(proposed_interval)));
    big_interval = 10*interval;
  }

  print_horizontal_line(xorigin,yorigin,xorigin+pagewidth,yorigin,left,right,BLACK);

  string = Genomicpos_commafmt(mincoord);
  printf("gsave\n");
  /* printf("0 setgray\n"); */
  printf("%.2f %.2f moveto\n",xorigin,yorigin);
  printf("0 5 rmoveto\n");
  printf("(%s) dup stringwidth pop neg 0 rmoveto show\n",string);
  printf("grestore\n");
  FREE(string);

  string = Genomicpos_commafmt(maxcoord);
  printf("gsave\n");
  /* printf("0 setgray\n"); */
  printf("%.2f %.2f moveto\n",xorigin+pagewidth,yorigin);
  printf("0 5 rmoveto\n");
  printf("(%s) show\n",string);
  printf("grestore\n");
  FREE(string);


  for (chrpos = interval * (int) (mincoord/interval); chrpos < maxcoord;
       chrpos += interval) {
    xpos = calc_xpos(chrpos,maxcoord,mincoord,left,right);

    if (chrpos % big_interval == 0) {

      print_vertical_line(xpos,yorigin,xpos,yorigin-8,BLACK);
      string = Genomicpos_commafmt(chrpos);
      printf("gsave\n");
      /* printf("0 setgray\n"); */
      printf("%.2f %.2f moveto\n",xpos,yorigin);
      printf("0 5 rmoveto\n");
      printf("(%s) dup stringwidth pop 2 div neg 0 rmoveto show\n",string);
      printf("grestore\n");
      FREE(string);

    } else {
      print_vertical_line(xpos,yorigin,xpos,yorigin-4,BLACK);
    }
  }

  return;
}


static void
scale_print (double xorigin, double yorigin, double yheight, double ydensity) {
  print_vertical_line(xorigin - 10,yorigin,xorigin - 10,yorigin + yheight,BLACK);
  printf("gsave\n");
  /* printf("0 setgray\n"); */
  printf("%.2f %.2f moveto\n",xorigin - 10,yorigin + yheight);
  printf("-5 -2 rmoveto\n");
  printf("(%d) dup stringwidth pop neg 0 rmoveto show\n",(int) rint(ydensity));
  printf("grestore\n");
  return;
}


static void
tally_print (long int *tally_matches, long int *tally_mismatches,
	     Genomicpos_T mincoord, Genomicpos_T maxcoord, double binstep,
	     double xorigin, double yorigin, double yheight) {
  Genomicpos_T chrpos;
  int nbins, bini, nextbini, bini_max;
  double bin_low, bin_high, marker, width, pct, maxpct;
  double *binx, total;
  double ymax = 0.0;
  double yfactor, height1, height2;
  double *maf;


  nbins = (int) ((maxcoord - mincoord + 1)/binstep + 1);
  binx = (double *) CALLOC(nbins,sizeof(double));
  maf = (double *) CALLOC(nbins,sizeof(double));

  for (chrpos = mincoord; chrpos <= maxcoord; chrpos++) {
    total = (double) (tally_matches[chrpos - mincoord] + tally_mismatches[chrpos - mincoord]);

    bin_low = (double) (chrpos - mincoord)/binstep;
    bin_high = (double) (chrpos + 1U - mincoord)/binstep;
    width = bin_high - bin_low;

    marker = bin_low;
    bini = (int) marker;
    nextbini = bini + 1;
    maxpct = 0.0;
    while (nextbini < bin_high) {
      pct = ((double) nextbini - marker)/width;
      if (pct > maxpct) {
	maxpct = pct;
	bini_max = bini;
      }
      binx[bini] += pct * total;

      marker = nextbini;
      bini++;
      nextbini++;
    }

    pct = (bin_high - marker)/width;
    if (pct > maxpct) {
      maxpct = pct;
      bini_max = bini;
    }
    binx[bini] += pct * total;

    if ((double) tally_mismatches[chrpos - mincoord]/total > maf[bini_max]) {
      maf[bini_max] = (double) tally_mismatches[chrpos - mincoord]/(double) total;
    }
  }


  for (bini = 0; bini < nbins; bini++) {
    if (binx[bini] > ymax) {
      ymax = binx[bini];
    }
  }

  if (ymax > 0) {
    scale_print(xorigin,yorigin,yheight,/*ydensity*/ymax/binstep);
    yfactor = yheight / (double) ymax;

    for (bini = 0; bini < nbins; bini++) {
      if (binx[bini] > 0) {
	height2 = binx[bini] * yfactor;
	height1 = maf[bini] * height2;
	print_vertical_line(xorigin + bini,yorigin,xorigin + bini,yorigin + height1,
			    DARK_BLUE);
	print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
			    GRAY);
	  
      }
    }
  }

  FREE(maf);
  FREE(binx);

  return;
}


static void
tally_print_lh (long int *tally_matches_low, long int *tally_mismatches_low,
		long int *tally_matches_high, long int *tally_mismatches_high,
		Genomicpos_T mincoord, Genomicpos_T maxcoord, double binstep,
		double xorigin, double yorigin, double yheight) {
  Genomicpos_T chrpos;
  int nbins, bini, nextbini;
  double bin_low, bin_high, marker, width, pct;
  double *binx_lowend, *binx_highend, total_lowend, total_highend;
  double ymax = 0.0, y;
  double yfactor, height1, height2;


  nbins = (int) ((maxcoord - mincoord + 1)/binstep + 1);
  binx_lowend = (double *) CALLOC(nbins,sizeof(double));
  binx_highend = (double *) CALLOC(nbins,sizeof(double));

  for (chrpos = mincoord; chrpos <= maxcoord; chrpos++) {
    total_lowend = (double) (tally_matches_low[chrpos - mincoord] + tally_mismatches_low[chrpos - mincoord]);
    total_highend = (double) (tally_matches_high[chrpos - mincoord] + tally_mismatches_high[chrpos - mincoord]);

    bin_low = (double) (chrpos - mincoord)/binstep;
    bin_high = (double) (chrpos + 1U - mincoord)/binstep;
    width = bin_high - bin_low;

    marker = bin_low;
    bini = (int) marker;
    nextbini = bini + 1;
    while (nextbini < bin_high) {
      pct = ((double) nextbini - marker)/width;
      binx_lowend[bini] += pct * total_lowend;
      binx_highend[bini] += pct * total_highend;

      marker = nextbini;
      bini++;
      nextbini++;
    }

    pct = (bin_high - marker)/width;
    binx_lowend[bini] += pct * total_lowend;
    binx_highend[bini] += pct * total_highend;
  }
  fprintf(stderr,"\n");

  for (bini = 0; bini < nbins; bini++) {
    if ((y = binx_lowend[bini] + binx_highend[bini]) > ymax) {
      ymax = y;
    }
  }

  if (ymax > 0) {
    print_vertical_line(xorigin - 10,yorigin,xorigin - 10,yorigin + yheight,BLACK);
    printf("gsave\n");
    /* printf("0 setgray\n"); */
    printf("%.2f %.2f moveto\n",xorigin - 10,yorigin + yheight);
    printf("-5 -2 rmoveto\n");
    printf("(%d) dup stringwidth pop neg 0 rmoveto show\n",(int) rint(ymax/binstep));
    printf("grestore\n");

    yfactor = yheight / (double) ymax;

    for (bini = 0; bini < nbins; bini++) {
      if (binx_lowend[bini] == 0 && binx_highend[bini] == 0) {
	/* Skip */
      } else {
	height1 = height2 = 0.0;
	if (monochromep == true) {
	  if (binx_lowend[bini] + binx_highend[bini] > 0) {
	    height2 = height1 + (binx_lowend[bini] + binx_highend[bini]) * yfactor;
	    print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
				DARK_BLUE);
	  }
	} else {
	  if (binx_lowend[bini] > 0) {
	    height2 = height1 + binx_lowend[bini] * yfactor;
	    print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
				BREWER_RED);
	    height1 = height2;
	  }
	  if (binx_highend[bini] > 0) {
	    height2 = height1 + binx_highend[bini] * yfactor;
	    print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
				BREWER_GREEN);
	    /* height1 = height2; */
	  }
	}
      }
    }
  }

  FREE(binx_highend);
  FREE(binx_lowend);

  return;
}


static void
value_print (double *values, Genomicpos_T mincoord, Genomicpos_T maxcoord,
	     double binstep, double xorigin, double ycenter, double yheight) {
  int nbins, bini, nextbini;
  double bin_low, bin_high, marker, width, pct;
  Genomicpos_T chrpos;

  double yfactor, height;
  double *binx, ymax = 0.0;

  nbins = (int) ((maxcoord - mincoord + 1)/binstep + 1);
  binx = (double *) CALLOC(nbins,sizeof(double));

  for (chrpos = mincoord; chrpos <= maxcoord; chrpos++) {
    bin_low = (double) (chrpos - mincoord)/binstep;
    bin_high = (double) (chrpos + 1U - mincoord)/binstep;
    width = bin_high - bin_low;

    marker = bin_low;
    bini = (int) marker;
    nextbini = bini + 1;
    while (nextbini < bin_high) {
      pct = ((double) nextbini - marker)/width;
      binx[bini] += pct * values[chrpos - mincoord];

      marker = nextbini;
      bini++;
      nextbini++;
    }

    pct = (bin_high - marker)/width;
    binx[bini] += pct * values[chrpos - mincoord];
  }

  for (bini = 0; bini < nbins; bini++) {
    if (binx[bini] > 0.0 && binx[bini] > ymax) {
      ymax = binx[bini];
    } else if (binx[bini] < 0.0 && -binx[bini] > ymax) {
      ymax = -binx[bini];
    }
  }

  if (ymax > 0) {
    yfactor = yheight / (double) ymax;

    for (bini = 0; bini < nbins; bini++) {
      if (binx[bini] != 0.0) {
	height = binx[bini] * yfactor;
	print_vertical_line(xorigin + bini,ycenter,xorigin + bini,ycenter + height,
			    GRAY);
      }
    }
  }

  FREE(binx);

  return;
}


static void
extents_print_one (long int *extents, Genomicpos_T chrlength,
		   Genomicpos_T mincoord, Genomicpos_T maxcoord,
		   double binstep, double xorigin, double yorigin, double yheight,
		   double red, double green, double blue) {
  Genomicpos_T chrpos;
  int nbins, bini, nextbini;
  double bin_low, bin_high, marker, width, pct;
  double *binx;
  double ymax = 0.0, y;
  double yfactor, height;

  nbins = (int) ((maxcoord - mincoord + 1)/binstep + 1);
  binx = (double *) CALLOC(nbins,sizeof(double));

  for (chrpos = mincoord; chrpos <= maxcoord; chrpos++) {
    bin_low = (double) (chrpos - mincoord)/binstep;
    bin_high = (double) (chrpos + 1U - mincoord)/binstep;
    width = bin_high - bin_low;

    marker = bin_low;
    bini = (int) marker;
    nextbini = bini + 1;
    while (nextbini < bin_high) {
      pct = ((double) nextbini - marker)/width;
      if (chrpos < chrlength) {
	binx[bini] += pct * extents[chrpos];
      }
      marker = nextbini;
      bini++;
      nextbini++;
    }

    pct = (bin_high - marker)/width;
    if (chrpos < chrlength) {
      binx[bini] += pct * extents[chrpos];
    }
  }
    

  for (bini = 0; bini < nbins; bini++) {
    if ((y = binx[bini]) > ymax) {
      ymax = y;
    }
  }

  if (ymax > 0) {
    scale_print(xorigin,yorigin,yheight,/*ydensity*/ymax/binstep);
    yfactor = yheight / (double) ymax;

    for (bini = 0; bini < nbins; bini++) {
      if (binx[bini] == 0) {
	/* Skip */
      } else {
	height = binx[bini] * yfactor;
	
	print_vertical_line(xorigin + bini,yorigin,xorigin + bini,yorigin + height,
			    red,green,blue);
      }
    }
  }
    
  FREE(binx);

  return;
}


static void
extents_print_dir (long int *fwd_extents, long int *rev_extents, long int *null_extents,
		   Genomicpos_T chrlength, Genomicpos_T mincoord, Genomicpos_T maxcoord,
		   double binstep, double xorigin, double yorigin, double yheight) {
  Genomicpos_T chrpos;
  int nbins, bini, nextbini;
  double bin_low, bin_high, marker, width, pct;
  double *binx_fwd, *binx_rev, *binx_null;
  double ymax = 0.0, y;
  double yfactor, height1, height2;
  

  nbins = (int) ((maxcoord - mincoord + 1)/binstep + 1);
  binx_fwd = (double *) CALLOC(nbins,sizeof(double));
  binx_rev = (double *) CALLOC(nbins,sizeof(double));
  binx_null = (double *) CALLOC(nbins,sizeof(double));

  for (chrpos = mincoord; chrpos <= maxcoord; chrpos++) {
    bin_low = (double) (chrpos - mincoord)/binstep;
    bin_high = (double) (chrpos + 1U - mincoord)/binstep;
    width = bin_high - bin_low;

    marker = bin_low;
    bini = (int) marker;
    nextbini = bini + 1;
    while (nextbini < bin_high) {
      pct = ((double) nextbini - marker)/width;
      if (fwd_extents && chrpos < chrlength) {
	binx_fwd[bini] += pct * fwd_extents[chrpos];
      }
      if (rev_extents && chrpos < chrlength) {
	binx_rev[bini] += pct * rev_extents[chrpos];
      }
      if (null_extents && chrpos < chrlength) {
	binx_null[bini] += pct * null_extents[chrpos];
      }

      marker = nextbini;
      bini++;
      nextbini++;
    }

    pct = (bin_high - marker)/width;
    if (fwd_extents && chrpos < chrlength) {
      binx_fwd[bini] += pct * fwd_extents[chrpos];
    }
    if (rev_extents && chrpos < chrlength) {
      binx_rev[bini] += pct * rev_extents[chrpos];
    }
    if (null_extents && chrpos < chrlength) {
      binx_null[bini] += pct * null_extents[chrpos];
    }
  }
    
  for (bini = 0; bini < nbins; bini++) {
    if ((y = binx_fwd[bini] + binx_rev[bini] + binx_null[bini]) > ymax) {
      ymax = y;
    }
  }

  if (ymax > 0) {
    scale_print(xorigin,yorigin,yheight,/*ydensity*/ymax/binstep);
    yfactor = yheight / (double) ymax;

    for (bini = 0; bini < nbins; bini++) {
      if (binx_fwd[bini] == 0 && binx_rev[bini] == 0 && binx_null[bini] == 0) {
	/* Skip */
      } else {
	height1 = height2 = 0.0;
	if (binx_fwd[bini] > 0) {
	  height2 = height1 + binx_fwd[bini] * yfactor;
	  print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
			      BLACK);
	  height1 = height2;
	}
	if (binx_null[bini] > 0) {
	  height2 = height1 + binx_null[bini] * yfactor;
	  print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
			      GRAY);
	  height1 = height2;
	}
	if (binx_rev[bini] > 0) {
	  height2 = height1 + binx_rev[bini] * yfactor;
	  print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
			      RED);
	  /* height1 = height2; */
	}
      }
    }
  }
    
  FREE(binx_null);
  FREE(binx_rev);
  FREE(binx_fwd);

  return;
}


static void
extents_print_crosshyb (long int *primary_extents, long int *crosshyb_extents,
			Genomicpos_T chrlength, Genomicpos_T mincoord, Genomicpos_T maxcoord,
			double binstep, double xorigin, double yorigin, double yheight) {
  Genomicpos_T chrpos;
  int nbins, bini, nextbini;
  double bin_low, bin_high, marker, width, pct;
  double *binx_primary, *binx_crosshyb;
  double ymax = 0.0, y;
  double yfactor, height1, height2;

  nbins = (int) ((maxcoord - mincoord + 1)/binstep + 1);
  binx_primary = (double *) CALLOC(nbins,sizeof(double));
  binx_crosshyb = (double *) CALLOC(nbins,sizeof(double));

  for (chrpos = mincoord; chrpos <= maxcoord; chrpos++) {
    bin_low = (double) (chrpos - mincoord)/binstep;
    bin_high = (double) (chrpos + 1U - mincoord)/binstep;
    width = bin_high - bin_low;

    marker = bin_low;
    bini = (int) marker;
    nextbini = bini + 1;
    while (nextbini < bin_high) {
      pct = ((double) nextbini - marker)/width;
      if (primary_extents && chrpos < chrlength) {
	binx_primary[bini] += pct * primary_extents[chrpos];
      }
      if (crosshyb_extents && chrpos < chrlength) {
	binx_crosshyb[bini] += pct * crosshyb_extents[chrpos];
      }

      marker = nextbini;
      bini++;
      nextbini++;
    }

    pct = (bin_high - marker)/width;
    if (primary_extents && chrpos < chrlength) {
      binx_primary[bini] += pct * primary_extents[chrpos];
    }
    if (crosshyb_extents && chrpos < chrlength) {
      binx_crosshyb[bini] += pct * crosshyb_extents[chrpos];
    }
  }
    
  for (bini = 0; bini < nbins; bini++) {
    if ((y = binx_primary[bini] + binx_crosshyb[bini]) > ymax) {
      ymax = y;
    }
  }

  if (ymax > 0) {
    scale_print(xorigin,yorigin,yheight,/*ydensity*/ymax/binstep);
    yfactor = yheight / (double) ymax;

    for (bini = 0; bini < nbins; bini++) {
      if (binx_primary[bini] == 0 && binx_crosshyb[bini] == 0) {
	/* Skip */
      } else {
	height1 = height2 = 0.0;
	if (binx_primary[bini] > 0) {
	  height2 = height1 + binx_primary[bini] * yfactor;
	  print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
			      DARK_BLUE);
	  height1 = height2;
	}
	if (binx_crosshyb[bini] > 0) {
	  height2 = height1 + binx_crosshyb[bini] * yfactor;
	  print_vertical_line(xorigin + bini,yorigin + height1,xorigin + bini,yorigin + height2,
			      BEIGE);
	  /* height1 = height2; */
	}
      }
    }
  }
    
  FREE(binx_crosshyb);
  FREE(binx_primary);

  return;
}



static void
splices_print (List_T splices, Genomicpos_T mincoord, Genomicpos_T maxcoord,
	       double left, double right, double yorigin, double yheight, double ydelta) {
  List_T p;
  double xpos1, xpos2, ypos;
  Splice_T splice;

  for (p = splices; p != NULL; p = List_next(p)) {
    splice = (Splice_T) List_head(p);
    if (Splice_level(splice) >= 0) {
      /* ypos = yorigin + yheight - ydelta * Splice_level(splice); */
      ypos = calc_ypos(/*lineindex*/Splice_level(splice) - 1,/*linesep*/ydelta,/*top*/yorigin + yheight,/*bottom*/0);
      xpos1 = calc_xpos(Splice_low(splice),maxcoord,mincoord,left,right);
      xpos2 = calc_xpos(Splice_high(splice),maxcoord,mincoord,left,right);

#if 0
      if (Splice_validp(splice) == false) {
	print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,GRAY);
      } else if (Splice_sign(splice) > 0) {
	print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,BLACK);
      } else {
	print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,RED);
      }
#else
      if (Splice_sign(splice) > 0) {
	print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,BLACK);
      } else if (Splice_sign(splice) < 0) {
	print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,RED);
      } else {
	print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,GRAY);
      }
#endif

      printf("gsave\n");
      /* printf("0 setgray\n"); */
      printf("%.2f %.2f moveto\n",xpos1,ypos);
      printf("0 5 rmoveto\n");
      printf("(%d) show\n",Splice_count(splice));
      printf("grestore\n");
    }
  }

  return;
}


/* Modified from revise_read in bamtally.c */
static void
print_seqdiffs (Genomicpos_T chrpos_low, Intlist_T types, Uintlist_T npositions,
		char *shortread, char *genomic, Genomicpos_T chrstart, Genomicpos_T chrend,
		double left, double right, double ypos1, double ypos2) {
  Genomicpos_T pos;
  char *p, *q;
  Intlist_T a;
  Uintlist_T b;
  unsigned int mlength;
  int type;
  double xpos, xpos1, xpos2;
  Genomicpos_T mincoord = chrstart, maxcoord = chrend;


  debug1(fprintf(stderr,"Printing read at %u\n",chrpos_low));

  pos = chrpos_low - 1U;		/* Bamread reports chrpos as 1-based */

  p = shortread;
  for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
    type = Intlist_head(a);
    mlength = Uintlist_head(b);
    if (type == 'S') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      p += mlength;
    } else if (type == 'H') {
      /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
      /* p += mlength; -- hard clip means query sequence is absent */
    } else if (type == 'N') {
      pos += mlength;

    } else if (type == 'P') {
      /* Phantom nucleotides that are inserted in the reference
	 without modifying the genomicpos.  Like a 'D' but leaves pos
	 unaltered. */

    } else if (type == 'I') {

      /* PRINT INSERTION AT POS */
      xpos1 = calc_xpos(pos,maxcoord,mincoord,left,right);
      xpos2 = calc_xpos(pos+mlength,maxcoord,mincoord,left,right);
      print_vertical_line(xpos1,ypos1,xpos2,ypos2,BLACK); /* Slopes upward right */

      p += mlength;

    } else if (type == 'D') {

      /* PRINT DELETION AT POS FOR MLENGTH (otherwise appears as a blank space) */
      xpos1 = calc_xpos(pos,maxcoord,mincoord,left,right);
      xpos2 = calc_xpos(pos+mlength,maxcoord,mincoord,left,right);
      print_vertical_line(xpos1,ypos2,xpos2,ypos1,BLACK); /* Slopes downward right */

      pos += mlength;

    } else if (type == 'M') {

      /* Assume we have genome from chrstart to chrend */
      debug1(fprintf(stderr,"Genomic pos is %u\n",pos));

      /* psave = p; qsave = q; */
      while (mlength-- > /* trimright */ 0) {
	if (pos < chrstart) {
	  debug1(fprintf(stderr,"Skipping pos %u\n",pos));
	} else if (pos >= chrend) {
	  debug1(fprintf(stderr,"Skipping pos %u\n",pos));
	} else {
	  q = &(genomic[pos - chrstart]);
	  debug1(fprintf(stderr,"Processing %c and %c at pos %u, mlength %u",*p,*q,pos+1U,mlength));
	  if (*p != toupper(*q)) {
	    debug1(fprintf(stderr," => difference"));
	    xpos = calc_xpos(pos,maxcoord,mincoord,left,right);
	    print_vertical_line(xpos,ypos1,xpos,ypos2,RED);
	  }
	  debug1(fprintf(stderr,"\n"));
	}

	p++;
	q++;
	pos++;
      }

    } else {
      fprintf(stderr,"Cannot parse type '%c'\n",type);
      exit(9);
    }
  }
  
  return;
}



/* mincoord is chrstart and maxcoord is chrend */
static void
bampairs_print (List_T bampairs, Genomicpos_T mincoord, Genomicpos_T maxcoord,
		double left, double right, double yorigin, double yheight, double ydelta,
		Genome_T genome, Genomicpos_T chroffset) {
  List_T p;
  double xpos1, xpos2, ypos, ypos1, ypos2;
  Bampair_T bampair;
  Bamline_T bamline;
  Uintlist_T chrpos_lows, chrpos_highs, a, b;
  Uintlist_T splice_lows, splice_highs;
  Intlist_T splice_signs, c;

  char *genomic;
  Genomicpos_T chrstart = mincoord, chrend = maxcoord;

  genomic = (char *) CALLOC(chrend - chrstart + 1 + 1,sizeof(char));
  Genome_fill_buffer_simple(genome,/*left*/chroffset+chrstart,chrend-chrstart+1,genomic);

  for (p = bampairs; p != NULL; p = List_next(p)) {
    bampair = (Bampair_T) List_head(p);
    if (Bampair_level(bampair) >= 0) {
      /* ypos = yorigin + yheight - ydelta * Bampair_level(bampair); */
      ypos = calc_ypos(/*lineindex*/Bampair_level(bampair) - 1,/*linesep*/ydelta,/*top*/yorigin + yheight,/*bottom*/0);
      xpos1 = calc_xpos(Bampair_chrpos_low(bampair),maxcoord,mincoord,left,right);
      xpos2 = calc_xpos(Bampair_chrpos_high(bampair),maxcoord,mincoord,left,right);

      print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,GRAY);

      ypos1 = ypos + 1.0;
      ypos2 = ypos - 1.0;

      Bampair_details(&chrpos_lows,&chrpos_highs,&splice_lows,&splice_highs,&splice_signs,bampair);

      /* Splices */
      a = splice_lows;
      b = splice_highs;
      c = splice_signs;
      while (a != NULL && b != NULL && c != NULL) {
	xpos1 = calc_xpos(Uintlist_head(a),maxcoord,mincoord,left,right);
	xpos2 = calc_xpos(Uintlist_head(b),maxcoord,mincoord,left,right);
	if (Intlist_head(c) > 0) {
	  print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,BLACK);
	} else if (Intlist_head(c) < 0) {
	  print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,RED);
	} else {
	  print_horizontal_line(xpos1,ypos,xpos2,ypos,left,right,GRAY);
	}
	a = Uintlist_next(a);
	b = Uintlist_next(b);
	c = Intlist_next(c);
      }
      Uintlist_free(&splice_lows);
      Uintlist_free(&splice_highs);
      Intlist_free(&splice_signs);

      /* Regions */
      a = chrpos_lows;
      b = chrpos_highs;
      while (a != NULL && b != NULL) {
	xpos1 = calc_xpos(Uintlist_head(a),maxcoord,mincoord,left,right);
	xpos2 = calc_xpos(Uintlist_head(b),maxcoord,mincoord,left,right);

	if (Bampair_primaryp(bampair) == true) {
	  print_rect(xpos1,ypos1,xpos2,ypos2,VERY_LIGHT_BLUE);
	} else if (Bampair_good_unique_p(bampair) == true) {
	  print_rect(xpos1,ypos1,xpos2,ypos2,LIGHT_BLUE);
	} else {
	  print_rect(xpos1,ypos1,xpos2,ypos2,BEIGE);
	}
	a = Uintlist_next(a);
	b = Uintlist_next(b);
      }
      Uintlist_free(&chrpos_lows);
      Uintlist_free(&chrpos_highs);

      if (print_seqdiffs_p == true) {
	/* Sequence differences */
	if ((bamline = Bampair_bamline_low(bampair)) != NULL) {
	  print_seqdiffs(Bamline_chrpos_low(bamline),
			 Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
			 Bamline_read(bamline),genomic,chrstart,chrend,left,right,
			 ypos1,ypos2);
	}
	
	if ((bamline = Bampair_bamline_high(bampair)) != NULL) {
	  print_seqdiffs(Bamline_chrpos_low(bamline),
			 Bamline_cigar_types(bamline),Bamline_cigar_npositions(bamline),
			 Bamline_read(bamline),genomic,chrstart,chrend,left,right,
			 ypos1,ypos2);
	}
      }
    }
  }

  FREE(genomic);

  return;
}



/************************************************************************
 *   Genes
 ************************************************************************/

#define GENEHEIGHT 4.0
#define GENES_YDELTA 14.0

#define COLOR_A 0.0,0.8,0.0	/* green */
#define COLOR_C 0.0,0.0,0.8	/* blue */
#define COLOR_G 1.0,0.733,0.0	/* orange */
#define COLOR_T 0.8,0.0,0.0	/* red */
#define COLOR_N 0.72,0.72,0.72	/* gray */


static IIT_T iit_for_sort;

static int
gene_compare (const void *x, const void *y) {
  int a = * (int *) x;
  int b = * (int *) y;
  unsigned int lowa, lowb;

  lowa = Interval_low(IIT_interval(iit_for_sort,a));
  lowb = Interval_low(IIT_interval(iit_for_sort,b));
  if (lowa < lowb) {
    return -1;
  } else if (lowa > lowb) {
    return +1;
  } else {
    return 0;
  }
}

#define WIDTHFORLABELS 20000001	/* 20 million nt */
#define WIDTHFOREXONS 10000001	/* 10 million nt */
#define PIXELSPERCHAR 6.0
#define PIXELSPERCHAR_ASCII 6.0 



static void
plot_exons (char *annotation, double ylevel, double xhigh, double xlow,
	    unsigned int maxcoord, unsigned int mincoord,
	    double left, double right, double red, double green, double blue) {
  char *p = annotation;
  unsigned int startpos, endpos;
  double xstart, xend;

  /* Skip first line */
  while (*p != '\0' && *p != '\n') {
    p++;
  }
  if (*p == '\n') {
    p++;
  }

  if (*p == '\0') {
    /* IIT file has no exon information */
    /* Solid gene */
    print_rect(xlow,ylevel,xhigh,ylevel+GENEHEIGHT,red,green,blue);
  } else {
    print_horizontal_line(xlow,ylevel+GENEHEIGHT/2.0,xhigh,ylevel+GENEHEIGHT/2.0,
			  left,right,red,green,blue);

    while (*p != '\0') {
      sscanf(p,"%u %u",&startpos,&endpos);
      xstart = calc_xpos(startpos,maxcoord,mincoord,left,right);
      xend = calc_xpos(endpos,maxcoord,mincoord,left,right);

      print_rect(xstart,ylevel,xend,ylevel+GENEHEIGHT,red,green,blue);

      /* Read to end of line */
      while (*p != '\0' && *p != '\n') {
	p++;
      }
      if (*p == '\n') {
	p++;
      }
    }
  }

  return;
}
    


static int
genes_print (IIT_T iit, char *chr, unsigned int mincoord, unsigned int maxcoord, int typeint,
	     unsigned int chrlength, double left, double right, double top, double bottom,
	     bool allgenesp) {
  int maxlevel_printed = -1, maxlevel = -1, *matches, nmatches, level, index, i, sign, type;
  Interval_T interval;
  unsigned int low, high, width, widthforlabels, widthforexons;
  char *label, *annotation, *restofheader, *typestring;
  bool donep, allocp, allocp2;

  double ylevel, xlow, xhigh, *rightmost, buffer;
  double red, green, blue;


  if (allgenesp == true) {
    if (typeint > 0) {
      matches = IIT_get_typed(&nmatches,iit,chr,1U,chrlength,typeint,/*sortp*/false);
    } else {
      matches = IIT_get(&nmatches,iit,chr,1U,chrlength,/*sortp*/false);
    }
    widthforlabels = 3*WIDTHFORLABELS;
    widthforexons = 3*WIDTHFOREXONS;
    fprintf(stderr,"%d genes found in entire chromosome\n",nmatches);
  } else {
    if (typeint > 0) {
      matches = IIT_get_typed(&nmatches,iit,chr,mincoord,maxcoord,typeint,/*sortp*/false);
    } else {
      matches = IIT_get(&nmatches,iit,chr,mincoord,maxcoord,/*sortp*/false);
    }
    widthforlabels = WIDTHFORLABELS;
    widthforexons = WIDTHFOREXONS;
    fprintf(stderr,"%d genes found in range\n",nmatches);
  }

  if (nmatches == 0) {
    return 0;
  }
  rightmost = (double *) CALLOC(nmatches,sizeof(double));
  for (i = 0; i < nmatches; i++) {
    rightmost[i] = -10000000.0;
  }

  width = maxcoord - (mincoord + 1U);


  iit_for_sort = iit;
  qsort(matches,nmatches,sizeof(int),gene_compare);

  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = IIT_interval(iit,index);
    annotation = IIT_annotation(&restofheader,iit,index,&allocp);
    low = Interval_low(interval);
    high = Interval_high(interval);
    if (allgenesp == false) {
      if (low < mincoord) {
	low = mincoord;
      }
      if (high > maxcoord) {
	high = maxcoord;
      }
    }
    xlow = calc_xpos(low,maxcoord,mincoord,left,right);
    xhigh = calc_xpos(high,maxcoord,mincoord,left,right);

    /* Find appropriate level */
    level = 0;
    donep = false;
    while (!donep) {
      if (level > maxlevel) {
	donep = true;
	maxlevel = level;
      } else if (rightmost[level] < xlow) {
	donep = true;
      } else {
	level++;
      }
    }
    
    if (width > widthforlabels) {
      rightmost[level] = xhigh + 3.0;
    } else if ((label = IIT_label(iit,index,&allocp2)) == NULL) {
      rightmost[level] = xhigh + 3.0;
    } else if (label[0] == '\0') {
      rightmost[level] = xhigh + 3.0;
      if (allocp2 == true) {
	FREE(label);
      }

    } else {
      buffer = strlen(label)*PIXELSPERCHAR;

      if (xlow + buffer > xhigh + 3.0) {
	rightmost[level] = xlow + buffer;
      } else {
	rightmost[level] = xhigh + 3.0;
      }
    }

    debug(fprintf(stderr,"%d: For match %d at %u..%u, assigned to level %d and altered rightmost at that level to %f\n",
		  i,index,low,high,level,rightmost[level]));

    ylevel = calc_ypos(level,GENES_YDELTA,top,bottom);

    if (xhigh < xlow) {
      fprintf(stderr,"low = %u => xlow %f\n",low,xlow);
      fprintf(stderr,"high = %u => xhigh %f\n",high,xhigh);
      abort();
    } else if (xlow < right && xhigh > left) {
      if (IIT_version(iit) < 2) {
	/* In version 1 IITs, FWD/REV types used to indicate sign */
	type = Interval_type(interval);
	typestring = IIT_typestring(iit,type);
	if (!strcmp(typestring,"FWD")) {
	  red = 0.0; green = 0.0, blue = 0.0;
	} else if (!strcmp(typestring,"REV")) {
	  red = 1.0; green = 0.0, blue = 0.0;
	} else {
	  red = 0.0; green = 1.0, blue = 0.0;
	}
      } else {
	/* IIT version 2 or greater */
	if ((sign = Interval_sign(interval)) > 0) {
	  red = 0.0; green = 0.0, blue = 0.0;
	} else if (sign < 0) {
	  red = 1.0; green = 0.0, blue = 0.0;
	} else {
	  red = 0.0; green = 1.0, blue = 0.0;
	}
      }

      if (level > maxlevel_printed) {
	maxlevel_printed = level;
      }

      if (width <= widthforexons) {
	plot_exons(annotation,ylevel,xhigh,xlow,maxcoord,mincoord,left,right,
		   red,green,blue);
      } else {
	/* Solid gene */
	print_rect(xlow,ylevel,xhigh,ylevel+GENEHEIGHT,red,green,blue);
      }
      debug(fprintf(stderr,"%d (%u..%u): level %d\n",i,low,high,level));

      if (width <= widthforlabels) {
	if ((label = IIT_label(iit,index,&allocp2)) != NULL && label[0] != '\0') {
	  printf("gsave\n");
	  /* printf("0 setgray\n"); */
	  printf("%.2f %.2f moveto\n",xlow,ylevel);
	  printf("0 5 rmoveto\n");
	  printf("(%s) show\n",label);
	  printf("grestore\n");
	}
	if (allocp2 == true) {
	  FREE(label);
	}
      }
    }

    if (allocp == true) {
      FREE(restofheader);
    }

  }

  FREE(rightmost);
  FREE(matches);

  return maxlevel_printed + 1;
}


/* Example: 55000000 56000000 0.75 */

static void
plot_fusions (FILE *fp, char *chr, unsigned int maxcoord, unsigned int mincoord,
	      double left, double right, double yorigin, double yheight,
	      double red, double green, double blue) {
  unsigned int startpos, endpos;
  double xstart, xend, yfrac, log_yfrac, count, intensity, log_intensity;
  char Buffer[1024], coord1[1024], coord2[1024], *p, *q;

  while (fgets(Buffer,1024,fp) != NULL) {
    if (sscanf(Buffer,"%s %s %lf",coord1,coord2,&count) < 3) {
      /* Skip */
    } else {
      p = coord1;
      while (*p != '\0' && *p != ':') {
	p++;
      }
      q = coord2;
      while (*q != '\0' && *q != ':') {
	q++;
      }
      if (*p == ':' && *q == ':' && !strncmp(coord1,chr,p-coord1) && !strncmp(coord2,chr,q-coord2)) {
	p++;
	q++;
	if (sscanf(p,"%u",&startpos) == 1 && sscanf(q,"%u",&endpos) == 1 &&
	    startpos >= mincoord && startpos <= maxcoord &&
	    endpos >= mincoord && endpos <= maxcoord) {
	  xstart = calc_xpos(startpos,maxcoord,mincoord,left,right);
	  xend = calc_xpos(endpos,maxcoord,mincoord,left,right);
	  if (endpos > startpos) {
	    yfrac = (double) (endpos - startpos)/(double) (maxcoord - mincoord);
	  } else {
	    yfrac = (double) (startpos - endpos)/(double) (maxcoord - mincoord);
	  }
	  log_yfrac = log(yfrac*E + 1.0 - yfrac);

	  intensity = count/400.0;
	  if (intensity > 1.0) {
	    intensity = 1.0;
	  }
	  log_intensity = log(intensity*E + 1.0 - intensity);
	  print_arc(xstart,xend,yorigin,/*yheight*/yheight*log_yfrac,
		    log_intensity*red + 1.0 - log_intensity,
		    log_intensity*green + 1.0 - log_intensity,
		    log_intensity*blue + 1.0 - log_intensity);
	  fprintf(stderr,"%s",Buffer);
	}
      }
    }
  }

  return;
}


#if 0
static void
compute_cum_int (long int *cum, long int *x, long int *y, Genomicpos_T chrlength) {
  Genomicpos_T chrpos;

  cum[0] = 0.0;
  for (chrpos = 1; chrpos <= chrlength; chrpos++) {
    cum[chrpos] = x[chrpos] + y[chrpos] + cum[chrpos - 1U];
  }
  return;
}
#endif

static void
compute_log_tally (double *log_tally, long int *x, long int *y, Genomicpos_T chrlength) {
  Genomicpos_T chrpos;

  /* For non-cum results, want < and not <= */
  for (chrpos = 0; chrpos < chrlength; chrpos++) {
    log_tally[chrpos] = log((double) (x[chrpos]+y[chrpos]+1));
  }
  
  return;
}

static void
compute_cum_double (double *cum, double *x, Genomicpos_T chrlength) {
  Genomicpos_T chrpos;

  cum[0] = 0.0;
  for (chrpos = 1; chrpos <= chrlength; chrpos++) {
    cum[chrpos] = x[chrpos] + cum[chrpos - 1U];
  }
  return;
}

static void
compute_diff (double *window_diff, double *cumlog_tally, Genomicpos_T chrlength) {
  Genomicpos_T chrpos;
  double sum_right, sum_left;

  for (chrpos = 30U + 1; chrpos <= chrlength - 30U - 1; chrpos++) {
    sum_right = cumlog_tally[chrpos + 30U - 1] - cumlog_tally[chrpos - 1];
    sum_left = cumlog_tally[chrpos - 1] - cumlog_tally[chrpos - 30U - 1];
    window_diff[chrpos] = sum_right - sum_left;
    /* window_sumx[(int) chrpos] - window_sumx[(int) (chrpos - 30U)]; */
#if 0
    if (window_diff[chrpos] != 0.0) {
      printf("Putting diff %.1f into chrpos %u\n",window_diff[chrpos],chrpos);
    }
#endif
  }

  return;
}


int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *fileroot = NULL, *mapdir = NULL;
  char *map_iitfile;
  IIT_T map_iit;
  FILE *fp;

  double binstep;
  long int *tally_matches = NULL, *tally_mismatches = NULL,
    *tally_matches_low = NULL, *tally_mismatches_low = NULL,
    *tally_matches_high = NULL, *tally_mismatches_high = NULL;
  List_T intervallist = NULL, labellist = NULL, datalist = NULL, p, t;
  double *window_diff, *cumlog_tally, *log_tally;
  int quality_counts_match[256], quality_counts_mismatch[256], i;

  Genomicpos_T chroffset = 0U, chrstart, chrend, chrlength;

  char *iitfile;
  IIT_T chromosome_iit;

  Genome_T genome = NULL;
  List_T bamfiles;
  Bamreader_T bamreader;
  int readlength, insertlength;
  Gstruct_T gstruct = NULL;
  List_T splices = NULL, bampairs = NULL;
  Bampair_T bampair;
  int nlevels;

  long int *fwd_extents = NULL, *rev_extents = NULL, *null_extents = NULL,
    *primary_extents = NULL, *crosshyb_extents = NULL;

  char *chr;
  Genomicpos_T genomicstart, genomiclength;
  bool revcomp;
  
  int xorigin, yorigin;
  double xfactor;
  int filei;

  int opt, c;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char **argstart;


  while ((opt = getopt_long(argc,argv,"D:d:a:t:f:q:n:C:U:P:M:m:p:s:",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else if (!strcmp(long_name,"sep-extents")) {
	separate_extents_p = true;
      } else if (!strcmp(long_name,"no-seqdiffs")) {
	print_seqdiffs_p = false;

      } else if (!strcmp(long_name,"depth")) {
	min_depth = atoi(optarg);
      } else if (!strcmp(long_name,"ignore-duplicates")) {
	ignore_duplicates_p = true;
      } else if (!strcmp(long_name,"allow-duplicates")) {
	ignore_duplicates_p = false;

      } else if (!strcmp(long_name,"read-group")) {
	desired_read_group = optarg;

#if 0
      } else if (!strcmp(long_name,"use-quality-const")) {
	quality_score_constant = atoi(optarg);
	if (quality_score_constant > MAX_QUALITY_SCORE) {
	  fprintf(stderr,"quality substitution score %d is > maximum %d\n",
		  quality_score_constant,MAX_QUALITY_SCORE);
	  exit(9);
	}
#endif
      } else if (!strcmp(long_name,"monochrome")) {
	monochromep = true;

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case 'a': arcfile = optarg; break;
    case 't': tracks = List_from_string(optarg); break;
    case 'f': fontsize = atoi(optarg); break;

    case 'M': user_mapdir = optarg; break;
    case 'm': map_iitfiles = List_from_string(optarg); break;

    case 'q': minimum_mapq = atoi(optarg); break;
    case 'n': maximum_nhits = atoi(optarg); break;

    case 'C':
      switch (atoi(optarg)) {
      case 0: need_concordant_p = false; break;
      case 1: need_concordant_p = true; break;
      default: fprintf(stderr,"Concordant mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'U':
      switch (atoi(optarg)) {
      case 0: need_unique_p = false; break;
      case 1: need_unique_p = true; break;
      default: fprintf(stderr,"Unique mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'P':
      switch (atoi(optarg)) {
      case 0: need_primary_p = false; break;
      case 1: need_primary_p = true; break;
      default: fprintf(stderr,"Primary mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'p': alloclength = strtoul(optarg,NULL,10); break;

    case 's': mincount = atoi(optarg); break;

    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d flag\n");
    print_program_usage();
    exit(9);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
    FREE(dbversion);
    FREE(dbroot);
  }

  if (tracks == NULL) {
    tracks = List_push(tracks,(void *) "tallies");
    tracks = List_push(tracks,(void *) "splices");
    tracks = List_push(tracks,(void *) "extents");
    tracks = List_push(tracks,(void *) "genes");
    tracks = List_push(tracks,(void *) "reads");
    tracks = List_reverse(tracks);
  }

  fprintf(stderr,"Starting to allocate memory for %u positions\n",alloclength);


  for (i = 0; i < 256; i++) {
    quality_counts_match[i] = 0;
    quality_counts_mismatch[i] = 0;
  }

  /* Need genome to determine wild-type, because "known gene" may not match reference genome */
  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);
  Maxent_hr_setup(Genome_blocks(genome));

  /* Start page */
  printf("%%!PS-Adobe-2.0\n");
  printf("%%%%Creator: bam_print program, written by Thomas Wu, Genentech, Inc.\n");
  printf("%%%%Pages: 1\n");
  printf("%%%%BoundingBox: 0 0 612 792\n");
  printf("%%%%PageOrientation: Portrait\n");
  printf("%%%%Orientation: Portrait\n");
  printf("%%%%DocumentFonts: Helvetica\n");
  printf("%%%%EndComments\n");
  printf("%%%%BeginProcSet\n");
  printf("%%%%EndProcSet\n");
  printf("%%%%EndProlog\n");
  printf("%%%%BeginSetup\n");
  printf("%%%%EndSetup\n");

  printf("0.5 setlinewidth\n");
  printf("/Helvetica findfont %d scalefont setfont\n",fontsize);

  print_arc_procedure();


  xorigin = xbase;
  yorigin = ybase + pageheight;

#if 0
  if (argc < 2) {
    /* We now allow printing of known genes, without any data */
    fprintf(stderr,"argc = %d\n",argc);
    fprintf(stderr,"Usage: bam_print -d <genome> [-t <tracks>] [-m <genes iit, comma delimited>] <coords> <bam file...>\n");
    exit(9);
  }
#endif

  /* Print arguments */
  /* printf("0 setgray\n"); */
  printf("gsave\n");
  printf("%d %d moveto\n",xbase,argsheight);
  /* printf("0 5 rmoveto\n"); */

  printf("(bam_print called with args");
  argstart = &(argv[-optind]);
  for (c = 1; c < argc + optind; c++) {
    printf(" %s",argstart[c]);
  }
  printf(") show\n");
  printf("grestore\n");


  Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
		       &chroffset,&chrlength,argv[0],genomesubdir,fileroot);
  fprintf(stderr,"GMAP index says chrlength is %u\n",chrlength);
  binstep = (double) (chrend - chrstart + 1)/(double) pagewidth;

  xfactor = (double) pagewidth / (double) (chrend - chrstart + 1);

  yorigin -= AXIS_YHEIGHT;
  axis_print(chrstart,chrend,xorigin,/*left*/xorigin,/*right*/xorigin+pagewidth,yorigin);
  yorigin -= SPACER_YHEIGHT;

  if (arcfile != NULL) {
    if ((fp = fopen(arcfile,"r")) == NULL) {
      fprintf(stderr,"Cannot open arc file %s\n",arcfile);
      exit(9);
    }
    yorigin -= 200;
    plot_fusions(fp,chromosome,/*maxcoord*/chrend,/*mincoord*/chrstart,
		 /*left*/xorigin,/*right*/xorigin+pagewidth,yorigin,/*yheight*/200,
		 0,0,0);
    fclose(fp);
    yorigin -= SPACER_YHEIGHT;
  }


  for (t = tracks; t != NULL; t = List_next(t)) {
    if (!strcmp(List_head(t),"genes")) {
      mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
      for (p = map_iitfiles; p != NULL; p = List_next(p)) {
	map_iitfile = (char *) List_head(p);
	if ((map_iit = IIT_read(map_iitfile,/*name*/NULL,true,/*divread*/READ_ONE,/*divstring*/chromosome,
				/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	  iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				    strlen(List_head(p))+1,sizeof(char));
	  sprintf(iitfile,"%s/%s",mapdir,map_iitfile);
	  if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_ONE,/*divstring*/chromosome,
				  /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	    fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
	  }
	  FREE(iitfile);
	}
	if (map_iit != NULL) {
	  nlevels = genes_print(map_iit,chromosome,chrstart,chrend,
				/*typeint*/0,chrlength,/*left*/xorigin,/*right*/xorigin+pagewidth,
				/*top*/yorigin,/*bottom*/yorigin,/*allgenesp*/false);
	  fprintf(stderr,"gene levels: %d\n",nlevels);
	  yorigin -= nlevels * GENES_YDELTA;
	  yorigin -= SPACER_YHEIGHT;
	  IIT_free(&map_iit);
	}
	FREE(map_iitfile);
      }
      List_free(&map_iitfiles);
      FREE(mapdir);
    }
  }

  for (t = tracks; t != NULL; t = List_next(t)) {
    if (!strcmp(List_head(t),"extents")) {
      need_gstruct_p = true;
    } else if (!strcmp(List_head(t),"splices")) {
      need_gstruct_p = true;
    }
  }

  if (need_gstruct_p == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);
  }

  for (filei = 1; filei < argc; filei++) {
    bamfiles = List_push(NULL,(void *) argv[filei]);
  }

  if (need_gstruct_p == true) {
    gstruct = Gstruct_bam_input(&readlength,&insertlength,bamfiles,
				chromosome,chrstart,chrend,
				/*ngoodhits_low_table*/NULL,/*ngoodhits_high_table*/NULL,
				/*genes_iit*/NULL,shortsplicedist,max_pairlength,genome,chromosome_iit,
				desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,need_unique_p,need_primary_p,
				ignore_duplicates_p,/*trust_sam_p*/true,need_canonical_p,
				/*bam_lacks_chr*/NULL,/*bam_lacks_chr_length*/0);

    if ((chr = Gstruct_next_chr(&splices,&fwd_extents,&rev_extents,&null_extents,
				&primary_extents,&crosshyb_extents,&chrlength,
				gstruct,max_exonlength,mincount_end_alt,
				minsupport,need_canonical_p,
				/*bam_lacks_chr*/NULL,/*bam_lacks_chr_length*/0)) == NULL) {
      fprintf(stderr,"No data in given region\n");
      exit(0);
    } else {
      fprintf(stderr,"BAM file says chrlength of chr %s is %u\n",chr,chrlength);
    }
  }

  for (filei = 1; filei < argc; filei++) {
    bamreader = Bamread_new(argv[filei]);
    for (t = tracks; t != NULL; t = List_next(t)) {
      if (!strcmp(List_head(t),"tallies")) {
	
	/* One tally */
	if (separate_tallies_p == false) {
	  Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
	  Bamtally_run(&tally_matches,&tally_mismatches,
		       &intervallist,&labellist,&datalist,
		       quality_counts_match,quality_counts_mismatch,
		       bamreader,genome,chromosome,chroffset,chrstart,chrend,/*map_iit*/NULL,
		       alloclength,/*resolve_low_table*/NULL,/*resolve_high_table*/NULL,
		       desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,need_concordant_p,
		       need_unique_p,need_primary_p,ignore_duplicates_p,
		       /*ignore_lowend_p*/false,/*ignore_highend_p*/false,
		       /*output_type*/OUTPUT_TALLY,/*blockp*/false,blocksize,
		       quality_score_adj,min_depth,/*variant_strands*/0,
		       /*genomic_diff_p*/false,/*signed_counts_p*/false,
		       /*ignore_query_Ns_p*/true,/*print_indels_p*/false,
		       /*print_totals_p*/false,/*print_cycles_p*/false,/*print_quality_scores_p*/false,
		       /*print_mapq_scores_p*/false,/*print_xs_scores_p*/false,/*want_genotypes_p*/false,
		       /*verbosep*/true,/*readlevel_p*/false,/*max_softclip*/0,/*print_noncovered_p*/false,
		       /*bamfile*/NULL);
	  Bamread_unlimit_region(bamreader);

	  yorigin -= TALLY_YHEIGHT;
	  tally_print(tally_matches,tally_mismatches,chrstart,chrend,binstep,
		      xorigin,yorigin,/*yheight*/TALLY_YHEIGHT);
	  yorigin -= SPACER_YHEIGHT;

	} else {
	  Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
	  Bamtally_run_lh(&tally_matches_low,&tally_mismatches_low,
			  &tally_matches_high,&tally_mismatches_high,
			  quality_counts_match,quality_counts_mismatch,
			  bamreader,genome,chromosome,chroffset,/*chrstart*/0U,/*chrend*/chrlength,
			  /*map_iit*/NULL,alloclength,
			  desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,need_concordant_p,
			  need_unique_p,need_primary_p,ignore_duplicates_p,
			  blocksize,quality_score_adj,min_depth,/*variant_strands*/0,/*genomic_diff_p*/false,
			  /*ignore_query_Ns_p*/true,/*verbosep*/true,/*readlevel_p*/false,
			  /*max_softclip*/0,/*print_noncovered_p*/false);
	  Bamread_unlimit_region(bamreader);

	  yorigin -= TALLY_YHEIGHT;
	  tally_print_lh(tally_matches_low,tally_mismatches_low,tally_matches_high,tally_mismatches_high,
			 chrstart,chrend,binstep,xorigin,yorigin,/*yheight*/TALLY_YHEIGHT);
	  yorigin -= SPACER_YHEIGHT;
	}

      } else if (!strcmp(List_head(t),"edges")) {
	/* Edges */
	log_tally = (double *) CALLOC(chrlength+1,sizeof(double));
	compute_log_tally(log_tally,tally_matches,tally_mismatches,chrlength);

	cumlog_tally = (double *) CALLOC(chrlength+1,sizeof(double));
	compute_cum_double(cumlog_tally,log_tally,chrlength);
	FREE(log_tally);

	window_diff = (double *) CALLOC(chrlength+1,sizeof(double));
	compute_diff(window_diff,cumlog_tally,chrlength);
	FREE(cumlog_tally);

	yorigin -= EDGES_YHEIGHT;
	value_print(window_diff,chrstart,chrend,binstep,
		    xorigin,/*ycenter*/yorigin+EDGES_YHEIGHT/2,/*yheight*/EDGES_YHEIGHT/2);
	yorigin -= SPACER_YHEIGHT;
	FREE(window_diff);

      } else if (!strcmp(List_head(t),"extents")) {
	yorigin -= EXTENTS_YHEIGHT;
	extents_print_crosshyb(primary_extents,crosshyb_extents,chrlength,
			       chrstart,chrend,binstep,xorigin,yorigin,/*yheight*/EXTENTS_YHEIGHT);
	yorigin -= SPACER_YHEIGHT;
      
	if (separate_extents_p == true) {
	  yorigin -= EXTENTS_YHEIGHT;
	  extents_print_one(fwd_extents,chrlength,chrstart,chrend,binstep,xorigin,yorigin,/*yheight*/EXTENTS_YHEIGHT,
			    BLACK);
	  yorigin -= SPACER_YHEIGHT;
	
	  yorigin -= EXTENTS_YHEIGHT;
	  extents_print_one(rev_extents,chrlength,chrstart,chrend,binstep,xorigin,yorigin,/*yheight*/EXTENTS_YHEIGHT,
			    RED);
	  yorigin -= SPACER_YHEIGHT;
	
	  yorigin -= EXTENTS_YHEIGHT;
	  extents_print_one(null_extents,chrlength,chrstart,chrend,binstep,xorigin,yorigin,/*yheight*/EXTENTS_YHEIGHT,
			    GRAY);
	  yorigin -= SPACER_YHEIGHT;

	} else {
	  yorigin -= EXTENTS_YHEIGHT;
	  extents_print_dir(fwd_extents,rev_extents,null_extents,chrlength,
			    chrstart,chrend,binstep,xorigin,yorigin,/*yheight*/EXTENTS_YHEIGHT);
	  yorigin -= SPACER_YHEIGHT;
	}

      } else if (!strcmp(List_head(t),"splices")) {
	nlevels = Splice_compute_levels(splices,chrstart,chrend,/*max_allowed_levels*/50,xfactor,mincount);
	fprintf(stderr,"splice levels: %d\n",nlevels);

	yorigin -= nlevels * SPLICES_YDELTA;
	splices_print(splices,chrstart,chrend,/*left*/xorigin,/*right*/xorigin+pagewidth,
		      yorigin,/*yheight*/nlevels * SPLICES_YDELTA,/*ydelta*/SPLICES_YDELTA);
	yorigin -= SPACER_YHEIGHT;

      } else if (!strcmp(List_head(t),"genes")) {
	/* Already handled above */

      } else if (!strcmp(List_head(t),"reads")) {
	Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
	bampairs = Bamread_all_pairs(bamreader,desired_read_group,minimum_mapq,good_unique_mapq,maximum_nhits,
				     need_unique_p,need_primary_p,ignore_duplicates_p,
				     need_concordant_p);
	Bamread_unlimit_region(bamreader);

	nlevels = Bampair_compute_levels(bampairs,chrstart,chrend,/*max_allowed_levels*/200,xfactor,
					 min_pairlength,/*only_internal_p*/true);
	fprintf(stderr,"read levels: %d\n",nlevels);

	yorigin -= SPACER_YHEIGHT;
	yorigin -= nlevels * RAWREADS_YDELTA;
	bampairs_print(bampairs,chrstart,chrend,/*left*/xorigin,/*right*/xorigin+pagewidth,
		       yorigin,/*yheight*/nlevels * RAWREADS_YDELTA,/*ydelta*/RAWREADS_YDELTA,
		       genome,chroffset);
	yorigin -= SPACER_YHEIGHT;

	for (p = bampairs; p != NULL; p = List_next(p)) {
	  bampair = (Bampair_T) List_head(p);
	  Bampair_free(&bampair);
	}
	List_free(&bampairs);

      } else {
	fprintf(stderr,"Don't recognize track %s.  Allowed values: tallies, splices, extents, genes, reads\n",
		(char *) List_head(t));
      }
    }

    FREE(tally_mismatches);
    FREE(tally_matches);

    FREE(tally_mismatches_high);
    FREE(tally_matches_high);
    FREE(tally_mismatches_low);
    FREE(tally_matches_low);

    if (gstruct != NULL) {
      FREE(fwd_extents);
      FREE(rev_extents);
      FREE(null_extents);
      FREE(primary_extents);
      FREE(crosshyb_extents);
      List_free(&splices);
      Gstruct_free(&gstruct);
    }

    Bamread_free(&bamreader);
  }

  List_free(&bamfiles);

  if (need_gstruct_p == true) {
    IIT_free(&chromosome_iit);
  }

  if (genome != NULL) {
    Genome_free(&genome);
  }

  FREE(fileroot);
  FREE(genomesubdir);
  List_free(&tracks);


#if 0
  /* Should have already been freed during printing */
  for (pos = 0; pos < alloclength; pos++) {
    for (p = mismatches_byshift[pos]; p != NULL; p = List_next(p)) {
      mismatch = (Mismatch_T) List_head(p);
      Mismatch_free(&mismatch);
    }
    List_free(&(mismatches_byshift[pos]));
  }
  FREE(mismatches_byshift);

  for (pos = 0; pos < alloclength; pos++) {
    for (p = mismatches_byquality[pos]; p != NULL; p = List_next(p)) {
      mismatch = (Mismatch_T) List_head(p);
      Mismatch_free(&mismatch);
    }
    List_free(&(mismatches_byquality[pos]));
  }
  FREE(mismatches_byquality);

  for (pos = 0; pos < alloclength; pos++) {
    for (p = matches_byshift[pos]; p != NULL; p = List_next(p)) {
      match = (Match_T) List_head(p);
      Match_free(&match);
    }
    List_free(&(matches_byshift[pos]));
  }
  FREE(matches_byshift);

  for (pos = 0; pos < alloclength; pos++) {
    for (p = matches_byquality[pos]; p != NULL; p = List_next(p)) {
      match = (Match_T) List_head(p);
      Match_free(&match);
    }
    List_free(&(matches_byquality[pos]));
  }
  FREE(matches_byquality);
#endif

  
  /* Finish page */
  printf("showpage\n");
  printf("%%%%Trailer\n");
  printf("%%%%EOF\n");

  return 0;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: bam_print [OPTIONS...] chromosome:range bamfile...  or\n\
       bam_print [OPTIONS...] chromosome: bamfile... \n\
\n\
where\n\
   range is startposition..endposition\n\
         or startposition+length (+ strand)\n\
\n\
Input options (must include -d)\n\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
  -M, --mapdir=STRING            Map file directory.  Program will for map file as given and also here.\n\
  -m, --map=STRING               Map file(s), comma-delimited\n\
\n\
Compute options\n\
  --read-group=STRING            Require alignments to have this read group in the RG field of\n\
                                   the BAM line\n\
  -q, --min-mapq=INT             Require alignments to have this mapping quality and higher\n\
                                   (default 0)\n\
  -n, --nhits=INT                Require alignments to have this many hits or fewer\n\
                                   (default: 1000000)\n\
  -C, --concordant=INT           Require alignments to be concordant (0=no [default], 1=yes)\n\
  -U, --unique=INT               Require alignments to be unique (0=no [default], 1=yes)\n\
  -P, --primary=INT              Require alignments to be primary (0=no [default], 1=yes)\n\
  --allow-duplicates             Allow alignments even if marked as duplicate (0x400) [default behavior]\n\
  --ignore-duplicates            Ignore alignments marked as duplicate (0x400)\n\
  --pairmax=INT                  Expected insert length (reserves memory for this amount, so\n\
                                   alignments longer than this value are discarded) [default=400000]\n\
\n\
Filtering of output (options may be combined)\n\
  --depth=INT                    Print only positions with this depth or more\n\
\n\
Output options\n\
  -t, --tracks=STRING            Comma-separated list of tracks to print.\n\
                                   Allowed values: tallies, splices, extents, genes, reads\n\
                                   genes will be taken from values for -m and -M\n\
  -s, --mincount=INT             Minimum count for splices (default 1)\n\
  -f, --fontsize=INT             Font size for text\n\
  --sep-extents                  Print fwd, rev, and null extents as separate tracks\n\
  --no-seqdiffs                  Do not print ticks at sequence differences in reads track\n\
\n\
");
    return;
}
