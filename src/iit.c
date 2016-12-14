#include <stdlib.h>
#include <gstruct/iit-read.h>

#include "iit.h"

static void R_IIT_free(SEXP iit_R) {
  IIT_T iit = R_ExternalPtrAddr(iit_R);
  IIT_free(&iit);
}

SEXP R_IIT_new(IIT_T iit) {
  SEXP iit_R;
  iit_R = R_MakeExternalPtr((void *) iit, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(iit_R, R_IIT_free);
  return iit_R;
}

/*
 IIT API:
 * IIT_read
 * IIT_annotation
 * IIT_interval(iit, which), Interval_low, Interval_length, Interval_sign

 Coordinate mapping (read divs, not labels):
   Any overlap:
   * IIT_get(&(*nmatches),*iit,*divstring,*coordstart,*coordend,sortp);
   By strand:
   * IIT_get_signed()
   Exact coordinate match:
   * IIT_get_exact_multiple()
 
 Can restrict by "type" of feature:
 * typeint = IIT_typeint(*iit,typestring)
 * IIT_get_typed(&(*nmatches),*iit,*divstring,*coordstart,*coordend,typeint,
                 sortp);
 * IIT_get_typed_signed()

 Restrict by label (read labels, not divs):
 * IIT_find(&(*nmatches),iit,query);

 Read standard annotation:
 * field = IIT_annotation(&restofheader,iit,index,&allocp);
 Or a specific field:
 * fieldint = IIT_fieldint(iit,fieldstring)
 * field = IIT_fieldvalue(iit,index,fieldint);

 Get labels:
 * IIT_label(iit, index, &allocp)
 
 Introspection: 
 Metadata:
 * IIT_name
 * IIT_version
 * IIT_total_nintervals
 Types:
 * IIT_ntypes(iit)
 * IIT_typestring (iit, typeint);
 Fields:
 * IIT_nfields(iit)
 * IIT_fieldstring(iit, fieldint);
*/

typedef struct IITMatches {
    IIT_T iit;
    int **subscripts;
    int *nsubscripts;
    int nqueries;
} IITMatches;

static IITMatches _new_IITMatches(IIT_T iit, int nqueries) {
    IITMatches matches;
    matches.iit = iit;
    matches.subscripts = (int **) R_alloc(sizeof(int *), nqueries);
    matches.nsubscripts = (int *) R_alloc(sizeof(int), nqueries);
    matches.nqueries = nqueries;
    return matches;
}

SEXP R_iit_open(SEXP iitfile_R, SEXP divread_R, SEXP labels_read_R) {
  char *iitfile = (char *) CHAR(asChar(iitfile_R));
  int divread = asLogical(divread_R) ? READ_ALL : READ_NONE;
  bool labels_read = asLogical(labels_read_R);
  IIT_T iit = IIT_read(iitfile, /*name*/NULL, /*readonlyp*/true,
                       divread, /*divstring*/NULL,
                       /*add_iit_p*/false, labels_read);
  return R_IIT_new(iit);
}

static IITMatches _iit_find(IIT_T iit, SEXP which_R) {
    IITMatches matches = _new_IITMatches(iit, length(which_R));
    for (int i = 0; i < length(which_R); i++) {
	matches.subscripts[i] =
	    IIT_find(matches.nsubscripts + i, iit,
		     (char *)CHAR(STRING_ELT(which_R, i)));
    }
    return matches;    
}

static IITMatches _iit_get_exact_multiple(IIT_T iit,
					  SEXP chr_R, int *start, int *end,
					  int type)
{
    IITMatches matches = _new_IITMatches(iit, length(chr_R));
    for (int i = 0; i < length(chr_R); i++) {
	matches.subscripts[i] =
	    IIT_get_exact_multiple(matches.nsubscripts + i, iit,
				   (char *)CHAR(STRING_ELT(chr_R, i)),
				   start[i], end[i], type);
    }
    return matches;
}

static IITMatches _iit_get_typed(IIT_T iit,
				 SEXP chr_R, int *start, int *end,
				 int type)
{
    IITMatches matches = _new_IITMatches(iit, length(chr_R));
    for (int i = 0; i < length(chr_R); i++) {
	matches.subscripts[i] =
	    IIT_get_typed(matches.nsubscripts + i, iit,
			  (char *)CHAR(STRING_ELT(chr_R, i)),
			  start[i], end[i], type,
			  /*sortp*/false);
    }
    return matches;
}

static IITMatches _iit_get(IIT_T iit, SEXP chr_R, int *start, int *end)
{
    IITMatches matches = _new_IITMatches(iit, length(chr_R));
    for (int i = 0; i < length(chr_R); i++) {
	matches.subscripts[i] =
	    IIT_get(matches.nsubscripts + i, iit,
		    (char *)CHAR(STRING_ELT(chr_R, i)),
		    start[i], end[i],
		    /*sortp*/false);
    }
    return matches;
}

static IITMatches _iit_get_typed_signed(IIT_T iit,
					SEXP chr_R, int *start, int *end,
					int type, int *sign)
{
    IITMatches matches = _new_IITMatches(iit, length(chr_R));
    for (int i = 0; i < length(chr_R); i++) {
	matches.subscripts[i] =
	    IIT_get_typed_signed(matches.nsubscripts + i, iit,
				 (char *)CHAR(STRING_ELT(chr_R, i)),
				 start[i], end[i], type, sign[i],
				 /*sortp*/false);
    }
    return matches;
}

static IITMatches _iit_get_signed(IIT_T iit,
				  SEXP chr_R, int *start, int *end,
				  int *sign)
{
    IITMatches matches = _new_IITMatches(iit, length(chr_R));
    for (int i = 0; i < length(chr_R); i++) {
	matches.subscripts[i] =
	    IIT_get_signed(matches.nsubscripts + i, iit,
			   (char *)CHAR(STRING_ELT(chr_R, i)),
			   start[i], end[i], sign[i],
			   /*sortp*/false);
    }
    return matches;
}

static IITMatches _iit_get_for_coords(IIT_T iit, SEXP which_R, SEXP type_R,
				      SEXP ignore_strand_R, SEXP exact_R)
{
    SEXP chr_R = VECTOR_ELT(which_R, 0);
    int *start = INTEGER(VECTOR_ELT(which_R, 1));
    int *end = INTEGER(VECTOR_ELT(which_R, 2));
    int *sign = INTEGER(VECTOR_ELT(which_R, 3));
    int type = type_R == R_NilValue ? 0 :
	IIT_typeint(iit, (char *)CHAR(asChar(type_R)));
    bool ignore_strand = asLogical(ignore_strand_R);
    bool exact = asLogical(exact_R);
    IITMatches matches;
    
    if (exact) {
	/* sign filtering happens in R */
	matches = _iit_get_exact_multiple(iit, chr_R, start, end, type);
    } else if (ignore_strand) {
	if (type > 0) {
	    matches = _iit_get_typed(iit, chr_R, start, end, type);
	} else {
	    matches = _iit_get(iit, chr_R, start, end);
	}
    } else {
	if (type > 0) {
	    matches = _iit_get_typed_signed(iit, chr_R, start, end, type, sign);
	} else {
	    matches = _iit_get_signed(iit, chr_R, start, end, sign);
	}
    }
    
    return matches;
}

static IITMatches _iit_get_for_labels(IIT_T iit, SEXP which_R) {
    return _iit_find(iit, which_R);
}

enum { CHR, START, WIDTH, STRAND, ANNO, ANS_LENGTH };

static SEXP _convert_matches(IITMatches matches, bool ret_ranges, SEXP fields_R)
{
    SEXP ans, chr_R, start_R, width_R, strand_R, anno_R;
    int nfields = fields_R == R_NilValue ? 1 : length(fields_R);
    int *fields;
    int nmatches = 0;
    IIT_T iit = matches.iit;

    for (int m = 0; m < matches.nqueries; m++) {
	nmatches += matches.nsubscripts[m];
    }

    PROTECT(ans = allocVector(VECSXP, ANS_LENGTH));
    if (ret_ranges) {
	chr_R = allocVector(STRSXP, nmatches);
	SET_VECTOR_ELT(ans, CHR, chr_R);
	start_R = allocVector(INTSXP, nmatches);
	SET_VECTOR_ELT(ans, START, start_R);
	width_R = allocVector(INTSXP, nmatches);
	SET_VECTOR_ELT(ans, WIDTH, width_R);
	strand_R = allocVector(INTSXP, nmatches);
	SET_VECTOR_ELT(ans, STRAND, strand_R);
    }
    anno_R = allocVector(VECSXP, nfields);
    SET_VECTOR_ELT(ans, ANNO, anno_R);

    for (int f = 0; f < nfields; f++) {
	SET_VECTOR_ELT(anno_R, f, allocVector(STRSXP, nmatches));
    }

    if (fields_R != R_NilValue) {
	fields = (int *)R_alloc(sizeof(int), nfields);
	for (int f = 0; f < nfields; f++) {
	    fields[f] = IIT_fieldint(iit, (char *)STRING_ELT(fields_R, f));
	}
    }
    
    for (int i = 0; i < nmatches; i++) {
	if (ret_ranges) {
	    Interval_T interval = IIT_interval(iit, i);
	    SET_STRING_ELT(chr_R, i, mkChar(IIT_divstring_from_index(iit, i)));
	    INTEGER(start_R)[i] = Interval_low(interval);
	    INTEGER(width_R)[i] = Interval_length(interval);
	    INTEGER(strand_R)[i] = Interval_sign(interval);
	}
	if (fields_R == R_NilValue) {
	    char *restofheader;
	    bool allocp;
	    SET_STRING_ELT(anno_R, i,
			   mkChar(IIT_annotation(&restofheader, iit, i,
						 &allocp)));
	    if (allocp == true) {
		free(restofheader);
	    }
	} else {
	    for (int f = 0; f < nfields; f++) {
		SET_STRING_ELT(VECTOR_ELT(anno_R, f), i,
			       mkChar(IIT_fieldvalue(iit, i, fields[f])));
	    }
	}
    }

    UNPROTECT(1);
    return ans;
}

SEXP R_iit_read(SEXP iit_R, SEXP which_R, SEXP type_R, SEXP fields_R,
		SEXP ignore_strand_R, SEXP exact_R, SEXP ret_ranges_R)
{
    IITMatches matches;
    IIT_T iit = R_ExternalPtrAddr(iit_R);
    bool ret_ranges = asLogical(ret_ranges_R);

    if (TYPEOF(which_R) == VECSXP) {
	matches = _iit_get_for_coords(iit, which_R, type_R,
				      ignore_strand_R, exact_R);
    } else {
	matches = _iit_get_for_labels(iit, which_R);
    }

    return _convert_matches(matches, ret_ranges, fields_R);
}

SEXP R_iit_typeNames(SEXP iit_R) {
    IIT_T iit = R_ExternalPtrAddr(iit_R);
    SEXP ans;
    PROTECT(ans = allocVector(STRSXP, IIT_ntypes(iit)));
    for (int i = 0; i < IIT_ntypes(iit); i++) {
	SET_STRING_ELT(ans, i, mkChar(IIT_typestring(iit, i)));
    }
    UNPROTECT(1);
    return ans;
}

SEXP R_iit_fieldNames(SEXP iit_R) {
    IIT_T iit = R_ExternalPtrAddr(iit_R);
    SEXP ans;
    PROTECT(ans = allocVector(STRSXP, IIT_nfields(iit)));
    for (int i = 0; i < IIT_nfields(iit); i++) {
	SET_STRING_ELT(ans, i, mkChar(IIT_fieldstring(iit, i)));
    }
    UNPROTECT(1);
    return ans;
}

SEXP R_iit_length(SEXP iit_R) {
    IIT_T iit = R_ExternalPtrAddr(iit_R);
    return ScalarInteger(IIT_total_nintervals(iit));
}

/* 'IIT_output_direct' is not included in the libgstruct binary */

/*
  #include <gstruct/iit-write.h>

  SEXP R_iit_write(SEXP tally_iit_R, SEXP iitfile_R) {
  IIT_T tally_iit = (IIT_T) R_ExternalPtrAddr(tally_iit_R);
  char *iitfile = (char *) CHAR(asChar(iitfile_R));
  IIT_output_direct(iitfile, tally_iit, IIT_LATEST_VERSION);
  return iitfile_R;
  }
*/

