/* $Id: dropgsw2.c $ */

/* copyright (c) 1996, 2014 by William R. Pearson and The Rector &
   Visitors of the University of Virginia */

/* Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing,
   software distributed under this License is distributed on an "AS
   IS" BASIS, WITHOUT WRRANTIES OR CONDITIONS OF ANY KIND, either
   express or implied.  See the License for the specific language
   governing permissions and limitations under the License. 
*/

/* 17-Aug-2006 - removed globals *sapp/last - alignment should be thread safe */

/* 12-Oct-2005 - converted to use a_res and aln for alignment coordinates */

/* 4-Nov-2004 - Diagonal Altivec Smith-Waterman included */

/* 14-May-2003 - modified to return alignment start at 0, rather than
   1, for begin:end alignments

   25-Feb-2003 - modified to support Altivec parallel Smith-Waterman

   22-Sep-2003 - removed Altivec support at request of Sencel lawyers
*/

/* this code uses an implementation of the Smith-Waterman algorithm
   designed by Phil Green, U. of Washington, that is 1.5 - 2X faster
   than my Miller and Myers implementation. */

/* the shortcuts used in this program prevent it from calculating scores
   that are less than the gap penalty for the first residue in a gap. As
   a result this code cannot be used with very large gap penalties, or
   with very short sequences, and probably should not be used with prss3.
*/

/* version 3.2 fixes a subtle bug that was encountered while running
   do_walign() interspersed with do_work().  This happens only with -m
   9 and pvcomplib.  The fix was to more explicitly zero-out ss[] at
   the beginning of do_work.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"

static char *verstr="7.2 Nov 2010";

#include "dropgsw2.h"
#include "smith_waterman_gpu.h"

#define DROP_INTERN
#include "drop_func.h"

struct swstr {int H, E;};

extern void init_karlin(const unsigned char *aa0, int n0, struct pstruct *ppst,
			double *aa0_f, double **kp);
extern int do_karlin(const unsigned char *aa1, int n1,
		     int **pam2, const struct pstruct *ppst,
		     double *aa0_f, double *kar_p, double *lambda, double *H);

extern int sw_walign (int **pam2p, int n0,
		      const unsigned char *aa1, int n1,
		      int q, int r,
		      struct swstr *ss,
		      struct a_res_str *a_res
		      );

extern struct a_res_str *
nsw_malign (int ***pam2p, int pam_ix, int n0,
	    const unsigned char *aa1, int n1,
	    int score_thresh, int max_res,
	    int gdelval, int ggapval, 
	    struct swstr *ss, 
	    struct a_res_str *cur_ares,
	    int (*fn_walign)
	    (
	     int **pam2p, int n0,
	     const unsigned char *aa1, int n1,
	     int q, int r,
	     struct swstr *ss,
	     struct a_res_str *a_res
	     ),
	    int do_rep
	    );

void
SIM(const unsigned char *A, /* seq1 indexed A[1..M] */
    const unsigned char *B, /* seq2 indexed B[1..N] */
    int M, int N,		/* len seq1, seq2 */
    struct pstruct *ppst,	/* parameters */
    int nseq,			/* nseq - number of different sequences */
    int mini_score,		/* cut-off score */
    int max_count,		/* number of alignments */
    struct a_res_str *a_res);	/* alignment result structure */

static int
FLOCAL_ALIGN(const unsigned char *aa0, const unsigned char *aa1,
	     int n0, int n1, int low, int up,
	     int **W, int GG,int HH, int MW,
	     struct f_struct *f_str);

extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);

int same_seq(const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1);

static int
prof_score(const unsigned char *aa1, int n0, int *pwaa_s);

/* initialize for Smith-Waterman optimal score */

void
init_work (unsigned char *aa0, int n0,
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
  int ip;
  int *pwaa_s, *pwaa_a;
  int e, f, i, j;
  struct f_struct *f_str;
  int **pam2p;
  struct swstr *ss;
  int nsq;

  if (ppst->ext_sq_set) {
    nsq = ppst->nsqx; ip = 1;
  }
  else {
    /* set for lower-case for memory mapped DBs with lower case encoding */
    nsq = ppst->nsqx; ip = 0;
  }

   /* allocate space for function globals */
   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   if((ppst->zsflag%10) == 6) {
     f_str->kar_p = NULL;
     init_karlin(aa0, n0, ppst, &f_str->aa0_f[0], &f_str->kar_p);
   }
  
   /* allocate space for the scoring arrays */
   if ((ss = (struct swstr *) calloc (n0+2, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, "cannot allocate ss array %3d\n", n0);
     exit (1);
   }
   ss++;

   ss[n0].H = -1;	/* this is used as a sentinel - normally H >= 0 */
   ss[n0].E = 1;
   f_str->ss = ss;

   /* initialize variable (-S) pam matrix */
   if ((f_str->waa_s= (int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate waa_s array %3d\n",nsq*n0);
     exit(1);
   }

   /* initialize pam2p[1] pointers */
   if ((f_str->pam2p[1]= (int **)calloc((n0+1),sizeof(int *))) == NULL) {
     fprintf(stderr,"cannot allocate pam2p[1] array %3d\n",n0);
     exit(1);
   }

   pam2p = f_str->pam2p[1];
   if ((pam2p[0]=(int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate pam2p[1][] array %3d\n",nsq*n0);
     exit(1);
   }

   for (i=1; i<n0; i++) {
     pam2p[i]= pam2p[0] + (i*(nsq+1));
   }

   /* initialize universal (alignment) matrix */
   if ((f_str->waa_a= (int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate waa_a struct %3d\n",nsq*n0);
     exit(1);
   }
   
   /* initialize pam2p[0] pointers */
   if ((f_str->pam2p[0]= (int **)calloc((n0+1),sizeof(int *))) == NULL) {
     fprintf(stderr,"cannot allocate pam2p[1] array %3d\n",n0);
     exit(1);
   }

   pam2p = f_str->pam2p[0];
   if ((pam2p[0]=(int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate pam2p[1][] array %3d\n",nsq*n0);
     exit(1);
   }

   for (i=1; i<n0; i++) {
     pam2p[i]= pam2p[0] + (i*(nsq+1));
   }

   /* 
      pwaa effectively has a sequence profile --
       pwaa[0..n0-1] has pam score for residue 0 (-BIGNUM)
       pwaa[n0..2n0-1] has pam scores for amino acid 1 (A)
       pwaa[2n0..3n0-1] has pam scores for amino acid 2 (R), ...

       thus: pwaa = f_str->waa_s + (*aa1p++)*n0; sets up pwaa so that
       *pwaa++ rapidly moves though the scores of the aa1p[] position
       without further indexing

       For a real sequence profile, pwaa[0..n0-1] vs ['A'] could have
       a different score in each position.
   */

   pwaa_s = f_str->waa_s;
   pwaa_a = f_str->waa_a;
   if (ppst->pam_pssm) {
     for (e = 0; e <=nsq; e++)	{	/* for each residue in the alphabet */
       for (f = 0; f < n0; f++) {	/* for each position in aa0 */
	 *pwaa_s++ = f_str->pam2p[ip][f][e] = ppst->pam2p[ip][f][e];
	 *pwaa_a++ = f_str->pam2p[0][f][e]  = ppst->pam2p[0][f][e];
       }
     }
   }
   else {	/* initialize scanning matrix */
     for (e = 0; e <=nsq; e++)	/* for each residue in the alphabet */
       for (f = 0; f < n0; f++)	{	/* for each position in aa0 */
	 *pwaa_s++ = f_str->pam2p[ip][f][e]= ppst->pam2[ip][aa0[f]][e];
	 *pwaa_a++ = f_str->pam2p[0][f][e] = ppst->pam2[0][aa0[f]][e];
       }
   }

    /* minimum allocation for alignment */
    f_str->max_res =  max(3*n0/2,MIN_RES);

    *f_arg = f_str;
}

void close_work (const unsigned char *aa0, int n0,
		 struct pstruct *ppst,
		 struct f_struct **f_arg)
{
  struct f_struct *f_str;

  f_str = *f_arg;

  if (f_str != NULL) {
    if (f_str->kar_p !=NULL) free(f_str->kar_p);
    f_str->ss--;
    free(f_str->ss);
    free(f_str->waa_a);
    free(f_str->pam2p[0][0]);
    free(f_str->pam2p[0]);
    free(f_str->waa_s);
    free(f_str->pam2p[1][0]);
    free(f_str->pam2p[1]);

    free(f_str);
    *f_arg = NULL;
  }
}


/* pstring1 is a message to the manager, currently 512 */
/*void get_param(struct pstruct *ppst,char *pstring1)*/
void    get_param (const struct pstruct *ppst,
		   char **pstring1, char *pstring2,
		   struct score_count_s *s_cnt_info)
{
  char pg_str[120];
  char psi_str[120];

strncpy(pg_str,"Smith-Waterman (PGopt)",sizeof(pg_str));

  if (ppst->pam_pssm) { strncpy(psi_str,"-PSI",sizeof(psi_str));}
  else { psi_str[0]='\0';}

  sprintf (pstring1[0], "%s (%s)", pg_str, verstr);
  sprintf (pstring1[1], 
#ifdef OLD_FASTA_GAP
	   "%s matrix%s (%d:%d)%s, gap-penalty: %d/%d",
#else
	   "%s matrix%s (%d:%d)%s, open/ext: %d/%d",
#endif
	   ppst->pam_name, psi_str, ppst->pam_h,ppst->pam_l, 
	   (ppst->ext_sq_set)?"xS":"\0", ppst->gdelval, ppst->ggapval);

   if (pstring2 != NULL) {
#ifdef OLD_FASTA_GAP
     sprintf(pstring2,"; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s%s (%d:%d)%s\n; pg_gap-pen: %d %d\n",
#else
     sprintf(pstring2,"; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s%s (%d:%d)%s\n; pg_open-ext: %d %d\n",
#endif
	     pg_str,verstr,ppst->pam_name,psi_str,ppst->pam_h,ppst->pam_l, 
	     (ppst->ext_sq_set)?"xS":"\0",ppst->gdelval,ppst->ggapval);
   }
}

void do_work (const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      int frame,
	      const struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, int shuff_flg, struct rstruct *rst,
	      struct score_count_s *sc_info)
{
  int     score;
  double lambda, H;
  int i;
  
#ifdef LALIGN
  if (same_seq(aa0,n0,aa1,n1)) {
    rst->score[0] = prof_score(aa1, n0, f_str->waa_s);
    return;
  }
#endif

  rst->alg_info = 0;
  rst->valid_stat = 1;
  sc_info->s_cnt[0]++;
  sc_info->tot_scores++;

  score = FLOCAL_ALIGN(aa0,aa1,n0,n1,0,0,
                       NULL,
#ifdef OLD_FASTA_GAP
                       -(ppst->gdelval - ppst->ggapval),
#else
                       -ppst->gdelval,
#endif
                       -ppst->ggapval,0,f_str);



  rst->score[0] = score;
  rst->score[1] = rst->score[2] = 0;

  if(((ppst->zsflag % 10) == 6) &&
     (do_karlin(aa1, n1, ppst->pam2[0], ppst,f_str->aa0_f, 
		f_str->kar_p, &lambda, &H)>0)) {
    rst->comp = 1.0/lambda;
    rst->H = H;

  }
  else {rst->comp = rst->H = -1.0;}

}


/*  Carrick notes:
    Smith-Waterman algorithm, written by Dr. Pearson.
    
    Plans:
      - Rename variables (for clarity).
      - parallelize.
*/
static int
FLOCAL_ALIGN(const unsigned char *aa0, const unsigned char *aa1,
	     int n0, int n1, int low, int up,
	     int **W, int GG,int HH, int MW,
	     struct f_struct *f_str) {
  return smith_waterman( aa0, aa1, n0, n1, low, up, W, GG, HH, MW, f_str );
}		/* here we should be all done */

void do_opt (const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1,
	     int frame,
	     struct pstruct *ppst, struct f_struct *f_str,
	     struct rstruct *rst)
{
}

struct a_res_str *
do_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, int repeat_thresh,
	   struct pstruct *ppst, 
	   struct f_struct *f_str, 
	   int *have_ares)
{
  int a_res_index;
  struct a_res_str *a_res, *tmp_a_res;


  *have_ares = 0x3;	/* set 0x2 bit to indicate local copy */

  if ((a_res = (struct a_res_str *)calloc(1, sizeof(struct a_res_str)))==NULL) {
    fprintf(stderr," [do_walign] Cannot allocate a_res");
    return NULL;
  }

#ifndef LALIGN
  a_res = nsw_malign(f_str->pam2p, (ppst->ext_sq_set ? 1 : 0), n0, aa1, n1,
		     repeat_thresh, f_str->max_res,
		     -ppst->gdelval, -ppst->ggapval,
		     f_str->ss, a_res,
		     &sw_walign, ppst->do_rep);

#else	/* LALIGN */
  if (!ppst->show_ident && same_seq(aa0, n0, aa1, n1)) ppst->nseq = 1;
  else ppst->nseq = 2;

  SIM(aa0-1, aa1-1, n0, n1, ppst, ppst->nseq, repeat_thresh, ppst->max_repeat, a_res);
#endif

  /* set a_res->index for alignments */

  a_res_index = 0;
  for (tmp_a_res=a_res; tmp_a_res; tmp_a_res = tmp_a_res->next) {
    tmp_a_res->index = a_res_index++;
  }

  return a_res;
}

/*
#define XTERNAL
#include "upam.h"

void
print_seq_prof(unsigned char *A, int M,
	       unsigned char *B, int N,
	       int **w, int iw, int dir) {
  char c_max;

  int i_max, j_max, i,j;

  char *c_dir="LRlr";

  for (i=1; i<=min(60,M); i++) {
    fprintf(stderr,"%c",aa[A[i]]);
  }
  fprintf(stderr, - %d\n,M);

  for (i=0; i<min(60,M); i++) {
    i_max = -1;
    for (j=1; j<21; j++) {
      if (w[iw+i][j]> i_max) {

	i_max = w[iw+i][j]; 
	j_max = j;
      }
    }
    fprintf(stderr,"%c",aa[j_max]);
  }
  fputc(':',stderr);

  for (i=1; i<=min(60,N); i++) {
    fprintf(stderr,"%c",aa[B[i]]);
  }

  fprintf(stderr," -%c: %d,%d\n",c_dir[dir],M,N);
}
*/

void
pre_cons(const unsigned char *aa1, int n1, int frame, struct f_struct *f_str) {

#ifdef TFAST
  f_str->n10 = aatran(aa1,f_str->aa1x,n1,frame);

#endif

}

/* aln_func_vals - set up aln.qlfact, qlrev, llfact, llmult, frame, llrev */
/* call from calcons, calc_id, calc_code */
void 
aln_func_vals(int frame, struct a_struct *aln) {

  aln->llfact = aln->llmult = aln->qlfact = 1;
  aln->llrev = 0;
  if (frame > 0) aln->qlrev = 1;
  else aln->qlrev = 0;
  aln->frame = 0;
}


/* calculate the 100% identical score */
int
prof_score(const unsigned char *aa1p, int n0, int *pwaa_s)
{
  int sum=0;

  while (*aa1p) {
    sum += pwaa_s[(*aa1p++)*n0];
    pwaa_s++;
  }
  return sum;
}


int same_seq(const unsigned char *aa0, int n0,
	     const unsigned char *aa1, int n1)
{
  const unsigned char *ap0, *ap1;
  int cnt=0;

  if (n0 != n1) return 0;

  ap0 = aa0;
  ap1 = aa1;
  
  while ( *ap0 && *ap0++ == *ap1++ ) {cnt++;}
  if (cnt != n0) return 0;
  return 1;
}

