#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"
#include "dropgsw2.h"


// inline void gpu_handle_error( cudaError_t err, const char* file, int line, int abort = 1 )
// {
// 	if (err != cudaSuccess)
// 	{
// 		fprintf (stderr, "gpu error %s, %s, %d\n", cudaGetErrorString (err), file, line);
// 		if (abort)
// 			exit (EXIT_FAILURE);
// 	}
// }
// #define gpu_err_chk(e) {gpu_handle_error( e, __FILE__, __LINE__ );}


// kernel function
// __global__
// void smith_waterman_kernel (const unsigned char *aal)
// {
// 	int i = blockIdx.x * blockDim.x + threadIdx.x;
// }

extern "C" {
  #include "smith_waterman_gpu.h"
}

struct swstr {int H, E;};

extern "C"
int smith_waterman( const unsigned char *aa0,
                    const unsigned char *aa1,
                    int n0,
                    int n1,
                    int low,
                    int up,
                    int **W,
                    int GG,
                    int HH,
                    int MW,
                    struct f_struct *f_str)   //defined in dropgsw.h
{
  // cudaError_t err;
  // const unsigned char *d_aa1;
  // const unsigned char *d_aa1p;
  const unsigned char *aa1p;
  // struct f_struct *d_f_str;
  // register int *d_pwaa;
  register int *pwaa;
  // register struct swstr *d_ssj;
  register struct swstr *ssj;
  // struct swstr *d_ss;
  struct swstr *ss;
  register int h, e, f, p;
  int temp, score;
  int gap_ext, n_gap_init;
  // int *d_score_pt;
  int *score_pt;

  ss = f_str->ss;
  ss[n0].H = -1;
  ss[n0].E = 1;

  n_gap_init = GG + HH;
  gap_ext = -HH;	/* GG, HH are both positive,
                gap_ext penalty should be negative */

  score = 0;
  for (h=0; h<n0; h++) {	  /* initialize 0th row */
    ss[h].H = ss[h].E = 0;

  }
  
  aa1p=aa1;

  // err = cudaMalloc ((void**) &d_aa1, n1 * sizeof(char));
  // gpu_err_chk(err);
  // err = cudaMalloc ((void**) &d_aa1p, n1 * sizeof(char));
  // gpu_err_chk(err);
  // err = cudaMalloc ((void**) &d_f_str, sizeof(f_struct));
  // gpu_err_chk(err);
  // err = cudaMalloc ((void**) &d_pwaa, sizeof(int));
  // gpu_err_chk(err);
  // err = cudaMalloc ((void**) &d_ssj, sizeof(swstr));
  // gpu_err_chk(err);
  // err = cudaMalloc ((void**) &d_ss, sizeof(swstr));
  // gpu_err_chk(err);
  // err = cudaMalloc ((void**) &d_score_pt, sizeof(int));
  // gpu_err_chk(err);

  // err = cudaMemcpy (d_aa1, h_aa1, n1 * sizeof(char), cudaMemcpyHostToDevice);
  // gpu_err_chk(err);
  // err = cudaMemcpy (d_aa1p, aa1p, n1 * sizeof(char), cudaMemcpyHostToDevice);
  // gpu_err_chk(err);
  // err = cudaMemcpy (d_f_str, f_str, sizeof(f_struct), cudaMemcpyHostToDevice);
  // gpu_err_chk(err);
  // err = cudaMemcpy (d_pwaa, pwaa, sizeof(int), cudaMemcpyHostToDevice);
  // gpu_err_chk(err);
  // err = cudaMemcpy (d_ssj, ssj, sizeof(swstr), cudaMemcpyHostToDevice);
  // gpu_err_chk(err);
  // err = cudaMemcpy (d_ss, ss, sizeof(swstr), cudaMemcpyHostToDevice);
  // gpu_err_chk(err);
  // err = cudaMemcpy (d_score_pt, score_pt, sizeof(int), cudaMemcpyHostToDevice);
  // gpu_err_chk(err);
  
  while (*aa1p) {		/* relies on d_aa1[n1]==0 for EOS flag */
    /* waa_s has the offsets for each residue in d_aa0 into pam2
  */
    /* waa_s has complexity (-S) dependent scores */
    pwaa = f_str->waa_s + (*aa1p++)*n0;
    ssj = ss;

    e = f = h = p = 0;
  zero_f:	/* in this section left-gap f==0, and is never examined */

    while (1) {	/* build until h > n_gap_init (f < 0 until h > n_gap_init) */
              /* bump through the pam[][]'s for each of the d_aa1[] matches to
              d_aa0[], because of the way *d_pwaa is set up */


      h = p + *pwaa++;		/* increment diag value */
      p = ssj->H;		/* get next diag value */
      if ((e = ssj->E) > 0 ) {	/* >0 from up-gap */
    if (p == -1) goto next_row;	/* done, -1=d_ss[n0].H sentinel */
    if (h < e) h = e;	/* up-gap better than diag */
    else 
      if (h > n_gap_init) {	/* we won't starting a new up-gap */
        e += gap_ext;	/* but we might be extending one */
        goto transition;	/* good h > n_gap_diag; scan f */
      }
    e += gap_ext;		/* up-gap decreased */
    ssj->E =  (e > 0) ?  e : 0;	/* set to 0 if < 0 */
    ssj++->H = h;		/* diag match updated */
      }
      else {			/* up-gap (->E) is 0 */

    if ( h > 0) {		/* diag > 0 */
      if (h > n_gap_init) {	/* we won't be starting a new up-gap */
        e = 0;		/* and we won't be extending one */
        goto transition;	/* good h > n_gap_diag; scan f */
      }
      ssj++->H = h;		/* update diag */
    }
    else ssj++->H = 0;	/* update diag to 0 */
      }
    }

    /* here h > n_gap_init and h > e, => the next f will be > 0 */
  transition:

#ifdef DEBUG
    if ( h > 10000) 
      fprintf(stderr,"h: %d d_ssj: %d\n",h, (int)(ssj-ss));
#endif
    if ( score < h ) score = h;	/* save best score, only when h > n_gap_init */

    temp = h - n_gap_init;	/* best score for starting a new gap */
    if ( f < temp ) f = temp;	/* start a left-gap? */
    if ( e < temp ) e = temp;	/* start an up-gap? */
    ssj->E = ( e > 0 ) ? e : 0;	/* update up-gap */

    ssj++->H = h;		/* update diag */
    e = 0;

    do {			/* stay here until f <= 0 */
      h = p + *pwaa++;		/* diag + match/mismatch */
      p = ssj->H;		/* save next (right) diag */

      if ( h < f ) h = f;	/* update diag using left gap */
      f += gap_ext;		/* update next left-gap */

      if ((e = ssj->E) > 0) {	/* good up gap */
    if (p == -1) goto next_row;	/* at the end of the row */
    if ( h < e ) h = e;	/* update diag using up-gap */
    else
      if ( h > n_gap_init ) {
        e += gap_ext;	/* update up gap */
        goto transition;	/* good diag > n_gap_init, restart */
      }
    e += gap_ext;		/* update up-gap */
    ssj->E = (e > 0) ? e : 0;	/* e must be >= 0 */
    ssj++->H = h;		/* update diag */
      }
      else {			/* up-gap <= 0 */
    if ( h > n_gap_init ) {
      e = 0;
      goto transition;	/* good diag > n_gap_init; restart */
    }
    ssj++->H = h;		/* update diag */
      }
    } while ( f > 0 );		/* while left gap f > 0  */
    goto zero_f;		/* otherwise, go to f==0 section */
  next_row:
    ;
  }		/* end while(*aap1) {} */

  return score;  
}



