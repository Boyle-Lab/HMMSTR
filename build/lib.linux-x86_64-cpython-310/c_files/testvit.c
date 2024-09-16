/*
 *  File: testvit.c
 *
 *  process viterbi inputs, run viterbi, and process outputs
 *
 *  The HMM structure and some codes are borrowed and modified from Kanungo's
 *  original HMM program.
 *  Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999. http://www.kanungo.com/software/software.html.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <glib.h> //added
#include "nrutil.h"
#include "hmm.h"
static char rcsid[] = "$Id: testvit.c,v 1.3 1998/02/23 07:39:07 kanungo Exp kanungo $";

int main (int argc, char **argv)
{
	int 	t, T,T_new,repeat_len,prefix_idx;
	HMM  	hmm;
	int	*O;	/* observation sequence O[1..T] */
	int	*q;	/* state sequence q[1..T] */
	double **delta;
	int	**psi;
	int **traceback_dir; //added, for Viterbi modification
	double 	proba, logproba;
	FILE	*fp, *traceback_del, *traceback, *delta_out;
	int i; //added
	int n; //added
	int offset; //added
	char *S; //added, input sequence as chars
	GList *list_log = NULL;
	GList *likelihoods_ints = NULL;
	GList *likelihoods_dec = NULL;



	if (argc != 6) { //changed to 4 so we can include T and O separately

		printf("%d",argc);
		printf("Usage error \n");
		printf("Usage: testvit <model.hmm> <obs.seq length> <obs.seq> \n");
		exit (1);
	}
	//open HMM file
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[1]);
		exit (1);
	}
	//read HMM file according to ReadHMM
	ReadHMM(fp, &hmm);
	fclose(fp);

	T = atoi(argv[2]); //read T from input
	S = argv[3];
	repeat_len =  atoi(argv[4]); // for calculation of RDn case
	prefix_idx = atoi(argv[5]); // for calculation of RDn case

	O = ivector(1,T);

	i=1;
	while (sscanf(S, "%d%n", &n, &offset)==1){ //added this to read string of sequence
		O[i] = n;
		S += offset;

		i++;
	}




	q = ivector(1,T);

	traceback_dir = imatrix(1, T, 1, hmm.N); //added
	delta = dmatrix(1, T, 1, hmm.N);
	psi = imatrix(1, T, 1, hmm.N);

	ViterbiLog(&hmm, T, O, delta, psi, traceback_dir,&list_log,&likelihoods_ints,&likelihoods_dec,&T_new, q, &logproba, repeat_len, prefix_idx);

	fprintf(stdout, "Viterbi  MLE log prob = %E\n", logproba);
	fprintf(stdout, "Optimal state sequence:\n");
	//PrintSequence(stdout, T, q);
	PrintPath(stdout, T_new, list_log);
	PrintLikelihoods(stdout, T_new, likelihoods_ints);
	PrintLikelihoods(stdout, T_new, likelihoods_dec);

	printf("------------------------------------\n");
	printf("The two log probabilites and optimal state sequences\n");
	printf("should identical (within numerical precision). \n");

	free_ivector(q, 1, T);
	free_ivector(O, 1, T);
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_imatrix(traceback_dir,1,T,1,hmm.N); //added
	g_list_free(list_log);
	g_list_free(likelihoods_ints);
	g_list_free(likelihoods_dec);
	free_dmatrix(delta, 1, T, 1, hmm.N);
	FreeHMM(&hmm);
}
