/*
 *  File: viterbi.c
 *
 *  viterbi step and hidden states labeling.
 *
 *  The HMM structure and some codes are borrowed and modified from Kanungo's
 *  original HMM program.
 *  Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999. http://www.kanungo.com/software/software.html.
 *
 */

#include <math.h>
#include "hmm.h"
#include "nrutil.h"
#include <glib.h>
static char rcsid[] = "$Id: viterbi.c,v 1.1 1999/05/06 05:25:37 kanungo Exp kanungo $";

#define VITHUGE  100000000000.0

void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi, int **traceback_dir,GList **list,
	int *T_new, int *q, double *pprob)
{
	int 	i, j;	/* state indices */
	int  	t;	/* time index */
	int blank_index = 5;
	int	maxvalind;
	double	maxval, val,val2,deletion,regular; //val2,deletion,and regular added for deletion cases

	/* 1. Initialization  */

	for (i = 1; i <= phmm->N; i++) {
		delta[1][i] = phmm->pi[i] * (phmm->B[i][O[1]]);
		psi[1][i] = 0;
		traceback_dir[1][i]=0; //added
	}

	/* 2. Recursion */

	for (t = 2; t <= T; t++) {
		for (j = 1; j <= phmm->N; j++) {
			maxval = 0.0;
			maxvalind = 1;
			val2 = 0.0;
			for (i = 1; i <= phmm->N; i++) {
				val = delta[t-1][i]*(phmm->A[i][j]);
				//added deletion case
				if ((i-2)%3 == 0 && (j-5)==i){

					val2 = delta[t][i]*(phmm->A[i][j]);
					if (val2 > val){
						val = val2;

					}
				}
				//end of added
				if (val > maxval) {
					maxval = val;
					maxvalind = i;
				}
				//added traceback code, we only want traceback to deletion nodes to happen here
				if (maxval == val2 && val2 !=0) {//track if we chose a deletion as our max value

					traceback_dir[t][j] = 1;
				}
        else{
					traceback_dir[t][j] = 0;
				}
				//end of added traceback code
			}
			deletion = maxval*(phmm->B[j][blank_index]);
			regular = maxval*(phmm->B[j][O[t]]);
			if(deletion > regular){
					delta[t][j] = deletion;
			}
			else{
				delta[t][j] = regular;
			}
			psi[t][j] = maxvalind;
		}
	}

	/* 3. Termination */

	*pprob = 0.0;
	q[T] = 1;
	for (i = 1; i <= phmm->N; i++) {
                if (delta[T][i] > *pprob) {
			*pprob = delta[T][i];
			q[T] = i;
		}
	}

	/* 4. Path (state sequence) backtracking */

/* 4. Modified to use a linked list so we can include additional states in path */
	*list = g_list_append(*list, GINT_TO_POINTER (q[T])); //head is equal to last state
	t=T;
	*T_new = 1; // I don't think this is necessary cause we're accessing the previous node in the list
	while (t >= 2){
		if (traceback_dir[t][GPOINTER_TO_INT(g_list_last(*list)->data)] == 1){
			*list = g_list_append(*list,GINT_TO_POINTER(psi[t][GPOINTER_TO_INT(g_list_last(*list)->data)]));
			*T_new+=1;
			continue;
		}
		else{
			*list = g_list_append(*list,GINT_TO_POINTER(psi[t][GPOINTER_TO_INT(g_list_last(*list)->data)]));
			t = t-1;
			*T_new+=1;
		}
	}

}
void ViterbiLog(HMM *phmm, int T, int *O, double **delta, int **psi,int **traceback_dir,GList **list,GList **likelihoods_ints,GList **likelihoods_dec,int *T_new,
        int *q, double *pprob, int repeat_len, int prefix_idx)
{
        int     i, j,l,p,k;   /* state indices */
        int     t;      /* time index */
				int blank_index=5; //index of "blank"
        int     maxvalind;
        double  maxval, val,val2,deletion,regular; //val2,deletion,and regular added to account for deltions
				double  **biot;
		int	RM1, RD1; // booleans for RDn case
		int continue_update, stop_early, next_del, curr_del; // for RDn updating loop

	/* 0. Preprocessing */
  //this section pre-converts everything to log, I don't think I need to change anything
	for (i = 1; i <= phmm->N; i++)
		phmm->pi[i] = log(phmm->pi[i]);
	for (i = 1; i <= phmm->N; i++)
		for (j = 1; j <= phmm->N; j++) {
			phmm->A[i][j] = log(phmm->A[i][j]);
		}

	biot = dmatrix(1, phmm->N, 1, T);
	for (i = 1; i <= phmm->N; i++)
		for (t = 1; t <= T; t++) {
			biot[i][t] = log(phmm->B[i][O[t]]);
		}
        /* 1. Initialization  */

        for (i = 1; i <= phmm->N; i++) {
                delta[1][i] = phmm->pi[i] + biot[i][1];
                psi[1][i] = 0;
								traceback_dir[1][i] = 0;
        }

        /* 2. Recursion */

        for (t = 2; t <= T; t++) {
                for (j = 1; j <= phmm->N; j++) {
                        maxval = -VITHUGE;
                        maxvalind = 1;
												val2 = -VITHUGE;
                        for (i = 1; i <= phmm->N; i++) {
								
								if ((i-2)%3 == 0 && (i !=phmm->N -1)){ //i is a deletion index, check the same level. this will not account for RDn yet
									if((j-5)==i || (j-3)==i || (j-1)==i){
										val2 = delta[t][i] + (phmm->A[i][j]);
										val=val2;
									}
								}
								else{ //i is not a deletion
									val = delta[t-1][i] + (phmm->A[i][j]);
								}
                                if (val > maxval) {
                                        maxval = val;
                                        maxvalind = i;
                                }
																if (maxval == val2 && val2 !=-VITHUGE) {//track if we chose a deletion as our max value
																	traceback_dir[t][j] = 1;
																}
		        									else{
																traceback_dir[t][j] = 0;
															}
															//end of added traceback code
														}
										deletion = maxval + log(phmm->B[j][blank_index]);
										regular = maxval + biot[j][t];
										if(deletion > regular){
											delta[t][j] = deletion;
										}
										else{
											delta[t][j] = regular;
										}
                    psi[t][j] = maxvalind;

                }
				//row t is filled, check for RDn being a better path retroactively
				continue_update = 1;
				stop_early = 0;
				//find the earliest deletion in RDn's path so we can stop before then, we don't want it to come all the way around and it shouldn't
				curr_del = prefix_idx+3*(repeat_len-1)+1;
				next_del = prefix_idx+3*(repeat_len-1)+1-3;
				while(psi[t][curr_del] == next_del){
					//check if we have traced back to the prefix
					if (next_del <= prefix_idx){
						break;
					}
					stop_early++; //increment to indicated we moved back in the path
					curr_del = next_del;
					next_del = curr_del - 3; // check next deletion
				}
				//check RM1, if RM1 is more likely from RDn than what it currently is, change the value (subtract out the emission when comparing the transition). No downstream affects because no value on level t is dependent on it
				if((delta[t][prefix_idx+3*(repeat_len-1)+1] + (phmm->A[prefix_idx+3*(repeat_len-1)+1][10])) > (delta[t][10] - biot[10][t])){
					//transitioning from RDn is better than previous recorded transition for RM1
					delta[t][10] = delta[t][prefix_idx+3*(repeat_len-1)+1] + (phmm->A[prefix_idx+3*(repeat_len-1)+1][10]) + biot[10][t]; //update delta
					psi[t][10] = prefix_idx+3*(repeat_len-1)+1; //update traceback to RDn
					traceback_dir[t][10] = 1; //update to indicate coming from a deletion
				}
				//check RD1 --> more complicated. If coming from RDn is more probabe than current path to RD1, change the value and all affected values UP TO RDn
				// We cannot allow this to also change RDn because it will cause an infinite loop.--> check if new value creates maximum path to: RM2, RI1, and/or RD2-->RD(n-1) (this will also require checking all M and I states from these deletion states)
				// if at anypoint (1) RD3 is chosen as the most likely previous state and/or (2) RD1 or following deletions are chosen as the NEW path, update both tracbacks
				// this is essentially going to be what we did above but restricted to only transitions that are done on row t 

				//if this is a homopolymer or a dinucleotide repeat, don't do any of this
				if (repeat_len == 1 || repeat_len == 2){
					continue;
				}
				//check RD1
				if ((delta[t][prefix_idx+3*(repeat_len-1)+1] + phmm->A[prefix_idx+3*(repeat_len-1)+1][8]) > (delta[t][8] - log(phmm->B[8][blank_index]))){
					//update RD1
					delta[t][8] = delta[t][prefix_idx+3*(repeat_len-1)+1] + phmm->A[prefix_idx+3*(repeat_len-1)+1][8] + log(phmm->B[8][blank_index]);
					psi[t][8] = prefix_idx+3*(repeat_len-1)+1;
					traceback_dir[t][8] = 1;

					// update related paths
					l = prefix_idx + 1; //current deletion index
					p = prefix_idx+2; //current insertion index
					k = prefix_idx + 6; //current next match index
					while (continue_update == 1){
						//check if we should be checking the next deletion's paths at all
						if (l == prefix_idx+3*(repeat_len-1)+1){
							break;
						}
						//check if next match came from current deletion
						if ((delta[t][l] + phmm->A[l][k]) > (delta[t][k] - biot[k][t])){
							//current deletion to next match is optimal, update path
							delta[t][k] = delta[t][l] + phmm->A[l][k] + biot[k][t];
							psi[t][k] = l;
							traceback_dir[t][k] = 1;
						}
						//check if current deletion to insertion is optimal
						if ((delta[t][l] + phmm->A[l][p]) > (delta[t][p] - biot[p][t])){
							delta[t][p] = delta[t][l] + phmm->A[l][p] + biot[p][t];
							psi[t][p] = l;
							traceback_dir[t][p] = 1;
						}
						//check if we want to update the next deletion state, if its the deletion before RDn or RD(n-1) is in RDns path, break
						if(l+stop_early*3+3 == (prefix_idx+3*(repeat_len-1)+1) && stop_early != 0 || l+3 == (prefix_idx+3*(repeat_len-1)+1)){ // this is the expanded version
							break;
						}
						//check if current deletion is optimal for next deletion
						if((delta[t][l] + phmm->A[l][l+3]) > (delta[t][l+3] - log(phmm->B[l+3][blank_index]))){
							delta[t][l+3] = delta[t][l] + phmm->A[l][l+3] + log(phmm->B[l+3][blank_index]);
							psi[t][l+3] = l;
							traceback_dir[t][l+3] = 1;

							//update indeces for next loop
							l = l+3;
							p=p+3;
							k=k+3;
						}
						else{
							//if the deletion path isn't optimal, no paths to update
							continue_update = 0;
						}
					
					}
				}
        }

        /* 3. Termination */

        *pprob = -VITHUGE;
        q[T] = 1;
        for (i = 1; i <= phmm->N; i++) {
                if (delta[T][i] > *pprob) {
                        *pprob = delta[T][i];
                        q[T] = i;
                }
        }


	/* 4. Path (state sequence) backtracking */

	*list = g_list_append(*list, GINT_TO_POINTER (q[T])); //head is equal to last state
	//adding a parallel list to record likelihoods, if this is inefficient then get rid of this
	*likelihoods_ints = g_list_append(*likelihoods_ints, GINT_TO_POINTER (*pprob));
	*likelihoods_dec = g_list_append(*likelihoods_dec, GINT_TO_POINTER (1000*(*pprob)-1000*((int)*pprob)));
	t=T;
	*T_new = 1; 
	while (t >= 2){
		if (traceback_dir[t][GPOINTER_TO_INT(g_list_last(*list)->data)] == 1){
			*list = g_list_append(*list,GINT_TO_POINTER(psi[t][GPOINTER_TO_INT(g_list_last(*list)->data)]));
			*likelihoods_ints = g_list_append(*likelihoods_ints,GINT_TO_POINTER(delta[t-1][GPOINTER_TO_INT(g_list_last(*list)->data)]));
			*likelihoods_dec = g_list_append(*likelihoods_dec,GINT_TO_POINTER(1000*delta[t-1][GPOINTER_TO_INT(g_list_last(*list)->data)]-1000*(int)delta[t-1][GPOINTER_TO_INT(g_list_last(*list)->data)]));
			*T_new+=1;
			continue;
		}
		else{
			*list = g_list_append(*list,GINT_TO_POINTER(psi[t][GPOINTER_TO_INT(g_list_last(*list)->data)]));
			*likelihoods_ints = g_list_append(*likelihoods_ints,GINT_TO_POINTER(delta[t-1][GPOINTER_TO_INT(g_list_last(*list)->data)]));
			*likelihoods_dec= g_list_append(*likelihoods_dec,GINT_TO_POINTER(1000*delta[t-1][GPOINTER_TO_INT(g_list_last(*list)->data)]-1000*(int)delta[t-1][GPOINTER_TO_INT(g_list_last(*list)->data)])); //I don't *think* that I am suppose to be indexing delta differently here, this seems wrong
			t = t-1;
			*T_new+=1;
		}
	}

}
