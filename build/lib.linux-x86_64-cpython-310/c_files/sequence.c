/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   22 February 1998
**      File:   sequence.c
**      Purpose: Routines for generating, reading and writing sequence of
**		observation symbols.
**      Organization: University of Maryland
**
**	Update:
**	Author: Tapas Kanungo
**	Purpose: To make calls to generic random number generators
**		and to change the seed everytime the software is executed.
**
**      $Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"

static char rcsid[] = "$Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $";

void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q)
{
        int     t = 1;
        int     q_t, o_t;
	hmmsetseed(seed);
        q[1] = GenInitalState(phmm);
        printf("%d",q[1]);
        O[1] = GenSymbol(phmm, q[1]);

        for (t = 2; t <= T; t++) {
                q[t] = GenNextState(phmm, q[t-1]);
                O[t] = GenSymbol(phmm, q[t]);
        }
}

int GenInitalState(HMM *phmm)
{
        double val, accum;
        int i, q_t;

        val = hmmgetrand();
        accum = 0.0;
        q_t = phmm->N;
        for (i = 1; i <= phmm->N; i++) {
                if (val < phmm->pi[i] + accum) {
                        q_t = i;
                        break;
                }
                else {
                                accum += phmm->pi[i];
                }
        }

        return q_t;
}

int GenNextState(HMM *phmm, int q_t)
{
        double val, accum;
        int j, q_next;

        val = hmmgetrand();
        accum = 0.0;
        q_next = phmm->N;
        for (j = 1; j <= phmm->N; j++) {
                if ( val < phmm->A[q_t][j] + accum ) {
                        q_next = j;
                        break;
                }
                else
                        accum += phmm->A[q_t][j];
        }

        return q_next;
}
int GenSymbol(HMM *phmm, int q_t)
{
        double val, accum;
        int j, o_t;

        val = hmmgetrand();
        accum = 0.0;
        o_t = phmm->M;
        for (j = 1; j <= phmm->M; j++) {
                if ( val < phmm->B[q_t][j] + accum ) {
                       o_t = j;
                       break;
                }
                else
                        accum += phmm->B[q_t][j];
        }

        return o_t;
}

void ReadSequence(FILE *fp, int *pT, int **pO)
{
        int *O;
        int i;

        fscanf(fp, "T= %d\n", pT);
        O = ivector(1,*pT);
        for (i=1; i <= *pT; i++)
                fscanf(fp,"%d", &O[i]);
        *pO = O;
}

void ReadSequences(FILE *fp, int *pN, int **pT, int ***pO) //STILL NEED TO FIGURE OUT HOW TO FREE THIS FROM MEMORY
{
  fscanf(fp, "N= %d\n", pN); //number of sequences
  int *O[*pN]; //jagged array of pointers to each sequence
  int *size; // array of sizes for each sequence

  size = ivector(0,*pN-1);
  //add sizes to size array and allocate memory to O accordingly
  for (int i=0; i < *pN; i++){
          fscanf(fp,"%d", &size[i]); //read each T and store in sizes
          //allocate memory in O
          O[i] = malloc(sizeof(int)*size[i]);
        }
  //add numbers to each sequence in O
  int k = 0;
  for (int i = 0; i < *pN; i++){
    int* p = O[i]; // start our pointer at the current sequence
    for (int j = 0; j<size[k]; j++){
      fscanf(fp,"%d", p);
      p++; //move pointer
    }
    fscanf(fp,"\n"); //should be at end of line after size[k] elements
    k++; // next index in size
  }
  *pT = size;
  *pO = O;
  //int* p = O[0];
  //printf("%d\n",*p);
  k = 0;
  int* p = O[1];
  printf("%d\n",*p);
  p++;
  printf("%d\n",*p);
  //built in check, remove when tested
  // Display elements in Jagged array
    for (int i = 0; i < *pN; i++) {

        int* p = O[i];
        for (int j = 0; j < size[k]; j++) {

            printf("%d ", *p);
            // move the pointer to the next element
            p++;
        }
        printf("\n");
        k++;
        // move the pointer to the next row
        //O[i]++;
    }
}
void PrintSequences(int *pN, int *size, int ***O){
  int k = 0;
  printf("in PrintSequences\n");
  printf("%d\n",size[0]);
  int* p = *O[0];
  printf("%d\n",*p);
  //built in check, remove when tested
  // Display elements in Jagged array
  for (int i = 0; i < *pN; i++) {

      int* p = *O[i];
      for (int j = 0; j < size[k]; j++) {

          printf("%d ", *p);
          // move the pointer to the next element
          p++;
      }
      printf("\n");
      k++;
      // move the pointer to the next row
      //O[i]++;
    }
}
void PrintSequence(FILE *fp, int T, int *O) //need to change this to print List
{
        int i;

        fprintf(fp, "T= %d\n", T);
        for (i=1; i <= T; i++)
                fprintf(fp,"%d ", O[i]);
	printf("\n");

}
void PrintPath(FILE *fp, int T_new, GList *list)
{
  fprintf(fp, "T= %d\n", T_new);
  GList *curr = g_list_last(list);
  while(curr!=NULL){
    fprintf(fp,"%d ", GPOINTER_TO_INT(curr->data));
    curr = g_list_previous(curr);
  }
  printf("\n");


}
void PrintLikelihoods(FILE *fp, int T_new, GList *likelihoods)
{
  GList *curr = g_list_last(likelihoods);
  while(curr!=NULL){
    fprintf(fp,"%d ", GPOINTER_TO_INT(curr->data));
    curr = g_list_previous(curr);
  }
  printf("\n");
}

void PrintAllOuts(FILE *fp, int T_new, GList *list,double ** delta, int N, int T){
  int i, j, k;
  for (i = 1; i <= N; i++) {
          for (j = 1; j <= T; j++) {
                  fprintf(fp, "%f ", delta[i][j] );
                }
        fprintf(fp, "\n");
      }
    fprintf(fp, "T= %d\n", T_new);
    GList *curr = g_list_last(list);
    while(curr!=NULL){
      fprintf(fp,"%d ", GPOINTER_TO_INT(curr->data));
      curr = g_list_previous(curr);
      }
    printf("\n");
}
