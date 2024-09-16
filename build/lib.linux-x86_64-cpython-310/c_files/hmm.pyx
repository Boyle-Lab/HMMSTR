# Import the converted pxd declarations
cimport c_files.hmm as chmm
import cython
from libc.stdio cimport FILE,stdout, sscanf, fclose, fopen, fprintf, stderr
from libc.stdlib cimport malloc, free, atoi
# Wrapper functions for the HMM functions
cdef void py_ReadHMM(FILE* fp, chmm.HMM* phmm):
    chmm.ReadHMM(fp, phmm)

cdef void py_PrintHMM(FILE* fp, chmm.HMM* phmm):
    chmm.PrintHMM(fp, phmm)

cdef void py_InitHMM(chmm.HMM* phmm, int N, int M, int seed):
    chmm.InitHMM(phmm, N, M, seed)

cdef void py_CopyHMM(chmm.HMM* phmm1, chmm.HMM* phmm2):
    chmm.CopyHMM(phmm1, phmm2)

cdef void py_FreeHMM(chmm.HMM* phmm):
    chmm.FreeHMM(phmm)

cdef void py_PrintDelta(FILE* fp, double** delta, chmm.HMM* phmm, int T):
    chmm.PrintDelta(fp, delta, phmm, T)

cdef void py_PrintTraceBack(FILE* fp, int** delta, chmm.HMM* phmm, int T):
    chmm.PrintTraceBack(fp, delta, phmm, T)

cdef void py_ReadSequence(FILE* fp, int* pT, int** pO):
    chmm.ReadSequence(fp, pT, pO)

cdef void py_ReadSequences(FILE* fp, int* pN, int** pT, int*** pO):
    chmm.ReadSequences(fp, pN, pT, pO)

cdef void py_PrintSequences(int* pN, int* size, int*** O):
    chmm.PrintSequences(pN, size, O)

cdef void py_PrintSequence(FILE* fp, int T, int* O):
    chmm.PrintSequence(fp, T, O)

cdef void py_PrintPath(FILE* fp, int T_new, chmm.GList* list):
    chmm.PrintPath(fp, T_new, list)

cdef void py_PrintLikelihoods(FILE* fp, int T_new, chmm.GList* likelihoods):
    chmm.PrintLikelihoods(fp, T_new, likelihoods)

cdef void py_GenSequenceArray(chmm.HMM* phmm, int seed, int T, int* O, int* q):
    chmm.GenSequenceArray(phmm, seed, T, O, q)

cdef int py_GenInitalState(chmm.HMM* phmm):
    return chmm.GenInitalState(phmm)

cdef int py_GenNextState(chmm.HMM* phmm, int q_t):
    return chmm.GenNextState(phmm, q_t)

cdef int py_GenSymbol(chmm.HMM* phmm, int q_t):
    return chmm.GenSymbol(phmm, q_t)

cdef void py_PrintAllOuts(FILE* fp, int T_new, chmm.GList* list, double** delta, int N, int T):
    chmm.PrintAllOuts(fp, T_new, list, delta, N, T)

cdef void py_Viterbi(chmm.HMM* phmm, int T, int* O, double** delta, int** psi, int** traceback_dir, chmm.GList** list, int* T_new, int* q, double* pprob):
    chmm.Viterbi(phmm, T, O, delta, psi, traceback_dir, list, T_new, q, pprob)

cdef void py_ViterbiLog(chmm.HMM* phmm, int T, int* O, double** delta, int** psi, int** traceback_dir, chmm.GList** list, chmm.GList** likelihoods_ints, chmm.GList** likelihoods_dec, int* T_new, int* q, double* pprob, int repeat_len, int prefix_idx):
    chmm.ViterbiLog(phmm, T, O, delta, psi, traceback_dir, list, likelihoods_ints, likelihoods_dec, T_new, q, pprob, repeat_len, prefix_idx)

cdef int py_hmmgetseed():
    return chmm.hmmgetseed()

cdef void py_hmmsetseed(int seed):
    chmm.hmmsetseed(seed)

cdef double py_hmmgetrand():
    return chmm.hmmgetrand()

cdef void py_free(void* ptr):
    free(ptr)
cdef void cy_main(bytes model_file, bytes obs_seq_length, bytes obs_seq,bytes repeat_len_str,bytes prefix_idx_str):
    cdef int t, T, T_new, repeat_len, prefix_idx
    cdef chmm.HMM hmm
    cdef int* O
    cdef int* q
    cdef double** delta
    cdef int** psi
    cdef int** traceback_dir
    cdef double proba, logproba
    cdef FILE *fp, *traceback_del, *traceback, *delta_out
    cdef int i, n, offset
    cdef char* S
    cdef chmm.GList *list_log = NULL
    cdef chmm.GList *likelihoods_ints = NULL
    cdef chmm.GList *likelihoods_dec = NULL

    #if argc != 6:
     #   print("Usage error")
      #  print("Usage: testvit <model.hmm> <obs.seq length> <obs.seq>")
      #  return

    # open HMM file
    fp = fopen(model_file, "r")
    if fp == NULL:
        fprintf(stderr, "Error: File %s not found\n", model_file)
        return

    # read HMM file according to ReadHMM
    chmm.ReadHMM(fp, &hmm)
    fclose(fp)

    T = atoi(obs_seq_length) # read T from input
    S = obs_seq
    repeat_len =  atoi(repeat_len_str)# for calculation of RDn case
    prefix_idx = atoi(prefix_idx_str)# for calculation of RDn case

    O = chmm.ivector(1, T)
    i = 1
    while (sscanf(S, "%d%n", &n, &offset) == 1):
        O[i] = n
        S += offset
        i += 1

    q = chmm.ivector(1, T)
    traceback_dir = chmm.imatrix(1, T, 1, hmm.N)
    delta = chmm.dmatrix(1, T, 1, hmm.N)
    psi = chmm.imatrix(1, T, 1, hmm.N)

    chmm.ViterbiLog(&hmm, T, O, delta, psi, traceback_dir, &list_log,&likelihoods_ints,&likelihoods_dec, &T_new, q, &logproba, repeat_len, prefix_idx)

    fprintf(stdout, "Viterbi  MLE log prob = %E\n", logproba)
    fprintf(stdout, "Optimal state sequence:\n")
    chmm.PrintPath(stdout, T_new, list_log)
    chmm.PrintLikelihoods(stdout, T_new, likelihoods_ints)
    chmm.PrintLikelihoods(stdout, T_new, likelihoods_dec)

    chmm.free_ivector(q, 1, T)
    chmm.free_ivector(O, 1, T)
    chmm.free_imatrix(psi, 1, T, 1, hmm.N)
    chmm.free_imatrix(traceback_dir, 1, T, 1, hmm.N)
    chmm.g_list_free(list_log);
    chmm.g_list_free(likelihoods_ints);
    chmm.g_list_free(likelihoods_dec);
    chmm.free_dmatrix(delta, 1, T, 1, hmm.N)
    chmm.py_FreeHMM(&hmm)

# Example usage of the HMM struct and functions
def example_usage():
    cdef chmm.HMM hmm
    cdef int T = 10
    cdef int* O_ptr = <int*>malloc((T+1) * sizeof(int))
    cdef int* p_ptr = <int*>malloc((T+1) * sizeof(int))
    # Initialize the HMM
    chmm.py_InitHMM(&hmm, 3, 2, 123)
    print("made it past initHMM")
    # Generate a sequence using the HMM
    chmm.py_GenSequenceArray(&hmm, 456, T, O_ptr, p_ptr)
    print("made it past py_GenSequenceArray")
    # Convert the sequence to a Python list
    O = [O_ptr[i] for i in range(1,T+1)]
    p = [p_ptr[i] for i in range(1,T+1)]
    print(p)
    print("made it past making O a list")
    # Print the generated sequence
    print("Generated Sequence:")
    print(O)

    # Clean up
    chmm.py_FreeHMM(&hmm)
    free(O_ptr)
    free(p_ptr)
    print("O freed")

# Run the example usage
#example_usage()
def run_main(model_file,obs_seq_length,obs_seq,repeat_len,prefix_idx):
    # ...
    model_file_bytes = model_file.encode('utf-8')
    obs_seq_length_bytes = obs_seq_length.encode('utf-8')
    obs_seq_bytes = obs_seq.encode('utf-8')
    repeat_len_bytes = repeat_len.encode('utf-8')
    prefix_idx_bytes = prefix_idx.encode('utf-8')

    cy_main(model_file_bytes, obs_seq_length_bytes, obs_seq_bytes, repeat_len_bytes, prefix_idx_bytes)

    #model_file = b'test_out/CANVAS_random_coord/test_outs_CANVAS.hmm'
    #obs_seq_length = b'831'
    #obs_seq = b'2 3 4 2 1 1 1 1 1 3 4 4 1 1 4 1 4 1 2 4 1 1 3 2 1 1 2 4 1 1 4 1 1 1 3 3 4 4 2 3 1 4 3 3 1 1 1 1 4 4 4 1 3 3 1 4 1 1 3 2 2 4 1 1 4 1 2 2 2 4 2 1 1 1 1 2 2 3 2 3 1 4 3 3 2 4 2 3 3 4 2 1 2 2 2 3 1 4 1 1 1 1 2 4 1 1 1 1 1 4 3 1 1 4 2 2 3 2 1 4 4 3 3 1 4 4 2 4 2 4 4 2 4 4 3 2 3 1 3 4 3 3 2 4 2 1 1 2 3 3 3 1 4 3 1 2 2 4 2 4 4 4 1 4 1 3 4 4 4 3 3 1 1 2 3 1 3 2 2 3 1 4 4 2 3 1 4 4 1 4 2 2 2 4 1 4 1 3 3 1 4 3 2 2 4 4 3 3 1 1 2 1 2 4 4 2 1 1 1 1 1 3 3 3 3 1 2 3 2 3 2 1 4 4 1 1 1 1 1 2 1 3 1 1 1 1 1 2 2 1 4 3 2 4 4 1 3 1 2 4 4 2 4 4 3 1 3 4 3 4 1 3 2 1 2 1 4 2 3 3 3 1 4 3 2 1 3 2 3 1 4 4 1 4 4 3 2 4 1 4 4 3 1 4 4 1 4 1 2 2 3 1 3 2 2 4 1 1 3 3 2 4 4 4 1 4 4 3 1 4 1 4 4 2 2 4 3 1 4 2 4 1 4 3 3 4 1 4 1 2 2 4 3 4 3 3 1 3 2 4 3 1 3 2 3 3 1 4 2 3 2 4 4 4 3 1 1 3 1 4 1 4 3 1 1 4 1 3 2 3 2 4 2 2 2 3 1 1 1 1 1 1 1 1 4 1 1 1 1 4 1 1 1 1 1 4 1 1 1 1 4 1 1 1 1 4 1 1 1 1 4 1 1 1 4 1 1 1 1 1 4 3 1 2 4 2 2 3 2 1 1 1 4 1 4 1 1 2 1 4 3 1 1 3 4 4 2 4 3 1 4 3 2 4 4 1 4 2 1 2 3 1 2 4 1 2 1 4 1 4 1 2 2 1 2 4 4 4 1 2 2 1 2 1 3 1 4 4 3 1 1 1 1 1 3 1 3 2 4 2 3 1 4 2 2 2 4 1 1 2 2 4 1 1 4 1 4 4 1 2 4 4 1 4 1 1 1 4 4 1 3 1 1 1 1 2 4 1 1 1 4 1 1 4 4 3 2 4 2 3 1 4 1 3 2 2 3 2 4 1 4 1 2 2 3 3 1 4 2 4 1 3 1 4 4 1 3 1 1 3 1 3 1 4 3 2 1 2 2 3 2 2 3 1 1 4 1 1 1 1 1 4 4 4 1 1 4 1 1 2 4 4 3 3 3 3 1 1 1 1 4 3 1 1 2 2 3 1 4 1 4 1 2 3 3 1 1 3 3 1 2 3 1 4 1 1 3 2 3 3 3 2 3 2 2 3 1 1 1 1 4 1 2 4 4 1 4 3 2 2 4 4 2 2 2 3 1 1 3 1 4 1 2 4 1 4 1 3 1 4 3 3 1 3 3 1 3 3 3 3 1 1 1 4 3 3 2 2 4 2 4 4 4 3 1 4 4 1 2 4 4 3 3 3 1 4 3 1 4 1 4 3 2 4 3 1 4 1 1 4 3 1 4 4 1 3 1 4 3 2 4 3 3 2 4 2 4 1 3 2 4 2 4 4 4 4 4 2 1 4 4 4 3 3 1 3 2 1 1 3 3 3 1 4 2 4 4 4 2 3 2 1 4 1 4 4 4 3 1 4 4 4 3 1 2 3 1 1 4 3 3'
    #repeat_len = b'5'
    #prefix_idx = b'91'

    #cy_main(model_file, obs_seq_length, obs_seq,repeat_len,prefix_idx)

    # ...
