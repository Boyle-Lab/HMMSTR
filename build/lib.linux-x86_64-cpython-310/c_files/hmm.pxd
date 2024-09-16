from libc.stdio cimport FILE
cdef extern from "glib.h":
    ctypedef struct _GList:
        void* data
        _GList* next
        _GList* prev

    ctypedef _GList* GList

cdef extern from "glist.h":
    # Doubly linked lists
    GList* g_list_alloc() 
    void g_list_free(GList* list)
    void g_list_free_1(GList* list)
    void g_list_free_full(GList* list, void* free_func)
    GList* g_list_append(GList* list, void* data)
    GList* g_list_prepend(GList* list, void* data)
    GList* g_list_insert(GList* list, void* data, int position)
    GList* g_list_insert_sorted(GList* list, void* data, void* func)
    GList* g_list_insert_sorted_with_data(GList* list, void* data, void* func, void* user_data)
    GList* g_list_insert_before(GList* list, GList* sibling, void* data)
    GList* g_list_insert_before_link(GList* list, GList* sibling, GList* link_)
    GList* g_list_concat(GList* list1, GList* list2)
    GList* g_list_remove(GList* list, void* data)
    GList* g_list_remove_all(GList* list, void* data)
    GList* g_list_remove_link(GList* list, GList* llink)
    GList* g_list_delete_link(GList* list, GList* link_)
    GList* g_list_reverse(GList* list)
    GList* g_list_copy(GList* list)
    GList* g_list_copy_deep(GList* list, void* func, void* user_data)
    GList* g_list_nth(GList* list, unsigned int n)
    GList* g_list_nth_prev(GList* list, unsigned int n)
    GList* g_list_find(GList* list, void* data)
    GList* g_list_find_custom(GList* list, void* data, void* func)
    int g_list_position(GList* list, GList* llink)
    int g_list_index(GList* list, void* data)
    GList* g_list_last(GList* list)
    GList* g_list_first(GList* list)
    unsigned int g_list_length(GList* list)
    void g_list_foreach(GList* list, void* func, void* user_data)
    GList* g_list_sort(GList* list, void* compare_func)
    GList* g_list_sort_with_data(GList* list, void* compare_func, void* user_data)
    void* g_list_nth_data(GList* list, unsigned int n)
    void g_clear_list(GList** list_ptr, void* destroy)

#define g_list_previous(list) ((list) ? (((GList *)(list))->prev) : NULL)
#define g_list_next(list) ((list) ? (((GList *)(list))->next) : NULL)
cdef extern from "nrutil.h":
    float *vector(int nl,int nh)
    float **matrix(int nrl,int nrh,int ncl,int nch)
    float **convert_matrix()
    double *dvector(int nl,int nh)
    double **dmatrix(int nrl,int nrh,int ncl,int nch)
    int *ivector(int nl, int nh)
    int **imatrix(int nrl,int nrh,int ncl,int nch)
    float **submatrix(float** a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl)
    void free_vector(float* v,int nl, int nh)
    void free_dvector(double* v, int nl, int nh)
    void free_ivector(int* v,int nl,int nh)
    void free_matrix(float** m,int nrl, int nrh, int ncl, int nch)
    void free_dmatrix(double** m, int nrl, int nrh, int ncl, int nch)
    void free_imatrix(int** m, int nrl, int nrh, int ncl, int nch)
    void free_submatrix(float** b, int nrl, int nrh, int ncl, int nch)
    void free_convert_matrix(float* a, int nrl, int nrh, int ncl, int nch)
    void nrerror(char errortext[])
cdef extern from "hmm.h":
    ctypedef struct HMM:
        int N
        int M
        double** A
        double** B
        double* pi

    void ReadHMM(FILE* fp, HMM* phmm)
    void PrintHMM(FILE* fp, HMM* phmm)
    void InitHMM(HMM* phmm, int N, int M, int seed)
    void CopyHMM(HMM* phmm1, HMM* phmm2)
    void FreeHMM(HMM* phmm)
    void PrintDelta(FILE* fp, double** delta, HMM* phmm, int T)
    void PrintTraceBack(FILE* fp, int** delta, HMM* phmm, int T)
    void ReadSequence(FILE* fp, int* pT, int** pO)
    void ReadSequences(FILE* fp, int* pN, int** pT, int*** pO)
    void PrintSequences(int* pN, int* size, int*** O)
    void PrintSequence(FILE* fp, int T, int* O)
    void PrintPath(FILE* fp, int T_new, GList* list)
    void PrintLikelihoods(FILE* fp, int T_new, GList* likelihoods)
    void GenSequenceArray(HMM* phmm, int seed, int T, int* O, int* q)
    int GenInitalState(HMM* phmm)
    int GenNextState(HMM* phmm, int q_t)
    int GenSymbol(HMM* phmm, int q_t)
    void PrintAllOuts(FILE* fp, int T_new, GList* list, double** delta, int N, int T)

    void Viterbi(HMM* phmm, int T, int* O, double** delta, int** psi, int** traceback_dir, GList** list, int* T_new, int* q, double* pprob)
    void ViterbiLog(HMM* phmm, int T, int* O, double** delta, int** psi, int** traceback_dir, GList** list, GList** likelihoods_ints, GList** likelihoods_dec, int* T_new, int* q, double* pprob, int repeat_len, int prefix_idx)

    int hmmgetseed()
    void hmmsetseed(int seed)
    double hmmgetrand()

    # Macros
    # Note: These will be inlined by the C compiler
    int MAX(int x, int y)
    int MIN(int x, int y)

# Wrapper functions for the HMM functions
cdef void py_ReadHMM(FILE* fp, HMM* phmm)
cdef void py_PrintHMM(FILE* fp, HMM* phmm)
cdef void py_InitHMM(HMM* phmm, int N, int M, int seed)
cdef void py_CopyHMM(HMM* phmm1, HMM* phmm2)
cdef void py_FreeHMM(HMM* phmm)
cdef void py_PrintDelta(FILE* fp, double** delta, HMM* phmm, int T)
cdef void py_PrintTraceBack(FILE* fp, int** delta, HMM* phmm, int T)
cdef void py_ReadSequence(FILE* fp, int* pT, int** pO)
cdef void py_ReadSequences(FILE* fp, int* pN, int** pT, int*** pO)
cdef void py_PrintSequences(int* pN, int* size, int*** O)

cdef void py_PrintSequence(FILE* fp, int T, int* O)
cdef void py_PrintPath(FILE* fp, int T_new, GList* list)
cdef void py_PrintLikelihoods(FILE* fp, int T_new, GList* likelihoods)
cdef void py_GenSequenceArray(HMM* phmm, int seed, int T, int* O, int* q)
cdef int py_GenInitalState(HMM* phmm)
cdef int py_GenNextState(HMM* phmm, int q_t)
cdef int py_GenSymbol(HMM* phmm, int q_t)
cdef void py_PrintAllOuts(FILE* fp, int T_new, GList* list, double** delta, int N, int T)

cdef void py_Viterbi(HMM* phmm, int T, int* O, double** delta, int** psi, int** traceback_dir, GList** list, int* T_new, int* q, double* pprob)

cdef void py_ViterbiLog(HMM* phmm, int T, int* O, double** delta, int** psi, int** traceback_dir, GList** list, GList** likelihoods_ints, GList** likelihoods_dec, int* T_new, int* q, double* pprob, int repeat_len, int prefix_idx)

cdef int py_hmmgetseed()
cdef void py_hmmsetseed(int seed)
cdef double py_hmmgetrand()
# Wrapper for the `free` function
cdef extern from "stdlib.h":
    void free(void* ptr)

cdef void py_free(void* ptr)
