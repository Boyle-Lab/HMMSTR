#include <Python.h>
#include "libsais.h"

#include <stdio.h>

#include <string.h>
#include "arrayobject.h"

//#include <stdio.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/* C vector utility functions */ 
PyArrayObject *pyvector(PyObject *objin);
int *pyvector_to_Carrayptrs(PyArrayObject *arrayin);
int  not_intvector(PyArrayObject *vec);

/* Vector Utility functions */
PyArrayObject *pyvector(PyObject *objin) 
{
    return (PyArrayObject *) PyArray_ContiguousFromObject(objin, NPY_INT, 1, 1); 
}

/* Create 1D Carray from PyArray */
int *pyvector_to_Carrayptrs(PyArrayObject *arrayin) 
{
    return (int *) PyArray_DATA(arrayin);  /* pointer to arrayin data as double */
}

int compare_int(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

/* Check that PyArrayObject is an int type and a vector */ 
int  not_intvector(PyArrayObject *vec)
{
    
    if (PyArray_DESCR(vec)->type_num != NPY_INT || PyArray_NDIM(vec) != 1)
    {  
        PyErr_SetString(PyExc_ValueError, "Array must be of type Int and 1 dimensional (n).");
        return 1;
    }
    return 0;
}

static PyObject *python_sais(PyObject *self, PyObject *arg)
{
    const unsigned char *T;
    Py_ssize_t n;

    
    PyArrayObject *SA_np, *LCP_np;
    int *SA, *LCP;

    if(!PyUnicode_Check(arg))
    {
        PyErr_SetString(PyExc_StopIteration, "String must be of type Unicode.");
        return NULL;
    }

    if(!PyUnicode_IS_ASCII(arg))
    {
        PyErr_SetString(PyExc_StopIteration, "String must be limited to ASCII characters.");
        return NULL;
    }
    T = (const unsigned char *) PyUnicode_AsUTF8AndSize(arg, &n);
    if(T == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to convert string to UTF-8.");
        return NULL;
    }
    
    
    npy_intp dims[2];
    dims[0] = n;

    SA_np = (PyArrayObject *) PyArray_ZEROS(1, dims, NPY_INT, 0);
    SA = pyvector_to_Carrayptrs(SA_np);
    int res = libsais(T, SA, n, 0, NULL);
    if (res < 0)
    {
        PyErr_SetString(PyExc_StopIteration, "Error occurred in SA-IS.");
        return NULL;
    }
    
    LCP_np = (PyArrayObject *) PyArray_ZEROS(1, dims, NPY_INT, 0);
    LCP = pyvector_to_Carrayptrs(LCP_np);
    
    int32_t *plcp = malloc(n * sizeof(int32_t));  
    if (plcp == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to allocate memory.");
        return NULL;
    }
    
    res = libsais_plcp(T, SA, plcp, n);
    if (res < 0)
    {
        PyErr_SetString(PyExc_StopIteration, "Error occurred in PLCP creation.");
        return NULL;
    }
    res = libsais_lcp(plcp, SA, LCP, n);
    if (res < 0)
    {
        PyErr_SetString(PyExc_StopIteration, "Error occurred in LCP creation.");
        return NULL;
    }
    free(plcp);
    return Py_BuildValue("NN", SA_np, LCP_np);
}

static PyObject *python_sais_int(PyObject *self, PyObject *args)
{
    PyArrayObject *T_np, *SA_np;
    int *T, *SA;
    int i, k;
    if (!PyArg_ParseTuple(args, "O!i", &PyArray_Type, &T_np, &k))
        return NULL;
    if (T_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "T cannot be None.");
        return NULL;
    }
    if (not_intvector(T_np))
        return NULL;
    if (k <= 0)
    {
        PyErr_SetString(PyExc_StopIteration, "Alphabet size k must be greater than 0.");
        return NULL;
    }
    T = pyvector_to_Carrayptrs(T_np);
    int n = PyArray_DIMS(T_np)[0];
    for (i = 0; i < n; i++)
        if (T[i] < 0 || T[i] >= k)
        {
            PyErr_SetString(PyExc_StopIteration, "Array elements must be >= 0 and < k (alphabet size).");
            return NULL;
        }
    npy_intp dims[2];
    dims[0] = n;
    SA_np = (PyArrayObject *) PyArray_ZEROS(1, dims, NPY_INT, 0);
    SA = pyvector_to_Carrayptrs(SA_np);
    
    int res = libsais_int(T, SA, n, k, 0);
    if (res < 0)
    {
        PyErr_SetString(PyExc_StopIteration, "Error occurred in SA-IS.");
        return NULL;
    }

    return Py_BuildValue("N", SA_np);
}

PyObject *python_min_string(PyObject *self, PyObject *arg)
{
    const unsigned char *s;
    Py_ssize_t n;
    if(!PyUnicode_Check(arg))
    {
        PyErr_SetString(PyExc_StopIteration, "String must be of type Unicode.");
        return NULL;
    }

    if(!PyUnicode_IS_ASCII(arg))
    {
        PyErr_SetString(PyExc_StopIteration, "String must be limited to ASCII characters.");
        return NULL;
    }
    s = (const unsigned char *) PyUnicode_AsUTF8AndSize(arg, &n);
    if(s == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to convert string to UTF-8.");
        return NULL;
    }
    /*Booth's lexicographically minimal string rotation algorithm*/

	int i, j, k = 0;
    
    int *f = malloc(2 * n * sizeof(int));
    if(f == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to allocate memory.");
        return NULL;
    }

	for(j= 0; j	< 2 * n; j++)
	{
		f[j] = -1;
	}

	for(j = 1; j < 2 * n; j++)
	{
		i = f[j - k - 1];
		while(i != -1 && s[j % n] != s[(k + i + 1) % n])
        {
            if(s[j % n] < s[(k + i + 1) % n])
            {
                k = j - i - 1;
            }
            i = f[i];
        }
        if(i == -1 && s[j % n] != s[(k + i + 1) % n])
        {
            if(s[j % n] < s[(k + i + 1) % n])
            {
                k = j;
            }
            f[j - k] = -1;
        }
        else
        {
            f[j - k] = i + 1;
        }

	}
	free(f);
	PyObject * result = PyUnicode_New(n, 255);
	Py_UCS1 * buf = PyUnicode_1BYTE_DATA(result);
	for(i = 0; i < n; i++)
	{
	    buf[i] = s[(i + k) % n];

	}
	//buf[n] = '\0';

    return result;

}

PyObject *python_max_suffix(PyObject *self, PyObject *arg)
{
    PyArrayObject *res_np;
    const unsigned char *s;
    Py_ssize_t n;
    if(!PyUnicode_Check(arg))
    {
        PyErr_SetString(PyExc_StopIteration, "String must be of type Unicode.");
        return NULL;
    }

    if(!PyUnicode_IS_ASCII(arg))
    {
        PyErr_SetString(PyExc_StopIteration, "String must be limited to ASCII characters.");
        return NULL;
    }
    s = (const unsigned char *) PyUnicode_AsUTF8AndSize(arg, &n);
    if(s == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to convert string to UTF-8.");
        return NULL;
    }
 
    int *res;
    int i;
    npy_intp dims[2];
    dims[0] = n;
    res_np = (PyArrayObject *) PyArray_ZEROS(1, dims, NPY_INT, 0);
    res = pyvector_to_Carrayptrs(res_np);

    int cur_max_size = 0;
    for(i = n - 1; i >= 0; i--)
    {
        if(s[i] == '$' || s[i] == '#')
        {
            cur_max_size = 0;
        }
        else
        {
            cur_max_size++;
        }
        res[i] = cur_max_size;
    }
    return Py_BuildValue("N", res_np);
}


int _is_kmer_repetitive(const unsigned char *s, int len)
{
    int repeat_len;
    int repetitive = 0;
    for(repeat_len = 1; repeat_len <= (len / 2); repeat_len++)
    {
        //check if the string can have full copies of the substring of length repeat_len
        if(len % repeat_len != 0)
        {
            continue;
        }
        int nunit = len / repeat_len;

        //check if the string is repetitive

        repetitive = repeat_len;
        for(int i = 1; i < nunit; i++)
        {
            if(memcmp(s, s + i * repeat_len, repeat_len) != 0)
            {
                repetitive = 0;
                break;
            }
        }
        if(repetitive)
        {
            break;
        }
    }
    return repetitive;
}

PyObject *python_kmer_repetitive(PyObject * self, PyObject *args)
{
    const unsigned char *s;
    Py_ssize_t n;
    if(!PyUnicode_Check(args))
    {
        PyErr_SetString(PyExc_StopIteration, "String must be of type Unicode.");
        return NULL;
    }

    if(!PyUnicode_IS_ASCII(args))
    {
        PyErr_SetString(PyExc_StopIteration, "String must be limited to ASCII characters.");
        return NULL;
    }
    s = (const unsigned char *) PyUnicode_AsUTF8AndSize(args, &n);
    if(s == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to convert string to UTF-8.");
        return NULL;
    }
    int result = _is_kmer_repetitive(s, n);

    return Py_BuildValue("i", result);
}


PyObject *python_kmer_mask_potential(PyObject *self, PyObject *args)
{
    PyArrayObject *SA_np;
    PyArrayObject *index_np;
    PyArrayObject *mask_np;
    int *SA;
    int *index;
    int *mask;

    int kmer_len, kmer_cnt, sa_idx;
    int *kmer_start_pos, *kmer_consecutive_done;
    if (!PyArg_ParseTuple(args, "O!O!O!iii", &PyArray_Type, &SA_np, &PyArray_Type, &mask_np, &PyArray_Type, &index_np, &kmer_len, &sa_idx, &kmer_cnt))
        return NULL;
    
    if (SA_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "SA array cannot be None.");
        return NULL;
    }
    if (mask_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Mask array cannot be None.");
        return NULL;
    }
    if (index_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Index array cannot be None.");
        return NULL;
    }

    //check if SA_np and mask_np have same length
    if(PyArray_DIM(SA_np, 0) != PyArray_DIM(mask_np, 0))
    {
        PyErr_SetString(PyExc_StopIteration, "SA and mask arrays must have same length.");
        return NULL;
    }

    SA = pyvector_to_Carrayptrs(SA_np); 
    mask = pyvector_to_Carrayptrs(mask_np);
    index = pyvector_to_Carrayptrs(index_np);


    kmer_start_pos = malloc(kmer_cnt * sizeof(int));
    kmer_consecutive_done = malloc(kmer_cnt * sizeof(int));
   
    if(kmer_start_pos == NULL || kmer_consecutive_done == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to allocate memory.");
        return NULL;
    }

    int i;
    int real_kmer_cnt = 0;
    int pos;
    //get all the kmer start positions
    for(i = 0; i < kmer_cnt; i++)
    {
        pos = SA[sa_idx+i];
        if(mask[pos] >= kmer_len)
        {
            kmer_start_pos[real_kmer_cnt] = pos;
            real_kmer_cnt++;
        }

        kmer_consecutive_done[i] = 0;
    }
    //sort the kmer start positions
    qsort(kmer_start_pos, real_kmer_cnt, sizeof(int), compare_int);

    int mask_size = 0;
    int max_indiv_count = 0;
    int start_indiv_cum_count = 0;
    int max_consecutive_size = 0;
    int index_pos = 0;
    int last_masked = -1;
    for(i = 0; i < real_kmer_cnt; i++)
    {
        //check if we have passed the index (next sequence segment)
        while(kmer_start_pos[i] >= index[index_pos])
        {
            max_indiv_count = MAX(i - start_indiv_cum_count, max_indiv_count);
            start_indiv_cum_count = i;
            index_pos++;
        }


        //if previous masked kmer overlaps, then this position becomes masked, and should not be counted
        if(last_masked == -1 || (kmer_start_pos[i] >= (kmer_start_pos[last_masked] + kmer_len)))
        {
            mask_size += kmer_len;
            last_masked = i;
        }

        if(kmer_consecutive_done[i] == 0)
        {
            int cur_max_consecutive_size = kmer_len;
            int j=i + 1;
            int last_elem =i;
            while (j < real_kmer_cnt && kmer_start_pos[j] <= (kmer_start_pos[last_elem] + kmer_len))
            {
                if (kmer_start_pos[j] == (kmer_start_pos[last_elem] + kmer_len))
                {
                    cur_max_consecutive_size += kmer_len;
                    last_elem = j;
                    kmer_consecutive_done[j] = 1; //dont start a future round here
                }
                j++;
            }
            max_consecutive_size = MAX(cur_max_consecutive_size, max_consecutive_size);
        }

    }
    max_indiv_count = MAX(i - start_indiv_cum_count, max_indiv_count);
    

    free(kmer_start_pos);
    free(kmer_consecutive_done);




    return Py_BuildValue("iii", mask_size, max_indiv_count, max_consecutive_size / kmer_len);

}


PyObject *python_kmer_mask(PyObject *self, PyObject *args)
{
    PyArrayObject *SA_np;
    PyArrayObject *mask_np;
    const unsigned char *T;
    const unsigned char *replace;
    PyObject * result;
    int i, n;
    int *SA;
    int *mask;
    int min_consecutive;
    int kmer_len, kmer_cnt, sa_idx;
    int *kmer_start_pos, *kmer_consecutive_done;
    if (!PyArg_ParseTuple(args, "sO!O!iiiis", &T, &PyArray_Type, &SA_np, &PyArray_Type, &mask_np, &kmer_len, &sa_idx, &kmer_cnt, &min_consecutive, &replace))
        return NULL;
    
    if(strlen((const char *) replace) != 1)
    {
        PyErr_SetString(PyExc_StopIteration, "Replacement string must be a single character.");
        return NULL;
    }

    if (SA_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "SA array cannot be None.");
        return NULL;
    }
    if (mask_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Mask array cannot be None.");
        return NULL;
    }

    //check if SA_np and mask_np have same length
    if(PyArray_DIM(SA_np, 0) != PyArray_DIM(mask_np, 0))
    {
        PyErr_SetString(PyExc_StopIteration, "SA and mask arrays must have same length.");
        return NULL;
    }

    SA = pyvector_to_Carrayptrs(SA_np); 
    mask = pyvector_to_Carrayptrs(mask_np);
    n = strlen((const char *)T);

    kmer_start_pos = malloc(kmer_cnt * sizeof(int));
    kmer_consecutive_done = malloc(kmer_cnt * sizeof(int));
    result = PyUnicode_New(n, 127);

    if(kmer_start_pos == NULL || kmer_consecutive_done == NULL || result == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to allocate memory.");
        return NULL;
    }

	Py_UCS1 * buf = PyUnicode_1BYTE_DATA(result);
	memcpy(buf, T, n);


    int real_kmer_cnt = 0;
    int pos;

    //first get all the kmer start positions
    for(i = 0; i < kmer_cnt; i++)
    {
        pos = SA[sa_idx+i];
        if(mask[pos] >= kmer_len)
        {
            kmer_start_pos[real_kmer_cnt] = pos;
            kmer_consecutive_done[real_kmer_cnt] = 0;
            real_kmer_cnt++;
        }

    }
    qsort(kmer_start_pos, real_kmer_cnt, sizeof(int), compare_int);
    i = 0;
    int marked_pos = 0;
    while(i < real_kmer_cnt)
    {
        if(kmer_consecutive_done[i] == 0)
        {
            int consecutive_count = 1;
            int j=i + 1;
            int last_elem =i;
            
            while (j < real_kmer_cnt && kmer_start_pos[j] <= (kmer_start_pos[last_elem] + kmer_len))
            {
                if (kmer_start_pos[j] == (kmer_start_pos[last_elem] + kmer_len))
                {
                    kmer_consecutive_done[j] = 1; //dont start a future round here
                    consecutive_count += 1;
                    last_elem = j;
                }
                j++;
            }

            if(consecutive_count >= min_consecutive)
            {
                memset(buf + kmer_start_pos[i], replace[0], consecutive_count * kmer_len);
                kmer_consecutive_done[i] = 2;
                marked_pos += 1;
                //mark all matches but also all inbetween matches as done
                j = i+1;
                while(j < real_kmer_cnt && kmer_start_pos[j] < (kmer_start_pos[i] + consecutive_count * kmer_len))
                {
                    if((kmer_start_pos[j] - kmer_start_pos[i]) % kmer_len == 0)
                    {
                        kmer_consecutive_done[j] = 2; //positions that are marked as masked
                        marked_pos += 1;
                    }
                    else
                    {
                        kmer_consecutive_done[j] = 1; //inbetween positions, dont start a future round here
                    }
                    j++;
                }
            }

        }

        i++;
    
    }
    
    npy_intp dims[2];
    dims[0] = marked_pos;

    PyArrayObject *res_np = (PyArrayObject *) PyArray_ZEROS(1, dims, NPY_INT, 0);
    int *res = pyvector_to_Carrayptrs(res_np);
    pos = 0;
    for(i = 0; i < real_kmer_cnt; i++)
    {
        if(kmer_consecutive_done[i] == 2)
        {
            res[pos++] = kmer_start_pos[i];
        }
    }


    free(kmer_start_pos);
    free(kmer_consecutive_done);

    return Py_BuildValue("NN", result, res_np);

}

struct kmer_count
{
    int sa_index;
    int count;
};

PyObject *python_kmer_count(PyObject *self, PyObject *args)
{

	//SA, LCP, index, max_kmer_size, min_kmer, max_kmer
    PyArrayObject *SA_np, *LCP_np, *max_kmer_size_np;
    int *SA, *LCP, *max_kmer_size;
    int min_kmer, max_kmer;
    int i, min_count;
    Py_ssize_t j, n;
    const unsigned char *T;
    if (!PyArg_ParseTuple(args, "sO!O!O!iii", &T, &PyArray_Type, &SA_np, &PyArray_Type, &LCP_np, &PyArray_Type, &max_kmer_size_np, &min_kmer, &max_kmer, &min_count))
        return NULL;
    if (SA_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "SA array cannot be None.");
        return NULL;
    }
    if (LCP_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "LCP cannot be None.");
        return NULL;
    }
    if (max_kmer_size_np == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Max_kmer_size array cannot be None.");
        return NULL;
    }
    if (not_intvector(SA_np))
        return NULL;
    if (not_intvector(LCP_np))
        return NULL;
    if (not_intvector(max_kmer_size_np))
        return NULL;
    n = PyArray_DIM(SA_np, 0);
    if (n != (Py_ssize_t) strlen((const char *) T))
    {
        PyErr_SetString(PyExc_StopIteration, "SA and T arrays must be of same length.");
        return NULL;
    }
    if (n != PyArray_DIM(LCP_np, 0))
    {
        PyErr_SetString(PyExc_StopIteration, "SA and LCP arrays must be of same length.");
        return NULL;
    }
    if (n != PyArray_DIM(max_kmer_size_np, 0))
    {
        PyErr_SetString(PyExc_StopIteration, "SA and max_kmer_size arrays must be of same length.");
        return NULL;
    }
    SA = pyvector_to_Carrayptrs(SA_np);
    LCP = pyvector_to_Carrayptrs(LCP_np);
    max_kmer_size = pyvector_to_Carrayptrs(max_kmer_size_np);
	
	//generate array of active kmers
	//kmer class: kmer_size, sa_index, sa_count
	//for through suffix array
	//- up to lcp value: increment count of active kmers
	//- above lcp value: finish existing kmer (add to list), start new kmer at index (for each kmer size up to length of string)
	//finish all remaining kmers
    //return list of kmers



    struct kmer_count *active_kmers;
    active_kmers = (struct kmer_count *) malloc(sizeof(struct kmer_count) * (max_kmer + 1));
    if (active_kmers == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to allocate memory for active kmers.");
        return NULL;
    }
    for (i = min_kmer; i <= max_kmer; i++)
    {
        active_kmers[i].sa_index = -1;
        active_kmers[i].count = 0;
    }

    PyObject *kmer_list = PyList_New(0);
    PyObject *r;
    for(i = 0; i < n; i++)
    {
        for(j = min_kmer; j <= max_kmer; j++)
        {
            //kmer cannot be longer than this length at this location, if there are active kmers, need to store them. 
            if (j > LCP[i])
            {
                if (active_kmers[j].sa_index != -1 && active_kmers[j].count >= min_count)
                {
                    if (_is_kmer_repetitive(T + SA[active_kmers[j].sa_index], j) == 0)
                    {
                        r = Py_BuildValue("iii", j, active_kmers[j].sa_index, active_kmers[j].count);
                        PyList_Append(kmer_list, r);
                    }
                }
                if(j > max_kmer_size[MIN(SA[i],n-1)])
                {
                    active_kmers[j].sa_index = -1;
                }
                else
                {
                    active_kmers[j].sa_index = i;
                    active_kmers[j].count = 1;
                }
            }
            else
            { 
                //if j> max_kmer_size[MIN(SA[i],n-1)], we still count on,
                // as we need to finish the kmer. Filter out later.
                active_kmers[j].count++;
            }
        }
        
    }
    for(j = min_kmer; j <= max_kmer; j++)
    {
        if (active_kmers[j].sa_index != -1 && active_kmers[j].count >= min_count)
        {
            if (_is_kmer_repetitive(T + SA[active_kmers[j].sa_index], j)==0)
            {
                r = Py_BuildValue("iii", j, active_kmers[j].sa_index, active_kmers[j].count);
                PyList_Append(kmer_list, r);
            }
        }
    }
    free(active_kmers);
    return kmer_list;
}

PyObject * python_kmer_mask_simple(PyObject * self, PyObject *args)
{
    const unsigned char *T;
    const unsigned char *replace;
    const unsigned char *kmer;
    int n;
    int kmer_size;

    if (!PyArg_ParseTuple(args, "sss", &T, &kmer, &replace))
        return NULL;

    n = strlen((const char *) T);
    kmer_size = strlen((const char *) kmer);
    if(strlen((const char *) replace) != 1)
    {
        PyErr_SetString(PyExc_StopIteration, "Replace string must be a single character.");
        return NULL;
    }

    PyObject *kmer_positions = PyList_New(0);
    PyObject *result = PyUnicode_New(n, 127);

    if (kmer_positions == NULL || result == NULL)
    {
        PyErr_SetString(PyExc_StopIteration, "Unable to allocate memory.");
        return NULL;
    }

	Py_UCS1 * buf = PyUnicode_1BYTE_DATA(result);
	memcpy(buf, T, n);
    int i = 0;
    

    while (i < n)
    {
        if (strncmp((const char *)T + i, (const char *) kmer, kmer_size) == 0)
        {
            PyList_Append(kmer_positions, Py_BuildValue("i", i));
            memset(buf + i, replace[0], kmer_size);
            i += kmer_size;
        }
        else
        {
            i++;
        }
    }
    return Py_BuildValue("NN", result, kmer_positions);

}



static PyMethodDef ModuleMethods[] = {
    {"sais",  python_sais, METH_O, "Construct a Suffix Array for a given string.\n:param string T : character string for which SA should be constructed.\n:returns ndarray SA : constructed suffix array."},
    {"sais_int",  python_sais_int, METH_VARARGS, "Construct a Suffix Array for a given NumPy integer array.\n:param ndarray T : int array for which SA should be constructed.\n:param k : alphabet size. All integers in T must be >= 0 and < k.\n:return ndarray SA : suffix array for T."},
    {"min_string",  python_min_string, METH_O, "Returns the least lexicographic cyclical rotation index of a string.\n:param s: str.\n:returns str: string rotated to the minimal lexicographic position;"},
    {"max_suffix",  python_max_suffix, METH_O, "Returns the maximum suffix length bounded by string separator $, # or string end.\n:param s: str.\n:returns A: array of max kmer size at each location of string s;"},
    {"kmer_count",  python_kmer_count, METH_VARARGS, "Constructs a list of kmers from a suffix and LCP array.\nFilters out repetitive k-mers. \n:param string T: character string for which SA was constructed.\n:param ndarray SA: sufffix array\n:param ndarray LCP: LCP array\n:param ndarray mask: masking array\n:param int k_min: minimal k-mer size\n:param int k_max: maximum k-mer size\n:param int min_count: minimal numer of occurences (possibly overlapping).\n:returns list kmers: list of tuples(kmer_len, sa_array_idx, kmer_cnt)."},
    {"kmer_mask_potential",  python_kmer_mask_potential, METH_VARARGS, "Determines masking potential of a k-mer in the seqeuence.\n:param ndarray SA: suffix array.\n:param ndarray mask: masking array.\n:param ndarray index: end positions of each substring.\n:param kmer_len: length of kmer.\n:param int sa_index: Index in suffix arra of kmer start.\n:param int kmer_cnt: number of kmer occurences in suffix array.\n:return tuple of (total_masked (excluding overlaps), max_count_in_indiv_seq (max kmer count per seq), max_continuous_masked)."},
    {"kmer_mask",  python_kmer_mask, METH_VARARGS, "Masks a k-mer in a seqeuence using a suffix array.\n:param str T: sequence.\n:param ndarray SA: suffix array.\n:param ndarray mask: masking array.\n:param kmer_len: length of kmer.\n:param int sa_index: Index in suffix arra of kmer start.\n:param int kmer_cnt: number of kmer occurences in suffix array.\n:param int min_consecutive: only mask regions with at least 'min_consecutive' k-mers'.\n:param str replace: character to replace with.\n:return tuple of (new_string, start of k-mer positions that were masked)."},
    {"kmer_mask_simple",  python_kmer_mask_simple, METH_VARARGS, "Masks a k-mer in a sequence.\n:param str T: sequence.\n:param str kmer: kmer to replace.\n:param str replace: character to replace with.\n:return tuple of (new_string, list of masked start positions)."},
    {"kmer_repetitive", python_kmer_repetitive, METH_O, "Determines if a kmer consists of a repetitive sub-kmer.\n:param s : kmer (only ASCII supported).\n:return int: 0 if not repetitive, otherwise length of minimal repetitive sub-kmer."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef mod =
{
    PyModuleDef_HEAD_INIT,
    "pylibsais",
    "Suffix Array and kmer library using the SA-IS algorithm",
    -1,
    ModuleMethods
};

PyMODINIT_FUNC PyInit_pylibsais(void)
{
    import_array(); // NumPy
    return PyModule_Create(&mod);
}
#else

PyMODINIT_FUNC initpylibsais(void)
{
    (void) Py_InitModule("pylibsais", ModuleMethods);
    import_array(); // NumPy
}

#endif
