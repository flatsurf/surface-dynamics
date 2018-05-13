from libc.stdint cimport uint64_t

cdef extern from "int_iet.h":
    ctypedef struct label:
        pass
    ctypedef struct int_iet:
        pass
    ctypedef int_iet int_iet_t[1]

    ctypedef struct li_vector_iterator:
        uint64_t * x
    ctypedef li_vector_iterator li_vector_iterator_t[1]

    # memory allocation
    void int_iet_init(int_iet_t t, unsigned int n)
    void int_iet_clear(int_iet_t t)

    # safety check
    int  int_iet_check(int_iet_t t)

    # set data
    void int_iet_set_labels_and_twin(int_iet_t t, int * labels, int * twin, int k)
    void int_iet_set_lengths(int_iet_t t, uint64_t * lengths)

    # output
    void int_iet_print(int_iet_t t)

    # number of cylinders
    int int_iet_num_cylinders(uint64_t * widths, int_iet_t t)

    # iterator through vectors of given sum n and length k
    int int_vector_first(uint64_t * x, int n, int k)
    int int_vector_next(uint64_t * x, int n, int k)

    #iteration through integer vectors of given sum and length
    void int_li_vector_init(li_vector_iterator_t t, uint64_t n, int kfree, int ktop, int kbot)
    void int_li_vector_clear(li_vector_iterator_t t);
    int int_li_vector_prefirst(li_vector_iterator_t t)
    int int_li_vector_first_or_next(li_vector_iterator_t t)
    void int_li_vector_info(li_vector_iterator_t t)

