from libc.stdint cimport uint64_t

cdef extern from "int_iet.h":
    ctypedef struct label:
        pass
    ctypedef struct int_iet:
        pass
    ctypedef int_iet int_iet_t[1]

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
    int int_iet_num_cylinders(int_iet_t t)
