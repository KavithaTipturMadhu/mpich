/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 */

#ifndef NETLOC_UTIL_H_INCLUDED
#define NETLOC_UTIL_H_INCLUDED

#include "netloc.h"

typedef enum {
    MPIR_NETLOC_NETWORK_TYPE__FAT_TREE,
    MPIR_NETLOC_NETWORK_TYPE__TORUS,
    MPIR_NETLOC_NETWORK_TYPE__INVALID,
} MPIR_Netloc_network_type;

typedef union {
    MPIR_Netloc_network_type type;

    union {
        struct {
            /* the levels at which the host and switch nodes are in
             * the network */
            int *node_levels;
        } fat_tree;
        struct {
            int dimension;
            int *geometry;
        } torus;
    } u;
} MPIR_Netloc_network_attributes;

int MPIR_Netloc_parse_topology(netloc_topology_t topology,
                               MPIR_Netloc_network_attributes * network_attr);

#endif /* NETLOC_UTIL_H_INCLUDED */
