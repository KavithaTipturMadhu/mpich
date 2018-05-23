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

typedef struct {
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

int MPIR_Netloc_parse_topology(const char *attributes_file, netloc_topology_t topology,
                               MPIR_Netloc_network_attributes * network_attr);

int MPIR_get_network_end_point(netloc_topology_t netloc_topology,
                               hwloc_topology_t hwloc_topology, hwloc_cpuset_t bindset,
                               netloc_node_t ** end_point);

/* Index of host node in the list of host nodes in the network topology object*/
int MPIR_get_host_index_in_network(netloc_topology_t netloc_topology,
                                   hwloc_topology_t hwloc_topology, hwloc_cpuset_t bindset,
                                   int *index);

int MPIR_get_switches_at_level(netloc_topology_t netloc_topology,
                               MPIR_Netloc_network_attributes attributes, int level,
                               netloc_node_t *** switches_at_level, int *switch_count);

#endif /* NETLOC_UTIL_H_INCLUDED */
