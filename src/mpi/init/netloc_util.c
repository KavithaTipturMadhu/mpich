/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2018 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 */

#include "mpiimpl.h"

#ifdef HAVE_NETLOC
#include "netloc_util.h"
#include "mpl.h"

static int get_fat_tree_attributes(netloc_topology_t topology,
                                   MPIR_Netloc_network_attributes * network_attr)
{
    int i, j;
    netloc_node_t **traversal_order;
    int traversal_start = 0, traversal_end = 0;
    int visited_count = 0;
    int *node_levels;

    /* indexed by node id, visited_node_list[i] > 1 indicates that the
     * node i has been visited */
    network_attr->u.fat_tree.node_levels = (int *) MPL_malloc(sizeof(int) * topology->num_nodes, MPL_MEM_OTHER);

    for (i = 0; i < topology->num_nodes; i++) {
        network_attr->u.fat_tree.node_levels[i] = -1;
    }

    /* Traversal order never exceeds the total number of nodes anyway */
    traversal_order = (netloc_node_t **) MPL_malloc(sizeof(unsigned) * topology->num_nodes, MPL_MEM_OTHER);

    /* Initialize host depths to zero */
    for (i = 0; i < topology->num_nodes; i++) {
        if (topology->nodes[i]->node_type == NETLOC_NODE_TYPE_HOST) {
            int num_edges;
            netloc_edge_t **edges;

            /* Mark the host node as visited */
            network_attr->u.fat_tree.node_levels[topology->nodes[i]->__uid__] = 0;
            visited_count++;

            /* Copy all parents without duplicates to the traversal list */
            netloc_get_all_edges(topology, topology->nodes[i], &num_edges, &edges);
            for (j = 0; j < num_edges; j++) {
                if (network_attr->u.fat_tree.node_levels[edges[j]->dest_node->__uid__] < 0) {
                    traversal_order[traversal_end++] = edges[j]->dest_node;
                    network_attr->u.fat_tree.node_levels[edges[j]->dest_node->__uid__] = 1;
                    visited_count++;
                }
            }
        }
    }

    while (visited_count < topology->num_nodes) {
        for (i = traversal_start; i < traversal_end; i++) {
            netloc_node_t *traversed_node_index = traversal_order[i];
            int num_edges;
            netloc_edge_t **edges;
            int depth = network_attr->u.fat_tree.node_levels[traversed_node_index->physical_id_int];

            /* find all nodes not visited with an edge from the current node */
            netloc_get_all_edges(topology, traversed_node_index, &num_edges, &edges);
            for (j = 0; j < num_edges; j++) {
                if (network_attr->u.fat_tree.node_levels[edges[j]->dest_node->__uid__] < 0 ||
                    (depth + 1) < network_attr->u.fat_tree.node_levels[edges[j]->dest_node->__uid__]) {
                    traversal_order[traversal_end++] = edges[j]->dest_node;
                    network_attr->u.fat_tree.node_levels[edges[j]->dest_node->physical_id_int] = depth + 1;
                    visited_count++;
                }
            }
            traversal_start++;
        }
    }
    return 0;
}

#undef FUNCNAME
#define FUNCNAME MPIR_Netloc_parse_topology
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIR_Netloc_parse_topology(netloc_topology_t topology,
                               MPIR_Netloc_network_attributes * network_attr)
{
    int i, j;
    unsigned *switch_node_depth;
    netloc_dt_lookup_table_t *host_nodes =
        (netloc_dt_lookup_table_t *) MPL_malloc(sizeof(netloc_dt_lookup_table_t), MPL_MEM_OTHER);
    int host_node_count = 0;
    struct netloc_dt_lookup_table_iterator *hti = NULL;
    int node_query_ret;
    int mpi_errno = MPI_SUCCESS;

    node_query_ret = netloc_get_all_host_nodes(topology, host_nodes);
    MPIR_ERR_CHKANDJUMP(node_query_ret, mpi_errno, MPI_ERR_OTHER, "**netloc_topo_load");

    host_node_count = netloc_lookup_table_size(*host_nodes);

    network_attr->type = MPIR_NETLOC_NETWORK_TYPE__FAT_TREE;
    hti = netloc_dt_lookup_table_iterator_t_construct(*host_nodes);
    while (!netloc_lookup_table_iterator_at_end(hti)) {
        netloc_node_t *node = (netloc_node_t *) netloc_lookup_table_iterator_next_entry(hti);
        int num_edges;
        netloc_edge_t **edges = (netloc_edge_t **) MPL_malloc(sizeof(netloc_edge_t *), MPL_MEM_OTHER);

        if (node == NULL)
            break;
        netloc_get_all_edges(topology, node, &num_edges, &edges);
        if (num_edges > 0) {
            /* Find the destination node */
            netloc_node_t *dest_node = NULL;
            for (i = 0; i < num_edges; i++) {
                if (edges[i]->dest_node->node_type != NETLOC_NODE_TYPE_SWITCH) {
                    network_attr->type = MPIR_NETLOC_NETWORK_TYPE__INVALID;
                    break;
                }
                if (dest_node != NULL && edges[i]->dest_node != dest_node) {
                    network_attr->type = MPIR_NETLOC_NETWORK_TYPE__INVALID;
                    break;
                }
            }
        }
    }

    switch (network_attr->type) {
    case MPIR_NETLOC_NETWORK_TYPE__FAT_TREE:
        get_fat_tree_attributes(MPIR_Process.netloc_topology, network_attr);
        break;
    case MPIR_NETLOC_NETWORK_TYPE__TORUS:
        assert(0);
        break;
    default:
        assert(0);
        break;
    }

fn_exit:
    return mpi_errno;

fn_fail:
    goto fn_exit;
}

#endif
