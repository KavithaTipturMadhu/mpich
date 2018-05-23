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
    int mpi_errno = MPI_SUCCESS;

    /* indexed by node id, visited_node_list[i] > 1 indicates that the
     * node i has been visited */
    network_attr->u.fat_tree.node_levels =
        (int *) MPL_malloc(sizeof(int) * topology->num_nodes, MPL_MEM_OTHER);

    for (i = 0; i < topology->num_nodes; i++) {
        network_attr->u.fat_tree.node_levels[i] = -1;
    }
    /* Traversal order never exceeds the total number of nodes anyway */
    traversal_order =
        (netloc_node_t **) MPL_malloc(sizeof(netloc_node_t *) * topology->num_nodes, MPL_MEM_OTHER);

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
                    /*Switch levels start from 1 */
                    network_attr->u.fat_tree.node_levels[edges[j]->dest_node->__uid__] = 1;
                    visited_count++;
                }
            }
        }
    }

    while (visited_count < topology->num_nodes) {
        netloc_node_t *traversed_node = traversal_order[traversal_start++];
        int num_edges;
        netloc_edge_t **edges;
        int depth = network_attr->u.fat_tree.node_levels[traversed_node->__uid__];

        /* find all nodes not visited with an edge from the current node */
        netloc_get_all_edges(topology, traversed_node, &num_edges, &edges);
        for (j = 0; j < num_edges; j++) {
            if (network_attr->u.fat_tree.node_levels[edges[j]->dest_node->__uid__] < 0 ||
                (depth + 1) < network_attr->u.fat_tree.node_levels[edges[j]->dest_node->__uid__]) {
                traversal_order[traversal_end++] = edges[j]->dest_node;
                network_attr->u.fat_tree.node_levels[edges[j]->dest_node->__uid__] = depth + 1;
                visited_count++;
            }
        }
    }

  fn_exit:
    if (traversal_order != NULL) {
        MPL_free(traversal_order);
    }
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}

static int get_network_type(netloc_topology_t netloc_topology,
                            MPIR_Netloc_network_attributes ** network_attr)
{
    int i, j, k, host_node_count;
    netloc_dt_lookup_table_t *host_nodes;
    struct netloc_dt_lookup_table_iterator *hti = NULL;
    int node_query_ret;
    int mpi_errno = MPI_SUCCESS;

    host_nodes =
        (netloc_dt_lookup_table_t *) MPL_malloc(sizeof(netloc_dt_lookup_table_t), MPL_MEM_OTHER);
    node_query_ret = netloc_get_all_host_nodes(netloc_topology, host_nodes);

    if (node_query_ret) {
        (*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__INVALID;
    } else {
        host_node_count = netloc_lookup_table_size(*host_nodes);
        (*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__FAT_TREE;
        hti = netloc_dt_lookup_table_iterator_t_construct(*host_nodes);
        netloc_node_t *start_host_node = NULL;

        while (!netloc_lookup_table_iterator_at_end(hti)) {
            netloc_node_t *node = (netloc_node_t *) netloc_lookup_table_iterator_next_entry(hti);
            int num_edges;
            netloc_edge_t **edges;

            if (node == NULL) {
                break;
            }
            if (start_host_node == NULL) {
                start_host_node = node;
            }
            netloc_get_all_edges(netloc_topology, node, &num_edges, &edges);
            if (num_edges > 0) {
                /* Find the destination node */
                netloc_node_t *dest_node = NULL;
                for (i = 0; i < num_edges; i++) {
                    if (edges[i]->dest_node->node_type != NETLOC_NODE_TYPE_SWITCH ||
                        (dest_node != NULL && edges[i]->dest_node != dest_node)) {
                        (*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__INVALID;
                        break;
                    }
                    dest_node = edges[i]->dest_node;
                }
            }
        }
        /* Check if the tree is a regular fat tree, not a clos network */
        if ((*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__FAT_TREE) {
            netloc_node_t **traversal_order =
                (netloc_node_t **) MPL_malloc(sizeof(netloc_node_t *) * netloc_topology->num_nodes,
                                              MPL_MEM_OTHER);
            int start_index = 0, end_index = 0;
            /*Check for cycles */
            MPIR_Assert(start_host_node != NULL);
            traversal_order[end_index++] = start_host_node;

            /*Traverse edges till you spot a cycle */
            while (end_index < netloc_topology->num_nodes) {
                netloc_node_t *node = traversal_order[start_index++];
                /* Add all neighbors to traversal list */
                netloc_edge_t **edges;
                int num_edges;
                netloc_get_all_edges(netloc_topology, &node, &num_edges, &edges);
                if (num_edges > 2) {
                    (*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__INVALID;
                    break;
                }
                for (j = 0; j < num_edges; j++) {
                    netloc_node_t *neighbor = edges[j]->dest_node;
                    int visited = 0;
                    /*Check if the neighbor has been visited before */
                    for (k = 0; k < end_index; k++) {
                        if (traversal_order[k]->__uid__ == neighbor->__uid__) {
                            int l;
                            netloc_edge_t **neighbor_edges;
                            int neighbor_num_edges;
                            netloc_get_all_edges(netloc_topology, &traversal_order[k],
                                                 &neighbor_num_edges, &neighbor_edges);
                            /*We need to do the following check since links are bidirectional */
                            for (l = 0; l < neighbor_num_edges; l++) {
                                if (neighbor_edges[l] == neighbor) {
                                    break;
                                }
                            }
                            if (l == neighbor_num_edges) {
                                visited = 1;
                                break;
                            }
                        }
                    }

                    if (visited) {
                        (*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__INVALID;
                        break;
                    }
                    traversal_order[end_index++] = neighbor;
                }
                if ((*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__INVALID) {
                    break;
                }
            }
            if (traversal_order != NULL) {
                MPL_free(traversal_order);
            }
        }

        if ((*network_attr)->type == MPIR_NETLOC_NETWORK_TYPE__INVALID) {
            (*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__TORUS;
            /*Check necessary condition for a torus i.e., outdegree of each node is the same */
            int num_edges = -1;
            int first_node = 1;
            hti = netloc_dt_lookup_table_iterator_t_construct(*host_nodes);

            while (!netloc_lookup_table_iterator_at_end(hti)) {
                netloc_node_t *node =
                    (netloc_node_t *) netloc_lookup_table_iterator_next_entry(hti);
                netloc_edge_t **edges;
                int num_edges_per_node;
                netloc_get_all_edges(netloc_topology, node, &num_edges_per_node, &edges);
                if (first_node) {
                    num_edges = num_edges_per_node;
                } else if (num_edges != num_edges_per_node) {
                    (*network_attr)->type = MPIR_NETLOC_NETWORK_TYPE__INVALID;
                    break;
                }
                first_node++;
            }
        }
    }

  fn_exit:
    if (host_nodes != NULL)
        MPL_free(host_nodes);
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}

#undef FUNCNAME
#define FUNCNAME MPIR_Netloc_parse_topology
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIR_Netloc_parse_topology(const char *attribute_file, netloc_topology_t netloc_topology,
                               MPIR_Netloc_network_attributes * network_attr)
{
    int mpi_errno = MPI_SUCCESS;
    int type_query_ret = get_network_type(netloc_topology, &network_attr);
    switch (network_attr->type) {
        case MPIR_NETLOC_NETWORK_TYPE__FAT_TREE:
            get_fat_tree_attributes(MPIR_Process.netloc_topology, network_attr);
            break;
        default:
            MPIR_ERR_CHKANDJUMP(type_query_ret, mpi_errno, MPI_ERR_OTHER,
                                "**netloc_topo_identification");
            break;
    }

  fn_exit:
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}

static void get_physical_address(hwloc_obj_t hwloc_obj, char **device_physical_address)
{
    char *physical_address = NULL;
    int i;
    int mpi_errno = MPI_SUCCESS;

    for (i = 0; i < hwloc_obj->infos_count; i++) {
        /*Check if node guid matches for infiniband networks */
        if (!strcmp(hwloc_obj->infos[i].name, "NodeGUID") ||
            !strcmp(hwloc_obj->infos[i].name, "Address")) {
            physical_address =
                (char *) MPL_malloc(sizeof(hwloc_obj->infos[i].value), MPL_MEM_OTHER);
            strcpy(physical_address, hwloc_obj->infos[i].value);
            break;
        }
    }
  fn_exit:
    *device_physical_address = physical_address;
    return mpi_errno;

  fn_fail:
    goto fn_exit;

}

#undef FUNCNAME
#define FUNCNAME MPIR_get_network_end_point
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIR_get_network_end_point(netloc_topology_t netloc_topology, hwloc_topology_t hwloc_topology,
                               hwloc_cpuset_t bindset, netloc_node_t ** end_point)
{
    hwloc_obj_t io_device = NULL;
    int mpi_errno = MPI_SUCCESS;
    int node_query_ret;
    netloc_dt_lookup_table_t *host_nodes;
    struct netloc_dt_lookup_table_iterator *hti = NULL;
    int i;

    node_query_ret = netloc_get_all_host_nodes(netloc_topology, host_nodes);
    if (node_query_ret) {
        /* No nodes found, which means that topology load failed */
        MPIR_ERR_CHKANDJUMP(node_query_ret, mpi_errno, MPI_ERR_OTHER, "**netloc_topo_load");
    }

    while ((io_device = hwloc_get_next_osdev(hwloc_topology, io_device))
           != NULL) {
        char *physical_address = NULL;
        /* Check if the physical id of the io device matches a node in netloc topology tree */
        get_physical_address(io_device, &physical_address);
        if (physical_address != NULL) {
            /* Find the node in netloc tree with the same physical id */
            hti = netloc_dt_lookup_table_iterator_t_construct(*host_nodes);
            while (!netloc_lookup_table_iterator_at_end(hti)) {
                netloc_node_t *node =
                    (netloc_node_t *) netloc_lookup_table_iterator_next_entry(hti);
                if (node == NULL) {
                    break;
                }
                if (strcmp(physical_address, node->physical_id)) {
                    *end_point = node;
                    break;
                }
            }
            MPL_free(physical_address);
        }
    }
  fn_exit:
    if (host_nodes != NULL) {
        MPL_free(host_nodes);
    }
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}

#undef FUNCNAME
#define FUNCNAME MPIR_get_host_index_in_network
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIR_get_host_index_in_network(netloc_topology_t netloc_topology,
                                   hwloc_topology_t hwloc_topology, hwloc_cpuset_t bindset,
                                   int *index)
{
    hwloc_obj_t io_device = NULL;
    int mpi_errno = MPI_SUCCESS;
    int node_query_ret;
    netloc_dt_lookup_table_t *host_nodes;
    struct netloc_dt_lookup_table_iterator *hti = NULL;
    int i, host_node_count;

    node_query_ret = netloc_get_all_host_nodes(netloc_topology, host_nodes);
    if (node_query_ret) {
        /* No nodes found, which means that topology load failed */
        MPIR_ERR_CHKANDJUMP(node_query_ret, mpi_errno, MPI_ERR_OTHER, "**netloc_topo_load");
    }

    while ((io_device = hwloc_get_next_osdev(hwloc_topology, io_device))
           != NULL) {
        char *physical_address = NULL;

        /* Check if the physical id of the io device matches a node in netloc topology tree */
        get_physical_address(io_device, &physical_address);

        /* Find the node in netloc tree with the same physical id */
        host_node_count = 0;
        hti = netloc_dt_lookup_table_iterator_t_construct(*host_nodes);
        if (physical_address != NULL) {
            while (!netloc_lookup_table_iterator_at_end(hti)) {
                netloc_node_t *node =
                    (netloc_node_t *) netloc_lookup_table_iterator_next_entry(hti);
                if (node == NULL) {
                    break;
                }
                if (strcmp(physical_address, node->physical_id)) {
                    *index = host_node_count;
                    break;
                }
                host_node_count++;
            }
            MPIR_Assert((*index) >= 0);

            MPL_free(physical_address);
        }
    }
  fn_exit:
    if (host_nodes != NULL) {
        MPL_free(host_nodes);
    }
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}

#undef FUNCNAME
#define FUNCNAME MPIR_get_switches_at_level
#undef FCNAME
#define FCNAME MPL_QUOTE(FUNCNAME)
int MPIR_get_switches_at_level(netloc_topology_t topology,
                               MPIR_Netloc_network_attributes attributes, int level,
                               netloc_node_t *** switches_at_level, int *switch_count)
{
    netloc_node_t **switches;
    netloc_dt_lookup_table_t *switch_nodes;
    struct netloc_dt_lookup_table_iterator *hti = NULL;
    int i = 0;
    int node_query_ret = netloc_get_all_switch_nodes(topology, switch_nodes);
    int mpi_errno = MPI_SUCCESS;

    switches = (netloc_node_t **) MPL_malloc(sizeof(netloc_node_t *), MPL_MEM_OTHER);

    hti = netloc_dt_lookup_table_iterator_t_construct(*switch_nodes);
    while (!netloc_lookup_table_iterator_at_end(hti)) {
        netloc_node_t *switch_node = (netloc_node_t *) netloc_lookup_table_iterator_next_entry(hti);
        if (switch_node == NULL) {
            break;
        }
        if (attributes.u.fat_tree.node_levels[switch_node->__uid__] = level) {
            switches =
                (netloc_node_t **) MPL_realloc(switches, sizeof(netloc_node_t *) * i,
                                               MPL_MEM_OTHER);
            switches[i++] = switch_node;
        }
    }

    if (i == 0) {
        MPIR_ERR_CHKANDJUMP(node_query_ret, mpi_errno, MPI_ERR_OTHER, "**netloc_topo_load");
    }

  fn_exit:
    *switches_at_level = switches;
    *switch_count = i;
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}
#endif
