//
// Created by Alexander Elder on 2/11/18.
//

#ifndef ALGORITHMS_GRAPH_H
#define ALGORITHMS_GRAPH_H

#include <iostream>
#include <list>
#include <set>
#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#define INF 0x3f3f3f3f

namespace elder {

    /**
     * The following code contains a general directed graph data structure as well as
     * a couple commonly used algorithms. I wrote this primarily to practice C++ and gain a better
     * understanding of combinatorial algorithms. I may add to this periodically.
     */

    struct d_vertex;
    struct d_edge;

    /**
     * Data structure for a directed vertex.
     */
    struct d_vertex {
        long index;
        long weight;
        std::list<d_edge *> edges_in;
        std::list<d_edge *> edges_out;
    };

    /**
     * Data structure for a directed edge.
     */
    struct d_edge {
        long length;
        d_vertex *origin;
        d_vertex *destination;
    };

    /**
     * Data structure for output from shortest paths algorithms.
     */
    struct shortest_paths {
        bool negative_cost_cycle;
        d_vertex *source;
        std::vector<long> path_map;
        std::vector<long> distances;

        void print_distances() {
            if (negative_cost_cycle) {
                std::cout << "Negative cost cycle detected.\n";
            } else {
                std::cout << "Shortest path distances from " << source->index << std::endl;
                for (int i = 0; i < distances.size(); ++i) {
                    std::cout << i << " : " << distances[i] << std::endl;
                }
            }
        };
    };


    /**
     * Class for directed graphs.
     */
    class d_graph {
        long V; // number of vertices
        long E; // number of edges
        bool negative_edges;
        bool enforce_positive_edges;

        /**
         * Add an edge to the graph.
         * @param u pointer to origin vertex
         * @param v pointer to destination vertex
         * @param w edge weight
         */
        void insert_edge(d_vertex *u, d_vertex *v, long w) {
            if (w < 0 && enforce_positive_edges) {
                std::cout << "WARNING: This graph has positive vertices enforced. Failed to add edge from "
                          << u->index << " to " << v->index << " of length " << w << std::endl;
                return;
            } else {
                if (!negative_edges && w < 0) {
                    negative_edges = true;
                }
                d_edge *e = new d_edge;
                e->origin = u;
                e->destination = v;
                e->length = w;

                u->edges_out.push_back(e);
                v->edges_in.push_back(e);

                E++;// increment the number of edges
            }
        }

    public:
        d_vertex *vertices;// dynamic array of vertex pointers
        /**
         * Constructor.
         * @param V : Number of vertices.
         */
        explicit d_graph(long V) : V(V) {
            enforce_positive_edges = false;
            E = 0;
            vertices = new d_vertex[V];
            for (long i = 0; i < V; ++i) {
                vertices[i].index = i;
            }
        }

        /**
         * Destructor. Deletes vertices and edges.
         */
        virtual ~d_graph() {
            for (int i = 0; i < V; ++i) {
                std::list<d_edge *>::iterator it;
                d_vertex *v = &vertices[i];
                for (it = v->edges_out.begin(); it != v->edges_out.end(); ++it) {
                    delete (*it);
                }
            }
            delete[] vertices;
        }

        /**
         * Public insert edge function. Avoids memory management.
         * @param u : index of origin
         * @param v : index of destination
         * @param w : edge weight
         */
        void insert_edge(long u, long v, long w) {
            insert_edge(&vertices[u], &vertices[v], w);
        }

        /**
         * Efficient implementation of Dijkstra's algorithm using priority queue data structure. Based off
         * geeksforgeeks.com's implementation. Finds all shortest paths to input vertex if the graph has only
         * positive edge weights.
         * @param source_index
         * @return
         */
        shortest_paths Dijkstra_shortest_path(long source_index) {
            // create priority que to allow O(nlogm) lookup of new vectors
            //    priority_queue< pair<long, long>, vector< pair<long, long> > > p;
            std::priority_queue<std::pair<long, d_vertex *>, std::vector<std::pair<long, d_vertex *> > > p;

            // initialize distances vector with all distances infinite
            std::vector<long> dist(V, INF);
            std::vector<long> parent(V, -1);

            // insert source into priority queue and initialize its distance to itself as zero
            d_vertex *source = &vertices[source_index];
            p.push(std::make_pair(0, source));
            dist[source_index] = 0;

            // greedy search
            while (!p.empty()) {
                // the minimum distance is stored as the first member of the pair while
                // the index is stored as the second member in the priority queue
                d_vertex *u = p.top().second;
                p.pop();//

                // iterate over the edges of the minimum index
                std::list<d_edge *>::iterator i;
                for (i = u->edges_out.begin(); i != u->edges_out.end(); ++i) {

                    // label and weight of the current adjacent vertex
                    d_vertex *v = (*i)->destination;
                    long weight = (*i)->length;

                    // update the distance to v if there is a shortened path to v through u
                    if (dist[v->index] > dist[u->index] + weight) {
                        dist[v->index] = dist[u->index] + weight;
                        p.push(std::make_pair(dist[v->index], v));
                        parent[v->index] = u->index;
                    }
                }
            }
            // return the distance from source of each node, associated by index
            shortest_paths sp;
            sp.source = source;
            sp.distances = dist;
            sp.negative_cost_cycle = false;
            sp.path_map = parent;

            return sp;
        }

        /**
         * Bellman Ford routing algorithm for all shortest paths to a vertex
         * given a graph with positive and negative edge weights.
         * @param source_index
         * @return
         */
        shortest_paths Bellman_Ford_shortest_paths(long source_index){
            long i,j, s = source_index, n = V+1; // encode memory for extra iteration here
            bool negative_cost_cycle = false;

            // initialize shortest paths data structure
            shortest_paths sp;
            sp.source = &vertices[source_index];
            sp.negative_cost_cycle = false;

            // initialize 2d array to store intermediate values. For first iteration iterate source as 0 and all other vectors
            // as INF
            std::vector<long> A_prev (V);
            std::vector<long> A_cur  (V);
            for(j=0; j<V; ++j){
                A_prev[j] = INF;
            }
            A_prev[s] = 0;

            // initialize a 2d array to store shortest paths
            std::vector<long> B_prev (V);
            std::vector<long> B_cur  (V);
            for(j=0; j<V; ++j){
                B_prev[j] = 0;
            }

            // main loop
            for(i=1; i<n-1; ++i){

                bool constant_iteration = true;
                for(j=0; j<V; ++j){
                    // the static case happens when there is no shorter path found
                    long static_case = A_prev[j];
                    d_vertex* v = &vertices[j];
                    std::list<d_edge*>::iterator w;
                    long test_case = INF;
                    long test_index= 0;
                    // search for shorter path
                    for(w = v->edges_in.begin(); w != v->edges_in.end(); ++w){
                        long prev_index = (*w)->origin->index;
                        if( A_prev[prev_index] + (*w)->length < test_case){
                            test_index= (*w)->origin->index;
                            test_case = A_prev[test_index] + (*w)->length;
                        }
                    }

                    // take min of test_case and degenerate_case and update the shortest paths list
                    if(test_case < static_case){
                        A_cur[j] = test_case;
                        B_cur[j] = test_index;
                    }
                    else{
                        A_cur[j] = static_case;
                        B_cur[j] = B_prev[j];
                    }
                    // halt check
                    if(A_cur[j] != A_prev[j]){
                        constant_iteration = false;
                    }
                }
                // check early halting condition
                if(constant_iteration){
                    break;
                }
                else{
                    // swap arrays for memory optimization
                    std::vector<long> temp = A_cur;
                    A_cur = A_prev;
                    A_prev = temp;

                    temp = B_cur;
                    B_cur = B_prev;
                    B_prev = temp;
                }
            }

            // if the loop runs to completion, check for negative cost cycles with an extra iteration
            if(i == n-1){
                for(j=0; j<V; ++j){
                    // store static case
                    long static_case = A_prev[j];
                    d_vertex* v = &vertices[j];
                    std::list<d_edge*>::iterator w;
                    long test_case = INF;
                    // search for shorter path
                    for(w = v->edges_in.begin(); w != v->edges_in.end(); ++w){
                        if( A_prev[(*w)->origin->index] + (*w)->length < test_case){
                            test_case = A_prev[(*w)->origin->index] + (*w)->length;
                        }
                    }
                    A_cur[j] = std::min(test_case, static_case);
                    //halt check
                    if(A_cur[j] != A_prev[j]){
                        negative_cost_cycle = true;
                    }
                }
                if(negative_cost_cycle){
                    sp.negative_cost_cycle = true;
                }
            }
            if(!negative_cost_cycle){
                // if no negative cost cycles, assign the distances and path_map values in the return data structure
                sp.distances = A_cur;
                sp.path_map  = B_cur;
            }

            // return shortest path data structure
            return sp;
        }

        /**
         * Johnson's algorithm for finding all shortest paths given a graph with positive and negative edge weights.
         * Complexity is O(V^2logV + VE) where V is the number of vertices and E is the number of edges.
         * @param verbose
         * @return
         */
        std::vector<shortest_paths> all_shortest_paths(bool verbose) {
            if(verbose){
                std::cout << "Computing all source shortest paths using Johnson's Algorithm.\n";
            }
            // initialize vector of shortest path data structures-for each vertex, log the shortest path to all other vertices
            std::vector<shortest_paths> sps (V);
            // Implementation of Johnson's all pairs shortest path algorithm
            // step 1: form expanded graph G' by adding a new vertex s that has edge length (s,v) = 0 for all v in G'
            d_graph G_prime (V+1);
            for(int i=0; i<V; ++i){// copy all of G's connections
                d_vertex *v1 = &vertices[i];
                for(auto edge:v1->edges_out){
                    G_prime.insert_edge(i, edge->destination->index, edge->length);
                }
            }
            for(int i=0; i<V; ++i){// add edges to new point
                G_prime.insert_edge(V, i, 0);
            }

            // step 2: run Bellman_ford shortest paths algorithm with new node as source (index V)
            if(verbose){
                std::cout << "Running Bellman_Ford_shortest_path on modified graph.\n";
            }
            shortest_paths sp_prime = G_prime.Bellman_Ford_shortest_paths(V);
            if(sp_prime.negative_cost_cycle){
                std::cout << "WARNING! all_shortest_paths: negative cost structure detected, returning empty vector.\n";
                for(int i=0; i<sps.size(); ++i){
                    sps[i].negative_cost_cycle = true;
                }
                return sps;
            }

            // step 3: compute new edge weights
            if(verbose){
                std::cout << "Reweighting edges.\n";
            }
            for(int i=0; i<V; ++i){
                d_vertex *v1 = &vertices[i];

                // update edge weights for edges originating at vertex i
                for(auto e1: v1->edges_out){
                    long p_u = sp_prime.distances[i];
                    long p_v = sp_prime.distances[e1->destination->index];
                    e1->length = e1->length + p_u - p_v;
                }
            }

            if(verbose){
                std::cout << "Computing shortest paths with Dijkstra's algorithm.\n";
            }
            // step 4: compute all edge weights using Dijkstra's algorithm
            for(int i=0; i<V; ++i){
                sps[i] = Dijkstra_shortest_path(i);

                if(verbose){
                    if(i == V/4){
                        std::cout << "    25% complete\n";
                    }
                    if(i == V/2){
                        std::cout << "    50% complete\n";
                    }
                    if(i == 3*V/4){
                        std::cout << "    75% complete\n";
                    }
                }

            }
            if(verbose){
                std::cout << "Shortest path computation complete. Reverting to original weights.\n";
            }
            // step 5: return shortest path distance
            for(int i=0; i<V; ++i){
                d_vertex *v1 = &vertices[i];

                // update edge weights for edges originating at vertex i
                for(auto e1: v1->edges_out){
                    long p_u = sp_prime.distances[i];
                    long p_v = sp_prime.distances[e1->destination->index];
                    e1->length = e1->length - p_u + p_v;
                }
            }
            // transform shortest paths distances
            for(int i=0; i<V; ++i){
                for(int j=0; j<V; ++j){
                    sps[j].distances[i] = sps[j].distances[i] - sp_prime.distances[j] + sp_prime.distances[i];
                }
            }

            return sps;

        }
    };

}

#endif //ALGORITHMS_GRAPH_H
