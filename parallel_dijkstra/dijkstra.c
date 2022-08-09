#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <omp.h>

static unsigned long int next = 1;

int my_rand(void) {
        return ((next = next * 1103515245 + 12345) % ((u_long) RAND_MAX + 1));
}

void my_srand(unsigned int seed) {
        next = seed;
}

struct Graph {
        int nNodes;
        int *nEdges;
        int **edges;
        int **w;
};

struct Graph *createRandomGraph(int nNodes, int nEdges, int seed) {
        my_srand(seed);

        struct Graph *graph = (struct Graph *) malloc(sizeof(struct Graph));
        graph->nNodes = nNodes;
        graph->nEdges = (int *) malloc(sizeof(int) * nNodes);
        graph->edges = (int **) malloc(sizeof(int *) * nNodes);
        graph->w = (int **) malloc(sizeof(int *) * nNodes);

        //execution time
        double i_exec_t, f_exec_t, exec_time;

        // Initial operation time
        i_exec_t = omp_get_wtime();

	int k, v;
        //#pragma omp parallel for firstprivate(nNodes)
        for (v = 0; v < nNodes; v++) {
                        graph->edges[v] = (int *) malloc(sizeof(int) * nNodes);
                        graph->w[v] = (int *) malloc(sizeof(int) * nNodes);
                        graph->nEdges[v] = 0;
        }

        // Final operation time
        f_exec_t = omp_get_wtime();
        exec_time = f_exec_t - i_exec_t;

        //print execution time and number of threads
        printf("\nCreate Random Graph operation time, for 1: %lf seconds\n", exec_time);

        // Initial operation time
        i_exec_t = omp_get_wtime();

        int source = 0;
        for (source = 0; source < nNodes; source++) {
                int nArestasVertice = (double) nEdges / nNodes * (0.5 + my_rand() / (double) RAND_MAX);
                //#pragma omp parallel for firstprivate(source) private(k) shared(nNodes)
                for (k = nArestasVertice; k >= 0; k--) {
                        int dest, w;
                        //#pragma omp critical
                        //{
                                dest = my_rand() % nNodes;
                                w = 1 + (my_rand() % 10);
                        //}
                        graph->edges[source][graph->nEdges[source]] = dest;
                        graph->w[source][graph->nEdges[source]++] = w;
                }
        }

        // Final operation time
        f_exec_t = omp_get_wtime();
        exec_time = f_exec_t - i_exec_t;

        //print execution time and number of threads
        printf("Create Random Graph operation time, for 2: %lf seconds\n\n", exec_time);

        return graph;
}

int *dijkstra(struct Graph *graph, int source) {
        int nNodes = graph->nNodes;
        int *visited = (int *) malloc(sizeof(int) * nNodes);
        int *distances = (int *) malloc(sizeof(int) * nNodes);
        int k, v;

        for (v = 0; v < nNodes; v++) {
                distances[v] = INT_MAX;
                visited[v] = 0;
        }
        distances[source] = 0;
        visited[source] = 1;

        //execution time
        double i_exec_t, f_exec_t, exec_time;

        // Initial operation time
        i_exec_t = omp_get_wtime();

        for (k = 0; k < graph->nEdges[source]; k++)
                distances[graph->edges[source][k]] = graph->w[source][k];

        // Final operation time
        f_exec_t = omp_get_wtime();
        exec_time = f_exec_t - i_exec_t;

        //print execution time and number of threads
        printf("Djikstra operation time, for 1: %lf seconds\n", exec_time);

        // Initial operation time
        i_exec_t = omp_get_wtime();

        for (v = 1; v < nNodes; v++) {
                int min = 0;
                int minValue = INT_MAX;
                #pragma omp parallel for shared(minValue,min)
                for (k = 0; k < nNodes; k++)
                        if (visited[k] == 0 && distances[k] < minValue) {
                                minValue = distances[k];
                                min = k;
                        }

		visited[min] = 1;

                #pragma omp parallel for
                for (k = 0; k < graph->nEdges[min]; k++) {
                        int dest = graph->edges[min][k];
                        if (distances[dest] > distances[min] + graph->w[min][k])
                                distances[dest] = distances[min] + graph->w[min][k];
                }
        }

        // Final operation time
        f_exec_t = omp_get_wtime();
        exec_time = f_exec_t - i_exec_t;

        //print execution time and number of threads
        printf("Djikstra operation time, for 2: %lf seconds\n\n", exec_time);

        free(visited);

        return distances;
}

int main(int argc, char ** argv) {
        int nNodes;
        int nEdges;
        int seed;

        //execution time
        double i_exec_t, f_exec_t, exec_time;

        if (argc == 4) {
                nNodes = atoi(argv[1]);
                nEdges = atoi(argv[2]);
                seed = atoi(argv[3]);
        }
        else if (argc == 5) {
                nNodes = atoi(argv[1]);
                nEdges = atoi(argv[2]);
                seed = atoi(argv[3]);
                omp_set_num_threads(atoi(argv[4]));
        }
         else {
                int aux = fscanf(stdin, "%d %d %d", &nNodes, &nEdges, &seed);
	}

        // Initial operation time
        i_exec_t = omp_get_wtime();

        nEdges = nNodes * nEdges;

        struct Graph *graph = createRandomGraph(nNodes, nEdges, seed);

        int *dist = dijkstra(graph, 0);

        double mean = 0;
        int v;
        for (v = 0; v < graph->nNodes; v++)
                mean += dist[v];

        fprintf(stdout, "Result: %.2f\n", mean / nNodes);

        // Final operation time
        f_exec_t = omp_get_wtime();
        exec_time = f_exec_t - i_exec_t;

        //print execution time and number of threads
        printf("\nOperation time total: %lf seconds\n\n", exec_time);
        printf("Threads: %d\n\n", omp_get_max_threads());
        return 0;
}