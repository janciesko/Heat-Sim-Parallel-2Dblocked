#include <mpi.h>
#include <cassert>
#include <iostream>
#include <cstring>
#include <math.h>

#include "common/heat.hpp"

#ifdef _OMPSS_2
#include <nanos6/debug.h>
#endif

#ifdef INTEROPERABILITY
#define DESIRED_THREAD_LEVEL (MPI_THREAD_MULTIPLE+1)
#else
#define DESIRED_THREAD_LEVEL (MPI_THREAD_MULTIPLE)
#endif

#define isSingleProcess false

void generateImage(const HeatConfiguration &conf, int rowBlocks, int colBlocks, int rowBlocksPerRank, int colBlocksPerRank);

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, DESIRED_THREAD_LEVEL, &provided);
	assert(provided == DESIRED_THREAD_LEVEL);
	
	int rank, rank_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	
	HeatConfiguration conf = readConfiguration(argc, argv);

	assert (rank_size ==  conf.processLayout.x * conf.processLayout.y);
	ProcessLayout rank2D {rank / conf.processLayout.y, rank % conf.processLayout.y};

	refineConfiguration(conf, conf.processLayout.x * BSX, conf.processLayout.y * BSY, isSingleProcess);
	if (!rank) printConfiguration(conf);
	
	conf.rowBlocks = conf.rows / BSX;
	conf.colBlocks = conf.cols / BSY;
	int rowBlocks = conf.rowBlocks;
	int colBlocks = conf.colBlocks;
	int rowBlocksPerRank = conf.rowBlocks / conf.processLayout.x;
	int colBlocksPerRank = conf.colBlocks / conf.processLayout.y;
	


	int err = initialize(conf, rowBlocksPerRank, colBlocksPerRank, rank2D);
	assert(!err);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Solve the problem
	double start = get_time();
	solve(conf.matrix, rowBlocksPerRank, colBlocksPerRank, conf, conf.halos_row, conf.halos_col, rank2D);
	double end = get_time();
	
	if (!rank) {
		long totalElements = (long)conf.rows * (long)conf.cols;
		double performance = totalElements * (long)conf.timesteps;
		performance = performance / (end - start);
		performance = performance / 1000000.0;
		
#ifdef _OMPSS_2
		int threads = nanos_get_num_cpus();
#else
		int threads = 1;
#endif
		
		fprintf(stdout, "rows, %d, cols, %d, rows_per_rank, %d, total, %ld, total_per_rank, %ld, bs, %d"
				" ,ranks, %d, threads, %d, timesteps, %d, time, %f, performance, %f\n",
				conf.rows, conf.cols, conf.rows / rank_size, totalElements, totalElements / rank_size,
				BSX, rank_size, threads, conf.timesteps, end - start, performance);
	}
	
	if (conf.generateImage) {
		generateImage(conf, rowBlocks, colBlocks, rowBlocksPerRank, colBlocksPerRank);
	}

	
	err = finalize(conf);
	assert(!err);
	
	MPI_Finalize();
	
	return 0;
}


void generateImage(const HeatConfiguration &conf, int rowBlocks, int colBlocks, int rowBlocksPerRank, int colBlocksPerRank)
{
    int rank, rank_size;
    block_t *auxMatrix = nullptr, *pivotedMatrix = nullptr ;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);

    if (!rank) {
            auxMatrix = (block_t *) malloc(rowBlocks * colBlocks * sizeof(block_t));
            if (auxMatrix == NULL) {
                    fprintf(stderr, "Memory cannot be allocated!\n");
                    exit(1);
            }

            initializeMatrix(conf, auxMatrix, rowBlocks, colBlocks);
    }

    int count = rowBlocksPerRank  * colBlocksPerRank * BSX * BSY;
    MPI_Gather(
            conf.matrix, count, MPI_DOUBLE,
            auxMatrix, count, MPI_DOUBLE,
            0, MPI_COMM_WORLD
    );

	// Pivot blocks for correct visualization (MPI gather returns other data to proc distro)
    if (!rank) {

		pivotedMatrix = (block_t *) malloc(rowBlocks * colBlocks * sizeof(block_t));
		if (pivotedMatrix == NULL) {
			fprintf(stderr, "Memory cannot be allocated!\n");
			exit(1);
		}

		int counter = 0;
		for(int i = 0; i < conf.processLayout.x; i ++) {
			for (int j = 0; j < conf.processLayout.y; j ++)	{
				for (int k = 0; k < rowBlocksPerRank; ++k )	{
					for (int l = 0; l < colBlocksPerRank; ++l){
						memcpy (pivotedMatrix[ i * (colBlocks * rowBlocksPerRank) + j * colBlocksPerRank + k * colBlocks + l],
								auxMatrix[counter++],
								sizeof(block_t));



					}
				}
			}
		}
    }

    if (!rank) {
            int err = writeImage(conf.imageFileName, pivotedMatrix, rowBlocks, colBlocks);
            assert(!err);

            free(auxMatrix);
            free(pivotedMatrix);
    }

}

