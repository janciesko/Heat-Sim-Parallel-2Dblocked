#include <cassert>
#include <iostream>

#include "common/heat.hpp"

#ifdef _OMPSS_2
#include <nanos6/debug.h>
#endif

#define isSingleProcess true

int main(int argc, char **argv)
{
	HeatConfiguration conf = readConfiguration(argc, argv);

	refineConfiguration(conf, BSX, BSY, isSingleProcess);
	printConfiguration(conf);
	
	conf.rowBlocks = conf.rows / BSX;
	conf.colBlocks = conf.cols / BSY;
	int rowBlocks = conf.rowBlocks;
	int colBlocks = conf.colBlocks;
	
	int err = initialize(conf, colBlocks, rowBlocks);
	assert(!err);

	
	// Solve the problem
	double start = get_time();
	double residual = solve(conf.matrix, rowBlocks, colBlocks, conf, conf.halos_row, conf.halos_col);
	double end = get_time();
	
	long totalElements = (long)conf.rows * (long)conf.cols;
	double performance = totalElements * (long)conf.timesteps;
	performance = performance / (end - start);
	performance = performance / 1000000.0;
	
#ifdef _OMPSS_2
	int threads = nanos_get_num_cpus();
#else
	int threads = 1;
#endif
	
	fprintf(stdout, "rows, %d, cols, %d, total, %ld, bs, %d, threads, %d, timesteps, %d, time, %f, performance, %f\n",
		conf.rows, conf.cols, totalElements, BSX, threads, conf.timesteps, end - start, performance);
	
	if (conf.generateImage) {
		err = writeImage(conf.imageFileName, conf.matrix, rowBlocks, colBlocks);
		assert(!err);
	}
	
	err = finalize(conf);
	assert(!err);
	
	return 0;
}
