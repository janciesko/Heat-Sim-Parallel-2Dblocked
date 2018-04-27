#ifndef HEAT_HPP
#define HEAT_HPP

#include <string>

#include "common/matrix.hpp"


template<class T = int>
inline T round(T a, T b)
{
	return ((a + b - 1) / b) * b;
}

struct HeatSource {
	float row;
	float col;
	float range;
	float temperature;
	
	HeatSource() :
		row(0.0),
		col(0.0),
		range(0.0),
		temperature(0.0)
	{
	}
};

struct ProcessLayout
{
	int x;
	int y;

	int getEast (ProcessLayout & pl){return x * pl.y + y - 1;}
	int getWest (ProcessLayout & pl){return x * pl.y + y + 1;}
	int getNorth(ProcessLayout & pl){return (x-1) * pl.y + y;}
	int getSouth(ProcessLayout & pl){return (x+1) * pl.y + y;}

	ProcessLayout(int a, int b) {x = a; y = b;}
	ProcessLayout(const ProcessLayout &pl) {x = pl.x; y = pl.y;}
};

struct HeatConfiguration {
	int timesteps;
	int rows;
	int cols;
	int rowBlocks;
	int colBlocks;
	block_t *matrix;
	col_t *halos_row[2];
	row_t *halos_col[2];
	bool isSingleProcess;
	int numHeatSources;
	HeatSource *heatSources;
	std::string confFileName;
	std::string imageFileName;
	bool generateImage;
	ProcessLayout processLayout;
	
	HeatConfiguration() :
		timesteps(0),
		rows(0),
		cols(0),
		rowBlocks(0),
		colBlocks(0),
		matrix(nullptr),
		halos_row{nullptr},
		halos_col{nullptr},
		isSingleProcess (true),
		numHeatSources(0),
		heatSources(nullptr),
		confFileName("heat.conf"),
		imageFileName("heat.ppm"),
		generateImage(false),
		processLayout{1,1}
	{
	}
};

int initialize(HeatConfiguration &conf, int rowBlocks, int colBlocks , ProcessLayout r = ProcessLayout (0,0));
int finalize(HeatConfiguration &conf);
int writeImage(std::string fileName, block_t *matrix, int rowBlocks, int colBlocks);
HeatConfiguration readConfiguration(int argc, char **argv);
void refineConfiguration(HeatConfiguration &conf, int rowValue, int colValue, bool isSingleProcess);
void printConfiguration(const HeatConfiguration &conf);
void initializeMatrix(const HeatConfiguration &conf, block_t *matrix, int rowBlocks, int colBlocks, ProcessLayout r = ProcessLayout (0,0));
void initializeHalos(const HeatConfiguration &conf, block_t *matrix, int rowBlocks, int colBlocks, ProcessLayout r = ProcessLayout (0,0));
double get_time();
double solve(block_t *matrix, int rowBlocks, int colBlocks, HeatConfiguration &conf,  row_t ** halos_row = nullptr, col_t ** halos_col = nullptr,  ProcessLayout r = ProcessLayout (0,0));

#endif // HEAT_HPP
