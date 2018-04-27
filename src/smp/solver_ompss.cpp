#include <algorithm>
#include <iostream>

#include "common/heat.hpp"


inline double solveBlock(block_t *matrix, row_t ** halo_row, col_t ** halo_col, int nbx, int nby, int bx, int by)
{
	block_t &targetBlock = matrix[bx*nby + by];
	const block_t &centerBlock = matrix[bx*nby + by];
	const block_t &topBlock    = matrix[(bx-1)*nby + by];
	const block_t &leftBlock   = matrix[bx*nby + (by-1)];
	const block_t &rightBlock  = matrix[bx*nby + (by+1)];
	const block_t &bottomBlock = matrix[(bx+1)*nby + by];

	const row_t &haloTop    = (bx == 0)    ? halo_row[top]   [by] : topBlock[BSX-1];
	const row_t &haloBottom = (bx == nbx-1)? halo_row[bottom][by] : bottomBlock[0];

	double sum = 0.0;
	for (int x = 0; x < BSX; ++x) {

		const row_t &topRow    = (x > 0)     ? centerBlock[x-1] : haloTop;
		const row_t &bottomRow = (x < BSX-1) ? centerBlock[x+1] : haloBottom;

		const double halo_left_element  = (by == 0)     ? halo_col[left] [bx][x] : leftBlock [x][BSY-1];
		const double halo_right_element = (by == nby-1) ? halo_col[right][bx][x] : rightBlock[x][0];

		for (int y = 0; y < BSY; ++y) {

			double leftElement  = (y > 0)     ? centerBlock[x][y-1] : halo_left_element;
			double rightElement = (y < BSY-1) ? centerBlock[x][y+1] : halo_right_element;

			double value = 0.25 * (topRow[y] + bottomRow[y] + leftElement + rightElement);
			double diff = value - targetBlock[x][y];
			sum += diff * diff;

			targetBlock[x][y] = value;
		}
	}

	return sum;
}

inline void gaussSeidelSolver(block_t * matrix, row_t ** halo_row, col_t ** halo_col, int nbx, int nby)
{
	double unew, diff, sum = 0.0;

	int bsX = BSX;
	int bsY = BSY;

	for (int bx = 0; bx < nbx; ++bx) {
		for (int by = 0; by < nby; ++by) {


			block_t * topBlock    = (bx == 0)     ? nullptr : & matrix[(bx-1)*nby + by];
			block_t * bottomBlock = (bx == nbx-1) ? nullptr : & matrix[(bx+1)*nby + by];
			block_t * leftBlock   = (by == 0)     ? nullptr : & matrix[bx*nby + (by-1)];
			block_t * rightBlock  = (by == nby-1) ? nullptr : & matrix[bx*nby + (by+1)];

			#pragma oss task label(gauss seidel) \
				in ([1]topBlock)              \
				in ([1]leftBlock)             \
				in ([1]rightBlock)            \
				in ([1]bottomBlock)           \
				inout(([nbx][nby]matrix)[bx][by])
			solveBlock(matrix, halo_row, halo_col, nbx, nby, bx, by);
		}
	}
}

double solve(block_t *matrix, int rowBlocks, int colBlocks, HeatConfiguration &conf,  row_t ** halos_row, col_t ** halos_col, ProcessLayout rank2D)
{
	double residual = 0.0;

	for (int t = 0; t < conf.timesteps; ++t) {
		gaussSeidelSolver(matrix, halos_row, halos_col, rowBlocks, colBlocks);
	}
	#pragma oss taskwait

	return residual;
}
