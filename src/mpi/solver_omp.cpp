#include <mpi.h>
#include <algorithm>

#include "common/heat.hpp"

inline void solveBlock(block_t *matrix, row_t ** halo_row, col_t ** halo_col, int nbx, int nby, int bx, int by, ProcessLayout rank2D, HeatConfiguration &conf)
{
	block_t &targetBlock = matrix[bx*nby + by];
	const block_t &centerBlock = matrix[bx*nby + by];
	const block_t &topBlock    = matrix[(bx-1)*nby + by];
	const block_t &leftBlock   = matrix[bx*nby + (by-1)];
	const block_t &rightBlock  = matrix[bx*nby + (by+1)];
	const block_t &bottomBlock = matrix[(bx+1)*nby + by];

	const row_t &halo_top    = (bx == 0)    ? halo_row[top]   [by] : topBlock[BSX-1];
	const row_t &halo_bottom = (bx == nbx-1)? halo_row[bottom][by] : bottomBlock[0];

	double sum = 0.0;
	for (int x = 0; x < BSX; ++x) {

		const row_t &topRow    = (x > 0)     ? centerBlock[x-1] : halo_top;
		const row_t &bottomRow = (x < BSX-1) ? centerBlock[x+1] : halo_bottom;

		const double halo_left_element  = (by == 0)     ? halo_col[left] [bx][x] : leftBlock [x][BSY-1];
		const double halo_right_element = (by == nby-1) ? halo_col[right][bx][x] : rightBlock[x][0];

		for (int y = 0; y < BSY; ++y) {

			double leftElement  = (y > 0)     ? centerBlock[x][y-1] : halo_left_element;
			double rightElement = (y < BSY-1) ? centerBlock[x][y+1] : halo_right_element;

			double value = 0.25 * (topRow[y] + bottomRow[y] + leftElement + rightElement);
			double diff = value - targetBlock[x][y];

			halo_col[left] [bx][x]  = (rank2D.y != 0                      && by == 0     && y == 0)     ? value : halo_col[left] [bx][x] ;
			halo_col[right][bx][x]  = (rank2D.y != conf.processLayout.x-1 && by == nby-1 && y == BSY-1) ? value : halo_col[right][bx][x] ;

			targetBlock[x][y] = value;
		}
	}
}

//A
inline void sendFirstComputeRow(block_t *matrix, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	for (int by = 0; by < nby; ++by) {

		MPI_Send(&matrix[by][0][0], BSY, MPI_DOUBLE, rank2D.getNorth(conf.processLayout), by, MPI_COMM_WORLD);
	}
}

//A
inline void receiveLowerBorder(row_t * halo, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	for (int by = 0; by < nby; ++by) {
		MPI_Recv(&halo[by], BSY, MPI_DOUBLE, rank2D.getSouth(conf.processLayout), by, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

//B
inline void sendLastComputeRow(block_t *matrix, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	for (int by = 0; by < nby; ++by) {
		MPI_Send(&matrix[(nbx-1)*nby + by][BSX-1][0], BSY, MPI_DOUBLE, rank2D.getSouth(conf.processLayout), by, MPI_COMM_WORLD);
	}
}

//B
inline void receiveUpperBorder(row_t * halo, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	for (int by = 0; by < nby; ++by) {
		MPI_Recv(&halo[by], BSY, MPI_DOUBLE, rank2D.getNorth(conf.processLayout), by, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

//C
inline void sendLeftBorder(col_t *halo, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	for (int bx = 0; bx < nbx; ++bx) {
		MPI_Send(&halo[bx], BSX, MPI_DOUBLE, rank2D.getEast(conf.processLayout), bx+nby, MPI_COMM_WORLD);
	}
}

//C
inline void receiveRightBorder(col_t * halo, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	for (int bx = 0; bx < nbx; ++bx) {
		MPI_Recv(&halo[bx], BSX, MPI_DOUBLE, rank2D.getWest(conf.processLayout), bx+nby, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

//D
inline void sendRightBorder(col_t * halo, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	for (int bx = 0; bx < nbx; ++bx) {
		MPI_Send(&halo[bx], BSX, MPI_DOUBLE, rank2D.getWest(conf.processLayout), bx+nby, MPI_COMM_WORLD);
	}
}

//D
inline void receiveLeftBorder(col_t * halo, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	for (int bx = 0; bx < nbx; ++bx) {
		MPI_Recv(&halo[bx], BSX, MPI_DOUBLE, rank2D.getEast(conf.processLayout), bx+nby, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

inline void solveGaussSeidel(block_t *matrix, row_t ** halo_row, col_t ** halo_col, int nbx, int nby, ProcessLayout rank2D, HeatConfiguration &conf)
{
	if(rank2D.x != 0) {
		sendFirstComputeRow(matrix, nbx, nby, rank2D, conf);          //A
		receiveUpperBorder(halo_row[top], nbx, nby, rank2D, conf);    //B
	}
	if(rank2D.y != 0) {
		sendLeftBorder(halo_col[left], nbx, nby, rank2D, conf);   //C
		receiveLeftBorder(halo_col[left], nbx, nby, rank2D, conf);//D

	}

	if(rank2D.x != conf.processLayout.x-1) {
		receiveLowerBorder(halo_row[bottom], nbx, nby, rank2D, conf); //A
	}

	if(rank2D.y != conf.processLayout.y-1) {
		receiveRightBorder(halo_col[right], nbx, nby, rank2D, conf);  //C

	}

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
			solveBlock(matrix, halo_row, halo_col, nbx, nby, bx, by, rank2D, conf);
		}
	}

	#pragma oss taskwait

	if(rank2D.x != conf.processLayout.x-1) {
		sendLastComputeRow(matrix, nbx, nby, rank2D, conf);         //B
	}

	if(rank2D.y != conf.processLayout.y-1) {
		sendRightBorder(halo_col[right], nbx, nby, rank2D, conf);   //D
	}
}

double solve(block_t *matrix, int rowBlocks, int colBlocks, HeatConfiguration &conf,  row_t ** halo_row, col_t ** halo_col, ProcessLayout rank2D )
{
	for (int t = 0; t < conf.timesteps; ++t) {
		solveGaussSeidel(matrix, halo_row, halo_col, rowBlocks, colBlocks, rank2D, conf);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return 0.0;
}

