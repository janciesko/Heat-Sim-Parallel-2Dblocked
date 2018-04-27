#ifndef MATRIX_HPP
#define MATRIX_HPP

#ifndef BSX
#define BSX 1024
#endif

#ifndef BSY
#define BSY BSX
#endif

#define left 0
#define right 1
#define top 0
#define bottom 1

// Definition of types
typedef double row_t[BSY];
typedef double col_t[BSX];
typedef row_t block_t[BSX];

// Useful functions for matrices
template <typename Func>
inline void traverseRow(block_t *matrix, int numColBlocks, int row, int startCol, int endCol, Func func)
{
	int rowBlock = row / BSX;
	for (int col = startCol; col < endCol; ++col) {
		int colBlock = col / BSY;
		block_t &block = matrix[rowBlock * numColBlocks + colBlock];
		func(row, col, block[row % BSX][col % BSY]);
	}
}

template <typename Func>
inline void traverseByRows(block_t *matrix, int rowBlocks, int colBlocks, Func func)
{
	int numRows = rowBlocks * BSX;
	int numCols = colBlocks * BSY;

	for (int x = 0; x < numRows; ++x) {
		traverseRow(matrix, colBlocks, x, 0, numCols, func);
	}
}


template <typename Func>
inline void copyByRows(block_t *source, block_t * target, int rowBlocks, int colBlocks, Func func)
{
	int numRows = rowBlocks * BSX;
		for (int x = 0; x < numRows; ++x) {
		func(target, source);
	}
}


template <typename Func>
inline void traverseRowHalo(col_t * row,int r, int startCol, int endCol, Func func)
{;
	for (int col = startCol; col < endCol; ++col) {
		func(r, col, row[col / BSY][col % BSY]);
	}
}

template <typename Func>
inline void traverseColHalo(col_t * column,int col, int startRow, int endRow, Func func)
{;
	for (int row = startRow; row < endRow; ++row) {
		func(row, col, column[row / BSX][row % BSX]);
	}
}

#endif // MATRIX_HPP
