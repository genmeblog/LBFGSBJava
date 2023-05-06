package org.generateme.lbfgsb;

public final class Matrix {
	public double[] mat;
	public int rows;
	public int cols;

	public Matrix(double[] mat, int rows, int cols) {
		this.mat = mat;
		this.rows = rows;
		this.cols = cols;
	}

	public Matrix(int rows, int cols) {
		this(new double[rows * cols], rows, cols);
	}

	public static Matrix resize(Matrix m, int rows, int cols) {
		if (m == null) {
			return new Matrix(rows, cols);
		} else {
			if (m.rows == rows && m.cols == cols) {
				return m;
			} else {
				return new Matrix(Vector.resize(m.mat, rows * cols), rows, cols);
			}
		}
	}

	public void setAll(double val) {
		Vector.setAll(mat, val);
	}
	
	public void set(int row, int col, double val) {
		mat[col * rows + row] = val;
	}

	public void setCol(int col, double[] v) {
		System.arraycopy(v, 0, mat, col * rows, rows);
	}

	public double[] getCol(int col) {
		double[] res = new double[rows];
		System.arraycopy(mat, col * rows, res, 0, rows);
		return res;
	}

	public void getCol(int col, double[] res) {
		System.arraycopy(mat, col * rows, res, 0, rows);
	}

	public double get(int row, int col) {
		return mat[col * rows + row];
	}

	public int index(int row, int col) {
		return col * rows + row;
	}

	public void setDiag(double val) {
		for (int i = 0; i < Math.min(rows, cols); i++) {
			set(i, i, val);
		}
	}

	public double colSquaredNorm(int col) {
		double res = 0.0;
		int idx = col * rows;
		for (int i = 0; i < rows; i++) {
			res += mat[idx + i] * mat[idx + i];
		}
		return res;
	}

	public double colNorm(int col) {
		return Math.sqrt(colSquaredNorm(col));
	}

	public double colDot(int col, double[] v) {
		double res = 0.0;
		int idx = col * rows;
		for (int i = 0; i < rows; i++) {
			res += mat[idx + i] * v[i];
		}

		return res;
	}

	public void mulv(double[] v, double[] res) {
		for (int row = 0; row < rows; row++) {
			double dot = 0.0;
			for (int col = 0; col < cols; col++) {
				dot += get(row, col) * v[col];
			}
			res[row] = dot;
		}
	}

	public void muls(double v) {
		for (int i = 0; i < mat.length; i++) {
			mat[i] *= v;
		}
	}
}
