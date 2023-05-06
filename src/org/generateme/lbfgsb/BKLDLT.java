// https://lbfgspp.statr.me/doc/BKLDLT_8h_source.html

package org.generateme.lbfgsb;

import static org.generateme.lbfgsb.Debug.*;

import java.util.ArrayList;

//Bunch-Kaufman LDLT decomposition
public final class BKLDLT {

	static final double alpha = (1.0 + Math.sqrt(17.0)) / 8.0;
	
	public enum Info {
		SUCCESSFUL, NOT_COMPUTED, NUMERICAL_ISSUE
	}

	// Compressed permutations
	class Pair {
		int a, b;

		public Pair(int a, int b) {
			this.a = a;
			this.b = b;
		}

		public String toString() {
			return "[" + a + ", " + b + "]";
		}
	}

	// Index reference, mutable
	class Index {
		int v;

		public Index(int v) {
			this.v = v;
		}
	}

	public static final void swap(double[] a, int i, int j) {
		double t = a[i];
		a[i] = a[j];
		a[j] = t;
	}

	public static final void swap_ranges(double a[], int start, int end, int target) {
		for (int i = start, j = target; i < end; i++, j++) {
			swap(a, i, j);
		}
	}

	public boolean computed = false;
	public Info info = Info.NOT_COMPUTED;

	public int n = 0;
	public double[] data; // lower triangular matrix, column-wise
	public int[] colptr; // indices of columns
	public int[] perm; // permutations
	public ArrayList<Pair> permc; // compressed permutations

	// index of value
	public int index(int i, int j) {
		return colptr[j] + (i - j);
	}

	// value from matrix
	public double coeff(int i, int j) {
		return data[colptr[j] + (i - j)];
	}

	// diagonal value
	public double diag_coeff(int i) {
		return data[colptr[i]];
	}

	// find and store column pointers
	private void compute_pointer() {
		colptr = new int[n];

		int head = 0;
		for (int i = 0; i < n; i++) {
			colptr[i] = head;
			head += n - i;
		}
	}

	// copy matrix data
	private void copy_data(Matrix mat) {
		for (int j = 0; j < n; j++) {
			int begin = mat.index(j, j);
			int len = n - j;
			System.arraycopy(mat.mat, begin, data, colptr[j], len);
		}
	}

	private void compress_permutation() {
		permc = new ArrayList<Pair>(n);

		for (int i = 0; i < n; i++) {
			int idx = perm[i] >= 0 ? perm[i] : -perm[i] - 1;
			if (idx != i) {
				permc.add(new Pair(i, idx));
			}
		}
	}

	public double find_lambda(int k, Index r) {
		int head = colptr[k];
		int end = colptr[k + 1];
		r.v = k + 1;
		double lambda = Math.abs(data[head + 1]);

		for (int ptr = head + 2; ptr < end; ptr++) {
			double abs_elem = Math.abs(data[ptr]);
			if (lambda < abs_elem) {
				lambda = abs_elem;
				r.v = k + (ptr - head);
			}
		}

		return lambda;
	}

	public double find_sigma(int k, int r, Index p) {
		double sigma = -1.0;

		if (r < n - 1) {
			sigma = find_lambda(r, p);
		}

		for (int j = k; j < r; j++) {
			double abs_elem = Math.abs(coeff(r, j));
			if (sigma < abs_elem) {
				sigma = abs_elem;
				p.v = j;
			}
		}

		return sigma;
	}

	public void pivoting_1x1(int k, int r) {
		if (k != r) {
			swap(data, colptr[k], colptr[r]);
			swap_ranges(data, index(r + 1, k), colptr[k + 1], index(r + 1, r));
			int src = index(k + 1, k);
			for (int j = k + 1; j < r; j++, src++) {
				swap(data, src, index(r, j));
			}
		}
		perm[k] = r;
	}

	public void pivoting_2x2(int k, int r, int p) {
		pivoting_1x1(k, p);
		pivoting_1x1(k + 1, r);
		swap(data, index(k + 1, k), index(r, k));
		perm[k] = -perm[k] - 1;
		perm[k + 1] = -perm[k + 1] - 1;
	}

	public void interchange_rows(int r1, int r2, int c1, int c2) {
		if (r1 != r2) {
			for (int j = c1; j <= c2; j++) {
				swap(data, index(r1, j), index(r2, j));
			}
		}
	}

	public boolean permutate_mat(int k) {
		Index r = new Index(k);
		Index p = new Index(k);

		double lambda = find_lambda(k, r);

		if (lambda > 0) {

			double abs_akk = Math.abs(diag_coeff(k));
			double alambda = alpha * lambda;

			if (abs_akk < alambda) {

				double sigma = find_sigma(k, r.v, p);

				if (sigma * abs_akk < alambda * lambda) {

					if (abs_akk >= alpha * sigma) {
						pivoting_1x1(k, r.v);
						interchange_rows(k, r.v, 0, k - 1);
						return true;
					} else {
						p = new Index(k);
						pivoting_2x2(k, r.v, p.v);
						interchange_rows(k, p.v, 0, k - 1);
						interchange_rows(k + 1, r.v, 0, k - 1);
						return false;
					}

				}
			}
		}

		return true;
	}

	public Info gaussian_elimination_1x1(int k) {

		double akk = diag_coeff(k);

		if (akk == 0.0) {
			return Info.NUMERICAL_ISSUE;
		}

		data[colptr[k]] = 1.0 / akk;

		int lptr = colptr[k] + 1;
		int ldim = n - k - 1;

//		MapVec l(lptr, ldim);
//    for (Index j = 0; j < ldim; j++)
//    {
//        MapVec(col_pointer(j + k + 1), ldim - j).noalias() -= (lptr[j] / akk) * l.tail(ldim - j);
//    }
		for (int j = 0; j < ldim; j++) {
			int col = colptr[j + k + 1];
			double v1 = data[lptr + j] / akk;
			for (int i = 0; i < ldim - j; i++) {
				double v = v1 * data[lptr + j + i]; // tail, idx of starting column = ldim-(ldim-j) = j
				data[col + i] -= v;
			}
		}

		for (int i = 0; i < ldim; i++) {
			data[lptr + i] /= akk;
		}

		return Info.SUCCESSFUL;
	}

	public Info gaussian_elimination_2x2(int k) {

		int e11 = colptr[k];
		int e21 = index(k + 1, k);
		int e22 = colptr[k + 1];

		double de11 = data[e11];
		double de21 = data[e21];
		double de22 = data[e22];
		double delta = de11 * de22 - de21 * de21;

		if (delta == 0.0) {
			return Info.NUMERICAL_ISSUE;
		}

		// inverse_inplace_2x2
		data[e11] = de22 / delta;
		data[e21] = -de21 / delta;
		data[e22] = de11 / delta;

		de11 = data[e11];
		de21 = data[e21];
		de22 = data[e22];

		int l1ptr = index(k + 2, k);
		int l2ptr = index(k + 2, k + 1);
		int ldim = n - k - 2;

		if (ldim != 0) {
			double[] X0 = new double[ldim];
			double[] X1 = new double[ldim];

			for (int i = 0; i < ldim; i++) {
				X0[i] = data[l1ptr + i] * de11 + data[l2ptr + i] * de21;
				X1[i] = data[l1ptr + i] * de21 + data[l2ptr + i] * de22;
			}

			for (int j = 0; j < ldim; j++) {
				int col = colptr[j + k + 2];
				for (int i = 0; i < ldim - j; i++) {
					double v = X0[j + i] * data[l1ptr + j] + X1[j + i] * data[l2ptr + j];
					data[col + i] -= v;
				}
			}

			for (int i = 0; i < ldim; i++) {
				data[l1ptr + i] = X0[i];
				data[l2ptr + i] = X1[i];
			}
		}

		return Info.SUCCESSFUL;
	}

	public void compute(Matrix mat) {
		if(DEBUG) debug('-', "compute BKLDLT");
		if(DEBUG) debug("mat: ", mat);
		
		n = mat.rows;
		perm = new int[n];

		// setLinSpaced
		for (int i = 0; i < n; i++) {
			perm[i] = i;
		}

		data = new double[(n * (n + 1)) / 2];
		compute_pointer();
		copy_data(mat);

		int k;
		for (k = 0; k < n - 1; k++) {
			boolean is_1x1 = permutate_mat(k);
			if (is_1x1) {
				info = gaussian_elimination_1x1(k);
			} else {
				info = gaussian_elimination_2x2(k);
				k++;
			}

			if (info != Info.SUCCESSFUL) {
				break;
			}
		}

		if (k == n - 1) {
			double akk = diag_coeff(k);
			if (akk == 0.0) {
				info = Info.NUMERICAL_ISSUE;
				data[colptr[k]] = Double.NaN;
			} else {
				data[colptr[k]] = 1.0 / data[colptr[k]];
			}
		}

		compress_permutation();

		computed = true;
	
		if(DEBUG) debug('-', "compute BKLDLT - end");
	}

	public void solve_inplace(double[] b) {
		if(DEBUG) debug('-', "solve BKLDLT");
		if(DEBUG) debug("b: ", b);

		for (Pair p : permc) {
			swap(b, p.a, p.b);
		}

		int end = perm[n - 1] < 0 ? n - 3 : n - 2;
		int i;
		for (i = 0; i <= end; i++) {
			int b1size = n - i - 1;
			int b2size = b1size - 1;
			if (perm[i] >= 0) {
				int l = index(i + 1, i);
				for (int k = 0; k < b1size; k++) {
					b[i + 1 + k] -= data[l + k] * b[i];
				}
			} else {
				int l1 = index(i + 2, i);
				int l2 = index(i + 2, i + 1);
				for (int k = 0; k < b2size; k++) {
					b[i + 2 + k] -= data[l1 + k] * b[i] + data[l2 + k] * b[i + 1];
				}
				i++;
			}
		}

		for (i = 0; i < n; i++) {
			double e11 = diag_coeff(i);

			if (perm[i] >= 0) {
				b[i] *= e11;
			} else {
				double e21 = coeff(i + 1, i);
				double e22 = diag_coeff(i + 1);
				double wi = b[i] * e11 + b[i + 1] * e21;
				b[i + 1] = b[i] * e21 + b[i + 1] * e22;
				b[i] = wi;
				i++;
			}
		}

		i = perm[n - 1] < 0 ? n - 3 : n - 2;
		for (; i >= 0; i--) {
			int ldim = n - i - 1;
			int l = index(i + 1, i);
			double dot1 = 0.0;
			for (int k = 0; k < ldim; k++) {
				dot1 += b[i + 1 + k] * data[l + k];
			}
			b[i] -= dot1;

			if (perm[i] < 0) {
				int l2 = index(i + 1, i - 1);
				double dot2 = 0.0;
				for (int k = 0; k < ldim; k++) {
					dot2 += b[i + 1 + k] * data[l2 + k];
				}
				b[i - 1] -= dot2;
				i--;
			}
		}

		for (i = permc.size() - 1; i >= 0; i--) {
			Pair p = permc.get(i);
			swap(b, p.a, p.b);
		}
		
		if(DEBUG) debug('-', "solve BKLDLT - end");
	}

	public double[] solve(double[] b) {
		double[] copy = b.clone();
		solve_inplace(copy);
		return copy;
	}

	public BKLDLT(Matrix mat) {
		compute(mat);
	}

	public BKLDLT() {}
}
