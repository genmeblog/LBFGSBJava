package org.generateme.lbfgsb;

import static org.generateme.lbfgsb.Debug.*;

import java.util.ArrayList;

public final class BFGSMat {
	public int m_m;
	public double m_theta;
	public Matrix m_s;
	public Matrix m_y;
	public double[] m_ys;
	public double[] m_alpha;
	public int m_ncorr;
	public int m_ptr;

	public Matrix permMinv;
	public BKLDLT permMsolver;

	public void reset(int n, int m) {
		m_m = m;
		m_theta = 1.0;
		m_s = Matrix.resize(m_s, n, m);
		m_y = Matrix.resize(m_y, n, m);
		m_ys = Vector.resize(m_ys, m);
		m_alpha = Vector.resize(m_alpha, m);
		m_ncorr = 0;
		m_ptr = m;

		permMinv = Matrix.resize(permMinv, 2 * m, 2 * m);
		permMinv.setAll(0.0);
		permMinv.setDiag(1.0); // sets diagonal to 1.0
		permMsolver = new BKLDLT();
	}

	public BFGSMat() {
	}

	public BFGSMat(int n, int m) {
		reset(n, m);
	}

	public void add_correction(double[] s, double[] y) {
		if (DEBUG) {
			debug('-', "add correction");
			debug("s: ", s);
			debug("y: ", y);
		}

		int loc = m_ptr % m_m;

		m_s.setCol(loc, s);
		m_y.setCol(loc, y);

		double ys = Vector.dot(s, y);
		m_ys[loc] = ys;

		m_theta = m_y.colSquaredNorm(loc) / ys;

		if (m_ncorr < m_m)
			m_ncorr++;

		m_ptr = loc + 1;

		permMinv.set(loc, loc, -ys);

		for (int i = 0; i < m_ncorr; i++) {
			double Ss = this.m_s.colDot(i, s);
			permMinv.set(m_m + loc, m_m + i, Ss);
			permMinv.set(m_m + i, m_m + loc, Ss);
		}

		int len = m_ncorr - 1;
		if (m_ncorr >= m_m) {
			for (int i = 0; i < m_m; i++) {
				permMinv.set(m_m + i, loc, 0.0);
			}
		}

		int yloc = (loc + m_m - 1) % m_m;
		int mloc = m_m + loc;
		for (int i = 0; i < len; i++) {
			permMinv.set(mloc, yloc, m_y.colDot(yloc, s));
			yloc = (yloc + m_m - 1) % m_m;
		}

		for (int i = 0; i < m_m; i++) {
			for (int j = 0; j < m_m; j++) {
				permMinv.set(m_m + i, m_m + j, permMinv.get(m_m + i, m_m + j) * m_theta);
			}
		}

		permMsolver.compute(permMinv);

		for (int i = 0; i < m_m; i++) {
			for (int j = 0; j < m_m; j++) {
				permMinv.set(m_m + i, m_m + j, permMinv.get(m_m + i, m_m + j) / m_theta);
			}
		}

		if (DEBUG)
			debug('-', "add correction - end");
	}

	public void apply_Wtv(double[] v, double[] res) {
		if (DEBUG) {
			debug('-', "apply_Wtv");
			debug("v:  ", v);
		}

		for (int i = 0; i < m_ncorr; i++)
			res[i] = m_y.colDot(i, v);

		for (int i = 0; i < m_ncorr; i++)
			res[i + m_ncorr] = m_theta * m_s.colDot(i, v);

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "apply_Wtv - end");
		}
	}

	public void apply_Mv(double[] v, double[] res) {
		if (DEBUG) {
			debug('-', "apply Mv");
			debug("v:  ", v);
		}

		if (m_ncorr < 1) {
			if (DEBUG)
				debug('-', "leaving apply_Mv, m_ncorr < 1");
			return;
		}

		double[] vpadding = new double[2 * m_m];
		for (int i = 0; i < m_ncorr; i++) {
			vpadding[i] = v[i];
			vpadding[m_m + i] = v[m_ncorr + i];
		}

		permMsolver.solve_inplace(vpadding);

		for (int i = 0; i < m_ncorr; i++) {
			res[i] = vpadding[i];
			res[i + m_ncorr] = vpadding[m_m + i];
		}

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "apply Mv - end");
		}
	}

	public double[] Wb(int b) {
		if (DEBUG) {
			debug('-', "Wb");
			debug("b: " + b);
		}

		double[] res = new double[2 * m_ncorr];

		for (int j = 0; j < m_ncorr; j++) {
			res[j] = m_y.get(b, j);
			res[m_ncorr + j] = m_theta * m_s.get(b, j);
		}

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "Wb - end");
		}

		return res;
	}

	public Matrix Wb(ArrayList<Integer> b) {
		if (DEBUG) {
			debug('-', "Wb");
			debug("b: " + b);
		}

		int nb = b.size();
		Matrix res = new Matrix(nb, 2 * m_ncorr);

		for (int j = 0; j < m_ncorr; j++) { // col
			for (int i = 0; i < nb; i++) { // row
				int row = b.get(i);
				res.set(i, j, m_y.get(row, j));
				res.set(i, j + m_ncorr, m_s.get(row, j));
			}
		}

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "Wb - end");
		}

		return res;
	}

	public boolean apply_WtPv(ArrayList<Integer> P_set, double[] v, double[] res) {
		return apply_WtPv(P_set, v, res, false);
	}

	public boolean apply_WtPv(ArrayList<Integer> P_set, double[] v, double[] res, boolean test_zero) {
		if (DEBUG) {
			debug('-', "apply_WtPv, test_zero=" + test_zero);
			debug("P_set: " + P_set);
			debug("v: ", v);
		}

		ArrayList<Integer> Pptr = P_set;
		double[] vptr = v;
		int nP = P_set.size();

		if (test_zero) {
			ArrayList<Integer> P_reduced = new ArrayList<Integer>(nP);
			ArrayList<Double> v_reduced = new ArrayList<Double>(nP);

			for (int i = 0; i < nP; i++) {
				if (vptr[i] != 0.0) {
					P_reduced.add(Pptr.get(i));
					v_reduced.add(v[i]);
				}
			}

			nP = P_reduced.size();
			Pptr = P_reduced;
			vptr = new double[nP];
			for (int i = 0; i < nP; i++) {
				vptr[i] = v_reduced.get(i);
			}
		}

		if (m_ncorr < 1 || nP < 1) {
			if (DEBUG)
				debug('-', "leaving apply_WtPv");
			Vector.setAll(res, 0.0);
			return false;
		}

		for (int j = 0; j < m_ncorr; j++) {
			double resy = 0.0;
			double ress = 0.0;
			for (int i = 0; i < nP; i++) {
				int row = Pptr.get(i);
				resy += m_y.get(row, j) * vptr[i];
				ress += m_s.get(row, j) * vptr[i];
			}
			res[j] = resy;
			res[m_ncorr + j] = m_theta * ress;
		}

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "apply_WtPv - end");
		}

		return true;
	}

	public boolean apply_PtWMv(ArrayList<Integer> P_set, double[] v, double[] res, double scale) {
		if (DEBUG) {
			debug('-', "apply_PtWMv, scale=" + scale);
			debug("P_set: " + P_set);
			debug("v: ", v);
		}

		int nP = P_set.size();
		Vector.setAll(res, 0.0);
		if (m_ncorr < 1 || nP < 1) {
			if (DEBUG)
				debug('-', "leaving apply_PTWMv, m_ncorr < 1 || np < 1");
			return false;
		}

		double[] Mv = new double[2 * m_ncorr];
		apply_Mv(v, Mv);
		for (int i = 0; i < m_ncorr; i++) {
			Mv[i + m_ncorr] *= m_theta;
		}

		for (int j = 0; j < m_ncorr; j++) {
			double Mvy = Mv[j];
			double Mvs = Mv[m_ncorr + j];
			for (int i = 0; i < nP; i++) {
				int row = P_set.get(i);
				res[i] += Mvy * m_y.get(row, j) + Mvs * m_s.get(row, j);
			}
		}

		for (int i = 0; i < res.length; i++) {
			res[i] *= scale;
		}

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "apply_PtWMv - end");
		}

		return true;
	}

	public boolean apply_PtWMv(Matrix WP, double[] v, double[] res, double scale) {
		if (DEBUG) {
			debug('-', "apply_PtWMv, scale=" + scale);
			debug("WP: ", WP);
			debug("v:", v);
		}

		int nP = WP.rows;

		if (m_ncorr < 1 || nP < 1) {
			if (DEBUG)
				debug('-', "leaving apply_PtWMv, m_ncorr < 1 || nP < 1");
			Vector.setAll(res, 0.0);
			return false;
		}

		double[] Mv = new double[2 * m_ncorr];
		apply_Mv(v, Mv);

		for (int i = 0; i < m_ncorr; i++) {
			Mv[i + m_ncorr] *= m_theta;
		}

		for (int i = 0; i < res.length; i++) {
			double dot = 0.0;
			for (int k = 0; k < Mv.length; k++) {
				dot += WP.get(i, k) * Mv[k];
			}
			res[i] = scale * dot;
		}

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "apply_PtWMv - end");
		}

		return true;
	}

	public void compute_FtBAb(Matrix WF, ArrayList<Integer> fv_set, ArrayList<Integer> newact_set, double[] Wd,
			double[] drt, double[] res) {

		if (DEBUG) {
			debug('-', "compute_FtBAb");
			debug("WF: ", WF);
			debug("fv_set: " + fv_set);
			debug("newact_set: " + newact_set);
			debug("Wd: ", Wd);
			debug("drt: ", drt);
		}

		int nact = newact_set.size();
		int nfree = WF.rows;

		if (m_ncorr < 1 || nact < 1 || nfree < 1) {
			if (DEBUG)
				debug('-', "leaving compute_FtBAb, m_ncorr < 1 || nact < 1 || nfree < 1");
			Vector.setAll(res, 0.0);
			return;
		}

		double[] rhs = new double[2 * m_ncorr];
		if (nact <= nfree) {
			double[] Ad = new double[nfree];
			for (int i = 0; i < nact; i++)
				Ad[i] = drt[newact_set.get(i)];
			apply_WtPv(newact_set, Ad, rhs);
		} else {
			double[] Fd = new double[nfree];
			for (int i = 0; i < nfree; i++)
				Fd[i] = drt[fv_set.get(i)];

			for (int i = 0; i < rhs.length; i++) {
				double dot = 0.0;
				for (int k = 0; k < nfree; k++) {
					dot += WF.get(k, i) * Fd[k];
				}
				rhs[i] = dot;
			}

			for (int i = 0; i < m_ncorr; i++) {
				rhs[i + m_ncorr] *= m_theta;
			}

			Vector.sub(Wd, rhs, rhs);
		}

		apply_PtWMv(WF, rhs, res, -1.0);

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "compute_FtBAb - end");
		}
	}

	public boolean apply_PtBQv(final Matrix WP, final ArrayList<Integer> Q_set, final double[] v, final double[] res) {
		return apply_PtBQv(WP, Q_set, v, res, false);
	}

	public boolean apply_PtBQv(final Matrix WP, final ArrayList<Integer> Q_set, final double[] v, final double[] res,
			boolean test_zero) {

		if (DEBUG) {
			debug('-', "PtBQv, test_zero=" + test_zero);
			debug("WP: ", WP);
			debug("Q_set: " + Q_set);
			debug("v: ", v);
		}

		int nP = WP.rows;
		int nQ = Q_set.size();

		if (m_ncorr < 1 || nP < 1 || nQ < 1) {
			if (DEBUG)
				debug('-', "leaving PtBQv, m_ncorr < 1 || nP < 1 || nQ < 1");

			Vector.setAll(res, 0.0);
			return false;
		}

		double[] WQtv = new double[2 * m_ncorr];
		boolean nonzero = apply_WtPv(Q_set, v, WQtv, test_zero);
		if (!nonzero) {
			if (DEBUG)
				debug('-', "leaving PtBQv, !nonzero");
			Vector.setAll(res, 0.0);
			return false;
		}

		double[] MWQtv = new double[2 * m_ncorr];
		apply_Mv(WQtv, MWQtv);

		for (int i = 0; i < m_ncorr; i++) {
			MWQtv[i + m_ncorr] *= m_theta;
		}

		for (int row = 0; row < WP.rows; row++) {
			double dot = 0.0;
			for (int col = 0; col < WP.cols; col++) {
				dot += WP.get(row, col) * MWQtv[col];
			}
			res[row] = -dot;
		}

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "PtBQv - end");
		}

		return true;
	}

	public void solve_PtBP(Matrix WP, double[] v, double[] res) {
		if (DEBUG) {
			debug('-', "solve_PtBP");
			debug("WP: ", WP);
			debug("v: ", v);
		}

		int nP = WP.rows;
		if (m_ncorr < 1 || nP < 1) {
			for (int i = 0; i < res.length; i++) {
				res[i] = v[i] / m_theta;
			}

			if (DEBUG) {
				debug("res: ", res);
				debug('-', "leaving PtBQv, m_ncorr < 1 || nP < 1");
			}

			return;
		}

		Matrix mid = new Matrix(2 * m_ncorr, 2 * m_ncorr);
		for (int j = 0; j < m_ncorr; j++) {
			for (int i = 0; i < m_ncorr - j; i++) {
				double dot = 0.0;
				for (int k = 0; k < nP; k++) {
					dot += WP.get(k, i + j) * WP.get(k, j);
				}

				mid.set(i + j, j, permMinv.get(i + j, j) - dot / m_theta);
			}
		}

		for (int i = 0; i < m_ncorr; i++) {
			for (int j = 0; j < m_ncorr; j++) {
				double dot = 0.0;
				for (int k = 0; k < nP; k++) {
					dot += WP.get(k, i + m_ncorr) * WP.get(k, j);
				}
				mid.set(i + m_ncorr, j, permMinv.get(i + m_m, j) - dot);
			}
		}

		for (int j = 0; j < m_ncorr; j++) {
			for (int i = 0; i < m_ncorr - j; i++) {
				double dot = 0.0;
				int right = WP.cols - (m_ncorr - j) + i;

				for (int k = 0; k < WP.rows; k++) {
					dot += WP.get(k, right) * WP.get(k, m_ncorr + j);
				}
				mid.set(i + m_ncorr + j, m_ncorr + j, m_theta * (permMinv.get(i + m_m + j, m_m + j) - dot));
			}
		}

		BKLDLT midsolver = new BKLDLT(mid);

		double[] WPv = new double[WP.cols];
		for (int i = 0; i < WP.cols; i++) {
			for (int j = 0; j < WP.rows; j++) {
				WPv[i] += WP.get(j, i) * v[j];
			}
		}

		for (int i = WPv.length - 1; i >= m_ncorr; i--) {
			WPv[i] *= m_theta;
		}

		midsolver.solve_inplace(WPv);

		for (int i = WPv.length - 1; i >= m_ncorr; i--) {
			WPv[i] *= m_theta;
		}

		double t2 = m_theta * m_theta;

		for (int i = 0; i < res.length; i++) {
			double dot = 0.0;
			for (int k = 0; k < WPv.length; k++) {
				dot += WP.get(i, k) * WPv[k];
			}
			res[i] = v[i] / m_theta + dot / t2;
		}

		if (DEBUG) {
			debug("res: ", res);
			debug('-', "PtBQv - end");
		}
	}
}
