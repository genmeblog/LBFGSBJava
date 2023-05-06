package org.generateme.lbfgsb;

import static org.generateme.lbfgsb.Debug.*;

import java.util.ArrayList;

public final class SubspaceMin {

	public final static double meps = -Math.ulp(1.0);

	public final static void subvec_assign(double[] v, ArrayList<Integer> ind, double[] rhs) {
		int nsub = ind.size();
		for (int i = 0; i < nsub; i++) {
			v[ind.get(i)] = rhs[i];
		}
	}

	public final static double[] subvec(final double[] v, final ArrayList<Integer> ind) {
		int nsub = ind.size();
		double[] res = new double[nsub];
		for (int i = 0; i < nsub; i++)
			res[i] = v[ind.get(i)];
		return res;
	}

	public final static boolean in_bounds(double[] x, double[] lb, double[] ub) {
		int n = x.length;
		for (int i = 0; i < n; i++) {
			if (x[i] < lb[i] || x[i] > ub[i])
				return false;
		}
		return true;
	}

	public final static boolean P_converged(ArrayList<Integer> yP_set, double[] vecy, double[] vecl, double[] vecu) {
		int nP = yP_set.size();
		for (int i = 0; i < nP; i++) {
			int coord = yP_set.get(i);
			if (vecy[coord] < vecl[coord] || vecy[coord] > vecu[coord])
				return false;
		}
		return true;
	}

	public final static boolean L_converged(ArrayList<Integer> yL_set, double[] lambda) {
		int nL = yL_set.size();
		for (int i = 0; i < nL; i++) {
			int coord = yL_set.get(i);
			if (lambda[coord] < 0.0)
				return false;
		}
		return true;
	}

	public final static boolean U_converged(ArrayList<Integer> yU_set, double[] mu) {
		int nU = yU_set.size();
		for (int i = 0; i < nU; i++) {
			int coord = yU_set.get(i);
			if (mu[coord] < 0.0)
				return false;
		}
		return true;
	}

	public final static void subspace_minimize(BFGSMat bfgs, double[] x0, double[] g, double[] lb, double[] ub,
			Cauchy cauchy, int maxit, double[] drt) {
		if (DEBUG)
			debug('-', "subspace minimize");

		double[] Wd = cauchy.vecc;

		Vector.sub(cauchy.xcp, x0, drt);

		int nfree = cauchy.fv_set.size();

		if (nfree < 1) {
			if (DEBUG)
				debug('-', "leaving subspace_minimize, nfree<1");
			return;
		}

		Matrix WF = bfgs.Wb(cauchy.fv_set);

		double[] vecc = new double[nfree];
		bfgs.compute_FtBAb(WF, cauchy.fv_set, cauchy.newact_set, Wd, drt, vecc);
		double[] vecl = new double[nfree];
		double[] vecu = new double[nfree];

		for (int i = 0; i < nfree; i++) {
			int coord = cauchy.fv_set.get(i);
			vecl[i] = lb[coord] - x0[coord];
			vecu[i] = ub[coord] - x0[coord];
			vecc[i] += g[coord];
		}

		double[] vecy = new double[nfree];
		double[] veccm = vecc.clone();
		for (int i = 0; i < veccm.length; i++) {
			veccm[i] *= -1.0;
		}

		bfgs.solve_PtBP(WF, veccm, vecy);

		if (in_bounds(vecy, vecl, vecu)) {
			subvec_assign(drt, cauchy.fv_set, vecy);
			if (DEBUG)
				debug('-', "leaving subspace_minimize, solution has been found");
			return;
		}

		double[] yfallback = vecy.clone();

		double[] lambda = new double[nfree];
		double[] mu = new double[nfree];

		ArrayList<Integer> L_set = new ArrayList<Integer>(nfree / 3);
		ArrayList<Integer> U_set = new ArrayList<Integer>(nfree / 3);
		ArrayList<Integer> yL_set = new ArrayList<Integer>(nfree / 3);
		ArrayList<Integer> yU_set = new ArrayList<Integer>(nfree / 3);
		ArrayList<Integer> P_set = new ArrayList<Integer>(nfree);
		ArrayList<Integer> yP_set = new ArrayList<Integer>(nfree);

		int k;
		for (k = 0; k < maxit; k++) {
			L_set.clear();
			U_set.clear();
			P_set.clear();
			yL_set.clear();
			yU_set.clear();
			yP_set.clear();

			for (int i = 0; i < nfree; i++) {
				int coord = cauchy.fv_set.get(i);
				double li = vecl[i];
				double ui = vecu[i];

				if ((vecy[i] < li) || (vecy[i] == li && lambda[i] >= 0.0)) {
					L_set.add(coord);
					yL_set.add(i);
					vecy[i] = li;
					mu[i] = 0.0;
				} else if ((vecy[i] > ui) || (vecy[i] == ui && mu[i] >= 0.0)) {
					U_set.add(coord);
					yU_set.add(i);
					vecy[i] = ui;
					lambda[i] = 0.0;
				} else {
					P_set.add(coord);
					yP_set.add(i);
					lambda[i] = 0.0;
					mu[i] = 0.0;
				}
			}

			Matrix WP = bfgs.Wb(P_set);
			int nP = P_set.size();

			if (nP > 0) {
				double[] rhs = subvec(vecc, yP_set);
				double[] lL = subvec(vecl, yL_set);
				double[] uU = subvec(vecu, yU_set);
				double[] tmp = new double[nP];
				boolean nonzero = bfgs.apply_PtBQv(WP, L_set, lL, tmp, true);
				if (nonzero)
					for (int i = 0; i < rhs.length; i++) {
						rhs[i] += tmp[i];
					}
				nonzero = bfgs.apply_PtBQv(WP, U_set, uU, tmp, true);
				if (nonzero)
					for (int i = 0; i < rhs.length; i++) {
						rhs[i] += tmp[i];
					}

				for (int i = 0; i < rhs.length; i++) {
					rhs[i] *= -1.0;
				}

				bfgs.solve_PtBP(WP, rhs, tmp);
				subvec_assign(vecy, yP_set, tmp);
			}

			int nL = L_set.size();
			int nU = U_set.size();
			double[] Fy = new double[2 * bfgs.m_ncorr];
			if (nL > 0 || nU > 0)
				bfgs.apply_WtPv(cauchy.fv_set, vecy, Fy);

			if (nL > 0) {
				double[] res = new double[L_set.size()];
				bfgs.apply_PtWMv(L_set, Fy, res, -1.0);
				double[] svc = subvec(vecc, yL_set);
				double[] svy = subvec(vecy, yL_set);
				for (int i = 0; i < res.length; i++) {
					res[i] += svc[i] + bfgs.m_theta * svy[i];
				}
				subvec_assign(lambda, yL_set, res);
			}

			if (nU > 0) {
				double[] negRes = new double[U_set.size()];
				bfgs.apply_PtWMv(U_set, Fy, negRes, -1.0);
				double[] svc = subvec(vecc, yU_set);
				double[] svy = subvec(vecy, yU_set);
				for (int i = 0; i < negRes.length; i++) {
					negRes[i] += (svc[i] + bfgs.m_theta * svy[i]);
					negRes[i] *= -1.0;
				}
				subvec_assign(mu, yU_set, negRes);
			}

			if (L_converged(yL_set, lambda) && U_converged(yU_set, mu) && P_converged(yP_set, vecy, vecl, vecu))
				break;
		}

		if (k >= maxit) {
			for (int i = 0; i < vecy.length; i++) {
				vecy[i] = Math.max(Math.min(vecy[i], vecu[i]), vecl[i]);
			}
			subvec_assign(drt, cauchy.fv_set, vecy);

			double dg = Vector.dot(drt, g);

			if (dg <= meps) {
				if (DEBUG)
					debug('-', "leaving subspace_minimize, projected");
				return;
			}

			for (int i = 0; i < vecy.length; i++) {
				vecy[i] = Math.max(Math.min(yfallback[i], vecu[i]), vecl[i]);
			}

			subvec_assign(drt, cauchy.fv_set, vecy);
			dg = Vector.dot(drt, g);

			if (dg <= meps) {
				if (DEBUG)
					debug('-', "leaving subspace_minimize, projected unconstrained");
				return;
			}

			// If still not, fall back to the unconstrained solution

			subvec_assign(drt, cauchy.fv_set, yfallback);

			if (DEBUG)
				debug('-', "leaving subspace_minimize, projected unconstrained");

			return;
		}

		subvec_assign(drt, cauchy.fv_set, vecy);

		if (DEBUG)
			debug('-', "leaving subspace_minimize, converged");

	}
}
