package org.generateme.lbfgsb;

import static org.generateme.lbfgsb.Debug.*;

import java.util.ArrayList;
import java.util.Collections;

public final class Cauchy {

	public double[] xcp;
	public double[] vecc;
	public ArrayList<Integer> newact_set;
	public ArrayList<Integer> fv_set;

	public static final double eps = Math.ulp(1.0);

	public static final int search_greater(double[] brk, ArrayList<Integer> ord, double t, int start) {
		int nord = ord.size();
		int i;
		for (i = start; i < nord; i++) {
			if (brk[ord.get(i)] > t)
				break;
		}
		return i;
	}

	public Cauchy(BFGSMat bfgs, double[] x0, double[] g, double[] lb, double[] ub) {
		if (DEBUG) {
			debug('=', "Cauchy");
			debug("x0: ", x0);
			debug(" g: ", g);
		}

		int n = x0.length;
		xcp = x0.clone();
		vecc = new double[2 * bfgs.m_ncorr];
		newact_set = new ArrayList<Integer>(n);
		fv_set = new ArrayList<Integer>(n);

		double[] brk = new double[n];
		double[] vecd = new double[n];

		ArrayList<Integer> ord = new ArrayList<Integer>(n);

		for (int i = 0; i < n; i++) {
			if (g[i] < 0.0) {
				brk[i] = (x0[i] - ub[i]) / g[i];
			} else if (g[i] > 0.0) {
				brk[i] = (x0[i] - lb[i]) / g[i];
			} else {
				brk[i] = Double.POSITIVE_INFINITY;
			}

			boolean iszero = (brk[i] == 0.0);
			vecd[i] = iszero ? 0.0 : -g[i];

			if (Double.isInfinite(brk[i])) {
				fv_set.add(i);
			} else if (!iszero) {
				ord.add(i);
			}
		}

		Collections.sort(ord, (a, b) -> Double.compare(brk[a], brk[b]));

		int nord = ord.size();
		int nfree = fv_set.size();

		if (DEBUG) {
			debug("    brk: ", brk);
			debug("   vecd: ", vecd);
			debug("   nord: " + nord);
			debug("    ord: " + ord);
			debug("  nfree: " + nfree);
			debug(" fv_set: " + fv_set);
		}

		if ((nfree < 1) && (nord < 1)) {
			if (DEBUG)
				debug('=', "leaving Cauchy, nfree < 1 && nord < 1");
			return;
		}

		double[] vecp = new double[2 * bfgs.m_ncorr];
		bfgs.apply_Wtv(vecd, vecp);
		double fp = -Vector.squaredNorm(vecd);

		double[] cache = new double[2 * bfgs.m_ncorr];
		bfgs.apply_Mv(vecp, cache);
		double fpp = -bfgs.m_theta * fp - Vector.dot(vecp, cache);

		double deltatmin = -fp / fpp;

		double il = 0.0;
		int b = 0;
		double iu = nord < 1 ? Double.POSITIVE_INFINITY : brk[ord.get(b)];
		double deltat = iu - il;

		boolean crossed_all = false;
		double[] wact = new double[2 * bfgs.m_ncorr];

		while (deltatmin >= deltat) {
			for (int i = 0; i < vecc.length; i++) {
				vecc[i] += deltat * vecp[i];
			}

			int act_begin = b;
			int act_end = search_greater(brk, ord, iu, b) - 1;

			if ((nfree == 0) && (act_end == nord - 1)) {
				for (int i = act_begin; i <= act_end; i++) {
					int act = ord.get(i);
					xcp[act] = vecd[act] > 0.0 ? ub[act] : lb[act];
					newact_set.add(act);
				}
				crossed_all = true;
				break;
			}

			fp += deltat * fpp;

			for (int i = act_begin; i <= act_end; i++) {
				int act = ord.get(i);
				xcp[act] = vecd[act] > 0.0 ? ub[act] : lb[act];
				double zact = xcp[act] - x0[act];
				double gact = g[act];
				double ggact = gact * gact;
				wact = bfgs.Wb(act);
				bfgs.apply_Mv(wact, cache);
				fp += ggact + bfgs.m_theta * gact * zact - gact * Vector.dot(cache, vecc);
				fpp -= bfgs.m_theta * ggact + 2 * gact * Vector.dot(cache, vecp) + ggact * Vector.dot(cache, wact);
				for (int j = 0; j < vecp.length; j++) {
					vecp[j] += gact * wact[j];
				}
				vecd[act] = 0.0;
				newact_set.add(act);
			}

			deltatmin = -fp / fpp;
			il = iu;
			b = act_end + 1;
			if (b >= nord) {
				break;
			}
			iu = brk[ord.get(b)];
			deltat = iu - il;
		}

		if (fpp < eps) {
			deltatmin = -fp / eps;
		}

		if (!crossed_all) {
			deltatmin = Math.max(deltatmin, 0.0);
			for (int j = 0; j < vecc.length; j++) {
				vecc[j] += deltatmin * vecp[j];
			}
			double tfinal = il + deltatmin;
			for (int i = 0; i < nfree; i++) {
				int coord = fv_set.get(i);
				xcp[coord] = x0[coord] + tfinal * vecd[coord];
			}
			for (int i = b; i < nord; i++) {
				int coord = ord.get(i);
				xcp[coord] = x0[coord] + tfinal * vecd[coord];
				fv_set.add(coord);
			}
		}

		if (DEBUG) {
			debug("vecc: ", vecc);
			debug("xcp: ", xcp);
			debug("newact_set: " + newact_set);
			debug("fv_set: " + fv_set);
			debug('=', "Cauchy - end");
		}
	}
}
