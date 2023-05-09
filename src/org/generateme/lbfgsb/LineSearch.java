package org.generateme.lbfgsb;

import static org.generateme.lbfgsb.Debug.*;

public final class LineSearch {

	static public class Bool {
		boolean b;

		public Bool() {
			this(false);
		}

		public Bool(boolean b) {
			this.b = b;
		}
	}

	public double fx;
	public double step;
	public double dg;

	public static final double eps = Math.ulp(1.0);

	public static final double quadratic_minimizer(double a, double b, double fa, double ga, double fb) {
		double ba = b - a;
		return a + 0.5 * ba * ba * ga / (fa - fb + ba * ga);
	}

	public static final double quadratic_minimizer(double a, double b, double ga, double gb) {
		return a + (b - a) * ga / (ga - gb);
	}

	public static final double cubic_minimizer(double a, double b, double fa, double fb, double ga, double gb,
			Bool exists) {
		double apb = a + b;
		double ba = b - a;
		double ba2 = ba * ba;
		double fba = fb - fa;
		double gba = gb - ga;

		double z3 = (ga + gb) * ba - 2.0 * fba;
		double z2 = 0.5 * (gba * ba2 - 3.0 * apb * z3);
		double z1 = fba * ba2 - apb * z2 - (a * apb + b * b) * z3;

		if (Math.abs(z3) < eps * Math.abs(z2) || Math.abs(z3) < eps * Math.abs(z1)) {
			// Minimizer exists if c2 > 0
			exists.b = (z2 * ba > 0.0);
			// Return the end point if the minimizer does not exist
			return exists.b ? (-0.5 * z1 / z2) : b;
		}

		// Now we can assume z3 > 0
		// The minimizer is a solution to the equation c1 + 2*c2 * x + 3*c3 * x^2 = 0
		// roots = -(z2/z3) / 3 (+-) sqrt((z2/z3)^2 - 3 * (z1/z3)) / 3
		//
		// Let u = z2/(3z3) and v = z1/z2
		// The minimizer exists if v/u <= 1
		double u = z2 / (3.0 * z3), v = z1 / z2;
		double vu = v / u;
		exists.b = (vu <= 1.0);
		if (!exists.b)
			return b;

		// We need to find a numerically stable way to compute the roots, as z3 may
		// still be small
		//
		// If |u| >= |v|, let w = 1 + sqrt(1-v/u), and then
		// r1 = -u * w, r2 = -v / w, r1 does not need to be the smaller one
		//
		// If |u| < |v|, we must have uv <= 0, and then
		// r = -u (+-) sqrt(delta), where
		// sqrt(delta) = sqrt(|u|) * sqrt(|v|) * sqrt(1-u/v)
		double r1 = 0.0, r2 = 0.0;
		if (Math.abs(u) >= Math.abs(v)) {
			double w = 1.0 + Math.sqrt(1.0 - vu);
			r1 = -u * w;
			r2 = -v / w;
		} else {
			double sqrtd = Math.sqrt(Math.abs(u)) * Math.sqrt(Math.abs(v)) * Math.sqrt(1 - u / v);
			r1 = -u - sqrtd;
			r2 = -u + sqrtd;
		}
		return (z3 * ba > 0.0) ? (Math.max(r1, r2)) : (Math.min(r1, r2));

	}

	public static final double deltal = 1.1;
	public static final double deltau = 0.66;

	public static final double step_selection(double al, double au, double at, double fl, double fu, double ft, double gl,
			double gu, double gt) {

		if (al == au)
			return al;

		if (Double.isInfinite(ft) || Double.isInfinite(gt))
			return (al + at) / 2.0;

		Bool ac_exists = new Bool();
		double ac = cubic_minimizer(al, at, fl, ft, gl, gt, ac_exists);
		double aq = quadratic_minimizer(al, at, fl, gl, ft);

		if (ft > fl) {
			if (!ac_exists.b)
				return aq;

			return (Math.abs(ac - al) < Math.abs(aq - al)) ? ac : ((aq + ac) / 2.0);

		}

		double as = quadratic_minimizer(al, at, gl, gt);
		if (gt * gl < 0.0)
			return (Math.abs(ac - at) >= Math.abs(as - at)) ? ac : as;

		if (Math.abs(gt) < Math.abs(gl)) {
			double res = ac_exists.b && (((ac - at) * (at - al)) > 0.0) && (Math.abs(ac - at) < Math.abs(as - at)) ? ac : as;
			return (at > al) ? Math.min(at + deltau * (au - at), res) : Math.max(at + deltau * (au - at), res);
		}

		if (Double.isInfinite(au) || Double.isInfinite(fu) || Double.isInfinite(gu))
			return at + deltal * (at - al);

		Bool ae_exists = new Bool();
		double ae = cubic_minimizer(at, au, ft, fu, gt, gu, ae_exists);
		return (at > al) ? Math.min(at + deltau * (au - at), ae) : Math.max(at + deltau * (au - at), ae);
	}

	public static final double delta = 1.1;

	public LineSearch(IGradFunction f, Parameters param, double[] xp, double[] drt, double step_max, double _step,
			double _fx, double[] grad, double _dg, double[] x) throws LBFGSBException {
		if (DEBUG) {
			debug('-', "line search");
			debug("      xp: ", xp);
			debug("       x: ", x);
			debug("      fx: " + _fx);
			debug("    grad: ", grad);
			debug("      dg: " + _dg);
			debug("    step: " + _step);
			debug("step_max: " + step_max);
			debug("     drt: ", drt);
		}

		fx = _fx;
		step = _step;
		dg = _dg;

		if (step <= 0.0)
			throw new LBFGSBException("step must be positive, step=" + step);
		if (step > step_max)
			throw new LBFGSBException("step exceeds max step, step=" + step + ", step_max=" + step_max);

		double fx_init = fx;
		double dg_init = dg;

		double test_decr = param.ftol * dg_init;
		double test_curv = -param.wolfe * dg_init;
		double I_lo = 0.0, I_hi = Double.POSITIVE_INFINITY;
		double fI_lo = 0.0, fI_hi = Double.POSITIVE_INFINITY;
		double gI_lo = (1.0 - param.ftol) * dg_init, gI_hi = Double.POSITIVE_INFINITY;

		for (int i = 0; i < x.length; i++) {
			x[i] = xp[i] + step * drt[i];
		}

		fx = f.eval(x, grad);
		dg = Vector.dot(grad, drt);

		if (DEBUG) {
			debug("-> before wolfe and loop");
			debug("test_decr: " + test_decr);
			debug("test_curv: " + test_curv);
			debug("        x: ", x);
			debug("       fx: " + fx);
			debug("     grad: ", grad);
			debug("       dg: " + dg);
			debug(
					"wolfe cond 1: " + fx + " <= " + (fx_init + step * test_decr) + " == " + (fx <= fx_init + step * test_decr));
			debug("wolfe cond 2: " + Math.abs(dg) + " <= " + test_curv + " == " + (Math.abs(dg) <= test_curv));
		}

		if ((fx <= fx_init + step * test_decr) && (Math.abs(dg) <= test_curv)) {
			if (DEBUG)
				debug('-', "leaving line search, criteria met");
			return;
		}

		if (DEBUG)
			debug("-> entering loop");

		int iter;
		for (iter = 0; iter < param.max_linesearch; iter++) {
			double ft = fx - fx_init - step * test_decr;
			double gt = dg - param.ftol * dg_init;

			if (DEBUG) {
				debug("iter: " + iter);
				debug("  ft: " + ft);
				debug("  gt: " + gt);
			}

			double new_step;
			if (ft > fI_lo) {
				new_step = step_selection(I_lo, I_hi, step, fI_lo, fI_hi, ft, gI_lo, gI_hi, gt);

				if (new_step <= param.min_step)
					new_step = (I_lo + step) / 2.0;

				I_hi = step;
				fI_hi = ft;
				gI_hi = gt;

				if (DEBUG) {
					debug("-- case 1, " + ft + " > " + fI_lo);
					debug("-- new_step: " + new_step);
				}
			} else if (gt * (I_lo - step) > 0.0) {
				new_step = Math.min(step_max, step + delta * (step - I_lo));

				I_lo = step;
				fI_lo = ft;
				gI_lo = gt;

				if (DEBUG) {
					debug("-- case 2, " + (gt * (I_lo - step)) + " > 0.0");
					debug("-- new_step: " + new_step);
				}
			} else {
				new_step = step_selection(I_lo, I_hi, step, fI_lo, fI_hi, ft, gI_lo, gI_hi, gt);

				I_hi = I_lo;
				fI_hi = fI_lo;
				gI_hi = gI_lo;

				I_lo = step;
				fI_lo = ft;
				gI_lo = gt;

				if (DEBUG) {
					debug("-- case 3");
					debug("-- new_step: " + new_step);
				}
			}

			if (step == step_max && new_step >= step_max) {
				if (DEBUG)
					debug('-', "leaving line search, maximum step size reached");
				return;
			}

			step = new_step;

			if (step < param.min_step)
				throw new LBFGSBException("the line search step became smaller than the minimum value allowed");

			if (step > param.max_step)
				throw new LBFGSBException("the line search step became larger than the maximum value allowed");

			for (int i = 0; i < x.length; i++) {
				x[i] = xp[i] + step * drt[i];
			}
			fx = f.eval(x, grad);
			dg = Vector.dot(grad, drt);

			if (DEBUG) {
				debug("     x: ", x);
				debug("    fx: " + fx);
				debug("  grad: ", grad);
				debug("    dg: " + dg);
				debug("  wolfe cond 1: " + fx + " <= " + (fx_init + step * test_decr) + " == "
						+ (fx <= fx_init + step * test_decr));
				debug("  wolfe cond 2: " + Math.abs(dg) + " <= " + test_curv + " == " + (Math.abs(dg) <= test_curv));
			}

			if ((fx <= fx_init + step * test_decr) && (Math.abs(dg) <= test_curv)) {
				if (DEBUG)
					debug('-', "leaving line search, criteria met (2)");
				return;
			}

			if (step >= step_max) {
				double ft2 = fx - fx_init - step * test_decr;
				if (ft2 <= fI_lo) {
					if (DEBUG)
						debug('-', "leaving line search, maximum step size reached (2)");
					return;
				}
			}
		}

		throw new LBFGSBException("the line search routine reached the maximum number of iterations");
	}
}
