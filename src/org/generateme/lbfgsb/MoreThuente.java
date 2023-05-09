package org.generateme.lbfgsb;

import static org.generateme.lbfgsb.Debug.DEBUG;
import static org.generateme.lbfgsb.Debug.debug;

// https://github.com/JuliaNLSolvers/LineSearches.jl/blob/master/src/morethuente.jl
// https://github.com/SurajGupta/r-source/blob/master/src/appl/lbfgsb.c#L2976

public class MoreThuente {

	static public enum RESULT {
		NONE, CONVERGED, OUT_RANGE, XTOL, STPMAX, STPMIN, MAX_ITERS
	}

	static public class Bool {
		boolean b;

		public Bool() {
			this(false);
		}

		public Bool(boolean b) {
			this.b = b;
		}
	}

	static public class CStep {
		public double stx, fx, dx, sty, fy, dy, stp, fp, dp;
		public boolean bracketed;
	}

	static public class PhiDPhi {
		IGradFunction f;
		double[] drt; // search direction (cauchy point - x)
		double[] xp; // previous x
		double[] x; // current x
		double[] grad; // gradient
		public double dg; // dot(grad,drt)

		public PhiDPhi(IGradFunction f, double[] x, double[] grad, double[] xp, double[] drt) {
			this.f = f;
			this.x = x;
			this.grad = grad;
			this.xp = xp;
			this.drt = drt;
		}

		public double evaluate(double alpha) {
			for (int i = 0; i < x.length; i++) {
				x[i] = xp[i] + alpha * drt[i];
			}

			double f = this.f.eval(x, grad);
			dg = Vector.dot(grad, drt);
			return f;
		}
	}

	public double fx;
	public double step;
	public double dg;
	public RESULT info;

	public static final double eps = Math.ulp(1.0);
	public static final int iterfinitemax = (int) (-(Math.log(eps)) / Math.log(2.0));

	public MoreThuente(IGradFunction fun, Parameters param, double[] xp, double[] drt, double step_max, double _step,
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

		if (dg >= 0.0)
			throw new LBFGSBException("the moving direction does not decrease the objective function value, dg=" + dg);

		double stp = step;
		double stpmin = param.min_step;
		double stpmax = step_max;

		PhiDPhi phidphi = new PhiDPhi(fun, x, grad, xp, drt);

		info = RESULT.NONE;
		Bool bracketed = new Bool();
		boolean stage1 = true;

		double finit = fx;
		double ginit = dg;
		double gtest = param.ftol * dg; // first wolfe condition
		double ctest = -param.wolfe * dg; // curvature test, second wolfe condition
		double width = stpmax - stpmin;
		double width1 = 2.0 * width;

		double stx = 0.0;
		double fx = finit;
		double gx = ginit;
		double sty = 0.0;
		double fy = finit;
		double gy = ginit;

		double stmin = 0.0;
		double stmax = stp + stp * 4.0;

		CStep cs = new CStep();

		double f = phidphi.evaluate(stp);
		dg = phidphi.dg;

		int iterfinite = 0;
		while ((Double.isInfinite(fx) || Double.isInfinite(dg)) && (iterfinite < iterfinitemax)) {
			stp = stp / 2.0;

			f = phidphi.evaluate(stp);
			dg = phidphi.dg;

			stx = stp * 7.0 / 8.0;
		}

		if (DEBUG) {
			debug('>', "entering loop");
			debug("     alpha: " + stp);
			debug("       stx: " + stx);
		}

		int iter = 0;
		while (true) {
			if (DEBUG) {
				debug("  line search iter:" + iter);
			}

			f = phidphi.evaluate(stp);
			dg = phidphi.dg;

			if (Math.abs(dg) < eps) {
				this.step = stp;
				this.fx = f;

				if (DEBUG) {
					debug("    step: " + stp);
					debug("      fx: " + f);
					debug("       x: ", x);
					debug("    grad: ", grad);
					debug('-', "leaving line search, dg = " + dg);
				}
				return;
			}

			double ftest = finit + stp * gtest;

			if (stage1 && f < ftest && dg >= 0.0) {
				stage1 = false;
			}

			if ((bracketed.b && (stp <= stmin || stp >= stmax))) {
				info = RESULT.OUT_RANGE;
			}
			if (stp == stpmax && f <= ftest && dg <= gtest) {
				info = RESULT.STPMAX;
			}
			if (stp == stpmin && (f > ftest || dg >= gtest)) {
				info = RESULT.STPMIN;
			}
			if (iter >= param.max_linesearch) {
				info = RESULT.MAX_ITERS;
			}
			if (bracketed.b && (stmax - stmin) <= param.xtol * stmax) {
				info = RESULT.XTOL;
			}
			if (f <= ftest && Math.abs(dg) <= ctest) {
				info = RESULT.CONVERGED;
			}

			if (info != RESULT.NONE) {
				this.step = stp;
				this.fx = f;

				if (DEBUG) {
					debug("    step: " + stp);
					debug("      fx: " + f);
					debug("       x: ", x);
					debug("    grad: ", grad);
					debug('-', "leaving line search, info = " + info);
				}
				return;
			}

			if (stage1 && f < fx && f > ftest) {
				double fm = f - stp * gtest;
				double fxm = fx - stx * gtest;
				double fym = fy - sty * gtest;
				double dgm = dg - gtest;
				double gxm = gx - gtest;
				double gym = gy - gtest;

				cs.stx = stx;
				cs.fx = fxm;
				cs.dx = gxm;
				cs.sty = sty;
				cs.fy = fym;
				cs.dy = gym;
				cs.stp = stp;
				cs.fp = fm;
				cs.dp = dgm;
				cs.bracketed = bracketed.b;

				cstep(cs, stmin, stmax);

				stx = cs.stx;
				fxm = cs.fx;
				gxm = cs.dx;
				sty = cs.sty;
				fym = cs.fy;
				gym = cs.dy;
				stp = cs.stp;
				fm = cs.fp;
				bracketed.b = cs.bracketed;

				fx = fxm + stx * gtest;
				fy = fym + sty * gtest;
				gx = gxm + gtest;
				gy = gym + gtest;
			} else {
				cs.stx = stx;
				cs.fx = fx;
				cs.dx = gx;
				cs.sty = sty;
				cs.fy = fy;
				cs.dy = gy;
				cs.stp = stp;
				cs.fp = f;
				cs.dp = dg;
				cs.bracketed = bracketed.b;

				cstep(cs, stmin, stmax);

				stx = cs.stx;
				fx = cs.fx;
				gx = cs.dx;
				sty = cs.sty;
				fy = cs.fy;
				gy = cs.dy;
				stp = cs.stp;
				f = cs.fp;
				dg = cs.dp;
				bracketed.b = cs.bracketed;
			}

			if (bracketed.b) {
				double adiff = Math.abs(sty - stx);
				if (adiff >= width1 * 0.66) {
					stp = stx + (sty - stx) / 2.0;
				}
				width1 = width;
				width = adiff;
			}

			if (bracketed.b) {
				stmin = Math.min(stx, sty);
				stmax = Math.max(stx, sty);
			} else {
				stmin = stp + 1.1 * (stp - stx);
				stmax = stp + 4.0 * (stp - stx);
			}

			stp = stp < stpmin ? stpmin : stp;
			stp = stp > stpmax ? stpmax : stp;

			if (DEBUG) {
				debug("  stmin: " + stmin);
				debug("  stmax: " + stmax);
				debug("  alpha: " + stp);
			}

			iter++;

			if ((bracketed.b && (stp <= stmin || stp >= stmax)) || (bracketed.b && ((stmax - stmin) <= (param.xtol * stmax)))
					|| iter >= param.max_linesearch) {
				if (DEBUG)
					debug("  fallback to stx: " + stx);
				stp = stx;
			}

		}
	}

	private void cstep(CStep cs, double stpmin, double stpmax) {
		double sgnd = cs.dp * (cs.dx / Math.abs(cs.dx));

		if (DEBUG) {
			debug("<< cstep");
			debug(" stx=" + cs.stx);
			debug("  fx=" + cs.fx);
			debug("  dx=" + cs.dx);
			debug(" sty=" + cs.sty);
			debug("  fy=" + cs.fy);
			debug("  dy=" + cs.dy);
			debug(" stp=" + cs.stp);
			debug("  fp=" + cs.fp);
			debug("  dp=" + cs.dp);
			debug("  br=" + cs.bracketed);
			debug(" min=" + stpmin);
			debug(" max=" + stpmax);
			debug("sgnd=" + sgnd);
		}

		double stpf;

		if (cs.fp > cs.fx) {

			if (DEBUG)
				debug("= Case 1");

			double theta = 3.0 * (cs.fx - cs.fp) / (cs.stp - cs.stx) + cs.dx + cs.dp;
			double s = Math.max(Math.max(Math.abs(theta), Math.abs(cs.dx)), Math.abs(cs.dp));
			double d1 = theta / s;
			double gamm = s * Math.sqrt(d1 * d1 - (cs.dx / s) * (cs.dp / s));
			if (cs.stp < cs.stx) {
				gamm = -gamm;
			}
			double p = gamm - cs.dx + theta;
			double q = gamm - cs.dx + gamm + cs.dp;
			double r = p / q;
			double stpc = cs.stx + r * (cs.stp - cs.stx);
			double stpq = cs.stx + cs.dx / ((cs.fx - cs.fp) / (cs.stp - cs.stx) + cs.dx) / 2.0 * (cs.stp - cs.stx);
			if (Math.abs(stpc - cs.stx) < Math.abs(stpq - cs.stx)) {
				stpf = stpc;
			} else {
				stpf = (stpc + stpq) / 2.0;
			}
			if (DEBUG)
				debug("= stpf: " + stpf);
			cs.bracketed = true;

		} else if (sgnd < 0.0) {
			if (DEBUG)
				debug("= Case 2");

			double theta = 3.0 * (cs.fx - cs.fp) / (cs.stp - cs.stx) + cs.dx + cs.dp;
			double s = Math.max(Math.max(Math.abs(theta), Math.abs(cs.dx)), Math.abs(cs.dp));
			double d1 = theta / s;
			double gamm = s * Math.sqrt(d1 * d1 - (cs.dx / s) * (cs.dp / s));
			if (cs.stp > cs.stx) {
				gamm = -gamm;
			}
			double p = gamm - cs.dp + theta;
			double q = gamm - cs.dp + gamm + cs.dx;
			double r = p / q;
			double stpc = cs.stp + r * (cs.stx - cs.stp);
			double stpq = cs.stp + (cs.dp / (cs.dp - cs.dx)) * (cs.stx - cs.stp);
			if (Math.abs(stpc - cs.stp) < Math.abs(stpq - cs.stp)) {
				stpf = stpc;
			} else {
				stpf = stpq;
			}
			if (DEBUG)
				debug("= stpf: " + stpf);
			cs.bracketed = true;

		} else if (Math.abs(cs.dp) < Math.abs(cs.dx)) {
			if (DEBUG)
				debug("= Case 3");

			double theta = 3.0 * (cs.fx - cs.fp) / (cs.stp - cs.stx) + cs.dx + cs.dp;
			double s = Math.max(Math.max(Math.abs(theta), Math.abs(cs.dx)), Math.abs(cs.dp));

			double d1 = theta / s;
			d1 = d1 * d1 - (cs.dx / s) * (cs.dp / s);
			double gamm = d1 <= 0.0 ? 0.0 : s * Math.sqrt(d1);

			if (cs.stp > cs.stx) {
				gamm = -gamm;
			}

			double p = gamm - cs.dp + theta;
			double q = gamm + (cs.dx - cs.dp) + gamm;
			double r = p / q;

			double stpc;
			if (r < 0.0 && gamm != 0.0) {
				stpc = cs.stp + r * (cs.stx - cs.stp);
			} else if (cs.stp > cs.stx) {
				stpc = stpmax;
			} else {
				stpc = stpmin;
			}

			double stpq = cs.stp + cs.dp / (cs.dp - cs.dx) * (cs.stx - cs.stp);
			if (cs.bracketed) {
				if (Math.abs(stpc - cs.stp) < Math.abs(stpq - cs.stp)) {
					stpf = stpc;
				} else {
					stpf = stpq;
				}
				d1 = cs.stp + (cs.sty - cs.stp) * 0.66;
				if (cs.stp > cs.stx) {
					stpf = Math.min(d1, stpf);
				} else {
					stpf = Math.max(d1, stpf);
				}
			} else {
				if (Math.abs(stpc - cs.stp) > Math.abs(stpq - cs.stp)) {
					stpf = stpc;
				} else {
					stpf = stpq;
				}
				stpf = Math.min(Math.max(stpmin, stpf), stpmax);
			}
			if (DEBUG)
				debug("= stpf: " + stpf);
		} else {
			if (DEBUG)
				debug("= Case 4");
			if (cs.bracketed) {
				double theta = 3.0 * (cs.fp - cs.fy) / (cs.sty - cs.stp) + cs.dy + cs.dp;
				double s = Math.max(Math.max(Math.abs(theta), Math.abs(cs.dy)), Math.abs(cs.dp));
				double d1 = theta / s;
				double gamm = s * Math.sqrt(d1 * d1 - (cs.dy / s) * (cs.dp / s));
				if (cs.stp > cs.sty) {
					gamm = -gamm;
				}
				double p = cs.stp - cs.dp + theta;
				double q = gamm - cs.dp + gamm + cs.dy;
				double r = p / q;
				stpf = cs.stp + r * (cs.sty - cs.stp);
			} else if (cs.stp > cs.stx) {
				stpf = stpmax;
			} else {
				stpf = stpmin;
			}

			if (DEBUG)
				debug("= stpf: " + stpf);
		}

		if (cs.fp > cs.fx) {
			cs.sty = cs.stp;
			cs.fy = cs.fp;
			cs.dy = cs.dp;
		} else {
			if (sgnd < 0.0) {
				cs.sty = cs.stx;
				cs.fy = cs.fx;
				cs.dy = cs.dx;
			}
			cs.stx = cs.stp;
			cs.fx = cs.fp;
			cs.dx = cs.dp;
		}

		cs.stp = stpf;

		if (DEBUG) {
			debug("= [2]alpha: " + cs.stp);
		}
	}
}
