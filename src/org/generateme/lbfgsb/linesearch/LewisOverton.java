package org.generateme.lbfgsb.linesearch;

import static org.generateme.lbfgsb.Debug.DEBUG;
import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;
import org.generateme.lbfgsb.Vector;

// https://github.com/ZJU-FAST-Lab/LBFGS-Lite/blob/master/include/lbfgs.hpp#L276
public class LewisOverton extends AbstractLineSearch {
	
	static public enum RESULT {
		NONE, CONVERGED, STPMAX, STPMIN, MAX_ITERS, ZERODG, STEPTOL
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
	
	public RESULT info;
	
	public static final double eps = Math.ulp(1.0);
	public static final int iterfinitemax = (int) (-(Math.log(eps)) / Math.log(2.0));

	private void finish(RESULT info, double f, double stp) {
		this.fx = f;
		this.step = stp;
		this.info = info;
		
		if (DEBUG) {
			debug("    step: " + stp);
			debug("      fx: " + f);
			debug('-', "leaving line search, dg = " + dg);
		}
	}
	
	public LewisOverton(IGradFunction fun, Parameters param, double[] xp, double[] drt, double step_max, double _step,
			double _fx, double[] grad, double _dg, double[] x) throws LBFGSBException {
		if (DEBUG) {
			debug('-', "LewisOverton line search");
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

		info = RESULT.NONE;

		boolean brackt = false;
		boolean touched = false;
		
		double finit = fx;
		double dgtest = param.ftol * dg;
		double dstest = param.wolfe * dg;
		double mu = 0.0;
		double nu = stpmax;
	
		PhiDPhi phidphi = new PhiDPhi(fun, x, grad, xp, drt);
		
		double f = phidphi.evaluate(stp);
		dg = phidphi.dg;

		int iterfinite = 0;
		while ((Double.isInfinite(fx) || Double.isInfinite(dg)) && (iterfinite < iterfinitemax)) {
			stp = stp / 2.0;

			f = phidphi.evaluate(stp);
			dg = phidphi.dg;
		}

		if (DEBUG) {
			debug('>', "entering loop");
			debug("       stp: " + stp);
		}

		int iter = 0;
		while (true) {
			if (DEBUG) {
				debug("  line search iter:" + iter);
			}

			f = phidphi.evaluate(stp);
			dg = phidphi.dg;

			if (Math.abs(dg) < eps) {
				finish(RESULT.ZERODG, f, stp);
				return;
			}

			if(f > finit + stp * dgtest) {
				mu = stp;
				brackt = true;
			} else {
				if(dg < dstest) {
					mu = stp;
				} else {
					finish(RESULT.CONVERGED, f, stp);
					return;
				}
			}
			
			if(iter >= param.max_linesearch) {
				finish(RESULT.MAX_ITERS, f, stp);
				return;
			}

			if(brackt && (nu - mu) < eps * nu) {
				finish(RESULT.STEPTOL, f, stp);
				return;
			}
			
			if(brackt) {
				stp = 0.5 * (mu + nu);
			} else {
				stp *= 2.0;
			}
			
			if(stp < stpmin) {
				finish(RESULT.STPMIN, f, stpmin);
				return;
			} else {
				if(touched) {
					finish(RESULT.STPMAX, f, stpmax);
					return;
				} else {
					touched = true;
					stp = stpmax;
				}
			}
			
			iter++;

		}
	}
}
