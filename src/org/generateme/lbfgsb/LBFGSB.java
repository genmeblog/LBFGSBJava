package org.generateme.lbfgsb;

import static org.generateme.lbfgsb.Debug.*;

import org.generateme.lbfgsb.linesearch.*;

public final class LBFGSB {
	public Parameters m_param;
	public BFGSMat m_bfgs;
	public double[] m_fx;
	public double[] m_xp;
	public double[] m_grad;
	public double[] m_gradp;
	public double[] m_drt;
	public double m_projgnorm;
	public double fx;
	public int k;
		
	public void reset(int n) {
		m_bfgs.reset(n, m_param.m);
		m_xp = Vector.resize(m_xp, n);
		m_grad = Vector.resize(m_grad,n);
		m_gradp = Vector.resize(m_gradp, n);
		m_drt = Vector.resize(m_drt, n);
		if(m_param.past > 0) {
			m_fx = Vector.resize(m_fx, m_param.past);
		}
	}
	
	public LBFGSB() { this(new Parameters()); }
	public LBFGSB(Parameters param) {
		this.m_param = param;
		this.m_bfgs = new BFGSMat();
	}
	
	public final static void force_bounds(double[] x, double[] lb, double[] ub) {
		for(int i=0;i<x.length;i++) {
			x[i] = Math.max(Math.min(x[i], ub[i]), lb[i]);
		}
	}
	
	public final static double proj_grad_norm(double[] x, double[] g, double[] lb, double[] ub) {
		double res = 0.0;
		
		for(int i=0; i<x.length; i++) {
			double proj = Math.max(Math.min(x[i]-g[i], ub[i]), lb[i]);
			res = Math.max(res, Math.abs(proj-x[i]));
		}
		
		return res;
	}
	
	public final static double max_step_size(double[] x0, double[] drt, double[] lb, double[] ub) {
		double step = Double.POSITIVE_INFINITY;
		
		for(int i=0; i<x0.length; i++) {
			if(drt[i] > 0.0) {
				step = Math.min(step, (ub[i]-x0[i])/drt[i]);
			} else if (drt[i] < 0) {
				step = Math.min(step, (lb[i]-x0[i])/drt[i]);
			}
		}
		
		return step;
	}
	
	public static final double eps = Math.ulp(1.0);
	
	// f -function, x - initial and final value, lb - lower bounds, ub - upper bounds 
	public double[] minimize(IGradFunction f, double[] in, double[] lb, double[] ub) throws LBFGSBException {
		if(DEBUG) debug('=', "entering minimization");
		
		double[] x = in.clone();
		
		int n = x.length;
		
		if(lb.length != n || ub.length != n)
			throw new LBFGSBException("lb and ub must have the same size as x");
		
		force_bounds(x,lb,ub);
		
		reset(n);
		int fpast = m_param.past;
		
		fx = f.eval(x, m_grad);
		m_projgnorm = proj_grad_norm(x, m_grad, lb, ub);
		
		if(DEBUG) debug("initial");
		if(DEBUG) debug("  fx:        " + fx);
		if(DEBUG) debug("  projgnorm: " + m_projgnorm);
		if(DEBUG) debug("  x:         ", x);
		if(DEBUG) debug("  grad:      ", m_grad);
		
		if(fpast > 0) {
			m_fx[0] = fx;
		}
	
		if(m_projgnorm <= m_param.epsilon || m_projgnorm <= m_param.epsilon_rel * Vector.norm(x)) {
			if(DEBUG) debug('=', "leaving minimization, projgnorm less than epsilon, projgnorm = " + m_projgnorm);
			return x;
		}
		
		// Cauchy stores xcp, vecc, newact_set and fv_set
		Cauchy cauchy = new Cauchy(m_bfgs, x, m_grad, lb, ub);
		
		Vector.sub(cauchy.xcp, x, m_drt);
//		Vector.normalize(m_drt);
		
		double[] vecs = new double[n];
		double[] vecy = new double[n];
		
		k = 1;
	
		for(;;) {
			if(DEBUG) debug('#', "K = " + k);
			
			m_xp = x.clone();
			m_gradp = m_grad.clone();
			
			double dg = Vector.dot(m_grad, m_drt);
			double step_max = max_step_size(x, m_drt, lb, ub);
			
			if(dg >= 0.0 || step_max <= m_param.min_step) {
				Vector.sub(cauchy.xcp, x, m_drt);
				m_bfgs.reset(n, m_param.m);
				dg = Vector.dot(m_grad, m_drt);
				step_max = max_step_size(x, m_drt, lb, ub);
			}
			
			step_max = Math.min(step_max, m_param.max_step);
			double step = Math.min(1.0, step_max);
			
			AbstractLineSearch ls;
			switch(m_param.linesearch) {
			case MORETHUENTE_LBFGSPP:
				ls = new LineSearch(f, m_param, m_xp, m_drt, step_max, step, fx, m_grad, dg, x);
				break;
			case MORETHUENTE_ORIG_STRONG:
				ls = new MoreThuente(f, m_param, m_xp, m_drt, step_max, step, fx, m_grad, dg, x, false);
				break;
			case LEWISOVERTON:
				ls = new LewisOverton(f, m_param, m_xp, m_drt, step_max, step, fx, m_grad, dg, x);
				break;
			case MORETHUENTE_ORIG_WEAK:
			default:
				ls = new MoreThuente(f, m_param, m_xp, m_drt, step_max, step, fx, m_grad, dg, x, true);
				break;
			}
			
			fx = ls.get_fx();
			step = ls.get_step();
			dg = ls.get_dg();
			
			m_projgnorm = proj_grad_norm(x,m_grad,lb,ub);

			if(DEBUG) debug("  fx:        " + fx);
			if(DEBUG) debug("  projgnorm: " + m_projgnorm);
			if(DEBUG) debug("  x:         ", x);
			if(DEBUG) debug("  grad:      ", m_grad);
			
			if(m_projgnorm <= m_param.epsilon || m_projgnorm <= m_param.epsilon_rel * Vector.norm(x)) {
				if(DEBUG) debug('=', "leaving minimization, projgnorm less than epsilons, projgnorm = " + m_projgnorm);
				return x;
			}
		
			if(fpast > 0 ) {
				double fxd = this.m_fx[k % fpast];
				if( k>=fpast && Math.abs(fxd - fx) <= m_param.delta * Math.max(Math.max(Math.abs(fx), Math.abs(fxd)), 1.0) ) {
					if(DEBUG) debug('=', "leaving minimization, past results less than delta");
					return x;
				}
				this.m_fx[k % fpast] = fx;
			}
			
			if(m_param.max_iterations != 0 && k>=m_param.max_iterations) {
				if(DEBUG) debug('=', "leaving minimization, max iterations reached");
				return x;
			}
		
			Vector.sub(x, m_xp, vecs);
			Vector.sub(m_grad, m_gradp, vecy);

			if(Vector.dot(vecs, vecy) > eps * Vector.squaredNorm(vecy)) {
				m_bfgs.add_correction(vecs, vecy);
			}
			
			force_bounds(x,lb,ub);			
			cauchy = new Cauchy(m_bfgs, x, m_grad, lb, ub);

			SubspaceMin.subspace_minimize(m_bfgs, x, m_grad, lb, ub, cauchy, m_param.max_submin, m_drt);
			
			k++;
		}

	}
}
