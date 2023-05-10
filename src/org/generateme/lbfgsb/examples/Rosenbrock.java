package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// ROSENBROCK
// https://www.sfu.ca/~ssurjano/rosen.html
// Global minimum: f(1,1,...,1) = 0;
public class Rosenbrock implements IGradFunction {
	private int n;

	public Rosenbrock(int n) {
		this.n = n;
	}

	public Rosenbrock() {
		this(5);
	}

	public boolean in_place_gradient() {
		return true;
	}

	public double evaluate(double[] x, double[] grad) {
		double fx = (x[0] - 1.0) * (x[0] - 1.0);
		grad[0] = 2.0 * (x[0] - 1.0) + 16.0 * (x[0] * x[0] - x[1]) * x[0];
		for (int i = 1; i < n; i++) {
			double v = x[i] - x[i - 1] * x[i - 1];
			fx += 4.0 * v*v;
			if (i == n - 1) {
				grad[i] = 8.0 * v;
			} else {
				grad[i] = 8.0 * v + 16.0 * (x[i] * x[i] - x[i + 1]) * x[i];
			}
		}
		return fx;
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);
		
		// converges to global minimum
		try {
			double[] res = lbfgsb.minimize(new Rosenbrock(), new double[] { 2, -4, 2, 4, -2 },
					new double[] { -5, -5, -5, -5, -5 }, new double[] { 10, 10, 10, 10, 10 });
			debug('!', "RESULT");
			debug("k = " + lbfgsb.k);
			debug("x = ", res);
			debug("fx = " + lbfgsb.fx);
			debug("grad = ", lbfgsb.m_grad);
		} catch (LBFGSBException e) {
			e.printStackTrace();
		}

		Rosenbrock f = new Rosenbrock();
		double[] g = new double[5];
 		debug("res=" + f.eval(new double[] {2,-4,2,4,-2}, g));
	
	}

}
