package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// Parabola, f(x)=2x^2-x+3
// Global minimum: f(0.25) = 2.875
public class Parabola implements IGradFunction {

	public double evaluate(double[] x, double[] grad) {
		double xx = x[0];
		grad[0] = 4 * xx - 1;
		return 2 * xx * xx - xx + 3;
	}

	public boolean in_place_gradient() {
		return true;
	}

	public static void main(String[] args) {

		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		// converges to global minimum
		try {
			double[] res = lbfgsb.minimize(new Parabola(), new double[] { -2 }, new double[] { -5 }, new double[] { 5 });
			debug('!', "RESULT");
			debug("k = " + lbfgsb.k);
			debug("x = ", res);
			debug("fx = " + lbfgsb.fx);
			debug("grad = ", lbfgsb.m_grad);
		} catch (LBFGSBException e) {
			e.printStackTrace();
		}
	}

}
