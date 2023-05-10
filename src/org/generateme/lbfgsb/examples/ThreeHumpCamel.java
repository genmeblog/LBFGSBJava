package org.generateme.lbfgsb.examples;

import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;

// THREE-HUMP CAMEL FUNCTION
// https://www.sfu.ca/~ssurjano/camel3.html
// Global minimum: f(0,0)=0 
public final class ThreeHumpCamel implements IGradFunction {

	public double evaluate(double[] in, double[] grad) {
		double x = in[0];
		double x2 = x * x;
		double x4 = x2 * x2;
		double y = in[1];

		grad[0] = 4 * x - 4.20 * x2 * x + x4 * x + y;
		grad[1] = x + 2 * y;

		return 2 * x2 - 1.05 * x4 + x4 * x2 / 6.0 + x * y + y * y;
	}

	public boolean in_place_gradient() {
		return true;
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;

		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		// converges to global minimum
		try {
//			double[] res = lbfgsb.minimize(new ThreeHumpCamel(), new double[] { -5, 5 }, new double[] { -5, -5 },
//					new double[] { 5, 5 });
			double[] res = lbfgsb.minimize(new ThreeHumpCamel(), new double[] { 2, 2 }, new double[] { 0, -5 },
					new double[] { 1, 1 });

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
