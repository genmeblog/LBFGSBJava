package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// https://www.sfu.ca/~ssurjano/booth.html
// f(1,3) = 0
public class Booth implements IGradFunction {
	public double evaluate(double[] in, double[] grad) {
		double x1 = in[0];
		double x2 = in[1];

		double a = x1 + 2.0 * x2 - 7.0;
		double b = 2.0 * x1 + x2 - 5.0;

		grad[0] = 2.0 * a + 4.0 * b;
		grad[1] = 4.0 * a + 2.0 * b;

		return a * a + b * b;
	}

	public boolean in_place_gradient() {
		return true;
	}
	
	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		try {
	//		double[] res = lbfgsb.minimize(new Booth(), new double[] { -10,0.1 }, new double[] { -10, -10 },
	//				new double[] { 10, 10 });
		double[] res = lbfgsb.minimize(new Booth(), new double[] { -1,-10 }, new double[] { -10, 3 },
					new double[] { 10, 3.1 });
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
