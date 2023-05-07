package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;
import org.generateme.lbfgsb.IGradFunction;

// CROSS-IN-TRAY
// https://www.sfu.ca/~ssurjano/crossit.html
// f(+/-1.3491, +/-1.3491) = -2.06261
public class CrossInTray implements IGradFunction {

	public double evaluate(double[] in) {
		double x1 = in[0];
		double x2 = in[1];

		return -0.0001 * Math.pow(
				Math.abs(Math.sin(x1) * Math.sin(x2) * Math.exp(Math.abs(100.0 - (Math.sqrt(x1 * x1 + x2 * x2) / Math.PI))))
						+ 1.0,
				0.1);
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		try {
			double[] res = lbfgsb.minimize(new CrossInTray(), new double[] { 3, -3 }, new double[] { -10, -10 },
					new double[] { 10, 10 });

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
