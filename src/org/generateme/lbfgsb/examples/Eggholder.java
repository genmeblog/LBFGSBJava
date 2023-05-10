package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// EGGHOLDER
// https://www.sfu.ca/~ssurjano/egg.html
// f(512,404.2319) = -959.6407
public class Eggholder implements IGradFunction {
	public double evaluate(double[] in) {
		double x1 = in[0];
		double x2 = in[1];

		double x2_47 = x2 + 47.0;
		return -x2_47 * Math.sin(Math.sqrt(Math.abs(x2_47 + x1 / 2.0))) - x1 * Math.sin(Math.sqrt(Math.abs(x1 - x2_47)));
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		try {
			double[] res = lbfgsb.minimize(new Eggholder(), new double[] { 500, 350 }, new double[] { -512, -512 },
					new double[] { 512, 512 });

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
