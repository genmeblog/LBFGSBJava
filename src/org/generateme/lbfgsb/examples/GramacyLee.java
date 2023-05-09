package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// Gramacy & Lee (2012)
// https://www.sfu.ca/~ssurjano/grlee12.html
public class GramacyLee implements IGradFunction {
	public double evaluate(double[] in) {
		double x = in[0];

		double v = (x-1)*(x-1);
		return Math.sin(10.0*Math.PI*x)/(2.0*x)+v*v;
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		param.max_linesearch = 1000;
		LBFGSB lbfgsb = new LBFGSB(param);

		// a lot of fails in line search
		
		try {
			double[] res = lbfgsb.minimize(new GramacyLee(), new double[] { 2.51 }, new double[] { 0.5 },
					new double[] { 2.5 });

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
