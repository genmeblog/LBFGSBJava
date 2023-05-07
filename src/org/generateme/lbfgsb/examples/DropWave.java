package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// DROP-WAVE
// https://www.sfu.ca/~ssurjano/drop.html
// f(0,0)=-1
public class DropWave implements IGradFunction {
	public double evaluate(double[] in) {
		double x1 = in[0];
		double x2 = in[1];

		double v = x1*x1+x2*x2;
		return - (1.0 + Math.cos(12.0*v)) / (0.5 * v + 2.0);
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		try {
			double[] res = lbfgsb.minimize(new DropWave(), new double[] { -0.1, 0.1 }, new double[] { -5, -5 },
					new double[] { 5, 5 });

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
