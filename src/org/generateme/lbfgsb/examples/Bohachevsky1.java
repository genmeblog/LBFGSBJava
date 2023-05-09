package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// https://www.sfu.ca/~ssurjano/boha.html
public class Bohachevsky1 implements IGradFunction {
	
	public double evaluate(double[] in, double[] g) {
		double x1 = in[0];
		double x2 = in[1];

		g[0] = 2 * x1 + 0.3 * Math.sin(3.0 * Math.PI * x1) * 3.0 * Math.PI;
		g[1] = 4 * x2 + 0.4 * Math.sin(4.0 * Math.PI * x2) * 4.0 * Math.PI;
		return x1 * x1 + 2.0 * x2 * x2 - 0.3 * Math.cos(3.0 * Math.PI * x1) - 0.4 * Math.cos(4.0 * Math.PI * x2) + 0.7;
	}

	public boolean in_place_gradient() {
		return true;
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;
				
		Parameters param = new Parameters();
		param.max_linesearch = 100;
		param.wolfe = 0.95;
		LBFGSB lbfgsb = new LBFGSB(param);

		try {
			double[] res = lbfgsb.minimize(new Bohachevsky1(), new double[] { 0,0.1 }, new double[] { -100, -100 },
					new double[] { 100, 100 });

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
