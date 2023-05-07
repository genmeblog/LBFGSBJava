package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// https://www.sfu.ca/~ssurjano/boha.html
public class Bohachevsky2 implements IGradFunction {
	
	public double evaluate(double[] in) {
		double x1 = in[0];
		double x2 = in[1];

		return x1 * x1 + 2.0 * x2 * x2 - 0.3 * Math.cos(3.0 * Math.PI * x1) * Math.cos(4.0 * Math.PI * x2) + 0.3;
	}



	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		param.max_linesearch = 100;
		param.wolfe = 0.95;
		LBFGSB lbfgsb = new LBFGSB(param);

		try {
			double[] res = lbfgsb.minimize(new Bohachevsky2(), new double[] { -100,-100 }, new double[] { -100, -100 },
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
