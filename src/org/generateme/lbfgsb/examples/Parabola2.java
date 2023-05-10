package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import java.util.Arrays;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

public class Parabola2 implements IGradFunction {

	public double evaluate(double[] in) {
		double x = in[0];
		double y = in[1];
		double t = x * x + y * y;
		if(x<1.0 || y<1.0) {
			debug('W', "outside " + Arrays.toString(in));
		}
		return  t + Math.sin(t) * Math.sin(t);
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;

		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		// converges to global minimum
		try {
			double[] res = lbfgsb.minimize(new Parabola2(), new double[] { 3, 3 }, new double[] { 1, 1 },
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
