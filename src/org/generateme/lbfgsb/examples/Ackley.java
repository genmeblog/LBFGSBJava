package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// ACKLEY
// https://www.sfu.ca/~ssurjano/ackley.html
// Global minimum: f(0,0,...,0) = 0
public class Ackley implements IGradFunction {
	private int d;
	private double a, b, c;

	public double evaluate(double[] in) {
		double sumx2 = 0.0;
		double sumcos = 0.0;
		for (double x : in) {
			sumx2 += x * x;
			sumcos += Math.cos(c * x);
		}

		return -a * Math.exp(-b * Math.sqrt(sumx2 / d)) - Math.exp(sumcos / d) + a + Math.E;

	}

	public Ackley() {
		this(5);
	}

	public Ackley(int d) {
		this(d, 20.0, 0.2, Math.PI * 2.0);
	}

	public Ackley(int d, double a, double b, double c) {
		this.d = d;
		this.a = a;
		this.b = b;
		this.c = c;
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		try {
//			double[] res = lbfgsb.minimize(new Ackley(), new double[] { 0.5, 0.2, -0.7, 0.7, -0.6 },
//					new double[] { -32, -32, -32, -32, -32 }, new double[] { 32, 32, 32, 32, 32 });
//			double[] res = lbfgsb.minimize(new Ackley(), new double[] { 0.65, -0.5, 0.0, 0.5, 0.6 },
//					new double[] { 0,-10,0,0,0 }, new double[] { 32, 32, 32, 32, 32 });
			double[] res = lbfgsb.minimize(new Ackley(10),
					new double[] { 0.65, -0.5, 0.0, 0.5, 0.2, 0.2, 0.2, -0.2, -0.2, -0.5 },
					new double[] { -32, -32, -32, -32, 0, -32, -32, -32, -32, -32 },
					new double[] { 32, 32, 32, 32, 32, 32, 32, 32, 32, 32 });

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
