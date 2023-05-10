package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

public class Michalewicz implements IGradFunction {
	public double m2;
	public int d;
	
	public Michalewicz(int d, double m) {
		this.m2 = 2.0*m;
		this.d = d;
	}
	
	public Michalewicz(int d) {
		this(d, 10.0);
	}
	
	public Michalewicz() {
		this(5);
	}
	
	public double evaluate(double[] in) {	
		double res = 0.0;
		for(int i=1;i<=d;i++) {
			double x = in[i-1];
			res += Math.sin(x)*Math.pow(Math.sin((i*x*x)/Math.PI), m2);
		}
		return -res;
	}
	
	
	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		double pi = Math.PI;
		
		// converges to global minimum
		try {
		//	double[] res = lbfgsb.minimize(new Michalewicz(), new double[] { 2,2,2,2,2 }, new double[] { 0,0,0,0,0 }, new double[] { pi,pi,pi,pi,pi });
			double[] res = lbfgsb.minimize(new Michalewicz(2), new double[] { 2,2}, new double[] { 0,0 }, new double[] { pi,pi });
			
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
