package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

public class Perm implements IGradFunction {
  private int d;
	private double[][] ij;
  private double[] jbeta;
  
	public Perm(double beta, int d) {
		this.d = d;
		
		ij = new double[d][d];
		jbeta = new double[d];
		
		for(int j=1; j<=d; j++) {
			jbeta[j-1] = j + beta;
			for(int i=1; i<=d; i++) {
				ij[i-1][j-1] = 1.0 / Math.pow(j, i);
			}
		}
	}

	public Perm(int d) { this(1.0, d); }
	public Perm() { this(5); }

	public double evaluate(double[] in) {
		double res = 0.0;
		for(int i=1; i<=d; i++) {
			double resj = 0.0;
			for(int j=1; j<=d; j++) {
				resj += jbeta[j-1] * (Math.pow(in[j-1], i) - ij[i-1][j-1]); 
			}
			res += resj*resj;
		}
		return res;
	}

	public static void main(String[] args) {

		Debug.DEBUG = true;
		
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		int d = 3;
		
		double[] z = new double[] { 0,0,0};
		Perm f = new Perm(d);
		debug("res = " + f.evaluate(z));
		
		// converges to global minimum
		try {
			double[] res = lbfgsb.minimize(new Perm(d), new double[] { 1,-1,-3 }, new double[] { -d,-d,-d  }, new double[] { d,d,d });
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
