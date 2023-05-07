package org.generateme.lbfgsb.examples;

import static org.generateme.lbfgsb.Debug.debug;

import org.generateme.lbfgsb.Debug;
import org.generateme.lbfgsb.IGradFunction;
import org.generateme.lbfgsb.LBFGSB;
import org.generateme.lbfgsb.LBFGSBException;
import org.generateme.lbfgsb.Parameters;

// BUKIN FUNCTION N. 6
// https://www.sfu.ca/~ssurjano/bukin6.html
// f(-10,1) = 0;
public class Buckin implements IGradFunction {
  public double evaluate(double[] in) {
  	double x1 = in[0];
  	double x2 = in[1];
  	
  	return 100.0 * Math.sqrt(Math.abs(x2 - 0.01*x1*x1)) + 0.01 * Math.abs(x1 + 10.0);
  }

  public static void main(String[] args) {

  	Debug.DEBUG = true;
  	
		Parameters param = new Parameters();
		LBFGSB lbfgsb = new LBFGSB(param);

		// FAILS
		
		try {
			double[] res = lbfgsb.minimize(new Buckin(),
					new double[] { -10.1, 1 },
					new double[] { -15, -3},
					new double[] { -5, 3});

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

