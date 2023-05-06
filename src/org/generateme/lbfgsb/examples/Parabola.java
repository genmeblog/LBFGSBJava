package org.generateme.lbfgsb.examples;

import org.generateme.lbfgsb.IGradFunction;

public class Parabola implements IGradFunction {

	@Override
	public double evaluate(double[] x, double[] grad) {
		double xx = x[0];
		
		grad[0] = 4*xx-1;
		
		return 2*xx*xx-xx+3;
	}

}
