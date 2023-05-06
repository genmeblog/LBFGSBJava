package org.generateme.lbfgsb.examples;

import org.generateme.lbfgsb.IGradFunction;

public final class Camel implements IGradFunction {
	public double evaluate(double[] in) {
		double x = in[0];
		double x2 = x*x;
		double x4 = x2*x2;
		double y = in[1];
		return 2*x2 - 1.05*x4 + x4*x2/6.0 + x*y + y*y;
	}
	public double evaluate(double[] in, double[] grad) {
		double x = in[0];
		double x2 = x*x;
		double x4 = x2*x2;
		double y = in[1];
		
		grad[0] = 4*x-4.20*x2*x+x4*x+y;
		grad[1] = x+2*y;
		
		return 2*x2 - 1.05*x4 + x4*x2/6.0 + x*y + y*y; 
	}
	public void gradient(double[] in, double[] grad) {
		double x = in[0];
		double x2 = x*x;
		double x4 = x2*x2;
		double y = in[1];
		
		grad[0] = 4*x-4.20*x2*x+x4*x+y;
		grad[1] = x+2*y; 	
	}
	public boolean in_place_gradient() {return true;}
	public Camel() {}
}
