package org.generateme.lbfgsb;

public interface IGradFunction {
	default double evaluate(double[] x, double[] grad) {
		return Double.NaN;
	}
	
	default double evaluate(double[] x) {
		return Double.NaN;
	}
	
	default void gradient(double[] x, double[] grad) {
		gradient(x,grad,1.0e-4);
	}
	
	// finite difference, symmetrical gradient, stores result in grad[]
	default void gradient(double[] x, double[] grad, double eps) {
		int n = grad.length;
		
		for(int i=0;i<n;i++) {
			double tmp = x[i];
			double x1 = tmp-eps;
			double x2 = tmp+eps;
			x[i] = x1;
			double y1 = evaluate(x);
			x[i] = x2;
			double y2 = evaluate(x);
			x[i] = tmp; // restore
			grad[i] = (y2-y1)/(2.0*eps);
		}
	}
	
	default boolean in_place_gradient() {
		return false;
	}

	default double eval(double[] x, double[] grad) {
		if(this.in_place_gradient()) {
			return this.evaluate(x, grad);
		} else {
			this.gradient(x,grad);
			return this.evaluate(x);	
		}
	}
	
}
