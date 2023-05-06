package org.generateme.lbfgsb.examples;

import org.generateme.lbfgsb.IGradFunction;

public class Rosenbrock implements IGradFunction {
	private int n;
	
	public Rosenbrock(int n) {
		this.n = n;
	}

	public double evaluate(double[] x, double[] grad) {
	    double fx = (x[0] - 1.0) * (x[0] - 1.0);
        grad[0] = 2 * (x[0] - 1) + 16 * (x[0] * x[0] - x[1]) * x[0];
        for(int i = 1; i < n; i++)
        {
            fx += 4 * Math.pow(x[i] - x[i - 1] * x[i - 1], 2);
            if(i == n - 1)
            {
                grad[i] = 8 * (x[i] - x[i - 1] * x[i - 1]);
            } else {
                grad[i] = 8 * (x[i] - x[i - 1] * x[i - 1]) + 16 * (x[i] * x[i] - x[i + 1]) * x[i];
            }
        }
        return fx;
	}

}
