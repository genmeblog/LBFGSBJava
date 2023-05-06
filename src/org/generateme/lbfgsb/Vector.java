package org.generateme.lbfgsb;

public final class Vector {
	
	// Eigen resize functionality
	public static final double[] resize(double[] a, int n) {
		 if(a == null) {
			 return new double[n];
		 } else {
			 if(a.length == n) {
				 return a;
			 } else {
				 double[] res = new double[n];
				 System.arraycopy(a, 0, res, 0, n);
				 return res;
			 }
		 }
	}
	
	public static final double dot(double[] a, double[] b) {
		double res = 0.0;
		for(int i=0;i<a.length; i++) {
			res += a[i] * b[i];
		}
		return res;
	}
	
	public static final double squaredNorm(double[] x) {
		double res = 0.0;
		for(int i=0;i<x.length;i++) {
			res += x[i]*x[i];
		}
		return res;
	}
	
	public static final double norm(double[] x) {
		return Math.sqrt(squaredNorm(x));
	}
	
	public static final void normalize(double[] x) {
		double n = norm(x);
		for(int i=0;i<x.length;i++) {
			x[i] /= n;
		}
	}	
	
	public static final void setAll(double[] x, double v) {
		for(int i=0;i<x.length;i++) {
			x[i] = v;
		}
	}
	
	public static final void sub(double[] a, double[] b, double[] res) {
		for(int i=0;i<res.length;i++) {
			res[i] = a[i]-b[i];
		}
	}
}
