package org.generateme.lbfgsb;

public final class Parameters {
	public int m = 6;
	public double epsilon = 1.0e-7;
	public double epsilon_rel = 1.0e-7;
	public int past = 3;
	public double delta = 1.0e-10;
	public int max_iterations = 1000;
	public int max_submin = 10;
	public int max_linesearch = 20;
	public double xtol = 1.0e-8; // MoreThuente
	public double min_step = 1.0e-20;
	public double max_step = 1.0e20;
	public double ftol = 1.0e-4;
	public double wolfe = 0.9;
}
