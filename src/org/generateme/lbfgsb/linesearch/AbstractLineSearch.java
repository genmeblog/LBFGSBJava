package org.generateme.lbfgsb.linesearch;

public abstract class AbstractLineSearch {
	public double fx;
	public double step;
	public double dg;
	
	public double get_fx() {return fx;}
	public double get_step() {return step;}
	public double get_dg() {return dg;}
}
