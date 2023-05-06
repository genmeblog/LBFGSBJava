package org.generateme.lbfgsb;

import static org.generateme.lbfgsb.Debug.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.generateme.lbfgsb.examples.*;

public class Core {
	
	public static void main(String[] args) {
//		System.out.println("Hello");

	//	Matrix mat = new Matrix(new double[] {1, 2, -3, -6, 10, 2, 11, 3, 7, 20, -3, 3, -5, 0,  2, -6, 7, 0, -20, 13, 10, 20, 2, 13, -1}, 5, 5);
	//	debug("Mat = ", mat);
		
//		BKLDLT bkldlt = new BKLDLT(mat);
		
//		System.out.println(Arrays.toString(bkldlt.data));
		
//		double[] res = bkldlt.solve(new double[] {1,20,10,30,50});
	
//		System.out.println(Arrays.toString(res));
	
//		BFGSMat bfgs = new BFGSMat(5,6);
	
//		bfgs.add_correction(new double[] { 0     ,      0 ,  -0.999987, -0.00357138,  0.00357138},
//				new double[] { 0 , 31.9996, -263.884 , 38.4302, 0.199895});
//		bfgs.add_correction(new double[] {  0        , 0 ,-0.210351, -0.996429,  0.175486},
//				new double[] {  0 , 6.73123, 0.645311, -261.783,  41.2326});
//		bfgs.add_correction(new double[] { 0    ,      0 ,-0.0640728,          0,  0.0260629},
//				new double[] { 0,  2.05033, -7.96426, 0.967843, 0.208503});
		
//		System.out.println(Arrays.toString(bfgs.permMinv.mat));
		
		Param param = new Param();
		param.max_iterations = 1000;
		LBFGSB lbfgsb = new LBFGSB(param);

		try {
			double[] x = new double[] {0,3,3,3,3};
			int k = lbfgsb.minimize(new Rosenbrock(5), x, 
					new double[] {0,1,1,1,Double.NEGATIVE_INFINITY}, 
					new double[] {1,4,4,4,4});
			debug('!', "RESULT");
			debug("k = " + k);
			debug("X = ", x);
			debug("fx = " + lbfgsb.fx);
			debug("grad = ", lbfgsb.m_grad);
		} catch (LBFGSBException e) {
			e.printStackTrace();
		}
	}

}
