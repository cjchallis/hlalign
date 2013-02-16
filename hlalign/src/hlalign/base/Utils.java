package hlalign.base;

import java.util.*;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Utils {
	
	public static final double log0 = Double.NEGATIVE_INFINITY;
	
	static double logAdd(double log_a, double log_b){
		if (log_a == log0)
			return log_b;
		else if (log_b == log0)
			return log_a;
		else if (log_a > log_b)
			return log_a + Math.log(1 + Math.exp(log_b - log_a));
		else
			return log_b + Math.log(1 + Math.exp(log_a - log_b));
	}
	
	static int sample(double[] logProbs)
	{
	  double u = generator.nextDouble();
	  // System.out.println("Random number: " + u);
	  // for(int i = 0; i < logProbs.length; i++)
	  //	  System.out.println(Math.exp(logProbs[i]));
	  double csum = 0;
	  int s = -1;
	  while(u > csum){
	    s++;  
	    csum += Math.exp(logProbs[s]);
	  }
	  return s;
	}
	
	/**
	 * The random number generator used throughout the program.
	 */
	public static Random generator = new Random(1);
	
	public static int[] copyOf(int[] array) {
		int len = array.length;
		int[] copy = new int[len];
		System.arraycopy(array, 0, copy, 0, len);
		return copy;
	}
	
	public static double[] copyOf(double[] array) {
		int len = array.length;
		double[] copy = new double[len];
		System.arraycopy(array, 0, copy, 0, len);
		return copy;
	}
	
	/**
	 * For an n X 3 coordinate matrix, calculate the 1 X 3 mean vector
	 * @param A - coordinate matrix
	 * @return mean vector
	 */	
	static RealVector meanVector(RealMatrix A){
		RealVector mean = new ArrayRealVector(new double[3]);
		for(int i = 0; i < 3; i ++){
			for(int j = 0; j < A.getColumn(0).length; j++)
				mean.addToEntry(i, A.getEntry(j, i));
			mean.setEntry(i, mean.getEntry(i) / A.getColumn(0).length);
		}
		return mean;
	}
	
}

