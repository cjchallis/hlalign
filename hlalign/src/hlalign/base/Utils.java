package hlalign.base;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;

public class Utils {
	
	public static final boolean DEBUG = false;
	
	public static double MIN_EDGE_LENGTH = .0001;
	
	public static final double log0 = Double.NEGATIVE_INFINITY;
	
	public static double logAdd(double log_a, double log_b){
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
	// public static Random generator = new Random(1);
	public static RandomGenerator generator = new Well19937c(1);
	
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
	
	/**
	 * During the burnin, the proposalWidthControlVariable for all continuous parameters
	 * is adjusted in order to ensure that the average acceptance rate is between 
	 * MIN_ACCEPTANCE and MAX_ACCEPTANCE where possible. 
	 * This is done by repeatedly multiplying the proposalWidthControlVariable
	 * by SPAN_MULTIPLIER until the acceptance falls within the desired range.
	 */
	public static final double SPAN_MULTIPLIER = 0.7;
	/**
	 * During the burnin, the proposalWidthControlVariable for all McmcMove objects
	 * is adjusted (if <code>McmcMove.autoTune=true</code>) in order to ensure that the average 
	 * acceptance rate is between MIN_ACCEPTANCE and MAX_ACCEPTANCE where possible. 
	 * This is done by repeatedly multiplying the proposalWidthControlVariable
	 * by SPAN_MULTIPLIER until the acceptance falls within the desired range.
	 */
	public static final double MIN_ACCEPTANCE = 0.2;
	// Put the minimum a bit higher than we want it to be
	// because as the parameters converge the acceptance rate
	// typically goes down
	/**
	 * During the burnin, the proposalWidthControlVariable for all continuous parameters
	 * is adjusted in order to ensure that the average acceptance rate is between 
	 * MIN_ACCEPTANCE and MAX_ACCEPTANCE where possible. 
	 * This is done by repeatedly multiplying the proposalWidthControlVariable
	 * by SPAN_MULTIPLIER until the acceptance falls within the desired range.
	 */
	public static final double MAX_ACCEPTANCE = 0.4;
	
	
	/** 
	 * @param x
	 * @param shape
	 * @param rate
	 * @return The unnormalised log density of Gamma(x | shape, rate)
	 */
	public static double logGammaDensity(double x, double shape, double rate) {
       if (x < 0) {
            return Double.NEGATIVE_INFINITY;
       }
       return (shape-1) * FastMath.log(x) - rate * x;
    }
}

