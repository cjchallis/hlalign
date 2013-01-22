package hlalign.base;

import java.util.*;

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
	
	static int sample(double[] probs)
	{
	  double u = generator.nextDouble();
	  double csum = 0;
	  int s = -1;
	  while(u > csum){
	    s++;  
	    csum += Math.exp(probs[s]);
	  }
	  return s;
	}
	
	/**
	 * The random number generator used throughout the program.
	 */
	public static Random generator = new Random(1);
	
}

