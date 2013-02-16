package hlalign.base;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.NormalDistribution;

public class OU {
	/** rate of structural drift */
	double sigma2;
	
	/** marginal variance, defined as sigma2 / (2 * theta)
	 * where theta is mean reversion coefficient */
	double tau;
	
	/** measurement error - structural differences
	 * not explained by diffusion process */
	double epsilon;
	
	public OU(double s, double t, double e){
		sigma2 = s;
		tau = t;
		epsilon = e;
	}

}
