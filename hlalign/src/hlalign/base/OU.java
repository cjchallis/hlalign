package hlalign.base;

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
	
	public double[] marginal(double[][] X){
		double[] pX = new double[X.length + 2];
		pX[0] = pX[1] = Utils.log0;
		
		NormalDistribution stationary = new NormalDistribution(0, tau);

		for(int i = 2; i < pX.length; i++)
			for(int j = 0; j < 3; j++)
				pX[i] += stationary.density(X[i-2][j]);
		
		return pX;		    
	}
	
	public double[][] joint(double[][] X, double[][] Y, double t){
		
		int lx = X.length + 2;
		int ly = Y.length + 2;
		
		double[][] pXY = new double[lx][ly];
		
		double var = tau + epsilon, covar = tau * Math.exp(-t * sigma2 / (2 * tau));
		
		double[][] covMat = new double[][] { {var, covar}, {covar, var} };

		MultiNormCholesky joint = new MultiNormCholesky(new double[2], covMat);
		
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < ly; j++)
				pXY[i][j] = Utils.log0;

		for(int i = 2; i < lx; i++)
			for(int j = 0; j < 2; j++)
				pXY[i][j] = Utils.log0;
			
		for(int i = 2; i < lx; i++)
			for(int j = 2; j < ly; j++)
				for(int k = 0; k < 3; k++)
					pXY[i][j] += joint.logDensity(new double[] {X[i-2][k], Y[j-2][k]});
		
		return pXY;

	}
}
