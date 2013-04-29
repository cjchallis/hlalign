package hlalign.base;

public class SubstitutionModel {
	public double[][] u;	// eigenvector matrix
	public double [] d;		// eigenvalues
	public double [][] v;	// inverse eigenvector matrix

	public double[] pi;

	int n = d.length;

	public char[] alphabet;

	public double[][] calcSubMatrix(double t){
		double exp_dtk;
		double [][] subMatrix = new double[n][n];
		for(int k = 0; k < n; k++){
			exp_dtk = Math.exp(d[k]*t);
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					subMatrix[i][j] += u[i][k] * exp_dtk * v[k][j];
				}
			}
		}
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				subMatrix[i][j] = Math.log(subMatrix[i][j]);
		return subMatrix;
	}
}
