package hlalign.base;

public class IndelModel {
	
	/** array containing the parameters of the model */
	public double[] params;
	
	/** number of states in pair HMM representation of model */
	public int states;
	
	/** transition matrix between states */
	public double[][] trans;
	
	/** contains a row for each state
	 *  X is ancestor, Y descendant
	 *  X  Y  val
	 *  1  1   3
	 *  1  0   2
	 *  0  1   1
	 *  0  0   0*/
	public int[][] emit;
	
	/** contains values in 'val' column above */
	public int[] emitVal;
	
	/** For heterogeneous models, gives rate class for each state */
	public int[] type;
	
	/** calculate @trans from @params 
	 * should be overridden by specific model types*/
	public void calcTrans(double t){}
	
	public void setParams(double[] pars){
		params = pars;
	}
}
