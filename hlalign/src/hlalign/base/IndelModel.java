package hlalign.base;

public class IndelModel {
	
	/** array containing the parameters of the model */
	double[] params;
	
	/** number of states in pair HMM representation of model */
	int states;
	
	/** transition matrix between states */
	double[][] trans;
	
	/** contains a column for each state
	 *  a column of 1 0 means that the state emits to the ancestor protein 
	 *  but not the descendant*/
	int[][] emitPattern;
	
	/** calculate @trans from @params */
	public void calcTrans(double t){}
	
}
