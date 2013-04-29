package hlalign.base;

public class IndelModel {
	
	public HMM two;
	public HMM three;
	public HMM four;
	public int types;
	
	/** array containing the parameters of the model */
	public double[] params;

	public void setParams(double[] pars){
		params = pars;
	}
	
	public class HMM {

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

		/** True for silent states other than start and end */
		public boolean[] silent;
		
		public HMM(int[][] em){
			emit = em;
			
			states = emit.length;
			
			emitVal = new int[states];
			for(int i = 0; i < states; i++)
				emitVal[i] = 2*emit[i][0] + emit[i][1];			
		
			/** only 1 type in TKF91 */
			type = new int[states];
			
			silent = new boolean[states];
			for(int i = 0; i < states; i++)		// not Start, not End, and nothing emitted
				silent[i] = (i != 0) & (i != (states-1)) & (emitVal[i] == 0);
			
		}

		/** calculate trans from params 
		 * should be overridden by specific model types*/
		public void calcTrans(double t){}

		/** calculate trans from params 
		 * should be overridden by specific model types*/
		public void calcTrans(double t1, double t2){}

		/** calculate trans from params 
		 * should be overridden by specific model types*/
		public void calcTrans(double t1, double t2, double t3){}

	}

}
