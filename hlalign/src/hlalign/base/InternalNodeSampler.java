package hlalign.base;

public class InternalNodeSampler {
	// There are now many states to keep track of
	// -x-- vs -x.. vs .x-- etc
	// and for different combinations of rate groups
	
	// Do each neighbor of internal node separately, parent seems easiest
		
	/** Dynamic programming matrix */
	public double [][] DPmat;
	
	/** HMM governing transitions */
	IndelModel indels;
	
	/** length of observations and number of states in pair HMM */
	int lx, S;
	
	/** an array for each state containing the states to which it is equivalent with respect to fixed 'through' alignments */
	int[][] equalStates;
	
	int[] origAlign;
	
	public boolean noEmit;
	public boolean forwardRun;
	
	public int[] parent(int[] alignwx, int[] alignwz, int[] alignxz, int[][] fourAlign){
		
		/* state format for hierarchical links:
		 *  rate class i = 0:I-1
		 *  0: Start
		 *  5i + 1: Match
		 *  5i + 2: Delete
		 *  5i + 3: Insert
		 *  5i + 4: Super Delete
		 *  5i + 5: Super Insert
		 *  5I + 1: Quiet
		 *  5I + 2: End
		*/
		
		// If want to work with some states forbidden, just have transition to them 0 but keep
		// same number of states
		
		/* first need to define equal states in a generic way for hierarchical links
		 * for parent, within each rate class we have:
		 * wx: {wx, w-, -x}
		 * w-: {wx, w-, -x}
		 * -x: {null, -x}
		 * w.: {w.}
		 * .x: {null, .x}
		*/
		
		equalStates[0] = new int[0];
		for(int i = 0; i < indels.types; i++){
				equalStates[5*i + 1] = new int[] {5*i + 1, 5*i + 2, 5*i + 3};
				equalStates[5*i + 2] = new int[] {5*i + 1, 5*i + 2, 5*i + 3};
				equalStates[5*i + 3] = new int[] {5*i + 3};
				equalStates[5*i + 4] = new int[] {5*i + 4};
				equalStates[5*i + 5] = new int[] {5*i + 5};
		}

		
		// initialize DP matrix
		for(int i = 0; i < lx; i++)
			for(int s = 0; s < S; s++)
				DPmat[i][s] = Utils.log0;
		
		
		
		// corresponds to start state
		DPmat[0][0] = 0;

		// DP algorithm
		for(int i = 0; i < fourAlign[0].length; i++){
			if(i > 0){
				for(int s : equalStates[origAlign[i]]){
					if(!indels.two.silent[s]){
						for(int t : equalStates[origAlign[i-1]])			// position to look back to depends on emission pattern of state in HMM
							DPmat[i][s] = Utils.logAdd( DPmat[i][s], indels.two.trans[t][s] + DPmat[i-1][t] );		    
						// DPmat[i][s] += emit(i, s);
					}
				}
			}
			for(int s : equalStates[origAlign[i]]){
				if(indels.two.silent[s]){
					for(int t : equalStates[origAlign[i]])
						DPmat[i][s] =  Utils.logAdd( DPmat[i][s], indels.two.trans[t][s] + DPmat[i][t] );
					DPmat[i][s] = DPmat[i][s] - Math.log(1 - Math.exp(indels.two.trans[s][s]));
				}
			}
		}
		forwardRun = true;
		
		
		
		return new int[5];
	}
}
