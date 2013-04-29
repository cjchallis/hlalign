package hlalign.base;

import java.util.ArrayList;

public class SingleDP {
	
	/** Dynamic programming matrix */
	public double [][] DPmat;
	
	/** HMM governing transitions */
	IndelModel.HMM hmm;
	
	/** Marginal probabilities for object X (ancestor) */
	double[] logMargProbsX;
	
	/** length of observations and number of states in pair HMM */
	int lx, S;
	
	/** an array for each state containing the states to which it is equivalent with respect to fixed 'through' alignments */
	int[][] equalStates;
	
	int[] origAlign;
	
	public boolean noEmit;
	public boolean forwardRun;
	
	public SingleDP(IndelModel.HMM h, int[] oldAlign, boolean silent){	
		hmm = h;
		lx = oldAlign.length;
		S = hmm.trans.length;
		noEmit = silent;
		forwardRun = false;
		DPmat = new double[lx][S];
		origAlign = oldAlign;
		logMargProbsX = new double[lx];
	}
	
	public void forward(){

		
		equalStates = new int[][] {
				{0, 12},
				{1, 12},
				{2, 12},
				{3, 12},
				{4, 5, 12},
				{4, 5, 12},
				{6, 11, 12},
				{7, 11, 12},
				{8, 10, 12},
				{9, 12},
				{8, 10, 12},
				{6, 7, 11, 12},
				{12},
				{13}
		};
		
		// initialize DP matrix
		for(int i = 0; i < lx; i++)
			for(int s = 0; s < S; s++)
				DPmat[i][s] = Utils.log0;
		
		// corresponds to start state
		DPmat[0][0] = 0;

		// DP algorithm
		for(int i = 0; i < lx; i++){
			if(i > 0){
				for(int s : equalStates[origAlign[i]]){
					if(!hmm.silent[s]){
						for(int t : equalStates[origAlign[i-1]])
							DPmat[i][s] = Utils.logAdd( DPmat[i][s], hmm.trans[t][s] + DPmat[i-1][t] );		    
						DPmat[i][s] += emit(i, s);
					}
				}
			}
			for(int s : equalStates[origAlign[i]]){
				if(hmm.silent[s]){
					for(int t : equalStates[origAlign[i]])
						DPmat[i][s] =  Utils.logAdd( DPmat[i][s], hmm.trans[t][s] + DPmat[i][t] );
					DPmat[i][s] = DPmat[i][s] - Math.log(1 - Math.exp(hmm.trans[s][s]));
				}
			}
		}
		forwardRun = true;
		
		System.out.println("DPmat");
		for(int s = 0; s < hmm.states; s++){
			for(int i = 0; i < 5; i++)
				System.out.print(DPmat[i][s]);
			System.out.println();
		}
	}
	

	/** Returns marginal probability over all alignments
	 * @return ml
	 */
	public double getMarginalLikelihood(){
		
		if(!forwardRun)
			System.out.println("Run forward first");
		
		double ml = Utils.log0;
		for(int s = 1; s < S; s++)  // sum starts at 1 to skip Start state
			ml = Utils.logAdd(ml, hmm.trans[s][S-1] + DPmat[lx-1][s]);  // S-1 is always End state
		return ml;
	}
	
	public double emit(int i, int s){
		double logp;
		switch (hmm.emitVal[s]){
		case 0:
			logp = 0;
			if(hmm.silent[s])
				logp = -Math.log(1 - Math.exp(hmm.trans[s][s]));
			break;
		case 1:
			logp = noEmit ? 0 : logMargProbsX[i];
			break;
		default:
			logp = Utils.log0;
		}
		return logp;
	}
	
	public int[] backward(){
		
		if(!forwardRun)
			System.out.println("Run forward first");
		
		double[] probs = new double[S];
		
		// state vector will always begin with match state corresponding to the start state
		int i = lx - 1;
		ArrayList<Integer> backAlign = new ArrayList<Integer>();
		
		double ml = getMarginalLikelihood();
		System.out.println("ML: " + ml);
		for(int k = 0; k < S; k++){
			probs[k] = hmm.trans[k][S-1] + DPmat[lx-1][k] - ml;
			System.out.println("Prob: " + Math.exp(probs[k]));
		}
		

		int s = Utils.sample(probs);
		backAlign.add(S-1);
		System.out.println("Sampled: " + s);
		backAlign.add(s);
		int j = 0;
		int t = 0;
		int inc = 0;
		
		while(i > 0){
			System.out.println("State: " + s);
			j = 0;
			inc = hmm.silent[s] ? 0 : 1;
			probs = new double[equalStates[origAlign[i-1]].length];
			/*System.out.println("DP - emit: " + (DPmat[i][s] - emit(i, s)));
			for(int k = 0; k < hmm.states; k++){
				System.out.println("DP: " + DPmat[i-1][k]);
				System.out.println("Tr: " + hmm.trans[k][s]);
			}*/
			
			for(int k : equalStates[origAlign[i-1]]){
				probs[j] = hmm.trans[k][s] + (DPmat[i - inc] [k]) - (DPmat[i][s] - emit(i, s));
				j++;
			}
			for(int k = 0; k < probs.length; k++)
				System.out.println("Prob: " + Math.exp(probs[k]));
			t = equalStates[origAlign[i-1]][Utils.sample(probs)];
			backAlign.add(t);
			i -= inc;
			s = t;
		}
		
		int[] align = new int[backAlign.size()];
		for(int k = 0; k < align.length; k++)
			align[k] = backAlign.get(align.length - k - 1);
		
		return align;
			
	}
}
