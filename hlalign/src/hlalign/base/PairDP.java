package hlalign.base;

import java.util.*;

/** Class for performing dynamic programming for a pair HMM
 * intended to be general to any pair HMM - transition probabilities
 * in an IndelModel object and emission probability arrays are required 
 * in the constructor
 * @author Challis
 */
public class PairDP {
	
	/** Dynamic programming matrix */
	public double [][][] DPmat;
	
	/** HMM governing transitions */
	IndelModel.HMM hmm;
	
	/** Joint probabilities of each pairwise combination of objects being 
	 * aligned (could be individual proteins or entire subtrees)
	 */
	double[][] logJointProbs;
	
	/** Marginal probabilities for object X (ancestor) */
	double[] logMargProbsX;
	
	/** Marginal probabilities for object Y (descendant) */
	double[] logMargProbsY;
	
	/** lengths of X and Y and number of states in pair HMM */
	int lx, ly, S;
		
	public boolean forwardRun;
	
	public PairDP(double[] pX, double[] pY, double[][] pXY, IndelModel.HMM h){
		logMargProbsX = pX;
		logMargProbsY = pY;
		logJointProbs = pXY;
		hmm = h;
		lx = pX.length;
		ly = pY.length;
		S = hmm.trans.length;
		forwardRun = false;
		DPmat = new double[lx][ly][S];

	}
	
	public void forward(){
		
		// initialize DP matrix
		for(int i = 0; i < lx; i++)
			for(int j = 0; j < ly; j++)
				for(int s = 0; s < S; s++)
					DPmat[i][j][s] = Utils.log0;
		
		// corresponds to start state
		DPmat[1][1][0] = 0;

		//DP algorithm
		for(int i = 1; i < lx; i++){
		  for(int j = 1 + ((i==1) ? 1 : 0); j < ly; j++){
		    for(int s = 0; s < S; s++){
		    	for(int t = 0; t < S; t++)					// position to look back to depends on emission pattern of state in HMM
		    		DPmat[i][j][s] = Utils.logAdd( DPmat[i][j][s], hmm.trans[t][s] + DPmat[i - hmm.emit[s][0]] [j - hmm.emit[s][1]] [t] );		    

		    	DPmat[i][j][s] += emit(i, j, s);
		    }
		  }
		}	
		forwardRun = true;
	}
	
	/** Returns marginal probability over all alignments
	 * @return ml
	 */
	public double getMarginalLikelihood(){
		
		if(!forwardRun)
			forward();
		
		double ml = Utils.log0;
		for(int s = 1; s < S; s++)  // sum starts at 1 to skip Start state
			ml = Utils.logAdd(ml, hmm.trans[s][S-1] + DPmat[lx-1][ly-1][s]);  // S-1 is always End state
		return ml;
	}
	
	public double emit(int i, int j, int s){
		double logp;
		switch (hmm.emitVal[s]){
		case 0:
			logp = 0;
			if(hmm.silent[s])
				logp = -Math.log(1 - Math.exp(hmm.trans[s][s]));
			break;
		case 1:
			logp = logMargProbsY[j];
			break;
		case 2:
			logp = logMargProbsX[i];
			break;
		case 3:
			logp = logJointProbs[i][j];
			break;
		default:
			logp = Utils.log0;
		}
		return logp;
	}
	
	public int[] backward(){
		
		if(!forwardRun)
			forward();
		
		double[] probs = new double[S];
		
		// state vector will always begin with match state corresponding to the start state
		int i = lx - 1;
		int j = ly - 1;
		ArrayList<Integer> backAlign = new ArrayList<Integer>();
		
		double ml = getMarginalLikelihood();
		for(int k = 0; k < S; k++)
			probs[k] = hmm.trans[k][S-1] + DPmat[lx-1][ly-1][k] - ml;
		

		int s = Utils.sample(probs);
		backAlign.add(S-1);
		backAlign.add(s);
		
		while( (i > 1) || (j > 1) ){
		//	System.out.println("State: " + s);
			for(int k = 0; k < S; k++)
				probs[k] = hmm.trans[k][s] + (DPmat[i - hmm.emit[s][0]] [j - hmm.emit[s][1]] [k]) - (DPmat[i][j][s] - emit(i, j, s));
			i -= hmm.emit[s][0];
			j -= hmm.emit[s][1];
			s = Utils.sample(probs);
			backAlign.add(s);
		}
		
		int[] align = new int[backAlign.size()];
		for(int k = 0; k < align.length; k++)
			align[k] = backAlign.get(align.length - k - 1);
		
		return align;
			
	}
	
	public void forward2(){
		
		// initialize DP matrix
		for(int i = 0; i < lx; i++)
			for(int j = 0; j < ly; j++)
				for(int s = 0; s < S; s++)
					DPmat[i][j][s] = Utils.log0;
		
		// corresponds to start state
		DPmat[1][1][0] = 0;

		//DP algorithm
		for(int i = 1; i < lx; i++){
		  for(int j = 1 + ((i==1) ? 1 : 0); j < ly; j++){
		    for(int s = 0; s < S-1; s++){
		    	if(!hmm.silent[s]){
		    		for(int t = 0; t < S-1; t++)					// position to look back to depends on emission pattern of state in HMM
		    			DPmat[i][j][s] = Utils.logAdd( DPmat[i][j][s], 
		    					hmm.trans[t][s] + DPmat[i - hmm.emit[s][0]] [j - hmm.emit[s][1]] [t] );		    
		    		DPmat[i][j][s] += emit(i, j, s);
		    	}
		    }
		    for(int s = 0; s < S-1; s++){
		    	if(hmm.silent[s]){
		    		for(int t = 0; t < S-1; t++)
		    			DPmat[i][j][s] =  Utils.logAdd( DPmat[i][j][s], 
		    					hmm.trans[t][s] + DPmat[i - hmm.emit[s][0]] [j - hmm.emit[s][1]] [t] );
		    		DPmat[i][j][s] = DPmat[i][j][s] - Math.log(1 - Math.exp(hmm.trans[s][s]));
		    	}
		    }
		  }
		}	
		forwardRun = true;
	}
	
}
