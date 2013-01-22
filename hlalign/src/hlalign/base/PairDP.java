package hlalign.base;

import java.util.*;

public class PairDP {
	
	public double [][][] DPmat;
	IndelModel id;
	double[][] logJointProbs;
	double[] logMargProbsX;
	double[] logMargProbsY;
	int lx, ly, S;
		
	public boolean forwardRun;
	
	public PairDP(double[] pX, double[] pY, double[][] pXY, IndelModel indels){
		logMargProbsX = pX;
		logMargProbsY = pY;
		logJointProbs = pXY;
		id = indels;
		lx = pX.length;
		ly = pY.length;
		S = id.trans.length;
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
		    	
		    	for(int t = 0; t < S; t++)
		    		DPmat[i][j][t] = Utils.logAdd( DPmat[i][j][t], id.trans[s][t] + DPmat[i - id.emit[t][0]] [j - id.emit[t][1]] [s] );		    
		    	
		    	DPmat[i][j][s] += emit(i, j, s);
		    	
		    }
		  }
		}	
		
		forwardRun = true;
	}
	/** Returns marginal probability over all alignments
	 * should not be called without first calling forward()
	 * @return
	 */
	public double getMarginalLikelihood(){
		
		if(!forwardRun)
			forward();
		
		double ml = Math.log(0);
		for(int s = 1; s < S; s++)  // sum starts at 1 to skip Start state
			ml = Utils.logAdd(ml, id.trans[s][S-1] + DPmat[lx-1][ly-1][s]);  // S-1 is always End state
		return ml;
	}
	
	public double emit(int i, int j, int s){
		double logp;
		switch (id.emitVal[s]){
		case 0:
			logp = 0;
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
			probs[k] = id.trans[k][S-1] + DPmat[lx-1][ly-1][k] - ml;

		int s = Utils.sample(probs);
		backAlign.add(S-1);
		backAlign.add(s);
		
		while( (i > 1) || (j > 1) ){
			
			for(int k = 0; k < S; k++)
				probs[k] = id.trans[k][s] + DPmat[i - id.emit[s][0]] [j - id.emit[s][1]] [k] - (DPmat[i][j][s] - emit(i, j, s));
			s = Utils.sample(probs);
			i -= id.emit[s][0];
			j -= id.emit[s][1];
			backAlign.add(s);
		}
		
		int[] align = new int[backAlign.size()];
		for(int k = 0; k < align.length; k++)
			align[k] = backAlign.get(align.length - k - 1);
		
		return align;
			
	}
	
}
