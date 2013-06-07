package hlalign.base;

import java.util.ArrayList;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import hlalign.indels.*;
import hlalign.io.*;
import hlalign.mcmc.*;

// import statalign.model.ext.plugins.StructAlign.MultiNormCholesky;

public class MCMC {
	public int B;
	public int N;
	
	double currLogLike;
	
	public SubstitutionPIPModel subModel;
	DataReader dataReader;
	String[] names;
	public double[][][] coords;
	public char[][] seqs;
	public int[][] seqsInt;
	public double[] sigma2;
	double[][] fullCovar;
	public Tree tree;
	Alignment[] align;
	// int[][] alignArray;
	Structure structure;
	public McmcMove[] moves;
	public boolean burnin = true;
	public ArrayList<Double> logLikeTrace;
	
	public MCMC(int b, int n, DataReader dr, String[] pnames){
		B = b;
		N = n;
		names = pnames;
		sigma2 = new double[b+n];
		//tree = new Tree[b+n];
		align = new Alignment[b+n];
		dataReader = dr;
		coords = dr.coords;
		seqs = dr.seqs;
		seqsInt = dr.seqsInt;
		logLikeTrace = new ArrayList<Double>(0);
		
		for(int i = 0; i < coords.length; i++){
			RealMatrix temp = new Array2DRowRealMatrix(coords[i]);
			RealVector mean = Utils.meanVector(temp);
			for(int j = 0; j < coords[i].length; j++)
				coords[i][j]= temp.getRowVector(j).subtract(mean).toArray();
		}
	}
	
	public void run(){
		initialize();	// set initial tree, alignment, parameters, and calculate log likelihood
		
		createMoves();
		
		deterministicScan();
		
		printSamples();
		
		finishUp();
		
		/*
		System.out.println("LL: " + structure.logLikelihood(tree));
		
		
		double[] pX = structure.marginal(coords[0]);
		double[] pY = structure.marginal(coords[1]);
		double[][] pXY = structure.joint(coords[0], coords[1], 1);
		
		TKF91 tkf91 = new TKF91(new double[]{.03, .033});
		tkf91.two.calcTrans(1);
		
		PairDP dp = new PairDP(pX, pY, pXY, tkf91.two);
		
		dp.forward();
		System.out.println("ML: " + dp.getMarginalLikelihood());
		int[] align = dp.backward();
		System.out.println("Single align:");
		for(int i = 0; i < align.length; i++)
			System.out.print(align[i]);
		
		int[][] alignInds = structure.makeAlignInds(tree.alignArray);
		structure.marginTree(coords, alignInds, 2, true, tree);
		structure.marginTree(coords, alignInds, 2, false, tree);
		structure.jointTree(coords, alignInds, 2, false, tree);
		
		double[] sums = new double[tkf91.three.states];
		tkf91.three.calcTrans(3, 1);
		
		for(int i = 0; i < tkf91.three.states; i++){
			for(int j = 0; j < tkf91.three.states; j++){
				sums[i] += Math.exp(tkf91.three.trans[i][j]);
			}
			System.out.println(sums[i]);
		}
		
		System.out.println("Four:");
		
		tkf91.four.calcTrans(1,2,3);
		// TODO Figure out where these lines should actually go
		tkf91.four.silent[12] = true;
		tkf91.four.emitVal = new int[tkf91.four.states];
		for(int i = 1; i < tkf91.four.states - 3; i++)
			tkf91.four.emitVal[i] = 1;


		System.out.println("TKF92 matrix sums:");
		TKF92 tkf92 = new TKF92(new double[] {.03, .033, .7});
		tkf92.two.calcTrans(1);
		
		sums = new double[tkf92.two.states];
		
		for(int i = 0; i < tkf92.two.states; i++){
			for(int j = 0; j < tkf92.two.states; j++){
				sums[i] += Math.exp(tkf92.two.trans[i][j]);
			}
			System.out.println(sums[i]);
		}
		
		
		for(int i = 0; i < tkf92.two.states; i++){
			for(int j = 0; j < tkf92.two.states; j++){
				System.out.print(Math.exp(tkf92.two.trans[i][j]) + "  ");
			}
			System.out.println();
		}
		*/


	}
	
	public boolean isParamChangeAccepted(double logPropPriorRatio, McmcMove move){
		double u = Math.log(Utils.generator.nextDouble());
		double r = currLogLike - move.getOldll() + logPropPriorRatio;
		return (r > u);
	}
	
	public double getLogLike(){
		return currLogLike;
	}
	public void setLogLike(double ll){
		currLogLike = ll;
	}
	
	public void initialize(){
		subModel = new SubstitutionPIPModel(dataReader.subQ, dataReader.e);
				
		structure = new Structure(coords, .1, 100, .1);
		tree = new Tree(names, coords, seqsInt, this);
		/* calc PIP ml */
		currLogLike = tree.calcML();
	}
	public void createMoves(){
		// 2 moves: eta and zeta
		// 2n - 2 moves: branch lengths
		moves = new McmcMove[2*names.length];
		
		GammaPrior etaPrior = new GammaPrior(1.5, .01);
		GammaProposal etaProp = new GammaProposal(1, 1);
		moves[0] = new EtaMove(this, etaPrior, etaProp, "eta");
		moves[0].proposalWidthControlVariable = 1;
		
		GammaPrior zetaPrior = new GammaPrior(1, 1);
		GammaProposal zetaProp = new GammaProposal(1, 1);
		moves[1] = new ZetaMove(this, zetaPrior, zetaProp, "zeta");
		moves[1].proposalWidthControlVariable = 1;
		
		for(int i = 2; i < moves.length; i++){
			GammaPrior edgePrior = new GammaPrior(1,1);
			GammaProposal edgeProp = new GammaProposal(1,1);
			moves[i] = new EdgeMove(this, i-2, edgePrior, edgeProp, "edge"+(i-2));
		}
			
	}
	
	public void deterministicScan(){
		// burn-in
		for(int i = 0; i < B; i++){
			for(int j = 0; j < moves.length; j++)
				moves[j].move(tree);
			if(i % Utils.CHECK_PROPOSAL_WIDTHS == 0)
				modifyProposalWidths();
		}
		
		for(int i = 0; i < N; i++)
			for(int j = 0; j < moves.length; j++)
				moves[j].move(tree);
	}
	
	public void printSamples(){
		System.out.println("Print samples:");
		System.out.println("N: " + N);
		for(int i = 0; i < moves[0].sample.size(); i++)
			System.out.println(moves[0].sample.get(i));
	}
	
	public void finishUp(){
		System.out.println("Acceptance rates:");
		for(int i = 0; i < moves.length; i++)
			System.out.println(moves[i].name + ": " + moves[i].acceptanceRate() );		
		
	}
	
	public void modifyProposalWidths() {
		for (McmcMove m : moves) {
			if (!m.autoTune) { continue; }
			if (m.proposalCount > Utils.MIN_SAMPLES_FOR_ACC_ESTIMATE) {
				if (m.acceptanceRate() < m.minAcceptance) {
					m.proposalWidthControlVariable *= m.spanMultiplier;
					m.proposalCount = 0;
					m.acceptanceCount = 0;
				}
				else if (m.acceptanceRate() > m.maxAcceptance) {
					m.proposalWidthControlVariable /= m.spanMultiplier;
					m.proposalCount = 0;
					m.acceptanceCount = 0;
				}
			}
		}
	}
}
