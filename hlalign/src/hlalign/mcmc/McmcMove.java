package hlalign.mcmc;

import hlalign.base.Utils;
import hlalign.mcmc.ParameterInterface;
import hlalign.mcmc.PriorDistribution;
import hlalign.base.MCMC;

public abstract class McmcMove {

	protected MCMC owner;
	public MCMC getOwner(){
		return owner;
	}
	
	protected PriorDistribution<? extends Object> prior;
	
	protected ParameterInterface param;
	public ParameterInterface getParam() {
		return param;
	}
	
	public int proposalCount = 0;
	public int acceptanceCount = 0;
	public boolean lastMoveAccepted = false;	
	public boolean moveProposed = false;

	
	public double proposalWidthControlVariable = 1.0;
	public double spanMultiplier = Utils.SPAN_MULTIPLIER;
	public double minAcceptance = Utils.MIN_ACCEPTANCE;
	public double maxAcceptance = Utils.MAX_ACCEPTANCE;
	public boolean autoTune = true;
	// TODO Add constructor fields for specifying the above two
	// variables.
	
	public String name;
	private long time = 0;
	public long getTime() {
		return time;
	}
	public double acceptanceRate() {
		return (double) acceptanceCount / (double) proposalCount;
	}
	
	public abstract void copyState(Object externalState);
	public abstract double proposal(Object externalState); 
	// Modifies variables and returns logProposalRatio
	
	public abstract double logPriorDensity(Object externalState);
	public abstract void updateLikelihood(Object externalState); 
	public abstract void restoreState(Object externalState);
	public void afterMove(Object externalState) { }

	public boolean isParamChangeAccepted(double logProposalRatio) {
		return getOwner().isParamChangeAccepted(logProposalRatio,this);
	}
	
	public void move(Object externalState) {
		if (Utils.DEBUG) System.out.println("Executing move:\""+name+"\"");
		time -= System.currentTimeMillis();
		proposalCount++;
		moveProposed = true;
		copyState(externalState);
		double logProposalRatio = -logPriorDensity(externalState);
		logProposalRatio += proposal(externalState); 
		logProposalRatio += logPriorDensity(externalState);
		if (logProposalRatio != Double.NEGATIVE_INFINITY) {
			updateLikelihood(externalState);
		}
		if(isParamChangeAccepted(logProposalRatio)) {
			acceptanceCount++;
			lastMoveAccepted = true;
		}
		else {
			lastMoveAccepted = false;
			restoreState(externalState);
		}
		afterMove(externalState);
		time += System.currentTimeMillis();
	}
}
