package hlalign.mcmc;

import java.util.ArrayList;
import hlalign.base.Tree;
import hlalign.base.MCMC;
import hlalign.mcmc.McmcMove;
import hlalign.mcmc.ParameterInterface;
import hlalign.mcmc.PriorDistribution;
import hlalign.mcmc.ProposalDistribution;

public abstract class ContinuousPositiveParameterMove extends McmcMove {
	
	protected Tree tree = null;

	private ProposalDistribution<Double> proposalDistribution;

	protected double oldpar;
	protected double minValue = 0.0;
	protected double maxValue = Double.POSITIVE_INFINITY;
	// If the proposal takes the parameter below/above these values, 
	// the move is rejected.
	
	public ContinuousPositiveParameterMove (MCMC m,
			ParameterInterface p, 
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		owner = m;
		param = p;
		prior = pr;
		name = n;
		proposalDistribution = prop;
		autoTune = true;
		sample = new ArrayList<Double>(0);
	}
	public void setMinValue(double x) {
		minValue = x;
	}
	public void setMaxValue(double x) {
		maxValue = x;
	}
	public void copyState(Object externalState) {
		// externalState is not used for this move - only present for more complicated
		// moves that require restoring more than a single parameter
		oldpar = param.get();
		oldll = owner.getLogLike();
	}

	public double proposal(Object externalState) {
		// externalState is not used for this move - only present for more complicated
		// moves that require restoring more than a single parameter
		proposalDistribution.updateProposal(proposalWidthControlVariable,param.get());
		double par = proposalDistribution.sample(); 
		param.set(par);
		sample.add(par);
		
		if (param.get() < minValue || param.get() > maxValue) {
			return(Double.NEGATIVE_INFINITY);
		}

		/** - log p(new | old) */
		double logProposalDensity = -proposalDistribution.logDensity(param.get());
		
		proposalDistribution.updateProposal(proposalWidthControlVariable,param.get());
		
		/** + log p(old | new) */
		logProposalDensity += proposalDistribution.logDensity(oldpar);
		
		return logProposalDensity;
	}
	
	public double logPriorDensity(Object externalState) {
		if (param.get() < minValue || param.get() > maxValue) {
			return(Double.NEGATIVE_INFINITY);
		}
		else {
			// return 0;
			double logDensity = prior.logDensityUnnormalised(param.get());
			// Since we're only using this in ratios, there's no
			// need to compute the normalising constant, which is good, 
			// because some priors may be improper.
			// NB be careful with this though -- an improper prior should
			// only be used if the posterior can be shown to be proper.
			return logDensity;
		}
	}
	public abstract void updateLikelihood(Object externalState);
	public void restoreState(Object externalState) {
		param.set(oldpar);
		//replace the most recent sample (proposed value) with oldpar
		sample.set(sample.size()-1, oldpar);
		owner.setLogLike(oldll);
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		proposalDistribution.updateProposal(proposalWidthControlVariable, currentParam);
	}
}
