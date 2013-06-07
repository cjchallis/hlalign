package hlalign.mcmc;

import hlalign.base.MCMC;
import hlalign.mcmc.ContinuousPositiveParameterMove;
import hlalign.mcmc.ParameterInterface;
import hlalign.mcmc.PriorDistribution;
import hlalign.mcmc.ProposalDistribution;

public abstract class IndelMove extends ContinuousPositiveParameterMove {
	
	class EtaInterface implements ParameterInterface {
		// eta = lambda / mu
		public double get() {
			return owner.subModel.l / owner.subModel.mu;
		}
		public void set(double eta) {
			double zeta = owner.subModel.l * owner.subModel.mu;
			// change l and mu, keeping zeta fixed
			owner.subModel.l = Math.sqrt(eta * zeta);
			owner.subModel.mu = Math.sqrt(zeta / eta);
			owner.subModel.updateDecomposition();
		}
	}
	class ZetaInterface implements ParameterInterface {
		// zeta = lambda * mu
		public double get() {
			return(owner.subModel.l * owner.subModel.mu);
		}
		public void set(double zeta) {
			double eta = owner.subModel.l / owner.subModel.mu;
			// change l and mu, keeping eta fixed
			owner.subModel.l = Math.sqrt(eta * zeta);
			owner.subModel.mu = Math.sqrt(zeta / eta);
		}
	}
	
	public IndelMove (MCMC m, 
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,null,pr,prop,n);		
	}

	@Override
	public void copyState(Object externalState) {
		super.copyState(externalState);
	}

	
	public void updateLikelihood(Object externalState) {
		owner.tree.calcML();
	}
	
	@Override
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
	}

}

