package hlalign.mcmc;

import hlalign.base.MCMC;
import hlalign.mcmc.PriorDistribution;
import hlalign.mcmc.ProposalDistribution;

public class EtaMove extends IndelMove{
	
	public EtaMove (MCMC m, 
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,pr,prop,n);	
		param = new EtaInterface();
	}
}
