package hlalign.mcmc;

import hlalign.base.MCMC;


public class ZetaMove extends IndelMove{
	
	public ZetaMove (MCMC m, 
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,pr,prop,n);	
		param = new ZetaInterface();
	}
}

