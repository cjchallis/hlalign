package hlalign.mcmc;

import hlalign.base.Utils;
import hlalign.base.MCMC;
import hlalign.mcmc.ContinuousPositiveParameterMove;
import hlalign.mcmc.ParameterInterface;
import hlalign.mcmc.PriorDistribution;
import hlalign.mcmc.ProposalDistribution;

public class EdgeMove extends ContinuousPositiveParameterMove {

	private int index;

	class EdgeInterface implements ParameterInterface {
		int ind;
		EdgeInterface (int i) {
			ind = i;
		}
		public double get() {
			return owner.tree.vertex[ind].edgeLength;
		}
		public void set(double x) {
			owner.tree.vertex[ind].edgeLength = x;
		}
	}
	
	public EdgeMove (MCMC m, 
			int edgeIndex,
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(m,null,pr,prop,n);
		index = edgeIndex;
		param = new EdgeInterface(index);
		//minValue = 0.01;
		minValue = Utils.MIN_EDGE_LENGTH;
	}
	public EdgeMove (MCMC m, 
			int edgeIndex,
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, 
			double propVar, String n) {
		super(m,null,pr,prop,n);
		index = edgeIndex;
		param = new EdgeInterface(index);
		//minValue = 0.01;
		minValue = Utils.MIN_EDGE_LENGTH;
		proposalWidthControlVariable = propVar;
	}
	void setEdgeIndex(int i) {
		index = i;
		param = new EdgeInterface(index);
	}
	int getEdgeIndex() {
		return index;
	}
	@Override
	public void copyState(Object externalState) {
		super.copyState(externalState);
	}
	/*@Override
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterEdgeLenChange(tree,tree.vertex[index],lastMoveAccepted);
	}*/
	@Override
	public double proposal(Object externalState) {
		double logProposalRatio = super.proposal(externalState);
		if (param.get() >= minValue && param.get() < maxValue) {
			owner.tree.vertex[index].edgeChangeUpdate();
		}
		return logProposalRatio;
	}
	public void updateLikelihood(Object externalState) {
		if (param.get() >= minValue && param.get() < maxValue) {
			owner.tree.calcML();
		}
	}
	@Override
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		owner.tree.vertex[index].edgeChangeUpdate();
	}
		
}
