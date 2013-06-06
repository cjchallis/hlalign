package hlalign.mcmc;

public interface PriorDistribution<T> {
	public abstract double logDensity(Double x);
	public abstract double logDensityUnnormalised(Double x);
}
