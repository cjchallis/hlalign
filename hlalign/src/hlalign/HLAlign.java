package hlalign;
import hlalign.io.DataReader;
import hlalign.base.MCMC;

public class HLAlign {
	
	public static void main(String[] args){
		System.out.println("Hierarchical Links Alignment v0.0");
		
		DataReader datareader = new DataReader();
		datareader.readAllCoords(args);
		
		MCMC mcmc = new MCMC(5, 5, datareader.coords, args);
		mcmc.run();
		System.out.println("Done");
	}
}
