package jackSolver;

public class ElementBasisSet{
	public static double[] zetas = {0.0,1.24,2.0925};//corresponds to the elements H=1,He=2, 0=non existant element unless you play mass effect
	public static int STOnG = 3;
	public static double[][] coef = {{1,0,0},{0.678914,0.430129,0},{0.444635,0.535328,0.154329}};
	public static double[][] expon = {{0.270950,0,0},{0.151623,0.851819,0},{0.109818,0.405771,2.22766}};
	
	GaussianOrbital[] orbitals;
	int Za;
	public ElementBasisSet(){		
	}
	
	public void setSTO3G(int Zain){
		Za=Zain;
		if(Za==0){
			orbitals=null;
		}else if(Za==1){//hydrogen
			orbitals = new GaussianOrbital[1];
			orbitals[0] = new GaussianOrbital(0, 0, 0, GaussianOrbital.expon[3], GaussianOrbital.coef[3], GaussianOrbital.zetas[Za]);
		}else if(Za==2){//helium
			orbitals = new GaussianOrbital[1];
			orbitals[0] = new GaussianOrbital(0, 0, 0, GaussianOrbital.expon[3], GaussianOrbital.coef[3], GaussianOrbital.zetas[Za]);
		}else if(Za==6){//lithium
			orbitals = new GaussianOrbital[5];
			orbitals[0] = new GaussianOrbital(0, 0, 0, GaussianOrbital.expon[3], GaussianOrbital.coef[3], GaussianOrbital.zetas[Za]);
			
		}
		
		
	}
	
}
