package jackSolver;





public class GaussianOrbital {
	/* 
	 * PSI=SIGMA(Di*psii(alphai))
	 * 
	 * psi=SIGMAoveri(X^l*Y^m*Z^n*EXP(-alphai*r^2))
	 * 
	 * 
	 * 
	 */
	
	//sto-ng data for hydrogen
	public static double[][] zetas = {{0.0},{1.24},{2.0925},{1,2}};//corresponds to the elements H=1,He=2, 0=non existant element unless you play mass effect
	public static int STOnG = 3;
	public static double[][] coef = {{1,0,0},{0.678914,0.430129,0},{0.444635,0.535328,0.154329}};
	public static double[][] expon = {{0.270950,0,0},{0.151623,0.851819,0},{0.109818,0.405771,2.22766}};
	
	
	public double[] lmn; //specifies the 
	public double alpha[]; //constains the alphas
	public double x,y,z; //location
	public double D[]; //contains the contraction coeficients
	public double Zeta;
	public GaussianOrbital(double l,double m, double n, double[] alphain, double[] Din,double xin, double yin, double zin){
		lmn = new double[3];
		lmn[0]=l;
		lmn[1]=m;
		lmn[2]=n;
		
		alpha=alphain;
		D=Din;
				
		x=xin;
		y=yin;
		z=zin;
		
	}
	
	public GaussianOrbital(double l,double m, double n, double[] alphain, double[] Din,double zetain){
		lmn = new double[3];
		lmn[0]=l;
		lmn[1]=m;
		lmn[2]=n;
		
		alpha=alphain;
		D=Din;
		Zeta=zetain;
	}
	
	public String toString(){
		return "GTO x: "+x+" y: "+y+" z: "+z;
	}
	
	public static GaussianOrbital[] generateOrbitalArray(Atom[] atoms,int STOnG) {
		GaussianOrbital[] orbitals = new GaussianOrbital[atoms.length];
		for(int i=0;i<orbitals.length;i++){
			//if(DEBUG)System.out.println("......Adding Orbital: "+i);
			orbitals[i]=generateorbital(atoms[i],STOnG);
			
		}
		return orbitals;
	}

	public static GaussianOrbital generateorbital(Atom atom,int STOnG) {
		
		// will eventually properly fill out the full contract orbitals with alphas and zetas
		
		/*
		 * STO-nG basis set data
		 * 
		 * PSI = sigma(Di*PSIi(alphai))
	
			COEF aka d for hydrogen
			
			d1s hydrogen sto 1g
			1.0D0,
			0.0D0,
			0.0
			
			d1s hydrogen sto 2g
			0.678914D0,
			0.430129D0,
			0.0D0,
			
			 
			 d1s sto 3g
			 0.444635D0,
			 0.535328D0,
			 0.154329D0,
			 
			 
			 EXPONS aka ALPHAS
			 
			 alpha 1s STO-1G
			 0.270950D0,
			 0.0D0
			 0.0D0,
			 
			 alpha 1s STO-2G
			 0.151623D0,
			 0.851819D0, 
			 0.0D0, 
			 
			 alpha 1s sto 3g
			 0.109818D0, 
			 0.405771D0, 
			 2.22760D0
		 */
		
		
		
		//double zeta1 = 2.0925;
		//double zeta2 = 1.24;  i thinkthis is for
		
		double[] zetas = {0.0,1.24,2.0925};//corresponds to the elements H=1,He=2, 0=non existant element unless you play mass effect
		//if(int STOnG = 1;
		double[][] coef = {{1,0,0},{0.678914,0.430129,0},{0.444635,0.535328,0.154329}};
		double[][] expon = {{0.270950,0,0},{0.151623,0.851819,0},{0.109818,0.405771,2.22766}};
		//double alpha =0.270950;//this is from appendix b sto-1g;
		
		double[] alphas = new double[STOnG];
		double[] Ds = new double[STOnG];
		
		for(int i=0;i<alphas.length;i++){
			//if(DEBUG)System.out.println("generating orbitals atom.element: "+atom.element+"  expon[STOnG][i]:  "+expon[STOnG-1][i]);
			alphas[i]=expon[STOnG-1][i]*zetas[atom.element]*zetas[atom.element];
			//if(DEBUG)System.out.println("generating orbitals alphas[i]:  "+alphas[i]);
			Ds[i]=coef[STOnG-1][i]*Math.pow(2.0*alphas[i]/Math.PI, 0.75);
		}
		
		
		return new GaussianOrbital(0, 0, 0, alphas,Ds, atom.x, atom.y, atom.z);
	}

	public static double calcR(GaussianOrbital a, GaussianOrbital b){
		System.out.println("CalcR a:"+a);
		double x = a.x-b.x;
		double y = a.y-b.y;
		double z = a.z-b.z;	
		double RaRbdifsqrd = x*x+y*y+z*z;
		
		return Math.sqrt(RaRbdifsqrd);
	}
	
	public static double calcRabsqrd(GaussianOrbital a, GaussianOrbital b){
		double x = a.x-b.x;
		double y = a.y-b.y;
		double z = a.z-b.z;	
		double RaRbdifsqrd = x*x+y*y+z*z;
		
		return RaRbdifsqrd;
	}
	
	public static double calcRpcsqrd(GaussianOrbital a, int Aalphai,GaussianOrbital b,  int Balphai,  double Xc, double Yc, double Zc){
		/*
		 * 
		 *   A
		 *   |
		 *   |
		 *   P-------------C
		 *   |    Rpc
		 *   |
		 *   |
		 *   |
		 *   B
		 *   
		 *   Rp=(alpha*Ra+beta*Rb)/(alpha+beta)
		 *   
		 *   
		 */
		double Xp =(a.alpha[Aalphai]*a.x+b.alpha[Balphai]*b.x)/(a.alpha[Aalphai]+b.alpha[Balphai]);
		double Yp =(a.alpha[Aalphai]*a.y+b.alpha[Balphai]*b.y)/(a.alpha[Aalphai]+b.alpha[Balphai]);
		double Zp =(a.alpha[Aalphai]*a.z+b.alpha[Balphai]*b.z)/(a.alpha[Aalphai]+b.alpha[Balphai]);
		
		double x = Xp-Xc;
		double y = Yp-Yc;
		double z = Zp-Zc;	
		double RaRbdifsqrd = x*x+y*y+z*z;
		
		return RaRbdifsqrd;
		
		
	}

	public static double calcRpqsqrd(GaussianOrbital[] orbitals, int a,
			int b, int c, int d, int alphai, int betai, int gammai, int deltai) {
		
		/*
		 *                 C
		 *   A             |
		 *   |             |
		 *   |             |
		 *   P-------------Q
		 *   |    Rpc      |
		 *   |             |
		 *   |             D
		 *   |
		 *   B
		 *   
		 *   Rp=(alpha*Ra+beta*Rb)/(alpha+beta)
		 *   
		 *   
		 */
		// TODO Auto-generated method stub
		Vec3D Rp= calcRp(orbitals[a],alphai,orbitals[b],betai);
		Vec3D Rq= calcRp(orbitals[c],gammai,orbitals[d],deltai);
		
		//System.out.println(Rp.toString()+Rq.toString());
		
		double Xrpc=Rp.x-Rq.x;
		double Yrpc=Rp.y-Rq.y;
		double Zrpc=Rp.z-Rq.z;
		
		double Rpqsqrd = Xrpc*Xrpc+Yrpc*Yrpc+Zrpc*Zrpc;
		
		return Rpqsqrd;
	}
	
	public static Vec3D calcRp(GaussianOrbital a, int Aalphai,GaussianOrbital b,  int Balphai){
		/*
		 * new center between 2 gaussians
		 * 
		 *   A
		 *   |
		 *   |
		 *   p
		 *   |
		 *   |
		 *   |
		 *   B
		 *   
		 * 
		 */
		double Xp =(a.alpha[Aalphai]*a.x+b.alpha[Balphai]*b.x)/(a.alpha[Aalphai]+b.alpha[Balphai]);
		double Yp =(a.alpha[Aalphai]*a.y+b.alpha[Balphai]*b.y)/(a.alpha[Aalphai]+b.alpha[Balphai]);
		double Zp =(a.alpha[Aalphai]*a.z+b.alpha[Balphai]*b.z)/(a.alpha[Aalphai]+b.alpha[Balphai]);
		Vec3D vec = new Vec3D(Xp,Yp,Zp);
		
		return vec;
		
	}
	
	
	
	

}
