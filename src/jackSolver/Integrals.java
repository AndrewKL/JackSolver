package jackSolver;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import com.sun.corba.se.spi.orbutil.threadpool.ThreadPool;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class Integrals {

	public static DenseDoubleMatrix2D GenerateSOverlapIntegrals(
			GaussianOrbital[] orbitals) {
		
		/*
		 * generates a KxK matrix with the overlap integrals
		 * 
		 * k=number of orbitals
		 */
		DenseDoubleMatrix2D Soverlapintegrals = new DenseDoubleMatrix2D(orbitals.length, orbitals.length);
		for(int i=0;i<Soverlapintegrals.rows();i++){
			for(int j=0;j<Soverlapintegrals.columns();j++){
				//if(DEBUG)System.out.println("calc S mu: "+i+" nu: "+j);
				
				Soverlapintegrals.set(i, j, Integrals.CalculateSintegral(orbitals,i,j));
			}
				
		}
		return Soverlapintegrals;
	}

	public static double CalculateSintegral(GaussianOrbital[] orbitals, int a,
			int b) {
		/*
		 *  S u,v = Integral(psi*u(1)*psiu(1)
		 *  
		 *  (A|B) = [pi/(alpha+beta)]^3/2 *EXP(-alpha*beta/(alpha+beta)* |Ra-Rb|^2) 
		 *  
		 *  alpha = first gaussian exponent
		 *  beta = second gaussian exponent
		 *  
		 *  scans through two contracted gaussian orbitals and calculates the integrals for the sets and sums them up 
		 *  
		 */
		double S = 0;
		for(int i=0;i<orbitals[a].alpha.length;i++){//scans through contracted gassian orbital a's simple orbitals
			for(int j=0;j<orbitals[b].alpha.length;j++){
				S += Integrals.calculateSIntegral(orbitals[a].alpha[i],
						orbitals[b].alpha[j],
						GaussianOrbital.calcRabsqrd(orbitals[a], orbitals[b]))
						* orbitals[a].D[i]*orbitals[b].D[j];
			}
			
		}
		
		
		return S;
	}

	public static double calculateSIntegral(double alpha, double beta,
			double RaRbdifsqrd) {
		/*
		 *  S u,v = Integral(psi*u(1)*psiu(1)
		 *  
		 *  (A|B) = [pi/(alpha+beta)]^3/2 *EXP(-alpha*beta/(alpha+beta)* |Ra-Rb|^2) 
		 *  
		 *  alpha = first gaussian exponent
		 *  beta = second gaussian exponent
		 *  
		 *  scans through two contracted gaussian orbitals and calculates the integrals for the sets and sums them up 
		 *  
		 */
		
		double powterm = Math.PI/(alpha+beta);
		
		double normalizationconst=Math.pow(powterm, 1.5);
		double insideexp = -(alpha*beta*RaRbdifsqrd)/(alpha+beta);
		double integral =Math.exp(insideexp);
		
		//if(DEBUG)System.out.println("calc s integral  alpha: "+alpha+"  beta: "+beta+"  RaRbdifsqrd: "+RaRbdifsqrd +"   "+(normalizationconst*integral));
		//if(DEBUG)System.out.println("                 powterm: "+powterm+" normalization const: "+ normalizationconst);
		//if(DEBUG)System.out.println("                 integralcomp: "+integral+"  insideexp: "+insideexp);
		return normalizationconst*integral;
		
	}

	public static double calcTkineticenergyintegral(GaussianOrbital[] orbitals,
			int a, int b) {
		
		//this scans through the contracted gaussian orbitals
		
		double T = 0; //kinetic energy sum
		for(int i=0;i<orbitals[a].alpha.length;i++){//scans through contracted gassian orbital a's simple orbitals
			for(int j=0;j<orbitals[b].alpha.length;j++){
				T += Integrals.calculateTkineticenergyIntegral(orbitals[a].alpha[i],
						orbitals[b].alpha[j],
						GaussianOrbital.calcRabsqrd(orbitals[a], orbitals[b]))
						* orbitals[a].D[i]*orbitals[b].D[j];
			}
			
		}
		return T;
	}

	public static double calculateTkineticenergyIntegral(double alpha, double beta,
			double RaRbsqrd) {
		double T = (alpha*beta)/(alpha+beta)*(3-2*alpha*beta*RaRbsqrd/(alpha+beta))*Math.pow(Math.PI/(alpha+beta), 1.5)*Math.exp(-alpha*beta*RaRbsqrd/(alpha+beta));
		
		return T;
	}

	public static double calcVnucintegral(GaussianOrbital[] orbitals, Atom[] atoms, int a, int b) {
		
		//double Rcpsqrd =0;
		 //Vnuc energy sum
		double Vnucsum =0;
		for(int k=0;k<atoms.length;k++){
			double Vnucsumatom = 0;
			for(int j=0;j<orbitals[b].alpha.length;j++){
				for(int i=0;i<orbitals[a].alpha.length;i++){//scans through contracted gassian orbital a's simple orbitals//scanning through atoms
					double Rcpsqrd = GaussianOrbital.calcRpcsqrd(orbitals[a], i, orbitals[b], j, atoms[k].x, atoms[k].y, atoms[k].z);
					double Rabsqrd = GaussianOrbital.calcRabsqrd(orbitals[a], orbitals[b]);
					double Vnuctemp = Integrals.calculateVnucIntegral(orbitals[a].alpha[i],
							orbitals[b].alpha[j],
							Rabsqrd,Rcpsqrd,
							atoms[k].element)
							* orbitals[a].D[i]*orbitals[b].D[j];
					//if(DEBUG)System.out.println("calcVnuc atom :"+k+"  a: "+(a+1)+"  b: "+(b+1)+" Vnuc: "+Vnuctemp);
					Vnucsumatom += Vnuctemp;
					
				}
				
			}
			//if(DEBUG)System.out.println(" Calc Vnuc atom k: "+k+" vnucsum atom: "+Vnucsumatom);
			Vnucsum += Vnucsumatom;
			
		}
		
		
		return Vnucsum;
	}

	public static double calculateVnucIntegral(double alpha, double beta,
			double Rabsqrd, double Rcpsqrd, int Zelement) {
				
		double V = 2*Math.PI/(alpha+beta)*Integrals.F0((alpha+beta)*Rcpsqrd)*Math.exp(-alpha*beta*Rabsqrd/(alpha+beta));
		double Vnuc=-V*Zelement;
		return Vnuc;
	}

	public static double calculateFourCenterIntegral(double alpha, double beta,
			double gamma, double delta, double Rabsqrd,
			double Rcdsqrd, double Rpqsqrd) {
		//this does the actual calculation of the integral
		
		return 2*Math.pow(Math.PI, 2.5)
				/((alpha+beta)*(gamma+delta)*Math.sqrt(alpha+beta+gamma+delta))
				*Math.exp(-(alpha*beta/(alpha+beta)*Rabsqrd)-(gamma*delta/(gamma+delta)*Rcdsqrd))
				*Integrals.F0((alpha+beta)*(gamma+delta)/(alpha+beta+gamma+delta)*Rpqsqrd);
		
	}

	public static double calculateFourCenterIntegrals(
			GaussianOrbital[] orbitals, int a, int b, int c, int d) {
		// TODO Auto-generated method stub
		double fourCtotal = 0; //kinetic energy sum
		for(int i=0;i<orbitals[a].alpha.length;i++){//scans through contracted gassian orbital a's simple gaussian orbitals
			for(int j=0;j<orbitals[b].alpha.length;j++){
				for(int k=0;k<orbitals[c].alpha.length;k++){
					for(int l=0;l<orbitals[d].alpha.length;l++){
						//System.out.println("4ci a: "+(a+1)+" b: "+(+b+1)+" c: "+(c+1)+" d: "+(1+d));
						//System.out.println(" i: "+i+" j "+j+" k "+k );
						//System.out.println(" alpha: "+orbitals[a].alpha[i]+" beta: "+orbitals[b].alpha[j]+" gamma: "+orbitals[c].alpha[k]+" delta: "+orbitals[d].alpha[l]);
						
						double Rpqsqrd = GaussianOrbital.calcRpqsqrd(orbitals,a,b,c,d,i,j,k,l);
						
						//System.out.println(" Rpqsrd: "+Rpqsqrd);
						//double Rabsqrd = GaussianOrbital.calcRabsqrd(orbitals[a], orbitals[b]);
						//System.out.println(" Rabsqrd: "+Rabsqrd);
						//double Rcdsqrd = GaussianOrbital.calcRabsqrd(orbitals[c], orbitals[d]);
						//System.out.println(" Rcdsqrd: "+Rcdsqrd);
						
						double fourCI = calculateFourCenterIntegral(orbitals[a].alpha[i],
								orbitals[b].alpha[j],
								orbitals[c].alpha[k],
								orbitals[d].alpha[l],
								GaussianOrbital.calcRabsqrd(orbitals[a], orbitals[b]),
								GaussianOrbital.calcRabsqrd(orbitals[c], orbitals[d]),
								Rpqsqrd);								
						
						//System.out.println(" 4ci: "+fourCI);
						double D = orbitals[a].D[i]*orbitals[b].D[j]*orbitals[c].D[k]*orbitals[d].D[l];
						//System.out.println(" D: "+D);
						double fourCsubtotal = fourCI*D;
						//System.out.println(" fourCsubtotal: "+fourCsubtotal);
						
						fourCtotal += fourCsubtotal;
					}
				}
				
			}
			
		}
		return fourCtotal;
	}

	public static double[][][][] calculate4centerintegrals(final GaussianOrbital[] orbitals) {
		
		
		final double[][][][] fourcenterintegrals = new double[orbitals.length][orbitals.length][orbitals.length][orbitals.length];
		
		ExecutorService pool =  Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		System.out.println("Threads available: "+Runtime.getRuntime().availableProcessors());
		//System.out.println("starting loop");
		for(int z=0;z<fourcenterintegrals.length;z++){
			final int i = z;
			pool.execute(new Runnable(){public void run(){
				//System.out.println("new thread started i: "+i);
				for(int j=0;j<fourcenterintegrals[i].length;j++){
					for(int k=0;k<fourcenterintegrals[i][j].length;k++){
						for(int l=0;l<fourcenterintegrals[i][j][k].length;l++){
							//System.out.println(" calcing 2ci");
							fourcenterintegrals[i][j][k][l]=calculateFourCenterIntegrals(orbitals,i,j,k,l);
						}
					}
				}
				//System.out.println("thread done");
				
		        
		    }});

			
			
		}
		pool.shutdown();
		//System.out.println("loop done");
		try {pool.awaitTermination(5, TimeUnit.MINUTES );// blocks till everything is done or 5 minutes have passed.
		} catch (InterruptedException e) {	e.printStackTrace();
		} 
		System.out.println("integrals calculated: "+(orbitals.length*orbitals.length*orbitals.length*orbitals.length));
		return fourcenterintegrals;
	}
	
	
	public static double F0(double d) {
		
		//calculates the f function for F0 only s type
		if(d<0.000001)return 1-(d/3);
		
		return Math.pow(Math.PI/d, 0.5)*ErrorFunction.erf(Math.pow(d, 0.5))/2;
	}

}
