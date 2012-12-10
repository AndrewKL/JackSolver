package jackSolver;


import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;

public class JackRHFSolver {
	
	/*
	 * JackRHFSolver
	 * 
	 * by Andrew Long 
	 * Advisor: Jason A.C. Clyburne 
	 * with additional assistance from Cory Pye.
	 * 
	 * based on the Fortran IV program by attila Szabo and Neil S. Ostlund
	 * in the book Modern Quantum Chemisty: introduction to advanced electronic 
	 * structure theory/ Attila Szabo, Neil S. Ostlund
	 *  
	 *  If you're just learning this stuff it is well worth your time to get the 
	 *  book. It's excellent.  The steps in the comments in the main part of this 
	 *  program match up with the steps on page 146 of the book.
	 *  
	 *  
	 */
	
	
	
	
	static boolean DEBUG =true;
	
	
	public static void main(String[] args){
		System.out.println("Check 1 2 3");
		System.out.println("calculation RHF for H-He+ ");
		Atom atom1 = new Atom(0,0,1.463200,1);
		Atom atom2 = new Atom(0,0,0,2);
		//Atom atom3 = new Atom(0,2,2,1);
		//Atom atom4 = new Atom(0,2,3,1);
		Atom[] atoms = {atom2,atom1};
		
		int numberofelectrons =2;
		System.out.println("Calculating H-He energy");
		
		HFSolver(atoms,numberofelectrons);
		
	}

	public static double HFSolver(Atom[] atoms, int numberofelectrons) {
		
		int STOnG = 3;
		
		
		
		
		//Writing Information about the input
		if(DEBUG)System.out.println("JQuantum RHFSolver by Andrew Long");
		
		if(numberofelectrons%2!=0){
			if(DEBUG)System.out.println("uneven number of electrons");
			return 0;
		}
		if(DEBUG)System.out.println("===List of Atoms=============");
		if(DEBUG){
			for(int i=0;i<atoms.length;i++){
				System.out.println(" Z: "+atoms[i].element+" x: "+atoms[i].x+" y: "+atoms[i].y+" z: "+atoms[i].z);
			}
		}
		if(DEBUG)System.out.println("...Generating Orbital Array");
		
		GaussianOrbital[] orbitals = generateOrbitalArray(atoms,STOnG);
		if(DEBUG)for(int i=0;i<orbitals.length;i++)System.out.println(orbitals[i]);
		
		//Here's where the real work begins.
		
		
		
		DoubleMatrix2D sOverlapMatrix, 
		T,
		Vnuc,
		fFockMatrix, 
		fPrimeMatrix,
		fPrimeDiagonalizedMatrix,
		cCoefficientsMatrix, 
		cPrimeMatrix,
		pDensityMatrix, 
		oldpDensityMatrix,
		xMatrix, 
		xAdjointMatrix,
		gMatrix,
		hCore;
		
		
		//step #2 calculate integrals
		
		sOverlapMatrix = GenerateSOverlapIntegrals(orbitals);
		if(DEBUG)System.out.println("===S Matrix==========");
		if(DEBUG)System.out.println(sOverlapMatrix);
		
		T = new DenseDoubleMatrix2D(orbitals.length, orbitals.length);
		Vnuc = new DenseDoubleMatrix2D(orbitals.length, orbitals.length);
		
		hCore = CalculateHCore(orbitals, atoms, T, Vnuc);
		if(DEBUG)System.out.println("===Hcore=============");
		if(DEBUG)System.out.println(hCore);
		if(DEBUG)System.out.println("===T=================");
		if(DEBUG)System.out.println(T);
		if(DEBUG)System.out.println("===Vnuc==============");
		if(DEBUG)System.out.println(Vnuc);
		
		//step #3 diagonalize S obtain  transform matrix X
		
		
		
		
		
		xMatrix = cannonicalOrthogonalization(sOverlapMatrix);
		if(DEBUG)System.out.println("===X matrix===========");
		if(DEBUG)System.out.println(xMatrix);
		
		
		//step #3 four integral debugging step  in the book they calculate these and put the in matrix
		
		//fourcenterintegrals = calculate4centerintegrals(orbitals);//this is now done in the test cases
		
		
		
		//step #4 estimate density matrix
		
		if(DEBUG)System.out.println("===P density matrix====");
		pDensityMatrix = EstimateDenstiyMatrix(orbitals,numberofelectrons);
		if(DEBUG)System.out.println(pDensityMatrix);
		
		
		//begin SCF loop
		Algebra alg = new Algebra();
		gMatrix = new DenseDoubleMatrix2D(orbitals.length,orbitals.length);
		fFockMatrix = new DenseDoubleMatrix2D(orbitals.length,orbitals.length);
		
		double currentElectronicEnergy=0;
				
		
		
		int MaxSCFiterations = 30;
		double ConvergenceCriteria = 0.0001;
		
		double currentConvergence = 1.0;
		int currentIteration = 0;
		while(currentIteration<MaxSCFiterations && currentConvergence>ConvergenceCriteria){
			currentIteration++;
			
			System.out.println("========================================");
			System.out.println("===SCF LOOP=============================");
			System.out.println(" current iteration: "+currentIteration);
			
			//step #5 calculate Matrix G
			
			gMatrix.assign(0);
			calcGMatrix(gMatrix, pDensityMatrix,orbitals);
			System.out.println("===G matrix===========");
			System.out.println(gMatrix);
			
			
			
			//step #6 calculate fock matrix F=Hcore+G
			calcFockMatrix(fFockMatrix,hCore,gMatrix);
			System.out.println("===fFockMatrix========");
			System.out.println(fFockMatrix);
			
			//step #6 and 1/3  calculate the current energy
			System.out.println("======================");
			currentElectronicEnergy = calculateElectronicEnergy(pDensityMatrix, hCore, fFockMatrix);
			System.out.println(" electronic energy: "+currentElectronicEnergy);
			
			//step #7 calc transformed fock matrix F'=X(+)FX
			xAdjointMatrix = adjoint(xMatrix);
			fPrimeMatrix=alg.mult(xAdjointMatrix, alg.mult(fFockMatrix, xMatrix));
			System.out.println("===fPrimeMatrix========");
			System.out.println(fPrimeMatrix);
			
			
			
			//step #8 diagonalize F' to get C' and epsilon
			// F'C'=C'epsilon
			
			System.out.println("===fPrimeDiagonalizedMatrix========");
			//EigenvalueDecomposition evd = new EigenvalueDecomposition(fPrimeMatrix);//aka diagonalize
			fPrimeDiagonalizedMatrix=diagonalize2x2MatrixgetEigenValues(fPrimeMatrix);
			//fPrimeDiagonalizedMatrix=evd.getD();
			System.out.println(fPrimeDiagonalizedMatrix);
			
			System.out.println("===cPrimeMatrix========");
			cPrimeMatrix=diagonalize2x2MatrixgetEigenVectors(fPrimeMatrix);
			
			
			//cPrimeMatrix=evd.getV();
			System.out.println(cPrimeMatrix);
			
			//step #9 calc C=XC'
			System.out.println("===cCoefficientsMatrix========");
			cCoefficientsMatrix = alg.mult(xMatrix, cPrimeMatrix);
			System.out.println(cCoefficientsMatrix);
			
			//step #10 calc new density matrix
			System.out.println("===new pDensityMatrix========");
			oldpDensityMatrix=pDensityMatrix.copy();
			calcNewDensityMatrix(pDensityMatrix,cCoefficientsMatrix,numberofelectrons);
			System.out.println(pDensityMatrix);
			
			//step #11 test convergence and calc expectation values and other quantities of interest
			
			currentConvergence=calcConvergence(pDensityMatrix,oldpDensityMatrix);
			System.out.println("=============================");
			System.out.println("currentConvergence : "+currentConvergence);
			
			
			
		}
		
		//clean up
		System.out.println("========================================");
		System.out.println("========================================");
		
		if(currentConvergence>ConvergenceCriteria){
			System.out.println("Convergence Criteria not met after max loops. Ending scf loop");
		}else if(currentConvergence<=ConvergenceCriteria){
			System.out.println("Convergence Criteria met. Ending scf loop");
			
		}
		System.out.println("Electronic Energy: "+currentElectronicEnergy);
		double nuclearRepulsionEnergy = calcNuclearRepulsionEnergy(atoms);
		System.out.println("Nuclear Repulsion Energy: "+nuclearRepulsionEnergy);
		double totalEnergy = currentElectronicEnergy+nuclearRepulsionEnergy;
		System.out.println("Total Energy: "+totalEnergy);
		
		
		System.out.println();
		System.out.println("===Mulliken Populaton========");
		System.out.println(alg.mult(pDensityMatrix, sOverlapMatrix));
		
		
		
		
		return totalEnergy;
	}
	

	
	public static double calcNuclearRepulsionEnergy(Atom[] atoms) {
		double sum = 0;
		for(int i=0;i<atoms.length;i++){
			for(int j=0;j<i;j++){
				//System.out.println(" NRE loop i"+i+" j"+j);
				sum += atoms[i].element*atoms[j].element/Atom.distance(atoms[i],atoms[j]);
				//System.out.print(" Za: "+atoms[i].element+" Zb: "+atoms[j].element+" dist: "+Atom.distance(atoms[i],atoms[j]));
				//System.out.println(" += "+(atoms[i].element*atoms[j].element/Atom.distance(atoms[i],atoms[j])));
			}
		}
		return sum;
	}

	public static double calcConvergence(DoubleMatrix2D pDensityMatrix,
			DoubleMatrix2D oldpDensityMatrix) {
		double delta = 0;
		for(int i=0;i<pDensityMatrix.rows();i++){
			for(int j=0;j<pDensityMatrix.columns();j++){
				
				delta += (pDensityMatrix.get(i,j)-oldpDensityMatrix.get(i, j))*
						(pDensityMatrix.get(i,j)-oldpDensityMatrix.get(i, j));
			}
		}
		delta=Math.sqrt(delta/(pDensityMatrix.rows()*pDensityMatrix.columns()));
		return delta;
	}

	public static void calcNewDensityMatrix(DoubleMatrix2D pDensityMatrix,
			DoubleMatrix2D cCoefficientsMatrix, int numberofelectrons) {
		
		pDensityMatrix.assign(0); //zero the density Matrix
		for(int i=0;i<pDensityMatrix.rows();i++){
			for(int j=0;j<pDensityMatrix.columns();j++){
				for(int n=0;n<(numberofelectrons/2);n++){
					//System.out.println(" i"+i+"j"+j+" n"+n+ " currentP "+pDensityMatrix.get(i,j)+" Cin"+cCoefficientsMatrix.get(i, n)+" Cjn"+cCoefficientsMatrix.get(j, n));
					pDensityMatrix.set(i, j, pDensityMatrix.get(i,j)+2*cCoefficientsMatrix.get(i, n)*cCoefficientsMatrix.get(j, n));
				}
			}
		}


		
	}

	private static double calculateElectronicEnergy(
			DoubleMatrix2D pDensityMatrix, DoubleMatrix2D hCore,
			DoubleMatrix2D fFockMatrix) {
		// TODO Auto-generated method stub
		double energy =0;
		for(int i = 0;i<fFockMatrix.rows();i++){
			for(int j=0;j<fFockMatrix.columns();j++){
				energy += 0.5*pDensityMatrix.get(i, j)*(hCore.get(i,j)+fFockMatrix.get(i,j));
			}
		}
		return energy;
	}

	private static void calcFockMatrix(DoubleMatrix2D fFockMatrix,
			DoubleMatrix2D hCore, DoubleMatrix2D gMatrix) {
		
		for(int i = 0;i<fFockMatrix.rows();i++){
			for(int j=0;j<fFockMatrix.columns();j++){
				fFockMatrix.set(i, j, hCore.get(i, j)+gMatrix.get(i, j));
			}
		}
		
	}

	public static void calcGMatrix(DoubleMatrix2D gMatrix,
			DoubleMatrix2D pDensityMatrix, GaussianOrbital[] orbitals) {
		/*
		 * Gmu,nu=SIGM[Aalpa to N/2]SIGMA[lamda to K]SIGMA[sigma to K] C[lamda,.....blah balah ablah
		 * 
		 *       =SIMGA[lamda to K][sigma to K](P[lamda,sigma]*[(mu nu|sigma lamda)-1/2(mu lamda|sigma nu)]
		 * 
		 * 
		 */
		
		for(int i = 0;i<pDensityMatrix.rows();i++){
			for(int j=0;j<pDensityMatrix.columns();j++){
				gMatrix.set(i, j, calcGelement(orbitals,i,j,pDensityMatrix));
			}
		}
		
	}

	public static double calcGelement(GaussianOrbital[] orbitals, int i,
			int j, DoubleMatrix2D pDensityMatrix) {
		// TODO Auto-generated method stub
		
		/*
		 * Gmu,nu=SIGM[Aalpa to N/2]SIGMA[lamda to K]SIGMA[sigma to K] C[lamda,.....blah balah ablah
		 * 
		 *       =SIMGA[lamda to K][sigma to K](P[lamda,sigma]*[(mu nu|sigma lamda)-1/2(mu lamda|sigma nu)]
		 * 
		 * 
		 */
		//System.out.println(" orb length"+orbitals.length);
		double sum=0;
		for(int k=0;k<orbitals.length;k++){
			for(int l=0;l<orbitals.length;l++){
				//System.out.println(" P[k][l]"+pDensityMatrix.get(k, l)+" 4CI: "+calculateFourCenterIntegrals(orbitals,i,j,k,l)+" 4CI"+calculateFourCenterIntegrals(orbitals,i,l,k,j));
				sum+=pDensityMatrix.get(k, l)*(calculateFourCenterIntegrals(orbitals,i,j,k,l)-(0.5*calculateFourCenterIntegrals(orbitals,i,l,k,j)));
			}
		}
		//System.out.println(" sum: "+sum);
		return sum;
	}

	public static DoubleMatrix2D cannonicalOrthogonalization(DoubleMatrix2D SMatrix) {
		//based on the book method for 2x2 matrices
		
		//X=Us^-1/2
		
		
		DoubleMatrix2D SdiagMatrix = diagonalize2x2MatrixgetEigenValues(SMatrix);
		DoubleMatrix2D UMatrix = diagonalize2x2MatrixgetEigenVectors(SMatrix);
		System.out.println("===SdiagMatrix===");
		System.out.println(SdiagMatrix);
		
		Algebra alg = new Algebra();
		
		inversesqrt(SdiagMatrix);
		DoubleMatrix2D XMatrix = alg.mult(UMatrix, SdiagMatrix);
		
		return XMatrix;
	}
	
	

	public static double[][][][] calculate4centerintegrals(
			GaussianOrbital[] orbitals) {
		
		//this is a test function
		double[][][][] fourcenterintegrals = new double[orbitals.length][orbitals.length][orbitals.length][orbitals.length];


		
		for(int i=0;i<fourcenterintegrals.length;i++){
			for(int j=0;j<fourcenterintegrals[i].length;j++){
				for(int k=0;k<fourcenterintegrals[i][j].length;k++){
					for(int l=0;l<fourcenterintegrals[i][j][k].length;l++){
						fourcenterintegrals[i][j][k][l]=calculateFourCenterIntegrals(orbitals,i,j,k,l);
					}
				}
			}
			
		}
		return fourcenterintegrals;
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
	
	

	public static double calculateFourCenterIntegral(double alpha, double beta,
			double gamma, double delta, double Rabsqrd,
			double Rcdsqrd, double Rpqsqrd) {
		//this does the actual calculation of the integral
		
		return 2*Math.pow(Math.PI, 2.5)
				/((alpha+beta)*(gamma+delta)*Math.sqrt(alpha+beta+gamma+delta))
				*Math.exp(-(alpha*beta/(alpha+beta)*Rabsqrd)-(gamma*delta/(gamma+delta)*Rcdsqrd))
				*F0((alpha+beta)*(gamma+delta)/(alpha+beta+gamma+delta)*Rpqsqrd);
		
	}

	/*public static DoubleMatrix2D symetricOrthogonalization(
			DoubleMatrix2D soverlapintegrals) {
		DoubleMatrix2D Sdiagonalized;
		DoubleMatrix2D UnitaryMatrixforS;
		//DoubleMatrix2D adjointunitarymatrix;
		
		
		Sdiagonalized = DiagonalizeMatrix(soverlapintegrals);
		if(DEBUG)System.out.println("===diagonalized S=====");
		if(DEBUG)System.out.println(Sdiagonalized);
		if(DEBUG)System.out.println("======================");
		
		
		
		
		
		UnitaryMatrixforS = Sdiagonalized;
		
		//adjointunitarymatrix = adjoint(UnitaryMatrixforS);
		
		if(DEBUG)System.out.println("===Sdiag^-1/2=========");
		inversesqrt(Sdiagonalized);
		if(DEBUG)System.out.println(Sdiagonalized);
		
		
		Algebra alg = new Algebra();
		return alg.mult(UnitaryMatrixforS, Sdiagonalized);
	}*/

	/*public static DoubleMatrix2D DiagonalizeMatrix(
			DoubleMatrix2D matrix) {
		if(matrix.rows()==2&&matrix.columns()==2){
			//to keep this consistent with the book in the 2x2 matrix example
			return diagonalize2x2MatrixgetEigenValues(matrix);
		}
		EigenvalueDecomposition evd = new EigenvalueDecomposition(matrix);
		
		if(DEBUG)System.out.println("===diagonalized S=====");
		if(DEBUG)System.out.println(evd);
		if(DEBUG)System.out.println("======================");
		//if(DEBUG)System.out.println(evd.getD);
		return evd.getD();

		
	}*/
	
	
	public static DoubleMatrix2D diagonalize2x2MatrixgetEigenValues(DoubleMatrix2D matrix){
		//based on the diagonalization method in the book
		
		DoubleMatrix2D eigenvalues = new DenseDoubleMatrix2D(2, 2);
		double theta;
		
		if(Math.abs(matrix.get(0, 0)-matrix.get(1,1))<0.0000000001) theta = Math.PI/4;//tests for if the matrix is symmetric
		else theta = 0.5*Math.atan(2*matrix.get(0, 1)/(matrix.get(0, 0)-matrix.get(1,1)));
		
		eigenvalues.set(0, 0, 
				matrix.get(0, 0)*Math.cos(theta)*Math.cos(theta)
				+matrix.get(1, 1)*Math.sin(theta)*Math.sin(theta)
				+matrix.get(0, 1)*Math.sin(2*theta));
		eigenvalues.set(1, 1, 
				matrix.get(0, 0)*Math.sin(theta)*Math.sin(theta)
				+matrix.get(1, 1)*Math.cos(theta)*Math.cos(theta)
				-matrix.get(0, 1)*Math.sin(2*theta));
		eigenvalues.set(1,0,0);
		eigenvalues.set(0,1,0);
		
		return eigenvalues;
	}
	
	
	public static DoubleMatrix2D diagonalize2x2MatrixgetEigenVectors(DoubleMatrix2D matrix){
		DoubleMatrix2D eigenvectors = new DenseDoubleMatrix2D(2, 2);
		double theta;
		System.out.println(" "+matrix.get(0, 1)+" "+matrix.get(0, 0)+" "+matrix.get(1,1));

		if(Math.abs(matrix.get(0, 0)-matrix.get(1,1))<0.0000001) theta = Math.PI/4;//tests for if the matrix is symmetric
		else theta = 0.5*Math.atan(2*matrix.get(0, 1)/(matrix.get(0, 0)-matrix.get(1,1)));
		
		System.out.println(" theta: "+theta);
		
		eigenvectors.set(0, 0, Math.cos(theta));
		eigenvectors.set(1,0,Math.sin(theta));
		eigenvectors.set(0,1,Math.sin(theta));
		eigenvectors.set(1,1,-Math.cos(theta));
		
		return eigenvectors;
	}
	

	public static DoubleMatrix2D adjoint(DoubleMatrix2D unitaryMatrixforS) {
		
		/*
		 *  adjoint def
		 *  
		 *  (Adagger)i,j = complexconjugate(Aj,i)
		 * 
		 */
		DoubleMatrix2D adjointmatrix = new DenseDoubleMatrix2D(unitaryMatrixforS.rows(), unitaryMatrixforS.columns());
		for(int i=0;i<unitaryMatrixforS.rows();i++){
			for(int j=0;j<unitaryMatrixforS.columns();j++){
				adjointmatrix.set(i, j, unitaryMatrixforS.get(j, i));
			}
		}
		return adjointmatrix;
	}

	public static void inversesqrt(DoubleMatrix2D sdiagonalized) {
		
		for(int i=0;i<sdiagonalized.columns();i++){
			sdiagonalized.set(i, i, 1/Math.sqrt(sdiagonalized.get(i, i)));
		}
	}

	public static DenseDoubleMatrix2D CalculateHCore(
			GaussianOrbital[] orbitals, Atom[] atoms, DoubleMatrix2D T,  DoubleMatrix2D Vnuc) {
		/*
		 * Hcore matrix = KxK matrix
		 * 
		 * Hcore[mu,nu]=T[mu,nu]+Vnucl[mu,nu] 
		 * 
		 * 
		 */
		
		DenseDoubleMatrix2D Hcorematrix = new DenseDoubleMatrix2D(orbitals.length, orbitals.length);
		
		//filling up matrix
		for(int i=0;i<Hcorematrix.rows();i++){
			for(int j=0;j<Hcorematrix.columns();j++){
				//if(DEBUG)System.out.println("calc S mu: "+i+" nu: "+j);
				//Hcore[mu,nu]=T[mu,nu]+Vnucl[mu,nu] 
				//if(DEBUG)System.out.println(" Hcore "+(i+1)+" "+(j+1)+"  T: "+calcTkineticenergyintegral(orbitals,i,j)+"  Vnuc: "+calcVnucintegral(orbitals,atoms,i,j));
				T.set(i, j, calcTkineticenergyintegral(orbitals,i,j));
				Vnuc.set(i, j, calcVnucintegral(orbitals,atoms,i,j));
				Hcorematrix.set(i, j, calcTkineticenergyintegral(orbitals,i,j)+calcVnucintegral(orbitals,atoms,i,j));
				
			}
				
		}
		return Hcorematrix;
		
		
		
		
		
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
					double Vnuctemp = calculateVnucIntegral(orbitals[a].alpha[i],
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
				
		double V = 2*Math.PI/(alpha+beta)*F0((alpha+beta)*Rcpsqrd)*Math.exp(-alpha*beta*Rabsqrd/(alpha+beta));
		double Vnuc=-V*Zelement;
		return Vnuc;
	}

	public static double F0(double d) {
		
		//calculates the f function for F0 only s type
		if(d<0.000001)return 1-(d/3);
		
		return Math.pow(Math.PI/d, 0.5)*ErrorFunction.erf(Math.pow(d, 0.5))/2;
	}
	
	/*public static double derf(double i){
		double P = 0.3275911;
		double[] A = {0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429}; 
		double t = 1/(1+P*i);
		double Tn = t;
		double poly =A[0]*Tn;
		for(int j=1;j<A.length;j++){
			Tn = Tn*t;
			poly += A[j]*Tn; 
		}
		
		
		return 1-poly*Math.exp(i*i);
		
	}*/

	public static double calcTkineticenergyintegral(GaussianOrbital[] orbitals,
			int a, int b) {
		
		//this scans through the contracted gaussian orbitals
		
		double T = 0; //kinetic energy sum
		for(int i=0;i<orbitals[a].alpha.length;i++){//scans through contracted gassian orbital a's simple orbitals
			for(int j=0;j<orbitals[b].alpha.length;j++){
				T += calculateTkineticenergyIntegral(orbitals[a].alpha[i],
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
			 
			 alpha 1s STO-3G
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

	

	public static DoubleMatrix2D EstimateDenstiyMatrix(
			GaussianOrbital[] orbitals, int numberofelectrons) {
		/*  Pu,v = 2*sumovera(Cu,a*Cv,a)
		 * 
		 *  Density matrix elements are formed from the expansion cooefficent matrix C
		 *  also sometimes known as the charge density bond order matrix
		 *  
		 *  for now this just puts in 0
		 */
		
		DoubleMatrix2D P = new DenseDoubleMatrix2D(orbitals.length, orbitals.length);
		P.assign(0.0);
		return P;
	}

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
				
				Soverlapintegrals.set(i, j, CalculateSintegral(orbitals,i,j));
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
				S += calculateSIntegral(orbitals[a].alpha[i],
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
	

	
	
	

}
