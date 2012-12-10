package jackSolver;


import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;

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
		
		HFSolver(atoms,numberofelectrons,2);
		
	}

	public static double HFSolver(Atom[] atoms, int numberofelectrons, int verbose) {
		
		/*
		 * verbose = 0 no text
		 *         = 1 basic info no matrices
		 *         = 2 matrices
		 *         =
		
		 * 
		 */
		
		
		int STOnG = 3;
		
		
		
		
		//Writing Information about the input
		if(verbose>0)System.out.println("JackRHFSolver by Andrew Long");
		
		if(numberofelectrons%2!=0){
			if(verbose>0)System.out.println("uneven number of electrons");
			return 0;
		}
		if(verbose>0)System.out.println("===List of Atoms=============");
		if(verbose>0){
			for(int i=0;i<atoms.length;i++){
				System.out.println(" Atom: "+i+" Z: "+atoms[i].element+" x: "+atoms[i].x+" y: "+atoms[i].y+" z: "+atoms[i].z);
			}
		}
		if(verbose>0)System.out.println("...Generating Orbital Array");
		
		GaussianOrbital[] orbitals = GaussianOrbital.generateOrbitalArray(atoms,STOnG);
		if(verbose>1)for(int i=0;i<orbitals.length;i++)System.out.println(orbitals[i]);
		
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
		
		sOverlapMatrix = Integrals.GenerateSOverlapIntegrals(orbitals);
		if(verbose>0)System.out.println("===S Matrix==========");
		if(verbose>1)System.out.println(sOverlapMatrix);
		
		T = new DenseDoubleMatrix2D(orbitals.length, orbitals.length);
		Vnuc = new DenseDoubleMatrix2D(orbitals.length, orbitals.length);
		
		hCore = CalculateHCore(orbitals, atoms, T, Vnuc);
		if(verbose>0)System.out.println("===Hcore=============");
		if(verbose>1)System.out.println(hCore);
		if(verbose>0)System.out.println("===T=================");
		if(verbose>1)System.out.println(T);
		if(verbose>0)System.out.println("===Vnuc==============");
		if(verbose>1)System.out.println(Vnuc);
		
		//step #3 diagonalize S obtain  transform matrix X
		
		
		
		
		
		xMatrix = cannonicalOrthogonalization(sOverlapMatrix);
		if(verbose>0)System.out.println("===X matrix===========");
		if(verbose>1)System.out.println(xMatrix);
		
		
		//step #3 four integral debugging step  in the book they calculate these and put the in matrix
		if(verbose>0)System.out.println("===four center integrals===");
		double[][][][]fourcenterintegrals = Integrals.calculate4centerintegrals(orbitals);//this is now done in the test cases
		
		
		
		//step #4 estimate density matrix
		
		if(DEBUG)System.out.println("===P density matrix====");
		pDensityMatrix = EstimateDenstiyMatrix(orbitals,numberofelectrons);
		if(verbose>1)System.out.println(pDensityMatrix);
		
		
		//begin SCF loop
		Algebra alg = new Algebra();
		gMatrix = new DenseDoubleMatrix2D(orbitals.length,orbitals.length);
		fFockMatrix = new DenseDoubleMatrix2D(orbitals.length,orbitals.length);
		
		double currentElectronicEnergy=0;
				
		
		
		int MaxSCFiterations = 500;
		double ConvergenceCriteria = 0.0001;
		
		double currentConvergence = 1.0;
		int currentIteration = 0;
		while(currentIteration<MaxSCFiterations && currentConvergence>ConvergenceCriteria){
			currentIteration++;
			
			if(verbose>1)System.out.println("========================================");
			if(verbose>1)System.out.println("===SCF LOOP=============================");
			if(verbose>1)System.out.println(" current iteration: "+currentIteration);
			
			//step #5 calculate Matrix G
			
			gMatrix.assign(0);
			calcGMatrix(gMatrix, pDensityMatrix,orbitals,fourcenterintegrals);
			if(verbose>1)System.out.println("===G matrix===========");
			if(verbose>1)System.out.println(gMatrix);
			
			
			
			//step #6 calculate fock matrix F=Hcore+G
			calcFockMatrix(fFockMatrix,hCore,gMatrix);
			if(verbose>1)System.out.println("===fFockMatrix========");
			if(verbose>1)System.out.println(fFockMatrix);
			
			//step #6 and 1/3  calculate the current energy
			System.out.println("======================");
			currentElectronicEnergy = calculateElectronicEnergy(pDensityMatrix, hCore, fFockMatrix);
			if(verbose>1)System.out.println(" electronic energy: "+currentElectronicEnergy);
			
			//step #7 calc transformed fock matrix F'=X(+)FX
			xAdjointMatrix = adjoint(xMatrix);
			fPrimeMatrix=alg.mult(xAdjointMatrix, alg.mult(fFockMatrix, xMatrix));
			if(verbose>1)System.out.println("===fPrimeMatrix========");
			if(verbose>1)System.out.println(fPrimeMatrix);
			
			
			
			//step #8 diagonalize F' to get C' and epsilon
			// F'C'=C'epsilon
			
			if(verbose>1)System.out.println("===fPrimeDiagonalizedMatrix========");
			//SingularValueDecomposition svd = new SingularValueDecomposition(fPrimeMatrix);
			EigenvalueDecomposition evd = new EigenvalueDecomposition(fPrimeMatrix);//aka diagonalize
			//fPrimeDiagonalizedMatrix=diagonalize2x2MatrixgetEigenValues(fPrimeMatrix);
			
			fPrimeDiagonalizedMatrix=evd.getD();
			cPrimeMatrix=evd.getV();
			//sortEigenvaluesAndVectors(fPrimeDiagonalizedMatrix,cPrimeMatrix);
			
			if(verbose>1)System.out.println(fPrimeDiagonalizedMatrix);
			
			if(verbose>1)System.out.println("===cPrimeMatrix========");
			//cPrimeMatrix=diagonalize2x2MatrixgetEigenVectors(fPrimeMatrix);
			
			//System.out.println("===old cPrimeMatrix========");
			
			if(verbose>1)System.out.println(cPrimeMatrix);
			//if(verbose>1)System.out.println(evd);
			
			
			
			//cPrimeMatrix=evd.getV();
			if(verbose>1)System.out.println(cPrimeMatrix);
			
			//step #9 calc C=XC'
			if(verbose>1)System.out.println("===cCoefficientsMatrix========");
			cCoefficientsMatrix = alg.mult(xMatrix, cPrimeMatrix);
			if(verbose>1)System.out.println(cCoefficientsMatrix);
			
			//step #10 calc new density matrix
			if(verbose>1)System.out.println("===new pDensityMatrix========");
			oldpDensityMatrix=pDensityMatrix.copy();
			calcNewDensityMatrix(pDensityMatrix,oldpDensityMatrix,cCoefficientsMatrix,numberofelectrons);
			if(verbose>1)System.out.println(pDensityMatrix);
			
			//step #11 test convergence and calc expectation values and other quantities of interest
			
			currentConvergence=calcConvergence(pDensityMatrix,oldpDensityMatrix);
			if(verbose>1)System.out.println("=============================");
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
		if(verbose>1)System.out.println("===Mulliken Populaton========");
		if(verbose>1)System.out.println(alg.mult(pDensityMatrix, sOverlapMatrix));
		
		
		
		
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

	public static void calcNewDensityMatrix(DoubleMatrix2D pDensityMatrix,DoubleMatrix2D oldpDensityMatrix, 
			DoubleMatrix2D cCoefficientsMatrix,  int numberofelectrons) {
		
		pDensityMatrix.assign(0); //zero the density Matrix
		for(int i=0;i<pDensityMatrix.rows();i++){
			for(int j=0;j<pDensityMatrix.columns();j++){
				for(int n=0;n<(numberofelectrons/2);n++){
					//System.out.println(" i"+i+"j"+j+" n"+n+ " currentP "+pDensityMatrix.get(i,j)+" Cin"+cCoefficientsMatrix.get(i, n)+" Cjn"+cCoefficientsMatrix.get(j, n));
					double newpValue = pDensityMatrix.get(i,j)+2*cCoefficientsMatrix.get(i, n)*cCoefficientsMatrix.get(j, n);
					double newaveragedpValue = 0.9*newpValue+0.1*oldpDensityMatrix.get(i,j);
					pDensityMatrix.set(i, j, newaveragedpValue);
				}
			}
		}


		
	}

	private static double calculateElectronicEnergy(
			DoubleMatrix2D pDensityMatrix, DoubleMatrix2D hCore,
			DoubleMatrix2D fFockMatrix) {
		
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
			DoubleMatrix2D pDensityMatrix, GaussianOrbital[] orbitals, double[][][][] fourcenterintegrals) {
		/*
		 * Gmu,nu=SIGM[Aalpa to N/2]SIGMA[lamda to K]SIGMA[sigma to K] C[lamda,.....blah balah ablah
		 * 
		 *       =SIMGA[lamda to K][sigma to K](P[lamda,sigma]*[(mu nu|sigma lamda)-1/2(mu lamda|sigma nu)]
		 * 
		 * 
		 */
		if(fourcenterintegrals==null){
			for(int i = 0;i<pDensityMatrix.rows();i++){
				for(int j=0;j<pDensityMatrix.columns();j++){
					gMatrix.set(i, j, calcGelement(orbitals,i,j,pDensityMatrix));
				}
			}
		}else{
			for(int i = 0;i<pDensityMatrix.rows();i++){
				for(int j=0;j<pDensityMatrix.columns();j++){
					gMatrix.set(i, j, calcGelement(fourcenterintegrals,i,j,pDensityMatrix));
				}
			}
		}
		
		
	}

	private static double calcGelement(double[][][][] fourcenterintegrals,
			int i, int j, DoubleMatrix2D pDensityMatrix) {
		/*
		 * Gmu,nu=SIGM[Aalpa to N/2]SIGMA[lamda to K]SIGMA[sigma to K] C[lamda,.....blah balah ablah
		 * 
		 *       =SIMGA[lamda to K][sigma to K](P[lamda,sigma]*[(mu nu|sigma lamda)-1/2(mu lamda|sigma nu)]
		 * 
		 * 
		 */
		//System.out.println(" orb length"+orbitals.length);
		double sum=0;
		for(int k=0;k<fourcenterintegrals.length;k++){
			for(int l=0;l<fourcenterintegrals.length;l++){
				//System.out.println(" P[k][l]"+pDensityMatrix.get(k, l)+" 4CI: "+calculateFourCenterIntegrals(orbitals,i,j,k,l)+" 4CI"+calculateFourCenterIntegrals(orbitals,i,l,k,j));
				sum+=pDensityMatrix.get(k, l)*(fourcenterintegrals[i][j][k][l]-(0.5*fourcenterintegrals[i][l][k][j]));
			}
		}
		//System.out.println(" sum: "+sum);
		return sum;
	}

	public static double calcGelement(GaussianOrbital[] orbitals, int i,
			int j, DoubleMatrix2D pDensityMatrix) {
		
		
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
				sum+=pDensityMatrix.get(k, l)*(Integrals.calculateFourCenterIntegrals(orbitals,i,j,k,l)-(0.5*Integrals.calculateFourCenterIntegrals(orbitals,i,l,k,j)));
			}
		}
		//System.out.println(" sum: "+sum);
		return sum;
	}

	
	public static DoubleMatrix2D cannonicalOrthogonalization(DoubleMatrix2D SMatrix) {
		//based on the book method for 2x2 matrices
		
		//X=Us^-1/2
		
		EigenvalueDecomposition svd = new EigenvalueDecomposition(SMatrix);
				
		DoubleMatrix2D SdiagMatrix = svd.getD();//diagonalize2x2MatrixgetEigenValues(SMatrix);
		DoubleMatrix2D UMatrix = svd.getV();//diagonalize2x2MatrixgetEigenVectors(SMatrix);
		
		sortEigenvaluesAndVectors(SdiagMatrix,UMatrix);
		
		Algebra alg = new Algebra();
						
		JackRHFSolver.inversesqrt(SdiagMatrix);
		DoubleMatrix2D XMatrix = alg.mult(UMatrix, SdiagMatrix);
		
		return XMatrix;
	}
	
	
	

	public static void sortEigenvaluesAndVectors(DoubleMatrix2D sdiagMatrix,
			DoubleMatrix2D uMatrix) {
		
		for(int i=0;i<sdiagMatrix.columns();i++){
			for(int j=0;j<sdiagMatrix.columns();j++){
				if(sdiagMatrix.get(i, i)>sdiagMatrix.get(j, j)){
					double temp=sdiagMatrix.get(i,i);
					sdiagMatrix.set(i,i,sdiagMatrix.get(j,j));
					sdiagMatrix.set(j, j, temp);
					
					exchangeColumnsOfMatrix(i,j,uMatrix);
				}
			}
		}
		// TODO Auto-generated method stub
		
	}

	public static void exchangeColumnsOfMatrix(int i, int j,
			DoubleMatrix2D matrix) {
		for(int k=0;k<matrix.rows();k++){
			double temp = matrix.get(k, i);
			matrix.set(k, i, matrix.get(k,j));
			matrix.set(k, j, temp);
		}
		
	}

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
				T.set(i, j, Integrals.calcTkineticenergyintegral(orbitals,i,j));
				Vnuc.set(i, j, Integrals.calcVnucintegral(orbitals,atoms,i,j));
				Hcorematrix.set(i, j, Integrals.calcTkineticenergyintegral(orbitals,i,j)+Integrals.calcVnucintegral(orbitals,atoms,i,j));
				
			}
				
		}
		return Hcorematrix;
		
		
		
		
		
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
	

	
	
	

}
