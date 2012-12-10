
package cern.colt.matrix.linalg; 
 
import cern.colt.matrix.DoubleFactory1D; 
import cern.colt.matrix.DoubleFactory2D; 
import cern.colt.matrix.DoubleMatrix1D; 
import cern.colt.matrix.DoubleMatrix2D; 
import cern.jet.math.Functions; 
 
/** 
 * @author <a href="razvan.surdulescu@post.harvard.edu">Razvan Surdulescu</a>, Copyright (c) 2003 
 */ 
public class LanczosDecomposition { 
    private static final double EPSILON = 1E-6; 
 
    private boolean m_log; 
 
    private DoubleMatrix1D m_eigenvalues; 
    private DoubleMatrix2D m_eigenvectors; 
 
    public LanczosDecomposition(DoubleMatrix2D matrix) { 
        this(matrix, matrix.columns(), false); 
    } 
 
    public LanczosDecomposition(DoubleMatrix2D matrix, int sought) { 
        this(matrix, sought, false); 
    } 
 
    public LanczosDecomposition(DoubleMatrix2D matrix, boolean log) { 
        this(matrix, matrix.columns(), log); 
    } 
 
    public LanczosDecomposition(DoubleMatrix2D matrix, int sought, boolean log) { 
        Property.DEFAULT.checkSquare(matrix); 
 
        if (!Property.DEFAULT.isSymmetric(matrix)) { 
            throw new IllegalArgumentException("Matrix must be symmetric: " + cern.colt.matrix.doublealgo.Formatter.shape(matrix)); 
        } 
 
        m_log = log; 
 
        lanczos(matrix, sought); 
    } 
 
    /** 
     * Returns the eigenvalues. 
     */ 
    public DoubleMatrix1D getRealEigenvalues() { 
        return m_eigenvalues; 
    } 
 
    /** 
     * Returns the eigenvector matrix. 
     */ 
    public DoubleMatrix2D getV() { 
        return m_eigenvectors; 
    } 
 
    /** 
     * Iterative Lanczos algorithm for finding the approximate eigenvalues and 
     * eigenvectors of a matrix. 
     * 
     * This implementation was pieced together using information from: 
     * <ul> 
     * <li>Lecture notes from <a href="http://www.cs.berkeley.edu/~demmel/cs267/lecture20/lecture20.html"> 
     * CS267</a> at Berkeley.</li> 
     * <li>Axel Ruhe's <a href="http://www.cs.utk.edu/~dongarra/etemplates/node103.html">description</a> 
     * of the Lanczos Algorithm.</li> 
     * </ul> 
     * 
     * @param input The matrix whose eigenvalues we wish to compute. 
     * @param sought The number of eigenvalues we wish to compute for the input matrix. 
     * This number cannot exceed the number of columns for the matrix. 
     */ 
    private void lanczos(DoubleMatrix2D input, int sought) { 
        int n = input.columns(); 
 
        if (sought > n) { 
            throw new IllegalArgumentException("Number of eigensolutions sought cannot exceed size of the matrix: " + 
                    sought + " > " + n); 
        } 
 
        // random vector 
        DoubleMatrix1D r = DoubleFactory1D.dense.make(n); 
        for (int i = 0; i < n; i++) { 
            r.set(i, Math.random()); 
        } 
        if (m_log) { 
            System.out.println("Random vector: " + r); 
        } 
 
        // initial size for algorithm data structures 
        // this may grow after a number of iterations 
        int size = sought; 
 
        // diagonal elements of T 
        DoubleMatrix1D a = DoubleFactory1D.dense.make(size); 
 
        // off-diagonal elements of T 
        DoubleMatrix1D b = DoubleFactory1D.dense.make(size); 
 
        // basis vectors for the Krylov subspace 
        DoubleMatrix2D v = DoubleFactory2D.dense.make(n, size); 
 
        // arrays used in the QL decomposition 
        double[] d = new double[size + 1]; 
        double[] e = new double[size + 1]; 
        double[][] z = new double[size + 1][size + 1]; 
 
        // which ritzvalues are converged? 
        boolean[] c = new boolean[size]; 
 
        // how many converged eigenvalues have been found? 
        int found = 0; 
 
        // algorithm iterations 
        int i = 0; 
 
        while (true) { 
            if (m_log) { 
                System.out.println("Iteration=" + i); 
            } 
 
            // v(i) = r / b(i-1) 
            if (i == 0) { 
                v.viewColumn(i).assign(r).assign(Functions.div(Math.sqrt(Algebra.DEFAULT.norm2(r)))); 
            } else if (b.get(i - 1) != 0) { 
                v.viewColumn(i).assign(r).assign(Functions.div(b.get(i - 1))); 
            } 
            if (m_log) { 
                System.out.println("v(" + i + ")=" + v.viewColumn(i)); 
            } 
 
            // r = A * v(i) 
            r.assign(Algebra.DEFAULT.mult(input, v.viewColumn(i))); 
            if (m_log) { 
                System.out.println("r = A * v(" + i + ") = " + r); 
            } 
 
            // r = r - b(i-1) * v(i-1) 
            if (i == 0) { 
                // v(i-1) = 0, so r is unchanged in this case 
            } else { 
                DoubleMatrix1D t1 = v.viewColumn(i - 1).copy().assign(Functions.mult(b.get(i - 1))); 
                r.assign(t1, Functions.minus); 
            } 
            if (m_log) { 
                System.out.println("r = r - b(" + (i - 1) + ") * v(" + (i - 1) + ") = " + r); 
            } 
 
            // a(i) = v(i)' * r 
            a.set(i, Algebra.DEFAULT.mult(v.viewColumn(i), r)); 
            if (m_log) { 
                System.out.println("a(" + i + ") = v(" + i + ")' * r = " + a); 
            } 
 
            // r = r - a(i)*v(i) 
            DoubleMatrix1D t2 = v.viewColumn(i).copy().assign(Functions.mult(a.get(i))); 
            r.assign(t2, Functions.minus); 
            if (m_log) { 
                System.out.println("r = r - a(" + i + ") * v(" + i + ") = " + r); 
            } 
 
            // TODO: re-orthogonalize if necessary 
 
            // b(i) = norm(r) 
            b.set(i, Math.sqrt(Algebra.DEFAULT.norm2(r))); 
            if (m_log) { 
                System.out.println("b(" + i + ") = norm(r) = " + b); 
            } 
 
            // prepare to compute the eigenvalues and eigenvectors 
            // of the tridiagonal matrix defined by a and b 
            System.arraycopy(a.toArray(), 0, d, 1, i + 1); 
            System.arraycopy(b.toArray(), 0, e, 2, i); 
 
            for (int p = 1; p <= i + 1; p++) { 
                for (int q = 1; q <= i + 1; q++) { 
                    if (p == q) { 
                        z[p][q] = 1; 
                    } else { 
                        z[p][q] = 0; 
                    } 
                } 
            } 
 
            // compute the eigenvalues and eigenvectors of the 
            // tridiagonal matrix 
            tqli(d, e, i + 1, z); 
 
            // count the number of converged eigenvalues 
            found = 0; 
            for (int j = 0; j <= i; j++) { 
                double ritzvalue = d[j + 1]; 
                double residual = Math.abs(b.get(i) * z[i + 1][j + 1]); 
 
                if (residual <= EPSILON) { 
                    c[j] = true; 
                    found++; 
                } else { 
                    c[j] = false; 
                } 
 
                if (m_log) { 
                    System.out.println("Ritz value[" + j + "]=" + ritzvalue + ", residual[" + j + "]=" + residual); 
                } 
            } 
 
            if (found >= sought) { 
                break; 
            } 
 
            i++; 
 
            if (i >= n) { 
                break; 
            } else if (i >= size) { 
                if (m_log) { 
                    System.out.println("Growing arrays from " + size + " to " + (2 * size)); 
                } 
 
                size = 2 * size; 
 
                a = grow(a, size); 
                b = grow(b, size); 
                v = grow(v, size); 
                d = grow(d, size + 1); 
                e = grow(e, size + 1); 
                z = grow(z, size + 1); 
                c = grow(c, size + 1); 
            } 
        } 
 
        m_eigenvalues = DoubleFactory1D.dense.make(found); 
        m_eigenvectors = DoubleFactory2D.dense.make(n, found); 
        DoubleMatrix2D ritzvectors = DoubleFactory2D.dense.make(z).viewPart(1, 1, size, size); 
 
        int index = 0; 
        for (int col = 0; col < i; col++) { 
            if (c[col]) { 
                m_eigenvalues.set(index, d[col + 1]); 
                m_eigenvectors.viewColumn(index).assign(Algebra.DEFAULT.mult(v, ritzvectors.viewColumn(col))); 
                index++; 
            } 
        } 
        if (m_log) { 
            System.out.println("Eigenvalues: " + m_eigenvalues); 
            System.out.println("Eigenvectors: " + m_eigenvectors); 
        } 
    } 
 
    private DoubleMatrix1D grow(DoubleMatrix1D matrix, int length) { 
        DoubleMatrix1D result = DoubleFactory1D.dense.make(length); 
        for (int index = 0; index < matrix.size(); index++) { 
            result.setQuick(index, matrix.getQuick(index)); 
        } 
        return result; 
    } 
 
    private DoubleMatrix2D grow(DoubleMatrix2D matrix, int columns) { 
        DoubleMatrix2D result = DoubleFactory2D.dense.make(matrix.rows(), columns); 
        for (int col = 0; col < matrix.columns(); col++) { 
            result.viewColumn(col).assign(matrix.viewColumn(col)); 
        } 
        return result; 
    } 
 
    private double[] grow(double[] array, int length) { 
        double[] result = new double[length]; 
        System.arraycopy(array, 0, result, 0, array.length); 
        return result; 
    } 
 
    private double[][] grow(double[][] array, int length) { 
        double[][] result = new double[length][length]; 
        for (int row = 0; row < array.length; row++) { 
            System.arraycopy(array[row], 0, result[row], 0, array.length); 
        } 
        return result; 
    } 
 
    private boolean[] grow(boolean[] array, int length) { 
        boolean[] result = new boolean[length]; 
        System.arraycopy(array, 0, result, 0, array.length); 
        return result; 
    } 
 
    /** 
     * Return the absolute value of a with the same sign as b. 
     */ 
    private double sign(double a, double b) { 
        return (b >= 0.0 ? Math.abs(a) : -Math.abs(a)); 
    } 
 
    /** 
     * Returns sqrt(a^2 + b^2) without under/overflow. 
     */ 
    private double pythag(double a, double b) { 
        double r; 
        if (Math.abs(a) > Math.abs(b)) { 
            r = b / a; 
            r = Math.abs(a) * Math.sqrt(1 + r * r); 
        } else if (b != 0) { 
            r = a / b; 
            r = Math.abs(b) * Math.sqrt(1 + r * r); 
        } else { 
            r = 0.0; 
        } 
        return r; 
    } 
 
    /** 
     * "Tridiagonal QL Implicit" routine for computing eigenvalues and eigenvectors of a symmetric, 
     * real, tridiagonal matrix. 
     * 
     * The routine works extremely well in practice. The number of iterations for the first few 
     * eigenvalues might be 4 or 5, say, but meanwhile the off-diagonal elements in the lower right-hand 
     * corner have been reduced too. The later eigenvalues are liberated with very little work. The 
     * average number of iterations per eigenvalue is typically 1.3 - 1.6. The operation count per 
     * iteration is O(n), with a fairly large effective coefficient, say, ~20n. The total operation count 
     * for the diagonalization is then ~20n * (1.3 - 1.6)n = ~30n^2. If the eigenvectors are required, 
     * there is an additional, much larger, workload of about 3n^3 operations. 
     * 
     * This implementation is taken directly from "Numerical Recipes in C" 
     * <a href="http://www.ulib.org/webRoot/Books/Numerical_Recipes/bookcpdf/c11-3.pdf">Chapter 11.3</a> 
     * and translated to Java. 
     * 
     * @param d [1..n] array. On input, it contains the diagonal elements of the tridiagonal matrix. 
     * On output, it contains the eigenvalues. 
     * @param e [1..n] array. On input, it contains the subdiagonal elements of the tridiagonal 
     * matrix, with e[1] arbitrary. On output, its contents are destroyed. 
     * @param n The size of all parameter arrays. 
     * @param z [1..n][1..n] array. On input, it contains the identity matrix. On output, the kth column 
     * of z returns the normalized eigenvector corresponding to d[k]. 
     */ 
    private void tqli(double d[], double e[], int n, double[][] z) { 
        int i; 
        // Convenient to renumber the elements of e. 
        for (i = 2; i <= n; i++) { 
            e[i - 1] = e[i]; 
        } 
        e[n] = 0.0; 
        for (int l = 1; l <= n; l++) { 
            int iter = 0; 
            int m; 
            do { 
                // Look for a single small subdiagonal element to split the matrix. 
                for (m = l; m <= n - 1; m++) { 
                    double dd = Math.abs(d[m]) + Math.abs(d[m + 1]); 
                    if (Math.abs(e[m]) + dd == dd) { 
                        break; 
                    } 
                } 
                if (m != l) { 
                    iter = iter + 1; 
                    if (iter == 30) { 
                        throw new RuntimeException("Too many iterations"); 
                    } 
                    // Form shift. 
                    double g = (d[l + 1] - d[l]) / (2.0 * e[l]); 
                    double r = pythag(g, 1.0); 
                    // This is dm / ks. 
                    g = d[m] - d[l] + e[l] / (g + sign(r, g)); 
                    double s, c; 
                    s = c = 1.0; 
                    double p = 0.0; 
                    // A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form. 
                    for (i = m - 1; i >= l; i--) { 
                        double f = s * e[i]; 
                        double b = c * e[i]; 
                        e[i + 1] = (r = pythag(f, g)); 
                        // Recover from underflow. 
                        if (r == 0.0) { 
                            d[i + 1] -= p; 
                            e[m] = 0.0; 
                            break; 
                        } 
                        s = f / r; 
                        c = g / r; 
                        g = d[i + 1] - p; 
                        r = (d[i] - g) * s + 2.0 * c * b; 
                        d[i + 1] = g + (p = s * r); 
                        g = c * r - b; 
                        // Form eigenvectors (optional). 
                        for (int k = 1; k <= n; k++) { 
                            f = z[k][i + 1]; 
                            z[k][i + 1] = s * z[k][i] + c * f; 
                            z[k][i] = c * z[k][i] - s * f; 
                        } 
                    } 
                    if (r == 0.0 && i >= l) { 
                        continue; 
                    } 
                    d[l] -= p; 
                    e[l] = g; 
                    e[m] = 0.0; 
                } 
            } while (m != l); 
        } 
    } 
    public String toString() {
    	StringBuffer buf = new StringBuffer();
    	String unknown = "Illegal operation or error: ";

    	buf.append("---------------------------------------------------------------------\n");
    	buf.append("LANCZOSDecomposition(A) --> D, real eigen values \n");
    	buf.append("---------------------------------------------------------------------\n");

    	buf.append("realEigenvalues = ");
    	try { buf.append(String.valueOf(this.getRealEigenvalues()));} 
    	catch (IllegalArgumentException exc) { buf.append(unknown+exc.getMessage()); }
    		
    	/*buf.append("\nimagEigenvalues = ");
    	try { buf.append(String.valueOf(this.getImagEigenvalues()));} 
    	catch (IllegalArgumentException exc) { buf.append(unknown+exc.getMessage()); }
    		
    	buf.append("\n\nD = ");
    	try { buf.append(String.valueOf(this.getD()));} 
    	catch (IllegalArgumentException exc) { buf.append(unknown+exc.getMessage()); }*/
    	
    	buf.append("\n\nV = ");
    	try { buf.append(String.valueOf(this.getV()));} 
    	catch (IllegalArgumentException exc) { buf.append(unknown+exc.getMessage()); }
    	
    	return buf.toString();
    }
} 

