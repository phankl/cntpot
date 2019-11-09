#include "linalg.h"

namespace linalg{
	/****************************************************
	QR Factorisation with Modified Gram-Schmidt algorithm
	
	Input: quadratic matrix a
	Output: orthogonal matrix q, upper triangular matrix r
	****************************************************/
	void qrFactorisation(
		const std::vector<std::vector<double>>& a, 
		std::vector<std::vector<double>>& q, 
		std::vector<std::vector<double>>& r)
	{
		int n = a.size();
		q.resize(n);
		r.resize(n);
		for(int i = 0; i < n; i++){
			q[i].resize(n);
			r[i].resize(n);			
		}

		std::vector<std::vector<double>> aTranspose(n, std::vector<double>(n, 0));
		std::vector<std::vector<double>> qTranspose(n, std::vector<double>(n, 0));
		std::vector<double> temp(n, 0);
		transpose(a, aTranspose);

		std::vector<std::vector<double>> u(n, std::vector<double>(n, 0));

		//Compute u vectors with modified Gram-Schmidt algorithm
		for(int i = 0; i < n; i++){
			u[i] = aTranspose[i];
		}
		for(int i = 0; i < n; i++){
			r[i][i] = norm(u[i]);
			if(r[i][i] != 0) scale(1.0/r[i][i],u[i],qTranspose[i]);
			for(int j = i+1; j < n; j++){
				r[i][j] = dot(qTranspose[i],u[j]);
				scale(r[i][j],qTranspose[i],temp);
				subtract(u[j],temp,u[j]);
			}
		}
		transpose(qTranspose,q);
	};

	/****************************************************
	QR Eigenvalue Algorithm
	
	Input: quadratic matrix a
	Output: list of eigenvalues, list of eigenvectors
	****************************************************/
	void qrEigenvalue(
		const std::vector<std::vector<double>>& a,
		std::vector<double>& eigenvalues,
		std::vector<std::vector<double>>& eigenvectors
		)
	{
		int n = a.size();
		std::vector<std::vector<double>> q(n,std::vector<double>(n,0));
		std::vector<std::vector<double>> r(n,std::vector<double>(n,0));
		std::vector<std::vector<double>> diff(n,std::vector<double>(n,0));
		std::vector<std::vector<double>> qTotal(n,std::vector<double>(n,0));

		unityMatrix(qTotal);

		std::vector<std::vector<double>> d = a;
		std::vector<std::vector<double>> temp = a;

		double error;
		double errorLimit = 1.0e-18;
		do{
			//Apply Wilkinson shift
			double diag1 = d[n-2][n-2];
			double diag2 = d[n-1][n-1];
			double b = d[n-2][n-1];
			double bSquared = b * b;
			double delta = 0.5 * (diag2 - diag1);
			
			int sign;
			if(delta > 0) sign = 1;
			else sign = -1;
			
			double shift;
			if(delta == 0) shift = diag2 - fabs(b);
			else shift = diag2 - bSquared/(delta + sign*sqrt(delta*delta + bSquared));

			std::vector<std::vector<double>> shiftMatrix(n,std::vector<double>(n,0));
			unityMatrix(shiftMatrix);
			scaleMatrix(shift,shiftMatrix,shiftMatrix);
			subtractMatrix(d,shiftMatrix,d);

			qrFactorisation(d,q,r);
			matrixProduct(r,q,d);
			addMatrix(d,shiftMatrix,d);
			matrixProduct(qTotal,q,qTotal);

			//subtractMatrix(d,temp,diff);
			//error = diagonalNormMatrix(diff) / n;
			//temp = d;

			//Check convergence
			for(int i = 0; i < n-1; i++){
				error = 2 * fabs(d[i][i+1]) / (fabs(d[i][i]) + fabs(d[i+1][i+1]));
				if(error < errorLimit){ 
					d[i][i+1] = 0;
					d[i+1][i] = 0;
				}
			}
			error = 2 * fabs(d[0][1]) / (fabs(d[0][0]) + fabs(d[1][1]));
		}
		while(error > errorLimit);

		eigenvalues.resize(n);
		eigenvectors.resize(n);
		for(int i = 0; i < n; i++){
			eigenvalues[i] = d[i][i];
			eigenvectors[i].resize(n);
		}
		transpose(qTotal,eigenvectors);
	}

	/****************************************************
	Quadratic Matrix Transpose
	
	Input: quadratic matrix a
	Output: transposed matrix b
	****************************************************/
	void transpose(
		const std::vector<std::vector<double>>& a,
		std::vector<std::vector<double>>& b
		)
	{
		int n = a.size();
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				b[i][j] = a[j][i];
			}
		}
	}

	/****************************************************
	Matrix Product
	
	Input: matrices a, b
	Output: matrix product c
	****************************************************/
	void matrixProduct(
		const std::vector<std::vector<double>>& a,
		const std::vector<std::vector<double>>& b,
		std::vector<std::vector<double>>& c
		)
	{
		int m = a.size();
		int n = b[0].size();
		int o = a[0].size();

		std::vector<std::vector<double>> d(m, std::vector<double>(n, 0));

		//Compute product
		for(int i = 0; i < m; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < o; k++){
					d[i][j] += a[i][k] * b[k][j];
				}
			}
		}

		c = d;
	}


	/****************************************************
	Vector Addition
	
	Input: vectors a, b
	Output: vector c
	****************************************************/
	void add(
		const std::vector<double>& a,
		const std::vector<double>& b,
		std::vector<double>& c
		)
	{
		int n = a.size();
		for(int i = 0; i < n; i++){
			c[i] = a[i] + b[i];
		}
	}

	/****************************************************
	Vector Subtraction
	
	Input: vectors a, b
	Output: vector c
	****************************************************/
	void subtract(
		const std::vector<double>& a,
		const std::vector<double>& b,
		std::vector<double>& c
		)
	{
		int n = a.size();
		for(int i = 0; i < n; i++){
			c[i] = a[i] - b[i];
		}
	}

	/****************************************************
	Normalise vector
	
	Input: vector a
	Output: normalised vector b
	****************************************************/
	void normalise(
		const std::vector<double>& a,
		std::vector<double>& b
		)
	{
		int n = a.size();
		double aNorm = norm(a);
		if(aNorm != 0) scale(1.0/aNorm,a,b);
		else b = a;
	}

	/****************************************************
	Scale vector
	
	Input: scale factor s, vector a
	Output: scaled vector b
	****************************************************/
	void scale(
		double s,
		const std::vector<double>& a,
		std::vector<double>& b
		)
	{
		int n = a.size();
		for(int i = 0; i < n; i++){
			b[i] = s * a[i];
		}
	}

	/****************************************************
	Vector Dot Product
	
	Input: vectors a, b
	Output: dot product
	****************************************************/
	double dot(
		const std::vector<double>& a,
		const std::vector<double>& b
		)
	{
		int n = a.size();
		double result = 0;
		for(int i = 0; i < n; i++){
			result += a[i]*b[i];
		}

		return result;
	}


	/****************************************************
	Vector Norm
	
	Input: vector a
	Output: vector norm
	****************************************************/
	double norm(
		const std::vector<double>& a
		)
	{
		return sqrt(dot(a,a));
	}

	/****************************************************
	Matrix Addition
	
	Input: matrices a, b
	Output: matrix c
	****************************************************/
	void addMatrix(
		const std::vector<std::vector<double>>& a,
		const std::vector<std::vector<double>>& b,
		std::vector<std::vector<double>>& c
		)
	{
		int n = a.size();
		int m = a[0].size();
		for(int i = 0; i < n; i++){
			for(int j = 0; j < m; j++){
				c[i][j] = a[i][j] + b[i][j];
			}
		}
	}

	/****************************************************
	Matrix Scalar Multiplication
	
	Input: scale factor s, matrix a
	Output: matrix b
	****************************************************/
	void scaleMatrix(
		double s,
		const std::vector<std::vector<double>>& a,
		std::vector<std::vector<double>>& b
		)
	{
		int n = a.size();
		int m = a[0].size();
		for(int i = 0; i < n; i++){
			for(int j = 0; j < m; j++){
				b[i][j] = s * a[i][j];
			}
		}
	}



	/****************************************************
	Matrix Subtraction
	
	Input: matrices a, b
	Output: matrix c
	****************************************************/
	void subtractMatrix(
		const std::vector<std::vector<double>>& a,
		const std::vector<std::vector<double>>& b,
		std::vector<std::vector<double>>& c
		)
	{
		int n = a.size();
		int m = a[0].size();
		for(int i = 0; i < n; i++){
			for(int j = 0; j < m; j++){
				c[i][j] = a[i][j] - b[i][j];
			}
		}
	}

	/****************************************************
	Unity Matrix
	
	Input: matrix a
	Output: unity matrix a
	****************************************************/
	void unityMatrix(
		std::vector<std::vector<double>>& a
		)
	{
		int n = a.size();
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				if(i == j) a[i][j] = 1;
				else a[i][j] = 0;
			}
		}
	}


	/****************************************************
	Matrix Dot Product
	
	Input: matrices a, b
	Output: dot product
	****************************************************/
	double dotMatrix(
		const std::vector<std::vector<double>>& a,
		const std::vector<std::vector<double>>& b
		)
	{
		int n = a.size();
		int m = a[0].size();
		double result = 0;
		for(int i = 0; i < n; i++){
			for(int j = 0; j < m; j++){
				result += a[i][j]*b[i][j];
			}
		}

		return result;
	}


	/****************************************************
	Frobenius Norm of a Matrix
	
	Input: matrix a
	Output: frobenius norm
	****************************************************/
	double normMatrix(
		const std::vector<std::vector<double>>& a
		)
	{
		return sqrt(dotMatrix(a,a));
	}

	/****************************************************
	Diagonal Norm of a Matrix
	
	Input: matrix a
	Output: frobenius norm
	****************************************************/
	double diagonalNormMatrix(
		const std::vector<std::vector<double>>& a
		)
	{
		double sum = 0;
		int n = a.size();
		for(int i = 0; i < n; i++){
			sum += a[i][i] * a[i][i];
		}

		return sqrt(sum);
	}

}
