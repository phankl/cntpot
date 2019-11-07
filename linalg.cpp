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
		transpose(a, aTranspose);

		std::vector<std::vector<double>> u(n, std::vector<double>(n, 0));

		//Compute u vectors with modified Gram-Schmidt algorithm
		for(int i = 0; i < n; i++){
			u[i] = aTranspose[i];
			for(int j = 0; j < i; j++){
				std::vector<double> projectionVector(n, 0);
				projection(u[i],u[j],projectionVector);
				subtract(u[i],projectionVector,u[i]);
			}
			normalise(u[i],u[i]);
		}

		//Get Q and R
		transpose(u,q);
		for(int i = 0; i < n; i++){
			for(int j = i; j < n; j++){
				r[i][j] = dot(u[i],aTranspose[j]);
			}
		}
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
		std::vector<std::vector<double>> d(n, std::vector<double>(n, 0));
		std::vector<std::vector<double>> q(n, std::vector<double>(n, 0));
		std::vector<std::vector<double>> r(n, std::vector<double>(n, 0));
		std::vector<std::vector<double>> diff(n, std::vector<double>(n, 0));

		qrFactorisation(a,q,r);
		std::vector<std::vector<double>> qTotal = q;
		std::vector<std::vector<double>> temp = a;

		double error = 1.0;
		double errorLimit = 1.0e-8;

		while(error > errorLimit){
			matrixProduct(r,q,d);
			qrFactorisation(d,q,r);

			matrixProduct(qTotal,q,qTotal);

			subtractMatrix(d,temp,diff);
			error = normMatrix(diff);
			temp = d;
		}

		eigenvalues.resize(n);
		eigenvectors.resize(n);
		for(int i = 0; i < n; i++){
			eigenvalues[i] = d[i][i];
			eigenvectors[i].resize(n);
		}
		transpose(qTotal,eigenvectors);
	}

	/****************************************************
	Gram-Schmidt Projection
	
	Input: vectors a, u
	Output: projection vector p
	****************************************************/
	void projection(
		const std::vector<double>& a,
		const std::vector<double>& u,	
		std::vector<double>& r
		)
	{
		int n = a.size();

		double uNorm = dot(u,u);
		double fraction; 
		if(uNorm == 0) fraction = 0;
		else fraction = dot(u,a) / uNorm;
		
		r.resize(n);
		for(int i = 0; i < n; i++){
			r[i] = fraction * u[i];
		}
	};

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
}
