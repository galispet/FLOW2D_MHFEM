#pragma once


#include <iostream>

#include <Eigen/Sparse>
#include <Eigen/Dense>

template<typename real>
class Vector;




template<typename real, unsigned m>
class SquareMatrix {

public:

	real * const data;
	unsigned const n = m;

	SquareMatrix() : data(new real[m*m])  {

		//this->data = new real[n*n];
		this->setZero();

	};
	SquareMatrix(SquareMatrix const & mat) : data(new real[m*m]) {

		for (unsigned i = 0; i < n*n; i++)
			this->data[i] = mat.data[i];

	};
	~SquareMatrix() {

		if (this->data)
			delete[] this->data;

	};


	inline real operator()(unsigned const i, unsigned const j) const {

		return this->data[j + m * i];

	};
	inline real & operator()(unsigned i, unsigned j) {

		return this->data[j + m * i];

	};
	inline SquareMatrix & operator=(SquareMatrix const & mat) {

		for (unsigned i = 0; i < n*m; i++)
			this->data[i] = mat.data[i];

		return * this;

	};
	inline SquareMatrix operator*(SquareMatrix const & mat) const {


		SquareMatrix<real, m> result;

		for (unsigned i = 0; i < m; i++) {
			for (unsigned j = 0; j < m; j++) {

				real s = 0.0;

				for (unsigned k = 0; k < m; k++)
					s += this->data[k + m * i] * mat.data[j + m * k];

				result.data[j + m * i] = s;

			}
		}

		return result;

	};


	inline void setZero() {

		for (unsigned i = 0; i < n*n; i++)
			this->data[i] = 0.0;

	};
	inline void inverseInPlace() {

		if (m == 2) {

			real const aa = this->data[0];
			real const bb = this->data[1];
			real const cc = this->data[2];
			real const dd = this->data[3];

			real const iiDet = (real) 1.0 / (aa*dd - bb * cc);

			this->data[0] = +iiDet * dd;
			this->data[1] = -iiDet * bb;
			this->data[2] = -iiDet * cc;
			this->data[3] = +iiDet * aa;

			return;

		}
		if (m == 3) {

			real const a = this->data[0];
			real const b = this->data[1];
			real const c = this->data[2];
			real const d = this->data[3];
			real const e = this->data[4];
			real const f = this->data[5];
			real const g = this->data[6];
			real const h = this->data[7];
			real const i = this->data[8];

			real const iDet = (real) 1.0 / (a*(e*i - f * h) - b * (d*i - f * g) + c * (d*h - e * g));

			this->data[0] = iDet * (e*i - f * h);
			this->data[1] = iDet * (c*h - b * i);
			this->data[2] = iDet * (b*f - c * e);
			this->data[3] = iDet * (f*g - d * i);
			this->data[4] = iDet * (a*i - c * g);
			this->data[5] = iDet * (c*d - a * f);
			this->data[6] = iDet * (d*h - e * g);
			this->data[7] = iDet * (b*g - a * h);
			this->data[8] = iDet * (a*e - b * d);

			return;

		}

		/*
		switch (m) {

		case 2: 

			real const aa = this->data[0];
			real const bb = this->data[1];
			real const cc = this->data[2];
			real const dd = this->data[3];

			real const iiDet = (real) 1.0 / (aa*dd - bb * cc);

			this->data[0] = +iiDet * dd;
			this->data[1] = -iiDet * bb;
			this->data[2] = -iiDet * cc;
			this->data[3] = +iiDet * aa;

			return;

		case 3:

			real const a = this->data[0];
			real const b = this->data[1];
			real const c = this->data[2];
			real const d = this->data[3];
			real const e = this->data[4];
			real const f = this->data[5];
			real const g = this->data[6];
			real const h = this->data[7];
			real const i = this->data[8];

			real const iDet = (real) 1.0 / (a*(e*i - f * h) - b * (d*i - f * g) + c * (d*h - e * g));

			this->data[0] = iDet * (e*i - f * h);
			this->data[1] = iDet * (c*h - b * i);
			this->data[2] = iDet * (b*f - c * e);
			this->data[3] = iDet * (f*g - d * i);
			this->data[4] = iDet * (a*i - c * g);
			this->data[5] = iDet * (c*d - a * f);
			this->data[6] = iDet * (d*h - e * g);
			this->data[7] = iDet * (b*g - a * h);
			this->data[8] = iDet * (a*e - b * d);

			return;

		}
		*/

	};
	inline SquareMatrix inverse() {


		SquareMatrix<real, m> result;

		if (m == 2) {

			real const aa = this->data[0];
			real const bb = this->data[1];
			real const cc = this->data[2];
			real const dd = this->data[3];

			real const iiDet = (real) 1.0 / (aa*dd - bb * cc);

			result.data[0] = +iiDet * dd;
			result.data[1] = -iiDet * bb;
			result.data[2] = -iiDet * cc;
			result.data[3] = +iiDet * aa;


			return result;

		}
		if (m == 3) {

			real const a = this->data[0];
			real const b = this->data[1];
			real const c = this->data[2];
			real const d = this->data[3];
			real const e = this->data[4];
			real const f = this->data[5];
			real const g = this->data[6];
			real const h = this->data[7];
			real const i = this->data[8];

			real const iDet = (real) 1.0 / (a*(e*i - f * h) - b * (d*i - f * g) + c * (d*h - e * g));

			result.data[0] = iDet * (e*i - f * h);
			result.data[1] = iDet * (c*h - b * i);
			result.data[2] = iDet * (b*f - c * e);
			result.data[3] = iDet * (f*g - d * i);
			result.data[4] = iDet * (a*i - c * g);
			result.data[5] = iDet * (c*d - a * f);
			result.data[6] = iDet * (d*h - e * g);
			result.data[7] = iDet * (b*g - a * h);
			result.data[8] = iDet * (a*e - b * d);

			return result;

		}

		return result;

	};


};



template<typename real>
class Matrix {

	template<unsigned N1, unsigned N2, unsigned N3>
	friend class CoeffMatrix;

public:

	real * data;
	unsigned numRows;
	unsigned numCols;

	Matrix() {
	
		this->data = NULL;
	
		this->numRows = 0;
		this->numCols = 0;
	
	};
	Matrix(unsigned rows, unsigned cols) {
	//
		this->data = new real[rows*cols];
		//
		this->numRows = rows;
		this->numCols = cols;
		//
		this->setZero();
		//
	};
	Matrix(Matrix const & mat) {


		unsigned const n = mat.numRows;
		unsigned const m = mat.numCols;

		this->data = new real[n*m];

		this->numRows = n;
		this->numCols = m;

		for (unsigned i = 0; i < n*m; i++)
			this->data[i] = mat.data[i];

	};
	virtual ~Matrix() {

		// Base class destructor is implemnted as virtual -> derived destructor is called also
		if (this->data)
			delete[] this->data;

	};


	inline real operator()(unsigned i, unsigned j) const {

		return this->data[j + numCols * i];

	};
	inline real & operator()(unsigned i, unsigned j) {

		return this->data[j + numCols * i];

	};
	inline Matrix & operator=(Matrix const & mat) {


		if (this == &mat)
			return *this;



		unsigned const n = mat.numRows;
		unsigned const m = mat.numCols;

		//delete[] this->data;
		//this->data = new real[n*m];

		this->numRows = n;
		this->numCols = m;

		for (unsigned i = 0; i < n*m; i++)
			this->data[i] = mat.data[i];

		return *this;

	};
	inline Matrix operator*(Matrix const & mat) const {


		unsigned const n = this->numRows;
		unsigned const m = this->numCols;

		unsigned const on = mat.numRows;
		unsigned const om = mat.numCols;

		Matrix result(n, om);


		for (unsigned i = 0; i < n; i++) {

			for (unsigned j = 0; j < om; j++) {

				real s = 0.0;

				for (unsigned k = 0; k < m; k++)
					s += this->data[k + m * i] * mat.data[j + om * k];

				result.data[j + om * i] = s;

			}
		}

		return result;

	};
	inline Vector<real> operator*(Vector<real> const & vec) const {


		unsigned const n = this->numRows;
		unsigned const m = this->numCols;

		Vector<real> result(n);


		for (unsigned i = 0; i < n; i++) {

			real s = 0.0;
			for (unsigned k = 0; k < m; k++)
				s += vec.data[k] * this->data[k + m * i];

			result.data[i] = s;

		}

		return result;

	};
	inline Matrix operator*(real const & a) const {


		unsigned const n = this->numRows;
		unsigned const m = this->numCols;

		Matrix result(n, m);

		for (unsigned i = 0; i < n*m; i++)
			result->data[i] = a * this->data[i];

		return result;

	};
	friend inline Matrix operator*(real const & a, Matrix const & mat) {


		unsigned const n = mat.numRows;
		unsigned const m = mat.numCols;

		Matrix result(n, m);

		for (unsigned i = 0; i < n*m; i++)
			result.data[i] = a * mat.data[i];

		return result;

	};
	inline Matrix operator/(real const & a) const {


		unsigned const n = this->numRows;
		unsigned const m = this->numCols;

		Matrix result(n, m);

		for (unsigned i = 0; i < n*m; i++)
			result->data[i] = this->data[i] / a;

		return result;

	};

	inline void setZero() {

		for (unsigned i = 0; i < this->numRows*this->numCols; i++)
			this->data[i] = 0.0;

	};
	inline Vector<real> getColumn(unsigned j) const {


		Vector<real> result(this->numRows);

		for (unsigned i = 0; i < this->numRows; i++)
			result.data[i] = this->data[j + i * this->numCols];

		return result;

	};
	inline void inverseInPlace() {

		if (this->numCols == 2 && this->numRows == 2) {

			real const a = this->data[0];
			real const b = this->data[1];
			real const c = this->data[2];
			real const d = this->data[3];

			real const iDet = (real) 1.0 / (a*d - b * c);

			this->data[0] =  iDet * d;
			this->data[1] = -iDet * b;
			this->data[2] = -iDet * c;
			this->data[3] =  iDet * a;
			
			return;

		}
		//else if (this->numCols == 3 && this->numRows == 3) {
		//
		//
		//	real const a = this->data[0];
		//	real const b = this->data[1];
		//	real const c = this->data[2];
		//	real const d = this->data[3];
		//	real const e = this->data[4];
		//	real const f = this->data[5];
		//	real const g = this->data[6];
		//	real const h = this->data[7];
		//	real const i = this->data[8];
		//
		//	real const iDet = (real) 1.0 / (a*(e*i - f * h) - b * (d*i - f * g) + c * (d*h - e * g));
		//
		//	this->data[0] = iDet * (e*i - f * h);
		//	this->data[1] = iDet * (c*h - b * i);
		//	this->data[2] = iDet * (b*f - c * e);
		//	this->data[3] = iDet * (f*g - d * i);
		//	this->data[4] = iDet * (a*i - c * g);
		//	this->data[5] = iDet * (c*d - a * f);
		//	this->data[6] = iDet * (d*h - e * g);
		//	this->data[7] = iDet * (b*g - a * h);
		//	this->data[8] = iDet * (a*e - b * d);
		//
		//	return;
		//
		//}

		//Eigen::Matrix<real, 8, 8> I = Eigen::Matrix<real, 8, 8>::Identity();
		//mat = mat.llt().solve(I);

		Eigen::Map<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix(this->data, this->numRows, this->numCols);
		matrix = matrix.inverse();
		
		for (unsigned i = 0; i < this->numRows*this->numCols; i++)
			this->data[i] = matrix.data()[i];
	
	};
	inline Matrix inverse() {


		unsigned const n = this->numRows;
		unsigned const m = this->numRows;

		Matrix result(n, m);

		if (n == 2 && m == 2) {

			real const a = this->data[0];
			real const b = this->data[1];
			real const c = this->data[2];
			real const d = this->data[3];

			real const idet = (real) 1.0 / (a*d - b * c);

			result.data[0] = idet * d;
			result.data[1] = -idet * b;
			result.data[2] = -idet * c;
			result.data[3] = idet * a;

			return result;

		}

		real * temp = new real[n*m];

		for (unsigned i = 0; i < n*m; i++)
			temp[i] = this->data[i];


		Eigen::Map<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix(temp, n, m);
		matrix = matrix.inverse();

		for (unsigned i = 0; i < n*m; i++)
			result.data[i] = matrix.data()[i];

		delete[] temp;

		return result;

	};
	inline void transposeInPlace() {

		for (unsigned i = 0; i < this->numRows; i++)
			for (unsigned j = i + 1; j < this->numCols; j++) {

				real const temp = data[j + this->numRows * i];
				this->data[j + this->numRows * i] = this->data[i + this->numRows * j];
				this->data[i + this->numRows * j] = temp;

			}

				

	};


	friend std::ostream & operator<<(std::ostream & ost, Matrix<real> const & m) {

		Eigen::MatrixXd temp(m.numRows, m.numCols);

		for (unsigned i = 0; i < m.numRows; i++)
			for (unsigned j = 0; j < m.numCols; j++)
				temp.coeffRef(i, j) = m(i, j);
		
		ost << temp << std::endl;

		return ost;

	};

};


template<typename real>
class Vector : public Matrix<real> {

private:

	using Matrix<real>::getColumn;

public:

	Vector() {

		this->data = NULL;

		this->numRows = 0;
		this->numCols = 0;

	};
	Vector(unsigned rows) {

		this->data = new real[rows];

		this->numRows = rows;
		this->numCols = 1;

		this->setZero();


	};
	Vector(Vector const & vec) {

		this->numRows = vec.numRows;

		this->data = new real[this->numRows];

		for (unsigned i = 0; i < this->numRows; i++)
			this->data[i] = vec.data[i];

	};
	//~Vector() {


	//	// Because Base class Matrix will clean it.
	//	// Base class destructor is implemnted as virtual -> derived destructor is called also, but we do not need it here

	//	//if (this->data)
	//	//	delete[] this->data;

	//};

	inline real operator()(unsigned i) const {

		return this->data[i];

	};
	inline real & operator()(unsigned i) {

		return this->data[i];

	};
	inline Vector & operator=(Vector const & vec) {


		if (this == &vec)
			return *this;

		this->numRows = vec.numRows;

		for (unsigned i = 0; i < vec.numRows; i++)
			this->data[i] = vec.data[i];

		return *this;

	};
	inline Vector operator*(Matrix<real> const & mat) const {


		unsigned const n = this->numRows;
		unsigned const om = mat->numCols;

		Vector result(n);


		for (unsigned i = 0; i < n; i++) {

			real s = 0.0;
			for (unsigned k = 0; k < n; k++)
				s += mat.data[k + om * i] * this->data[k];

			result.data[i] = s;

		}

		return result;

	};
	inline Vector operator*(real const & a) const {


		unsigned const n = this->numRows;

		Vector result(n);

		for (unsigned i = 0; i < n; i++)
			result.data[i] = a * this->data[i];

		return result;

	};
	friend inline Vector operator*(real const & a, Vector const & vec) {


		unsigned const n = vec.numRows;

		Vector result(n);

		for (unsigned i = 0; i < n; i++)
			result.data[i] = a * vec.data[i];

		return result;

	};
	inline Vector operator/(real const & a) const {


		unsigned const n = this->numRows;

		Vector result(n);

		for (unsigned i = 0; i < n; i++)
			result.data[i] = this->data[i] / a;

		return result;

	};

	inline real norm() const {

		real s = 0.0;

		for (unsigned k = 0; k < this->numRows; k++)
			s += (real) (this->data[k] * this->data[k]);

		return (real) sqrt(s);

	};

};

template<typename real>
inline real dot(Vector<real> const & v, Vector<real> const & w) {


	real s = 0.0;

	for (unsigned k = 0; k < v.numRows; k++)
		s += v.data[k] * w.data[k];

	return s;

};
template<typename real>
inline real dot(Vector<real> & v, Vector<real> & w) {


	real s = 0.0;

	for (unsigned k = 0; k < v.numRows; k++)
		s += v.data[k] * w.data[k];

	return s;

};