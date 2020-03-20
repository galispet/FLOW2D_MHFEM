#pragma once

#include "matrix.h"

#include <iostream>


#define MEMORY_ALIGNEMENT 64

template<unsigned N1, unsigned N2, unsigned N3>
class CoeffMatrix3D {

	template<typename real>
	friend class Matrix;

	template<unsigned N1, unsigned N2, unsigned N3>
	friend class CoeffMatrixOfElement;

public:

	CoeffMatrix3D() {

		data = NULL;
		length = 0;

	};
	~CoeffMatrix3D() {

		if (data)
			_mm_free(data);
			//delete[] data;

	};

	void setNumberOfElements(unsigned const numElements) {

		length = numElements * (N1 * N2 * N3);
		//data = new double[length];

		data = (double *)_mm_malloc(numElements * (N1 * N2 * N3) * sizeof(double), MEMORY_ALIGNEMENT);

	};
	void setZero() {

		for (unsigned i = 0; i < length; i++)
			data[i] = 0.0;

	};

	inline double & setCoeff(unsigned const elementIndex, unsigned const n1, unsigned const n2, unsigned const n3) {

		return data[n3 + N3 * (n2 + N2 * (n1 + N1 * elementIndex))];

	};
	inline double operator()(unsigned const elementIndex, unsigned const n1, unsigned const n2, unsigned const n3) const {

		return data[n3 + N3 * (n2 + N2 * (n1 + N1 * elementIndex))];

	};


	// Operators
	void operator=(CoeffMatrix3D const & m) {

		for (unsigned i = 0; i < length; i++)
			this->data[i] = m.data[i];

	};

	
	// There could be problem, which matrix to extract if (n1,n2) or (n1,n3) or (n2,n3) ... define, which degree of freedom to hold const -> maybe put extra template argument as usnigned
	//template<typename real>
	//Matrix<real> getMatrix(unsigned const elementIndex, unsigned const n1);


private:

	double * data;
	unsigned length;

};

template<unsigned N1, unsigned N2>
class CoeffMatrix2D {

	template<typename real>
	friend class Matrix;

public:

	CoeffMatrix2D() {

		data = NULL;
		length = 0;

	};
	~CoeffMatrix2D() {

		if (data)
			_mm_free(data);
			//delete[] data;

	};

	void setNumberOfElements(unsigned const numElements) {

		length = numElements * (N1 * N2);
		//data = new double[length];

		data = (double *)_mm_malloc(numElements * (N1 * N2) * sizeof(double), MEMORY_ALIGNEMENT);

	};
	void setZero() {

		for (unsigned i = 0; i < length; i++)
			data[i] = 0.0;

	};

	inline double & setCoeff(unsigned const elementIndex, unsigned const n1, unsigned const n2) {

		return data[n2 + N2 * (n1 + N1 * elementIndex)];

	};
	inline double operator()(unsigned const elementIndex, unsigned const n1, unsigned const n2) const {

		return data[n2 + N2 * (n1 + N1 * elementIndex)];

	};

	void operator=(CoeffMatrix2D const & m) {

		for (unsigned i = 0; i < length; i++)
			this->data[i] = m.data[i];

	};


private:

	double * data;
	unsigned length;

};

template<unsigned N1, unsigned numElements = 0>
class CoeffMatrix1D {

	template<typename real>
	friend class Matrix;

public:


	CoeffMatrix1D() {

		data = NULL;
		length = 0;

		if (numElements) {

			length = numElements * (N1);
			data = new double[length];

		}

	};
	~CoeffMatrix1D() {

		if (data)
			_mm_free(data);

			//delete[] data;

	};

	void setNumberOfElements(unsigned const numElements) {

		length = numElements * (N1);
		//data = new double[length];

		data = (double *)_mm_malloc(numElements * (N1) * sizeof(double), MEMORY_ALIGNEMENT);

	};
	void setZero() {

		for (unsigned i = 0; i < length; i++)
			data[i] = 0.0;

	};

	inline double & setCoeff(unsigned const elementIndex, unsigned const n1) {

		return data[n1 + N1 * elementIndex];

	};
	inline double operator()(unsigned const elementIndex, unsigned const n1) const {

		return data[n1 + N1 * elementIndex];

	};

	void operator=(CoeffMatrix1D const & m) {

		for (unsigned i = 0; i < length; i++)
			this->data[i] = m.data[i];

	};


private:

	double * data;
	unsigned length;

};


//
//template<unsigned N1, unsigned N2, unsigned N3>
//class Struct3  {
//
//	template<unsigned N1, unsigned N2, unsigned N3>
//	class SoACoeffMatrix3D;
//
//public:
//
//	double * const Data1;
//	double * const Data2;
//	double * const Data3;
//
//	/*Struct3() : Data1((double *)_mm_malloc(N1 * sizeof(double), MEMORY_ALIGNEMENT)),
//				Data2((double *)_mm_malloc(N2 * sizeof(double), MEMORY_ALIGNEMENT)),
//				Data3((double *)_mm_malloc(N3 * sizeof(double), MEMORY_ALIGNEMENT))
//	{};*/
//	Struct3() : Data1(new double[N1]),
//		Data2(new double[N2]),
//		Data3(new double[N3])
//	{};
//	~Struct3() {
//
//		_mm_free(Data1);
//		_mm_free(Data2);
//		_mm_free(Data3);
//
//	};
//
//};


template<unsigned N1, unsigned N2, unsigned N3>
struct Struct3 {


	//double * const Data1 = (double *)_mm_malloc(N1 * sizeof(double), MEMORY_ALIGNEMENT);
	//double * const Data2 = (double *)_mm_malloc(N2 * sizeof(double), MEMORY_ALIGNEMENT);
	//double * const Data3 = (double *)_mm_malloc(N3 * sizeof(double), MEMORY_ALIGNEMENT);

	//double * const Data1 = new double[N1];
	//double * const Data2 = new double[N2];
	//double * const Data3 = new double[N3];

	double Data1[N1];
	double Data2[N2];
	double Data3[N3];

	~Struct3() {

		//_mm_free(Data1);
		//_mm_free(Data2);
		//_mm_free(Data3);

	};

};


template<unsigned N1, unsigned N2, unsigned N3>
class SoACoeffMatrix3D {

	template<unsigned N1, unsigned N2, unsigned N3>
	friend class CoeffMatrixOfElement;

public:

	SoACoeffMatrix3D() {

	};
	~SoACoeffMatrix3D() {

		//if (data)
		//	_mm_free(data);
		////delete[] data;

	};

	
	void setNumberOfElements(unsigned const numElements) {

		length = numElements * (N1 * N2 * N3);
		data = new Struct3<N1, N2, N3>[numElements];
		Maybe this all is not needed. Data are aleradz good aligned (array of SoA) in container i currently have .. check it
		double * a = &(data->Data1[0]);
		double * b = &(data->Data1[0]);
		double * c = &(data->Data3[0]);

		//__m256d reg = _mm256_loadu_pd(a);

		int differenceInBytesAB = (b - a) * sizeof(double);
		int differenceInBytesAC = (c - a) * sizeof(double);
		int differenceInBytesBC = (c - b) * sizeof(double);

	};

	//void setZero() {
	//
	//	for (unsigned i = 0; i < length; i++)
	//		data[i] = 0.0;
	//
	//};
	//
	//inline double & setCoeff(unsigned const elementIndex, unsigned const n1, unsigned const n2, unsigned const n3) {
	//
	//	return data[n3 + N3 * (n2 + N2 * (n1 + N1 * elementIndex))];
	//
	//};
	//inline double operator()(unsigned const elementIndex, unsigned const n1, unsigned const n2, unsigned const n3) const {
	//
	//	return data[n3 + N3 * (n2 + N2 * (n1 + N1 * elementIndex))];
	//
	//};

private:

	Struct3<N1, N2, N3> * data = NULL;

	unsigned length = 0;

};
//
//template<unsigned N1, unsigned N2>
//class SoACoeffMatrix2D {
//
//	template<typename real>
//	friend class Matrix;
//
//public:
//
//	CoeffMatrix2D() {
//
//		data = NULL;
//		length = 0;
//
//	};
//	~CoeffMatrix2D() {
//
//		if (data)
//			_mm_free(data);
//		//delete[] data;
//
//	};
//
//	void setNumberOfElements(unsigned const numElements) {
//
//		length = numElements * (N1 * N2);
//		//data = new double[length];
//
//		data = (double *)_mm_malloc(numElements * (N1 * N2) * sizeof(double), MEMORY_ALIGNEMENT);
//
//	};
//	void setZero() {
//
//		for (unsigned i = 0; i < length; i++)
//			data[i] = 0.0;
//
//	};
//
//	inline double & setCoeff(unsigned const elementIndex, unsigned const n1, unsigned const n2) {
//
//		return data[n2 + N2 * (n1 + N1 * elementIndex)];
//
//	};
//	inline double operator()(unsigned const elementIndex, unsigned const n1, unsigned const n2) const {
//
//		return data[n2 + N2 * (n1 + N1 * elementIndex)];
//
//	};
//
//	void operator=(CoeffMatrix2D const & m) {
//
//		for (unsigned i = 0; i < length; i++)
//			this->data[i] = m.data[i];
//
//	};
//
//
//private:
//
//	double * data;
//	unsigned length;
//
//};
//
//template<unsigned N1, unsigned numElements = 0>
//class SoACoeffMatrix1D {
//
//	template<typename real>
//	friend class Matrix;
//
//public:
//
//
//	CoeffMatrix1D() {
//
//		data = NULL;
//		length = 0;
//
//		if (numElements) {
//
//			length = numElements * (N1);
//			data = new double[length];
//
//		}
//
//	};
//	~CoeffMatrix1D() {
//
//		if (data)
//			_mm_free(data);
//
//		//delete[] data;
//
//	};
//
//	void setNumberOfElements(unsigned const numElements) {
//
//		length = numElements * (N1);
//		//data = new double[length];
//
//		data = (double *)_mm_malloc(numElements * (N1) * sizeof(double), MEMORY_ALIGNEMENT);
//
//	};
//	void setZero() {
//
//		for (unsigned i = 0; i < length; i++)
//			data[i] = 0.0;
//
//	};
//
//	inline double & setCoeff(unsigned const elementIndex, unsigned const n1) {
//
//		return data[n1 + N1 * elementIndex];
//
//	};
//	inline double operator()(unsigned const elementIndex, unsigned const n1) const {
//
//		return data[n1 + N1 * elementIndex];
//
//	};
//
//	void operator=(CoeffMatrix1D const & m) {
//
//		for (unsigned i = 0; i < length; i++)
//			this->data[i] = m.data[i];
//
//	};
//
//
//private:
//
//	double * data;
//	unsigned length;
//
//};























template<unsigned N1, unsigned N2, unsigned N3>
class SingleCoeffMatrix3D {

public:

	SingleCoeffMatrix3D() : length(N1 * N2 * N3), data(new double[N1 * N2 * N3]) {

	};
	~SingleCoeffMatrix3D() {

		delete[] data;

	};

	void setZero() {

		for (unsigned i = 0; i < length; i++)
			data[i] = 0.0;

	};

	inline double & setCoeff(unsigned const n1, unsigned const n2, unsigned const n3) {

		return data[n3 + N3 * (n2 + N2 * n1)];

	};
	inline double operator()(unsigned const n1, unsigned const n2, unsigned const n3) const {

		return data[n3 + N3 * (n2 + N2 * n1)];

	};

	void operator=(SingleCoeffMatrix3D const & m) {

		for (unsigned i = 0; i < length; i++)
			this->data[i] = m.data[i];

	};


private:

	double * const data;
	unsigned const length;

};

template<unsigned N1, unsigned N2>
class SingleCoeffMatrix2D {

public:

	SingleCoeffMatrix2D() : length(N1 * N2), data(new double[N1 * N2]) {

	};
	~SingleCoeffMatrix2D() {

		delete[] data;

	};

	void setZero() {

		for (unsigned i = 0; i < length; i++)
			data[i] = 0.0;

	};

	inline double & setCoeff(unsigned const n1, unsigned const n2) {

		return data[n2 + N2 * n1];

	};
	inline double operator()(unsigned const n1, unsigned const n2) const {

		return data[n2 + N2 * n1];

	};

	void operator=(SingleCoeffMatrix2D const & m) {

		for (unsigned i = 0; i < length; i++)
			this->data[i] = m.data[i];

	};


private:

	double * const data;
	unsigned const length;

};







template<unsigned N1, unsigned N2, unsigned N3>
class CoeffMatrixOfElement : public CoeffMatrix<N1, N2, N3> {

	//template<typename real>
	//friend class Matrix;

public:

	CoeffMatrixOfElement(CoeffMatrix3D<N1, N2, N3> & other, unsigned const elementIndex);
	~CoeffMatrixOfElement();


	inline double & setCoeff(unsigned const n1, unsigned const n2, unsigned const n3);
	inline double operator()(unsigned const n1, unsigned const n2, unsigned const n3) const;



private:

	double * data;

};


template<unsigned N1, unsigned N2, unsigned N3>
CoeffMatrixOfElement<N1, N2, N3>::CoeffMatrixOfElement(CoeffMatrix3D<N1, N2, N3> & other, unsigned const elementIndex) {

	data = other.data + elementIndex * (N1*N2*N3);

};
template<unsigned N1, unsigned N2, unsigned N3>
CoeffMatrixOfElement<N1, N2, N3>::~CoeffMatrixOfElement() {

	// Hope it doesnt NULL other.data :)) I don't have such a profound knowledge of c++
	data = NULL;

};
template<unsigned N1, unsigned N2, unsigned N3>
inline double & CoeffMatrixOfElement<N1, N2, N3>::setCoeff(unsigned const n1, unsigned const n2, unsigned const n3) {

	return data[n3 + N3 * (n2 + N2 * n1)];

};
template<unsigned N1, unsigned N2, unsigned N3>
inline double CoeffMatrixOfElement<N1, N2, N3>::operator()(unsigned const n1, unsigned const n2, unsigned const n3) const {

	return data[n3 + N3 * (n2 + N2 * n1)];

};









/*
template<unsigned N1, unsigned N2, unsigned N3>
template<typename real>
Matrix<real> CoeffMatrix<N1, N2, N3>::getMatrix(unsigned const elementIndex, unsigned const n1) {


	Matrix<real> result(N2, N3);

	unsigned const beginRow = 0;
	unsigned const beginCol = 0;

	double * const tempData = &this->data[beginCol + N3 * (beginRow + N2 * (n1 + N1 * elementIndex))];
	//double * const tempData = coeffMat.data[n3 + N3 * (n2 + N2 * (n1 + N1 * elementIndex))];

	for (unsigned i = 0; i < N2*N3; i++)
		result.data[i] = tempData[i];

	return result;




	//Matrix<real> JF2 = affineMappingMatrixJF.getMatrix<real>(k_index, 0);

	//std::cout << JF(0, 0) << "  " << JF2(0, 0) << std::endl;
	//std::cout << JF(0, 1) << "  " << JF2(0, 1) << std::endl;
	//std::cout << JF(1, 0) << "  " << JF2(1, 0) << std::endl;
	//std::cout << JF(1, 1) << "  " << JF2(1, 1) << std::endl;
	//std::cout << std::endl;

};
*/

/*

//template<class T, unsigned ... RestDimension>
//class MultiDimensionalArray;
//
//template<class T, unsigned FirstDimension>
//class MultiDimensionalArray<T, FirstDimension> {
//
//	typedef T type[FirstDimension];
//	type data;
//	T & operator[](unsigned i) { return data[i]; }
//
//};
//
//template<class T, unsigned FirstDimension, unsigned ... RestDimension>
//class MultiDimensionalArray<T, FirstDimension, RestDimension...> {
//
//	typedef typename MultiDimensionalArray<T, RestDimension ...>::type OneDimensionDownArrayT;
//	typedef OneDimensionDownArrayT type[FirstDimension];
//	type data;
//	OneDimensionDownArrayT & operator[](unsigned i) { return data[i]; }
//
//};
//
//
//
//
//
//template<class T, unsigned ... RestDimension>
//class MultiDimensionalArray;
//
//template<class T, unsigned FirstDimension>
//class MultiDimensionalArray<T, FirstDimension> {
//
//
//
//	T * data;
//	T & operator[](unsigned i) { return data[i]; }
//
//};


//template<class T, unsigned FirstDimension, unsigned ... RestDimension>
//class MultiDimensionalArray<T, FirstDimension, RestDimension...> {
//
//public:
//
//	MultiDimensionalArray() {
//
//		unsigned length = FirstDimension + sum(RestDimension);
//
//		data = new T[length];
//
//	};
//	~MultiDimensionalArray() {
//
//
//		delete[] data;
//
//	}
//
//
//private:
//
//	T * data;
//
//
//	template<unsigned K>
//	unsigned sum(K i) { return i; };
//
//	template<unsigned K, unsigned ... Rest>
//	unsigned sum(K i, Rest ... rest) { return i + sum(rest...); }
//
//};


//template<class T, unsigned ... RestDimension>
//class MultiDimensionalArray;
//
//template<class T, unsigned FirstDimension>
//class MultiDimensionalArray<T, FirstDimension> {
//
//	typedef T type[FirstDimension];
//	type data;
//	T & operator[](unsigned i) { return data[i]; }
//
//};
//
//template<class T, unsigned FirstDimension, unsigned ... RestDimension>
//class MultiDimensionalArray<T, FirstDimension, RestDimension...> {
//
//public:
//
//	MultiDimensionalArray() {
//
//		unsigned length = sizeof(T)
//
//	};
//	~MultiDimensionalArray() {
//
//
//		delete[] data;
//
//	}
//
//
//private:
//
//	T * data;
//
//
//	template<unsigned K>
//	unsigned sum(K i) { return i; };
//
//	template<unsigned K, unsigned ... Rest>
//	unsigned sum(K i, Rest ... rest) { return i + sum(rest...); }
//
//};
*/