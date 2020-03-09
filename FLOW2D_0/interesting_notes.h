#pragma once

/*

-	float a[] means the same thing as float* a,
	the [] syntax is a hint to the caller of the function
	that it can work on arrays rather than single variables.

-	you need to use (&array) to clarify to the compiler that
	you want a reference to an array, rather than the
	(invalid) array of references int & array[100];.
	"array of references" - which isn't legal.

-	void foo(int * x);
	void foo(int x[100]);
	void foo(int x[]);
	These three are different ways of declaring the same function.
	They're all treated as taking an int * parameter, you can
	pass any size array to them.

*/

/*

-	const reference to an object T (T const & ) -> const is redunant, reerence is always cosnt (it will never be reference to anz=ythig else)

*/

/*

-	can't do  the following templated call of constructor


template<unsigned N1>
class CoeffMatrix1D {


	template<unsigned const numElements>
	CoeffMatrix1D() {

		length = numElements * (N1);
		data = new double[length];

	};

	OR

	template<unsigned const numElements>
	CoeffMatrix1D(unsigned const numElements) {

		length = numElements * (N1);
		data = new double[length];

	};

}

CoeffMatrix1D<3> QuadraturePointsAndWeightsOnReferenceTriangle<4>;

-	is it possible somehow?

*/
