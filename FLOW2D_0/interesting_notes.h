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

