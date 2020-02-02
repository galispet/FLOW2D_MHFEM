#pragma once

/**


	1.  - Big issue !!! 
		- make the 'straight_walk' algorithm robust. When considering large num of vertices - it's the fastest !!!
		- instead of orientation test (at initializing step only?) use 't->neighbor_ccw_cw'
		- Implement Orthogonal walk - may be useful
		- Implement Visibility walk - isn't it lawson which i already have ????
		- Make search a point in triangle more robust
		- try to use point in the center of the triangle?
		- try to parallelize, e.g. each thread with different algorithm or from different starting triangles
		- implement walk on edges : probably pointers to prev and next edges will be needed

	2.	- Why, in some cases, is 'in_circle' test not symmetric? 
		- i.e.   incircle(t,ov) != incircle(ot,v)
		- could it be, because of bounding triangle vertices?

	3.	- Error tolerance should be proportional to the length of edge , area of triangle , ...

	4.	- Make Point on edge more less precise i.e. point will fall on edge if e.g. dist(v,e) < 1.e-3 * length(e)

	5.	- Precomputing determinants can be made with additional data struture, that will be deleted after triangulation => no additional memory is consumed

	6.  - Implement 'infinite bounding triangle'
		- My idiotic solution : if the input vertice is outside bounding triangle -> make the triangle inflate. We have to check coordinates of each inserted point :(

8. Speed up thing consedering 'orient' function to inline everywhere where it can work
	also try to reduce function calling, figure out where it can be eliminated atc.


11. Verify my triangulation with matlab's -> generate same points and compare 


13. MEMORY MANAGEMENT -> allocate all memory at the beggining of triangulation
	vector.reserve doesnt allocate memory
	maybe i have to implement some memory pool which i will work with ???? 



17. When inserting constrained edge, keep the pointer to it. When the edge is refined by inserting point to it, so we can keep track to its subsegments


18. Barycentric refinment (subdivision)

19. Try to normalize given PSLG onto [0,1] x [0,1] for better accuracy. new_x = (old_x-min_x)/d_max, new_y = (old_y-min_y)/d_max, d_max=max{x_max-x_min,y_max-y_min}






*/





// Walking algorithm for constraints -> when hit constraint (hole=no triangle) then move along the constraintsuntil an triangle which is interscted with a ray is found

//'insert_constraint_intersection' make vertices more distributed so they do not form clusters and therefore create good quality triangles

// make laplacian smoothing similiar to the book of 'Jiri Blazek - Computational Fluid Dynamics_ Principles and Applications (2015, Butterworth-Heinemann)' p. 394/388

// insert vertex into edge -> if the redunandt triangles are deleted in first place before refining, there can occur a situtantion, when will be created only 2 triangles and not four
// therefore there must be ensure that i am not calling a method on NULL >]]


// !!!!!!!!!!  Deletion of redundant triangles is connected to super triangle only ! wont work for concave region !!! spred the virus also for the triangles connected to super triangles -> it will work :)
//				this way maybe the convex hull / rather points on the boundary can also be obtained


// get rid of pointers and make primitives verteces and adjacent triangle made of indeces to their respective vertices






/* ------------- CHANGES ------------- */
// 9.7. : 
// - Dumped the super edges from inserting into 'edge_tr' container. In remove_super_triangle, therefore no need to get rid of them explicitly (they would in the first three positions)

// !!! Create adjacency list of triangles to each vertex. Do it when deleting holes and outside triangles, so we can easily rotary traverse around every Vertex so we do not need
//	   to change direction when encounter NULL !!