#pragma once


enum class NUMBERING { CM, RCM, RedBlue, Xmajor, Ymajor };

enum class SMOOTHING { LAPLACIAN };

enum class MESH_QUALITY { MINIMUM_ANGLE, INCIRCLE_RADIUS, SHORTEST_EDGE };


enum class INSERT_CONSTRAINT { ENFORCE, FLIP, HALVING, INTERSECTION };

enum class REFINEMENT { RUPPERT, CHEW, CHEW2 };
enum class REFINEMENT_PRIORITY { WORST, BEST };


enum class LOCATION { ON_VERTEX, ON_EDGE, IN_TRIANGLE, NOT_VISIBLE };

enum class V_MARKER { FREE, CONSTRAINED };			// Vertex marker
enum class E_MARKER { NONE, DIRICHLET, NEUMANN};	// Edge marker
enum class T_MARKER { NONE, OUTSIDE, INSIDE };		// Triangle marker

enum class POLYGON_TYPE { OPEN, CLOSED, HOLE };				// Constraint broken line is 'open', or forms a 'closed' polygon



		/*****************************************************************************/
		/*                                                                           */
		/*    - Exact predicates : geometric predicates such as						 */
		/*					      'in_circle', 'in_triangle', 'orientation', etc...  */
		/*						  are computed exactly (robustly)					 */
		/*																			 */
		/*	  - Exact construct  : geometric constructs such as						 */
		/*						  'edge_intersection', 'barycentric_coordinates',    */
		/*						  are computed exactly (robustly)					 */
		/*						 : more functions will be preseneted with automatic  */
		/*						   grid generation, mesh smoothing etc..	         */
		/*																			 */
		/*																			 */
		/*****************************************************************************/
enum class GEOMETRIC_KERNEL {

	EXACT_PREDICATES_EXACT_CONSTRUCT,
	EXACT_PREDICATES_INEXACT_CONSTRUCT,
	INEXACT_PREDICATES_EXACT_CONSTRUCT,
	INEXACT_PREDICATES_INEXACT_CONSTRUCT

};


	/*****************************************************************************/
	/*                                                                           */
	/*    - Check if two vertices have 'almost_equal' coordinates                */
	/*		  : if yes, one of these vertices is not inserted into triangulation */	
	/*		  : it may cause problems when inserting constraints ?				 */
	/*				- have to think about that more								 */
	/*																			 */
	/*****************************************************************************/
enum class CHECK_KERNEL { 
	
	CHECK_NONE, 
	CHECK_DUPLICATE 

};
