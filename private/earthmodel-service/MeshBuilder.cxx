#include <math.h>
#include <vector>
#include <array>
#include <memory>
#include <tuple>
#include <algorithm>

#include "earthmodel-service/MeshBuilder.h"

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace earthmodel {
namespace Mesh {

#define INSIDE 0
#define OUTSIDE 1

void Voxel::AddPoint(Point const & point) {
	if(n_points == 0) {
		for(unsigned int i=0; i<3; ++i) {
			min_extent[i] = point[i];
			max_extent[i] = point[i];
		}
	} else {
		for(unsigned int i=0; i<3; ++i) {
			min_extent[i] = std::min(min_extent[i], point[i]);
			max_extent[i] = std::max(max_extent[i], point[i]);
		}
	}
	n_points += 1;
}

bool Voxel::Contains(Voxel const & box) const {
	return  min_extent[0] <= box.min_extent[0] and max_extent[0] >= box.max_extent[0]
		and min_extent[1] <= box.min_extent[1] and max_extent[1] >= box.max_extent[1]
		and min_extent[2] <= box.min_extent[2] and max_extent[2] >= box.max_extent[2];
}

double Voxel::SurfaceArea() const {
	Point length;
	for(unsigned int i=0; i<3; ++i) {
		length[i] = std::abs(max_extent[i] - min_extent[i]);
	}
	return 2.0 * (length[0] * (length[1] + length[2]) + length[1] * length[2]);
}

void Voxel::Split(AxisAlignedPlane const & plane, Voxel & VL, Voxel & VR) const {
	int axis = int(plane.axis);
	VL = (*this);
	VL.depth += 1;
	VR = VL;
	VL.max_extent[axis] = plane.position;
	VR.min_extent[axis] = plane.position;
}

double Voxel::EmptyVoxelBias(int NTrianglesLeft, int NTrianglesRight) {
	if(NTrianglesLeft == 0 or NTrianglesRight == 0) {
		return 0.8;
	} else {
		return 1.0;
	}
}

double Voxel::VoxelSAHSplitCost(double probability_left, double probability_right, int num_left, int num_right, double traversal_cost, double intersection_cost) {
	return Voxel::EmptyVoxelBias(num_left, num_right) * (
			traversal_cost + intersection_cost
			* (probability_left * num_left + probability_right * num_right)
			);
}

std::tuple<double, PlanarEventSide> Voxel::VoxelSAHSplitCost(AxisAlignedPlane const & p, int num_left, int num_right, int num_planar, double traversal_cost, double intersection_cost) const {
	Voxel VL, VR;
	Split(p, VL, VR); // Make two voxels from one
	double SA_tot = this->SurfaceArea();
	double SA_L = VL.SurfaceArea();
	double SA_R = VR.SurfaceArea();
	double PL = SA_L / SA_tot; // Probability of hitting voxel is proportional to surface area
	double PR = SA_R / SA_tot;

	// Check what is more expensive: planar events on left side, or planar events on right side
	double CpL = Voxel::VoxelSAHSplitCost(PL, PR, num_left + num_planar, num_right, traversal_cost, intersection_cost);
	double CpR = Voxel::VoxelSAHSplitCost(PL, PR, num_left, num_right + num_planar, traversal_cost, intersection_cost);
	if(CpL < CpR) {
		return std::tuple<double, PlanarEventSide>(CpL, PlanarEventSide::LEFT);
	} else {
		return std::tuple<double, PlanarEventSide>(CpR, PlanarEventSide::RIGHT);
	}
}

std::tuple<AxisAlignedPlane, PlanarEventSide, double> Voxel::FindSplitPlane(int Num_tris, std::vector<Event> const & events, double traversal_cost, double intersection_cost) const {
	// Find the lowest cost plane across all 3 dimensions
	// Relevant plane positions are when a triangle swaps sides from one voxel to another
	// Assume events vector is already sorted
	std::array<int, 3> NL = {0, 0, 0}; // Number of triangles left of the plane
	std::array<int, 3> NP = {0, 0, 0}; // Number of triangles in the plane
	std::array<int, 3> NR = {Num_tris, Num_tris, Num_tris}; // Number of triangles right of the plane
	AxisAlignedPlane best_plane;
	PlanarEventSide best_side;
	double best_cost = 0;
	bool have_cost = false;
	// Sweep across the possible planes, keeping track of how many triangles are in each voxel
	// Loop does 3 dimensions simultaneously ==> Triangle counts are sepearate for each dimension
	for(unsigned int i=0; i<events.size(); ++i) {
		AxisAlignedPlane plane;
		plane.axis = events[i].axis;
		plane.position = events[i].position;
		int dim = int(plane.axis);
		int p_start = 0;
		int p_end = 0;
		int p_planar = 0;
		while(i < events.size() and events[i].axis == plane.axis and events[i].position == plane.position and events[i].type == EventType::END) {
			++p_end; ++i;
		}
		while(i < events.size() and events[i].axis == plane.axis and events[i].position == plane.position and events[i].type == EventType::PLANAR) {
			++p_planar; ++i;
		}
		while(i < events.size() and events[i].axis == plane.axis and events[i].position == plane.position and events[i].type == EventType::START) {
			++p_start; ++i;
		}

		NP[dim] = p_planar;
		NR[dim] -= p_planar + p_end;

		std::tuple<double, PlanarEventSide> cost_side =
			this->VoxelSAHSplitCost(plane, NL[dim], NR[dim], NP[dim], traversal_cost, intersection_cost);
		if((not have_cost) or std::get<0>(cost_side) < best_cost) {
			have_cost = true;
			best_cost = std::get<0>(cost_side);
			best_side = std::get<1>(cost_side);
			best_plane = plane;
		}

		NL[dim] += p_planar + p_start;
		NP[dim] = 0;
	}
	return std::tuple<AxisAlignedPlane, PlanarEventSide, double>(best_plane, best_side, best_cost);
}

bool Voxel::Intersects(TData const & triangle) const {
	Point lengths = subtract(this->max_extent, this->min_extent);

	TData scaled_tri = triangle;

	for(unsigned int dim=0; dim<3; ++dim) {
	lengths[dim] = std::abs(lengths[dim]);
		for(unsigned int p_i=0; p_i<3; ++p_i) {
			scaled_tri[p_i][dim] -= this->min_extent[dim];
			scaled_tri[p_i][dim] /= lengths[dim];
		}
	}

	return t_c_intersection(scaled_tri) == INSIDE;
}

bool Voxel::Intersects(Voxel const & voxel) const {
    bool status = true;

    // AABBs cannot intersect if they are separated along any dimension
    for(int i=0; i<3; ++i) {
        double min1 = this->min_extent[i];
        double max1 = this->max_extent[i];
        double min2 = voxel.min_extent[i];
        double max2 = voxel.max_extent[i];
        status &= not (max1 < min2 || min1 > max2);
    }  // END for all dimensions

    return status;
}

PolygonData Voxel::Clip(TData const & tri) const {

	// Use two polygons with pointers for 'back-buffer'-like swapping
	const int MAX_VERTS = 6;
	PolygonData poly[2] = {PolygonData(MAX_VERTS), PolygonData(MAX_VERTS)};
	PolygonData* currentPoly = &poly[0];
	PolygonData* prevPoly = &poly[1];

	// First check if the triangle is contained in the bbox, if not we are empty
	Voxel triBox;
	triBox.AddPoint(tri[0]);
	triBox.AddPoint(tri[1]);
	triBox.AddPoint(tri[2]);

	if(not this->Intersects(triBox)) {
		return *currentPoly;  // note: currentPoly is empty
	}

	// Set up the initial polygon
	currentPoly->push_back(tri[0]);
	currentPoly->push_back(tri[1]);
	currentPoly->push_back(tri[2]);

	// If all the vertices are contained, we have no work to do
	if(this->Contains(triBox)) {
		return *currentPoly;  // Note current poly has the same verts as tri
	}

	// Loop through the planes of the bbox and clip the vertices
	for(int dim = 0; dim < 3; ++dim) {
		// Optimization note: we should be able to save some work based on
		// the clipping plane and the triangle's bounding box

		if(triBox.max_extent[dim] > this->min_extent[dim]) {
            Mesh::swap(prevPoly, currentPoly);
			clipAxisPlane(prevPoly, currentPoly, 2 * dim + 0, this->min_extent[dim]);
		}
		if(triBox.min_extent[dim] < this->max_extent[dim]) {
            Mesh::swap(prevPoly, currentPoly);
			clipAxisPlane(prevPoly, currentPoly, 2 * dim + 1, this->max_extent[dim]);
		}
	}

	return *currentPoly;
}



///////////////
// Fast Box-Triangle intersection check
///////////////

#define EPS 10e-5
#define SIGN3( A ) \
  (((A)[0] < EPS) ? 4 : 0 | ((A)[0] > -EPS) ? 32 : 0 | \
   ((A)[1] < EPS) ? 2 : 0 | ((A)[1] > -EPS) ? 16 : 0 | \
   ((A)[2] < EPS) ? 1 : 0 | ((A)[2] > -EPS) ? 8 : 0)

#define CROSS( A, B, C ) { \
  (C)[0] =  (A)[1] * (B)[2] - (A)[2] * (B)[1]; \
  (C)[1] = -(A)[0] * (B)[2] + (A)[2] * (B)[0]; \
  (C)[2] =  (A)[0] * (B)[1] - (A)[1] * (B)[0]; \
   }
#define SUB( A, B, C ) { \
  (C)[0] =  (A)[0] - (B)[0]; \
  (C)[1] =  (A)[1] - (B)[1]; \
  (C)[2] =  (A)[2] - (B)[2]; \
   }
#define LERP( A, B, C) ((B)+(A)*((C)-(B)))
#define MIN3(a,b,c) ((((a)<(b))&&((a)<(c))) ? (a) : (((b)<(c)) ? (b) : (c)))
#define MAX3(a,b,c) ((((a)>(b))&&((a)>(c))) ? (a) : (((b)>(c)) ? (b) : (c)))

/*___________________________________________________________________________*/

/* Which of the six face-plane(s) is point P outside of? */

long face_plane(Point p) {
	long outcode;
	outcode = 0;
	if (p[0] >  .5) outcode |= 0x01;
	if (p[0] < -.5) outcode |= 0x02;
	if (p[1] >  .5) outcode |= 0x04;
	if (p[1] < -.5) outcode |= 0x08;
	if (p[2] >  .5) outcode |= 0x10;
	if (p[2] < -.5) outcode |= 0x20;
	return(outcode);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Which of the twelve edge plane(s) is point P outside of? */

long bevel_2d(Point p) {
	long outcode;

	outcode = 0;
	if ( p[0] + p[1] > 1.0) outcode |= 0x001;
	if ( p[0] - p[1] > 1.0) outcode |= 0x002;
	if (-p[0] + p[1] > 1.0) outcode |= 0x004;
	if (-p[0] - p[1] > 1.0) outcode |= 0x008;
	if ( p[0] + p[2] > 1.0) outcode |= 0x010;
	if ( p[0] - p[2] > 1.0) outcode |= 0x020;
	if (-p[0] + p[2] > 1.0) outcode |= 0x040;
	if (-p[0] - p[2] > 1.0) outcode |= 0x080;
	if ( p[1] + p[2] > 1.0) outcode |= 0x100;
	if ( p[1] - p[2] > 1.0) outcode |= 0x200;
	if (-p[1] + p[2] > 1.0) outcode |= 0x400;
	if (-p[1] - p[2] > 1.0) outcode |= 0x800;
	return(outcode);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Which of the eight corner plane(s) is point P outside of? */

long bevel_3d(Point p) {
	long outcode;

   outcode = 0;
   if (( p[0] + p[1] + p[2]) > 1.5) outcode |= 0x01;
   if (( p[0] + p[1] - p[2]) > 1.5) outcode |= 0x02;
   if (( p[0] - p[1] + p[2]) > 1.5) outcode |= 0x04;
   if (( p[0] - p[1] - p[2]) > 1.5) outcode |= 0x08;
   if ((-p[0] + p[1] + p[2]) > 1.5) outcode |= 0x10;
   if ((-p[0] + p[1] - p[2]) > 1.5) outcode |= 0x20;
   if ((-p[0] - p[1] + p[2]) > 1.5) outcode |= 0x40;
   if ((-p[0] - p[1] - p[2]) > 1.5) outcode |= 0x80;
   return(outcode);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Test the point "alpha" of the way from P1 to P2 */
/* See if it is on a face of the cube              */
/* Consider only faces in "mask"                   */

long check_point(Point p1, Point p2, float alpha, long mask) {
	Point plane_point;

	plane_point[0] = LERP(alpha, p1[0], p2[0]);
	plane_point[1] = LERP(alpha, p1[1], p2[1]);
	plane_point[2] = LERP(alpha, p1[2], p2[2]);
	return(face_plane(plane_point) & mask);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Compute intersection of P1 --> P2 line segment with face planes */
/* Then test intersection point to see if it is on cube face       */
/* Consider only face planes in "outcode_diff"                     */
/* Note: Zero bits in "outcode_diff" means face line is outside of */

long check_line(Point p1, Point p2, long outcode_diff) {
	if ((0x01 & outcode_diff) != 0)
		if (check_point(p1,p2,( 0.5f-p1[0])/(p2[0]-p1[0]),0x3e) == INSIDE)
			return(INSIDE);
	if ((0x02 & outcode_diff) != 0)
		if (check_point(p1,p2,(-0.5f-p1[0])/(p2[0]-p1[0]),0x3d) == INSIDE)
            return(INSIDE);
	if ((0x04 & outcode_diff) != 0)
		if (check_point(p1,p2,( 0.5f-p1[1])/(p2[1]-p1[1]),0x3b) == INSIDE)
            return(INSIDE);
	if ((0x08 & outcode_diff) != 0)
		if (check_point(p1,p2,(-0.5f-p1[1])/(p2[1]-p1[1]),0x37) == INSIDE)
            return(INSIDE);
	if ((0x10 & outcode_diff) != 0)
		if (check_point(p1,p2,( 0.5f-p1[2])/(p2[2]-p1[2]),0x2f) == INSIDE)
            return(INSIDE);
	if ((0x20 & outcode_diff) != 0)
		if (check_point(p1,p2,(-0.5f-p1[2])/(p2[2]-p1[2]),0x1f) == INSIDE)
            return(INSIDE);
	return(OUTSIDE);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Test if 3D point is inside 3D triangle */

long point_triangle_intersection(Point p, TData t) {
	long sign12,sign23,sign31;
	Point vect12,vect23,vect31,vect1h,vect2h,vect3h;
	Point cross12_1p,cross23_2p,cross31_3p;

	/* First, a quick bounding-box test:                               */
	/* If P is outside triangle bbox, there cannot be an intersection. */

	if (p[0] > MAX3(t[0][0], t[1][0], t[2][0])) return(OUTSIDE);
	if (p[1] > MAX3(t[0][1], t[1][1], t[2][1])) return(OUTSIDE);
	if (p[2] > MAX3(t[0][2], t[1][2], t[2][2])) return(OUTSIDE);
	if (p[0] < MIN3(t[0][0], t[1][0], t[2][0])) return(OUTSIDE);
	if (p[1] < MIN3(t[0][1], t[1][1], t[2][1])) return(OUTSIDE);
	if (p[2] < MIN3(t[0][2], t[1][2], t[2][2])) return(OUTSIDE);

	/* For each triangle side, make a vector out of it by subtracting vertexes; */
	/* make another vector from one vertex to point P.                          */
	/* The crossproduct of these two vectors is orthogonal to both and the      */
	/* signs of its X,Y,Z components indicate whether P was to the inside or    */
	/* to the outside of this triangle side.                                    */

	SUB(t[0], t[1], vect12)
		SUB(t[0],    p, vect1h);
	CROSS(vect12, vect1h, cross12_1p)
		sign12 = SIGN3(cross12_1p);      /* Extract X,Y,Z signs as 0..7 or 0...63 integer */

	SUB(t[1], t[2], vect23)
		SUB(t[1],    p, vect2h);
	CROSS(vect23, vect2h, cross23_2p)
		sign23 = SIGN3(cross23_2p);

	SUB(t[2], t[0], vect31)
		SUB(t[2],    p, vect3h);
	CROSS(vect31, vect3h, cross31_3p)
		sign31 = SIGN3(cross31_3p);

	/* If all three crossproduct vectors agree in their component signs,  */
	/* then the point must be inside all three.                           */
	/* P cannot be OUTSIDE all three sides simultaneously.                */

	/* this is the old test; with the revised SIGN3() macro, the test
	 * needs to be revised. */
	return ((sign12 & sign23 & sign31) == 0) ? OUTSIDE : INSIDE;
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/**********************************************/
/* This is the main algorithm procedure.      */
/* Triangle t is compared with a unit cube,   */
/* centered on the origin.                    */
/* It returns INSIDE (0) or OUTSIDE(1) if t   */
/* intersects or does not intersect the cube. */
/**********************************************/

IntersectionResut t_c_intersection(TData t) {
    long v1_test,v2_test,v3_test;
    float d,denom;
    Point vect12,vect13,norm;
    Point hitpp,hitpn,hitnp,hitnn;

    /* First compare all three vertexes with all six face-planes */
    /* If any vertex is inside the cube, return immediately!     */

    if ((v1_test = face_plane(t[0])) == INSIDE) return(INSIDE);
    if ((v2_test = face_plane(t[1])) == INSIDE) return(INSIDE);
    if ((v3_test = face_plane(t[2])) == INSIDE) return(INSIDE);

    /* If all three vertexes were outside of one or more face-planes, */
    /* return immediately with a trivial rejection!                   */

    if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);

    /* Now do the same trivial rejection test for the 12 edge planes */

    v1_test |= bevel_2d(t[0]) << 8;
    v2_test |= bevel_2d(t[1]) << 8;
    v3_test |= bevel_2d(t[2]) << 8;
    if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);

    /* Now do the same trivial rejection test for the 8 corner planes */

    v1_test |= bevel_3d(t[0]) << 24;
    v2_test |= bevel_3d(t[1]) << 24;
    v3_test |= bevel_3d(t[2]) << 24;
    if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);

    /* If vertex 1 and 2, as a pair, cannot be trivially rejected */
    /* by the above tests, then see if the v1-->v2 triangle edge  */
    /* intersects the cube.  Do the same for v1-->v3 and v2-->v3. */
    /* Pass to the intersection algorithm the "OR" of the outcode */
    /* bits, so that only those cube faces which are spanned by   */
    /* each triangle edge need be tested.                         */

    if ((v1_test & v2_test) == 0)
        if (check_line(t[0],t[1],v1_test|v2_test) == INSIDE)
            return(INSIDE);
    if ((v1_test & v3_test) == 0)
        if (check_line(t[0],t[2],v1_test|v3_test) == INSIDE)
            return(INSIDE);
    if ((v2_test & v3_test) == 0)
        if (check_line(t[1],t[2],v2_test|v3_test) == INSIDE)
            return(INSIDE);

    /* By now, we know that the triangle is not off to any side,     */
    /* and that its sides do not penetrate the cube.  We must now    */
    /* test for the cube intersecting the interior of the triangle.  */
    /* We do this by looking for intersections between the cube      */
    /* diagonals and the triangle...first finding the intersection   */
    /* of the four diagonals with the plane of the triangle, and     */
    /* then if that intersection is inside the cube, pursuing        */
    /* whether the intersection point is inside the triangle itself. */

    /* To find plane of the triangle, first perform crossproduct on  */
    /* two triangle side vectors to compute the normal vector.       */

    SUB(t[0],t[1],vect12);
    SUB(t[0],t[2],vect13);
    CROSS(vect12,vect13,norm);

    /* The normal vector "Vnorm" X,Y,Z components are the coefficients */
    /* of the triangles AX + BY + CZ + D = 0 plane equation.  If we   */
    /* solve the plane equation for X=Y=Z (a diagonal), we get        */
    /* -D/(A+B+C) as a metric of the distance from cube center to the */
    /* diagonal/plane intersection.  If this is between -0.5 and 0.5, */
    /* the intersection is inside the cube.  If so, we continue by    */
    /* doing a point/triangle intersection.                           */
    /* Do this for all four diagonals.                                */

    d = norm[0] * t[0][0] + norm[1] * t[0][1] + norm[2] * t[0][2];

    /* if one of the diagonals is parallel to the plane, the other will intersect the plane */
    if(fabs(denom=(norm[0] + norm[1] + norm[2]))>EPS) {
        /* skip parallel diagonals to the plane; division by 0 can occur */
        hitpp[0] = hitpp[1] = hitpp[2] = d / denom;
        if (fabs(hitpp[0]) <= 0.5)
            if (point_triangle_intersection(hitpp,t) == INSIDE)
                return(INSIDE);
    }
    if(fabs(denom=(norm[0] + norm[1] - norm[2]))>EPS) {
        hitpn[2] = -(hitpn[0] = hitpn[1] = d / denom);
        if (fabs(hitpn[0]) <= 0.5)
            if (point_triangle_intersection(hitpn,t) == INSIDE)
                return(INSIDE);
    }
    if(fabs(denom=(norm[0] - norm[1] + norm[2]))>EPS) {
        hitnp[1] = -(hitnp[0] = hitnp[2] = d / denom);
        if (fabs(hitnp[0]) <= 0.5)
            if (point_triangle_intersection(hitnp,t) == INSIDE)
                return(INSIDE);
    }
    if(fabs(denom=(norm[0] - norm[1] - norm[2]))>EPS) {
        hitnn[1] = hitnn[2] = -(hitnn[0] = d / denom);
        if (fabs(hitnn[0]) <= 0.5)
            if (point_triangle_intersection(hitnn,t) == INSIDE)
                return(INSIDE);
    }

    /* No edge touched the cube; no cube diagonal touched the triangle. */
    /* We're done...there was no intersection.                          */

    return(OUTSIDE);
}

#undef INSIDE
#undef OUTSIDE
#undef MAX3
#undef MIN3
#undef LERP
#undef SUB
#undef CROSS
#undef SIGN3
#undef EPS

///////////////
// Clip triangle to voxel
///////////////

double dot(Point const & a, Point const & b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Point mul(Point const & a, double x) {
    Point ret = a;
    ret[0] *= x;
    ret[1] *= x;
    ret[2] *= x;
    return ret;
}

Point add(Point const & a, Point const & b) {
    Point ret = a;
    ret[0] += b[0];
    ret[1] += b[1];
    ret[2] += b[2];
    return ret;
}

Point subtract(Point const & a, Point const & b) {
    Point ret = a;
    ret[0] -= b[0];
    ret[1] -= b[1];
    ret[2] -= b[2];
    return ret;
}


bool isEven(int i) {
	return i % 2 == 0;
}

OrientationResult classifyPointAxisPlane(Point const & pt,
		int index,
		double val,
		const double eps) {
	// Note: we are exploiting the fact that the planes are axis aligned
	// So the dot product is +/- the given coordinate.
	// In general, we would need to call distance(pt, plane) here
	double dist = isEven(index) ? val - pt[index / 2] : pt[index / 2] - val;

	if(dist > eps) {
		return OrientationResult::ON_POSITIVE_SIDE;
	}
	if(dist < -eps) {
		return OrientationResult::ON_NEGATIVE_SIDE;
	}

	return OrientationResult::ON_BOUNDARY;
}

Point findIntersectionPoint(Point const & a,
		Point const & b,
		int index,
		double val) {
	// Need to find a parameter t for the point pt, such that,
	// * 0 <= t <= 1
	// * pt = a + t * (b-a)
	// * pt[ index/2]  == val
	double t = (val - a[index / 2]) / (b[index / 2] - a[index / 2]);
	Point ret = add(a, mul(subtract(b, a), t));
	return ret;
}

void clipAxisPlane(PolygonData const * prevPoly,
		PolygonData * currentPoly,
		int index,
		double val)
{

	currentPoly->clear();
	int numVerts = prevPoly->size();
	if(numVerts == 0) {
		return;
	}
	// Initialize point a with the last vertex of the polygon
	const Point* a = &(*prevPoly)[numVerts - 1];
	int aSide = classifyPointAxisPlane(*a, index, val);
	for(int i = 0; i < numVerts; ++i) {
		const Point* b = &(*prevPoly)[i];
		int bSide = classifyPointAxisPlane(*b, index, val);
		switch(bSide) {
            case OrientationResult::ON_POSITIVE_SIDE:
				if(aSide == OrientationResult::ON_NEGATIVE_SIDE) {
					currentPoly->push_back(findIntersectionPoint(*a, *b, index, val));
				}
				break;
			case ON_BOUNDARY:
				if(aSide == OrientationResult::ON_NEGATIVE_SIDE) {
					currentPoly->push_back(*b);
				}
				break;
			case OrientationResult::ON_NEGATIVE_SIDE:
				switch(aSide) {
                    case OrientationResult::ON_POSITIVE_SIDE:
						currentPoly->push_back(findIntersectionPoint(*a, *b, index, val));
						currentPoly->push_back(*b);
						break;
					case ON_BOUNDARY:
						currentPoly->push_back(*a);
						currentPoly->push_back(*b);
						break;
					case OrientationResult::ON_NEGATIVE_SIDE:
						currentPoly->push_back(*b);
						break;
				}
				break;
		}
		// swap a and b
		a = b;
		aSide = bSide;
	}
}

///////////////
// Split events
///////////////

void ClassifyEventLeftRightBoth(std::vector<Event> & E, AxisAlignedPlane const & p, PlanarEventSide side) {
    for(unsigned int i=0; i<E.size(); ++i) {
        Event & e = E[i];
        E[i].side = EventPlaneSide::BOTH;
    }
    EventPlaneSide default_planar_side;
    if(side == PlanarEventSide::LEFT) {
        default_planar_side = EventPlaneSide::LEFT;
    } else {
        default_planar_side = EventPlaneSide::RIGHT;
    }
    for(unsigned int i=0; i<E.size(); ++i) {
        Event & e = E[i];
        if(e.type == EventType::END
                and e.axis == p.axis
                and e.position <= p.position) {
            e.side = EventPlaneSide::LEFT;
        } else if (e.type == EventType::START
                and e.axis == p.axis
                and e.position >= p.position) {
            e.side = EventPlaneSide::RIGHT
        } else if (e.type == EventType::PLANAR
                and e.axis == p.axis) {
            if(e.position < p.position) {
                e.side = EventPlaneSide::LEFT;
            } else if(e.position > p.position) {
                e.side = EventPlaneSide::RIGHT;
            } else if(e.position == p.position) {
                e.side = default_planar_side;
            }
        }
    }
}

void AddPlanarEvent(std::vector<Event> & events, Voxel const & tri_box, Axis axis, Triangle triangle) {
    Event e;
    int dim = int(axis);
    e.position = tri_box.min_extent[dim];
    e.axis = axis;
    e.type = EventType::PLANAR;
    e.triangle = triangle;
    e.side = EventPlaneSide::BOTH;
    events.push_back(e);
}

void AddStartEndEvents(std::vector<Event> & events, Voxel const & tri_box, Axis axis, Triangle triangle) {
    Event e;
    int dim = int(axis);
    e.position = tri_box.min_extent[dim];
    e.axis = axis;
    e.type = EventType::START;
    e.triangle = triangle;
    e.side = EventPlaneSide::BOTH;
    events.push_back(e);
    Event e_end;
    e.position = tri_box.max_extent[dim];
    e.type = EventType::END;
    events.push_back(e);
}

void GenerateNonClippedTriangleVoxelEvents(std::vector<Event> & events, TData const & tri_data, Triangle triangle) {
    Voxel box;
    for(unsigned int point_i=0; point_i<3; ++point_i) {
        box.AddPoint(tri_data[point_i]);
    }

    for(unsigned int dim_i=0; dim_i<3; ++dim_i) {
        if(box.min_extent[dim_i] == box.max_extent[dim_i]) {
            // Planar event
            AddPlanarEvent(events, box, Axis(dim_i), triangle);
        } else {
            // Start and End events
            AddStartEndEvents(events, box, Axis(dim_i), triangle);
        }
    }
}

void GenerateClippedTriangleVoxelEvents(std::vector<Event> & events, TData const & tri_data, Triangle triangle, Voxel const & voxel_box) {
    PolygonData p = voxel_box.Clip(tri_data);

    Voxel box;
    for(unsigned int point_i=0; point_i<p.size(); ++point_i) {
        box.AddPoint(p[point_i]);
    }

    for(unsigned int dim_i=0; dim_i<3; ++dim_i) {
        if(box.min_extent[dim_i] == box.max_extent[dim_i]) {
            // Planar event
            AddPlanarEvent(events, box, Axis(dim_i), triangle);
        } else {
            // Start and End events
            AddStartEndEvents(events, box, Axis(dim_i), triangle);
        }
    }
}

void GeneratePlaneEvents(std::vector<Event> & events_L, std::vector<Event> & events_R, std::vector<TData> const & triangle_data, std::vector<Triangle> const & intersecting_tris, Voxel const & voxel, AxisAlignedPlane const & plane) {
    Voxel VL;
    Voxel VR;
    voxel.SplitBox(plane, VL, VR);

    for(unsigned int i=0; i<triangles.size()) {
        GenerateClippedTriangleVoxelEvents(events_L, triangle_data[intersecting_tris[i]], intersecting_tris[i], VL);
        GenerateClippedTriangleVoxelEvents(events_R, triangle_data[intersecting_tris[i]], intersecting_tris[i], VR);
    }
}

int TauEventType(EventType etype) {
    return int(etype);
}


bool EventCompare(Event const & a, Event const & b) {
    // Is a < b
    return (a.position < b.position)
        or (a.position == b.position and TauEventType(a) < TauEventType(b));
}

void SplitEventsByPlane(std::vector<Event> const & events,
        std::vector<TData> const & triangle_data,
        Voxel const & voxel,
        AxisAlignedPlane const & plane,
        std::vector<Event> & EL,
        std::vector<Event> & ER,
        std::vector<Triangle> & TL,
        std::vector<Triangle> & TR,
        PlanarEventSide const & side) {
    std::vector<Event> ELtemp;
    std::vector<Event> ERtemp;
    std::vector<Event> EBLtemp;
    std::vector<Event> EBRtemp;
    ClassifyEventLeftRightBoth(events, plane, side);
    std::vector<Triangle> intersecting_tris;
    for(unsigned int i=0; i<events.size(); ++i) {
        Event & e = events[i];
        if(e.side == EventPlaneSide::BOTH) {
            intersecting_tris.push_back(e.triangle);
        } else if (e.side == EventPlaneSide::LEFT) {
            ELtemp.push_back(e);
        } else if (e.side == EventPlaneSide::RIGHT) {
            ERtemp.push_back(e);
        }
    }
    GeneratePlaneEvents(EBLtemp, EBRtemp, triangle_data, intersecting_tris, voxel, plane);
    std::sort(EBLtemp.begin(), EBLtemp.end(), EventCompare);
    std::sort(EBRtemp.begin(), EBRtemp.end(), EventCompare);
    std::merge(ELtemp.begin(), ELtemp.end(), EBLtemp.begin(), EBLtemp.end(), EL.begin(), EventCompare);
    std::merge(ERtemp.begin(), ERtemp.end(), EBRtemp.begin(), EBRtemp.end(), ER.begin(), EventCompare);
	for(unsigned int i=0; ++i; i<EL.size()) {
        Event const & e = EL[i];
        if(e.axis != plane.axis)
            continue;
        TL.append(e.triangle);
    }
	for(unsigned int i=0; ++i; i<ER.size()) {
        Event const & e = ER[i];
        if(e.axis != plane.axis)
            continue;
        TR.append(e.triangle);
    }
}

///////////////
// Build KD Tree
///////////////

std::shared_ptr<KDNode> RecBuild(std::vector<TData> const & triangle_data, std::vector<Triangle> const & T, Voxel & V, std::vector<Event> const & events, double traversal_cost, double intersection_cost, int max_depth) {
    std::tuple<AxisAlignedPlane, PlanarEventSide, double> plane_side_cost = V.FindSplitPlane(T.size(), events, traversal_cost, intersection_cost);
    double cost = std::get<2>(plane_side_cost);
    double termination_cost = intersection_cost * T.size();

    if(cost > termination_cost or V.depth >= max_depth) {
        return std::make_shared<KDNode>(V, T);
    }

    AxisAlignedPlane const & plane = std::get<0>(plane_side_cost);
    PlanarEventSide const & side = std::get<1>(plane_side_cost);

    std::vector<Event> EL, ER;
    std::vector<Triangle> TL, TR;
    SplitEventsByPlane(events, triangle_data, V, plane, EL, ER, TL, TR, side);
    Voxel VL, VR;
    V.SplitBox(plane, VL, VR);
    return std::make_shared<KDNode>(V, RecBuild(triangle_data, TL, VL, EL, traversal_cost, intersection_cost, max_depth), RecBuild(triangle_data, TR, VR, ER, traversal_cost, intersection_cost, max_depth));
}

std::shared_ptr<KDNode> BuildKDTree(std::vector<TData> const & triangle_data, double traversal_cost, double intersection_cost, int max_depth) {
    // Generate all events and compute world bounding box
    Voxel V;
    std::vector<Event> events;
    for(unsigned int i=0; i<triangles.size(); ++i) {
        GenerateNonClippedTriangleVoxelEvents(events, triangles[i], i);
        box.addPoint(triangles[i][0]);
        box.addPoint(triangles[i][1]);
        box.addPoint(triangles[i][2]);
    }
    V.depth = 0;

    // Sort events once
    std::sort(std::begin(events), std::end(events), EventCompare);

    // Start with all triangle IDs
    std::vector<int> T(triangles.size());
    std::iota(std::begin(T), std::end(T), 0); // Fill T with 0, 1, ..., triangles.size()

    // Recursively build the tree
    return RecBuild(triangles, T, V, events, traversal_cost, intersection_cost);
}

} // namespace Mesh
} // namespace earthmodel
