#pragma once
#ifndef SIREN_MeshBuilder_H
#define SIREN_MeshBuilder_H

#include <set>
#include <map>
#include <array>
#include <tuple>
#include <memory>
#include <vector>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/serialization/array.h"

namespace siren {
namespace geometry {
namespace Mesh {

typedef int                   Vertex;
typedef int                   VertexID;
typedef int                   EdgeID;
typedef int                   TriangleID;
typedef std::array<Vertex, 2> Edge;
typedef std::array<Vertex, 3> Triangle;
typedef std::array<double, 3> VData;
typedef std::array<VData,  2> EData;
typedef std::array<VData,  3> TData;
typedef VData Point;

typedef std::vector<Point> PolygonData;

struct VAttribute {
    VData data;
    std::set<Edge> eset;
    std::set<Triangle> tset;
    bool operator==(VAttribute const & other) const;
    bool operator!=(VAttribute const & other) const;
    bool operator<(VAttribute const & other) const;
};

struct EAttribute {
    EData data;
    std::set<Triangle> tset;
    bool operator==(EAttribute const & other) const;
    bool operator!=(EAttribute const & other) const;
    bool operator<(EAttribute const & other) const;
};

struct TAttribute {
    TData data;
    bool operator==(TAttribute const & other) const;
    bool operator!=(TAttribute const & other) const;
    bool operator<(TAttribute const & other) const;
};

typedef std::pair<Vertex,   VAttribute> VPair;
typedef std::pair<Edge,     EAttribute> EPair;
typedef std::pair<Triangle, TAttribute> TPair;
typedef std::vector<VAttribute>        VMap;
typedef std::map<Edge,     EAttribute> EMap;
typedef std::map<Triangle, TAttribute> TMap;


struct TMesh {
    VMap vmap;
    EMap emap;
    TMap tmap;
    bool operator==(TMesh const & other) const;
    bool operator!=(TMesh const & other) const;
    bool operator<(TMesh const & other) const;
};

enum class Axis {
    X=0, Y=1, Z=2
};

enum class PlanarEventSide {
    LEFT,
    RIGHT
};

enum class EventPlaneSide {
    LEFT,
    RIGHT,
    BOTH
};

enum class EventType {
    END = 0,
    PLANAR = 1,
    START = 2
};

struct AxisAlignedPlane {
    Axis axis;
    double position;
};

struct Event {
    Axis axis;
    double position;
    EventType type;
    TriangleID triangle;
};

struct Voxel {
    int depth = 0;
    int n_points = 0;
    Point min_extent;
    Point max_extent;
    void AddPoint(Point const & point);
    bool Contains(Voxel const & box) const;
	double SurfaceArea() const;
    void Split(AxisAlignedPlane const & plane, Voxel & VL, Voxel & VR) const;
    static double EmptyVoxelBias(int NTrianglesLeft, int NTrianglesRight);
    static double VoxelSAHSplitCost(double probability_left, double probability_right, int num_left, int num_right, double traversal_cost, double intersection_cost);
    std::tuple<double, PlanarEventSide> VoxelSAHSplitCost(AxisAlignedPlane const & p, int num_left, int num_right, int num_planar, double traversal_cost, double intersection_cost) const;
    std::tuple<AxisAlignedPlane, PlanarEventSide, double> FindSplitPlane(int Num_tris, std::vector<Event> const & events, double traversal_cost, double intersection_cost) const;
    bool Intersects(TData const & triangle) const;
    bool Intersects(Voxel const & Voxel) const;
    PolygonData Clip(TData const & triangle) const;
};

long face_plane(Point p);
long bevel_2d(Point p);
long bevel_3d(Point p);
long check_point(Point p1, Point p2, float alpha, long mask);
long check_line(Point p1, Point p2, long outcode_diff);
long point_triangle_intersection(Point p, TData t);
long t_c_intersection(TData t);

enum class OrientationResult {
    ON_BOUNDARY,      /*!< primitive is on the boundary of a primitive      */
    ON_POSITIVE_SIDE, /*!< primitive is on the positive side of a primitive */
    ON_NEGATIVE_SIDE  /*!< primitive is on the negative side of a primitive */
};

double dot(Point const & a, Point const & b);
Point mul(Point const & a, double x);
Point add(Point const & a, Point const & b);
Point subtract(Point const & a, Point const & b);
bool isEvent(int i);
OrientationResult classifyPointAxisPlane(Point const & pt, int index, double val, const double eps = 1e-8);
Point findIntersectionPoint(Point const & a, Point const & b, int index, double val);
void clipAxisPlane(PolygonData const * prevPoly, PolygonData * currentPoly, int index, double val);

std::vector<EventPlaneSide> ClassifyEventLeftRightBoth(std::vector<Event> & E, AxisAlignedPlane const & p, PlanarEventSide side);
void AddPlanarEvent(std::vector<Event> & events, Voxel const & tri_box, AxisAlignedPlane axis, int tri_id);
void AddStartEndEvents(std::vector<Event> & events, Voxel const & tri_box, AxisAlignedPlane axis, int tri_id);
void GenerateNonClippedTriangleVoxelEvents(std::vector<Event> & events, TData const & tri_data, TriangleID triangle);
void GenerateClippedTriangleVoxelEvents(std::vector<Event> & events, TData const & tri_data, TriangleID triangle, Voxel const & voxel_box);
void GeneratePlaneEvents(std::vector<Event> & events_L, std::vector<Event> & events_R, std::vector<TData> const & triangle_data, std::vector<TriangleID> const & intersecting_tris, Voxel const & voxel, AxisAlignedPlane const & plane);
int TauEventType(EventType etype);
bool EventCompare(Event const & a, Event const & b);
void SplitEventsByPlane(std::vector<Event> const & events, std::vector<TData> const & triangle_data, Voxel const & voxel, AxisAlignedPlane const & plane, std::vector<Event> & EL, std::vector<Event> & ER, PlanarEventSide const & side);


struct KDNode {
    bool terminal = false;
    Voxel V;
    std::vector<TriangleID> T;
    std::shared_ptr<KDNode> left;
    std::shared_ptr<KDNode> right;
    KDNode(Voxel const & voxel, std::vector<TriangleID> const & triangles) :
        terminal(true), V(voxel), T(triangles) {};
    KDNode(Voxel const & voxel, std::shared_ptr<KDNode> l, std::shared_ptr<KDNode> r) :
        terminal(false), V(voxel), left(l), right(r) {};
};

std::shared_ptr<KDNode> RecBuild(std::vector<TData> const & triangle_data, std::vector<TriangleID> const & T, Voxel & V, std::vector<Event> const & events, double traversal_cost, double intersection_cost, int max_depth);
std::shared_ptr<KDNode> BuildKDTree(std::vector<TData> const & triangles_data, double traversal_cost, double intersection_cost, int max_depth);


} // namespace Mesh
} // namespace geometry
} // namespace siren

#endif // SIREN_MeshBuilder_H

