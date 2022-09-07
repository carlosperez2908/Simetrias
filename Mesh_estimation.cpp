#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <fstream>
#include <cassert>
#include <math.h>

#include <CGAL/property_map.h>
#include <CGAL/Random.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/Ridges.h>
#include <CGAL/Umbilics.h>
#include <CGAL/Plane_3.h>
#include <CGAL/intersections.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <boost/iterator/counting_iterator.hpp>

#if defined(CGAL_USE_BOOST_PROGRAM_OPTIONS) && ! defined(DONT_USE_BOOST_PROGRAM_OPTIONS)
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#endif

using namespace std;

#include "PolyhedralSurf.h"
#include "PolyhedralSurf_operations.h"
#include "PolyhedralSurf_rings.h"

//Kernel of the PolyhedralSurf
typedef double                DFT;
typedef CGAL::Simple_cartesian<DFT>  Data_Kernel;
typedef Data_Kernel::Point_3  DPoint;
typedef Data_Kernel::Vector_3 DVector;

//HDS
typedef PolyhedralSurf::Vertex_handle Vertex_handle;
typedef PolyhedralSurf::Vertex Vertex;
typedef PolyhedralSurf::Halfedge_handle Halfedge_handle;
typedef PolyhedralSurf::Halfedge Halfedge;
typedef PolyhedralSurf::Vertex_iterator Vertex_iterator;
typedef PolyhedralSurf::Facet_handle Facet_handle;
typedef PolyhedralSurf::Facet Facet;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Line_3 Line_3;
typedef Kernel::Intersect_3 Intersect_3;

//definition of a non-mutable lvalue property map,
//with the get function as a friend function to give it
//access to the private member
class PointPropertyMap{
	const std::vector<Point_3>& points;
public:
	typedef Point_3 value_type;
	typedef const value_type& reference;
	typedef std::size_t key_type;
	typedef boost::lvalue_property_map_tag category;
	PointPropertyMap(const std::vector<Point_3>& pts) :points(pts){}
	reference operator[](key_type k) const { return points[k]; }
	friend reference get(const PointPropertyMap& ppmap, key_type i)
	{
		return ppmap[i];
	}
};


typedef CGAL::Random_points_in_cube_3<Point_3>									Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>											Traits_base;
typedef CGAL::Search_traits_adapter<std::size_t, PointPropertyMap, Traits_base> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>								K_neighbor_search;
typedef K_neighbor_search::Tree													Tree;
typedef Tree::Splitter															Splitter;
typedef K_neighbor_search::Distance												Distance;

typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;

//typedef CGAL::Linear_algebraCd<DFT> LA;
//typedef LA::Matrix LMatrix; typedef LA::Vector LVector; typedef LA::Matrix::Identity Identity; 
//LMatrix identity = LMatrix(3, Identity());

struct Hedge_cmp{
  bool operator()(Halfedge_handle a,  Halfedge_handle b) const{
    return &*a < &*b;
  }
};

struct Facet_cmp{
  bool operator()(Facet_handle a, Facet_handle b) const{
    return &*a < &*b;
  }
};

const double EPSILON = 0.66;  //0.33;
const double CLUSTER_EPSILON = 0.75;
const double KERNEL_BW = 0.75; //2.25
const string filePath = "30_original.off";
const double MM_CURVATURE_THRESHOLD = 0.66;
const double RT_CURVATURE = 0.95;
std::ofstream out_verbose_rt, out_verbose_rt2, out_verbose_rt3;
std::ofstream out_verbose_fx, out_verbose_fx2, out_verbose_fx3;
std::ofstream out_verbose_cr, out_verbose_cr2, out_verbose_cr3;

double diagonal_bbox;
bool detectReflectional = true;
bool detectRotational = false;
bool verification = false;

//Vertex property map, with std::map
typedef std::map<Vertex*, int> Vertex2int_map_type;
typedef boost::associative_property_map< Vertex2int_map_type > Vertex_PM_type;
typedef T_PolyhedralSurf_rings<PolyhedralSurf, Vertex_PM_type > Poly_rings;

//Hedge property map, with enriched Halfedge with its length
// typedef HEdge_PM<PolyhedralSurf> Hedge_PM_type;
// typedef T_PolyhedralSurf_hedge_ops<PolyhedralSurf, Hedge_PM_type> Poly_hedge_ops;
//Hedge property map, with std::map
typedef std::map<Halfedge_handle, double, Hedge_cmp> Hedge2double_map_type;
typedef boost::associative_property_map<Hedge2double_map_type> Hedge_PM_type;
typedef T_PolyhedralSurf_hedge_ops<PolyhedralSurf, Hedge_PM_type> Poly_hedge_ops;

// //Facet property map with enriched Facet with its normal
// typedef Facet_PM<PolyhedralSurf> Facet_PM_type;
// typedef T_PolyhedralSurf_facet_ops<PolyhedralSurf, Facet_PM_type> Poly_facet_ops;
//Facet property map, with std::map
typedef std::map<Facet_handle, Vector_3, Facet_cmp> Facet2normal_map_type;
typedef boost::associative_property_map<Facet2normal_map_type> Facet_PM_type;
typedef T_PolyhedralSurf_facet_ops<PolyhedralSurf, Facet_PM_type> Poly_facet_ops;

typedef double LFT;
typedef CGAL::Simple_cartesian<LFT>     Local_Kernel;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;

typedef CGAL::Monge_via_jet_fitting<Kernel>    Monge_via_jet_fitting;
typedef Monge_via_jet_fitting::Monge_form      Monge_form;

// default parameter values and global variables
unsigned int d_fitting = 4; //2
unsigned int d_monge = 4; //2
unsigned int nb_rings = 0;//seek min # of rings to get the required #pts
unsigned int nb_points_to_use = 0;//
bool verbose = true;
bool verbose_fx = false;
bool verbose_cr = false;
unsigned int min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;

typedef std::vector<Point_with_normal>                       pwnVector;	
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// In Efficient_RANSAC_traits the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<Kernel, pwnVector, Point_map, Normal_map> TraitsRansac;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<TraitsRansac> Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<TraitsRansac>            Plane;

//gather points around the vertex v using rings on the
//polyhedralsurf. the collection of points resorts to 3 alternatives:
// 1. the exact number of points to be used
// 2. the exact number of rings to be used
// 3. nothing is specified
void gather_fitting_points(Vertex* v,
			   std::vector<DPoint> &in_points,
			   Vertex_PM_type& vpm)
{
  //container to collect vertices of v on the PolyhedralSurf
  std::vector<Vertex*> gathered;
  //initialize
  in_points.clear();

  //OPTION -p nb_points_to_use, with nb_points_to_use != 0. Collect
  //enough rings and discard some points of the last collected ring to
  //get the exact "nb_points_to_use"
  if ( nb_points_to_use != 0 ) {
    Poly_rings::collect_enough_rings(v, nb_points_to_use, gathered, vpm);
    if ( gathered.size() > nb_points_to_use ) gathered.resize(nb_points_to_use);
  }
  else { // nb_points_to_use=0, this is the default and the option -p is not considered;
    // then option -a nb_rings is checked. If nb_rings=0, collect
    // enough rings to get the min_nb_points required for the fitting
    // else collect the nb_rings required
    if ( nb_rings == 0 )
      Poly_rings::collect_enough_rings(v, min_nb_points, gathered, vpm);
    else Poly_rings::collect_i_rings(v, nb_rings, gathered, vpm);
  }

  //store the gathered points
  std::vector<Vertex*>::iterator
    itb = gathered.begin(), ite = gathered.end();
  CGAL_For_all(itb,ite) in_points.push_back((*itb)->point());
}

#pragma region TransformationPairs 
struct PointSt
{
	DPoint p;
	Vector_3 normal;
	Vector_3 kmin;
	Vector_3 kmax;
	double k1;
	double k2;
};

struct ReflectionPair
{
	int i;
	int j;
	double PlaneA;
	double PlaneB;
	double PlaneC;
	double PlaneD;
	Point_with_normal PointN;
};

struct RotTransPair
{
	int i;
	int j;
	double Scale;
	double Rx;
	double Ry;
	double Rz;
	double Tx;
	double Ty;
	double Tz;
};

struct ConsRotPair
{
	int i1;
	int j1;
	int i2;
	int j2;
	Kernel::Vector_3 n;
	Point_3 p;
};
#pragma endregion 

//These vectors contain the same number of elements
std::vector<PointSt> pointsSampled;
std::vector<Point_3> points; //Points of sample P

std::vector<PointSt> uPointSampledList;
std::vector<Point_3> uPointList; //Points of sample M
std::vector<int> uIndices;

struct ContinuousRotationTransform {
	DPoint operator()(const Vertex& v) const {
		Point_3 vertex(v.point().x(), v.point().y(), v.point().z());
		Kernel::Point_3 lPoint(-0.46, 0.0, -6.97);
		Kernel::Vector_3 lVector(0.0, 1.0, 0.0);
		Kernel::Line_3 line(lPoint, lVector);

		Point_3 proj = line.projection(vertex);
		DPoint dPoint(2.0 * proj.x() - vertex.x(), 2.0 * proj.y() - vertex.y(), 2.0 * proj.z() - vertex.z());
		return dPoint;
	}
};

struct ReflectionTransform {
	DPoint operator()(const Vertex& v) const {
		Point_3 vertex(v.point().x(), v.point().y(), v.point().z());
		Plane_3 plane(1.0, 0.0, 0.0, 0.0);

		Point_3 proj = plane.projection(vertex);
		DPoint dPoint(2.0 * proj.x() - vertex.x(), 2.0 * proj.y() - vertex.y(), 2.0 * proj.z() - vertex.z());
		return dPoint;
	}
};

#pragma region Clustering
struct Cluster {
	RotTransPair mode;
	std::vector<RotTransPair> original_points;
	std::vector<RotTransPair> shifted_points;
};

struct ClusterFX {
	ReflectionPair mode;
	std::vector<ReflectionPair> original_points;
	std::vector<ReflectionPair> shifted_points;
	ReflectionPair transformation;
};

struct ClusterCR {
	ConsRotPair mode;
	std::vector<ConsRotPair> original_points;
	std::vector<ConsRotPair> shifted_points;
	ConsRotPair transformation;
};

struct PointClusterFX {
	std::vector<int> points;
	std::vector<int> reflectedPoints;
};

struct PointClusterRT {
	std::vector<int> points;
	std::vector<int> rotatedPoints;
};

struct PointClusterCR {
	std::vector<int> points;
	std::vector<int> reflectedPoints;
};
#pragma endregion

#pragma region Utilities
Kernel::Vector_3 CorrectDirection(Kernel::Vector_3 nVector)
{
	if (abs(nVector.x()) > abs(nVector.y()))
	{
		if (abs(nVector.x()) > abs(nVector.z()))
		{
			if (nVector.x() < 0) { nVector = -nVector; }
		}
		else
			if (nVector.z() < 0) { nVector = -nVector; }
	}
	else
	{
		if (abs(nVector.y()) > abs(nVector.z()))
		{
			if (nVector.y() < 0) { nVector = -nVector; }
		}
		else
			if (nVector.z() < 0) { nVector = -nVector; }
	}
	return nVector;
}

Plane_3 averagePlane(std::vector<Plane_3> vPlanes)
{
	Point_3 origin(0, 0, 0);
	Point_3 pAccumulator(0, 0, 0);
	Kernel::Vector_3 vAccumulator(0, 0, 0);
	double pX = 0.0;
	double pY = 0.0;
	double pZ = 0.0;
	int cSize = vPlanes.size();

	for (int id = 0; id < cSize; id++)
	{
		Plane_3 plane = vPlanes[id];
		Kernel::Vector_3 vDirection = CorrectDirection(plane.orthogonal_vector());
		Point_3 pt = plane.projection(origin);
		Kernel::Vector_3 vPt = pt - origin;

		pAccumulator = pAccumulator + vPt / cSize;
		vAccumulator = vAccumulator + vDirection / cSize;
	}
	vAccumulator = vAccumulator / sqrt(vAccumulator.squared_length());

	Plane_3 pAverage(pAccumulator, vAccumulator);
	return pAverage;
}

RotTransPair averageTransform(std::vector<RotTransPair> vRotations)
{
	double rX = 0.0;
	double rY = 0.0;
	double rZ = 0.0;

	double tX = 0.0;
	double tY = 0.0;
	double tZ = 0.0;

	int cSize = vRotations.size();
	for (int id = 0; id < cSize; id++)
	{
		rX += vRotations[id].Rx;
		rY += vRotations[id].Ry;
		rZ += vRotations[id].Rz;

		tX += vRotations[id].Tx;
		tY += vRotations[id].Ty;
		tZ += vRotations[id].Tz;
	}
	RotTransPair rtAverage;
	rtAverage.Rx = rX / cSize;
	rtAverage.Ry = rY / cSize;
	rtAverage.Rz = rZ / cSize;
	rtAverage.Tx = tX / cSize;
	rtAverage.Ty = tY / cSize;
	rtAverage.Tz = tZ / cSize;

	return rtAverage;
}

Eigen::Matrix3d alignmentMatrix(Vector_3 pINormal, Vector_3 pJNormal, Vector_3 rAxis)
{
	Eigen::Matrix3d mIdentity = Eigen::Matrix3d::Identity();
	double cosAngle = pINormal * pJNormal;	
	Vector_3 rAxisCP = CGAL::cross_product(pINormal, pJNormal);
	double sinAngle = sqrt(rAxisCP.squared_length());

	if (rAxis * rAxisCP < 0)
		sinAngle = -sinAngle;
	
	rAxisCP = rAxisCP / sinAngle;

	Eigen::Matrix3d crossProduct;
	crossProduct(0, 0) = 0; crossProduct(0, 1) = -rAxisCP.z(); crossProduct(0, 2) = rAxisCP.y();
	crossProduct(1, 0) = rAxisCP.z(); crossProduct(1, 1) = 0; crossProduct(1, 2) = -rAxisCP.x();
	crossProduct(2, 0) = -rAxisCP.y(); crossProduct(2, 1) = rAxisCP.x(); crossProduct(2, 2) = 0;

	Eigen::Vector3d vrAxis(rAxisCP.x(), rAxisCP.y(), rAxisCP.z());
	Eigen::RowVector3d rrAxis = vrAxis.transpose();

	Eigen::Matrix3d tensorProduct = vrAxis * rrAxis;
	//Eigen::Matrix3d tensorProduct = crossProduct * crossProduct;	
	Eigen::Matrix3d mAlign = cosAngle * mIdentity + sinAngle * crossProduct + (1 - cosAngle) * tensorProduct;
	//Eigen::Matrix3d mAlign = mIdentity + crossProduct + (1 - cosAngle) * tensorProduct / sinAngle2;

	//double det = mAlign.determinant();	
	//Eigen::Vector3d cI(pINormal.x(), pINormal.y(), pINormal.z());
	//Eigen::Vector3d cIa = mAlign * cI;
	//Vector_3 pIa(cIa.x(), cIa.y(), cIa.z());
	//double validation = pIa * pJNormal;
	return mAlign;
}

double euclidean_distance(const vector<double> &point_a, const vector<double> &point_b)
{
	double total = 0;
	for (int i = 0; i<point_a.size(); i++) {
		total += (point_a[i] - point_b[i]) * (point_a[i] - point_b[i]);
	}
	return sqrt(total);
}
#pragma endregion

#pragma region Kernel
double gaussian_kernel(double distance, double kernel_bandwidth)
{
	double temp = exp(-1.0 / 2.0 * (distance * distance) / (kernel_bandwidth * kernel_bandwidth));
	return temp;
}

double epanechnikov_kernel(double distance, double kernel_bandwidth)
{
	double temp = std::max(0.0, 1.0 - (distance * distance) / (kernel_bandwidth * kernel_bandwidth));
	return temp;
}
#pragma endregion

#pragma region Reflectional Clustering
double distanceMetricFX(const ReflectionPair &pointA, const ReflectionPair &pointB)
{
	const double B = 3.0 / diagonal_bbox; //0.0025; //25; // 3
	const double A = 3.0; //50; //30 * sqrt(3.0)
	
	Plane_3 planeA(pointA.PlaneA, pointA.PlaneB, pointA.PlaneC, pointA.PlaneD);
	Plane_3 planeB(pointB.PlaneA, pointB.PlaneB, pointB.PlaneC, pointB.PlaneD);
	double cosineN = (planeA.orthogonal_vector() * planeB.orthogonal_vector()) / 
		sqrt(planeA.orthogonal_vector().squared_length() * planeB.orthogonal_vector().squared_length());
	double angle = acos(abs(cosineN));

	Point_3 origin(0, 0, 0);
	double distance = sqrt(CGAL::squared_distance(planeA.projection(origin), planeB.projection(origin)));

	double total = B * distance * angle;
	return total;
}

ReflectionPair shiftPointFX(const ReflectionPair &point, const vector<ReflectionPair> &points, double kernel_bandwidth)
{
	ReflectionPair shifted_point = point;
	shifted_point.PlaneA = 0;
	shifted_point.PlaneB = 0;
	shifted_point.PlaneC = 0;
	shifted_point.PlaneD = 0;
	double total_weight = 0;
	for (int i = 0; i < points.size(); i++) {
		ReflectionPair temp_point = points[i];
		double distance = distanceMetricFX(point, temp_point);
		
		if (distance < kernel_bandwidth)
		{
			double weight = epanechnikov_kernel(distance, kernel_bandwidth); //gaussian_kernel
			shifted_point.PlaneA += temp_point.PlaneA * weight;
			shifted_point.PlaneB += temp_point.PlaneB * weight;
			shifted_point.PlaneC += temp_point.PlaneC * weight;
			shifted_point.PlaneD += temp_point.PlaneD * weight;
			total_weight += weight;
		}
	}
	if (total_weight > 0)
	{
		shifted_point.PlaneA /= total_weight;
		shifted_point.PlaneB /= total_weight;
		shifted_point.PlaneC /= total_weight;
		shifted_point.PlaneD /= total_weight;
		return shifted_point;
	}
	else
		return point;
}

vector<ReflectionPair> meanshiftFX(const vector<ReflectionPair> &points, double kernel_bandwidth)
{
	vector<bool> stop_moving(points.size(), false);
	vector<ReflectionPair> shifted_points = points;
	double max_shift_distance;

	do {
		max_shift_distance = 0;
		for (int i = 0; i < shifted_points.size(); i++)
		{
			if (!stop_moving[i])
			{
				ReflectionPair point_new = shiftPointFX(shifted_points[i], points, kernel_bandwidth);
				double shift_distance = distanceMetricFX(point_new, shifted_points[i]);
				if (shift_distance > max_shift_distance) {
					max_shift_distance = shift_distance;
				}
				if (shift_distance <= EPSILON) {
					stop_moving[i] = true;
				}
				shifted_points[i] = point_new;
			}
		}
	} while (max_shift_distance > EPSILON);

	return shifted_points;
}

vector<ClusterFX> clusterFX(const vector<ReflectionPair> &points, const vector<ReflectionPair> &shifted_points)
{
	vector<ClusterFX> clusters;
	for (int i = 0; i < shifted_points.size(); i++) {
		int c = 0;
		for (; c < clusters.size(); c++) {
			if (distanceMetricFX(shifted_points[i], clusters[c].mode) <= CLUSTER_EPSILON) {
				break;
			}
		}
		if (c == clusters.size()) {
			ClusterFX clus;
			clus.mode = shifted_points[i];
			clusters.push_back(clus);
		}
		clusters[c].original_points.push_back(points[i]);
		clusters[c].shifted_points.push_back(shifted_points[i]);
	}

	return clusters;
}

vector<ClusterFX> clusterFX(const vector<ReflectionPair> &points, double kernel_bandwidth)
{
	vector<ReflectionPair> shifted_points = meanshiftFX(points, kernel_bandwidth);
	return clusterFX(points, shifted_points);
}
#pragma endregion

#pragma region Rotational-Traslational Clustering
double distanceMetricRT(const RotTransPair &pointA, const RotTransPair &pointB)
{
	//const double B1 = 0.5 * pow(1 / 10.0, 2);
	const double B2 = 4.0; // *pow(1 / M_PI, 2);
	const double B3 = 0.0016; //pow(1 / 200.0, 2); //0.5 * diagonal_bbox;

	double total = //B1 * pow(pointA.Scale - pointB.Scale, 2) +
		B2 * (pow(pointA.Rx - pointB.Rx, 2) + 4 * pow(pointA.Ry - pointB.Ry, 2) + pow(pointA.Rz - pointB.Rz, 2));// +
		//B3 * (pow(pointA.Tx - pointB.Tx, 2) + pow(pointA.Ty - pointB.Ty, 2) + pow(pointA.Tz - pointB.Tz, 2));
	return sqrt(total);
}

RotTransPair shiftPointRT(const RotTransPair &point, const vector<RotTransPair> &points, double kernel_bandwidth)
{
	RotTransPair shifted_point = point;
	shifted_point.Rx = 0;
	shifted_point.Ry = 0;
	shifted_point.Rz = 0;
	shifted_point.Tx = 0;
	shifted_point.Ty = 0;
	shifted_point.Tz = 0;
	//shifted_point.Scale = 0;
	
	double total_weight = 0;
	for (int i = 0; i < points.size(); i++) {
		RotTransPair temp_point = points[i];
		double distance = distanceMetricRT(point, temp_point);
		
		if (distance < kernel_bandwidth)
		{
			double weight = epanechnikov_kernel(distance, kernel_bandwidth); //gaussian_kernel
			shifted_point.Rx += temp_point.Rx * weight;
			shifted_point.Ry += temp_point.Ry * weight;
			shifted_point.Rz += temp_point.Rz * weight;
			shifted_point.Tx += temp_point.Tx * weight;
			shifted_point.Ty += temp_point.Ty * weight;
			shifted_point.Tz += temp_point.Tz * weight;
			//shifted_point.Scale += temp_point.Scale * weight;
			total_weight += weight;
		}
	}
	
	shifted_point.Rx /= total_weight;
	shifted_point.Ry /= total_weight;
	shifted_point.Rz /= total_weight;
	shifted_point.Tx /= total_weight;
	shifted_point.Ty /= total_weight;
	shifted_point.Tz /= total_weight;
	//shifted_point.Scale /= total_weight;
	
	return shifted_point;
}

vector<RotTransPair> meanshiftRT(const vector<RotTransPair> &points, double kernel_bandwidth)
{
	vector<bool> stop_moving(points.size(), false);
	vector<RotTransPair> shifted_points = points;
	double max_shift_distance;

	do {
		max_shift_distance = 0;
		for (int i = 0; i < shifted_points.size(); i++)
		{
			if (!stop_moving[i]) 
			{
				RotTransPair point_new = shiftPointRT(shifted_points[i], points, kernel_bandwidth);
				double shift_distance = distanceMetricRT(point_new, shifted_points[i]);
				if (shift_distance > max_shift_distance) {
					max_shift_distance = shift_distance;
				}
				if (shift_distance <= EPSILON) {
					stop_moving[i] = true;
				}
				shifted_points[i] = point_new;
			}
		}		
	} while (max_shift_distance > EPSILON);

	return shifted_points;
}

vector<Cluster> clusterRT(const vector<RotTransPair> &points, const vector<RotTransPair> &shifted_points)
{
	vector<Cluster> clusters;
	for (int i = 0; i < shifted_points.size(); i++) {
		int c = 0;
		for (; c < clusters.size(); c++) {
			if (distanceMetricRT(shifted_points[i], clusters[c].mode) <= CLUSTER_EPSILON) {
				break;
			}
		}
		if (c == clusters.size()) {
			Cluster clus;
			clus.mode = shifted_points[i];
			clusters.push_back(clus);
		}
		clusters[c].original_points.push_back(points[i]);
		clusters[c].shifted_points.push_back(shifted_points[i]);
	}

	return clusters;
}

vector<Cluster> clusterRT(const vector<RotTransPair> &points, double kernel_bandwidth)
{
	vector<RotTransPair> shifted_points = meanshiftRT(points, kernel_bandwidth);
	return clusterRT(points, shifted_points);
}
#pragma endregion

#pragma region Continous Rotational Clustering
double distanceMetricCR(const ConsRotPair &pointA, const ConsRotPair &pointB)
{
	double cosine = (pointA.n * pointB.n) / sqrt(pointA.n.squared_length() * pointB.n.squared_length());
	double angle = acos(abs(cosine));
	Line_3 lA(pointA.p, pointA.n);
	Line_3 lB(pointB.p, pointB.n);
	double distance = sqrt(CGAL::squared_distance(pointA.p, pointB.p));
	//double distance = sqrt(CGAL::squared_distance(lA, lB));

	const double B = 1 / diagonal_bbox;
	const double A = 3.0;
	double total = A * angle + B * distance;
	return total;
}

ConsRotPair shiftPointCR(const ConsRotPair &point, const vector<ConsRotPair> &points, double kernel_bandwidth)
{
	double nx = 0;
	double ny = 0;
	double nz = 0;
	double px = 0;
	double py = 0;
	double pz = 0;
	
	double total_weight = 0;
	for (int i = 0; i < points.size(); i++) {
		ConsRotPair temp_point = points[i];
		double distance = distanceMetricCR(point, temp_point);

		if (distance < kernel_bandwidth)
		{
			double weight = epanechnikov_kernel(distance, kernel_bandwidth); //gaussian_kernel
			nx += temp_point.n.x() * weight;
			ny += temp_point.n.y() * weight;
			nz += temp_point.n.z() * weight;

			px += temp_point.p.x() * weight;
			py += temp_point.p.y() * weight;
			pz += temp_point.p.z() * weight;
			total_weight += weight;
		}
	}
	ConsRotPair shifted_point;

	if (total_weight > 0)
	{
		nx /= total_weight;
		ny /= total_weight;
		nz /= total_weight;

		px /= total_weight;
		py /= total_weight;
		pz /= total_weight;

		Kernel::Vector_3 n(nx, ny, nz);
		n = n / sqrt(n.squared_length());
		Point_3 p(px, py, pz);
		shifted_point.n = n;
		shifted_point.p = p;
	}
	else
		shifted_point = point;
	
	return shifted_point;
}

vector<ConsRotPair> meanshiftCR(const vector<ConsRotPair> &points, double kernel_bandwidth)
{
	vector<bool> stop_moving(points.size(), false);
	vector<ConsRotPair> shifted_points = points;
	double max_shift_distance;

	do {
		max_shift_distance = 0;
		for (int i = 0; i < shifted_points.size(); i++)
		{
			if (!stop_moving[i])
			{
				ConsRotPair point_new = shiftPointCR(shifted_points[i], points, kernel_bandwidth);
				double shift_distance = distanceMetricCR(point_new, shifted_points[i]);
				if (shift_distance > max_shift_distance) {
					max_shift_distance = shift_distance;
				}
				if (shift_distance <= EPSILON) {
					stop_moving[i] = true;
				}
				shifted_points[i] = point_new;
			}
		}
	} while (max_shift_distance > EPSILON);

	return shifted_points;
}

vector<ClusterCR> clusterCR(const vector<ConsRotPair> &points, const vector<ConsRotPair> &shifted_points)
{
	vector<ClusterCR> clusters;
	for (int i = 0; i < shifted_points.size(); i++) {
		int c = 0;
		for (; c < clusters.size(); c++) {
			if (distanceMetricCR(shifted_points[i], clusters[c].mode) <= CLUSTER_EPSILON) {
				break;
			}
		}
		if (c == clusters.size()) {
			ClusterCR clus;
			clus.mode = shifted_points[i];
			clusters.push_back(clus);
		}
		clusters[c].original_points.push_back(points[i]);
		clusters[c].shifted_points.push_back(shifted_points[i]);
	}

	return clusters;
}

vector<ClusterCR> clusterCR(const vector<ConsRotPair> &points, double kernel_bandwidth)
{
	vector<ConsRotPair> shifted_points = meanshiftCR(points, kernel_bandwidth);
	return clusterCR(points, shifted_points);
}
#pragma endregion

#pragma region Verification
Point_3 reflectPointPlane(PointSt point, ReflectionPair plane)
{
	Eigen::Vector3d vNeighbor(point.p.x(), point.p.y(), point.p.z());
	Eigen::Vector3d vPlaneNormal(plane.PlaneA, plane.PlaneB, plane.PlaneC);

	double k = 2.0 * (plane.PlaneA * vNeighbor.x() + plane.PlaneB * vNeighbor.y() + plane.PlaneC * vNeighbor.z() + plane.PlaneD);
	Eigen::Vector3d pTrNeighbor = vNeighbor - k * vPlaneNormal;

	//We must find the closests points in the sample, that must be located within some threshold.
	Point_3 qTrPoint(pTrNeighbor.x(), pTrNeighbor.y(), pTrNeighbor.z());
	return qTrPoint;
}

Point_3 reflectPointLine(PointSt point, ConsRotPair line)
{
	Point_3 vertex(point.p.x(), point.p.y(), point.p.z());
	Kernel::Line_3 linea(line.p, line.n);
	Point_3 proj = linea.projection(vertex);
	Point_3 qTrPoint(2.0 * proj.x() - vertex.x(), 2.0 * proj.y() - vertex.y(), 2.0 * proj.z() - vertex.z());
	return qTrPoint;
}

std::vector<int> SearchNeighbors(int ixPoint, bool subSample)
{
	//points
	PointPropertyMap ppmap(uPointList);
	// Insert number_of_data_points in the tree
	Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(uPointList.size()), Splitter(), Traits(ppmap));
	Distance dist(ppmap);
		
	int K = 8;
	double distanceThreshold = 3.0;
	std::vector<int> neighbors;	
	Point_3 queryI = uPointList[ixPoint];

	// search K nearest neighbours
	K_neighbor_search search(tree, queryI, K, 0, true, dist);
	K_neighbor_search::iterator it = search.begin();
	bool thresholdExceeded = false;

	while (it != search.end() && !thresholdExceeded)
	{
		if (subSample)
		{
			if (dist.inverse_of_transformed_distance(it->second) < distanceThreshold)
				neighbors.push_back(it->first);
			else
				thresholdExceeded = true;
		}
		else
		{
			if (ixPoint != it->first && dist.inverse_of_transformed_distance(it->second) < distanceThreshold)
				neighbors.push_back(it->first);
			else
				thresholdExceeded = true;
		}
		it++;
	}

	return neighbors;
}

int VerifyPointInClusterFx(int ixPoint, ReflectionPair rPair)
{
	PointPropertyMap ppmap(uPointList);
	// Insert number_of_data_points in the tree
	Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(uPointList.size()), Splitter(), Traits(ppmap));
	Distance dist(ppmap);

	int K = 2;
	double distanceThreshold = 3.0; // 5.0; // 1.25;
	double rAngleThreshold = 0.25;
	//PointSt pointSt = pointsSampled[ixPoint];
	PointSt pointSt = uPointSampledList[ixPoint];

	//We must find the closests points in the sample, that must be located within some threshold.
	Point_3 reflectedPoint = reflectPointPlane(pointSt, rPair);
	K_neighbor_search search(tree, reflectedPoint, K, 0, true, dist);
	K_neighbor_search::iterator it = search.begin();
		
	Vector_3 vtPlaneNormal(rPair.PlaneA, rPair.PlaneB, rPair.PlaneC);
	
	bool found = false;
	int ixReflected = -1;
	//Verify if the point's directions have coherence respect with the converted point.
	while (it != search.end() && !found)
	{
		//if (it->first != rPair.j)
		//{
			double distanceSeedFX = dist.inverse_of_transformed_distance(it->second);

			if ((distanceSeedFX < distanceThreshold)
				//&& (abs((pointSt.normal + uPointSampledList[it->first].normal) * vtPlaneNormal) < rAngleThreshold)
				)
			{
				ixReflected = it->first;
				found = true;
			}
		//}
		it++;
	}
	return ixReflected;
}

int VerifyPointInClusterCr(int ixPoint, ConsRotPair rPair)
{
	PointPropertyMap ppmap(uPointList);
	// Insert number_of_data_points in the tree
	Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(uPointList.size()), Splitter(), Traits(ppmap));
	Distance dist(ppmap);

	int K = 2;
	double distanceThreshold = 3.0; // 5.0; // 1.25;
	double rAngleThreshold = 0.25;
	//PointSt pointSt = pointsSampled[ixPoint];
	PointSt pointSt = uPointSampledList[ixPoint];

	//We must find the closests points in the sample, that must be located within some threshold.
	Point_3 reflectedPoint = reflectPointLine(pointSt, rPair);
	K_neighbor_search search(tree, reflectedPoint, K, 0, true, dist);
	K_neighbor_search::iterator it = search.begin();

	Vector_3 vtLineNormal(rPair.n.x(), rPair.n.y(), rPair.n.z());

	bool found = false;
	int ixReflected = -1;
	//Verify if the point's directions have coherence respect with the converted point.
	while (it != search.end() && !found)
	{
		//if (it->first != rPair.j1)
		//{
			double distanceSeedCR = dist.inverse_of_transformed_distance(it->second);

			if ((distanceSeedCR < distanceThreshold)
				//&& (abs((pointSt.normal + uPointSampledList[it->first].normal) * vtLineNormal) < rAngleThreshold)
				)
			{
				ixReflected = it->first;
				found = true;
			}
		//}
		it++;
	}
	return ixReflected;
}

int VerifyPointInClusterRt(int ixPoint, RotTransPair rPair)
{
	PointPropertyMap ppmap(points);
	// Insert number_of_data_points in the tree
	Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(points.size()), Splitter(), Traits(ppmap));
	Distance dist(ppmap);

	int K = 9;
	double distanceThreshold = 2.5; // 1.25;
	double rAngleThreshold = 0.75;
	
	//PointSt pointSt = pointsSampled[ixPoint];
	PointSt pointSt = uPointSampledList[ixPoint];

	//We must find the closests points in the sample, that must be located within some threshold.
	Eigen::Vector3d traslation(rPair.Tx, rPair.Ty, rPair.Tz);
	Eigen::Vector3d vPointSt(pointSt.p.x(), pointSt.p.y(), pointSt.p.z());
	Eigen::Vector3d vPointKmin(pointSt.kmin.x(), pointSt.kmin.y(), pointSt.kmin.z());
	Eigen::Vector3d vPointKmax(pointSt.kmax.x(), pointSt.kmax.y(), pointSt.kmax.z());
	Eigen::Vector3d vPointNormal(pointSt.normal.x(), pointSt.normal.y(), pointSt.normal.z());

	Eigen::Matrix3d rotationX;
	rotationX(0, 0) = 1;	rotationX(0, 1) = 0;				rotationX(0, 2) = 0;
	rotationX(1, 0) = 0;	rotationX(1, 1) = cosf(rPair.Rx);	rotationX(1, 2) = -sinf(rPair.Rx);
	rotationX(2, 0) = 0;	rotationX(2, 1) = sinf(rPair.Rx);	rotationX(2, 2) = cosf(rPair.Rx);
	
	Eigen::Matrix3d rotationY;
	rotationY(0, 0) = cosf(rPair.Ry);	rotationY(0, 1) = 0;	rotationY(0, 2) = sinf(rPair.Ry);
	rotationY(1, 0) = 0;				rotationY(1, 1) = 1;	rotationY(1, 2) = 0;
	rotationY(2, 0) = -sinf(rPair.Ry);	rotationY(2, 1) = 0;	rotationY(2, 2) = cosf(rPair.Ry);
	
	Eigen::Matrix3d rotationZ;
	rotationZ(0, 0) = cosf(rPair.Rz);	rotationZ(0, 1) = -sinf(rPair.Rz);	rotationZ(0, 2) = 0;
	rotationZ(1, 0) = sinf(rPair.Rz);	rotationZ(1, 1) = cosf(rPair.Rz);	rotationZ(1, 2) = 0;
	rotationZ(2, 0) = 0;				rotationZ(2, 1) = 0;				rotationZ(2, 2) = 1;

	Eigen::Matrix3d rotationMatrix = rotationZ * rotationY * rotationX;
	Eigen::Vector3d vTransformedPoint = rotationMatrix * vPointSt + traslation;
	Eigen::Vector3d vTransformedNormal = rotationMatrix * vPointNormal;
	Eigen::Vector3d vTransformedKmin = rotationMatrix * vPointKmin;
	Eigen::Vector3d vTransformedKmax = rotationMatrix * vPointKmax;

	Point_3 transformedPoint(vTransformedPoint.x(), vTransformedPoint.y(), vTransformedPoint.z());
	Vector_3 transformedNormal(vTransformedNormal.x(), vTransformedNormal.y(), vTransformedNormal.z());
	Vector_3 transformedKmin(vTransformedKmin.x(), vTransformedKmin.y(), vTransformedKmin.z());
	Vector_3 transformedKmax(vTransformedKmax.x(), vTransformedKmax.y(), vTransformedKmax.z());
	K_neighbor_search search(tree, transformedPoint, K, 0, true, dist);
	K_neighbor_search::iterator it = search.begin();
	
	bool found = false;
	int ixReflected = -1;
	//Verify if the point's directions have coherence respect with the converted point.
	while (it != search.end() && !found)
	{
		if (it->first != rPair.j)
		{
			double distanceSeedRN = dist.inverse_of_transformed_distance(it->second);

			//if (distanceSeedRN < distanceThreshold
			//	&& (abs(transformedNormal * pointsSampled[it->first].normal) > rAngleThreshold)
			//	&& (abs(transformedKmin * pointsSampled[it->first].kmin) > rAngleThreshold)
			//	&& (abs(transformedKmax * pointsSampled[it->first].kmax) > rAngleThreshold))

			if (distanceSeedRN < distanceThreshold)				
			{
				if ((abs(transformedNormal * uPointSampledList[it->first].normal) > rAngleThreshold)
					&& (abs(transformedKmin * uPointSampledList[it->first].kmin) > rAngleThreshold)
					&& (abs(transformedKmax * uPointSampledList[it->first].kmax) > rAngleThreshold))
				{
					ixReflected = it->first;
					found = true;
				}
			}
			else
				break;
		}
		it++;
	}
	return ixReflected;
}
#pragma endregion 

#pragma region Detection

void RefineRotational(std::vector<Cluster> clustersRt)
{
	std::vector<Cluster> clustersRtText;
	std::vector<PointClusterRT> patchedClustersRT;

	for (int ix = 0; ix < clustersRt.size(); ix++)
	{
		if (clustersRt[ix].original_points.size() > 3)
		{
			//pointsSampled
			clustersRtText.push_back(clustersRt[ix]);
			//std::cout << "Cluster " << ix << " mode: " << clustersFx[ix].mode.i << " (" << points[clustersFx[ix].mode.i] << ")";

			PointSt point_mode = pointsSampled[clustersRt[ix].mode.i];
			//ReflectionPair rPair = clustersFx[ix].original_points[0]; //Only the first reflection Pair?

			std::vector<int> pointsCluster;
			std::vector<int> pointsClusterRt;
			//pointsCluster.resize(clustersFx[ix].original_points.size());

			//pointsCluster.push_back(clustersFx[ix].mode.i);
			for (int s = 0; s < clustersRt[ix].original_points.size(); s++)
			{
				//Add the points included in the cluster from the seed.
				RotTransPair rPair = clustersRt[ix].original_points[s];

				std::vector<int> pcNeighborhood;
				std::vector<int> pcNeighborhoodRx;
				pcNeighborhood.push_back(uIndices[rPair.i]); //Replace with index in the uList
				pcNeighborhoodRx.push_back(uIndices[rPair.j]);

				for (int k = 0; k < pcNeighborhood.size(); k++)
				{
					bool pointClusterIncluded = false;
					std::vector<int> neighbors; // = SearchNeighbors(pcNeighborhood[k]);

					for (int nI = 0; nI < neighbors.size(); nI++)
					{
						bool alreadyIncluded = false;

						//Check if the point has not been yet included.
						int nP = 0, nC = 0;
						while (nP < pcNeighborhood.size() && !alreadyIncluded)
						{
							if (pcNeighborhood[nP] == neighbors[nI])
								alreadyIncluded = true;
							nP++;
						}
						while (nC < pointsCluster.size() && !alreadyIncluded)
						{
							if (pointsCluster[nC] == neighbors[nI])
								alreadyIncluded = true;
							nC++;
						}

						int rfPoint = VerifyPointInClusterRt(neighbors[nI], rPair);
						if (!alreadyIncluded && rfPoint != -1)
						{
							pcNeighborhood.push_back(neighbors[nI]);
							pcNeighborhoodRx.push_back(rfPoint);
						}
					}
				}

				if (pcNeighborhood.size() > 1)
					for (int p = 0; p < pcNeighborhood.size(); p++)
					{
						pointsCluster.push_back(pcNeighborhood[p]);
						pointsClusterRt.push_back(pcNeighborhoodRx[p]);
					}
			}

			PointClusterRT pClusterRT;
			pClusterRT.points = pointsCluster;
			pClusterRT.rotatedPoints = pointsClusterRt;
			patchedClustersRT.push_back(pClusterRT);
		}
	}

	std::cout << std::setprecision(3) << std::fixed;
	std::cout << "Number of patches: " << patchedClustersRT.size() << endl;

	for (int ix = 0; ix < patchedClustersRT.size(); ix++)
	{
		if (patchedClustersRT[ix].points.size() > 0)
		{
			RotTransPair mode1 = clustersRtText[ix].mode;
			std::cout << "Cluster " << ix << " size: " << patchedClustersRT[ix].points.size() << " mode" <<
				" x: " << mode1.Rx << " y: " << mode1.Ry << " z: " << mode1.Rz << endl;
		}
	}

	for (int ix = 0; ix < clustersRt.size(); ix++)
	{
		if (clustersRt[ix].original_points.size() > 3)
			clustersRtText.push_back(clustersRt[ix]);
	}
}

void DetectionDiscreteRotational(std::vector<RotTransPair> vRotTransPair)
{
	const clock_t cluster_begin = clock();
	
	bool refine = false;
	std::vector<RotTransPair> shifted_pointsRT = meanshiftRT(vRotTransPair, KERNEL_BW);

	if (verbose)
	{
		for (int k = 0; k < shifted_pointsRT.size(); k++)
		{
			for (int m = 0; m < shifted_pointsRT.size(); m++)
			{
				double oDistance = distanceMetricRT(vRotTransPair[k], vRotTransPair[m]);
				double sDistance = distanceMetricRT(shifted_pointsRT[k], shifted_pointsRT[m]);

				if (m < shifted_pointsRT.size() - 1)
				{
					out_verbose_rt << k << "-" << m << "|";
					out_verbose_rt2 << oDistance << "|";
					out_verbose_rt3 << sDistance << "|";
				}
				else
				{
					out_verbose_rt << k << "-" << m;
					out_verbose_rt2 << oDistance;
					out_verbose_rt3 << sDistance;
				}
			}
			out_verbose_rt2 << std::endl;
			out_verbose_rt3 << std::endl;
		}
	}

	vector<Cluster> clustersRt = clusterRT(vRotTransPair, shifted_pointsRT);
	std::cout << "Number of clusters: " << clustersRt.size() << endl;	
	
	const clock_t cluster_end = clock();
	std::cout << "Time of cluster: " << float(cluster_end - cluster_begin) / CLOCKS_PER_SEC << " seconds." << endl;

	for (int ix = 0; ix < clustersRt.size(); ix++)
	{
		int cSize = clustersRt[ix].original_points.size();
		RotTransPair modeRt = clustersRt[ix].mode;
		std::cout << "Cluster " << ix << " size: " << cSize << " mode" << " x: " << modeRt.Rx << " y: " << modeRt.Ry << " z: " << modeRt.Rz << endl;
				
		if (cSize > 4)
		{
			std::vector<RotTransPair> vTransformations;
			RotTransPair avRotation = averageTransform(clustersRt[ix].original_points);
			std::cout << "        " << ix << " average: " << avRotation.Rx << "x " << avRotation.Ry << "y " << avRotation.Rz << "z " << endl;
			int nCount = 0;

			for (int id = 0; id < cSize; id++)
			{
				double distance = distanceMetricRT(avRotation, clustersRt[ix].original_points[id]);
				
				if (distance < CLUSTER_EPSILON)
				{
					vTransformations.push_back(clustersRt[ix].original_points[id]);
					nCount++;
				}
			}

			RotTransPair avRotationN = averageTransform(vTransformations);
			if (nCount > 0)
				std::cout << "        " << ix << " size: " << nCount << " outlier: " << avRotationN.Rx << "x " << avRotationN.Ry << "y " << avRotationN.Rz << "z " << endl;
			else
				std::cout << "        No average available" << endl;
		}
	}

	if (refine)
		RefineRotational(clustersRt);
}

void RefineReflectional(ClusterFX clusterFx, ReflectionPair modeFxCenter)
{
	std::vector<int> pointsCluster;
	std::vector<int> pointsClusterFx;
	PointPropertyMap ppmap(uPointList);
	// Insert number_of_data_points in the tree
	Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(uPointList.size()), Splitter(), Traits(ppmap));
	Distance dist(ppmap);

	for (int s = 0; s < clusterFx.original_points.size(); s++)
	{
		//Add the points included in the cluster from the seed.
		ReflectionPair initialRPair = clusterFx.original_points[s];

		std::vector<int> pcNeighborhood;
		std::vector<int> pcNeighborhoodRx;

		int k = 0;
		bool firstPoint = true;
		int sampleIx = uIndices[initialRPair.i];
		int kTree = 8;
		double distanceThreshold = 3.0;
		
		while (firstPoint || k < pcNeighborhood.size())
		{
			std::vector<int> neighbors;

			if (k == 0)
			{
				//neighbors = SearchNeighbors(sampleIx, true);
				Point_3 queryI = uPointList[sampleIx];
				// search K nearest neighbours
				K_neighbor_search search(tree, queryI, kTree, 0, true, dist);
				K_neighbor_search::iterator it = search.begin();
				bool thresholdExceeded = false;

				while (it != search.end() && !thresholdExceeded)
				{					
					if (dist.inverse_of_transformed_distance(it->second) < distanceThreshold)
						neighbors.push_back(it->first);
					else
						thresholdExceeded = true;
					it++;
				}
			}
			else
			{
				//neighbors = SearchNeighbors(pcNeighborhood[k], false);
				Point_3 queryI = uPointList[pcNeighborhood[k]];
				// search K nearest neighbours
				K_neighbor_search search(tree, queryI, kTree, 0, true, dist);
				K_neighbor_search::iterator it = search.begin();
				bool thresholdExceeded = false;

				while (it != search.end() && !thresholdExceeded)
				{					
					if (pcNeighborhood[k] != it->first && dist.inverse_of_transformed_distance(it->second) < distanceThreshold)
						neighbors.push_back(it->first);
					else
						thresholdExceeded = true;					
					it++;
				}
			}
			
			for (int nI = 0; nI < neighbors.size(); nI++)
			{
				bool neighborhoodIncluded = false;
				//Check if the point has not been yet included.
				int nP = 0, nC = 0;
				while (nP < pcNeighborhood.size() && !neighborhoodIncluded)
				{
					if (neighbors[nI] == pcNeighborhood[nP])
						neighborhoodIncluded = true;
					nP++;
				}

				if (!neighborhoodIncluded)
				{
					bool pointClusterIncluded = false;
					int nC = 0;
					while (nC < pointsCluster.size() && !pointClusterIncluded)
					{
						if (pointsCluster[nC] == neighbors[nI])
							pointClusterIncluded = true;
						nC++;
					}

					//int rfPoint = VerifyPointInClusterFx(neighbors[nI], modeFxCenter);
					int K = 2;
					double distanceThreshold = 3.0; // 5.0; // 1.25;
					double rAngleThreshold = 0.25;
					//PointSt pointSt = pointsSampled[ixPoint];
					PointSt pointSt = uPointSampledList[neighbors[nI]];

					//We must find the closests points in the sample, that must be located within some threshold.
					Point_3 reflectedPoint = reflectPointPlane(pointSt, modeFxCenter);
					K_neighbor_search search(tree, reflectedPoint, K, 0, true, dist);
					K_neighbor_search::iterator it = search.begin();

					Vector_3 vtPlaneNormal(modeFxCenter.PlaneA, modeFxCenter.PlaneB, modeFxCenter.PlaneC);

					bool found = false;
					int rfPoint = -1;
					//Verify if the point's directions have coherence respect with the converted point.
					while (it != search.end() && !found)
					{
						double distanceSeedFX = dist.inverse_of_transformed_distance(it->second);

						if ((distanceSeedFX < distanceThreshold))
						{
							rfPoint = it->first;
							found = true;
						}
						it++;
					}

					if (!pointClusterIncluded && rfPoint != -1)
					{
						pcNeighborhood.push_back(neighbors[nI]);
						pcNeighborhoodRx.push_back(rfPoint);
					}
				}
			}
			
			k++;
			firstPoint = false;
		}

		for (int p = 0; p < pcNeighborhood.size(); p++)
		{
			pointsCluster.push_back(pcNeighborhood[p]);
			pointsClusterFx.push_back(pcNeighborhoodRx[p]);
		}
	}	
	std::cout << "         Corrected size: " << pointsCluster.size() << endl;
}

void DetectionReflectional(std::vector<ReflectionPair> vReflectionPair)
{
	const clock_t cluster_begin = clock();
	bool refine = true;
	std::vector<ReflectionPair> shifted_pointsFX = meanshiftFX(vReflectionPair, KERNEL_BW);

	if (verbose_fx)
	{
		for (int k = 0; k < shifted_pointsFX.size(); k++)
		{
			for (int m = 0; m < shifted_pointsFX.size(); m++)
			{
				double oDistance = distanceMetricFX(vReflectionPair[k], vReflectionPair[m]);
				double sDistance = distanceMetricFX(shifted_pointsFX[k], shifted_pointsFX[m]);

				if (m < shifted_pointsFX.size() - 1)
				{
					//out_verbose << shifted_pointsFX[k].i << "-" << shifted_pointsFX[k].j << "|";
					out_verbose_fx << k << "-" << m << "|";
					out_verbose_fx2 << oDistance << "|";
					out_verbose_fx3 << sDistance << "|";
				}
				else
				{
					out_verbose_fx << k << "-" << m;
					out_verbose_fx2 << oDistance;
					out_verbose_fx3 << sDistance;
				}
			}
			out_verbose_fx << std::endl;
			out_verbose_fx2 << std::endl;
			out_verbose_fx3 << std::endl;
		}
	}
	vector<ClusterFX> clustersFx = clusterFX(vReflectionPair, shifted_pointsFX);
	vector<ClusterFX> clustersFxText;
	std::vector<PointClusterFX> patchedClustersFX;

	std::cout << "Number of clusters: " << clustersFx.size() << endl;
	const clock_t cluster_end = clock();
	std::cout << "Time of cluster: " << float(cluster_end - cluster_begin) / CLOCKS_PER_SEC << " seconds." << endl;

	int sampleSize = 0;
	for (int ix = 0; ix < clustersFx.size(); ix++)
		sampleSize += clustersFx[ix].original_points.size();

	for (int ix = 0; ix < clustersFx.size(); ix++)
	{
		int cSize = clustersFx[ix].original_points.size();

		if (cSize > 0.05 * sampleSize)
		{
			ReflectionPair modeFx = clustersFx[ix].original_points[0];

			int idx = -1;
			double distance = 100000.0;

			for (int id = 0; id < cSize; id++)
			{
				double sDistance = 0;

				for (int id_aux = 0; id_aux < cSize; id_aux++)
					sDistance += distanceMetricFX(clustersFx[ix].original_points[id], clustersFx[ix].original_points[id_aux]);

				if (distance > sDistance)
				{
					distance = sDistance;
					idx = id;
				}
			}

			ReflectionPair modeFxCenter = clustersFx[ix].original_points[idx];
			std::cout << "Cluster " << ix << " size: " << cSize << " center: " << modeFxCenter.PlaneA << "x " << modeFxCenter.PlaneB << "y " <<
				modeFxCenter.PlaneC << "z " << modeFxCenter.PlaneD << endl;
			
			Kernel::Vector_3 vCenter();
			Kernel::Vector_3 sumDistance(0.0, 0.0, 0.0);
			double sSquaredError = 0.0;

			for (int k = 0; k < clustersFx[ix].original_points.size(); k++)
			{
				DPoint pJ = pointsSampled[clustersFx[ix].original_points[k].j].p;
				Point_3 pComparison(pJ.x(), pJ.y(), pJ.z());
				Point_3 pReflected = reflectPointPlane(pointsSampled[clustersFx[ix].original_points[k].i], modeFxCenter);
				
				Kernel::Vector_3 vDistance = pComparison - pReflected;
				sumDistance = sumDistance + vDistance;
				sSquaredError += vDistance.squared_length();
			}
			sumDistance = sumDistance / clustersFx[ix].original_points.size();
			sSquaredError = sqrt(sSquaredError / clustersFx[ix].original_points.size());

			if (refine)
				RefineReflectional(clustersFx[ix], modeFxCenter);
		}
	}
}

void RefineContinuousRotational(ClusterCR clusterCr, ConsRotPair modeCrCenter)
{
	std::vector<int> pointsCluster;
	std::vector<int> pointsClusterCr;
	PointPropertyMap ppmap(uPointList);
	// Insert number_of_data_points in the tree
	Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(uPointList.size()), Splitter(), Traits(ppmap));
	Distance dist(ppmap);

	for (int s = 0; s < clusterCr.original_points.size(); s++)
	{
		//Add the points included in the cluster from the seed.
		ConsRotPair initialRPair = clusterCr.original_points[s];

		std::vector<int> pcNeighborhood;
		std::vector<int> pcNeighborhoodCr;

		int k = 0;
		int originPoint = 0;
		std::vector<int> sampleIx;
		sampleIx.resize(4);
		sampleIx[0] = uIndices[initialRPair.i1];
		sampleIx[1] = uIndices[initialRPair.i2];
		sampleIx[2] = uIndices[initialRPair.j1];
		sampleIx[3] = uIndices[initialRPair.j2];
		int kTree = 8;
		double distanceThreshold = 3.0;

		while (originPoint  < 4 || k < pcNeighborhood.size())
		{
			std::vector<int> neighbors;

			if (originPoint < 4)
			{
				//neighbors = SearchNeighbors(sampleIx[originPoint], true);
				Point_3 queryI = uPointList[sampleIx[originPoint]];
				// search K nearest neighbours
				K_neighbor_search search(tree, queryI, kTree, 0, true, dist);
				K_neighbor_search::iterator it = search.begin();
				bool thresholdExceeded = false;

				while (it != search.end() && !thresholdExceeded)
				{					
					if (dist.inverse_of_transformed_distance(it->second) < distanceThreshold)
						neighbors.push_back(it->first);
					else
						thresholdExceeded = true;					
					it++;
				}
			}				
			else				
			{
				//neighbors = SearchNeighbors(pcNeighborhood[k], false);
				Point_3 queryI = uPointList[pcNeighborhood[k]];
				// search K nearest neighbours
				K_neighbor_search search(tree, queryI, kTree, 0, true, dist);
				K_neighbor_search::iterator it = search.begin();
				bool thresholdExceeded = false;

				while (it != search.end() && !thresholdExceeded)
				{					
					if (pcNeighborhood[k] != it->first && dist.inverse_of_transformed_distance(it->second) < distanceThreshold)
						neighbors.push_back(it->first);
					else
						thresholdExceeded = true;
					it++;
				}
			}

			for (int nI = 0; nI < neighbors.size(); nI++)
			{
				bool neighborhoodIncluded = false;
				//Check if the point has not been yet included.
				int nP = 0, nC = 0;
				while (nP < pcNeighborhood.size() && !neighborhoodIncluded)
				{
					if (neighbors[nI] == pcNeighborhood[nP])
						neighborhoodIncluded = true;
					nP++;
				}

				if (!neighborhoodIncluded)
				{
					bool pointClusterIncluded = false;
					int nC = 0;
					while (nC < pointsCluster.size() && !pointClusterIncluded)
					{
						if (pointsCluster[nC] == neighbors[nI])
							pointClusterIncluded = true;
						nC++;
					}

					//int rfPoint = VerifyPointInClusterCr(neighbors[nI], modeCrCenter);
					int K = 2;
					double distanceThreshold = 3.0; // 5.0; // 1.25;
					double rAngleThreshold = 0.25;
					//PointSt pointSt = pointsSampled[ixPoint];
					PointSt pointSt = uPointSampledList[neighbors[nI]];

					//We must find the closests points in the sample, that must be located within some threshold.
					Point_3 reflectedPoint = reflectPointLine(pointSt, modeCrCenter);
					K_neighbor_search search(tree, reflectedPoint, K, 0, true, dist);
					K_neighbor_search::iterator it = search.begin();

					Vector_3 vtLineNormal(modeCrCenter.n.x(), modeCrCenter.n.y(), modeCrCenter.n.z());

					bool found = false;
					int rfPoint = -1;
					//Verify if the point's directions have coherence respect with the converted point.
					while (it != search.end() && !found)
					{
						double distanceSeedCR = dist.inverse_of_transformed_distance(it->second);

						if ((distanceSeedCR < distanceThreshold))
						{
							rfPoint = it->first;
							found = true;
						}
						it++;
					}

					if (!pointClusterIncluded && rfPoint != -1)
					{
						pcNeighborhood.push_back(neighbors[nI]);
						pcNeighborhoodCr.push_back(rfPoint);
					}
				}
			}

			if (originPoint < 4)
				originPoint++;
			else
				k++;
		}

		for (int p = 0; p < pcNeighborhood.size(); p++)
		{
			pointsCluster.push_back(pcNeighborhood[p]);
			pointsClusterCr.push_back(pcNeighborhoodCr[p]);
		}
	}
	std::cout << "         Corrected size: " << pointsCluster.size() << endl;
}

void DetectContinuousRotational(std::vector<ConsRotPair> vConsRotPair)
{
	const clock_t cluster_begin = clock();
	std::vector<ConsRotPair> shifted_pointsCR = meanshiftCR(vConsRotPair, KERNEL_BW);
 	vector<ClusterCR> clustersCR = clusterCR(vConsRotPair, shifted_pointsCR);
	vector<ClusterCR> clustersCRText;
	std::vector<PointClusterFX> patchedClustersCR;
	bool refine = true;

	std::cout << "Number of clusters: " << clustersCR.size() << endl;
	const clock_t cluster_end = clock();
	std::cout << "Time of cluster: " << float(cluster_end - cluster_begin) / CLOCKS_PER_SEC << " seconds." << endl;

	int sampleSize = 0;
	for (int ix = 0; ix < clustersCR.size(); ix++)
		sampleSize += clustersCR[ix].original_points.size();

	for (int ix = 0; ix < clustersCR.size(); ix++)
	{
		int cSize = clustersCR[ix].original_points.size();

		if (cSize > 0.05 * sampleSize)
		{
			ConsRotPair modeCR = clustersCR[ix].original_points[0];
			//std::cout << "Cluster " << ix << " size: " << cSize << " mode: n: " << modeCR.n << " p: " << modeCR.p << endl;

			int idx = -1;
			double distance = 100000.0;

			for (int id = 0; id < cSize; id++)
			{
				double sDistance = 0;

				for (int id_aux = 0; id_aux < cSize; id_aux++)
					sDistance += distanceMetricCR(clustersCR[ix].original_points[id], clustersCR[ix].original_points[id_aux]);

				if (distance > sDistance)
				{
					distance = sDistance;
					idx = id;
				}
			}

			ConsRotPair modeCrCenter = clustersCR[ix].original_points[idx];
			std::cout << "Cluster " << ix << " size: " << cSize << " center: n : " << modeCrCenter.n << " p : " << modeCrCenter.p << endl;

			if (refine)
				RefineContinuousRotational(clustersCR[ix], modeCrCenter);
		}
	}

	if (verbose_cr)
	{
		for (int k = 0; k < shifted_pointsCR.size(); k++)
		{
			for (int m = 0; m < shifted_pointsCR.size(); m++)
			{
				double oDistance = distanceMetricCR(vConsRotPair[k], vConsRotPair[m]);
				double sDistance = distanceMetricCR(shifted_pointsCR[k], shifted_pointsCR[m]);

				if (m < shifted_pointsCR.size() - 1)
				{
					out_verbose_cr << k << "-" << m << "|";
					out_verbose_cr2 << sDistance << "|";
					out_verbose_cr3 << oDistance << "|";
				}
				else
				{
					out_verbose_cr << k << "-" << m;
					out_verbose_cr2 << sDistance;
					out_verbose_cr3 << oDistance;
				}
			}
			out_verbose_cr << std::endl;
			out_verbose_cr2 << std::endl;
			out_verbose_cr3 << std::endl;
		}
	}
}

#pragma endregion

#pragma region Base
std::vector<ReflectionPair> BaseReflection()
{
	std::vector<ReflectionPair> vReflectionPair;
	int sample_size = pointsSampled.size();
	//int sampleSizeM = (int)(sample_size / 2) + 1;

	//for (int i = 0; i < sampleSizeM; i++)
	for (int i = 0; i < sample_size - 1; i++)
	{
		PointSt pI = pointsSampled[i];

		//for (int j = sampleSizeM; j < sample_size; j++)
		for (int j = i + 1; j < sample_size; j++)
		{
			PointSt pJ = pointsSampled[j];
			double ratio_1 = abs(pI.k1 / pJ.k1);
			double ratio_2 = abs(pI.k2 / pJ.k2);

			//uniform scaling, rotation, and translation leave the ratio of principal curvatures unchanged
			if (ratio_1 > RT_CURVATURE && ratio_1 > 1 / RT_CURVATURE &&
				ratio_2 > RT_CURVATURE && ratio_2 > 1 / RT_CURVATURE &&
				ratio_1 / ratio_2 > RT_CURVATURE && ratio_2 / ratio_1 > RT_CURVATURE)
			{
				Vector_3 reflexNormal = pI.p - pJ.p;
				reflexNormal = reflexNormal / sqrt(reflexNormal.squared_length());

				Vector_3 normalBalance = pI.normal + pJ.normal;
				normalBalance = normalBalance / sqrt(normalBalance.squared_length());

				//Evaluate normal inconsistency
				if (abs(normalBalance * reflexNormal) < 1 - RT_CURVATURE)
				{
					if (abs(reflexNormal.x()) > abs(reflexNormal.y()))
					{
						if (abs(reflexNormal.x()) > abs(reflexNormal.z()))
						{
							if (reflexNormal.x() < 0) { reflexNormal = -reflexNormal; }
						}
						else
							if (reflexNormal.z() < 0) { reflexNormal = -reflexNormal; }
					}
					else
					{
						if (abs(reflexNormal.y()) > abs(reflexNormal.z()))
						{
							if (reflexNormal.y() < 0) { reflexNormal = -reflexNormal; }
						}
						else
							if (reflexNormal.z() < 0) { reflexNormal = -reflexNormal; }
					}

					Vector_3 middlePoint = Vector_3(pI.p.x() + pJ.p.x(), pI.p.y() + pJ.p.y(), pI.p.z() + pJ.p.z());
					middlePoint = middlePoint / 2.0;

					ReflectionPair cPlaneReflection;
					cPlaneReflection.i = i;
					cPlaneReflection.j = j;
					cPlaneReflection.PlaneA = reflexNormal.x();
					cPlaneReflection.PlaneB = reflexNormal.y();
					cPlaneReflection.PlaneC = reflexNormal.z();
					cPlaneReflection.PlaneD = -reflexNormal * middlePoint;

					Point_with_normal a;
					a.first = Point_3(middlePoint.x(), middlePoint.y(), middlePoint.z());
					a.second = Kernel::Vector_3(reflexNormal.x(), reflexNormal.y(), reflexNormal.z());
					cPlaneReflection.PointN = a;
					vReflectionPair.push_back(cPlaneReflection);
				}
			}
		}
	}
	return vReflectionPair;
}

std::vector<ConsRotPair> BaseContinuousRotation(std::vector<ReflectionPair> vReflectionPair)
{
	std::vector<ConsRotPair> vConsRotPair;
	for (int m = 0; m < vReflectionPair.size() - 1; m++)
	{
		for (int n = m + 1; n < vReflectionPair.size(); n++)
		{
			//Compare pairs
			ReflectionPair rpM = vReflectionPair[m];
			ReflectionPair rpN = vReflectionPair[n];

			Plane_3 planeM(rpM.PlaneA, rpM.PlaneB, rpM.PlaneC, rpM.PlaneD);
			Plane_3 planeN(rpN.PlaneA, rpN.PlaneB, rpN.PlaneC, rpN.PlaneD);
			CGAL::cpp11::result_of<Intersect_3(Plane_3, Plane_3)>::type result = intersection(planeM, planeN);

			Point_3 origin(0, 0, 0);
			if (const Line_3* s = boost::get<Line_3>(&*result))
			{
				PointSt pI1 = pointsSampled[rpM.i];
				PointSt pJ1 = pointsSampled[rpM.j];
				PointSt pI2 = pointsSampled[rpN.i];
				PointSt pJ2 = pointsSampled[rpN.j];

				if (((pI1.k1 / pI2.k1 > RT_CURVATURE) && (pI2.k1 / pI1.k1 > RT_CURVATURE) && (pI1.k2 / pI2.k2 > RT_CURVATURE) && (pI2.k2 / pI1.k2 > RT_CURVATURE)) ||
					((pI1.k1 / pJ2.k1 > RT_CURVATURE) && (pJ2.k1 / pI1.k1 > RT_CURVATURE) && (pI1.k2 / pJ2.k2 > RT_CURVATURE) && (pJ2.k2 / pI1.k2 > RT_CURVATURE)))
				{
					Point_3 pI1p(pI1.p.x(), pI1.p.y(), pI1.p.z());
					Point_3 pJ1p(pJ1.p.x(), pJ1.p.y(), pJ1.p.z());
					double distance1I = sqrt(CGAL::squared_distance(pI1p, *s));
					//double distance1J = CGAL::squared_distance(pJ1p, *s);
					Point_3 pI2p(pI2.p.x(), pI2.p.y(), pI2.p.z());
					Point_3 pJ2p(pJ2.p.x(), pJ2.p.y(), pJ2.p.z());
					double distance2I = sqrt(CGAL::squared_distance(pI2p, *s));
					//double distance2J = CGAL::squared_distance(pJ2p, *s);
					
					Point_3 proj1 = s->projection(pI1p);
					Point_3 proj2 = s->projection(pI2p);

					if ((abs(distance2I - distance1I) < 1.00) && abs(sqrt(CGAL::squared_distance(proj1, proj2))< 1.00))
					{
						/*
						Vector_3 normalBalance1 = pI1.normal + pJ1.normal;
						normalBalance1 = normalBalance1 / sqrt(normalBalance1.squared_length());
						Point_3 pM(pI1.p.x() / 2 + pJ1.p.x() / 2, pI1.p.y() / 2 + pJ1.p.y() / 2, pI1.p.z() / 2 + pJ1.p.z() / 2);

						Vector_3 normalBalance2 = pI2.normal + pJ2.normal;
						normalBalance2 = normalBalance2 / sqrt(normalBalance2.squared_length());
						Point_3 pN(pI2.p.x() / 2 + pJ2.p.x() / 2, pI2.p.y() / 2 + pJ2.p.y() / 2, pI2.p.z() / 2 + pJ2.p.z() / 2);
						Vector_3 difference1(pI1.p.x() - pJ1.p.x(), pI1.p.y() - pJ1.p.y(), pI1.p.z() - pJ1.p.z());
						Vector_3 difference2(pI2.p.x() - pJ2.p.x(), pI2.p.y() - pJ2.p.y(), pI2.p.z() - pJ2.p.z());
						//Vector_3 product = CGAL::cross_product(difference1, difference2);
						//product = product / sqrt(product.squared_length());

						double angle = abs(normalBalance1 * normalBalance2);
						double distance1 = sqrt(CGAL::squared_distance(pM, *s)); //planeN));
						double distance2 = sqrt(CGAL::squared_distance(pN, *s)); //planeM));
						*/
						//if (abs(distance1 - distance2) < 1.25) //(angle > RT_CURVATURE)
						//{
							Point_3 point = s->projection(origin);
							Kernel::Vector_3 vectorL = s->to_vector();
							vectorL = vectorL / sqrt(vectorL.squared_length());
							vectorL = CorrectDirection(vectorL);

							ConsRotPair crPair;
							crPair.i1 = rpM.i;
							crPair.j1 = rpM.j;
							crPair.i2 = rpN.i;
							crPair.j2 = rpN.j;
							crPair.n = vectorL;
							crPair.p = point;
							vConsRotPair.push_back(crPair);
						//}
					}
				}
			}
		}
	}
	return vConsRotPair;
}

std::vector<RotTransPair> BaseDiscreteRotation()
{
	std::vector<RotTransPair> vRotTransPair;
	int sample_size = pointsSampled.size();
	int sampleSizeM = sample_size / 2;

	//Find rotational-traslational symmetry
	for (int i = 0; i < sampleSizeM; i++)
	{
		PointSt pI = pointsSampled[i];
		Eigen::Vector3d pI_p(pI.p.x(), pI.p.y(), pI.p.z());

		for (int j = sampleSizeM; j < sample_size; j++)
		{
			PointSt pJ = pointsSampled[j];
			Eigen::Vector3d pJ_p(pJ.p.x(), pJ.p.y(), pJ.p.z());

			double ratio_1 = abs(pI.k1 / pJ.k1);
			double ratio_2 = abs(pI.k2 / pJ.k2);

			//uniform scaling, rotation, and translation leave the ratio of principal curvatures unchanged
			if (ratio_1 / ratio_2 > RT_CURVATURE && ratio_2 / ratio_1 > RT_CURVATURE &&
				ratio_1 > RT_CURVATURE &&  ratio_2 > RT_CURVATURE && ratio_1 < 1 / RT_CURVATURE &&  ratio_2 < 1 / RT_CURVATURE)
			{
				double cAngle1 = pI.kmin * pI.kmax; //0
				double cAngle2 = pJ.kmin * pJ.kmax;//0

				Vector_3 rAxis = CGAL::cross_product(pI.normal, pJ.normal);
				Eigen::Matrix3d mAlign = alignmentMatrix(pI.normal, pJ.normal, rAxis);

				Eigen::Vector3d pIn(pI.normal.x(), pI.normal.y(), pI.normal.z());
				Eigen::Vector3d pJn(pJ.normal.x(), pJ.normal.y(), pJ.normal.z());

				//Validations
				//Curvatures * normals
				//double cCNiMin = pI.normal * pI.kmin;//0
				//double cCNiMax = pI.normal * pI.kmax;//0
				//double cCNjMin = pJ.normal * pJ.kmin;//0
				//double cCNjMax = pJ.normal * pJ.kmax;//0
				Eigen::Vector3d c1I(pI.kmin.x(), pI.kmin.y(), pI.kmin.z());
				Eigen::Vector3d c2I(pI.kmax.x(), pI.kmax.y(), pI.kmax.z());
				Eigen::Vector3d c1J(pJ.kmin.x(), pJ.kmin.y(), pJ.kmin.z());
				Eigen::Vector3d c2J(pJ.kmax.x(), pJ.kmax.y(), pJ.kmax.z());

				Eigen::Vector3d c1Ia = mAlign * c1I;
				Eigen::Vector3d c2Ia = mAlign * c2I;
				Vector_3 c1IAlign = Vector_3(c1Ia.x(), c1Ia.y(), c1Ia.z());
				Vector_3 c2IAlign = Vector_3(c2Ia.x(), c2Ia.y(), c2Ia.z());

				Eigen::Vector3d xpJn = mAlign * pIn;
				Vector_3 vpJn(xpJn.x(), xpJn.y(), xpJn.z());

				//After first rotation:
				//double cC1 = pJ.normal * c1IAlign;
				//double cC2 = pJ.normal * c2IAlign;
				//double cC3 = c1IAlign * c2IAlign;
				//double cC01 = vpJn * c1IAlign;
				//double cC02 = vpJn * c2IAlign;
				//Vector_3 vAxis = CGAL::cross_product(c1IAlign, c2IAlign);

				Eigen::Matrix3d mAlignC1;
				Vector_3 jKmin = pJ.kmin;
				Vector_3 jKmax = pJ.kmax;
				//double cC01x = vpJn * jKmin;
				Vector_3 rAxisN = CGAL::cross_product(c1IAlign, jKmin);
				//double cNP = rAxisN * pJ.normal;
				//double cAngleC = c1IAlign * pJ.kmin;

				//if (acos(cAngleC) < acos(-cAngleC))
				mAlignC1 = alignmentMatrix(c1IAlign, jKmin, pJ.normal); // pJ.normal);
				//else
				//	mAlignC1 = alignmentMatrix(-c1IAlign, jKmin, rAxisN); // pJ.normal);

				Eigen::Vector3d pJKnMin = mAlignC1 * c1Ia;
				Eigen::Vector3d pJKnMax = mAlignC1 * c2Ia;
				//if (acos(c1IAlign * pJ.kmin) >= acos(-c1IAlign * pJ.kmin))
				//	pJKn = -pJKn;

				Vector_3 cpJKMin = Vector_3(pJKnMin.x(), pJKnMin.y(), pJKnMin.z());
				Vector_3 cpJKMax = Vector_3(pJKnMax.x(), pJKnMax.y(), pJKnMax.z());

				double productMin = cpJKMin * jKmin;
				double productMax = cpJKMax * jKmax;
				if (productMin < 0.99)
					int val = -1;

				Eigen::Matrix3d mAlignFinal = mAlignC1 * mAlign;
				double scale_ij = (ratio_1 + ratio_2) / 2;
				Eigen::Vector3d mTraslation = pJ_p - mAlignFinal * pI_p;

				/*
				//Reference: Rotation Matrix to Euler angles (ZYX). web.mit.edu/2.05/www/Handout/HO2.PDF
				double Ry = M_PI / 2;
				double Rz = 0;
				double Rx = atan2(mAlignFinal(0, 1), mAlignFinal(1, 1));
				if (mAlignFinal(0, 0) != 0.0 || mAlignFinal(1, 0) != 0.0)
				{
				Ry = atan2(-mAlignFinal(2, 0), sqrt(pow(mAlignFinal(0, 0), 2) + pow(mAlignFinal(1, 0), 2)));
				Rz = atan2(mAlignFinal(1, 0), mAlignFinal(0, 0));
				Rx = atan2(mAlignFinal(2, 1), mAlignFinal(2, 2));
				}
				*/

				double Rx = 0;
				double Ry = 0;
				double Rz = 0;
				//Reference: www.geometrictools.com/Documentation/EulerAngles.pdf
				double sinThreshold = 0.999;

				if (mAlignFinal(1, 2) < sinThreshold)
				{
					if (mAlignFinal(1, 2) > -sinThreshold)
					{
						Rx = asin(-mAlignFinal(1, 2));
						Ry = atan2(mAlignFinal(0, 2), mAlignFinal(2, 2));
						Rz = atan2(mAlignFinal(1, 0), mAlignFinal(1, 1));
					}
					else
					{
						Rx = M_PI / 2;
						Ry = -atan2(-mAlignFinal(0, 1), mAlignFinal(0, 0));
						Rz = 0;
					}
				}
				else
				{
					Rx = -M_PI / 2;
					Ry = atan2(-mAlignFinal(0, 1), mAlignFinal(0, 0));
					Rz = 0;
				}

				Eigen::Matrix3d rotationX;
				rotationX(0, 0) = 1;	rotationX(0, 1) = 0;		rotationX(0, 2) = 0;
				rotationX(1, 0) = 0;	rotationX(1, 1) = cosf(Rx);	rotationX(1, 2) = -sinf(Rx);
				rotationX(2, 0) = 0;	rotationX(2, 1) = sinf(Rx);	rotationX(2, 2) = cosf(Rx);

				Eigen::Matrix3d rotationY;
				rotationY(0, 0) = cosf(Ry);		rotationY(0, 1) = 0;	rotationY(0, 2) = sinf(Ry);
				rotationY(1, 0) = 0;			rotationY(1, 1) = 1;	rotationY(1, 2) = 0;
				rotationY(2, 0) = -sinf(Ry);	rotationY(2, 1) = 0;	rotationY(2, 2) = cosf(Ry);

				Eigen::Matrix3d rotationZ;
				rotationZ(0, 0) = cosf(Rz);	rotationZ(0, 1) = -sinf(Rz);	rotationZ(0, 2) = 0;
				rotationZ(1, 0) = sinf(Rz);	rotationZ(1, 1) = cosf(Rz);		rotationZ(1, 2) = 0;
				rotationZ(2, 0) = 0;		rotationZ(2, 1) = 0;			rotationZ(2, 2) = 1;

				Eigen::Matrix3d vTransformedPoint_01 = rotationY * rotationX * rotationZ; //rotationZ * rotationY * rotationX;

				RotTransPair rtsTransformation;
				rtsTransformation.i = i;
				rtsTransformation.j = j;
				rtsTransformation.Scale = 1;
				rtsTransformation.Rx = Rx;
				rtsTransformation.Ry = Ry;
				rtsTransformation.Rz = Rz;
				rtsTransformation.Tx = mTraslation.x();
				rtsTransformation.Ty = mTraslation.y();
				rtsTransformation.Tz = mTraslation.z();

				vRotTransPair.push_back(rtsTransformation);
			}
		}
	}
	return vRotTransPair;
}

#pragma endregion

///////////////MAIN///////////////////////////////////////////////////////
#if defined(CGAL_USE_BOOST_PROGRAM_OPTIONS) && ! defined(DONT_USE_BOOST_PROGRAM_OPTIONS)
int main(int argc, char *argv[])
#else
int main()
#endif
{
  string if_name_string;
  string if_name; //input file name
  string w_if_name;  //as above, but / replaced by _
  string res4openGL_fname;
  string verbose_fname;
  string verbose_fname2;
  string verbose_fname3;
  std::ofstream out_4ogl; //out_verbose, out_verbose2;

	try 
	{
		#if defined(CGAL_USE_BOOST_PROGRAM_OPTIONS) && ! defined(DONT_USE_BOOST_PROGRAM_OPTIONS)
			po::options_description desc("Allowed options");
			//C:\dev\CGAL - 4.9\build - example\Jet_fitting_3\Debug\data
			desc.add_options()
			  ("help,h", "produce help message.")
			  ("input-file,f", po::value<string>(&if_name_string)->default_value(filePath),
			   "name of the input off file")
			  ("degree-jet,d", po::value<unsigned int>(&d_fitting)->default_value(4),
			   "degree of the jet, 1 <= degre-jet <= 4")
			  ("degree-monge,m", po::value<unsigned int>(&d_monge)->default_value(4),
			   "degree of the Monge rep, 1 <= degree-monge <= degree-jet")
			  ("nb-rings,a", po::value<unsigned int>(&nb_rings)->default_value(0),
			   "number of rings to collect neighbors. 0 means collect enough rings to make appro possible a>=1 fixes the nb of rings to be collected")
			  ("nb-points,p", po::value<unsigned int>(&nb_points_to_use)->default_value(0),
			   "number of neighbors to use.  0 means this option is not considered, this is the default p>=1 fixes the nb of points to be used")
			  ("verbose,v", po::value<bool>(&verbose)->default_value(true),
			   "verbose output on text file")
			  ;

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);

			if (vm.count("help")) {
			  cout << desc << "\n";
			  return 1;
			}
		#else
			std::cerr << "Command-line options require Boost.ProgramOptions" << std::endl;
			if_name_string = "Ark_HM_10_lo-0.off";
			d_fitting = 2;
			d_monge = 2;
			nb_rings = 0;
			nb_points_to_use = 0;
			verbose = false;
		#endif
	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}
	catch(...) {
		cerr << "Exception of unknown type!\n";
	}

	const clock_t begin_time = clock();

	//modify global variables which are fct of options:
	min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;
	if (nb_points_to_use < min_nb_points && nb_points_to_use != 0)
	{
		std::cerr << "the nb of points asked is not enough to perform the fitting" << std::endl; exit(0);
	}

	//prepare output file names
	//--------------------------
	std::cerr << "if_name_string" << if_name_string  << std::endl;
	if_name = if_name_string;

	w_if_name = if_name;
	for(unsigned int i=0; i<w_if_name.size(); i++)
		if (w_if_name[i] == '/') w_if_name[i]='_';
	//cerr << if_name << '\n';
	//cerr << w_if_name << '\n';

	res4openGL_fname = w_if_name + ".4ogl.txt";
	//std::cerr << "res4openGL_fname" << res4openGL_fname  << std::endl;
	out_4ogl.open(res4openGL_fname.c_str(), std::ios::out);
	assert(out_4ogl.good());
	//if verbose only...
	if(verbose){
		//Rotation
		verbose_fname  = w_if_name + "_rt.txt";
		out_verbose_rt.open(verbose_fname.c_str(), std::ios::out);
		assert(out_verbose_rt.good());
		CGAL::set_pretty_mode(out_verbose_rt);

		verbose_fname2 = w_if_name + "_rt2.txt";
		out_verbose_rt2.open(verbose_fname2.c_str(), std::ios::out);
		assert(out_verbose_rt2.good());
		CGAL::set_pretty_mode(out_verbose_rt2);

		verbose_fname3 = w_if_name + "_rt3.txt";
		out_verbose_rt3.open(verbose_fname3.c_str(), std::ios::out);
		assert(out_verbose_rt3.good());
		CGAL::set_pretty_mode(out_verbose_rt3);
		
		//Discrete Reflection
		verbose_fname = w_if_name + "_fx.txt";
		out_verbose_fx.open(verbose_fname.c_str(), std::ios::out);
		assert(out_verbose_fx.good());
		CGAL::set_pretty_mode(out_verbose_fx);

		verbose_fname2 = w_if_name + "_fx2.txt";
		out_verbose_fx2.open(verbose_fname2.c_str(), std::ios::out);
		assert(out_verbose_fx2.good());
		CGAL::set_pretty_mode(out_verbose_fx2);

		verbose_fname3 = w_if_name + "_fx3.txt";
		out_verbose_fx3.open(verbose_fname3.c_str(), std::ios::out);
		assert(out_verbose_fx3.good());
		CGAL::set_pretty_mode(out_verbose_fx3);

		///Continuous Rotation
		verbose_fname = w_if_name + "_cr.txt";
		out_verbose_cr.open(verbose_fname.c_str(), std::ios::out);
		assert(out_verbose_cr.good());
		CGAL::set_pretty_mode(out_verbose_cr);

		verbose_fname2 = w_if_name + "_cr2.txt";
		out_verbose_cr2.open(verbose_fname2.c_str(), std::ios::out);
		assert(out_verbose_cr2.good());
		CGAL::set_pretty_mode(out_verbose_cr2);

		verbose_fname3 = w_if_name + "_cr3.txt";
		out_verbose_cr3.open(verbose_fname3.c_str(), std::ios::out);
		assert(out_verbose_cr3.good());
		CGAL::set_pretty_mode(out_verbose_cr3);
	}
	unsigned int nb_vertices_considered = 0;//count vertices for verbose

	//load the model from <mesh.off>
	//------------------------------
	PolyhedralSurf P;
	std::ifstream stream(if_name.c_str());
	stream >> P;
	std::cout << "loadMesh...  "<< "Polysurf with " << P.size_of_vertices()
		<< " vertices and " << P.size_of_facets()
		<< " facets. " << std::endl;
	
	//exit if not enough points in the model
	//if (min_nb_points > P.size_of_vertices())    exit(0);
		
	const clock_t read_time = clock();
	std::cout << "Time of reading: " << float(read_time - begin_time) / CLOCKS_PER_SEC << " seconds." << endl;

	CGAL::Random random = CGAL::Random::Random();
	int total_size = num_vertices(P);
	//	int total_size = P.size_of_vertices();
	
	//create property maps
	//-----------------------------
	//Vertex, using a std::map
	Vertex2int_map_type vertex2props;
	Vertex_PM_type vpm(vertex2props);

	//Hedge, with enriched hedge
	//HEdgePM_type hepm = get_hepm(boost::edge_weight_t(), P);
	//Hedge, using a std::map
	Hedge2double_map_type hedge2props;
	//Hedge_PM_type hepm(hedge2props); //CEPS

	//Facet PM, with enriched Facet
	//FacetPM_type fpm = get_fpm(boost::vertex_attribute_t(), P);
	//Facet PM, with std::map
	Facet2normal_map_type facet2props;
	Facet_PM_type fpm(facet2props);

	//initialize Polyhedral data : length of edges, normal of facets
	//Poly_hedge_ops::compute_edges_length(P, hepm); //CEPS
	Poly_facet_ops::compute_facets_normals(P, fpm);

	CGAL::Simple_cartesian<DFT>::Iso_cuboid_3 bc3 = CGAL::bounding_box(P.points_begin(), P.points_end());
	diagonal_bbox = 0.5 * sqrt(squared_distance(bc3.min(), bc3.max()));

	double length_x = bc3.max().x() - bc3.min().x();
	double length_y = bc3.max().y() - bc3.min().y();
	double length_z = bc3.max().z() - bc3.min().z();
	double area = 0.125 * (length_x * length_y + length_x * length_z + length_y * length_z);

	if (detectReflectional || detectRotational)
	{
		//MAIN LOOP: perform calculation for each vertex
		//----------------------------------------------
		std::vector<DPoint> in_points;  //container for data points
		Vertex_iterator vitb, vite;

		//initialize the tag of all vertices to -1
		vitb = P.vertices_begin(); vite = P.vertices_end();
		CGAL_For_all(vitb, vite) put(vpm, &(*vitb), -1);
		vitb = P.vertices_begin(); vite = P.vertices_end();

		std::vector<ReflectionPair> vReflectionPair;
		std::vector<RotTransPair> vRotTransPair;
		std::vector<CGAL::Simple_cartesian<DFT>::Point_3> points_3;

		//const int PRESIZE = 200;
		//const int SAMPLESIZE = 2;

		//First sampling
		int preSample_size = 0.64 * area; //pow(diagonal_bbox, 2);// (int)(sqrt(total_size) * 2);
		vector<int> preSample_indices;
		preSample_indices.resize(preSample_size);
		for (int i = 0; i < preSample_size; i++)
		{
			int cociente = (int)total_size / preSample_size;
			int index = random.uniform_int(cociente * i, min(cociente * (i + 1) - 1, total_size - 1));
			preSample_indices[i] = index;
		}

		uPointList.resize(preSample_size);
		uPointSampledList.resize(preSample_size);
		points_3.resize(preSample_size);

		int indice = 0;
		int preSampleIndex = 0;
		int secondSampleIndex = 0;

		//Determining the size of the sample
		int sample_size = 0.09 * area; //pow(diagonal_bbox, 2);
		vector<int> sample_indices;
		vector<int> validation_indices;

		if (sample_size > preSample_size)
			sample_size = preSample_size;

		sample_indices.resize(sample_size);
		validation_indices.resize(sample_size);

		for (int i = 0; i < sample_size; i++)
		{
			int cociente = (int)preSample_size / sample_size;
			int index = random.uniform_int(cociente * i, min(cociente * (i + 1) - 1, total_size - 1));
			sample_indices[i] = preSample_indices[index]; //index
			validation_indices[i] = index;
		}

		for (; vitb != vite; vitb++)
		{
			//initialize
			Vertex* v = &(*vitb);

			//PreSample
			if (preSampleIndex < preSample_size && indice == preSample_indices[preSampleIndex])
			{
				DPoint pPointSampled;
				pPointSampled = v->point();				
				points_3[preSampleIndex] = v->point();

				Point_3 uPoint(pPointSampled.x(), pPointSampled.y(), pPointSampled.z());
				uPointList[preSampleIndex] = uPoint;

				PointSt pointSampled;
				pointSampled.p = pPointSampled;
				uPointSampledList[preSampleIndex] = pointSampled;

				preSampleIndex++;
			}

			in_points.clear();
			My_Monge_form monge_form;
			//PreSample
			if (secondSampleIndex < sample_size && indice == sample_indices[secondSampleIndex])
			{
				//gather points around the vertex using rings
				gather_fitting_points(v, in_points, vpm);

				//skip if the nb of points is to small
				if (in_points.size() < min_nb_points)
				{
					std::cerr << "not enough pts for fitting this vertex" << in_points.size() << std::endl;
					continue;
				}

				// perform the fitting
				My_Monge_via_jet_fitting monge_fit;
				monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);
				//switch min-max ppal curv/dir wrt the mesh orientation
				const DVector normal_mesh = Poly_facet_ops::compute_vertex_average_unit_normal(v, fpm);
				monge_form.comply_wrt_given_normal(normal_mesh);

				PointSt pointSampled;
				pointSampled.p = v->point();
				pointSampled.kmin = monge_form.minimal_principal_direction();
				pointSampled.kmax = monge_form.maximal_principal_direction();
				//Vector_3 normal = monge_form.normal_direction();
				Vector_3 normalVector = CGAL::cross_product(pointSampled.kmin, pointSampled.kmax);
				pointSampled.normal = normalVector / sqrt(normalVector.squared_length());
				pointSampled.k1 = monge_form.principal_curvatures(1);
				pointSampled.k2 = monge_form.principal_curvatures(0);
				points_3[secondSampleIndex] = v->point();

				if (abs(pointSampled.k1) + abs(pointSampled.k2) > 0.05)
				{
					/*
					Point_3 uPoint(pointSampled.p.x(), pointSampled.p.y(), pointSampled.p.z());
					uPointList[preSampleIndex] = uPoint;
					uPointSampledList[preSampleIndex] = pointSampled;
					*/
					pointsSampled.push_back(pointSampled); //[sampleIndex]
					Point_3 point_3d(pointSampled.p.x(), pointSampled.p.y(), pointSampled.p.z());
					points.push_back(point_3d);
					uIndices.push_back(validation_indices[secondSampleIndex]); //uIndices.push_back(indice);
				}

				secondSampleIndex++;
			}
			indice++;
		} //all vertices processed
		
		//CGAL::Simple_cartesian<DFT>::Iso_cuboid_3 bc3 = CGAL::bounding_box(points_3.begin(), points_3.end());
		//diagonal_bbox = 0.5 * sqrt(squared_distance(bc3.min(), bc3.max()));
		Point_3 center(bc3.min().x() / 2 + bc3.max().x() / 2, bc3.min().y() / 2 + bc3.max().y() / 2, bc3.min().z() / 2 + bc3.max().z() / 2);

		std::cout << "        Center: (" << center.x() << ", " << center.y() << ", " << center.z() << ")" << endl;
		const clock_t preSample_time = clock();
		std::cout << "Time of presample: " << float(preSample_time - read_time) / CLOCKS_PER_SEC << " seconds." << endl;
		std::cout << "Presample: " << preSample_size << endl;

		/*
		for (; vitb != vite; vitb++)
		{
			//initialize
			Vertex* v = &(*vitb);
			in_points.clear();
			My_Monge_form monge_form;

			//PreSample
			if (preSampleIndex < preSample_size && indice == preSample_indices[preSampleIndex])
			{
				//gather points around the vertex using rings
				gather_fitting_points(v, in_points, vpm);

				//skip if the nb of points is to small
				if (in_points.size() < min_nb_points)
				{
					std::cerr << "not enough pts for fitting this vertex" << in_points.size() << std::endl;
					continue;
				}

				// perform the fitting
				My_Monge_via_jet_fitting monge_fit;
				monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);
				//switch min-max ppal curv/dir wrt the mesh orientation
				const DVector normal_mesh = Poly_facet_ops::compute_vertex_average_unit_normal(v, fpm);
				monge_form.comply_wrt_given_normal(normal_mesh);

				PointSt pointSampled;
				pointSampled.p = v->point();
				pointSampled.kmin = monge_form.minimal_principal_direction();
				pointSampled.kmax = monge_form.maximal_principal_direction();
				//Vector_3 normal = monge_form.normal_direction();
				Vector_3 normalVector = CGAL::cross_product(pointSampled.kmin, pointSampled.kmax);
				pointSampled.normal = normalVector / sqrt(normalVector.squared_length());
				pointSampled.k1 = monge_form.principal_curvatures(1);
				pointSampled.k2 = monge_form.principal_curvatures(0);
				points_3[preSampleIndex] = v->point();

				Point_3 uPoint(pointSampled.p.x(), pointSampled.p.y(), pointSampled.p.z());
				uPointList[preSampleIndex] = uPoint;
				uPointSampledList[preSampleIndex] = pointSampled;
				preSampleIndex++;
			}
			
			indice++;
		} //all vertices processed
		*/
				
		//vitb = P.vertices_begin(); vite = P.vertices_end();
		/*
		for (int sampleIndex = 0; sampleIndex < sample_size; sampleIndex++)
		{
			indice = sample_indices[sampleIndex];
			PointSt pointSampled = uPointSampledList[indice];

			//if (abs(pointSampled.k1 / pointSampled.k2) < MM_CURVATURE_THRESHOLD || abs(pointSampled.k2 / pointSampled.k1) < MM_CURVATURE_THRESHOLD)
			if (abs(pointSampled.k1) + abs(pointSampled.k2) > 0.05)
			{
				pointsSampled.push_back(pointSampled); //[sampleIndex]

				Point_3 point_3d(pointSampled.p.x(), pointSampled.p.y(), pointSampled.p.z());
				points.push_back(point_3d);
				uIndices.push_back(indice);
			}
		}
		*/
		/*
		preSampleIndex = 0;
		indice = 0;

		for (; vitb != vite; vitb++)
		{
			//initialize
			Vertex* v = &(*vitb);
			in_points.clear();
			My_Monge_form monge_form;

			//PreSample
			if (preSampleIndex < sample_size && indice == sample_indices[preSampleIndex])
			{
				//gather points around the vertex using rings
				gather_fitting_points(v, in_points, vpm);

				//skip if the nb of points is to small
				if (in_points.size() < min_nb_points)
				{
					std::cerr << "not enough pts for fitting this vertex" << in_points.size() << std::endl;
					continue;
				}

				// perform the fitting
				My_Monge_via_jet_fitting monge_fit;
				monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);
				//switch min-max ppal curv/dir wrt the mesh orientation
				const DVector normal_mesh = Poly_facet_ops::compute_vertex_average_unit_normal(v, fpm);
				monge_form.comply_wrt_given_normal(normal_mesh);

				PointSt pointSampled;
				pointSampled.p = v->point();
				pointSampled.kmin = monge_form.minimal_principal_direction();
				pointSampled.kmax = monge_form.maximal_principal_direction();
				//Vector_3 normal = monge_form.normal_direction();
				Vector_3 normalVector = CGAL::cross_product(pointSampled.kmin, pointSampled.kmax);
				pointSampled.normal = normalVector / sqrt(normalVector.squared_length());
				pointSampled.k1 = monge_form.principal_curvatures(1);
				pointSampled.k2 = monge_form.principal_curvatures(0);
				points_3[preSampleIndex] = v->point();

				if (abs(pointSampled.k1) + abs(pointSampled.k2) > 0.05)
				{
					//Point_3 uPoint(pointSampled.p.x(), pointSampled.p.y(), pointSampled.p.z());
					//uPointList[preSampleIndex] = uPoint;
					//uPointSampledList[preSampleIndex] = pointSampled;
					pointsSampled.push_back(pointSampled); //[sampleIndex]
					Point_3 point_3d(pointSampled.p.x(), pointSampled.p.y(), pointSampled.p.z());
					points.push_back(point_3d);
					uIndices.push_back(validation_indices[preSampleIndex]); //uIndices.push_back(indice);
				}

				preSampleIndex++;
			}

			indice++;
		} //all vertices processed
		*/

		/*
		for (int m = 0; m < pointsSampled.size(); m++)
		{
		if (verbose)
		out_verbose2 << "Point " << m << ": p(" <<
		pointsSampled[m].p.x() << ", " << pointsSampled[m].p.y() << ", " << pointsSampled[m].p.z() << ")" <<
		"n: (" << pointsSampled[m].normal.x() << ", " << pointsSampled[m].normal.y() << ", " << pointsSampled[m].normal.z() << ") " <<
		"c1: (" << pointsSampled[m].kmin.x() << ", " << pointsSampled[m].kmin.y() << ", " << pointsSampled[m].kmin.z() << ") " <<
		"c2: (" << pointsSampled[m].kmax.x() << ", " << pointsSampled[m].kmax.y() << ", " << pointsSampled[m].kmax.z() << ") " <<
		std::endl;
		}
		*/

		const clock_t curvature_time = clock();
		std::cout << "Time of curvature: " << float(curvature_time - preSample_time) / CLOCKS_PER_SEC << " seconds." << endl;
		std::cout << "Sample: " << sample_size << " filtered: " << uIndices.size() << endl;

		//Find reflectional symmetry
		sample_size = pointsSampled.size();
		int sampleSizeM = sample_size / 2;
		std::vector<ConsRotPair> vConsRotPair;

		if (detectReflectional)
		{
			vReflectionPair = BaseReflection();
			vConsRotPair = BaseContinuousRotation(vReflectionPair);
		}

		if (detectRotational)
			vRotTransPair = BaseDiscreteRotation();

		const unsigned int K = 9;
		PointPropertyMap ppmap(points);
		// Insert number_of_data_points in the tree
		Tree tree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(points.size()), Splitter(), Traits(ppmap));
		Distance tr_dist(ppmap);

		//vector<Cluster> clusters = clusterReflect(vRotTransPair, kernel_bandwidth);	
		//vector<Cluster> clusters = clusterRT(vRotTransPair, kernel_bandwidth);

		std::cout << std::setprecision(3) << std::fixed;

		if (detectRotational)
		{
			std::cout << "Discrete Rotations: " << vRotTransPair.size() << " samples" << endl;
			DetectionDiscreteRotational(vRotTransPair);
		}

		if (detectReflectional)
		{
			std::cout << "Reflection: " << vReflectionPair.size() << " samples" << endl;
			DetectionReflectional(vReflectionPair);

			std::cout << "Continuous Rotation: " << vConsRotPair.size() << " samples" << endl;
			DetectContinuousRotational(vConsRotPair);

			/*
			for (int ix = 0; ix < clustersFx.size(); ix++)
			{
			if ((clustersFx[ix].original_points.size() > 3) && (refine))
			{
			//Addition points through ICP
			ReflectionPair tReflection = clustersFx[ix].transformation;

			for (int s = 0; s < clustersFx[ix].original_points.size(); s++)
			{
			//Add the points included in the cluster from the seed.
			ReflectionPair rPair = clustersFx[ix].original_points[s];
			std::vector<int> pcNeighborhood;
			std::vector<int> pcNeighborhoodRx;
			pcNeighborhood.push_back(rPair.i);
			pcNeighborhoodRx.push_back(rPair.j);

			for (int k = 0; k < pcNeighborhood.size(); k++)
			{
			PointSt pointSt = uPointSampledList[pcNeighborhood[k]];
			//We must find the closests points in the sample, that must be located within some threshold.
			Point_3 reflectedPoint = reflectPointPlane(pointSt, tReflection);

			DPoint rxPoint = uPointSampledList[pcNeighborhoodRx[k]].p;
			Point_3 transformedPoint(rxPoint.x(), rxPoint.y(), rxPoint.z());

			double distance = sqrt(CGAL::squared_distance(reflectedPoint, transformedPoint));
			}
			}
			}
			}
			*/
		}

		const clock_t sampling_time = clock();
		std::cout << "Time of computing: " << float(sampling_time - curvature_time) / CLOCKS_PER_SEC << " seconds." << endl;
	}

	bool refined = false;
	
	if (refined)
	{
		std::vector<DPoint> pts;                    // smooth the old vertices
		pts.reserve(P.size_of_vertices());  // get intermediate space for the new points
		std::transform(P.vertices_begin(), P.vertices_end(), std::back_inserter(pts), ReflectionTransform());
		std::copy(pts.begin(), pts.end(), P.points_begin());

		std::ofstream refined_off("refined.off");
		refined_off << P;
		refined_off.close();
		std::cout << "Saved file." << endl;
	}

	//cleanup filenames
	//------------------

	out_4ogl.close();

	if (verbose) {
		out_verbose_rt.close();
		out_verbose_rt2.close();

		out_verbose_fx.close();
		out_verbose_fx2.close();

		out_verbose_cr.close();
		out_verbose_cr2.close();
	}
	return 0;
}