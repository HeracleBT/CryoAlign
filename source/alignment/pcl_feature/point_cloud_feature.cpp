#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/keypoints/iss_3d.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/3dsc.h>
#include <pcl/search/kdtree.h>
#include <pcl/console/time.h>
#include <pcl/filters/random_sample.h>
#include <pcl/registration/ia_ransac.h>
//#include <boost/thread/thread.hpp>
#include <time.h> // 时间
#include <pcl/common/transforms.h>
#include <pcl/features/fpfh.h>
#include <pcl/features/pfh.h>
#include <pcl/features/shot.h>
#include <pcl/features/rops_estimation.h>
#include <pcl/surface/gp3.h>
#include <pcl/features/usc.h>
// #include <pcl/features/shot_omp.h>

#include <pcl/registration/ia_fpcs.h> // 4PCS算法


#include <pcl/filters/filter.h>


#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <pcl/keypoints/sift_keypoint.h>   // SIFT关键点相关
#include <pcl/pcl_macros.h>
#include <pcl/registration/correspondence_estimation.h>
#include <pcl/registration/correspondence_rejection_sample_consensus.h>
#include <pcl/registration/correspondence_rejection_distance.h>
#include <fstream>

#include <pcl/common/common_headers.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace std;

#define MAXBUFSIZE  ((int) 1e6)

typedef pcl::PointXYZ PointT;
typedef pcl::Normal PointNT;
typedef pcl::SHOT352 FeatureT;

// compute SHOT features
void compute_shot_keypoint(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::Normal>::Ptr cloud_normal, pcl::PointCloud<pcl::PointXYZ>::Ptr keypoints, float radius)
{   
    pcl::SHOTEstimation<pcl::PointXYZ, pcl::Normal, pcl::SHOT352> shot;
    shot.setRadiusSearch(radius);
    shot.setInputCloud(keypoints);
    shot.setSearchSurface(cloud);
    shot.setInputNormals(cloud_normal);
    pcl::PointCloud<FeatureT>::Ptr features(new pcl::PointCloud<FeatureT>);
    shot.compute(*features);

	for (int i = 0; i < features->size(); i++) 
	{
		pcl::SHOT352 descriptor = features->points[i];
		cout << descriptor << endl;
	}

    // ofstream write;
	// write.open(output_dir);
	// for (int i = 0; i < features->size(); i++) 
	// {
	// 	pcl::SHOT352 descriptor = features->points[i];
	// 	write << descriptor;
	// 	write << "\n";
	// }
	// write.close();

}

Eigen::MatrixXd readMatrix(const char *filename)
    {
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
        {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
    };


int main(int argc, char** argv)
{

	pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
    pcl::io::loadPCDFile(argv[1], *cloud);

    // cout << "read pcd" << endl;

    pcl::PointCloud<pcl::Normal>::Ptr cloud_normal(new pcl::PointCloud<pcl::Normal>);
    cloud_normal->width = cloud->width;
    cloud_normal->height = cloud->height;
    cloud_normal->is_dense = true;
    cloud_normal->points.resize(cloud_normal->width * cloud_normal->height);

    Eigen::MatrixXd normals = readMatrix(argv[2]);
    for (int i = 0; i < cloud->points.size(); i++)
    {
        cloud_normal->points[i].normal_x = normals(i, 0);
        cloud_normal->points[i].normal_y = normals(i, 1);
        cloud_normal->points[i].normal_z = normals(i, 2);
    }
    pcl::PointCloud<PointT>::Ptr keypoints(new pcl::PointCloud<PointT>);
    pcl::io::loadPCDFile(argv[3], *keypoints);
    float radius = atof(argv[4]);
    compute_shot_keypoint(cloud, cloud_normal, keypoints, radius);

    return 0;
}

