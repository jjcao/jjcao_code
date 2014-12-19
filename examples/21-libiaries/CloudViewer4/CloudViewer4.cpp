/*=============================================================================
 * CloudViewer4.cpp
 * A straight forward, simple point cloud visualization.
 * This function views a labeled point cloud with user-defined colors.
 * Syntax:
 *   [] = CloudViewer4(input_cloud, labels, color, param)
 * Where
 *   This function has no output parameters.
 *   input_cloud: The input point cloud represented by an N x 3 matrix.
 *   labels: The label of the points represented by an N x 1 vector starting from 1. 
 *   color: The color of the regions represented by an M x 3 matrix.
 *   param: The rendering parameters represented by a 1 x 4 vector.
 *     param[1]: The size of rendered points.
 *     param[2], param[3], param[4]: The RGB values of the background.
 *
 * This is a MEX-file for MATLAB.
 * Baochang Han hanbc@mail.dlut.edu.cn
 * 2013-03-05
 *=============================================================================*/


// MEX includes
#include <mex.h>

// PCL includes
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>

// Additional dependencies
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libmex.lib")
#pragma comment(lib, "libmat.lib")

#ifdef _DEBUG
#pragma comment(lib, "pcl_visualization_debug.lib")
#pragma comment(lib, "pcl_common_debug.lib")
#else
#pragma comment(lib, "pcl_visualization_release.lib")
#pragma comment(lib, "pcl_common_release.lib")
#endif // _DEBUG

// Global varibles for viewer parameters
static double point_size = 2.0;
static double background[3] = {1.0, 1.0, 1.0};

// Set the parameters
static void viewerSetting(pcl::visualization::PCLVisualizer& viewer)
{
	// Set the background color
	viewer.setBackgroundColor(background[0], background[1], background[2]);

	// Set the point size
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, point_size);

	// Remove the coordinate system
	viewer.removeCoordinateSystem();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	// Get the data from prhs
	if (nrhs != 3 && nrhs != 4)
	{
		mexErrMsgTxt("Three or four input arguments required.");
	}

	if (mxGetN(prhs[0]) != 3)
	{
		mexErrMsgTxt("The input point cloud must be an N x 3 matrix.");
	}

	if (mxGetN(prhs[1]) != 1)
	{
		mexErrMsgTxt("The label of the input point cloud must be an N x 1 vector.");
	}

	if (mxGetM(prhs[0]) != mxGetM(prhs[1]))
	{
		mexErrMsgTxt("The rows of input_cloud and labels must be the same."); 
	}

	if (mxGetN(prhs[2]) != 3)
	{
		mexErrMsgTxt("The color table must be an N x 3 matrix.");
	}

	int num = mxGetM(prhs[0]);
	double *pt = mxGetPr(prhs[0]);
	double *labels = mxGetPr(prhs[1]);
	double *color = mxGetPr(prhs[2]);

	if (nrhs == 4)
	{
		int temp = mxGetN(prhs[3]);
		if (mxGetM(prhs[3]) != 1 || (temp != 1 && temp != 4))
		{
			mexErrMsgTxt("The parameters must be a 1 x 1 or 1 x 4 vector.");
		}

		double *param = mxGetPr(prhs[3]);
		point_size = *param;

		if (temp == 4)
		{
			background[0] = *(param + 1);
			background[1] = *(param + 2);
			background[2] = *(param + 3);
		}	
	}

	// Generate the color table
	int num_label = 0;
	for (int i = 0; i < num; i++)
	{
		int temp_label = (int)(*(labels + i));
		num_label = (temp_label > num_label) ? (temp_label) : (num_label);
	}

	if (num_label != mxGetM(prhs[2]))
	{
		mexErrMsgTxt("The number of color must be equal to the number of regions.");
	}

	// Construct the input point cloud
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>(num, 1));
	for (unsigned int i = 0; i < cloud->points.size(); i++)
	{
		cloud->points[i].x = *(pt + i);
		cloud->points[i].y = *(pt + num + i);
		cloud->points[i].z = *(pt + 2 * num + i);

		int temp_label = (int)(*(labels + i)) - 1;
		cloud->points[i].r = pcl_lrint(*(color + temp_label) * 255.0);
		cloud->points[i].g = pcl_lrint(*(color + num_label + temp_label) * 255.0);
		cloud->points[i].b = pcl_lrint(*(color + 2 * num_label + temp_label) * 255.0);
	}

	// View the point cloud
	pcl::visualization::CloudViewer viewer("Cloud Viewer");
	viewer.showCloud(cloud);
	viewer.runOnVisualizationThreadOnce (viewerSetting);

	while (!viewer.wasStopped())
	{
		// No code
	}

	// Postprocess
	if (nlhs != 0)
	{
		mexErrMsgTxt("No output argument required.");
	}

	return;
}