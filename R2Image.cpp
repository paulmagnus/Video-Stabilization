// Source file for image class

// Include files

#define _USE_MATH_DEFINES
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <float.h>
#include <string>
#include <sstream>

// comment to remove print statements
#define PRINT_PROGRESS_REPORT

void print(const char* format, ...) {
#ifdef PRINT_PROGRESS_REPORT
	va_list argptr;
	va_start(argptr, format);
	vfprintf(stderr, format, argptr);
	va_end(argptr);
#endif
}

using namespace std;

typedef struct HarrisContainer {
	int posx;
	int posy;
	float harrisScore;

	HarrisContainer(int x, int y, float harris) {
		posx = x;
		posy = y;
		harrisScore = harris;
	}
};

struct Pair {
	vec2 a, b;
	Pair() : a(), b() {};
	Pair(vec2 & a, vec2 & b) : a(a), b(b) {};
};

mat3 getHomographyMatrix(const std::vector<Pair> & pairs) {
	double** linEquations = dmatrix(1, pairs.size() * 2, 1, 9);

	for (int i = 0; i < pairs.size(); i++) {
		Pair pair = pairs[i];

		linEquations[2 * i + 1][1] = 0;
		linEquations[2 * i + 1][2] = 0;
		linEquations[2 * i + 1][3] = 0;
		linEquations[2 * i + 1][4] = -pair.a.x;
		linEquations[2 * i + 1][5] = -pair.a.y;
		linEquations[2 * i + 1][6] = -1;
		linEquations[2 * i + 1][7] = pair.b.y * pair.a.x;
		linEquations[2 * i + 1][8] = pair.b.y * pair.a.y;
		linEquations[2 * i + 1][9] = pair.b.y;

		linEquations[2 * i + 2][1] = pair.a.x;
		linEquations[2 * i + 2][2] = pair.a.y;
		linEquations[2 * i + 2][3] = 1;
		linEquations[2 * i + 2][4] = 0;
		linEquations[2 * i + 2][5] = 0;
		linEquations[2 * i + 2][6] = 0;
		linEquations[2 * i + 2][7] = -pair.a.x * pair.b.x;
		linEquations[2 * i + 2][8] = -pair.a.x * pair.b.y;
		linEquations[2 * i + 2][9] = -pair.b.x;
	}

	double singularValues[10] = { 0,0,0,0,0,0,0,0,0,0 };
	double** nullspaceMatrix = dmatrix(1, 9, 1, 9);

	svdcmp(linEquations, pairs.size() * 2, 9, singularValues, nullspaceMatrix);

	// find smallest singular value
	int smallestIndex = 1;
	for (int i = 2; i < 10; i++)
		if (singularValues[i] < singularValues[smallestIndex])
			smallestIndex = i;

	return mat3(nullspaceMatrix[1][smallestIndex], nullspaceMatrix[2][smallestIndex], nullspaceMatrix[3][smallestIndex],
		nullspaceMatrix[4][smallestIndex], nullspaceMatrix[5][smallestIndex], nullspaceMatrix[6][smallestIndex],
		nullspaceMatrix[7][smallestIndex], nullspaceMatrix[8][smallestIndex], nullspaceMatrix[9][smallestIndex]);
}

bool isElement(std::vector<unsigned int> v, unsigned int e) {
	for (int i = 0; i < v.size(); i++)
		if (v[i] == e)
			return true;
	return false;
}

std::vector<unsigned int> getRandoms(int nRand, int mod) {
	std::vector<unsigned int> rands;
	while (rands.size() < nRand) {
		unsigned int r = rand() % mod;
		if (!isElement(rands, r))
			rands.push_back(r);
	}
	return rands;
}

bool harris_sort_function(const HarrisContainer &left, const HarrisContainer &right) {
	return left.harrisScore > right.harrisScore;
}

/* A simple helper function to allow for easy edge detection */
int bound(int val, int max) {
	if (val < 0)
		return 0;
	if (val >= max)
		return max - 1;
	return val;
}

////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
	: pixels(NULL),
	npixels(0),
	width(0),
	height(0)
{
}



R2Image::
R2Image(const char *filename)
	: pixels(NULL),
	npixels(0),
	width(0),
	height(0)
{
	// Read image
	Read(filename);
}



R2Image::
R2Image(int width, int height)
	: pixels(NULL),
	npixels(width * height),
	width(width),
	height(height)
{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
	: pixels(NULL),
	npixels(width * height),
	width(width),
	height(height)
{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels
	for (int i = 0; i < npixels; i++)
		pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
	: pixels(NULL),
	npixels(image.npixels),
	width(image.width),
	height(image.height)

{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels
	for (int i = 0; i < npixels; i++)
		pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
	// Free image pixels
	if (pixels) delete[] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
	// Delete previous pixels
	if (pixels) { delete[] pixels; pixels = NULL; }

	// Reset width and height
	npixels = image.npixels;
	width = image.width;
	height = image.height;

	// Allocate new pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels
	for (int i = 0; i < npixels; i++)
		pixels[i] = image.pixels[i];

	// Return image
	return *this;
}


void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	R2Point p1(1.2, 3.5);
	R2Point p2(2.1, 2.2);
	R2Point p3(0.2, 1.6);
	R2Point p4(0.0, 0.5);
	R2Point p5(-0.2, 4.2);

	// build the 5x6 matrix of equations
	double** linEquations = dmatrix(1, 5, 1, 6);

	linEquations[1][1] = p1[0] * p1[0];
	linEquations[1][2] = p1[0] * p1[1];
	linEquations[1][3] = p1[1] * p1[1];
	linEquations[1][4] = p1[0];
	linEquations[1][5] = p1[1];
	linEquations[1][6] = 1.0;

	linEquations[2][1] = p2[0] * p2[0];
	linEquations[2][2] = p2[0] * p2[1];
	linEquations[2][3] = p2[1] * p2[1];
	linEquations[2][4] = p2[0];
	linEquations[2][5] = p2[1];
	linEquations[2][6] = 1.0;

	linEquations[3][1] = p3[0] * p3[0];
	linEquations[3][2] = p3[0] * p3[1];
	linEquations[3][3] = p3[1] * p3[1];
	linEquations[3][4] = p3[0];
	linEquations[3][5] = p3[1];
	linEquations[3][6] = 1.0;

	linEquations[4][1] = p4[0] * p4[0];
	linEquations[4][2] = p4[0] * p4[1];
	linEquations[4][3] = p4[1] * p4[1];
	linEquations[4][4] = p4[0];
	linEquations[4][5] = p4[1];
	linEquations[4][6] = 1.0;

	linEquations[5][1] = p5[0] * p5[0];
	linEquations[5][2] = p5[0] * p5[1];
	linEquations[5][3] = p5[1] * p5[1];
	linEquations[5][4] = p5[0];
	linEquations[5][5] = p5[1];
	linEquations[5][6] = 1.0;

	printf("\n Fitting a conic to five points:\n");
	printf("Point #1: %f,%f\n", p1[0], p1[1]);
	printf("Point #2: %f,%f\n", p2[0], p2[1]);
	printf("Point #3: %f,%f\n", p3[0], p3[1]);
	printf("Point #4: %f,%f\n", p4[0], p4[1]);
	printf("Point #5: %f,%f\n", p5[0], p5[1]);

	// compute the SVD
	double singularValues[7]; // 1..6
	double** nullspaceMatrix = dmatrix(1, 6, 1, 6);
	svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// get the result
	printf("\n Singular values: %f, %f, %f, %f, %f, %f\n", singularValues[1], singularValues[2], singularValues[3], singularValues[4], singularValues[5], singularValues[6]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for (int i = 2; i < 7; i++) if (singularValues[i] < singularValues[smallestIndex]) smallestIndex = i;

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n", nullspaceMatrix[1][smallestIndex], nullspaceMatrix[2][smallestIndex], nullspaceMatrix[3][smallestIndex], nullspaceMatrix[4][smallestIndex], nullspaceMatrix[5][smallestIndex], nullspaceMatrix[6][smallestIndex]);

	// make sure the solution is correct:
	printf("Equation #1 result: %f\n", p1[0] * p1[0] * nullspaceMatrix[1][smallestIndex] +
		p1[0] * p1[1] * nullspaceMatrix[2][smallestIndex] +
		p1[1] * p1[1] * nullspaceMatrix[3][smallestIndex] +
		p1[0] * nullspaceMatrix[4][smallestIndex] +
		p1[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	printf("Equation #2 result: %f\n", p2[0] * p2[0] * nullspaceMatrix[1][smallestIndex] +
		p2[0] * p2[1] * nullspaceMatrix[2][smallestIndex] +
		p2[1] * p2[1] * nullspaceMatrix[3][smallestIndex] +
		p2[0] * nullspaceMatrix[4][smallestIndex] +
		p2[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	printf("Equation #3 result: %f\n", p3[0] * p3[0] * nullspaceMatrix[1][smallestIndex] +
		p3[0] * p3[1] * nullspaceMatrix[2][smallestIndex] +
		p3[1] * p3[1] * nullspaceMatrix[3][smallestIndex] +
		p3[0] * nullspaceMatrix[4][smallestIndex] +
		p3[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	printf("Equation #4 result: %f\n", p4[0] * p4[0] * nullspaceMatrix[1][smallestIndex] +
		p4[0] * p4[1] * nullspaceMatrix[2][smallestIndex] +
		p4[1] * p4[1] * nullspaceMatrix[3][smallestIndex] +
		p4[0] * nullspaceMatrix[4][smallestIndex] +
		p4[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	printf("Equation #5 result: %f\n", p5[0] * p5[0] * nullspaceMatrix[1][smallestIndex] +
		p5[0] * p5[1] * nullspaceMatrix[2][smallestIndex] +
		p5[1] * p5[1] * nullspaceMatrix[3][smallestIndex] +
		p5[0] * nullspaceMatrix[4][smallestIndex] +
		p5[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	R2Point test_point(0.34, -2.8);

	printf("A point off the conic: %f\n", test_point[0] * test_point[0] * nullspaceMatrix[1][smallestIndex] +
		test_point[0] * test_point[1] * nullspaceMatrix[2][smallestIndex] +
		test_point[1] * test_point[1] * nullspaceMatrix[3][smallestIndex] +
		test_point[0] * nullspaceMatrix[4][smallestIndex] +
		test_point[1] * nullspaceMatrix[5][smallestIndex] +
		nullspaceMatrix[6][smallestIndex]);

	return;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
	// Brighten the image by multiplying each pixel component by the factor.
	// This is implemented for you as an example of how to access and set pixels
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			Pixel(i, j) *= factor;
			Pixel(i, j).Clamp();
		}
	}
}

float pixelAverage(R2Pixel pixel) {
#pragma warning(push)
#pragma warning(disable: 4244)
	return sqrtf(pixel.Blue() * pixel.Blue()
		+ pixel.Green() * pixel.Green()
		+ pixel.Red() * pixel.Red());
#pragma warning(pop)
}

void R2Image::
SobelX(void)
{
	// Apply the Sobel oprator to the image in X direction

	float kernel[3][3] = { {-0.25, 0, 0.25},
						{-0.5, 0, 0.5},
						{-0.25, 0, 0.25} };

	R2Image tempImage(*this);
	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++) {
			float sum = 0;
			for (int m = -1; m < 2; m++)
				for (int n = -1; n < 2; n++)
					sum += pixelAverage(tempImage.Pixel(bound(x + m, width), bound(y + n, height))) * kernel[n + 1][m + 1];
			Pixel(x, y) = R2Pixel(sum, sum, sum, 1);
		}
}

void R2Image::
SobelY(void)
{
	// Apply the Sobel oprator to the image in Y direction

	float kernel[3][3] = { {-0.25, -0.5, -0.25},
						{0, 0, 0},
						{0.25, 0.5, 0.25} };

	R2Image tempImage(*this);
	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++) {
			float sum = 0;
			for (int m = -1; m < 2; m++)
				for (int n = -1; n < 2; n++)
					sum += pixelAverage(tempImage.Pixel(bound(x + m, width), bound(y + n, height))) * kernel[n + 1][m + 1];
			Pixel(x, y) = R2Pixel(sum, sum, sum, 1);
		}
}

void R2Image::
LoG(void)
{
	// Apply the LoG oprator to the image

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
	// Changes the saturation of an image
	// Find a formula that changes the saturation without affecting the image brightness

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}

/* A simple helper function for the blur method so that I can computer the gaussian for any sigma*/
float computeGaussian(float sigma, float position) {
#pragma warning(suppress: 4244)
	return pow(M_E, -(float)(position * position) / (2 * sigma * sigma)) / sqrt(2 * M_PI * sigma * sigma);
}

// Linear filtering ////////////////////////////////////////////////

void R2Image::
Blur(double sigma)
{
	// Gaussian blur of the image. Separable solution is preferred
	std::vector<float> kernel(0);

	float sum = 0;
	// compute the gaussian values for the kernel
	for (int i = 0; i < ceil(6 * sigma + 1); i++) {
#pragma warning(suppress: 4244)
		float gauss = computeGaussian(sigma, i - 3 * sigma);
		kernel.push_back(gauss);
		sum += gauss;
	}

	for (unsigned int i = 0; i < kernel.size(); i++)
		kernel[i] = kernel[i] / sum;

	R2Image tempImage(width, height);

	// blur in the x direction
	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++) {
			R2Pixel sum(0, 0, 0, 1);
			for (unsigned int i = 0; i < kernel.size(); i++) {
#pragma warning(suppress: 4244)
				sum += Pixel(bound(ceil(x + i - 3 * sigma), width), y) * kernel.at(i);
			}
			tempImage.Pixel(x, y) = sum;
		}

	// blur in the y direction
	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++) {
			R2Pixel sum(0, 0, 0, 1);
			for (unsigned int i = 0; i < kernel.size(); i++) {
#pragma warning(suppress: 4244)
				sum += tempImage.Pixel(x, bound(ceil(y + i - 3 * sigma), height)) * kernel.at(i);
			}
			Pixel(x, y) = sum;
		}
}

void R2Image::
Harris(double sigma)
{
	// Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges

	print("Calculating base images...");

	R2Image Img1(*this);
	Img1.SobelX();

	R2Image Img2(*this);
	Img2.SobelY();

	R2Image Img3(Img2);

	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++) {
			Img3.Pixel(x, y) = Img3.Pixel(x, y) * Img1.Pixel(x, y);
			Img2.Pixel(x, y) = Img2.Pixel(x, y) * Img2.Pixel(x, y);
			Img1.Pixel(x, y) = Img1.Pixel(x, y) * Img1.Pixel(x, y);
		}
	print("Done\n");

	print("Blurring Sobel X squared...");
	Img1.Blur(sigma);
	print("Done\n");

	print("Blurring Sobel Y squared...");
	Img2.Blur(sigma);
	print("Done\n");

	print("Blurring Sobel X * Sobel Y...");
	Img3.Blur(sigma);
	print("Done\n");

	print("Calculating final pixels");
	for (int x = 0; x < Img1.width; x++)
		for (int y = 0; y < Img1.height; y++) {
			Pixel(x, y) = 4 * (Img1.Pixel(x, y) * Img2.Pixel(x, y) - Img3.Pixel(x, y) * Img3.Pixel(x, y) - 0.04 * ((Img1.Pixel(x, y) + Img2.Pixel(x, y)) * (Img1.Pixel(x, y) + Img2.Pixel(x, y))));
			Pixel(x, y) = Pixel(x, y) + R2Pixel(0.5, 0.5, 0.5, 1);
			Pixel(x, y).Clamp();
		}
	print("Done\n");
}

float getSSD(int xRange, int yRange, R2Image& img, int x1, int y1, R2Image& other, int x2, int y2) {
	float sum = 0;
	for (int dx = -1 * xRange; dx < xRange; dx++) {
		for (int dy = -1 * yRange; dy < yRange; dy++) {
			R2Pixel diff = img.Pixel(bound(x1 + dx, img.Width()), bound(y1 + dy, img.Height()))
				- other.Pixel(bound(x2 + dx, other.Width()), bound(y2 + dy, other.Height()));
			diff = diff * diff;
			sum += diff.Red() + diff.Green() + diff.Blue();
		}
	}
	return sum;
}

void R2Image::line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
	/*
	if (x0 > 3 && x0 < width - 3 && y0>3 && y0 < height - 3)
	{
		for (int x = x0 - 3; x <= x0 + 3; x++)
		{
			for (int y = y0 - 3; y <= y0 + 3; y++)
			{
				Pixel(x, y).Reset(r, g, b, 1.0);
			}
		}
	}
	*/

	if (x0 > x1)
	{
		int x = y1;
		y1 = y0;
		y0 = x;

		x = x1;
		x1 = x0;
		x0 = x;
	}
	int deltax = x1 - x0;
	int deltay = y1 - y0;
	float error = 0;
	float deltaerr = 0.0;
	if (deltax != 0) deltaerr = fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
																		// note that this division needs to be done in a way that preserves the fractional part
	int y = y0;
	for (int x = x0; x <= x1; x++)
	{
		Pixel(x, y).Reset(r, g, b, 1.0);
		error = error + deltaerr;
		if (error >= 0.5)
		{
			if (deltay > 0) y = y + 1;
			else y = y - 1;

			error = error - 1.0;
		}
	}
}

void R2Image::TranslationMatch(double sigma, int numFeatures, R2Image* other) {

	/* Find features in this image */
	R2Image harris(*this);
	harris.Harris(sigma);
	print("Done computing harris\n");

	std::vector<vec2> featurePoints(0);
	print("Obtain %d feature points\n", numFeatures);
	while (featurePoints.size() < numFeatures) {
		if (featurePoints.size() % 10 == 0) {
			print("\33[2K\r");
			print("%d features obtained", featurePoints.size());
		}

		// find next highest harris value
		vec2 pixel;
		float harrisVal = 0;
		for (int x = 3 * sigma; x < width - 3 * sigma; x++)	// do not consider edge features
			for (int y = 3 * sigma; y < height - 3 * sigma; y++) {	// do not consider edge features
				float ave = pixelAverage(harris.Pixel(x, y));
				if (ave > harrisVal) {
					harrisVal = ave;
					pixel.x = x;
					pixel.y = y;
				}
			}
		featurePoints.push_back(pixel);

		// no pixel within 10 pixels should be considered further
		for (int dx = -10; dx <= 10; dx++)
			for (int dy = -10; dy <= 10; dy++)
				harris.Pixel(bound(pixel.x + dx, width), bound(pixel.y + dy, height)) = R2Pixel(0, 0, 0, 1);
	}
	print("\33[2K\r");
	print("%d features obtained\n", featurePoints.size());
	print("Collected all feature points\n\n");

	print("Searching for translated features\n");
	/* Search for the same features in the other image */
	std::vector<vec2> otherFeatures(0);
	for (unsigned int i = 0; i < featurePoints.size(); i++) {
		print("\r%u features found", otherFeatures.size());
		float bestSSD = FLT_MAX;
		vec2 loc(0, 0);
		// search near original feature location
		for (int dx = -width * 0.1; dx < width * 0.1; dx++) {
			for (int dy = -height * 0.1; dy < height * 0.1; dy++) {
				float ssd = getSSD(3 * sigma, 3 * sigma, *this, featurePoints[i].x, featurePoints[i].y, *other, featurePoints[i].x + dx, featurePoints[i].y + dy);
				if (ssd <= bestSSD) {
					bestSSD = ssd;
					loc.x = bound(featurePoints[i].x + dx, width);
					loc.y = bound(featurePoints[i].y + dy, height);
				}
			}
		}

		if (loc.x == -1 || loc.y == -1) {
			print("\nError in SSD calculation");
			exit(1);
		}

		otherFeatures.push_back(loc);
	}

	if (featurePoints.size() != otherFeatures.size()) {
		print("\nError in finding translated features");
		exit(1);
	}

	print("\nDetermining false matches");

	std::vector<unsigned int> tracks;
	// do N (50) trials
	for (int n = 0; n < numFeatures / 10; n++) {
		std::vector<unsigned int> inlierTracks;

		// randomly select a single track
		int index = rand() % numFeatures;
		vec2 randomTrack = otherFeatures[index] - featurePoints[index];

		// check all features
		for (unsigned int i = 0; i < numFeatures; i++) {
			// project features according to random track
			vec2 projectedFeature = featurePoints[i] + randomTrack;

			// determine distance between actual motion and projected motion
			double dist = (projectedFeature - otherFeatures[i]).length();

			// if within distance threshold add to inlierTracks
			if (dist <= 4)
				inlierTracks.push_back(i);
		}
		if (inlierTracks.size() > tracks.size())
			tracks = inlierTracks;

		print("Tracks: %u\n", tracks.size());
	}
	printf("Done\n");

	// draw feature changes
	print("\nDrawing...");
	int index = 0;
	for (unsigned int i = 0; i < featurePoints.size(); i++)
		if (index < tracks.size() && tracks[index] == i) {
			line(featurePoints[i].x, otherFeatures[i].x, featurePoints[i].y, otherFeatures[i].y, 0, 1, 0);
			index++;
		}
		else
			line(featurePoints[i].x, otherFeatures[i].x, featurePoints[i].y, otherFeatures[i].y, 1, 0, 0);
	print("\nDone\n");
}

void R2Image::Features(double sigma, int numFeatures) {
	R2Image harris(*this);
	harris.Harris(sigma);
	print("Done computing harris\n");

	std::vector<vec2> featurePoints(0);
	print("Obtain %d feature points\n", numFeatures);
	while (featurePoints.size() < numFeatures) {
		if (featurePoints.size() % 10 == 0) {
			print("\33[2K\r");
			print("%d features obtained", featurePoints.size());
		}

		// find next highest harris value
		vec2 pixel;
		float harrisVal = 0;
		for (int x = 3 * sigma; x < width - 3 * sigma; x++)	// do not consider edge features
			for (int y = 3 * sigma; y < height - 3 * sigma; y++) {	// do not consider edge features
				float ave = pixelAverage(harris.Pixel(x, y));
				if (ave > harrisVal) {
					harrisVal = ave;
					pixel.x = x;
					pixel.y = y;
				}
			}
		featurePoints.push_back(pixel);

		// no pixel within 10 pixels should be considered further
		for (int dx = -10; dx <= 10; dx++)
			for (int dy = -10; dy <= 10; dy++)
				harris.Pixel(bound(pixel.x + dx, width), bound(pixel.y + dy, height)) = R2Pixel(0, 0, 0, 1);
	}
	print("\33[2K\r");
	print("%d features obtained\n", featurePoints.size());
	print("Collected all feature points\n");

	//	featurePoints[0].print();

	print("Drawing feature points...");
	for (unsigned int i = 0; i < featurePoints.size(); i++) {
		int x = featurePoints[i].x - 3 * sigma;
		for (int j = -3 * sigma; j < 3 * sigma; j++)
			Pixel(x, featurePoints[i].y + j) = R2Pixel(1, 0, 0, 1);
		x = featurePoints[i].x + 3 * sigma;
		for (int j = -3 * sigma; j < 3 * sigma; j++)
			Pixel(x, featurePoints[i].y + j) = R2Pixel(1, 0, 0, 1);
		int y = featurePoints[i].y - 3 * sigma;
		for (int j = -3 * sigma; j < 3 * sigma; j++)
			Pixel(featurePoints[i].x + j, y) = R2Pixel(1, 0, 0, 1);
		y = featurePoints[i].y + 3 * sigma;
		for (int j = -3 * sigma; j < 3 * sigma; j++)
			Pixel(featurePoints[i].x + j, y) = R2Pixel(1, 0, 0, 1);
	}
	print("Done\n");
}


void R2Image::HomographyMatch(double sigma, int numFeatures, R2Image* other) {
	/* Find features in this image */
	R2Image harris(*this);
	harris.Harris(sigma);
	print("Done computing harris\n");

	std::vector<vec2> featurePoints(0);
	print("Obtain %d feature points\n", numFeatures);
	while (featurePoints.size() < numFeatures) {
		if (featurePoints.size() % 10 == 0) {
			print("\33[2K\r");
			print("%d features obtained", featurePoints.size());
		}

		// find next highest harris value
		vec2 pixel;
		float harrisVal = 0;
		for (int x = 3 * sigma; x < width - 3 * sigma; x++)	// do not consider edge features
			for (int y = 3 * sigma; y < height - 3 * sigma; y++) {	// do not consider edge features
				float ave = pixelAverage(harris.Pixel(x, y));
				if (ave > harrisVal) {
					harrisVal = ave;
					pixel.x = x;
					pixel.y = y;
				}
			}
		featurePoints.push_back(pixel);

		// no pixel within 10 pixels should be considered further
		for (int dx = -10; dx <= 10; dx++)
			for (int dy = -10; dy <= 10; dy++)
				harris.Pixel(bound(pixel.x + dx, width), bound(pixel.y + dy, height)) = R2Pixel(0, 0, 0, 1);
	}
	print("\33[2K\r");
	print("%d features obtained\n", featurePoints.size());
	print("Collected all feature points\n\n");

	print("Searching for transformed features\n");
	/* Search for the same features in the other image */
	std::vector<vec2> otherFeatures(0);
	for (unsigned int i = 0; i < featurePoints.size(); i++) {
		print("\r%u features found", otherFeatures.size());
		float bestSSD = FLT_MAX;
		vec2 loc(0, 0);
		// search near original feature location
		for (int dx = -width * 0.1; dx < width * 0.1; dx++) {
			for (int dy = -height * 0.1; dy < height * 0.1; dy++) {
				float ssd = getSSD(3 * sigma, 3 * sigma, *this, featurePoints[i].x, featurePoints[i].y, *other, featurePoints[i].x + dx, featurePoints[i].y + dy);
				if (ssd <= bestSSD) {
					bestSSD = ssd;
					loc.x = bound(featurePoints[i].x + dx, width);
					loc.y = bound(featurePoints[i].y + dy, height);
				}
			}
		}

		if (loc.x == -1 || loc.y == -1) {
			print("\nError in SSD calculation");
			exit(1);
		}

		otherFeatures.push_back(loc);
	}

	if (featurePoints.size() != otherFeatures.size()) {
		print("\nError in finding transformed features");
		exit(1);
	}

	print("\nDetermining false matches\n");

	std::vector<unsigned int> tracks;

	for (int n = 0; n < numFeatures / 10; n++) {
		std::vector<unsigned int> inlierTracks;

		// randomly select 4 tracks
		std::vector<unsigned int> indecies = getRandoms(4, numFeatures);

		std::vector<Pair> pairs;
		for (int i = 0; i < 4; i++) {
			pairs.push_back(Pair(featurePoints[indecies[i]], otherFeatures[indecies[i]]));
		}


		mat3 transformationMatrix = getHomographyMatrix(pairs);

		// check all features
		for (unsigned int i = 0; i < numFeatures; i++) {
			// project features according to transformation matrix
			vec3 feature(featurePoints[i].x, featurePoints[i].y, 1);
			vec3 projectedFeatureH = transformationMatrix * feature;

			if (projectedFeatureH.z != 0) {
				vec2 projectedFeature(projectedFeatureH.x / projectedFeatureH.z, projectedFeatureH.y / projectedFeatureH.z);

				// determine distance between actual motion and projected motion
				double dist = (projectedFeature - otherFeatures[i]).length();

				// if within distance threshold, add to inlierTracks
				if (dist <= 4)
					inlierTracks.push_back(i);
			}
		}
		if (inlierTracks.size() > tracks.size())
			tracks = inlierTracks;
	}
	print("Done\n");

	// draw feature changes
	print("\nDrawing...");
	int index = 0;
	for (unsigned int i = 0; i < featurePoints.size(); i++)
		if (index < tracks.size() && tracks[index] == i) {
			line(featurePoints[i].x, otherFeatures[i].x, featurePoints[i].y, otherFeatures[i].y, 0, 1, 0);
			index++;
		}
		else
			line(featurePoints[i].x, otherFeatures[i].x, featurePoints[i].y, otherFeatures[i].y, 1, 0, 0);
	print("\nDone\n");
}

void R2Image::
Sharpen()
{
	// Sharpen an image using a linear filter. Use a kernel of your choosing.

	int kernel[3][3] = { {-1, -1, -1},
					{-1, 8, -1},
					{-1, -1, -1} };

	R2Image tempImage(*this);
	for (int x = 1; x < width - 1; x++)
		for (int y = 1; y < height - 1; y++) {
			R2Pixel sum(0, 0, 0, 1);
			for (int m = -1; m < 2; m++)
				for (int n = -1; n < 2; n++)
					sum += tempImage.Pixel(x + m, y + n) * kernel[n + 1][m + 1];
			Pixel(x, y) = sum / 8 + Pixel(x, y);
			Pixel(x, y).Clamp();
		}
}


void R2Image::
blendOtherImageTranslated(R2Image * otherImage, double sigma, int numFeatures)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage"
	// into this image with a 50% opacity.
	TranslationMatch(sigma, numFeatures, otherImage);
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage, double sigma, int numFeatures)
{
	R2Image harris(*this);
	harris.Harris(sigma);
	print("Done computing harris\n");

	std::vector<vec2> featurePoints(0);
	print("Obtain %d feature points\n", numFeatures);
	while (featurePoints.size() < numFeatures) {
		if (featurePoints.size() % 10 == 0) {
			print("\33[2K\r");
			print("%d features obtained", featurePoints.size());
		}

		// find next highest harris value
		vec2 pixel;
		float harrisVal = 0;
		for (int x = 3 * sigma; x < width - 3 * sigma; x++)	// do not consider edge features
			for (int y = 3 * sigma; y < height - 3 * sigma; y++) {	// do not consider edge features
				float ave = pixelAverage(harris.Pixel(x, y));
				if (ave > harrisVal) {
					harrisVal = ave;
					pixel.x = x;
					pixel.y = y;
				}
			}
		featurePoints.push_back(pixel);

		// no pixel within 10 pixels should be considered further
		for (int dx = -10; dx <= 10; dx++)
			for (int dy = -10; dy <= 10; dy++)
				harris.Pixel(bound(pixel.x + dx, width), bound(pixel.y + dy, height)) = R2Pixel(0, 0, 0, 1);
	}
	print("\33[2K\r");
	print("%d features obtained\n", featurePoints.size());
	print("Collected all feature points\n\n");

	print("Searching for transformed features\n");
	/* Search for the same features in the other image */
	std::vector<vec2> otherFeatures(0);
	for (unsigned int i = 0; i < featurePoints.size(); i++) {
		print("\r%u features found", otherFeatures.size());
		float bestSSD = FLT_MAX;
		vec2 loc(0, 0);
		// search near original feature location
		for (int dx = -width * 0.1; dx < width * 0.1; dx++) {
			for (int dy = -height * 0.1; dy < height * 0.1; dy++) {
				float ssd = getSSD(3 * sigma, 3 * sigma, *this, featurePoints[i].x, featurePoints[i].y, *otherImage, featurePoints[i].x + dx, featurePoints[i].y + dy);
				if (ssd <= bestSSD) {
					bestSSD = ssd;
					loc.x = bound(featurePoints[i].x + dx, width);
					loc.y = bound(featurePoints[i].y + dy, height);
				}
			}
		}

		if (loc.x == -1 || loc.y == -1) {
			print("\nError in SSD calculation");
			exit(1);
		}

		otherFeatures.push_back(loc);
	}

	if (featurePoints.size() != otherFeatures.size()) {
		print("\nError in finding transformed features");
		exit(1);
	}

	print("\nDetermining false matches\n");

	std::vector<unsigned int> tracks;

	mat3 bestTransformMatrix(0, 0, 0, 0, 0, 0, 0, 0, 0);

	for (int n = 0; n < numFeatures / 10; n++) {
		std::vector<unsigned int> inlierTracks;

		// randomly select 4 tracks
		std::vector<unsigned int> indecies = getRandoms(4, numFeatures);

		std::vector<Pair> pairs;
		for (int i = 0; i < 4; i++) {
			pairs.push_back(Pair(featurePoints[indecies[i]], otherFeatures[indecies[i]]));
		}


		mat3 transformationMatrix = getHomographyMatrix(pairs);

		// check all features
		for (unsigned int i = 0; i < numFeatures; i++) {
			// project features according to transformation matrix
			vec3 feature(featurePoints[i].x, featurePoints[i].y, 1);
			vec3 projectedFeatureH = transformationMatrix * feature;

			if (projectedFeatureH.z != 0) {
				vec2 projectedFeature(projectedFeatureH.x / projectedFeatureH.z, projectedFeatureH.y / projectedFeatureH.z);

				// determine distance between actual motion and projected motion
				double dist = (projectedFeature - otherFeatures[i]).length();

				// if within distance threshold, add to inlierTracks
				if (dist <= 4)
					inlierTracks.push_back(i);
			}
		}
		if (inlierTracks.size() > tracks.size()) {
			tracks = inlierTracks;
			bestTransformMatrix = transformationMatrix;
		}
	}
	print("Done\n");

	print("Blending images...");
	for (double x = 0; x < width; x++) {
		for (double y = 0; y < height; y++) {
			vec3 projectedPoint = bestTransformMatrix * vec3(x, y, 1);
			if (projectedPoint.z == 0)
				print("Error: w value was 0");
			else {
				int projectedX = round(projectedPoint.x / projectedPoint.z);
				int projectedY = round(projectedPoint.y / projectedPoint.z);

				if (projectedX >= 0 && projectedX < width && projectedY >= 0 && projectedY < height) {
					// projected point is on the other image
					Pixel(x, y) = Pixel(x, y) * 0.5 + otherImage->Pixel(projectedX, projectedY) * 0.5;
				}
			}
		}
	}
	print("Done\n");
}

std::string
padZeroNumber(int num) {
	/*
	stringstream ss;

	ss << num;
	string ret;
	ss >> ret;
	*/

	std::string ret = std::to_string(num);

	int str_length = ret.length();
	for (int i = 0; i < 5 - str_length; i++) {
		ret = "0" + ret;
	}
	return ret;
}

// update 12/8
int R2Image::transformFrames(const std::vector<R2Image*> & images,
	const std::vector<mat3*> & transformationMatricies,
	const std::string & folderName) {

	assert(images.size() == transformationMatricies.size());

	for (int i = 0; i < images.size(); i++) {
		std::string current = padZeroNumber(i);
		std::string curFileName = folderName + "/" + current + ".jpg";
		if (!images[i]->transformFrame(transformationMatricies[i], curFileName.c_str()))
			return 0;	// failure
	}
	// success
	return 1;
}

// update 12/8
int R2Image::transformFrame(mat3 * transformMatrix, const char* filename) {
	R2Image tmpImage(*this);
	for (double x = 0; x < width; x++) {
		for (double y = 0; y < height; y++) {
			vec3 projectedPoint = *transformMatrix * vec3(x, y, 1);
			if (projectedPoint.z == 0)
				print("Error: w value was 0");
			else {
				int projectedX = round(projectedPoint.x / projectedPoint.z);
				int projectedY = round(projectedPoint.y / projectedPoint.z);

				if (projectedX >= 0 && projectedX < width && projectedY >= 0 && projectedY < height) {
					// projected point is on the other image
					tmpImage.Pixel(x, y) = this->Pixel(projectedX, projectedY);
				}
				else
					tmpImage.SetPixel(x, y, R2Pixel(0, 0, 0, 1));
			}
		}
	}
	return tmpImage.Write(filename);
}

bool pointIsValid(vec2 point, int width, int height) {
	return point.x >= 0 && point.x < width && point.y >= 0 && point.y < height;
}

void R2Image::stabilize(std::string dirName, int numImage, int numFeatures, double sigma) {
	if (numFeatures < 4) {
		fprintf(stderr, "Number of features must be at least 4 you gave %d features", numFeatures);
		exit(1);
	}

	std::vector<track> tracks = findTrack(dirName, numImage, numFeatures, sigma);

	for (track t : tracks) {
		for (int i = 0; i < t.size() - 1; i++) {
			vec2 point1 = t.getPosition(i);
			vec2 point2 = t.getPosition(i + 1);
			if(pointIsValid(point1, width, height) && pointIsValid(point2, width, height))
				line(point1.x, point2.x, point1.y, point2.y, 1, 0, 0);
		}
	}
}


struct feature {
	vec2 point;
	int track;

	feature() : point(), track(0) {}
	feature(vec2 p, int t) : point(p), track(t) {}
};

struct list {
	std::vector<feature> v;
	list() : v(0) {}
	list(std::vector<feature> vec) : v(vec) {}
	feature at(int i) {
		return v.at(i);
	}
	size_t size() {
		return v.size();
	}
};

bool contains(list v, const vec2 & o) {
	for (int i = 0; i < v.size(); i++)
		if (v.at(i).point.x == o.x && v.at(i).point.y == o.y)
			return true;
	return false;
}

std::vector<vec2> toVec2List(const std::vector<feature> & list) {
	std::vector<vec2> r(0);
	for (feature f : list)
		r.push_back(f.point);
	return r;
}

struct fPair {
	feature a, b;
	fPair(feature & a, feature & b) : a(a), b(b) {};
};

std::vector<Pair> toVec2List(const std::vector<fPair> & list) {
	std::vector<Pair> r(0);
	for (fPair f : list)
		r.push_back(Pair(f.a.point, f.b.point));
	return r;
}

std::vector<track>
R2Image::findTrack(std::string dirName, int numImages, int numFeatures, double sigma) {

	std::vector<track> tracks(numFeatures);

	std::vector<feature> featurePoints(0);
	int nextFeature = 0;

	for (int i = 1; i < numImages; i++) {
		std::string current = padZeroNumber(i);
		std::string curFileName = dirName + "/" + current + ".jpg";
		std::string next = padZeroNumber(i + 1);
		std::string nextFileName = dirName + "/" + next + ".jpg";
		R2Image curImage(curFileName.c_str());
		*this = curImage;
		R2Image other(nextFileName.c_str());
		print("\nMatching %s to %s\n", curFileName.c_str(), nextFileName.c_str());

		// other features are numFeatures of feature points on the second frame

		/*
		The following section should be updated to check for translation only (see MatchTranslation above),
		much of the code from MatchTranslation may be able to be copied right into here just as we did with
		MatchHomography before.
		It should also be updated to save the dx and dy for each frame rather than the feature positions.
		We could just save the dx and dy as a single "track" since a track object can hold a pair of floats for each frame.
		*/
		if (featurePoints.size() < numFeatures) {
			/* Find features in this image */
			R2Image harris(*this);
			harris.Harris(sigma);
			print("Done computing harris\n");


			print("Obtain %d feature points\n", numFeatures);
			while (featurePoints.size() < numFeatures) {
				if (featurePoints.size() % 10 == 0) {
					print("\33[2K\r");
					print("%d features obtained", featurePoints.size());
				}

				// find next highest harris value
				vec2 pixel;
				float harrisVal = 0;
				for (int x = 3 * sigma; x < width - 3 * sigma; x++)	// do not consider edge features
					for (int y = 3 * sigma; y < height - 3 * sigma; y++) {	// do not consider edge features
						float ave = pixelAverage(harris.Pixel(x, y));
						if (ave > harrisVal) {
							harrisVal = ave;
							pixel.x = x;
							pixel.y = y;
						}
					}
				if (!contains(list(featurePoints), pixel)) {
					featurePoints.push_back(feature(pixel, nextFeature));
					nextFeature++;
				}

				// no pixel within 10 pixels should be considered further
				for (int dx = -10; dx <= 10; dx++)
					for (int dy = -10; dy <= 10; dy++)
						harris.Pixel(bound(pixel.x + dx, width), bound(pixel.y + dy, height)) = R2Pixel(0, 0, 0, 1);
			}
			print("\33[2K\r");
			print("%d features obtained\n", featurePoints.size());
			print("Collected all feature points\n\n");
		}

		print("Searching for transformed features\n");
		/* Search for the same features in the other image */
		std::vector<feature> otherFeatures(0);
		for (unsigned int j = 0; j < featurePoints.size(); j++) {
			print("\r%u features found", otherFeatures.size());
			float bestSSD = FLT_MAX;
			vec2 loc(0, 0);
			// search near original feature location
			for (int dx = -width * 0.1; dx < width * 0.1; dx++) {
				for (int dy = -height * 0.1; dy < height * 0.1; dy++) {
					float ssd = getSSD(3 * sigma, 3 * sigma,
						*this, featurePoints[j].point.x, featurePoints[j].point.y, other,
						featurePoints[j].point.x + dx, featurePoints[j].point.y + dy);

					if (ssd <= bestSSD) {
						bestSSD = ssd;
						loc.x = bound(featurePoints[j].point.x + dx, width);
						loc.y = bound(featurePoints[j].point.y + dy, height);
					}
				}
			}

			if (loc.x == -1 || loc.y == -1) {
				print("\nError in SSD calculation");
				exit(1);
			}

			otherFeatures.push_back(feature(loc, featurePoints[j].track));
		}
		print("\r%u features found", otherFeatures.size());

		if (featurePoints.size() != otherFeatures.size()) {
			print("\nError in finding transformed features");
			exit(1);
		}

		// verify features through RANSAC
		print("\nDetermining false matches\n");

		std::vector<unsigned int> paths;

		mat3 bestTransformMatrix(0, 0, 0, 0, 0, 0, 0, 0, 0);

		for (int n = 0; n < numFeatures / 10; n++) {
			std::vector<unsigned int> inlierPaths;

			// randomly select 4 tracks
			std::vector<unsigned int> indecies = getRandoms(4, numFeatures);

			std::vector<fPair> pairs;
			for (int j = 0; j < 4; j++) {
				pairs.push_back(fPair(featurePoints[indecies[j]], otherFeatures[indecies[j]]));
			}


			mat3 transformationMatrix = getHomographyMatrix(toVec2List(pairs));

			// check all features
			for (unsigned int j = 0; j < numFeatures; j++) {
				// project features according to transformation matrix
				vec3 feature(featurePoints[j].point.x, featurePoints[j].point.y, 1);
				vec3 projectedFeatureH = transformationMatrix * feature;

				if (projectedFeatureH.z != 0) {
					vec2 projectedFeature(projectedFeatureH.x / projectedFeatureH.z, projectedFeatureH.y / projectedFeatureH.z);

					// determine distance between actual motion and projected motion
					double dist = (projectedFeature - otherFeatures[j].point).length();

					// if within distance threshold, add to inlierTracks
					if (dist <= 4)
						inlierPaths.push_back(j);
				}
			}
			if (inlierPaths.size() > paths.size()) {
				paths = inlierPaths;
				bestTransformMatrix = transformationMatrix;
			}
		}
		print("Done\n");


		print("Removing false matches...");
		// set features that are not vaild based on RANSAC
		int index = 0;
		for (int j = 0; j < featurePoints.size(); j++)
			if (index < paths.size() && paths[index] == j) {
				index++;
			}
			else
				otherFeatures[j].point.x = -1;
		print("Done\n");

		// first time, add first image's features to list
		if (tracks[0].size() == 0) {
			for (int j = 0; j < featurePoints.size(); j++) {
				tracks[featurePoints[j].track].addPosition(featurePoints[j].point);
			}
		}

		// add other features to list
		for (int j = 0; j < otherFeatures.size(); j++) {
			if(featurePoints[j].track < tracks.size())
				tracks[featurePoints[j].track].addPosition(otherFeatures[j].point);
			else {
				tracks.push_back(track());
				assert(featurePoints[j].track < tracks.size());
				for (int k = 0; k < i; k++)
					tracks[featurePoints[j].track].addPosition(vec2(-1, -1));
				tracks[featurePoints[j].track].addPosition(otherFeatures[j].point);
			}
		}

		// remove bad tracks
		for (int j = 0; j < otherFeatures.size(); j++) {
			if (!tracks[otherFeatures[j].track].isValid()) {
				otherFeatures.erase(otherFeatures.begin() + j);
				j--;
			}
		}

		featurePoints = otherFeatures;
	}
	return tracks;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
	// Initialize everything
	if (pixels) { delete[] pixels; pixels = NULL; }
	npixels = width = height = 0;

	// Parse input filename extension
	char *input_extension;
	if (!(input_extension = (char*)strrchr(filename, '.'))) {
		fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
		return 0;
	}

	// Read file of appropriate type
	if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
	else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
	else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
	else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

	// Should never get here
	fprintf(stderr, "Unrecognized image file extension");
	return 0;
}



int R2Image::
Write(const char *filename) const
{
	// Parse input filename extension
	char *input_extension;
	if (!(input_extension = (char*)strrchr(filename, '.'))) {
		fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
		return 0;
	}

	// Write file of appropriate type
	if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
	else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
	else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
	else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

	// Should never get here
	fprintf(stderr, "Unrecognized image file extension");
	return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
	unsigned short int bfType;
	unsigned int bfSize;
	unsigned short int bfReserved1;
	unsigned short int bfReserved2;
	unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
	unsigned int biSize;
	int biWidth;
	int biHeight;
	unsigned short int biPlanes;
	unsigned short int biBitCount;
	unsigned int biCompression;
	unsigned int biSizeImage;
	int biXPelsPerMeter;
	int biYPelsPerMeter;
	unsigned int biClrUsed;
	unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
	unsigned char rgbtBlue;
	unsigned char rgbtGreen;
	unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
	unsigned char rgbBlue;
	unsigned char rgbGreen;
	unsigned char rgbRed;
	unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
	// Read a unsigned short int from a file in little endian format
	unsigned short int lsb, msb;
	lsb = getc(fp);
	msb = getc(fp);
	return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
	// Write a unsigned short int to a file in little endian format
	unsigned char lsb = (unsigned char)(x & 0x00FF); putc(lsb, fp);
	unsigned char msb = (unsigned char)(x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
	// Read a unsigned int word from a file in little endian format
	unsigned int b1 = getc(fp);
	unsigned int b2 = getc(fp);
	unsigned int b3 = getc(fp);
	unsigned int b4 = getc(fp);
	return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
	// Write a unsigned int to a file in little endian format
	unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
	unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
	unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
	unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
	// Read a int word from a file in little endian format
	int b1 = getc(fp);
	int b2 = getc(fp);
	int b3 = getc(fp);
	int b4 = getc(fp);
	return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
	// Write a int to a file in little endian format
	char b1 = (x & 0x000000FF); putc(b1, fp);
	char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
	char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
	char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	/* Read file header */
	BITMAPFILEHEADER bmfh;
	bmfh.bfType = WordReadLE(fp);
	bmfh.bfSize = DWordReadLE(fp);
	bmfh.bfReserved1 = WordReadLE(fp);
	bmfh.bfReserved2 = WordReadLE(fp);
	bmfh.bfOffBits = DWordReadLE(fp);

	/* Check file header */
	assert(bmfh.bfType == BMP_BF_TYPE);
	/* ignore bmfh.bfSize */
	/* ignore bmfh.bfReserved1 */
	/* ignore bmfh.bfReserved2 */
	assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

	/* Read info header */
	BITMAPINFOHEADER bmih;
	bmih.biSize = DWordReadLE(fp);
	bmih.biWidth = LongReadLE(fp);
	bmih.biHeight = LongReadLE(fp);
	bmih.biPlanes = WordReadLE(fp);
	bmih.biBitCount = WordReadLE(fp);
	bmih.biCompression = DWordReadLE(fp);
	bmih.biSizeImage = DWordReadLE(fp);
	bmih.biXPelsPerMeter = LongReadLE(fp);
	bmih.biYPelsPerMeter = LongReadLE(fp);
	bmih.biClrUsed = DWordReadLE(fp);
	bmih.biClrImportant = DWordReadLE(fp);

	// Check info header
	assert(bmih.biSize == BMP_BI_SIZE);
	assert(bmih.biWidth > 0);
	assert(bmih.biHeight > 0);
	assert(bmih.biPlanes == 1);
	assert(bmih.biBitCount == 24);  /* RGB */
	assert(bmih.biCompression == BI_RGB);   /* RGB */
	int lineLength = bmih.biWidth * 3;  /* RGB */
	if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
	assert(bmih.biSizeImage == (unsigned int)lineLength * (unsigned int)bmih.biHeight);

	// Assign width, height, and number of pixels
	width = bmih.biWidth;
	height = bmih.biHeight;
	npixels = width * height;

	// Allocate unsigned char buffer for reading pixels
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = bmih.biSizeImage;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Read buffer
	fseek(fp, (long)bmfh.bfOffBits, SEEK_SET);
	if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
		fprintf(stderr, "Error while reading BMP file %s", filename);
		return 0;
	}

	// Close file
	fclose(fp);

	// Allocate pixels for image
	pixels = new R2Pixel[width * height];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Assign pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			double b = (double)*(p++) / 255;
			double g = (double)*(p++) / 255;
			double r = (double)*(p++) / 255;
			R2Pixel pixel(r, g, b, 1);
			SetPixel(i, j, pixel);
		}
	}

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return success
	return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
	// Open file
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	// Compute number of bytes in row
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

	// Write file header
	BITMAPFILEHEADER bmfh;
	bmfh.bfType = BMP_BF_TYPE;
	bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
	bmfh.bfReserved1 = 0;
	bmfh.bfReserved2 = 0;
	bmfh.bfOffBits = BMP_BF_OFF_BITS;
	WordWriteLE(bmfh.bfType, fp);
	DWordWriteLE(bmfh.bfSize, fp);
	WordWriteLE(bmfh.bfReserved1, fp);
	WordWriteLE(bmfh.bfReserved2, fp);
	DWordWriteLE(bmfh.bfOffBits, fp);

	// Write info header
	BITMAPINFOHEADER bmih;
	bmih.biSize = BMP_BI_SIZE;
	bmih.biWidth = width;
	bmih.biHeight = height;
	bmih.biPlanes = 1;
	bmih.biBitCount = 24;       /* RGB */
	bmih.biCompression = BI_RGB;    /* RGB */
	bmih.biSizeImage = rowsize * (unsigned int)bmih.biHeight;  /* RGB */
	bmih.biXPelsPerMeter = 2925;
	bmih.biYPelsPerMeter = 2925;
	bmih.biClrUsed = 0;
	bmih.biClrImportant = 0;
	DWordWriteLE(bmih.biSize, fp);
	LongWriteLE(bmih.biWidth, fp);
	LongWriteLE(bmih.biHeight, fp);
	WordWriteLE(bmih.biPlanes, fp);
	WordWriteLE(bmih.biBitCount, fp);
	DWordWriteLE(bmih.biCompression, fp);
	DWordWriteLE(bmih.biSizeImage, fp);
	LongWriteLE(bmih.biXPelsPerMeter, fp);
	LongWriteLE(bmih.biYPelsPerMeter, fp);
	DWordWriteLE(bmih.biClrUsed, fp);
	DWordWriteLE(bmih.biClrImportant, fp);

	// Write image, swapping blue and red in each pixel
	int pad = rowsize - width * 3;
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			const R2Pixel& pixel = (*this)[i][j];
			double r = 255.0 * pixel.Red();
			double g = 255.0 * pixel.Green();
			double b = 255.0 * pixel.Blue();
			if (r >= 255) r = 255;
			if (g >= 255) g = 255;
			if (b >= 255) b = 255;
			fputc((unsigned char)b, fp);
			fputc((unsigned char)g, fp);
			fputc((unsigned char)r, fp);
		}

		// Pad row
		for (int i = 0; i < pad; i++) fputc(0, fp);
	}

	// Close file
	fclose(fp);

	// Return success
	return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	// Read PPM file magic identifier
	char buffer[128];
	if (!fgets(buffer, 128, fp)) {
		fprintf(stderr, "Unable to read magic id in PPM file");
		fclose(fp);
		return 0;
	}

	// skip comments
	int c = getc(fp);
	while (c == '#') {
		while (c != '\n') c = getc(fp);
		c = getc(fp);
	}
	ungetc(c, fp);

	// Read width and height
	if (fscanf(fp, "%d%d", &width, &height) != 2) {
		fprintf(stderr, "Unable to read width and height in PPM file");
		fclose(fp);
		return 0;
	}

	// Read max value
	double max_value;
	if (fscanf(fp, "%lf", &max_value) != 1) {
		fprintf(stderr, "Unable to read max_value in PPM file");
		fclose(fp);
		return 0;
	}

	// Allocate image pixels
	pixels = new R2Pixel[width * height];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for PPM file");
		fclose(fp);
		return 0;
	}

	// Check if raw or ascii file
	if (!strcmp(buffer, "P6\n")) {
		// Read up to one character of whitespace (\n) after max_value
		int c = getc(fp);
		if (!isspace(c)) putc(c, fp);

		// Read raw image data
		// First ppm pixel is top-left, so read in opposite scan-line order
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				double r = (double)getc(fp) / max_value;
				double g = (double)getc(fp) / max_value;
				double b = (double)getc(fp) / max_value;
				R2Pixel pixel(r, g, b, 1);
				SetPixel(i, j, pixel);
			}
		}
	}
	else {
		// Read asci image data
		// First ppm pixel is top-left, so read in opposite scan-line order
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				// Read pixel values
				int red, green, blue;
				if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
					fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
					fclose(fp);
					return 0;
				}

				// Assign pixel values
				double r = (double)red / max_value;
				double g = (double)green / max_value;
				double b = (double)blue / max_value;
				R2Pixel pixel(r, g, b, 1);
				SetPixel(i, j, pixel);
			}
		}
	}

	// Close file
	fclose(fp);

	// Return success
	return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
	// Check type
	if (ascii) {
		// Open file
		FILE *fp = fopen(filename, "w");
		if (!fp) {
			fprintf(stderr, "Unable to open image file: %s", filename);
			return 0;
		}

		// Print PPM image file
		// First ppm pixel is top-left, so write in opposite scan-line order
		fprintf(fp, "P3\n");
		fprintf(fp, "%d %d\n", width, height);
		fprintf(fp, "255\n");
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				const R2Pixel& p = (*this)[i][j];
				int r = (int)(255 * p.Red());
				int g = (int)(255 * p.Green());
				int b = (int)(255 * p.Blue());
				fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
				if (((i + 1) % 4) == 0) fprintf(fp, "\n");
			}
			if ((width % 4) != 0) fprintf(fp, "\n");
		}
		fprintf(fp, "\n");

		// Close file
		fclose(fp);
	}
	else {
		// Open file
		FILE *fp = fopen(filename, "wb");
		if (!fp) {
			fprintf(stderr, "Unable to open image file: %s", filename);
			return 0;
		}

		// Print PPM image file
		// First ppm pixel is top-left, so write in opposite scan-line order
		fprintf(fp, "P6\n");
		fprintf(fp, "%d %d\n", width, height);
		fprintf(fp, "255\n");
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				const R2Pixel& p = (*this)[i][j];
				int r = (int)(255 * p.Red());
				int g = (int)(255 * p.Green());
				int b = (int)(255 * p.Blue());
				fprintf(fp, "%c%c%c", r, g, b);
			}
		}

		// Close file
		fclose(fp);
	}

	// Return success
	return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
extern "C" {
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
};
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	// Initialize decompression info
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, fp);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);

	// Remember image attributes
	width = cinfo.output_width;
	height = cinfo.output_height;
	npixels = width * height;
	int ncomponents = cinfo.output_components;

	// Allocate pixels for image
	pixels = new R2Pixel[npixels];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Allocate unsigned char buffer for reading image
	int rowsize = ncomponents * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = rowsize * height;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
		fclose(fp);
		return 0;
	}

	// Read scan lines
	// First jpeg pixel is top-left, so read pixels in opposite scan-line order
	while (cinfo.output_scanline < cinfo.output_height) {
		int scanline = cinfo.output_height - cinfo.output_scanline - 1;
		unsigned char *row_pointer = &buffer[scanline * rowsize];
		jpeg_read_scanlines(&cinfo, &row_pointer, 1);
	}

	// Free everything
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);

	// Close file
	fclose(fp);

	// Assign pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			double r, g, b, a;
			if (ncomponents == 1) {
				r = g = b = (double)*(p++) / 255;
				a = 1;
			}
			else if (ncomponents == 1) {
				r = g = b = (double)*(p++) / 255;
				a = 1;
				p++;
			}
			else if (ncomponents == 3) {
				r = (double)*(p++) / 255;
				g = (double)*(p++) / 255;
				b = (double)*(p++) / 255;
				a = 1;
			}
			else if (ncomponents == 4) {
				r = (double)*(p++) / 255;
				g = (double)*(p++) / 255;
				b = (double)*(p++) / 255;
				a = (double)*(p++) / 255;
			}
			else {
				fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
				return 0;
			}
			R2Pixel pixel(r, g, b, a);
			SetPixel(i, j, pixel);
		}
	}

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return success
	return 1;
#else
	fprintf(stderr, "JPEG not supported");
	return 0;
#endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
	// Open file
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s", filename);
		return 0;
	}

	// Initialize compression info
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, fp);
	cinfo.image_width = width; 	/* image width and height, in pixels */
	cinfo.image_height = height;
	cinfo.input_components = 3;		/* # of color components per pixel */
	cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
	cinfo.dct_method = JDCT_ISLOW;
	jpeg_set_defaults(&cinfo);
	cinfo.optimize_coding = TRUE;
	jpeg_set_quality(&cinfo, 95, TRUE);
	jpeg_start_compress(&cinfo, TRUE);

	// Allocate unsigned char buffer for reading image
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = rowsize * height;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
		fclose(fp);
		return 0;
	}

	// Fill buffer with pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			const R2Pixel& pixel = (*this)[i][j];
			int r = (int)(255 * pixel.Red());
			int g = (int)(255 * pixel.Green());
			int b = (int)(255 * pixel.Blue());
			if (r > 255) r = 255;
			if (g > 255) g = 255;
			if (b > 255) b = 255;
			*(p++) = r;
			*(p++) = g;
			*(p++) = b;
		}
	}



	// Output scan lines
	// First jpeg pixel is top-left, so write in opposite scan-line order
	while (cinfo.next_scanline < cinfo.image_height) {
		int scanline = cinfo.image_height - cinfo.next_scanline - 1;
		unsigned char *row_pointer = &buffer[scanline * rowsize];
		jpeg_write_scanlines(&cinfo, &row_pointer, 1);
	}

	// Free everything
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	// Close file
	fclose(fp);

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return number of bytes written
	return 1;
#else
	fprintf(stderr, "JPEG not supported");
	return 0;
#endif
}
