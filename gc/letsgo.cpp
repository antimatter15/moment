//////////////////////////////////////////////////////////////////////////////
// wolo
/////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "GCoptimization.h"



#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


using namespace std;

const int VID_WIDTH = 299;
const int VID_HEIGHT = 168;
const int VID_FRAMES = 40;
const int CHANNELS = 3;
const int CSTATIC = .5;

//all possible start times
int start_time_at_index(int i){
	return i*4;
}

int start_time_index(int start_time){
	return start_time/4;
}

//all possible periods
int period_at_index(int i){
	return (i+1)*4;
}

int period_index(int period){
	return period/4-1;
}


int stage1_cur_period = 12;

//video frame=time, height=rows=ycoords, width=cols=xcoords, channels=rgb
int video [VID_FRAMES][VID_HEIGHT][VID_WIDTH][CHANNELS];

//[pixel][period] = optimal start time
int biglookupthing [VID_WIDTH * VID_HEIGHT][VID_FRAMES/4];

void setOptimalStartTimeIndexForPeriodIndex(int pixel, int period_index, int start_time_index){
	biglookupthing[pixel][period_index] = start_time_index;
}

int getOptimalStartTimeForPeriodIndex(int pixel, int period_index){
	return start_time_at_index(biglookupthing[pixel][period_index]);
}

struct ForDataFn{
	int numLab;
	int *data;
};


int smoothFn(int p1, int p2, int l1, int l2)
{
	if ( (l1-l2)*(l1-l2) <= 4 ) return((l1-l2)*(l1-l2));
	else return(4);
}

int dataFn(int p, int l, void *data)
{
	ForDataFn *myData = (ForDataFn *) data;
	int numLab = myData->numLab;
	
	return( myData->data[p*numLab+l] );
}

int getTemporalChromaticNorm(int x, int y, int t1, int t2){

	// printf("Got here...\n");
	// printf("%d %d\n", (t1 + VID_FRAMES) % VID_FRAMES, (t2 + VID_FRAMES) % VID_FRAMES);

	int r = pow(video[(t1 + VID_FRAMES) % VID_FRAMES][y][x][0] - video[(t2 + VID_FRAMES) % VID_FRAMES][y][x][0],2);
	int g = pow(video[(t1 + VID_FRAMES) % VID_FRAMES][y][x][1] - video[(t2 + VID_FRAMES) % VID_FRAMES][y][x][1],2);
	int b = pow(video[(t1 + VID_FRAMES) % VID_FRAMES][y][x][2] - video[(t2 + VID_FRAMES) % VID_FRAMES][y][x][2],2);

	int full = r+g+b;


	// printf("full: %d = %d+%d+%d \n", full,r,g,b);

	return full;
}

int temporalNeighborDiff(int x, int y, int t){
	//return threenorm(video[x[0]][x[1]][t], video[x[0]][x[1]][t+1]); 
	return (int) round(sqrt(getTemporalChromaticNorm(x, y, t, t+1)));
}
 
void printVec(vector<int>& vec){
	int size =  vec.size();
	printf("size %d", size);
	// cout<<"\nThe given array is ... "<<endl;
    for (int i = 0; i < size; i++) {
    	printf("ded tyehsd %d %d\n", vec[i], i);
        // cout<<vec[i]<<" ";
    }
}

int simpleFuckingMedian(vector<int>& list) {
    std::sort (list.begin(), list.end());
    return list[list.size() / 2];
}

int MutuallyAssuredDestruction(vector<int>& list) {
	int size = list.size();
	int median = simpleFuckingMedian(list);
	vector<int> deviations(size);
	for(int i = 0; i < size; i++){
		deviations[i] = abs(list[i] - median);
	}
	return simpleFuckingMedian(deviations);
}



double gaussianLookup(double distanceSquared, double sigma){
	return 1.0 / (sigma * sqrt(2.0 * 3.14159)) * exp(-distanceSquared / (2 * sigma * sigma));
}

int getIntensityAtPosition(int x, int y, int t){
	t = (t + VID_FRAMES) % VID_FRAMES;
	x = max(0, min(VID_WIDTH, x));
	y = max(0, min(VID_HEIGHT, y));
	return (video[t][y][x][0] * 212 + video[t][y][x][1] *  715 + video[t][y][x][2] * 72) / 1000;
}


int gaussianDifferenceNorm(int x0, int y0, int t0){
	int spatialRadius = 10;
	int temporalRadius = 10;
	int totesMcGoats = 0;
	for(int dx = -spatialRadius; dx < spatialRadius; dx++){
		for(int dy = -spatialRadius; dy < spatialRadius; dy++){
			for(int dt = -temporalRadius; dt < temporalRadius; dt++){
				totesMcGoats +=
					gaussianLookup(dt*dt,1.2)*
					pow(
						gaussianLookup(dx * dx + dy * dy + dt * dt, 0.9)
							* getIntensityAtPosition(x0 + dx, y0 + dy, t0 + dt) -
						gaussianLookup(dx * dx + dy * dy + (dt + 1) * (dt + 1), 0.9)
							* getIntensityAtPosition(x0 + dy, y0 + dy, t0 + dt + 1)
					, 2);
			}
		}
	}
	return sqrt(totesMcGoats);
}



int staticEnergy(int x, int y){
	vector<int> spatioTemporalGaussianNorms(VID_FRAMES);
	for (int t = 0; t < VID_FRAMES; ++t)
	{
		spatioTemporalGaussianNorms[t] = gaussianDifferenceNorm(x, y, t);
	}

	return CSTATIC*min(1,100*MutuallyAssuredDestruction(spatioTemporalGaussianNorms));
}






//temporal energy at this pixel w/ this period + start_time
int getJiggawatts(int x, int y, int start_time, int period){
	// printVec(pos);

	int staticE = 0;
	if (period==1)
	{
		staticE = staticEnergy(x,y);
	}

	vector<int> neighbors(period);

	for (int t = 0; t < period; t++){
		// printf("getjigg: getting temporalNeighborDiff: %d, %d, %d \n",x, y, start_time + t);
		neighbors[t] = temporalNeighborDiff(x, y, start_time + t);
		// printf("getjigg: temporalNeighborDiff: %d \n", neighbors[t]);
	}

	// printVec(neighbors);
	int magic = 1 + 400 * MutuallyAssuredDestruction(neighbors);
	// printf("magic: %d\n",magic);

	return staticE+(getTemporalChromaticNorm(x, y, start_time    , start_time + period    ) +
			getTemporalChromaticNorm(x, y, start_time - 1, start_time + period - 1)) / magic;
}

// http://stackoverflow.com/a/4229930/205784
int gcd(int a, int b){
    for (;;){
        if (a == 0) return b;
        b %= a;
        if (b == 0) return a;
        a %= b;
    }
}

int lcm(int a, int b){
    int temp = gcd(a, b);
    return temp ? (a / temp * b) : 0;
}

//time start time pixel number
int fie(int t, int sx, int px) {
       int phi = sx - (sx % px) + (t % px);
       if( (t % px) < (sx % px)) {
               phi += px;
       }
       return phi;
}



int getSpatialChromaticDiff(int x1, int y1, int x2, int y2, int t){

	return (int) sqrt( pow(video[(t + VID_FRAMES) % VID_FRAMES][y1][x1][0] - video[(t + VID_FRAMES) % VID_FRAMES][y2][x2][0],2)+
				     + pow(video[(t + VID_FRAMES) % VID_FRAMES][y1][x1][1] - video[(t + VID_FRAMES) % VID_FRAMES][y2][x2][1],2)+
				     + pow(video[(t + VID_FRAMES) % VID_FRAMES][y1][x1][2] - video[(t + VID_FRAMES) % VID_FRAMES][y2][x2][2],2));
}


int invlambda(int x1, int y1, int x2, int y2, int T) {

	vector<int> spatialAdjacentPairs(T);

	for (int t = 0; t < T; t++){
		spatialAdjacentPairs[t] = getSpatialChromaticDiff(x1, y1, x2, y2, t);
	}

	return 1 + 100 * MutuallyAssuredDestruction( spatialAdjacentPairs );
}

//xy, start time, period
int muskEnergy(int x1, int y1, int x2, int y2, int s1, int s2, int p1, int p2){
	int psi = 0;
	int T = lcm(p1, p2);
	for(int t = 0; t < T; t++){
		int fi_1 = fie(t, s1, p1),
			fi_2 = fie(t, s2, p2);
		psi+=getTemporalChromaticNorm(x1, y1,fi_1,fi_2);
		psi+=getTemporalChromaticNorm(x2, y2,fi_1,fi_2);
	}

	// printf("psi: %d T: %d\n", psi, T);
	psi = psi / T; // todo non-int arith
	// printf("\nmusk %d, %d\n", psi, invlambda(x1, y1, x2, y2, T));
	return psi / invlambda(x1, y1, x2, y2, T);
}



int muskEnergyWrapperGivenPeriod(int p1, int p2, int l1, int l2){
	return muskEnergy( p1 % VID_WIDTH, p1/VID_WIDTH,  p2 % VID_WIDTH, p2/VID_WIDTH, l1, l2, stage1_cur_period, stage1_cur_period);
}



////////////////////////////////////////////////////////////////////////////////
//given a global period, what are the best start times?
//executed in a for loop over each period
void RockNRollGivenPeriod(int width, int height, int num_pixels, int period)
{
	// printf("got here");

	stage1_cur_period = period;


	int num_labels = VID_FRAMES/4;

	int *result = new int[num_pixels];   // stores result of optimization

	// first set up the array for data costs
	int *data = new int[num_pixels*num_labels];
	// printf("data\n");
	for ( int i = 0; i < num_pixels; i++ ){
		for (int l = 0; l < num_labels; l++ ){ //l is the index of the start time

			// printf("getting jiggawatts %d, %d, %d, %d\n", i/width, i % width, l*4, stage1_cur_period);
			data[i*num_labels+l] = getJiggawatts(i/width, i % width, start_time_at_index(l), period);
			// printf("%d: %d \n", i * num_labels, data[i*num_labels+l]);
		}
	}

	printf("calculated jiggawatts\n");

	try{
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(width,height,num_labels);

		// set up the needed data to pass to function for the data costs
		ForDataFn toFn;
		toFn.data = data;
		toFn.numLab = num_labels;
		printf("setting data cost\n");
		gc->setDataCost(&dataFn,&toFn);

		printf("setting smooth cost\n");
		// smoothness comes from function pointer
		gc->setSmoothCost(&muskEnergyWrapperGivenPeriod);
		printf("smoove\n");
		printf("\nBefore optimization energy is %d\n",gc->compute_energy());
		printf("doing expansion\n");
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		printf("has optimized\n");
		printf("\nAfter optimization energy is %d\n",gc->compute_energy());

		for ( int  pixel = 0; pixel < num_pixels; pixel++ ){
			result[pixel] = gc->whatLabel(pixel);
			setOptimalStartTimeIndexForPeriodIndex(pixel, period_index(period), result[pixel]);
			// printf(" %d", result[i]);
		}

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

	delete [] result;
	delete [] data;

}


int muskEnergyWrapperStage2(int p1, int p2, int l1, int l2){
	// printf("called musk energy stage 2\n");

	return muskEnergy( p1 % VID_WIDTH, p1/VID_WIDTH,  p2 % VID_WIDTH, p2/VID_WIDTH, getOptimalStartTimeForPeriodIndex(p1,l1), getOptimalStartTimeForPeriodIndex(p2,l2), period_at_index(l1), period_at_index(l2));
}

////////////////////////////////////////////////////////////////////////////////
//given the best start times for each pixel for each period, what are the best per pixel periods?
void RockNRollGivenStartTimes(int width,int height,int num_pixels)
{

	int num_labels = VID_FRAMES/4-1;

	int *result = new int[num_pixels];   // stores result of optimization

	// first set up the array for data costs
	// l is the period of the pixel
	printf("starting rnr given start times\n");

	int *data = new int[num_pixels*num_labels];
	for ( int i = 0; i < num_pixels; i++ ){
		for (int l = 0; l < num_labels; l++ ){ //l is the period in the optimal period / start time pair
			// std::vector<int> pos(2);
			// int a[] = {num_pixels/width, num_pixels % width};
			// vector<int> pos(a, a + 2);
			// pos.push_back(num_pixels/width);
			// pos.push_back(num_pixels % width);
			data[i*num_labels+l] = getJiggawatts(num_pixels/width, num_pixels % width, getOptimalStartTimeForPeriodIndex(i,l), period_at_index(l));
		}
	}
	printf("finished populating data array\n");

	try{
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(width,height,num_labels);

		// set up the needed data to pass to function for the data costs
		ForDataFn toFn;
		toFn.data = data;
		toFn.numLab = num_labels;

		gc->setDataCost(&dataFn,&toFn);

		// smoothness comes from function pointer
		gc->setSmoothCost(&muskEnergyWrapperStage2);

		printf("\nBefore _final_ optimization energy is %d",gc->compute_energy());
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		printf("\nAfter optimization energy is %d",gc->compute_energy());

		for ( int  i = 0; i < num_pixels; i++ ){
			result[i] = gc->whatLabel(i);
			printf("%d ", result[i]*4);			
		}

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

	delete [] result;
	delete [] data;

}


int main(int argc, char **argv)
{


	// int width = 10;
	// int height = 5;
	// int num_pixels = width*height;
	// int num_labels = 7;
	printf("starting\n");

	FILE *fileptr;
    unsigned char *buffer;
    long filelen;

    fileptr = fopen("../stab/nine.bin", "rb");  // Open the file in binary mode
    fseek(fileptr, 0, SEEK_END);          // Jump to the end of the file
    filelen = ftell(fileptr);             // Get the current byte offset in the file
    rewind(fileptr);                      // Jump back to the beginning of the file

    buffer = (unsigned char *)malloc((filelen+1)*sizeof(unsigned char)); // Enough memory for file + \0
    fread(buffer, filelen, 1, fileptr); // Read in the entire file
    fclose(fileptr); // Close the file

    printf("length of file %d %d\n", filelen, VID_WIDTH * VID_HEIGHT * VID_FRAMES * CHANNELS);
    int pos = 0;
	for(int i = 0; i < VID_FRAMES; i++){
		for(int j = 0; j < VID_HEIGHT; j++){
			for(int k = 0; k < VID_WIDTH; k++){
				for(int l = 0; l < CHANNELS; l++){
					// read one byte and stick it in 
					video[i][j][k][l] = buffer[pos];
					if(pos<10){
						printf("yo %d\n", buffer[pos]);
					}
					pos++;
				}
			}
		}
	}
	printf("finished reading file\n");

	for (int i = 4; i < VID_FRAMES; i += 4) {
		//fill in the appropriate entries in biglookup
		printf("started %d %d\n",i,VID_FRAMES);
		RockNRollGivenPeriod(VID_WIDTH, VID_HEIGHT, VID_WIDTH*VID_HEIGHT, i);
		printf("finished %d %d\n\n",i,VID_FRAMES);
	}

	printf("Stage 1 complete\n");

	RockNRollGivenStartTimes(VID_WIDTH, VID_HEIGHT, VID_WIDTH*VID_HEIGHT);

	printf("\n  Finished %d (%d) clock per sec %d\n",clock()/CLOCKS_PER_SEC,clock(),CLOCKS_PER_SEC);

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////