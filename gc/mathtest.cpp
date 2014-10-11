#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

int findMedian(vector<int> vec){
	int median;
    size_t size = vec.size();
    median = vec[(size/2)];
    return median;
}

int findMedianOfMedians(vector<vector<int> > values){
    vector<int> medians;
    for (int i = 0; i < values.size(); i++) {
        int m = findMedian(values[i]);
        medians.push_back(m);
    }
    return findMedian(medians);
}

int selectionByMedianOfMedians(const vector<int> values, int k){
	//    Divide the list into n/5 lists of 5 elements each
    vector<vector<int> > vec2D;
    int count = 0;
    while (count != values.size()) {
        int countRow = 0;
        vector<int> row;
        while ((countRow < 5) && (count < values.size())) {
            row.push_back(values[count]);
            count++;
            countRow++;
        }
        vec2D.push_back(row);
    }
    int m = findMedianOfMedians(vec2D);
	//    Partition the list into unique elements larger than 'm' (call this sublist L1) and
	//    those smaller them 'm' (call this sublist L2)
    vector<int> L1, L2;
    for (int i = 0; i < vec2D.size(); i++) {
        for (int j = 0; j < vec2D[i].size(); j++) {
            if (vec2D[i][j] > m) {
                L1.push_back(vec2D[i][j]);
            }else if (vec2D[i][j] < m){
                L2.push_back(vec2D[i][j]);
            }
        }
    }
	if ((k - 1) == L1.size()) {
        return m;
    }else if (k <= L1.size()) {
        return selectionByMedianOfMedians(L1, k);
    }else if (k > (L1.size() + 1)){
        return selectionByMedianOfMedians(L2, k-((int)L1.size())-1);
    }
}



int median(int intlist[]){
	return intlist[sizeof(intlist)/sizeof(*intlist)/2];
}

int mad(int intlist[]){
	int med = median(intlist);
	int length = sizeof(intlist)/sizeof(*intlist);
	int* absolute_deviations = new int[length];
	for(int i = 0;i<length; i++){
		absolute_deviations[i] = abs(intlist[i] - med);
	}
	return median(absolute_deviations);
}

int threeselfdot(int v1[], int v2[]){
	return pow(v1[0] - v2[0],2) 
	     + pow(v1[1] - v2[1],2)
	     + pow(v1[2] - v2[2],2);
}

int threenorm(int v1[], int v2[]){
	return (int) sqrt(threeselfdot(v1,v2));	
}

// http://stackoverflow.com/a/4229930/205784
int gcd(int a, int b)
{
    for (;;)
    {
        if (a == 0) return b;
        b %= a;
        if (b == 0) return a;
        a %= b;
    }
}

int lcm(int a, int b)
{
    int temp = gcd(a, b);
    return temp ? (a / temp * b) : 0;
}



void printVec(vector<int>& vec){
	cout<<"\nThe given array is : "<<endl;
    for (int i = 0; i < vec.size(); i++) {
        cout<<vec[i]<<" ";
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


int main(int argc, char **argv){

    

    // char blah[][10] = (char[][10])buffer;



	int list[] = {1, 2, 3, 4, 5, 6, 7, 8, 8, 3, 2};
	int len = sizeof(list)/sizeof(*list);
	printf("\n  length %d", len);

	printf("\n  median of stuff %d", median(list));
	printf("\n  MAD %d", median(list));
	printf("\n  median of stuff %d", median(list));

	printf("\n  hashtag yolo");

	int values[] = {2, 3, 5, 4, 1, 12, 11, 13, 16, 7, 8, 6, 10, 9, 17, 15, 19, 20, 18, 23, 21, 22, 25, 24, 14, 2, 3, 5, 4, 1, 12, 11, 13, 16, 7, 8, 6, 10, 9, 17, 15, 19, 20, 18, 23, 21, 22, 25, 24, 14, 2, 3, 5, 4, 1, 12, 11, 13, 16, 7, 8, 6, 10, 9, 17, 15, 19, 20, 18, 23, 21, 22, 25, 24, 14, 2, 3, 5, 4, 1, 12, 11, 13, 16, 7, 8, 6, 10, 9, 17, 15, 19, 20, 18, 23, 21, 22, 25, 24, 14};


    vector<int> vec(values, values + 43);

    printVec(vec);

    printf("\n CALCULATING MEDIAN OF ABSOLUTE DEVIATIONS %d \n", MutuallyAssuredDestruction(vec));

	printVec(vec);    


    // cout<<"The given array is : "<<endl;
    // for (int i = 0; i < vec.size(); i++) {
    //     cout<<vec[i]<<" ";
    // }

    // selectionByMedianOfMedians(vec, 8);
}