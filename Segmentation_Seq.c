#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "cv.h"
#include "highgui.h"
//####################Global Variables######################
IplImage* orgImg ;
//###################Declarations###########################
void buildWeightsMatrix(int*,int,int); 
int calColumn(int,int,int);
void calWeight(int,int,int,int[],int,int*);
void constructSegmentedImage(int*,int,int);
int* countRegions(int*,int);
int* countSort(int,int*,int);
int getGrayLevel(int,int);
void initMatrix(int*,int,int);
int indexOf(int,int[]);
int regionOf(int,int*);
int* segmentation(int*,int,int,int);
void unionRegions(int,int,int*);
//###################Functions##############################
//Convert the image matrix to a weighted graph matrix
void buildWeightsMatrix(int*weightsMatrix,int width,int height)
{		
    int i,j;
    int mainPix;
    int pixNum=height*width;
    int matrixWidth=4;
    int max=256;
    int row=0;
	
    initMatrix(weightsMatrix,matrixWidth*pixNum,max);
	
  	for(i = 0;i < height; i++)
	{
        for(j = 0;j < width; j++)
		{
			mainPix=getGrayLevel(i,j);
			if(i==0)
			{
				if(j==0)
				{
					int skip[]={1,-1};
					calWeight(i,j,row,skip,mainPix,weightsMatrix);
				}
				else
					if(j==width-1)
					{
						int skip[]={0,3,-1};
						calWeight(i,j,row,skip,mainPix,weightsMatrix);
					}
					else
					{
						int skip[]={-1};
						calWeight(i,j,row,skip,mainPix,weightsMatrix);
					}
			}
			else
				if(i==height-1)
				{
					if(j==0)
					{
						int skip[]={1,2,3,-1};
						calWeight(i,j,row,skip,mainPix,weightsMatrix);
					}
					else
						if(j==width-1)
						{
							int skip[]={0,1,2,3,-1};
							calWeight(i,j,row,skip,mainPix,weightsMatrix);
						}
						else
						{
							int skip[]={1,2,3,-1};
							calWeight(i,j,row,skip,mainPix,weightsMatrix);
						}
				}	
				else
					if(j==0)
					{
						int skip[]={1,-1};
						calWeight(i,j,row,skip,mainPix,weightsMatrix);
					}
					else
						if(j==width-1)
						{	
							int skip[]={0,3,-1};
							calWeight(i,j,row,skip,mainPix,weightsMatrix);
						}
					else
					{
						int skip[]={-1};
						calWeight(i,j,row,skip,mainPix,weightsMatrix);
					}
			row++;
		}		
	}
}
//convert j from weightsMatrix actual indecies to pixNum
int calColumn(int j,int i,int width)
{
	if(j==0)
		return i+1;
	if(j==1)
		return i+width-1;
	if(j==2)
		return i+width;
	return i+width+1;	
}
//Calculate the weight of the edges between the main pixel & its next 4 neighbour
void calWeight(int i,int j,int row,int skip[],int mainPix,int*weightsMatrix) 
{
	int col;
	int matrixWidth=4;
	int index=row*matrixWidth;
	for(col=0;col<matrixWidth;col++)
	{
		if(col==0)
			j++;
		if(col==1)
		{
			i++;
			j=j-2;
		}
		if((indexOf(col,skip))>=0)
			continue;
		index=index+col;
		weightsMatrix[index]=abs(mainPix-getGrayLevel(i,j));
		j++;
	}
}
//Construct & display the final segmented image
void constructSegmentedImage(int* regions,int width,int height)
{
    CvSize imgSize; 
	imgSize.width = width; 
    imgSize.height = height;
    IplImage* segImg = cvCreateImage( imgSize, 8, 1 );
	char* imgName;
	int i,j,index;
    int pixs=width*height;  
    int* count=countRegions(regions,pixs);
    
    for(index=0;index<pixs;index++)
		if(count[index]>1000)
		{
			for ( i = 0; i <height ; i++ )
			{
				for ( j = 0; j < width; j++ )
				{
					if(regions[i*width+j]==index)
						((uchar*)(segImg->imageData + i * segImg->widthStep)) [ j ]=getGrayLevel(i,j);
					//else
						//((uchar*)(segImg->imageData + i * segImg->widthStep)) [ j ]=255;
				}
			}
			char* imgPath=malloc(1000*sizeof(char));
			char imgNum[100];
			char imgFldr[]="./SEQ/seg";
			sprintf(imgNum,"%d",index);
			strcpy(imgPath,imgFldr);
			strcat(imgPath,imgNum);
			strcat(imgPath,".jpg");
			cvSaveImage(imgPath,segImg,0);
			cvSet(segImg,cvScalar(256,256,256,256),NULL);
			free(imgPath);
			//break;
	}
	free(count);
	cvReleaseImage(&segImg);
}
//count the size of each region
int* countRegions(int* regions,int pixs)
{
	int i;
	int* count=malloc(pixs*sizeof(int));
		
	initMatrix(count,pixs,0);
	
	for(i = 0 ; i < pixs ; i++)
		count [ regions [ i ] ] = count [ regions[ i ] ] + 1;
		
	return count;
}
//Counting Sort Algorithm to sort the edges according to their weights
int* countSort(int max,int* matrix,int size)
{
	int* sortedIndexes=malloc(size*sizeof(int));
	int* count=malloc((max)*sizeof(int));
	int i;
	
	initMatrix(count,max+1,0);
   
	for(i = 0 ; i < size ; i++)
		count [ matrix [ i ] ] = count [ matrix[ i ] ] + 1;
   	  
	for( i = 1 ; i < max+1 ; i++)
		count[ i ] = count[ i ] + count[ i - 1];  
   
	for(i = size-1;i>=0;i--)
	{
		sortedIndexes [count[matrix[i]]-1] = i ;
		count[matrix[i]]  = count[matrix[i]]-1;
	}
	
	free(count);
	return sortedIndexes;
}
//get the gray level of a pixel
int getGrayLevel(int i,int j)
{
	if( orgImg == NULL ) 
        fprintf( stderr, "Cannot load file %s!\n", "./images/T.jpg" );
	if(j==-1)
		return ((int)((uchar*)(orgImg->imageData))[i]);
	else
		return ((int)((uchar*)(orgImg->imageData + i*orgImg->widthStep))[j]);
}
//init a matrix with specific value
void initMatrix(int* matrix,int size,int initValue)
{
	int i;
	for(i = 0; i < size; i++)
    {
    	if(initValue==-1)
    		matrix[i] = i;
    	else
        	matrix[i] = initValue;      
    }
}
//check if the num is in the array & return its index if exisits
int indexOf(int num,int arr[])
{
	int i=0;
	while((arr[i]!=-1))
	{			
		if(arr[i]==num)
			return i;
		i++;
	}
	return -1;
}
//get the set containing v
int regionOf(int v,int* regions)
{
    while(regions[v] != v)
       v = regions[v];        
	return v;
}
//Segmentation part
int* segmentation(int* weightsMatrix,int size,int width,int max)
{
    int  i,j,ri,rj,r,c,index;
	int matrixWidth=4;
    int* regions=malloc(size*sizeof(int));
    int* sortedIndexes=countSort(max,weightsMatrix,size*matrixWidth); 
   	 
   	initMatrix(regions,size,-1);
       
  	for(index = 0; index < size*matrixWidth; index++)
    {
    	r=sortedIndexes[index]/matrixWidth;
    	c=sortedIndexes[index]%matrixWidth;
    	if(weightsMatrix[r*matrixWidth+c]==max) 
    		break;      
	    i = r;
	    j = calColumn(c,r,width); 
	    if(j>=size)
	    	continue;
	    ri = regionOf(i, regions);
    	rj = regionOf(j, regions);
    	if ((ri != rj))
	        unionRegions(ri,rj,regions);  
    }
    free(sortedIndexes);
    return regions;
}
//Union different regions
void unionRegions(int i,int j,int* regions)
{
    if(j > i)
        regions[j]= i;   
    else
        regions[i]= j;
}
//////////////////////////////////////////////////////////////////////// 
int main()
{
	int width,height;
    int pixNum;
	int max=256;
	int matrixWidth=4;
	int* regions;
	int* weightsMatrix ;
	
    orgImg = cvLoadImage( "./images/Webcam_grayscale.jpg",0);
	if( orgImg == NULL ) 
        fprintf( stderr, "Cannot load file %s!\n", "./images/TvBig.jpg" );
        
	width=orgImg->width;
	height=orgImg->height;
	pixNum=height*width;
	weightsMatrix=malloc(matrixWidth * pixNum * sizeof(int));
	
    buildWeightsMatrix(weightsMatrix,width,height);
    regions=segmentation(weightsMatrix,pixNum,width,max);
    free(weightsMatrix);
  	constructSegmentedImage(regions,width,height);
  	
  	cvReleaseImage(&orgImg);
  	free(regions);
    return 0;
}
