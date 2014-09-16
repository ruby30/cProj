#include<stdio.h>
#include<stdlib.h>
#include "cv.h"
#include "highgui.h"
#include "mpi.h"
//####################Global Variables######################
IplImage* orgImg ;
//###################Declarations###########################
int* buildWeightsMatrix(int[2],int,int); 
int calColumn(int,int,int);
void calWeight(int,int,int,int[],int,int*);
void constructSegmentedImage(int*,int,int);
int* countRegions(int*,int);
int* countSort(int,int*,int);
int getGrayLevel(int,int);
int* imageToGraph(int,int,int,int);
void initCountsDispls(int*,int*,int,int,int,int,int);
void	initMatrix(int*,int,int);
void initRegions(int*,int,int);
int indexOf(int,int[]);
int* partialSegmentation(int*,int,int*,int,int);
int regionOf(int,int*,int);
void replaceRegion(int*,int,int,int,int);
void resolveXBorderEdges(int*,int*,int,int,int);
int* segmentation(int*,int,int,int,int);
void unionRegions(int,int,int*,int);
//###################Functions##############################
//Convert the image matrix to a weighted graph matrix
int* buildWeightsMatrix(int fromTo[2],int width,int height)
{
    int i,j;
    int mainPix;
    int row=0;
    int max=256;
    int matrixWidth=4;
    int pixNum;
    int fromRow=fromTo[0];
    int toRow;
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if(fromTo[1]==height)
		toRow=fromTo[1];
	else
		toRow=fromTo[1]-1;
	pixNum=width*(toRow-fromRow);

    int* subWeights=malloc(pixNum * matrixWidth * sizeof(int));
	//init weightsMatrix
	initMatrix(subWeights,pixNum*matrixWidth,max);
	//Set subWeights
   for(i = fromRow;i < toRow; i++)
	{
        for(j = 0;j < width; j++)
		{
			mainPix=getGrayLevel(i,j);
			if(i==0)
			{
				if(j==0)
				{
					int skip[]={1,-1};
					calWeight(i,j,row,skip,mainPix,subWeights);
				}
				else
					if(j==width-1)
					{
						int skip[]={0,3,-1};
						calWeight(i,j,row,skip,mainPix,subWeights);
					}
					else
					{
						int skip[]={-1};
						calWeight(i,j,row,skip,mainPix,subWeights);
					}
			}
			else
				if(i==height-1)
				{
					if(j==0)
					{
						int skip[]={1,2,3,-1};
						calWeight(i,j,row,skip,mainPix,subWeights);
					}
					else
						if(j==width-1)
						{
							int skip[]={0,1,2,3,-1};
							calWeight(i,j,row,skip,mainPix,subWeights);
						}
						else
						{
							int skip[]={1,2,3,-1};
							calWeight(i,j,row,skip,mainPix,subWeights);
						}
				}	
				else
					if(j==0)
					{
						int skip[]={1,-1};
						calWeight(i,j,row,skip,mainPix,subWeights);
					}
					else
						if(j==width-1)
						{	
							int skip[]={0,3,-1};
							calWeight(i,j,row,skip,mainPix,subWeights);
						}
					else
					{
						int skip[]={-1};
						calWeight(i,j,row,skip,mainPix,subWeights);
					}
			row++;
		}		
	}
	return subWeights;
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
void calWeight(int i,int j,int row,int skip[],int mainPix,int* weights) 
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
		weights[index]=abs(mainPix-getGrayLevel(i,j));
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
					else
						((uchar*)(segImg->imageData + i * segImg->widthStep)) [ j ]=255;
				}
			}
			char* imgPath=malloc(1000*sizeof(char));
			char imgNum[100];
			char imgFldr[]="./HPC/seg";
			sprintf(imgNum,"%d",index);
			strcpy(imgPath,imgFldr);
			strcat(imgPath,imgNum);
			strcat(imgPath,".jpg");
			cvSaveImage(imgPath,segImg,0);
			cvSet(segImg,cvScalar(256,256,256,256),NULL);
			free(imgPath);
	}
	free(count);
	cvReleaseImage(&segImg);
}
//count the size of each region
int* countRegions(int* regions,int pixs)
{
	int i,index;
	int* count=malloc(pixs*sizeof(int));
		
	initMatrix(count,pixs,0);
	
	for(i = 0 ; i < pixs ; i++)
	{
		if(regions[i]==i)
		index=regionOf(i, regions,0);
		count [ index ] = count [ index ] + 1;
	}
		
	return count;
}
//Counting Sort Algorithm to sort the edges according to their weights
int* countSort(int max,int* matrix,int rows)
{
	int i;
	int cols=4;
	int size=rows*cols;
	int* count=malloc((max+1)*sizeof(int));
	int* sortedIndexes=malloc(size*sizeof(int));
	
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
	if(j==-1)
		return ((int)((uchar*)(orgImg->imageData))[i]);
	else
		return ((int)((uchar*)(orgImg->imageData + i*orgImg->widthStep))[j]);
}
//Convert image to a weighted graph
int* imageToGraph(int width,int height,int rank,int size)
{	
	int i;
	int root=0;
	int subRows=height/size;	
	int* subWeights;
	int imgIndecies[size][2];	
	int fromTo[2];
	imgIndecies[0][0]=0;//subRows-1
	imgIndecies[0][1]=subRows;	
	for(i=1;i<size-1;i++)
	{
		imgIndecies[i][0]=(i*subRows)-1;//subRows
		imgIndecies[i][1]=(i+1)*subRows;
	}
	imgIndecies[size-1][0]=((size-1)*subRows)-1;//subCol+1+height%size
	imgIndecies[size-1][1]=height;
	
	MPI_Scatter(imgIndecies,2, MPI_INT,fromTo,2,MPI_INT, root, MPI_COMM_WORLD);	
	subWeights=buildWeightsMatrix(fromTo,width,height);
	
	return subWeights;
}
//Set Counts & displs values
void initCountsDispls(int* counts,int* displs,int subRows,int width,int matrixWidth,int mod,int size)
{
	MPI_Barrier(MPI_COMM_WORLD);
	int i;
	for(i=0;i<size;i++)
	{
		if(i==0)
			counts[i]=width*(subRows-1)*matrixWidth;
		else
			if(i==size-1)
				counts[i]=width*(subRows+1+mod)*matrixWidth;
			else
				counts[i]=width*subRows*matrixWidth;			
	}
	displs[0]=0;
	for(i=1;i<size;i++)
	{
		displs[i]=displs[i-1]+counts[i-1];
	}
}
//init a matrix with initValue
void	initMatrix(int* matrix,int size,int initValue)
{
	int i;
	for(i=0;i<size;i++)
	{
		matrix[i]=initValue;
	}
}
//init regions for first use
void initRegions(int* regions,int pixs,int startIndex)
{
	int i;
	for(i = 0; i < pixs; i++)
	{
	    regions[i] = startIndex;
	    startIndex++;
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
//Segmentation part
int* partialSegmentation(int* subWeights,int subSize,int* delayed,int startIndex,int width)
{
    int  i,j,ri,rj,r,c,index;
    int* sortedIndexes;
    int max=256;
    int matrixWidth=4;
	int delayedWidth=2;
    int pixs=subSize/matrixWidth;
    startIndex=startIndex/matrixWidth;
    int delayCount=0;
    int* regions=malloc(pixs * sizeof(int));

    initRegions(regions,pixs,startIndex);
    sortedIndexes=countSort(max,subWeights,pixs); 
	
   for(index = 0; index < subSize; index++)
    {
    	r=sortedIndexes[index]/matrixWidth;
    	c=sortedIndexes[index]%matrixWidth;
    	if(subWeights[r*matrixWidth+c]==max) 
    		break;            
	    i = r;
	    j = calColumn(c,r,width); 
	    if(j>=pixs)//cross-border edge
	    {
	    	delayed[delayCount*delayedWidth]=i;
	    	delayed[delayCount*delayedWidth+1]=j;
	    	delayCount++; 	
	    }
	    else
	    {
			ri = regionOf(i, regions,startIndex);
			rj = regionOf(j, regions,startIndex);
			if ((ri != rj) )
			    unionRegions(ri,rj,regions,startIndex); 
	    } 
    }
    
    free(sortedIndexes);
    return regions;
}
//get the set containing v
int regionOf(int v,int *regions,int startIndex)
{
	int index=v;
	v+=startIndex;
    while(regions[index] != v)
    {
		v = regions[index]; 
		index=v-startIndex;     
     }  
	return v;
}
//replace region index in specific interval
void replaceRegion(int* regions,int oldRegion,int newRegion,int from,int to)
{
	int i;
	for(i=from;i<to;i++)
	{
		if(regions[i]==oldRegion)
			regions[i]=newRegion;
	}
}
//calculate & union xBorder edges
void resolveXBorderEdges(int* delayed,int* regions,int pixs,int width,int size)
{
	int index,counter,i,j,blockNum;
	int delayedWidth=2;
	int delayedPerPix=3;
	int delayBoundary=delayedPerPix*width*delayedWidth;
	int regionBoundary=pixs/size;	
	
	for(blockNum=0;blockNum<size-1;blockNum++)
	{
		index=blockNum*delayBoundary;
		while(delayed[index]!=-1)
		{
			counter=blockNum*regionBoundary;
			i=delayed[index]+counter;
	   	 	j = delayed[index+1]+counter;
	   	 	regions[j]=regions[i];
			index+=delayedWidth;
		}
	}
}
//Big segmentation scatter weightsMatrix(pixNums * 4) & gather regions(pixNum * 4)
int* segmentation(int* subWeights,int width,int pixNums,int rank,int size)
{
    int root=0;    
	int matrixWidth=4;
	int delayedPerPix=3;
	int delayedWidth=2;
	int regionsWidth=1;
	int height=pixNums/width;
	int subRows=height/size;
	int localDelayedSize=delayedPerPix*delayedWidth*width; 
	int* localRegions;	
	int* displs=malloc(size*sizeof(int));
	int* counts=malloc(size*sizeof(int));	
	int* localDelayed=malloc(localDelayedSize*sizeof(int));	
	int* regions;	
	int* delayedQ;
	
	if(rank==root)
	{
		regions=malloc(pixNums*sizeof(int));
		delayedQ=malloc(size*localDelayedSize*sizeof(int));
	}
	
	initMatrix(localDelayed,localDelayedSize,-1);
	initCountsDispls(counts,displs,subRows,width,matrixWidth,pixNums%size,size);	
	localRegions=partialSegmentation(subWeights,counts[rank],localDelayed,displs[rank],width);	
	initCountsDispls(counts,displs,subRows,width,regionsWidth,pixNums%size,size);
	
	MPI_Gatherv(localRegions,counts[rank], MPI_INT,regions,counts,displs,MPI_INT,root,MPI_COMM_WORLD);
	MPI_Gather(localDelayed,localDelayedSize, MPI_INT,delayedQ,localDelayedSize,MPI_INT,root,MPI_COMM_WORLD);	
	
	if(rank==root)
		resolveXBorderEdges(delayedQ,regions,pixNums,width,size);
		
	free(localDelayed);
	free(localRegions);	
	free(displs);
	free(counts);
	return regions;
}
//Union ij from different sets
void unionRegions(int ri,int rj,int *regions,int startIndex)
{
    if(rj > ri)
        regions[rj-startIndex] = ri;
    else
        regions[ri-startIndex] = rj;
}
//////////////////////////////////////////////////////////////////////// 
int main(int argc, char **argv)
{
	int width,height;
    int pixNum;
	int rank,size;
	int root=0;
	int* subWeights;
	int* regions;
	
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );

    orgImg = cvLoadImage( "./images/pic1.jpg",0);
	if( orgImg == NULL ) 
        fprintf( stderr, "Cannot load file %s!\n", "./images/man.jpg" );   
    
	width=orgImg->width;
	height=orgImg->height;
	pixNum=height*width;
	
	subWeights=imageToGraph(width,height,rank,size);
    regions=segmentation(subWeights,width,pixNum,rank,size);
  	free(subWeights);
    if(rank==root)
    {
    	constructSegmentedImage(regions,width,height);
    	free(regions);
    }
    
    cvReleaseImage(&orgImg);
    MPI_Finalize();    
    
    return 0;
}
