/*I called it Binary Merge Sort because thread number should be the power of 2.
exp: 1,2,4,8,16 and so on.
*/

#include <stdio.h>
#include <math.h>
#include <time.h> // for CPU time
#include <sys/time.h> //for gettimeofday
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define max_len 1000000
#define LENGTH 40

//necessary function
int* AllocateSubArrayStartEndIndex(int len,int numberOfThreads);

int* AllocateSubArrayStartEndIndexWithShift(int len,int numberOfThreads);

int* AllocateSubArrayLen(int len,int numberOfThreads);

int* AllocateSubArrayLenWithShift(int len,int numberOfThreads);

int ReturnCountInside(int value,int segment);

int ReturnEndingCount(int value,int segment);

void GenerateDisplacementAndShifCountForScatter(int* displ, int* sendCount, int numberOfThreads,int len, int* startEndIndex);




//Sorting
void Sort(double* array, int sizeOfArray);
bool CheckIfAllSorted(double* endBuffer, double* headerBuffer, int numberOfThreads );

//Injection arrays
void InjectArrayAccordingToSortedHeadElement1(double* b, double* temp_b,int injectQueue,int* sortPathIndexArray, int* displ, int* sendCount, int* startEndIndexForInjections,int* startEndIndex);
void InjectArrayAccordingToSortedHeadElement(double* b, double* temp_b,int injectQueue,int* sortPathIndexArray, int* displ, int* sendCount, int* startEndIndexForInjections);
void MakeSpaceInsidePointerArray(double* temp_b, int arrayLen,int startIndex, int spaceAmount);
void InjectArray(double* b, double* temp_b,int injectQueue, int* displ, int* sendCount, int* startEndIndex);
//Printing function
void PrintDoubleArray(double* array, int sizeOfArray);
void PrintIntArray(int* array, int sizeOfArray);
void CopyDoubleArray(double* save, double* copy, int sizeOfArray );
void CopyIntArray(int* save, int* copy, int sizeOfArray );
void SortElements(double* b, double* temp_b, int currentIndex, int numberofThreads, int* displ, int* sendCount);
double* MergeArrays(double* a ,double* b, int arrayLen1,int arrayLen2);
void MergeArraysOnThread(double* subArray ,int arrayStartIndex1, int arrayStartIndex2,int arrayCount1,int arrayCount2);
void GenerateMergeArrays(double* b, double* temp_b,int injectQueue, int* displ, int* sendCount, int* startEndIndex);
void PushMerge(double* b, double* temp_b, int merge_index, int* startEndIndex, int* displ, int* sendCount);
double* GetSegmentFromGlobalArray(double* globalArray, int startIndex, int sendCount);
int powerOfTwo(int n);
int* GetSegmentOfIntegerArray(int* array, int lenOfArray);
void SetIntArrayZero(int* array, int lenOfArray);

int main(int argc, char *argv[]){
/*  start the main function */

    //  MPI Paralel variables
    int currentThreadIndex,numberOfThreads;
    int hostIndex=0;
    
    //variables for reading double elements from file
    int i=0,lengthOfDoubleArray,ind[max_len+1],j,cur,prev;
    //double b[max_len+1],c[max_len+1],new,cnew,time;
    static double b[max_len];
     
    double* temp_b; //Temporary array for injecting the segments
    
    /*The program will read and sort and do some manipulations according to this bool variable in while loop
    if all is sorted the processing will be false and program will be ended...*/
    bool processing=true;
    int threadingCounter=0; //will count how many while loop will be runned
    int endMainLoopCounter=1; //TESTING
    
    
    /*Communication data*/
    int* startEndIndex;

    int* displ;
    int* sendCount;

    int* displ_merging;
    int* sendCount_merging;



    double* sortingHeadBuffer;
    double* sortingEndBuffer;

    double *headerBuffer;
    double *endBuffer;  

    //time variables
    clock_t cpu_before_reading,cpu_after_reading,cpu_after_calculation_end;
    struct timeval wall_before_reading,wall_after_reading,wall_after_calculation_end;
    double wall_time_for_reading_and_comm_process,wall_time_for_calculation_process,total_wall_time_process;

    double* wall_time_for_calculations;




    //Threading count
    int powerOftwo=1;
    int currentDispatchNumber=0;

    //File name and file opening variables
    char name[LENGTH]="100k.txt",line[LENGTH];
    FILE *fp;
    //time variables
    clock_t cpu0,cpu1,cpu2,cpu3; // clock_t defined in <time.h> and <sys/types.h> as int
    struct timeval time0, time1,time2,time3; // for wall clock in s and us
    double  dtime12,dtime03; // for wall clock in s (real number)

    //TIME CONTROLLER*************************************************************
    clock_t cputime; /* clock_t defined in <time.h> and <sys/types.h> as int */
    double dtime;

    struct timeval start, end;


    cpu0 = clock();    // assign initial CPU time (IN CPU CLOCKS)
    gettimeofday(&time0, NULL); // returns structure with time in s and us (microseconds)
   
    /*Initialize MPI*/
    MPI_Init(&argc,&argv); // initialise MPI
    MPI_Comm_size(MPI_COMM_WORLD,&numberOfThreads); // return total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD,&currentThreadIndex); // return number of this processor, me=0..nproc-1

    /*Allocate necessary data which is going to be used in all threads.
     start end index data according to data len and numberOfThreads*/
    startEndIndex=malloc((numberOfThreads+1)*sizeof(int));

    /*Those variables are necessary for MPI_Scatterv..*/
    displ=malloc((numberOfThreads)*sizeof(int));
    sendCount=malloc((numberOfThreads)*sizeof(int));    

    displ_merging = malloc((numberOfThreads)*sizeof(int));  
    sendCount_merging = malloc((numberOfThreads)*sizeof(int)); 

    



    MPI_Barrier(MPI_COMM_WORLD); //wait for allocation

    /* Read file and generate necessary data
     for sending and receiving data from MPI threads*/
    if(currentThreadIndex == hostIndex){

      //Before reading time
      cpu_before_reading = clock();
      gettimeofday(&wall_before_reading, NULL); 
      wall_time_for_reading_and_comm_process = 0;

      fp=fopen(name,"r");
      printf("\n file <%s> has opened in hostIndex[%i]\n ",name,currentThreadIndex);
      while(1){ //1 serves as true, i.e. condition which is always true
          if(fgets(line, LENGTH,fp)==NULL)break; // finish reading when empty line is read
          if(sscanf(line, "%lf",&b[i])==-1) break; // finish reading after error
          i++;
      }
      lengthOfDoubleArray=i;
      
      //start calculating time
      cputime = clock();    // assign initial CPU time (IN CPU CLOCKS)
      gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)


      //Allocate temporary array for copying it and inserting segments
      temp_b = malloc(sizeof(double)*(lengthOfDoubleArray));


      //Generate necessary data
      startEndIndex = AllocateSubArrayStartEndIndex(lengthOfDoubleArray,numberOfThreads); 
      //generate data with shif depending on len of array

      //startEndIndex = AllocateSubArrayStartEndIndexWithShift(lengthOfDoubleArray,numberOfThreads); //Generate data with shift
      //PrintIntArray(startEndIndex,numberOfThreads+1);

      //Generate data for MPI_Scatterv 
      GenerateDisplacementAndShifCountForScatter(displ,sendCount,numberOfThreads,lengthOfDoubleArray,startEndIndex);

      //PrintIntArray(displ,numberOfThreads);
      //PrintIntArray(sendCount,numberOfThreads);
      
      CopyIntArray(displ_merging,displ,numberOfThreads);
      CopyIntArray(sendCount_merging,sendCount,numberOfThreads);
      


      /*Chek the number of threads*/
      if( powerOfTwo(numberOfThreads)==1 ){
        

        //printf("\n**** [%i] \n", powerOfTwo(numberOfThreads));


        powerOftwo =  (ceil(log2(numberOfThreads)))+1;
      
        currentDispatchNumber = numberOfThreads;
      
        printf("\n Total Number of Threads [%i] is [%i] power of 2, we can safely proceed the process. \n",numberOfThreads,powerOftwo-1);

        /*In this case powerOftwo is process number*/ 
        /* 1 read data on host thread and dispatch it, then sort it. */
        /* Rest will be doing merging till last merge will be done on host index*/
              

      }
      else{
        int status;
         printf("\n Total Number of Threads [%i] is NOT power of 2, terminating process... \n",numberOfThreads);
        exit(status);
      }

    }//end of file reading in host Thread

    MPI_Barrier(MPI_COMM_WORLD); //wait for reading file
  
    MPI_Bcast(&(processing), 1, MPI_C_BOOL, hostIndex, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD); //wait to broeadcast



    MPI_Comm_size(MPI_COMM_WORLD,&numberOfThreads); // return total number of processors

    //Time after reading and communication
        cpu_after_reading = clock();
        gettimeofday(&wall_after_reading, NULL); 
    /*sorting process has started and will be done in this loop*/



    while(processing==true){

        MPI_Bcast(&(threadingCounter), 1, MPI_INT, hostIndex, MPI_COMM_WORLD); //Broadcast threading counter every iteration of loop
        MPI_Bcast(&(currentDispatchNumber), 1, MPI_INT, hostIndex, MPI_COMM_WORLD);


        MPI_Bcast(&(startEndIndex[0]), numberOfThreads+1, MPI_INT, hostIndex, MPI_COMM_WORLD); //it gives memory leak problem on thread# 5

        MPI_Barrier(MPI_COMM_WORLD); //wait to broadcast
      
        //For the MPI_Scatterv
        MPI_Bcast(&(sendCount[0]), numberOfThreads, MPI_INT, hostIndex, MPI_COMM_WORLD); /* we can broadcast smaller data no problem*/
        MPI_Bcast(&(displ[0]), numberOfThreads, MPI_INT, hostIndex, MPI_COMM_WORLD);
        //This is for merging
        MPI_Bcast(&(displ_merging[0]), numberOfThreads, MPI_INT, hostIndex, MPI_COMM_WORLD); /* we can broadcast smaller data no problem*/
        MPI_Bcast(&(sendCount_merging[0]), numberOfThreads, MPI_INT, hostIndex, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD); //wait to broadcast

        int lenOfLocalSubBuffer=sendCount[currentThreadIndex];


        /*Generate local sub buffer for each thread*/
        double *localSubBuffer=malloc(sizeof(double) * (lenOfLocalSubBuffer));
        

        MPI_Scatterv(&b, sendCount ,displ, MPI_DOUBLE, localSubBuffer, lenOfLocalSubBuffer , MPI_DOUBLE, hostIndex, MPI_COMM_WORLD);

   
        
        MPI_Barrier(MPI_COMM_WORLD);//Wait for scatter segment of array

         
    
    wall_time_for_reading_and_comm_process = ((wall_after_reading.tv_sec  - wall_before_reading.tv_sec)+(wall_after_reading.tv_usec - wall_before_reading.tv_usec)/1e6);;

        if(threadingCounter==0){
            //This has to be done once.
            
            Sort(localSubBuffer,lenOfLocalSubBuffer); /*Sort segment of array*/

            //printf("\n [%i] \n",threadingCounter);

            printf("\nSorting on Thread [%i] has finished...  \n",currentThreadIndex);

            currentDispatchNumber=currentDispatchNumber/2;
            

        }else{
            

             if(currentThreadIndex<currentDispatchNumber){
                
              int lenOfLocalMergeSubBuffer = lenOfLocalSubBuffer;  //get rid of this

              double* mergeSubThreadArray = malloc(sizeof(double) * (lenOfLocalMergeSubBuffer));


              int lenOfFirstMergeArray=sendCount_merging[currentThreadIndex*2];
              int firtMergeArrayStartIndex=displ_merging[currentThreadIndex*2]-displ[currentThreadIndex];
              int lenOfSecondMergeArray=sendCount_merging[currentThreadIndex*2+1];
              int secondMergeArrayStartIndex=displ_merging[currentThreadIndex*2+1]-displ[currentThreadIndex];



              //printf("\n first array start [%i] count [%i] \n",firtMergeArrayStartIndex,lenOfFirstMergeArray);
              //printf("\n second array start [%i] count [%i] \n",secondMergeArrayStartIndex,lenOfSecondMergeArray);
              //printf("\n len of sub array [%i] \n",lenOfLocalMergeSubBuffer);

       
             
              MergeArraysOnThread(localSubBuffer,firtMergeArrayStartIndex,secondMergeArrayStartIndex,lenOfFirstMergeArray,lenOfSecondMergeArray);
              printf("\nCurrent SubStep Merging Process on threadNumber [%i] has finished... \n",currentThreadIndex);

              }


            /* Maybe put this one to the loop checker*/
            if(numberOfThreads==1){

            }
            else{
                currentDispatchNumber=currentDispatchNumber/2;

            }

        }

        MPI_Barrier(MPI_COMM_WORLD);//Wait for sorting.




       // MPI_Gatherv(localSubBuffer, sendCount[threadingCounter] , MPI_DOUBLE, &b, sendCount ,displ, MPI_DOUBLE, hostIndex, MPI_COMM_WORLD);
        MPI_Gatherv(localSubBuffer, lenOfLocalSubBuffer , MPI_DOUBLE, &b, sendCount ,displ, MPI_DOUBLE, hostIndex, MPI_COMM_WORLD);

        //MPI_Gatherv(localSubBuffer, testSendCount[threadingCounter] , MPI_DOUBLE, &b, testSendCount ,testDisplacement, MPI_DOUBLE, hostIndex, MPI_COMM_WORLD);

        //change sendCount 

        MPI_Barrier(MPI_COMM_WORLD); 

        free(localSubBuffer);

        if(currentThreadIndex==hostIndex){

       
         
          printf("\n****** CURRENT THREADING PROCESS STEP NUMBER  [%i] ******\n",threadingCounter);
          if(threadingCounter==0){

            printf("\n Sorting Process on step number [%i] has finished ******\n",threadingCounter);

          }
          else{

            printf("\n Merging Process on step number [%i] has finished ******\n",threadingCounter);
          }



          threadingCounter++;
          

          CopyIntArray(displ_merging,displ,numberOfThreads);
          CopyIntArray(sendCount_merging,sendCount,numberOfThreads);
      



           //Change displacement count array
          int* sendCountCorrection =  malloc(sizeof(int) * (currentDispatchNumber));


          int lastNullCount=0;
          for(int i=0;i<currentDispatchNumber;i++){
              sendCountCorrection[i] = sendCount[i*2]+sendCount[i*2+1];
              lastNullCount =  sendCountCorrection[i];
          }

          for(int i=0;i<numberOfThreads;i++){
              if(i<currentDispatchNumber){
                sendCount[i] = sendCountCorrection[i];
              }else{
                sendCount[i] =  lastNullCount;
              }
          }

          

          //Correct the displacement
          int* displCorrection =  malloc(sizeof(int) * (currentDispatchNumber));
          
          int lastNullDisp=0;
          displCorrection[0]=0;
          lastNullDisp =  lengthOfDoubleArray;
          for(int i=1;i<currentDispatchNumber;i++){
              displCorrection[i] = displCorrection[i-1]+sendCount[i-1];
              
          }

          for(int i=0;i<numberOfThreads;i++){
              if(i<currentDispatchNumber){
              displ[i] = displCorrection[i];
              }else{
              displ[i] =  lastNullDisp;
              }
          }



        //If sorted cut the loop and the sorting 
        if(threadingCounter==numberOfThreads || currentDispatchNumber==0){
          

          //Time after reading and communication
          cpu_after_calculation_end = clock();
          gettimeofday(&wall_after_calculation_end, NULL); 
    
          wall_time_for_calculation_process = ((wall_after_calculation_end.tv_sec  - wall_after_reading.tv_sec)+(wall_after_calculation_end.tv_usec - wall_after_reading.tv_usec)/1e6);
          //wall_time_for_calculations[currentThreadIndex]=wall_time_for_calculation_process;

          total_wall_time_process = ((wall_after_calculation_end.tv_sec  - wall_before_reading.tv_sec)+(wall_after_calculation_end.tv_usec - wall_before_reading.tv_usec)/1e6);;


          cputime= clock()-cputime;      // calculate  cpu clock time as difference of times after-before
          gettimeofday(&end, NULL);
          dtime = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
        
          //printf("Elapsed wall time: %f\n", dtime);
          //printf("Elapsed CPU  time: %f\n", (float) cputime/CLOCKS_PER_SEC);



          printf("\n------------------------------------------");
        printf("\nElapsed Wall time for reading: %f\n",  wall_time_for_reading_and_comm_process);
        printf("Elapsed Wall time calculations and comm : %f\n", wall_time_for_calculation_process);
        printf("Elapsed total Wall time: %f", total_wall_time_process);
        printf("\n------------------------------------------\n");

        //New CPU time definition
        clock_t total_reading_and_comm_time = cpu_after_reading-cpu_before_reading;
        clock_t total_time_for_calculations = cpu_after_calculation_end - cpu_after_reading;
        clock_t total_cpu_time  = cpu_after_calculation_end-cpu_before_reading;

        printf("\n------------------------------------------");
        printf("\nElapsed CPU time for reading: %f\n", (float) total_reading_and_comm_time/CLOCKS_PER_SEC);
        printf("Elapsed CPU time calculation and comm: %f\n", (float) total_time_for_calculations/CLOCKS_PER_SEC);
        printf("Elapsed total CPU time: %f", (float) total_cpu_time/CLOCKS_PER_SEC);
        printf("\n------------------------------------------\n");



          fp=fopen("sorted_data.txt","w");   
          for(i=0;i<lengthOfDoubleArray;i++){fprintf(fp,"%lf\n",b[i]);}
          fclose(fp);

          printf("\n file <%s> has closed in hostIndex[%i]\n ",name,currentThreadIndex);
          printf("Sorting process has finished by threading count [%i], Finalizing the MPI..\n",threadingCounter);
          processing=false;
          } 
          //kill the process
        }//host thread process condition





    }//Main looping processs

    //MPI_Finalize();
    return 0;
}/*END OF MAIN FUNCTION------------------------------------------------------------------------------------------------*/



int* AllocateSubArrayStartEndIndex(int len,int numberOfThreads){
  int *subArrayLen;
  
  subArrayLen = malloc((numberOfThreads)*sizeof(int));

  subArrayLen = AllocateSubArrayLen(len,numberOfThreads);

  // PrintIntArray(subArrayLen,numberOfThreads);
    //Allocate index 
    int* startEndIndex;
    startEndIndex=malloc((numberOfThreads+1)*sizeof(int));

    startEndIndex[0]=0;
   
    for(int i=0;i<numberOfThreads+1;i++){
    if(i==0){
      startEndIndex[i+1] = subArrayLen[i];
     
      }
    else{
      startEndIndex[i+1] = startEndIndex[i]+subArrayLen[i];
    
      }
    }
  return  startEndIndex;
}

int* AllocateSubArrayStartEndIndexWithShift(int len,int numberOfThreads){
  int *subArrayLen;
  
  subArrayLen = malloc((numberOfThreads)*sizeof(int));

  subArrayLen = AllocateSubArrayLenWithShift(len,numberOfThreads);


     //PrintIntArray(subArrayLen,numberOfThreads);
    //Allocate index 
    int* startEndIndex;
    startEndIndex=malloc((numberOfThreads+1)*sizeof(int));

    startEndIndex[0]=0;
   
    for(int i=0;i<numberOfThreads+1;i++){
    if(i==0){
      startEndIndex[i+1] = subArrayLen[i];
      }
    else{
      startEndIndex[i+1] = startEndIndex[i]+subArrayLen[i];
    
      }
    }
  return  startEndIndex;
}

int* AllocateSubArrayLen(int len,int numberOfThreads){
  int *subArrayLen;
  
  subArrayLen = malloc((numberOfThreads)*sizeof(int));

  for(int i=0;i<numberOfThreads;i++){
    if(i==0){
      int currentSubArrayLen = ReturnCountInside(len,numberOfThreads)+ReturnEndingCount(len,numberOfThreads);
      //printf("\n subbufferLenght %i \n",currentSubArrayLen);
      subArrayLen[i] =currentSubArrayLen; 
      }
    else{
      int currentSubArrayLen = ReturnCountInside(len,numberOfThreads);
      //printf("\n subbufferLenght %i \n",currentSubArrayLen);
      subArrayLen[i] = currentSubArrayLen;
      }
    }

  return subArrayLen;
}

int* AllocateSubArrayLenWithShift(int len,int numberOfThreads){

  int *subArrayLen;
  int endingCount=0;
  
  subArrayLen = malloc((numberOfThreads)*sizeof(int));

  if(ReturnEndingCount(len,numberOfThreads)==0){
    endingCount=1;
  }else{
    endingCount=ReturnEndingCount(len,numberOfThreads);
  }

  for(int i=0;i<numberOfThreads;i++){
    if(i==0){
      int currentSubArrayLen = ReturnCountInside(len,numberOfThreads)+endingCount;
      //printf("\n subbufferLenght %i \n",currentSubArrayLen);
      subArrayLen[i] =currentSubArrayLen; 
      }
    else if(i==numberOfThreads-1){
      int currentSubArrayLen = ReturnCountInside(len,numberOfThreads);
      //printf("\n subbufferLenght %i \n",currentSubArrayLen);
      subArrayLen[i] = currentSubArrayLen-1;

    }
    else{
      int currentSubArrayLen = ReturnCountInside(len,numberOfThreads);
      //printf("\n subbufferLenght %i \n",currentSubArrayLen);
      subArrayLen[i] = currentSubArrayLen;
      }
    }

  return subArrayLen;
}

int ReturnCountInside(int value,int segment){
    
  int count=0;

  count = (int) value/segment;
  return count;
}

int ReturnEndingCount(int value,int segment){
  int endingCount=0;

  if(value%segment!=0){

    endingCount = value-segment*ReturnCountInside(value,segment);

  }

  return endingCount;
}

void PrintDoubleArray(double* array, int sizeOfArray){
  printf("\n Double array \n");
    for(int i=0;i<sizeOfArray;i++){
    printf(" %f ", array[i]);
  }

}
//Sorting function
void PrintIntArray(int* array, int sizeOfArray){
  printf("\n int array *********");
    for(int i=0;i<sizeOfArray;i++){
    printf("\n [%i] %i \n", i,array[i]);
  }

}

void GenerateDisplacementAndShifCountForScatter(int* displ, int* sendCount, int numberOfThreads,int len, int* startEndIndex){
    
    for(int i=0;i<numberOfThreads;i++){
      displ[i] = startEndIndex[i];
      sendCount[i] = startEndIndex[i+1]-startEndIndex[i];
    }
}

void Sort(double* array, int sizeOfArray){

  for(int i=0;i<sizeOfArray-1;i++){
    
    for(int j=0;j<sizeOfArray-i-1;j++){

        if(array[j]>array[j+1]){
          double temporary = array[j];
          array[j]=array[j+1];
          array[j+1] = temporary;
        }
        //PrintArray(array,sizeOfArray);
    } 
  }
}

//Injecting segments inside total array
void InjectArrayAccordingToSortedHeadElement1(double* b, double* temp_b,int injectQueue,int* sortPathIndexArray, int* displ, int* sendCount, int* startEndIndexForInjections,int* startEndIndex){
    //Sorted array queue
    //int intejctingIndex=sortPathIndexArray[injectQueue];
    int intejctingIndex=injectQueue;
    int startIndex=startEndIndexForInjections[0];


   //int endIndex = startEndIndexForInjections[intejctingIndex];

    int endIndex = startEndIndex[injectQueue];
    printf("\n insertion array start[%i] end[%i] index \n",startIndex,endIndex);

    int injectedArrayStart = displ[intejctingIndex];
    int injectedArrayEnd = (displ[intejctingIndex]+sendCount[intejctingIndex]);

    printf("\n injected  array start[%i] end[%i] index \n",injectedArrayStart,injectedArrayEnd);

    int currentStartCounter1=0; //Counter for starting injecting in array in previous array

    bool stop = false;

    //new Array is temp_b with start 0 and end 

  

     for(int externalCounter=startIndex;externalCounter<endIndex && !stop;externalCounter++){
        //swap places
        for(int internalCounter=displ[intejctingIndex];internalCounter<(displ[intejctingIndex]+sendCount[intejctingIndex]) && !stop;internalCounter++){
          if(temp_b[externalCounter]>b[internalCounter]){
              printf("\n bi [%f] binternal [%f]\n",temp_b[externalCounter],b[internalCounter]);
            currentStartCounter1 = externalCounter;
            stop=true;
            break;
          }else{
             currentStartCounter1 = externalCounter;
          }
        }
      }


      MakeSpaceInsidePointerArray(temp_b,endIndex,currentStartCounter1,sendCount[intejctingIndex]);

      //add with counter
      int memoryCounter=0;
      for(int internalCounter=displ[intejctingIndex];internalCounter<(displ[intejctingIndex]+sendCount[intejctingIndex]);internalCounter++){
          *(&(temp_b[currentStartCounter1])+memoryCounter) = b[internalCounter]; //Iterate over memory 
           memoryCounter++;
      }
}

void InjectArrayAccordingToSortedHeadElement(double* b, double* temp_b,int injectQueue,int* sortPathIndexArray, int* displ, int* sendCount, int* startEndIndexForInjections){
    //Sorted array queue
   // int intejctingIndex=sortPathIndexArray[injectQueue];
    int intejctingIndex=injectQueue;
    int startIndex=startEndIndexForInjections[0];


    int endIndex = startEndIndexForInjections[intejctingIndex+1];

    printf("\n insertion array start[%i] end[%i] index \n",startIndex,endIndex);

    int injectedArrayStart = displ[intejctingIndex];
    int injectedArrayEnd = (displ[intejctingIndex]+sendCount[intejctingIndex]);

    printf("\n injected  array start[%i] end[%i] index \n",injectedArrayStart,injectedArrayEnd);

    int currentStartCounter1=0; //Counter for starting injecting in array in previous array

    bool stop = false;

     for(int externalCounter=startIndex;externalCounter<endIndex && !stop;externalCounter++){
        //swap places
        for(int internalCounter=displ[intejctingIndex];internalCounter<(displ[intejctingIndex]+sendCount[intejctingIndex]) && !stop;internalCounter++){
          if(b[externalCounter]>b[internalCounter]){
              printf("\n bi [%f] binternal [%f]\n",b[externalCounter],b[internalCounter]);
            currentStartCounter1 = externalCounter;
            stop=true;
            break;
          }
        }
      } 

      if(stop==false){
        currentStartCounter1=displ[0]+sendCount[0]-1;
      }
      //Refil the new array according to memory indexing
      int simpleMemoryCounter=0;
      
      for(int i=displ[0];i<(currentStartCounter1);i++){
          *(&(temp_b[0])+simpleMemoryCounter) = b[i];
          simpleMemoryCounter++;
      }
      
      for(int internalCounter=displ[intejctingIndex];internalCounter<(displ[intejctingIndex]+sendCount[intejctingIndex]);internalCounter++){
          *(&(temp_b[0])+simpleMemoryCounter) = b[internalCounter]; //Iterate over memory 
           simpleMemoryCounter++;
      }
      
     

      for(int i=currentStartCounter1;i<sendCount[0];i++){
      
         *(&(temp_b[0])+simpleMemoryCounter) = b[i];
         simpleMemoryCounter++;
      }


}

void InjectArray(double* b, double* temp_b, int injectQueue, int* displ, int* sendCount, int* startEndIndex){
    //Sorted array queue
   // int intejctingIndex=sortPathIndexArray[injectQueue];
    int intejctingIndex=injectQueue;

    int startIndex=startEndIndex[0];
    int endIndex = startEndIndex[intejctingIndex];

    
    
   // printf("\n insertion array start[%i] end[%i] index \n",startIndex,endIndex);

    int injectedArrayStart = displ[intejctingIndex];
    int injectedArrayEnd = (displ[intejctingIndex]+sendCount[intejctingIndex]);

    //int injectedArrayStart = startEndIndex[intejctingIndex];
    //int injectedArrayEnd = startEndIndex[intejctingIndex+1];

    //printf("\n injected  array start[%i] end[%i] index \n",injectedArrayStart,injectedArrayEnd);
    
    int currentInjectionStartCounter=0; //Counter for starting injecting in array in previous array

    bool stop = false;

    for(int externalCounter=startIndex;externalCounter<endIndex && !stop;externalCounter++){
        //swap places
        for(int internalCounter=injectedArrayStart;internalCounter<injectedArrayEnd && !stop;internalCounter++){
          if(b[externalCounter]>b[internalCounter]){
            //printf("\n bi [%f] binternal [%f]\n",b[externalCounter],b[internalCounter]);
            currentInjectionStartCounter = externalCounter;
            stop=true;
            break;
          }
        }
      } 

    if(stop==false){
    currentInjectionStartCounter=endIndex;
    }
    // CopyArraySegmment(); 
    if(injectQueue==1){
      CopyDoubleArray(temp_b,b,endIndex);
    }
   
    MakeSpaceInsidePointerArray(temp_b,endIndex,currentInjectionStartCounter,sendCount[intejctingIndex]);
    
    //add with counter
    int memoryCounter=0;
    for(int internalCounter=displ[intejctingIndex];internalCounter<(displ[intejctingIndex]+sendCount[intejctingIndex]);internalCounter++){
        *(&(temp_b[currentInjectionStartCounter])+memoryCounter) = b[internalCounter]; //Iterate over memory 
        memoryCounter++;
    }


}

void MakeSpaceInsidePointerArray(double* temp_b, int arrayLen,int startIndex, int spaceAmount){

  //startIndex+spaceAmount<arrayLen

  for(int i=arrayLen;i>=startIndex;i--){

    double tempValue = temp_b[i];
    *(&(temp_b[0])+i+spaceAmount) = tempValue;

    if(i<=startIndex+spaceAmount && i >=startIndex){ /*We can play whith condition*/
      temp_b[i]=0;
    }

  }

}

void CopyDoubleArray(double* save, double* copy, int sizeOfArray ){
     for(int j=0;j<sizeOfArray;j++){
      save[j] = copy[j];
    } 
}

void CopyIntArray(int* save, int* copy, int sizeOfArray ){
     for(int j=0;j<sizeOfArray;j++){
      save[j] = copy[j];
    } 
}



bool CheckIfAllSorted(double* endBuffer, double* headerBuffer, int numberOfThreads ){
  
  bool sorted=false;
  int counter=0;

  //PrintDoubleArray(endBuffer,numberOfThreads);
  //PrintDoubleArray(headerBuffer,numberOfThreads);

  for(int i=0;i<numberOfThreads-1;i++){
    if(endBuffer[i]<headerBuffer[i+1] ){
      counter++;
    //  printf("\n sorted [%i] endbuffer [%f] headerBuffer [%f]",i,endBuffer[i],headerBuffer[i+1]);
    }
  }

  if(counter==numberOfThreads-1){
    sorted=true;
  }

  //printf(sorted ? "true" : "false");
  return sorted;
}



void SortElements(double* b, double* temp_b, int currentIndex, int numberofThreads, int* displ, int* sendCount){

    double* currentElementsToSort = malloc(sizeof(double) * (numberofThreads));
    
    //avoid if first index has gap

    for(int i=0;i<numberofThreads;i++){
        
        currentElementsToSort[i] = b[displ[i]+currentIndex];
    }

    //PrintDoubleArray(currentElementsToSort,numberofThreads);
    
    Sort(currentElementsToSort,numberofThreads);

    //PrintDoubleArray(currentElementsToSort,numberofThreads);

    for(int i=0;i<numberofThreads;i++){
        
        temp_b[currentIndex*numberofThreads+i] = currentElementsToSort[i];
    }

   

}

double* MergeArrays(double* arr1 ,double* arr2, int arrayLen1,int arrayLen2){

   // double* mergedArrays = malloc(sizeof(double) * (arrayLen1+arrayLen2));
    

    //declare an empty array (array III) of size n1+n2
    //int size = sizeof(arr1)/sizeof(arr1[0]) + sizeof(arr2)/sizeof(arr2[0]);
    int size = arrayLen1+arrayLen2;

    //printf(" Mergin started with  size [%i]",size);

    double* arr3 = malloc(sizeof(double) * (size));

    //Sorted Merge

    int i=0, j=0, k=0;
    while(i<arrayLen1 && j<arrayLen2){
      if(arr1[i]< arr2[j]){
            arr3[k] = arr1[i];
            k++;
            i++;
        }else{
          arr3[k]=arr2[j];
          k++;
          j++;
        }
    }

    while(i<arrayLen1){
      arr3[k] = arr1[i];
      k++;
      i++;
    }

    while(j<arrayLen2){
      arr3[k] = arr2[j];
      k++;
      j++;
    }
    
    

    //PrintDoubleArray(arr3,size);
    return arr3;
}


void PushMerge(double* b, double* temp_b, int merge_index, int* startEndIndex, int* displ, int* sendCount){

  int startPush = startEndIndex[merge_index];
  int endPush =  startEndIndex[merge_index+1];


  double* mergeArray1 = malloc(sizeof(double) * (sendCount[merge_index]));
  mergeArray1 = GetSegmentFromGlobalArray(b,displ[merge_index],sendCount[merge_index]);
  

  double* mergeArray2 = malloc(sizeof(double) * (sendCount[merge_index+1]));
  mergeArray2 = GetSegmentFromGlobalArray(b,displ[merge_index+1],sendCount[merge_index+1]);


  printf("\n Merge %i %i \n",merge_index,merge_index+1);

  printf("\n Array Lens %i %i \n",sendCount[merge_index],sendCount[merge_index+1]);
  PrintDoubleArray(mergeArray1,sendCount[merge_index]);
  PrintDoubleArray(mergeArray2,sendCount[merge_index+1]);

  MergeArrays(mergeArray1 ,mergeArray2, sendCount[merge_index] ,sendCount[merge_index+1]);
  
}

double* GetSegmentFromGlobalArray(double* globalArray, int startIndex, int sendCount){
  
  double* segmnetArray = malloc(sizeof(double) * (sendCount));

  int counter=0;
  for(int i=startIndex;i<startIndex+sendCount;i++){
    segmnetArray[counter] = globalArray[i]; 
    counter++;
  }

  return segmnetArray;
}

int powerOfTwo(int n){
        //if(n==0) { return 0; }
        //while(n != 1)
        //{   
        //  printf("\n** [%i] \n",n);
        //    n = n/2;  

        //    if(n%2 != 0 && n != 1){ return 0; }
       // }

       // return 1; 
    if (n == 0)
      return 0;
    while (n != 1) {
      if (n % 2 != 0)
        return 0;
      n = n / 2;
      }
    return 1;
  }

int* GetSegmentOfIntegerArray(int* array, int lenOfArray){
    int* returnArray = malloc(sizeof(int) * (lenOfArray));
    for(int i=0;i<lenOfArray;i++){
      returnArray[i] = array[i];
    }

    return returnArray;
}

void SetIntArrayZero(int* array, int lenOfArray){

  for(int i=0;i<lenOfArray;i++){
    array[i]=0;
  }
}


void MergeArraysOnThread(double* subArray , int arrayStartIndex1, int arrayStartIndex2,int arrayCount1,int arrayCount2){

    //allocate first and to arrays
    double* firstArray = malloc(sizeof(double)*(arrayCount1));
    double* secondArray = malloc(sizeof(double)*(arrayCount2));

    double* testArray = malloc(sizeof(double)*(arrayCount2+arrayCount1));

    //we have to add loop for this 
    //we can optimize it also

    for(int i=0;i<arrayCount1;i++){
      firstArray[i]=subArray[arrayStartIndex1+i];
    }

     for(int i=0;i<arrayCount2;i++){
      secondArray[i]=subArray[arrayStartIndex2+i];
    }
    


    testArray = MergeArrays(firstArray,secondArray,arrayCount1,arrayCount2);

    for(int i=0;i<arrayCount1+arrayCount2;i++){
      
                subArray[i]=testArray[i];
                //localSubBuffer[i]= localSubBuffer[i]+1;
              } //test
    
    // This one should work wihout this but it is not working workout on this

    //ree(firstArray);
    //free(secondArray);

    free(testArray);
    free(firstArray);
    free(secondArray);

   // printf("\n First element first array [%f]\n",firstArray[0]);
   // printf("\n last element first array [%f]\n",firstArray[arrayCount1-1]);
   // printf("\n First element second array [%f]\n",secondArray[0]);
   // printf("\n last element second array [%f]\n",secondArray[arrayCount2-1]);

    // deallocate

}




