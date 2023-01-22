//Random data generator 
#include <stdio.h>
#include <math.h>
#include <time.h> // for CPU time
#include <sys/time.h> //for gettimeofday
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define LENGTH 40



int  main(){
  
    long long int counter=0;

    double current_number=0;

    double last_number=0;

    bool sorted=true;


    char name[LENGTH]="1mout.txt",line[LENGTH];

    FILE *fp;

   
    
    printf("Input NAME of sorted data file\n");
    scanf("%s",name); // reading of filename from keyboard 

    fp=fopen(name,"r");
      printf("\n sorted data file <%s> has opened for testing purposes \n ",name);

      while(1){ 
          if(fgets(line, LENGTH,fp)==NULL) break; 
          if(sscanf(line, "%lf",&current_number)==-1) break; 

          if(counter==0){
            last_number=current_number;
            counter++;
          }else{

            if(current_number<last_number){
                //not sorted
                printf("\n Not sorted at line [%i]\n",counter+1);
                sorted=false;
                break;

            }else{
                last_number=current_number;
                counter++;
            }
          }
      }


    if(sorted==true){
        printf("\n data in file <%s> is sorted..\n",name);
    }



    fclose(fp);


    return 0;
}