// Write this notation in 20221115
// This program is for reading a bit file, and showing the data
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <NIDAQmx.h>
#include <stdlib.h> /*system is defined in this header */

//------------ Define some const parameters ---------------
#define if_formalFileName 1
#define TARGET_SAMPLE_NUM 85000000 
// #define TARGET_SAMPLE_NUM 1000*6


#define bitcheck(byte,nbit) ((byte) &   (1<<(nbit)))

void main() {
    int TARGET_BYTE_NUM = ceil((double)TARGET_SAMPLE_NUM/8);
    int TARGET_BIT_RESIDUAL = TARGET_SAMPLE_NUM - floor((double)TARGET_SAMPLE_NUM/8)*8;
    
    static uint8_t array[TARGET_SAMPLE_NUM];
    long sum=0;
    double mean=0;
    int sample_count = 0;
    int high_count = 0;
    int high_count_continous = 0;    
    int trial_marker_count = 0;
    int other_marker_count = 0;
    int end_marker_count = 0;
        
    FILE *f, *fp;
    char fileName[100];
    int tempLen = 0;
    
    if(if_formalFileName == 1) {
        fp=popen("find -wholename \"./MKdata/*PTB*MKdata\" | sort -r | head -n 1", "r");
        fgets(fileName, sizeof(fileName), fp); 
        pclose(fp);
                
        //strcpy(fileName, fileName+2);
        tempLen = strlen(fileName);
        fileName[tempLen-1] = '\0';
        //printf("tempLen=%d\n", tempLen);
        
        
//         strcpy(fileName, "./MKdata/2022-11-15-2P-marker-1-fromMarkerCount0.MKdata");
//         strcpy(fileName, "./MKdata/2022-11-21-PTB-marker-11-fromMarkerCount0.MKdata");

        //printf("File name is %s!\n", fileName+9);                     
        
        f = fopen(fileName, "rb+");
    }
    else {
        f = fopen("array.data", "rb+");
    }
    
    fseek(f, 0, SEEK_END);// file pointer go to end
    int fileSize = ftell(f);
    fseek(f, 0, SEEK_SET);// file pointer go to begin
    
    fread(array, sizeof(uint8_t), fileSize, f);
    
    int ACTUAL_SAMPLE_NUM = fileSize*8;
    int ACTUAL_BYTE_NUM = fileSize;
    printf("TARGET_SAMPLE_NUM = %d\n",TARGET_SAMPLE_NUM);
    printf("ACTUAL_SAMPLE_NUM = %d\n",ACTUAL_SAMPLE_NUM);
    
    for (int i=0;i<ACTUAL_BYTE_NUM;i++) {
        for (int j=0;j<8;j++) {
            //if (i!=TARGET_BYTE_NUM-1 || TARGET_BIT_RESIDUAL==0 || j<TARGET_BIT_RESIDUAL) {
            if (i>=0) {
                if (bitcheck(array[i], j)) {
                    high_count_continous++;
                    high_count++;
                }
                else {
                    if (high_count_continous > 0) {
                        printf("There are %d high samples! (sample_count %d to sample_count %d).\n",\
                                high_count_continous, sample_count-high_count_continous, sample_count-1);
                        
                        //if (high_count_continous == 20) {
                        if (high_count_continous > 60 && high_count_continous <= 100) {
                        	trial_marker_count++;                        	
                        }
                        if (high_count_continous <= 60) {
                        	other_marker_count++;                        	
                        }      
                        if (high_count_continous > 150) {
                        	end_marker_count++;                        	
                        }                          
                        high_count_continous = 0;
                    }
                }
                sample_count++;
            }
            else {
                break;
            }
        }
    }
    printf("high_count = %d.\n", high_count);
    printf("trial_marker_count = %d.\n", trial_marker_count);
    printf("other_marker_count = %d.\n", other_marker_count);
    printf("end_marker_count = %d.\n", end_marker_count);
    printf("total_marker_count = %d.\n", trial_marker_count+other_marker_count+end_marker_count);
    
    
    if(if_formalFileName == 1) {
        printf("File name is %s!\n", fileName+9); 
    }
    
    fclose(f);
}
