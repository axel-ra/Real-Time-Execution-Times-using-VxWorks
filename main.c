/* 
 * 
 * 
 * Disclosures/References:
 * [1] Declaration and memory allocation for double pointer. 
 * 	   Credit: https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/#
 * 	   
 * [2] https://www.geeksforgeeks.org/c-rand-function/
 * 
 * [3] I got the C-style code (Pseudocode?) for the IJK, JKI, KIJ matrix multiplication
 * 	   algorithms from the following location:
 * 	   https://courses.engr.illinois.edu/cs232/sp2009/lectures/X18.pdf
 * 	   I just had to take care of the matrixes.
 * 
 * [4] I got the algorithm for the Strassen's matrix multiplication from the following video
 * 	   https://www.youtube.com/watch?v=0oJyNmEbS4w
 * 	   I put the algorithm into C code and basically adapted the algorithm from the video
 * 	   to the C code I used in my program.
 * 	   
 * [5] A classmate mentioned the sysTimestamp() function to me after I asked him something.
 * 	   I ended up using the sysTimestamp() function instead of tickGet().
 * 
 * See report document for other references.
 * 
 */

#include <taskLib.h>
#include <stdio.h>
#include <kernelLib.h>

#include "helper.h" // another file I made with helper functions

// Idea of how to declare the locations in the array is from sample file for this project
unsigned int exec_times[] = {0, 0, 0, 0}; 

// VxWorks libraries
#include <taskHookLib.h>

// C libraries
#include <unistd.h>
#include <stdint.h>
#include <sysLib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

struct data {
	char* action;
	unsigned int timeStamp;
};

// Four arrays to store information from the hook function
// for each of the tasks (each of the algorithms).
struct data timedata0[10000];
struct data timedata1[10000];
struct data timedata2[10000];
struct data timedata3[10000];

// Indexes used to access the arrays created above
int timedata0_index = 0;
int timedata1_index = 0;
int timedata2_index = 0;
int timedata3_index = 0;

TASK_ID id0, id1, id2, id3, guest_id;

// declaring the pointers for matrix a and b (the matrices to be multiplied)
double** _a;
double** _b;

void switchHook (WIND_TCB *pOldTcb, WIND_TCB *pNewTcb){
	
	WIND_TCB* tTcb0 = taskTcb(id0);
	WIND_TCB* tTcb1 = taskTcb(id1);
	WIND_TCB* tTcb2 = taskTcb(id2);
	WIND_TCB* tTcb3 = taskTcb(id3);
	
	// if the old task and the new task are different, proceed
	
	if (pOldTcb != pNewTcb){
		
		
		if (pOldTcb == tTcb0 || pOldTcb == tTcb1 || pOldTcb == tTcb2 || pOldTcb == tTcb3 ||
			pNewTcb == tTcb0 || pNewTcb == tTcb1 || pNewTcb == tTcb2 || pNewTcb == tTcb3) {
				
				struct data tmpOld, tmpNew;
				
				if ((pOldTcb == tTcb0 || pOldTcb == tTcb1 || pOldTcb == tTcb2 || pOldTcb == tTcb3) && (pNewTcb == tTcb0 || pNewTcb == tTcb1 || pNewTcb == tTcb2 || pNewTcb == tTcb3))
				{
					// Both old task and new task are from Task0, Task1, Task2, and Task3
					
					// Summarized explanation:
					// The code checks which task was the old task, gets the timestamp using sysTimestamp(),
					// and stores the information in the array of data (user-defined struct data)
					// timedata0, timedata1, timedata2, or timedata3 -depending on the which is the matched task.
					
					if (pOldTcb == tTcb0){
						tmpOld.action = "old task: task0 -IJK task has been preempted-";
						tmpOld.timeStamp = sysTimestamp();
						if (timedata0_index < 10000){
							timedata0[timedata0_index] = tmpOld;
							timedata0_index++;
						}
						
					}
					else if (pOldTcb == tTcb1){
						tmpOld.action = "old task: task1 -JKI task has been preempted-";
						tmpOld.timeStamp = sysTimestamp();
						if (timedata1_index < 10000){
							timedata1[timedata1_index] = tmpOld;
							timedata1_index++;
						}
						
					}
					else if (pOldTcb == tTcb2){
						tmpOld.action = "old task: task2 -KIJ task has been preempted-";
						tmpOld.timeStamp = sysTimestamp();
						if (timedata2_index < 10000){
							timedata2[timedata2_index] = tmpOld;
							timedata2_index++;
						}
						
					}
					else if (pOldTcb == tTcb3){
						tmpOld.action = "old task: task3 -S. Matrix Multiplication task has been preempted-";
						tmpOld.timeStamp = sysTimestamp();
						if (timedata3_index < 10000){
							timedata3[timedata3_index] = tmpOld;
							timedata3_index++;
						}
						
					}
					
					// Check which task is the new task:
					
					// Summarized explanation:
					// The process is the same than in the previous chain of
					// if statements, but now the code is checking which out of 
					// Task0, Task1, Task2, and Task3 is the new task.
					
					if (pNewTcb == tTcb0){
						tmpNew.action = "new task: task0 -IJK task scheduled to run-";
						tmpNew.timeStamp = sysTimestamp();
						if (timedata0_index < 10000){
							timedata0[timedata0_index] = tmpNew;
							timedata0_index++;
						}
						
					}
					else if (pNewTcb == tTcb1){
						tmpNew.action = "new task: task1 -JKI task scheduled to run-";
						tmpNew.timeStamp = sysTimestamp();
						if (timedata1_index < 10000){
							timedata1[timedata1_index] = tmpNew;
							timedata1_index++;
						}
						
					}
					else if (pNewTcb == tTcb2){
						tmpNew.action = "new task: task2 -KIJ task scheduled to run-";
						tmpNew.timeStamp = sysTimestamp();
						if (timedata2_index < 10000){
							timedata2[timedata2_index] = tmpNew;
							timedata2_index++;
						}
						
					}
					else if (pNewTcb == tTcb3){
						tmpNew.action = "new task: task3 -S. Matrix Multiplication task scheduled to run-";
						tmpNew.timeStamp = sysTimestamp();
						if (timedata3_index < 10000){
							timedata3[timedata3_index] = tmpNew;
							timedata3_index++;
						}
						
					}
					
				}
				else if (pOldTcb == tTcb0 || pOldTcb == tTcb1 || pOldTcb == tTcb2 || pOldTcb == tTcb3){
					
					// Checking which task out of Task0, Task1, Task2, and Task3 was the old task,
					// and getting and writing the information accordingly
					
					if (pOldTcb == tTcb0){
						tmpOld.action = "old task: task0 -IJK task has been preempted-, new task is non-matrixmultiplication task";
						tmpOld.timeStamp = sysTimestamp();
						if (timedata0_index < 10000){
							timedata0[timedata0_index] = tmpOld;
							timedata0_index++;
						}
						
					}
					else if (pOldTcb == tTcb1){
						tmpOld.action = "old task: task1 -JKI task has been preempted-, new task is non-matrixmultiplication task";
						tmpOld.timeStamp = sysTimestamp();
						if (timedata1_index < 10000){
							timedata1[timedata1_index] = tmpOld;
							timedata1_index++;
						}
						
					}
					else if (pOldTcb == tTcb2){
						tmpOld.action = "old task: task2 -KIJ task has been preempted-, new task is non-matrixmultiplication task";
						tmpOld.timeStamp = sysTimestamp();
						if (timedata2_index < 10000){
							timedata2[timedata2_index] = tmpOld;
							timedata2_index++;
						}
						
					}
					else if (pOldTcb == tTcb3){
						tmpOld.action = "old task: task3 -S. Matrix Multiplication task has been preempted-, new task is non-matrixmultiplication task";
						tmpOld.timeStamp = sysTimestamp();
						if (timedata3_index < 10000){
							timedata3[timedata3_index] = tmpOld;
							timedata3_index++;
						}
						
					}
					
				}
				else if (pNewTcb == tTcb0 || pNewTcb == tTcb1 || pNewTcb == tTcb2 || pNewTcb == tTcb3){

					// Checking which task out of Task0, Task1, Task2, and Task3 is the new task,
					// and getting and writing the information accordingly
					
					if (pNewTcb == tTcb0){
						tmpNew.action = "new task: task0 -IJK task scheduled to run-,  old task is non-matrixmultiplication task";
						tmpNew.timeStamp = sysTimestamp();
						if (timedata0_index < 10000){
							timedata0[timedata0_index] = tmpNew;
							timedata0_index++;
						}
						
					}
					else if (pNewTcb == tTcb1){
						tmpNew.action = "new task: task1 -JKI task scheduled to run-, old task is non-matrixmultiplication task";
						tmpNew.timeStamp = sysTimestamp();
						if (timedata1_index < 10000){
							timedata1[timedata1_index] = tmpNew;
							timedata1_index++;
						}
						
					}
					else if (pNewTcb == tTcb2){
						tmpNew.action = "new task: task2 -KIJ task scheduled to run-, old task is non-matrixmultiplication task";
						tmpNew.timeStamp = sysTimestamp();
						if (timedata2_index < 10000){
							timedata2[timedata2_index] = tmpNew;
							timedata2_index++;
						}
						
					}
					else if (pNewTcb == tTcb3){
						tmpNew.action = "new task: task3 -S. Matrix Multiplication task scheduled to run-, old task is non-matrixmultiplication task";
						tmpNew.timeStamp = sysTimestamp();
						if (timedata3_index < 10000){
							timedata3[timedata3_index] = tmpNew;
							timedata3_index++;
						}
						
					}
					
					
				}
			
			}
		
	}
	
	
	return;
}

void initializeArrays(void){
	
	// initialize all locations of the timeStamp arrays (from struct "data") 
	// in the timedata arrays to 0
	for (int i = 0; i < 10000; i++){
		timedata0[i].timeStamp = 0;
		timedata1[i].timeStamp = 0;
		timedata2[i].timeStamp = 0;
		timedata3[i].timeStamp = 0;
	}
	
}

bool isQEmpty(int q){
	// return true if the array with the timestamps for a given
	// task q is "empty", return false if it is not. 
	// ("empty" means all locations in the array timeStamp 
	// have the value 0 in them).
	
	if (q == 0){ 		// q = 0 for task0
		
		for (int i = 0; i < 10000; i++){
			if (timedata0[i].timeStamp != 0)
				return false;
		}
		return true;
		
	}
	else if (q == 1){ 	// q = 1 for task1
		
		for (int i = 0; i < 10000; i++){
			if (timedata1[i].timeStamp != 0)
				return false;
		}
		return true;
		
	}
	else if (q == 2){ 	// q = 2 for task2
		
		for (int i = 0; i < 10000; i++){
			if (timedata2[i].timeStamp != 0)
				return false;
		}
		return true;
		
	}
	else if (q == 3){ 	// q = 3 for task3
		
		for (int i = 0; i < 10000; i++){
			if (timedata3[i].timeStamp != 0)
				return false;
		}
		return true;
		
	}
	else {
		
		return (ERROR); // return error if invalid task number is passed to isQEmpty
	
	}
}

void guest( void ){
	
	// while task0, task1, task2, and task3 have not finished, guest task keeps rescheduling itself 
	while (taskIdVerify(id0) == OK || taskIdVerify(id1) == OK || taskIdVerify(id2) == OK || taskIdVerify(id3) == OK){
		taskDelay(10);
	}
	
	// Task0, Task1, Task2, and Task3 (All the tasks performing matrix multiplication)
	// have ended execution, proceed to calculate the execution times:
	
    unsigned int l = 0;
    unsigned int r = 0;
    unsigned int diff = 0;
	
	// if array with timestamps for task0 is not empty
	if (!isQEmpty(0)){
		
	    for (int i = 0; i < 10000; i++){
	    	if (timedata0[i].timeStamp != 0){
	    		if (strstr(timedata0[i].action, "new task: task0") != NULL)
	    			l = timedata0[i].timeStamp;
	    		else if (strstr(timedata0[i].action, "old task: task0") != NULL)
	    			r = timedata0[i].timeStamp;
	    	
	    		if (l != 0 && r != 0){
	    			//printf("subtracting %d - %d\n", (unsigned int) r, (unsigned int) l);
	    			diff = r - l;
	    			//printf("adding difference %d\n", (unsigned int) diff);
	    			exec_times[0] += diff;
	    			l = 0;
	    			r = 0;
	    		}
	    	}
	    }
	}
	
	// if array with timestamps for task1 is not empty
	if (!isQEmpty(1)){
	    
	    for (int i = 0; i < 10000; i++){
	    	if (timedata1[i].timeStamp != 0){
	    		if (strstr(timedata1[i].action, "new task: task1") != NULL)
	    			l = timedata1[i].timeStamp;
	    		else if (strstr(timedata1[i].action, "old task: task1") != NULL)
	    			r = timedata1[i].timeStamp;
	    	
	    		if (l != 0 && r != 0){
	    			//printf("subtracting %d - %d\n", (unsigned int) r, (unsigned int) l);
	    			diff = r - l;
	    			//printf("adding difference %d\n", (unsigned int) diff);
	    			exec_times[1] += diff;
	    			l = 0;
	    			r = 0;
	    		}
	    	}
	    }
	}
	
	// if array with timestamps for task2 is not empty
	if (!isQEmpty(2)){
	    
	    for (int i = 0; i < 10000; i++){
	    	if (timedata2[i].timeStamp != 0){
	    		if (strstr(timedata2[i].action, "new task: task2") != NULL)
	    			l = timedata2[i].timeStamp;
	    		else if (strstr(timedata2[i].action, "old task: task2") != NULL)
	    			r = timedata2[i].timeStamp;
	    	
	    		if (l != 0 && r != 0){
	    			//printf("subtracting %d - %d\n", (unsigned int) r, (unsigned int) l);
	    			diff = r - l;
	    			//printf("adding difference %d\n", (unsigned int) diff);
	    			exec_times[2] += diff;
	    			l = 0;
	    			r = 0;
	    		}
	    	}
	    }
	}
	
	// if array with timestamps for task3 is not empty
	if (!isQEmpty(3)){
		
	    for (int i = 0; i < 10000; i++){
	    	if (timedata3[i].timeStamp != 0){
	    		if (strstr(timedata3[i].action, "new task: task3") != NULL)
	    			l = timedata3[i].timeStamp;
	    		else if (strstr(timedata3[i].action, "old task: task3") != NULL)
	    			r = timedata3[i].timeStamp;
	    	
	    		if (l != 0 && r != 0){
	    			//printf("subtracting %d - %d\n", (unsigned int) r, (unsigned int) l);
	    			diff = r - l;
	    			//printf("adding difference %d\n", (unsigned int) diff);
	    			exec_times[3] += diff;
	    			l = 0;
	    			r = 0;
	    		}
	    	}
	    }
	}

	// Printing the results:
	
	unsigned int systimefreq = sysTimestampFreq();
	
	// calculations for task0
	double et0 = ((double) exec_times[0]) / ((double) systimefreq);
	double ms0 = et0*( (double) 1000);
	
	// calculations for task1
	double et1 = ((double) exec_times[1]) / ((double) systimefreq);
	double ms1 = et1*( (double) 1000);
	
	// calculations for task2
	double et2 = ((double) exec_times[2]) / ((double) systimefreq);
	double ms2 = et2*( (double) 1000);
	
	// calculations for task3
	double et3 = ((double) exec_times[3]) / ((double) systimefreq);
	double ms3 = et3*( (double) 1000);
	
	//printf("Sysclockrate is %d\n", (unsigned int) systimefreq);
		
	printf("Execution times are:\n");
		
	// task 0
	printf("IJK Algorithm:\n");
	printf("%.25f s\n", (double) et0);
	printf("%.25f ms\n", ms0);
		
	// task 1
	printf("JKI Algorithm:\n");
	printf("%.25f s\n", (double) et1);
	printf("%.25f ms\n", ms1);
		
	// task 2
	printf("KIJ Algorithm:\n");
	printf("%.25f s\n", (double) et2);
	printf("%.25f ms\n", ms2);
		
	// task 3
	printf("Strassen's Algorithm:\n");
	printf("%.25f s\n", (double) et3);
	printf("%.25f ms\n", ms3);
		
	printf("End of program.");
	
	
}

void IJK(int tnum, int _n){
    double** c = pointerCGenerator(_n);

    double sum;

    for (int i = 0; i < _n; i++){
        for (int j = 0; j < _n; j++){
            sum = 0.0;
            for (int k = 0; k < _n; k++){
                sum += _a[i][k] * _b[k][j];
            }
            c[i][j] = sum;
        }
    }
    


}

void JKI(int tnum, int _n){
    double** c = pointerCGenerator(_n);

    double res;

    for (int j = 0; j < _n; j++){
        for (int k = 0; k < _n; k++){
            res = _b[k][j];
            for (int i = 0; i < _n; i++)
                c[i][j] += _a[i][k] * res;
        }
    }
    

    
}

void KIJ(int tnum, int _n){
    double** c = pointerCGenerator(_n);

    double res;

    for (int k = 0; k < _n; k++){
        for (int i = 0; i < _n; i++){
            res = _a[i][k];
            for (int j = 0; j < _n; j++)
                c[i][j] += res * _b[k][j];
        }
    }
    


}

// functions for Strassen matrix multiplication
double** MAddition(double** m1, double** m2, int n){
    double** res = (double**) malloc(n*sizeof(double*)); 		// Credit: [1]

    for (int i = 0; i < n;i++) 							 		// Credit: [1]
        res[i] = (double*) malloc(n*sizeof(double));			// Credit: [1]
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            res[i][j] = m1[i][j] + m2[i][j];
    }

    return res;
}

double** SMM(double** a, double** b, int n){
    if (n <= 2){
        double** res = (double**) malloc(n*sizeof(double*));	// Credit: [1]

        for (int i = 0; i < n;i++)								// Credit: [1]
            res[i] = (double*) malloc(n*sizeof(double));		// Credit: [1]

        res[0][0] = (a[0][0] * b[0][0]) + (a[0][1] * b[1][0]);
        res[0][1] = (a[0][0] * b[0][1]) + (a[0][1] * b[1][1]);
        res[1][0] = (a[1][0] * b[0][0]) + (a[1][1] * b[1][0]);
        res[1][1] = (a[1][0] * b[0][1]) + (a[1][1] * b[1][1]);

        return res;
    }
    else {
        int mid = n/2;

        double** A11;
        double** A12; 
        double** A21;
        double** A22;
        double** B11;
        double** B12;
        double** B21;
        double** B22;

        double** C11;
        double** C12; 
        double** C21;
        double** C22;

        A11 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            A11[i] = (double*) malloc(mid*sizeof(double));

        A12 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            A12[i] = (double*) malloc(mid*sizeof(double));

        A21 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            A21[i] = (double*) malloc(mid*sizeof(double));

        A22 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            A22[i] = (double*) malloc(mid*sizeof(double));

        B11 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            B11[i] = (double*) malloc(mid*sizeof(double));

        B12 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            B12[i] = (double*) malloc(mid*sizeof(double));

        B21 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            B21[i] = (double*) malloc(mid*sizeof(double));

        B22 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            B22[i] = (double*) malloc(mid*sizeof(double));
        
        // initialize A11, A12, A21, A22, B11, B12, B21, B22
        for (int i = 0; i<mid; i++){
            for (int j = 0; j<mid;j++){
                A11[i][j] = a[i][j];
                B11[i][j] = b[i][j];
            }
        }

        for (int i = 0; i<mid; i++){
            for (int j = 0; j<mid;j++){
                A12[i][j] = a[i][j + mid];
                B12[i][j] = b[i][j + mid];
            }
        }

        for (int i = 0; i<mid; i++){
            for (int j = 0; j<mid;j++){
                A21[i][j] = a[i + mid][j];
                B21[i][j] = b[i + mid][j];
            }
        }

        for (int i = 0; i<mid; i++){
            for (int j = 0; j<mid;j++){
                A22[i][j] = a[i + mid][j + mid];
                B22[i][j] = b[i + mid][j + mid];
            }
        }



        C11 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            C11[i] = (double*) malloc(mid*sizeof(double));

        C12 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            C12[i] = (double*) malloc(mid*sizeof(double));

        C21 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            C21[i] = (double*) malloc(mid*sizeof(double));

        C22 = (double**) malloc(mid*sizeof(double*));
        for (int i = 0; i < mid;i++)
            C22[i] = (double*) malloc(mid*sizeof(double));

        C11 = MAddition(SMM(A11, B11, mid),SMM(A12, B21, mid), mid);
        C12 = MAddition(SMM(A11, B12, mid),SMM(A12, B22, mid), mid);
        C21 = MAddition(SMM(A21, B11, mid),SMM(A22, B21, mid), mid);
        C22 = MAddition(SMM(A21, B12, mid),SMM(A22, B22, mid), mid);
      

        double** finalAnswer = (double**) malloc(n*sizeof(double*));
        
        for (int i = 0; i < n;i++)
            finalAnswer[i] = (double*) malloc(n*sizeof(double));

        for (int i = 0; i < mid; i++){
            for (int j = 0; j < mid; j++)
                finalAnswer[i][j] = C11[i][j];
        }

        for (int i = 0; i < mid; i++){
            for (int j = mid; j < n; j++)
                finalAnswer[i][j] = C12[i][j - mid];
        }

        for (int i = mid; i < n; i++){
            for (int j = 0; j < mid; j++)
                finalAnswer[i][j] = C21[i - mid][j];
        }

        for (int i = mid; i < n; i++){
            for (int j = mid; j < n; j++)
                finalAnswer[i][j] = C22[i - mid][j - mid];
        }

        return finalAnswer;
        
    }
}


void Strassen(int tnum, int _n){

    double** result = SMM(_a, _b, _n);
   
}

void CreateTasks(void)
{	
		// ******
		// Code below from sample file shown by TA Thomas in the video
	
        kernelTimeSlice(1);

        cpuset_t affinity;
        CPUSET_ZERO(affinity);
        CPUSET_SET(affinity, 0);
        taskCpuAffinitySet(taskIdSelf(), affinity);
        
        // ******

        int r = rand() % (7 - 3 + 1) + 3; 				// Credit: [2]
        
        // For simplicity, I am only generating matrices which size is a power of 2,
        // between 2exp_3 (2 to the 3rd) and 2exp_7.
        // (The program would seem to halt with a 2exp_8 matrix on vxworks).
        
        // Computing 2exp_r (2 to the r)
        int n = 1;
        for (int i = 1; i <= r; i++)
        	n = n*2;
        
        printf("Multiplying a %d-by-%d matrix\n", (int) n, (int) n);
        
        _a = matrixGenerator(n);
        _b = matrixGenerator(n);
        

        // Create hook
        if ((taskSwitchHookAdd((FUNCPTR)switchHook)) != OK){
         	printf("Error creating the switch hook\n");
        }
        
        // Initialize arrays where the timestamps will be added
        initializeArrays();
        
        // Code from initial sample file (I think the name was cmk.c) and adapted for
        // my program:
        
        id0 = taskSpawn("Task0", 210, 0, 4096, (FUNCPTR) IJK, 0, n, 0, 0, 0, 0, 0, 0, 0, 0);
        id1 = taskSpawn("Task1", 210, 0, 4096, (FUNCPTR) JKI, 1, n, 0, 0, 0, 0, 0, 0, 0, 0);
        id2 = taskSpawn("Task2", 210, 0, 4096, (FUNCPTR) KIJ, 2, n, 0, 0, 0, 0, 0, 0, 0, 0);
        id3 = taskSpawn("Task3", 230, 0, 4096, (FUNCPTR) Strassen, 3, n, 0, 0, 0, 0, 0, 0, 0, 0);
        guest_id = taskSpawn("GuestTask", 250, 0, 4096, (FUNCPTR) guest, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        
}
