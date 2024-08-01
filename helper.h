/* 
 * Student: Axel Rolando Alvarenga Munoz
 * UH ID: 2075248
 * COSC 4331 Homework 2 - Helper Functions file
 * 
 * 
 * Disclosures/References:
 * [1] Declaration and memory allocation for double pointer. 
 * 	   Credit: https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/#
 * 	   
 * 
 *  See report document for other references.
 *  
 */

double** matrixGenerator(double _n){
    double** matrix = (double**) malloc(_n*sizeof(double*));				// Credit: [1]

    for (int i = 0; i < _n;i++)												// Credit: [1]
        matrix[i] = (double*) malloc(_n*sizeof(double));					// Credit: [1]

    for (int i = 0; i < _n;i++){
        for (int j = 0; j < _n;j++)
            matrix[i][j] = rand();
    }
    
    return matrix;
}

double** pointerCGenerator (int _n) {
    double** matrix = (double**) malloc(_n*sizeof(double*));				// Credit: [1]

    for (int i = 0; i < _n;i++)												// Credit: [1]
        matrix[i] = (double*) malloc(_n*sizeof(double));					// Credit: [1]

    for (int i = 0; i < _n;i++){
        for (int j = 0; j < _n;j++)
            matrix[i][j] = 0;
    }

    return matrix;
}
