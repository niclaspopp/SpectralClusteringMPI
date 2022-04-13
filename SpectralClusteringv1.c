#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>

#define MATSIZE 1024
#define k 4 // number of clusters
#define d 4 // number of dimensions for k-means
#define MASTER 0               // taskid of fMaster
#define FROM_MASTER 1          // setting a message type
#define FROM_WORKER 2          // setting a message type 

// ----------------------------------------------------------------
// Functions for Initialisation
// ----------------------------------------------------------------
double **readmatrix(size_t *rows, size_t *cols, const char *filename){
    if(rows == NULL || cols == NULL || filename == NULL)
        return NULL;

    *rows = 0;
    *cols = 0;

    FILE *fp = fopen(filename, "r");

    if(fp == NULL)
    {
        fprintf(stderr, "could not open %s: %s\n", filename, strerror(errno));
        return NULL;
    }

    double **matrix = NULL, **tmp;

    char line[1024];

    while(fgets(line, sizeof line, fp))
    {
        if(*cols == 0)
        {
            // determine the size of the columns based on
            // the first row
            char *scan = line;
            double dummy;
            int offset = 0;
            while(sscanf(scan, "%lf%n", &dummy, &offset) == 1)
            {
                scan += offset;
                (*cols)++;
            }
        }

        tmp = realloc(matrix, (*rows + 1) * sizeof *matrix);

        if(tmp == NULL)
        {
            fclose(fp);
            return matrix; // return all you've parsed so far
        }

        matrix = tmp;

        matrix[*rows] = calloc(*cols, sizeof *matrix[*rows]);

        if(matrix[*rows] == NULL)
        {
            fclose(fp);
            if(*rows == 0) // failed in the first row, free everything
            {
                fclose(fp);
                free(matrix);
                return NULL;
            }

            return matrix; // return all you've parsed so far
        }

        int offset = 0;
        char *scan = line;
        for(size_t j = 0; j < *cols; ++j)
        {
            if(sscanf(scan, "%lf%n", matrix[*rows] + j, &offset) == 1)
                scan += offset;
            else
                matrix[*rows][j] = 0; // could not read, set cell to 0
        }

        // incrementing rows
        (*rows)++;
    }

    fclose(fp);

    return matrix;
}

float* vectorize(const int num_elements, double M[MATSIZE][k]){
    float *vec = (float *)malloc(sizeof(float) * num_elements);
    for (int i=0; i<MATSIZE; i++){
        for (int j=0; j<d; j++){
            vec[i*d+j]=M[i][j];
        }
    }
    return vec;
}

double** alloc_matrix(int m, int n)
{
    double**  rows = NULL;
    double*   data = NULL;
    int i;

    assert( (m >= 0) && (n >= 0) );

    if ((m == 0)||(n == 0))
    {
        return rows;
    }

    rows = (double**) malloc ( m *     sizeof(double*) );
    data = (double*)  malloc ( m * n * sizeof(double)  );

    if (!rows||!data){
        fprintf ( stderr, "Out of memory.\n" );
        exit    ( 1 );
    }

    for ( i = 0; i < m; i++ ){
        rows[i] = &data[i * n];
    }

    return rows;
}

double* vectorize_alloc(double **M){
    double *vec = (double *)malloc(sizeof(double) * (MATSIZE*MATSIZE));
    for (int i=0; i<MATSIZE; i++){
        for (int j=0; j<MATSIZE; j++){
            vec[i*MATSIZE+j]=M[i][j];
        }
    }
    return vec;
}

// ----------------------------------------------------------------
// Functions for Odd-Even Sort
// ----------------------------------------------------------------
// Helper function to exchange two elements
void exchange(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

//split the list into two lists, smaller and larger than pivot
int split (double arr[], int startIndex, int EndIndex)
{
    double pivot = arr[EndIndex];    // pivot
    int i = (startIndex - 1);  // Index of smaller element

    for (int j = startIndex; j <= EndIndex - 1; j++)
    {
        if (arr[j] <= pivot)
        {
            i++;    // increase the index of the smaller element
            exchange(&arr[i], &arr[j]);
        }
    }
    exchange(&arr[i + 1], &arr[EndIndex]);
    return (i + 1);
}

//startIndex: starting index, EndIndex: end index
void QuickSort(double arr[], int startIndex, int EndIndex)
{
    if (startIndex < EndIndex)
    {
        // si
        int si = split(arr, startIndex, EndIndex);
        // Separately sort elements before split and after the split
        QuickSort(arr, startIndex, si - 1);
        QuickSort(arr, si + 1, EndIndex);
    }
}

// Sort the larger half of the array after each merging step
void SortLargerHalf(double *x, double *a, int nx, int na){
  double temp[nx];
  int i = nx-1;
  int j = na-1;
  for (int n = nx-1; n >= 0; --n) {
    if (x[i] >= a[j]) {
      temp[n] = x[i];
      i--;
    }
    else {
      temp[n] = a[j];
      j--;
    }
  }
  for (int n = 0; n < nx; ++n) {
    x[n] = temp[n];
  }
  return;
}

// Sort the smaller half of the array after each merging step
void SortSmallerHalf(double *x, double *a, int nx, int na){
  double temp[nx];
  int i = 0;
  int j = 0;

  for (int n = 0; n < nx; ++n) {
    if (x[i] < a[j]) {
      temp[n] = x[i];
      i++;
    }else{
      temp[n] = a[j];
      j++;
    }
 
  }

  for (int n = 0; n < nx; ++n) {
    x[n] = temp[n];
  }
  
  return;
}

// ----------------------------------------------------------------
// Functions for k-means
// ----------------------------------------------------------------
// Distance**2 between d-vectors pointed to by v1, v2.
float distance2(const float *v1, const float *v2, const int de) {
  float dist = 0.0;
  for (int i=0; i<de; i++) {
    float diff = v1[i] - v2[i];
    dist += diff * diff;
  }
  return dist;
}

// Assign a site to the correct cluster by computing its distances to each cluster centroid.
int assign_site(const float* site, float* centroids, const int ka, const int de) {
  int best_cluster = 0;
  float best_dist = distance2(site, centroids, de);
  float* centroid = centroids + d;
  for (int c = 1; c < ka; c++, centroid += de) {
    float dist = distance2(site, centroid, de);
    if (dist < best_dist) {
      best_cluster = c;
      best_dist = dist;
    }
  }
  return best_cluster;
}

// Add a site (vector) into a sum of sites (vector).
void add_site(const float * site, float * sum, const int de) {
  for (int i=0; i<de; i++) {
    sum[i] += site[i];
  }
}

// Print the centroids one per line.
void print_centroids(float * centroids, const int ka, const int de) {
  float *p = centroids;
  printf("Centroids:\n");
  for (int i = 0; i<ka; i++) {
    for (int j = 0; j<de; j++, p++) {
      printf("%f ", *p);
    }
    printf("\n");
  }
}

// ----------------------------------------------------------------
// MAIN
// ----------------------------------------------------------------
int main(int argc, char *argv[]){
    
    int numtasks,taskid,numworkers,source,dest,mtype,rows,averow,extra,offset,m;
    int i,j,rc,u,q,l,K;
    //double A[MATSIZE][MATSIZE],ATEMP[MATSIZE][MATSIZE],Aj[MATSIZE][MATSIZE],R[MATSIZE][MATSIZE],Q[MATSIZE][MATSIZE];
    double **A, **ATEMP, **Aj, **R, **Q;
    double qg[MATSIZE],rg[MATSIZE],r[MATSIZE], eigvals[MATSIZE];
    double t,normsq;
    double *Avec;
    double startQR,startSort,startKM,endQR,endSort,endKM; // used for timing

    A = alloc_matrix(MATSIZE,MATSIZE);
    ATEMP = alloc_matrix(MATSIZE,MATSIZE);
    Aj = alloc_matrix(MATSIZE,MATSIZE);
    R = alloc_matrix(MATSIZE,MATSIZE);
    Q = alloc_matrix(MATSIZE,MATSIZE);
    


    MPI_Init(NULL,NULL);
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    K = MATSIZE/numtasks;
    int localsize = K*MATSIZE;
    double qglocal[K],rglocal[K];
    double *aj;
    aj = malloc(sizeof(double)*localsize);
    
    if (numtasks < 2 ) {
        printf("Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }
    if (MATSIZE%numtasks != 0) {
        printf("Number of processes has to divide number of nodes. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }
    numworkers = numtasks-1;

    // ----------------------------------------------------------------
    // initialization
    // ----------------------------------------------------------------
    
    if (taskid == MASTER){
        printf("Spectral clustering with %d tasks.\n",numtasks);
        // printf("Reading matix...\n");
        size_t cols_r, rows_r;
        double **matrix = readmatrix(&rows_r, &cols_r, "Laplacian_test.dat");
        for (i=0; i<MATSIZE; i++){
            eigvals[i]=1;
            for (l=0; l<MATSIZE; l++){
                A[i][j] = matrix[i][j];
                ATEMP[i][j] = A[i][j];
                R[i][j] = 0;
                Q[i][j] = 0;
            }
        }
    }

    
    // ----------------------------------------------------------------
    // QR algorithm with Double Gram Schmidt
    // ----------------------------------------------------------------
    // Start QR timing
    MPI_Barrier(MPI_COMM_WORLD); 
    startQR = MPI_Wtime();
    for(u=0; u<2; u++){
        // Get candidate for q vector
        if (taskid == MASTER)
        {
            // Initialize Aj
            for (i=0; i<MATSIZE; i++){
                for (j=0; j<MATSIZE; j++){
                Aj[i][j]=ATEMP[i][j];
                }
            }
            // get q vector
            for (j = 0; j < MATSIZE; j++ ){
                qg[j] = Aj[j][i];
            }
            normsq = 0;
        }
        
        // normalize q vector, calculcate norm in parallel
        MPI_Scatter(qg,K, MPI_DOUBLE, qglocal, K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double localnorm = 0;
        for (j = 0; j < K; j++ ){
                localnorm = localnorm + qglocal[j]*qglocal[j];
                
        }
        MPI_Reduce(&localnorm, &normsq, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        if (taskid == MASTER)
        {
            for (j = 0; j < MATSIZE; j++ ){
                qg[j] = qg[j]/sqrt(normsq);
            }
        }
        
        // calculate r
        // first broadcast q
        MPI_Bcast(&qg, MATSIZE, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        // Scatter Aj
        if (taskid == MASTER)
        {
            Avec = vectorize_alloc(Aj);  
        }
        MPI_Scatter(Avec, K*MATSIZE, MPI_DOUBLE, aj, K*MATSIZE, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        // Calculate r vector (distributed)
        for (j = 0; j < K; j++ ){
            rglocal[j] = 0;
            for (l = 0; l < MATSIZE; l++ ){
                rglocal[j] = rglocal[j] + aj[j*MATSIZE + l]*qg[l];
            }
        }
        // gather r
        //printf("Process %d, my value = %f.\n", taskid, rglocal[0]);
        // Determine root's rank
        int root_rank = 0;
    
        // Get my rank
        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
        // Define my value
        double my_value = rglocal[0] * 100;
    
        if(my_rank == root_rank)
        {
            double rbuffer[MATSIZE];
            MPI_Gather(rglocal, K, MPI_DOUBLE, rbuffer, K, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
            for (l = 0; l < MATSIZE; l++ ){
                rg[l]=rbuffer[l];
            } 
        }
        else
        {
            MPI_Gather(rglocal, K, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
        }
        // go on with Gram Schmidt
        if(taskid == MASTER)
        {
                // Update Q and R
                for (j = 0; j < MATSIZE; j++ ){
                    Q[j][i] = qg[j];
                    R[i][j] = rg[j];
                }
                // Update Aj
                for (j=0; j<MATSIZE; j++){
                    for (l=0; l<MATSIZE; l++){
                        Aj[l][j]=Aj[l][j]-qg[l]*rg[j];
                    }
                }
                // Update i-th eignvalue
                eigvals[i]=eigvals[i]*R[i][i];
        }

        // Gram-Schidt finished --------------------------------------------------
        // Parallelized matrix multiplication-------------------------------------
        if(taskid == MASTER){
            averow = MATSIZE/numworkers;
            extra = MATSIZE%numworkers;
            offset = 0;
            mtype = FROM_MASTER;
            
            for (dest=1; dest<=numworkers; dest++)
            {
                rows = (dest <= extra) ? averow+1 : averow;   	
                // printf("Sending %d rows to task %d offset=%d\n",rows,dest,offset);
                MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                MPI_Send(&(R[offset][0]), rows*MATSIZE, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
                MPI_Send(&Q[0][0], MATSIZE*MATSIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
                offset = offset + rows;
            }
            mtype = FROM_WORKER;
            for (i=1; i<=numworkers; i++){
                source = i;
                MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&ATEMP[offset][0], rows*MATSIZE, MPI_DOUBLE, source, mtype, 
                MPI_COMM_WORLD, &status);
                // printf("Received results from task %d\n",source);
            }
        }
        // Worker task
        if (taskid > MASTER)
        {
            
            mtype = FROM_MASTER;
            MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&(R[0][0]), rows*MATSIZE, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&Q[0][0], MATSIZE*MATSIZE, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
            for (q=0; q<MATSIZE; q++){
                for (i=0; i<rows; i++){
                    ATEMP[i][q] = 0.0;
                    for (j=0; j<MATSIZE; j++){
                        ATEMP[i][q] = ATEMP[i][q] + R[i][j] * Q[j][q];
                    }
                }
            }
            mtype = FROM_WORKER;
            MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
            MPI_Send(&ATEMP[0][0], rows*MATSIZE, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
            
        }
        
    }
    
    // stop QR timing
    MPI_Barrier(MPI_COMM_WORLD); 
    endQR = MPI_Wtime();
    if (taskid == MASTER){
        printf("QR xecution time: %e\n",endQR-startQR);
        free(Avec);
    }

    // ----------------------------------------------------------------
    // Odd-Even Sort of the eigenvalues to get the correct eigenvectors
    // ----------------------------------------------------------------
    // Start sorting timing
    MPI_Barrier(MPI_COMM_WORLD); 
    startSort = MPI_Wtime();
    int P, p, NN, I, step, odd;
    double start,end; // used for timing
    bool evenprocess, evenphase;
    double *x;
    double *a;
    double *send;
    double *receive;
    int *indices;
    int *x_inds;
    int tag = 42;

    // number of ranks: P
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    // current rank: p
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

    // Find problem size N from command line 
    NN = MATSIZE;

    // local size
    //I = (NN+P-p-1)/P;
    I = MATSIZE/P;
    
    // initialization of random number generator
    srandom(p+1);

    a = malloc(sizeof(x) * I +1);
    x = malloc(sizeof(x) * I);
    x_inds = malloc(sizeof(x) * I);

    // Populate index array
    
    if (taskid == MASTER){
        indices = malloc(sizeof(x) * MATSIZE);
        send = malloc(sizeof(x) * MATSIZE);
        for (j=0; j<MATSIZE; j++){
            indices[j]=j;
            send[j]=eigvals[j];
        }       
    }

    
    // Root sends each process its share of the eigenvalues and the indices
    MPI_Scatter(send,I, MPI_DOUBLE, x, I, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(indices,I, MPI_INT, x_inds, I, MPI_INT, 0, MPI_COMM_WORLD);

    // start with local sorting
    QuickSort(x, 0, I-1);

    // check if process is odd or even
    odd = P % 2;
    evenprocess = ((p % 2) == 0);

    // start with even phase
    evenphase = 1;

    // perform odd-even sorting
    
    for (step = 0; step < P; step++) {
        if (evenphase == evenprocess) {
            if (p != P-1) {
                MPI_Send(x, I, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
                MPI_Recv(a, I, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
                SortSmallerHalf(x, a, I, I);
            }
        }
        else {
            if (p != 0) {
                MPI_Recv(a, I, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
                MPI_Send(x, I, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
                SortLargerHalf(x, a, I, I);
            }
        }
        evenphase = !evenphase;
    }

    MPI_Gather(x,I, MPI_DOUBLE, send, I, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // stop sorting timing
    MPI_Barrier(MPI_COMM_WORLD); 
    endSort = MPI_Wtime();
    if (taskid == MASTER)
        printf("Sorting execution time: %e\n",endSort-startSort);

    // ----------------------------------------------------------------
    // k-means
    // ----------------------------------------------------------------
    // Start k-means timing
    MPI_Barrier(MPI_COMM_WORLD); 
    startKM = MPI_Wtime();

    int sites_per_proc = MATSIZE/numtasks;
    
    // Alloc memory for sites, sums, counts, centroids, labales
    float* sites;  
    sites = malloc(sites_per_proc * d * sizeof(float));
    float* sums;
    sums = malloc(k * d * sizeof(float));
    int* counts;
    counts = malloc(k * sizeof(int));
    float* centroids;
    centroids = malloc(k * d * sizeof(float));
    int* labels;
    labels = malloc(sites_per_proc * sizeof(int));

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Data structures maintained only in root process.
    float* all_sites = NULL;
    float* grand_sums = NULL;
    int* grand_counts = NULL;
    int* all_labels;
    
    if (rank == 0) {
        // Update Aj
        for (i=0; i<MATSIZE; i++){
            for (j=0; j<MATSIZE; j++){
                Aj[i][j]=ATEMP[i][j];
            }
        }
        // Get U Matrix form trace version of Rayleigh-Ritz
        double U[MATSIZE][k];
        for (i=0; i<MATSIZE; i++){
            for (j=0; j<k; j++){
                U[i][j]=Aj[i][MATSIZE-k+j];
            }
        }

        all_sites = vectorize(d * sites_per_proc * nprocs,U);

        // Take the first k sites as the initial cluster centroids.
        for (int i = 0; i < k * d; i++) {
            centroids[i] = all_sites[i]; 
        }
        //print_centroids(centroids, k, d);
        grand_sums = malloc(k * d * sizeof(float));
        grand_counts = malloc(k * sizeof(int));
        all_labels = malloc(nprocs * sites_per_proc * sizeof(int));
    }

    // Root sends each process its share of sites.
    MPI_Scatter(all_sites,d*sites_per_proc, MPI_FLOAT, sites, d*sites_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    float norm = 1.0;  // Will tell us if centroids have moved.
  
    while (norm > 0.00001) { // While they've moved...

        // Broadcast the current cluster centroids to all processes.
        MPI_Bcast(centroids, k*d, MPI_FLOAT,0, MPI_COMM_WORLD);

        // Each process reinitializes its cluster accumulators.
        for (int i = 0; i < k*d; i++) sums[i] = 0.0;
        for (int i = 0; i < k; i++) counts[i] = 0;

        // Find the closest centroid to each site and assign to cluster.
        float* site = sites;
        for (int i = 0; i < sites_per_proc; i++, site += d) {
        int cluster = assign_site(site, centroids, k, d);
        // Record the assignment of the site to the cluster.
        counts[cluster]++;
        add_site(site, &sums[cluster*d], d);
        }

        // Gather and sum at root all cluster sums for individual processes.
        MPI_Reduce(sums, grand_sums, k * d, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(counts, grand_counts, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
        // Root process computes new centroids by dividing sums per cluster
        // by count per cluster.
        for (int i = 0; i<k; i++) {
            for (int j = 0; j<d; j++) {
                int dij = d*i + j;
                grand_sums[dij] /= grand_counts[i];
            }
        }
        // Have the centroids changed much?
        norm = distance2(grand_sums, centroids, d*k);
        
        //printf("norm: %f\n",norm);
        // Copy new centroids from grand_sums into centroids.
        for (int i=0; i<k*d; i++) {
            centroids[i] = grand_sums[i];
        }
        //print_centroids(centroids,k,d);
        }
        // Broadcast the norm.  All processes will use this in the loop test.
        MPI_Bcast(&norm, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    // Now centroids are fixed, so compute a final label for each site.
    float* site = sites;
    for (int i = 0; i < sites_per_proc; i++, site += d) {
        labels[i] = assign_site(site, centroids, k, d);
    }

    // Gather all labels into root process.
    MPI_Gather(labels, sites_per_proc, MPI_INT, all_labels, sites_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
    // stop k-means timing
    MPI_Barrier(MPI_COMM_WORLD); 
    endKM = MPI_Wtime();
    if (taskid == MASTER)
        printf("k-means execution time: %e\n",endKM-startKM);

    MPI_Finalize();
    free(A);free(Aj);free(ATEMP);free(Q);free(R);
    return 0;
}