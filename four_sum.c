


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <math.h>



long long K,L,M,N;

#define MAX_PAIRS 100000

long long agg, nagg, value;

//First 1000 Primes To Serve As Coefficients
long long primes[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919
};






#define HASH_SIZE 100003  // Large prime number for hash table size

// Hash table node
typedef struct HashNode {
    long long key;
    struct HashNode* next;
} HashNode;

// Hash table array
HashNode* hash_table[HASH_SIZE] = {NULL};

// Hash function
unsigned int hash_function(long long key) {
    return (unsigned int)(key % HASH_SIZE);
}

// Insert a value into the hash table
void insert_hash(long long key) {
    unsigned int index = hash_function(key);

    // Check if key already exists
    HashNode* node = hash_table[index];
    while (node) {
        if (node->key == key) {
            return;  // Already exists, no need to insert again
        }
        node = node->next;
    }

    // Create new hash node
    HashNode* new_node = (HashNode*)malloc(sizeof(HashNode));
    new_node->key = key;
    new_node->next = hash_table[index];
    hash_table[index] = new_node;
}

// Search for a value in the hash table
int search_hash(long long key) {
    unsigned int index = hash_function(key);
    HashNode* node = hash_table[index];
    while (node) {
        if (node->key == key) return 1;  // Found
        node = node->next;
    }
    return 0;  // Not found
}

// Clear the hash table (free memory)
void clear_hash_table() {
    for (int i = 0; i < HASH_SIZE; i++) {
        HashNode* node = hash_table[i];
        while (node) {
            HashNode* temp = node;
            node = node->next;
            free(temp);
        }
        hash_table[i] = NULL;
    }
}




// 4SUM function with IDs
int four_sum_with_ids(long long **A_rows, long long A_counts, 
    long long **B_rows, long long B_counts, 
    long long **C_rows, long long C_counts,
    long long **D_rows, long long D_counts
    , int target) {

        //To Mimizize hash FILLING (create an upper and lower bound)
        long long min_ab = target - C_rows[C_counts - 1][1] - D_rows[D_counts - 1][1];
        long long max_ab = target - C_rows[0][1] - D_rows[0][1];

        printf("In");


        // for (int a = 0; a < A_counts; a++){
        //     for (int b = 0; b < B_counts; b++){
        //         for (int c = 0; c < C_counts; c++){
        //             for (int d = 0; d < D_counts; d++){
        //                 if (A_rows[a][1] + B_rows[b][1] + C_rows[c][1] + D_rows[d][1] == target){
        //                     printf("    SOLUTION: %lld %lld %lld %lld  ", A_rows[a][0], B_rows[b][0], C_rows[c][0], D_rows[d][0]);
        //                 }
        //             }
        //         }
        //     }
        // }

        // Step 1: Store all (A + B) sums in the hash table
    for (int a = 0; a < A_counts; a++) {
        for (int b = 0; b < B_counts; b++) {
            long long sum_ab = A_rows[a][1] + B_rows[b][1];
            if ((sum_ab >= min_ab) & (sum_ab <= max_ab)){
                insert_hash(sum_ab);
            }
        }
    }

    // Step 2: Check if (target - (C + D)) exists in the hash table
    for (int c = 0; c < C_counts; c++) {
        for (int d = 0; d < D_counts; d++) {
            long long sum_cd = C_rows[c][1] + D_rows[d][1];
            long long complement = target - sum_cd;

            // If complement exists, print solution
            if (search_hash(complement)) {
                printf("SOLUTION FOUND (*): %lld + %lld + %lld + %lld = %d\n", 
                       C_rows[c][1], D_rows[d][1], complement, sum_cd, target);
            }
        }
    }

    // Clean up the hash table
    clear_hash_table();



 
}



















//GOOD
int compare(const void *a, const void *b) {
    long long val1 = (*(long long**)a)[1];
    long long val2 = (*(long long**)b)[1];

    if (val1 < val2) return -1;
    if (val1 > val2) return 1;
    return 0;

}








void solve(long long ****A_rows, long long **A_counts, 
            long long ****B_rows, long long **B_counts, 
            long long ****C_rows, long long **C_counts,
            long long ****D_rows, long long **D_counts) {
     
    long long goal_agg = 0;
    for (int x = 0; x < K; x++){
        goal_agg += (primes[x]*L);
    }    
    printf("\nGOAL: %lld\n", goal_agg);
    
    long long goal_mod = goal_agg % M;

    long long an,bn,cn,dn;
    long long **ar;
    long long **br;
    long long **cr;
    long long **dr;  

    printf("\nCoverage:\n");
    for (int y = 0; y < M; y++){
        printf("\n---------");
        for (int z = 0; z < M; z++){
            printf("%lld, %lld, %lld, %lld\n",A_counts[y][z],B_counts[y][z],C_counts[y][z],D_counts[y][z]);
        }
    }

    printf("\nExact\n");
    for (int x = 0; x < M; x++){
        for (int z = 0; z < M; z++){
            for (int y = 0; y < A_counts[x][z]; y++){
                printf("A: %lld %lld %lld   %lld %lld\n", x, z, y, A_rows[x][z][y][0],A_rows[x][z][y][1]);
            }
        }
    }
    for (int x = 0; x < M; x++){
        for (int z = 0; z < M; z++){
            for (int y = 0; y < B_counts[x][z]; y++){
                printf("B: %lld %lld %lld   %lld %lld\n", x, z, y, B_rows[x][z][y][0],B_rows[x][z][y][1]);
            }
        }
    }
    for (int x = 0; x < M; x++){
        for (int z = 0; z < M; z++){
            for (int y = 0; y < C_counts[x][z]; y++){
                printf("C: %lld %lld %lld   %lld %lld\n", x, z, y, C_rows[x][z][y][0],C_rows[x][z][y][1]);
            }
        }
    }
    for (int x = 0; x < M; x++){
        for (int z = 0; z < M; z++){
            for (int y = 0; y < D_counts[x][z]; y++){
                printf("D: %lld %lld %lld   %lld %lld\n", x, z, y, D_rows[x][z][y][0],D_rows[x][z][y][1]);
            }
        }
    }



    int theoretical = 0;
    int actual = 0;

    for (long long a1 = 0; a1 < M; a1++){
        for (long long b1 = 0; b1 < M; b1++){
            for (long long c1 = 0; c1 < M; c1++){
                for (long long d1 = 0; d1 < M; d1++){
                    
                    if (((a1 + b1 + c1 + d1) % M) == goal_agg % M){

                        for (long long a2 = 0; a2 < M; a2++){
                            for (long long b2 = 0; b2 < M; b2++){
                                for (long long c2 = 0; c2 < M; c2++){
                                    for (long long d2 = 0; d2 < M; d2++){


                                        if (1==1){
                                            if (((a2 + b2 + c2 + d2) % M) == goal_agg % M){
                                                actual += 1;
                                                if ((A_counts[a1][a2] > 0) & (B_counts[b1][b2] > 0) &
                                                    (C_counts[c1][c2] > 0) & (D_counts[d1][d2] > 0)){
                                                        printf("ABOUT %lld %lld %lld", a1, b1, c1, d1);
                                                four_sum_with_ids(A_rows[a1][a2], A_counts[a1][a2],
                                                                  B_rows[b1][b2], B_counts[b1][b2],
                                                                  C_rows[c1][c2], C_counts[c1][c2],
                                                                  D_rows[d1][d2], D_counts[d1][d2],
                                                                  goal_agg);
                                                
                                                    }

                                                
                                            }
                                        }
                                        theoretical += 1;

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("%lld %lld", theoretical, actual);
  


    // int f = 0;
    // int g = 0;

    // // // For every A section
    // for (long long am = 0; am < M; am++) {
    //     an = A_counts[am];
    //     ar = A_rows[am];  
 
    //    long long mod_rem = (goal_mod - am + M) % M;

    //     // Find Bm, Cm such that remainders line up
    //     for (long long bm = 0; bm < M; bm++) {
    //         bn = B_counts[bm];
    //         br = B_rows[bm];  
    //         for (long long cm = 0; cm < M; cm++) {
    //             cn = C_counts[cm];
    //             cr = C_rows[cm];  

    //             for (long long dm = 0; dm < M; dm++){
    //                 dn = D_counts[dm];
    //                 dr = D_rows[dm];  

    //                 if (((bm + cm + dm) % M) == mod_rem) {
    //                     f += 1;
    //                     if ((an > 0) & (bn > 0) & (cn > 0) & (dn > 0)){
                            
    //                         four_sum_with_ids(ar, an, br, bn, cr, cn, dr, dn, goal_agg);
    //                     } 
    //                 } 
    //                 g += 1;
    //             }
    //         }
    //     }
    // }

    // printf("LENGTH: %d  %d", f, g);
 }








//GOOD
void setup(long long *****rows, long long ***counts, const char *filename){

    //Long Term Storage
    *rows = (long long ****)malloc(M*sizeof(long long ***));
    *counts = (long long **)malloc(M*sizeof(long long *));

    for (int i = 0; i < M; i++) {
        (*counts)[i] = (long long *)calloc(M, sizeof(long long));
        (*rows)[i] = (long long ***)malloc(M * sizeof(long long **));
    }

    //Temp Storage
    long long **temp = (long long **)malloc(N * sizeof(long long*));
    for (int i = 0; i < N; ++i){
        temp[i] = (long long *)malloc(3 * sizeof(long long));
    }

    long long **group_pos = (long long **)malloc(M * sizeof(long long *));
    for (int i = 0; i < M; i++){
        group_pos[i] = (long long *)malloc(M * sizeof(long long));
    }

    //File
    FILE *file = fopen(filename, "r");

    printf("\nOpening File:");

    //Move File To Temp
    for (long long n = 0; n < N; n++) {

        //Dispay
        if (n % N/100 == 0){printf("P");}

        printf(" A");

        //Determine Agg
        agg = 0;
        nagg = 0;
        for (long long k = 0; k < K; k++) {
            fscanf(file, "%lld", &value); 
            agg += (primes[k] * value); 
            nagg += (primes[K-k-1] * value);
        }

        //Set Data
        temp[n][0] = n; 
        temp[n][1] = agg;  
        temp[n][2] = nagg; 
        
        (*counts)[agg % M][nagg % M] += 1; 

    }
    fclose(file);  

    

    //Size Long Term
    for (int m = 0; m < M; m++){
        for (int n = 0; n < M; n++){
            (*rows)[m][n] = (long long **)malloc((*counts)[m][n] * sizeof(long long *));
        }
        
    }
    

    //Fill Long Term
    int group, ngroup, size;
    for (int n = 0; n < N; n++) {
        
        //Display
        if (n % N/50 == 0){printf("P");}

        printf(" 1");
        
        group = temp[n][1] % M;
        ngroup = temp[n][2] % M;
        size = group_pos[group][ngroup];

        printf("2:_%lld_%lld_%lld_",size, group, ngroup);

        (*rows)[group][ngroup][size] = (long long *)malloc(3 * sizeof(long long));
        (*rows)[group][ngroup][size][0] = temp[n][0];  
        (*rows)[group][ngroup][size][1] = temp[n][1];  

        

        group_pos[group][ngroup] += 1;

        printf("4\n");


    }

    printf(" A");

   

    //TODO: requires for range cutoffs for ab (worth it?????)
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < M; n++){
            qsort((*rows)[m][n], (*counts)[m][n], sizeof(long long *), compare);
        }
        
    }

    printf(" B");


    for (int x = 0; x < N; x++){
        free(temp[x]);
        free(group_pos[x]);
    }
    free(temp);
    free(group_pos);
}

void cleanup(long long *****rows, long long ***counts) {
    if (*rows) {
        for (int i = 0; i < M; i++) {
            if ((*rows)[i]) {
                for (int j = 0; j < M; j++) {
                    if ((*rows)[i][j]) {
                        for (int k = 0; k < (*counts)[i][j]; k++) {
                            free((*rows)[i][j][k]);  // Free each row entry
                        }
                        free((*rows)[i][j]);  // Free row container
                    }
                }
                free((*rows)[i]);  // Free outer array
            }
        }
        free(*rows);  // Free the top-level pointer
        *rows = NULL;
    }

    if (*counts) {
        for (int i = 0; i < M; i++) {
            free((*counts)[i]);  // Free counts array for each row
        }
        free(*counts);
        *counts = NULL;
    }
}



int main(int argc, char *argv[]) {
    srand(time(NULL));

    K = atoll(argv[1]);
    L = atoll(argv[2]);
    M = atoll(argv[3]);
    N = atoll(argv[4]);

    long long ****A_rows;
    long long ****B_rows;
    long long ****C_rows;
    long long ****D_rows;

    long long **A_counts;
    long long **B_counts;
    long long **C_counts;
    long long **D_counts;

    setup(&A_rows, &A_counts, "A.txt");
    setup(&B_rows, &B_counts, "B.txt");
    setup(&C_rows, &C_counts, "C.txt");
    setup(&D_rows, &D_counts, "D.txt");

    solve(A_rows, A_counts, B_rows, B_counts, C_rows, C_counts, D_rows, D_counts);

    cleanup(&A_rows, &A_counts);
    cleanup(&B_rows, &B_counts);
    cleanup(&C_rows, &C_counts);
    cleanup(&D_rows, &D_counts);
}

