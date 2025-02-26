#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

long long K, L, N;



// Function to generate a random float between 0 and 1
float rand_float() {
    return (float)rand() / (float)RAND_MAX;
}

// Function to sort an array of floats
void sort_floats(float arr[], int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                float temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

// Split function (Generate Row)
void split(int lamb, int m, int values[]) {
    float breaks[K];
    for (int i = 0; i < m - 1; i++) {
        breaks[i] = rand_float() * lamb;
    }
    sort_floats(breaks, m - 1);

    values[0] = round(breaks[0]);
    for (int i = 1; i < m - 1; i++) {
        values[i] = round(breaks[i] - breaks[i-1]);
    }
    values[m-1] = round(lamb - breaks[m-2]);

    // Ensure no zeros and adjust to maintain sum
    for (int i = 0; i < m; i++) {
        if (values[i] == 0) {
            values[i] = 1;
            int j = (i + 1) % m;
            while (values[j] <= 1) {
                j = (j + 1) % m;
            }
            values[j]--;
        }
    }
}

// Goal split function (Generate 3 Solution Rows)
void goal_split(int array[3][K]) {
    int row_sum = K * L / 3;

    for (int j = 0; j < K; j++) {
        int values[3] = {1, 1, 1};
        int remaining = L - 3;
        for (int _ = 0; _ < remaining; _++) {
            values[rand() % 3]++;
        }
        for (int i = 0; i < 3; i++) {
            array[i][j] = values[i];
        }
    }

    for (int i = 0; i < 2; i++) {
        int diff = row_sum;
        for (int j = 0; j < K; j++) {
            diff -= array[i][j];
        }
        while (diff != 0) {
            if (diff > 0) {
                int col = rand() % K;
                if (array[i][col] < L - 1 && array[2][col] > 1) {
                    array[i][col]++;
                    array[2][col]--;
                    diff--;
                }
            } else {
                int col = rand() % K;
                if (array[i][col] > 1 && array[2][col] < L - 1) {
                    array[i][col]--;
                    array[2][col]++;
                    diff++;
                }
            }
        }
    }
}

// Generate files function
void generate_files() {

    //Generates Rows To Be Solution
    int sol_rows[3];
    for (int i = 0; i < 3; i++) {
        sol_rows[i] = rand() % N;
    }

    //Generates Solution
    int sol_split[3][K];
    goal_split(sol_split);

    //Files To Create
    FILE *files[3];
    char *filenames[] = {"A.txt", "B.txt", "C.txt"};

    for (int i = 0; i < 3; i++) {

        //Opening File
        printf("\n BUILDING: %s: ", filenames[i]);
        files[i] = fopen(filenames[i], "w");
        if (files[i] == NULL) {
            printf("Error opening file %s\n", filenames[i]);
            return;
        }

        //Adding Rows
        for (int row = 0; row < N; row++) {

            //Progress
            if (row % N/100 == 0){printf("P");}

            //Adding To File
            if (row == sol_rows[i]) {
                for (int j = 0; j < K; j++) {
                    fprintf(files[i], "%d ", sol_split[i][j]);
                }
            } else {
                int values[K];
                split(L, K, values);
                for (int j = 0; j < K; j++) {
                    fprintf(files[i], "%d ", values[j]);
                }
            }
            fprintf(files[i], "\n");
        }

        fclose(files[i]);
    }

    printf("\n Solution Rows: %d %d %d\n", sol_rows[0], sol_rows[1], sol_rows[2]);
}


int main(int argc, char *argv[]) {
    srand(time(NULL));

    K = atoll(argv[1]);
    L = atoll(argv[2]);
    N = atoll(argv[3]);

    generate_files();

    return 0;
}