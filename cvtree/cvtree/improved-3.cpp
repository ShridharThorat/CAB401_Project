#include "cvtree_parallel.hpp"
using namespace std;

int number_bacteria;
char** bacteria_name;
long M, M1, M2;

void Init()
{
    M2 = 1;
    for (int i=0; i<LEN-2; i++)    // M2 = AA_NUMBER ^ (LEN-2);
        M2 *= AA_NUMBER;
    M1 = M2 * AA_NUMBER;        // M1 = AA_NUMBER ^ (LEN-1);
    M  = M1 *AA_NUMBER;            // M  = AA_NUMBER ^ (LEN);
}

void InitVectors(long* vector, long* second, long one_l[], long indexs, long total, long total_l, long complement)
{
    vector = new long [M];
    second = new long [M1];
    //        one_l = new long [AA_NUMBER];
    memset(vector, 0, M * sizeof(long));
    memset(second, 0, M1 * sizeof(long));
    memset(one_l, 0, AA_NUMBER * sizeof(long));
    total = 0;
    total_l = 0;
    complement = 0;
}
void init_buffer(char* buffer, long* second, long& complement, long& indexs, long one_l[], long& total_l)
{
    complement++;
    indexs = 0;
    for (int i=0; i<LEN-1; i++)
    {
        short enc = encode(buffer[i]);
        one_l[enc]++;
        total_l++;
        indexs = indexs * AA_NUMBER + enc;
    }
    second[indexs]++;
}
void cont_buffer(char ch, long* vector, long* second, long one_l[], long& indexs, long& total, long& total_l, long& complement)
{
    short enc = encode(ch);
    one_l[enc]++;
    total_l++;
    long index = indexs * AA_NUMBER + enc;
    vector[index]++;
    total++;
    indexs = (indexs % M2) * AA_NUMBER + enc;
    second[indexs]++;
}

class Bacteria_2
{
public:
    long count;
    double* tv;
    long *ti;
    
    Bacteria_2(char* filename)
    {
        // Variables
        long* vector;
        long* second;
        long one_l[AA_NUMBER];
        //    long* one_l;
        long indexs = 0;
        long total;
        long total_l;
        long complement;
        
        // InitVectors
        vector = new long [M];
        second = new long [M1];
        //        one_l = new long [AA_NUMBER];
        memset(vector, 0, M * sizeof(long));
        memset(second, 0, M1 * sizeof(long));
        memset(one_l, 0, AA_NUMBER * sizeof(long));
        total = 0;
        total_l = 0;
        complement = 0;
//        InitVectors(vector, second, one_l, indexs, total, total_l, complement);
        
        // Get the file from directory
        string directory = "/Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree/cvtree/data/";
        string filePath = directory.append(filename);
        
        FILE* bacteria_file = fopen(filePath.c_str(), "r");
        if(bacteria_file == NULL)
        {
            fprintf(stderr, "Error: failed to open file %s\n", filename);
            exit(1);
        }
        // Read the file
        char ch;
        while ((ch = fgetc(bacteria_file)) != EOF)
        {
            if (ch == '>')
            {
                while (fgetc(bacteria_file) != '\n'); // skip rest of line
                
                char buffer[LEN-1];
                fread(buffer, sizeof(char), LEN-1, bacteria_file);
//                init_buffer(buffer);
                init_buffer(buffer, second, complement, indexs, one_l, total_l);
            }
            else if (ch != '\n' && ch != '\t' && ch != '\r')
                cont_buffer(ch, vector, second, one_l, indexs, total, total_l, complement);
        }
        
        // Calculations
        // not entirely sure but looks like `stochastic`
        long total_plus_complement = total + complement;
        double total_div_2 = total * 0.5;
        int i_mod_aa_number = 0;
        int i_div_aa_number = 0;
        long i_mod_M1 = 0;
        long i_div_M1 = 0;
        
        double one_l_div_total[AA_NUMBER];
        for (int i=0; i<AA_NUMBER; i++)
            one_l_div_total[i] = (double)one_l[i] / total_l;
        
        double* second_div_total = new double[M1];
        for (int i=0; i<M1; i++)
            second_div_total[i] = (double)second[i] / total_plus_complement;
        
        count = 0;
        double* t = new double[M];
        
        for(long i=0; i<M; i++)
        {
            double p1 = second_div_total[i_div_aa_number];
            double p2 = one_l_div_total[i_mod_aa_number];
            double p3 = second_div_total[i_mod_M1];
            double p4 = one_l_div_total[i_div_M1];
            double stochastic =  (p1 * p2 + p3 * p4) * total_div_2;
            
            if (i_mod_aa_number == AA_NUMBER-1)
            {
                i_mod_aa_number = 0;
                i_div_aa_number++;
            }
            else
                i_mod_aa_number++;
            
            if (i_mod_M1 == M1-1)
            {
                i_mod_M1 = 0;
                i_div_M1++;
            }
            else
                i_mod_M1++;
            
            if (stochastic > EPSILON)
            {
                t[i] = (vector[i] - stochastic) / stochastic;
                count++;
            }
            else
                t[i] = 0;
        }
        
        // Reset values
        delete[] second_div_total;
        delete vector;
        delete second;
        //        delete one_l;
        
        
        tv = new double[count];
        ti = new long[count];
        
        int pos = 0;
        for (long i=0; i<M; i++)
        {
            if (t[i] != 0)
            {
                tv[pos] = t[i];
                ti[pos] = i;
                pos++;
            }
        }
        delete[] t;
        
        fclose (bacteria_file);
    }

};

double CompareBacteria_2(Bacteria_2* b1, Bacteria_2* b2)
{
    double correlation = 0;
    double vector_len1=0;
    double vector_len2=0;
    long p1 = 0;
    long p2 = 0;
    while (p1 < b1->count && p2 < b2->count)
    {
        long n1 = b1->ti[p1];
        long n2 = b2->ti[p2];
        if (n1 < n2)
        {
            double t1 = b1->tv[p1];
            vector_len1 += (t1 * t1);
            p1++;
        }
        else if (n2 < n1)
        {
            double t2 = b2->tv[p2];
            p2++;
            vector_len2 += (t2 * t2);
        }
        else
        {
            double t1 = b1->tv[p1++];
            double t2 = b2->tv[p2++];
            vector_len1 += (t1 * t1);
            vector_len2 += (t2 * t2);
            correlation += t1 * t2;
        }
    }
    while (p1 < b1->count)
    {
        double t1 = b1->tv[p1++];
        vector_len1 += (t1 * t1);
    }
    while (p2 < b2->count)
    {
        double t2 = b2->tv[p2++];
        vector_len2 += (t2 * t2);
    }

    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

class Bacteria
{
private:
    long* vector;
    long* second;
    long one_l[AA_NUMBER];
    //    long* one_l;
    long indexs;
    long total;
    long total_l;
    long complement;
    
    void InitVectors()
    {
        vector = new long [M];
        second = new long [M1];
        //        one_l = new long [AA_NUMBER];
        memset(vector, 0, M * sizeof(long));
        memset(second, 0, M1 * sizeof(long));
        memset(one_l, 0, AA_NUMBER * sizeof(long));
        total = 0;
        total_l = 0;
        complement = 0;
    }
    
    void init_buffer(char* buffer)
    {
        complement++;
        indexs = 0;
        for (int i=0; i<LEN-1; i++)
        {
            short enc = encode(buffer[i]);
            one_l[enc]++;
            total_l++;
            indexs = indexs * AA_NUMBER + enc;
        }
        second[indexs]++;
    }
    
    void cont_buffer(char ch)
    {
        short enc = encode(ch);
        one_l[enc]++;
        total_l++;
        long index = indexs * AA_NUMBER + enc;
        vector[index]++;
        total++;
        indexs = (indexs % M2) * AA_NUMBER + enc;
        second[indexs]++;
    }
    
public:
    long count;
    double* tv;
    long *ti;
    
    void Calculate(char* filename)
    {
        string directory = "/Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree/cvtree/data/";
        string filePath = directory.append(filename);
        
        FILE* bacteria_file = fopen(filePath.c_str(), "r");
        if(bacteria_file == NULL)
        {
            fprintf(stderr, "Error: failed to open file %s\n", filename);
            exit(1);
        }
        InitVectors();
        
        char ch;
        while ((ch = fgetc(bacteria_file)) != EOF)
        {
            if (ch == '>')
            {
                while (fgetc(bacteria_file) != '\n'); // skip rest of line
                
                char buffer[LEN-1];
                fread(buffer, sizeof(char), LEN-1, bacteria_file);
                init_buffer(buffer);
            }
            else if (ch != '\n' && ch != '\t' && ch != '\r')
                cont_buffer(ch);
        }
        
        // not entirely sure but looks like `stochastic`
        long total_plus_complement = total + complement;
        double total_div_2 = total * 0.5;
        int i_mod_aa_number = 0;
        int i_div_aa_number = 0;
        long i_mod_M1 = 0;
        long i_div_M1 = 0;
        
        double one_l_div_total[AA_NUMBER];
        for (int i=0; i<AA_NUMBER; i++)
            one_l_div_total[i] = (double)one_l[i] / total_l;
        
        double* second_div_total = new double[M1];
        for (int i=0; i<M1; i++)
            second_div_total[i] = (double)second[i] / total_plus_complement;
        
        count = 0;
        double* t = new double[M];
        
        // Para-04: Parallelise in a straightforward way
        omp_set_num_threads(8);
        #pragma omp parallel
        {
            long local_i_mod_aa_number = i_mod_aa_number;
            long local_i_div_aa_number = i_div_aa_number;
            long local_i_mod_M1 = i_mod_M1;
            long local_i_div_M1 = i_div_M1;
            int local_count = 0;

            #pragma omp for
            for (long i = 0; i < M; i++) {
                double p1 = second_div_total[local_i_div_aa_number];
                double p2 = one_l_div_total[local_i_mod_aa_number];
                double p3 = second_div_total[local_i_mod_M1];
                double p4 = one_l_div_total[local_i_div_M1];
                double stochastic = (p1 * p2 + p3 * p4) * total_div_2;

                if (local_i_mod_aa_number == AA_NUMBER - 1) {
                    local_i_mod_aa_number = 0;
                    local_i_div_aa_number++;
                } else {
                    local_i_mod_aa_number++;
                }

                if (local_i_mod_M1 == M1 - 1) {
                    local_i_mod_M1 = 0;
                    local_i_div_M1++;
                } else {
                    local_i_mod_M1++;
                }

                if (stochastic > EPSILON) {
                    t[i] = (vector[i] - stochastic) / stochastic;
                    local_count++;
                } else {
                    t[i] = 0;
                }
            }

            // Aggregate local counts if needed
            #pragma omp atomic
            count += local_count;
        }
        
        // Reset values
        delete[] second_div_total;
        delete vector;
        delete second;
        //        delete one_l;
        
        
        tv = new double[count];
        ti = new long[count];
        
        int pos = 0;
        // Para-05
        #pragma omp parallel for num_threads(8)
        for (long i=0; i<M; i++)
        {
            if (t[i] != 0)
            {
                tv[pos] = t[i];
                ti[pos] = i;
                pos++;
            }
        }
        delete[] t;
        
        fclose (bacteria_file);
    }
    
    Bacteria(char* filename)
    {
        Calculate(filename);
    }
};

void ReadInputFile(const char* input_name)
{
    FILE* input_file = fopen(input_name, "r");
    if(input_file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
        exit(1);
    }

    
    fscanf(input_file, "%d", &number_bacteria);
    bacteria_name = new char*[number_bacteria];

    
    for(long i=0;i<number_bacteria;i++)
    {
        bacteria_name[i] = new char[20];
        fscanf(input_file, "%s", bacteria_name[i]);
        strcat(bacteria_name[i], ".faa");
    }
    fclose(input_file);
}

double CompareBacteria(Bacteria* b1, Bacteria* b2)
{
    double correlation = 0;
    double vector_len1=0;
    double vector_len2=0;
    long p1 = 0;
    long p2 = 0;
    while (p1 < b1->count && p2 < b2->count)
    {
        long n1 = b1->ti[p1];
        long n2 = b2->ti[p2];
        if (n1 < n2)
        {
            double t1 = b1->tv[p1];
            vector_len1 += (t1 * t1);
            p1++;
        }
        else if (n2 < n1)
        {
            double t2 = b2->tv[p2];
            p2++;
            vector_len2 += (t2 * t2);
        }
        else
        {
            double t1 = b1->tv[p1++];
            double t2 = b2->tv[p2++];
            vector_len1 += (t1 * t1);
            vector_len2 += (t2 * t2);
            correlation += t1 * t2;
        }
    }
    while (p1 < b1->count)
    {
        double t1 = b1->tv[p1++];
        vector_len1 += (t1 * t1);
    }
    while (p2 < b2->count)
    {
        double t2 = b2->tv[p2++];
        vector_len2 += (t2 * t2);
    }

    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

// Doesn't work
double CompareBacteria_parallel(Bacteria* b1, Bacteria* b2)
{
    double correlation = 0;
    double vector_len1=0;
    double vector_len2=0;
    long p1 = 0;
    long p2 = 0;
    
    int num_threads = 4;
    
    // Calculate the size of each part
    long b1_part_size = b1->count / num_threads;
    int b1_remainder = b1->count % num_threads;
    
    long b2_part_size = b2->count / num_threads;
    int b2_remainder = b2->count % num_threads;

    // Create a vector to store the parts
    vector<long> b1_count_thread(num_threads, b1_part_size);
    vector<long> b2_count_thread(num_threads, b2_part_size);

    // Distribute the remainder equally among the first 'remainder' parts
    for (int i = 0; i < b1_remainder; i++) {
        b1_count_thread[i]++;
    }
    for (int i = 0; i < b2_remainder; i++) {
        b2_count_thread[i]++;
    }
            
    // Separating variables into `num_thread` partitions
    long p1_thread[num_threads];
    long p2_thread[num_threads];
    double vector_len1_thread[num_threads];
    double vector_len2_thread[num_threads];
    double correlation_thread[num_threads];
    
    omp_set_num_threads(num_threads);
    // Run the entire section in parallel for `num_threads` threads
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        // Initialise values for this thread
        p1_thread[thread_id] = 0; p2_thread[thread_id] = 0;
        vector_len1_thread[thread_id] = 0; vector_len2_thread[thread_id] = 0;
        correlation_thread[thread_id] = 0;
        while(p1_thread[thread_id] < b1_count_thread[thread_id] && p2_thread[thread_id] < b2_count_thread[thread_id])
        {
            long n1 = b1->ti[p1_thread[thread_id]];
            long n2 = b2->ti[p2_thread[thread_id]];
            
            if (n1 < n2)
            {
                double t1 = b1->tv[p1_thread[thread_id]];
                #pragma omp critical
                {
                    vector_len1_thread[thread_id] += (t1 * t1);
                    p1_thread[thread_id]++;
                }
            }
            else if (n2 < n1)
            {
                double t2 = b2->tv[p2_thread[thread_id]];
                #pragma omp critical
                {
                    vector_len2_thread[thread_id] += (t2 * t2);
                    p2_thread[thread_id]++;
                }
            }
            else
            {
                #pragma omp critical
                {
                    double t1 = b1->tv[p1_thread[thread_id]++];
                    double t2 = b2->tv[p2_thread[thread_id]++];
                    vector_len1_thread[thread_id] += (t1 * t1);
                    vector_len2_thread[thread_id] += (t2 * t2);
                    correlation_thread[thread_id] += t1 * t2;
                }
            }
        }
        while(p1_thread[thread_id] < b1_count_thread[thread_id])
        {
            #pragma omp critical
            {
                double t1 = b1->tv[p1_thread[thread_id]++];
                vector_len1_thread[thread_id] += (t1 * t1);
            }
        }
        while(p2_thread[thread_id] < b2_count_thread[thread_id])
        {
            #pragma omp critical
            {
                double t2 = b2->tv[p2_thread[thread_id]++];
                vector_len2_thread[thread_id] += (t2 * t2);
            }
        }
    }
    
    for(int i =0; i<num_threads; i++)
    {
        correlation += correlation_thread[i];
        vector_len1 += vector_len1_thread[i];
        vector_len2 += vector_len2_thread[i];
    }

    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

double CompareBacteria_parallel_2(int num_threads, Bacteria* b1, Bacteria* b2)
{
    double correlation = 0;
    double vector_len1 = 0;
    double vector_len2 = 0;
    long p1 = 0;
    long p2 = 0;

    #pragma omp parallel num_threads(num_threads) shared(correlation, vector_len1, vector_len2, p1, p2)
    {
        #pragma omp for reduction(+:correlation) reduction(+:vector_len1) reduction(+:vector_len2)
        for (int i = 0; i < num_threads; i++) {
            long local_p1 = p1;
            long local_p2 = p2;
            double local_vector_len1 = 0;
            double local_vector_len2 = 0;
            double local_correlation = 0;

            while (local_p1 < b1->count && local_p2 < b2->count)
            {
                long n1 = b1->ti[local_p1];
                long n2 = b2->ti[local_p2];

                if (n1 < n2)
                {
                    double t1 = b1->tv[local_p1];
                    local_vector_len1 += (t1 * t1);
                    local_p1++;
                }
                else if (n2 < n1)
                {
                    double t2 = b2->tv[local_p2];
                    local_p2++;
                    local_vector_len2 += (t2 * t2);
                }
                else
                {
                    double t1 = b1->tv[local_p1++];
                    double t2 = b2->tv[local_p2++];
                    local_vector_len1 += (t1 * t1);
                    local_vector_len2 += (t2 * t2);
                    local_correlation += t1 * t2;
                }
            }

            while (local_p1 < b1->count)
            {
                double t1 = b1->tv[local_p1++];
                local_vector_len1 += (t1 * t1);
            }
            while (local_p2 < b2->count)
            {
                double t2 = b2->tv[local_p2++];
                local_vector_len2 += (t2 * t2);
            }

            // Update the shared variables with thread-local results
            #pragma omp critical
            {
                correlation += local_correlation;
                vector_len1 += local_vector_len1;
                vector_len2 += local_vector_len2;
                p1 = local_p1;
                p2 = local_p2;
            }
        }
    }

    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

// Doesn't work
double CompareBacteria_parallel_3(int num_threads, Bacteria* b1, Bacteria* b2)
{
    double correlation = 0;
    double vector_len1 = 0;
    double vector_len2 = 0;
    long p1 = 0;
    long p2 = 0;
    
    // Shared variables for reduction
    double local_correlation = 0;
    double local_vector_len1 = 0;
    double local_vector_len2 = 0;

    #pragma omp parallel num_threads(num_threads) private(local_correlation, local_vector_len1, local_vector_len2) shared(correlation, vector_len1, vector_len2, p1, p2)
    {
        #pragma omp for
        for (int i = 0; i < num_threads; i++) {
            long local_p1 = 0;
            long local_p2 = 0;

            double thread_local_correlation = 0;
            double thread_local_vector_len1 = 0;
            double thread_local_vector_len2 = 0;

            while (local_p1 < b1->count && local_p2 < b2->count)
            {
                long n1 = b1->ti[local_p1];
                long n2 = b2->ti[local_p2];

                if (n1 < n2)
                {
                    double t1 = b1->tv[local_p1];
                    thread_local_vector_len1 += (t1 * t1);
                    local_p1++;
                }
                else if (n2 < n1)
                {
                    double t2 = b2->tv[local_p2];
                    thread_local_vector_len2 += (t2 * t2);
                    local_p2++;
                }
                else
                {
                    double t1 = b1->tv[local_p1++];
                    double t2 = b2->tv[local_p2++];
                    thread_local_vector_len1 += (t1 * t1);
                    thread_local_vector_len2 += (t2 * t2);
                    thread_local_correlation += t1 * t2;
                }
            }

            // Accumulate thread-local results
            #pragma omp critical
            {
                local_correlation += thread_local_correlation;
                local_vector_len1 += thread_local_vector_len1;
                local_vector_len2 += thread_local_vector_len2;
                p1 += local_p1;
                p2 += local_p2;
            }
        }
    }

    // Combine the results from all threads
    #pragma omp master
    {
        correlation = local_correlation;
        vector_len1 = local_vector_len1;
        vector_len2 = local_vector_len2;
    }

    // Ensure all threads have completed before calculating the final result
//    #pragma omp barrier

    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}


void CompareAllBacteria()
{
    int num_threads = 1;
    // Para-02-version-1-start: A bit more work, partition loading into independent sections
    // each nested for loop is its own, independent
//    int partitions = 3; // totalThreads - 1: Can be maximum of 7
//    vector<thread> threads;
//    
//    Bacteria** b = new Bacteria*[number_bacteria];
//    int start = 0;
//    int end = 0;
//    for(int i=0; i<partitions; ++i)
//    {
//        start = i*(number_bacteria/partitions);
//        end = (i+1) * (number_bacteria/partitions);
//        
//        threads.emplace_back([start,end, b]()
//        {
//            for(int j = start; j <end; j++)
//            {
//                printf("load %d of %d\n", j+1, number_bacteria);
//                b[j] = new Bacteria(bacteria_name[j]);
//            }
//        });
//    }
//    // Wait for threads to finish: prevents errors, bugs and 'kills' threads
//    // Failing to join the threads can lead to issues like accessing resources that are being modified by the
//    //      threads or not properly releasing thread-related resources, causing errors and unexpected behavior.
//    // When you join the threads, you are essentially waiting for each thread to finish its work before continuing
//    //      with the rest of the program. This ensures that all threads have completed their tasks before you
//    //      proceed with any further actions that depend on the results of those threads.
//    for (auto& thread : threads)
//    {
//        thread.join();
//    }
//    
//    // Do left over sequentially
//    if(end!=number_bacteria)
//    {
//        for(int k=end; k<number_bacteria;k++)
//        {
//            printf("leftover load %d of %d\n", k+1, number_bacteria);
//            b[k] = new Bacteria(bacteria_name[k]);
//        }
//    }
    // Para-02-version-1-end
        
    // num_threads = 01: 32 seconds
    // num_threads = 02: 21 seconds
    // num_threads = 03: 16 seconds
    // num_threads = 04: 15 seconds
    // num_threads = 05: 14 seconds
    // num_threads = 06: 13 seconds
    // num_threads = 07: 14 seconds
    // num_threads = 08: 14 seconds
    // num_threads = 09: 13 seconds
    // num_threads = 10: 15 seconds
    // num_threads = 11: 15 seconds
    // num_threads = 12: 15 seconds
    // num_threads = 13: 15 seconds
    // num_threads = 14: 15 seconds
    // num_threads = 15: 15 seconds
    // num_threads = 16: 15 seconds
        
    Bacteria** b = new Bacteria*[number_bacteria];
//    Bacteria_2** b_2 = new Bacteria_2*[number_bacteria];

//    Bacteria Timing
//    Sequential
//    Loading Bacteria: 20 seconds
//    CompareBacteria: 14 seconds
//    Printing Compare: 0 seconds
//    time elapsed: 34 seconds
//
//    Parallel: 8 Threads
//    Loading Bacteria: 9 seconds
//    CompareBacteria: 4 seconds
//    Printing Compare: 0 seconds
//    time elapsed: 13 seconds
    
//    Bacteria_2 Timing
//    Sequential
//    Loading Bacteria: 21 seconds
//    CompareBacteria: 13 seconds
//    Printing Compare: 0 seconds
//    time elapsed: 34 seconds
//
//    Parallel: 8 Threads
//    Loading Bacteria: 8 seconds
//    CompareBacteria: 5 seconds
//    Printing Compare: 0 seconds
//    time elapsed: 13 seconds
    
    time_t t1 = time(NULL);
    // Para-02-version-2: Same speedup compared to version-1
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i<number_bacteria; i++)
    {
        printf("load %d of %d RUNNING ON THREAD: %d \n", i + 1, number_bacteria, omp_get_thread_num()+1);
        b[i] = new Bacteria(bacteria_name[i]);
//        b_2[i] = new Bacteria_2(bacteria_name[i]);
    }
    time_t t2 = time(NULL);
    printf("Loading Bacteria: %ld seconds\n", t2 - t1);
    
    // Para-01: Straightforward parallelisation
    time_t t3 = time(NULL);
    vector<tuple<int, int, double>> correlations;
//    #pragma omp parallel for num_threads(num_threads/2)
    for(int i=0; i<number_bacteria-1; i++)
    {
        // 13 seconds with Para-01 here & Para-02
        #pragma omp parallel for num_threads(2)
        for(int j=i+1; j<number_bacteria; j++)
        {
            double correlation = CompareBacteria(b[i], b[j]);
//            double correlation = CompareBacteria_2(b_2[i], b_2[j]);
//            double correlation = CompareBacteria_parallel_2(2, b[i], b[j]);
            correlations.push_back(make_tuple(i, j, correlation));
        }
    }
    time_t t4 = time(NULL);
    printf("CompareBacteria: %ld seconds\n", t4 - t3);
    
    time_t t5 = time(NULL);
    // Sort the correlations vector based on i and j before printing
    sort(correlations.begin(), correlations.end(),[](const tuple<int, int, double>& a, const std::tuple<int, int, double>& b) {
        int ai, aj, bi, bj;
        tie(ai, aj, std::ignore) = a;
        tie(bi, bj, std::ignore) = b;
        if (ai != bi) {
            return ai < bi;
        } else {
            return aj < bj;
        }
    });
    
    // Print the results sequentially
    for (const auto& tuple : correlations) {
        int i, j;
        double correlation;
        tie(i, j, correlation) = tuple;
        printf("%2d %2d -> %.20lf\n", i, j, correlation);
    }
    time_t t6 = time(NULL);
    printf("Printing Compare: %ld seconds\n", t6 - t5);
}

void Delay()
{
    // Delay to open profiler
    cout << "Start of delayed method" << endl;

    // Delay for 10 seconds
    this_thread::sleep_for(chrono::seconds(10));

    // Continue with the work after the delay
    cout << "End of delayed method" << endl;
}

int main(int argc,char * argv[])
{
//    Delay();
    time_t t1 = time(NULL);
    Init();
    ReadInputFile("/Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree/cvtree/list.txt");
    CompareAllBacteria();

    time_t t2 = time(NULL);
    printf("time elapsed: %ld seconds\n", t2 - t1);
    return 0;
}

