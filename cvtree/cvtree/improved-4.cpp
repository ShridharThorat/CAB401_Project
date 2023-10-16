////
////  improved-4.cpp
////  cvtree
////
////  Created by Shridhar Thorat on 3/10/2023.
////
//
//#include "cvtree_parallel.hpp"
//using namespace std;
//using namespace std::chrono;
//
//int number_bacteria;
//char** bacteria_name;
//long M, M1, M2;
//
//void Init()
//{
//    M2 = 1;
//    for (int i=0; i<LEN-2; i++)    // M2 = AA_NUMBER ^ (LEN-2);
//        M2 *= AA_NUMBER;
//    M1 = M2 * AA_NUMBER;        // M1 = AA_NUMBER ^ (LEN-1);
//    M  = M1 *AA_NUMBER;            // M  = AA_NUMBER ^ (LEN);
//}
//
//class Bacteria
//{
//private:
//    long* vector;
//    long* second;
//    long* one_l;
//    long indexs;
//    long total;
//    long total_l;
//    long complement;
//    
//    void InitVectors()
//    {
//        vector = new long [M];
//        second = new long [M1];
//        one_l = new long [AA_NUMBER];
//        memset(vector, 0, M * sizeof(long));
//        memset(second, 0, M1 * sizeof(long));
//        memset(one_l, 0, AA_NUMBER * sizeof(long));
//        total = 0;
//        total_l = 0;
//        complement = 0;
//    }
//    
//    void init_buffer(char* buffer)
//    {
//        complement++;
//        indexs = 0;
//        for (int i=0; i<LEN-1; i++)
//        {
////            cout << buffer[i];
//            short enc = encode(buffer[i]);
//            one_l[enc]++;
//            total_l++;
//            indexs = indexs * AA_NUMBER + enc;
//        }
//        second[indexs]++;
//    }
//    
//    void cont_buffer(char ch)
//    {
////        cout << ch;
//        short enc = encode(ch);
//        one_l[enc]++;
//        total_l++;
//        long index = indexs * AA_NUMBER + enc;
//        vector[index]++;
//        total++;
////        indexs = (indexs % M2) * AA_NUMBER + enc;
//        indexs = indexs % M2;
//        indexs = indexs*AA_NUMBER;
//        indexs = indexs + enc;
//        second[indexs]++;
//    }
//    
//public:
//    long count;
//    double* tv;
//    long *ti;
//    
//    Bacteria(int num_threads, char* filename)
//    {
//        string directory = "/Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree/cvtree/data/";
//        string filePath = directory.append(filename);
//        
//        FILE* bacteria_file = fopen(filePath.c_str(), "r");
//        if(bacteria_file == NULL)
//        {
//            fprintf(stderr, "Error: failed to open file %s\n", filename);
//            exit(1);
//        }
//        InitVectors();
//        
//        char ch;
//        while ((ch = fgetc(bacteria_file)) != EOF)
//        {
//            if (ch == '>')
//            {
//                while (fgetc(bacteria_file) != '\n'); // skip rest of line
//                
//                char buffer[LEN-1];
//                fread(buffer, sizeof(char), LEN-1, bacteria_file);
//                init_buffer(buffer);
//            }
//            else if (ch != '\n' && ch != '\t' && ch != '\r')
//                cont_buffer(ch);
//        }
//        
//        // not entirely sure but looks like `stochastic`
//        long total_plus_complement = total + complement;
//        double total_div_2 = total * 0.5;
//        int i_mod_aa_number = 0;
//        int i_div_aa_number = 0;
//        long i_mod_M1 = 0;
//        long i_div_M1 = 0;
//
//        double one_l_div_total[AA_NUMBER];
//        for (int i=0; i<AA_NUMBER; i++)
//            one_l_div_total[i] = (double)one_l[i] / total_l;
//
//        double* second_div_total = new double[M1];
//        for (int i=0; i<M1; i++)
//            second_div_total[i] = (double)second[i] / total_plus_complement;
//
//        count = 0;
//        double* t = new double[M];
//
//        for(long i=0; i<M; i++)
//        {
//            // stochastic_compute. See tidy.cpp
//            double p1 = second_div_total[i_div_aa_number];
//            double p2 = one_l_div_total[i_mod_aa_number];
//            double p3 = second_div_total[i_mod_M1];
//            double p4 = one_l_div_total[i_div_M1];
//            double stochastic =  (p1 * p2 + p3 * p4) * total_div_2;
//
//            if (i_mod_aa_number == AA_NUMBER-1)
//            {
//                i_mod_aa_number = 0;
//                i_div_aa_number++;
//            }
//            else
//                i_mod_aa_number++;
//
//            if (i_mod_M1 == M1-1)
//            {
//                i_mod_M1 = 0;
//                i_div_M1++;
//            }
//            else
//                i_mod_M1++;
//
//            if (stochastic > EPSILON)
//            {
//                // i is the index of the current k-mer (6-mer)
//                // if t[i] is +ve, we see the k-mer more than we expect
//                // if t[i] is -ve, we see the k-mer less than we expected
//                // TODO: NOTE: t represents the ratio between the amount of times we see it (vector[i]) and how many times we expect to see it (stochastic). We can call it the deviation.
//                // Could also say TODO: NOTE: normalized difference between the observed count and the expected count of a k-mer in the corresponding bacteria sequence.
//                t[i] = (vector[i] - stochastic) / stochastic; //
//                count++; // count represents the number of k-mers that
//            }
//            else
//                t[i] = 0; // we don't care when we see less than we expect
//        }
//        
//        // Reset values
//        delete[] second_div_total;
//        delete vector;
//        delete second;
//        
//        tv = new double[count];
//        ti = new long[count];
//        
//        int pos = 0;
//        // i is the k-mer, and pos is the index of the significant k-mers
//        for (long i=0; i<M; i++)
//        {
//            if (t[i] != 0) // i.e for significant deviations only
//            {
//                tv[pos] = t[i]; // TODO: NOTE: tv stores the deviation value for k-mer at index pos
//                ti[pos] = i; // TODO: NOTE: ti stores the actual k-mer at the index
//                pos++;
//            }
//        }
//        delete[] t;
//        
//        fclose (bacteria_file);
//    }
//};
//
//void ReadInputFile(const char* input_name)
//{
//    FILE* input_file = fopen(input_name, "r");
//    if(input_file == NULL)
//    {
//        fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
//        exit(1);
//    }
//    
//    fscanf(input_file, "%d", &number_bacteria);
//    bacteria_name = new char*[number_bacteria];
//
//    for(long i=0;i<number_bacteria;i++)
//    {
//        bacteria_name[i] = new char[20];
//        fscanf(input_file, "%s", bacteria_name[i]);
//        strcat(bacteria_name[i], ".faa");
//    }
//    fclose(input_file);
//}
//
//double CompareBacteria(Bacteria* b1, Bacteria* b2)
//{
////    auto start = high_resolution_clock::now();
//    double correlation = 0;
//    double vector_len1=0;
//    double vector_len2=0;
//    long p1 = 0;
//    long p2 = 0;
//
//    // Seq-01: Change while loops into for loops because for loops are better
////    for (long i = 0; p1 < b1->count && p2 < b2->count; i++)
//    while(p1 < b1->count && p2 < b2->count)
//    {
//        long n1 = b1->ti[p1];
//        long n2 = b2->ti[p2];
//        if (n1 < n2) // if the significant 6-mer in bacteria 2 is before the 6-mer in bacteria 1
//        {
//            double t1 = b1->tv[p1];
//            vector_len1 += (t1 * t1); // TODO: NOTE: squared lengths (norms) of vector for each 
//            p1++;
//        }
//        else if (n2 < n1)
//        {
//            double t2 = b2->tv[p2];
//            p2++;
//            vector_len2 += (t2 * t2);
//        }
//        else // If the 6-mer's in both b1 and b2 are the same
//        {
//            double t1 = b1->tv[p1++];
//            double t2 = b2->tv[p2++];
//            vector_len1 += (t1 * t1);
//            vector_len2 += (t2 * t2);
//            correlation += t1 * t2;
//        }
//    }
//
//    // Calculate vector_len1 and vector_len2 for the remaining elements
//    while (p1 < b1->count)
//    {
//        double t1 = b1->tv[p1++];
//        vector_len1 += (t1 * t1);
//    }
//    while (p2 < b2->count)
//    {
//        double t2 = b2->tv[p2++];
//        vector_len2 += (t2 * t2);
//    }
//    
//    
////    printf("<%f, %f> \n",vector_len1,vector_len2);
//    // TODO: NOTE: Normalising so all values are between -1 and 1.
//    double val = correlation / (sqrt(vector_len1) * sqrt(vector_len2));
////    auto stop = high_resolution_clock::now();
////
////    auto duration = duration_cast<microseconds>(stop - start);
////    cout << "Time taken by function: "<< duration.count() << " microseconds" << endl;
//    return val;
//}
//
//
//void CompareAllBacteria()
//{
//    os_log_t loadLog = os_log_create("cvtree", "loadingLogger");
//    os_log_t compareLog = os_log_create("cvtree", "compareLogger");
//        
//    int num_threads = 2;
//    Bacteria** b = new Bacteria*[number_bacteria];
//        
//    time_t t1 = time(NULL);
//    // Para-02-version-2: Same speedup compared to version-1
//    os_signpost_interval_begin(loadLog,1, "Load Bacteria");
//    #pragma omp parallel for num_threads(num_threads)
//    for (int i = 0; i<number_bacteria; i++)
//    {
//        printf("load %d of %d on thread: %d \n", i + 1, number_bacteria, omp_get_thread_num()+1);
//        b[i] = new Bacteria(num_threads, bacteria_name[i]);
//    }
//    os_signpost_interval_end(loadLog,1, "Load Bacteria");
//    time_t t2 = time(NULL);
//    printf("Loading Bacteria: %ld seconds\n", t2 - t1);
//    
//    
//    // Para-01: Straightforward parallelisation
//    time_t t3 = time(NULL);
//    vector<tuple<int, int, double>> correlations;
//    
//    os_signpost_interval_begin(compareLog,2, "Compare Bacteria");
//    for(int i=0; i<number_bacteria-1; i++)
//    {
//        #pragma omp parallel for num_threads(num_threads)
//        for(int j=i+1; j<number_bacteria; j++)
//        {
//            double correlation = CompareBacteria(b[i], b[j]);
//            
//            // Push the result into a thread-local vector
//            correlations.push_back(make_tuple(i, j, correlation));
//        }
//    }
//    os_signpost_interval_end(compareLog,2, "Compare Bacteria");
//    time_t t4 = time(NULL);
//    printf("CompareBacteria: %ld seconds\n", t4 - t3);
//    
//    time_t t5 = time(NULL);
//    // Sort the correlations vector based on i and j before printing
//    sort(correlations.begin(), correlations.end(),[](const tuple<int, int, double>& a, const std::tuple<int, int, double>& b) {
//        int ai, aj, bi, bj;
//        tie(ai, aj, std::ignore) = a;
//        tie(bi, bj, std::ignore) = b;
//        if (ai != bi) {
//            return ai < bi;
//        } else {
//            return aj < bj;
//        }
//    });
//    
//    // Print the results sequentially
//    for (const auto& tuple : correlations) {
//        int i, j;
//        double correlation;
//        tie(i, j, correlation) = tuple;
//        printf("%2d %2d -> %.20lf\n", i, j, correlation);
//    }
//    
//    time_t t6 = time(NULL);
//    printf("Printing Compare: %ld seconds\n", t6 - t5);
//}
//
//int main()
//{
//    Init();
//    ReadInputFile("list.txt");
//    CompareAllBacteria();
//    return 0;
//}
