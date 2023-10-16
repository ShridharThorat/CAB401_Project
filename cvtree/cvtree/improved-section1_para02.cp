//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <time.h>
//#include <math.h>
//#include <string>
//#include <cstring>
//#include <chrono>
//#include <iostream>
//#include <fstream>
//#include <os/log.h>
//#include <os/signpost.h>
//#include <vector>
//#include <numeric>
//using namespace std;
//
//int number_bacteria;
//char** bacteria_name;
//long M, M1, M2;
//short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
//#define encode(ch) code[ch-'A']
//#define LEN 6
//#define AA_NUMBER 20
//#define EPSILON 1e-010
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
//void InitVectors_threads(int num_threads)
//{
//    
//}
//
//class Bacteria
//{
//private:
//    long* vector;
//    long* second;
//    long one_l[AA_NUMBER];
//    long indexs;
//    long total;
//    long total_l;
//    long complement;
//
//    void InitVectors()
//    {
//        vector = new long [M];
//        second = new long [M1];
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
//        short enc = encode(ch);
////        cout << "ch '" << ch <<"': enc: "<< enc << endl;
//        one_l[enc]++;
//        total_l++;
//        long index = indexs * AA_NUMBER + enc;
////        cout << "index: " << index << endl;
//        vector[index]++;
//        total++;
//        indexs = (indexs % M2) * AA_NUMBER + enc;
////        cout << "indexs: " << indexs << endl;
//        second[indexs]++;
//    }
//
//public:
//    long count;
//    double* tv;
//    long *ti;
//
//    Bacteria(char* filename)
//    {
//        // Add folder directory
//        string directory = "/Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree/cvtree/data/";
//        string filePath = directory.append(filename);
//
//        // Get the file
//        FILE* bacteria_file = fopen(filePath.c_str(), "r");
//        if(bacteria_file == NULL)
//        {
//            fprintf(stderr, "Error: failed to open file %s\n", filename);
//            exit(1);
//        }
//        // Initialise Vectors to 0
//        InitVectors();
//
//        // Calculate k-mer, (k-1)-mer and 1-mer for entire bacteria
//        // TODO: the while, else-if can be parallelised but required 1) refactoring the sequential to something that can be better parallelised. Then 2) Create N threads - for each other character and do con_buffer simultaneously
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
//            double p1 = second_div_total[i_div_aa_number];
//            double p2 = one_l_div_total[i_mod_aa_number];
//            double p3 = second_div_total[i_mod_M1];
//            double p4 = one_l_div_total[i_div_M1];
//            double stochastic =  (p1 * p2 + p3 * p4) * total_div_2;
//
////            i_mod_aa_number = (i_mod_aa_number + 1) % AA_NUMBER;
////            if (i_mod_aa_number == 0) {
////                i_div_aa_number++;
////            }
////
////            i_mod_M1 = (i_mod_M1 + 1) % M1;
////            if (i_mod_M1 == 0) {
////                i_div_M1++;
////            }
////
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
//                t[i] = (vector[i] - stochastic) / stochastic;
//                count++;
//            }
//            else
//                t[i] = 0;
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
//        for (long i=0; i<M; i++)
//        {
//            if (t[i] != 0)
//            {
//                tv[pos] = t[i];
//                ti[pos] = i;
//                pos++;
//            }
//        }
//        delete[] t;
//
//        fclose (bacteria_file);
//    }
//};
//
//void ReadInputFile(const char* input_name) // was const originally but not in the video
//{
//    // Open the file
//    FILE* input_file = fopen(input_name, "r");
//    if(input_file == NULL)
//    {
//        fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
//        exit(1);
//    }
//
//    // Read the first line of the input file this is a number telling how many files we have
//    fscanf(input_file, "%d", &number_bacteria);
//    bacteria_name = new char*[number_bacteria];
//
//    // Read every string one by one into our bacteria_name global variable
//    for(long i=0;i<number_bacteria;i++)
//    {
//        bacteria_name[i] = new char[20];
//        fscanf(input_file, "%s", bacteria_name[i]);
//        strcat(bacteria_name[i], ".faa"); // concatenates the file name with .faa
//    }
//    fclose(input_file);
//}
//
//double CompareBacteria(Bacteria* b1, Bacteria* b2)
//{
//    double correlation = 0;
//    double vector_len1=0;
//    double vector_len2=0;
//    long p1 = 0;
//    long p2 = 0;
//    while (p1 < b1->count && p2 < b2->count)
//    {
//        long n1 = b1->ti[p1];
//        long n2 = b2->ti[p2];
//        if (n1 < n2)
//        {
//            double t1 = b1->tv[p1];
//            vector_len1 += (t1 * t1);
//            p1++;
//        }
//        else if (n2 < n1)
//        {
//            double t2 = b2->tv[p2];
//            p2++;
//            vector_len2 += (t2 * t2);
//        }
//        else
//        {
//            double t1 = b1->tv[p1++];
//            double t2 = b2->tv[p2++];
//            vector_len1 += (t1 * t1);
//            vector_len2 += (t2 * t2);
//            correlation += t1 * t2;
//        }
//    }
//    while (p1 < b1->count)
//    {
////        long n1 = b1->ti[p1];
//        double t1 = b1->tv[p1++];
//        vector_len1 += (t1 * t1);
//    }
//    while (p2 < b2->count)
//    {
////        long n2 = b2->ti[p2];
//        double t2 = b2->tv[p2++];
//        vector_len2 += (t2 * t2);
//    }
//
//    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
//}
//
//void CompareAllBacteria()
//{
//    int num_threads = 8;
//    //    os_log_t loadLog = os_log_create("cvtree", "loadingLogger");
//    //    os_log_t compareLog = os_log_create("cvtree", "compareLogger");
//    
//    Bacteria** b = new Bacteria*[number_bacteria];
//    auto t1 = std::chrono::high_resolution_clock::now();
//    //    os_signpost_interval_begin(loadLog,1, "Load Bacteria");
//    #pragma omp parallel for num_threads(num_threads)
//    for(int i=0; i<number_bacteria; i++)
//    {
//        printf("load %d of %d\n", i+1, number_bacteria);
//        b[i] = new Bacteria(bacteria_name[i]);
//    }
//    //    os_signpost_interval_end(loadLog,1, "Load Bacteria");
//    auto t2 = std::chrono::high_resolution_clock::now();
//    auto load = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1);
//    printf("Loading Bacteria: %f seconds\n", load.count()/1000.0);
////    v.push_back(load.count()/1000.0);
//    
//    auto t3 = std::chrono::high_resolution_clock::now();
//    //    os_signpost_interval_begin(compareLog,2, "Compare Bacteria");
//    for(int i=0; i<number_bacteria-1; i++)
//    {
//        for(int j=i+1; j<number_bacteria; j++)
//        {
//            printf("%2d %2d -> ", i, j);
//            double correlation = CompareBacteria(b[i], b[j]);
//            printf("%.20lf\n", correlation);
//        }
//    }
//    //    os_signpost_interval_end(compareLog,2, "Compare Bacteria");
//    auto t4 = std::chrono::high_resolution_clock::now();
//    auto compare =std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
//    printf("CompareBacteria: %f seconds\n", compare.count()/1000.0);
//
//    printf("\n");
//    
////    printf("Average load %f", average(loading_times));
////    printf("Average total %f", average(total_times));
//}
//
//int main(int argc,char * argv[])
//{
////    int trials = 2;
////    std::vector<float> total_times;
////    std::vector<float> loading_times;
////    for(int i=0;i<trials;i++)
////    {
////        printf("\nTrial %d\n",i);
//    auto t1 = std::chrono::high_resolution_clock::now();
//    
//    Init();
//    ReadInputFile("/Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree/cvtree/list.txt");
//    CompareAllBacteria();
////    CompareAllBacteria(loading_times);
//    
//    auto t2 = std::chrono::high_resolution_clock::now();
//    auto total =std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
//    printf("time elapsed: %f seconds\n", total.count()/1000.0);
////        total_times.push_back(total.count()/1000.0);
////    }
////    float avg_total = reduce(total_times.begin(), total_times.end()) / total_times.size();
////    float avg_load = reduce(loading_times.begin(), loading_times.end()) / loading_times.size();
////    cout << "avg total: " << avg_total << endl;
////    cout << "avg load: " << avg_load << endl;
//    return 0;
//}
