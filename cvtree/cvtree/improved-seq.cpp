//
//  improved-seq.cpp
//  improved-seq
//
//  Created by Shridhar Thorat on 11/10/2023.
//  ISXX - standard for XXth Improvement to the Sequential

#include "cvtree_parallel.hpp"

using namespace std;

int number_bacteria;
char** bacteria_name;
long M, M1, M2;

void ReadInputFile(const char* input_name) // was const originally but not in the video
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
        strcat(bacteria_name[i], ".faa"); // concatenates the file name with .faa
    }
    fclose(input_file);
}

void Init()
{
    M2 = 1;
    for (int i=0; i<LEN-2; i++)
        M2 *= AA_NUMBER;
    M1 = M2 * AA_NUMBER;
    M  = M1 *AA_NUMBER;
}

class Bacteria
{
private:
    long* vector;
    long* second;
    long one_l[AA_NUMBER];
    long indexs;
    long total;
    long total_l;
    long complement;
    
    void InitVectors()
    {
        vector = new long [M];
        second = new long [M1];
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
        total_l += LEN-1; // IS01: O(1) rather than O(LEN-2)
        for (int i=0; i<LEN-1; i++)
        {
            short enc = encode(buffer[i]);
            one_l[enc]++;
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

    Bacteria(char* filename)
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
                while (fgetc(bacteria_file) != '\n');

                char buffer[LEN-1];
                fread(buffer, sizeof(char), LEN-1, bacteria_file);
                init_buffer(buffer);
            }
            else if (ch != '\n' && ch != '\t' && ch != '\r')
                cont_buffer(ch);
        }
        
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
            second_div_total[i] = (double)second[i] / (total_plus_complement);

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

/// Calculates the starting index for bacteria j in an array that stores values for every comparison
/// IS02
int start(int j)
{
    int start = 0;
    for(int l=0; l<j;l++)
    {
        start += number_bacteria-1-l;
    }
    return start-1;
}


void CompareAllBacteria()
{
    Bacteria** b = new Bacteria*[number_bacteria];
    int* comparisons = new int[820]; // track comparisons, instead of just writing 1, just write the correlation value since complexity doesn't increase
    int comparison_i = 0;

    for(int i=0; i<number_bacteria; i++)
    {
        printf("load %d of %d\n", i, number_bacteria-1);
        b[i] = new Bacteria(bacteria_name[i]);
        
        // IS02: Reduced nested for loop to single loop
        for(int j=0; j<i; j++)
        {
            comparison_i = start(j) + i - j;
            printf("%2d %2d -> ", j, i);
            double correlation = CompareBacteria(b[i], b[j]);
            printf("%2d %2d -> %.20lf\n", j, i, correlation);
            comparisons[comparison_i] = correlation;
        }
    }
    
}

int main(int argc,char * argv[])
{
    time_t t1 = time(NULL);

    Init();
    ReadInputFile("/Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree/cvtree/list.txt");
    CompareAllBacteria();

    time_t t2 = time(NULL);
    printf("time elapsed: %ld seconds\n", (t2) - (t1));
    return 0;
}

