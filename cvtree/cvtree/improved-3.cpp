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

void CompareAllBacteria()
{
    // Para-02-version-1-start: A bit more work, partition loading into independent sections
    // each nested for loop is its own, independent
    int partitions = 3; // totalThreads - 1: Can be maximum of 7
    vector<thread> threads;
    
    int num_threads = 6;
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
    
    // Para-02-version-2: Same speedup compared to version-1
    Bacteria** b = new Bacteria*[number_bacteria];
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i<number_bacteria; i++)
    {
        printf("load %d of %d RUNNING ON THREAD: %d \n", i + 1, number_bacteria, omp_get_thread_num());
        b[i] = new Bacteria(bacteria_name[i]);
    }
    
    
    // Para-01: Straightforward parallelisation
    #pragma omp parallel for num_threads(num_threads)
    for(int i=0; i<number_bacteria-1; i++)
    {
        // 13 seconds with Para-01 here & Para-02
        for(int j=i+1; j<number_bacteria; j++)
        {
            double correlation = CompareBacteria(b[i], b[j]);
            printf("%2d %2d -> %.20lf\n", i, j, correlation);
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
    printf("time elapsed: %ld seconds\n", t2 - t1);
    return 0;
}
