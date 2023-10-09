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
//using namespace std;
//
//int number_bacteria;
//char** bacteria_name;
//long M, M1, M2;
//short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
//#define encode(ch)		code[ch-'A']
//#define LEN				6
//#define AA_NUMBER		20 // how many amino acids we can have
//#define	EPSILON			1e-010
//
//void Init()
//{
//    // Think of base 20 instead of 2 -- like in binary
//	M2 = 1;
//    // The number of possible (k-2)-mers     4-mers
//	for (int i=0; i<LEN-2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
//		M2 *= AA_NUMBER;
//    // The number of possible (k-1)-mers     5-mers
//	M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
//    // The number of possible k-mers         6-mers
//	M  = M1 *AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
//}
//
//class Bacteria
//{
//private:
//	long* second;
//	long one_l[AA_NUMBER];
//	long indexs;
//	long total;
//	long total_l;
//	long complement;
//
//    // Initialises the count values for every single kmer to 0
//	void InitVectors()
//	{
//		vector = new long [M];  // each element is a k-mer
//		second = new long [M1]; // each element is a (k-1)-mer
//		memset(vector, 0, M * sizeof(long));    // sets the vector to 0
//		memset(second, 0, M1 * sizeof(long));
//		memset(one_l, 0, AA_NUMBER * sizeof(long)); // Each 1-mer. eg. number of times we see acid A, etc
//		total = 0;
//		total_l = 0;
//		complement = 0;
//	}
//
//    // creates a buffer of length 5
//    // initialises the indexs that we can modify as we carry on
//	void init_buffer(char* buffer)
//	{
//		complement++;
//		indexs = 0;
//		for (int i=0; i<LEN-1; i++)
//		{
//			short enc = encode(buffer[i]);
//			one_l[enc]++;
//			total_l++;
//			indexs = indexs * AA_NUMBER + enc;
//		}
//		second[indexs]++; // TODO: NOTE: See how this only appends for (k-1)-mer and not the k-mer
//	}
//
//    /// Continues incrementing the number of k, k-1 and 1-mer by doing a 'trick' with the index since numbering is base 2-
//    /// - Parameter ch: a character that is the last char in a k-mer
//	void cont_buffer(char ch)
//	{
//		short enc = encode(ch); // enc is encoding a number of the alphabet between 1 and 19 (base 20 encoding)
//		one_l[enc]++;
//		total_l++;
//        // This is a trick, we multiply our previous index + the enc for our kth character
//		long index = indexs * AA_NUMBER + enc;
//		vector[index]++; // TODO: NOTE: increment the vector to the next k-mer
//		total++;
//        // we can use the index for our k-mer to create a reduced index for the 5-mer
//		indexs = (indexs % M2) * AA_NUMBER + enc; // MOD to base of M2 to reduce calculation
//		second[indexs]++;
//	}
//
//public:
//	long* vector;
//
//	Bacteria(char* filename)
//	{
//        // Add folder directory
//        string directory = "data/";
//        string filePath = directory.append(filename);
//
//        // Get the file
//        FILE* bacteria_file = fopen(filePath.c_str(), "r");
//        if(bacteria_file == NULL)
//        {
//            fprintf(stderr, "Error: failed to open file %s\n", filename);
//            exit(1);
//        }
//		InitVectors(); // Vectors of 0's
//
//        // Parse the bacteria_file char-by-char
//		char ch;
//		while ((ch = fgetc(bacteria_file)) != EOF)
//		{
//			if (ch == '>') // if we're reading a new gene
//			{
//				while (fgetc(bacteria_file) != '\n'); // continue until
//
//
//                // TODO: NOTE: This effectivey only reads the first 5 characters of the gene, we do a trick in cont_buffer that will let us skip k-mers each time
//				char buffer[LEN-1];
//				fread(buffer, sizeof(char), LEN-1, bacteria_file);
//                // TODO: NOTE: A lot of work happens here
//				init_buffer(buffer); // First 5 character of each gene
//			}
//			else if (ch != '\n') // not at the end of the file
//
//				cont_buffer(ch); // Every other 5-char in the gene
//		}
//		fclose (bacteria_file);
//	}
//
//    // Resets a bacteria for future calculations
//	~Bacteria()
//	{
//		delete vector;
//		delete second;
//	}
//
//    // TODO: NOTE: returns the amount of times we expect to see kmer i
//    // this is based on (k-1)-mer and 1-mer
//    // TODO: NOTE: Example
//    // line 740 `AcMNPV.faa`
//    // MNRFFR = 15/ 10,000 (for example):   10,000 is the amount of k-mers
//    // MNRFF occured at rate p1 and R occurred at rate p2
//    // so we can calculate the probability of R coming after MNRFF
//    // We can also do get p of M and p of NRFFR and calculate that
//	double stochastic_compute(long i)
//	{
//        // First way eg: MNRFF and R
//		double p1 = (double)second[i / AA_NUMBER] / (total + complement); // probability of (k-1)-mer
//		double p2 = (double) one_l[i % AA_NUMBER] / total_l;             // probability of 1-mer
//        // Second way eg: M and NRFFR
//        double p3 = (double)second[i % M1] / (total + complement);
//		double p4 = (double) one_l[i / M1] / total_l;
//		return total * (p1*p2 + p3*p4) / 2; // Average it
//	}
//};
//
//// Reads an input file that has the number of files as the first line
//// and file names (without .faa) in the rest of the file
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
//	bacteria_name = new char*[number_bacteria];
//
//    // Read every string one by one into our bacteria_name global variable
//	for(long i=0;i<number_bacteria;i++)
//	{
//        bacteria_name[i] = new char[20];
//        fscanf(input_file, "%s", bacteria_name[i]);
//        strcat(bacteria_name[i], ".faa"); // concatenates the file name with .faa
//	}
//	fclose(input_file);
//}
//
//// NOTE: Can parallelise here
//double CompareBacteria(Bacteria* b1, Bacteria* b2)
//{
//	double correlation = 0;
//	double vector_len1=0;
//	double vector_len2=0;
//
//	for(long i=0; i<M; i++) // for every k-mer
//	{
//        // NOTE: How many times we expect to see a kmer in bacteria 1
//		double stochastic1 = b1->stochastic_compute(i);
//		double t1;
//		if (stochastic1 > EPSILON)
//            // if t1 is +ve, we see the k-mer more than we expect
//            // if t2 is -ve, we see the k-mer less than we expected
//			t1 = (b1->vector[i] - stochastic1) / stochastic1;
//		else
//			t1=0;
//		vector_len1 += (t1 * t1);
//
//		double stochastic2 = b2->stochastic_compute(i);
//        double t2;
//		if (stochastic2 > EPSILON)
//			t2 = (b2->vector[i] - stochastic2) / stochastic2;
//		else
//			t2 = 0;
//		vector_len2 += (t2 * t2);
//
//		correlation = correlation + t1 * t2; // // for every k-mer we add the corr score (could be 0 or -ve)
//	}
//
//	return correlation / (sqrt(vector_len1) * sqrt(vector_len2)); // normalise
//}
//
//// TODO: Profile this method -- very similar to matrix multiplication where you can simply have a thread for each comparison. Though you don't want to compare a bacteria to one that's already been compared to so be careful
//void CompareAllBacteria()
//{
//    for(int i=0; i<number_bacteria-1; i++) // For each bacteria
//	{
//		Bacteria* b1 = new Bacteria(bacteria_name[i]); // create a bacteria (from a file)
//
//		for(int j=i+1; j<number_bacteria; j++) // Compare against kmers of all others
//		{
//			Bacteria* b2 = new Bacteria(bacteria_name[j]);
//			double correlation = CompareBacteria(b1, b2); // TODO: Parallelisation can be done here
//			printf("%03d %03d -> %.10lf\n", i, j, correlation);
//			delete b2;
//		}
//		delete b1;
//	}
////    original_data.close();
//}
//
//int main(int argc,char * argv[])
//{
//	time_t t1 = time(NULL);
//
//	Init(); // Intialises constant variables
//    ReadInputFile("list.txt");
//    // Important: step 1 is to profile this
//	CompareAllBacteria();
//
//	time_t t2 = time(NULL);
//	printf("time elapsed: %ld seconds\n", t2 - t1);
//	return 0;
//}
