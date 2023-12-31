// cd /Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree_project/cvtree_project
// chpl -o cvtree improved.chpl;
// clear; chpl -o cvtree improved.chpl; ./cvtree
use Math;
use Time;
use Random;
use IO;
use IO.FormattedIO;
use Path;

var number_bacteria: int;
var bacteria_domain = {0..-1};
var bacteria_name: [bacteria_domain] string;
var M, M1, M2: int;
const LEN: int = 6;
const AA_NUMBER: int = 20; // number of amino acids
const EPSILON: real= 1e-10; //0.0000000001;

var D: domain(1) = {0..25};
var amino_acids: [D] string = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'];
var amino_vals: [D] int = [0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3];

proc encode(amino_acid: string): int{
    // TODO: forall for data-thread, coforall for thread
    var amino_index: int;
    var found = amino_acids.find(amino_acid, amino_index);
    if found then return amino_vals[amino_index];
    else return -1;
}

proc Init()
{
    M2 = 1;
    for i in 1..(LEN-2) do {
        M2 *= AA_NUMBER; // M2 = AA_NUMBER ^ (LEN-2);
    }
    // The number of (k-1)-mers     5-mers
    M1 = M2* AA_NUMBER; // M1 = AA_NUMBER ^ (LEN-1);
    // The number of k-mers         6-mers
    M = M1 * AA_NUMBER; // M  = AA_NUMBER ^ (LEN);
}

class Bacteria
{
    // TODO: Note No need for InitVectors since these variables are pre initialised in Chapel
    var lastVectorIndex: int = M-1;
    var lastSecondIndex: int = M1-1;
    var lastOneIndex: int = AA_NUMBER-1;
    var vector: [0..lastVectorIndex] int;
    var second: [0..lastSecondIndex] int;
    var one_l: [0..lastOneIndex] int;
    var indexs, total, total_l, complement: int;
    
    var count: int;
    var tv_ti_dom = {0..-1};
    var tv: [tv_ti_dom] real;
    var ti: [tv_ti_dom] int;
    
    proc InitVectors()
    {
        this.vector = 0;
        this.second = 0;
        this.one_l = 0;
        this.indexs = 0;
        this.total = 0;
        this.total_l = 0;
        this.complement = 0;
    }
    
    proc init_buffer(buffer: string)
    {
        this.complement += 1;
        this.indexs = 0;
        for i in 0..(LEN-2) do
        {
            var enc = encode(buffer[i]);
            this.one_l[enc] += 1;
            // writeln("enc ", i, ": is ", enc); // debugging
            this.total_l += 1;
            this.indexs = indexs*AA_NUMBER + enc;
            // writeln("indexs ", i, ": is ", indexs); // debugging
        }
        this.second[indexs] += 1;
    }
    
    proc cont_buffer(ch: string)
    {
        var enc = encode(ch);
//        writeln("ch '", ch, "': enc: ", enc); // debugging
        this.one_l[enc] += 1;
        this.total_l += 1;
        var index_: int = indexs * AA_NUMBER + enc; // Note: original has 'index' but this is reserved in chapel
//        writeln("index: ", index_); // debugging
        this.vector[index_] += 1;
        this.total += 1;
        this.indexs = (indexs % M2) * AA_NUMBER + enc;
//        writeln("indexs: ", indexs); // debugging
        this.second[indexs] += 1;
    }
    
    proc init()
    {
        this.complete();
        this.InitVectors();
    }
    
    // Chapel's styling for constructors
    proc init(filename: string)
    {
        this.complete(); // Allows calling InitVectors later
        this.InitVectors();
        this.readBacteria(filename);
        
        // not entirely sure but looks like `stochastic`
        var total_plus_complement: int;
        total_plus_complement = this.total + this.complement;
        
        var total_div_2: real;
        total_div_2 = this.total * 0.5;
        
        var i_mod_aa_number: int;
        var i_div_aa_number: int;
        var i_mod_M1: int;
        var i_div_M1: int;
        
        var one_l_div_total: [0..this.lastOneIndex] real;
        for i in 0..this.lastOneIndex do
        {
            one_l_div_total[i] = this.one_l[i]:real / this.total_l;
        }
        
        var second_div_total: [0..this.lastSecondIndex] real;
        for i in 0..this.lastSecondIndex
        {
            second_div_total[i] = this.second[i]:real / total_plus_complement;
        }
        
        var t: [0..this.lastVectorIndex] real;
        
        for i in 0..this.lastVectorIndex do
        {
            var p1: real(64) = second_div_total[i_div_aa_number];
            var p2: real(64) = one_l_div_total[i_mod_aa_number];
            var p3: real(64) = second_div_total[i_mod_M1];
            var p4: real(64) = one_l_div_total[i_div_M1];
            
            var stochastic: real = (p1 * p2 + p3 * p4) * total_div_2;
            
            // For ones
            if(i_mod_aa_number == AA_NUMBER-1) then
            {
                i_mod_aa_number = 0;
                i_div_aa_number += 1;
            }
            else { i_mod_aa_number += 1; }
            
            // For (k-1)-mers
            if (i_mod_M1 == M1-1)
            {
                i_mod_M1 = 0;
                i_div_M1 += 1;
            }
            else { i_mod_M1 += 1; }
            
            if(stochastic > EPSILON)
            {
                t[i] = (this.vector[i] - stochastic) / stochastic;
                count += 1;
            }
            else { t[i] = 0; }
        }
        
        // Reset values TODO: Figure out how to do this in chapel
        
        tv_ti_dom = {0..(this.count-1)}; // implicitly resize ti and tv
        
        var pos: int = 0;
        for i in 0..lastVectorIndex do
        {
            if t[i] != 0 then
            {
                this.tv[pos] = t[i];
                this.ti[pos] = i;
                pos += 1;
            }
        }
//        delete t;

    }
    
    proc readBacteria(filename: string)
    {
        try{
            var directory: string = "data/";
            var filePath: string = directory + filename;
            var bacteria_file = open(filePath, ioMode.r);
            
            var bacteria_file_reader = bacteria_file.reader();
            var line: string;
            var ch: string;
//            writeln("Reading file: ", filename); // debugging
            while (bacteria_file_reader.readString(ch,1))
            {
                if ch == ">" then
                {
                    bacteria_file_reader.readString(ch,1);
                    while ( ch != '\n')
                    {
                        bacteria_file_reader.readString(ch,1);
                    }
                    var buffer_string: string;
                    bacteria_file_reader.readString(buffer_string, LEN-1); // (k-1)-mer long
                    // writeln("buffer_string: ", buffer_string); // debugging
                    init_buffer(buffer_string);
                    
                }
                else if ch != "\n" && ch != "\t" && ch != "\r"  then
                {
                    cont_buffer(ch);
                }
            }
            bacteria_file.close();
//            writeln("FINISHED: ", filename); // debugging
        } catch e: EofError {
            writeln("r is at EOF");
        } catch e: UnexpectedEofError {
            writeln("unable to read an 'int'");
        } catch e: SystemError {
            writeln("system error in IO implementation: ", e);
        } catch e: Error {
            writeln("something else went wrong...");
        }
        
    }
    
    proc get_count(): int { return this.count; }
//    proc get_ti(): [tv_ti_dom] real {return this.ti;}
    proc get_ti(i: int): int {return this.ti[i]; }
//    proc get_tv(): [tv_ti_dom] real {return this.tv;}
    proc get_tv(i: int): real {return this.tv[i]; }
}

proc ReadInputFile(input_name: string)
{
    try
    {
        var input_file = open(input_name, ioMode.r);
        var input_file_reader = input_file.reader();
        var line: string;
        
        // Get number of bacteria and resize bacteria_name
        number_bacteria = input_file_reader.readLine(): int;
        writeln("number_bacteria: ", number_bacteria); // debugging
        var end = number_bacteria - 1;
        bacteria_domain = {0..end}; // implicty resizing of size of bacteria_name

//        forall i in 0..end do
        for i in 0..end do
        {
//            var line: string;
            input_file_reader.readLine(line);
            line = line.strip('\n'); // safety
            line = line.strip('\t');
            line = line.strip('\r'); // Note: this is the real problem
            bacteria_name[i] = ''.join([line,".faa"]);
            writeln("bacteria_name[", i, "]: ", bacteria_name[i]);
        }
        input_file.close();
        // while (input_file_reader.readLine(line))
        // {
        //     if line != "\n"
        //     {
        //         line = line.strip('\n'); // safety
        //         line = line.strip('\t');
        //         line = line.strip('\r'); // Note: this is the real problem
        //         bacteria_name[i] = ''.join([line,".faa"]);
        //         i += 1;
        //     }
        // }
        
        
    } catch e: EofError {
        writeln("Error: failed to open file %s (Hint: check your working directory)\n");
    } catch e: UnexpectedEofError {
        writeln("unable to read an 'int'");
    } catch e: SystemError {
        writeln("system error in IO implementation: ", e);
    } catch e: Error {
        writeln("something else went wrong...");
    }
}

proc CompareBacteria(b1: Bacteria, b2: Bacteria): real
{
    writeln("b1 count: ", b1.get_count());
    writeln("b2 count: ", b2.get_count());
    var correlation: real = 0.0;
    var vector_len1: real = 0.0;
    var vector_len2: real = 0.0;
    var p1: int = 0;
    var p2: int = 0;
    
    var b1_count: int = b1.get_count();
    var b2_count: int = b2.get_count();
    
    var i: int = 0;
    while p1 < b1_count && p2 < b2_count do
    {
        if i==0 then {writeln("first area"); i+=1;}
        var n1: int = b1.get_ti(p1);
        var n2: int = b2.get_ti(p2);
        
        if n1 < n2 then
        {
            var t1: real = b1.get_tv(p1);
            vector_len1 = vector_len1 + (t1 * t1);
            p1 += 1;
        }
        else if n1 > n2 then
        {
            var t2: real = b2.get_tv(p2);
            p2 += 1;
            vector_len2 = vector_len2 + (t2 * t2);
        }
        else
        {
            var t1: real = b1.get_tv(p1);
            var t2: real = b2.get_tv(p2);
            vector_len1 += (t1 * t1);
            vector_len2 += (t2 * t2);
            correlation += (t1 * t2);
            p1 += 1;
            p2 += 1;
        }
    }
    i=0;
    while (p1 < b1_count)
    {
        if i==0 then {writeln("second area"); i+=1;}
        var t1: real = b1.get_tv(p1);
        vector_len1 += (t1 * t1);
        p1 += 1;
    }
    i=0;
    while (p2 < b2_count)
    {
        if i==0 then {writeln("third area"); i+=1;}
        var t2: real = b2.get_tv(p2);
        vector_len2 += (t2 * t2);
        p2 += 1;
    }
    
    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

proc CompareAllBacteria()
{
    var b: [0..number_bacteria-1] Bacteria;
    var correlation: real;
    
    // Load Bacteria
    for i in 0..number_bacteria-1 do
    {
        writeln("load ", i+1, " of ", number_bacteria);
        b[i] = new Bacteria(bacteria_name[i]);
    }
    
    // Comparison
    for i in 0..number_bacteria-2 do
    {
        for j in i+1..number_bacteria-1 do
        {
            correlation = CompareBacteria(b[i], b[j]);
            writeln(i, " -> ", j, ": ", correlation);
        }
    }
}

proc main(){
    var temp: real = 0.03363286346937661986;
    writef("%.20dr\n", temp);
//    var t: stopwatch;
//    t.start();
//    Init();
//    t.stop();
//    writeln("Initialise Vectors: ", t.elapsed()*1000, " milliseconds");
//    t.clear();
//
//    t.start();
//    ReadInputFile("list.txt");
//    t.stop();
//    writeln("Read Input File: ", t.elapsed()*1000, " milliseconds");
//    t.clear();
//
//    t.start();
//    CompareAllBacteria();
//    t.stop();
//    writeln("time elapsed: ",t.elapsed(), " seconds\n",  " seconds");
}

