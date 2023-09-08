// cd /Users/Shridhar/Library/Developer/Xcode/XCodeProjects/CAB401_Project/cvtree_project/cvtree_project
// chpl -o tidy_chpl tidy.chpl; ./tidy_chpl
use Math;
use Time;
use Random;

var number_bacteria: int;
var bacteria_name: string;
var M, M1, M2: int;
const LEN: int = 6;
const AA_NUMBER: int = 20; // number of amino acids
const EPSILON: real= 1e-010; //0.0000000001;

var D: domain(1) = {0..25};
var amino_acids: [D] string = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'];
var amino_vals: [D] int = [0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3];

proc encode(amino: string): int{
    // TODO: forall for data-thread, coforall for thread
    for i in D.dim(0) do{
        if amino_acids[i] == amino then
            return amino_vals[i];
    }
    return -1;
}

proc Init()
{
    M2 = 1;
    for i in 1..LEN-2 do {
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
    var second: [0..M1-1] int;
    var one_l: [0..AA_NUMBER-1] int;
    var indexs, total, total_l, complement: int;
    var vector: [0..M-1] int;
    
    proc init_buffer(buffer: [0..LEN-2] string)
    {
        complement += 1;
        indexs = 0;
        for i in 0..LEN-2 do{
            enc = encode(buffer[i]);
            one_l[enc] += 1;
            total_l += 1;
            indexs = indexs * AA_NUMBER + enc;
        }
        second[indexs] += 1;
    }
    proc cont_buffer(ch: string)
    {
        enc = encode(ch);
        one_l[enc] += 1;
        total_l += 1;
        index = (indexs * AA_NUMBER) + enc;
        vector[index]  += 1;
        total
    }
}

proc main(){
    Init();
//    var t: stopwatch;
//    t.start();
//    for i in 0..LEN-3 do {
//        M2 *= AA_NUMBER; // M2 = AA_NUMBER ^ (LEN-2);
//    }
//    t.stop();
//    writeln(t.elapsed(), " milliseconds");
//    t.clear();
//
//    t.start();
//    M2 = AA_NUMBER ** 2; // M2 = AA_NUMBER ^ (LEN-2);
//    t.stop();
//    writeln(t.elapsed(), " milliseconds");
//    t.clear();
    
}



