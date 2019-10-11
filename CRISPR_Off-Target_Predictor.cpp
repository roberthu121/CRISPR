/**
 * File: CRISPR_Off-Target_Predictor.cpp
 * ----------------
 * Robert Hu
 *
 * Predicts the sites of potential off-target edits by CRISPR/Cas9
 * editing given a guide RNA strand and the genome to be tested. Also
 * allows for the user to pick their PAM sequence of choice. Uses a
 * BLAST-like algorithim to shorten the run time of the program by selecting
 * a key word (the last few nucleotides at the 3' end) before checking
 * for potential sites with a number of mismatches under the restriction.
 * The data is then printed to show all the off-target sequences and their
 * position and chromosome
 */

#include <iostream>
#include <fstream>
#include <string>
#include "console.h"
#include "random.h"
#include "simpio.h"
#include "vector.h"
#include "strlib.h"
#include "map.h"
#include "filelib.h"
using namespace std;

/**
 * Struct that stores the index, the sequence, and the chromosome of an off-target location
 */
struct offTarget {
    int index;
    string sequence;
    string chromosome;
};

/**
 * Generates a random nucleotide (A, T, G, or C)
 */
char getRandomNucleotide() {
    char nucleotide;
    int num = randomInteger(1, 4);
    switch (num)
    {
    case 1:
        nucleotide = 'A';
        break;
    case 2:
        nucleotide = 'T';
        break;
    case 3:
        nucleotide = 'G';
        break;
    case 4:
        nucleotide = 'C';
        break;
    }
    return nucleotide;
}

/**
 * Prompts the user for a length and then generates a random sequence of DNA
 * of that length
 */
string generateRandomSequence() {
    int lengthOfDNA = getInteger("Enter the length of a random DNA sequence: ");
    string DNA = "";
    for (int i = 0; i < lengthOfDNA; i++) {
        DNA += getRandomNucleotide();
    }
    cout << "Random Sequence Generated." << endl;
    return DNA;
}

/**
 * Prompts the user to pick a Protospacer Adjacent Motif sequence, which is
 * right behing the 3' end of the crDNA
 */
string promptUserPAM() {
    string PAM;
    cout << "1. SpCas9 from Streptococcus pyogenes: 5'-NGG-3' " << endl;
    cout << "2. SpCas9 from Streptococcus pyogenes: 5'-NRG-3' (R = A or G) " << endl;
    cout << "3. StCas9 from Streptococcus thermophilus: 5'-NNAGAAW-3' (W = A or T) " << endl;
    cout << "4. NmCas9 from Neisseria meningitidis: 5'-NNNNGMTT-3' (M = A or C) " << endl;
    cout << "5. SaCas9 from Staphylococcus aureus: 5'-NNGRRT-'3 (R = A or G) " << endl;
    cout << "6. CjCas9 from Campylobacter jejuni: 5'-NNNVRYAC-3' (V = G or C or A, R = A or G, Y = C or T) " << endl;
    cout << "7. CjCas9 from Campylobacter jejuni: 5'-NNNNRYAC-3' (R = A or G, Y = C or T) " << endl;
    cout << "8. AsCpf1 from Acidaminococcus or LbCpf1 from Lachnospiraceae: 5'-TTTN-3' " << endl;
    cout << "9. AsCpf1 from Acidaminococcus or LbCpf1 from Lachnospiraceae: 5'-TTTV-3' (V = G or C or A) " << endl;
    cout << "10. SpCas9 from Streptococcus pasteurianus: 5'-NNGTGA-3' " << endl;
    cout << "11. FnCpf1 from Francisella: 5'-TTN-3' " << endl;
    cout << "12. SaCas9 from Staphylococcus aureus: 5'-NNNRRT-'3 (R = A or G) " << endl;
    int num = getInteger("Choose the number of one of the PAM (Protospacer Adjacent Motif) sequences above: ");
    switch (num) {
    case 1:
        PAM = "NGG";
        break;
    case 2:
        PAM = "NRG";
        break;
    case 3:
        PAM = "NNAGAAW";
        break;
    case 4:
        PAM = "NNNNGMTT";
        break;
    case 5:
        PAM = "NNGRRT";
        break;
    case 6:
        PAM = "NNNVRYAC";
        break;
    case 7:
        PAM = "NNNNRYAC";
        break;
    case 8:
        PAM = "TTTN";
        break;
    case 9:
        PAM = "TTTV";
        break;
    case 10:
        PAM = "NNGTGA";
        break;
    case 11:
        PAM = "TTN";
        break;
    case 12:
        PAM = "NNNRRT";
        break;
    }
    return PAM;
}

/**
 * Checks to see whether the guide RNA sequence the user enters
 * is valid (contains only the nucleotides A, U, G, and C
 */
bool validcrRNA(const string crRNA, const int lengthOfcrRNA) {
    if (crRNA.size() != lengthOfcrRNA) return false;
    for (int i = 0; i < lengthOfcrRNA; i++) {
        if (!isalpha(crRNA[i])) return false;
        else {
            char c = toUpperCase(crRNA[i]);
            if (c != 'A' && c != 'U' && c != 'G' && c != 'C') return false;
        }
    }
    return true;
}

/**
 * Changes the guide RNA into a CRISPR DNA sequence by swapping
 * all appearances of a U with a T
 */
string crRNATocrDNA (const string crRNA) {
    string crDNA = "";
    for (int i = 0; i < crRNA.size(); i++) {
        if (crRNA[i] == 'U') {
            crDNA += 'T';
        } else {
            crDNA += crRNA[i];
        }
    }
    return crDNA;
}

/**
 * Prompts the user to enter a guide RNA sequence of around 20 nucleotides long
 * and checks the validity of the sequence before changing it to CRISPR DNA
 *
 * Prompts the user to restrict the number of mismatches that the search will allow
 */
string promptUsercrDNA (int& numMismatches) {
    int lengthOfcrRNA = getInteger("Enter the length of the guide RNA sequence (between 15 - 25 nucleotides): ");
    while (lengthOfcrRNA < 15 || lengthOfcrRNA > 25) {
        lengthOfcrRNA = getInteger("Not in the correct range of 15 - 25, please re-enter the length: ");
    }
    string crRNA = getLine("Enter the 5' to 3' " + integerToString(lengthOfcrRNA) + "-nucleotide guide RNA sequence to be edited with CRISPR: ");
    while (!validcrRNA(crRNA, lengthOfcrRNA)) {
        if (crRNA.size() != lengthOfcrRNA) {
            crRNA = getLine("Sequence entered is not " + integerToString(lengthOfcrRNA)
                            + " nucleotides long, please re-enter the sequence: ");
        } else {
            crRNA = getLine("The sequence entered is invalid, please re-enter: ");
        }
    }
    crRNA = toUpperCase(crRNA);
    numMismatches = getInteger("Enter the number of mismatches allowed (between 0 - 9): ");
    while (numMismatches < 0 || numMismatches > 9) {
        numMismatches = getInteger("Not in the correct range of 0 - 9, please re-enter the number of mismatches: ");
    }
    return crRNATocrDNA(crRNA);
}

/**
 * Reads the file using and ifstream and concatenates line by line of the files
 * into the string that represents the DNA strand
 */
void readFile(string& DNA, const string filename) {
    ifstream infile;
    infile.open(filename.c_str());
    string line;
    getline(infile, line);
    while (getline(infile, line)) {
        DNA += line;
    }
    infile.close();
}

/**
 * Reads the files of all human chromosomes and stores them in a map
 * that has the chromosome as the key
 */
void inputHumanChromosomes(Map<string, string>& humanChromosomes) {
    for (int i = 1; i < 25; i++) {
        string DNA;
        string chromosome = integerToString(i);
        if (i == 23) {
            chromosome = "X";
        } else if (i == 24) {
            chromosome = "Y";
        }
        string filename = "human_chromosome_" + chromosome + ".txt";
        readFile(DNA, filename);
        humanChromosomes.add(chromosome, DNA);
        cout << "Chromosome " + chromosome  + " inputted." << endl;
    }
}

/**
 * Prompts the user to either import a file, import the entire human genome, or generate
 * a random sequence of DNA
 */
string promptUserDNA(Map<string, string>& humanChromosomes, string& chromosome) {
    string DNA = "";
    bool import = getYesOrNo("Do you want to import a sequence? ");
    if (import) {
        bool human = getYesOrNo("Do you want to import the human genome? ");
        if (human) {
            inputHumanChromosomes(humanChromosomes);
        } else {
            string filename = getLine("Enter the name of the file for the chromosome: ");
            while (!isFile(filename)) {
                filename = getLine("Invalid file name, please re-enter: ");
            }
            chromosome = getLine("Chromosome name: ");
            readFile(DNA, filename);
        }
    } else {
        DNA = generateRandomSequence();
    }
    return DNA;
}

/**
 * Checks the validity of the potential PAM sequence in the DNA where a PAM sequence
 * is expected and returns the actual PAM sequence (which replaces the non-exclusive
 * filler nucleotides with the actual nucleotides that reside there)
 */
bool validPAM(string& currentPAM, const string PAM, const string potentialPAM) {
    for (int i = 0; i < PAM.length(); i++) {
        if (potentialPAM[i] == PAM[i]) {
            currentPAM += potentialPAM[i];
        } else if (PAM[i] == 'N') {
            if (potentialPAM[i] == 'A' || potentialPAM[i] == 'T' || potentialPAM[i] == 'G' || potentialPAM[i] == 'C') {
                currentPAM += potentialPAM[i];
            } else return false;
        } else if (PAM[i] == 'R') {
            if (potentialPAM[i] == 'A' || potentialPAM[i] == 'G') {
                currentPAM += potentialPAM[i];
            } else return false;
        } else if (PAM[i] == 'W') {
            if (potentialPAM[i] == 'A' || potentialPAM[i] == 'T') {
                currentPAM += potentialPAM[i];
            } else return false;
        } else if (PAM[i] == 'V') {
            if (potentialPAM[i] == 'G' || potentialPAM[i] == 'C' || potentialPAM[i] == 'A') {
                currentPAM += potentialPAM[i];
            } else return false;
        } else if (PAM[i] == 'Y') {
            if (potentialPAM[i] == 'C' || potentialPAM[i] == 'T') {
                currentPAM += potentialPAM[i];
            } else return false;
        } else return false;
    }
    return true;
}

/**
 * Recursively compares the remainder of the guide sequence with the DNA at
 * matching positions, storing the sequence as an offTarget struct in the
 * case that it fulfills the entire length underneath the designated amount
 * of mismatches allowed
 */
void compareSequence(Map<int, Vector<offTarget> >& offTargets, const string remainder, const string remainderDNA, string sequence,
                     int index, int count, int mismatches, const int numMismatches, const string key, const string chromosome) {
    // Base case where the number of mismatches is great than the amount allowed
    if (mismatches > numMismatches) return;
    // Base case where the entire remainder sequence is examined within
    // the constraints of the number of mismatches
    if (count == 0) {
        offTarget current = {index + 1, sequence + key, chromosome};
        if (!offTargets.containsKey(mismatches)) {
            Vector<offTarget> offSequences;
            offSequences.add(current);
            offTargets.add(mismatches, offSequences);
        } else {
            offTargets[mismatches].add(current);
        }
        return;
    }
    // If the nucleotide is different, the number of mismatches is incremented and
    // the differing nucleotide is indicated by being lowercase
    if (remainderDNA[count - 1] != remainder[count - 1]) {
        sequence = toLowerCase(remainderDNA[count - 1]) + sequence;
        compareSequence(offTargets, remainder, remainderDNA, sequence, index - 1, count - 1, mismatches + 1, numMismatches, key, chromosome);
    } else {
        sequence = remainderDNA[count - 1] + sequence;
        compareSequence(offTargets, remainder, remainderDNA, sequence, index - 1, count - 1, mismatches, numMismatches, key, chromosome);
    }
}

/**
 * Analyzes the sequence of DNA it receives using a BLAST-like algorithm,
 * first searching for key sequences that are exactly the same as the
 * end of the crDNA near the 3' end, then checking for the prescence
 * of a PAM sequence immediately adjacent to the end. Successful indices
 * are mapped to their corresponding PAM sequence before all these sequences
 * corresponding to the indices are compared recursively against the crDNA
 */
void analyzeSequence(Map<int, Vector<offTarget> >& offTargets, const string DNA, const string crDNA, const string PAM, const int numMismatches, const string chromosome) {
    Map<int, string> indices;
    int sizeOfKey = 3;
    string key = crDNA.substr(crDNA.size() - sizeOfKey, sizeOfKey);
    string remainder = crDNA.substr(0, crDNA.size() - sizeOfKey);
    for (int i = 0; i < DNA.size(); i++) {
        if (DNA[i] == key[0]) {
            if (DNA.substr(i, sizeOfKey) == key) {
                string currentPAM = "";
                if (validPAM(currentPAM, PAM, DNA.substr(i + sizeOfKey, PAM.size()))) {
                    indices.add(i, currentPAM);
                }
            }
        }
    }
    Vector<int> keys = indices.keys();
    for (int i = 0; i < keys.size(); i++) {
        compareSequence(offTargets, remainder, DNA.substr(keys[i] - remainder.size(), remainder.size()), "", keys[i] - 1, remainder.size(), 0, numMismatches, key + indices[keys[i]], chromosome);
    }
}

/**
 * Displays the data separated by the number of mismatches, labeled
 * with their chromosome and position number
 */
void showData(const Map<int, Vector<offTarget> >& offTargets) {
    Vector<int> keys = offTargets.keys();
    for (int i = 0; i < keys.size(); i++) {
        cout << integerToString(keys[i]) + " mismatches " + "(" + integerToString(offTargets[keys[i]].size()) + ")" << endl;
        cout << "--------------------------------------------------------" << endl;
        for (int j = 0; j < offTargets[keys[i]].size(); j++) {
            cout << "Chromosome " + offTargets[keys[i]][j].chromosome + " Position " + integerToString(offTargets[keys[i]][j].index) + ": " + offTargets[keys[i]][j].sequence << endl;
        }
        cout << endl;
    }
}

/**
 * Entry point to the program.
 */
int main() {
    Map<int, Vector<offTarget> > offTargets;
    int numMismatches;
    Map<string, string> humanChromosomes;
    string chromosome;
    string PAM = promptUserPAM();
    string crDNA = promptUsercrDNA(numMismatches);
    string DNA = promptUserDNA(humanChromosomes, chromosome);
    cout << "Predicting off-target edits... " << endl;
    // Checks to analyze either the human genome or another sequence
    if (humanChromosomes.size() > 0) {
        for (int i = 1; i < 25; i++) {
            string chromosome = integerToString(i);
            if (i == 23) {
                chromosome = "X";
            } else if (i == 24) {
                chromosome = "Y";
            }
            analyzeSequence(offTargets, humanChromosomes[chromosome], crDNA, PAM, numMismatches, chromosome);
            cout << "Chromosome " + chromosome + " sequenced." << endl;
        }
    } else {
        analyzeSequence(offTargets, DNA, crDNA, PAM, numMismatches, chromosome);
    }
    getLine("Press enter to continue");
    showData(offTargets);
    return 0;
}
