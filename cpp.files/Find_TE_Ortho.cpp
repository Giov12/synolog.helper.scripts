#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <tuple>
#include <zlib.h>
#include <cstring>
#include <fstream>
#include <getopt.h>

using std::string;
using std::vector;
using std::cout;
using std::unordered_map;
using std::tuple;
using std::fstream;

// class to hold TE information (location & putative orthologs)
class TE{
    public:
        string family, subfamily;
        int start, end;
        bool assigned = false; // set true when/if assigned
        vector<int> ortho_indices; // indices of the (co-) orthologs
        char orientation;

        void assign(int ortho_index){
            this->assigned = true; // constant update
            ortho_indices.push_back(ortho_index);
        }
        
        TE(string &family, string &subfamily, int start, int end, char orientation){
            this->family = family;
            this->subfamily = family;
            this->start = start;
            this->end = end;
            this->orientation = orientation;
        }

        ~TE(void){
            ortho_indices.~vector(); // destroy contents if present
        }
};

// will loop through the memberships file to get the orthologous loci pairs
class Ortho_Pair{
    public:
        // initialize members
        string spA, spB, spA_chrom, spB_chrom, geneID_A, geneID_B, cluster_id;
        int start_A, end_A, start_B, end_B;

        Ortho_Pair(string &spA, string &spB, string &spA_chrom, string &spB_chrom,
                   string &geneID_A, string &geneID_B, string &cluster_id, int start_A,
                   int end_A, int start_B, int end_B){
                    // assign values
                    this->spA = spA;
                    this->spB = spB;
                    this->spA_chrom = spA_chrom;
                    this->spB_chrom = spB_chrom;
                    this->geneID_A = geneID_A;
                    this->geneID_B = geneID_B;
                    this->cluster_id = cluster_id;
                    this->start_A = start_A;
                    this->end_A = end_A;
                    this->start_B = start_B;
                    this->end_B = end_B;
                   }
        ~Ortho_Pair(void){};
};

// returns a line for an uncompressed file
int readline(fstream &fh, string &line){

    bool eol = false;

    do {
        if (getline(fh, line)){
            eol = true;
        }
    } while (!fh.eof() && !eol);

    if (fh.eof()){
        return 0;
    }

    return 1;
}

// returns a line for a gzipped file
int readgzline(gzFile &fh, char *line, int *size){

    char buffer[*size];
    size_t buffer_len, line_len = 0;
    bool eol = false; // end of line

    memset(line, 0, *size); // fill will 0s

    do{
        if (gzgets(fh, buffer, *size) == NULL){
            break; // done with file
        }
        buffer_len = std::strlen(buffer);
        if (buffer_len > 0 && buffer[buffer_len - 1] == '\n'){
            // reached end of line
            eol = true;
            std::strcat(line, buffer);
        }
        else{
            *size *= 2; // double size to read in
            line_len += buffer_len;
            line = (char *) realloc(line, *size);
            std::strcat(line, buffer); // merge for now
        }
    } while (!gzeof(fh) && !eol);

    if (gzeof(fh)){
        // done with file
        return 0;
    }
    return 1;
}

// split a line to the desired fields
void make_fields(char *line, const char delimiter, vector<string> &fields){
    
    string field = ""; //make empty field
    size_t line_len = std::strlen(line);

    for (int i = 0; i < line_len; i++){
        if (line[i] == delimiter){
            fields.push_back(field);
            field = "";
        }
        else{
            field += line[i];
        }
    }
    // add the last field
    fields.push_back(field);
}

void get_species_ids(const char *mem_tsv, string &spA, string &spB){
    /*
        parse membership fille name to get Species IDs
    */
    string mem_tsv_str = string(mem_tsv);
    string basename = mem_tsv_str.substr(mem_tsv_str.find_last_of("/\\") + 1);
    int index = basename.find("_consyn_cluster_membership.tsv");
    basename = basename.substr(0, index);

    for (int i = 0; i < basename.length(); i++){
        if (basename[i] == '-'){
            spA = basename.substr(0, i - 1);
            spB = basename.substr(i + 1);
            break;
        }
    }
}

void parse_chromsomes_tsv(unordered_map<string, int> &spA_lengths, unordered_map<string, int> &spB_lengths, 
                        string &spA, string &spB, const char *chroms_tsv, const char *mem_tsv, int min_len){

    // first get species IDs
    get_species_ids(mem_tsv, spA, spB);

    cout << "Using " << spA << " as the first species\nUsing " << spB << " as the second species\n";

    // now to go through the file
    vector<string> fields;
    size_t path_len = std::strlen(chroms_tsv);
    int buffersize = 1024;
    char *line = new char[buffersize];

    if (chroms_tsv[path_len - 3] == '.' && chroms_tsv[path_len - 2] == 'g' && chroms_tsv[path_len - 1] == 'z'){
        gzFile gz = gzopen(chroms_tsv, "rb");
        if (gz == NULL){
            std::cerr << "Error: Could not open file " << chroms_tsv << '\n';
            EXIT_FAILURE;
        }
        while (readgzline(gz, line, &buffersize)){
            if (line[0] != '#'){
                make_fields(line, '\t', fields);
                if (fields[0] == spA && std::stoi(fields[2]) >= min_len){
                    spA_lengths[fields[1]] = std::stoi(fields[2]);
                }
                else if (fields[0] == spB && stoi(fields[2]) >= min_len){
                    spB_lengths[fields[1]] = std::stoi(fields[2]);
                }
                fields.clear();
            }
        }
        delete []line;
    }
    else{
        fstream fh;
        delete []line; // using a str instead
        string line;
        fh.open(chroms_tsv, std::ios::in);
        while (readline(fh, line)){
            if (line[0] != '#'){
                char *tmp = new char[line.length() + 1];
                std::strcpy(tmp, line.c_str());
                make_fields(tmp, '\t', fields);
                delete []tmp;
                if (fields[0] == spA && std::stoi(fields[2]) >= min_len){
                    spA_lengths[fields[1]] = std::stoi(fields[2]);
                }
                else if (fields[0] == spB && stoi(fields[2]) >= min_len){
                    spB_lengths[fields[1]] = std::stoi(fields[2]);
                }
                fields.clear();
            }
        }
    }


}

void help(void){

    cout << "./Find_TE_Ortho.cpp\n" << "-m/-memberships_tsv spp1-spp2_consyn_cluster_membership.tsv [required]\n"
        << "-g/-species1_gff species1.gff [required]\n-G/-species2_gff species2.gff [required]\n" <<
        "-c/-chromosomes_tsv chromosomes.tsv [required]\n-M/-minimum_len int [default = 1000000]\n" <<
        "-d/-max_distance int [default = 1000000]\n";
    EXIT_SUCCESS;
}

bool check_if_null(char * file, const char filetype[]){

    if (file == nullptr){
        cout << filetype << " is required\n";
        return true;
    }

    return false;
}

int main(int argc, char * argv[]){

    // some defaults
    uint max_distance = 1000000, min_len = 1000000;
    char * chromosomes_tsv = nullptr;
    char * species1_gff = nullptr;
    char * species2_gff = nullptr;
    char * membership_tsv = nullptr;

    // get the command line options
    int c;
    while (1){
        static struct option long_options[] = {
            {"help", no_argument, NULL, 'h'},
            {"memberships_tsv", required_argument, NULL, 'm'},
            {"species1_gff", required_argument, NULL, 'g'},
            {"species2_gff", required_argument, NULL, 'G'},
            {"chromosomes_tsv", required_argument, NULL, 'c'},
            {"minimum_len", optional_argument, NULL, 'M'},
            {"max_distance", optional_argument, NULL, 'd'},
            {0, 0, 0, 0}
        };

        // according to the docs, getopt_long stores the option index here
        int option_index = 0;

        c = getopt_long(argc, argv, "hm:g:G:c:M:d:", long_options, &option_index);

        if (c == -1){
            break;
        }
        switch (c){
        case 'h': help(); break;
        case 'm': {
            const char * membership_tsv = optarg;
            break; }
        case 'g': {
            const char * species1_gff = optarg;
            break; }
        case 'G': {
            const char * species2_gff = optarg;
            break; }
        case 'c': {
            const char * chromosomes_tsv = optarg;
            break; }
        case 'M': {
            min_len = (uint)std::stoi(optarg);
            break; }
        case 'd': {
            max_distance = (uint)std::stoi(optarg);
            break; }
        default: help(); break;
        }
    }

    // check input files
    if (check_if_null(membership_tsv, "Membership file")) { return 1;}
    if (check_if_null(species1_gff, "Species 1 gff")) { return 1;}
    if (check_if_null(species2_gff, "Species 2 gff")) { return 1;}
    if (check_if_null(chromosomes_tsv, "chromosomes file")) { return 1;}

    // get the species IDs and chromosomes
    unordered_map<string, int> spA_lengths, spB_lengths;
    string spA, spB;

    parse_chromsomes_tsv(spA_lengths, spB_lengths, spA, spB, chromosomes_tsv, membership_tsv, min_len);

    return 0;
}

