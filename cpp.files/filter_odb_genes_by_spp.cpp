#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <unordered_set>
#include <vector>
#include <zlib.h>

using std::string;
using std::unordered_set;
using std::vector;
using std::cerr;
using std::cout;
using std::fstream;
using std::to_string;

//
// this utility is design to be a simple filter for odb12v2_OG2genes.tab.gz
// to filter down to a set of specified species
//

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

int
parse_tabular(string &line, vector<string> &parts){

    int start  = 0, end = 0;

    while (end < line.size()){
        if (line[end] == '\t'){
            parts.emplace_back(line.substr(start, end - start));
            start = end + 1;
        }
        end++;
    }

    if (start < line.size()){
        parts.emplace_back(line.substr(start));
    }

    if (parts.size() < 3){
        cerr << "Error: Malformed line\n" << line << '\n';
        exit(1);
    }

    return 0;
}


int 
load_species(const string &infile, unordered_set<string> &species){

    //
    // this function will load a list of species into a set
    //

    fstream fh(infile);

    if (!fh.is_open()){
        cerr << "Error: Unable to open " << infile << '\n';
        exit(1);
    }

    string spp;
    while (std::getline(fh, spp)){
        if (spp.empty() || spp[0] == '#'){
            continue;
        }
        species.insert(spp);
    }

    fh.close();

    cerr << "Loaded " << species.size() << " target species\n";

    return 0;
}

string
get_gzline(gzFile fh, bool &eof){
    //
    // construct a string that reaches the '\n' character
    //
    string line;
    const int buff_size = 8192;
    char buffer[buff_size];
    bool chars_read = false; // were characters read

    while (true){
        char *read_chars = gzgets(fh, buffer, buff_size);

        if (read_chars == NULL){
            break; // reach the end of the file stream
        }
        chars_read = true;
        line      += buffer;
        if (!line.empty() && line.back() == '\n'){
            break;
        }
    }

    eof = !chars_read; // will be true if no characters read
    return line;
}

int
filter_genes(const string &genes_file, unordered_set<string> &species){

    //
    // filter the odb12v2_genes.tab.gz to only records with species
    // from the target set
    //

    gzFile fh = gzopen(genes_file.c_str(), "rb");
    if (fh == NULL){
        cerr << "Error: could not open " << genes_file << '\n';
        exit(1);
    }

    string outname = "genes." + to_string(species.size()) + ".spp.odb12v2_genes.tab.gz";
    gzFile ofh = gzopen(outname.c_str(), "wb");
    if (ofh == NULL){
        cerr << "Error: could not open " << outname << '\n';
        exit(1);
    }

    const int buff_size = 8192;
    char buffer[buff_size];
    vector<string> parts;
    string line;
    bool eof;
    long tot = 0;

    while (true){
        eof  = false;
        line = get_gzline(fh, eof);

        if (eof){
            break; // end of file
        }

        // strip new line characters
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }
        // it was an empty line
        if (line.empty()){
            continue;
        }

        parts.clear();
        parse_tabular(line, parts);

        const string &spp = parts[1]; 

        if (species.count(spp) == 0){
            continue;
        }

        gzputs(ofh, line.c_str());
        gzputs(ofh, "\n");
        tot++;
    }

    gzclose(fh);
    gzclose(ofh);

    cerr << "Wrote " << tot << " records\n";

    return 0;
}

void
help(){
    cerr << "./filter_genes_by_species \n"
         << "-s species_list [single column list of species IDs]\n"
         << "-g genes_file [odb12v2_genes.tab.gz file]\n";
    exit(0);
}

int main(int argc, char *argv[]){

    string spp_file, genes_file;
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }
    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-s" && i + 1 < argc){
            spp_file = string(argv[i + 1]);
        }
        else if (arg == "-g" && i + 1 < argc){
            genes_file = string(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (spp_file.empty() && genes_file.empty()){
        help();
    }

    if (!file_exists(spp_file) || !file_exists(genes_file)){
        cerr << "One or more input files not found\n";
        exit(1);
    }

    // grab the species
    unordered_set<string> species;
    load_species(spp_file, species);
    filter_genes(genes_file,species);

    return 0;
}