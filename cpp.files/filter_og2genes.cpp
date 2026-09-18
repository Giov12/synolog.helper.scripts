//
// parse a bam file to find reads bridging 2 seperate contigs
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <string>
#include <unordered_set>
#include <vector>
#include <cstdio>
#include <cctype>
#include <zlib.h>

// make sure to compile with -lz

using std::string;
using std::unordered_set;
using std::vector;
using std::fstream;
using std::iostream;
using std::cout;
using std::cerr;
using std::to_string;

typedef unsigned int uint;

int
parse_tabular(string &line, vector<string> &parts){

    int start  = 0, end = 0;
    int colons = 0;

    while (end < line.size()){
        if (line[end] == '\t'){
            parts.emplace_back(line.substr(start, end - start));
            start = end + 1;
        }
        else if (line[end] == ':'){
            colons++;
        }
        end++;
    }

    if (start < line.size()){
        parts.emplace_back(line.substr(start));
    }

    if (colons != 1 || parts.size() != 2){
        cerr << "Invalid line found\n" << line;
        exit(1);
    }

    return 0;
}

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}



int parse_targets(string &infile, unordered_set<string> &targets, bool is_og){
    //
    // parse the targets file and add the records to the target set
    //

    fstream file(infile);

    if (!file.is_open()){
        cerr << "Error: Could not open " << infile << '\n';
        exit(1);
    }

    string target;
    while (std::getline(file, target)){
        targets.insert(target);
    }

    file.close();

    string t = is_og ? "groups" : "species";

    cerr << "Loaded " << targets.size() << " target " << t << '\n';

    return 0;
}

int
filter_orthogroups(string og2genes_file, unordered_set<string> &og_ids, unordered_set<string> &spp_ids){
    //
    // this is the main work horse for this function
    //

    int spp_count = spp_ids.size(), og_count = og_ids.size();

    gzFile fh = gzopen(og2genes_file.c_str(), "rb");

    if (fh == NULL){
        cerr << "Error: could not open " << og2genes_file << '\n';
        exit(1);
    }

    string outname = "filtered." + to_string(og_count) + ".OGs." + to_string(spp_count) +
                     ".spp.odb12v2_OG2genes.tab.gz";

    gzFile ofh = gzopen(outname.c_str(), "wb");
    
    if (ofh == NULL){
        cerr << "Error: could not open " << outname << '\n';
        exit(1);
    }

    // create the buffer
    const int buff_size = 8192;
    char buffer[buff_size];
    long tot = 0; // keep count

    // start parsing
    vector<string> parts;
    while (gzgets(fh, buffer, buff_size) != NULL){
        string line(buffer);

        if (!line.empty() && line.back() == '\n'){
            // strip new line char
            line.pop_back();
        }
        if (line.empty()){
            continue;
        }
        // clear then fill
        parts.clear();
        parse_tabular(line, parts);

        const string og_id = parts[0];

        // first checkpoint
        if (og_ids.count(og_id) == 0){
            continue;
        }

        // second checkpoint
        const string spp_id = parts[1].substr(0, parts[1].find(':'));

        if (spp_ids.count(spp_id) == 0){
            continue;
        }

        // write this record
        gzprintf(ofh, "%s\t%s\n", og_id.c_str(), parts[1].c_str());
        tot++;

    }

    // close out
    gzclose(fh);
    gzclose(ofh);

    cerr << "Retained " << tot << " records\n";

    return 0;
}

void
help(){
    //
    // help message
    //
    cerr << "./filter_og2_genes \n"
         << "-s species_list [single column list of species IDs]\n"
         << "-g og_list [single column list of OG to select]\n"
         << "-o [odb12v2_OG2genes.tab.gz file]\n";
    exit(0);

}

int
main(int argc, char *args[]){

    string og_file, spp_file, og2genes_file;
    unordered_set<string> spp_ids, og_ids;

    for (int i = 1; i < argc; i++){
        string arg = args[i];
        if (arg == "-s" && i + 1 < argc){
            spp_file = string(args[i + 1]);
        }
        else if (arg == "-g" && i + 1 < argc){
            og_file = string(args[i + 1]);
        }
        else if (arg == "-o" && i + 1 < argc){
            og2genes_file = string(args[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }

    if (og_file.size() == 0 || !file_exists(og_file)){
        cerr << "Unable to use orthogroup list " << og_file << '\n';
        exit(1);
    }
    else if (spp_file.size() == 0 || !file_exists(spp_file)){
        cerr << "Unable to use species list " << spp_file << '\n';
        exit(1);
    }
    else if (og2genes_file.size() == 0 || !file_exists(og2genes_file)){
        cerr << "Unable to use odb12v2_OG2genes.tab.gz file provided as " << og2genes_file << '\n';
        exit(1);
    }

    //
    // step 1: load all the targets
    //
    parse_targets(og_file, og_ids, true);
    parse_targets(spp_file, spp_ids, false);

    //
    // step 2: write the output file as we go
    //
    filter_orthogroups(og2genes_file, og_ids, spp_ids);

    return 0;
}