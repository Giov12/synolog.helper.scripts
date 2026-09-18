#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

using std::string;
using std::unordered_map;
using std::vector;
using std::cerr;
using std::cout;

//
// prioritize to ensembl
// & fallback to uniprot
//
struct GeneInfo {
    string ensembl;
    string uniprot;
};

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

    if (start < line.size() && end - start > 1){
        parts.emplace_back(line.substr(start));
    }

    return 0;
}

string
strip_version(const string &ensembl_field){
    
    // empty field
    if (ensembl_field.empty()){
        return ""; 
    }

    // find the first colon
    size_t idx = ensembl_field.find(';');
    string first = (idx == string::npos) ? ensembl_field : ensembl_field.substr(0, idx);

    // remove any version tags
    size_t jdx = first.find('.');

    return (jdx == string::npos) ? first : first.substr(0, jdx);
}

int 
load_genes(const string &genes_file, unordered_map<string, GeneInfo> &gene_map){

    //
    // this function will create a mapping
    // where a gene ID from orthoDB -> {ensembl,uniprot}
    //

    gzFile fh = gzopen(genes_file.c_str(), "rb");
    if (fh == NULL){
        cerr << "Error: could not open " << genes_file << '\n';
        exit(1);
    }

    //
    // buffer & vector to store fields
    //
    const int buff_size = 8192;
    char buffer[buff_size];
    vector<string> parts;
    long count = 0;

    while (gzgets(fh, buffer, buff_size) != NULL){
        string line(buffer);
        
        // strip \n character
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }
        // start with a clean slate
        parts.clear();
        parse_tabular(line, parts);

        // verify
        if (parts.size() < 6){
            cerr << "Error: Malformed line\n" << line << '\n';
            exit(1); 
        }


        const string &orthodb_id = parts[0];
        const string &uniprot    = parts[4];
        const string &ensembl    = parts[5];

        gene_map[orthodb_id] = { strip_version(ensembl), uniprot };
        count++;
    }
    gzclose(fh);
    cerr << "Loaded " << count << " gene records\n";

    return 0;
}

void join_and_write(const string &og2genes_file, const string &outname,
                     unordered_map<string, GeneInfo> &gene_map){

    gzFile fh = gzopen(og2genes_file.c_str(), "rb");
    if (fh == NULL){
        cerr << "Error: could not open " << og2genes_file << '\n';
        exit(1);
    }
    gzbuffer(fh, 1 << 20);

    gzFile ofh = gzopen(outname.c_str(), "wb");
    if (ofh == NULL){
        cerr << "Error: could not open " << outname << '\n';
        exit(1);
    }

    gzprintf(ofh, "Orthogroup\tGeneID\tSpecies\n");

    const int buff_size = 8192;
    char buffer[buff_size];
    vector<string> parts;
    long written = 0, missing = 0;

    while (gzgets(fh, buffer, buff_size) != NULL){
        string line(buffer);
        if (!line.empty() && line.back() == '\n') line.pop_back();
        if (line.empty()) continue;

        split_tab(line, parts);
        const string &og_id = parts[0];
        const string &orthodb_gene_id = parts[1];

        size_t colon = orthodb_gene_id.find(':');
        string species_id = orthodb_gene_id.substr(0, colon);

        auto it = gene_map.find(orthodb_gene_id);
        if (it == gene_map.end()){
            missing++;
            continue;
        }

        // prefer Ensembl gene id; fall back to UniProt if Ensembl is empty
        const string &gene_id = it->second.ensembl.empty() ? it->second.uniprot : it->second.ensembl;

        if (gene_id.empty()){
            missing++;
            continue;
        }

        gzprintf(ofh, "%s\t%s\t%s\n", og_id.c_str(), gene_id.c_str(), species_id.c_str());
        written++;
    }

    gzclose(fh);
    gzclose(ofh);

    cerr << "Wrote " << written << " records, " << missing << " skipped (no gene id found)\n";
}

void
help(){
    cerr << "Usage: ./make_orthodb_refOGs -G filtered.odb12v2_OG2genes.tab.gz -g filtered.odb12v2_genes.tab.gz\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string og2genes_file, genes_file, outname;
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }
    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-G" && i + 1 < argc){
            og2genes_file = string(argv[i + 1]);
        }
        else if (arg == "-g" && i + 1 < argc){
            genes_file = string(argv[i + 1]);
        }
        else if (arg == "-o" && i + 1 < argc){
            outname = string(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (og2genes_file.empty() && genes_file.empty()){
        help();
    }

    if (!file_exists(og2genes_file) || !file_exists(genes_file)){
        cerr << "One or more input files not found\n";
        exit(1);
    }

    // populate the gene map
    unordered_map<string, GeneInfo> gene_map;
    load_genes(genes_file, gene_map);
    join_and_write(og2genes_file, outname, gene_map);

    return 0;
}