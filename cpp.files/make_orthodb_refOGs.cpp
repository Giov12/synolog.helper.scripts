#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

using std::string;
using std::fstream;
using std::unordered_map;
using std::vector;
using std::cerr;
using std::cout;

typedef unsigned int uint;

//
// hols different identifiers
// for a specific gene including ensembl
// protein, uniprot, and gene_id
//
struct GeneInfo {
    string ensembl;
    string uniprot;
    string protein_id;
    string gene_id;
};

struct SppNames {
    string orthdb_id;
    string new_id;
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

int
create_spp_map(const string &infile, vector<SppNames> &spp_names){
    //
    // this will create a mapping for
    // orthoDB species ID -> new org ID
    //

    if (infile.empty()){
        return 1; // nothing to do here
    }

    fstream fh(infile);

    if (!fh.is_open()){
        cerr << "Error: Unable to open " << infile << '\n';
        exit(1);
    }

    string line;
    while (std::getline(fh, line)){
        if (line.empty() || line[0] == '#'){
            continue;
        }
        if (line.back() == '\n'){
            line.pop_back();
        }

        SppNames sn;

        size_t idx = line.find('\t');
        if (idx == string::npos){
            continue;
        }
        sn.orthdb_id = line.substr(0, idx);
        sn.new_id    = line.substr(idx + 1);
        spp_names.push_back(sn);
    }

    fh.close();

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

string
extract_geneid(const string &synomoms){
    //
    // helper function to grab the gene ids from the 4th column of
    // odb12v2_genes.tab.gz
    //

    if (synomoms.empty()){
        return "";
    }

    size_t start = 0, next = string::npos, length = synomoms.size();
    string part;

    // GeneID is completely full of digits. Need to check if this the ID
    bool all_digits;

    // iterate over a ';' delimited string
    while (start <= length){
        next = synomoms.find(';', start);
        part = next == string::npos ? synomoms.substr(start) : synomoms.substr(start, next - start);
        
        // assume true unless proven wrong
        all_digits = !part.empty();
        if (all_digits){
            for (uint i = 0; i < part.size(); i++){
                if (!isdigit((unsigned char)part[i])){
                    all_digits = false;
                    break;
                }
            }
        }
        // we found it
        if (all_digits){
            return part;
        }
        // if not, exit if we reached the end
        if (next == string::npos){
            break;
        }

        start = next + 1;
    }

    return "";
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

    // objects used for parsing
    vector<string> parts;
    string line;
    bool   eof;
    long count = 0;

    while (true){
        line = get_gzline(fh, eof);
        if (eof){
            break; // reach the end of the line
        }
        
        // strip new line characters
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
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
        const string &protein_id = parts[2];
        const string &synonyms   = parts[3];
        const string &uniprot    = parts[4];
        const string &ensembl    = parts[5];

        gene_map[orthodb_id] = { strip_version(ensembl), uniprot, protein_id, extract_geneid(synonyms)};
        count++;
    }
    gzclose(fh);
    cerr << "Loaded " << count << " gene records\n";

    return 0;
}

int
join_and_write(const string &og2genes_file, const string &outname,
            unordered_map<string, GeneInfo> &gene_map, vector<SppNames> &spp_names){

    //
    // create a single tsv file with all the orthogroups but use
    // ensembl or uniprot id's instead of orthodb ids
    //

    gzFile fh = gzopen(og2genes_file.c_str(), "rb");
    if (fh == NULL){
        cerr << "Error: could not open " << og2genes_file << '\n';
        exit(1);
    }

    gzFile ofh = gzopen(outname.c_str(), "wb");
    if (ofh == NULL){
        cerr << "Error: could not open " << outname << '\n';
        exit(1);
    }

    gzFile ofh2 = gzopen("missing_ids.tsv.gz", "wb");
    if (ofh2 == NULL){
        cerr << "Error: could not open " << outname << '\n';
        exit(1);
    }

    gzprintf(ofh, "#Orthogroup\tGeneID\tSpecies\tSource(s)\n");

    const int buff_size = 8192;
    char buffer[buff_size];
    vector<string> parts;
    string line;
    bool   eof;

    //
    // several counters
    // number of elements/genes written & number of
    // genes not found in the map & number of genes
    // that lack a usable ID
    //

    long written = 0, not_in_map = 0, no_gene_id = 0;

    while (true){
        line = get_gzline(fh, eof);

        if (eof){
            break;
        }

        // remove last line character
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }
        parts.clear();
        parse_tabular(line, parts);

        const string &og_id           = parts[0]; // orthogroup ID
        const string &orthodb_gene_id = parts[1]; // orthoDB gene ID

        size_t idx = orthodb_gene_id.find(':');
        string species_id = orthodb_gene_id.substr(0, idx);

        // check if we need to update the name
        for (uint i = 0; i < spp_names.size(); i++){
            if (species_id == spp_names[i].orthdb_id){
                species_id = spp_names[i].new_id;
                break;
            }
        }

        auto it = gene_map.find(orthodb_gene_id);
        // is this gene ID in our mapping
        if (it == gene_map.end()){
            not_in_map++;
            continue;
        }

        // down the priority change
        string gene_id, id_source;
        const GeneInfo &info = it->second;
        if (!info.gene_id.empty()){
            gene_id   = info.gene_id;
            id_source = "geneid";
        }
        if (!info.protein_id.empty()){
            gene_id   += ',' + info.protein_id;
            id_source += ",protein_id";
        }
        else if (!info.ensembl.empty()){
            gene_id   += ',' + info.ensembl;
            id_source += ",ensembl";
        }
        else if (!info.uniprot.empty()){
            gene_id   += ',' + info.uniprot;
            id_source += ",uniprot";
        }
        //
        // this gene won't be trackable
        //
        if (gene_id.empty()){
            gzprintf(ofh2, "%s\t%s\t%s\n", og_id.c_str(), orthodb_gene_id.c_str(), species_id.c_str());
            no_gene_id++;
            continue;
        }

        if (gene_id[0] == ','){
            gene_id   = gene_id.substr(1); // no GeneID was found
            id_source = id_source.substr(1);
        }

        gzprintf(ofh, "%s\t%s\t%s\t%s\n", og_id.c_str(), gene_id.c_str(), species_id.c_str(), id_source.c_str());
        written++;
    }

    gzclose(fh);
    gzclose(ofh);
    gzclose(ofh2);

    cerr << "Wrote " << written << " records\n";
    cerr << not_in_map << " not found in gene map (possible species-list mismatch)\n";
    cerr << no_gene_id << " found but had no Ensembl/UniProt id\n";

    return 0;
}

void
help(){
    cerr << "Usage: ./make_orthodb_refOGs -G filtered.odb12v2_OG2genes.tab.gz -g filtered.odb12v2_genes.tab.gz "
         << "-o output.tsv.gz [optional] -s spp_map.tsv [optional]\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string og2genes_file, genes_file, outname, spp_map_file;
    
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
        else if (arg == "-s" && i + 1 < argc){
            spp_map_file = string(argv[i + 1]);
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

    // collect any species ID mappings
    vector<SppNames> spp_names;
    create_spp_map(spp_map_file, spp_names);

    //
    // create an output name if needed
    //
    if (outname.empty()){
        outname = "OrthoDB.RefOGs.tsv.gz";
    }

    join_and_write(og2genes_file, outname, gene_map, spp_names);

    return 0;
}