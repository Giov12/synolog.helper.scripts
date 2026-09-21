#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include <glob.h>

using std::string;
using std::fstream;
using std::unordered_map;
using std::vector;
using std::cerr;
using std::cout;

typedef unsigned int uint;

struct Entry {
    string key;
    string value;
};

struct GeneMap {
    unordered_map<string, string> gene_id_map;
    unordered_map<string, string> protein_map;
};

// global map of org -> protein & GeneID -> gene_id
unordered_map<string, GeneMap> gene_map;

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

bool
is_directory(const string &path){
    struct stat buffer;

    if (stat(path.c_str(), &buffer) != 0){
        return false;  // path doesn't exist 
    }

    return S_ISDIR(buffer.st_mode);
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
parse_tabular(string &line, vector<string> &parts){

    //
    // parse a '\t' delimited line
    //

    int start  = 0, end = 0;

    //
    // start from an empty vector
    //
    parts.clear();

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

    return 0;
}

int
get_element_info(string &element_srcs, string &element_ids, vector<Entry> &parts){

    //
    // parse both columns that are ',' delimited from the orthogroups file
    //


    //
    // start from an empty vector
    //
    parts.clear();

    vector<string> srcs, ids;

    string part;

    // parse one column at a time
    // split element_srcs on ','
    size_t start = 0, next;
    while (start <= element_srcs.size()){
        next = element_srcs.find(',', start);
        part = next == string::npos ? element_srcs.substr(start) : element_srcs.substr(start, next - start);
        srcs.push_back(part);
        if (next == string::npos){
            break;
        }
        start = next + 1;
    }

    // split element_ids on ','
    start = 0;
    while (start <= element_ids.size()){
        next = element_ids.find(',', start);
        part = next == string::npos ? element_ids.substr(start) : element_ids.substr(start, next - start);
        ids.push_back(part);
        if (next == string::npos){
            break;
        }
        start = next + 1;
    }

    if (srcs.size() != ids.size()){
        cerr << "Error: mismatched counts between IDSource (" << srcs.size()
             << ") and GeneID (" << ids.size() << ") fields: "
             << element_srcs << " / " << element_ids << '\n';
        exit(1);
    }

    for (size_t i = 0; i < srcs.size(); i++){
        Entry e;
        e.key   = srcs[i];
        e.value = ids[i];
        parts.push_back(e);
    }
    return 0;
}

int
parse_attributes(string &attributes, vector<Entry> &parts){

    //
    // parse a ';' delimited string
    //


    if (attributes.empty()){
        return 1;
    }

    size_t start = 0, next = string::npos, length = attributes.size();
    string part;

    // iterate over a ';' delimited string
    while (start <= length){
        next = attributes.find(';', start);
        part = next == string::npos ? attributes.substr(start) : attributes.substr(start, next - start);
        
        // strip whitespace
        uint i = 0;
        while (i < part.size() && part[i] == ' '){
            i++;
        }

       part = part.substr(i);

       if (!part.empty()){
            size_t idx = part.find(' '); // find if we have a key value pair
            if (idx != string::npos){
                Entry e;
                string key   = part.substr(0, idx);
                string value = part.substr(idx + 1);
                // remove qoutes
                if (value.size() >= 2 && value[0] == '"' && value.back() == '"'){
                    value = value.substr(1, value.size() - 2);
                }
                if (key == "db_xref" && value.size() >= 7 && value.substr(0, 7) == "GeneID:"){
                    key   = "GeneID";
                    idx   = value.find(':');
                    value = value.substr(idx + 1); 
                }
                e.key   = key;
                e.value = value;
                parts.push_back(e);
            }
       }
        // we reached the end
        if (next == string::npos){
            break;
        }
        start = next + 1;
    }

    return 0;
}

vector<string>
get_annotations(const string &ann_dir){
    //
    // return a list containg the paths of all
    // the annotations at this directory
    //

    const string pattern = ann_dir.back() == '/' ? 
                          ann_dir + "*.gtf.gz" : ann_dir + "/*.gtf.gz";
    glob_t glob_result;
    
    int return_value = glob(pattern.c_str(), GLOB_TILDE, nullptr, &glob_result);

    // did we fail to find anything
    if (return_value != 0){
        globfree(&glob_result); // empty this countainer
        cerr << "No *.gtf.files found at " << ann_dir << '\n';
        exit(1);
    }

    //
    // create a vector with all the annotations
    //
    vector<string> gtfs;
    for (uint i = 0; i < glob_result.gl_pathc; i++){
        gtfs.push_back(string(glob_result.gl_pathv[i]));
    }

    globfree(&glob_result);
    return gtfs;

}

string
get_spp_id(const string &ann){
    //
    // get the species synolog id from this annotation file
    // spp.def.gtf.gz <- file naming scheme
    //

    size_t idx = ann.find_last_of('/');
    
    // get basename
    string bname = idx == string::npos ? ann : ann.substr(idx + 1);
    idx = bname.find('.'); // get position of first '.'
    
    return idx == string::npos ? bname : bname.substr(0, idx);
}

GeneMap
parse_annotation(const string &ann){
    //
    // return a protein to ensembl mapping
    //

    gzFile fh = gzopen(ann.c_str(), "rb");

    if (fh == NULL){
        cerr << "Error: could not open " << ann << '\n';
        exit(1);
    }

    //
    // create the two mappings
    //
    GeneMap org_map;

    //
    // create the objects we need to store info
    //
    vector<string> parts;
    vector<Entry> entries;
    string line, prot_id, gene_id, GeneID; 
    bool eof;
    uint missing_prots = 0;

    while (true){
        line = get_gzline(fh, eof);

        if (eof){
            break; // end of file
        }
        if (line.empty() || line[0] == '#'){
            continue; // skip comment & empty lines
        }
        parse_tabular(line, parts);
        if (parts.size() < 9){
            cerr << "Malformed line in " << ann << '\n' << line;
            exit(1);
        }
        if (parts[2] == "CDS" || parts[2] == "cds" || parts[2] == "gene"){
            // parse through everything
            entries.clear();
            prot_id.clear();
            GeneID.clear();
            gene_id.clear();
            parse_attributes(parts[8], entries);
            for (uint i = 0; i < entries.size(); i++){
                Entry &e = entries[i];
                if (e.key == "protein_id"){
                    prot_id = e.value;
                }
                else if (e.key == "GeneID"){
                    GeneID = e.value;
                }
                else if (e.key == "gene_id"){
                    gene_id = e.value;
                }
            }
            if (!gene_id.empty()){
                if (parts[2] != "gene" && prot_id.empty()){
                    // cerr << "Unable to find protein_id in " << ann << '\n' << line;
                    missing_prots++;
                }
                else if (!prot_id.empty()) {
                   org_map.protein_map[prot_id] = gene_id; 
                }    
                if (!GeneID.empty()){
                    org_map.gene_id_map[GeneID] = gene_id;
                }
            }
            else {
                cerr << "No gene_id found in " << ann << " for record\n" << line;
                exit(1);
            }
        } // end of cds parse
    } // end of parsing

    gzclose(fh);

    cerr << "Unable to find protein_ids for " << missing_prots << " records in " << ann << '\n';
    cerr << "Recovered protein_ids for " << org_map.protein_map.size() << " records in " << ann << '\n';

    return org_map;
}

int
create_gene_map(const string &ann_dir){
    //
    // function orchestrate getting & parsing the annotations
    //
    
    vector<string> anns = get_annotations(ann_dir);

    string spp;
    for (uint i = 0; i < anns.size(); i++){
        spp = get_spp_id(anns[i]);
        gene_map[spp] = parse_annotation(anns[i]);
    }

    return 0;
}

string
make_output_name(const string &infile){
    //
    // add a prefix to this file
    //
    
    size_t idx = infile.find_last_of('/');

    string bname = idx == string::npos ? infile : infile.substr(idx + 1);

    return "renamed." + bname;
}

int
rename_orthologs(const string &orthogroups_file){

   //
   // create a new file with gene_ids used
   // to reflect synolog's gene naming scheme
   //

    gzFile fh = gzopen(orthogroups_file.c_str(), "rb");
    if (fh == NULL){
        cerr << "Error: could not open " << orthogroups_file << '\n';
        exit(1);
    }

    string outname = make_output_name(orthogroups_file);
    gzFile ofh = gzopen(outname.c_str(), "wb");
    if (ofh == NULL){
        cerr << "Error: could not open " << outname << '\n';
        exit(1);
    }

    vector<string> parts;
    vector<Entry> entries;
    string line, element_ids, element_srcs, spp;
    string grp_id, gene_id, prot_id, src;
    bool   eof;
    uint bad = 0;

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

        if (line[0] == '#'){
            gzputs(ofh, line.c_str());
            continue;
        }


        // remove last line character
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')){
            line.pop_back();
        }

        if (line.empty()){
            continue;
        }
        parse_tabular(line, parts);

        grp_id       = parts[0];
        spp          = parts[2];
        element_ids  = parts[1];
        element_srcs = parts[3];

        get_element_info(element_srcs, element_ids, entries);
        gene_id.clear();
        prot_id.clear();
        src.clear();

        for (uint i = 0; i < entries.size(); i++){
            Entry &e = entries[i];
            string element_src = e.key;
            string element_id  = e.value;
            
            if (gene_id.empty() && element_src == "geneid"){
                gene_id = gene_map[spp].gene_id_map[element_id];
            }
            if (gene_id.empty() && element_src == "protein_id"){
                gene_id = gene_map[spp].protein_map[element_id];
            }
            if (element_src == "protein_id"){
                prot_id = element_id;
            }
            if (!gene_id.empty()){
                src = "renamed";
                break;
            }
        }
        if (gene_id.empty()){
            // cerr << "Error: Unable to find gene_id for " << element_ids << " in " << spp << '\n';
            // exit(1);
            gene_id = prot_id.empty() ? element_ids : prot_id;
            src = element_srcs;
            bad++;
        }
        else{
            written++;
        }
        // else {
        //     cerr << "Error: unrecognized element source " << element_src 
        //          << " for " << element_id << " in " << spp << '\n';
        //     exit(1);
        // }

        gzprintf(ofh, "%s\t%s\t%s\t%s\n", grp_id.c_str(), gene_id.c_str(), spp.c_str(), src.c_str());
    }

    cerr << "Failed to get IDs for " << bad << " records\n";
    cerr << "Renamed " << written << " records\n";

    gzclose(fh);
    gzclose(ofh);

    return 0;
}

void
help(){
    cerr << "Usage: ./rename_odb_genes -o orthodb_orthogroups.tsv.gz -a /path/to/gtf_annotations/\n";
    exit(1);
}

int main(int argc, char *argv[]){

    string orthogroups_file, ann_dir;
    
    // expect at least 2 inputs
    if (argc < 3){
        help();
    }

    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if (arg == "-o" && i + 1 < argc){
            orthogroups_file = string(argv[i + 1]);
        }
        else if (arg == "-a" && i + 1 < argc){
            ann_dir = string(argv[i + 1]);
        }
        else if (arg == "-h"){
            help();
        }
    }
    if (orthogroups_file.empty() && ann_dir.empty()){
        help();
    }

    if (!file_exists(orthogroups_file)){
        cerr << "Unable to find " << orthogroups_file << '\n';
        exit(1);
    }
    if (!is_directory(ann_dir)){
        cerr << ann_dir << " is not a directory\n";
        exit(1);
    }

    // populate the gene map
    create_gene_map(ann_dir);

    // now rewrite
    rename_orthologs(orthogroups_file);

    return 0;
}