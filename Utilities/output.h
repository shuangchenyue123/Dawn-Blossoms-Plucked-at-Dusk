// Created by zhangshengbin on 2021/11/29.
#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "DBscan.h"
#include "other_function.h"
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>


class Output {
public:
    Output(OthoHGT &obj_ortho, DBscan &obj_cluster, map<string, string> &prot_name_sequence_map):_obj_cluster(obj_cluster), _obj_ortho(obj_ortho), _prot_name_sequence_map(prot_name_sequence_map){}

    inline void out_event_result();
    inline void out_HGT_percent_result();
    inline void out_ortho_protein_parallel(int nthreads, string &path);
    inline void out_ortho_protein(my_iter begin, my_iter end, const string &final_path);

    ~Output(){}

private:
    OthoHGT &_obj_ortho;
    DBscan &_obj_cluster;
    map<string, string> &_prot_name_sequence_map;
};



inline void Output::out_event_result(){
    ofstream out_e("Events.tsv");
    if(! out_e)
        handle_error("Error: can't open output file: Events.tsv");
    else{
        // 这些时间的意思都是距离今天多少年
        int pos = 0;
        out_e << "gene_present_time(mya)\thorizontal_transfer_time(mya)\tevent_confidence\torthogroup\tcluster" << endl;
        for(auto &it:_obj_cluster.p_h_time_vec){
            out_e << it[0] << "\t" << it[1] << "\t" << _obj_cluster.event_confidence_vec[pos] << "\t" << _obj_cluster.event_ortho_vec[pos] << "\t" << "cluster" <<_obj_cluster.event_cluster_vec[pos] << endl;
            pos++;
        }
    }
    out_e.close();
    cout << "Output Events finish." << endl;
}



inline void Output::out_HGT_percent_result() {
    ofstream out_e("genome_HGT.tsv");
    if(! out_e)
        handle_error("Error: can't open output file: genome_HGT.tsv");
    else{
        out_e << "genome\tHGT_gene_number_min\tHGT_gene_number_max" << endl;
        for(auto &i:_obj_ortho.genome_HGT_num){
            out_e << change_name_back(i.first) << '\t' << i.second[0] << '\t' << i.second[1] << endl;
        }
    }
    out_e.close();
    cout << "Output genome HGT gene number finish." << endl;
}




inline void Output::out_ortho_protein_parallel(int nthreads, string &path){
    vector<thread> threads(nthreads);
    const string final_path = path + "Every_ortho_sequences/";
    const string command = "mkdir -p " + final_path;
    int mkdir_result = system(command.c_str());
    if(mkdir_result == -1)
        handle_error("Unable to make a new folder for the ortho protein sequences");
    if(_obj_ortho.HGTortho_rtree.size() < nthreads)
        nthreads = _obj_ortho.HGTortho_rtree.size();
    int grainsize = _obj_ortho.HGTortho_rtree.size()/nthreads + 1;
    map<string, vector<string> > order_change_map;
    move_map_order(_obj_ortho.HGTortho_rtree, order_change_map, nthreads, grainsize);

    auto temp_iter = order_change_map.begin();
    auto work_iter_b = temp_iter;
    for(int i = 0; i < grainsize; i++)
        temp_iter++;
    auto work_iter_e = temp_iter;

    for(auto it = threads.begin(); it != threads.end() - 1; ++it) {
        // 原理就是把 map 中的数据拆开了，一部分一部分进行处理
        // 等号右边创建了一个thread对象，将其拷贝初始化给线程池中的某个thread
        *it = thread(&Output::out_ortho_protein, this, work_iter_b, work_iter_e, final_path);
        auto temp_iter2 = work_iter_b;
        work_iter_b = work_iter_e;
        for (int i = 0; i < grainsize; i++)
            temp_iter2++;
        work_iter_e = temp_iter2;
    }
    threads.back() = thread(&Output::out_ortho_protein, this, work_iter_b, order_change_map.end(), final_path);
    for(auto &i:threads)
        i.join();
    cout << "Output protein sequence file finish." << endl;
}


inline void Output::out_ortho_protein(my_iter begin, my_iter end, const string &final_path){
    for(auto it = begin; it != end; ++it){
        ofstream out_p(final_path + it->first + "_protein_summary.fasta");
        for(auto &prot:_obj_ortho.ortho_prot[it->first]){
            auto tar_iter = _prot_name_sequence_map.find(prot);
            if(tar_iter == _prot_name_sequence_map.end()){
                cerr << "Protein name doesn't find in the summary protein file: " << prot << endl;
                exit(255);
            }
            out_p << ">" << change_name_back(prot) << endl;
            out_p << tar_iter->second;
        }
        out_p.close();
    }
}

#endif
