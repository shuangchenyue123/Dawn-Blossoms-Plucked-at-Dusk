// Created by zhangshengbin 2021.11.21

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include "other_function.h"
#include "sliding_mode.h"
#include "output.h"
#define WHITE "\033[37m"      /* Magenta */
#define GREEN   "\033[36m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define RESET   "\033[0m"
using namespace std;


void flower(){
    cout << WHITE << "                                        7LI   ..   iv7" << endl;
    cout << "                                ..     XjvR7 2X7r:Lr:i." << endl;
    cout << "                               jJri:  .2irLsL7i.r7::::" << endl;
    cout << "                               i7i:ir  7i7v1:::ir:.:::.rQM.:i." << endl;
    cout << "                          i:    r::.::..:7vs:.:7::.:irPd5SZgBg" << endl;
    cout << "                          QUv7i  r::::ir.:vvi.ir::.ivJ5SXKIjR" << endl;
    cout << "                   i2vr:   ii::i.iri:::7i.rLi:7i:.:KJ15IS12M" << endl;
    cout << "                   :K2Y77r:.ii::::rri:irv:i7iir::iS11S2S127   .iJBvsda" << endl;
    cout << "                      .rrii::::::::rriii7i:i:r7:iUu5SISu5r..uQQBBLdas" << endl;
    cout << "                        .:ii:::::i::i:::::iiir7iYEqKS5IurrEDMQ2isda" << endl;
    cout << "                 M5sLv7i::irrriii::::" << YELLOW << "::rrrri::::r7" << WHITE << "KbXJs7udDRvsdas" << endl;
    cout << "                 :irrrriiirr77r:ii:" << YELLOW << ".i7i:::::::::::.." << WHITE << "rvuqdgrasdsa" << endl;
    cout << "                   1s..irrrrrrri::" << YELLOW << ".iri::   .riiii:." << WHITE << " vPKP7" << endl;
    cout << "                            .iiii:" << YELLOW << ".i7i:i::.rv7ri:."<< WHITE << "  rbPq." << endl;
    cout << "                           ..irri. " << YELLOW << ".riiiir7ri:.." << WHITE << "  .i5I7vsgg5r7:1sadsa" << endl;
    cout << "                      7XQQQRJv777i. " << YELLOW << "::.." << WHITE << "        .Xgb7rr7LvvL1IPZ." << endl;
    cout << "                   vBBRdsv777vYsu1jvi" << YELLOW << "....:...." << WHITE << "rKDZEPjvi  ......." << endl;
    cout << "                   :..vYrr7YJY7r: .i:7ZEDgQgvbEPL:.YQQBJ" << endl;
    cout << "                     UB2u7:       ::7dbbggi iMEXIq.   .." << endl;
    cout << "                                 2ivZDgEr    rQgZP5" << endl;
    cout << "                                S71MQ:         :iii" << endl;
    cout << "                                BBBj" << GREEN << endl;
    cout << "                                             ." << endl;
    cout << "                                             ." << endl;
    cout << "                                             ." << endl;
    cout << "                                             ." << endl;
    cout << "                                             ." << endl;
    cout << "                                             ." << endl;
    cout << "                                            .." << endl;
    cout << "                                            .   :i:." << endl;
    cout << "                                            : :ui" << endl;
    cout << "                                           .. r" << endl;
    cout << "                                            i" << endl;
    cout << "                                            i" << endl;
    cout << "                                           r." << endl;
    cout << "                                           7" << endl;
    cout << "                                          r." << endl;
    cout << "                                         .Y" << RESET << endl;
}


void usage(){
    printf("Usage:\n");
    printf("\tDBPD --io orthogroup_file --it species_tree_file [options]\n");
}

void print_help(){

    cout << endl;
    // 程序介绍部分
    cout << "DBPD(Dawn Blossoms Plucked at Duck) - a HGT gene detection tool" << endl;
    cout << "When use this program please confirm the leaf node names of the specie tree are correspond\n"
            "to the names in the orthogroups file" << endl;
    cout << endl;

    flower();

    // 大体用法部分
    usage();
    cout << endl;

    // 普通选项部分
    cout << "GENERAL OPTIONS:" << endl;
    cout << "-h, --help                  Print help information" << endl;
    cout << "--io FILE                   Input orthogroups file: xxx.proteinortho.tsv file provided by\n"
            "                            proteinortho or Orthogroups.tsv file provided by orthofinder2" << endl;
    cout << "--it FILE                   Input species tree file" << endl;
    cout << "--ip FILE                   Input the protein sequences of all genomes in one file, our\n"
            "                            program will produce every orthogroup\'s protein sequences in\n"
            "                            fasta format" << endl;
    cout << "--min_th NUM                The threshold for searching the distribution depth of every\n"
            "                            orthogroup on the tree, if you open the sliding mode of our\n"
            "                            program, this threshold will be used as the starting value in\n"
            "                            the sliding process, default is 0.6" << endl;
    cout << "-t, --threads NUM           Number of threads to use, default is 1.0" << endl;
    cout << "-p, --paint                 Output events scatter diagram, default is false" << endl;
    cout << "-v, --version               Display version number" << endl;
    cout << endl;

    //高级选项部分
    // 树置信度部分
    cout << "TREE CONFIDENCE OPTIONS:" << endl;
    cout << "-b, --bootstrap             The input phylogenetic tree contains bootstrap test, default\n"
            "                            is false" << endl;
    cout << "-s, --sh                    The input phylogenetic tree contains sh test, default is false" << endl;
    cout << endl;

    // 树时间标注部分
    cout << "DATING TREE OPTIONS:" << endl;
    cout << "--time_iqtree2              The input time tree is produced by iqtree2 or the same format" << endl;
    cout << "--time_treePL               The input time tree is produced by treePL or the same format" << endl;
    cout << endl;


    // 滑行阈值部分
    cout << "SLDING MODE OPTIONS:" << endl;
    cout << "--slding_mode               Open the sliding threshold mode" << endl;
    cout << endl;

    // DBPD时间、距离聚类部分
    cout << "EVENTS CLUSTERING OPTIONS:" << endl;
    cout << "-M, --Minpts NUM            Threshold number of neighbors for DBscan to confirm a cluster\n"
            "                            center, default is 4.0" << endl;
    cout << "-e, --eps NUM               The distance threshold for DBscan to find a neighbor, default\n"
            "                            is 100.0" << endl;
    cout << endl;


    // 尾部信息
    cout << "DBPD version v1.2.0 is developed by Zhang Shengbin in Shandong university, if you have\n"
            "any questions or advice, Please contact me by sending email to 202012632@mail.sdu.edu.cn" << endl;
    cout << endl;
}

struct option longopts[] = {
        {"io", required_argument, NULL, 1},
        {"it", required_argument, NULL, 2},
        {"ip", required_argument, NULL, 3},
        {"min_th",required_argument, NULL, 4},
        {"paint", required_argument, NULL, 'p'},
        {"eps", required_argument, NULL, 'e'},
        {"Minpts", required_argument, NULL, 'M'},
        {"threads", required_argument, NULL, 't'},
        {"bootstrap", no_argument, NULL, 'b'},
        {"sh", no_argument, NULL, 's'},
        {"sliding_mode", no_argument, NULL, 6},
        {"time_iqtree2", no_argument, NULL, 7},
        {"time_treePL", no_argument, NULL, 8},
        {"version", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        {0, 0, 0, 0}
};

int main(int argc, char *argv[]){
    int c;//这个专门用来返回 getopt 的结果

    // 下面是相应的参数和 getopt 结果分开定义，更清晰
    double eps = 100, Minpts = 4, cpu = 1.0, min_t = 0.6;
    const string version_name = "DBPD(Dawn Blossoms Plucked at Duck) v1.2.0";
    string tree_root = "WHOLE", error_input, orthofile, treefile, proteinfile;
    bool orthofinder = false, bootstrap = false, sh = false, sl_mode = false, time_iq = false, time_tPL = false, paint = false;
    map<string, string> prot_name_sequence_map;

    // 表示我们自定义输入参数的报错信息
    opterr = 0;
    // 下面是对输入的信息进行拆解
    while ((c = getopt_long(argc, argv, ":e:t:bpsvM:hc", longopts, NULL)) != -1){
        switch (c) {
            case 'v':
                cout << version_name << endl;
                return 0;
            case 'b':
                bootstrap = true;
                break;
            case 's':
                sh = true;
                break;
            case 'e':
                eps = atof(optarg);
                break;
            case 't':
                cpu = atoi(optarg);
                break;
            case 'M':
                Minpts = atoi(optarg);
                break;
            case 'h':
                print_help();
                return 0;
            case 'p':
                paint = true;
                break;
            case 1:
                orthofile = optarg;
                break;
            case 2:
                treefile = optarg;
                break;
            case 3:
                proteinfile = optarg;
                break;
            case 4:
                min_t = atof(optarg);
                break;
            case 6:
                sl_mode = true;
                break;
            case 7:
                time_iq = true;
                break;
            case 8:
                time_tPL = true;
                break;
            case '?':
                // 如果是短选项错误 optopt 里是短选项，如果是长选项其值是0
                if(optopt){
                    printf("Unrecognized option: '%c'\n", optopt);
                    cerr << "To learn more information: run 'DBPD -h'" << endl;
                    exit(255);
                }
                else{
                    // 从前面往后一直找，找到这个发现解析不来了 optind 中存储的是下一个要找的位置
                    printf("Unrecognized option: \"%s\"\n", argv[optind-1]);
                    cerr << "To learn more information: run 'DBPD -h'" << endl;
                    exit(255);
                }

            case ':':
                if (optopt) {
                    printf("Option: \"%s\" needs a parameter\n", argv[optind-1]);
                    cerr << "To learn more information: run 'DBPD -h'" << endl;
                    exit(255);
                }
                else {
                    printf("Option: '%c' needs a parameter\n", optopt);
                    cerr << "To learn more information: run 'DBPD -h'" << endl;
                    exit(255);
                }
        }
    }

    // 下面对拆解的信息进行处理
    // 先看软件按运行必须的两个文件是否都输入了，都有首先证明运行是不成问题的
    if(orthofile.empty() or treefile.empty()){
        usage();
        handle_error("The orthogroups file and species tree file are necessary, to show all available options: run 'DBPD -h'");

    }

    // 检查时间文件，两个开关不能同时打开
    if(time_iq and time_tPL){
        handle_error("The time tree can only be produced by one software, to learn more information: run 'DBPD -h'");
    }

    if(!time_iq and !time_tPL){
        handle_error("Our program can only work on time tree with all tips dating 0 mya (million year ago), to learn more information: run 'DBPD -h'");
    }

    // 滑动阈值打开后，max_t和min_t的大小不对
    if(sl_mode = true){
        if(min_t > 1){
            handle_error("min_t must be a value smaller than 1, to learn more information: run 'DBPD -h'");
        }
    }

    // sh 和 bootstrap 值一个都没有是不被我们许可的
    if(!sh and ! bootstrap){
        handle_error("Please make sure that there is at least one parameter to evaluate the confidence of the tree, to learn more information: run 'DBPD -h'");
    }


    cout << "======================================Welcome to use DBPD======================================" << endl;
    cout << version_name << " is a HGT detection tool." << endl;
    flower();


    // 下面是树文件的处理
    Parse_Tree tree(treefile, bootstrap, sh, time_iq);
    tree.read_tree();
    tree.parse_tree(tree._tree_text, tree_root);
    tree.judge_tree_root();
    auto nodes_num = tree.parent_child_map.size();
    tree.find_tree_leaf(tree_root);
    tree.make_tree_attribution_map(tree_root);
    if(!time_iq){
        tree.oat = tree.tree_evolution_map[tree.leaf_name_to_code_map[tree.tree_leaf_map[tree_root][0]]];
        cout << "The oldest ancestor is " << tree.oat << " million years ago." << endl;
    }
    auto leafs_num = tree.tree_leaf_map["WHOLE"].size();
    cout << "We detect " << nodes_num << " internal nodes in your species tree." << endl;
    cout << "We detect " << leafs_num << " leaves in your species tree." << endl;

    cout << "======================================Parse tree finish======================================" << endl;


    sliding_mode sl_matrix(tree, min_t);
    // 开启滑动模型
    if(sl_mode){
        cout << "You are using the sliding mode, we are carrying out preparations." << endl;
        sl_matrix.calculate_every_node_threshold();
        cout << "The threshold matrix has been finished." << endl;
    }

    // 下面是聚类文件的处理-水平转移推断
    OthoHGT ortho(orthofile ,tree, sl_mode, min_t, sl_matrix);
    ortho.read_file();
    cout << "======================================Read orthogroups finish======================================" << endl;
    cout << "Now we are working on the HGT infer." << endl;
    ortho.statistic_HGT_genes_paraller(cpu);
    cout << endl;
    cout << "======================================HGT infer finish======================================" << endl;


    // 对基因产生时间和HGT时间聚类
    DBscan dbscan(Minpts, eps, time_tPL, ortho, tree);
    dbscan.DBscan_clsuter();
    cout << endl;
    cout << "======================================DBscan cluster finish======================================" << endl;


    // 对结果进行输出
    Output out(ortho, dbscan, prot_name_sequence_map);
    out.out_event_result();
    out.out_HGT_percent_result();
    if(!proteinfile.empty()){
        string path = "";
        auto position = proteinfile.find_last_of("/");
        if(position != proteinfile.npos)
            path = proteinfile.substr(0,position+1);
        read_summary_protein_file(proteinfile, prot_name_sequence_map);
        out.out_ortho_protein_parallel(cpu, path);
    }
    cout << "======================================Output result finish======================================" << endl;

    //调用R语言脚本绘图
    if(paint) {
        system("./paint.R");
        cout << "======================================Paint finish======================================" << endl;
    }
}

