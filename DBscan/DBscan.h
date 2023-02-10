// Created by zhangshengbin on 2021/11/26.
#ifndef _DB_H
#define _DB_H


#include <numeric>
#include "ortho_class.h"
#include <set>
#include <unordered_set>

using namespace std;
typedef vector<int> pl;

class DBscan {

public:

    DBscan(double Minpts, double eps, bool time_t, OthoHGT &obj_ortho, Parse_Tree &obj_tree):_Minpts(Minpts), _eps(eps), _eps_pow(eps*eps), _time_t(time_t), _obj_ortho(obj_ortho), _obj_tree(obj_tree){}

    // DBscan 的核心就在这里了
    void DBscan_clsuter();

    ~DBscan(){
        for(auto &i:_kd_value_map){
            delete i.first;
        }
    }

    vector<vector<double> > p_h_time_vec;
    vector<string> event_ortho_vec;
    vector<double> event_confidence_vec;
    vector<int> event_cluster_vec;

private:
    vector<int> _event_pos_vec;

    unordered_map<int, int> _point_present_num_map;
    map<vector<double>, int> _point_present_pos_map;
    unordered_map<int, int> _pos_cluster_map;

    // 一个父平面下面可能有1-2个子平面的
    map<pl*, vector<pl*> > _kd_tree_map;
    // 这个都是一一对应的关系
    map<pl*, pl* > _kd_tree_reverse_map;
    map<pl*, vector<double>>_kd_value_map;
    pl* _root_plane;

    map<int,vector<int>> _cluster_map;
    unordered_set<int> _cluster_finish_set;
    vector<int> _core_ortho_vec;

    OthoHGT &_obj_ortho;
    Parse_Tree &_obj_tree;

    bool _time_t;
    double _Minpts, _eps, _eps_pow;
    int _core_judge_num = 0;

    // 这个函数和 DBscan 其实没关系，只是构建 kd 树的时候可以十分方便给出各个事件的置信度
    inline double _give_point_confidence(vector<string> &tree_vec);

    // 构建 kd 树
    inline void _make_kd_tree_map();
    // 配合上一个函数进行 kd 树的构建，把基因组分开
    void _classification(vector<int> &next_vec_index, double median, vector<int> &father_plane, string type = "p_t", double hl = 0);

    // 对代查找的的数据在 kd 树上进行分类，找到其最邻近的叶子节点
    pl* _find_leaf(int event, pl* every_plane);
    // 代查找的数据在 kd 树上进行检索，找到满足条件的所有节点
    void _kd_tree_search(int event, pl* now_plane, vector<int> &neighbor_list, vector<pl*> &inspect_plane_vec, bool core, bool fflag = true);
    // 检查兄弟节点是否满足条件，配合上一个函数完成 kd 树搜索
    void _inspect_borther_node(int event, pl* father_plane, vector<int> &neighbor_list, vector<pl*> &inspect_plane_vec, bool core);

    // 寻找一个核心所有能连接的点，配合上一个函数完成 Dbscan 算法
    void _find_all_neighor_of_one_core(int event, int cluster_num, vector<int> &judge_vec, double &process, const double total);

};

inline double DBscan::_give_point_confidence(vector<string> &tree_vec){
    double SUM = 0;
    for(auto &tree:tree_vec){
        double test = _obj_tree.tree_test_map[tree];
        SUM += test;
    }
    double MEAN = SUM/tree_vec.size();
    return MEAN;
}

// 这个注意求的是叶子节点到目标节点的距离
// 文件是不同的时间格式处理是不一样的
inline void DBscan::_make_kd_tree_map(){
    vector <double> temp_present_vec;
    int pos = 0;
    for(auto &ortho:_obj_ortho.HGTortho_ptree){
        // 所有叶子节点时间都是0，所以我们随便挑一个就可以
        double p_t, h_t;
        // 对出现时间进行处理
        string arbitrarily_p_n = ortho.second[0];
        if(_time_t) {
            p_t = _obj_tree.oat - _obj_tree.tree_evolution_map[arbitrarily_p_n];
            if (p_t < 0.0001)
                p_t = 0;
        }
        else
            p_t = _obj_tree.search_iqtree2_time(arbitrarily_p_n);

        // 对水平转移时间进行处理
        for(auto &receive_tree:_obj_ortho.HGTortho_rtree[ortho.first]){
            if(_time_t) {
                h_t = _obj_tree.oat - _obj_tree.tree_evolution_map[receive_tree];
                if(h_t < 0.0001)
                    h_t = 0;
            }
            else
                h_t = _obj_tree.search_iqtree2_time(receive_tree);
            vector<string> temp_h_t;
            temp_h_t.push_back(receive_tree);
            vector<double> temp{p_t, h_t};


            // 如果说有两个节点，本来是一致的，但是因为减法的计算精度问题这里认为不一致，导致的后果是
            // kd树上多一个点，计算距离的时候它需要多算一下，距离计算因为阈值很高，肯定不受精度的影响
            // 总而言之，精度对于DBscan聚类来说没有影响
            // 另外就是输出的时候这个点和本来与它一样的点坐标上可能差一点点，但是这个对于后续都没有影响，所以无所谓了
            auto target = _point_present_pos_map.find(temp);
            if (target == _point_present_pos_map.end()) {
                _event_pos_vec.push_back(pos);
                _point_present_num_map[pos] = 1;
                _point_present_pos_map[temp] = pos;
                temp_present_vec.push_back(p_t);
            }
            else
                _point_present_num_map[target->second]++;
            p_h_time_vec.push_back(temp);
            event_ortho_vec.push_back(ortho.first);
            double confi = (_give_point_confidence(temp_h_t) + _give_point_confidence(ortho.second))/2;
            event_confidence_vec.push_back(confi);
            pos++;
        }
    }

    // 处理一下没有方向的水平基因转移
    // 按照事件的思想，no_direction_HGT_vec里的每一个min_key其实都只能算是一个事件，这一事件的特点就是 p 和 h 事件一致
    // 而且不区分哪颗树代表 p，哪棵树代表 h，总的置信度也是大家求平均就好了
    for(auto &no_direction_eve:_obj_ortho.no_direction_HGT_vec){
        string random_tree = no_direction_eve[0];
        double time;

        if(_time_t) {
            time = _obj_tree.oat - _obj_tree.tree_evolution_map[random_tree];
            if(time < 0.0001)
                time = 0;
        }
        else
            time = _obj_tree.search_iqtree2_time(random_tree);


        vector<double> temp{time, time};

        auto target = _point_present_pos_map.find(temp);

        if(target == _point_present_pos_map.end()) {
            _event_pos_vec.push_back(pos);
            _point_present_num_map[pos] = 1;
            _point_present_pos_map[temp] = pos;
            temp_present_vec.push_back(time);
        }
        else
            _point_present_num_map[target->second]++;

        p_h_time_vec.push_back(temp);
        // 加入的随机5字符就是为了这里能够对应起来
        event_ortho_vec.push_back(_obj_ortho.no_direction_ortho_map[no_direction_eve]);
        // 这里把最后加的字符标志删掉，后面要对里面的树进行遍历，这个要碍事了
        no_direction_eve.pop_back();
        double confi = _give_point_confidence(no_direction_eve);
        event_confidence_vec.push_back(confi);
        pos++;
    }

    // 这里用 0 代表第一次，既不是高的情况，也不是低的情况
    // 之后用 1 代表大，-1 代表小

    double median_p_t = get_median(temp_present_vec);
    vector<int> temp;
    _classification(_event_pos_vec, median_p_t, temp, "p_t", 0);
}


#endif
