// Created by zhangshengbin on 2021/11/26.
#include "DBscan.h"
#include <algorithm>

void DBscan::_classification(vector<int> &next_vec_index, double median, vector<int> &father_plane, string type, double hl){
    vector<double> h_vec, l_vec;
    vector<int> next_h_v, next_l_v;
    pl* classification_plane= new pl;

    classification_plane->push_back(hl);

    // 对所有的 HGT 事件进行循环分发
    for(auto &event:next_vec_index){
        if(type == "p_t") {
            // 这个地方其实不怎么存在浮点数做差的问题，因为media值
            if (abs(p_h_time_vec[event][0] - median) < 0.0001) {
                classification_plane->push_back(event);
            }else if (p_h_time_vec[event][0] < median) {
                l_vec.push_back(p_h_time_vec[event][1]);
                next_l_v.push_back(event);
            }else if(p_h_time_vec[event][0] > median) {
                h_vec.push_back(p_h_time_vec[event][1]);
                next_h_v.push_back(event);
            }
        }else {
            if (abs(p_h_time_vec[event][1] - median) < 0.0001) {
                classification_plane->push_back(event);
            }else if (p_h_time_vec[event][1] < median) {
                l_vec.push_back(p_h_time_vec[event][0]);
                next_l_v.push_back(event);
            }else if(p_h_time_vec[event][1] > median) {
                h_vec.push_back(p_h_time_vec[event][0]);
                next_h_v.push_back(event);
            }
        }
    }
    if(type == "p_t"){
        // 0 代表 p_t，1代表h_t
        if(hl == 0)
            _root_plane = classification_plane;
        vector<double> temp{median, 0};
        _kd_value_map[classification_plane] = temp;
        if(hl != 0){
            _kd_tree_map[&father_plane].push_back(classification_plane);
            _kd_tree_reverse_map[classification_plane] = &father_plane;
        }
        double h_vec_m = get_median(h_vec);
        double l_vec_m = get_median(l_vec);
        if(!next_h_v.empty())
            _classification(next_h_v, h_vec_m, *classification_plane, "h_t", 1);
        if(!next_l_v.empty())
            _classification(next_l_v, l_vec_m, *classification_plane, "h_t", -1);
    }
    else{
        // 0 代表 p_t，1代表h_t
        vector<double> temp{median, 1};
        _kd_value_map[classification_plane] = temp;
        _kd_tree_map[&father_plane].push_back(classification_plane);
        _kd_tree_reverse_map[classification_plane] = &father_plane;
        double h_vec_m = get_median(h_vec);
        double l_vec_m = get_median(l_vec);
        if(!next_h_v.empty())
            _classification(next_h_v, h_vec_m, *classification_plane, "p_t", 1);
        if(!next_l_v.empty())
            _classification(next_l_v, l_vec_m, *classification_plane, "p_t", -1);
    }
}

// cluster_finish_vec 是已被分类的 event 的集合
// inspect_plane_vec 里面是对每一次search来说已经检查过的平面, 在改变过之后的位置添加，因为每一个核心都会重置，所以这样就可以了
void DBscan::_kd_tree_search(int event, pl* now_plane, vector<int> &neighbor_list, vector<pl*> &inspect_plane_vec, bool core, bool fflag){

    // 寻找聚类核的时候满足条件就可以提前撤退了
    if(_core_judge_num > _Minpts)
        return;

    // 对于每一个聚类来说，只要是检查过了的平面，不再重复检查，inspect_dic 针对每一个聚类次次清空
    inspect_plane_vec.push_back(now_plane);
    // 这里就不仅仅是查父节点了，父节点查完，兄弟节点一查到底，然后返回来继续向根节点进发，查到根节点了就结束了，如果需要根节点的兄弟节点会在结束之前被全查一遍
    // 这里也分情况讨论了，具体讨论如下：
    // 这个节点是否是刚刚开始的节点，也就是叶节点，如果是叶子节点不存在检查相邻空间的问题
    // 这个节点是x轴节点还是y轴节点
    // 这个节点是不是根节点
    // 这个节点是否 “被” 过线了
    if(_kd_value_map[now_plane][1] == 0){
        if(abs(p_h_time_vec[event][0] - _kd_value_map[now_plane][0]) <= _eps){
            // 这个地方要来计算最终的距离了
            // 还没被分类，这个点(是点不是面)就查，分类了就不查了，继续看下一个点，这一层上的点查完了，按规矩取下一层查
            for(auto pt_now_event = (now_plane->begin() + 1); pt_now_event != now_plane->end(); pt_now_event++){
                if(_cluster_finish_set.find(*pt_now_event) != _cluster_finish_set.end()){
                    swap(*pt_now_event, now_plane->back());
                    now_plane->pop_back();
                    pt_now_event--;
                    continue;
                }
                double temp_d1 = p_h_time_vec[event][0]-p_h_time_vec[*pt_now_event][0], temp_d2 = p_h_time_vec[event][1]-p_h_time_vec[*pt_now_event][1], distance = temp_d1*temp_d1 + temp_d2*temp_d2;
                if(distance < _eps_pow){
                    neighbor_list.push_back(*pt_now_event);
                    if(core)
                        _core_judge_num += _point_present_num_map[*pt_now_event];
                    else{
                        swap(*pt_now_event, now_plane->back());
                        now_plane->pop_back();
                        pt_now_event--;
                    }
                }
            }

            if(fflag or (_kd_tree_map[now_plane].size() == 1)) {
                if (now_plane != _root_plane)
                    _kd_tree_search(event, _kd_tree_reverse_map[now_plane], neighbor_list, inspect_plane_vec, core,false);
            }else{
                _inspect_borther_node(event, now_plane, neighbor_list, inspect_plane_vec, core);
                if(now_plane != _root_plane)
                    _kd_tree_search(event, _kd_tree_reverse_map[now_plane], neighbor_list, inspect_plane_vec, core, false);
            }
        }else
            if(now_plane != _root_plane)
                _kd_tree_search(event, _kd_tree_reverse_map[now_plane], neighbor_list, inspect_plane_vec, core, false);

    }else{
        if(abs(p_h_time_vec[event][1] - _kd_value_map[now_plane][0]) <= _eps){
            // 这个地方要来计算最终的距离了
            // 还没被分类，这个点(是点不是面)就查，分类了就不查了，继续看下一个点，这一层上的点查完了，按规矩取下一层查
            for(auto pt_now_event = now_plane->begin() + 1; pt_now_event != now_plane->end(); pt_now_event++){
                if(_cluster_finish_set.find(*pt_now_event) != _cluster_finish_set.end()){
                    swap(*pt_now_event, now_plane->back());
                    now_plane->pop_back();
                    pt_now_event--;
                    continue;
                }
                double temp_d1 = p_h_time_vec[event][0]-p_h_time_vec[*pt_now_event][0], temp_d2 = p_h_time_vec[event][1]-p_h_time_vec[*pt_now_event][1], distance = temp_d1*temp_d1 + temp_d2*temp_d2;
                if(distance < _eps_pow){
                    neighbor_list.push_back(*pt_now_event);
                    if(core)
                        _core_judge_num += _point_present_num_map[*pt_now_event];
                    else{
                        swap(*pt_now_event, now_plane->back());
                        now_plane->pop_back();
                        pt_now_event--;
                    }
                }
            }
            if(fflag or (_kd_tree_map[now_plane].size() == 1)) {
                if (now_plane != _root_plane)
                    _kd_tree_search(event, _kd_tree_reverse_map[now_plane], neighbor_list, inspect_plane_vec, core,false);
            }else{
                _inspect_borther_node(event, now_plane, neighbor_list, inspect_plane_vec, core);
                if(now_plane != _root_plane)
                    _kd_tree_search(event, _kd_tree_reverse_map[now_plane], neighbor_list, inspect_plane_vec, core, false);
            }
        }else
            if(now_plane != _root_plane)
                    _kd_tree_search(event, _kd_tree_reverse_map[now_plane], neighbor_list, inspect_plane_vec, core, false);
    }
}

void DBscan::_inspect_borther_node(int event, pl* father_plane, vector<int> &neighbor_list, vector<pl*> &inspect_plane_vec, bool core) {
    // 寻找聚类核的时候满足条件提前结束
    if(_core_judge_num > _Minpts)
        return;

    if (_kd_tree_map.count(father_plane) != 0) {
        for (auto &child_p:_kd_tree_map[father_plane]) {
            if (find(inspect_plane_vec.begin(), inspect_plane_vec.end(), child_p) == inspect_plane_vec.end()) {
                if (_kd_value_map[child_p][1] == 0) {
                    if (abs(p_h_time_vec[event][0] - _kd_value_map[child_p][0]) <= _eps) {
                        for (auto now_event = child_p->begin() + 1; now_event != child_p->end(); now_event++) {
                            if(_cluster_finish_set.find(*now_event) != _cluster_finish_set.end()){
                                swap(*now_event, child_p->back());
                                child_p->pop_back();
                                now_event --;
                            }
                            double temp_d1 = p_h_time_vec[event][0]-p_h_time_vec[*now_event][0], temp_d2 = p_h_time_vec[event][1]-p_h_time_vec[*now_event][1], distance = temp_d1*temp_d1 + temp_d2*temp_d2;
                            if (distance < _eps_pow) {
                                neighbor_list.push_back(*now_event);
                                if(core)
                                    _core_judge_num += _point_present_num_map[*now_event];
                                else{
                                    swap(*now_event, child_p->back());
                                    child_p->pop_back();
                                    now_event --;
                                }
                            }
                        }
                    }
                } else {
                    if (abs(p_h_time_vec[event][1] - _kd_value_map[child_p][0]) <= _eps) {
                        for (auto now_event = child_p->begin() + 1; now_event != child_p->end(); now_event++) {
                            if(_cluster_finish_set.find(*now_event) != _cluster_finish_set.end()){
                                swap(*now_event, child_p->back());
                                child_p->pop_back();
                                now_event --;
                            }
                            double temp_d1 = p_h_time_vec[event][0]-p_h_time_vec[*now_event][0], temp_d2 = p_h_time_vec[event][1]-p_h_time_vec[*now_event][1], distance = temp_d1*temp_d1 + temp_d2*temp_d2;
                            if (distance < _eps_pow) {
                                neighbor_list.push_back(*now_event);
                                if(core)
                                    _core_judge_num += _point_present_num_map[*now_event];
                                else{
                                    swap(*now_event, child_p->back());
                                    child_p->pop_back();
                                    now_event --;
                                }
                            }
                        }
                    }
                }
                _inspect_borther_node(event, child_p, neighbor_list, inspect_plane_vec, core);
            }
        }
    }
}

pl* DBscan::_find_leaf(int event, pl* every_plane){
    // 这个代表是叶子节点，可以返回了
    if(!_kd_tree_map.count(every_plane))
        return every_plane;
    // 不是叶子节点，那么继续往下找
    vector<pl*> child_plane_vec = _kd_tree_map[every_plane];

    pl* temp_h_eve = nullptr;
    pl* temp_l_eve = nullptr;

    for(auto child_plane:child_plane_vec){
       if(child_plane->front() == 1)
           temp_h_eve = child_plane;
       else
           temp_l_eve = child_plane;
    }

    if(_kd_value_map[every_plane][1] == 0){
        if(p_h_time_vec[event][0] >= _kd_value_map[every_plane][0]){
            if(temp_h_eve != nullptr)
                return _find_leaf(event, temp_h_eve);
            else
                return temp_l_eve;
        }
        else {
            if (temp_l_eve != nullptr)
                return _find_leaf(event, temp_l_eve);
            else
                return temp_h_eve;
        }
    }
    else {
        if (p_h_time_vec[event][1] >= _kd_value_map[every_plane][0]) {
            if (temp_h_eve != nullptr)
                return _find_leaf(event, temp_h_eve);
            else
                return temp_l_eve;
        }
        else {
            if (temp_l_eve != nullptr)
                return _find_leaf(event, temp_l_eve);
            else
                return temp_h_eve;
        }
    }
}

void DBscan::DBscan_clsuter(){

    _make_kd_tree_map();
    cout << "Make kd tree finish." << endl;

    // 寻找到聚类的核心点
    srand((unsigned)time(NULL));

    for(auto &eve_event:_event_pos_vec){
        pl* search_start_leaf = _find_leaf(eve_event, _root_plane);
        vector<int> neighbor_list;
        vector<pl*> inspect_plane_vec;

        _kd_tree_search(eve_event, search_start_leaf, neighbor_list, inspect_plane_vec, true);
        if(_core_judge_num > _Minpts)
            _core_ortho_vec.push_back(eve_event);
        _core_judge_num = 0;
    }
    cout << "Find cluster core finish." << endl;
    cout << "Now we are working on the DBscan clustering." << endl;

    int cluster_num = 1;
    double process = 0, total = _core_ortho_vec.size();
    vector<int> judge_vec = _core_ortho_vec;
    cout << "We find "<< total << " cluster core." << endl;
    while (!_core_ortho_vec.empty()){
        process ++;
        cout << '\r' << "now: " << process/total;
        int random = rand() % _core_ortho_vec.size();
        auto event = _core_ortho_vec[random];
        _cluster_map[cluster_num].push_back(event);

        // 一开始的这个核心点在树上也是存在的，这里要屏蔽一下它自己，防止又在邻居里出现，会导致二次删除
        // 在树上遇到直接给它删了就可以了
        // 树上其它的点互相之间是没有重复的，所以不需要添加到_cluster_finish_set里
        // 如果是某个点邻域点的直接被删除了，不是的本来就应该留着

        _cluster_finish_set.insert(event);

        swap(*(_core_ortho_vec.begin() + random), _core_ortho_vec.back());
        _core_ortho_vec.pop_back();

        swap(*find(_event_pos_vec.begin(), _event_pos_vec.end(), event), _event_pos_vec.back());
        _event_pos_vec.pop_back();

        _find_all_neighor_of_one_core(event, cluster_num, judge_vec, process, total);
        cluster_num++;
    }

    for(auto &event:_event_pos_vec)
        _cluster_map[0].push_back(event);

    for(auto &event:_cluster_map){
        for(auto &pos:event.second){
            _pos_cluster_map[pos] = event.first;
        }
    }

    for(auto &event:p_h_time_vec){
        auto pos = _point_present_pos_map[event];
        event_cluster_vec.push_back(_pos_cluster_map[pos]);
    }
}


void DBscan::_find_all_neighor_of_one_core(int event, int cluster_num, vector<int> &judge_vec, double &process, const double total){
    pl* search_start_leaf = _find_leaf(event, _root_plane);
    vector<int> neighbor_list;
    vector<pl*> inspect_plane_vec;
    _kd_tree_search(event, search_start_leaf, neighbor_list, inspect_plane_vec, false);

    if(!neighbor_list.empty())
    {
        for(auto &other_event:neighbor_list){
            _cluster_map[cluster_num].push_back(other_event);
            swap(*find(_event_pos_vec.begin(), _event_pos_vec.end(), other_event), _event_pos_vec.back());
            _event_pos_vec.pop_back();
        }
        for(auto &another_event:neighbor_list){
            auto iter = find(judge_vec.begin(), judge_vec.end(), another_event);
            if(iter != judge_vec.end()){

                swap(*find(_core_ortho_vec.begin(), _core_ortho_vec.end(), another_event), _core_ortho_vec.back());
                _core_ortho_vec.pop_back();

                process ++;
                cout << '\r' << "now: " << process/total;
                _find_all_neighor_of_one_core(another_event, cluster_num, judge_vec, process, total);
            }
        }
    }
}