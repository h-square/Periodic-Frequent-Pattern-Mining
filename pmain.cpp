#include <cassert>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <utility>
#include <bits/stdc++.h>
#include <fstream>
#include <omp.h>

#include "fptree.hpp"

#define pb push_back

using namespace std;

FPNode::FPNode(const Item& item, const std::shared_ptr<FPNode>& parent) :
    item( item ), frequency( 1 ), node_link( nullptr ), parent( parent ), children()
{
}

FPTree::FPTree(const vector<int> tids,const std::vector<Transaction>& transactions, const int minimum_support_threshold, const int maximum_periodicity) :
    root( std::make_shared<FPNode>( Item{}, nullptr ) ), header_table(),
    minimum_support_threshold( minimum_support_threshold ),
    maximum_periodicity(maximum_periodicity)
{
    int max_threads = 2;
    omp_set_num_threads(max_threads);
    int number_of_transactions = tids.size();
    vector<map<Item, int>> partial_supports(max_threads);
    vector<map<Item, set<int>>> partial_tid_list(max_threads);
    set<Item> set_items;

    //distribute the TDB 
    int i=0;
    #pragma omp parallel default(shared) private(i)
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(static,1)
        for(i=0;i<number_of_transactions;i++){
            const Transaction& transaction = transactions[i];
            for(const Item& item : transaction){
                partial_supports[thread_id][item]++;
                partial_tid_list[thread_id][item].insert(tids[i]);
                set_items.insert(item);
            }
        }
    }
    /*
    for(int j=0;j<partial_supports.size();j++){
         cout<<j<<endl;
         for(auto it=partial_tid_list[j].begin();it!=partial_tid_list[j].end();it++){
             cout<<it->first<<" ";
             set<int> temp = it->second;
             for(int te:temp){
                 cout<<te<<" ";
             }
             cout<<endl;
         }
    }
    */

    map<Item, int> global_support;
    map<Item, set<int>> global_tid_list;
    vector<Item> items(set_items.begin(),set_items.end());
    set<Item> set_fitems;
    
    //compute partial-supports and partial tid-lists 
    //filter with minSup
    i=0;
    #pragma omp parallel default(shared) private(i)
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(static,1)
        for(i=0;i<items.size();i++){
            for(int j=0;j<max_threads;j++){
                if(partial_supports[j][items[i]]){
                    global_support[items[i]]+=partial_supports[j][items[i]];
                    global_tid_list[items[i]].insert(partial_tid_list[j][items[i]].begin(),partial_tid_list[j][items[i]].end());
                }
            }
            if(global_support[items[i]] >= minimum_support_threshold) {
                #pragma omp critical
                set_fitems.insert(items[i]);
            }
        }
    }

    map<Item, int> frequency_by_item;
    vector<Item> fitems(set_fitems.begin(),set_fitems.end());
        
    //compute period and filter with maxPer
    i = 0;
    #pragma omp parallel default(shared) private(i)
    {
        #pragma omp for schedule(static,1)
        for(i = 0; i < fitems.size(); i++) {
            int periodicity = -1;
            int lasttid = 0;
            for(auto it = global_tid_list[fitems[i]].begin(); it != global_tid_list[fitems[i]].end(); it++) {
                periodicity = max(*it - lasttid, periodicity);
                lasttid = *it;
            }
            periodicity = max(periodicity, (int)tids.size() - lasttid);
            if(periodicity <= maximum_periodicity) {
                frequency_by_item[fitems[i]] = global_support[fitems[i]];
            }
        }
    }
    
    /*
    for(auto it=frequency_by_item.begin();it!=frequency_by_item.end();it++){
        cout<<it->first<<" ";
        set<int> temp = global_tid_list[it->first];
        for(int te:temp){
            cout<<te<<" ";
        }
        cout<<endl;
    }
    */

    struct frequency_comparator
    {
        bool operator()(const std::pair<Item, uint64_t> &lhs, const std::pair<Item, uint64_t> &rhs) const
        {
            return std::tie(lhs.second, lhs.first) > std::tie(rhs.second, rhs.first);
        }
    };
    std::set<std::pair<Item, int>, frequency_comparator> items_ordered_by_frequency(frequency_by_item.cbegin(), frequency_by_item.cend());
}

int main(int argc, char* argv[]){

    if(argc!=4){
        cout<<"Invalid Arguments";
        return 0;
    }

    ifstream fin;
    fin.open(argv[1]);
    
    const int min_sup = stoi(argv[2]);
    const int max_per = stoi(argv[3]);

    vector<Transaction> transactions;
    int len=0;
    while(fin){
        string line;
        getline(fin,line);
        Transaction s;
        stringstream ss(line);
        while(ss>>line){
            s.push_back(line);
        }
        if(s.size())
            transactions.push_back(s);
        len++;
    }
    fin.close();

    if(!len){
        cout<<"Empty Database";
        return 0;
    }

    len-=1;
    vector<int> tids;
    for(int i=1;i<=len;i++){
        tids.pb(i);
    }

    // cout<<transactions.size()<<endl;
    // for(int i=0;i<transactions.size();i++){
    //     for(int j=0;j<transactions[i].size();j++){
    //         cout<<transactions[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }

    const FPTree fptree{ tids, transactions, min_sup, max_per };
}