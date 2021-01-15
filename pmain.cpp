/*
    Author: Rag Patel, Husain Hirani
*/

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cstdint>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <fstream>
#include <sstream>
#include <omp.h>
#include<sys/time.h>

#include "fptree.hpp"

#define pb push_back

using namespace std;

FPNode::FPNode(const Item& item, const std::shared_ptr<FPNode>& parent) :
    item( item ), frequency( 1 ), node_link( nullptr ), parent( parent ), children(), tid_list()
{
}

FPTree::FPTree(const vector<int> tids, const std::vector<Transaction>& transactions, const int minimum_support_threshold, const int maximum_periodicity) :
    root( std::make_shared<FPNode>( Item{}, nullptr ) ), header_table(), minimum_support_threshold( minimum_support_threshold ), maximum_periodicity(maximum_periodicity)
{}

FPTree::FPTree(const int minimum_support_threshold, const int maximum_periodicity) :
    root( std::make_shared<FPNode>( Item{}, nullptr ) ), header_table(), minimum_support_threshold( minimum_support_threshold ), maximum_periodicity(maximum_periodicity) 
{}

FPTree::FPTree(const vector<int> tids,const std::vector<Transaction>& transactions, const int minimum_support_threshold, const int maximum_periodicity, const int max_threads) :
    root( std::make_shared<FPNode>( Item{}, nullptr ) ), header_table(), minimum_support_threshold( minimum_support_threshold ), maximum_periodicity(maximum_periodicity)
    {
    omp_set_num_threads(max_threads);
    int number_of_transactions = tids.size();
    vector<map<Item, int>> partial_supports(max_threads);
    vector<map<Item, set<int>>> partial_tid_list(max_threads);
    vector<set<Item>> partial_items(max_threads);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //distribute the TDB 
    int i=0;
    #pragma omp parallel default(shared) private(i)
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(static)
        for(i=0;i<number_of_transactions;i++){
            const Transaction& transaction = transactions[i];
            for(const Item& item : transaction){
                partial_supports[thread_id][item]++;
                partial_tid_list[thread_id][item].insert(tids[i]);
                partial_items[thread_id].insert(item);
            }
        }
    }
    set<Item> set_items;
    for(i = 0; i < max_threads; i++) {
        set_items.insert(partial_items[i].begin(), partial_items[i].end());
    }

    cout << "1-itemset\n";
    for(auto it : set_items) {
        cout << it << " ";
    }
    cout << endl;
    for(int j=0;j<partial_supports.size();j++){
         cout<<"tid: " << j << endl;
         for(auto it=partial_tid_list[j].begin();it!=partial_tid_list[j].end();it++){
             set<int> temp = it->second;
             cout<<it->first<<" : "<< partial_supports[j][it->first] << " : ";
             for(int te:temp){
                 cout<<te<<" ";
             }
             cout<<endl;
         }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //compute global-supports, global tid-lists and filter with minSup
    map<Item, int> global_supports;
    map<Item, set<int>> global_tid_lists;
    vector<Item> items(set_items.begin(),set_items.end()); 
    set<Item> set_fitems;
    i=0;
    #pragma omp parallel default(shared) private(i)
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(static,1) collapse(1)
        for(i=0;i<items.size();i++){
            int global_support = 0;
            set<int> global_tid_list;
            for(int j=0;j<max_threads;j++){
                if(partial_supports[j].find(items[i]) != partial_supports[j].end()) {
                    global_support = global_support + partial_supports[j][items[i]];
                    global_tid_list.insert(partial_tid_list[j][items[i]].begin(),partial_tid_list[j][items[i]].end());
                }
            }
            if(global_support >= minimum_support_threshold) {
                #pragma omp critical 
                {
                    set_fitems.insert(items[i]);
                    global_supports[items[i]] = global_support;
                    global_tid_lists[items[i]].insert(global_tid_list.begin(), global_tid_list.end()); 
                }
            }
        }
    }

    cout << "Frequent-Items: \n";
    for(auto it:set_fitems) {
        cout << it << " ";
    }
    cout << endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    //compute periodicity and filter with maxPer
    set<Item> fpitems;
    map<Item, int> frequency_by_item;
    vector<Item> fitems(set_fitems.begin(),set_fitems.end());
    i = 0;
    #pragma omp parallel default(shared) private(i)
    {
        #pragma omp for schedule(dynamic,1)
        for(i = 0; i < fitems.size(); i++) {
            int periodicity = -1;
            int lasttid = 0;
            for(auto it = global_tid_lists[fitems[i]].begin(); it != global_tid_lists[fitems[i]].end(); it++) {
                periodicity = max(*it - lasttid, periodicity);
                lasttid = *it;
            }
            periodicity = max(periodicity, (int)tids.size() - lasttid);
            if(periodicity <= maximum_periodicity) {
                #pragma omp critical
                {
                    fpitems.insert(fitems[i]);
                }
                
            }
        }
    }
    for(auto it : fpitems) {
        frequency_by_item[it] = global_supports[it];
    }
    cout <<"Periodic-Frequent Items: \n";
    for(auto it=frequency_by_item.begin();it!=frequency_by_item.end();it++){
        cout<<it->first<<" : ";
        set<int> temp = global_tid_lists[it->first];
        for(int te:temp){
            cout<<te<<" ";
        }
        cout<<endl;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //sort the frequent-items by frequency
    struct frequency_comparator
    {
        bool operator()(const std::pair<Item, uint64_t> &lhs, const std::pair<Item, uint64_t> &rhs) const
        {
            return std::tie(lhs.second, lhs.first) > std::tie(rhs.second, rhs.first);
        }
    };
    std::set<std::pair<Item, int>, frequency_comparator> items_ordered_by_frequency(frequency_by_item.cbegin(), frequency_by_item.cend());

    cout << "Ordered By Frequency: ";
    for(auto te:items_ordered_by_frequency) {
        cout << te.first << " ";
    }
    cout << endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //construct partial PF-trees
    vector< std::shared_ptr<FPTree> > fptrees(max_threads, nullptr);
    for(i = 0; i < max_threads; i++) fptrees[i] = make_shared<FPTree>(minimum_support_threshold, maximum_periodicity);
    i = 0;
    #pragma omp parallel default(shared) private(i)
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(static)
        for(i = 0; i < number_of_transactions; i++) {
            const Transaction& transaction = transactions[i];
            auto curr_fpnode = fptrees[thread_id]->root;
            auto& curr_header_table = fptrees[thread_id]->header_table;
            
            for ( const auto& pair : items_ordered_by_frequency ) {
                const Item& item = pair.first;
                //check if item is present in current transaction
                if ( std::find( transaction.cbegin(), transaction.cend(), item ) != transaction.cend() ) {
                    
                    const auto it = std::find_if(
                        curr_fpnode->children.cbegin(), curr_fpnode->children.cend(),  [item](const std::shared_ptr<FPNode>& fpnode) {
                            return fpnode->item == item; 
                    } );

                    if ( it == curr_fpnode->children.cend() ) {
                        // the child doesn't exist, create a new node
                        const auto curr_fpnode_new_child = std::make_shared<FPNode>( item, curr_fpnode );

                        // add the new node to the tree
                        curr_fpnode->children.push_back( curr_fpnode_new_child );

                        // update the node-link structure
                        if ( curr_header_table.count( curr_fpnode_new_child->item ) ) {
                            auto prev_fpnode = curr_header_table[curr_fpnode_new_child->item];
                            while ( prev_fpnode->node_link ) { prev_fpnode = prev_fpnode->node_link; }
                            prev_fpnode->node_link = curr_fpnode_new_child;
                        }
                        else {
                            curr_header_table[curr_fpnode_new_child->item] = curr_fpnode_new_child;
                        }

                        // advance to the next node of the current transaction
                        curr_fpnode = curr_fpnode_new_child;
                    }
                    else {
                    // the child exist, increment its frequency
                    auto curr_fpnode_child = *it;
                    ++curr_fpnode_child->frequency;

                    // advance to the next node of the current transaction
                    curr_fpnode = curr_fpnode_child;
                    }
                }            
            }

            if(curr_fpnode != fptrees[thread_id]->root)
                curr_fpnode->tid_list.insert(tids[i]);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //merging local PF-trees
    vector<Item> pfitems(items_ordered_by_frequency.size());
    for(auto it : items_ordered_by_frequency) {
        pfitems.push_back(it.first);
    }
    set<Item> vis;
    i = 0;
    
    #pragma omp parallel default(shared) private(i)
    {
        #pragma omp for schedule(dynamic)
        for(i = 0; i < pfitems.size(); i++) {
            Item curr_item = pfitems[i];
            shared_ptr<FPNode> last_ptr=nullptr;
            for(int j=0;j<max_threads;j++){
                if(fptrees[j]->header_table.find(curr_item)==fptrees[j]->header_table.end()){
                    continue;
                }
                if(last_ptr){
                    last_ptr->node_link = fptrees[j]->header_table[curr_item];
                    last_ptr = last_ptr->node_link;
                }
                else{
                    last_ptr = fptrees[j]->header_table[curr_item];
                    header_table[curr_item] = fptrees[j]->header_table[curr_item];
                }
                if(((last_ptr->parent.lock())->parent.lock())==nullptr && vis.find(last_ptr->item)==vis.end()){
                    vis.insert(last_ptr->item);
                    root->children.push_back(last_ptr);
                    last_ptr->parent = root;
                }
                while(last_ptr->node_link){
                    last_ptr = last_ptr->node_link;
                    if(((last_ptr->parent.lock())->parent.lock())==nullptr && vis.find(last_ptr->item)==vis.end()){
                        vis.insert(last_ptr->item);
                        root->children.push_back(last_ptr);
                        last_ptr->parent = root;
                    }
                }
            }    
        }
    }
}

void get_walltime_(double* wcTime) 
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	*wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);
}

void get_walltime(double* wcTime) 
{
	get_walltime_(wcTime);
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

    //const FPTree fptree{ tids, transactions, min_sup, max_per };
    const FPTree parallelfptree{ tids, transactions, min_sup, max_per, 2};
}