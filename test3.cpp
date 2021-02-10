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
int TDB_SIZE = 0;
using namespace std;

double phase_times[5] = {};

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

FPNode::FPNode(const Item& item, const std::shared_ptr<FPNode>& parent) :
    item( item ), frequency( 1 ), node_link( nullptr ), parent( parent ), children(), tid_list()
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
	double start = 0, end = 0;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //distribute the TDB 
    int i=0;
	get_walltime(&start);
    #pragma omp parallel default(shared) private(i) num_threads(max_threads)
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(dynamic)
        for(i=0;i<number_of_transactions;i++){
            const Transaction& transaction = transactions[i];
            for(const Item& item : transaction){
                partial_supports[thread_id][item]++;
                partial_tid_list[thread_id][item].insert(tids[i]);
                partial_items[thread_id].insert(item);
            }
        }
    }
	get_walltime(&end);
	phase_times[0] += (end - start);
    set<Item> set_items;
    for(i = 0; i < max_threads; i++) {
        set_items.insert(partial_items[i].begin(), partial_items[i].end());
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //compute global-supports, global tid-lists and filter with minSup
    map<Item, int> global_supports;
    map<Item, set<int>> global_tid_lists;
    vector<Item> items(set_items.begin(),set_items.end()); 
    set<Item> set_fitems;
    i=0;
	get_walltime(&start);
    #pragma omp parallel default(shared) private(i) num_threads(max_threads)
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
	get_walltime(&end);
	phase_times[1] += (end - start);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    //compute periodicity and filter with maxPer
    set<Item> fpitems;
    map<Item, int> frequency_by_item;
    vector<Item> fitems(set_fitems.begin(),set_fitems.end());
    i = 0;
	get_walltime(&start);
    #pragma omp parallel default(shared) private(i) num_threads(max_threads)
    {
        #pragma omp for schedule(dynamic,1)
        for(i = 0; i < fitems.size(); i++) {
            int periodicity = -1;
            int lasttid = 0;
            for(auto it = global_tid_lists[fitems[i]].begin(); it != global_tid_lists[fitems[i]].end(); it++) {
                periodicity = max(*it - lasttid, periodicity);
                lasttid = *it;
            }
            periodicity = max(periodicity, (int)TDB_SIZE - lasttid);
            if(periodicity <= maximum_periodicity) {
                #pragma omp critical
                {
                    fpitems.insert(fitems[i]);
                }
                
            }
        }
    }
    get_walltime(&end);
	phase_times[2] += (end - start);
	for(auto it : fpitems) {
        frequency_by_item[it] = global_supports[it];
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
    for(const auto& pair : items_ordered_by_frequency) {
        items_with_frequency.push_back(pair);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //construct partial PF-trees
    vector< std::shared_ptr<FPTree> > fptrees(max_threads, nullptr);
    for(i = 0; i < max_threads; i++) fptrees[i] = make_shared<FPTree>(minimum_support_threshold, maximum_periodicity);
    i = 0;
	get_walltime(&start);
    #pragma omp parallel default(shared) private(i) num_threads(max_threads)
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
	get_walltime(&end);
	phase_times[3] += (end - start);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //merging local PF-trees
    set<Item> vis;
    i = 0;
	get_walltime(&start);
    for(auto it:items_ordered_by_frequency){
        Item curr_item = it.first;
        shared_ptr<FPNode> last_ptr = nullptr;
        for(int j=0;j<max_threads;j++){
            if(fptrees[j]->header_table.find(curr_item)==fptrees[j]->header_table.end()){
                continue;
            }
            if(last_ptr != nullptr){
                last_ptr->node_link = fptrees[j]->header_table[curr_item];
                last_ptr = last_ptr->node_link;
            }
            else{
                last_ptr = fptrees[j]->header_table[curr_item];
                header_table[curr_item] = fptrees[j]->header_table[curr_item];
            }
            if(((last_ptr->parent.lock())->parent.lock())==nullptr){
                root->children.push_back(last_ptr);
                last_ptr->parent = root;
            }
            while(last_ptr->node_link){
                last_ptr = last_ptr->node_link;
                if(((last_ptr->parent.lock())->parent.lock())==nullptr){
                    root->children.push_back(last_ptr);
                    last_ptr->parent = root;
                }
            }
        }
    }
	get_walltime(&end);
	phase_times[4] += (end - start);
}

int main(int argc, char* argv[]){

    if(argc!=7){
        cout<<"Invalid Arguments";
        return 0;
    }

    ifstream fin;
    fin.open(argv[1]);
    
    const double min_sup_percentage = double(stof(argv[2]));
    const double max_per_percentage = double(stof(argv[3]));
	const int max_threads = stoi(argv[4]);
    const int max_runs_parallel = stoi(argv[5]);
    const int start_thread = stoi(argv[6]);

    int flag = 1;
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
        tids.push_back(i);
    }
    TDB_SIZE = len;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    int min_sup = int(0.01 * min_sup_percentage * len);
    int max_per = int(0.01 * max_per_percentage * len);
    cout << "TDB length: " << len << endl;
    cout << "minSup: " << min_sup_percentage << "%, " << min_sup << endl;
    cout << "maxPer: " << max_per_percentage << "%, " << max_per << endl;
    cout << "available threads: " << omp_get_max_threads() << endl; 

    for(int p = start_thread; p <= max_threads; p++) {
		for(int i = 0; i < 5; i++) {
			phase_times[i] = 0.0;
		}
		cout << "#Threads = " << p << endl;
		double start = 0.0, end = 0.0, total_time = 0.0;
		for(int run = 0; run < max_runs_parallel; run++) {
			start = omp_get_wtime();
			//////////////////////////////////////////////////////////////////////
			const FPTree parallelfptree{ tids, transactions, min_sup, max_per, p};
			//////////////////////////////////////////////////////////////////////
			end = omp_get_wtime();
			total_time += end - start;
		}
		total_time = total_time / max_runs_parallel;
		cout << "total-time = " << total_time << endl; 
		cout << "t = " <<  endl;		
		for(int i = 0; i < 5; i++) {
			printf(" %3.4lf ", phase_times[i]/max_runs_parallel);
		}
		cout << endl << "% = " << endl;
		for(int i = 0; i < 5; i++) {
			printf(" %3.4lf ", 100*phase_times[i]/(max_runs_parallel * total_time));
		}
		cout << endl;
	}

	return 0;
}