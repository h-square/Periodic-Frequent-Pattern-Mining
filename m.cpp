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

FPNode::FPNode(const Item& item, const std::shared_ptr<FPNode>& parent) :
    item( item ), frequency( 1 ), node_link( nullptr ), parent( parent ), children(), tid_list()
{
}

FPTree::FPTree(const vector<int> tids, const std::vector<Transaction>& transactions, const int minimum_support_threshold, const int maximum_periodicity) :
    root( std::make_shared<FPNode>( Item{}, nullptr ) ), header_table(), minimum_support_threshold( minimum_support_threshold ), maximum_periodicity(maximum_periodicity)
{
     // scan the transactions counting the frequency of each item
    map<Item, int> frequency_by_item;
    map<Item, set<int>> tids_by_item;
    int i=0;
    for ( const Transaction& transaction : transactions ) {
        for ( const Item& item : transaction ) {
            ++frequency_by_item[item];
            tids_by_item[item].insert(tids[i]);
        }
        i++;
    }

    // keep only items which have a frequency greater or equal than the minimum support threshold
    for ( auto it = frequency_by_item.cbegin(); it != frequency_by_item.cend(); ) {
        const int item_frequency = (*it).second;
        if ( item_frequency < minimum_support_threshold ) {
            tids_by_item.erase(it->first);
            frequency_by_item.erase( it++ );
        }
        else { ++it; }
    }

    map<Item,int> period_by_item;
    for(auto it = tids_by_item.cbegin(); it!=tids_by_item.cend();it++){
        set<int> temp = it->second;
        int periodicity=-1;
        int lasttid = 0;
        for(auto itr = temp.begin(); itr != temp.end();itr++){
            periodicity = max(*itr - lasttid, periodicity);
            lasttid = *itr;
        }
        periodicity = max(periodicity, (int)TDB_SIZE-lasttid);
        period_by_item[it->first]=periodicity;
    }

    for ( auto it = period_by_item.cbegin(); it != period_by_item.cend(); ) {
        const int item_period = (*it).second;
        if ( item_period > maximum_periodicity ) {
            frequency_by_item.erase(it->first);
            tids_by_item.erase(it->first);
            period_by_item.erase( it++ );
        }
        else { ++it; }
    }

    // order items by decreasing frequency
    struct frequency_comparator
    {
        bool operator()(const std::pair<Item, uint64_t> &lhs, const std::pair<Item, uint64_t> &rhs) const
        {
            return std::tie(lhs.second, lhs.first) > std::tie(rhs.second, rhs.first);
        }
    };
    std::set<std::pair<Item, int>, frequency_comparator> items_ordered_by_frequency(frequency_by_item.cbegin(), frequency_by_item.cend());

    for(const auto& pair : items_ordered_by_frequency){
        const Item& item = pair.first;
        items_with_frequency.push_back(pair);
    }

    // start tree construction
	i=0;
    // scan the transactions again
    for ( const Transaction& transaction : transactions ) {
        auto curr_fpnode = root;
        // select and sort the frequent items in transaction according to the order of items_ordered_by_frequency
        for ( const auto& pair : items_ordered_by_frequency ) {
            const Item& item = pair.first;
            // check if item is contained in the current transaction
            if ( std::find( transaction.cbegin(), transaction.cend(), item ) != transaction.cend() ) {
                // insert item in the tree

                // check if curr_fpnode has a child curr_fpnode_child such that curr_fpnode_child.item = item
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
                    if ( header_table.count( curr_fpnode_new_child->item ) ) {
                        auto prev_fpnode = header_table[curr_fpnode_new_child->item];
                        while ( prev_fpnode->node_link ) { prev_fpnode = prev_fpnode->node_link; }
                        prev_fpnode->node_link = curr_fpnode_new_child;
                    }
                    else {
                        header_table[curr_fpnode_new_child->item] = curr_fpnode_new_child;
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
        if(curr_fpnode != root)
            curr_fpnode->tid_list.insert(tids[i]);
        i++;
    }
}

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

    // cout << "1-itemset\n";
    // for(auto it : set_items) {
    //     cout << it << " ";
    // }
    // cout << endl;
    // for(int j=0;j<partial_supports.size();j++){
    //      cout<<"tid: " << j << endl;
    //      for(auto it=partial_tid_list[j].begin();it!=partial_tid_list[j].end();it++){
    //          set<int> temp = it->second;
    //          cout<<it->first<<" : "<< partial_supports[j][it->first] << " : ";
    //          for(int te:temp){
    //              cout<<te<<" ";
    //          }
    //          cout<<endl;
    //      }
    // }

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

    // cout << "Frequent-Items: \n";
    // for(auto it:set_fitems) {
    //     cout << it << " ";
    // }
    // cout << endl;
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
            periodicity = max(periodicity, (int)TDB_SIZE - lasttid);
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
    // cout <<"Periodic-Frequent Items: \n";
    // for(auto it=frequency_by_item.begin();it!=frequency_by_item.end();it++){
    //     cout<<it->first<<" : ";
    //     set<int> temp = global_tid_lists[it->first];
    //     for(int te:temp){
    //         cout<<te<<" ";
    //     }
    //     cout<<endl;
    // }
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
    // cout << "Ordered By Frequency: ";
    // for(auto te:items_ordered_by_frequency) {
    //     cout << te.first << " ";
    // }
    // cout << endl;
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
    set<Item> vis;
    i = 0;
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
}

bool FPTree::empty() const
{
    assert( root );
    return root->children.size() == 0;
}

bool contains_single_path(const std::shared_ptr<FPNode>& fpnode)
{
    assert( fpnode );
    if ( fpnode->children.size() == 0 ) { return true; }
    if ( fpnode->children.size() > 1 ) { return false; }
    return contains_single_path( fpnode->children.front() );
}

bool contains_single_path(const FPTree& fptree)
{
    return fptree.empty() || contains_single_path( fptree.root );
}

std::set<Pattern> parallel_fptree_growth(const FPTree& fptree, const int max_threads)
{
    if ( fptree.empty() ) { return {}; }

    if ( contains_single_path( fptree ) ) {
        // generate all possible combinations of the items in the tree

        std::set<Pattern> single_path_patterns;

        // for each node in the tree
        assert( fptree.root->children.size() == 1 );
        auto curr_fpnode = fptree.root->children.front();
        while ( curr_fpnode ) {
            const Item& curr_fpnode_item = curr_fpnode->item;
            const int curr_fpnode_frequency = curr_fpnode->frequency;
            set<int> curr_fpnode_tids = curr_fpnode->tid_list;

            // add a pattern formed only by the item of the current node
            Pattern new_pattern{ { curr_fpnode_item }, curr_fpnode_tids };
            single_path_patterns.insert( new_pattern );

            // create a new pattern by adding the item of the current node to each pattern generated until now
            for ( const Pattern& pattern : single_path_patterns ) {
                Pattern new_pattern{ pattern };
                new_pattern.first.insert( curr_fpnode_item );
                //assert( curr_fpnode_frequency <= pattern.second );
                new_pattern.second = curr_fpnode_tids;

                single_path_patterns.insert( new_pattern );
            }

            // advance to the next node until the end of the tree
            assert( curr_fpnode->children.size() <= 1 );
            if ( curr_fpnode->children.size() == 1 ) { curr_fpnode = curr_fpnode->children.front(); }
            else { curr_fpnode = nullptr; }
        }
        return single_path_patterns;
    }
    else {
        // generate conditional fptrees for each different item in the fptree, then join the results

        std::set<Pattern> multi_path_patterns;

        // for each item in the FP-list 
        for (int i = fptree.items_with_frequency.size() - 1; i >= 0; i-- ) {
            const Item& curr_item = fptree.items_with_frequency[i].first;
            
            std::vector<TransformedPrefixPath> conditional_pattern_base;

            // for each path in the header_table (relative to the current item)
            auto ht = fptree.header_table;
            auto path_starting_fpnode = ht[curr_item];
    
            while ( path_starting_fpnode ){

                set<int> path_starting_fpnode_tids = path_starting_fpnode->tid_list;

                auto curr_path_fpnode = path_starting_fpnode->parent.lock();
                if ( curr_path_fpnode->parent.lock() ) {
                    TransformedPrefixPath transformed_prefix_path{ {}, path_starting_fpnode_tids };

                    while ( curr_path_fpnode->parent.lock() ) {
                        transformed_prefix_path.first.push_back( curr_path_fpnode->item );
                        curr_path_fpnode = curr_path_fpnode->parent.lock();
                    }
                    conditional_pattern_base.push_back( transformed_prefix_path );
                }

                path_starting_fpnode = path_starting_fpnode->node_link;
            }
            
            // generate the transactions that represent the conditional pattern base
            std::vector<Transaction> conditional_fptree_transactions;
            vector<int> conditional_fptree_tids;

            for ( const TransformedPrefixPath& transformed_prefix_path : conditional_pattern_base ) {
                const std::vector<Item>& transformed_prefix_path_items = transformed_prefix_path.first;
                //const uint64_t transformed_prefix_path_items_frequency = transformed_prefix_path.second;
                set<int> transformed_prefix_path_items_tids = transformed_prefix_path.second;

                Transaction transaction = transformed_prefix_path_items;

                // add the same transaction transformed_prefix_path_items_frequency times
                for ( auto it = transformed_prefix_path_items_tids.begin(); it != transformed_prefix_path_items_tids.end(); it++ ) {
                    conditional_fptree_tids.push_back(*it);
                    conditional_fptree_transactions.push_back( transaction );
                }
            }
            
            // int j=0;
            // for(auto transaction:conditional_fptree_transactions){
            //     for(Item it:transaction){
            //         cout<<it<<" ";
            //     }
            //     cout<<" Tids "<<conditional_fptree_tids[j++]<<endl;
            // }
            // break;
            
            
            // build the conditional fptree relative to the current item with the transactions just generated
            const FPTree conditional_fptree( conditional_fptree_tids, conditional_fptree_transactions, fptree.minimum_support_threshold, fptree.maximum_periodicity, max_threads);
            // call recursively fptree_growth on the conditional fptree (empty fptree: no patterns)

            std::set<Pattern> conditional_patterns = parallel_fptree_growth( conditional_fptree, max_threads);

            // construct patterns relative to the current item using both the current item and the conditional patterns
            std::set<Pattern> curr_item_patterns;

            // the first pattern is made only by the current item
            // compute the frequency of this pattern by summing the frequency of the nodes which have the same item (follow the node links)
            int curr_item_frequency = 0;
            set<int> curr_item_tids;
            auto fpnode = ht[curr_item];
            while ( fpnode ) {
                curr_item_frequency += fpnode->frequency;
                curr_item_tids.insert(fpnode->tid_list.begin(),fpnode->tid_list.end());
                fpnode = fpnode->node_link;
            }
            // add the pattern as a result
            Pattern pattern{ { curr_item }, curr_item_tids };
            curr_item_patterns.insert( pattern );

            // the next patterns are generated by adding the current item to each conditional pattern
            for ( const Pattern& pattern : conditional_patterns ) {
                Pattern new_pattern{ pattern };
                new_pattern.first.insert( curr_item );
                //assert( curr_item_frequency >= pattern.second );
                new_pattern.second = pattern.second;

                curr_item_patterns.insert( { new_pattern } );
            }

            // join the patterns generated by the current item with all the other items of the fptree
            multi_path_patterns.insert( curr_item_patterns.cbegin(), curr_item_patterns.cend() );

            auto leaf_fpnode = ht[curr_item];
            while(leaf_fpnode){
                set<int> leaf_fpnode_tids = leaf_fpnode->tid_list;
                auto parent = leaf_fpnode->parent.lock();
                if(parent != nullptr){
                    parent->tid_list.insert(leaf_fpnode_tids.begin(), leaf_fpnode_tids.end());
                }
                
                leaf_fpnode = leaf_fpnode->node_link;

                for(int j=0;j<parent->children.size();j++){
                    if(parent->children[j]->item==curr_item){
                        parent->children.erase(parent->children.begin()+j);
                        break;
                    }
                }
            }   
        }
        return multi_path_patterns;
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

    if(argc!=5){
        cout<<"Invalid Arguments";
        return 0;
    }

    ifstream fin;
    fin.open(argv[1]);
    
    const int min_sup = stoi(argv[2]);
    const int max_per = stoi(argv[3]);
	const int max_thds = stoi(argv[4]);
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

    //const FPTree fptree{ tids, transactions, min_sup, max_per };
    const FPTree parallelfptree{ tids, transactions, min_sup, max_per, max_thds};
    const std::set<Pattern> patterns = parallel_fptree_growth( parallelfptree, max_thds);
    cout << patterns.size();
}