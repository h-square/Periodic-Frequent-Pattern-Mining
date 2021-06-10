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
#include <sys/time.h>
#include "fptree.hpp"
using namespace std;

int TDB_SIZE = 0;
double total_time_fptree = 0;
double phase_time[4]={0};

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

//constructor for item and it's parent
FPNode::FPNode(const Item& item, const std::shared_ptr<FPNode>& parent) :
    item( item ), frequency( 1 ), node_link( nullptr ), parent( parent ), children(), tid_list()
{
}

FPTree::FPTree(const vector<int> tids, const std::vector<Transaction>& transactions, const int minimum_support_threshold, const int maximum_periodicity) :
    root( std::make_shared<FPNode>( Item{}, nullptr ) ), header_table(), minimum_support_threshold( minimum_support_threshold ), maximum_periodicity(maximum_periodicity)
{
    
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
    double st,en;
    int i=0;
    get_walltime(&st);
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
    get_walltime(&en);
    phase_time[0]+=(en-st);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //compute global-supports, global tid-lists and filter with minSup
    
    map<Item, int> global_supports;
    map<Item, set<int>> global_tid_lists;
    vector<Item> items(set_items.begin(),set_items.end()); 
    set<Item> set_fitems;
    i=0;
    get_walltime(&st);
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
    get_walltime(&en);
    phase_time[1]+=(en-st);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    //compute periodicity and filter with maxPer
    set<Item> fpitems;
    map<Item, int> frequency_by_item;
    vector<Item> fitems(set_fitems.begin(),set_fitems.end());
    i = 0;
    get_walltime(&st);
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
    get_walltime(&en);
    phase_time[2]+=(en-st);
    for(auto it : fpitems) {
        frequency_by_item[it] = global_supports[it];
    }
    
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
    //construct partial PF-trees
    
    vector< std::shared_ptr<FPTree> > fptrees(max_threads, nullptr);
    for(i = 0; i < max_threads; i++) fptrees[i] = make_shared<FPTree>(minimum_support_threshold, maximum_periodicity);
    i = 0;
    get_walltime(&st);
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
    get_walltime(&en);
    phase_time[3]+=(en-st);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //merging local PF-trees
    set<Item> vis;
    vector<shared_ptr<FPNode>> partial_root[max_threads];
    i=0;
    #pragma omp parallel default(shared) private(i) num_threads(max_threads)
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(static)
        for(i=0;i<items_with_frequency.size();i++){
            Item curr_item = items_with_frequency[i].first;
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
                    #pragma omp critical
                    {
                        header_table[curr_item] = fptrees[j]->header_table[curr_item];
                    }
                }
                if(((last_ptr->parent.lock())->parent.lock())==nullptr){
                    partial_root[thread_id].push_back(last_ptr);
                    //root->children.push_back(last_ptr);
                    last_ptr->parent = root;
                }
                while(last_ptr->node_link){
                    last_ptr = last_ptr->node_link;
                    if(((last_ptr->parent.lock())->parent.lock())==nullptr){
                        partial_root[thread_id].push_back(last_ptr);
                        //root->children.push_back(last_ptr);
                        last_ptr->parent = root;
                    }
                }
            }
        }
    }
    for(int i=0;i<max_threads;i++){
        for(int j=0;j<partial_root[i].size();j++){
            root->children.push_back(partial_root[i][j]);
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

std::set<onlyPattern> parallel_fptree_growth(const FPTree& fptree, const int max_threads)
{
    if ( fptree.empty() ) { return {}; }

    if ( contains_single_path( fptree ) ) {
        // generate all possible combinations of the items in the tree

        std::set<onlyPattern> single_path_patterns;

        // for each node in the tree
        assert( fptree.root->children.size() == 1 );
        auto curr_fpnode = fptree.root->children.front();
        while ( curr_fpnode ) {
            const Item& curr_fpnode_item = curr_fpnode->item;
            const int curr_fpnode_frequency = curr_fpnode->frequency;
            set<int> curr_fpnode_tids = curr_fpnode->tid_list;

            // add a pattern formed only by the item of the current node
            onlyPattern new_pattern{ { curr_fpnode_item } };
            single_path_patterns.insert( new_pattern );

            // create a new pattern by adding the item of the current node to each pattern generated until now
            for ( const onlyPattern& onlypattern : single_path_patterns ) {
                onlyPattern new_pattern{ onlypattern };
                new_pattern.insert( curr_fpnode_item );
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

        std::set<onlyPattern> multi_path_patterns;

        // for each item in the FP-list 
        for (int i = fptree.items_with_frequency.size() - 1; i >= 0; i-- ) {
            const Item& curr_item = fptree.items_with_frequency[i].first;
            
            vector<TransformedPrefixPath> conditional_pattern_base;
            vector<TransformedPrefixPath> partial_conditional_pattern_base[max_threads];

            // for each path in the header_table (relative to the current item)
            auto ht = fptree.header_table;
            auto path_starting_fpnode = ht[curr_item];

            vector<shared_ptr<FPNode>> item_fpnodes;
            while(path_starting_fpnode){
                item_fpnodes.push_back(path_starting_fpnode);
                path_starting_fpnode = path_starting_fpnode->node_link;
            }

            int iterator=0;
            #pragma omp parallel default(shared) private(iterator)
            {
                int thread_id = omp_get_thread_num();
                #pragma omp for schedule(static)
                for(iterator=0;iterator<item_fpnodes.size();iterator++){
                    auto curr_path_starting_fpnode = item_fpnodes[iterator];    
                    set<int> path_starting_fpnode_tids = curr_path_starting_fpnode->tid_list;

                    auto curr_path_fpnode = curr_path_starting_fpnode->parent.lock();
                    if ( curr_path_fpnode->parent.lock() ) {
                        TransformedPrefixPath transformed_prefix_path{ {}, path_starting_fpnode_tids };

                        while ( curr_path_fpnode->parent.lock() ) {
                            transformed_prefix_path.first.push_back( curr_path_fpnode->item );
                            curr_path_fpnode = curr_path_fpnode->parent.lock();
                        }
                        partial_conditional_pattern_base[thread_id].push_back( transformed_prefix_path );
                    }
                }
            }

            for(int i=0;i<max_threads;i++){
                for(int j=0;j<partial_conditional_pattern_base[i].size();j++){
                    conditional_pattern_base.push_back(partial_conditional_pattern_base[i][j]);
                }
            }
            
            // generate the transactions that represent the conditional pattern base
            vector<Transaction> conditional_fptree_transactions;
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
            
            double start,end;
            get_walltime(&start);
            const FPTree conditional_fptree( conditional_fptree_tids, conditional_fptree_transactions, fptree.minimum_support_threshold, fptree.maximum_periodicity, max_threads);
            get_walltime(&end);
            total_time_fptree += ((double) (end - start)) ;
            // call recursively fptree_growth on the conditional fptree (empty fptree: no patterns)

            std::set<onlyPattern> conditional_patterns = parallel_fptree_growth( conditional_fptree, max_threads);

            // construct patterns relative to the current item using both the current item and the conditional patterns
            std::set<onlyPattern> curr_item_patterns;

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
            onlyPattern onlypattern{ {curr_item} };
            curr_item_patterns.insert( onlypattern );

            // the next patterns are generated by adding the current item to each conditional pattern
            for ( const onlyPattern& onlypattern : conditional_patterns ) {
                onlyPattern new_pattern{ onlypattern };
                new_pattern.insert( curr_item );
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

int main(int argc, char* argv[]){

    if(argc!=6){
        cout<<"Invalid Arguments";
        return 0;
    }

    ifstream fin;
    fin.open(argv[1]);
    
    const double min_sup_percentage = double(stof(argv[2]));
    const double max_per_percentage = double(stof(argv[3]));
    const int max_threads = stoi(argv[4]);
    const int start_thread = stoi(argv[5]);

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

    cout<<"Number of transactions: "<<len<<endl;

    int min_sup = int(0.01 * min_sup_percentage * len);
    int max_per = int(0.01 * max_per_percentage * len);

    cout<<"Minimum Support: "<<min_sup<<endl;
    cout<<"Maximum Periodicity: "<<max_per<<endl;
    
    //////////////////////////////////////////////////////
    //For graphs
    vector<int> numberOfThreads;
    vector<double> totalRunTime,treeConstructionTime,miningRunTime;
    //////////////////////////////////////////////////////
    
    for(int p=start_thread;p<=max_threads;p++){
        cout<<"For "<<p<<" threads"<<endl;
        double start,end,st,en;
        total_time_fptree=0;
        for(int i=0;i<4;i++){
            phase_time[i]=0;
        }
        get_walltime(&start);

        get_walltime(&st);
        const FPTree parallelfptree{ tids, transactions, min_sup, max_per, p};
        get_walltime(&en);
        total_time_fptree+=(en-st);
        const std::set<onlyPattern> patterns = parallel_fptree_growth( parallelfptree, p);
        get_walltime(&end);
        double total_time = ((double) (end - start));

        cout<<"Number of patterns "<<patterns.size()<<endl;

        cout<<"Total Time taken "<<total_time<<endl;
        totalRunTime.push_back(total_time);
        cout<<"Construction: "<<total_time_fptree<<" "<<"Mining: "<<total_time-total_time_fptree<<endl;;
        treeConstructionTime.push_back(total_time_fptree);
        miningRunTime.push_back(total_time-total_time_fptree);
        cout<<"All phase timings for "<<p<<" threads: [ ";
        for(int i=0;i<4;i++){
            cout<<phase_time[i]<<", ";
        }
        cout<<" ] "<<endl<<endl;
    }

    cout<<"For Graphs: "<<endl;
    cout<<"Number of threads: [";
    for(int t:numberOfThreads){
        cout<<t<<", ";
    }
    cout<<"]"<<endl;
    cout<<"Total timing: [";
    for(double t:totalRunTime){
        cout<<t<<", ";
    }
    cout<<"]"<<endl;
    cout<<"Constructor time: [";
    for(double t:treeConstructionTime){
        cout<<t<<", ";
    }
    cout<<"]"<<endl;
    cout<<"Mining Time: [";
    for(double t:treeConstructionTime){
        cout<<t<<", ";
    }
    cout<<"]"<<endl;
}