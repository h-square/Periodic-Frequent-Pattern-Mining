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

#include "fptree.hpp"

#define pb push_back
int TDB_SIZE = 0;
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
        periodicity = max(periodicity, (int)tids.size()-lasttid);
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

    // for(auto it=frequency_by_item.begin();it!=frequency_by_item.end();it++){
    //     cout<<it->first<<" "<<frequency_by_item[it->first]<<" "<<period_by_item[it->first]<<endl;
    // }

    for(const auto& pair : items_ordered_by_frequency){
        const Item& item = pair.first;
        items_with_frequency.pb(pair);
        //cout<<pair.first<<" "<<pair.second<<" "<<period_by_item[pair.first]<<endl;
    }

    // start tree construction

    i=0;
    // scan the transactions again
    for ( const Transaction& transaction : transactions ) {
        auto curr_fpnode = root;
        // select and sort the frequent items in transaction according to the order of items_ordered_by_frequency
        for ( const auto& pair : items_ordered_by_frequency ) {
            const Item& item = pair.first;
            //cout<<item<<" ";
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
        if(curr_fpnode)
            curr_fpnode->tid_list.insert(tids[i]);
        //cout<<endl;
        i++;
    }

    // auto cur_node = root;
    // queue<shared_ptr<FPNode>> q;
    // map<shared_ptr<FPNode>,int> lev;
    // q.push(cur_node);
    // lev[cur_node]=0;
    // int prevLevel=-1;
    // while(q.size()){
    //     auto t = q.front();
    //     q.pop();
    //     if(lev[t]!=prevLevel){
    //         prevLevel=lev[t];
    //         cout<<endl;
    //     }
    //     if(t == nullptr)
    //         cout<<0<<" ";
    //     else{
    //         cout<<t->item<<" ( ";
    //         for(auto it=t->tid_list.begin();it!=t->tid_list.end();it++){
    //             cout<<*it<<" ";
    //         }
    //         cout<<") ";
    //     }
    //     for(auto v: t->children){
    //         if(lev.find(v)==lev.end()){
    //             q.push(v);
    //             lev[v] = lev[t]+1;
    //         }
    //     }
    // }
    // cout<<endl;
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

std::set<onlyPattern> fptree_growth(const FPTree& fptree)
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
                    conditional_fptree_tids.pb(*it);
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
            const FPTree conditional_fptree( conditional_fptree_tids, conditional_fptree_transactions, fptree.minimum_support_threshold, fptree.maximum_periodicity);
            // call recursively fptree_growth on the conditional fptree (empty fptree: no patterns)

            std::set<onlyPattern> conditional_patterns = fptree_growth( conditional_fptree);

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
    // test.pb({"1","2","5"});
    // test.pb({"2","4"});
    // test.pb({"2","3"});
    // test.pb({"1","2","4"});
    // test.pb({"1","3"});
    // test.pb({"2","3"});
    // test.pb({"1","3"});
    // test.pb({"1","2","3","5"});
    // test.pb({"1","2","3"});

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
    TDB_SIZE = len;
    // cout<<transactions.size()<<endl;
    // for(int i=0;i<transactions.size();i++){
    //     for(int j=0;j<transactions[i].size();j++){
    //         cout<<transactions[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }

    const FPTree fptree{ tids, transactions, min_sup, max_per };

    const std::set<onlyPattern> patterns = fptree_growth( fptree );

    /*auto it=patterns.begin();
    while(it!=patterns.end()){
        set<Item> s=(*it).first;
        auto itr=s.begin();
        while(itr!=s.end()){
            cout<<(*itr)<<" ";
            itr++;
        }
        cout<<"Tids: ";
        for(auto itr=it->second.begin();itr!=it->second.end();itr++){
            cout<<*itr<<" ";
        }
        cout<<endl;
        it++;
    }*/
    cout << patterns.size();
}