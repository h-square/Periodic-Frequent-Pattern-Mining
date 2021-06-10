#ifndef PFTREE_HPP
#define PFTREE_HPP

#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <utility>

using namespace std;

using Item = string;
using Transaction = vector<Item>;
using TransformedPrefixPath = pair<vector<Item>, set<int>>;
using Pattern = set<Item>;

struct PFNode {
    const Item item;
    int frequency;
    shared_ptr<PFNode> node_link;
    weak_ptr<PFNode> parent;
    vector<shared_ptr<PFNode>> children;
    set<int> tid_list;

    PFNode(const Item&, const shared_ptr<PFNode>&);
};

struct PFTree {
    shared_ptr<PFNode> root;
    map<Item, shared_ptr<PFNode>> header_table;
    const int minimum_support_threshold;
    const int maximum_periodicity;
    vector<pair<Item,int>> pfitems_by_frequency;

    // constructor for serial PF-tree 
    PFTree(const vector<int>, const vector<Transaction>&, const int, const int);

    // constructor for parallel PF-tree 
    PFTree(const vector<int>, const vector<Transaction>&, const int, const int, const int);
    
    PFTree(const int, const int);
    bool empty() const;
};

#endif  // PFTREE_HPP