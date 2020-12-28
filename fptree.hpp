#ifndef FPTREE_HPP
#define FPTREE_HPP

#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <utility>

using namespace std;

using Item = std::string;
using Transaction = std::vector<Item>;
using TransformedPrefixPath = std::pair<std::vector<Item>, std::set<int>>;
using Pattern = std::pair<std::set<Item>, set<int>>;


struct FPNode {
    const Item item;
    int frequency;
    std::shared_ptr<FPNode> node_link;
    std::weak_ptr<FPNode> parent;
    std::vector<std::shared_ptr<FPNode>> children;
    std::set<int> tid_list;

    FPNode(const Item&, const std::shared_ptr<FPNode>&);
};

struct FPTree {
    std::shared_ptr<FPNode> root;
    std::map<Item, std::shared_ptr<FPNode>> header_table;
    const int minimum_support_threshold;
    const int maximum_periodicity;
    vector<pair<Item,int>> items_with_frequency;

    FPTree(const std::vector<int>, const std::vector<Transaction>&, const int, const int);

    bool empty() const;
};


std::set<Pattern> fptree_growth(const FPTree&);


#endif  // FPTREE_HPP